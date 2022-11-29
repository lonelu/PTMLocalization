using EngineLayer;
using EngineLayer.GlycoSearch;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace PTMLocalization
{
    /**
     * Whole PSM table container with support for getting specific column indices/etc
     */
    public class MSFragger_PSMTable
    {
        public MSFragger_PSMTable(string FilePath)
        {
            this.FilePath = FilePath;
            string[] allLines = File.ReadAllLines(FilePath);
            Headers = allLines[0].Split("\t");
            PSMdata = new List<string>(allLines);
            PSMdata.RemoveAt(0);    // remove header line

            // init headers
            DeltaMassCol = GetColumnIndex("Delta Mass");
            AssignedModCol = GetColumnIndex("Assigned Modifications");
            PeptideCol = GetColumnIndex("Peptide");
            ModifiedPeptideCol = GetColumnIndex("Modified Peptide");
            SpectrumCol = GetColumnIndex("Spectrum");
            ObservedModCol = GetColumnIndex("Observed Modifications");
            int calMZCol = GetColumnIndex("Calibrated Observed M/Z");
            if (calMZCol != -1)
            {
                PrecursorMZCol = calMZCol;
            }
            else
            {
                PrecursorMZCol = GetColumnIndex("Observed M/Z");
            }
            CalcMZCol = GetColumnIndex("Calculated M/Z");
            CalcPeptideMassCol = GetColumnIndex("Calculated Peptide Mass");
            ChargeCol = GetColumnIndex("Charge");
        }

        public string FilePath { get; set; }
        public string[] Headers { get; set; }
        public List<string> PSMdata { get; }

        public int DeltaMassCol { get; set; }
        public int AssignedModCol { get; set; }
        public int PeptideCol { get; set; }
        public int ModifiedPeptideCol { get; set; }
        public int SpectrumCol { get; set; }
        public int ObservedModCol { get; set; }
        public int PrecursorMZCol { get; set; }
        public int CalcMZCol { get; set; }
        public int CalcPeptideMassCol { get; set; }
        public int ChargeCol { get; set; }

        public int GetColumnIndex(string columnName)
        {
            for (int i = 0; i < Headers.Length; i++)
            {
                if (Headers[i].Equals(columnName))
                {
                    return i;
                }
            }
            return -1;
        }

        public static string GetRawFile(string spectrum)
        {
            string[] splits = spectrum.Split(".");
            return splits[0];
        }

        public static int GetScanNum(string spectrum)
        {
            string[] splits = spectrum.Split(".");
            return Int32.Parse(splits[1]);
        }
        public static int GetScanCharge(string spectrum)
        {
            string[] splits = spectrum.Split(".");
            return Int32.Parse(splits[splits.Length - 1]);
        }

        /**
         * Parse the scan-pairing information table from FragPipe.
         * Format: scan1 \t scan2 \n
         */
        public static Dictionary<int, int> ParseScanPairTable(string filepath)
        {
            Dictionary<int, int> scanPairs = new();
            string[] allLines = File.ReadAllLines(filepath);
            foreach (string line in allLines)
            {
                string[] splits = line.Split('\t');
                try
                {
                    scanPairs.Add(Int32.Parse(splits[0]), Int32.Parse(splits[1]));
                } catch (ArgumentException)
                {
                    // multiple paired scans for this parent scan. Only take the first one
                    continue;
                }
            }
            return scanPairs;
        }

        /**
         * Add the known modifications to peptide sequence to translate to a PeptideWithSetModifications object
         */
        public static PeptideWithSetModifications GetMSFraggerPeptide(string baseSequence, Dictionary<int, string> modPositions, Dictionary<string, Proteomics.Modification> modDefinitions)
        {
            StringBuilder sb = new();
            for (int i = 0; i < baseSequence.Length; i++)
            {
                // add base sequence residue
                sb.Append(baseSequence[i]);
                if (modPositions.ContainsKey(i))
                {
                    // there is a mod at this position, place it as well
                    sb.Append(modPositions[i]);
                }
            }
            string peptideWithMods = sb.ToString();
            PeptideWithSetModifications finalPeptide = new PeptideWithSetModifications(peptideWithMods, modDefinitions);
            return finalPeptide;
        }

        /**
         * Insert the provided OPair output into the PSM line at the set insertion point, OR, overwrite the existing OPair
         * output if present. 
         */
        public string editPSMLine(GlycoSpectralMatch gsm, List<string> existingPSMline, Glycan[] globalGlycans, bool overwritePrevious, bool writeToAssignedMods)
        {
            int insertIndex = this.ObservedModCol + 1;

            if (!overwritePrevious)
            {
                existingPSMline.InsertRange(insertIndex, gsm.localizerOutput.Split("\t"));
            }
            else
            {
                // previous scan data exists - overwrite it [NOTE: debugging mode only - not for use with writeToAssignedMods]
                string[] localizerOutputEntries = gsm.localizerOutput.Split("\t");
                for (int i = 0; i < localizerOutputEntries.Length; i++)     // assumes same length localizer output as before
                {
                    existingPSMline[insertIndex + i] = localizerOutputEntries[i];
                }
            }

            if (gsm.LocalizedGlycan == null)
            {
                return String.Join("\t", existingPSMline);
            }

            // also edit Assigned Modifications (and modified peptide) if requested
            if (writeToAssignedMods)
            {               
                string peptide = existingPSMline[PeptideCol];
                List<string> newAssignedMods = new List<string>();
                // read existing assigned mods
                string prevAssignedModsStr = existingPSMline[AssignedModCol];
                if (prevAssignedModsStr.Length > 0) {
                    string[] prevAssignedMods = prevAssignedModsStr.Split(", ");
                    for (int i = 0; i < prevAssignedMods.Length; i++)
                    {
                        newAssignedMods.Add(prevAssignedMods[i]);
                    }
                }

                // add new assigned modification for each O-Pair mod
                double opairGlycanMass = 0;
                bool hasUnlocalizedGlycs = false;
                byte[] unlocalizedGlycanComp = gsm.getTotalKind();
                foreach (var localizedGlyc in gsm.LocalizedGlycan)
                {
                    if (localizedGlyc.IsLocalized)
                    {
                        unlocalizedGlycanComp = Glycan.subtractKind(unlocalizedGlycanComp, globalGlycans[localizedGlyc.GlycanID].Kind);     // keep track of what glycan(s) are left unlocalized

                        double mass = EngineLayer.GlycanBox.GlobalOGlycans[localizedGlyc.GlycanID].Mass * 0.00001;    // mass is saved as int * 10e5 in GlycanBox
                        opairGlycanMass += mass;
                        var peptide_site = localizedGlyc.ModSite - 1;
                        string aa = peptide.Substring(peptide_site - 1, 1);
                        var comp = EngineLayer.GlycanBox.GlobalOGlycans[localizedGlyc.GlycanID].Composition;
                        newAssignedMods.Add(string.Format("{0}{1}({2:0.0000})", peptide_site, aa, mass));

                        // edit modified peptide col
                        string originalModPeptide = existingPSMline[this.ModifiedPeptideCol];
                        if (originalModPeptide.Length == 0)
                        {
                            originalModPeptide = peptide;   // no modified peptide is written if no previous mods. Start with base sequence
                        }
                        int residueIndex = -1;
                        for (int i = 0; i < originalModPeptide.Length; i++)
                        {
                            if (originalModPeptide[i] >= 'A' && originalModPeptide[i] <= 'Z')
                            {
                                // this is an actual residue, increment residue counter
                                residueIndex++;
                            }
                            if (residueIndex == peptide_site)
                            {
                                // Found the glycan location - insert new mod and stop looking
                                double aaMass = Proteomics.AminoAcidPolymer.Residue.GetResidue(aa).MonoisotopicMass;
                                string modToInsert = String.Format("[{0:0}]", mass + aaMass);
                                string newModPep = originalModPeptide.Insert(i, modToInsert);
                                existingPSMline[this.ModifiedPeptideCol] = newModPep;
                                break;
                            }
                        }
                    }
                    else
                    {
                        // level 3 or otherwise unlocalized glycan. NOTE: there may be multiple entries for each possible equivalent site, so do not read the mass here as it may not be correct
                        hasUnlocalizedGlycs = true;
                    }
                }
                // if fewer glycosites than glycans (e.g., N-glycan/other glycan present), may have remaining localized glycans: check for them
                foreach (byte b in unlocalizedGlycanComp)
                {
                    if (b > 0)
                    {
                        hasUnlocalizedGlycs = true;
                        break;
                    }
                }
                // handle unlocalized glycan(s) writing to assigned mods
                if (hasUnlocalizedGlycs)
                {
                    newAssignedMods.Sort();
                    // at least one unlocalized glycan present. Determine total unlocalized mass
                    double unlocMass = Glycan.GetMass(unlocalizedGlycanComp) * 0.00001;
                    opairGlycanMass += unlocMass;
                    // add total unlocalized mass as an assigned mod
                    newAssignedMods.Insert(0, string.Format("{0:0.0000}", unlocMass));      // put the unlocalized mod first
                } else
                {
                    newAssignedMods.Sort();
                }

                // write final mods back to assigned mods column
                ModComparer mc = new ModComparer();
                newAssignedMods.Sort(mc);
                string finalAssignedMods = string.Join(",", newAssignedMods);
                existingPSMline[this.AssignedModCol] = finalAssignedMods;

                // edit delta mass, calculated peptide mass and m/z, and 
                double originalDeltaMass = Double.Parse(existingPSMline[this.DeltaMassCol]);
                existingPSMline[this.DeltaMassCol] = String.Format("{0:0.0000}", originalDeltaMass - opairGlycanMass);

                double prevCalcMZ = Double.Parse(existingPSMline[this.CalcMZCol]);
                double newCalcMZ = prevCalcMZ + (opairGlycanMass / double.Parse(existingPSMline[this.ChargeCol]));      // previous m/z already has proton mass(es) handled, so just add glycan mass/charge
                existingPSMline[this.CalcMZCol] = String.Format("{0:0.0000}", newCalcMZ);

                double prevCalcPeptideMass = Double.Parse(existingPSMline[this.CalcPeptideMassCol]);
                existingPSMline[this.CalcPeptideMassCol] = String.Format("{0:0.0000}", prevCalcPeptideMass + opairGlycanMass);
                
            }

            return String.Join("\t", existingPSMline);
        }

        // check if this PSM table has previously had OPair results written by looking for O-Pair column headers
        public bool hasPreviousOPairResults()
        {
            return GetColumnIndex("O-Pair Score") != -1;     // -1 means not found, no previous O-Pair results
        }

        
    }

    public class ModComparer : IComparer<string>
    {
        public ModComparer() { }
        public int Compare(string mod1, string mod2)
        {
            int pos1 = int.Parse(Regex.Match(mod1, "[0-9]+").Value);
            int pos2 = int.Parse(Regex.Match(mod2, "[0-9]+").Value);

            if (pos1 < pos2)
            {
                return -1;
            }
            else if (pos1 == pos2)
            {
                return 0;
            }
            else
            {
                return 1;
            }
        }
    }
}
