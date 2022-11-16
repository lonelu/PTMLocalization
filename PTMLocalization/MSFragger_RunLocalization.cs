using System;
using System.Collections.Generic;
using System.Linq;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using MzLibUtil;
using EngineLayer.GlycoSearch;
using IO.Mgf;
using IO.MzML;
using MassSpectrometry;
using EngineLayer;
using System.IO;

namespace PTMLocalization
{
    public class MSFragger_RunLocalization
    {
        private Tolerance ProductMassTolerance;
        private Tolerance PrecursorMassTolerance;
        private readonly string psmFile;
        private readonly string scanpairFile;
        private readonly string rawfilesDirectory;
        private readonly string o_glycan_database;
        private int maxOGlycansPerPeptide;
        private int[] isotopes;
        
        public static readonly double AveragineIsotopeMass = 1.00235;
        public static readonly string OPAIR_HEADERS = "OPair Score\tNumber of Glycans\tTotal Glycan Composition\tGlycan Site Composition(s)\tConfidence Level\tSite Probabilities\t138/144 Ratio\tPaired Scan Num";
        //public static readonly string OPAIR_HEADERS = "OPair Score\tNumber of Glycans\tTotal Glycan Composition\tGlycan Site Composition(s)\tConfidence Level\tSite Probabilities";
        public static readonly int OUTPUT_LENGTH = OPAIR_HEADERS.Split("\t").Length;
        public static readonly string EMPTY_OUTPUT = String.Join("\t", new string[OUTPUT_LENGTH]);
        public static readonly string EMPTY_OUTPUT_WITH_PAIRED_SCAN = String.Join("\t", new string[OUTPUT_LENGTH - 1]);

        private static readonly string CAL_INPUT_EXTENSION = "_calibrated.mzML";
        private static readonly string CAL_INPUT_EXTENSION_FALLBACK = "_calibrated.MGF";

        private static readonly double OXO138 = 138.05495;
        private static readonly double OXO144 = 144.06552;

        public MSFragger_RunLocalization(string _psmFile, string _scanpairFile, string _rawfilesDirectory, string _o_glycan_database_path, int _maxOGlycansPerPeptide, Tolerance _PrecursorMassTolerance, Tolerance _ProductMassTolerance, int[] _isotopes)
        {
            psmFile = _psmFile;
            scanpairFile = _scanpairFile;
            rawfilesDirectory = _rawfilesDirectory;
            o_glycan_database = _o_glycan_database_path;
            PrecursorMassTolerance = _PrecursorMassTolerance;
            ProductMassTolerance = _ProductMassTolerance;
            maxOGlycansPerPeptide = _maxOGlycansPerPeptide;
            isotopes = _isotopes;

            Setup();
        }

        /**
         * Load glycan database and initialize glycan boxes given database and max number of O-glycans per peptide provided for the search.
         */
        public void Setup()
        {
            GlobalVariables.SetUpGlobalVariables();
            GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(o_glycan_database, true, true).ToArray();
            GlycanBox.GlobalOGlycanMods = GlycanBox.BuildGlobalOGlycanMods(GlycanBox.GlobalOGlycans).ToArray();
            GlycanBox.OGlycanBoxes = GlycanBox.BuildGlycanBoxes(maxOGlycansPerPeptide, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods).OrderBy(p => p.Mass).ToArray();
        }


        /**
         * Main method to localize from MSFragger outputs. Loads the MSFragger PSM and scan-pair tables, parses PSM information
         * to determine which scans to localize and the input peptides/glycan masses for those scans, runs localization,
         * and writes results back to the PSM table. 
         */
        public void Localize()
        {
            Dictionary<int, int> scanPairs = new();
            if (scanpairFile != null)
            {
                // single scanPair file passed directly - use. Otherwise, will read below during PSM reading
                scanPairs = MSFragger_PSMTable.ParseScanPairTable(scanpairFile);
            }

            // read PSM table to prepare to pass scans to localizer
            MSFragger_PSMTable PSMtable = new(psmFile);
            string currentRawfile = "";
            string currentRawfileBase = "";
            MsDataFile currentMsDataFile = null;
            Dictionary<int, MsDataScan> msScans = new();
            List<string> output = new();
            bool overwritePrevious = false;
            if (PSMtable.Headers.Contains("OPair Score"))
            {
                // overwriting previous OP results, don't add new columns
                output.Add(String.Join("\t", PSMtable.Headers));
                overwritePrevious = true;
            } 
            else
            {
                output.Add(PSMtable.editPSMLine(OPAIR_HEADERS, PSMtable.Headers.ToList(), overwritePrevious));
            }

            Dictionary<string, bool> mgfNotFoundWarnings = new ();
            Dictionary<string, bool> mzmlNotFoundWarnings = new();
            // timers and etc
            System.Diagnostics.Stopwatch timer = new System.Diagnostics.Stopwatch();
            System.Diagnostics.Stopwatch totalTimer = new System.Diagnostics.Stopwatch();
            timer.Start();
            totalTimer.Start();
            int psmCount = 0;

            foreach (string PSMline in PSMtable.PSMdata)
            {
                if (timer.ElapsedMilliseconds > 10000)
                {
                    Console.Write("\r\tprogress: {0}/{1} PSMs processed in {2:0.0} s", psmCount, PSMtable.PSMdata.Count, totalTimer.ElapsedMilliseconds * 0.001);
                    timer.Restart();
                }
                psmCount++;

                // Parse PSM information
                string[] lineSplits = PSMline.Split("\t");
                string spectrumString = lineSplits[PSMtable.SpectrumCol];
                string peptide = lineSplits[PSMtable.PeptideCol];
                double deltaMass = Double.Parse(lineSplits[PSMtable.DeltaMassCol]);
                string assignedMods = lineSplits[PSMtable.AssignedModCol];
                double precursorMZ = Double.Parse(lineSplits[PSMtable.PrecursorMZCol]);
                if (deltaMass < 3.5 && deltaMass > -1.5)
                {
                    // unmodified peptide - no localization
                    output.Add(PSMtable.editPSMLine(EMPTY_OUTPUT, lineSplits.ToList(), overwritePrevious));
                    continue;
                }

                // Load scan (and rawfile if necessary)
                int scanNum = MSFragger_PSMTable.GetScanNum(spectrumString);
                string rawfileBase = MSFragger_PSMTable.GetRawFile(spectrumString);
                string rawfileName = rawfileBase + CAL_INPUT_EXTENSION;
                string scanpairName = rawfileBase + ".pairs";
                int precursorCharge = MSFragger_PSMTable.GetScanCharge(spectrumString);

                if (!rawfileBase.Equals(currentRawfileBase))
                {
                    // new rawfile: load scan data
                    string spectraFile = Path.Combine(rawfilesDirectory, rawfileName);
                    FilteringParams filter = new FilteringParams();
                    if (File.Exists(spectraFile))
                    {
                        try
                        {
                            currentMsDataFile = Mzml.LoadAllStaticData(spectraFile, filter, searchForCorrectMs1PrecursorScan:false);
                        } 
                        catch
                        {
                            currentMsDataFile = null;
                        }
                    } 
                    else if (File.Exists(Path.Combine(rawfilesDirectory, rawfileBase + CAL_INPUT_EXTENSION_FALLBACK)))
                    {
                        // support legacy "_calibrated.MGF" if present and "_calibrated.mzML" is not
                        spectraFile = Path.Combine(rawfilesDirectory, rawfileBase + CAL_INPUT_EXTENSION_FALLBACK);
                        filter = new FilteringParams();
                        Console.WriteLine("Calibrated mzML not found, using calibrated MGF for file {0}", rawfileBase);
                        try { 
                            currentMsDataFile = Mgf.LoadAllStaticData(spectraFile, filter);
                        }
                        catch
                        {
                            currentMsDataFile = null;
                        }
                    }
                    
                    if (currentMsDataFile is null)
                    {
                        // no calibrated file (mzML/MGF) found (or loading failed), try falling back to raw mzML
                        rawfileName = rawfileBase + ".mzML";
                        spectraFile = Path.Combine(rawfilesDirectory, rawfileName);
                        if (File.Exists(spectraFile))
                        {
                            currentMsDataFile = Mzml.LoadAllStaticData(spectraFile, filter, searchForCorrectMs1PrecursorScan:false);
                            // warn user that MGF wasn't found for this raw file if we haven't already done so
                            if (!mgfNotFoundWarnings.ContainsKey(rawfileBase))
                            {
                                Console.WriteLine("Warning: calibrated file not found or not loaded for raw file {0}, uncalibrated mzML will be used", rawfileBase);
                                mgfNotFoundWarnings.Add(rawfileBase, true);
                            }
                        } 
                        else 
                        {
                            if (!mzmlNotFoundWarnings.ContainsKey(rawfileBase))
                            {
                                Console.WriteLine("Warning: no MGF or mzML found for file {0}, PSMs from this file will NOT be localized!", rawfileBase);
                                mzmlNotFoundWarnings.Add(rawfileBase, true);
                            }
                            output.Add(PSMtable.editPSMLine(EMPTY_OUTPUT, lineSplits.ToList(), overwritePrevious));
                            continue;
                        }
                    }
                    List<MsDataScan> allScans = currentMsDataFile.GetAllScansList();
                    // save scans by scan number, since MGF file from MSFragger does not guarantee all scans will be saved
                    msScans = new();
                    foreach (MsDataScan currentScan in allScans)
                    {
                        msScans[currentScan.OneBasedScanNumber] = currentScan;
                    }

                    // load scan pair file if not done already
                    if (scanpairFile == null)
                    {
                        string pairsFile = Path.Combine(rawfilesDirectory, scanpairName);
                        scanPairs = MSFragger_PSMTable.ParseScanPairTable(pairsFile);
                    }
                    currentRawfile = rawfileName;
                    currentRawfileBase = rawfileBase;
                }

                // retrieve the right child scan
                int childScanNum;
                if (scanPairs.ContainsKey(scanNum))
                {
                    childScanNum = scanPairs[scanNum];
                }
                else
                {
                    // no paired child scan found - do not attempt localization
                    output.Add(PSMtable.editPSMLine("No paired scan" + EMPTY_OUTPUT, lineSplits.ToList(), overwritePrevious));
                    continue;
                }

                // retrieve the child scan
                MsDataScan ms2Scan;
                try
                {
                    ms2Scan = msScans[childScanNum];
                }
                catch (KeyNotFoundException)
                {
                    Console.Out.WriteLine(String.Format("Error: MS2 scan {0} not found in the MGF file. No localization performed", childScanNum));
                    output.Add(PSMtable.editPSMLine(EMPTY_OUTPUT, lineSplits.ToList(), overwritePrevious));
                    continue;
                }

                // finalize spectrum for search
                IsotopicEnvelope[] neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2Scan, 4, 3);
                var scan = new Ms2ScanWithSpecificMass(ms2Scan, precursorMZ, precursorCharge, currentRawfile, 4, 3, neutralExperimentalFragments);

                // initialize peptide with all non-glyco mods
                PeptideWithSetModifications peptideWithMods = getPeptideWithMSFraggerMods(assignedMods, peptide);

                // get oxo 138/144 ratio from the HCD scan
                MsDataScan hcdDataScan = msScans[scanNum];
                var hcdScan = new Ms2ScanWithSpecificMass(hcdDataScan, precursorMZ, precursorCharge, currentRawfile, 4, 3, neutralExperimentalFragments);
                double ratio = hcdScan.computeOxoRatio(OXO138, OXO144, ProductMassTolerance.Value);

                // finally, run localizer
                string localizerOutput = LocalizeOGlyc(scan, peptideWithMods, deltaMass, ratio);

                // write info back to PSM table
                output.Add(PSMtable.editPSMLine(localizerOutput, lineSplits.ToList(), overwritePrevious));
            }
            // write output to PSM table
            File.WriteAllLines(psmFile, output.ToArray());
            totalTimer.Stop();
            Console.WriteLine("\nFinished localization in {0:0.0} s", totalTimer.ElapsedMilliseconds * 0.001);
        }


        /**
         * Single scan O-glycan localization, given a peptide and glycan mass from MSFragger search
         */
        public string LocalizeOGlyc(Ms2ScanWithSpecificMass ms2Scan, PeptideWithSetModifications peptide, double glycanDeltaMass, double oxoRatio)
        {
            // generate peptide fragments
            List<Product> products = new List<Product>();
            peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);

            // generate glycan modification sites
            List<int> n_modPos = new List<int>();
            var o_modPos = GlycoPeptides.GetPossibleModSites(peptide, new string[] { "S", "T" }).OrderBy(p => p).ToList();

            int[] _modPos;
            string[] modMotifs;
            GlycoPeptides.GetModPosMotif(GlycoType.OGlycoPep, n_modPos.ToArray(), o_modPos.ToArray(), out _modPos, out modMotifs);

            // Get possible O-glycan boxes from delta mass
            List<LocalizationGraph> graphs = new List<LocalizationGraph>();
            foreach (int isotope in isotopes)
            {
                // include peptide mass for PPM calculation between precursor and glycan delta mass
                double currentDeltaMass = glycanDeltaMass - (isotope * AveragineIsotopeMass) + peptide.MonoisotopicMass;
                double maxMassErrorDa = currentDeltaMass * 0.000001 * PrecursorMassTolerance.Value;
                var possibleGlycanMassLow = currentDeltaMass - maxMassErrorDa - peptide.MonoisotopicMass;   // subtract peptide mass to get final glycan-only mass
                var possibleGlycanMassHigh = currentDeltaMass + maxMassErrorDa - peptide.MonoisotopicMass;
                int iDLow = GlycoPeptides.BinarySearchGetIndex(GlycanBox.OGlycanBoxes.Select(p => p.Mass).ToArray(), possibleGlycanMassLow);

                while (iDLow < GlycanBox.OGlycanBoxes.Length && GlycanBox.OGlycanBoxes[iDLow].Mass <= possibleGlycanMassHigh)
                {
                    LocalizationGraph localizationGraph = new LocalizationGraph(_modPos, modMotifs, GlycanBox.OGlycanBoxes[iDLow], GlycanBox.OGlycanBoxes[iDLow].ChildGlycanBoxes, iDLow);
                    LocalizationGraph.LocalizeMod(localizationGraph, ms2Scan, ProductMassTolerance,
                        products.Where(v => v.ProductType == ProductType.c || v.ProductType == ProductType.zDot).ToList(),
                        GlycoPeptides.GetLocalFragmentGlycan, GlycoPeptides.GetUnlocalFragmentGlycan);
                    graphs.Add(localizationGraph);
                    iDLow++;
                }
            }

            // save results
            if (graphs.Count == 0)
            {
                // no matching glycan boxes found to this delta mass, no localization can be done! 
                return "No match to glycan delta mass" + EMPTY_OUTPUT_WITH_PAIRED_SCAN + string.Format("\t{0}", ms2Scan.TheScan.OneBasedScanNumber);      // keep the psm line the same length as all other outputs
            }
            var best_graph = graphs.OrderByDescending(p => p.TotalScore).First();
            var gsm = new GlycoSpectralMatch(new List<LocalizationGraph> { best_graph }, GlycoType.OGlycoPep);
            gsm.ScanInfo_p = ms2Scan.TheScan.MassSpectrum.Size * ProductMassTolerance.GetRange(1000).Width / ms2Scan.TheScan.MassSpectrum.Range.Width;
            gsm.Thero_n = products.Count();
            gsm.oxoRatio = oxoRatio;
            GlycoSite.GlycoLocalizationCalculation(gsm, gsm.GlycanType, DissociationType.HCD, DissociationType.EThcD);

            return gsm.WriteLine(null) + string.Format("\t{0}", ms2Scan.TheScan.OneBasedScanNumber);
        }

        /**
         * Set up the peptide container with MSFragger assigned modifications included. 
         */
        private PeptideWithSetModifications getPeptideWithMSFraggerMods(string assignedMods, string peptide)
        {
            PeptideWithSetModifications peptideWithMods;
            if (assignedMods.Length > 0)
            {
                // MSFragger mods are present: generate a peptide with them included
                string[] assignedModSplits = assignedMods.Split(",");
                Dictionary<string, Modification> modDefinitionDict = new();
                Dictionary<int, string> modsForPeptideSeqDict = new();
                string modType = "MSFragger";

                foreach (string split in assignedModSplits)
                {
                    string[] massSplits = split.Split("(");
                    string residueType;
                    int position;
                    // handle terminal mods as well
                    if (split.Contains("N-term"))
                    {
                        // terminal mod, position and residue type are different (encoding as the terminal residue, not the terminus, because I don't know the syntax for the terminus)
                        residueType = peptide[0].ToString();
                        position = 0;
                    } 
                    else if (split.Contains("C-term"))
                    {
                        // terminal mod, position and residue type are different (encoding as the terminal residue, not the terminus, because I don't know the syntax for the terminus)
                        residueType = peptide[peptide.Length - 1].ToString();
                        position = peptide.Length - 1;
                    }
                    else
                    {
                        residueType = massSplits[0][massSplits[0].Length - 1].ToString();
                        position = Int32.Parse(massSplits[0].Substring(0, massSplits[0].Length - 1)) - 1;   // MSFragger mod positions are 1-indexed, convert to 0-index
                    }
                    
                    double mass = double.Parse(massSplits[1].Replace(")", ""));
                    string modName = "(" + massSplits[1];
                    // Generate a Modification of type "MSFragger" with ID "(mass)" 
                    ModificationMotif.TryGetMotif(residueType, out ModificationMotif motif2);
                    Modification mod = new Modification(_originalId: modName, _modificationType: modType, _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: mass);
                    modDefinitionDict.TryAdd(mod.IdWithMotif, mod);     // ignore if this mod definition is already present
                    modsForPeptideSeqDict[position] = string.Format("[{0}:{1}]", modType, mod.IdWithMotif);
                }
                peptideWithMods = MSFragger_PSMTable.GetMSFraggerPeptide(peptide, modsForPeptideSeqDict, modDefinitionDict);
            }
            else
            {
                // unmodified peptide
                peptideWithMods = new PeptideWithSetModifications(peptide, new Dictionary<string, Modification>());
            }
            return peptideWithMods;
        }
    }
}
