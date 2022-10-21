using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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
        }

        public string FilePath { get; set; }
        public string[] Headers { get; set; }
        public List<string> PSMdata { get; }

        public int DeltaMassCol { get; set; }
        public int AssignedModCol { get; set; }
        public int PeptideCol { get; set; }
        public int SpectrumCol { get; set; }
        public int ObservedModCol { get; set; }
        public int PrecursorMZCol { get; set; }

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
        public string editPSMLine(string localizerOutput, List<string> existingPSMline, bool overwritePrevious)
        {
            int insertIndex = GetColumnIndex("Observed Modifications") + 1;

            if (!overwritePrevious)
            {
                existingPSMline.InsertRange(insertIndex, localizerOutput.Split("\t"));
            }
            else
            {
                // previous scan data exists - overwrite it
                string[] localizerOutputEntries = localizerOutput.Split("\t");
                for (int i = 0; i < localizerOutputEntries.Length; i++)     // assumes same length localizer output as before
                {
                    existingPSMline[insertIndex + i] = localizerOutputEntries[i];
                }
            }
            return String.Join("\t", existingPSMline);
        }
    }
}
