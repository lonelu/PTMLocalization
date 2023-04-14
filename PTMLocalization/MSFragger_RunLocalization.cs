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
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Threading;

namespace PTMLocalization
{
    public class MSFragger_RunLocalization
    {
        private Tolerance ProductMassTolerance;
        private Tolerance PrecursorMassTolerance;
        private readonly string psmFile;
        private readonly string scanpairFile;
        private string rawfilesDirectory;
        private readonly string lcmsFileListPath;
        private readonly string o_glycan_database;
        private int maxOGlycansPerPeptide;
        private int[] isotopes;
        private bool filterOxonium;

        private Dictionary<string, string> lcmsPaths;

        public static readonly double AveragineIsotopeMass = 1.00235;
        public static readonly string OPAIR_HEADERS = "O-Pair Score\tNumber of Glycans\tTotal Glycan Composition\tGlycan Site Composition(s)\tConfidence Level\tSite Probabilities\t138/144 Ratio\tHas N-Glyc Sequon\tPaired Scan Num";
        //public static readonly string OPAIR_HEADERS = "O-Pair Score\tNumber of Glycans\tTotal Glycan Composition\tGlycan Site Composition(s)\tConfidence Level\tSite Probabilities";
        public static readonly int OUTPUT_LENGTH = OPAIR_HEADERS.Split("\t").Length;
        public static readonly string EMPTY_OUTPUT = String.Join("\t", new string[OUTPUT_LENGTH]);
        public static readonly string EMPTY_OUTPUT_WITH_PAIRED_SCAN = String.Join("\t", new string[OUTPUT_LENGTH - 1]);

        private static readonly string CAL_INPUT_EXTENSION = "_calibrated.mzML";
        private static readonly string CAL_INPUT_EXTENSION_FALLBACK = "_uncalibrated.mzML";
        private static readonly string CAL_INPUT_EXTENSION_FALLBACK2 = "_calibrated.MGF";

        private static readonly double OXO138 = 138.05495;
        private static readonly double OXO144 = 144.06552;
        private static readonly Regex NglycMotifRegex = new Regex("N[^P][ST]", RegexOptions.Compiled);

        public MSFragger_RunLocalization(string _psmFile, string _scanpairFile, string _rawfilesDirectory, string lcmsFileList, string _o_glycan_database_path, int _maxOGlycansPerPeptide, Tolerance _PrecursorMassTolerance, Tolerance _ProductMassTolerance, int[] _isotopes, bool _filterOxonium)
        {
            psmFile = _psmFile;
            scanpairFile = _scanpairFile;
            rawfilesDirectory = _rawfilesDirectory;
            lcmsFileListPath = lcmsFileList;
            o_glycan_database = _o_glycan_database_path;
            PrecursorMassTolerance = _PrecursorMassTolerance;
            ProductMassTolerance = _ProductMassTolerance;
            maxOGlycansPerPeptide = _maxOGlycansPerPeptide;
            isotopes = _isotopes;
            filterOxonium = _filterOxonium;

            Setup();
        }

        /**
         * Load glycan database and initialize glycan boxes given database and max number of O-glycans per peptide provided for the search.
         */
        public void Setup()
        {
            System.Diagnostics.Stopwatch timer = new();
            Console.Write("Loading Glycan Database...");
            timer.Start();
            GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(o_glycan_database, true, true).ToArray();
            GlycanBox.GlobalOGlycanMods = GlycanBox.BuildGlobalOGlycanMods(GlycanBox.GlobalOGlycans).ToArray();
            GlycanBox.OGlycanBoxes = GlycanBox.BuildGlycanBoxes(maxOGlycansPerPeptide, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods).OrderBy(p => p.Mass).ToArray();

            timer.Stop();
            Console.WriteLine(String.Format("done in {0:0.0}s\n\tLoaded {1} glycans, {2} glycan boxes", timer.ElapsedMilliseconds * 0.001, GlycanBox.GlobalOGlycans.Length, GlycanBox.OGlycanBoxes.Length));
        }

        private MsDataFile CheckAndLoadData(string rawfileBase, string rawfileName, Dictionary<string, bool> mzmlNotFoundWarnings, bool usingLcmsFilePath)
        {
            // new rawfile: load scan data
            string spectraFile = null;
            bool fileFound = false;
            if (lcmsPaths.Count > 0)
            {
                // read from FragPipe file list
                fileFound = lcmsPaths.TryGetValue(rawfileBase, out spectraFile);
            } 
            if (!fileFound && rawfilesDirectory != null)
            {
                // fall back to read from rawfile directory if provided
                spectraFile = Path.Combine(rawfilesDirectory, rawfileName);
            }

            if (spectraFile == null)
            {
                Console.WriteLine(String.Format("Error: could not find spectra file from FragPipe list for file {0}", rawfileBase));
                return null;
            } 

            // init directory path if reading from LCMS file list rather than from raw directory
            if (usingLcmsFilePath)
            {
                rawfilesDirectory = Path.GetDirectoryName(spectraFile);
            }

            FilteringParams filter = new FilteringParams();
            if (File.Exists(spectraFile))
            {
                try
                {
                    MsDataFile dataFile = Mzml.LoadAllStaticData(spectraFile, filter, searchForCorrectMs1PrecursorScan: false);
                    return dataFile;
                }
                catch (Exception ex)
                {
                    Console.Write("\nError loading calibrated mzML file, loading base mzML...");
                    //Console.WriteLine(String.Format("Error loading mzML file {0}", spectraFile));
                    //Console.Write(ex.ToString());
                }
            }
            // fallback attempts to locate the raw file if not using a filelist from FragPipe
            else if (File.Exists(Path.Combine(rawfilesDirectory, rawfileBase + CAL_INPUT_EXTENSION_FALLBACK)))
            {
                // support "_uncalibrated.mzML" if present and "_calibrated.mzML" is not
                spectraFile = Path.Combine(rawfilesDirectory, rawfileBase + CAL_INPUT_EXTENSION_FALLBACK);
                filter = new FilteringParams();
                try
                {
                    return Mgf.LoadAllStaticData(spectraFile, filter);
                }
                catch
                {
                    return null;
                }
            }
            else if (File.Exists(Path.Combine(rawfilesDirectory, rawfileBase + CAL_INPUT_EXTENSION_FALLBACK2)))
            {
                // support legacy "_calibrated.MGF" if present and other options are not
                spectraFile = Path.Combine(rawfilesDirectory, rawfileBase + CAL_INPUT_EXTENSION_FALLBACK2);
                filter = new FilteringParams();
                try
                {
                    return Mgf.LoadAllStaticData(spectraFile, filter);
                }
                catch
                {
                    return null;
                }
            }
            
            // no calibrated file (mzML/MGF) found (or loading failed), try falling back to raw mzML
            rawfileName = rawfileBase + ".mzML";
            spectraFile = Path.Combine(rawfilesDirectory, rawfileName);
            if (File.Exists(spectraFile))
            {
                try
                {
                    return Mzml.LoadAllStaticData(spectraFile, filter, searchForCorrectMs1PrecursorScan: false);
                }
                catch (Exception ex)
                {
                    Console.WriteLine("\nError loading base mzML file! No Localization can be done - please check the data file and try again.");
                    Console.WriteLine(ex.ToString());
                    return null;
                }
            }
            else
            {
                if (!mzmlNotFoundWarnings.ContainsKey(rawfileBase))
                {
                    Console.WriteLine("Warning: no mzML found for file {0}, PSMs from this file will NOT be localized!", rawfileBase);
                    mzmlNotFoundWarnings.Add(rawfileBase, true);
                }
                return null;
            }
        }

        /**
         * Reads FragPipe-generated list of raw file paths to use. Paths may be 
         * [filename]_calibrated.mzML, 
         * [filename]_uncalibrated.mzML, or 
         * [filename].mzML
         * Returns dictionary of base filename: full path for all files
         */
        private bool parseLCMSfilePaths()
        {
            lcmsPaths = new();
            try
            {
                string[] allLines = File.ReadAllLines(lcmsFileListPath);
                foreach (string line in allLines)
                {
                    // remove possible calibration signifiers
                    string strippedLine = line.Replace("\n", "").Trim();
                    string path = Path.GetFullPath(strippedLine);
                    string filename = Path.GetFileName(strippedLine);
                    string basename = Path.GetFileNameWithoutExtension(strippedLine);
                    basename = basename.Replace("_uncalibrated", "");
                    basename = basename.Replace("_calibrated", "");

                    if (filename.Contains("_calibrated"))
                    {
                        // calibrated mzML passed from FragPipe: look for calibrated, then uncalibrated, then original mzML
                        if (File.Exists(strippedLine))
                        {
                            lcmsPaths.Add(basename, path);
                        } 
                        else
                        {
                            // check uncal mzML
                            string uncalMzmlPath = path.Replace("_calibrated.mzML", "_uncalibrated.mzML");
                            if (File.Exists(uncalMzmlPath))
                            {
                                lcmsPaths.Add(basename, uncalMzmlPath);
                            }
                            else
                            {
                                // check original mzML
                                string originalMzmlPath = path.Replace("_calibrated.mzML", ".mzML");
                                if (File.Exists(originalMzmlPath))
                                {
                                    lcmsPaths.Add(basename, originalMzmlPath);
                                } 
                                else
                                {
                                    Console.WriteLine(String.Format("Error: no mzML file found for input {0}", strippedLine));
                                    return false;
                                }
                            }
                        }
                    } 
                    else
                    {
                        // non-calibrated file passed, read it exactly
                        if (File.Exists(strippedLine))
                        {
                            lcmsPaths.Add(basename, path);
                        }
                        else
                        {
                            Console.WriteLine(String.Format("Error: no mzML file found for input {0}", strippedLine));
                            return false;
                        }
                    }
                }
                return true;
            } catch
            {
                return false;
            }
        }


        /**
         * Main method to localize from MSFragger outputs. Loads the MSFragger PSM and scan-pair tables, parses PSM information
         * to determine which scans to localize and the input peptides/glycan masses for those scans, runs localization,
         * and writes results back to the PSM table. 
         */
        public int Localize()
        {
            List<string> allOutput = new();
            Dictionary<int, int> scanPairs = new();
            bool usingLcmsFilePath = false;
            if (scanpairFile != null)
            {
                // single scanPair file passed directly - use. Otherwise, will read below during PSM reading
                scanPairs = MSFragger_PSMTable.ParseScanPairTable(scanpairFile);
            }

            if (lcmsFileListPath != null)
            {
                bool parseSuccess = parseLCMSfilePaths();
                if (!parseSuccess && rawfilesDirectory == null)
                {
                    Console.WriteLine("Error: could not load FragPipe LCMS files list. Cannot localize");
                    return 1;
                }
                usingLcmsFilePath = true;
            }

            // read PSM table to prepare to pass scans to localizer
            MSFragger_PSMTable PSMtable = new(psmFile);

            bool overwritePrevious = false;
            if (PSMtable.Headers.Contains("O-Pair Score"))
            {
                Console.WriteLine("Error: PSM table contains previous O-Pair results. Overwriting existing results is not supported. Please try again with a fresh PSM table.");
                return 1;
            } 
            else
            {
                // write headers to PSM table
                GlycoSpectralMatch emptyMatch = new GlycoSpectralMatch();
                emptyMatch.localizerOutput = OPAIR_HEADERS;
                allOutput.Add(PSMtable.editPSMLine(emptyMatch, PSMtable.Headers.ToList(), GlycanBox.GlobalOGlycans, overwritePrevious, false));
            }

            Dictionary<string, bool> mzmlNotFoundWarnings = new();
            // timers and etc
            System.Diagnostics.Stopwatch timer = new();
            System.Diagnostics.Stopwatch totalTimer = new();
            totalTimer.Start();

            foreach (KeyValuePair<string, List<string>> rawfileEntry in PSMtable.RawfilePSMDict)
            {
                // load raw file and pairs table
                string rawfileName = rawfileEntry.Key + CAL_INPUT_EXTENSION;
                timer.Start();
                Console.Write(String.Format("\tLoading MS file {0}...", rawfileEntry.Key));
                MsDataFile currentMsDataFile = CheckAndLoadData(rawfileEntry.Key, rawfileName, mzmlNotFoundWarnings, usingLcmsFilePath);
                Console.Write(String.Format(" {0:0.0}s\n", timer.ElapsedMilliseconds * 0.001));
                timer.Reset();
                if (currentMsDataFile == null)
                {
                    // No raw data found! Write empty output for these PSMs
                    GlycoSpectralMatch emptyGSM = new GlycoSpectralMatch();
                    emptyGSM.localizerOutput = EMPTY_OUTPUT;
                    List<string> output = new();
                    foreach (string psmLine in rawfileEntry.Value)
                    {
                        string[] lineSplits = psmLine.Split('\t');
                        allOutput.Add(PSMtable.editPSMLine(emptyGSM, lineSplits.ToList(), GlycanBox.GlobalOGlycans, overwritePrevious, false));
                    }
                }
                else
                {
                    timer.Start();
                    List<MsDataScan> allScans = currentMsDataFile.GetAllScansList();
                    // Dictionary by scan number for lookup
                    Dictionary<int, MsDataScan> dataScansDict = new();
                    foreach (MsDataScan scan in allScans)
                    {
                        dataScansDict.Add(scan.OneBasedScanNumber, scan);
                    }
                    string scanpairName = rawfileEntry.Key + ".pairs";
                    string pairsFile = Path.Combine(rawfilesDirectory, scanpairName);
                    scanPairs = MSFragger_PSMTable.ParseScanPairTable(pairsFile);
                    Dictionary<int, string> scanDict = PSMtable.GetScanDict(rawfileEntry.Key);
                    Dictionary<int, string> output = new ();

                    // timer for PSM progress
                    int analyzedCount = 0;
                    double pct = 0;
                    Timer threadedTimer = new Timer(_ =>
                    {
                        Console.Write($"\x000D\t\t[progress: {analyzedCount}/{scanDict.Count} ({pct:0.0}%)]");
                    }, null, TimeSpan.Zero, TimeSpan.FromSeconds(10));

                    // localize all PSMs in parallel
                    ParallelOptions options = new();
                    options.MaxDegreeOfParallelism = 1;
                    Parallel.ForEach(scanDict, options, psmEntry =>
                    {
                        string psmOutput = LocalizePSM(dataScansDict, scanPairs, psmEntry.Value, PSMtable, overwritePrevious, rawfileName);
                        lock (output) // to synchronize access to the dictionary
                        {
                            output.Add(psmEntry.Key, psmOutput);
                            analyzedCount++;
                            pct = 100 * analyzedCount / (double) scanDict.Count;
                        }
                    });
                    Console.Write($"\x000D\t\t[progress: {analyzedCount}/{scanDict.Count} ({pct}%)]");
                    threadedTimer.Dispose();

                    // save results back to output list
                    List<string> sortedPSMs = output.OrderBy(x => x.Key)
                                                    .Select(x => x.Value)
                                                    .ToList();
                    allOutput.AddRange(sortedPSMs);
                    Console.WriteLine(String.Format("...done. {0} PSMs analyzed in {1:0.0}s", sortedPSMs.Count, timer.ElapsedMilliseconds * 0.001));
                    timer.Reset();
                }
            }

            File.WriteAllLines(psmFile, allOutput.ToArray());
            totalTimer.Stop();
            Console.WriteLine("\nFinished localization in {0:0.0} s", totalTimer.ElapsedMilliseconds * 0.001);
            return 0;
        }

        private string LocalizePSM(Dictionary<int, MsDataScan> dataScansDict, Dictionary<int, int> scanPairs, string PSMline, MSFragger_PSMTable PSMtable, bool overwritePrevious, string currentRawfile)
        {
            int scanNum = MSFragger_PSMTable.GetScanNum(PSMline);

            // Parse PSM information
            string[] lineSplits = PSMline.Split("\t");
            string spectrumString = lineSplits[PSMtable.SpectrumCol];
            string peptide = lineSplits[PSMtable.PeptideCol];
            double deltaMass = Double.Parse(lineSplits[PSMtable.DeltaMassCol]);
            string assignedMods = lineSplits[PSMtable.AssignedModCol];
            double precursorMZ = Double.Parse(lineSplits[PSMtable.PrecursorMZCol]);
            int precursorCharge = MSFragger_PSMTable.GetScanCharge(spectrumString);

            if (deltaMass < 3.5 && deltaMass > -1.5)
            {
                // unmodified peptide - no localization
                GlycoSpectralMatch emptyGSM = new GlycoSpectralMatch();
                emptyGSM.localizerOutput = EMPTY_OUTPUT;
                return PSMtable.editPSMLine(emptyGSM, lineSplits.ToList(), GlycanBox.GlobalOGlycans, overwritePrevious, false);
            }

            // retrieve the child scan num
            int childScanNum;
            if (scanPairs.ContainsKey(scanNum))
            {
                childScanNum = scanPairs[scanNum];
            }
            else
            {
                // no paired child scan found - do not attempt localization
                GlycoSpectralMatch emptyGSM = new GlycoSpectralMatch();
                emptyGSM.localizerOutput = "No paired scan" + EMPTY_OUTPUT;
                return PSMtable.editPSMLine(emptyGSM, lineSplits.ToList(), GlycanBox.GlobalOGlycans, overwritePrevious, false);
            }

            // retrieve the child scan
            MsDataScan ms2Scan;
            try
            {
                ms2Scan = dataScansDict[childScanNum];
            }
            catch (KeyNotFoundException)
            {
                Console.Out.WriteLine(String.Format("Error: MS2 scan {0} not found in the spectrum file. No localization performed", childScanNum));
                GlycoSpectralMatch emptyGSM = new GlycoSpectralMatch();
                emptyGSM.localizerOutput = EMPTY_OUTPUT;
                return PSMtable.editPSMLine(emptyGSM, lineSplits.ToList(), GlycanBox.GlobalOGlycans, overwritePrevious, false);
            } 

            // finalize spectrum for search
            IsotopicEnvelope[] neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2Scan, 4, 3);
            var scan = new Ms2ScanWithSpecificMass(ms2Scan, precursorMZ, precursorCharge, currentRawfile, 4, 3, neutralExperimentalFragments);

            // initialize peptide with all non-glyco mods
            PeptideWithSetModifications peptideWithMods = getPeptideWithMSFraggerMods(assignedMods, peptide);

            // get oxo 138/144 ratio from the HCD scan
            Dictionary<FilterRule, bool> oxoniumsToFilter = new();
            double ratio;
            try
            {
                MsDataScan hcdDataScan = dataScansDict[scanNum];
                var hcdScan = new Ms2ScanWithSpecificMass(hcdDataScan, precursorMZ, precursorCharge, currentRawfile, 4, 3, neutralExperimentalFragments);
                ratio = hcdScan.ComputeOxoRatio(OXO138, OXO144, ProductMassTolerance.Value);
                // todo: add parameter
                oxoniumsToFilter = hcdScan.FindOxoniums(GlobalVariables.OxoniumFilters, ProductMassTolerance.Value);
            }
            catch (KeyNotFoundException)
            {
                Console.WriteLine(string.Format("Error: HCD scan {0} not found in file {1}, could not calculate oxonium ratio or use oxonium filter", scanNum, currentRawfile));
                ratio = -1;
            }

            // finally, run localizer
            GlycoSpectralMatch gsm = LocalizeOGlyc(scan, peptideWithMods, deltaMass, ratio, oxoniumsToFilter);
            return PSMtable.editPSMLine(gsm, lineSplits.ToList(), GlycanBox.GlobalOGlycans, overwritePrevious, true);

        }

        /**
         * Single scan O-glycan localization, given a peptide and glycan mass from MSFragger search
         */
        public GlycoSpectralMatch LocalizeOGlyc(Ms2ScanWithSpecificMass ms2Scan, PeptideWithSetModifications peptide, double glycanDeltaMass, double oxoRatio, Dictionary<FilterRule, bool> oxoniumBools)
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
            bool oxoFilteredOut = false;
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
                    // filter based on oxonium ions if requested
                    bool oxoFilter = CheckOxoniumFilter(iDLow, oxoniumBools);
                    if (!oxoFilter)
                    {
                        iDLow++;
                        oxoFilteredOut = true;
                        continue;
                    }

                    // only consider possibilities with enough sites on the peptide to accomodate all glycans
                    if (_modPos.Length >= GlycanBox.OGlycanBoxes[iDLow].OGlycanCount)
                    {
                        LocalizationGraph localizationGraph = new LocalizationGraph(_modPos, modMotifs, GlycanBox.OGlycanBoxes[iDLow], GlycanBox.OGlycanBoxes[iDLow].ChildGlycanBoxes, iDLow);
                        LocalizationGraph.LocalizeMod(localizationGraph, ms2Scan, ProductMassTolerance,
                            products.Where(v => v.ProductType == ProductType.c || v.ProductType == ProductType.zDot).ToList(),
                            GlycoPeptides.GetLocalFragmentGlycan, GlycoPeptides.GetUnlocalFragmentGlycan);
                        graphs.Add(localizationGraph);
                    }
                    iDLow++;
                }
            }

            // save results
            GlycoSpectralMatch gsm;
            if (graphs.Count == 0)
            {
                // no matching glycan boxes found to this delta mass, no localization can be done! 
                gsm = new GlycoSpectralMatch();
                // keep the psm line the same length as all other outputs
                if (oxoFilteredOut)
                {
                    // filtered out a possible match - note that
                    gsm.localizerOutput = "No match after oxonium filtering" + EMPTY_OUTPUT_WITH_PAIRED_SCAN + string.Format("\t{0}", ms2Scan.TheScan.OneBasedScanNumber);
                }
                else
                {
                    gsm.localizerOutput = "No match to glycan delta mass" + EMPTY_OUTPUT_WITH_PAIRED_SCAN + string.Format("\t{0}", ms2Scan.TheScan.OneBasedScanNumber);
                }
                    return gsm;      
            }
            var best_graph = graphs.OrderByDescending(p => p.TotalScore).First();
            gsm = new GlycoSpectralMatch(new List<LocalizationGraph> { best_graph }, GlycoType.OGlycoPep);
            gsm.ScanInfo_p = ms2Scan.TheScan.MassSpectrum.Size * ProductMassTolerance.GetRange(1000).Width / ms2Scan.TheScan.MassSpectrum.Range.Width;
            gsm.Thero_n = products.Count();
            gsm.oxoRatio = oxoRatio;
            gsm.NGlycanMotifExist = NglycMotifRegex.IsMatch(peptide.BaseSequence);
            GlycoSite.GlycoLocalizationCalculation(gsm, gsm.GlycanType, DissociationType.HCD, DissociationType.EThcD);
            gsm.localizerOutput = gsm.WriteLine(null) + string.Format("\t{0}", ms2Scan.TheScan.OneBasedScanNumber);

            return gsm;
        }

        /**
         * Check if a particular glycan has the required oxonium ion(s) given its composition. 
         * Returns false if the glycan is NOT allowed based on oxonium, true if it is allowed. 
         * 
         * Requires a dictionary of FilterRule: bool indicating whether the required ion(s) have
         * been found for a given rule for this spectrum. 
         */
        public bool CheckOxoniumFilter(int iDLow, Dictionary<FilterRule, bool> oxoniumBools)
        {
            foreach (FilterRule rule in GlobalVariables.OxoniumFilters)
            {
                bool matchedKind = true;
                foreach (KeyValuePair<byte, int> requiredMonosaccharide in rule.MonosaccharidesRequired)
                {
                    if (GlycanBox.OGlycanBoxes[iDLow].Kind[requiredMonosaccharide.Key] < requiredMonosaccharide.Value)
                    {
                        matchedKind = false;
                    }
                }
                // glycan has the monosaccharide(s) of interest to this filter: make sure oxonium(s) found
                if (matchedKind && !oxoniumBools[rule])
                {
                    return false;
                }
            }
            return true;
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
                    string locationRestriction;
                    // handle terminal mods as well
                    if (split.Contains("N-term"))
                    {
                        // N-terminal mod. Set position to -1 and residue type to "X"
                        residueType = "X";
                        position = -1;
                        locationRestriction = "Peptide N-terminal.";
                    } 
                    else if (split.Contains("C-term"))
                    {
                        // C-terminal mod. Set position to peptide length and residue type to "X"
                        residueType = "X";
                        position = peptide.Length;
                        locationRestriction = "Peptide C-terminal.";
                    }
                    else
                    {
                        residueType = massSplits[0][massSplits[0].Length - 1].ToString();
                        position = Int32.Parse(massSplits[0].Substring(0, massSplits[0].Length - 1)) - 1;   // MSFragger mod positions are 1-indexed, convert to 0-index
                        locationRestriction = "Anywhere.";
                    }

                    double mass = double.Parse(massSplits[1].Replace(")", ""));
                    string modName = "(" + massSplits[1];
                    // Generate a Modification of type "MSFragger" with ID "(mass)" 
                    ModificationMotif.TryGetMotif(residueType, out ModificationMotif motif2);
                    Modification mod = new Modification(_originalId: modName, _modificationType: modType, _target: motif2, _locationRestriction: locationRestriction, _monoisotopicMass: mass);
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
