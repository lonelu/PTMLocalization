using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MzLibUtil;
using EngineLayer.GlycoSearch;
using IO.Mgf;
using System.Diagnostics;
using PTMLocalization;
using System;

namespace PTMLocalizationTest
{
    [TestFixture]
    public class Tests
    {

        [OneTimeSetUp]
        public static void Setup()
        {
            GlobalVariables.SetUpGlobalVariables(null, null);
            GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanLocations.Where(p => p.Contains("OGlycan.gdb")).First(), true, true).ToArray();
            GlycanBox.GlobalOGlycanMods = GlycanBox.BuildGlobalOGlycanMods(GlycanBox.GlobalOGlycans).ToArray();
            GlycanBox.OGlycanBoxes = GlycanBox.BuildGlycanBoxes(3, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods).OrderBy(p => p.Mass).ToArray();
        }


        [Test]
        public static void OGlycoTest_LocalizeMod()
        {
            //This test proves that the LocalizeMod method is the same as LocalizeOGlycan. 

            //Get glycanBox
            var boxId = 19;
            var glycanBox = GlycanBox.OGlycanBoxes[boxId];

            //Get unmodified peptide, products, allPossible modPos and all boxes.
            Protein prot = new Protein(sequence: "TTGSLEPSSGASGPQVSSVK", accession: "", isDecoy: false);
            PeptideWithSetModifications peptide = new PeptideWithSetModifications(prot, null, 1, 20, CleavageSpecificity.Full, "TTGSLEPSSGASGPQVSSVK", 0, new Dictionary<int, Modification>(), 0, null);

            List<Product> products = new List<Product>();
            peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);

            int[] modPos = GlycoPeptides.GetPossibleModSites(peptide, new string[] { "S", "T" }).OrderBy(v => v).ToArray();
            var boxes = GlycanBox.BuildChildGlycanBoxes(glycanBox.ModIds, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods).ToArray();
            Assert.That(boxes.Count() == 6);

            //Get Unlocal Fragment
            var unlocalCost = GlycoPeptides.GetUnlocalFragment(products, modPos, glycanBox);
            Assert.That(unlocalCost.Count == 5); //Basicly, the unlocal are c/z ions that don't localize glycosylation. 

            //Get scan
            Tolerance ProductMassTolerance = new PpmTolerance(20);

            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_4565.mgf");
            FilteringParams filter = new FilteringParams();
            var myMSDataFile = Mgf.LoadAllStaticData(spectraFile, filter);

            var msNScans = myMSDataFile.GetAllScansList().Where(x => x.MsnOrder > 1).ToArray();
            var ms2Scan = msNScans.Where(p => p.MsnOrder == 2).ToArray()[0];


            double precursorMZ = ms2Scan.SelectedIonMonoisotopicGuessMz.Value;
            //double precursorMZ = ms2Scan.SelectedIonMZ.Value;
            int precursorCharge = ms2Scan.SelectedIonChargeStateGuess.Value;

            IsotopicEnvelope[] neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2Scan, 4, 3);
            var scan = new Ms2ScanWithSpecificMass(ms2Scan, precursorMZ, precursorCharge, spectraFile, 4, 3, neutralExperimentalFragments);

            ////Known peptideWithMod match.
            //var peptideWithMod = GlycoPeptides.GlyGetTheoreticalPeptide(new int[3] { 10, 2, 3 }, peptide, glycanBox, GlycanBox.GlobalOGlycanMods);
            //Assert.That(peptideWithMod.FullSequence == "T[O-Glycosylation:H1N1 on X]T[O-Glycosylation:H1N1 on X]GSLEPSS[O-Glycosylation:N1 on X]GASGPQVSSVK");
            //List<Product> knownProducts = GlycoPeptides.OGlyGetTheoreticalFragments(DissociationType.EThcD, peptide, peptideWithMod);
            //var matchedKnownFragmentIons = MetaMorpheusEngine.MatchFragmentIons(scans.First(), knownProducts, commonParameters);

            //Graph Localization
            var boxMotifs = new string[] { "S/T", "S/T", "S/T", "S/T", "S/T", "S/T", "S/T", "S/T" };

            //LocalizeMod
            LocalizationGraph localizationGraph0 = new LocalizationGraph(modPos, boxMotifs, glycanBox, boxes, boxId);
            //LocalizationGraph.LocalizeMod(localizationGraph0, scans.First(), commonParameters.ProductMassTolerance, 
            //    products.Where(v => v.ProductType == ProductType.c || v.ProductType == ProductType.zDot).ToList(),
            //    GlycoPeptides.GetLocalFragment, GlycoPeptides.GetUnlocalFragment);
            LocalizationGraph.LocalizeMod(localizationGraph0, scan, ProductMassTolerance,
    products.Where(v => v.ProductType == ProductType.c || v.ProductType == ProductType.zDot).ToList(),
    GlycoPeptides.GetLocalFragmentGlycan, GlycoPeptides.GetUnlocalFragmentGlycan);
            var allPaths0 = LocalizationGraph.GetAllHighestScorePaths(localizationGraph0.array, localizationGraph0.ChildModBoxes);
            var knowPath0 = new int[8] { 2, 4, 4, 4, 5, 5, 5, 5 };
            Assert.That(Enumerable.SequenceEqual(knowPath0, allPaths0[0]));

            Assert.That(27.097639419729177 < localizationGraph0.TotalScore + 0.00001 && 27.097639419729177 > localizationGraph0.TotalScore - 0.00001);


            var gsm = new GlycoSpectralMatch(new List<LocalizationGraph> { localizationGraph0 }, GlycoType.OGlycoPep);
            gsm.ScanInfo_p = scan.TheScan.MassSpectrum.Size * ProductMassTolerance.GetRange(1000).Width / scan.TheScan.MassSpectrum.Range.Width;
            gsm.Thero_n = products.Count();
            GlycoSite.GlycoLocalizationCalculation(gsm, gsm.GlycanType, DissociationType.HCD, DissociationType.EThcD);

            var output = gsm.WriteLine(null);
        }


        [Test]
        public static void OGlycoTest_LocalizeMod2()
        {
            //In this test, the peptide contain Common mod.

            //Get unmodified peptide, products, allPossible modPos.

            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            Modification mod = new Modification(_originalId: "Deamidation on N", _modificationType: "Common Artifact", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 0.984016);

            PeptideWithSetModifications peptide = new PeptideWithSetModifications("STN[Common Artifact:Deamidation on N]ASTVPFR",
                new Dictionary<string, Modification> { { "Deamidation on N", mod } });

            List<Product> products = new List<Product>();
            peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);

            int[] modPos = GlycoPeptides.GetPossibleModSites(peptide, new string[] { "S", "T" }).OrderBy(v => v).ToArray();

            //Get scan
            Tolerance PrecursorTolerance = new PpmTolerance(10);
            Tolerance ProductMassTolerance = new PpmTolerance(20);

            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_5467.mgf"); //The Scans in the file are a pair.
            FilteringParams filter = new FilteringParams();
            var myMSDataFile = Mgf.LoadAllStaticData(spectraFile, filter);

            var msNScans = myMSDataFile.GetAllScansList().Where(x => x.MsnOrder > 1).ToArray();
            var ms2Scan = msNScans.Where(p => p.MsnOrder == 2).ToArray()[1];

            double precursorMZ = ms2Scan.SelectedIonMonoisotopicGuessMz.Value;
            //double precursorMZ = ms2Scan.SelectedIonMZ.Value;
            int precursorCharge = ms2Scan.SelectedIonChargeStateGuess.Value;

            IsotopicEnvelope[] neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2Scan, 4, 3);
            var scan = new Ms2ScanWithSpecificMass(ms2Scan, precursorMZ, precursorCharge, spectraFile, 4, 3, neutralExperimentalFragments);


            List<int> n_modPos = new List<int>();
            //var n_modPos = GlycoPeptides.GetPossibleModSites(peptide, new string[] { "Nxt", "Nxs" }).OrderBy(p => p).ToList();
            var o_modPos = GlycoPeptides.GetPossibleModSites(peptide, new string[] { "S", "T" }).OrderBy(p => p).ToList();

            int[] _modPos;
            string[] modMotifs;
            GlycoPeptides.GetModPosMotif(GlycoType.OGlycoPep, n_modPos.ToArray(), o_modPos.ToArray(), out _modPos, out modMotifs);

            var possibleGlycanMassLow = PrecursorTolerance.GetMinimumValue(scan.PrecursorMass) - peptide.MonoisotopicMass;
            var possibleGlycanMassHigh = PrecursorTolerance.GetMaximumValue(scan.PrecursorMass) - peptide.MonoisotopicMass;
            int iDLow = GlycoPeptides.BinarySearchGetIndex(GlycanBox.OGlycanBoxes.Select(p => p.Mass).ToArray(), possibleGlycanMassLow);

            List<LocalizationGraph> graphs = new List<LocalizationGraph>();
            while (iDLow < GlycanBox.OGlycanBoxes.Length && GlycanBox.OGlycanBoxes[iDLow].Mass <= possibleGlycanMassHigh)
            {
                LocalizationGraph localizationGraph = new LocalizationGraph(_modPos, modMotifs, GlycanBox.OGlycanBoxes[iDLow], GlycanBox.OGlycanBoxes[iDLow].ChildGlycanBoxes, iDLow);
                LocalizationGraph.LocalizeMod(localizationGraph, scan, ProductMassTolerance,
                    products.Where(v => v.ProductType == ProductType.c || v.ProductType == ProductType.zDot).ToList(),
                    GlycoPeptides.GetLocalFragmentGlycan, GlycoPeptides.GetUnlocalFragmentGlycan);
                graphs.Add(localizationGraph);
                iDLow++;
            }
            var best_graph = graphs.OrderByDescending(p => p.TotalScore).First();
            var gsm = new GlycoSpectralMatch(new List<LocalizationGraph> { best_graph }, GlycoType.OGlycoPep);
            gsm.ScanInfo_p = scan.TheScan.MassSpectrum.Size * ProductMassTolerance.GetRange(1000).Width / scan.TheScan.MassSpectrum.Range.Width;
            gsm.Thero_n = products.Count();
            GlycoSite.GlycoLocalizationCalculation(gsm, gsm.GlycanType, DissociationType.HCD, DissociationType.EThcD);

            var output = gsm.WriteLine(null);
        }

        

        [Test]
        public static void OGlycoTest_FragPipe_FullRun()
        {
            // load test parameters, spectrum files, and MSFragger PSM table
            Tolerance ProductMassTolerance = new PpmTolerance(10);
            Tolerance PrecursorMassTolerance = new PpmTolerance(30);
            string psmFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\psm.tsv");
            string scanpairFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\_scan-pairs.tsv");
            string rawfileDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData");
            string glycoDatabase = GlobalVariables.OGlycanLocations.Where(p => p.Contains("OGlycan.gdb")).First();
            int maxNumGlycans = 3;
            int[] isotopes = { 0, 1, 2 };
            bool filterOxonium = false;

            var localizer = new MSFragger_RunLocalization(psmFile, scanpairFile, rawfileDirectory, null, glycoDatabase, maxNumGlycans, PrecursorMassTolerance, ProductMassTolerance, isotopes, filterOxonium);
            localizer.Localize();
        }

        [Test]
        // original single scan MSFragger test
        public static void OGlycoTest_FragPipeInput()
        {
            //int scanNum, string peptideSequence, double deltaMass

            // Test method for moving towards MSFragger/FragPipe inputs/output

            //Get scan
            Tolerance ProductMassTolerance = new PpmTolerance(20);

            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_4565.mgf");
            //string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_calibrated.mgf");

            FilteringParams filter = new FilteringParams();
            var myMSDataFile = Mgf.LoadAllStaticData(spectraFile, filter);

            var msNScans = myMSDataFile.GetAllScansList().Where(x => x.MsnOrder > 1).ToArray();
            var ms2Scan = msNScans.Where(p => p.MsnOrder == 2).ToArray()[0];

            double precursorMZ = ms2Scan.SelectedIonMonoisotopicGuessMz.Value;
            //double precursorMZ = ms2Scan.SelectedIonMZ.Value;
            int precursorCharge = ms2Scan.SelectedIonChargeStateGuess.Value;

            IsotopicEnvelope[] neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2Scan, 4, 3);
            var scan = new Ms2ScanWithSpecificMass(ms2Scan, precursorMZ, precursorCharge, spectraFile, 4, 3, neutralExperimentalFragments);


            //Get unmodified peptide, products
            Protein prot = new Protein(sequence: "TTGSLEPSSGASGPQVSSVK", accession: "", isDecoy: false);
            PeptideWithSetModifications peptide = new PeptideWithSetModifications(prot, null, 1, 20, CleavageSpecificity.Full, "TTGSLEPSSGASGPQVSSVK", 0, new Dictionary<int, Modification>(), 0, null);

            List<Product> products = new List<Product>();
            peptide.Fragment(DissociationType.ETD, FragmentationTerminus.Both, products);
            int[] modPos = GlycoPeptides.GetPossibleModSites(peptide, new string[] { "S", "T" }).OrderBy(v => v).ToArray();

            double deltaMass = scan.PrecursorMass - peptide.MonoisotopicMass;
            //Get glycanBox (loop to get all different glycans searched)
            double bestScore = 0;
            GlycoSpectralMatch bestGSM = new GlycoSpectralMatch();
            for (int boxId = 0; boxId < GlycanBox.OGlycanBoxes.Length; boxId++)
            {
                var glycanBox = GlycanBox.OGlycanBoxes[boxId];

                // only search glycans within the delta mass window
                // TODO: change to accurate mass matching with isotopes
                if (-1.5 < (glycanBox.Mass - deltaMass) && (glycanBox.Mass - deltaMass) < 3.5)
                {
                    var boxes = GlycanBox.BuildChildGlycanBoxes(glycanBox.ModIds, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods).ToArray();
                    //Assert.That(boxes.Count() == 6);

                    //Get Unlocal Fragment
                    var unlocalCost = GlycoPeptides.GetUnlocalFragment(products, modPos, glycanBox);
                    //Assert.That(unlocalCost.Count == 5); //Basicly, the unlocal are c/z ions that don't localize glycosylation. 

                    //Graph Localization
                    var boxMotifs = new string[] { "S/T", "S/T", "S/T", "S/T", "S/T", "S/T", "S/T", "S/T" };

                    //LocalizeMod
                    LocalizationGraph localizationGraph0 = new LocalizationGraph(modPos, boxMotifs, glycanBox, boxes, boxId);

                    LocalizationGraph.LocalizeMod(localizationGraph0, scan, ProductMassTolerance,
            products.Where(v => v.ProductType == ProductType.c || v.ProductType == ProductType.zDot).ToList(),
            GlycoPeptides.GetLocalFragmentGlycan, GlycoPeptides.GetUnlocalFragmentGlycan);
                    var allPaths0 = LocalizationGraph.GetAllHighestScorePaths(localizationGraph0.array, localizationGraph0.ChildModBoxes);
                    var knowPath0 = new int[8] { 2, 4, 4, 4, 5, 5, 5, 5 };
                    //Assert.That(Enumerable.SequenceEqual(knowPath0, allPaths0[0]));

                    //Assert.That(27.097639419729177 < localizationGraph0.TotalScore + 0.00001 && 27.097639419729177 > localizationGraph0.TotalScore - 0.00001);

                    if (localizationGraph0.TotalScore > bestScore)
                    {
                        // save info
                        bestScore = localizationGraph0.TotalScore;
                        var gsm = new GlycoSpectralMatch(new List<LocalizationGraph> { localizationGraph0 }, GlycoType.OGlycoPep);
                        gsm.ScanInfo_p = scan.TheScan.MassSpectrum.Size * ProductMassTolerance.GetRange(1000).Width / scan.TheScan.MassSpectrum.Range.Width;
                        gsm.Thero_n = products.Count();
                        GlycoSite.GlycoLocalizationCalculation(gsm, gsm.GlycanType, DissociationType.HCD, DissociationType.EThcD);
                        bestGSM = gsm;
                    }
                }
            }
            var output = bestGSM.WriteLine(null);
            Debug.WriteLine(output);
        }
    }
}