using Chemistry;
using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;
using MzLibUtil;
using EngineLayer.GlycoSearch;
using MathNet.Numerics.LinearRegression;
using IO.Mgf;

namespace PTMLocalizationTest
{
    [TestFixture]
    public class Tests
    {
        private static GlycanBox[] OGlycanBoxes { get; set; }

        [OneTimeSetUp]
        public static void Setup()
        {
            GlobalVariables.SetUpGlobalVariables();
            GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanLocations.Where(p => p.Contains("OGlycan.gdb")).First(), true, true).ToArray();
            GlycanBox.GlobalOGlycanMods = GlycanBox.BuildGlobalOGlycanMods(GlycanBox.GlobalOGlycans).ToArray();
            OGlycanBoxes = GlycanBox.BuildGlycanBoxes(3, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods).OrderBy(p => p.Mass).ToArray();
        }

        [Test]
        public static void OGlycoTest_LocalizeMod()
        {
            //This test proves that the LocalizeMod method is the same as LocalizeOGlycan. 

            //Get glycanBox
            var glycanBox = OGlycanBoxes[19];

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
            LocalizationGraph localizationGraph = new LocalizationGraph(modPos, boxMotifs, glycanBox, boxes, -1);

            //LocalizeMod
            LocalizationGraph localizationGraph0 = new LocalizationGraph(modPos, boxMotifs, glycanBox, boxes, -1);
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

        }
    }
}