using MzLibUtil;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using EngineLayer;
using MassSpectrometry;
using Chemistry;

namespace EngineLayer.GlycoSearch
{
    public class RunLocalization
    {
        public static readonly double ToleranceForMassDifferentiation = 1e-9;
        private readonly int OxoniumIon204Index = 9; //Check Glycan.AllOxoniumIons

        private GlycoSearchType GlycoSearchType;
        private readonly bool MixedGlycanAllowed;
        private readonly int _maxOGlycanNum;
        private readonly int _maxNGlycanNum;

        private readonly string _oglycanDatabase;
        private readonly string _nglycanDatabase;
        public RunLocalization(string oglycanDatabase, string nglycanDatabase, GlycoSearchType glycoSearchType, bool mixedGlycanAllowed, int maxOGlycanNum, int maxNGlycanNum)
        {
            this.GlycoSearchType = glycoSearchType;
            this.MixedGlycanAllowed = mixedGlycanAllowed;
            this._maxOGlycanNum = maxOGlycanNum;
            this._maxNGlycanNum = maxNGlycanNum;

            this._oglycanDatabase = oglycanDatabase;
            this._nglycanDatabase = nglycanDatabase;

            if (glycoSearchType == GlycoSearchType.OGlycanSearch)
            {
                GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanLocations.Where(p => System.IO.Path.GetFileName(p) == _oglycanDatabase).First(), true, true).ToArray();
                GlycanBox.GlobalOGlycanMods = GlycanBox.BuildGlobalOGlycanMods(GlycanBox.GlobalOGlycans).ToArray();
                GlycanBox.OGlycanBoxes = GlycanBox.BuildGlycanBoxes(_maxOGlycanNum, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods).OrderBy(p => p.Mass).ToArray();
            }
            else if (glycoSearchType == GlycoSearchType.NGlycanSearch)
            {
                GlycanBox.GlobalNGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanLocations.Where(p => System.IO.Path.GetFileName(p) == _nglycanDatabase).First(), true, false).OrderBy(p => p.Mass).ToArray();
                GlycanBox.GlobalNGlycanMods = GlycanBox.BuildGlobalNGlycanMods(GlycanBox.GlobalNGlycans).ToArray();
                GlycanBox.NGlycanBoxes = GlycanBox.BuildGlycanBoxes(_maxNGlycanNum, GlycanBox.GlobalNGlycans, GlycanBox.GlobalNGlycanMods).OrderBy(p => p.Mass).ToArray();
                //TO THINK: Glycan Decoy database.

            }
            else if (glycoSearchType == GlycoSearchType.N_O_GlycanSearch)
            {
                GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanLocations.Where(p => System.IO.Path.GetFileName(p) == _oglycanDatabase).First(), true, true).ToArray();
                GlycanBox.GlobalOGlycanMods = GlycanBox.BuildGlobalOGlycanMods(GlycanBox.GlobalOGlycans).ToArray();
                GlycanBox.OGlycanBoxes = GlycanBox.BuildGlycanBoxes(_maxOGlycanNum, GlycanBox.GlobalOGlycans, GlycanBox.GlobalOGlycanMods).OrderBy(p => p.Mass).ToArray();

                GlycanBox.GlobalNGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanLocations.Where(p => System.IO.Path.GetFileName(p) == _nglycanDatabase).First(), true, false).ToArray();
                GlycanBox.GlobalNGlycanMods = GlycanBox.BuildGlobalNGlycanMods(GlycanBox.GlobalNGlycans).ToArray();
                GlycanBox.NGlycanBoxes = GlycanBox.BuildGlycanBoxes(_maxNGlycanNum, GlycanBox.GlobalNGlycans, GlycanBox.GlobalNGlycanMods).OrderBy(p => p.Mass).ToArray();

                if (MixedGlycanAllowed)
                {
                    GlycanBox.GlobalMixedGlycans = GlycanBox.GlobalOGlycans.Concat(GlycanBox.GlobalNGlycans).ToArray();
                    GlycanBox.GlobalMixedGlycanMods = GlycanBox.GlobalOGlycanMods.Concat(GlycanBox.GlobalNGlycanMods).ToArray();
                    GlycanBox.MixedModBoxes = GlycanBox.BuildMixedGlycanBoxes(GlycanBox.GlobalMixedGlycans, GlycanBox.GlobalMixedGlycanMods, _maxOGlycanNum, _maxNGlycanNum, GlycanBox.GlobalOGlycans.Length, GlycanBox.GlobalNGlycans.Length).OrderBy(p => p.Mass).ToArray();
                }
            }
        }

        public List<LocalizationGraph> RunLocalizationGraph(Ms2ScanWithSpecificMass theScan, List<Product> mainProducts, 
            List<Product> childProducts, int[] modPos, string[] modMotifs, int[] nPos, List<int> glycanBoxIds, 
            GlycanBox[] globalBoxes, Tolerance ProductMassTolerance)
        {
            List<LocalizationGraph> localizationGraphs = new List<LocalizationGraph>();

            double bestLocalizedScore = 0;

            foreach (var iDLow in glycanBoxIds)
            {
                LocalizationGraph localizationGraph = new LocalizationGraph(modPos, modMotifs, globalBoxes[iDLow], globalBoxes[iDLow].ChildGlycanBoxes, iDLow);

                LocalizationGraph.LocalizeMod(localizationGraph, theScan, ProductMassTolerance, mainProducts, GlycoPeptides.GetLocalFragmentGlycan, GlycoPeptides.GetUnlocalFragmentGlycan);

                if (theScan.ChildScans.Count > 0)
                {
                    foreach (var childScan in theScan.ChildScans)
                    {
                        LocalizationGraph.LocalizeMod(localizationGraph, childScan, ProductMassTolerance, childProducts, GlycoPeptides.GetLocalFragmentGlycan, GlycoPeptides.GetUnlocalFragmentGlycan);
                    }
                }

                double currentLocalizationScore = localizationGraph.TotalScore;

                if (currentLocalizationScore > bestLocalizedScore)
                {
                    bestLocalizedScore = currentLocalizationScore;
                    localizationGraphs.Clear();
                    localizationGraphs.Add(localizationGraph);
                }
                else if ((currentLocalizationScore <= bestLocalizedScore + 0.00000001 && currentLocalizationScore >= bestLocalizedScore - 0.00000001))
                {
                    localizationGraphs.Add(localizationGraph);
                }

            }

            return localizationGraphs;
        }

        public void ExtractGraphInfo(Ms2ScanWithSpecificMass theScan, PeptideWithSetModifications theScanBestPeptide, 
            List<LocalizationGraph> localizationGraphs, int[] NPos, GlycoType GType, DissociationType dissociationType, DissociationType MS2ChildScanDissociationType,
            Tolerance ProductMassTolerance, Glycan[] globalglycans, Modification[] globalmods)
        {
            var firstPath = LocalizationGraph.GetFirstPath(localizationGraphs[0].array, localizationGraphs[0].ChildModBoxes);
            var route = LocalizationGraph.GetLocalizedPath(localizationGraphs[0], firstPath);

            var peptideWithMod = GlycoPeptides.GlyGetTheoreticalPeptide(route, theScanBestPeptide, globalmods);

            Glycan[] glycans = new Glycan[localizationGraphs.First().ModBox.ModCount];
            for (int i = 0; i < localizationGraphs.First().ModBox.ModCount; i++)
            {
                glycans[i] = globalglycans[((GlycanBox)localizationGraphs.First().ModBox).ModIds[i]];
            }

            List<int> npos = route.Mods.Select(p => p.ModSite).ToArray().Intersect(NPos).ToList();

            var fragmentsForEachGlycoPeptide = GlycoPeptides.GlyGetTheoreticalFragments(GType, dissociationType, theScanBestPeptide, peptideWithMod, npos, glycans);

            var p = theScan.TheScan.MassSpectrum.Size * ProductMassTolerance.GetRange(1000).Width / theScan.TheScan.MassSpectrum.Range.Width;

            int n = fragmentsForEachGlycoPeptide.Where(v => v.ProductType != ProductType.D && v.ProductType != ProductType.Ycore && v.ProductType != ProductType.Y).Count();

            foreach (var childScan in theScan.ChildScans)
            {
                //People always use CID with low res. This special code works for Nic Scott's data. (High-HCD, High-EThcD, Low-CID, High-secHCD)
                if (childScan.TheScan.DissociationType == DissociationType.CID)
                {
                    continue;
                }
                var childFragments = GlycoPeptides.GlyGetTheoreticalFragments(GType, MS2ChildScanDissociationType, theScanBestPeptide, peptideWithMod, npos, glycans);

                n += childFragments.Where(v => v.ProductType != ProductType.D && v.ProductType != ProductType.Ycore && v.ProductType != ProductType.Y).Count();

                p += childScan.TheScan.MassSpectrum.Size * ProductMassTolerance.GetRange(1000).Width / childScan.TheScan.MassSpectrum.Range.Width;

            }



        }
            
    }
}
