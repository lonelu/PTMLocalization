using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System;
using System.Collections.Generic;
using System.Linq;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using MzLibUtil;
using EngineLayer.GlycoSearch;
using IO.Mgf;
using MassSpectrometry;
using EngineLayer;
using System.IO;

namespace PTMLocalization
{
    public class Task
    {
        public void run_msfragger(int productPpmTol, int PrecursorPpmTol, string psmFile, string scanpairFile, string rawfileDirectory, string glycoDatabase, int maxNumGlycans)
        {
            Tolerance ProductMassTolerance = new PpmTolerance(productPpmTol);
            Tolerance PrecursorMassTolerance = new PpmTolerance(PrecursorPpmTol);
            string _glycoDatabase = GlobalVariables.OGlycanLocations.Where(p => p.Contains(glycoDatabase)).First();

            //Tolerance ProductMassTolerance = new PpmTolerance(10);
            //Tolerance PrecursorMassTolerance = new PpmTolerance(30);
            //string psmFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\psm.tsv");
            //string scanpairFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\_scan-pairs.tsv");
            //string rawfileDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData");
            //string glycoDatabase = GlobalVariables.OGlycanLocations.Where(p => p.Contains("OGlycan.gdb")).First();
            //int maxNumGlycans = 3;
            int[] isotopes = { 0, 1, 2 };

            var localizer = new MSFragger_RunLocalization(psmFile, scanpairFile, rawfileDirectory, _glycoDatabase, maxNumGlycans, PrecursorMassTolerance, ProductMassTolerance, isotopes);
            localizer.Localize();
        }

    }
}
