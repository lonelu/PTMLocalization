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
        public int run_msfragger(double productPpmTol, double PrecursorPpmTol, string psmFile, string scanpairFile, string rawfileDirectory, string glycoDatabase, int maxNumGlycans, int minIsotopeError, int maxIsotopeError)
        {
            Tolerance ProductMassTolerance = new PpmTolerance(productPpmTol);
            Tolerance PrecursorMassTolerance = new PpmTolerance(PrecursorPpmTol);
            if (psmFile == null)
            {
                Console.WriteLine("No PSM file specified, exiting.");
                return 2;
            }


            if (!File.Exists(glycoDatabase))
            {
                // read internal databases if not passed a full path
                try
                {
                    glycoDatabase = GlobalVariables.OGlycanLocations.Where(p => p.Contains(glycoDatabase)).First();
                }
                catch (Exception ex)
                {
                    // file not found - warn user and exit
                    Console.WriteLine("No Glycan Database found at {0}. Exiting.", glycoDatabase);
                    return 2;
                }

            }

            int[] isotopes = new int[maxIsotopeError - minIsotopeError + 1];
            for (int i = 0; i < isotopes.Length; i++)
            {
                isotopes[i] = minIsotopeError + i;
            }

            var localizer = new MSFragger_RunLocalization(psmFile, scanpairFile, rawfileDirectory, glycoDatabase, maxNumGlycans, PrecursorMassTolerance, ProductMassTolerance, isotopes);
            int returnCode = localizer.Localize();
            return returnCode;
        }

    }
}
