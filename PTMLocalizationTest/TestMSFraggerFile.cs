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
using IO.MzML;

namespace PTMLocalizationTest
{
    [TestFixture]
    public class TestMSFraggerFile
    {
        [Test]
        public static void TestMSFraggerMZlib()
        {
            string origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "2019_09_16_StcEmix_35trig_EThcD25_rep1_4565_msfraggered.mzML");
            FilteringParams filter = new(200, 0.01, 1, null, false, false, true);
            bool lookForPrecursorScans = false;
            MsDataFile myFile = Mzml.LoadAllStaticData(origDataFile, filter, 1, lookForPrecursorScans);


        }
    }
}
