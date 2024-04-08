using Chemistry;
using Easy.Common.Extensions;
using MassSpectrometry;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using UsefulProteomicsDatabases;

namespace EngineLayer
{
    public static class GlobalVariables
    {
        // tranlsations for default FragPipe monosaccharides names
        public static Dictionary<string, char> fragpipeGlycsToSymbols = new Dictionary<string, char>
        {
            { "HexNAc", 'N' },
            { "Hex", 'H' },
            { "Fuc", 'F' },
            { "NeuAc", 'A' },
            { "NeuGc", 'G' },
            { "Phospho", 'P' },
            { "Sulfo", 'S' },
            { "Na", 'Y' },
            { "Acetyl", 'C' },
            { "Xylose", 'X' },
            { "Succinyl", 'U' },
            { "Formyl", 'M' },
        };


        public static string DataDir { get; private set; }

        public static Monosaccharide[] Monosaccharides { get; private set; }

        public static List<string> OGlycanLocations { get; private set; }

        public static List<string> NGlycanLocations { get; private set; }

        public static List<FilterRule> OxoniumFilters {get; private set; }

        public static HashSet<char> usedSymbols { get; private set; }

        public static void SetUpGlobalVariables(string? monosaccharidePath, string? oxoniumPath)
        {
            Loaders.LoadElements();

            SetUpDataDirectory();

            usedSymbols = new HashSet<char>();
            usedSymbols.AddRange(fragpipeGlycsToSymbols.Values);
            LoadGlycans(monosaccharidePath, oxoniumPath);
        }

        private static void SetUpDataDirectory()
        {
            // get data directory
            var pathToProgramFiles = Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles);
            if (!String.IsNullOrWhiteSpace(pathToProgramFiles) && AppDomain.CurrentDomain.BaseDirectory.Contains(pathToProgramFiles)
                && !AppDomain.CurrentDomain.BaseDirectory.Contains("Jenkins"))
            {
                DataDir = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "MetaMorpheus");
            }
            else
            {
                DataDir = AppDomain.CurrentDomain.BaseDirectory;
            }

        }

        private static void LoadGlycans(string? monosaccharidePath, string? oxoniumPath)
        {
            monosaccharidePath ??= Path.Combine(DataDir, @"Glycan_Mods", @"Monosaccharide.tsv");    // use default if not provided
            Monosaccharides = Monosaccharide.LoadMonosaccharide(monosaccharidePath);

            oxoniumPath ??= Path.Combine(DataDir, @"Glycan_Mods", @"OxoniumFilter.tsv");            // use default if not provided
            OxoniumFilters = FilterRule.LoadOxoniumFilters(oxoniumPath);

            OGlycanLocations = new List<string>();
            NGlycanLocations = new List<string>();

            foreach (var glycanFile in Directory.GetFiles(Path.Combine(DataDir, @"Glycan_Mods", @"OGlycan")))
            {
                OGlycanLocations.Add(glycanFile);
            }

            foreach (var glycanFile in Directory.GetFiles(Path.Combine(DataDir, @"Glycan_Mods", @"NGlycan")))
            {
                NGlycanLocations.Add(glycanFile);
            }
        }

    }
}
