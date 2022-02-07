﻿using Chemistry;
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
        public static string DataDir { get; private set; }

        public static Monosaccharide[] Monosaccharides { get; private set; }

        public static List<string> OGlycanLocations { get; private set; }

        public static List<string> NGlycanLocations { get; private set; }

        public static void SetUpGlobalVariables()
        {
            Loaders.LoadElements();

            SetUpDataDirectory();

            LoadGlycans();
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

        private static void LoadGlycans()
        {
            string monosaccharidePath = Path.Combine(DataDir, @"Glycan_Mods", @"Monosaccharide.tsv");
            Monosaccharides = Monosaccharide.LoadMonosaccharide(monosaccharidePath);

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
