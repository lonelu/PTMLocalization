using CommandLine;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;


namespace CMD
{
    public class CmdSettings
    {

        [Option('v', Default = VerbosityType.normal, HelpText = "[Optional] Determines how much text is written. Options are no output ('none'), minimal output and errors  ('minimal'), or normal ('normal')")]
        public VerbosityType Verbosity { get; set; }
        public enum VerbosityType { none, minimal, normal };

        [Option('b', Default = 10, HelpText = "[Optional] productPpmTol ")]
        public double productPpmTol { get; set; }

        [Option('c', Default = 30, HelpText = "[Optional] precursorPpmTol ")]
        public double precursorPpmTol { get; set; }

        [Option('g', Default = "OGlycan.gdb", HelpText = "[Optional] glycoDatabase ")]
        public string glycoDatabase { get; set; }

        [Option('n', Default = 3, HelpText = "[Optional] maxNumGlycans ")]
        public int maxNumGlycans { get; set; }

        [Option('i', Default = 0, HelpText = "[Optional] minIsotopeError ")]
        public int minIsotopeError { get; set; }

        [Option('j', Default = 2, HelpText = "[Optional] maxIsotopeError ")]
        public int maxIsotopeError { get; set; }

        [Option('d', HelpText = "[Optional] rawfileDirectory")]
        public string rawfileDirectory { get; set; }

        [Option('s', HelpText = "[Optional] psmFile")]
        public string psmFile { get; set; }

        [Option('p', HelpText = "[Optional] scanpairFile. Searches raw file directory if not provided ")]
        public string scanpairFile { get; set; }

        [Option('o', HelpText = "[Optional] Output folder")]
        public string outputFolder { get; set; }

        public override string ToString()
        {
            string x = "";
            x += " -b " + productPpmTol.ToString();

            x += " -c " + precursorPpmTol.ToString();

            x += " -g " + glycoDatabase;
            x += " -n " + maxNumGlycans.ToString();
            x += " -i " + minIsotopeError.ToString(); 
            x += " -j " + maxIsotopeError.ToString();
            x += " -d " + rawfileDirectory;
            x += " -s " + psmFile;
            x += " -p " + scanpairFile;
            x += " -o " + outputFolder;

            return x;
        }
    }
}
