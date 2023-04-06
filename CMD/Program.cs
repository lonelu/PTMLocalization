using CommandLine;
using CommandLine.Text;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using EngineLayer;
using PTMLocalization;

namespace CMD
{
    public static class Program
    {
        private static CmdSettings CommandLineSettings;

        public static int Main(string[] args)
        {
            int errorCode = 0;     

            var parser = new Parser(with => with.HelpWriter = null);
            var parserResult = parser.ParseArguments<CmdSettings>(args);

            parserResult
              .WithParsed<CmdSettings>(options => errorCode = Run(options))
              .WithNotParsed(errs => errorCode = DisplayHelp(parserResult, errs));

            return errorCode;
        }

        private static int Run(CmdSettings settings)
        {
            int errorCode = 0;

            if (settings.Verbosity == CmdSettings.VerbosityType.minimal || settings.Verbosity == CmdSettings.VerbosityType.normal)
            {
                Console.WriteLine("Welcome to PTM Localization.");
            }

            if (settings.Verbosity == CmdSettings.VerbosityType.minimal || settings.Verbosity == CmdSettings.VerbosityType.normal)
            {
                Console.WriteLine("PTMLocalization version: 1.0.");
            }

            // check settings and trim input paths
            bool pathsOK = CheckParamsAndPaths(settings);
            if (!pathsOK)
            {
                Console.WriteLine("Please fix the provided file path(s) and try again");
                errorCode = 1;
                return errorCode;
            }

            GlobalVariables.SetUpGlobalVariables();

            try
            {
                var task = new Task();
                errorCode = task.run_msfragger(settings.productPpmTol, settings.precursorPpmTol, settings.psmFile, settings.scanpairFile, settings.rawfileDirectory, settings.lcmsFilesList, settings.glycoDatabase, settings.maxNumGlycans, settings.minIsotopeError, settings.maxIsotopeError, settings.oxoFilter);
                if (errorCode == 0)
                {
                    Console.WriteLine("Run finished.");
                }
            }
            catch (Exception e)
            {
                while (e.InnerException != null)
                {
                    e = e.InnerException;
                }

                var message = "Run failed, Exception: " + e.Message;

                if (settings.Verbosity == CmdSettings.VerbosityType.minimal || settings.Verbosity == CmdSettings.VerbosityType.normal)
                {
                    Console.WriteLine(message);
                    Console.WriteLine(e.StackTrace);
                    Console.WriteLine(settings.ToString());
                }
                errorCode = 4;
            }

            return errorCode;
        }

        private static bool CheckParamsAndPaths(CmdSettings settings)
        {
            // trim file/directory paths (if provided) and confirm all specified files exist
            if (settings.psmFile != null)
            {
                settings.psmFile = settings.psmFile.Trim();
                if (!File.Exists(settings.psmFile))
                {
                    Console.WriteLine("Error: PSM file not found at {0}", settings.psmFile);
                    return false;
                }
            }
            if (settings.scanpairFile != null)
            {
                settings.scanpairFile = settings.scanpairFile.Trim();
                if (!File.Exists(settings.scanpairFile))
                {
                    Console.WriteLine("Error: Scan pair file not found at {0}", settings.scanpairFile);
                    return false;
                }
            }

            // raw file checking 
            bool rawfilesFound = false;
            if (settings.lcmsFilesList != null)
            {
                settings.lcmsFilesList = settings.lcmsFilesList.Trim();
                if (!File.Exists(settings.lcmsFilesList))
                {
                    Console.WriteLine("Error: lcms file not found at {0}", settings.lcmsFilesList);
                } 
                else
                {
                    rawfilesFound = true;
                }
            }
            if (settings.rawfileDirectory != null)
            {
                settings.rawfileDirectory = settings.rawfileDirectory.Trim();
                if (!Directory.Exists(settings.rawfileDirectory))
                {
                    if (!rawfilesFound)
                    {
                        // warn if an lcms file list has not been provided
                        Console.WriteLine("Error: rawfile directory not found at {0}", settings.rawfileDirectory);
                    }
                }
                else
                {
                    rawfilesFound = true;
                }
            }
            if (!rawfilesFound)
            {
                return false;
            }

            // glyco database checking
            if (settings.glycoDatabase != null)
            {
                settings.glycoDatabase = settings.glycoDatabase.Trim();
            }
            return true;
        }

        public static int DisplayHelp<T>(ParserResult<T> result, IEnumerable<Error> errs)
        {
            Console.WriteLine("Welcome to PTMLocation.");

            int errorCode = 0;

            var helpText = HelpText.AutoBuild(result, h =>
            {
                h.AdditionalNewLineAfterOption = false;
                h.Copyright = "";
                return HelpText.DefaultParsingErrorsHandler(result, h);
            }, e => e);

            helpText.MaximumDisplayWidth = 300;

            helpText.AddPostOptionsLine("Example usage (Windows): ");
            //helpText.AddPostOptionsLine("CMD.exe -d C:\\ExampleDatabase.fasta -s C:\\ExampleSpectra.mzML -t C:\\ExampleTask.toml");
            helpText.AddPostOptionsLine(Environment.NewLine);

            helpText.AddPostOptionsLine("Example usage (Linux): ");
            //helpText.AddPostOptionsLine("dotnet CMD.dll -d home/mydata/ExampleDatabase.fasta -s home/mydata/ExampleSpectra.mzML -t home/mydata/ExampleTask.toml");
            helpText.AddPostOptionsLine(Environment.NewLine);

            Console.WriteLine(helpText);

            if (errs.Any(x => x.Tag != ErrorType.HelpRequestedError))
            {
                errorCode = 1;
            }

            return errorCode;
        }

    }
}
