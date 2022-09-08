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

            GlobalVariables.SetUpGlobalVariables();

            //var task = new Task();
            //task.run_msfragger(settings.productPpmTol, settings.precursorPpmTol, settings.psmFile, settings.scanpairFile, settings.rawfileDirectory, settings.glycoDatabase, settings.maxNumGlycans);

            try
            {
                //var a = new PTMLocalization.RunLocalization();
                //a.test_run(settings.OutputFolder);
                var task = new Task();
                task.run_msfragger(settings.productPpmTol, settings.precursorPpmTol, settings.psmFile, settings.scanpairFile, settings.rawfileDirectory, settings.glycoDatabase, settings.maxNumGlycans);
                Console.WriteLine("Run finished.");
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
