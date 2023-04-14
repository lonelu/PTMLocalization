using Easy.Common.Extensions;
using EngineLayer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class FilterRule
    {
        public FilterRule(List<double> massList, Dictionary<byte, int> rules) 
        { 
            Masses = massList;
            MonosaccharidesRequired = rules;
        }

        public List<double> Masses { get; }

        public Dictionary<byte, int> MonosaccharidesRequired { get; }

        /**
        * Load a list of rules for oxonium ion-based filtering from a tsv file. 
        * File format: 
        * mass1,mass2,...\t ID1 \t minCount1 \t ID2 \t minCount2 (etc)
        */
        public static List<FilterRule> LoadOxoniumFilters(string filePath)
        {
            List<FilterRule> rules = new();
            using (StreamReader lines = new StreamReader(filePath))
            {
                bool firstLine = true;
                while (lines.Peek() != -1)
                {
                    string line = lines.ReadLine();

                    if (firstLine)
                    {
                        //Skip the first line
                        firstLine = false;
                        continue;
                    }

                    var splits = line.Replace('"', ' ').Split('\t');

                    // Parse mass(es)
                    var massSplits = splits[0].Split(',');
                    List<double> masses = new();
                    foreach (string mass in massSplits)
                    {
                        masses.Add(double.Parse(mass.Trim()));
                    }

                    // Parse monosaccharide rules
                    Dictionary<byte, int> ruleDict = new();
                    int i = 1;
                    while (i < splits.Length)
                    {
                        if (splits[i].IsNullOrEmpty())
                        {
                            i += 2;
                            continue;
                        }
                        byte id = byte.Parse(splits[i].Trim());
                        i++;
                        int count = int.Parse(splits[i].Trim());
                        i++;
                        ruleDict[id] = count;
                    }

                    FilterRule rule = new FilterRule(masses, ruleDict);
                    rules.Add(rule);
                }
            }
            return rules;
        }
    }
}
