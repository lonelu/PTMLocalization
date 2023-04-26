using Chemistry;
using MassSpectrometry;
using Nett;
using System;
using System.Collections.Generic;
using System.Linq;
using MzLibUtil;

namespace EngineLayer
{
    public class Ms2ScanWithSpecificMass
    {
        public Ms2ScanWithSpecificMass(MsDataScan mzLibScan, double precursorMonoisotopicPeakMz, int precursorCharge, string fullFilePath, double DeconvolutionMassTolerance, double DeconvolutionIntensityRatio, IsotopicEnvelope[] neutralExperimentalFragments = null)
        {
            PrecursorMonoisotopicPeakMz = precursorMonoisotopicPeakMz;
            PrecursorCharge = precursorCharge;
            PrecursorMass = PrecursorMonoisotopicPeakMz.ToMass(precursorCharge);
            FullFilePath = fullFilePath;
            ChildScans = new List<Ms2ScanWithSpecificMass>();
            NativeId = mzLibScan.NativeId;

            TheScan = mzLibScan;

            ExperimentalFragments = neutralExperimentalFragments ?? GetNeutralExperimentalFragments(mzLibScan, DeconvolutionMassTolerance, DeconvolutionIntensityRatio);
            
            if (ExperimentalFragments != null && ExperimentalFragments.Any())
            {
                DeconvolutedMonoisotopicMasses = ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray();
            }
            else
            {
                DeconvolutedMonoisotopicMasses = new double[0];
            }
            BasePeakIntensity = GetBasePeakIntensity();
        }

        /**
         * Compute oxonium intensity ratio (e.g. 138/144) for this scan. 
         */
        public double ComputeOxoRatio(double numeratorMass, double denominatorMass, Tolerance productTolerance)
        {
            GetClosestExperimentalFragmentMzWithinTol(numeratorMass, productTolerance, out double? numerator);
            GetClosestExperimentalFragmentMzWithinTol(denominatorMass, productTolerance, out double? denominator); 
            double ratio = (double)(numerator / denominator);
            return ratio;
        }

        /**
         * Determine if a given set of oxonium ions are found in the scan, within the provided PPM tolerance. 
         * For each provided FilterRule, checks that at least one required ion is found (true if so). 
         * Returns a dictionary of rule: oxoniums found or not 
         * minIntensity param is the base peak intensity pre-multiplied by minimum rel intensity user param
         */
        public Dictionary<FilterRule, bool> FindOxoniums(List<FilterRule> filterRules, Tolerance oxoTolerance, double minIntensity)
        {
            Dictionary<FilterRule, bool> outputBools = new();
            foreach (FilterRule rule in filterRules) 
            {
                outputBools[rule] = false;
                double totalIntensity = 0;
                foreach (double mz in rule.Masses)
                {
                    GetClosestExperimentalFragmentMzWithinTol(mz, oxoTolerance, out double? intensity);
                    if (intensity > 0)
                    {
                        totalIntensity += (double) intensity;
                    }
                }
                if (totalIntensity >= minIntensity)
                {
                    outputBools[rule] = true;
                }
            }
            return outputBools;
        }


        public MsDataScan TheScan { get; }
        public double PrecursorMonoisotopicPeakMz { get; }
        public double PrecursorMass { get; }
        public int PrecursorCharge { get; }
        public string FullFilePath { get; }
        public IsotopicEnvelope[] ExperimentalFragments { get; private set; }
        public List<Ms2ScanWithSpecificMass> ChildScans { get; set; } // MS2/MS3 scans that are children of this MS2 scan
        private double[] DeconvolutedMonoisotopicMasses;
        public string NativeId { get; } 

        public int OneBasedScanNumber => TheScan.OneBasedScanNumber;

        public int? OneBasedPrecursorScanNumber => TheScan.OneBasedPrecursorScanNumber;

        public double RetentionTime => TheScan.RetentionTime;

        public int NumPeaks => TheScan.MassSpectrum.Size;

        public double TotalIonCurrent => TheScan.TotalIonCurrent;

        public double BasePeakIntensity { get; set; }

        public static IsotopicEnvelope[] GetNeutralExperimentalFragments(MsDataScan scan, double DeconvolutionMassTolerance, double DeconvolutionIntensityRatio, bool AssumeOrphanPeaksAreZ1Fragments = true)
        {
            int minZ = 1;
            int maxZ = 10;

            var neutralExperimentalFragmentMasses = scan.MassSpectrum.Deconvolute(scan.MassSpectrum.Range,
                minZ, maxZ, DeconvolutionMassTolerance, DeconvolutionIntensityRatio).ToList();

            if (AssumeOrphanPeaksAreZ1Fragments)
            {
                HashSet<double> alreadyClaimedMzs = new HashSet<double>(neutralExperimentalFragmentMasses
                    .SelectMany(p => p.Peaks.Select(v => Chemistry.ClassExtensions.RoundedDouble(v.mz).Value)));

                for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
                {
                    double mz = scan.MassSpectrum.XArray[i];
                    double intensity = scan.MassSpectrum.YArray[i];

                    if (!alreadyClaimedMzs.Contains(Chemistry.ClassExtensions.RoundedDouble(mz).Value))
                    {
                        neutralExperimentalFragmentMasses.Add(new IsotopicEnvelope(
                            new List<(double mz, double intensity)> { (mz, intensity) },
                            mz.ToMass(1), 1, intensity, 0, 0));
                    }
                }
            }

            return neutralExperimentalFragmentMasses.OrderBy(p => p.MonoisotopicMass).ToArray();
        }

        public IsotopicEnvelope GetClosestExperimentalIsotopicEnvelope(double theoreticalNeutralMass)
        {
            if (DeconvolutedMonoisotopicMasses.Length == 0)
            {
                return null;
            }
            return ExperimentalFragments[GetClosestFragmentMass(theoreticalNeutralMass)];
        }

        public int GetClosestFragmentMass(double mass)
        {
            int index = Array.BinarySearch(DeconvolutedMonoisotopicMasses, mass);
            if (index >= 0)
            {
                return index;
            }
            index = ~index;

            if (index == DeconvolutedMonoisotopicMasses.Length)
            {
                return index - 1;
            }
            if (index == 0 || mass - DeconvolutedMonoisotopicMasses[index - 1] > DeconvolutedMonoisotopicMasses[index] - mass)
            {
                return index;
            }

            return index - 1;
        }

        //look for IsotopicEnvelopes which are in the range of acceptable mass 
        public IsotopicEnvelope[] GetClosestExperimentalIsotopicEnvelopeList(double minimumMass, double maxMass)
        {

            if (DeconvolutedMonoisotopicMasses.Length == 0)
            {
                return null;
            }

            //if no mass is in this range, then return null
            if (DeconvolutedMonoisotopicMasses[0] > maxMass || DeconvolutedMonoisotopicMasses.Last() < minimumMass)
            {
                return null;
            }

            int startIndex = GetClosestFragmentMass(minimumMass);
            int endIndex = GetClosestFragmentMass(maxMass);

            //the index we get from GetClosestFragmentMass is the closest mass, while the acceptable mass we need is between minimumMass and maxMass
            //so the startIndex mass is supposed to be larger than minimumMass, if not , then startIndex increases by 1;
            //the endIndex mass is supposed to be smaller than maxMass, if not , then endIndex decreases by 1;
            if (DeconvolutedMonoisotopicMasses[startIndex]<minimumMass)
            {
                startIndex = startIndex+1;
            }
            if(DeconvolutedMonoisotopicMasses[endIndex] > maxMass)
            {
                endIndex = endIndex - 1;
            }
            int length = endIndex - startIndex + 1;

            if (length < 1)
            {
                return null;
            }
            IsotopicEnvelope[] isotopicEnvelopes = ExperimentalFragments.Skip(startIndex).Take(length).ToArray();
            return isotopicEnvelopes;
        }

        /**
         * Like GetClosestExperimentalFragment but restricted to a given tolerance (ppm) around the theoretical m/z. Returns null
         * if no peaks found or if the closest peak is outside the tolerance. 
         * Returns intensity of 0 for peaks not found. 
         */
        public double? GetClosestExperimentalFragmentMzWithinTol(double theoreticalMz, Tolerance productTolerance, out double? intensity)
        {
            double mass = GetClosestExperimentalFragmentMz(theoreticalMz, out intensity);
            if (productTolerance.Within(mass, theoreticalMz)) 
            {
                return mass;
            }
            else
            {
                intensity = 0.0;
                return null;
            }
        }

        public double GetClosestExperimentalFragmentMz(double theoreticalMz, out double? intensity)
        {
            if (TheScan.MassSpectrum.XArray.Length == 0)
            {
                intensity = null;
                return 0;
            }
            intensity = TheScan.MassSpectrum.YArray[GetClosestFragmentMzIndex(theoreticalMz).Value];
            return TheScan.MassSpectrum.XArray[GetClosestFragmentMzIndex(theoreticalMz).Value];
        }

        private int? GetClosestFragmentMzIndex(double mz)
        {
            if (TheScan.MassSpectrum.XArray.Length == 0)
            {
                return null;
            }
            int index = Array.BinarySearch(TheScan.MassSpectrum.XArray, mz);
            if (index >= 0)
            {
                return index;
            }
            index = ~index;

            if (index >= TheScan.MassSpectrum.XArray.Length)
            {
                return index - 1;
            }
            if (index == 0)
            {
                return index;
            }

            if (mz - TheScan.MassSpectrum.XArray[index - 1] > TheScan.MassSpectrum.XArray[index] - mz)
            {
                return index;
            }
            return index - 1;

        }
        private double GetBasePeakIntensity()
        {
            if (TheScan.MassSpectrum.YArray.Length == 0)
            {
                return 0;
            }
            return TheScan.MassSpectrum.YArray.Max();
        }
    }
}