﻿using System.Collections.Generic;


namespace EngineLayer.GlycoSearch
{
    public class AdjNode
    {
        //AdjNode -> Adjactent node is used to build graph matrix for localizaiton. Each node in graph matrix contain Sources, max cost, current cost, etc.
        //For more about localization graph matrix, please check LocalizaitonGraph.cs. 
        public AdjNode(int i, int j, int modPos, ModBox box)
        {
            PosX = i;
            PosY = j;
            ModPos = modPos;
            ModBox = box;
        }

        public int PosX { get; }
        public int PosY { get; }

        public int ModPos { get; }
        public ModBox ModBox { get; }

        //sources are represented by index. Only track ones with highest cummulative cost
        public HashSet<int> CummulativeSources { get; set; } = new HashSet<int>();

        public double CummulativeCost { get; set; }

        public double CurrentCost { get; set; }

        public HashSet<int> AllSources { get; set; } = new HashSet<int>();
    }
}
