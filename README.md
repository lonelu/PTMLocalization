# PTMLocalization. 
The program could be used for O-Glycopeptide ETD/EThcD spectrum O-Glycan localization. 


# Glycan Databases

Users can add their own database for special purpose.

You can follow the exist glycan database files to create the new glycan database files. 

The database file should follow the following format:

    Composition based: HexNAc(1)Hex(1)Fuc(1). This is compitable with Byonic
    Structure based: (N(H)(F)). This is compitable with pGlyco2.

The program currently support the following monosaccharides. 

{"Hex"},{"HexNAc"},{"NeuAc"},{"NeuGc"},{"Fuc"},{"Phospho"},{"Sulfo"},{"Na"},{"Ac"},{"Xylose" }

| Name      | Symbol |
| ----------- | ----------- |
| Hex | H |
| HexNAc | N |
| NeuAc | A |
| NeuGc | G |
| Fuc | F |
| Phospho | P |
| Sulfo | S |
| Na | Y |
| Ac | C |
| Xylose | X |
| SuccinylHex | U |
| Formylation | N |

(If there are special requirement, please contact the developers. We have a new update in the near future to support the adding or any type of sugar.)

# Output

All potential glycan localizations: for unlocalized O-Glycans, show all potential the glycan on potential glycosites.
AllSiteSpecificLocalizationProbability: provide the site probability for every possible site.

Those are mainly for level 2 and level 3 ids. For example:

| Name      | Info |
| ----------- | ----------- |
| Peptide Base Sequence | STNASTVPFRNPDENSR |
| GlycanLocalizationLevel | Level2 |
| Localized Glycans with Peptide Site Specific Probability | [6,H1N1,0.994] |
| Localized Glycans with Protein Site Specific Probability | [223,H1N1,0.994] |
| All potential glycan localizations | {@10[2-1,6-1]}{@10[1-1,6-1]} |

* {@10[2-1,6-1]} corresponding to {GlycanGroupId 10, [Glycosite 2, glycanId 1], [Glycosite 6, glycanId 1]}

* AllSiteSpecificLocalizationProbability:{@1[1,0.500]}{@2[1,0.500]}{@5[1,0.006]}{@6[1,0.994]}{@16[1,0.000]}

* {@1[1,0.500]} corresponding to {Site 1, [glycanId 1, probability 0.5]}

For this peptide there are two H1N1 modified. We report H1N1 on S6 which could be localized with probability 0.994, the other glycan could be at either S1 or T2 based on the probabiliy and thus it is a leval 2 id.
