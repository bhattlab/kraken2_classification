import os
import pandas as pd
from itertools import compress
import argparse



'''
    parser = argparse.ArgumentParser(description="Convert a kraken2 \
        standard report into the mpa-style  .")
    parser.add_argument('-i',
                        required=True, 
                        action='store',
                        dest='input_report',
                        help='Input kraken2 report.')
    parser.add_argument('-o',
                        required=True, 
                        action='store',
                        dest='output',
                        help='Output mpa style report')
    parser.add_argument('--no_corrections',
                        required=False, 
                        action='store_true',
                        dest='no_corrections',
                        help='Some taxonomy levels \
                        are bugged in the mpa report because \
                        this script doesnt use the full ncbi taxonomy. \
                        Normally these are corrected based on a database \
                        built by comparing reports. This option allows \
                        you to skip the corrections.')
    args = parser.parse_args()
'''
def main():

    taxonomy_levels = ["D","P","C","O","F","G","S"]
    taxonomy_lower = [t.lower() for t in taxonomy_levels]

    report = pd.read_csv(snakemake.input[0], delimiter='\t', header=None)

    # filter to rows at taxonomy_levels
    report_filtered = report.loc[report[3].isin(taxonomy_levels)]
    # filter taxonomy names
    report[5] = [a.strip() for a in report[5]]

    # tax_dict keeps track of where we are in the taxonomy
    tax_dict = {t:"" for t in taxonomy_levels}

    # build report strings
    mpa_strings = []
    for i in range(report.shape[0]):
        # reset tax dict appropriately
        # what an awful way to do this
        if report.iloc[i, 3][0] == "D":
            tax_dict.update({"P":"", "C":"", "O":"", "F":"", "G":"", "S":""})
        elif report.iloc[i, 3][0] == "P":
            tax_dict.update({"C":"", "O":"", "F":"", "G":"", "S":""})
        elif report.iloc[i, 3][0] == "C":
            tax_dict.update({"O":"", "F":"", "G":"", "S":""})
        elif report.iloc[i, 3][0] == "O":
            tax_dict.update({"F":"", "G":"", "S":""})
        elif report.iloc[i, 3][0] == "F":
            tax_dict.update({"G":"", "S":""})
        elif report.iloc[i, 3][0] == "G":
            tax_dict.update({"S":""})
        elif report.iloc[i, 3][0] == "S":
            tax_dict.update({})
        else:
            continue

        if report.iloc[i, 3] in taxonomy_levels:
            tax_dict[report.iloc[i, 3]] = report.iloc[i, 5]
            mpa_strings.append(tax_dict_to_string(tax_dict))

    # correct names if desired
    #if not args.no_corrections:
    mpa_strings = correct_mpa_strings(mpa_strings)

    # build data frame of report
    mpa_data = {"col1":mpa_strings, "col2":list(report_filtered[1])}  
    mpa_df = pd.DataFrame(mpa_data)

    # write out
    mpa_df.to_csv(snakemake.output[0], sep='\t', header=False, index=False)

def tax_dict_to_string(tax_dict):
    taxonomy_levels = ["D","P","C","O","F","G","S"]
    taxonomy_lower = [t.lower() for t in taxonomy_levels]
    level_list = [tax_dict[t] for t in taxonomy_levels if tax_dict[t] != ""]
    add_letters = list(compress(taxonomy_lower, [tax_dict[t] != "" for t in taxonomy_levels]))

    mpa_string = ""
    for letter, level in zip(add_letters, level_list):
        # print(letter)
        # print(level)
        if letter != "d": 
            mpa_string = mpa_string + "|" + letter + "__" + level
        else: 
            mpa_string = mpa_string + letter + "__" + level

    return(mpa_string)

def correct_mpa_strings(mpa_strings):
    # database of strings to correct
    # shouldn't store this in this script...
    correction_dict =  {
    "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Candidatus Hamiltonella|s__secondary endosymbiont of Ctenarytaina eucalypti" : "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|s__secondary endosymbiont of Ctenarytaina eucalypti",
    "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Candidatus Hamiltonella|s__secondary endosymbiont of Heteropsylla cubana" : "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|s__secondary endosymbiont of Heteropsylla cubana",
    "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Ruthia|s__Calyptogena okutanii thioautotrophic gill symbiont" : "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Calyptogena okutanii thioautotrophic gill symbiont",
    "d__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|o__Pelagibacterales|f__Pelagibacteraceae|g__Candidatus Fonsibacter" : "d__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|o__Pelagibacterales|g__Candidatus Fonsibacter",
    "d__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|o__Pelagibacterales|f__Pelagibacteraceae|g__Candidatus Fonsibacter|s__Candidatus Fonsibacter ubiquis" : "d__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|o__Pelagibacterales|g__Candidatus Fonsibacter|s__Candidatus Fonsibacter ubiquis",
    "d__Bacteria|p__Proteobacteria|c__Deltaproteobacteria|o__Desulfurellales|f__Candidatus Desulfofervidaceae" : "d__Bacteria|p__Proteobacteria|c__Deltaproteobacteria|f__Candidatus Desulfofervidaceae",
    "d__Bacteria|p__Proteobacteria|c__Deltaproteobacteria|o__Desulfurellales|f__Candidatus Desulfofervidaceae|g__Candidatus Desulfofervidus" : "d__Bacteria|p__Proteobacteria|c__Deltaproteobacteria|f__Candidatus Desulfofervidaceae|g__Candidatus Desulfofervidus",
    "d__Bacteria|p__Proteobacteria|c__Deltaproteobacteria|o__Desulfurellales|f__Candidatus Desulfofervidaceae|g__Candidatus Desulfofervidus|s__Candidatus Desulfofervidus auxilii" : "d__Bacteria|p__Proteobacteria|c__Deltaproteobacteria|f__Candidatus Desulfofervidaceae|g__Candidatus Desulfofervidus|s__Candidatus Desulfofervidus auxilii",
    "d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|g__Mogibacterium|s__[Eubacterium] minutum" : "d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|s__[Eubacterium] minutum",
    "d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|g__Mogibacterium|s__[Eubacterium] sulci" : "d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|s__[Eubacterium] sulci",
    "d__Bacteria|p__Cyanobacteria|c__Gloeobacteria|o__Chroococcidiopsidales" : "d__Bacteria|p__Cyanobacteria|o__Chroococcidiopsidales",
    "d__Bacteria|p__Cyanobacteria|c__Gloeobacteria|o__Chroococcidiopsidales|f__Chroococcidiopsidaceae" : "d__Bacteria|p__Cyanobacteria|o__Chroococcidiopsidales|f__Chroococcidiopsidaceae",
    "d__Bacteria|p__Cyanobacteria|c__Gloeobacteria|o__Chroococcidiopsidales|f__Chroococcidiopsidaceae|g__Chroococcidiopsis" : "d__Bacteria|p__Cyanobacteria|o__Chroococcidiopsidales|f__Chroococcidiopsidaceae|g__Chroococcidiopsis",
    "d__Bacteria|p__Cyanobacteria|c__Gloeobacteria|o__Chroococcidiopsidales|f__Chroococcidiopsidaceae|g__Chroococcidiopsis|s__Chroococcidiopsis thermalis" : "d__Bacteria|p__Cyanobacteria|o__Chroococcidiopsidales|f__Chroococcidiopsidaceae|g__Chroococcidiopsis|s__Chroococcidiopsis thermalis",
    "d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|o__Dehalococcoidales|f__Dehalococcoidaceae|g__Dehalogenimonas" : "d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|g__Dehalogenimonas",
    "d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|o__Dehalococcoidales|f__Dehalococcoidaceae|g__Dehalogenimonas|s__Dehalogenimonas formicexedens" : "d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|g__Dehalogenimonas|s__Dehalogenimonas formicexedens",
    "d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|o__Dehalococcoidales|f__Dehalococcoidaceae|g__Dehalogenimonas|s__Dehalogenimonas sp. WBC-2" : "d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|g__Dehalogenimonas|s__Dehalogenimonas sp. WBC-2",
    "d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|o__Dehalococcoidales|f__Dehalococcoidaceae|g__Dehalogenimonas|s__Dehalogenimonas lykanthroporepellens" : "d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|g__Dehalogenimonas|s__Dehalogenimonas lykanthroporepellens",
    "d__Bacteria|p__Bacteroidetes|c__Chitinophagia|o__Bacteroidetes Order II. Incertae sedis" : "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis",
    "d__Bacteria|p__Bacteroidetes|c__Chitinophagia|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae" : "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae",
    "d__Bacteria|p__Bacteroidetes|c__Chitinophagia|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salinibacter" : "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salinibacter",
    "d__Bacteria|p__Bacteroidetes|c__Chitinophagia|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salinibacter|s__Salinibacter ruber" : "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salinibacter|s__Salinibacter ruber",
    "d__Bacteria|p__Bacteroidetes|c__Chitinophagia|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|s__Rhodothermaceae bacterium" : "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|s__Rhodothermaceae bacterium",
    "d__Bacteria|p__Bacteroidetes|c__Chitinophagia|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|s__Rhodothermaceae bacterium RA" : "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|s__Rhodothermaceae bacterium RA",
    "d__Bacteria|p__Bacteroidetes|c__Chitinophagia|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Rhodothermus" : "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Rhodothermus",
    "d__Bacteria|p__Bacteroidetes|c__Chitinophagia|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Rhodothermus|s__Rhodothermus marinus" : "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Rhodothermus|s__Rhodothermus marinus",
    "d__Bacteria|p__Candidatus Saccharibacteria|g__Candidatus Saccharimonas|s__Candidatus Saccharibacteria oral taxon TM7x" : "d__Bacteria|p__Candidatus Saccharibacteria|s__Candidatus Saccharibacteria oral taxon TM7x",
    "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Candidatus Regiella|s__secondary endosymbiont of Ctenarytaina eucalypti": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|s__secondary endosymbiont of Ctenarytaina eucalypti",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Candidatus Regiella|s__secondary endosymbiont of Heteropsylla cubana": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|s__secondary endosymbiont of Heteropsylla cubana",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thiodiazotropha|s__Solemya velum gill symbiont": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Solemya velum gill symbiont",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thiodiazotropha|s__endosymbiont of Galathealinum brachiosum": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__endosymbiont of Galathealinum brachiosum",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thiodiazotropha|s__Solemya velesiana gill symbiont": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Solemya velesiana gill symbiont",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thiodiazotropha|s__Bathymodiolus thermophilus thioautotrophic gill symbiont": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Bathymodiolus thermophilus thioautotrophic gill symbiont",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thiodiazotropha|s__enodsymbiont of Escarpia spicata": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__enodsymbiont of Escarpia spicata",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thiodiazotropha|s__Solemya elarraichensis gill symbiont": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Solemya elarraichensis gill symbiont",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Ruthia|s__Bathymodiolus azoricus thioautotrophic gill symbiont": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Bathymodiolus azoricus thioautotrophic gill symbiont",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Ruthia|s__Solemya pervernicosa gill symbiont": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Solemya pervernicosa gill symbiont",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Ruthia|s__endosymbiont of Ridgeia piscesae": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__endosymbiont of Ridgeia piscesae",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Ruthia|s__Bathymodiolus septemdierum thioautotrophic gill symbiont": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Bathymodiolus septemdierum thioautotrophic gill symbiont",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Endoriftia|s__endosymbiont of Tevnia jerichonana": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__endosymbiont of Tevnia jerichonana",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Endoriftia|s__endosymbiont of Riftia pachyptila": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__endosymbiont of Riftia pachyptila",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Endoriftia|s__endosymbiont of Seepiophila jonesi": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__endosymbiont of Seepiophila jonesi",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Endoriftia|s__endosymbiont of Lamellibrachia luymesi": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__endosymbiont of Lamellibrachia luymesi",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Endoriftia|s__Gammaproteobacteria bacterium LUC13013_P4": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Gammaproteobacteria bacterium LUC13013_P4",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Endoriftia|s__Gammaproteobacteria bacterium LUC13017_P7": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Gammaproteobacteria bacterium LUC13017_P7",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Endoriftia|s__Gammaproteobacteria bacterium LUC14_002_19_P1": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Gammaproteobacteria bacterium LUC14_002_19_P1",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Endoriftia|s__Gammaproteobacteria bacterium LUC004_P11": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Gammaproteobacteria bacterium LUC004_P11",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Endoriftia|s__Gammaproteobacteria bacterium LUC13016_P6": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Gammaproteobacteria bacterium LUC13016_P6",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Endoriftia|s__Gammaproteobacteria bacterium LUC006_P13": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Gammaproteobacteria bacterium LUC006_P13",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Endoriftia|s__Gammaproteobacteria bacterium LUC13018_P8": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Gammaproteobacteria bacterium LUC13018_P8",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Endoriftia|s__Gammaproteobacteria bacterium LUC005_P12": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Gammaproteobacteria bacterium LUC005_P12",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Endoriftia|s__Gammaproteobacteria bacterium LUC003_P10": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Gammaproteobacteria bacterium LUC003_P10",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Endoriftia|s__Gammaproteobacteria bacterium LUC001_P9": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Gammaproteobacteria bacterium LUC001_P9",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Endoriftia|s__Gammaproteobacteria bacterium LUC13015_P5": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|s__Gammaproteobacteria bacterium LUC13015_P5",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Sedimenticola": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Sedimenticola",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Sedimenticola|s__Sedimenticola selenatireducens": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Sedimenticola|s__Sedimenticola selenatireducens",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Sedimenticola|s__Sedimenticola thiotaurini": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Sedimenticola|s__Sedimenticola thiotaurini",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Sedimenticola|s__Sedimenticola sp.": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Sedimenticola|s__Sedimenticola sp.",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Wohlfahrtiimonas": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Wohlfahrtiimonas",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Wohlfahrtiimonas|s__Wohlfahrtiimonas chitiniclastica": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Wohlfahrtiimonas|s__Wohlfahrtiimonas chitiniclastica",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Wohlfahrtiimonas|s__Wohlfahrtiimonas populi": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Wohlfahrtiimonas|s__Wohlfahrtiimonas populi",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Wohlfahrtiimonas|s__Wohlfahrtiimonas larvae": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Wohlfahrtiimonas|s__Wohlfahrtiimonas larvae",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Wohlfahrtiimonas|s__Wohlfahrtiimonas sp. G9077": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Wohlfahrtiimonas|s__Wohlfahrtiimonas sp. G9077",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Gallaecimonas": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Gallaecimonas",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Gallaecimonas|s__Gallaecimonas pentaromativorans": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Gallaecimonas|s__Gallaecimonas pentaromativorans",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Gallaecimonas|s__Gallaecimonas sp. HK-28": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Gallaecimonas|s__Gallaecimonas sp. HK-28",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Gallaecimonas|s__Gallaecimonas xiamenensis": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Gallaecimonas|s__Gallaecimonas xiamenensis",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Candidatus Thioglobus": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thioglobus",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Candidatus Thioglobus|s__Candidatus Thioglobus singularis": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thioglobus|s__Candidatus Thioglobus singularis",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Candidatus Thioglobus|s__Candidatus Thioglobus perditus": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thioglobus|s__Candidatus Thioglobus perditus",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Candidatus Thioglobus|s__Candidatus Thioglobus autotrophicus": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thioglobus|s__Candidatus Thioglobus autotrophicus",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Candidatus Thioglobus|s__Candidatus Thioglobus sp. MED-G25": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thioglobus|s__Candidatus Thioglobus sp. MED-G25",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Candidatus Thioglobus|s__Candidatus Thioglobus sp. REDSEA-S14_B12": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thioglobus|s__Candidatus Thioglobus sp. REDSEA-S14_B12",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Candidatus Thioglobus|s__uncultured SUP05 cluster bacterium": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thioglobus|s__uncultured SUP05 cluster bacterium",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Candidatus Thioglobus|s__Candidatus Thioglobus sp. REDSEA-S12_B1": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thioglobus|s__Candidatus Thioglobus sp. REDSEA-S12_B1",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Candidatus Thioglobus|s__Candidatus Thioglobus sp. REDSEA-S03_B1": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thioglobus|s__Candidatus Thioglobus sp. REDSEA-S03_B1",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Candidatus Thioglobus|s__Candidatus Thioglobus sp.": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thioglobus|s__Candidatus Thioglobus sp.",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Candidatus Thioglobus|s__Candidatus Thioglobus sp. TMED218": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thioglobus|s__Candidatus Thioglobus sp. TMED218",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Candidatus Thioglobus|s__Candidatus Thioglobus sp. MED-G23": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Candidatus Thioglobus|s__Candidatus Thioglobus sp. MED-G23",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Ignatzschineria": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Ignatzschineria",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Ignatzschineria|s__Ignatzschineria cameli": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Ignatzschineria|s__Ignatzschineria cameli",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Ignatzschineria|s__Ignatzschineria larvae": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Ignatzschineria|s__Ignatzschineria larvae",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Ignatzschineria|s__Ignatzschineria ureiclastica": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Ignatzschineria|s__Ignatzschineria ureiclastica",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Ignatzschineria|s__Ignatzschineria indica": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Ignatzschineria|s__Ignatzschineria indica",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Ignatzschineria|s__Ignatzschineria sp. F8392": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Ignatzschineria|s__Ignatzschineria sp. F8392",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Pseudohongiella": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Pseudohongiella",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Pseudohongiella|s__Pseudohongiella acticola": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Pseudohongiella|s__Pseudohongiella acticola",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Pseudohongiella|s__Pseudohongiella spirulinae": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Pseudohongiella|s__Pseudohongiella spirulinae",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Pseudohongiella|s__Pseudohongiella nitratireducens": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Pseudohongiella|s__Pseudohongiella nitratireducens",
"d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|f__Candidatus Competibacteraceae|g__Pseudohongiella|s__Pseudohongiella sp.": "d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|g__Pseudohongiella|s__Pseudohongiella sp.",
"d__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|o__Rhizobiales|f__Salinarimonadaceae|g__Salinarimonas|s__Salinarimonadaceae bacterium HL-109": "d__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|o__Rhizobiales|f__Salinarimonadaceae|s__Salinarimonadaceae bacterium HL-109",
"d__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|g__Candidatus Micropelagos|s__PS1 clade bacterium": "d__Bacteria|p__Proteobacteria|c__Alphaproteobacteria|s__PS1 clade bacterium",
"d__Bacteria|p__Proteobacteria|c__Zetaproteobacteria|o__Mariprofundales|f__Mariprofundaceae|g__Ghiorsea": "d__Bacteria|p__Proteobacteria|c__Zetaproteobacteria|g__Ghiorsea",
"d__Bacteria|p__Proteobacteria|c__Zetaproteobacteria|o__Mariprofundales|f__Mariprofundaceae|g__Ghiorsea|s__Ghiorsea bivora": "d__Bacteria|p__Proteobacteria|c__Zetaproteobacteria|g__Ghiorsea|s__Ghiorsea bivora",
"d__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Micrococcales|f__Microbacteriaceae|g__Rhodoluna|s__Microbacteriaceae bacterium MWH-Ta3": "d__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Micrococcales|f__Microbacteriaceae|s__Microbacteriaceae bacterium MWH-Ta3",
"d__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|g__Candidatus Carbobacillus|s__[Flavobacterium] thermophilum": "d__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|s__[Flavobacterium] thermophilum",
"d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|g__Monoglobus|s__[Bacteroides] pectinophilus": "d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|s__[Bacteroides] pectinophilus",
"d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|g__Anaerovorax|s__[Eubacterium] infirmum": "d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|s__[Eubacterium] infirmum",
"d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|g__Mobilibacterium|s__[Eubacterium] minutum": "d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|s__[Eubacterium] minutum",
"d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|g__Mobilibacterium|s__[Eubacterium] sulci": "d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|s__[Eubacterium] sulci",
"d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|g__Mobilibacterium|s__[Eubacterium] nodatum": "d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|s__[Eubacterium] nodatum",
"d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|g__Mobilibacterium|s__[Eubacterium] brachy": "d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|s__[Eubacterium] brachy",
"d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|g__Ileibacterium|s__[Eubacterium] saphenum": "d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiales Family XIII. Incertae Sedis|s__[Eubacterium] saphenum",
"d__Bacteria|p__Firmicutes|c__Tissierellia|o__Tissierellales|f__Tissierellaceae|g__Anaerosalibacter|s__[Clostridium] ultunense": "d__Bacteria|p__Firmicutes|c__Tissierellia|o__Tissierellales|f__Tissierellaceae|s__[Clostridium] ultunense",
"d__Bacteria|p__Firmicutes|c__Tissierellia|o__Tissierellales|f__Gottschalkiaceae|g__Sporanaerobacter": "d__Bacteria|p__Firmicutes|c__Tissierellia|o__Tissierellales|g__Sporanaerobacter",
"d__Bacteria|p__Firmicutes|c__Tissierellia|o__Tissierellales|f__Gottschalkiaceae|g__Sporanaerobacter|s__Sporanaerobacter acetigenes": "d__Bacteria|p__Firmicutes|c__Tissierellia|o__Tissierellales|g__Sporanaerobacter|s__Sporanaerobacter acetigenes",
"d__Bacteria|p__Firmicutes|c__Tissierellia|o__Tissierellales|f__Gottschalkiaceae|g__Sporanaerobacter|s__Sporanaerobacter sp. PP17-6a": "d__Bacteria|p__Firmicutes|c__Tissierellia|o__Tissierellales|g__Sporanaerobacter|s__Sporanaerobacter sp. PP17-6a",
"d__Bacteria|p__Firmicutes|c__Tissierellia|g__Sedimentibacter|s__[Bacteroides] coagulans": "d__Bacteria|p__Firmicutes|c__Tissierellia|s__[Bacteroides] coagulans",
"d__Bacteria|p__Cyanobacteria|c__Gloeobacteria|o__Gloeoemargaritales": "d__Bacteria|p__Cyanobacteria|o__Gloeoemargaritales",
"d__Bacteria|p__Cyanobacteria|c__Gloeobacteria|o__Gloeoemargaritales|f__Gloeomargaritaceae": "d__Bacteria|p__Cyanobacteria|o__Gloeoemargaritales|f__Gloeomargaritaceae",
"d__Bacteria|p__Cyanobacteria|c__Gloeobacteria|o__Gloeoemargaritales|f__Gloeomargaritaceae|g__Gloeomargarita": "d__Bacteria|p__Cyanobacteria|o__Gloeoemargaritales|f__Gloeomargaritaceae|g__Gloeomargarita",
"d__Bacteria|p__Cyanobacteria|c__Gloeobacteria|o__Gloeoemargaritales|f__Gloeomargaritaceae|g__Gloeomargarita|s__Gloeomargarita lithophora": "d__Bacteria|p__Cyanobacteria|o__Gloeoemargaritales|f__Gloeomargaritaceae|g__Gloeomargarita|s__Gloeomargarita lithophora",
"d__Bacteria|p__Candidatus Margulisbacteria|c__Candidatus Sericytochromatia": "d__Bacteria|c__Candidatus Sericytochromatia",
"d__Bacteria|p__Candidatus Margulisbacteria|c__Candidatus Sericytochromatia|s__Candidatus Sericytochromatia bacterium GL2-53 LSPB_72": "d__Bacteria|c__Candidatus Sericytochromatia|s__Candidatus Sericytochromatia bacterium GL2-53 LSPB_72",
"d__Bacteria|p__Candidatus Margulisbacteria|c__Candidatus Sericytochromatia|s__Candidatus Sericytochromatia bacterium S15B-MN24 CBMW_12": "d__Bacteria|c__Candidatus Sericytochromatia|s__Candidatus Sericytochromatia bacterium S15B-MN24 CBMW_12",
"d__Bacteria|p__Candidatus Margulisbacteria|c__Candidatus Sericytochromatia|s__Candidatus Sericytochromatia bacterium S15B-MN24 RAAC_196": "d__Bacteria|c__Candidatus Sericytochromatia|s__Candidatus Sericytochromatia bacterium S15B-MN24 RAAC_196",
"d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|o__Dehalococcoidales|g__Dehalogenimonas": "d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|g__Dehalogenimonas",
"d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|o__Dehalococcoidales|g__Dehalogenimonas|s__Dehalogenimonas sp. GP": "d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|g__Dehalogenimonas|s__Dehalogenimonas sp. GP",
"d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|o__Dehalococcoidales|g__Dehalogenimonas|s__Dehalogenimonas alkenigignens": "d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|g__Dehalogenimonas|s__Dehalogenimonas alkenigignens",
"d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|o__Dehalococcoidales|g__Dehalogenimonas|s__Dehalogenimonas formicexedens": "d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|g__Dehalogenimonas|s__Dehalogenimonas formicexedens",
"d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|o__Dehalococcoidales|g__Dehalogenimonas|s__Dehalogenimonas sp. WBC-2": "d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|g__Dehalogenimonas|s__Dehalogenimonas sp. WBC-2",
"d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|o__Dehalococcoidales|g__Dehalogenimonas|s__Dehalogenimonas lykanthroporepellens": "d__Bacteria|p__Chloroflexi|c__Dehalococcoidia|g__Dehalogenimonas|s__Dehalogenimonas lykanthroporepellens",
"d__Bacteria|p__Chloroflexi|c__Thermomicrobia|g__Thermorudis|s__Thermomicrobia bacterium": "d__Bacteria|p__Chloroflexi|c__Thermomicrobia|s__Thermomicrobia bacterium",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salinibacter": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salinibacter",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salinibacter|s__Salinibacter ruber": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salinibacter|s__Salinibacter ruber",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salinibacter|s__Salinibacter altiplanensis": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salinibacter|s__Salinibacter altiplanensis",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salinibacter|s__Salinibacter sp. 10B": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salinibacter|s__Salinibacter sp. 10B",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salinibacter|s__Salinibacter sp. J07SB67": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salinibacter|s__Salinibacter sp. J07SB67",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|s__Rhodothermaceae bacterium": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|s__Rhodothermaceae bacterium",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|s__Rhodothermaceae bacterium bin80": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|s__Rhodothermaceae bacterium bin80",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|s__Rhodothermaceae bacterium RA": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|s__Rhodothermaceae bacterium RA",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|s__Rhodothermaceae bacterium Tc-Br11_B2g6_7": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|s__Rhodothermaceae bacterium Tc-Br11_B2g6_7",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|s__Rhodothermaceae bacterium TMED105": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|s__Rhodothermaceae bacterium TMED105",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Rhodothermus": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Rhodothermus",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Rhodothermus|s__Rhodothermus marinus": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Rhodothermus|s__Rhodothermus marinus",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Rhodothermus|s__Rhodothermus profundi": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Rhodothermus|s__Rhodothermus profundi",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Longibacter": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Longibacter",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Longibacter|s__Longibacter salinarum": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Longibacter|s__Longibacter salinarum",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Longimonas": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Longimonas",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Longimonas|s__Longimonas halophila": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Longimonas|s__Longimonas halophila",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salisaeta": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salisaeta",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salisaeta|s__Salisaeta longa": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salisaeta|s__Salisaeta longa",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salisaeta|s__Bacteroidetes bacterium UBA7368": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|s__Bacteroidetes bacterium UBA7368",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salisaeta|s__Bacteroidetes bacterium UBA2364": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|s__Bacteroidetes bacterium UBA2364",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salisaeta|s__Bacteroidetes bacterium UBA2024": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|s__Bacteroidetes bacterium UBA2024",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salisaeta|s__Bacteroidetes bacterium UBA2412": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|s__Bacteroidetes bacterium UBA2412",
"d__Bacteria|p__Bacteroidetes|c__Saprospiria|o__Bacteroidetes Order II. Incertae sedis|f__Rhodothermaceae|g__Salisaeta|s__Bacteroidetes bacterium UBA5586": "d__Bacteria|p__Bacteroidetes|o__Bacteroidetes Order II. Incertae sedis|s__Bacteroidetes bacterium UBA5586",
"d__Bacteria|p__Candidatus Kryptonia|g__Candidatus Aegiribacteria": "d__Bacteria|g__Candidatus Aegiribacteria",
"d__Bacteria|p__Candidatus Kryptonia|g__Candidatus Aegiribacteria|s__Candidatus Aegiribacteria bacterium MLS_C": "d__Bacteria|g__Candidatus Aegiribacteria|s__Candidatus Aegiribacteria bacterium MLS_C",
"d__Bacteria|p__Verrucomicrobia|c__Spartobacteria|o__Chthoniobacterales|g__Terrimicrobium": "d__Bacteria|p__Verrucomicrobia|c__Spartobacteria|g__Terrimicrobium",
"d__Bacteria|p__Verrucomicrobia|c__Spartobacteria|o__Chthoniobacterales|g__Terrimicrobium|s__Terrimicrobium sacchariphilum": "d__Bacteria|p__Verrucomicrobia|c__Spartobacteria|g__Terrimicrobium|s__Terrimicrobium sacchariphilum",
"d__Bacteria|p__Verrucomicrobia|c__Spartobacteria|o__Chthoniobacterales|g__Candidatus Xiphinematobacter": "d__Bacteria|p__Verrucomicrobia|c__Spartobacteria|g__Candidatus Xiphinematobacter",
"d__Bacteria|p__Verrucomicrobia|c__Spartobacteria|o__Chthoniobacterales|g__Candidatus Xiphinematobacter|s__Candidatus Xiphinematobacter sp. Idaho Grape": "d__Bacteria|p__Verrucomicrobia|c__Spartobacteria|g__Candidatus Xiphinematobacter|s__Candidatus Xiphinematobacter sp. Idaho Grape",
"d__Bacteria|p__Candidatus Saccharibacteria|g__Candidatus Saccharimonas|s__candidate division TM7 genomosp. GTL1": "d__Bacteria|p__Candidatus Saccharibacteria|s__candidate division TM7 genomosp. GTL1",
"d__Bacteria|p__Candidatus Riflebacteria|c__Candidatus Ozemobacteria|f__Candidatus Ozemobacteraceae|g__Candidatus Ozemobacter|s__Candidatus Riflebacteria bacterium HGW-Riflebacteria-2": "d__Bacteria|p__Candidatus Riflebacteria|s__Candidatus Riflebacteria bacterium HGW-Riflebacteria-2",
"d__Bacteria|p__Candidatus Riflebacteria|c__Candidatus Ozemobacteria|f__Candidatus Ozemobacteraceae|g__Candidatus Ozemobacter|s__Candidatus Riflebacteria bacterium GWC2_50_8": "d__Bacteria|p__Candidatus Riflebacteria|s__Candidatus Riflebacteria bacterium GWC2_50_8",
"d__Bacteria|p__Candidatus Riflebacteria|c__Candidatus Ozemobacteria|f__Candidatus Ozemobacteraceae|g__Candidatus Ozemobacter|s__Candidatus Riflebacteria bacterium HGW-Riflebacteria-1": "d__Bacteria|p__Candidatus Riflebacteria|s__Candidatus Riflebacteria bacterium HGW-Riflebacteria-1",
"d__Bacteria|p__Candidatus Riflebacteria|c__Candidatus Ozemobacteria|f__Candidatus Ozemobacteraceae|g__Candidatus Ozemobacter|s__Candidatus Riflebacteria bacterium": "d__Bacteria|p__Candidatus Riflebacteria|s__Candidatus Riflebacteria bacterium",
"d__Bacteria|p__Candidatus Riflebacteria|c__Candidatus Ozemobacteria|f__Candidatus Ozemobacteraceae|g__Candidatus Ozemobacter|s__Candidatus Riflebacteria bacterium RBG_13_59_9": "d__Bacteria|p__Candidatus Riflebacteria|s__Candidatus Riflebacteria bacterium RBG_13_59_9",
"d__Archaea|p__Euryarchaeota|c__Methanomicrobia|o__Methanosarcinales|g__Candidatus Syntrophoarchaeum|s__Methanosarcinales archaeon UBA261": "d__Archaea|p__Euryarchaeota|c__Methanomicrobia|o__Methanosarcinales|s__Methanosarcinales archaeon UBA261",
"d__Archaea|p__Euryarchaeota|c__Methanomicrobia|o__Methanosarcinales|g__Candidatus Syntrophoarchaeum|s__Methanosarcinales archaeon UBA203": "d__Archaea|p__Euryarchaeota|c__Methanomicrobia|o__Methanosarcinales|s__Methanosarcinales archaeon UBA203",
"d__Archaea|p__Euryarchaeota|c__Methanomicrobia|o__Methanosarcinales|g__Candidatus Syntrophoarchaeum|s__Methanosarcinales archaeon Methan_03": "d__Archaea|p__Euryarchaeota|c__Methanomicrobia|o__Methanosarcinales|s__Methanosarcinales archaeon Methan_03",
"d__Archaea|p__Euryarchaeota|c__Methanomicrobia|o__Methanosarcinales|g__Candidatus Syntrophoarchaeum|s__Methanosarcinales archaeon Methan_02": "d__Archaea|p__Euryarchaeota|c__Methanomicrobia|o__Methanosarcinales|s__Methanosarcinales archaeon Methan_02",
"d__Archaea|p__Euryarchaeota|c__Methanomicrobia|o__Methanosarcinales|g__Candidatus Syntrophoarchaeum|s__Methanosarcinales archaeon Methan_01": "d__Archaea|p__Euryarchaeota|c__Methanomicrobia|o__Methanosarcinales|s__Methanosarcinales archaeon Methan_01",
"d__Archaea|p__Thaumarchaeota|c__Nitrososphaeria|o__Cenarchaeales": "d__Archaea|p__Thaumarchaeota|o__Cenarchaeales",
"d__Archaea|p__Thaumarchaeota|c__Nitrososphaeria|o__Cenarchaeales|f__Cenarchaeaceae": "d__Archaea|p__Thaumarchaeota|o__Cenarchaeales|f__Cenarchaeaceae",
"d__Archaea|p__Thaumarchaeota|c__Nitrososphaeria|o__Cenarchaeales|f__Cenarchaeaceae|g__Cenarchaeum": "d__Archaea|p__Thaumarchaeota|o__Cenarchaeales|f__Cenarchaeaceae|g__Cenarchaeum",
"d__Archaea|p__Thaumarchaeota|c__Nitrososphaeria|o__Cenarchaeales|f__Cenarchaeaceae|g__Cenarchaeum|s__Cenarchaeum symbiosum": "d__Archaea|p__Thaumarchaeota|o__Cenarchaeales|f__Cenarchaeaceae|g__Cenarchaeum|s__Cenarchaeum symbiosum",
"d__Archaea|p__Candidatus Korarchaeota|g__Candidatus Korarchaeum|s__Candidatus Korarchaeota archaeon NZ13-K": "d__Archaea|p__Candidatus Korarchaeota|s__Candidatus Korarchaeota archaeon NZ13-K",
    }

    mpa_strings_corrected = [correction_dict[s] if s in correction_dict.keys() else s for s in mpa_strings]
    return(mpa_strings_corrected)

if __name__ == '__main__':
    main()