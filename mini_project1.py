def determine_data_type(value):
    """
    The function takes a string input and determines its data type to be either a float, int, or string. 
    """
    # BEGIN SOLUTION
    '''try:
        if float(value)==int(float(value)):
            return int
        else:
            return float
    except ValueError:
        return str'''
    try:
        if int(value) == float(value):
            return int
        else:
            return float
    except ValueError:
        try:
            float(value)
            if 'e' in value or '.' in value:
                return float
            else:
                return str
        except ValueError:
            return str
    # END SOLUTION


def determine_data_type_of_list(values):
    """
    Write a function whose input is a list of strings. 
    This function determines the correct data type of all the elements in the list. 
    For example, ['1', '2', '3'] is int, ['1.1', '2.2', '3.3'] is float, ['1.1', '2', '3.3'] 
    is also float, and ['1.1', '234String', '3.3'] is str. 
    The purpose of this function to figure out what to cast an array of strings to. 
    Some lists might be all ints, in which case the data type is int. 
    Some might be a mixture of ints and floats, in which case the data type will be a float. 
    Some lists might be a mixture of ints, floats, and strings, 
    in which case the data type of the list will be a string.
    NOTE: This function should use "determine_data_type" function you coded previously

    """
    # BEGIN SOLUTION
    c1=0
    c2=0
    c3=0
    for i in values:
        check_type=determine_data_type(i)
        if check_type==int:
            c1+= 1
        elif check_type==float:
            c2+= 1
        else:
            c3+= 1
    if c3>0:
        return str
    elif c2>0:
        return float
    else:
        return int
    # END SOLUTION


def format_sample_fields(format_field, sample_field):
    """
    Write a function whose inputs are a format field and sample field. 
    The format field looks like format_field = 'GT:AD:DP:GQ:PGT:PID:PL' and 
    the sample field looks like

    sample_field = {'XG102': '1/1:0,76:76:99:1|1:48306945_C_G:3353,229,0',
                'XG103': '1/1:0,52:52:99:.:.:1517,156,0',
                'XG104': '0/1:34,38:72:99:.:.:938,0,796',
                'XG202': '1/1:0,76:76:99:1|1:48306945_C_G:3353,229,0',
                'XG203': '1/1:0,52:52:99:.:.:1517,156,0',
                'XG204': '0/1:34,38:72:99:.:.:938,0,796',
                'XG302': '1/1:0,76:76:99:1|1:48306945_C_G:3353,229,0',
                'XG303': '1/1:0,52:52:99:.:.:1517,156,0',
                'XG304': '0/1:34,38:72:99:.:.:938,0,796',
                'XG402': '1/1:0,76:76:99:1|1:48306945_C_G:3353,229,0',
                'XG403': '1/1:0,52:52:99:.:.:1517,156,0',
                'XG404': '0/1:34,38:72:99:.:.:938,0,796'}

    Transform the inputs such that the output looks like this:

    output = {
        'XG102': {'AD': '0,76',
            'DP': '76',
            'GQ': '99',
            'GT': '1/1',
            'PGT': '1|1',
            'PID': '48306945_C_G',
            'PL': '3353,229,0'},
        'XG103': {'AD': '0,52',
                'DP': '52',
                'GQ': '99',
                'GT': '1/1',
                'PGT': '.',
                'PID': '.',
                'PL': '1517,156,0'},
        'XG104': {'AD': '34,38',
                'DP': '72',
                'GQ': '99',
                'GT': '0/1',
                'PGT': '.',
                'PID': '.',
                'PL': '938,0,796'},
        'XG202': {'AD': '0,76',
                'DP': '76',
                'GQ': '99',
                'GT': '1/1',
                'PGT': '1|1',
                'PID': '48306945_C_G',
                'PL': '3353,229,0'},
        'XG203': {'AD': '0,52',
                'DP': '52',
                'GQ': '99',
                'GT': '1/1',
                'PGT': '.',
                'PID': '.',
                'PL': '1517,156,0'},
        'XG204': {'AD': '34,38',
                'DP': '72',
                'GQ': '99',
                'GT': '0/1',
                'PGT': '.',
                'PID': '.',
                'PL': '938,0,796'},
        'XG302': {'AD': '0,76',
                'DP': '76',
                'GQ': '99',
                'GT': '1/1',
                'PGT': '1|1',
                'PID': '48306945_C_G',
                'PL': '3353,229,0'},
        'XG303': {'AD': '0,52',
                'DP': '52',
                'GQ': '99',
                'GT': '1/1',
                'PGT': '.',
                'PID': '.',
                'PL': '1517,156,0'},
        'XG304': {'AD': '34,38',
                'DP': '72',
                'GQ': '99',
                'GT': '0/1',
                'PGT': '.',
                'PID': '.',
                'PL': '938,0,796'},
        'XG402': {'AD': '0,76',
                'DP': '76',
                'GQ': '99',
                'GT': '1/1',
                'PGT': '1|1',
                'PID': '48306945_C_G',
                'PL': '3353,229,0'},
        'XG403': {'AD': '0,52',
                'DP': '52',
                'GQ': '99',
                'GT': '1/1',
                'PGT': '.',
                'PID': '.',
                'PL': '1517,156,0'},
        'XG404': {'AD': '34,38',
                'DP': '72',
                'GQ': '99',
                'GT': '0/1',
                'PGT': '.',
                'PID': '.',
                'PL': '938,0,796'}}
    """

    # BEGIN SOLUTION
    res1=format_field.split(':')
    l=[]
    m=[]
    for i in sample_field.keys():
        l.append(i)
    for i in sample_field.values():
        m.append(i)
    q=[]
    for i in range(len(m)):
        q.append(m[i].split(':'))
    p=[]
    for i in range(len(q)):
        var_dict = dict(zip(res1, q[i]))
        p.append(var_dict)
    d={}
    for i in range(len(l)):
        d[l[i]]=p[i]
    return d

    # END SOLUTION


def create_dict_from_line(header, line):
    """
    Given the header and a single line, transform them into dictionary as described above. 
    Header and line input are provided in this cell. 

    Write a function whose inputs are a list containing the vcf header and a variant line. 
    The function should return a dictionary using the header as keys and the variant line as values.
     The function should use the format_sample_fields you wrote previously to format the sample fields. 
     The output of the first line looks like this:

    {'ALT': 'G',
    'CHROM': '4',
    'FILTER': 'PASS',
    'ID': '.',
    'INFO': 'AC=1;AF=0.167;AN=6;BaseQRankSum=-2.542;ClippingRankSum=0;DP=180;ExcessHet=3.0103;FS=0;MLEAC=1;MLEAF=0.167;MQ=52.77;MQRankSum=-4.631;QD=0.39;ReadPosRankSum=1.45;SOR=0.758;VQSLOD=-8.209;culprit=MQ;ANNOVAR_DATE=2018-04-16;Func.refGene=intergenic;Gene.refGene=IL2,IL21;GeneDetail.refGene=dist=38536,dist=117597;ExonicFunc.refGene=.;AAChange.refGene=.;Func.ensGene=intergenic;Gene.ensGene=ENSG00000109471,ENSG00000138684;GeneDetail.ensGene=dist=38306,dist=117597;ExonicFunc.ensGene=.;AAChange.ensGene=.;cytoBand=4q27;gwasCatalog=.;tfbsConsSites=.;wgRna=.;targetScanS=.;Gene_symbol=.;OXPHOS_Complex=.;Ensembl_Gene_ID=.;Ensembl_Protein_ID=.;Uniprot_Name=.;Uniprot_ID=.;NCBI_Gene_ID=.;NCBI_Protein_ID=.;Gene_pos=.;AA_pos=.;AA_sub=.;Codon_sub=.;dbSNP_ID=.;PhyloP_46V=.;PhastCons_46V=.;PhyloP_100V=.;PhastCons_100V=.;SiteVar=.;PolyPhen2_prediction=.;PolyPhen2_score=.;SIFT_prediction=.;SIFT_score=.;FatHmm_prediction=.;FatHmm_score=.;PROVEAN_prediction=.;PROVEAN_score=.;MutAss_prediction=.;MutAss_score=.;EFIN_Swiss_Prot_Score=.;EFIN_Swiss_Prot_Prediction=.;EFIN_HumDiv_Score=.;EFIN_HumDiv_Prediction=.;CADD_score=.;CADD_Phred_score=.;CADD_prediction=.;Carol_prediction=.;Carol_score=.;Condel_score=.;Condel_pred=.;COVEC_WMV=.;COVEC_WMV_prediction=.;PolyPhen2_score_transf=.;PolyPhen2_pred_transf=.;SIFT_score_transf=.;SIFT_pred_transf=.;MutAss_score_transf=.;MutAss_pred_transf=.;Perc_coevo_Sites=.;Mean_MI_score=.;COSMIC_ID=.;Tumor_site=.;Examined_samples=.;Mutation_frequency=.;US=.;Status=.;Associated_disease=.;Presence_in_TD=.;Class_predicted=.;Prob_N=.;Prob_P=.;SIFT_score=.;SIFT_converted_rankscore=.;SIFT_pred=.;Polyphen2_HDIV_score=.;Polyphen2_HDIV_rankscore=.;Polyphen2_HDIV_pred=.;Polyphen2_HVAR_score=.;Polyphen2_HVAR_rankscore=.;Polyphen2_HVAR_pred=.;LRT_score=.;LRT_converted_rankscore=.;LRT_pred=.;MutationTaster_score=.;MutationTaster_converted_rankscore=.;MutationTaster_pred=.;MutationAssessor_score=.;MutationAssessor_score_rankscore=.;MutationAssessor_pred=.;FATHMM_score=.;FATHMM_converted_rankscore=.;FATHMM_pred=.;PROVEAN_score=.;PROVEAN_converted_rankscore=.;PROVEAN_pred=.;VEST3_score=.;VEST3_rankscore=.;MetaSVM_score=.;MetaSVM_rankscore=.;MetaSVM_pred=.;MetaLR_score=.;MetaLR_rankscore=.;MetaLR_pred=.;M-CAP_score=.;M-CAP_rankscore=.;M-CAP_pred=.;CADD_raw=.;CADD_raw_rankscore=.;CADD_phred=.;DANN_score=.;DANN_rankscore=.;fathmm-MKL_coding_score=.;fathmm-MKL_coding_rankscore=.;fathmm-MKL_coding_pred=.;Eigen_coding_or_noncoding=.;Eigen-raw=.;Eigen-PC-raw=.;GenoCanyon_score=.;GenoCanyon_score_rankscore=.;integrated_fitCons_score=.;integrated_fitCons_score_rankscore=.;integrated_confidence_value=.;GERP++_RS=.;GERP++_RS_rankscore=.;phyloP100way_vertebrate=.;phyloP100way_vertebrate_rankscore=.;phyloP20way_mammalian=.;phyloP20way_mammalian_rankscore=.;phastCons100way_vertebrate=.;phastCons100way_vertebrate_rankscore=.;phastCons20way_mammalian=.;phastCons20way_mammalian_rankscore=.;SiPhy_29way_logOdds=.;SiPhy_29way_logOdds_rankscore=.;Interpro_domain=.;GTEx_V6_gene=.;GTEx_V6_tissue=.;esp6500siv2_all=.;esp6500siv2_aa=.;esp6500siv2_ea=.;ExAC_ALL=.;ExAC_AFR=.;ExAC_AMR=.;ExAC_EAS=.;ExAC_FIN=.;ExAC_NFE=.;ExAC_OTH=.;ExAC_SAS=.;ExAC_nontcga_ALL=.;ExAC_nontcga_AFR=.;ExAC_nontcga_AMR=.;ExAC_nontcga_EAS=.;ExAC_nontcga_FIN=.;ExAC_nontcga_NFE=.;ExAC_nontcga_OTH=.;ExAC_nontcga_SAS=.;ExAC_nonpsych_ALL=.;ExAC_nonpsych_AFR=.;ExAC_nonpsych_AMR=.;ExAC_nonpsych_EAS=.;ExAC_nonpsych_FIN=.;ExAC_nonpsych_NFE=.;ExAC_nonpsych_OTH=.;ExAC_nonpsych_SAS=.;1000g2015aug_all=.;1000g2015aug_afr=.;1000g2015aug_amr=.;1000g2015aug_eur=.;1000g2015aug_sas=.;CLNALLELEID=.;CLNDN=.;CLNDISDB=.;CLNREVSTAT=.;CLNSIG=.;dbscSNV_ADA_SCORE=.;dbscSNV_RF_SCORE=.;snp138NonFlagged=.;avsnp150=.;CADD13_RawScore=0.015973;CADD13_PHRED=2.741;Eigen=-0.3239;REVEL=.;MCAP=.;Interpro_domain=.;ICGC_Id=.;ICGC_Occurrence=.;gnomAD_genome_ALL=0.0003;gnomAD_genome_AFR=0.0001;gnomAD_genome_AMR=0;gnomAD_genome_ASJ=0;gnomAD_genome_EAS=0.0007;gnomAD_genome_FIN=0.0009;gnomAD_genome_NFE=0.0002;gnomAD_genome_OTH=0.0011;gerp++gt2=.;cosmic70=.;InterVar_automated=.;PVS1=.;PS1=.;PS2=.;PS3=.;PS4=.;PM1=.;PM2=.;PM3=.;PM4=.;PM5=.;PM6=.;PP1=.;PP2=.;PP3=.;PP4=.;PP5=.;BA1=.;BS1=.;BS2=.;BS3=.;BS4=.;BP1=.;BP2=.;BP3=.;BP4=.;BP5=.;BP6=.;BP7=.;Kaviar_AF=.;Kaviar_AC=.;Kaviar_AN=.;ALLELE_END',
    'POS': '123416186',
    'QUAL': '23.25',
    'REF': 'A',
    'SAMPLE': {'XG102': {'AD': '51,8',
                      'DP': '59',
                      'GQ': '32',
                      'GT': '0/1',
                      'PL': '32,0,1388'},
            'XG103': {'AD': '47,0',
                      'DP': '47',
                      'GQ': '99',
                      'GT': '0/0',
                      'PL': '0,114,1353'},
            'XG104': {'AD': '74,0',
                      'DP': '74',
                      'GQ': '51',
                      'GT': '0/0',
                      'PL': '0,51,1827'},
            'XG202': {'AD': '51,8',
                      'DP': '59',
                      'GQ': '32',
                      'GT': '0/1',
                      'PL': '32,0,1388'},
            'XG203': {'AD': '47,0',
                      'DP': '47',
                      'GQ': '99',
                      'GT': '0/0',
                      'PL': '0,114,1353'},
            'XG204': {'AD': '74,0',
                      'DP': '74',
                      'GQ': '51',
                      'GT': '0/0',
                      'PL': '0,51,1827'},
            'XG302': {'AD': '51,8',
                      'DP': '59',
                      'GQ': '32',
                      'GT': '0/1',
                      'PL': '32,0,1388'},
            'XG303': {'AD': '47,0',
                      'DP': '47',
                      'GQ': '99',
                      'GT': '0/0',
                      'PL': '0,114,1353'},
            'XG304': {'AD': '74,0',
                      'DP': '74',
                      'GQ': '51',
                      'GT': '0/0',
                      'PL': '0,51,1827'},
            'XG402': {'AD': '51,8',
                      'DP': '59',
                      'GQ': '32',
                      'GT': '0/1',
                      'PL': '32,0,1388'},
            'XG403': {'AD': '47,0',
                      'DP': '47',
                      'GQ': '99',
                      'GT': '0/0',
                      'PL': '0,114,1353'},
            'XG404': {'AD': '74,0',
                      'DP': '74',
                      'GQ': '51',
                      'GT': '0/0',
                      'PL': '0,51,1827'}}}
    """
    # BEGIN SOLUTION
    d=dict(zip(header,line.split('\t')))
    cropdict= dict(list(d.items())[-12:])
    cropdict1= dict(list(d.items())[0:8])
    cropdict1.update({'SAMPLE': format_sample_fields(d["FORMAT"], cropdict)})
    return cropdict1
     #print(cropdict)
    #format_sample_fields(d['FORMAT'], cropdict)
    '''ff='GT:AD:DP:GQ:PL'
    res1=ff.split(':')
    #print(res1)
    l=[]
    m=[]
    for i in cropdict.keys():
        l.append(i)
    for i in cropdict.values():
        m.append(i)
    #print(l)
    #print(m)
    q=[]
    for i in range(len(m)):
        q.append(m[i].split(':'))
    #print(q)

    #print(new_list)
    p=[]
    for i in range(len(q)):
        var_dict = dict(zip(res1, q[i]))
        p.append(var_dict)
    #print(p)
    d={}
    for i in range(len(l)):
        d[l[i]]=p[i]
    #print(d)'''

    # END SOLUTION


def read_vcf_file(filename):
    """
    Write a function whose input is a filename for a vcf. 
    The function reads the vcf file one variant at a time and transforms it 
    into a dictionary using the create_dict_from_line function. 
    It returns a list containing all the variant dictionaries. 
    NOTE: Your function should be able to handle multiple lines.
    """
    # BEGIN SOLUTION
    header = None
    res= []
    with open(filename, encoding="utf-8") as file:
        for line in file:
            if not line.strip():
                continue
            if line[0]=='#' and line[1]=='#':
                continue
            if not header:
                header =line[1:].strip().split('\t')
                continue
            #resheader= '\t'.join([str(elem) for elem in header])
            line =line.strip().split('\t')
            resline= '\t'.join([str(elem) for elem in line])
            res.append(create_dict_from_line(header, resline))
    return res
    # END SOLUTION


def extract_info_field(data):
    """
    Write a function that extracts the info field from the data dictionary that was 
    created in the previous part. The function should return all the info field dictionaries as list. 
    """
    # BEGIN SOLUTION
    l=[]
    for i in data:
        if i.get('INFO'):
            l.append(i.get('INFO'))
    return l
    # END SOLUTION


def create_dictionary_of_info_field_values(data):
    """
    You now need to figure out that data types for each of the info fields. Begin by writing a function that first takes the info fields and turns them into a dictionary. Make sure to skip any fields that do not have a value or are missing a value.

    Note: only return keys that have values! 
    """

    # BEGIN SOLUTION
    from collections import defaultdict
    d =defaultdict(list)
    for ele in data:
        p=ele.split(';')
        #print(p)
        for i in p:
            if '=' not in i:
                continue
            key,value = i.split('=',1)
            if value=='.':
                continue
            if value not in d[key]:
                d[key].append(value)
    return d
    '''from collections import defaultdict
    result = defaultdict(list)
    for key, value in d.items():
        result[key].append(value)
    for key, value in d1.items():
        result[key].append(value)
    z=dict(result)
    result = {}
    new_dict = {}
    for key, value in d.items():
        unique_values = []
        for v in value:
            if v not in unique_values:
                unique_values.append(v)
        new_dict[key] = unique_values
    return new_dict
    # END SOLUTION'''




def determine_data_type_of_info_fields(data):
    """
    Write a function whose input is the output from create_dictionary_of_info_field_values 
    and uses the previously written function determine_data_type_of_list to determine 
    the data type of each of the info fields. The output is a dictionary whose 
    keys are the name of the info fields and values are the data type. 
    """
    # BEGIN SOLUTION
    '''d={}
    for k,v in data.items():
        c1=c2=c3=0
        for j in v:
            check_type=determine_data_type(j)
            if check_type==int:
                c1+= 1
            elif check_type==float:
                c2+= 1
            else:
                c3+= 1
        if c3>0:
            d[k]=str
        elif c2>0:
            d[k]=float
        else:
            d[k]=int
    return d'''
    d={}
    for k,v in data.items():
        check_type=determine_data_type_of_list(v)
        d[k]=check_type
    return d
    # END SOLUTION


def format_data(data, info_field_data_type):
    """
    Write a function whose first input is the data from read_vcf_file and 
    the second input is the output from determine_data_type_of_info_fields. 
    The function converts the info field into a dictionary and uses the data types 
    that you determined to cast each field into the correct data type. 
    Make sure to convert the POS to int and QUAL to float. 
    The output will look something like this (I have removed most of the fields):

    The output will look something like this

    {
            "ALT": "G",
            "CHROM": "4",
            "FILTER": "PASS",
            "ID": ".",
            "INFO": {

                "Gene.ensGene": "ENSG00000109471,ENSG00000138684",
                "Gene.refGene": "IL2,IL21",
                "GeneDetail.ensGene": "dist=38306,dist=117597",
                "GeneDetail.refGene": "dist=38536,dist=117597"
            },
            "POS": 123416186,
            "QUAL" :23.25,
            "REF": "A",
            "SAMPLE": {
                "XG102": {
                    "AD": "51,8",
                    "DP": "59",
                    "GQ": "32",
                    "GT": "0/1",
                    "PL": "32,0,1388"
                }
        }

    Additional hints: The function in part 9 takes in two inputs. 
    input #1 is all the data read from lab1_data.vcf and converted into a 
    list of dictionaries where each dictionary corresponds to a line in the vcf file. 
    input #2 is a dictionary that tells you what the data type of each of the info field is.

    The purpose of part 9 is update each of the fields in "data" input so 
    that the data type matches what you have determined it to be previously.
    POS is an integer and QUAL is a float. For the info fields, you have already 
    created a dictionary called info_field_check_typethat contains the information 
    for the data type of each of the info field. Now use this to cast the info field 
    into the correct data types.

    And the info field goes from being a string to a nested dictionary.

    NOTE: You can only test this function in the last part! There are not tests for it    

    """
    # BEGIN SOLUTION
    '''l=[]
    for dic in data:
        d={}
        dd={}
        for k,v in dic.items():
            if k=="POS":
                d[k]=int(v)
            elif k=="QUAL":
                d[k]=float(v)
            elif k=='INFO':
                p=v.split(';')
                for i in p:
                    if '=' in i:
                        key_value = i.split('=',1)
                        if key_value[1]!='.':
                            if key_value[1].isdigit():
                                dd[key_value[0]] = int(key_value[1])
                            elif key_value[1] and key_value[1].replace('.', '', 1).replace('-','',1).replace('e','2',1).isdigit():
                                dd[key_value[0]] = float(key_value[1])
                            else:
                                dd[key_value[0]] = key_value[1]
                d[k]=dd
            else:
                d[k]=v
        l.append(d)
    new_dict={}
    new_dict.update({'var':l})
    return new_dict'''
    l=[]
    for line in data:   
        line["POS"] = int(line["POS"])
        line["QUAL"] = float(line["QUAL"])
        list_comp=[str(ele) for ele in line["INFO"].split(';')]
        info_dict = create_dictionary_of_info_field_values(list_comp)
        line["INFO"] = {}
        for key, value in info_dict.items():
            if key in info_field_data_type:
                if info_field_data_type[key] == int:
                    value = int(float(value[0]))
                elif info_field_data_type[key] == float:
                    value = float(value[0])
                elif info_field_data_type[key] == str:
                    value = str(value[0])
                line["INFO"][key] = value
        l.append(line)  
    d={}
    d.update({'var':l})
    return d
    # END SOLUTION


def save_data_as_json(data, filename):
    """
    Write a function whose inputs are a Python dictionary and filename. 
    The function will saves the dictionary as a json file using the filename given. 
    Use the json library. 
    Use these options to correctly format your JSON -- 
    sort_keys=True, indent=4, separators=(',', ': '), ensure_ascii=False. 
    Use this function to save your parsed data as a json file.
    """
    # BEGIN SOLUTION
    import json
    with open(filename,'w',encoding='utf-8') as f:
        json.dump(data,f,sort_keys=True, indent=4, separators=(',', ': '), ensure_ascii=False)    

    # END SOLUTION'''


def load_data_from_json(filename):
    """
    Write a function whose input is a filename for a json file. 
    The function should use the filename to read the JSON file in 
    which you saved your final parsed data. 
    """
    # BEGIN SOLUTION
    import json
    with open(filename,encoding='utf-8') as f:
        res = json.load(f)
    return res

    # END SOLUTION


def find_variant(CHROM, REF, ALT, POS, filename):
    """
    Write a function whose inputs are CHROM, REF, ALT, POS, and filename. 
    Using these inputs, the function should load a JSON file using the given 
    filename and return a list of variants that match the given CHROM, REF, ALT, and POS. 
    """
    # BEGIN SOLUTION
    data=load_data_from_json(filename)
    res=data['var']
    l=[]
    for dic in res:
        if dic['CHROM'] == CHROM and dic['REF'] == REF and dic['ALT'] == ALT and dic['POS'] == POS:
            l.append(dic)
    return l
    # END SOLUTION


def pull_basic_and_predictor_fields(filename):
    """
    Load mini_project1_data.json and pull out all the variants that have a 
    """
    # BEGIN SOLUTION
    data=load_data_from_json(filename)

    resdic= []
    for var in data['var']:
        p = var.get('INFO')
        if any(key in p for key in ['FATHMM_pred', 'LRT_pred', 'MetaLR_pred', 'MetaSVM_pred', 'MutationAssessor_pred', 'MutationTaster_pred', 'PROVEAN_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'SIFT_pred', 'fathmm-MKL_coding_pred']):
            basic_info = {'CHROM': var['CHROM'],'POS': var['POS'],'REF': var['REF'],'ALT': var['ALT']}
            predictor_info = {}
            for predictor in ['FATHMM_pred', 'LRT_pred', 'MetaLR_pred', 'MetaSVM_pred', 'MutationAssessor_pred', 'MutationTaster_pred', 'PROVEAN_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'SIFT_pred', 'fathmm-MKL_coding_pred']:
                if predictor in p:
                    value = p[predictor]
                    if value!=None or value!='.':
                        predictor_info[predictor] = value
            output = {**basic_info, **predictor_info}
            resdic.append(output)
    l=[]
    for dic in resdic:
        sum=0
        if dic.get('FATHMM_pred')=='T':
            sum+=0
        elif dic.get('FATHMM_pred')=='D':
            sum+=1
        else:
            sum+=0
        if dic.get('LRT_pred')=='D':
            sum+=1
        elif dic.get('LRT_pred')=='N':
            sum+=0
        elif dic.get('LRT_pred')=='U':
            sum+=0
        else:
            sum+=0
        if dic.get('MetaLR_pred')=='T':
            sum+=0
        elif dic.get('MetaLR_pred')=='D':
            sum+=1
        else:
            sum+=0
        if dic.get('MetaSVM_pred')=='T':
            sum+=0
        elif dic.get('MetaSVM_pred')=='D':
            sum+=1
        else:
            sum+0
        if dic.get('MutationAssessor_pred')=='H':
            sum+=1
        elif dic.get('MutationAssessor_pred')=='N':
            sum+=0
        elif dic.get('MutationAssessor_pred')=='L':
            sum+=0.25
        elif dic.get('MutationAssessor_pred')=='M':
            sum+=0.5
        else:
            sum+=0
        if dic.get('MutationTaster_pred')=='D':
            sum+=1
        elif dic.get('MutationTaster_pred')=='P':
            sum+=0
        elif dic.get('MutationTaster_pred')=='N':
            sum+=0
        elif dic.get('MutationTaster_pred')=='A':
            sum+=1
        else:
            sum+=0
        if dic.get('PROVEAN_pred')=='D':
            sum+=1
        elif dic.get('PROVEAN_pred')=='N':
            sum+=0
        else:
            sum+=0
        if dic.get('Polyphen2_HDIV_pred')=='D':
            sum+=1
        elif dic.get('Polyphen2_HDIV_pred')=='B':
            sum+=0
        elif dic.get('Polyphen2_HDIV_pred')=='P':
            sum+=0.5
        else:
            sum+=0
        if dic.get('Polyphen2_HVAR_pred')=='D':
            sum+=1
        elif dic.get('Polyphen2_HVAR_pred')=='B':
            sum+=0
        elif dic.get('Polyphen2_HVAR_pred')=='P':
            sum+=0.5
        else:
            sum+=0
        if dic.get('SIFT_pred')=='D':
            sum+=1
        elif dic.get('SIFT_pred')=='T':
            sum+=0
        else:
            sum+=0
        l.append(sum)
    #print(l)
    i=0
    for dic in resdic:
        dic.update({'sum_predictor_values':l[i]})
        i+=1
    #from pprint import pprint
    #pprint(resdic)
    z=[]
    for dic in resdic:
        dic.pop('fathmm-MKL_coding_pred')
        z.append(dic)
    return z
    # from pprint import pprint
    # pprint(z)
    # END SOLUTION

def pull_basic_and_predictor_fields_gzip(filename):
    # BEGIN SOLUTION
    import gzip
    resdic= []
    with gzip.open(filename, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            info = fields[7]

            p = {}
            for i in info.split(';'):
                if '=' not in i:
                    continue
                key,value = i.split('=',1)
                if value=='.' or value==None:
                    continue
                else:
                    p[key] = value

            if any(key in p for key in ['FATHMM_pred', 'LRT_pred', 'MetaLR_pred', 'MetaSVM_pred', 'MutationAssessor_pred', 'MutationTaster_pred', 'PROVEAN_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'SIFT_pred','fathmm-MKL_coding_pred']):
                basic_info = {'CHROM': fields[0],'POS': fields[1],'REF': fields[3],'ALT': fields[4]}
                predictor_info = {}
                for predictor in ['FATHMM_pred', 'LRT_pred', 'MetaLR_pred', 'MetaSVM_pred', 'MutationAssessor_pred', 'MutationTaster_pred', 'PROVEAN_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'SIFT_pred','fathmm-MKL_coding_pred']:
                    if predictor in p:
                        value1 = p[predictor]
                        if value1:
                            predictor_info[predictor] = value1
                        else:
                            continue
                output = {**basic_info, **predictor_info}
                resdic.append(output)
    l=[]
    for dic in resdic:
        sum=0
        if dic.get('FATHMM_pred')=='T':
            sum+=0
        elif dic.get('FATHMM_pred')=='D':
            sum+=1
        else:
            sum+=0
        if dic.get('LRT_pred')=='D':
            sum+=1
        elif dic.get('LRT_pred')=='N':
            sum+=0
        elif dic.get('LRT_pred')=='U':
            sum+=0
        else:
            sum+=0
        if dic.get('MetaLR_pred')=='T':
            sum+=0
        elif dic.get('MetaLR_pred')=='D':
            sum+=1
        else:
            sum+=0
        if dic.get('MetaSVM_pred')=='T':
            sum+=0
        elif dic.get('MetaSVM_pred')=='D':
            sum+=1
        else:
            sum+0
        if dic.get('MutationAssessor_pred')=='H':
            sum+=1
        elif dic.get('MutationAssessor_pred')=='N':
            sum+=0
        elif dic.get('MutationAssessor_pred')=='L':
            sum+=0.25
        elif dic.get('MutationAssessor_pred')=='M':
            sum+=0.5
        else:
            sum+=0
        if dic.get('MutationTaster_pred')=='D':
            sum+=1
        elif dic.get('MutationTaster_pred')=='P':
            sum+=0
        elif dic.get('MutationTaster_pred')=='N':
            sum+=0
        elif dic.get('MutationTaster_pred')=='A':
            sum+=1
        else:
            sum+=0
        if dic.get('PROVEAN_pred')=='D':
            sum+=1
        elif dic.get('PROVEAN_pred')=='N':
            sum+=0
        else:
            sum+=0
        if dic.get('Polyphen2_HDIV_pred')=='D':
            sum+=1
        elif dic.get('Polyphen2_HDIV_pred')=='B':
            sum+=0
        elif dic.get('Polyphen2_HDIV_pred')=='P':
            sum+=0.5
        else:
            sum+=0
        if dic.get('Polyphen2_HVAR_pred')=='D':
            sum+=1
        elif dic.get('Polyphen2_HVAR_pred')=='B':
            sum+=0
        elif dic.get('Polyphen2_HVAR_pred')=='P':
            sum+=0.5
        else:
            sum+=0
        if dic.get('SIFT_pred')=='D':
            sum+=1
        elif dic.get('SIFT_pred')=='T':
            sum+=0
        else:
            sum+=0
        l.append(sum)
    #print(l)
    i=0
    for dic in resdic:
        dic.update({'sum_predictor_values':l[i]})
        i+=1
    from pprint import pprint
    #pprint(resdic)
    z=[]
    for dic in resdic:
        dic.pop('fathmm-MKL_coding_pred')
        z.append(dic)
    final=[]
    for dic in z:
        p=dic.keys()
        if len(p)>5:
            final.append(dic)


    import json
    with open('mini_project1_gzip.json', 'w') as f:
        json.dump(final,f,indent=2,sort_keys=True)

    return final

            
    # END SOLUTION

def return_all_non_zero_sum_predictor_values():
    # BEGIN SOLUTION
    import json
    with open('mini_project1_gzip.json') as f:
        data = json.load(f)
    res = []
    for var in data:
        if var['sum_predictor_values'] != 0:
            res.append(var)
    with open('sum_predictor_values_gt_zero.json', 'w') as f:
        json.dump(res,f,indent=2,sort_keys=True)
    # END SOLUTION
              

