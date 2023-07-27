import numpy as np

import methods._01_ORF_code as ORF_code
import methods._02_CTDcode as CTDcode
import methods._03_Fickettcode as Fickettcode
import methods._04_kmer_counts_molmap as kmer_counts
# import methods._05_Hexamercode as Hexamercode
import methods._06_proparcoder as proparcoder
import methods._07_GCcounts as GCcounts
import methods._08_edpfeature as edpfeature
import methods._09_StopCodon as StopCodon
import methods._10_onehot as onehot
import methods._15_SparseEncoding as SparseEncod
import methods._18_SStructure as SStructure
import methods._19_CTDcoder_st_ph as CTDcoderstph

import sys
import os
ac_path = os.path.join(os.path.dirname(__file__), 'repDNA')
sys.path.append(ac_path)
import ac as RNA_ac
import psenac as RNA_psenac
from rpy2.robjects import pandas2ri
import rpy2.robjects as robjects

import pandas as pd
# from Bio import SeqIO

dictMe = {
    'Open reading frame (1D)': '1',
    'Entropy density of transcript (1D)': '2_1',
    'Global descriptor (1D)': '2_2',
    'K-mer (1D)': '2_3',
    'Codon related (1D)': '3',
    'Pseudo protein related (1D)': '6',
    'Guanine-cytosine related (1D)': '7',
    'Sequence-intrinsic Features:One-hot encoding (2D)': '13',
    'Sparse encoding (2D)': '15',
    'Structure-based Features:One-hot encoding (2D)':'16',
    'Nucleotide related (1D)':'17',
    'Secondary structure (1D)':'18',
    'EIIP based spectrum (1D)':'18_1',
    
    'Solubility lipoaffinity (1D)':'19_1',
    'Partition coefficient (1D)':'19_101',
    'Polarizability refractivity (1D)':'19_2',
    'Hydrogen bond related (1D)':'19_3',
    'Topological indice (1D)':'19_4',
    'Molecular fingerprint (1D)':'19_5'
}

from rpy2.rinterface import StrSexpVector
def tuple_str(tpl):
    res = StrSexpVector(tpl)
    return res

def switch_meth(fun, textPath):
    if fun == '1':
        
        ORF = ORF_code.ORF_count(textPath).get_ORF()
        ORFedp = edpfeature.EDPcoder(textPath).getEDP_orf()
        hexmerORF = SStructure.makeORFEucDist(textPath)

        rownames_list = [hexmerORF.rownames[i] for i in range(len(hexmerORF.rownames))]

        pd_hexmerORF = pd.DataFrame(hexmerORF)
        pd_hexmerORF = pd_hexmerORF.T

        ORFedp_index = pd.DataFrame(ORFedp.index)
        ORFedp_index.index= ORFedp_index.iloc[:,0]

        pd_hexmerORF.index = rownames_list
        pd_hexmerORF_com = pd.merge(ORFedp_index, pd_hexmerORF,  left_index=True, how='left',right_index=True, sort=False)

        pd_hexmerORF_com01 = pd_hexmerORF.fillna(0)
        pd_hexmerORF_com02 = pd_hexmerORF_com01
        # print('pd_hexmerORF_com02.shape'.format(pd_hexmerORF_com02.shape))
        pd_hexmerORF_com02.columns = ['EucDist.LNC_orf', 'EucDist.PCT_orf', 'EucDist.Ratio_orf', 'LogDist.LNC_orf', 'LogDist.PCT_orf', 'LogDist.Ratio_orf', 'Hexamer.Score_orf']

        T1 = pd.concat([ORF, ORFedp], axis=1, join='inner')
        T1 = pd.concat([T1, pd_hexmerORF_com02], axis=1, join='inner')
        return T1
    elif fun == '2_1':
        EDP = edpfeature.EDPcoder(textPath).getEDP()        
        return EDP
    elif fun == '2_2':
        CTD = CTDcode.CTDcoder(textPath).get_ctd()
        return CTD
    elif fun == '2_3':
        Kmer1 = kmer_counts.BasicCounter(textPath, int(1)).get_counts()
        Kmer2 = kmer_counts.BasicCounter(textPath, int(2)).get_counts()
        Kmer3 = kmer_counts.BasicCounter(textPath, int(3)).get_counts()
        # Kmer4 = kmer_counts.BasicCounter(textPath, int(4)).get_counts()

        T1 = pd.concat([Kmer1, Kmer2], axis=1, join='inner')
        T1 = pd.concat([T1, Kmer3], axis=1, join='inner')
        return T1
    elif fun == '3':
        Fickett = Fickettcode.Fickettcoder(textPath).get_fickett()
        StopCod = StopCodon.get_stop(textPath)
        return pd.concat([Fickett, StopCod], axis=1, join='inner')
    elif fun == '6':
        return proparcoder.ProtPar(textPath).get_protper()
    elif fun == '7':
        return GCcounts.GCconder(textPath).get_gc()
    elif fun == '13':
        return onehot.Onehot(textPath).get_onehot(1000) #2D 方法
    elif fun == '15':
        return SparseEncod.get_encoding(textPath, 1000)
    elif fun == '16':
        seqname1, rnaresult1 = onehot.Onehot(textPath).get_onehot(1000)
        seqname2, rnaresult2 = SparseEncod.get_encoding(textPath, 1000)
        rnaresult3 = np.concatenate((rnaresult1, rnaresult2), axis=2)
        return seqname1, rnaresult3
    elif fun == '17':
        dac = RNA_ac.rna_dac(textPath,2)   
        dcc = RNA_ac.rna_dcc(textPath,2)
        psenac = RNA_psenac.rna_pc_psednc(textPath)
        SCPseDNC = RNA_psenac.rna_SCPseDNC(textPath)


        T1 = pd.concat([dac, dcc], axis=1, join='inner')
        T1 = pd.concat([T1, psenac], axis=1, join='inner')
        T1 = pd.concat([T1, SCPseDNC], axis=1, join='inner')

        return T1
    elif fun == '18':

        UTR_len = edpfeature.EDPcoder(textPath).getUTR_len()
        index_nms = UTR_len.index
        SStruc = SStructure.extract_SSfeatures(textPath)
        SStruc_T = pandas2ri.rpy2py(SStruc)
        # print(SStruc_T)
        SStruc_T.index = index_nms
        SStruc_T.columns = ["SLDLD: Structural logarithm distance to lncRNA of acguD", "SLDPD: Structural logarithm distance to pcRNA of acguD", "SLDRD: Structural logarithm distance acguD ratio", "SLDLN: Structural logarithm distance to lncRNA of acguACGU", "SLDPN: Structural logarithm distance to pcRNA of acguACGU", "SLDRN: Structural logarithm distance acguACGU ratio","SDMFE: Secondary structural minimum free energy", "SFPUS: Secondary structural UP frequency paired-unpaired"]

        return SStruc_T
    elif fun == '18_1':
        tran_len = edpfeature.EDPcoder(textPath).get_tran_len()
        index_nms = tran_len.index
        
        SStruc = SStructure.makeEIIP(textPath)        
        SStruc = pd.DataFrame(SStruc)
        SStruc_T = SStruc.T
        SStruc_T.index = index_nms
        SStruc_T.columns = ['EipSP: Electron-ion interaction pseudopotential signal peak','EipAP: Electron-ion interaction pseudopotential average power','EiSNR: Electron-ion interaction pseudopotential signal/noise ratio','EiPS0: Electron-ion interaction pseudopotential spectrum 0','EiPS1: Electron-ion interaction pseudopotential spectrum 0.25','EiPS2: Electron-ion interaction pseudopotential spectrum 0.5','EiPS3: Electron-ion interaction pseudopotential spectrum 0.75','EiPS4: Electron-ion interaction pseudopotential spectrum 1']


        return SStruc_T
    elif fun == '19_1':
        ACGT_encode = ['1000']
        Solubility = CTDcoderstph.CTDcoder(textPath,ACGT_encode).get_ctd()
        
        
        ACGT_encode_3 = ['1020']
        Solubility_3 = CTDcoderstph.CTDcoder_3class(textPath,ACGT_encode_3).get_ctd()
        T1 = pd.concat([Solubility, Solubility_3], axis=1, join='inner')
        # Kmers
        for i in range(1,4,1):
            Kmer_Solubility = kmer_counts.BasicCounter_2(textPath, int(i), encode = ACGT_encode).get_counts() 
            Kmer_Solubility_3 = kmer_counts.BasicCounter_3(textPath, int(i), encode = ACGT_encode_3).get_counts() 
            
            T1 = pd.concat([T1, Kmer_Solubility], axis=1, join='inner')
            T1 = pd.concat([T1, Kmer_Solubility_3], axis=1, join='inner')
        return T1
    elif fun == '19_101':
        ACGT_encode = ['0100']
        Partition = CTDcoderstph.CTDcoder(textPath,ACGT_encode).get_ctd()
        
        ACGT_encode_3 = ['0102']
        Partition_3 = CTDcoderstph.CTDcoder_3class(textPath,ACGT_encode_3).get_ctd()
        T1 = pd.concat([Partition, Partition_3], axis=1, join='inner')
        # Kmers
        for i in range(1,4,1):
            Kmer_Partition = kmer_counts.BasicCounter_2(textPath, int(i), encode = ACGT_encode).get_counts() 
            Kmer_Partition_3 = kmer_counts.BasicCounter_3(textPath, int(i), encode = ACGT_encode_3).get_counts() 
       
            T1 = pd.concat([T1, Kmer_Partition], axis=1, join='inner')
            T1 = pd.concat([T1, Kmer_Partition_3], axis=1, join='inner')
        
        return T1
    elif fun == '19_2':
        ACGT_encode = ['0101']
        T1 = CTDcoderstph.CTDcoder(textPath,ACGT_encode).get_ctd()

        # Kmers
        for i in range(1,4,1):
            Kmer_Polarizability = kmer_counts.BasicCounter_2(textPath, int(i), encode = ACGT_encode).get_counts() 
            T1 = pd.concat([T1, Kmer_Polarizability], axis=1, join='inner')
      
        return T1
    elif fun == '19_3':
        ACGT_encode = ['0010', '0110']
        Hydrogen = CTDcoderstph.CTDcoder(textPath,ACGT_encode).get_ctd()

        ACGT_encode_3 = ['1200', '0120']
        Hydrogen_3 = CTDcoderstph.CTDcoder_3class(textPath,ACGT_encode_3).get_ctd()
        T1 = pd.concat([Hydrogen, Hydrogen_3], axis=1, join='inner')

        ACGT_encode_00 = ['0010']
        ACGT_encode_01 = ['0110']
        # Kmers
        for i in range(1,4,1):
            Kmer_Hydrogen = kmer_counts.BasicCounter_2(textPath, int(i), encode = ACGT_encode_00).get_counts() 
            Kmer_Hydrogen_01 = kmer_counts.BasicCounter_2(textPath, int(i), encode = ACGT_encode_01).get_counts() 

            Kmer_Hydrogen_3 = kmer_counts.BasicCounter_3(textPath, int(i), encode = [ACGT_encode_3[0]]).get_counts()
            Kmer_Hydrogen_3_1 = kmer_counts.BasicCounter_3(textPath, int(i), encode = [ACGT_encode_3[1]]).get_counts()
       
            T1 = pd.concat([T1, Kmer_Hydrogen], axis=1, join='inner')
            T1 = pd.concat([T1, Kmer_Hydrogen_01], axis=1, join='inner')
            T1 = pd.concat([T1, Kmer_Hydrogen_3], axis=1, join='inner')
            T1 = pd.concat([T1, Kmer_Hydrogen_3_1], axis=1, join='inner')
        
        return T1
    elif fun == '19_4':
        ACGT_encode = ['0001']
        Topological = CTDcoderstph.CTDcoder(textPath,ACGT_encode).get_ctd()

        ACGT_encode_3 = ['1002', '0021']
        Topological_3 = CTDcoderstph.CTDcoder_3class(textPath,ACGT_encode_3).get_ctd()
        T1 = pd.concat([Topological, Topological_3], axis=1, join='inner')
        # Kmers
        for i in range(1,4,1):
            Kmer_Topological = kmer_counts.BasicCounter_2(textPath, int(i), encode = ACGT_encode).get_counts() 
            Kmer_Topological_3 = kmer_counts.BasicCounter_3(textPath, int(i), encode = [ACGT_encode_3[0]]).get_counts() 
            Kmer_Topological_3_1 = kmer_counts.BasicCounter_3(textPath, int(i), encode = [ACGT_encode_3[1]]).get_counts() 
           
            T1 = pd.concat([T1, Kmer_Topological], axis=1, join='inner')
            T1 = pd.concat([T1, Kmer_Topological_3], axis=1, join='inner')
            T1 = pd.concat([T1, Kmer_Topological_3_1], axis=1, join='inner')
        
        return T1
    elif fun == '19_5':
        ACGT_encode = ['0011']
        T1 = CTDcoderstph.CTDcoder(textPath,ACGT_encode).get_ctd()

        # Kmers
        for i in range(1,4,1):
            Kmer_Molecular = kmer_counts.BasicCounter_2(textPath, int(i), encode = ACGT_encode).get_counts() 

            T1 = pd.concat([T1, Kmer_Molecular], axis=1, join='inner')
        return T1
    else:
        return None


# new feature ========================================================================

list_new_featrures = ['SHDCA: Strong H-Bond donors composition of A','SHDAB: Strong H-Bond donors transition between A and B','SHDA0: Strong H-Bond donors distribution of 0.00A','SHDA1: Strong H-Bond donors distribution of 0.25A','SHDA2: Strong H-Bond donors distribution of 0.50A','SHDA3: Strong H-Bond donors distribution of 0.75A','SHDA4: Strong H-Bond donors distribution of 1.00A','MLFCA: Linear free energy composition of A','MLFCB: Linear free energy composition of B','MLFAB: Linear free energy transition between A and B','MLFA0: Linear free energy distribution of 0.00A','MLFA1: Linear free energy distribution of 0.25A','MLFA2: Linear free energy distribution of 0.50A','MLFA3: Linear free energy distribution of 0.75A','MLFA4: Linear free energy distribution of 1.00A','MLFB0: Linear free energy distribution of 0.00B','MLFB1: Linear free energy distribution of 0.25B','MLFB2: Linear free energy distribution of 0.50B','MLFB3: Linear free energy distribution of 0.75B','MLFB4: Linear free energy distribution of 1.00B','SHAAB: Strong H-Bond acceptors_3 transition between A and B','SHAAC: Strong H-Bond acceptors_3 transition between A and C','PHBAB: Potential Hydrogen Bonds_3 transition between A and B','PHBAC: Potential Hydrogen Bonds_3 transition between A and C','SHD: Strong H-Bond donors0','MLF: Linear free energy0','MLF: Linear free energy1','SHD: Strong H-Bond donors00','SHD: Strong H-Bond donors01','SHD: Strong H-Bond donors10','MLF: Linear free energy00','MLF: Linear free energy01','MLF: Linear free energy10','MLF: Linear free energy11','SHA: Strong H-Bond acceptors_301','SHA: Strong H-Bond acceptors_302','SHA: Strong H-Bond acceptors_310','SHA: Strong H-Bond acceptors_320','PHB: Potential Hydrogen Bonds_301','PHB: Potential Hydrogen Bonds_302','PHB: Potential Hydrogen Bonds_310','PHB: Potential Hydrogen Bonds_320','SHD: Strong H-Bond donors000','SHD: Strong H-Bond donors001','SHD: Strong H-Bond donors010','SHD: Strong H-Bond donors011','SHD: Strong H-Bond donors100','SHD: Strong H-Bond donors101','SHD: Strong H-Bond donors110','MLF: Linear free energy000','MLF: Linear free energy001','MLF: Linear free energy010','MLF: Linear free energy011','MLF: Linear free energy100','MLF: Linear free energy101','MLF: Linear free energy110','MLF: Linear free energy111','SHA: Strong H-Bond acceptors_3001','SHA: Strong H-Bond acceptors_3002','SHA: Strong H-Bond acceptors_3010','SHA: Strong H-Bond acceptors_3011','SHA: Strong H-Bond acceptors_3012','SHA: Strong H-Bond acceptors_3020','SHA: Strong H-Bond acceptors_3021','SHA: Strong H-Bond acceptors_3022','SHA: Strong H-Bond acceptors_3100','SHA: Strong H-Bond acceptors_3101','SHA: Strong H-Bond acceptors_3102','SHA: Strong H-Bond acceptors_3110','SHA: Strong H-Bond acceptors_3120','SHA: Strong H-Bond acceptors_3200','SHA: Strong H-Bond acceptors_3201','SHA: Strong H-Bond acceptors_3202','SHA: Strong H-Bond acceptors_3210','SHA: Strong H-Bond acceptors_3220','PHB: Potential Hydrogen Bonds_3001','PHB: Potential Hydrogen Bonds_3002','PHB: Potential Hydrogen Bonds_3010','PHB: Potential Hydrogen Bonds_3011','PHB: Potential Hydrogen Bonds_3012','PHB: Potential Hydrogen Bonds_3020','PHB: Potential Hydrogen Bonds_3021','PHB: Potential Hydrogen Bonds_3022','PHB: Potential Hydrogen Bonds_3100','PHB: Potential Hydrogen Bonds_3101','PHB: Potential Hydrogen Bonds_3102','PHB: Potential Hydrogen Bonds_3110','PHB: Potential Hydrogen Bonds_3120','PHB: Potential Hydrogen Bonds_3200','PHB: Potential Hydrogen Bonds_3201','PHB: Potential Hydrogen Bonds_3202','PHB: Potential Hydrogen Bonds_3210','PHB: Potential Hydrogen Bonds_3220','CNHCA: NH- count composition of A','CNHCB: NH- count composition of B','CNHAB: NH- count transition between A and B','CNHA0: NH- count distribution of 0.00A','CNHA1: NH- count distribution of 0.25A','CNHA2: NH- count distribution of 0.50A','CNHA3: NH- count distribution of 0.75A','CNHA4: NH- count distribution of 1.00A','CNHB0: NH- count distribution of 0.00B','CNHB1: NH- count distribution of 0.25B','CNHB2: NH- count distribution of 0.50B','CNHB3: NH- count distribution of 0.75B','CNHB4: NH- count distribution of 1.00B','CNH: NH- count0','CNH: NH- count1','CNH: NH- count00','CNH: NH- count01','CNH: NH- count10','CNH: NH- count11','CNH: NH- count000','CNH: NH- count001','CNH: NH- count010','CNH: NH- count011','CNH: NH- count100','CNH: NH- count101','CNH: NH- count110','CNH: NH- count111','MReCA: Molar refractivity composition of A','MReCB: Molar refractivity composition of B','MReAB: Molar refractivity transition between A and B','MReA0: Molar refractivity distribution of 0.00A','MReA1: Molar refractivity distribution of 0.25A','MReA2: Molar refractivity distribution of 0.50A','MReA3: Molar refractivity distribution of 0.75A','MReA4: Molar refractivity distribution of 1.00A','MReB0: Molar refractivity distribution of 0.00B','MReB1: Molar refractivity distribution of 0.25B','MReB2: Molar refractivity distribution of 0.50B','MReB3: Molar refractivity distribution of 0.75B','MReB4: Molar refractivity distribution of 1.00B','MRe: Molar refractivity0','MRe: Molar refractivity1','MRe: Molar refractivity00','MRe: Molar refractivity01','MRe: Molar refractivity10','MRe: Molar refractivity11','MRe: Molar refractivity000','MRe: Molar refractivity001','MRe: Molar refractivity010','MRe: Molar refractivity011','MRe: Molar refractivity100','MRe: Molar refractivity101','MRe: Molar refractivity110','MRe: Molar refractivity111','PSNCA: Primary or secondary nitrogens composition of A','PSNAB: Primary or secondary nitrogens transition between A and B','PSNA0: Primary or secondary nitrogens distribution of 0.00A','PSNA1: Primary or secondary nitrogens distribution of 0.25A','PSNA2: Primary or secondary nitrogens distribution of 0.50A','PSNA3: Primary or secondary nitrogens distribution of 0.75A','PSNA4: Primary or secondary nitrogens distribution of 1.00A','SLFAB: Sum of path lengths starting from oxygens_3 transition between A and B','SLFAC: Sum of path lengths starting from oxygens_3 transition between A and C','TPSAB: Topological polar surface area_3 transition between A and B','TPSAC: Topological polar surface area_3 transition between A and C','PSN: Primary or secondary nitrogens0','PSN: Primary or secondary nitrogens00','PSN: Primary or secondary nitrogens01','PSN: Primary or secondary nitrogens10','SLF: Sum of path lengths starting from oxygens_301','SLF: Sum of path lengths starting from oxygens_302','SLF: Sum of path lengths starting from oxygens_310','SLF: Sum of path lengths starting from oxygens_320','TPS: Topological polar surface area_301','TPS: Topological polar surface area_302','TPS: Topological polar surface area_310','TPS: Topological polar surface area_320','PSN: Primary or secondary nitrogens000','PSN: Primary or secondary nitrogens001','PSN: Primary or secondary nitrogens010','PSN: Primary or secondary nitrogens011','PSN: Primary or secondary nitrogens100','PSN: Primary or secondary nitrogens101','PSN: Primary or secondary nitrogens110','SLF: Sum of path lengths starting from oxygens_3001','SLF: Sum of path lengths starting from oxygens_3002','SLF: Sum of path lengths starting from oxygens_3010','SLF: Sum of path lengths starting from oxygens_3011','SLF: Sum of path lengths starting from oxygens_3012','SLF: Sum of path lengths starting from oxygens_3020','SLF: Sum of path lengths starting from oxygens_3021','SLF: Sum of path lengths starting from oxygens_3022','SLF: Sum of path lengths starting from oxygens_3100','SLF: Sum of path lengths starting from oxygens_3101','SLF: Sum of path lengths starting from oxygens_3102','SLF: Sum of path lengths starting from oxygens_3110','SLF: Sum of path lengths starting from oxygens_3120','SLF: Sum of path lengths starting from oxygens_3200','SLF: Sum of path lengths starting from oxygens_3201','SLF: Sum of path lengths starting from oxygens_3202','SLF: Sum of path lengths starting from oxygens_3210','SLF: Sum of path lengths starting from oxygens_3220','TPS: Topological polar surface area_3001','TPS: Topological polar surface area_3002','TPS: Topological polar surface area_3010','TPS: Topological polar surface area_3011','TPS: Topological polar surface area_3012','TPS: Topological polar surface area_3020','TPS: Topological polar surface area_3021','TPS: Topological polar surface area_3022','TPS: Topological polar surface area_3100','TPS: Topological polar surface area_3101','TPS: Topological polar surface area_3102','TPS: Topological polar surface area_3110','TPS: Topological polar surface area_3120','TPS: Topological polar surface area_3200','TPS: Topological polar surface area_3201','TPS: Topological polar surface area_3202','TPS: Topological polar surface area_3210','TPS: Topological polar surface area_3220','HPCCA: Gas-hexadecane PC composition of A','HPCAB: Gas-hexadecane PC transition between A and B','HPCA0: Gas-hexadecane PC distribution of 0.00A','HPCA1: Gas-hexadecane PC distribution of 0.25A','HPCA2: Gas-hexadecane PC distribution of 0.50A','HPCA3: Gas-hexadecane PC distribution of 0.75A','HPCA4: Gas-hexadecane PC distribution of 1.00A','HPCAB: Gas-hexadecane PC_3 transition between A and B','HPCAC: Gas-hexadecane PC_3 transition between A and C','HPC: Gas-hexadecane PC0','HPC: Gas-hexadecane PC00','HPC: Gas-hexadecane PC01','HPC: Gas-hexadecane PC10','HPC: Gas-hexadecane PC_301','HPC: Gas-hexadecane PC_302','HPC: Gas-hexadecane PC_310','HPC: Gas-hexadecane PC_320','HPC: Gas-hexadecane PC000','HPC: Gas-hexadecane PC001','HPC: Gas-hexadecane PC010','HPC: Gas-hexadecane PC011','HPC: Gas-hexadecane PC100','HPC: Gas-hexadecane PC101','HPC: Gas-hexadecane PC110','HPC: Gas-hexadecane PC_3001','HPC: Gas-hexadecane PC_3002','HPC: Gas-hexadecane PC_3010','HPC: Gas-hexadecane PC_3011','HPC: Gas-hexadecane PC_3012','HPC: Gas-hexadecane PC_3020','HPC: Gas-hexadecane PC_3021','HPC: Gas-hexadecane PC_3022','HPC: Gas-hexadecane PC_3100','HPC: Gas-hexadecane PC_3101','HPC: Gas-hexadecane PC_3102','HPC: Gas-hexadecane PC_3110','HPC: Gas-hexadecane PC_3120','HPC: Gas-hexadecane PC_3200','HPC: Gas-hexadecane PC_3201','HPC: Gas-hexadecane PC_3202','HPC: Gas-hexadecane PC_3210','HPC: Gas-hexadecane PC_3220','LFICA: Lipoaffinity index composition of A','LFIAB: Lipoaffinity index transition between A and B','LFIA0: Lipoaffinity index distribution of 0.00A','LFIA1: Lipoaffinity index distribution of 0.25A','LFIA2: Lipoaffinity index distribution of 0.50A','LFIA3: Lipoaffinity index distribution of 0.75A','LFIA4: Lipoaffinity index distribution of 1.00A','LFIAB: Lipoaffinity index_3 transition between A and B','LFIAC: Lipoaffinity index_3 transition between A and C','LFI: Lipoaffinity index0','LFI: Lipoaffinity index00','LFI: Lipoaffinity index01','LFI: Lipoaffinity index10','LFI: Lipoaffinity index_301','LFI: Lipoaffinity index_302','LFI: Lipoaffinity index_310','LFI: Lipoaffinity index_320','LFI: Lipoaffinity index000','LFI: Lipoaffinity index001','LFI: Lipoaffinity index010','LFI: Lipoaffinity index011','LFI: Lipoaffinity index100','LFI: Lipoaffinity index101','LFI: Lipoaffinity index110','LFI: Lipoaffinity index_3001','LFI: Lipoaffinity index_3002','LFI: Lipoaffinity index_3010','LFI: Lipoaffinity index_3011','LFI: Lipoaffinity index_3012','LFI: Lipoaffinity index_3020','LFI: Lipoaffinity index_3021','LFI: Lipoaffinity index_3022','LFI: Lipoaffinity index_3100','LFI: Lipoaffinity index_3101','LFI: Lipoaffinity index_3102','LFI: Lipoaffinity index_3110','LFI: Lipoaffinity index_3120','LFI: Lipoaffinity index_3200','LFI: Lipoaffinity index_3201','LFI: Lipoaffinity index_3202','LFI: Lipoaffinity index_3210','LFI: Lipoaffinity index_3220']


