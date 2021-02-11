# -*- coding: utf-8 -*-
"""
Spyder Editor

ADN -> ARN -> PROTEINA

Que busco Hacer?
1 use ADN and PROTEIN SEQUENCE DATA
2 TRANSLATE ADN to ARN
3 OBTENER ARN SEQUENCE


Pasos para Programar algoritmo:
1 descargar 'ADN SEQUENCE AND PROTEIN SEQUENCE DATA' 
2 import DNA data into Python
3 crear algoritmos que traduce DNA a aminoacids
4 revisar el traduccion sea correcta (con los resultados del profesor :)

/
TRABAJANDO PASO 1
#segmento a leer desde [21 a 938] -> la lectura correcta segun su pagina web
EL ADN  => https://www.ncbi.nlm.nih.gov/nuccore/NM_207618.2 => descarga FASTA

GGTCAGAAAAAGCCCTCTCCATGTCTACTCACGATACATCCCTGAAAACCACTGAGGAAGTGGCTTTTCA
GATCATCTTGCTTTGCCAGTTTGGGGTTGGGACTTTTGCCAATGTATTTCTCTTTGTCTATAATTTCTCT
CCAATCTCGACTGGTTCTAAACAGAGGCCCAGACAAGTGATTTTAAGACACATGGCTGTGGCCAATGCCT
TAACTCTCTTCCTCACTATATTTCCAAACAACATGATGACTTTTGCTCCAATTATTCCTCAAACTGACCT
CAAATGTAAATTAGAATTCTTCACTCGCCTCGTGGCAAGAAGCACAAACTTGTGTTCAACTTGTGTTCTG
AGTATCCATCAGTTTGTCACACTTGTTCCTGTTAATTCAGGTAAAGGAATACTCAGAGCAAGTGTCACAA
ACATGGCAAGTTATTCTTGTTACAGTTGTTGGTTCTTCAGTGTCTTAAATAACATCTACATTCCAATTAA
GGTCACTGGTCCACAGTTAACAGACAATAACAATAACTCTAAAAGCAAGTTGTTCTGTTCCACTTCTGAT
TTCAGTGTAGGCATTGTCTTCTTGAGGTTTGCCCATGATGCCACATTCATGAGCATCATGGTCTGGACCA
GTGTCTCCATGGTACTTCTCCTCCATAGACATTGTCAGAGAATGCAGTACATATTCACTCTCAATCAGGA
CCCCAGGGGCCAAGCAGAGACCACAGCAACCCATACTATCCTGATGCTGGTAGTCACATTTGTTGGCTTT
TATCTTCTAAGTCTTATTTGTATCATCTTTTACACCTATTTTATATATTCTCATCATTCCCTGAGGCATT
GCAATGACATTTTGGTTTCGGGTTTCCCTACAATTTCTCCTTTACTGTTGACCTTCAGAGACCCTAAGGG
TCCTTGTTCTGTGTTCTTCAACTGTTGAAAGCCAGAGTCACTAAAAATGCCAAACACAGAAGACAGCTTT
GCTAATACCATTAAATACTTTATTCCATAAATATGTTTTTAAAAGCTTGTATGAACAAGGTATGGTGCTC
ACTGCTATACTTATAAAAGAGTAAGGTTATAATCACTTGTTGATATGAAAAGATTTCTGGTTGGAATCTG
ATTGAAACAGTGAGTTATTCACCACCCTCCATTCTCT


Su secuencia de proteínas => https://www.ncbi.nlm.nih.gov/nuccore/NM_207618.2 => descarga CDS

MSTHDTSLKTTEEVAFQIILLCQFGVGTFANVFLFVYNFSPIST
GSKQRPRQVILRHMAVANALTLFLTIFPNNMMTFAPIIPQTDLKCKLEFFTRLVARST
NLCSTCVLSIHQFVTLVPVNSGKGILRASVTNMASYSCYSCWFFSVLNNIYIPIKVTG
PQLTDNNNNSKSKLFCSTSDFSVGIVFLRFAHDATFMSIMVWTSVSMVLLLHRHCQRM
QYIFTLNQDPRGQAETTATHTILMLVVTFVGFYLLSLICIIFYTYFIYSHHSLRHCND
ILVSGFPTISPLLLTFRDPKGPCSVFFNC

TRABAJANDO PASO 2

"""


def leer_squence(file):
    """lee una secuencia y la limpia
    """
    if 1 == 1:
        with open(file, 'r') as f:  #funcion with open :D
            sequence = f.read()
        # remover espacios que salen en la consola
        sequence = sequence.replace('\n', '')
        sequence = sequence.replace('\r', '')
    return sequence


"""TRABAJANDO PASO 3"""


def adn_arn(sequence):
    """
    Traducir una cadena que contiene una secuencia de nucleótidos(adn) 
    en una cadena que contiene el correspondiente secuencias de aminoácidos.(arn)
    los nucleótidos se traducen en tripletes utilizando el diccionario o tablas :)

    ver esto con help(transcription)


    Parameters
    ----------
    sequence : TYPE
        DESCRIPTION.

    Returns
    -------
    protein : TYPE
        DESCRIPTION.
    """

    # tabla para identificar scuencia
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }

    #ALGORITMOS :D
    if len(sequence) % 3 == 0:
        protein=""
        for i in range(0, len(sequence), 3):
            codon   = sequence[i : i+3]
            protein += table[codon]
        return protein


"""TRABAJANDO PASO 4 resultado"""

ad = './adn.txt'
pr = "./protein.txt"


leido = leer_squence(ad)
res = adn_arn(leido[20:935])   # mi arn obtenido

arnsegunWeb = leer_squence(pr) # arn segun la web

if res == arnsegunWeb:
    print('la secuncia de NM_207618.2 esta comprobada por el NCBI')

#res = adn_arn(sec[20:398])
# print(res)
#res = leer_squence(proteina)
# print(res)
