#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Comparison of the number of barcodes between real variants and false one.
    BND variants are not considered.
"""

import pysam, xlsxwriter


class Variant:
    '''
        Splits a string (vcf format) to create a variant.
        
        line -- string
    '''
    def __init__(self,line):
        line = line.split()
        self.chrom = line[0]
        self.pos = int(line[1])
        self.id = line[2]
        self.ref = line[3]
        self.alt = line[4]
        self.qual = line[5]
        self.filter = line[6]
        self.info = line[7]

    def get_end(self):
        '''
            Returns the end position of the variant.
        '''
        end = (self.info.split(';')[0]).split('=')[1]
        return int(end)

    def get_svtype(self):
        '''
            Returns the type of the variant.
        '''
        svtype= (self.info.split(';')[1]).split('=')[1]
        if len(svtype) != 3:
            svtype = "BND"
        return svtype

def trueSV(file):
    '''
        Returns a list of variants within a file.
        Registers first chromosome, beginning and end position.

        file -- file
    '''
    truth = []
    with open(file,"r") as filin:
        for line in filin:
            sv = line.split()
            truth.append((sv[0],int(sv[1]),int(sv[3])))
    return truth

    
def isValid(variant,L,m):
    '''
        Returns True if a Variant is valid, else False.

        variant -- Variant object
        L -- list of real variants
        m -- int
    '''
    for v in L:
        if variant.chrom==v[0] and abs(v[1]-variant.pos)<=m and abs(v[2]-variant.get_end())<=m:
            return True
    return False


def sortSV(vcf,bam,truth,margin):
    '''
        Creates results.xlsx containing the number of barcodes for real
        variants and false one.

        vcf -- vcf file with variants
        bam -- bam file with reads mapping in the genome reference
        truth -- file with real variants
        margin -- boolean
    '''
    row = row_bis = cpt = 0
    realSV = trueSV(truth)
    samfile = pysam.AlignmentFile(bam,"rb")
    workbook = xlsxwriter.Workbook('results.xlsx')
    worksheet = workbook.add_worksheet()
    with open(vcf,"r") as filin:
        for line in filin:
            if not line.startswith('#'):
                all_bx = [] # contains all barcodes for a variant
                v = Variant(line)
                if v.get_svtype() != "BND":
                    # finds all barcodes for a variant :
                    end = v.get_end()
                    for read in samfile.fetch(v.chrom,v.pos,end):
                        if read.has_tag('BX'):
                            bx = read.get_tag('BX')
                            all_bx.append(bx)
                    m = 100 if margin else 0
                    # variant is valid :
                    if isValid(v,realSV,m):
                        worksheet.write(row,0,v.chrom+":"+str(v.pos)+"-"+str(end))
                        worksheet.write(row,1,len(all_bx))
                        row += 1
                    # variant is not valid :
                    else:
                        worksheet.write(row_bis,3,v.chrom+":"+str(v.pos)+"-"+str(end))
                        worksheet.write(row_bis,4,len(all_bx))
                        row_bis += 1
                    cpt += 1
    workbook.close()
    samfile.close()

####################################################

#sortSV("candidateSV.vcf","possorted_bam.bam","Truth",True)
sortSV("candidateSV.vcf","possorted_bam.bam","Truth",False)

