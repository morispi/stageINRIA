#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Comparison of the number of barcodes between real variants and false one.
    BND variants are not considered.
"""

########## RAJOUTER EXCEPTION POUR LES GET


import pysam, xlsxwriter


class Variant:
    '''
        Splits a string (vcf format) to create a variant.
        
        line -- string
    '''
    def __init__(self,line):
        line = line.split()
        self.chrom = line[0]
        self.pos = self.get_pos(line)
        self.id = line[2]
        self.ref = line[3]
        self.alt = line[4]
        self.qual = line[5]
        self.filter = line[6]
        self.info = self.createDict(line)
        
    def createDict(self,line):
        '''
            Creates a dictionnary via all the info.
        '''
        info = line[7].split(";")
        d = {}
        for field in info:
            e = field.split('=')
            if len(e) == 1: # exception for keywords without value
                d[e[0]] = ""
            else:
                d[e[0]] = e[1] 
        return d
        
    def get_pos(self,line):
        '''
            Returns the variant's start position.
        '''
        if self.get_svtype() == "INS":
            return int(line[1]) - len(self.info["LEFT_SVINSSEQ"])
        return int(line[1])

    def get_end(self):
        '''
            Returns the variant's end position.
        '''
        if self.get_svtype() == "INS":
            if "RIGHT_SVINSSEQ" not in self.info:
                return self.pos + int(self.info["SVLEN"])
            else:
                return self.pos + len(self.info["RIGHT_SVINSSEQ"])
        else:
            return int(self.info["END"])

    def get_svtype(self):
        '''
            Returns the variant's type.
        '''
        return self.info["SVTYPE"]
    
    def get_svlen(self):
        '''
            Returns variant's length.
        '''
        if "SVLEN" in self.info:
            return int(self.info["SVLEN"])

def trueSV(file):
    '''
        Returns a list of variants within a file.
        Registers first chromosome, start and end position.

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
    row = row_bis = 0
    realSV = trueSV(truth)
    samfile = pysam.AlignmentFile(bam,"rb")
    workbook = xlsxwriter.Workbook('results.xlsx')
    worksheet = workbook.add_worksheet()
    with open(vcf,"r") as filin:
        # skips file's head :
        line = filin.readline()
        while line.startswith('#'):
            line = filin.readline()
        # for each variant :
        while line != '':
            all_bx = set() # contains all barcodes for the variant
            v = Variant(line)
            if v.get_svtype() != "BND":
                # finds all barcodes for a variant :
                end = v.get_end()
                for read in samfile.fetch(v.chrom,v.pos,end):
                    if read.has_tag('BX'):
                        bx = read.get_tag('BX')
                        all_bx.add(bx)
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
            line = filin.readline()
    workbook.close()
    samfile.close()

####################################################

#sortSV("candidateSV.vcf","possorted_bam.bam","Truth",True)
sortSV("candidateSV.vcf","possorted_bam.bam","Truth",False)

