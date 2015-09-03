#############################################################################################
##
##Digital_RADs_SingleLineFasta.py
##Version 1.0: Sept 2015
##Created by Michael Sorenson and Jeffrey DaCosta
##Copyright (c) 2012 Boston University. All rights reserved.
##
##Usage with 1 enzyme:
##python3 Digital_RADs_SingleLineFasta.py infile outfile 1 enzyme_seq bases
##
##Usage with 2 enzymes:
##python3 Digital_RADs_SingleLineFasta.py infile outfile 2 enzyme1_seq enzyme2_seq window_start window_end
##
#############################################################################################

import sys, os, re, subprocess

#define base composition function
def basecomp(seq):
    Acount = seq.count('A')
    Ccount = seq.count('C')
    Gcount = seq.count('G')
    Tcount = seq.count('T')
    Ncount = seq.count('N')

    totalGC = Ccount+Gcount
    if len(seq)-Ncount == 0:
        percentGC = 'na'
    else:
        percentGC = round(totalGC/(len(seq)-Ncount),3)
    return(percentGC)


#define complement sequence function
def complement(seq):
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    complseq = [complement[base] for base in seq]
    return complseq

#define reverse-complement sequence function
def reverse_complement(seq):
    seq = list(seq)
    seq.reverse()
    return ''.join(complement(seq))

print('\n********************************************************************************************\n'+
      'Digital_RADs_SingleLineFasta.py\n'+
      'Version 1.0: Sept 2015\n\n'+
      'Created by Michael Sorenson and Jeffrey DaCosta\n'+
      'Copyright (c) 2012 Boston University. All rights reserved.\n\n'+
      'Usage with 1 enzyme:\n'+
      'python3 Digital_RADs.py infile outfile 1 enzyme_seq bases\n\n'+
      'Usage with 2 enzymes:\n'+
      'python3 Digital_RADs.py infile outfile 2 enzyme1_seq enzyme2_seq window_start window_end\n'+
      '\nSee Digital_RADs.readme file for more details.\n'+
      '********************************************************************************************\n')

#define input file: one fasta file with 1-N chromosomes/contigs
infile_name = sys.argv[1]
if '.fa' in infile_name:
    infile = open(infile_name,'r')
elif '.fasta' in infile_name:
    infile = open(infile_name,'r')
elif '.fas' in infile_name:
    infile = open(infile_name,'r')
else:
    print('ERROR: input filename must end with ".fa", ".fas", or ".fasta"\n')
    quit()

#define output file
outfile_name = sys.argv[2]
outf = open(outfile_name,'w')

#define number of enzymes
enz_num = int(sys.argv[3])

if enz_num == 1:
    #define enzyme #1 sequence
    enz1 = sys.argv[4]
    #define window size
    bases = sys.argv[5]
    #print header of output file
    outf.write('contig\tposition\tdownstream\tdown_GC\tupstream\tup_GC\n')

    Enz1Sites = 0
    
    for line in infile:
        if line[0] != '>':
            print('\nERROR: infile does not appear to be in fasta format')
            print('Check that first line begins with ">"\n')
            quit()
        else:
            contigname = line[1:].strip('\n')
            print('Processing '+contigname)
            seq = infile.readline()
            seq = seq.strip('\n')

            sitelist=[]
            position=0
            while position < len(seq):
                site=seq.find(enz1,position)
                if site == -1:
                    break
                Enz1Sites += 1
                position=site+1
                sitelist.append(site)

            sitelist.sort()

            for i in sitelist:
                downfrag = seq[int(i):int(i)+int(bases)]
                upfrag = seq[int(i)+int(len(enz1))-int(bases):int(i)+int(len(enz1))]
                outf.write(contigname+'\t'+str(i+1)+'\t'+downfrag+'\t'+str(basecomp(downfrag))+'\t'+reverse_complement(upfrag)+'\t'+str(basecomp(upfrag))+'\n')

    print('\nFound '+str(Enz1Sites)+' '+enz1+' cut sites\n'+
          '\nFinished!\n\n')
    infile.close()
    outf.close()

elif enz_num == 2:
    #define enzyme #1 sequence
    enz1 = sys.argv[4]
    #define enzyme #2 sequence
    enz2 = sys.argv[5]
    #define beginning of size window
    win_start = int(sys.argv[6])
    #define end of size window
    win_end = int(sys.argv[7])
    
    #print header of output file
    outf.write('contig\tposition\tlength\tdirection\tsequence\tGC\n')

    #set counts to zero
    double1frags = 0
    double2frags = 0
    Enz1Sites = 0
    Enz2Sites = 0
    RADloci = 0
    up_count = 0
    down_count = 0
    size_distr=[0] * 1001

    for line in infile:
        if line[0] != '>':
            print('\nERROR: infile does not appear to be in fasta format')
            print('Check that first line begins with ">"\n')
            quit()
        else:
            contigname = line[1:].strip('\n')
            print('Processing '+contigname)
            seq = infile.readline()
            seq = seq.strip('\n')

            sitelist=[]
            position=0
            while position < len(seq):
                site=seq.find(enz1,position)
                if site == -1:
                    break
                Enz1Sites += 1
                position=site+1
                sitelist.append([site,1])

            position=0
            while position < len(seq):
                site=seq.find(enz2,position)
                if site == -1:
                    break
                Enz2Sites += 1
                position=site+1
                sitelist.append([site,2])

            sitelist.sort()
            for i in range(len(sitelist)-1):
                if sitelist[i][1]!=sitelist[i+1][1]:
                    if sitelist[i][1]==1:
                        size = (sitelist[i+1][0]+len(enz2))-sitelist[i][0]
                        if size >= win_start and size <= win_end:
                            RADloci += 1
                            down_count += 1
                            fullfrag = seq[sitelist[i][0]:sitelist[i+1][0]+len(enz2)]
                            outf.write(contigname+'\t'+str(sitelist[i][0]+1)+'\t'+str(size)+'\t1\t'+fullfrag+'\t'+str(basecomp(fullfrag))+'\n')
                        if size < 1000:
                            size_distr[size]+=1
                        else:
                            size_distr[1000]+=1
                    elif sitelist[i+1][1]==1:
                        size = (sitelist[i+1][0]+len(enz1))-sitelist[i][0]
                        if size >= win_start and size <= win_end:
                            RADloci += 1
                            up_count += 1
                            tag = seq[sitelist[i+1][0]-(size-len(enz1)):sitelist[i+1][0]+len(enz1)]
                            revtag = reverse_complement(tag)
                            outf.write(contigname+'\t'+str(sitelist[i+1][0]+len(enz1))+'\t'+str(size)+'\t-1\t'+revtag+'\t'+str(basecomp(revtag))+'\n')
                        if size < 1000:
                            size_distr[size]+=1
                        else:
                            size_distr[1000]+=1
                else:
                    if sitelist[i][1]==sitelist[i+1][1] and sitelist[i][1]==1:
                        size = (sitelist[i+1][0]+len(enz1))-sitelist[i][0]
                        if size >= win_start and size <= win_end:
                            double1frags += 1
                    else:
                        if sitelist[i][1]==sitelist[i+1][1] and sitelist[i][1]==2:
                            size = (sitelist[i+1][0]+len(enz1))-sitelist[i][0]
                            if size >= win_start and size <= win_end:
                                double2frags += 1

    infile = open(infile_name,'r')
    if '.fa' in infile_name:
        histf = open(infile_name.replace('.fa','_ddRAD.hist'),'w')
    elif '.fasta' in infile_name:
        histf = open(infile_name.replace('.fasta','_ddRAD.hist'),'w')
    elif '.fas' in infile_name:
        histf = open(infile_name.replace('.fas','_ddRAD.hist'),'w')

    base_count = 0
    freq_count = 0
    target = 25
    for i in size_distr:
        base_count += 1
        freq_count += i
        if base_count == target:
            histf.write(str(target-25)+'-'+str(target)+'\t'+str(freq_count)+'\n')
            freq_count = 0
            target = target + 25
        if base_count == 1001:
            break

    print('\nFound '+str(RADloci)+' ddRAD loci in '+str(win_start)+'-'+str(win_end)+' bp window')
    print('\t'+str(up_count)+' with '+enz2+' upstream of '+enz1)
    print('\t'+str(down_count)+' with '+enz2+' downstream of '+enz1+'\n')
    print('Found '+str(Enz1Sites)+' '+enz1+' cut sites')
    print('Found '+str(Enz2Sites)+' '+enz2+' cut sites')
    print('Found '+str(double1frags)+' double '+enz1+' in '+str(win_start)+'-'+str(win_end)+' bp window')
    print('Found '+str(double2frags)+' double '+enz2+' in '+str(win_start)+'-'+str(win_end)+' bp window')
    print('\nFinished!\n')

    infile.close()
    outf.close()

else:
    print('\nError: Number of enzymes must be 1 or 2\n')
