import math
import os
import json
import matplotlib


class Matrix():

    def __init__(self, ref_seqpath, fas_seqpath) -> None:
        
        self.reference_path = ref_seqpath
        self.query_path = fas_seqpath
        self.reference = ''
        self.query = ''

        # default 2000 bases per slice. (EX: 5,000,000 bp ref / 2,000 = 2,500 slices)
        # The index of a slice is ((slice number - 1) * 2000) - 1
        self.reference_slicehandler = []

        # index matches index of slice in slicehandler
        # value is the score of the matching slice
        self.query_scorehandler = []

        # default 24 bases per slice
        # slice's position in query string = ((slice number - 1) * 24) - 1
        self.query_slicehandler = []

        # index of even value: index of last common. odd: index of starting common.
        # Represents Gaps between the Query and the Ref. Could be a polymorphism
        # Could be an intron. Who knows.
        # EX: [2, 5, 7, 11]
        # ->  |||--|||---||||||
        self.ref_gaps = []
        
        # indexes of imperfect matches. Each query slice is aligned as best as possible
        # but if it isn't a perfect score, the next best is chosen and the differences
        # noted. 
        self.query_diffs = []  
        pass
    
    def build_Strings(self):
        with open(self.reference_path) as fileA, open(self.query_path) as fileB:
            self.reference = "".join([line.strip() for line in fileA if not line.startswith(">")])
            self.query = "".join([line.strip() for line in fileB if not line.startswith(">")])
        self.reference = self.reference.replace("N", "")
        return 0

    def build_QSlices(self, window=24):
        sizeofQuery = len(self.query)
        remainder = sizeofQuery % window

        countofSlices = int((sizeofQuery - remainder) / window)
        self.query_slicehandler.append(0)
        for i in range(2, (countofSlices + 1)):
            slicepos = ((i - 1) * window) 
            self.query_slicehandler.append(slicepos)
        
        self.query_slicehandler.append(len(self.query)-remainder)

        return 0

    def build_RSlices(self, window=2000):
        if int(len(self.reference) / window) <= 2:
            print("WINDOW of", window, "BASES is TOO SMALL FOR THE REFERENCE. LOWERING WINDOW SIZE BY FACTOR of 10.")
            window = int(window / 10)
        sizeofRef = len(self.reference)
        remainder = sizeofRef % window

        countofSlices = int((sizeofRef - remainder) / window)
        self.reference_slicehandler.append(0)

        for i in range(2, (countofSlices + 1)):
            slicepos = ((i - 1) * window)
            self.reference_slicehandler.append(slicepos)
        
        self.reference_slicehandler.append(len(self.reference)-remainder)

        return 0
    
    def run_slice(self, Rslice, Qslice, threshold):
        score = 0

        # Score = X of matches / Total matches
        # EX: 1/18 = 0.055. 16/18 = 0.889
        # The acceptable min score is threhold / total
        # default: 18/24 = 0.667

        startpoint = int(len(Qslice)/2)
        section_score = []
        best = 0

        for j in range (startpoint-1, len(Rslice)):

            # TEST THROTTLE. Remove for full run
            #if j == (startpoint + 5):
            #    break
            matches = 0
            count = 0
            streak = 0

            TEMPA = []
            TEMPB = []
            
            for i in range (1, len(Qslice)):
                #Don't be a range problem
                if j - count < 0 or count >= 23:
                    break
                #Reverse compare (tail to head) the query to
                #the section of the ref slice.
                #print(Qslice[i], "vs", Qslice[-i])
                TEMPA.append(Qslice[-i])
                TEMPB.append(Rslice[j-count])
                if Qslice[-i] == Rslice[j-count]:
                    matches += 1
                    streak += 1
                else:
                    streak = streak
                count = count + 1
              
            ################
            TEMPA.reverse()
            TEMPB.reverse()
            print("\n\t", ''.join(TEMPA))
            print("\t", ''.join(TEMPB))

            ################        
            score = matches / count
            if streak >= threshold:
                if score > best:
                    best = score
                    if count < 23:
                        best_qind = len(Qslice) - j - 2
                        best_rind = 0
                    else:
                        best_qind = 0
                        best_rind = j - (count - 1)
                    section_score = [best, best_qind, best_rind]

        return section_score

    def score_slices(self, refslice_nos, queryslice_no, threshold=18):

        all_scores = []
        
        # fetch the appropriate query slice. Queryslice_no is the id number, not the pos of the slice.
        query_slice = self.query[self.query_slicehandler[queryslice_no]:self.query_slicehandler[queryslice_no+1]
                                 ]
        # refslice_nose should be a list of refslice id numbers ex: [0, 1, 2, 3]
        # fetch the appropriate ref slice given the "number id" of the slice wanted. the handler gives the indexes
        for slice in refslice_nos:

            ref_slice = self.reference[self.reference_slicehandler[slice]:self.reference_slicehandler[slice+1]]
            print(ref_slice)
            result = self.run_slice(ref_slice, query_slice, threshold)
            if not result:
                continue
            elif result[0] > 0.82:
                all_scores.append([slice, result])
                break
            else:
                all_scores.append([slice, result])

        print("QUERY", query_slice)
        print("REFERENCE", ref_slice)

        return all_scores



# MAIN DRIVER

refpath = './Tools/testa.fasta'
querypath = './Tools/testb.fasta'
tabby = Matrix(refpath, querypath)

tabby.build_Strings()
tabby.build_QSlices()
tabby.build_RSlices()


print("\n////////////\n")

totalslices = len(tabby.reference_slicehandler) - 1
totalqueries = len(tabby.query_slicehandler) - 1

# Handle 0 - k-1 slices
ref_list = []
query_list = []
for i in range(0, (totalslices-2)):
    # THROTTLE CODE
    #if i >= int(totalslices/4):
    #    break
    ref_list.append(i)

for h in range(0, (totalqueries-1)):
    # THROTTLE CODE
    #if i >= int(totalslices/4):
    #    break
    query_list.append(h)

for each in query_list:
    # THROTTLE CODE
    # if each >= 5:
    #     break
    final_score =  tabby.score_slices(ref_list, each)
    tabby.query_scorehandler.append([each, final_score])

for i in tabby.query_scorehandler:
    print("Query", i[0], ": ", i[1])

# Handle k slice (last slice)
