import math
import os
import json
import matplotlib


class Matrix():

    def __init__(self, ref_seqpath, fas_seqpath, r_win = 2000, q_win = 24) -> None:
        
        self.reference_path = ref_seqpath
        self.query_path = fas_seqpath
        self.reference = ''
        self.query = ''

        self.ref_window = r_win
        self.query_window = q_win

        # default 2000 bases per slice. (EX: 5,000,000 bp ref / 2,000 = 2,500 slices)
        # slice handler index = slice id
        # The index of a slice in the mother is ((slice number - 1) * 2000) - 1
        self.reference_slicehandler = []

        # index matches index of slice in slicehandler
        # value is the score of the matching slice
        self.query_scorehandler = []

        # default 24 bases per slice
        # slice's position in query string = ((slice number - 1) * 24) - 1
        # slice id = index, slice value = pos in the mother 
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
    
    def run_slice(self, Rid, Qid, Rslice, Qslice, threshold):
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
            #print("\n\t", ''.join(TEMPA))
            #print("\t", ''.join(TEMPB))

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
                    # final qind = the start pos of the query slice within the query
                    # final rind = the start pos of the match location with the reference
                    final_qind = self.query_slicehandler[Qid] + best_qind
                   
                    final_rind = self.reference_slicehandler[Rid] + best_rind - 1

                    section_score = [best, final_qind, final_rind]

        return section_score

    def score_slices(self, refslice_nos, queryslice_no, threshold=18):

        all_scores = [0, [0, 0, 0]]
        
        # fetch the appropriate query slice. Queryslice_no is the id number, not the pos of the slice.
        query_slice = self.query[self.query_slicehandler[queryslice_no]:self.query_slicehandler[queryslice_no+1]]


        # refslice_nose should be a list of refslice id numbers ex: [0, 1, 2, 3]
        # fetch the appropriate ref slice given the "number id" of the slice wanted. the handler gives the indexes
        for slice in refslice_nos:
            ref_slice = self.reference[self.reference_slicehandler[slice]:self.reference_slicehandler[slice+1]]
            result = self.run_slice(slice, queryslice_no, ref_slice, query_slice, threshold)
            if not result:
                continue

            # And record your best scores
            elif result[0] > 0.97:
                all_scores = [slice, result]
                break
            else:
                if all_scores[1][0] < result[0]:
                    all_scores = [slice, result]

        #print("QUERY", query_slice)
        #print("REFERENCE", ref_slice)


        return all_scores
    

    # Gaps: [Gap Id, Internal/Sister (0, 1), gap_length,
            # [gap_start pos on ref (first mismatch), query_start_ind_of_diff]
                # [0, 0, 25, [964, 125]]
    def seek_internalgap(self, query_id, gap_id):
        query_select = self.query_scorehandler[query_id]
        query_string = self.query[query_select[1][1][1]:query_select[1][1][1]+24]
        ref_string = self.reference[query_select[1][1][2]:query_select[1][1][2]+24]
        ref_breaks = []
        breaks = []
        for i in range(0, len(query_string)):
            if query_string[i] == ref_string[i]:
                if not breaks:
                    continue
                else:
                    ref_breaks.append(breaks)
                    breaks = []
            else:
                breaks.append(i)

        for each in ref_breaks:
            ref_gap_ind = query_select[1][1][2] + each[-1]
            query_gap_ind = query_select[1][1][1] + each[-1]
            self.ref_gaps.append([gap_id, 0, len(each), [ref_gap_ind, query_gap_ind]])
            gap_id += 1

        return gap_id
    
    def seek_sistergap(self, query_id, gap_id):
        query_select = self.query_scorehandler[query_id]
        sister_select = self.query_scorehandler[query_id+1]

        break_start = query_select[1][1][2] + (self.query_window - query_select[1][0])
        break_end = sister_select[1][1][2] - 1
        print("Q", query_id, "From ", break_start, "to Q", query_id+1, break_end)
        gap_length = int(break_end-break_start)+1
        if gap_length <= 0:
            # ERROR!
            return gap_id
        self.ref_gaps.append([gap_id, 1, gap_length, [break_start, break_end]])
        gap_id += 1

        return gap_id
    
    def seek_allgaps(self):
        gap_id = 0
        for query in range(0, len(self.query_scorehandler)-1):
            if (self.query_scorehandler[query][1][1][0] < 1.0):
                gap_id = self.seek_internalgap(query, gap_id)
            if (self.query_scorehandler[query][1][1][2]+24) != self.query_scorehandler[query+1][1][1][2]:
                gap_id = self.seek_sistergap(query, gap_id)
                continue
            else:
                continue

            

        return
    
    def consolidate_gaps(self):
        groups = []
        for gap in range(0, len(self.ref_gaps)-1):
            overlap = self.ref_gaps[gap][3][0] + self.ref_gaps[gap][2]
            if overlap >= self.ref_gaps[gap+1][3][0]-1:
                gap_length = (self.ref_gaps[gap+1][3][0] - self.ref_gaps[gap][3][0]) + self.ref_gaps[gap+1][2]
                grouped_gap = [self.ref_gaps[gap][0], 3, gap_length, [self.ref_gaps[gap][3][0], self.ref_gaps[gap][3][0]+gap_length]]
                groups.append(grouped_gap)
        for id in groups:
            for par in self.ref_gaps:
               
                if not par:
                    continue
                elif id[0] == par[0]:
                    self.ref_gaps[self.ref_gaps.index(par)+1] = []
                    self.ref_gaps[self.ref_gaps.index(par)] = id
                    break
        self.ref_gaps.remove([])

        return

    def build_table(self):
        #each entry in the query_scorehandler should look like:
        # [ query id no, 
        # [[ref_slice id no, 
        # [match score (1.0 means perfect match), 
                # query_slice_start_ind (0 means whole query slice was aligned), 
                        # ref_slice_start_ind (the pos where the match alignment starts) ]]]
        # [0, [[5, [1.0, 0, 515]]]]
        

        """ We need to use each query in the scorehandler to "fill the puzzle" 
            Imagine the reference as a long sequence. Using the markers provided by the scorehandler, 
            make a table that highlights where the whole query matches the reference. If a query does 
            not have a perfect 1.0 score, there is a gap present. Could be an intron, could be a 
            polymorphism. We'll need to investigate deeper for those queries.
            If a query and its neighbor are not aligned to each other (the end of the query is the start
            of the neighbor) then there is a gap. Same issue as with those imperfect scores.
            -> a gap of larger than like 12 bases is probably not a polymorphism. For this program 
                we elect to ignore polymorphisms larger than SNPs (1 base)
        """



        return 0
    
######################################################################
#################         MAIN DRIVER            #####################
######################################################################

refpath = './Tools/testa.fasta'
querypath = './Tools/testb.fasta'



tabby = Matrix(refpath, querypath)

tabby.build_Strings()
tabby.build_QSlices()
tabby.build_RSlices()


totalslices = len(tabby.reference_slicehandler) - 1
totalqueries = len(tabby.query_slicehandler) - 1

# Handle 0 - k-1 slices
ref_list = []
query_list = []
for i in range(0, (totalslices)):
    # THROTTLE CODE
    #if i >= int(totalslices/4):
    #    break
    ref_list.append(i)

for h in range(0, (totalqueries)):
    # THROTTLE CODE
    #if i >= int(totalslices/4):
    #    break
    query_list.append(h)

for each in query_list:
    # THROTTLE CODE
    # if each >= 20:
    #    break
   final_score =  tabby.score_slices(ref_list, each)
   tabby.query_scorehandler.append([each, final_score])
   

for i in tabby.query_scorehandler:    
    print("Query", i[0], ": ", i[1])

calico = tabby.query_scorehandler[15]
print(calico)
TEMPC = tabby.query[calico[1][1][1]:calico[1][1][1]+24]
print(TEMPC)
print(tabby.reference[calico[1][1][2]:calico[1][1][2]+24])

tabby.seek_allgaps()
tabby.consolidate_gaps()

print(tabby.ref_gaps)
#tabby.seek_internalgap(33, 0)

# Handle k slice (last slice)