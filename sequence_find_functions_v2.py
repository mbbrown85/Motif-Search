import re
import os
import datetime
import json
import sys
from Bio.Seq import Seq

today = datetime.datetime.today()
today_date = f"{today.year}_{today.month}_{today.day}"

def create_folder(work_path):
    if os.path.exists(work_path) == False:
        print("Cannot find program file.")
    input_fold = os.path.join(work_path, "Input_files")
    temp_fold = os.path.join(input_fold, "tmp_files")
    output_fold = os.path.join(work_dir, "Output_files")
    for folder in input_fold, temp_fold, output_fold:
        if os.path.exists(folder):
            break
        else:
            os.makedirs(folder)
        
work_dir = os.path.dirname(os.path.abspath("sequence_find_functions_v2.py"))
create_folder(work_path)

class Meta():
    '''
    Creates an object which contains all the "meta" data for the run, including the motif 
    to be searched for, the distance from the start codon to be searched in, how many 
    instances of the motif are required, and what maximum distances apart from each other 
    they must be.
    '''
    def __init__(self):
        self.meta_dict = {}

    def get_parameters(self):        
        while True:
            motifs = input("Please input your motif.  If multiple, separate with commas:\n").upper()
            mot_list = motifs.split(',')
            check_list = []
            for m in mot_list:
                bases = set(m)
                if all(b in "AGTC" for b in bases) and len(bases) != 0:
                    check_list.append(m)
                else:
                    print("Invalid entry. Only use characters 'A', 'G', 'T', and 'C'.\n")
            break
        self.meta_dict["Motifs"] = check_list
        while True:
            try:
                start_distance = int(input("Please enter the distance from the start codon you wish to search in:\n"))
                self.meta_dict["Start Distance"] = start_distance
                break
            except ValueError:
                print("Please enter a whole number using only digits.\n")
        while True:
            try:
                motif_num = int(input("Please enter the minimum number of motifs present:\n"))
                self.meta_dict["Motif Number"] = motif_num
                break
            except ValueError:
                print("Please enter a whole number using only digits.\n")        
        while True:
            motif_range = input("""Please input the maximum basepair range that the motifs must be in.
                                If multiple ranges are desired, please separate them with commas.\n""")
            motif_list = motif_range.split(',')
            try:
                motif_range_list = [int(i) for i in motif_list]
                self.meta_dict["Motif Range"] = motif_range_list
                break
            except ValueError:
                print("Please use only whole numbers using only digits.\n")
        return self.meta_dict
    

class FileNames():
    def __init__(self, folder):
        self.folder = folder

    def get_file_names(self, type=''):
        self.filenames = []
        for file in os.listdir(self.folder):
            file_list = [file for file in os.listdir(self.folder) if os.path.isfile(os.path.join(self.folder, file)) == True]
        for file in file_list:
            if file.endswith(type):
                self.filenames.append(file)
            else:
                print(f"{file} is not a valid {type} file.")
        return self.filenames
    

class ReadFile():
    def __init__(self, filepath):
        self.path = filepath

    def read_file(self):
        with open(self.path) as f:
            self.contents = f.readlines()
        return self.contents

    def create_seq_dict(self):
        self.seq_dict = {}
        for line in self.contents:
            if line.startswith('>'):
                seq_name = line        
                self.seq_dict[seq_name] = []
            else:
                self.seq_dict[seq_name].append(line.rstrip('\n'))
        for key, value in self.seq_dict.items():
            self.seq_dict[key] = ''.join(self.seq_dict[key])
        return self.seq_dict
    
    def read_json(self):
        with open(self.path) as f:
            json_load = json.load(f)
            #print(json_load)
        return json_load
    

class WriteFile():
    def __init__(self, filepath):
        self.path = filepath
    
    def write_seq_file(self, key, value):
        with open(self.path, 'w') as wf:
            wf.write(key)
            wf.write(value)        
    
    def write_summary_file(self, meta_dict, results_dict):
        with open(self.path, 'w') as wf:
            results = list()
            results.append(meta_dict)
            results.append(results_dict)
            result_json = json.dumps(results)
            wf.write(result_json)


class ReadWrite():
    def __init__(self, input_fold, output_fold):
        self.input_fold = input_fold
        self.output_fold = output_fold
    
    def create_output_folder(self, today_date):
        i = 1
        while True:
            try:
                self.today_path = os.path.join(self.output_fold, f"{today_date}_{i}")
                os.mkdir(self.today_path)
                break
            except FileExistsError:
                i += 1
        return self.today_path

    def clear_tmp_folder(self):
        for del_file in os.listdir(temp_fold):
            try:
                del_path = os.path.join(temp_fold, del_file)
                os.remove(del_path)
            except:
                print(f"Error: {sys.exc_info}")

    def split_seq_files(self):
        filelist = FileNames(self.input_fold).get_file_names()
        for file in filelist:
            path = os.path.join(self.input_fold, file)
            rfile = ReadFile(path)
            rfile.read_file()
            seq_dict = rfile.create_seq_dict()
            for key, value in seq_dict.items():
                name = key.replace(':','_').replace(' ','_').replace('>','').replace('\n','')
                wpath = os.path.join(self.output_fold, name+'.txt')
                if os.path.exists(wpath) == False:
                    wfile = WriteFile(wpath)
                    wfile.write_seq_file(key, value)
                else:
                    os.remove(wpath)
                    wfile = WriteFile(wpath)
                    wfile.write_seq_file(key, value)    

    def create_reverse_files(self):
        filelist = FileNames(self.output_fold).get_file_names()
        for file in filelist:
            path = os.path.join(self.output_fold, file)
            rfile = ReadFile(path)
            rfile.read_file()
            seq_dict = rfile.create_seq_dict()
            sequence = Seq(list(seq_dict.values())[0])
            rev_seq = sequence.reverse_complement()
            new_name = file.strip('.txt')+"_reverse.txt"
            for key, value in seq_dict.items():
                wpath = os.path.join(self.output_fold, new_name)
                wfile = WriteFile(wpath)
                wfile.write_seq_file(key, str(rev_seq))

    def screen_seq_files(self, meta_dict):
        filelist = FileNames(self.input_fold).get_file_names()
        for file in filelist:
            path = os.path.join(self.input_fold, file)
            rfile = ReadFile(path)
            rfile.read_file()
            seq_dict = rfile.create_seq_dict()
            seq_name = list(seq_dict.keys())[0]
            sequence = list(seq_dict.values())[0]
            current_seq = Sequence(seq_name, sequence, meta_dict, file)
            current_seq.find_start()
            current_seq.find_motif_loc()
            results_dict = current_seq.screen_summary_dict()
               
            wpath = os.path.join(self.today_path, file)
            wfile = WriteFile(wpath)
            wfile.write_summary_file(meta_dict, results_dict)

    def collate_summary_files(self, meta_dict):
        filelist = FileNames(self.input_fold).get_file_names()
        range_count_dict = {}
        range_sum_dict = {}
        for mrange in meta_dict['Motif Range']:
            range_count_dict[mrange] = 0
            range_sum_dict[mrange] = []
        for file in filelist:
            path = os.path.join(self.input_fold, file)
            rfile = ReadFile(path)
            json_load = rfile.read_json()[1]
            for name, data in json_load.items():
                for mrange, motifs in data.items():
                    #print(f"The number of motifs for {name} in a range of {mrange} is {len(count)}")
                    mrange_key = int(mrange)
                    range_count_dict[mrange_key] = range_count_dict[mrange_key]+len(motifs)
                    if motifs is not {}:
                        range_sum_dict[mrange_key].append(motifs)
        print(range_count_dict)
        path = os.path.join(self.input_fold, "_Summary_file.txt")
        wfile = WriteFile(path)
        wfile.write_summary_file(meta_dict, range_sum_dict)
                        

class Sequence():
    '''
    Takes a sequence from a file in FASTA format and identifies every potential start codon (ATG)
    present within the sequence. Currently does not differentiate between real start codons and
    ATG which code for Met. For each ATG, then finds all the selected motifs with a certain 
    distance upstream of the start. Then loops through the list of motif locations and determines
    whether the selected number of motifs is present within the selected basepair range.  Returns
    a multilevel dictionary containing the sequence name, the basepair range, the location of a 
    given start codon, the locations of any motifs which meet the selected criteria, and the DNA 
    sequence 1500 bp downstream of the ATG.    
    '''
    def __init__(self, seq_name, sequence, meta_dict, filename):
        self.seq_name = seq_name
        self.sequence = sequence
        self.meta_dict = meta_dict
        self.filename = filename

    def find_start(self):
        seq_iter = re.finditer("ATG", self.sequence)
        self.starts = [m.start(0) for m in seq_iter]
        
    def find_motif_loc(self):
        self.motif_dict = {}
        if self.starts == "No ATG sequences found.":
            pass
        else:
            self.motif_dict[self.seq_name] = {}
            for motif_range in self.meta_dict["Motif Range"]:
                self.motif_dict[self.seq_name][motif_range] = {}
                for atg_loc in self.starts:
                    region = self.sequence[atg_loc-self.meta_dict['Start Distance']:atg_loc]
                    gene_seq = self.sequence[atg_loc:atg_loc+1500]
                    if region == []:
                        motif_loc = "The ATG is too close to the start of the sequence."
                        self.motif_dict[self.seq_name][motif_range][atg_loc] = motif_loc
                    else:
                        patterns = self.meta_dict["Motifs"]
                        re_patterns = '|'.join(re.escape(pattern) for pattern in patterns)
                        comp_pattern = re.compile(re_patterns)
                        motifs = re.findall(comp_pattern, region)
                        if len(motifs) == 0:
                            motif_loc = "No motifs found."
                            self.motif_dict[self.seq_name][motif_range][atg_loc] = motif_loc
                        elif len(motifs) < self.meta_dict["Motif Number"]:
                            motif_loc = "Number of motifs found lower than desired criteria."
                            self.motif_dict[self.seq_name][motif_range][atg_loc] = motif_loc
                        elif len(motifs) >= self.meta_dict["Motif Number"]:
                            motif_loc = [m.start(0) for m in re.finditer(comp_pattern, region)]
                            #self.motif_dict[motif_range][atg_loc] = []
                            for i in range(self.meta_dict['Motif Number']):
                                try:
                                    start_motif = motif_loc[i]
                                    end_motif = motif_loc[i+self.meta_dict["Motif Number"]-1]
                                    distance = end_motif-start_motif
                                    #print(f"Distance is {distance}.")
                                    if distance <= motif_range:
                                        motifs_within_range = []
                                        for j in range(self.meta_dict["Motif Number"]):
                                            motifs_within_range.append(motif_loc[i+j])
                                        motif_set = list(set(motifs_within_range))
                                        self.motif_dict[self.seq_name][motif_range][atg_loc] = [motif_set, gene_seq]
                                        #self.motif_dict[motif_range][atg_loc].append(gene_seq)
                                        #print(self.motif_dict[self.seq_name][motif_range][atg_loc])
                                        #print(type(self.motif_dict[self.seq_name][motif_range][atg_loc]))
                                    else:
                                        break
                                except IndexError:
                                    break     

    def screen_summary_dict(self):
        screen_dict = {}
        for seq_name, motif_ranges in self.motif_dict.items():
            screen_dict[seq_name] = {}
            for m_range, atg_locs in motif_ranges.items():
                screen_dict[seq_name][m_range] = {}
                for atg_loc, motifs in atg_locs.items():
                    if type(motifs[0]) != str:
                        screen_dict[seq_name][m_range][atg_loc] = motifs
                    else:
                        pass
        #print(f"Screen dictionary is {screen_dict}")
        results_dict = {}
        for seq_name, motif_range in screen_dict.items():
            #print(f"Seq is {seq_name}")
            results_dict[seq_name] = {}
            for m_range, atg_locs in motif_range.items():
                seen = set()
                results_dict[seq_name][m_range] = {}
                #print(f"Range is {m_range}")
                for atg_loc, motifs_seq in atg_locs.items():
                    #print(f"Motifs for {atg_loc} are at {motifs_seq}")
                    try:
                        motifs = motifs_seq[0]
                        gene_seq = motifs_seq[1]
                        adj_motifs = []
                        for i in motifs:
                            adj_value = int(atg_loc)-(meta_dict['Start Distance']-i)
                            adj_motifs.append(adj_value)
                        adj_tuple = tuple(sorted(adj_motifs))
                        if adj_tuple in seen:
                            #print(f"Motif at locations {adj_tuple} already seen.")
                            pass
                        else:
                            results_dict[seq_name][m_range][atg_loc] = [adj_motifs, gene_seq]
                            #print(f"Motifs seen for {atg_loc} at {adj_motifs}")
                            #print(results_dict[seq_name][m_range][atg_loc])
                            seen.add(adj_tuple)  
                    except:
                        print("Error:", sys.exc_info())
        #print(f"Items of Results dict are {results_dict.items()}")
        return results_dict


if __name__ == "__main__":
    meta_dict = Meta().get_parameters()
    seq_files = ReadWrite(input_fold, temp_fold)
    seq_files.clear_tmp_folder()
    seq_files.split_seq_files()
    seq_files.create_reverse_files()
    temp_files = ReadWrite(temp_fold, output_fold)
    today_path = temp_files.create_output_folder(today_date)
    temp_files.screen_seq_files(meta_dict)
    output_files = ReadWrite(today_path, today_path)
    output_files.collate_summary_files(meta_dict)
