import re, os, pandas as pd, concurrent.futures, subprocess, sys
from turtle import down
from datetime import datetime
from subscripts.covipipe_utilities import covipipe_housekeeper as hk
from subscripts.downstream.pipeline_report import copy_files_parallel
from subscripts.src.modules import (
    Module, #base pipeline wrapper class
    module_data, #module metadata 
    pipeline_path #path to pipeline home folder
)



###############################################################
# Defining module classes to extend basic module functionality
###############################################################


class Covid_assembly(Module):
    '''Class extends Module and implements pipeline-specific file processing methods'''

    def __init__(self, *args, **kwargs):
        super(Covid_assembly, self).__init__(*args, **kwargs) #extending parent init method
        self.latest_processing_id = None #latest processing id read from file
        self.processing_id_dict = {} #to map processing ids to original ids for the current batch of samples
        self.renamed_result_file_map = {} #to map file paths annotated with processing ids and with run id removed to the raw file paths
        self.backward_result_file_map = {} #to map file restored annotated paths to specific annotated path
        self.renamed_fastq_file_map = {} #to store original and illumina-formatted fastq file name
        self.backward_fastq_file_map = {} #storing reversed renamed_fastq_file_map to restore original fastq file names after processing

    def fill_input_dict(self, check_integrity:bool=False): #extends Module class
        '''
        Extends fill_input_dict method of the Module class to run integrity check on fastq files before adding them to the input dictionary.
        The integrity check runs by-default. Set check_integrity = False to avoid it.
        '''
        super(Covid_assembly, self).fill_input_dict() #running fill_input_dict of Module class - input_dict is filled by all fastq files found in folders and subfolders of input dir
        
        if check_integrity:
            fastq_total_count = len(self.input_dict["001.fastq.gz"]) #count total fastq found
            print(f"Performing fastq file integrity check:\n{fastq_total_count} total fastq files.") #announcing to the terminal
            failed_log_list = [] #to store data about fastq that failed integrity test
            with concurrent.futures.ProcessPoolExecutor() as executor: #using multiprocessing to speed-up (one check takes about 2.3 sec which translates to about 30 min for 384 nextseq batch using single thread)
                    results = [executor.submit(hk.check_fastq_integrity, fastq) for fastq in self.input_dict["001.fastq.gz"]] #submitting function calls to different processes
                    processed_count = 0 #counter for processed samples
                    for f in concurrent.futures.as_completed(results): #collecting results from different processes
                        result = f.result() #accessing results
                        if not result[0]:
                            failed_log_list.append(f'{result[1]} FAILED\n') #save result as failed to write to log later
                            if 'R1' in result[1]: #if read 1 is corrupted
                                self.input_dict["001.fastq.gz"].remove(result[1]) #remove from further processing
                                self.input_dict["001.fastq.gz"].remove(result[1].replace("R1", "R2"))
                            else: #if read 2 is corrupted
                                self.input_dict["001.fastq.gz"].remove(result[1]) #remove from further processing
                                self.input_dict["001.fastq.gz"].remove(result[1].replace("R2", "R1"))
                        processed_count += 1 #counting processed files
                        hk.printProgressBar(processed_count, fastq_total_count, prefix = 'Progress:', suffix = 'Complete', length = 50)
            with open(f'{self.output_path}input_integrity_check.log', "w+") as logfile: #loggin failed results
                for record in failed_log_list: logfile.write(f"{record}\n")
            print(f"Fastq integrity check complete!\n") #report completion
            if len(failed_log_list) > 0: #if some files have failed
                print(f"Please check {self.output_path}input_integrity_check.log for details on files that failed.")
            else: #if no files failed
                print(f"All fastq files are intact!")
                            


    def map_fastq_to_illumina(self):
        '''Generates illumina-format fastq file names for every fastq.gz file found in self.input_path and stores in self.renamed_fastq_file_map and self.backward_fastq_file_map'''
        fastq_path_list = hk.parse_folder(folder_pth_str=self.input_path, file_fmt_str='_[1,2].fastq.gz')
        fastq_total_count = len(fastq_path_list)
        if fastq_total_count > 0: #if non-illumina formatted samples present
            with concurrent.futures.ProcessPoolExecutor() as executor:
                print(f"Generating illumina format fastq file names:\n{fastq_total_count} total fastq files.") #announcing to the terminal
                results = [executor.submit(hk.name_formatter, path) for path in fastq_path_list] #submitting function calls to different processes
                processed_count = 0 #counter for processed samples
                for f in concurrent.futures.as_completed(results):
                    result = f.result()
                    self.renamed_fastq_file_map[result[0]] = result[1]
                    self.backward_fastq_file_map[result[1]] = result[0]
                    processed_count += 1
                    hk.printProgressBar(processed_count, fastq_total_count, prefix = 'Progress:', suffix = 'Complete', length = 50)
            print(f"Illumina name generation complete!\n") #report completion

    def convert_fastq_names(self, to_illumina:bool=True):
        '''
        Changes fastq file names for every fastq.gz file in self.input_path, using self.renamed_fastq_file_map. 
        If to_illumina set to False, performs renaming using self.backward_fastq_file_map to restore original file names.
        Creates renaming log files in the self.output_path, indicating succesful and failed renaming attempts.
        '''

        if to_illumina: #preprocessing fastq files
            renaming_logs = []
            fastq_total_count = len(self.renamed_fastq_file_map)
            if fastq_total_count > 0:
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    print(f"Renaming fastq files to illumina format names:\n{fastq_total_count} total fastq files.") #announcing to the terminal
                    results = [executor.submit(hk.renamer, path, self.renamed_fastq_file_map[path]) for path in self.renamed_fastq_file_map] #submitting function calls to different processes
                    processed_count = 0

                    for f in concurrent.futures.as_completed(results):
                        result = f.result()
                        renaming_logs.append(result)
                        processed_count += 1
                        hk.printProgressBar(processed_count, fastq_total_count, prefix = 'Progress:', suffix = 'Complete', length = 50)
                print(f"Fastq renaming to Illumina format finished!\n") #report completion
                with open(f'{self.output_path}fastq_forward_renaming.log', "w+") as rename_log: #log renaming process
                    for record in renaming_logs: rename_log.write(record)
                
        else: #restoring names of fastq files to original
            renaming_logs = []
            fastq_total_count = len(self.backward_fastq_file_map)
            if fastq_total_count > 0:
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    print(f"Restoring original (non-illumina) fastq file names:\n{fastq_total_count} total fastq files.") #announcing to the terminal
                    results = [executor.submit(hk.renamer, path, self.backward_fastq_file_map[path]) for path in self.backward_fastq_file_map] #submitting function calls to different processes
                    processed_count = 0

                    for f in concurrent.futures.as_completed(results):
                        result = f.result()
                        renaming_logs.append(result)
                        processed_count += 1
                        hk.printProgressBar(processed_count, fastq_total_count, prefix = 'Progress:', suffix = 'Complete', length = 50)
                print(f"Fastq renaming to original format finished!\n") #report completion
                with open(f'{self.output_path}fastq_backward_renaming.log', "w+") as rename_log: #log renaming process
                    for record in renaming_logs: rename_log.write(record)
                

    def restore_annotated_results(self):
        '''Checks if output directory contains annotations log file. If it is there, reads the contents and restores original name of all files according to the log.'''
        if os.path.isfile(f'{self.output_path}result_annotation.log') and os.stat(f'{self.output_path}result_annotation.log').st_size > 0:
            df = pd.read_table(f'{self.output_path}result_annotation.log', sep=" ", header=None)
            fastq_total_count = len(df)
            processed_count = 0
            with concurrent.futures.ProcessPoolExecutor() as executor:
                print(f"Removing processing id annotation from output file names for reprocessing:\n{fastq_total_count} total output files.") #announcing to the terminal
                results = [] 
                for i, path in enumerate(df[0]):#submitting function calls to different processes
                    self.backward_result_file_map[str(path)] = str(df[1][i])
                    results.append(executor.submit(hk.renamer, str(df[1][i]), str(path)))

                for _ in concurrent.futures.as_completed(results):
                    processed_count += 1
                    hk.printProgressBar(processed_count, fastq_total_count, prefix = 'Progress:', suffix = 'Complete', length = 50)
            print(f"Finished removing processing ids for reprocessing!\n") #report completion


    def compute_processing_ids(self):
        '''Calculates a range of processing ids for current batch of samples and stores it in self.processing_id_dict'''
        self.latest_processing_id = hk.read_plaintext_file(self.config_file['latest_id_file'])[0] #reading latest id
        processing_id_list = [self.latest_processing_id[:3]+str(int(self.latest_processing_id[3:])+i+1).zfill(6) for i in range(len(self.sample_sheet['sample_id']))] #generating processing ids
        self.processing_id_dict = {ids[0]:f'{ids[1]}_{ids[0]}' for ids in zip(list(self.sample_sheet['sample_id']),processing_id_list)} #generating numbered ids
        if not self.backward_result_file_map:
            hk.overwrite_plaintext_file(self.config_file['latest_id_file'],processing_id_list[-1]) #updating processing id after range of ids was used 


    def fill_output_file_maps(self):
        '''Saves original-formatted file path pair in self.renamed_result_file_map'''
        for file_path in self.config_file["assembly_target_files"]: #generating new file path
            for id in self.processing_id_dict:
                if id in file_path:
                    new_path = re.sub(id, self.processing_id_dict[id], file_path) #adding numbered id
                    new_path = re.sub(r"(_S[0-9]{3}|_S[0-9]{2}|_S[0-9]{1})", "", new_path) #removing illumina run sample number
                    if 'qualimap' in file_path: #sample id in directory name and directory should be named
                        old_path = os.path.dirname(file_path)
                        new_path = os.path.dirname(new_path)
                        self.renamed_result_file_map[old_path] = new_path
                    else:
                        self.renamed_result_file_map[file_path] = new_path #storing for forward renaming
                    break


    def annotate_processed_files(self):
        '''Renames processed files using self.renamed_result_file_map and creates a renaming log file in the output directory'''
        if self.backward_result_file_map:
            total_count = len(self.backward_result_file_map)
            processed_count = 0
            annotation_logs = []
            with concurrent.futures.ProcessPoolExecutor() as executor:
                print(f"Reannotating output files with processing ids:\n{total_count} total output files.") #announcing to the terminal
                results = [executor.submit(hk.renamer, path, self.backward_result_file_map[path]) for path in self.backward_result_file_map] #submitting function calls to different processes
                for f in concurrent.futures.as_completed(results):
                    result = f.result()
                    annotation_logs.append(result)
                    processed_count += 1
                    hk.printProgressBar(processed_count, total_count, prefix = 'Progress:', suffix = 'Complete', length = 50)
            print(f"Result annotation finished!\n") #report completion
            with open(f'{self.output_path}result_annotation.log', 'w+') as ann_log: #opening log file for writing  
                for log in annotation_logs: ann_log.write(log)         
        else:
            total_count = len(self.renamed_result_file_map)
            processed_count = 0
            annotation_logs = []
            with concurrent.futures.ProcessPoolExecutor() as executor:
                print(f"Annotating output files with processing ids:\n{total_count} total output files.") #announcing to the terminal
                results = [executor.submit(hk.renamer, path, self.renamed_result_file_map[path]) for path in self.renamed_result_file_map] #submitting function calls to different processes
                for f in concurrent.futures.as_completed(results):
                    result = f.result()
                    annotation_logs.append(result)
                    processed_count += 1
                    hk.printProgressBar(processed_count, total_count, prefix = 'Progress:', suffix = 'Complete', length = 50)

            with open(f'{self.output_path}result_annotation.log', 'w+') as ann_log: #opening log file for writing  
                for log in annotation_logs: ann_log.write(log)


    def switch_sample_ids(self, new_to_old:bool = False):
        '''
        Replaces sample_id column in self.sample_sheet with ids that contain processing number and lack run number.
        If new_to_old is set to True, performs the reverse operation.
        '''
        if new_to_old:
            temp_dict = {re.sub(r"(_S[0-9]{3}|_S[0-9]{2}|_S[0-9]{1})", "", self.processing_id_dict[id]):id for id in self.processing_id_dict} #removing run ids
            self.sample_sheet = hk.map_replace_column(self.sample_sheet, temp_dict, 'sample_id', 'sample_id')
        else:
            temp_dict = {id:re.sub(r"(_S[0-9]{3}|_S[0-9]{2}|_S[0-9]{1})", "", self.processing_id_dict[id]) for id in self.processing_id_dict} #removing run ids
            self.sample_sheet = hk.map_replace_column(self.sample_sheet, temp_dict, 'sample_id', 'sample_id')


class Covid_downstream():
    '''Class implements pipeline-specific downstream processing methods'''

    def __init__(self, start_date:str, end_date:str, skip_pango:bool, skip_heatmap:bool, skip_db_update:bool, skip_tessy:bool, update_share:bool):

        #command-line arguments
        self.start_date = start_date
        self.end_date = end_date
        self.skip_pango = skip_pango
        self.skip_heatmap = skip_heatmap
        self.skip_db_update = skip_db_update
        self.skip_tessy = skip_tessy
        self.update_share = update_share

        #static paths
        self.report_folder_path = "/mnt/beegfs2/home/groups/nmrl/virus_analysis/cov_analysis/reports/"
        self.mutation_file_folder = "/mnt/beegfs2/home/groups/nmrl/virus_analysis/cov_analysis/mutation_heatmap/mutation_files/"
        self.pipe_report_script = "/mnt/beegfs2/home/groups/nmrl/virus_analysis/cov_analysis/SARS-CoV2_assembly/subscripts/downstream/pipeline_report.py"
        self.mutstat_report_script = "/mnt/beegfs2/home/groups/nmrl/virus_analysis/cov_analysis/SARS-CoV2_assembly/subscripts/downstream/mutstat_report.py"
        self.tessy_report_script = "/mnt/beegfs2/home/groups/nmrl/virus_analysis/cov_analysis/SARS-CoV2_assembly/subscripts/downstream/generate_tessy_report.py"
        self.heatmap_update_script = "/mnt/beegfs2/home/groups/nmrl/virus_analysis/cov_analysis/SARS-CoV2_assembly/subscripts/downstream/update_heatmap_job.sh"
        self.weekly_report_script = "/mnt/beegfs2/home/groups/nmrl/virus_analysis/cov_analysis/SARS-CoV2_assembly/subscripts/downstream/run_papermill.sh"
        self.summary_file_script = "/mnt/beegfs2/home/groups/nmrl/virus_analysis/cov_analysis/SARS-CoV2_assembly/subscripts/downstream/update_database_file.py"
        self.update_share_script = "/mnt/beegfs2/home/groups/nmrl/virus_analysis/cov_analysis/SARS-CoV2_assembly/subscripts/downstream/copy_to_c19_share.sh"
        self.log_folder_path = "/mnt/beegfs2/home/groups/nmrl/virus_analysis/cov_analysis/SARS-CoV2_assembly/covipipe_job_logs/"
        self.papermill_path = "/mnt/beegfs2/home/groups/nmrl/virus_analysis/cov_analysis/SARS-CoV2_assembly/tools/rbase_env/bin/papermill"
        self.share_update_file_path = "/mnt/beegfs2/home/groups/nmrl/virus_analysis/cov_analysis/SARS-CoV2_assembly/resources/downstream/c19_share_update.txt"

        #operational variables
        self.current_summary_file = None
        self.pipeline_report_path = None
        self.filtered_report_path = None
        self.mutstat_report_path = None


    def generate_pipeline_report(self):
        '''
        Runs subscripts/downstream/pipeline_report.py script with arguments passed to the main script (date range or sample_id list).
        Saves path to the report file to self.pipeline_report_path and path to the report folder to self.report_folder_path.
        '''
        log_path = f'{self.log_folder_path}{datetime.now().strftime("%Y-%m-%d-%H-%M")}_covipipe_report.log'
        print('Generating pipeline report.')
        try:
            with open(log_path, 'w+') as log_file:
                log_file.write(self.pipe_report_script+" -d1 "+self.start_date+" -d2 "+self.end_date)
                command = [f'{self.pipe_report_script}', "-d1", self.start_date, "-d2", self.end_date]
                if self.skip_pango:
                    print('Pangolin typing skipped, assuming report files exist.')
                    command.append('-s')
                    subprocess.check_call(command, stdout=log_file, stderr=log_file)
                else: 
                    subprocess.check_call(command, stdout=log_file, stderr=log_file)
        except subprocess.CalledProcessError:
            sys.exit(f'Pipeline report generation failed\nSee {log_path} for details')
        self.report_folder_path = f"{self.report_folder_path}report_{datetime.today().date().strftime('%Y-%m-%d')}_{self.start_date}_{self.end_date}"
        self.pipeline_report_path = f"{self.report_folder_path}/pipeline_report.csv"


    def filter_pipeline_report(self):
        '''
        Reads report file from self.pipeline_report_path, filters the dataframe and saves as new file under report_folder_path.
        Saves path to the filtered report file to self.filtered_report_path.
        '''
        print('Filtering pipeline report.')
        df = pd.read_csv(self.pipeline_report_path)
        df = df[(df['genome_N_percentage'] <= 10) & (df['lineage'] != "Unassigned")]
        self.filtered_report_path = f'{self.report_folder_path}/filtered_pipeline_report.csv'
        df.to_csv(self.filtered_report_path, header=True, index=False)


    def generate_mutstat_report(self):
        '''
        Runs subscripts/downstream/mutstat_report.py script with arguments passed to the main script (date range or sample_id list).
        Saves path to the mutstat report to self.mutstat_report_path.
        '''
        print('Generating mutation-by-gene count report.')
        log_path = f'{self.log_folder_path}{datetime.now().strftime("%Y-%m-%d-%H-%M")}_mutstat_report.log'
        command = [self.mutstat_report_script, "-l", self.filtered_report_path, "-o", self.report_folder_path]
        try:
            with open(log_path, 'w+') as log_file:
                log_file.write(" ".join(command))
                subprocess.check_call(command, stdout=log_file, stderr=log_file)
        except subprocess.CalledProcessError:
            sys.exit(f'Mutstat report generation failed\nSee {log_path} for details')
        self.mutstat_report_path = f"{self.report_folder_path}/mutstat_report.csv"


    def update_summary_file(self):
        '''
        Runs subscripts/downstream/update_database_file.py script with arguments passed to the main script (date range or sample_id list).
        '''
        if not self.skip_db_update:
            log_path = f'{self.log_folder_path}{datetime.now().strftime("%Y-%m-%d-%H-%M")}_summary_update.log'
            print(f'Updating summary file.')
            command = [f'{self.summary_file_script}', "-p", self.pipeline_report_path]
            try:
                with open(log_path, 'w+') as log_file:
                    log_file.write(" ".join(command))
                    subprocess.check_call(command, stdout=log_file, stderr=log_file)
            except subprocess.CalledProcessError:
                sys.exit(f'Summary file update failed\nSee {log_path} for details')
        else:
            print(f'Summary file update skipped, assuming the file is up-to-date.')
        self.current_summary_file = [f'/mnt/beegfs2/home/groups/nmrl/virus_analysis/cov_analysis/analysis_history/{file}' for file in os.listdir('/mnt/beegfs2/home/groups/nmrl/virus_analysis/cov_analysis/analysis_history/') if 'summary_file' in file][0]


    def copy_mutation_files(self):
        '''
        Copies all .ann.csv files from self.report_folder_path/source_files/ to self.mutation_file_folder.
        '''
        path_list = [f'{self.report_folder_path}/source_files/{file}' for file in os.listdir(f'{self.report_folder_path}/source_files') if ".ann.csv" in file]
        print(f'Copying annotated variant files to {self.mutation_file_folder}')
        copy_files_parallel(path_list, self.mutation_file_folder, progress_bar=False)


    def update_mut_heatmap(self):
        '''
        Runs subscripts/downstream/update_heatmap_job.sh script with arguments passed to the main script (date range or sample_id list).
        '''
        if not self.skip_heatmap:
            print('Submitting job to HPC to update heatmap plot.')
            log_path = f'{self.log_folder_path}{datetime.now().strftime("%Y-%m-%d-%H-%M")}_heatmap_update.log'
            command = ['qsub', "-o", log_path, "-e", log_path, self.heatmap_update_script]
            try:
                with open(log_path, 'w+') as log_file:
                    log_file.write(" ".join(command))
                    subprocess.check_call(command, stdout=log_file, stderr=log_file)
            except subprocess.CalledProcessError:
                sys.exit(f'Summary file update failed\nSee {log_path} for details')
        else:
            print(f'Heatmap update skipped, assuming the heatmap is up-to-date.')


    def generate_weekly_report(self):
        '''
        Runs subscripts/downstream/report.py script with arguments passed to the main script (date range or sample_id list).
        '''
        print('Generating weekly report.')
        log_path = f'{self.log_folder_path}{datetime.now().strftime("%Y-%m-%d-%H-%M")}_weekly_report.log'
        command = [self.weekly_report_script, self.current_summary_file, self.filtered_report_path, self.mutstat_report_path, self.report_folder_path+"/", self.filtered_report_path, f"{self.report_folder_path}/report_notebook.ipynb"]
        try:
            with open(log_path, 'w+') as log_file:
                log_file.write(" ".join(command))
                subprocess.check_call(command, stdout=log_file, stderr=log_file)
        except subprocess.CalledProcessError:
            sys.exit(f'Weekly report generation failed\nSee {log_path} for details')


    def generate_tessy_report(self):
        '''
        Runs subscripts/downstream/generate_tessy_report.py script with arguments passed to the main script (date range or sample_id list).
        '''
        if not self.skip_tessy:
            print('Generating TESSY report.')
            log_path = f'{self.log_folder_path}{datetime.now().strftime("%Y-%m-%d-%H-%M")}_tessy_report.log'
            command = [self.tessy_report_script, self.filtered_report_path]
            try:
                with open(log_path, 'w+') as log_file:
                    log_file.write(" ".join(command))
                    subprocess.check_call(command, stdout=log_file, stderr=log_file)
            except subprocess.CalledProcessError:
                sys.exit(f'TESSY report generation failed\nSee {log_path} for details')
        else:
            print(f'TESSY report generation skipped.')


    def update_covidshare(self):
        '''
        Updates covidshare by copying bam, fasta and vcf files for all samples included in self.filtered_report_path file.
        '''
        if self.update_share:
            log_path = f'{self.log_folder_path}{datetime.now().strftime("%Y-%m-%d-%H-%M")}_update_share.log'
            os.system(f'''
            echo "Preparing list of samples."
            cat {self.filtered_report_path} | while IFS="," read a b c; do echo ${{b}}_${{a}} ; done > {self.share_update_file_path} 2>> {log_path}
            ''')
            print('Submitting job to HPC to update c19_data_share contents.')
            command = ['qsub', "-o", log_path, "-e", log_path, self.update_share_script]
            try:
                with open(log_path, 'w+') as log_file:
                    log_file.write(" ".join(command))
                    subprocess.check_call(command, stdout=log_file, stderr=log_file)
            except subprocess.CalledProcessError:
                sys.exit(f'Covidshare update failed\nSee {log_path} for details')


################################################
# Defining wrapper functions to call from main
################################################


def run_all(args, num_jobs):
    '''Wrapper function to run all modules sequentially.'''
    assembly = Covid_assembly(
        module_name='assembly',
        input_path=args.input,
        module_config=args.config,
        output_path=args.output_dir,
        run_mode=None,
        dry_run=args.dry_run,
        force_all=args.force_all,
        rule_graph=args.rule_graph,
        pack_output=None,
        unpack_output=None,
        job_name=module_data['assembly']['job_name'],
        patterns=module_data['assembly']['patterns'],
        targets=module_data['assembly']['targets'],
        requests=module_data['assembly']['requests'],
        snakefile_path=module_data['snakefiles']['assembly'],
        cluster_config_path=module_data['cluster_config']
    )

    downstream = Covid_downstream(
        start_date=args.start_date,
        end_date=args.end_date,
        skip_pango = args.skip_pango,
        skip_heatmap = args.skip_heatmap,
        skip_db_update = args.skip_db_update,
        skip_tessy = args.skip_tessy,
        update_share = args.update_share
    )

    if assembly.input_path:
        #Running assembly
        assembly.fill_input_dict()
        assembly.fill_sample_sheet()
        if assembly.unfold_output: assembly.unfold_output()
        assembly.make_output_dir()
        assembly.write_sample_sheet()
        assembly.fill_target_list()
        assembly.add_module_targets()
        assembly.add_output_dir()
        assembly.write_module_config()
        assembly.files_to_wd()
        try:
            assembly.run_module(job_count=num_jobs)
        except Exception as e:
            assembly.clear_working_directory() #to avoid manually moving files back to input
            raise e
        assembly.check_module_output()
        assembly.write_sample_sheet()
        assembly.clear_working_directory()
    else:
        print(f'Skipping assembly as path to a folder containing fastq files (-i argument) was not supplied.')


    if (downstream.start_date and downstream.end_date):
    #Running downstream
        downstream.generate_pipeline_report()
        downstream.filter_pipeline_report()
        downstream.generate_mutstat_report()
        downstream.update_summary_file()
        downstream.copy_mutation_files()
        downstream.update_mut_heatmap()
        downstream.generate_weekly_report()
        downstream.generate_tessy_report()
        downstream.update_covidshare()
    else:
        sys.exit(f'Skipping downstream as start-end dates (-d1; -d2) arguments were not supplied')
    

def run_assembly(args, num_jobs):
    '''Wrapper function to run only assembly module.'''
    assembly = Covid_assembly(
        module_name='assembly',
        input_path=args.input,
        module_config=args.config,
        output_path=args.output_dir,
        run_mode=None,
        dry_run=args.dry_run,
        force_all=args.force_all,
        rule_graph=args.rule_graph,
        pack_output=None,
        unpack_output=None,
        job_name=module_data['assembly']['job_name'],
        patterns=module_data['assembly']['patterns'],
        targets=module_data['assembly']['targets'],
        requests=module_data['assembly']['requests'],
        snakefile_path=module_data['snakefiles']['assembly'],
        cluster_config_path=module_data['cluster_config']
    )
    if assembly.input_path:
        assembly.make_output_dir()
        assembly.map_fastq_to_illumina()
        assembly.convert_fastq_names()
        assembly.restore_annotated_results()
        assembly.fill_input_dict(check_integrity=not args.skip_integrity_check)
        assembly.fill_sample_sheet()
        assembly.write_sample_sheet()
        assembly.fill_target_list()
        assembly.add_module_targets()
        assembly.add_output_dir()
        assembly.write_module_config()
        assembly.files_to_wd()
        try:
            assembly.run_module(job_count=num_jobs)
        except Exception as e:
            assembly.clear_working_directory()
            raise e
        assembly.check_module_output()
        assembly.compute_processing_ids()
        assembly.fill_output_file_maps()
        assembly.annotate_processed_files()
        assembly.switch_sample_ids()
        assembly.write_sample_sheet()
        assembly.clear_working_directory()
        assembly.convert_fastq_names(to_illumina=False)
    else:
        sys.exit(f'Path to a folder containing fastq files (-i argument) must be supplied to run assembly module.')
    
    #Housekeeping
    hk.asign_perm_rec(path_to_folder=assembly.output_path)
    hk.asign_perm_rec(f"{pipeline_path}/covipipe_job_logs/")
    hk.name_job_logs(pipeline_name='covipipe')
    if args.clean_job_logs:
        hk.remove_old_files(f"{pipeline_path}/covipipe_job_logs/")


def run_downstream(args):
    '''Wrapper function to run only downstream module.'''
    downstream = Covid_downstream(
        start_date=args.start_date,
        end_date=args.end_date,
        skip_pango = args.skip_pango,
        skip_heatmap = args.skip_heatmap,
        skip_db_update = args.skip_db_update,
        skip_tessy = args.skip_tessy,
        update_share = args.update_share
    )

    if (downstream.start_date and downstream.end_date):
    #Running downstream
        downstream.generate_pipeline_report()
        downstream.filter_pipeline_report()
        downstream.generate_mutstat_report()
        downstream.update_summary_file()
        downstream.copy_mutation_files()
        downstream.update_mut_heatmap()
        downstream.generate_weekly_report()
        downstream.generate_tessy_report()
        downstream.update_covidshare()
    else:
        sys.exit(f'Skipping downstream as start-end dates (-d1; -d2) arguments were not supplied')
    