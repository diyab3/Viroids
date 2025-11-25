import sys
import os

usage = 'Usage: ' + sys.argv[0] + ' <lane info file>'

if len(sys.argv) != 2:
    print(usage)
    exit();

lane_info_file = sys.argv[1]

with open(lane_info_file) as F:
    for line in F:
        sra_acc,name = line.strip().split()
        fastq_command = f'mv {sra_acc}.fastq {name}.fastq'
        
        """
        # Create a variable for the original log directory path
        original_dir = f'/home/bio/balakrid/viroids/data/Marquez-Molins_2023/{sra_acc}.log'
        # Create the new name for the log directory
        new_dir_name = "{name}.log"
        # Get the parent directory (of the log directory)
        parent_dir = os.path.dirname(original_dir)
        # Get the base name of the log directory
        current_dir_name = os.path.basename(original_dir)

        # Get a list of all the files in the log directory
        files = os.listdir(original_dir)
        # Sort the list of log files
        files.sort()
        
        print("# File rename commands:")
        # Create the new base name for all of our log files
        new_base_name = f'{name}'
        
        # For each file in our list
        for i, filename in enumerate(files, start=1):
            # Create a variable for the filepath
            old_path = os.path.join(original_dir, filename)

            # store the index of .log in the filename
            idx = filename.find(".log")
            # store everything from .log onwards in the filename
            log_part = filename[idx:]
            # create the new filename: new base name + .log onwards from the old filename
            new_filename = f"{new_base_name}{log_part}"
            
            # get the old filepath by joining the log directory path with the log file's filename
            old_file_path = os.path.join(current_dir_name, filename)
            # get what will be the new filepath by joining the log directory path with the log file's new name
            new_path = os.path.join(current_dir_name, new_filename)
            os.rename(old_path, new_path)
            
            # this will be our linux command for renaming the log file
            print(f"mv '{old_file_path}' '{new_file_path}'")
        
        print("\n# Directory rename command:")
        print(f"mv '{current_dir_name}' '{new_dir_name}'")
        """
        print(fastq_command)

