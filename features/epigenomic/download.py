import os
import subprocess
import sys

def run_script_in_folders(folders):
    for folder in folders:
        # Get absolute path of the folder
        abs_folder = os.path.abspath(folder)
        
        # Check if the folder exists
        if os.path.isdir(abs_folder):
            script_path = os.path.join(abs_folder, 'download_bigbed.sh')
            if os.path.isfile(script_path):
                print(f"Running {script_path} in folder: {abs_folder}")
                try:
                    # Change the current working directory to the folder
                    os.chdir(abs_folder)
                    
                    # Run the shell script using subprocess
                    subprocess.run(['bash', script_path], check=True)
                    print(f"Successfully ran {script_path}")
                except subprocess.CalledProcessError as e:
                    print(f"Error occurred while running {script_path}: {e}")
                except Exception as e:
                    print(f"An unexpected error occurred: {e}")
            else:
                print(f"{script_path} not found in folder: {abs_folder}")
        else:
            print(f"Folder {abs_folder} does not exist")

# Example usage
if __name__ == "__main__":
    # Check if exactly one command line argument is provided
    if len(sys.argv) != 2:
        print("Usage: python script.py <base_directory>")
        sys.exit(1)

    # Get the base directory from command line arguments
    base_directory = sys.argv[1]
    # Ensure the provided path is an absolute path
    base_directory = os.path.abspath(base_directory)

    # Collect all subdirectories in the base directory into a list
    folders_to_process = [os.path.join(base_directory, name) for name in os.listdir(base_directory) if os.path.isdir(os.path.join(base_directory, name))]

    # Call the function with the list of folders to process
    run_script_in_folders(folders_to_process)
