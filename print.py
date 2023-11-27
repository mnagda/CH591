import subprocess
import numpy as np

# List of values for AA
AA_values = np.linspace(0.01, 0.025, 20)  # Add more values as needed

# Empty lists to store data for each file
check_data = []
test_data = []
sample_data = []

for AA in AA_values:
    # Modify the user_input.py file with the current AA value
    with open('user_input.py', 'r') as file:
        data = file.readlines()

    # Update the AA value
    for i in range(len(data)):
        if data[i].startswith('AA ='):
            data[i] = f'AA = {AA}\n'

    # Write the modified data back to user_input.py
    with open('user_input.py', 'w') as file:
        file.writelines(data)

    # Run each .py file and capture the output
    files_to_run = ['test.py']  # List of .py files
    for file_to_run, data_list in zip(files_to_run, [check_data, test_data, sample_data]):
        output = subprocess.check_output(['python', file_to_run], universal_newlines=True)

        # Process the output data and append to respective lists
        output_lines = output.strip().split('\n')
        output_values = [float(val) for val in output_lines[-1].split()]  # Assuming the last line contains c1 c2 c3
        data_list.append([AA] + output_values)

# Write the collected data to separate files for each script
def write_data_to_file(file_name, data):
    with open(file_name, 'w') as file:
        for row in data:
            file.write(' '.join(map(str, row)) + '\n')

#write_data_to_file('check_data.txt', check_data)
write_data_to_file('test_data.txt', test_data)
#write_data_to_file('sample_data.txt', sample_data)
