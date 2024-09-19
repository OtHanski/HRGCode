from tkinter import filedialog

def open_file():
    file_path = filedialog.askopenfilename()
    return file_path

# Read the two column data
def read_data(file_path):
    with open(file_path, 'r') as file:
        data = file.readlines()
    return data

# Multiply the first column by 10
def process_data(data):
    processed_data = []
    for line in data:
        x, y = line.split()
        x = float(x) * 250
        processed_data.append((x, y))
    return processed_data

# Write the processed data to a new file
def write_data(processed_data, filename='processed_data.txt'):
    with open(filename, 'w') as file:
        for x, y in processed_data:
            file.write(f'{x} {y}\n')

# Main function
def main():
    file_path = open_file()
    data = read_data(file_path)
    processed_data = process_data(data)
    write_data(processed_data,filename = file_path)

if __name__ == '__main__':
    main()