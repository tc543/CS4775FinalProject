
input_file_path = 'human_chromosome.txt'  
output_file_path = 'extracted_content.txt' 

try:
    with open(input_file_path, 'r', encoding='utf-8') as infile:
        content = infile.read(60000)  
    
    with open(output_file_path, 'w', encoding='utf-8') as outfile:
        outfile.write(content)


except FileNotFoundError:
    print(f"wrong file name")
except Exception as e:
    print(f"bruh what?")
