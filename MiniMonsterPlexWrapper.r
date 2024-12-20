# Load the glue library for string interpolation
library(glue)

# Prompt the user to enter the output folder
output.folder =""
while (output.folder== ""){
	output.folder = readline(prompt = "Enter the output folder. The folder must not exist: ")

	if (dir.exists(output.folder)) {
		print("The folder already exists. Please enter a different folder.")
		output.folder = ""
	}
}

# Prompt the user to enter the metadata.csv file
metadata.file = ""
while (metadata.file == "") {
	metadata.file = readline(prompt="Enter the metadata.csv file: ") 
	if (metadata.file == "") {
		print("The metadata file is required.")
	}
}
# Prompt the user to enter the input folder
input.folder = ""
while (input.folder == "") {
	input.folder = readline(prompt="Enter the input folder: ")

	if (input.folder == "") {
		print("The input folder is required.")
	}
}
# Prompt the user to enter the space-separated isolate list (optional)
isolate.list = readline(prompt="Enter the space separated isolate list. Optional push enter if you want to skip: ")

# If the isolate list is empty, set it to NULL; otherwise, split it into a list
if (isolate.list == "") {
	isolate.list = NULL
} else {
	isolate.list = unlist(strsplit(isolate.list, " "))
}

# Prompt the user to enter the isolate list file (optional)
isolate.file = readline(prompt="Enter the isolate list file. Optional push enter if you want to skip: ")

# If the isolate file is empty, set it to NULL
if (isolate.file == "") {
isolate.file = NULL
}

# Prompt the user to enter the space-separated host list (optional)
host.list = readline(prompt="Enter the space seperated host list. Optional push enter if you want to skip: ")

# If the host list is empty, set it to NULL; otherwise, split it into a list
if (host.list == "") {
host.list = NULL
} else {
host.list = unlist(strsplit(host.list, " "))
}

# Prompt the user to enter the host list file (optional)
host.file = readline(prompt="Enter the host list file. Optional push enter if you want to skip: ")

# If the host file is empty, set it to NULL
if (host.file == "") {
	host.file = NULL
}

# Construct the command options based on the provided isolate and host lists/files
if (is.null(isolate.list)) {
	isolate.list.command = " "
} else {
	isolate.list.command = glue(' -i {isolate.list}')
}
if (is.null(isolate.file)) {
	isolate.file.command = " "
} else {
	isolate.file.command = glue(' -il {isolate.file}')
}
if (is.null(host.list)) {
	host.list.command = " "
} else {
	host.list.command = glue(' -hf {host.list}')
}
if (is.null(host.file)) {
	host.file.command = " "
} else {
	host.file.command = glue(' -hfl {host.file}')
}

# Construct the final command to run the MiniMonsterPlex.py script with the provided options
command = glue('python MiniMonsterPlex.py -o {output.folder} -m {metadata.file} -f {input.folder} {isolate.list.command} {isolate.file.command} {host.list.command} {host.file.command}')

# Execute the constructed command
system(command)
