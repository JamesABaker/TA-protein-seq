# p005methods.rb
# gets and chomp


#For custom files
'''
puts "Name your first file?"
STDOUT.flush
firstfile = gets.chomp
'''

#In the pipeline

firstfile = "C_terminal_single_TRANSMEM.txt"
secondfile = "S1_kalbfleisch_list_uniprot_excluding_predicted_proteins.txt"

firstlist = File.open(firstfile)

f_lines = firstlist.read.split("\n")

puts "These are the IDs in the first file."
puts f_lines

a = f_lines

secondlist = File.open(secondfile)

g_lines = secondlist.read.split("\n")

b = g_lines

puts "These are the IDs in the first file."
puts g_lines

puts "Here is some sort of sorting..."

new_protein_list = a - b

puts new_protein_list
