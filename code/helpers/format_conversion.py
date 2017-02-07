file = # add path
new_file=file.replace(".txt","_part1.txt")
write_file=open(new_file,mode="w")
for line in open(file, mode="r"):
    items = line.split(" ")
    x = 2
    while x < len(items)-1:
        write_file.write(items[0]+" "+items[1]+" "+items[x]+" "+items[x+1]+"\n")
        x += 2

write_file.close()
