import argparse

NEW_FILE="corrected_graph.txt"

def give_output(file, cpath):
    colorList=[]
    for line in open(cpath,mode="r"):
        colorList=line.strip().split(" ")
        print "Color list:",
        print colorList
    write_file=open(NEW_FILE,mode="w")
    for line in open(file, mode="r"):
        items = line.strip().split(" ")
        x = 2
        probNotTaking=1
        for x in range(2,len(items)-1,2):
            if not len(colorList)==0:
                if items[x+1] in colorList:
                    probNotTaking *= (1-float(items[x]))
            else:
                 probNotTaking *= (1-float(items[x]))
        probNotTaking=1-probNotTaking
        write_file.write(items[0]+" "+items[1]+" "+str(probNotTaking)+"\n")

    write_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fetch the specified number of authors from dblp.xml')
    parser.add_argument('--file', default=None, help='Path to the dblp.xml file')
    parser.add_argument('--cpath', default=None, help='Path to the color file')
    give_output(**parser.parse_args().__dict__)