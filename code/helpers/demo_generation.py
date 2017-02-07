import argparse
from lxml import etree
from collections import Counter


# FILE="/Users/bhavyachandra/Desktop/FY/Query/sample.xml"
# FILE="/Users/bhavyachandra/Desktop/FY/Query/dblp.xml"

class DemoParse:

    def __init__(self,file,threshold=5):
        self.file=file
        self.CONFERENCES=["sigmod","vldb","pods","icde","cikm","icdt"]
        self.authors=Counter()
        self.threshold=threshold
        self.current_count=0

    def checkConference(self, booktitle, elem):
        try:
            for conference in self.CONFERENCES:
                if conference in booktitle.text.lower():
                    for auth in elem.findall("author"):
                        self.authors[auth.text]+=1
                        if(self.authors[auth.text]==self.threshold):
                            self.current_count = self.current_count + 1
                    return True
                else:
                    return False
        except AttributeError:
            return False

    def parseXMLDOC(self,count):
            for event, elem in etree.iterparse(self.file,resolve_entities=False,recover=True,remove_blank_text=True):
               if elem.tag != "inproceedings":
                   continue

               booktitle=elem.find("booktitle")
               self.checkConference(booktitle,elem)
               elem.clear()

               if self.current_count >= count:
                   print ("Breaking out")
                   break

def give_output(path,number,min_count):
    demo=DemoParse(path,min_count)
    demo.parseXMLDOC(number)
    print (demo.authors)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fetch the specified number of authors from dblp.xml')
    parser.add_argument('--path', default=None, help='Path to the dblp.xml file')
    parser.add_argument('--number', type=int, default=10, help='The minimum number of authors required')
    parser.add_argument('--min_count', type=int, default=5, help='The minimum articles required to be a valid author')
    give_output(**parser.parse_args().__dict__)