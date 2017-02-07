from collections import OrderedDict

class dblp_graph_label_translation:

    def __init__(self,color_file, author_file):
        self.color_file_path=color_file
        self.author_file_path=author_file
        self.author_list=dict()
        self.original=[]
        self.color_list=dict()
        self.author_number_list=dict()
        self.color_number_list=dict()

    def author_to_number(self,dict_of_author):
        """
        creates dictionary containing mapping of author to number

        :param dict_of_author: A dict of author names vs value
        :return: nothing
        """
        self.original=sorted(dict_of_author.items(), key=lambda x:x[1], reverse=True)
        list_of_author=self._translate_author_dict(dict_of_author)
        new_author_list=self._translate_names(list_of_author).split(" ")
        for line in open(self.author_file_path,mode='r'):
            entry=line.split(" ")
            if entry[0].lower() in new_author_list:
                self.author_list[entry[0]]=int(entry[1])
            if len(self.author_list) == len(new_author_list):
                break
        print ("Unable to find the following:")
        temp3 = [x for x in new_author_list for y in self.author_list if x.lower()==y.lower()]
        print (len(set(new_author_list)-set(temp3)))
        print (set(new_author_list)-set(temp3))

    def number_to_author(self,list_of_numbers,sep=" "):
        """
        creates dictionary containing mapping of number to author

        :param list_of_numbers: A list of numbers sepearated by space or comma
        :return: dictionary
        """
        numbers = list_of_numbers.split(sep)
        for line in open(self.author_file_path,mode='r'):
            entry=line.split(" ")
            if entry[1].strip() in numbers:
                self.author_number_list[int(entry[1].strip())]=entry[0]
            if len(self.author_number_list) == len(numbers):
                break

    def _translate_names(self,list_of_author):
        """

        :param list_of_author:A list of author names
        :return: modified author names with underscore insertion
        """
        list_of_new_names=""
        for full_name in list_of_author:
            name_part=full_name.strip().split(" ")
            new_full_name=""
            for name in name_part:
                new_full_name+=name.strip()+"_"
            list_of_new_names+=new_full_name[:-1].lower()+" "
        return list_of_new_names

    def _translate_author_dict(self,auth_dict):
        list_of_auth=[]
        for key in auth_dict.keys():
            list_of_auth.append(key)
        return list_of_auth


    def color_to_number(self,list_of_color):
        """
        creates dictionary containing mapping of author to number
        :param list_of_color: A string of colour names, seperated by commas
        :return: nothing
        """
        new_color_list=list_of_color.split(" ")
        for line in open(self.color_file_path,mode='r'):
            entry=line.split(" ")
            if entry[0].lower() in new_color_list:
                self.color_list[entry[0]]=int(entry[1])
            if len(self.color_list) == len(new_color_list):
                break

    def number_to_color(self,list_of_numbers,sep=" "):
        """

        :param list_of_numbers: A list of numbers sepearated by space or comma
        :return: dictionary
        """
        numbers = list_of_numbers.split(sep)
        for line in open(self.color_file_path,mode='r'):
            entry=line.split(" ")
            if entry[1].strip() in numbers:
                self.color_number_list[int(entry[1].strip())]=entry[0]
            if len(self.color_number_list) == len(numbers):
                break

    def pretty_print(self,auth_to_num=True,num_to_auth=False,num_to_color=True,color_to_num=False):
        # author_list
        if auth_to_num:
            count=1
            pos_str=""
            for (key,value) in self.original:
                key_t=key.split(" ")
                temp=""
                for k in key_t:
                    temp+=k+"_"
                if temp[:-1] in self.author_list.keys():
                    print str(count)+". "+key+":"+str(self.author_list[temp[:-1]])
                    count+=1
                    pos_str+=str(self.author_list[temp[:-1]])+" "

            print ("List of only numbers")
            print (pos_str)

        # author_number_list
        if num_to_auth:
            for key,value in self.author_number_list.items():
                print str(key)+":"+value

        # num_to_color
        if num_to_color:
            for key,value in self.color_number_list.items():
                print value+":"+str(key)


        # color list
        if color_to_num:
            for key,value in self.color_list.items():
                print key+":"+str(value)

PATHC= # Add color file path
PATH= # Add author file path
lok=dblp_graph_label_translation (PATHC,PATH)

# lok.author_to_number()
# print(lok.author_list)
# lok.number_to_author("8737 16181 16527 16220 88415")
print(lok.author_number_list)
lok.number_to_color("197 64 31 43 215 3 19 279 63 113")
# print (lok.color_number_list)
# lok.color_to_number("datacenter system")
# print(lok.color_list)
lok.pretty_print()

