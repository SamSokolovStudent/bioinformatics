import re
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from tkinter.messagebox import showinfo

# Regex consensus patterns stored as global var as they remain unchanged
# and regex should not be recompiled on every use.
his_pattern = \
    re.compile(".[ST]G[LIVMFYW]{3}[GN][A-Z]{2}T[LIVM][A-Z]T[A-Z]{2}H.")
ser_pattern = \
    re.compile(".T[A-Z]{2}[GC][NQ]SGS[A-Z][LIVM][FY].")


class GeneBankData:
    """GeneBankData Class is used to store a dictionary/map
    to automatically retrieve and assign any data in a modular way in
    the init of the object."""

    def __init__(self, dictionary_map):
        self.__protein_id = dictionary_map["protein_id"]
        self.__protein_sequence = dictionary_map["translation"]
        self.__protein_product = dictionary_map["product"]
        self.__cds_range = dictionary_map["cds_range"]
        self.__matched_pattern = dictionary_map["matched_pattern"]
        self.__db_xref = dictionary_map["db_xref"]
        self.__strand = None

    def get_protein_id(self):
        """Returns the protein ID from the object.
        :return __protein_id"""
        return self.__protein_id

    def get_protein_sequence(self):
        """Returns the protein ID from the object.
        :return __protein_sequence"""
        return self.__protein_sequence

    def get_protein_product(self):
        """Returns the protein ID from the object.
        :return __protein_product"""
        return self.__protein_product

    def get_cds_range(self):
        """Returns the protein ID from the object.
        :return __cds_range"""
        return self.__cds_range

    def get_strand(self):
        """Returns the strand from the object.
        :return __strand"""
        return self.__strand

    def get_pattern(self):
        """Returns the pattern from the object.
        :return __pattern"""
        return self.__matched_pattern

    def set_strand(self, strand):
        """Determines the type of the strand and sets the correct value.
        :param strand"""
        if strand == "+":
            self.__strand = "Forward Strand"
        elif strand == "-":
            self.__strand = "Reverse Strand"
        else:
            self.__strand = "Unknown Strand"


def gbff_processor(filename, object_dictionary):
    """Opens the file, within scope of with open, processes the file
    to be converted into usable data through further functions.
    :param filename
    :param object_dictionary
    :return None"""
    with open(filename, encoding="utf8") as file:
        feature_checkpoint = False
        cds_buffer = []
        last_location = None
        while True:
            line = file.readline()
            if not line:
                break
            if line.startswith("FEATURES"):
                # checks to see if "FEATURES" has been passed as this
                # can be ignored.
                feature_checkpoint = True
            if feature_checkpoint:
                location_ = location(line)
                if location_ == "CDS":
                    # split on spaces and create list, get last element
                    # from list as that is of interest to set CDS_range
                    split_list = line.split(" ")
                    cds_range = split_list[-1]
                if location_ != "segment":
                    last_location = location_
                    if cds_buffer:
                        process_cds(cds_buffer, cds_range, object_dictionary)
                        cds_buffer = []
                elif last_location == "CDS":
                    cds_buffer.append(line)


def location(line):
    """Determines the position of the current line inside gbff file.
    :param line
    :returns source, gene, CDS, segment"""
    # Line strip can NOT be used in this scenario as "gene" is used in
    # associated segments and causes errors by prematurely cutting off.
    # This functions as an exception handler.
    if line.startswith("     source"):
        return "source"
    elif line.startswith("     gene"):
        return "gene"
    elif line.startswith("     CDS"):
        return "CDS"
    else:
        return "segment"


def process_cds(cds_buffer, cds_range, object_dictionary):
    """Uses list buffer and range and filters out unwanted elements,
    processes desired data into dictionaries.
    :param cds_buffer
    :param cds_range
    :param object_dictionary"""
    processed_cds_map = {}
    # set default values to None for use in sanity check.
    current_key, current_value = None, None
    for element in cds_buffer:
        element = element.strip()
        element = element.replace('"', "")
        if element.startswith("/"):
            current_key, current_value = element.replace("/", "")\
                .split("=")
        # sanity check / exception handle to verify if value is set.
        elif current_value:
            # if current value does not start with "/", implies that
            # next line is part of previous value. therefore adds it
            # to 'last' current value + 'current' current value.
            current_value += element
        # sanity check / exception handle to verify if value is set.
        if current_key:
            # if current key has a value, overrides, meaning that after
            # writing last value, will override with new value in case
            # of multi-line value.
            processed_cds_map[current_key] = current_value
    processed_cds_map["matched_pattern"] = []
    if his_pattern.search(processed_cds_map["translation"]):
        processed_cds_map["matched_pattern"].append("histidine")
    if ser_pattern.search(processed_cds_map["translation"]):
        processed_cds_map["matched_pattern"].append("serine")
    # If a RegEx pattern could NOT be found, does not create object out
    # of the current data.
    if not processed_cds_map["matched_pattern"]:
        return
    processed_cds_map["cds_range"] = cds_range
    object_generator(processed_cds_map, object_dictionary)


def object_generator(processed_cds_map, object_dictionary):
    """Generates and object and updates dictionary from main for later
    reference.
    :param processed_cds_map
    :param object_dictionary"""
    object_dictionary[processed_cds_map["db_xref"]]\
        = GeneBankData(processed_cds_map)


def gff_processor(filename, object_dictionary):
    """Reads out the .GFF file and creates a buffer of the data inside.
    :param filename
    :param object_dictionary"""
    with open(filename, encoding="utf8") as file:
        species_checkpoint = False
        while True:
            line = file.readline()
            # breaks if at end of file, indicated in file by triple #.
            if not line or line.startswith("###"):
                break
            if line.startswith("##species"):
                species_checkpoint = True
            elif species_checkpoint:
                gff_tab_buffer = line.split("\t")
                gff_data_organizer(gff_tab_buffer, object_dictionary)


def gff_data_organizer(gff_tab_buffer, object_dictionary):
    """Uses the database cross reference or 'Dxref' to link the
    two files, get the associated data desired and write it to the
    correct object."""
    # Passes the first line as it is not useful for data gathering.
    if gff_tab_buffer[2] == "region":
        return
    attribute_string = gff_tab_buffer[8]
    attribute_list = attribute_string.split(";")
    gff_attribute_map = {}
    for element in attribute_list:
        tag, value = element.split("=")
        gff_attribute_map[tag] = value
    # sanity check to catch possible KeyError exception.
    if gff_attribute_map["Dbxref"]:
        dbxref_value = gff_attribute_map["Dbxref"]
        dbxref_list = dbxref_value.split(",")
        # Use list comprehension in order to only pass on values
        # that contain "Gene ID:"
        try:
            id_ = [string for string in dbxref_list if "GeneID:" in string][0]
            if id_ in object_dictionary:
                object_dictionary[id_].set_strand(gff_tab_buffer[6])
        # In case of dbxref_list not having element containing "GeneID:"
        # print error message and exit program.
        except IndexError:
            print("No GeneID reference could be found.\n"
                  "Is .gff file correctly formatted?")
            exit()


def display_data(object_dictionary):
    """"Shows the desired data from the first 5 objectives on the test
    in a GUI through the generated object dictionary.
    :param object_dictionary"""
    root = tk.Tk()
    root.title("Protein GUI")
    window_width = round(root.winfo_screenwidth()/1.9)
    window_height = round(root.winfo_screenheight()/1.9)
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    center_x = int(screen_width / 2 - window_width / 2)
    center_y = int(screen_height / 2 - window_height / 2)
    root.geometry(f"{window_width}x{window_height}+{center_x}+{center_y}")
    root.iconbitmap("favicon.ico")
    columns = ("Product", "Protein_ID", "Active Site", "Strand", "Sequence")
    tree = ttk.Treeview(root, columns=columns, show='headings')
    tree.heading("Product", text="Product")
    tree.heading("Protein_ID", text="Protein_ID")
    tree.heading("Active Site", text="Active Site")
    tree.heading("Strand", text="Strand")
    tree.heading("Sequence", text="Sequence")
    data_list = []
    for i in object_dictionary:
        data = ([])
        data.append(object_dictionary[i].get_protein_product())
        data.append(object_dictionary[i].get_protein_id())
        data.append(object_dictionary[i].get_pattern())
        data.append(object_dictionary[i].get_strand())
        data.append(object_dictionary[i].get_protein_sequence())
        data_list.append(data)
    for data in data_list:
        tree.insert('', tk.END, values=data)

    # nested definition needed to pass tree.
    def item_selected(event):
        """"Allows the user to click the 'Treeview' to open up a display
        that shows the user the sequence of the protein.
        Note: Event is needed for functionality even though it is seen
        as 'not used'.
        :param event"""
        print(event)
        for selected_item in tree.selection():
            item = tree.item(selected_item)
            record = item['values'][4]
            # show a message
            showinfo(title='Sequence', message=record)
    tree.bind('<<TreeviewSelect>>', item_selected)
    tree.grid(row=0, column=0, sticky='nsew')
    scrollbar = ttk.Scrollbar(root, orient=tk.VERTICAL, command=tree.yview)
    tree.configure(yscroll=scrollbar.set)
    scrollbar.grid(row=0, column=1, sticky='ns')
    root.mainloop()


def display_strands(object_dictionary):
    """""This would display a nice looking graph if matplotlib
    wasn't completely fucking atrocious to use.
    :param object_dictionary"""
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    plt.title("Strands")
    plt.xlabel("Different strands")
    plt.ylabel("Values")
    value_list = []
    for key in object_dictionary:
        value = object_dictionary[key].get_strand()
        value_list.append(value)
    strands = ["Forward", "Reverse"]
    amount = len(value_list)
    ax.bar(strands, amount)
    plt.show()


def main():
    database_file = "GCF_000013425.1_ASM1342v1_genomic.gbff"
    gene_file = "GCF_000013425.1_ASM1342v1_genomic.gff"
    # Make a dictionaries/map to store genes that match the RegEx
    # of the consensus pattern.
    object_dictionary = {}
    gbff_processor(database_file, object_dictionary)
    gff_processor(gene_file, object_dictionary)
    display_strands(object_dictionary)
    display_data(object_dictionary)


if __name__ == '__main__':
    main()
