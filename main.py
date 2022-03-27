import os
from tkinter import *
from tkinter import filedialog
from tkinter import ttk
from ttkthemes import ThemedTk


class CommuniqueApp(object):
    def __init__(self, master: Tk) -> None:
        master.title("Communique")

        mainframe = ttk.Frame(master, padding="3 3 12 12")
        mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
        master.columnconfigure(0, weight=1)
        master.rowconfigure(0, weight=1)
        master.geometry("500x1200")  # Set window size

        # Sequence to optimize
        self.sequence_label_frame = ttk.Labelframe(mainframe, text="Sequence to Optimize")
        self.sequence_label_frame.pack(fill="both", expand="yes")
        upload_sequence = ttk.Button(self.sequence_label_frame, text="Upload Sequence", command=self.upload_sequence)
        upload_sequence.pack(side=TOP)
        # TODO - consider another widget with white background for the path
        self.sequence_path_label = ttk.Label(self.sequence_label_frame)
        self.sequence_path_label.pack(side=TOP)

        # Wanted hosts
        self.wanted_hosts_frame = ttk.Labelframe(mainframe, text="Wanted Hosts")
        self.wanted_hosts_frame.pack(fill="both", expand="yes")
        upload_wanted_hosts = ttk.Button(self.wanted_hosts_frame,
                                         text="Upload .gb files",
                                         command=self.upload_hosts_files)
        upload_wanted_hosts.pack(side=TOP)

        self.wanted_hosts_grid = Frame(self.wanted_hosts_frame)
        self.wanted_hosts_count = 0
        self.input_dict = {"organisms": {}}

        # TODO - after running "optimize", we need to regenerate the dictionary for the modules code

    def upload_sequence(self):
        sequence_file = filedialog.askopenfilename(filetypes=[("Fasta files", "*.fa",)])  # TODO - add *.fasta files
        self.sequence_path_label.config(text=sequence_file)
        # TODO - add Label or entry for the file name
        self.input_dict["sequence"] = sequence_file

    def upload_hosts_files(self):
        # TODO - need to add fltering and correct handling of organisms for various edge cases (together with unwanted hosts)
        organisms = self.input_dict["organisms"]
        hosts_files = filedialog.askopenfilename(filetypes=[("Genebank files", "*.gb")], multiple=True)
        if len(hosts_files) > 10:
            # TODO - convert to alert box
            print("Cannot upload more than 10 files. Please remove some files and try again")

        ttk.Label(self.wanted_hosts_grid, text="host name").grid(column=0, row=0)
        ttk.Label(self.wanted_hosts_grid, text="genome path").grid(column=1, row=0)
        ttk.Label(self.wanted_hosts_grid, text="optimization priority").grid(column=2, row=0)
        ttk.Label(self.wanted_hosts_grid, text="expression levels path (optional)").grid(column=3, row=0)

        initial_row = 1
        if self.wanted_hosts_count > 0:
            # TODO - handle duplicates
            initial_row = self.wanted_hosts_count + 1

        for index, genome_path in enumerate(hosts_files):
            row = initial_row + index

            host_name_var = StringVar()
            host_name_var.set(os.path.basename(genome_path))
            ttk.Entry(self.wanted_hosts_grid, textvariable=host_name_var).grid(column=0, row=row)

            genome_path_var = StringVar()
            genome_path_var.set(genome_path)
            # TODO - make the path scrollable + add option to change the path (and update path and name, accordingly)
            ttk.Entry(self.wanted_hosts_grid, textvariable=genome_path_var).grid(column=1, row=row)
            # ttk.Label(self.wanted_hosts_grid, text=genome_path).grid(column=1, row=row)
            optimization_priority = IntVar()
            optimization_priority.set(50)
            ttk.Spinbox(self.wanted_hosts_grid, from_=1, to=100, textvariable=optimization_priority).grid(column=2,
                                                                                                          row=row)

            expression_level_button = ttk.Button(self.wanted_hosts_grid, text="upload")
            expression_level_button.grid(column=3, row=row)
            expression_level_button.bind("<Button-1>", self.upload_expression)

            remove_button = ttk.Button(self.wanted_hosts_grid, text="remove")
            remove_button.grid(column=4, row=row)
            remove_button.bind("<Button-1>", self.remove_host)

            organism = {
                "host_name": host_name_var,
                "genome_path": genome_path_var,
                "optimized": True,
                "expression_csv": None,
                "optimization_priority": optimization_priority,
            }
            organisms[genome_path] = organism

        for child in self.wanted_hosts_grid.winfo_children():
            child.grid_configure(padx=5, pady=5)
        self.wanted_hosts_grid.pack()

        self.wanted_hosts_count += len(hosts_files)

    def upload_expression(self, event) -> None:
        expression_file_name = filedialog.askopenfilename(filetypes=[("csv", "*.csv")])
        upload_widget = event.widget
        host_row = upload_widget.grid_info()["row"]
        # TODO - define consts for column numbers
        genome_path = self.wanted_hosts_grid.grid_slaves(row=host_row, column=1)[0]
        host_entry = self.input_dict["organisms"][genome_path]
        host_entry["expression_csv"] = expression_file_name

    def remove_host(self, event) -> None:
        number_of_columns = 5
        remove_widget = event.widget
        row_to_remove = remove_widget.grid_info()["row"]
        remove_widget.grid_remove()
        for i in range(row_to_remove, self.wanted_hosts_count+1):
            if i == self.wanted_hosts_count:
                for j in range(number_of_columns-1):
                    self.wanted_hosts_grid.grid_slaves(row=i, column=j)[0].grid_remove()
            else:
                for j in range(number_of_columns):
                    self.wanted_hosts_grid.grid_slaves(row=i, column=j)[0].grid(row=i-1, column=j)
        self.wanted_hosts_count -= 1


if __name__ == '__main__':
    # root = ThemedTk(theme="arc")      # TODO - very very slow. Try to find another reasonable style
    root = Tk()
    CommuniqueApp(root)
    root.mainloop()
