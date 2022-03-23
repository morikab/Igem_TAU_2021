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
        self.hosts = []

        # self.feet = StringVar()
        # feet_entry = ttk.Entry(mainframe, width=7, textvariable=self.feet)
        # feet_entry.grid(column=2, row=1, sticky=(W, E))
        #
        # self.meters = StringVar()
        # ttk.Label(mainframe, textvariable=self.meters).grid(column=2, row=2, sticky=(W, E))
        #
        # ttk.Button(mainframe, text="Calculate", command=self.calculate).grid(column=3, row=3, sticky=W)
        #
        # ttk.Label(mainframe, text="feet").grid(column=3, row=1, sticky=W)
        # ttk.Label(mainframe, text="is equivalent to").grid(column=1, row=2, sticky=E)
        # ttk.Label(mainframe, text="meters").grid(column=3, row=2, sticky=W)
        #
        # for child in mainframe.winfo_children():
        #     child.grid_configure(padx=5, pady=5)
        #
        # feet_entry.focus()
        # root.bind("<Return>", self.calculate)
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

        for index, genome_path in enumerate(hosts_files):
            row = index + 1

            host_name = StringVar()
            host_name.set(os.path.basename(genome_path))
            ttk.Entry(self.wanted_hosts_grid, textvariable=host_name).grid(column=0, row=row)
            # host_name_entry.insert(0, os.path.basename(genome_path))

            # TODO - make it collapsible/fixed size
            ttk.Label(self.wanted_hosts_grid, text=genome_path).grid(column=1, row=row)
            optimization_priority = IntVar()
            optimization_priority.set(50)
            ttk.Spinbox(self.wanted_hosts_grid, from_=1, to=100, textvariable=optimization_priority).grid(column=2,
                                                                                                          row=row)
            ttk.Button(self.wanted_hosts_grid, text="remove", command=self.remove_host).grid(column=4, row=row)

            organism = {
                "host_name": host_name,
                "optimized": True,
                "expression_csv": None,  # TODO - fix this
                "optimization_priority": optimization_priority,
            }
            organisms[genome_path] = organism

        for child in self.wanted_hosts_grid.winfo_children():
            child.grid_configure(padx=5, pady=5)

        self.wanted_hosts_grid.pack()

    def remove_host(self):
        # TODO - add logic
        pass

    @staticmethod
    def clear_frame(frame: Frame) -> None:
        # destroy all widgets from frame
        for widget in frame.winfo_children():
            widget.destroy()


if __name__ == '__main__':
    # root = ThemedTk(theme="arc")      # TODO - very very slow. Try to find another reasonable style
    root = Tk()
    CommuniqueApp(root)
    root.mainloop()
