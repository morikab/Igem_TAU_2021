import os
import typing
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
        master.geometry("700x1200")  # Set window size

        # Scroll Bar
        self.scrollbar = ttk.Scrollbar(mainframe, orient="vertical")
        self.scrollbar.pack(side="right", fill="y")

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

        # Optimize Button
        self.optimize_button = ttk.Button(mainframe, text="Optimize", command=self.optimize)
        self.optimize_button.pack(side=BOTTOM, pady=20)

        # User Input Parameters
        self.organisms = {}
        self.sequence = None

    @property
    def wanted_hosts_count(self) -> int:
        return len({key: value for key, value in self.organisms.items() if value["optimized"]})

    def upload_sequence(self):
        # TODO - add *.fasta files
        sequence_file_name = filedialog.askopenfilename(filetypes=[("Fasta files", "*.fa",)])
        self.sequence_path_label.config(text=sequence_file_name)
        # TODO - add Label or entry for the file name
        self.sequence = sequence_file_name

    def upload_hosts_files(self):
        # TODO - throw error if using the same file twice (same/different groups)
        hosts_files = filedialog.askopenfilename(filetypes=[("Genebank files", "*.gb")], multiple=True)
        if len(hosts_files) > 10:   # TODO - is really necessary?
            # TODO - convert to alert box
            print("Cannot upload more than 10 files. Please remove some files and try again")

        ttk.Label(self.wanted_hosts_grid, text="host name").grid(column=0, row=0)
        ttk.Label(self.wanted_hosts_grid, text="genome path").grid(column=1, row=0)
        ttk.Label(self.wanted_hosts_grid, text="optimization priority").grid(column=2, row=0)
        ttk.Label(self.wanted_hosts_grid, text="expression levels path (optional)").grid(column=3, row=0)

        initial_row = 1
        if self.wanted_hosts_count > 0:
            initial_row = self.wanted_hosts_count + 1

        for index, genome_path in enumerate(hosts_files):
            row = initial_row + index

            host_name_var = StringVar()
            host_name_var.set(os.path.basename(genome_path))
            ttk.Entry(self.wanted_hosts_grid, textvariable=host_name_var).grid(column=0, row=row)

            genome_path_var = StringVar()
            genome_path_var.set(genome_path)
            # TODO - make the path scrollable
            ttk.Entry(self.wanted_hosts_grid, textvariable=genome_path_var).grid(column=1, row=row)
            optimization_priority_var = IntVar()
            optimization_priority_var.set(50)
            # TODO - recheck from and to parameters here.
            ttk.Spinbox(self.wanted_hosts_grid, from_=1, to=100, textvariable=optimization_priority_var).grid(column=2,
                                                                                                              row=row)

            self.create_upload_expression_frame(row=row)

            remove_button = ttk.Button(self.wanted_hosts_grid, text="remove")
            remove_button.grid(column=4, row=row)
            remove_button.bind("<Button-1>", self.remove_host)

            organism = {
                "host_name": host_name_var,
                "genome_path": genome_path_var,
                "optimized": True,
                "expression_csv": None,
                "optimization_priority": optimization_priority_var,
            }
            self.organisms[genome_path] = organism

        for child in self.wanted_hosts_grid.winfo_children():
            child.grid_configure(padx=5, pady=5)
        self.wanted_hosts_grid.pack()

    def upload_expression(self, event) -> None:
        expression_file_name = filedialog.askopenfilename(filetypes=[("csv", "*.csv")])
        upload_widget = event.widget
        host_row = upload_widget.master.grid_info()["row"]

        upload_widget.master.destroy()
        self.create_upload_expression_frame(row=host_row, expression_file_name=expression_file_name)

    def create_upload_expression_frame(self, row: int, expression_file_name: typing.Optional[str] = None):
        expression_frame = Frame(self.wanted_hosts_grid)
        expression_frame.grid(row=row, column=3)

        expression_level_button = ttk.Button(expression_frame, text="upload")
        expression_level_button.bind("<Button-1>", self.upload_expression)
        expression_level_button.pack()

        if expression_file_name is not None:
            expression_path_var = StringVar()
            expression_path_var.set(expression_file_name)
            expression_path_entry = ttk.Entry(expression_frame, textvariable=expression_path_var)
            # expression_label = ttk.Label(expression_frame, text=expression_file_name)
            expression_path_entry.pack()
            host_info = self.get_host_info(row=row)
            host_info["expression_csv"] = expression_path_var

    def get_host_info(self, row: int) -> typing.Dict:
        # TODO - define consts for column numbers
        genome_path = self.wanted_hosts_grid.grid_slaves(row=row, column=1)[-1].get()
        return self.organisms[genome_path]

    def remove_host(self, event) -> None:
        number_of_columns = 5
        remove_widget = event.widget
        row_to_remove = remove_widget.grid_info()["row"]
        genome_path = self.wanted_hosts_grid.grid_slaves(row=row_to_remove, column=1)[-1].get()
        remove_widget.grid_remove()
        for j in range(number_of_columns - 1):
            # TODO - check using destroy instead
            self.wanted_hosts_grid.grid_slaves(row=row_to_remove, column=j)[0].grid_remove()

        for i in range(row_to_remove+1, self.wanted_hosts_count+1):
            for j in range(number_of_columns):
                self.wanted_hosts_grid.grid_slaves(row=i, column=j)[0].grid(row=i-1, column=j)

        self.organisms.pop(genome_path)
        self.wanted_hosts_grid.pack()

        if not self.organisms:
            self.wanted_hosts_grid.destroy()
            self.wanted_hosts_grid = Frame(self.wanted_hosts_frame)

    def optimize(self) -> None:
        # TODO - create dict, run modules and display the final window
        print("Optimize")


if __name__ == "__main__":
    # root = ThemedTk(theme="arc")      # TODO - very very slow. Try to find another reasonable style
    root = Tk()
    CommuniqueApp(root)
    root.mainloop()
