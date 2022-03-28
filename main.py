import os
import typing
from tkinter import *
from tkinter import filedialog
from tkinter import ttk
from ttkthemes import ThemedTk


class CommuniqueApp(object):
    def __init__(self, master: Tk) -> None:
        master.title("Communique")

        self.mainframe = ttk.Frame(master, padding="3 3 12 12")
        self.mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
        master.columnconfigure(0, weight=1)
        master.rowconfigure(0, weight=1)
        master.geometry("700x1000")  # Set window size

        # Scroll Bar
        # TODO - fix it to make it scrollable for real
        self.scrollbar = ttk.Scrollbar(self.mainframe, orient="vertical")
        self.scrollbar.pack(side="right", fill="y")

        #########################
        # Sequence to optimize
        #########################
        self.sequence_label_frame = ttk.Labelframe(self.mainframe, text="Sequence to Optimize")
        self.sequence_label_frame.pack(fill="both", expand="yes")
        upload_sequence = ttk.Button(self.sequence_label_frame, text="Upload Sequence", command=self.upload_sequence)
        upload_sequence.pack(side=TOP)
        # TODO - consider another widget with white background for the path
        self.sequence_path_label = ttk.Label(self.sequence_label_frame)
        self.sequence_path_label.pack(side=TOP)

        #########################
        # Wanted hosts
        #########################
        self.wanted_hosts_frame = ttk.Labelframe(self.mainframe, text="Wanted Hosts")
        self.wanted_hosts_frame.pack(fill="both", expand="yes")
        upload_wanted_hosts = ttk.Button(self.wanted_hosts_frame,
                                         text="Upload .gb files",
                                         command=self.upload_wanted_hosts_files)
        upload_wanted_hosts.pack(side=TOP)
        self.wanted_hosts_grid = Frame(self.wanted_hosts_frame)

        #########################
        # Unwanted hosts
        #########################
        self.unwanted_hosts_frame = ttk.Labelframe(self.mainframe, text="Unwanted Hosts")
        self.unwanted_hosts_frame.pack(fill="both", expand="yes")
        upload_unwanted_hosts = ttk.Button(self.unwanted_hosts_frame,
                                           text="Upload .gb files",
                                           command=self.upload_unwanted_hosts_files)
        upload_unwanted_hosts.pack(side=TOP)
        self.unwanted_hosts_grid = Frame(self.unwanted_hosts_frame)

        # TODO - think again about the location of the buttons (right/left/center)
        bottom_frame = ttk.Frame(self.mainframe)
        bottom_frame.pack(side=BOTTOM, pady=20)

        # Advanced Options
        self.options_button = ttk.Button(bottom_frame, text="Advanced Options", command=self.advanced_options)
        self.options_button.grid(row=0, column=0)

        # Optimize Button
        self.optimize_button = ttk.Button(bottom_frame, text="Optimize", command=self.optimize)
        self.optimize_button.grid(row=0, column=1)

        # User Input Parameters
        self.organisms = {}
        self.sequence = None

        self.tuning_parameter = IntVar()
        self.tuning_parameter.set(50)       # TODO - move to constant

    def get_hosts_count(self, is_optimized: bool) -> int:
        return len({key: value for key, value in self.organisms.items() if value["optimized"] == is_optimized})

    def upload_sequence(self):
        # TODO - add *.fasta files
        sequence_file_name = filedialog.askopenfilename(filetypes=[("Fasta files", "*.fa",)])
        self.sequence_path_label.config(text=sequence_file_name)
        # TODO - add Label or entry for the file name
        self.sequence = sequence_file_name

    def upload_wanted_hosts_files(self) -> None:
        self.upload_hosts_files(grid=self.wanted_hosts_grid, is_optimized=True)

    def upload_unwanted_hosts_files(self) -> None:
        self.upload_hosts_files(grid=self.unwanted_hosts_grid, is_optimized=False)

    def upload_hosts_files(self, grid: Frame, is_optimized: bool) -> None:
        # TODO - throw error if using the same file twice (same/different groups)
        hosts_files = filedialog.askopenfilename(filetypes=[("Genebank files", "*.gb")], multiple=True)
        if len(hosts_files) > 10:   # TODO - is really necessary?
            # TODO - convert to alert box
            print("Cannot upload more than 10 files. Please remove some files and try again")

        ttk.Label(grid, text="host name").grid(column=0, row=0)
        ttk.Label(grid, text="genome path").grid(column=1, row=0)
        ttk.Label(grid, text="optimization priority").grid(column=2, row=0)
        ttk.Label(grid, text="expression levels (optional)").grid(column=3, row=0)

        initial_row = self.get_hosts_count(is_optimized) + 1

        for index, genome_path in enumerate(hosts_files):
            row = initial_row + index

            host_name_var = StringVar()
            host_name_var.set(os.path.basename(genome_path))
            ttk.Entry(grid, textvariable=host_name_var).grid(column=0, row=row)

            genome_path_var = StringVar()
            genome_path_var.set(genome_path)
            # TODO - make the path scrollable
            ttk.Entry(grid, textvariable=genome_path_var).grid(column=1, row=row)
            optimization_priority_var = IntVar()
            optimization_priority_var.set(50)
            ttk.Spinbox(grid, from_=1, to=100, textvariable=optimization_priority_var).grid(column=2, row=row)

            self.create_upload_expression_frame(grid=grid, row=row)

            remove_button = ttk.Button(grid, text="remove")
            remove_button.grid(column=4, row=row)
            remove_button.bind("<Button-1>", self.remove_wanted_host if is_optimized else self.remove_unwanted_host)

            organism = {
                "host_name": host_name_var,
                "genome_path": genome_path_var,
                "optimized": is_optimized,
                "expression_csv": None,
                "optimization_priority": optimization_priority_var,
            }
            self.organisms[genome_path] = organism

        for child in grid.winfo_children():
            child.grid_configure(padx=5, pady=5)
        grid.pack()

    def upload_expression(self, event) -> None:
        expression_file_name = filedialog.askopenfilename(filetypes=[("csv", "*.csv")])
        upload_widget = event.widget
        grid = upload_widget.master
        host_row = grid.grid_info()["row"]

        upload_widget.master.destroy()
        self.create_upload_expression_frame(grid=grid, row=host_row, expression_file_name=expression_file_name)

    def create_upload_expression_frame(
            self,
            grid: Frame,
            row: int,
            expression_file_name: typing.Optional[str] = None,
    ) -> None:
        expression_frame = Frame(grid)
        expression_frame.grid(row=row, column=3)

        expression_level_button = ttk.Button(expression_frame, text="upload")
        expression_level_button.bind("<Button-1>", self.upload_expression)
        expression_level_button.pack()

        if expression_file_name is not None:
            expression_path_var = StringVar()
            expression_path_var.set(expression_file_name)
            expression_path_entry = ttk.Entry(expression_frame, textvariable=expression_path_var)
            expression_path_entry.pack()

            host_info = self.get_host_info(grid=grid, row=row)
            host_info["expression_csv"] = expression_path_var

    def get_host_info(self, grid: Frame, row: int) -> typing.Dict:
        # TODO - define consts for column numbers
        genome_path = grid.grid_slaves(row=row, column=1)[-1].get()
        return self.organisms[genome_path]

    def remove_wanted_host(self, event) -> None:
        self.remove_host(event=event, is_optimized=True)

    def remove_unwanted_host(self, event) -> None:
        self.remove_host(event=event, is_optimized=False)

    def remove_host(self, event, is_optimized: bool) -> None:
        number_of_columns = 5
        remove_widget = event.widget
        grid = remove_widget.master
        row_to_remove = remove_widget.grid_info()["row"]
        genome_path = grid.grid_slaves(row=row_to_remove, column=1)[-1].get()
        remove_widget.grid_remove()
        for j in range(number_of_columns - 1):
            # TODO - check using destroy instead
            grid.grid_slaves(row=row_to_remove, column=j)[0].grid_remove()

        for i in range(row_to_remove+1, self.get_hosts_count(is_optimized)+1):
            for j in range(number_of_columns):
                grid.grid_slaves(row=i, column=j)[0].grid(row=i-1, column=j)

        self.organisms.pop(genome_path)
        grid.pack()

        if self.get_hosts_count(is_optimized) == 0:
            grid.destroy()

            if is_optimized:
                self.wanted_hosts_grid = Frame(self.wanted_hosts_frame)
            else:
                self.unwanted_hosts_grid = Frame(self.unwanted_hosts_frame)

    def advanced_options(self) -> None:
        options_window = Toplevel(self.mainframe)
        options_window.title("Advanced Options")
        options_window.geometry("300x200")

        options_frame = ttk.Frame(options_window, padding="3 3 12 12")
        ttk.Label(options_frame, text="Tuning Parameter: ").grid(row=0, column=0)
        ttk.Spinbox(options_frame, from_=1, to=100, textvariable=self.tuning_parameter).grid(row=0, column=1)
        # TODO - add clustering number (talk with Liyam)

        options_frame.pack()

    def optimize(self) -> None:
        # TODO - create dict, run modules and display the final window
        print("Optimize")


if __name__ == "__main__":
    # root = ThemedTk(theme="arc")      # TODO - very very slow. Try to find another reasonable style
    root = Tk()
    CommuniqueApp(root)
    root.mainloop()
