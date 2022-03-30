import typing
from pathlib import Path
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import ttk

from modules.main import run_modules


class CommuniqueApp(object):
    MAX_HOSTS_COUNT = 10
    INITIAL_WIDTH = 700
    INITIAL_HEIGHT = 1000
    DEFAULT_TUNING_PARAMETER_VALUE = 50
    DEFAULT_CLUSTERS_COUNT_VALUE = 2
    DEFAULT_PRIORITY_VALUE = 50
    HOST_NAME_COLUMN_INDEX = 0
    GENOME_PATH_COLUMN_INDEX = 1
    HOST_PRIORITY_COLUMN_INDEX = 2
    EXPRESSION_LEVEL_COLUMN_INDEX = 3
    REMOVE_HOST_COLUMN_INDEX = 4

    def __init__(self, master: tk.Tk) -> None:
        self.organisms = {}
        self.sequence = None
        self.tuning_parameter = None
        self.clusters_count = None

        # Widgets
        self.master = master
        self.mainframe = None
        self.sequence_path_label = None
        self.wanted_hosts_frame = None
        self.wanted_hosts_grid = None
        self.unwanted_hosts_frame = None
        self.unwanted_hosts_grid = None
        self.bottom_frame = None
        self.results_frame = None

        self.initialize_master_widget()
        self.prepare_for_new_run()

    def initialize_master_widget(self) -> None:
        self.master.title("Communique")
        self.master.geometry(F"{self.INITIAL_WIDTH}x{self.INITIAL_HEIGHT}")

    def prepare_for_new_run(self) -> None:
        canvas = tk.Canvas(self.master, width=self.INITIAL_WIDTH, height=self.INITIAL_HEIGHT)
        scrollbar = ttk.Scrollbar(self.master, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        canvas.pack(fill="both", expand=True)

        self.mainframe = ttk.Frame(canvas)
        window = canvas.create_window((0, 0), window=self.mainframe, anchor="nw")
        self.mainframe.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(window, width=e.width))
        canvas.bind_all("<MouseWheel>", lambda e: canvas.yview_scroll(-1 * e.delta//120, "units"))
        self.master.after_idle(canvas.yview_moveto, 0)

        #########################
        # Sequence to optimize
        #########################
        sequence_label_frame = ttk.Labelframe(self.mainframe, text="Sequence to Optimize")
        sequence_label_frame.pack(fill="both", expand="yes")
        upload_sequence = ttk.Button(sequence_label_frame, text="Upload Sequence", command=self.upload_sequence)
        upload_sequence.pack(side=tk.TOP)
        self.sequence_path_label = ttk.Label(sequence_label_frame)
        self.sequence_path_label.pack(side=tk.TOP, pady=5)

        #########################
        # Wanted hosts
        #########################
        self.wanted_hosts_frame = ttk.Labelframe(self.mainframe, text="Wanted Hosts")
        self.wanted_hosts_frame.pack(fill="both", expand="yes")
        upload_wanted_hosts = ttk.Button(self.wanted_hosts_frame,
                                         text="Upload .gb files",
                                         command=self.upload_wanted_hosts_files)
        upload_wanted_hosts.pack(side=tk.TOP, pady=5)
        self.wanted_hosts_grid = tk.Frame(self.wanted_hosts_frame)

        #########################
        # Unwanted hosts
        #########################
        self.unwanted_hosts_frame = ttk.Labelframe(self.mainframe, text="Unwanted Hosts")
        self.unwanted_hosts_frame.pack(fill="both", expand="yes")
        upload_unwanted_hosts = ttk.Button(self.unwanted_hosts_frame,
                                           text="Upload .gb files",
                                           command=self.upload_unwanted_hosts_files)
        upload_unwanted_hosts.pack(side=tk.TOP, pady=5)
        self.unwanted_hosts_grid = tk.Frame(self.unwanted_hosts_frame)

        #########################
        # Bottom Frame
        #########################
        self.bottom_frame = ttk.Frame(self.mainframe)
        self.bottom_frame.pack(side=tk.BOTTOM, pady=20)
        # Advanced Options
        options_button = ttk.Button(self.bottom_frame, text="Advanced Options", command=self.advanced_options)
        options_button.grid(row=0, column=0)
        # Optimize Button
        optimize_button = ttk.Button(self.bottom_frame, text="Optimize", command=self.optimize)
        optimize_button.grid(row=0, column=1)

        # User Input Parameters
        self.initialize_user_input_parameters()

    def initialize_user_input_parameters(self) -> None:
        self.organisms = {}
        self.sequence = None
        self.tuning_parameter = tk.IntVar()
        self.tuning_parameter.set(self.DEFAULT_TUNING_PARAMETER_VALUE)
        self.clusters_count = tk.IntVar()
        self.clusters_count.set(self.DEFAULT_CLUSTERS_COUNT_VALUE)

    def get_hosts_count(self, is_optimized: bool) -> int:
        return len({key: value for key, value in self.organisms.items() if value["optimized"] == is_optimized})

    def upload_sequence(self) -> None:
        sequence_file_name = filedialog.askopenfilename(filetypes=[("Fasta files", ".fa .fasta")])
        self.sequence_path_label.config(text=sequence_file_name)
        self.sequence = sequence_file_name

    def upload_wanted_hosts_files(self) -> None:
        self.upload_hosts_files(grid=self.wanted_hosts_grid, is_optimized=True)

    def upload_unwanted_hosts_files(self) -> None:
        self.upload_hosts_files(grid=self.unwanted_hosts_grid, is_optimized=False)

    def upload_hosts_files(self, grid: tk.Frame, is_optimized: bool) -> None:
        hosts_files = filedialog.askopenfilename(filetypes=[("Genebank files", ".gb")], multiple=True)
        if not self.validate_uploaded_files(hosts_files):
            return

        headers_row = 0
        ttk.Label(grid, text="host name").grid(column=self.HOST_NAME_COLUMN_INDEX, row=headers_row)
        ttk.Label(grid, text="genome path").grid(column=self.GENOME_PATH_COLUMN_INDEX, row=headers_row)
        ttk.Label(grid, text="optimization priority").grid(column=self.HOST_PRIORITY_COLUMN_INDEX, row=headers_row)
        ttk.Label(grid, text="expression levels (optional)").grid(column=self.EXPRESSION_LEVEL_COLUMN_INDEX,
                                                                  row=headers_row)

        initial_row = self.get_hosts_count(is_optimized) + 1

        for index, genome_path in enumerate(hosts_files):
            row = initial_row + index

            host_name_var = tk.StringVar()
            host_name = Path(genome_path).stem
            host_name_var.set(host_name)
            ttk.Entry(grid, textvariable=host_name_var).grid(column=self.HOST_NAME_COLUMN_INDEX, row=row)

            genome_path_var = tk.StringVar()
            genome_path_var.set(genome_path)
            ttk.Entry(grid, textvariable=genome_path_var).grid(column=self.GENOME_PATH_COLUMN_INDEX, row=row)
            optimization_priority_var = tk.IntVar()
            optimization_priority_var.set(self.DEFAULT_PRIORITY_VALUE)
            ttk.Spinbox(grid, from_=1, to=100, textvariable=optimization_priority_var).grid(
                column=self.HOST_PRIORITY_COLUMN_INDEX,
                row=row,
            )

            self.create_upload_expression_frame(grid=grid, row=row)

            remove_button = ttk.Button(grid, text="remove")
            remove_button.grid(column=self.REMOVE_HOST_COLUMN_INDEX, row=row)
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

    def validate_uploaded_files(self, hosts_files: typing.Sequence[str]) -> bool:
        if len(hosts_files) > self.MAX_HOSTS_COUNT:
            messagebox.showerror(title="Error",
                                 message=F"Please upload at most {self.MAX_HOSTS_COUNT} files.")
            return False

        for host in hosts_files:
            if host in self.organisms.keys():
                messagebox.showerror(title="Error",
                                     message=F"Host file {host} already uploaded.")
                return False

        return True

    def upload_host_expression_file(self, event) -> None:
        expression_file_name = filedialog.askopenfilename(filetypes=[("csv", " .csv")])
        upload_widget = event.widget
        grid = upload_widget.master.master
        host_row = upload_widget.master.grid_info()["row"]

        upload_widget.master.destroy()
        self.create_upload_expression_frame(grid=grid, row=host_row, expression_file_name=expression_file_name)

    def create_upload_expression_frame(
            self,
            grid: tk.Frame,
            row: int,
            expression_file_name: typing.Optional[str] = None,
    ) -> None:
        expression_frame = tk.Frame(grid)
        expression_frame.grid(row=row, column=self.EXPRESSION_LEVEL_COLUMN_INDEX)

        expression_level_button = ttk.Button(expression_frame, text="upload")
        expression_level_button.bind("<Button-1>", self.upload_host_expression_file)
        expression_level_button.pack()

        if expression_file_name is not None:
            expression_path_var = tk.StringVar()
            expression_path_var.set(expression_file_name)
            expression_path_entry = ttk.Entry(expression_frame, textvariable=expression_path_var)
            expression_path_entry.pack()

            host_info = self.get_host_info(grid=grid, row=row)
            host_info["expression_csv"] = expression_path_var

    def get_host_info(self, grid: tk.Frame, row: int) -> typing.Dict:
        genome_path = grid.grid_slaves(row=row, column=self.GENOME_PATH_COLUMN_INDEX)[-1].get()
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
        genome_path = grid.grid_slaves(row=row_to_remove, column=self.GENOME_PATH_COLUMN_INDEX)[-1].get()
        remove_widget.grid_remove()
        for j in range(number_of_columns - 1):
            # TODO - check using destroy instead
            grid.grid_slaves(row=row_to_remove, column=j)[0].destroy()
            # grid.grid_slaves(row=row_to_remove, column=j)[0].grid_remove()

        for i in range(row_to_remove+1, self.get_hosts_count(is_optimized)+1):
            for j in range(number_of_columns):
                grid.grid_slaves(row=i, column=j)[0].grid(row=i-1, column=j)

        self.organisms.pop(genome_path)
        grid.pack()

        if self.get_hosts_count(is_optimized) == 0:
            grid.destroy()
            if is_optimized:
                self.wanted_hosts_grid = tk.Frame(self.wanted_hosts_frame)
            else:
                self.unwanted_hosts_grid = tk.Frame(self.unwanted_hosts_frame)

    def advanced_options(self) -> None:
        options_window = tk.Toplevel(self.mainframe)
        options_window.title("Advanced Options")
        options_window.geometry("300x200")

        options_frame = ttk.Frame(options_window, padding="5 5 12 12")
        ttk.Label(options_frame, text="Tuning Parameter: ").grid(row=0, column=0)
        ttk.Spinbox(options_frame, from_=1, to=100, textvariable=self.tuning_parameter).grid(row=0, column=1)

        ttk.Label(options_frame, text="Clusters Count: ").grid(row=1, column=0)
        # TODO - what is the max clusters count?
        ttk.Spinbox(options_frame, from_=2, to=10, textvariable=self.clusters_count).grid(row=1, column=1)
        # TODO - add option to configure the optimization method in advanced options

        for child in options_frame.winfo_children():
            child.grid_configure(padx=5, pady=5)
        options_frame.pack()

    def optimize(self) -> None:
        if not self.validate_user_input():
            return

        user_input = self.generate_user_input()

        # TODO - alert on errors (if there are any)
        user_output, zip_file_path = run_modules(user_input_dict=user_input)

        #########################
        # Results widgets
        #########################
        self.bottom_frame.destroy()

        self.results_frame = ttk.Labelframe(self.mainframe, text="Results")
        self.results_frame.pack(fill="both", expand="yes")
        results_grid = ttk.Frame(self.results_frame)
        results_grid.pack(side=tk.TOP, pady=5)

        self.bottom_frame = ttk.Frame(self.mainframe)
        self.bottom_frame.pack(side=tk.BOTTOM, pady=20)

        # New Run
        ttk.Label(results_grid, text="Optimized sequence:").grid(row=0, column=0)
        optimized_sequence = tk.Text(results_grid)
        optimized_sequence.grid(row=0, column=1)
        optimized_sequence.insert(tk.END, user_output["final_sequence"])

        new_run_button = ttk.Button(self.bottom_frame, text="New run", command=self.recreate)
        new_run_button.grid(row=0, column=0)

    def generate_user_input(self) -> typing.Dict:
        user_input = {
            "sequence": self.sequence,
            "tuning_param": self.tuning_parameter.get() / 100,
            "clusters_count": self.clusters_count.get(),
        }
        input_organisms = {}
        for organism in self.organisms.values():
            organism_name = organism["host_name"].get()
            input_organisms[organism_name] = {
                "genome_path": organism["genome_path"].get(),
                "expression_csv": organism["expression_csv"].get() if organism["expression_csv"] else None,
                "optimized": organism["optimized"],
                "optimization_priority": organism["optimization_priority"].get(),
            }
        user_input["organisms"] = input_organisms
        return user_input

    def validate_user_input(self) -> bool:
        if self.sequence is None:
            messagebox.showerror(title="Error",
                                 message=F"Please choose a sequence file for optimization.")
            return False
        if self.get_hosts_count(is_optimized=True) == 0:
            messagebox.showerror(title="Error",
                                 message=F"Please choose at least one wanted host file.")
            return False
        if self.get_hosts_count(is_optimized=False) == 0:
            messagebox.showerror(title="Error",
                                 message=F"Please choose at least one unwanted host file.")
            return False

        return True

    def recreate(self) -> None:
        # TODO - should we recreate all the params, just run again the modules (add two buttons for both?)
        for widget in self.master.winfo_children():
            widget.destroy()
        self.prepare_for_new_run()


if __name__ == "__main__":
    root = tk.Tk()
    CommuniqueApp(root)
    root.mainloop()
