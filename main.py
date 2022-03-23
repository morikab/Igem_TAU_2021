from tkinter import *
from tkinter import filedialog
from tkinter import ttk


class CommuniqueApp(object):
    def __init__(self, root: Tk) -> None:
        root.title("Communique")

        mainframe = ttk.Frame(root, padding="3 3 12 12")
        mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        root.geometry("500x1200")  # Set window size

        # Sequence to optimize
        self.sequence_label_frame = ttk.Labelframe(mainframe, text="Sequence to Optimize")
        self.sequence_label_frame.pack(fill="both", expand="yes")    # TODO - explore options
        upload_sequence = ttk.Button(self.sequence_label_frame, text="Upload Sequence", command=self.upload_sequence)
        upload_sequence.pack(side=LEFT)
        # TODO - consider another widget with white background for the path
        self.sequence_path_label = ttk.Label(self.sequence_label_frame)
        self.sequence_path_label.pack(side=LEFT)

        # Wanted hosts
        self.wanted_hosts_frame = ttk.Labelframe(mainframe, text="Wanted Hosts")
        self.wanted_hosts_frame.pack(fill="both", expand="yes")  # TODO - explore options

        # Unwanted hosts
        self.unwanted_hosts_frame = ttk.Labelframe(mainframe, text="Unwanted Hosts")
        self.unwanted_hosts_frame.pack(fill="both", expand="yes")  # TODO - explore options

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

    def upload_sequence(self):
        sequence_file = filedialog.askopenfile()
        self.sequence_path_label.config(text=sequence_file.name)
        # TODO - display file name

    # def calculate(self, *args):
    #     try:
    #         value = float(self.feet.get())
    #         self.meters.set(int(0.3048 * value * 10000.0 + 0.5)/10000.0)
    #     except ValueError:
    #         pass


if __name__ == '__main__':
    root = Tk()
    CommuniqueApp(root)
    root.mainloop()
