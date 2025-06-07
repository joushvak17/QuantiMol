import tkinter as tk
import loguru
from inputframe import InputFrame
from outputframe import OutputFrame
from tableframe import TableFrame


class Interface(tk.Tk):
    # TODO:
    # 1 - Figure out what size geometry to use and if it should be resizable
    # 2 - Figure out how to set the icon and path correctly
    def __init__(
        self,
        screenName=None,
        baseName=None,
        className="Tk",
        useTk=True,
        sync=False,
        use=None,
    ):
        super().__init__(screenName, baseName, className, useTk, sync, use)
        self.title("QuantiMol")

        self.geometry("1250x800")
        self.resizable(False, False)

        self.configure(bg="#0d1117")

        iframe = InputFrame(self)
        oframe = OutputFrame(self)
        tframe = TableFrame(self)
        iframe.pack(side="left", fill="both", expand=False, padx=15, pady=15)
        tframe.pack(side="top", fill="both", expand=True, padx=15, pady=15)
        oframe.pack(side="top", fill="both", expand=True, padx=15, pady=15)


if __name__ == "__main__":
    loguru.logger.info("Starting QuantiMol Interface")
    root = Interface()
    root.mainloop()
    loguru.logger.info("QuantiMol Interface closed")
