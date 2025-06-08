import customtkinter as ctk
import loguru

from pathlib import Path
from frames.inputframe import InputFrame
from frames.outputframe import OutputFrame
from frames.tableframe import TableFrame


class Interface(ctk.CTk):
    def __init__(self, fg_color="#0d1117", **kwargs):
        super().__init__(fg_color, **kwargs)
        self.title("QuantiMol")
        self.geometry("1000x600")
        self.resizable(False, False)

        try:
            icon_path = Path(__file__).parent / "images" / "icon.ico"
            self.iconbitmap(icon_path)
        except Exception:
            loguru.logger.warning("Icon file not found, using default icon.")

        iframe = InputFrame(self)
        tframe = TableFrame(self)
        oframe = OutputFrame(self)
        iframe.pack(side="left", fill="both", expand=False, padx=15, pady=15)
        tframe.pack(side="top", fill="both", expand=False, padx=15, pady=15)
        oframe.pack(side="top", fill="both", expand=True, padx=15, pady=15)


if __name__ == "__main__":
    loguru.logger.info("Starting QuantiMol Interface")
    root = Interface()
    root.mainloop()
    loguru.logger.info("QuantiMol Interface closed")
