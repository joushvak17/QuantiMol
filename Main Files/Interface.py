import tkinter as tk
import os
import sys

from InputFrame import InputFrame
from OutputFrame import OutputFrame
from TableFrame import TableFrame

try:
    from ctypes import windll

    windll.shcore.SetProcessDpiAwareness(1)
except:
    pass


class Interface(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('EGFR QSAR Modeling Software')
        self.geometry('1400x900')
        self.resizable(False, False)
        self.configure(bg='#0F0F0F')

        script_dir = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
        icon_path = os.path.join(script_dir, 'Image Files', 'Icon.ico')

        # self.iconbitmap("Image Files/Icon.ico")
        self.iconbitmap(icon_path)

        table_frame = TableFrame(self)
        output_frame = OutputFrame(self)
        input_frame = InputFrame(self, table_frame, output_frame)
        output_frame.set_input_frame(input_frame)

        input_frame.pack(side='left', fill='both', expand=False, padx=15, pady=15)
        table_frame.pack(side='top', fill='both', expand=True, padx=15, pady=15)
        output_frame.pack(side='top', fill='both', expand=False, padx=15, pady=15)


if __name__ == '__main__':
    root = Interface()
    root.mainloop()
