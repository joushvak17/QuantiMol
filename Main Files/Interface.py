import tkinter as tk

try:
    from ctypes import windll

    windll.shcore.SetProcessDpiAwareness(1)
except:
    pass


class Interface(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('EGFR QSAR Modeling Software')
        self.geometry('1400x800')
        self.resizable(False, False)
        self.configure(bg='#0F0F0F')
        self.iconbitmap("Image Files/icon.ico")


root = Interface()
root.mainloop()
