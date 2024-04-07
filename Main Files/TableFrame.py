import tkinter as tk
from tkinter import ttk

import pandas as pd


class TableFrame(tk.Frame):
    def __init__(self, container, **kwargs):
        super().__init__(container, **kwargs)

        self.configure(width=700, height=350)

        vsb = tk.Scrollbar(self, orient='vertical', borderwidth=0, relief='flat')
        vsb.pack(side='left', fill='y')

        hsb = tk.Scrollbar(self, orient='horizontal', borderwidth=0, relief='flat')
        hsb.pack(side='bottom', fill='x')

        self.tree = ttk.Treeview(self, yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.tree.pack(fill='both', expand=True)

        vsb.config(command=self.tree.yview)
        hsb.config(command=self.tree.xview)

    def show_dataframe(self, df: pd.DataFrame):
        self.tree['columns'] = list(df.columns)
        self.tree['show'] = 'headings'

        for col in self.tree.get_children():
            self.tree.delete(col)

        for i, col_name in enumerate(df.columns):
            max_length = max(df[col_name].map(str).map(len).max(), len(str(col_name)))

            column_width = max_length * 8

            self.tree.column(i, width=column_width, stretch=tk.NO)
            self.tree.heading(i, text=col_name)

        for row in df.itertuples():
            self.tree.insert('', 'end', values=row[1:])
