import os
import tkinter as tk
from multiprocessing import cpu_count
from tkinter import filedialog
from tkinter import messagebox
from tkinter import ttk

import customtkinter as ctk
import numpy as np
import pandas as pd
from PIL import Image
from joblib import load

from DataOperations import compute_descriptors, descriptors_calculation


class InputFrame(tk.Frame):

    def on_tree_select(self, event):
        selected = self.treelist.selection()
        if selected:
            self.compute_button.configure(state=tk.NORMAL)
            self.delete_button.configure(state=tk.NORMAL)
        else:
            self.compute_button.configure(state=tk.DISABLED)
            self.delete_button.configure(state=tk.DISABLED)

    def upload_data(self):
        file_path = filedialog.askopenfilename()
        if file_path and file_path.endswith(".csv"):
            self.file_count += 1
            self.treelist.insert('', 'end', values=(self.file_count,
                                                    os.path.basename(file_path),
                                                    file_path))
        else:
            messagebox.showinfo("Error", "Please upload a .csv file")

    def delete_data(self):
        selected_items = self.treelist.selection()
        all_items = self.treelist.get_children('')
        remaining_data = []

        for item in all_items:
            if item in selected_items:
                self.file_count -= 1
            else:
                remaining_data.append(self.treelist.item(item)['values'][1:])

        self.treelist.delete(*all_items)

        for i, data in enumerate(remaining_data, start=1):
            self.treelist.insert('', 'end', values=(i, *data))

        remaining_after_delete = self.treelist.get_children()
        if remaining_after_delete:
            self.on_tree_select(None)
        else:
            self.compute_button.configure(state=tk.DISABLED)
            self.delete_button.configure(state=tk.DISABLED)

    def export_data(self):
        if self.computed_df is not None:
            filename = filedialog.asksaveasfilename(defaultextension='.csv',
                                                    filetypes=(('CSV files', '*.csv'),
                                                               ('All files', '*.*')))

            if filename:
                self.computed_df.to_csv(filename, index=False)
        else:
            tk.messagebox.showwarning('Data error', 'No data available for export!')

    def compute_data(self):
        selected_items = self.treelist.selection()
        num_processes = cpu_count()

        for item in selected_items:
            values = self.treelist.item(item, 'values')

            file_path = values[2]

            df = pd.read_csv(file_path)

            df.rename(columns={df.columns[0]: 'Smiles'}, inplace=True)
            descriptors_df = descriptors_calculation(df['Smiles'].tolist(), compute_descriptors, num_processes)

            model = load("MLModels/XGBClassifierEGFR.joblib")
            predictions = model.predict(descriptors_df)
            descriptors_df['Predicted Activity'] = predictions
            descriptors_df['Predicted Activity'] = descriptors_df['Predicted Activity'].map({0: 'Inactive',
                                                                                             1: 'Active'})

            num_smiles = len(df['Smiles'])
            df['Smile No.'] = np.arange(1, num_smiles + 1)

            self.computed_df = pd.concat([df['Smile No.'], descriptors_df], axis=1)
            self.output_frame.smiles_list = df.iloc[:, 0].tolist()
            self.output_frame.load_smiles()

        self.table_frame.show_dataframe(self.computed_df)

        self.compute_button.configure(state=tk.DISABLED)
        self.upload_button.configure(state=tk.DISABLED)
        self.export_button.configure(state=tk.NORMAL)
        self.visual_button.configure(state=tk.NORMAL)
        # self.analysis_button.configure(state=tk.NORMAL)

    def visualize_data(self):
        self.av_label.destroy()
        self.v_label = ctk.CTkLabel(self.av_frame,
                                    text='Enter the type of plot you would like to visualize:\n(Options '
                                         'include pariwise scatter plot, bar plot, and histogram plot)')
        self.v_label.pack(side='top', pady=10)

        self.v_entry = ctk.CTkEntry(self.av_frame, width=200)
        self.v_entry.pack(side='top', pady=10)

        self.v_button = ctk.CTkButton(self.av_frame,
                                      text='Enter',
                                      corner_radius=64,
                                      fg_color='#0F0F0F',
                                      border_color='#F6B17A',
                                      border_width=2,
                                      command=self.select_data)
        self.v_button.pack(side='top', pady=10)

    def select_data(self):
        entered_value = self.v_entry.get().lower()
        plot_types = ['bar plot', 'histogram plot', 'histogram', 'scatter plot', 'pairwise scatter plot']

        if entered_value in plot_types:
            if entered_value == 'scatter plot' or 'pairwise scatter plot':
                self.v_entry.delete(0, 'end')
                self.v_label.configure(text='Enter the two columns you want to compare against the target activity:')
                self.v_button.configure(text='Visualize')
            elif entered_value == 'bar plot' or 'histogram' or 'histogram plot':
                self.v_entry.delete(0, 'end')
        else:
            tk.messagebox.showerror('Invalid Input', 'Invalid plot type. '
                                                     'Please enter a valid type or check for spelling mistakes.')

    def analyze_data(self):
        self.av_label.destroy()
        self.a_label = ctk.CTkLabel(self.av_frame,
                                    text='Enter the type of analysis you would like to perform:')
        self.a_label.pack(side='top')

        self.a_entry = ctk.CTkEntry(self.av_frame, width=200)
        self.a_entry.pack(side='top')

        self.a_button = ctk.CTkButton(self.av_frame,
                                      text='Enter',
                                      corner_radius=64,
                                      fg_color='#0F0F0F',
                                      border_color='#F6B17A',
                                      border_width=2)
        self.a_button.pack(side='top')

    def __init__(self, container, table_frame=None, output_frame=None, **kwargs):
        super().__init__(container, **kwargs)
        self.a_button = None
        self.v_button = None
        self.a_entry = None
        self.v_entry = None
        self.a_label = None
        self.v_label = None
        self.computed_df = None
        self.table_frame = table_frame
        self.output_frame = output_frame
        self.configure(bg='#222831',
                       width=600,
                       height=700,
                       highlightcolor='white',
                       highlightbackground='white',
                       highlightthickness=1)

        upload_delete_frame = tk.Frame(self, bg='#222831')
        upload_delete_frame.pack(side='top', fill='x')

        udf_sep = ttk.Separator(self, orient='horizontal')
        udf_sep.pack(side='top', fill='x')

        upload_delete_frame.columnconfigure(0, weight=1)
        upload_delete_frame.columnconfigure(1, weight=1)

        upload_img = Image.open("Image Files/Upload.png")
        upload_img = upload_img.resize((50, 50))
        upload_img = ctk.CTkImage(upload_img)

        self.upload_button = ctk.CTkButton(upload_delete_frame,
                                           command=self.upload_data,
                                           corner_radius=64,
                                           fg_color='#0F0F0F',
                                           border_color='#F6B17A',
                                           border_width=2,
                                           text='Upload Data',
                                           image=upload_img)
        self.upload_button.grid(row=0, column=0, padx=20, pady=10)

        delete_img = Image.open("Image Files/Delete.png")
        delete_img = delete_img.resize((50, 50))
        delete_img = ctk.CTkImage(delete_img)

        self.delete_button = ctk.CTkButton(upload_delete_frame,
                                           command=self.delete_data,
                                           corner_radius=64,
                                           fg_color='#0F0F0F',
                                           border_color='#F6B17A',
                                           border_width=2,
                                           text='Delete Data',
                                           state=tk.DISABLED,
                                           image=delete_img)
        self.delete_button.grid(row=0, column=1, padx=20, pady=10)

        treelist_frame = tk.Frame(self, bg='#222831')
        treelist_frame.pack(side='top', fill='x', padx=10, pady=20)

        treelistframe_sep = ttk.Separator(self, orient='horizontal')
        treelistframe_sep.pack(side='top', fill='x')

        scrollbar = tk.Scrollbar(treelist_frame, borderwidth=0, relief='flat')
        scrollbar.pack(side='left', fill='y')

        self.file_count = 0

        self.treelist = ttk.Treeview(treelist_frame,
                                     columns=('File Index', 'File Name'),
                                     show='headings',
                                     yscrollcommand=scrollbar.set,
                                     selectmode='extended')
        self.treelist.bind('<<TreeviewSelect>>', self.on_tree_select)

        self.treelist.heading('File Index', text='File Index')
        self.treelist.heading('File Name', text='File Name')

        self.treelist.column('File Index', width=100, stretch=tk.NO, anchor='center')
        self.treelist.column('File Name', stretch=tk.YES, anchor='center')

        self.treelist.pack(side='left', fill='both', expand=True)

        scrollbar.configure(command=self.treelist.yview)

        compute_export_frame = tk.Frame(self, bg='#222831')
        compute_export_frame.pack(side='top', fill='x')

        compute_export_frame.columnconfigure(0, weight=1)
        compute_export_frame.columnconfigure(1, weight=1)

        compute_img = Image.open("Image Files/Compute.png")
        compute_img = compute_img.resize((50, 50))
        compute_img = ctk.CTkImage(compute_img)

        self.compute_button = ctk.CTkButton(compute_export_frame,
                                            corner_radius=64,
                                            fg_color='#0F0F0F',
                                            border_color='#F6B17A',
                                            border_width=2,
                                            text='Compute Data',
                                            state=tk.DISABLED,
                                            command=self.compute_data,
                                            image=compute_img)
        self.compute_button.grid(row=0, column=0, padx=20, pady=10)

        export_img = Image.open("Image Files/CSV.png")
        export_img = export_img.resize((50, 50))
        export_img = ctk.CTkImage(export_img)

        self.export_button = ctk.CTkButton(compute_export_frame,
                                           corner_radius=64, fg_color='#0F0F0F',
                                           border_color='#F6B17A',
                                           border_width=2,
                                           text='Export Data',
                                           state=tk.DISABLED,
                                           command=self.export_data,
                                           image=export_img)
        self.export_button.grid(row=0, column=1, padx=20, pady=10)

        self.analysis_visual_frame = tk.Frame(self, bg='#222831')
        self.analysis_visual_frame.pack(side='top', fill='x')

        avf_sep = ttk.Separator(self, orient='horizontal')
        avf_sep.pack(side='top', fill='x')

        self.analysis_visual_frame.columnconfigure(0, weight=1)
        self.analysis_visual_frame.columnconfigure(1, weight=1)

        analysis_img = Image.open("Image Files/Analysis.png")
        analysis_img = analysis_img.resize((50, 50))
        analysis_img = ctk.CTkImage(analysis_img)

        self.analysis_button = ctk.CTkButton(self.analysis_visual_frame,
                                             corner_radius=64,
                                             fg_color='#0F0F0F',
                                             border_color='#F6B17A',
                                             border_width=2,
                                             text='Data Analysis',
                                             state=tk.DISABLED,
                                             image=analysis_img,
                                             command=self.analyze_data)
        self.analysis_button.grid(row=0, column=0, padx=(30, 25), pady=10)

        visual_img = Image.open("Image Files/Visualization.png")
        visual_img = visual_img.resize((50, 50))
        visual_img = ctk.CTkImage(visual_img)

        self.visual_button = ctk.CTkButton(self.analysis_visual_frame,
                                           corner_radius=64,
                                           fg_color='#0F0F0F',
                                           border_color='#F6B17A',
                                           border_width=2,
                                           text='Data Visualization',
                                           state=tk.DISABLED,
                                           image=visual_img,
                                           command=self.visualize_data)
        self.visual_button.grid(row=0, column=1, padx=(30, 20), pady=10)

        self.av_frame = tk.Frame(self, bg='#222831')
        self.av_frame.pack(side='top', fill='both', expand=True)
        self.av_label = tk.Label(self.av_frame,
                                 bg='#222831',
                                 text='Data Analysis/Visualization Frame',
                                 font=('Default Font Family', 10),
                                 fg='#EEEEEE')
        self.av_label.pack(side='left', fill='both', expand=True)
