import os
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import ttk

import customtkinter as ctk
import numpy as np
import pandas as pd
from joblib import load
from PIL import Image
from rdkit import Chem
from rdkit.Chem import (rdMolDescriptors as
                        rdmd,
                        GraphDescriptors,
                        Descriptors,
                        FindMolChiralCenters)


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
        def is_csv(file_path):
            _, extension = os.path.splitext(file_path)
            return extension.lower() == '.csv'

        file_path = filedialog.askopenfilename()
        if file_path:
            if is_csv(file_path):
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

        for item in selected_items:
            values = self.treelist.item(item, 'values')

            file_path = values[2]

            df = pd.read_csv(file_path)

            def descriptors_calculation(smiles_list):
                mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

                wiener_res = []
                amat_list = [Chem.GetDistanceMatrix(mol) for mol in mol_list]
                for i, mol in enumerate(mol_list):
                    res = 0
                    amat = amat_list[i]
                    num_atoms = mol.GetNumAtoms()
                    for j in range(num_atoms):
                        for k in range(j + 1, num_atoms):
                            res += amat[j][k]
                    wiener_res.append(res)

                descriptor_data = {
                    'Molecular Weight': [Descriptors.ExactMolWt(mol) for mol in mol_list],
                    'Number of Rotatable Bonds': [rdmd.CalcNumRotatableBonds(mol) for mol in mol_list],
                    'Number of Atoms': [mol.GetNumAtoms() for mol in mol_list],
                    'Number of Bonds': [mol.GetNumBonds() for mol in mol_list],
                    'Count of Chiral Centers': [len(FindMolChiralCenters(mol, includeUnassigned=True)) for mol in
                                                mol_list],
                    'Number of Rings': [rdmd.CalcNumRings(mol) for mol in mol_list],
                    'Number of Aromatic Rings': [rdmd.CalcNumAromaticRings(mol) for mol in mol_list],
                    'Number of Hydrogen Bond Donors': [rdmd.CalcNumHBD(mol) for mol in mol_list],
                    'Number of Hydrogen Bond Acceptors': [rdmd.CalcNumHBA(mol) for mol in mol_list],
                    'Balaban J Index': [GraphDescriptors.BalabanJ(mol) for mol in mol_list],
                    'Wiener Index': wiener_res,
                    'LogP': [rdmd.CalcCrippenDescriptors(mol)[0] for mol in mol_list],
                    'TPSA': [rdmd.CalcTPSA(mol) for mol in mol_list],
                }

                descriptor_values = pd.DataFrame(descriptor_data)

                return descriptor_values.round(2)

            df.rename(columns={df.columns[0]: 'Smiles'}, inplace=True)
            descriptors_df = descriptors_calculation(df['Smiles'].tolist())

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

        self.table_frame.show_dataframe(self.computed_df)

        self.update_idletasks()
        self.compute_button.configure(state=tk.DISABLED)
        self.upload_button.configure(state=tk.DISABLED)
        self.export_button.configure(state=tk.NORMAL)

    def __init__(self, container, table_frame=None, output_frame=None, **kwargs):
        super().__init__(container, **kwargs)
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

        analysis_visual_frame = tk.Frame(self, bg='#222831')
        analysis_visual_frame.pack(side='top', fill='x')

        avf_sep = ttk.Separator(self, orient='horizontal')
        avf_sep.pack(side='top', fill='x')

        analysis_visual_frame.columnconfigure(0, weight=1)
        analysis_visual_frame.columnconfigure(1, weight=1)

        analysis_img = Image.open("Image Files/Analysis.png")
        analysis_img = analysis_img.resize((50, 50))
        analysis_img = ctk.CTkImage(analysis_img)

        self.analysis_button = ctk.CTkButton(analysis_visual_frame,
                                             corner_radius=64,
                                             fg_color='#0F0F0F',
                                             border_color='#F6B17A',
                                             border_width=2,
                                             text='Data Analysis',
                                             state=tk.DISABLED,
                                             image=analysis_img)
        self.analysis_button.grid(row=0, column=0, padx=(30, 25), pady=10)

        visual_img = Image.open("Image Files/Visualization.png")
        visual_img = visual_img.resize((50, 50))
        visual_img = ctk.CTkImage(visual_img)

        self.visual_button = ctk.CTkButton(analysis_visual_frame,
                                           corner_radius=64,
                                           fg_color='#0F0F0F',
                                           border_color='#F6B17A',
                                           border_width=2,
                                           text='Data Visualization',
                                           state=tk.DISABLED,
                                           image=visual_img)
        self.visual_button.grid(row=0, column=1, padx=(30, 20), pady=10)

        av_frame = tk.Frame(self, bg='#222831')
        av_frame.pack(side='top', fill='both', expand=True)
        self.av_label = tk.Label(av_frame,
                                 bg='#222831',
                                 text='Data Analysis/Visualization Frame',
                                 font=('Default Font Family', 10),
                                 fg='#EEEEEE')
        self.av_label.pack(side='left', fill='both', expand=True)
