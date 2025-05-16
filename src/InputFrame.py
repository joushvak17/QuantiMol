import os
import sys
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import ttk

import customtkinter as ctk
import numpy as np
import pandas as pd
import lightgbm
from molfeat.trans.fp import FPVecTransformer
from PIL import Image
from joblib import load
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from DataOperations import descriptors_calculation


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
            self.treelist.insert(
                "",
                "end",
                values=(self.file_count, os.path.basename(file_path), file_path),
            )
        else:
            messagebox.showinfo("Error", "Please upload a .csv file")

    def delete_data(self):
        selected_items = self.treelist.selection()
        all_items = self.treelist.get_children("")
        remaining_data = []

        for item in all_items:
            if item in selected_items:
                self.file_count -= 1
            else:
                remaining_data.append(self.treelist.item(item)["values"][1:])

        self.treelist.delete(*all_items)

        for i, data in enumerate(remaining_data, start=1):
            self.treelist.insert("", "end", values=(i, *data))

        remaining_after_delete = self.treelist.get_children()
        if remaining_after_delete:
            self.on_tree_select(None)
        else:
            self.compute_button.configure(state=tk.DISABLED)
            self.delete_button.configure(state=tk.DISABLED)

    def export_data(self):
        if self.computed_df is not None:
            filename = filedialog.asksaveasfilename(
                defaultextension=".csv",
                filetypes=(("CSV files", "*.csv"), ("All files", "*.*")),
            )

            if filename:
                self.computed_df.to_csv(filename, index=False)
        else:
            tk.messagebox.showwarning("Data error", "No data available for export!")

    def compute_data(self):
        self.compute_button.configure(state=tk.DISABLED)
        self.upload_button.configure(state=tk.DISABLED)
        self.delete_button.configure(state=tk.DISABLED)
        selected_items = self.treelist.selection()

        for item in selected_items:
            values = self.treelist.item(item, "values")

            file_path = values[2]

            df = pd.read_csv(file_path)

            df.rename(columns={df.columns[0]: "Smiles"}, inplace=True)
            # TODO: Uncomment the following lines to use molfeat featurizer
            transformer = FPVecTransformer(kind="desc2D")
            features = transformer(df["Smiles"])
            features_df = pd.DataFrame(features, index=df["Smiles"].index)
            descriptors_df = descriptors_calculation(df["Smiles"].tolist())
            model_df = pd.concat([descriptors_df, features_df], axis=1)
            # TODO: Uncomment the following lines to use lightgbm model
            descriptors_df = model_df

            # model_path = os.path.join(self.script_dir, 'MLModels', 'lgb_92_model.pkl')
            # model = load(model_path)
            # predictions = model.predict(model_df)
            # descriptors_df['Predicted Activity'] = predictions.round(0).astype(int)
            # descriptors_df['Predicted Activity'] = descriptors_df['Predicted Activity'].map({0: 'Inactive',
            #                                                                                  1: 'Active'})

            num_smiles = len(df["Smiles"])
            df["Smile No."] = np.arange(1, num_smiles + 1)

            self.computed_df = pd.concat([df["Smile No."], descriptors_df], axis=1)
            self.output_frame.smiles_list = df.iloc[:, 0].tolist()
            self.output_frame.load_smiles()

        self.table_frame.show_dataframe(self.computed_df)

        self.export_button.configure(state=tk.NORMAL)
        self.visual_button.configure(state=tk.NORMAL)
        self.analysis_button.configure(state=tk.NORMAL)

    def analyze_data(self):
        self.av_label.destroy()
        self.analysis_button.configure(state=tk.DISABLED)
        self.visual_button.configure(state=tk.DISABLED)
        self.a_label = ctk.CTkLabel(
            self.av_frame,
            text="Enter the type of analysis you would like to perform:\n(Only "
            "option as of right now is data summary.\n If "
            "you wish to do data summary and column,\n"
            ' enter "data summary, column name")',
        )
        self.a_entry = ctk.CTkEntry(self.av_frame, width=200)

        self.a_button = ctk.CTkButton(
            self.av_frame,
            text="Enter",
            corner_radius=64,
            fg_color="#0F0F0F",
            border_color="#F6B17A",
            border_width=2,
            command=self.v_select_data,
        )

        self.a_exit_button = ctk.CTkButton(
            self.av_frame,
            text="Exit",
            corner_radius=64,
            fg_color="#0F0F0F",
            border_color="#F6B17A",
            border_width=2,
            command=self.exit_analysis_data,
        )

        self.update_idletasks()
        self.a_label.pack(side="top", pady=10)
        self.a_entry.pack(side="top", pady=10)
        self.a_entry.bind("<Return>", lambda event: self.a_select_data())
        self.a_button.pack(side="top", pady=10)
        self.a_exit_button.pack(side="top", pady=10)

    def exit_analysis_data(self):
        if hasattr(self, "a_label"):
            self.a_label.pack_forget()
        if hasattr(self, "a_entry"):
            self.a_entry.pack_forget()
        if hasattr(self, "a_button"):
            self.a_button.pack_forget()
        if hasattr(self, "a_exit_button"):
            self.a_exit_button.pack_forget()
        if getattr(self, "tree_frame", None):
            self.tree_frame.pack_forget()

        self.update_idletasks()
        self.av_label = tk.Label(
            self.av_frame,
            bg="#222831",
            text="Data Analysis/Visualization Frame",
            font=("Default Font Family", 10),
            fg="#EEEEEE",
        )
        self.av_label.pack(side="left", fill="both", expand=True)
        self.analysis_button.configure(state=tk.NORMAL)
        self.visual_button.configure(state=tk.NORMAL)

    def a_select_data(self):
        entered_values = self.a_entry.get().lower().split(",")
        entered_values = [x.strip() for x in entered_values]

        self.lower_df = self.computed_df.copy()
        self.lower_df.columns = self.lower_df.columns.str.lower()
        columns = [col for col in self.lower_df.columns.tolist() if col != "smile no."]

        self.a_label.destroy()
        self.a_entry.destroy()
        self.a_button.destroy()

        self.tree_frame = tk.Frame(self.av_frame, width=480, height=320)
        self.tree_frame.pack_propagate(0)

        tree_scrollbar_x = tk.Scrollbar(self.tree_frame, orient="horizontal")

        tree = ttk.Treeview(
            self.tree_frame, xscrollcommand=tree_scrollbar_x.set, show="headings"
        )

        tree_scrollbar_x.config(command=tree.xview)

        self.tree_frame.pack()
        tree_scrollbar_x.pack(side="bottom", fill="x")
        tree.pack(fill="both", expand=True)

        if len(entered_values) == 1 and entered_values[0] == "data summary":
            data_summary = self.lower_df[columns].describe().round(2)
            tree["columns"] = ["summary"] + list(data_summary.columns)
            for i in tree["columns"]:
                tree.column(i, anchor="center", width=270)
                tree.heading(i, text=i, anchor="center")

            for index, row in data_summary.iterrows():
                tree.insert("", "end", values=[index] + row.to_list())
        elif len(entered_values) >= 2 and entered_values[0] == "data summary":
            for column_name in entered_values[1:]:
                column_name = column_name.lower()
                if column_name in columns:
                    column_summary = self.lower_df[column_name].describe().round(2)
                    tree["columns"] = ["summary", column_name]
                    for i in tree["columns"]:
                        tree.column(i, anchor="center", width=200)
                        tree.heading(i, text=i, anchor="center")

                    for index, value in column_summary.items():
                        tree.insert("", "end", values=[index, value])
                else:
                    tk.messagebox.showwarning(
                        "Data error", "No data available for analysis."
                    )
                    self.a_entry.delete(0, "end")
        else:
            tk.messagebox.showwarning(
                "Invalid Input", "Please enter a valid data analysis input."
            )
            self.a_entry.delete(0, "end")

    def visualize_data(self):
        self.av_label.destroy()
        self.visual_button.configure(state=tk.DISABLED)
        self.analysis_button.configure(state=tk.DISABLED)
        self.v_label = ctk.CTkLabel(
            self.av_frame,
            text="Enter the type of plot you would like to visualize:\n(Options "
            "include pairwise scatter plot, box plot, and histogram plot)",
        )

        self.v_entry = ctk.CTkEntry(self.av_frame, width=200)

        self.v_button = ctk.CTkButton(
            self.av_frame,
            text="Enter",
            corner_radius=64,
            fg_color="#0F0F0F",
            border_color="#F6B17A",
            border_width=2,
            command=self.v_select_data,
        )

        self.v_exit_button = ctk.CTkButton(
            self.av_frame,
            text="Exit",
            corner_radius=64,
            fg_color="#0F0F0F",
            border_color="#F6B17A",
            border_width=2,
            command=self.exit_visualize_data,
        )

        self.update_idletasks()
        self.v_label.pack(side="top", pady=10)
        self.v_entry.pack(side="top", pady=10)
        self.v_entry.bind("<Return>", lambda event: self.v_select_data())
        self.v_button.pack(side="top", pady=10)
        self.v_exit_button.pack(side="top", pady=10)

    def exit_visualize_data(self):
        self.after_cancel(self.visualize_data)
        self.after_cancel(self.v_select_data)
        self.after_cancel(self.hist)
        self.after_cancel(self.scatter)
        self.after_cancel(self.box)

        if hasattr(self, "v_label"):
            self.v_label.destroy()
        if hasattr(self, "v_entry"):
            self.v_entry.destroy()
        if hasattr(self, "v_button"):
            self.v_button.destroy()
        if hasattr(self, "v_exit_button"):
            self.v_exit_button.destroy()

        self.update_idletasks()
        self.av_label = tk.Label(
            self.av_frame,
            bg="#222831",
            text="Data Analysis/Visualization Frame",
            font=("Default Font Family", 10),
            fg="#EEEEEE",
        )
        self.av_label.pack(side="left", fill="both", expand=True)
        self.visual_button.configure(state=tk.NORMAL)
        self.analysis_button.configure(state=tk.NORMAL)

    def v_select_data(self):
        entered_value = self.v_entry.get().lower()
        plot_types = [
            "box plot",
            "boxplot",
            "histogram plot",
            "histogram",
            "scatter plot",
            "pairwise scatter plot",
            "scatter",
        ]

        if entered_value in plot_types:
            self.v_button.configure(state=tk.DISABLED)
            if (
                entered_value == "scatter plot"
                or entered_value == "pairwise scatter plot"
            ):
                self.v_button.configure(command=lambda: None)
                self.v_entry.delete(0, "end")
                self.v_entry.configure(width=300)
                self.v_entry.unbind("<Return>")
                self.v_entry.bind("<Return>", lambda event: self.scatter())
                self.v_label.configure(
                    text="Enter the two columns you want to compare against the target "
                    "activity:\n(Seperate the two columns by comma)"
                )
                self.v_button.configure(
                    text="Visualize", state=tk.NORMAL, command=self.scatter
                )
            elif entered_value == "histogram" or entered_value == "histogram plot":
                self.v_button.configure(command=lambda: None)
                self.v_entry.delete(0, "end")
                self.v_entry.configure(width=200)
                self.v_entry.unbind("<Return>")
                self.v_entry.bind("<Return>", lambda event: self.hist())
                self.v_label.configure(
                    text="Enter the column you want to visualize a histogram with:"
                )
                self.v_button.configure(
                    text="Visualize", state=tk.NORMAL, command=self.hist
                )
            else:
                self.v_button.configure(command=lambda: None)
                self.v_entry.delete(0, "end")
                self.v_entry.configure(width=200)
                self.v_entry.unbind("<Return>")
                self.v_entry.bind("<Return>", lambda event: self.box())
                self.v_label.configure(
                    text="Enter the column you want to visualize a box plot with:"
                )
                self.v_button.configure(
                    text="Visualize", state=tk.NORMAL, command=self.box
                )
        else:
            tk.messagebox.showerror(
                "Invalid Input",
                "Invalid plot type. "
                "Please enter a valid type or check for spelling mistakes.",
            )
            self.v_entry.delete(0, "end")

    def hist(self):
        entered_value = self.v_entry.get().lower().strip()

        self.temp_df = self.computed_df.copy()
        self.temp_df.columns = map(str.lower, self.temp_df.columns)

        if entered_value not in self.temp_df:
            tk.messagebox.showerror(
                "Invalid Input", "The columns specified were not found in the data."
            )
            self.v_entry.delete(0, "end")
            return

        fig = Figure(figsize=(7, 4), dpi=100)
        ax = fig.add_subplot(111)
        fig.subplots_adjust(bottom=0.2)
        ax.hist(self.temp_df[entered_value], bins=30, edgecolor="black")
        ax.set_xlabel(entered_value)
        ax.set_ylabel("Frequency")
        ax.set_title(f"Histogram of {entered_value.upper()}")

        hist_win = tk.Toplevel(self)
        hist_win.wm_title("Histogram")
        icon_path = os.path.join(self.script_dir, "Image Files", "Icon.ico")
        hist_win.iconbitmap(icon_path)

        canvas = FigureCanvasTkAgg(fig, master=hist_win)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

        toolbar = NavigationToolbar2Tk(canvas, hist_win)
        toolbar.update()

        def on_resize(event):
            width, height = event.width, event.height
            dpi = fig.get_dpi()
            fig.set_size_inches(width / dpi, height / dpi)
            canvas.draw()

        hist_win.bind("<Configure>", on_resize)

    def box(self):
        entered_value = self.v_entry.get().lower().strip()

        self.temp_df = self.computed_df.copy()
        self.temp_df.columns = map(str.lower, self.temp_df.columns)

        if entered_value not in self.temp_df.columns:
            tk.messagebox.showerror(
                "Invalid Input", "The column specified was not found in the data."
            )
            self.v_entry.delete(0, "end")
            return

        fig = Figure(figsize=(7, 4), dpi=100)
        ax = fig.add_subplot(111)
        fig.subplots_adjust(bottom=0.2)
        ax.boxplot(
            self.temp_df[entered_value],
            patch_artist=True,
            boxprops=dict(edgecolor="black"),
        )
        ax.set_xlabel(entered_value)
        ax.set_ylabel("Distribution")
        ax.set_title(f"Boxplot of {entered_value.upper()}")

        box_win = tk.Toplevel(self)
        box_win.wm_title("Boxplot")
        icon_path = os.path.join(self.script_dir, "Image Files", "Icon.ico")
        box_win.iconbitmap(icon_path)

        canvas = FigureCanvasTkAgg(fig, master=box_win)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

        toolbar = NavigationToolbar2Tk(canvas, box_win)
        toolbar.update()

        def on_resize(event):
            width, height = event.width, event.height
            dpi = fig.get_dpi()
            fig.set_size_inches(width / dpi, height / dpi)
            canvas.draw()

        box_win.bind("<Configure>", on_resize)

    def scatter(self):
        entered_columns = self.v_entry.get().split(",")
        if len(entered_columns) == 2:
            column_1 = entered_columns[0].lower().strip()
            column_2 = entered_columns[1].lower().strip()

            self.temp_df = self.computed_df.copy()
            self.temp_df.columns = map(str.lower, self.temp_df.columns)

            if column_1 not in self.temp_df or column_2 not in self.temp_df:
                tk.messagebox.showerror(
                    "Invalid Input",
                    "One or both of the columns specified were not found in the data.",
                )
                self.v_entry.delete(0, "end")
                return

            fig = Figure(figsize=(7, 4), dpi=100)
            ax = fig.add_subplot(111)
            fig.subplots_adjust(bottom=0.3)
            ax.scatter(self.temp_df[column_1], self.temp_df[column_2])
            ax.set_xlabel(column_1.upper())
            ax.set_ylabel(column_2.upper())
            ax.set_title(f"Scatter Plot of {column_1.upper()} vs {column_2.upper()}")

            scatter_win = tk.Toplevel(self)
            scatter_win.wm_title("Scatter Plot")
            icon_path = os.path.join(self.script_dir, "Image Files", "Icon.ico")
            scatter_win.iconbitmap(icon_path)

            canvas = FigureCanvasTkAgg(fig, master=scatter_win)
            canvas.draw()
            canvas.get_tk_widget().pack(fill="both", expand=True)

            toolbar = NavigationToolbar2Tk(canvas, scatter_win)
            toolbar.update()

            def on_resize(event):
                width, height = event.width, event.height
                dpi = fig.get_dpi()
                fig.set_size_inches(width / dpi, height / dpi)
                canvas.draw()

            scatter_win.bind("<Configure>", on_resize)

        else:
            tk.messagebox.showerror(
                "Invalid Input",
                "Please enter exactly two column names "
                "separated by a comma. Check for spelling mistakes.",
            )
            self.v_entry.delete(0, "end")

    def __init__(self, container, table_frame=None, output_frame=None, **kwargs):
        super().__init__(container, **kwargs)
        self.lower_df = None
        self.tree_frame = None
        self.a_exit_button = None
        self.v_exit_button = None
        self.temp_df = None
        self.a_button = None
        self.v_button = None
        self.a_entry = None
        self.v_entry = None
        self.a_label = None
        self.v_label = None
        self.computed_df = None
        self.table_frame = table_frame
        self.output_frame = output_frame
        self.configure(
            bg="#222831",
            width=600,
            height=700,
            highlightcolor="white",
            highlightbackground="white",
            highlightthickness=1,
        )
        self.script_dir = getattr(
            sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__))
        )

        upload_delete_frame = tk.Frame(self, bg="#222831")
        upload_delete_frame.pack(side="top", fill="x")

        udf_sep = ttk.Separator(self, orient="horizontal")
        udf_sep.pack(side="top", fill="x")

        upload_delete_frame.columnconfigure(0, weight=1)
        upload_delete_frame.columnconfigure(1, weight=1)

        upload_img = os.path.join(self.script_dir, "Image Files", "Upload.png")
        upload_img = Image.open(upload_img)
        upload_img = upload_img.resize((50, 50))
        upload_img = ctk.CTkImage(upload_img)

        self.upload_button = ctk.CTkButton(
            upload_delete_frame,
            command=self.upload_data,
            corner_radius=64,
            fg_color="#0F0F0F",
            border_color="#F6B17A",
            border_width=2,
            text="Upload Data",
            image=upload_img,
        )
        self.upload_button.grid(row=0, column=0, padx=20, pady=10)

        delete_img = os.path.join(self.script_dir, "Image Files", "Delete.png")
        delete_img = Image.open(delete_img)
        delete_img = delete_img.resize((50, 50))
        delete_img = ctk.CTkImage(delete_img)

        self.delete_button = ctk.CTkButton(
            upload_delete_frame,
            command=self.delete_data,
            corner_radius=64,
            fg_color="#0F0F0F",
            border_color="#F6B17A",
            border_width=2,
            text="Delete Data",
            state=tk.DISABLED,
            image=delete_img,
        )
        self.delete_button.grid(row=0, column=1, padx=20, pady=10)

        treelist_frame = tk.Frame(self, bg="#222831")
        treelist_frame.pack(side="top", fill="x", padx=10, pady=20)

        treelistframe_sep = ttk.Separator(self, orient="horizontal")
        treelistframe_sep.pack(side="top", fill="x")

        scrollbar = tk.Scrollbar(treelist_frame, borderwidth=0, relief="flat")
        scrollbar.pack(side="left", fill="y")

        self.file_count = 0

        self.treelist = ttk.Treeview(
            treelist_frame,
            columns=("File Index", "File Name"),
            show="headings",
            yscrollcommand=scrollbar.set,
            selectmode="extended",
        )
        self.treelist.bind("<<TreeviewSelect>>", self.on_tree_select)

        self.treelist.heading("File Index", text="File Index")
        self.treelist.heading("File Name", text="File Name")

        self.treelist.column("File Index", width=100, stretch=tk.NO, anchor="center")
        self.treelist.column("File Name", stretch=tk.YES, anchor="center")

        self.treelist.pack(side="left", fill="both", expand=True)

        scrollbar.configure(command=self.treelist.yview)

        compute_export_frame = tk.Frame(self, bg="#222831")
        compute_export_frame.pack(side="top", fill="x")

        compute_export_frame.columnconfigure(0, weight=1)
        compute_export_frame.columnconfigure(1, weight=1)

        compute_img = os.path.join(self.script_dir, "Image Files", "Compute.png")
        compute_img = Image.open(compute_img)
        compute_img = compute_img.resize((50, 50))
        compute_img = ctk.CTkImage(compute_img)

        self.compute_button = ctk.CTkButton(
            compute_export_frame,
            corner_radius=64,
            fg_color="#0F0F0F",
            border_color="#F6B17A",
            border_width=2,
            text="Compute Data",
            state=tk.DISABLED,
            command=self.compute_data,
            image=compute_img,
        )
        self.compute_button.grid(row=0, column=0, padx=20, pady=10)

        export_img = os.path.join(self.script_dir, "Image Files", "CSV.png")
        export_img = Image.open(export_img)
        export_img = export_img.resize((50, 50))
        export_img = ctk.CTkImage(export_img)

        self.export_button = ctk.CTkButton(
            compute_export_frame,
            corner_radius=64,
            fg_color="#0F0F0F",
            border_color="#F6B17A",
            border_width=2,
            text="Export Data",
            state=tk.DISABLED,
            command=self.export_data,
            image=export_img,
        )
        self.export_button.grid(row=0, column=1, padx=20, pady=10)

        self.analysis_visual_frame = tk.Frame(self, bg="#222831")
        self.analysis_visual_frame.pack(side="top", fill="x")

        avf_sep = ttk.Separator(self, orient="horizontal")
        avf_sep.pack(side="top", fill="x")

        self.analysis_visual_frame.columnconfigure(0, weight=1)
        self.analysis_visual_frame.columnconfigure(1, weight=1)

        analysis_img = os.path.join(self.script_dir, "Image Files", "Analysis.png")
        analysis_img = Image.open(analysis_img)
        analysis_img = analysis_img.resize((50, 50))
        analysis_img = ctk.CTkImage(analysis_img)

        self.analysis_button = ctk.CTkButton(
            self.analysis_visual_frame,
            corner_radius=64,
            fg_color="#0F0F0F",
            border_color="#F6B17A",
            border_width=2,
            text="Data Analysis",
            state=tk.DISABLED,
            image=analysis_img,
            command=self.analyze_data,
        )
        self.analysis_button.grid(row=0, column=0, padx=(30, 25), pady=10)

        visual_img = os.path.join(self.script_dir, "Image Files", "Visualization.png")
        visual_img = Image.open(visual_img)
        visual_img = visual_img.resize((50, 50))
        visual_img = ctk.CTkImage(visual_img)

        self.visual_button = ctk.CTkButton(
            self.analysis_visual_frame,
            corner_radius=64,
            fg_color="#0F0F0F",
            border_color="#F6B17A",
            border_width=2,
            text="Data Visualization",
            state=tk.DISABLED,
            image=visual_img,
            command=self.visualize_data,
        )
        self.visual_button.grid(row=0, column=1, padx=(30, 20), pady=10)

        self.av_frame = tk.Frame(self, bg="#222831")
        self.av_frame.pack(side="top", fill="both", expand=True)
        self.av_label = tk.Label(
            self.av_frame,
            bg="#222831",
            text="Data Analysis/Visualization Frame",
            font=("Default Font Family", 10),
            fg="#EEEEEE",
        )
        self.av_label.pack(side="left", fill="both", expand=True)
