import tkinter as tk
from tkinter import messagebox, ttk
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageTk
import io
import customtkinter as ctk


class OutputFrame(tk.Frame):
    def __init__(self, container):
        super().__init__(container)
        self.configure(bg='#222831',
                       width=700,
                       height=350,
                       highlightcolor='white',
                       highlightbackground='white',
                       highlightthickness=1)
        self.pack_propagate(0)
        self.input_frame = None

        button_frame = tk.Frame(self, bg='#222831')
        button_frame.pack(side='top', fill='x')

        top_sep = ttk.Separator(self, orient='horizontal')
        top_sep.pack(side='top', expand=True, fill='x')

        previous_img = Image.open("Image Files/Previous.png")
        previous_img = previous_img.resize((50, 50))
        previous_img = ctk.CTkImage(previous_img)

        self.btn_previous = ctk.CTkButton(button_frame,
                                          text='Previous',
                                          command=self.prev_image,
                                          state=tk.DISABLED,
                                          corner_radius=64,
                                          fg_color='#0F0F0F',
                                          border_color='#F6B17A',
                                          border_width=2,
                                          image=previous_img)
        self.btn_previous.pack(side='left', padx=5, pady=10)

        next_img = Image.open("Image Files/Next.png")
        next_img = next_img.resize((50, 50))
        next_img = ctk.CTkImage(next_img)

        self.btn_next = ctk.CTkButton(button_frame,
                                      text='Next',
                                      command=self.next_image,
                                      state=tk.DISABLED,
                                      corner_radius=64,
                                      fg_color='#0F0F0F',
                                      border_color="#F6B17A",
                                      border_width=2,
                                      image=next_img)
        self.btn_next.pack(side='left', padx=5, pady=10)

        self.label = ctk.CTkLabel(button_frame, text='Smile No.')
        self.label.pack(side='left', padx=10)

        self.entry = ctk.CTkEntry(button_frame)
        self.entry.pack(side='left')
        self.entry.bind("<Return>", lambda event: self.search_image())

        search_img = Image.open('Image Files/search.png')
        search_img = search_img.resize((50, 50))
        search_img = ctk.CTkImage(search_img)

        self.btn_search = ctk.CTkButton(button_frame,
                                        text='Search',
                                        command=self.search_image,
                                        state=tk.DISABLED,
                                        corner_radius=64,
                                        fg_color='#0F0F0F',
                                        border_color='#F6B17A',
                                        border_width=2,
                                        image=search_img)
        self.btn_search.pack(side='left', padx=5, pady=10)

        image_frame = tk.Frame(self, width=300, height=300)
        image_frame.pack_propagate(0)
        image_frame.pack(side='left', padx=(35, 10), pady=(20, 20))
        self.image_label = tk.Label(image_frame,
                                    bg='#222831',
                                    text='2D Molecular Image Frame',
                                    font=('Default Font Family', 10),
                                    fg='#EEEEEE')
        self.image_label.pack(side='left', fill='both', expand=True)

        mid_sep = ttk.Separator(self, orient='vertical')
        mid_sep.pack(side='left', fill='y', padx=(35, 5))

        self.info_text = ctk.CTkTextbox(self, bg_color='#222831', state='disabled')
        self.info_text.pack(expand=True, fill='x', padx=10)

        self.smiles_list = []
        self.position = 0

    def search_image(self):
        self.info_text.configure(state='normal')
        input_value = self.entry.get()
        if not input_value.isdigit():
            messagebox.showerror('Error', 'Please enter a valid number!')
            self.entry.delete(0, 'end')
            return
        index = int(input_value) - 1
        if index < 0 or index >= len(self.smiles_list):
            messagebox.showerror('Error', 'Index out of bounds!')
            self.entry.delete(0, 'end')
            return

        self.position = index
        self.show_smiles(self.smiles_list[index])

        if index <= 0:
            self.btn_previous.configure(state=tk.DISABLED)
        else:
            self.btn_previous.configure(state=tk.NORMAL)

        if index >= len(self.smiles_list) - 1:
            self.btn_next.configure(state=tk.DISABLED)
        else:
            self.btn_next.configure(state=tk.NORMAL)

        if self.input_frame:
            text_content = self.info_text.get('1.0', 'end-1c')
            if text_content.strip():
                self.info_text.delete('1.0', 'end')
            descriptor_values = self.input_frame.computed_df
            row_values = descriptor_values.iloc[index]
            for descriptor, value in row_values.items():
                line = f'{descriptor}: {value}\n'
                self.info_text.insert('end', line)
            self.info_text.configure(state="disabled")

        self.entry.delete(0, 'end')

    def load_smiles(self):
        self.info_text.configure(state='normal')
        if len(self.smiles_list) > 0:
            self.show_smiles(self.smiles_list[0])

            if self.input_frame:
                descriptor_values = self.input_frame.computed_df
                row_values = descriptor_values.iloc[0]
                for descriptor, value in row_values.items():
                    line = f'{descriptor}: {value}\n'
                    self.info_text.insert('end', line)

            self.btn_previous.configure(state=tk.DISABLED)

            if len(self.smiles_list) > 1:
                self.btn_next.configure(state=tk.NORMAL)
            else:
                self.btn_next.configure(state=tk.DISABLED)

            self.btn_search.configure(state=tk.NORMAL)

    def prev_image(self):
        self.info_text.configure(state='normal')
        if self.position > 0:
            self.position -= 1
            self.show_smiles(self.smiles_list[self.position])

            if self.input_frame:
                text_content = self.info_text.get('1.0', 'end-1c')
                if text_content.strip():
                    self.info_text.delete('1.0', 'end')
                descriptor_values = self.input_frame.computed_df
                row_values = descriptor_values.iloc[self.position]
                for descriptor, value in row_values.items():
                    line = f'{descriptor}: {value}\n'
                    self.info_text.insert('end', line)
                self.info_text.configure(state="disabled")

            self.btn_next.configure(state=tk.NORMAL)

        if self.position == 0:
            self.btn_previous.configure(state=tk.DISABLED)

    def next_image(self):
        self.info_text.configure(state='normal')
        if self.position < len(self.smiles_list) - 1:
            self.position += 1
            self.show_smiles(self.smiles_list[self.position])

            if self.input_frame:
                text_content = self.info_text.get('1.0', 'end-1c')
                if text_content.strip():
                    self.info_text.delete('1.0', 'end')
                descriptor_values = self.input_frame.computed_df
                row_values = descriptor_values.iloc[self.position]
                for descriptor, value in row_values.items():
                    line = f'{descriptor}: {value}\n'
                    self.info_text.insert('end', line)
                self.info_text.configure(state="disabled")

            self.btn_previous.configure(state=tk.NORMAL)

        if self.position == len(self.smiles_list) - 1:
            self.btn_next.configure(state=tk.DISABLED)

    def show_smiles(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        img = Draw.MolToImage(mol)
        img_byte_arr = io.BytesIO()
        img.save(img_byte_arr, format='PNG')
        img_byte_arr = img_byte_arr.getvalue()
        photo = ImageTk.PhotoImage(Image.open(io.BytesIO(img_byte_arr)))
        self.image_label.configure(image=photo, text='')
        self.image_label.image = photo

    def set_input_frame(self, input_frame):
        self.input_frame = input_frame