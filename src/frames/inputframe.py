import tkinter as tk
import customtkinter as ctk
from tkinter import ttk


class InputFrame(tk.Frame):
    def __init__(
        self,
        master=None,
        cnf=None,
        *,
        background=None,
        bd=0,
        bg="#161b22",
        border=0,
        borderwidth=0,
        class_="Frame",
        colormap="",
        container=False,
        cursor="",
        height=None,
        highlightbackground="white",
        highlightcolor="white",
        highlightthickness=1,
        name=None,
        padx=0,
        pady=0,
        relief="flat",
        takefocus=0,
        visual="",
        width=350,
    ):
        super().__init__(
            master,
            cnf,
            background=background,
            bd=bd,
            bg=bg,
            border=border,
            borderwidth=borderwidth,
            class_=class_,
            colormap=colormap,
            container=container,
            cursor=cursor,
            height=height,
            highlightbackground=highlightbackground,
            highlightcolor=highlightcolor,
            highlightthickness=highlightthickness,
            name=name,
            padx=padx,
            pady=pady,
            relief=relief,
            takefocus=takefocus,
            visual=visual,
            width=width,
        )

        self.propagate(False)

        # NOTE: Data Input Frame - This frame is used to upload and delete data.
        data_input_frame = tk.Frame(self, bg=bg)
        data_input_frame.pack(side="top", fill="x")
        data_input_frame.columnconfigure(0, weight=1)
        data_input_frame.columnconfigure(1, weight=1)
        self.upload_button = ctk.CTkButton(
            data_input_frame,
            height=10,
            fg_color="#0F0F0F",
            border_color="#40916c",
            border_width=2,
            text="Upload Data",
        )
        self.upload_button.grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        self.delete_button = ctk.CTkButton(
            data_input_frame,
            height=10,
            fg_color="#0F0F0F",
            border_color="#40916c",
            border_width=2,
            text="Delete Data",
            state="disabled",
        )
        self.delete_button.grid(row=0, column=1, padx=10, pady=10, sticky="ew")
        dif_sep = ttk.Separator(self, orient="horizontal")
        dif_sep.pack(side="top", fill="x")
