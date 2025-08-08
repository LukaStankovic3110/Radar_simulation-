import tkinter as tk
from tkinter import ttk

# Kreiraj glavni prozor
root = tk.Tk()
root.title("Target Configuration")

# Unos parametara
param_label = tk.Label(root, text="Unesite parametar:")
param_label.pack()

param_entry = tk.Entry(root)
param_entry.pack()

# Funkcija koja se poziva na klik OK dugmeta
def on_ok_click():
    param_value = param_entry.get()
    print(f"Uneti parametar: {param_value}")
    selected_targets = target_combobox.get()
    print(f"Broj targeta: {selected_targets}")

# OK dugme
ok_button = tk.Button(root, text="OK", command=on_ok_click)
ok_button.pack()

# Padajući meni za izbor broja targeta
target_label = tk.Label(root, text="Izaberite broj targeta:")
target_label.pack()

target_combobox = ttk.Combobox(root, values=[1, 2, 3, 4, 5])
target_combobox.pack()

# Pokreni GUI
root.mainloop()
