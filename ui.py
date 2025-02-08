import tkinter as tk
from tkinter import ttk

def update_speeds():
    selected_modules = [module for module, var in module_vars.items() if var.get()]
    for module in selected_modules:
        print(f"Updating speed for {module} to {speed_entries[module].get()} m/s")

# Create the main window
root = tk.Tk()
root.title("NEXAAS Wind")

# Target velocity display
target_velocity_label = tk.Label(root, text="Target velocity (m/s)", font=('Arial', 14))
target_velocity_label.grid(row=0, column=0, padx=10, pady=10)

target_velocity_value = tk.Label(root, text="0.00", font=('Arial', 14))
target_velocity_value.grid(row=0, column=1, padx=10, pady=10)

# Module selection
module_frame = tk.LabelFrame(root, text="Modules", padx=10, pady=10)
module_frame.grid(row=1, column=0, columnspan=2, padx=10, pady=10)

modules = ["Module 1", "Module 2", "Module 3", "Module 4"]
module_vars = {}
speed_entries = {}

for i, module in enumerate(modules):
    module_vars[module] = tk.BooleanVar()
    chk = tk.Checkbutton(module_frame, text=module, variable=module_vars[module])
    chk.grid(row=i, column=0, sticky='w')
    
    speed_entry = tk.Entry(module_frame, width=10)
    speed_entry.grid(row=i, column=1, padx=10)
    speed_entries[module] = speed_entry

# Buttons for module selection
select_all_button = tk.Button(module_frame, text="Select All", command=lambda: [var.set(True) for var in module_vars.values()])
select_all_button.grid(row=len(modules), column=0, pady=5)

deselect_all_button = tk.Button(module_frame, text="Deselect All", command=lambda: [var.set(False) for var in module_vars.values()])
deselect_all_button.grid(row=len(modules), column=1, pady=5)

# Gradient settings
gradient_frame = tk.LabelFrame(root, text="Gradient", padx=10, pady=10)
gradient_frame.grid(row=2, column=0, columnspan=2, padx=10, pady=10)

gradient_type_label = tk.Label(gradient_frame, text="Gradient Type")
gradient_type_label.grid(row=0, column=0, sticky='w')

gradient_type_combobox = ttk.Combobox(gradient_frame, values=["Linear", "Radial"])
gradient_type_combobox.grid(row=0, column=1, padx=10, pady=5)
gradient_type_combobox.current(0)

direction_x_label = tk.Label(gradient_frame, text="Direction X:")
direction_x_label.grid(row=1, column=0, sticky='w')

direction_x_entry = tk.Entry(gradient_frame, width=10)
direction_x_entry.grid(row=1, column=1, padx=10, pady=5)
direction_x_entry.insert(0, "0.00")

direction_y_label = tk.Label(gradient_frame, text="Y:")
direction_y_label.grid(row=2, column=0, sticky='w')

direction_y_entry = tk.Entry(gradient_frame, width=10)
direction_y_entry.grid(row=2, column=1, padx=10, pady=5)
direction_y_entry.insert(0, "0.00")

# Update button
update_button = tk.Button(root, text="Update Module Speeds", command=update_speeds)
update_button.grid(row=3, column=0, columnspan=2, pady=10)

# Start the Tkinter event loop
root.mainloop()