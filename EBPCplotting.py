import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from concurrent.futures import ThreadPoolExecutor
import matplotlib.colors as mcolors 
import matplotlib.dates as mdates
import numpy as np

def process_species(species):
    species_data = data[data['current_sci_name '] == species]
    status_at_dates = {}
    
    for date in all_dates:
        status_at_date = species_data[species_data['date_effective'] <= date].sort_values('date_effective').iloc[-1, [3, 4]] if not species_data[species_data['date_effective'] <= date].empty else None
        status_at_dates[date] = status_at_date
    
    return species, status_at_dates

data_file_path = r"data"  # Replace with your data file path
date_file_path = r"date"  # Replace with your all_dates file path
output_dir_path = r"output"  # Replace with your desired output directory path

print("Reading data from file...")
data = pd.read_csv(data_file_path)
data['date_effective'] = pd.to_datetime(data['date_effective'], format="%d/%m/%Y")

print("Reading all_dates from file...")
all_dates = pd.read_csv(date_file_path, header=None, nrows=1)
all_dates = pd.to_datetime(all_dates.stack().values, format="%d/%m/%Y")

print("Finding unique species...")
unique_species = data['current_sci_name '].unique()

# Only plot the first x amount of species
unique_species = unique_species[:200]

species_status = {}

print("Processing species data...")
with ThreadPoolExecutor() as executor:
    results = executor.map(process_species, unique_species)
    for species, status_at_dates in results:
        print(f"  Processed {species}...")
        species_status[species] = status_at_dates

# Color coding function for bars based on conservation status
def get_bar_color(status):
    color_map = {
        'Extinct': 'black',
        'Extinct in the wild': '#660066', #dark purple
        'Critically Endangered': '#ff0000', #red
        'Endangered': '#ed7727', # orange
        'Vulnerable': '#ffd966', # light yellow
        'Conservation Dependent': '#305496', # blue
    }
    return color_map.get(status, 'black')

def find_continuous_segments(species_data, all_dates):
    segments = []
    segment_start = None
    prev_status = None

    for idx, date in enumerate(all_dates[:-1]):
        next_date = all_dates[idx + 1]
        if species_data['date_effective'].isin(pd.date_range(date, next_date, closed='left')).any():
            current_status = species_data[species_data['date_effective'] <= date].sort_values('date_effective', ascending=False).iloc[0]['conservation_staus']
            if segment_start is None:
                segment_start = date
                prev_status = current_status
            elif prev_status != current_status:
                segments.append((segment_start, date))
                segment_start = date
                prev_status = current_status
        else:
            if segment_start is not None:
                segments.append((segment_start, date))
                segment_start = None

    if segment_start is not None:
        segments.append((segment_start, all_dates[-1]))

    return segments

# Plotting
species_index = {species: idx for idx, species in enumerate(unique_species)}

num_species = len(unique_species)
species_per_page = 40
num_pages = int(np.ceil(num_species / species_per_page))

with PdfPages(output_dir_path + '\\multipage_plot.pdf') as pdf:
    #for page in range(num_pages - 1, num_pages):  # Iterate only through the last page
    for page in range(num_pages):  # Iterate through all the pages
        species_subset = unique_species[page * species_per_page:(page + 1) * species_per_page]
        species_index = {species: idx for idx, species in enumerate(species_subset)}

        # Calculate the number of empty bars needed on the last page
        if page == num_pages - 1:
            num_empty_bars = (species_per_page * num_pages) - num_species
        else:
            num_empty_bars = 0

        print(f"Plotting page {page + 1}...")
        fig, ax = plt.subplots(figsize=(8.27, 11.69), dpi=300)  # Set figsize to A4 dimensions and DPI to 300

        min_date = data['date_effective'].min() - pd.Timedelta(days=1)
        max_date = data['date_effective'].max() + pd.Timedelta(days=1)

        # Add empty bars to fill the remaining spaces on the last page
        for i in range(num_empty_bars):
            empty_bar_index = i
            ax.broken_barh([(min_date, max_date - min_date)], (empty_bar_index - 0.4, 0.8), facecolor='none', linewidth=0)


        for species in species_subset:
            species_data = data[data['current_sci_name '] == species]
            segments = find_continuous_segments(species_data, all_dates)

            for start_date, end_date in segments:
                duration = end_date - start_date
                status = species_data[species_data['date_effective'] <= start_date].sort_values('date_effective', ascending=False).iloc[0]['conservation_staus']

                bar_color = get_bar_color(status)
                ax.broken_barh([(start_date, duration)], (species_index[species] - 0.4 + num_empty_bars, 0.8), color=bar_color, linewidth=0)

        # Add black solid lines between each species
        for i in range(len(species_subset) - 1):
            ax.axhline(i + 0.5 + num_empty_bars, color='black', linestyle='solid')

        ax.set_xlabel("Date")
        #ax.set_ylabel("Listed Scientific Name and Conservation Status")

        # Center the title at the top of the page
        title = f"Species Conservation Status Over Time (Species {page * species_per_page + 1} to {min((page + 1) * species_per_page, num_species)})"
        ax.set_title(title, y=1.02, fontsize=12)
        fig.subplots_adjust(top=0.93)

        ax.set_yticks(list(range(num_empty_bars, len(species_subset) + num_empty_bars)))
        ax.set_yticklabels([f"{data.loc[data['current_sci_name '] == species, 'listed_sci_name'].iloc[0]}" for species in species_subset], fontsize=6)
        
        ax.set_xlim(min_date, max_date)  # Set x-axis limits to the earliest and latest date_effective
        plt.xticks(rotation=45)

        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

print("Plotting complete.")
