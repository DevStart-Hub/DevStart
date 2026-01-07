import pandas as pd
import numpy as np
import tables
from datetime import datetime
import os

def convert_csv_to_hdf5(input_csv, output_hdf5, 
                        framerate=300, 
                        screen_size=[1920, 1080],
                        raw_format=False):
    """
    Convert CSV to DeToX HDF5 format with correct Root-Level metadata.
    """
    # Read original CSV
    df = pd.read_csv(input_csv, index_col=0)
    
    # =========================================================================
    # 1. Create GAZE DataFrame (with events merged inline)
    # =========================================================================
    gaze_df = pd.DataFrame()
    
    # TimeStamp normalization (start from 0)
    gaze_df['TimeStamp'] = (df['time'] - df['time'].iloc[0]).astype('float64')
    
    # Coordinates: Convert bottom-left origin to center origin (PsychoPy 'pix')
    gaze_df['Left_X'] = (df['L_X'] - screen_size[0]/2).astype('float64')
    gaze_df['Left_Y'] = (df['L_Y'] - screen_size[1]/2).astype('float64')
    gaze_df['Left_Validity'] = df['L_V'].map({True: 1, False: 0}).astype('int64')
    gaze_df['Left_Pupil'] = df['L_P'].astype('float64')
    gaze_df['Left_Pupil_Validity'] = gaze_df['Left_Validity']
    
    gaze_df['Right_X'] = (df['R_X'] - screen_size[0]/2).astype('float64')
    gaze_df['Right_Y'] = (df['R_Y'] - screen_size[1]/2).astype('float64')
    gaze_df['Right_Validity'] = df['R_V'].map({True: 1, False: 0}).astype('int64')
    gaze_df['Right_Pupil'] = df['R_P'].astype('float64')
    gaze_df['Right_Pupil_Validity'] = gaze_df['Right_Validity']
    
    # Events (inline)
    gaze_df['Events'] = df['Event'].fillna('').astype('string')
    
    # =========================================================================
    # 2. Create EVENTS DataFrame (separate table)
    # =========================================================================
    events_mask = df['Event'].notna()
    events_df = pd.DataFrame()
    events_df['TimeStamp'] = (df.loc[events_mask, 'time'] - df['time'].iloc[0]).astype('float64')
    events_df['Events'] = df.loc[events_mask, 'Event'].astype('string')
    
    # Add small negative noise to events so they appear just before the sample
    np.random.seed(42)
    noise = np.random.uniform(-5, 0, size=len(events_df))
    events_df['TimeStamp'] = events_df['TimeStamp'] + noise
    
    # =========================================================================
    # 3. Save to HDF5 with DeToX Metadata
    # =========================================================================
    
    # Prepare fixed-width strings for HDF5
    gaze_df['Events'] = gaze_df['Events'].astype('S50')
    events_df['Events'] = events_df['Events'].astype('S50')
    
    # Convert to structured arrays
    gaze_array = gaze_df.to_records(index=False)
    events_array = events_df.to_records(index=False)
    
    # Open File
    with tables.open_file(output_hdf5, mode='w') as f:
        
        # --- A. Root-Level Metadata (This is what DeToX checks) ---
        # Session Info
        f.root._v_attrs.filename = os.path.basename(output_hdf5)
        f.root._v_attrs.collection_date = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        
        # Hardware Info (As requested)
        f.root._v_attrs.eyetracker_model = "Tobii TX300"
        f.root._v_attrs.eyetracker_serial = "0000000000"  # Repeated 0s
        f.root._v_attrs.illumination_mode = "default"
        
        # Recording Settings
        f.root._v_attrs.framerate = int(framerate)
        
        # Display Configuration
        # Note: DeToX expects these as strings representation of python objects
        f.root._v_attrs.screen_size = str(list(screen_size))
        f.root._v_attrs.window_units = "pix"  
        
        # Data Format Settings
        f.root._v_attrs.raw_format = str(raw_format)
        f.root._v_attrs.coordinate_units = "pix"
        f.root._v_attrs.relative_timestamps = str(True) # Normalized time
        
        # --- B. Create Data Tables ---
        f.create_table(f.root, 'gaze', obj=gaze_array, title='Gaze data samples')
        f.create_table(f.root, 'events', obj=events_array, title='Event markers')

    print(f"✓ Converted {input_csv} -> {output_hdf5}")
    print(f"  Metadata applied: Tobii TX300, 300Hz, pix units, relative time")

    return gaze_df, events_df

if __name__ == "__main__":
    # Example usage
    input_file = r"C:\Users\tomma\OneDrive - Birkbeck, University of London\Personal\Workshop\Chieti_2025\WorkshopFiles\data\Child1.csv"
    output_file = r"C:\Users\tomma\OneDrive - Birkbeck, University of London\Personal\Workshop\Chieti_2025\WorkshopFiles\data\Child1.h5"
    
    convert_csv_to_hdf5(input_file, output_file)