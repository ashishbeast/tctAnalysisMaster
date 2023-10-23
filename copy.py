#!/usr/bin/env python3
import os, shutil

wafers = [1, 9, 12]
areas = {6:'A', 25:'B', 100:'C'}
types = {1:'0', 2:'1', 3:'2'}
fileTypes = {'':'1_VS_IR_a_Noise.tct', '_av_256':'1_VS_IR_a.tct', '_NoGain_av_256':'0_VS_IR_a.tct'}

for w in wafers:
    for a in areas:
        for t in types:
            for ft in fileTypes:
                source = f'w{w}Data/w_{w}_type_{t}_area_{a}{ft}'
                target = f'data/slapp_1MIP_W_{w}_pad_{areas[a]}_type_{types[t]}_gain_{fileTypes[ft]}'
                if os.path.exists(source):
                    try:
                        shutil.copy(source, target)
                    except Exception as e:
                        print(f"Error copying {source} to {target}: {e}")
                        continue
                else:
                    source = f'w{w}Data/w_{w}_type_{t}_area_{a}_Hole_4{ft}'
                    try:
                        shutil.copy(source, target)
                    except FileNotFoundError:
                        source = f'w{w}Data/w_{w}_type_{t}_area_{a}_Hole_1{ft}'
                        if os.path.exists(source):
                            shutil.copy(source, target)
                        else:
                            continue
                    except Exception as e:
                        print(f"Error copying {source} to {target}: {e}")
                        continue
                        
                    
