#!/usr/bin/env python
# coding: utf-8

# # AdhE initial screen assay
# Based on C:\Users\d89659k\Documents\Research\Ctherm CBP project\adhE mutant characterization - Angel Pech 2022\AdhE assay 7-6-2023 PROSS ALDH NADPH\AdhE init screen OT2 protocol 7-6-2023 OT2.ipynb

# # Protocol overview
# * Turn on plate reader and get protocol ready
# * Place labware onto OT2 robot (remember this protocol requires 2 boxes of p300 tips)
# * Load protocol on OT2 and calibrate if needed
# * Prepare stock solutions
# * Prepare buffer
# * Manually prepare 12-well reservoir plate
# * Add purified AdhE enzymes to 1st column of deepwell plate
# 
# # OT2 overview
# * Prepare AdhE dilution series
# * Make NADH and NADPH standard curves
# * Add buffer to wells
# * Heat plate to 50C so that when room temp. buffer is added, final temp is close to 37C
# * Add enzyme to wells to start the reaction
# * Manually apply sealing film and load the plate into the platereader

# Remember to include imports in protocol file for the robot
import logging
from opentrons import protocol_api, simulate, types
from opentrons.types import Location, Point # for making point offsets
import opentrons.execute
import pandas as pd # for getting data from Excel file
import itertools
import os
import io # for Excel sharing violation workaround
import math
#from IPython.display import display # this import is only necessary when converting to a standalone python protocol

# Set up logging
ot_logger = logging.getLogger('opentrons')
ot_logger.setLevel(logging.WARNING) # for opentrons, ignore all but warnings

FORMAT = "[%(name)s:%(funcName)s():%(lineno)s]%(levelname)s: %(message)s"
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # for this Jupyter Notebook

#logging.debug('test')

# Check to make sure I'm running the right version of the Opentrons library
logger.info(f'Opentrons version is: {opentrons.__version__}')

metadata = {'apiLevel': '2.12',
            'author'  : 'Dan Olson'}

# For displaying well status
pd.options.display.max_columns = 25 # to enable display of a complete 384 well plate with 24 columns

# Set up constants
DEFAULT_WELL_BOTTOM_CLEARANCE = 1.0 # mm above well bottom
ABOVE_LIQUID_HEIGHT = 8.0 # mm above well bottom
ASSAY_PLATE_WELL_VOLUME = 60 # ul volume of assay plate
ASSAY_PLATE_TOUCH_TIP_Z_OFFSET = -3.3 # mm below top of well. This should be 2 mm above the liquid when they're filled with 60 ul liquid
experiment_setup_excel_file_local = 'AdhE assay setup 7-6-2023.xlsx'
experiment_setup_excel_file_remote = r'/var/lib/jupyter/notebooks/' + experiment_setup_excel_file_local

# choose the right path for the file
# Note 7-10-2023
# When the Opentrons app runs this test, it does so from the directory of the local computer
# C:\Users\d89659k\AppData\Local\Programs\Opentrons
if os.path.exists(experiment_setup_excel_file_local):
    experiment_setup_excel_file = experiment_setup_excel_file_local
    logger.info(f'Excel setup file found in local directory: {experiment_setup_excel_file_local}')
elif os.path.exists(experiment_setup_excel_file_remote):
    experiment_setup_excel_file = experiment_setup_excel_file_remote
    logger.info(f'Excel setup file found in remote directory: {experiment_setup_excel_file_remote}')
else:
    #logger.warning(f'Excel setup file not found in either local or remote directories. Current directory: {os.getcwd()}')
    #experiment_setup_excel_file = experiment_setup_excel_file_remote
    raise ValueError(f'Excel setup file not found in either local {experiment_setup_excel_file_local} or remote {experiment_setup_excel_file_remote} directories. Current directory: {os.getcwd()}')




############### PROTOCOL ACTIONS ######################
def run(protocol: protocol_api.ProtocolContext):
    ####################### SET UP CUSTOM FUNCTIONS WITHIN THE RUN BLOCK #######################
    def blink_rail_lights(num_blinks = 3):
        """
        Blink the rail lights to indicate that user interaction is needed
        """
        for i in range(num_blinks):
            protocol.set_rail_lights(False)
            protocol.delay(seconds = 0.5)
            protocol.set_rail_lights(True)
            protocol.delay(seconds = 0.5)
   
    def reset_well_volumes():
        """Reset volumes in all wells to zero, for purposes of well volume tracking.
        
        Args:
            all_volumes : nested dictionary with labware and well_volume dictionaries
            
        Return:
            all_volumes with volume values set to zero.
        """
        result_dict = {}
        for k, val in all_volumes.items():
            name = k
            plate = val[0]
            plate_volumes = dict(zip(plate.wells_by_name(), itertools.repeat(0)))
    
            all_volumes[name][1] = plate_volumes
    
    def display_plate_status():
        """
        Display the current volume of each well in a given plate
        assume that all_volumes is a dictionary in the global namespace
    
        Args:
            all_volumes : nested dictionary with well volumes
        """
        for key, value in all_volumes.items():
        
            plate = value[0]
            well_volumes = value[1]
            assert len(plate.wells()) == len(well_volumes), f'plate and well_volumes have a different number of wells'
    
            print(f'Well volume for {key}')
            plate.rows_by_name().keys()
            plate.columns_by_name().keys()
            df = pd.DataFrame(plate.wells_by_name().keys(), columns = ['well'])
            df['row'] = df['well'].str.extract(r'([A-Z])')
            df['col'] = df['well'].str.extract(r'([0-9]{1,2})').astype(int)
            df['vol'] = df['well'].map(well_volumes)
    
            piv = pd.pivot_table(df, index = 'row', columns = 'col', values = 'vol', aggfunc = 'mean' )
            s = piv.style.format("{:.0f}") # 
            # display(s)
    
    def indices_from_well_name(well_name):
        """
        Given an alphanumeric well name, return the row and column index, starting from 0
        The first character must be a letter
        The next characters must be numbers
        """
        row_idx = ord(well_name[:1].upper())-65 # convert first letter of well ID to a row index, starting from zero
        col_idx = int(well_name[1:]) - 1
        return(row_idx, col_idx)
    
    def well_id_from_name(well_name, rows_in_plate = 16):
        """
        Given an alphanumeric well name, return the well ID
        Parameters:
           well_name: Name of the well, i.e. B15
           rows_in_plate: int, 16 for a 384 well plate
        """
        row, col = indices_from_well_name(well_name)
        
        return(col*16+row)
        
    def well_name_from_indices(row, col):
        """
        Given an row and col index, starting from 0
        Return a well name
        """
        row_letter = chr(row + 65)
        return(row_letter + str(col+1))
    
    def offset_well(well_name, col_offset = 0, row_offset = 0):
        """
        Given a starting well name, offset it by a certain number of rows and columns
        """
        row_num, col_num = indices_from_well_name(well_name)
        return(well_name_from_indices(row_num + row_offset, col_num + col_offset))
               
    def get_384_well_dispense_height(current_volume, volume_to_add):
        """
        Given a well on the 384 well assay_plate, and a volume that is being added,
        determine what height to use for dispensing to keep the pipette tip at the top of the liquid
        """
        
        well_depth = 11.43 # mm, from Opentrons library https://labware.opentrons.com/corning_384_wellplate_112ul_flat
        well_volume = 112 # ul, from Opentrons library https://labware.opentrons.com/corning_384_wellplate_112ul_flat
        EXTRA_Z_OFFSET = 0.8 # mm, extra height to add to keep pipette tip just above the liquid
        volume_conversion_factor = well_depth/well_volume # mm per ul volume
        
        #### need to fix this DO 4-13-2022
        #current_volume = assay_plate_well_volumes[well_id]
        
        if (current_volume + volume_to_add) > 100:
            logger.warning('You are planning to add more than 100 ul to a well, this volume is too large')
        
        height = (current_volume + volume_to_add) * volume_conversion_factor + EXTRA_Z_OFFSET 
        
        if height < 1:
            height = 1.0 # keep the dispense height at least 1 mm above the bottom of the well
        return height  
    
    def get_well_edge_offset(pipette, well, height_above_bottom):
        """
        Given a pipette and a well, provide the offset needed to pipette on the edge of the well
        
        Args:
            pipette : opentrons pipette object
            well : opentrons well object
            height_above_bottom : height above the bottom of the well, in mm 
        """
        EXTRA_OFFSET = 0.2 # additional offset to ensure tip is positively engaged with side of labware
        
        # calculate height below the top of the well, that is the key parameter for determining pipette offset
        height_below_well_top = well.depth - height_above_bottom
        assert height_below_well_top > 0, f"height_above_bottom is too high, pipette tip above well top"
        
        # set the tip width of the pipette, to allow for pipetting onto the sides of wells
        # the first parameter is the tip diameter (mm)
        # the second parameter is the tip slope parameter (mm dia per mm tip length)
        tip_dict = {'p20_single_gen2' : [0.90, 0.083],
                    'p300_multi_gen2' : [1.07, 0.093]
                   }
        try:            
            PIPETTE_TIP_WIDTH = tip_dict[pipette.name][0] + height_below_well_top * tip_dict[pipette.name][1]
        except KeyError as e:
            print(f'No tip width defined for the {pipette.name}')
        
        # calculate well diameter
        if isinstance(well, opentrons.protocol_api.labware.Well): # check if the destination is a well object
            well_diameter = well.width # note, this converts the well to a location object
            logger.debug(f'destination type is well: {well}')
        elif isinstance(well, opentrons.types.Location): # check if the destination is a location object
            well_diameter = well.labware.as_well().width
            logger.debug(f'destination type is location: {well}')
        else:
            raise RuntimeError(f'Destination type {type(well)} not allowed, must be a Well or Location object')
        
        offset = (well_diameter - PIPETTE_TIP_WIDTH)/2 + EXTRA_OFFSET
        return offset
    
    def p20_std_curve(component, buffer, start_wells, dilutions = 14, blanks = 2):
        """Set up a standard curve using the p20 single channel pipette on a 384 well plate called assay_plate. 
        
        Each dilution is 2-fold
        The first well is a 1:3 dilution of the component well concentration
        The default volume is 20 ul
        The final volume is 60 ul
        
        Args:
          component : a point object for the source of the starting component
          buffer : a point object for the source of the buffer component
          start_wells : list of opentrons well objects
          dilutions : number of dilutions in the dilution series (each dilution in a separate row)
          blanks : number of blank wells at the end of the dilution series
        """
        protocol.home() # home the robot at the beginning of the standard curve to improve XY accuracy
        assert dilutions > 1, f"A dilution series must have more than 1 dilution"
            
        
        for start_well in start_wells:
            protocol.comment(f'  ### SETTING UP STANDARD CURVE in well {start_well.well_name} ###')
            logger.info(f'  ### SETTING UP STANDARD CURVE in well {start_well.well_name} ###')
            (stdcrv_rowidx, stdcrv_colidx) = indices_from_well_name(start_well.well_name)
            logger.debug(f'Starting column:{stdcrv_colidx}, Starting row:{stdcrv_rowidx}')
    
            # define the wells to add buffer to (all dilutions except the first, and all blanks)
            stdcrv_wells = assay_plate.columns()[stdcrv_colidx][stdcrv_rowidx:stdcrv_rowidx + dilutions + blanks] 
    
            # add buffer to all wells, except first ones
            transfer_volume = 20
            dispense_height = 2
            source = buffer
            destination = [well for well in stdcrv_wells[1:]] # all stdcrv_wells except the first
            logger.debug(f"Buffer wells:{[well.well_name for well in stdcrv_wells[1:]]}")
            custom_transfer(pipette = p20, 
                            transfer_volume = transfer_volume, 
                            source = source, 
                            destination_list = destination, 
                            dispense_height = dispense_height,
                            offset_direction = 'left',
                            do_blowout = False,
                            drop_tip = True
                           ) # perform the transfer
    
            # add component to first well
            # define the wells
            protocol.comment(f'      Adding component from {component} to first well of assay plate')
            logger.info(f'      Adding component from {component} to first well of assay plate')
            transfer_volume = 20
            source = component
            dispense_height = 2
            destination = [start_well]
            custom_transfer(pipette = p20, 
                            transfer_volume = transfer_volume, 
                            source = source, 
                            destination_list = destination, 
                            dispense_height = dispense_height,
                            offset_direction = 'left',
                            do_blowout = False,
                            drop_tip = False
                           ) # perform the transfer
    
            # add component to second well (which already has buffer in it)
            dispense_height = 4
            destination = [assay_plate.rows()[stdcrv_rowidx+1][stdcrv_colidx]]
            custom_transfer(pipette = p20, 
                    transfer_volume = transfer_volume, 
                    source = source, 
                    destination_list = destination, 
                    dispense_height = dispense_height,
                    offset_direction = 'left', 
                    do_blowout = False,
                    drop_tip = False,
                    prewet_tip = False        
                   ) # perform the transfer
            
            # make serial dilution of component (usually NADH or NADPH)
            protocol.comment('      Making serial dilution of component')
            logger.info('      Making serial dilution of component')
            serial_dilution(pipette = p20, transfer_volume = transfer_volume, mix_volume = transfer_volume, 
                            wells = stdcrv_wells[1:dilutions], mix_steps = 5, dispense_height = 4, 
                            blowout_height = 4, mix_before = True, do_blowout = False)
            
        # add buffer to any columns that are part of the standard curve
        buffer_wells = [] # empty list to hold wells
        for row in reversed(range(2)): # for a 384 well plate, we only address rows A and B
            for start_well in start_wells:
                assert start_well.well_name[0] == 'A', 'The starting well must be in Row A'
                (stdcrv_rowidx, stdcrv_colidx) = indices_from_well_name(start_well.well_name)
                buffer_wells.append(assay_plate.columns()[stdcrv_colidx][row])
    
        # convert buffer wells to locations
        buffer_well_locations = [well for well in buffer_wells] 
        logger.debug(f"Buffer wells:{[well.well_name for well in buffer_wells]}")
        
        # perform transfer
        # custom_transfer deals with volume tracking automatically
        custom_transfer(pipette = p300m,
                        transfer_volume = 40,
                        source = buffer,
                        destination_list = buffer_well_locations,
                        dispense_height = 6,
                        offset_direction = 'right',
                        prewet_tip = True,
                        do_touch_tip = False,
                        do_blowout = False,
                        transfer_type = 'distribute',
                        drop_tip = True,
                        evenly_split_transfers = True, 
                       )
        
    def serial_dilution(pipette, transfer_volume, mix_volume, wells, mix_steps = 5, dispense_height = 1, 
                        blowout_height = 10, mix_before = False, do_blowout = True, always_get_new_tip = False):
        """Perform serial dilution
        
        Args:
            pipette : opentrons pipette object
            transfer volume : float
            mix_volume : float
            wells : list of opentrons well objects
                aspirate from the first well, serial transfer into remaining wells
            mix_steps : int
                number of times to mix
            dispense_height : float
                set value > 1 to allow for axial mixing
            blowout_height : float
            do_blowout : boolean, false disables blowout steps
            always_get_new_tip: boolean, true forces tip discard step between each dilution
        """
        # pick up a tip if the pipette doesn't have one
        if not pipette.has_tip:
            logger.info(f"Pick up new pipette tip")
            pipette.pick_up_tip()
            
        if mix_before:
            # mix
            logger.info(f'Mixing before serial transfer')
            for j in range(mix_steps):
                logger.debug(f'...mix step {j}')
                pipette.aspirate(mix_volume, wells[0], rate = 1.0)
                pipette.dispense(mix_volume, wells[0].bottom(z = dispense_height), rate = 10.0)
                protocol.delay(seconds = 0.5) # brief delay for better mixing
            if (do_blowout):
                logger.info(f"Blowing out pipette")
                pipette.blow_out(wells[0].bottom(z = blowout_height))
            pipette.touch_tip(speed = 20.0) # touch tip to knock off droplets, use minimum speed
            
        
        for i in range(len(wells)-1):
            logger.info(f'Serial transfer from well {wells[i]} to well {wells[i+1]}')
            pipette.aspirate(transfer_volume, wells[i], rate = 1.0)
            pipette.dispense(transfer_volume, wells[i+1], rate = 10.0)
            adjust_well_volume(transfer_volume, wells[i], wells[i+1], pipette.channels)
    
            # mix
            for j in range(mix_steps):
                logger.debug(f'...mix step {j}')
                pipette.aspirate(mix_volume, wells[i+1], rate = 10.0)
                pipette.dispense(mix_volume, wells[i+1].bottom(z = dispense_height), rate = 10.0)
                protocol.delay(seconds = 0.5) # brief delay for better mixing
            if (do_blowout):
                logger.info(f"Blowing out pipette")
                pipette.blow_out(wells[i+1].bottom(z = blowout_height))
            pipette.touch_tip(speed = 20.0) # touch tip to knock off droplets, use minimum speed
            if (always_get_new_tip):
                pipette.drop_tip()
                logger.info(f"Pick up new pipette tip")
                pipette.pick_up_tip()
    
        # discard extra volume from the final wells
        logger.debug(f"Discarding extra volume from well {wells[i+1]}")
        pipette.aspirate(transfer_volume, wells[i+1])
        adjust_well_volume(transfer_volume, wells[i+1], trash, pipette.channels)
        pipette.drop_tip() # drop tip in the trash
         
    def custom_transfer(pipette, 
                        transfer_volume, 
                        source, 
                        destination_list, 
                        source_blowout_height = 25, 
                        dest_blowout_height = 8,
                        dispense_height = 1,
                        offset_direction = 'center', 
                        prewet_tip = True,
                        drop_tip = True, 
                        do_blowout = True,
                        do_touch_tip = True,
                        transfer_type = 'single',
                        dispense_rate = 1.0, # aspirate and dispense rate,
                        inspect_tips = False, 
                        evenly_split_transfers = False, # when True, this will split a set of 18 dispenses into 2 groups of 9 
                       ):
        """
        Custom transfer routine for the p20 or p300m pipettes for transferring liquid from one source
        well to one or more destination wells via a series of individual transfers.
        
        Args:
          pipette: Opentrons pipette object
          transfer_volume: volume to be transferred in ul
          source: source location object or well object
          destination_list: list of destination location objects or well objects
          source_blowout_height: source_blowout height in mm above bottom of well
          dest_blowout_height: destination blowout height in mm above bottom of well
          dispense_height: dispense height above the bottom of the well, in mm
          offset_direction: location within well to pipette. Can be center, left, right, top, or bottom.
          prewet_tip: boolean, true causes pipette to aspirate and dispense once to prewet the pipette tip
          drop_tip: boolean, true causes pipette tip to be dropped at the end of the transfer
          do_blowout: boolean, false suppresses all blowout steps
          do_touch_tip: boolean, false suppresses all touch_tip steps
          transfer_type:
              'single' -- one aspirate for each dispense step
              'distribute' -- aspirate a large volume and then distribute to multiple wells, refilling pipette as needed
        """       
                    
        # check to make sure the pipetting volume is in a range that is allowed by the pipette
        assert (transfer_volume <= pipette.max_volume) & (transfer_volume >= pipette.min_volume), 'Transfer volume out of range' 
                    
        # get point location if source object is a well
        if isinstance(source, opentrons.protocol_api.labware.Well):
            source = source.bottom(z=1) # return well bottom location
            logger.debug(f'Transferring {transfer_volume} ul from {source.labware}')
        
        # pick up a tip if the pipette doesn't have one
        if not pipette.has_tip:
            pipette.pick_up_tip()
        
        # prewet tip
        if (prewet_tip):
            pipette.aspirate(pipette.max_volume, source)
            pipette.dispense(pipette.max_volume, source)
            if (do_blowout):
                pipette.blow_out(source.labware.as_well().bottom(z = source_blowout_height)) # blow out above aspirate location
            pipette.touch_tip()
        
        if transfer_type == 'single':
            for well in destination_list:
                # perform transfer
                logger.debug(f'Performing SINGLE transfer from well {source.labware} to well {well}')
                logger.debug(f'Pipette {pipette} current volume is {pipette.current_volume}')
                pipette.aspirate(transfer_volume, source)
                pipette.touch_tip()
                offset_dispense(pipette, 
                                well, 
                                transfer_volume, 
                                dispense_height = dispense_height,
                                offset_direction = offset_direction, 
                                do_blowout = do_blowout, 
                                do_touch_tip = do_touch_tip,
                                dispense_rate = dispense_rate,
                                inspect_tips = inspect_tips,
                               )
    
                # adjust well volumes
                adjust_well_volume(transfer_volume, source, well, pipette.channels)
    
            if drop_tip:
                logger.debug('Dropping tips')
                pipette.drop_tip()
                
        elif transfer_type == 'distribute':
            logger.debug('############# PERFORMING DISTRIBUTE TRANSFER ################')
            working_destination_list = destination_list.copy() # list of wells that still need to be distributed into
            
            # calculate the maximum number of wells that can be filled with one aspirate
            carryover = 0.1 # fraction of pipette volume to use as carryover
            max_wells_per_aspirate = int((pipette.max_volume*(1-carryover))/transfer_volume)
            logger.info(f'max_wells_per_aspirate = {max_wells_per_aspirate}')
            
            if (evenly_split_transfers == True):
            # re-calculate the maximum number of wells to fill with one aspirate to ensure that the aspirate volume is split evenly
                # calculate the number of wells per aspirate if the aspirations are divided into equal sized chunks
                total_wells = len(working_destination_list) # total number of wells that will be distributed into
                num_transfer_groups = math.ceil(total_wells/max_wells_per_aspirate)
                max_wells_per_aspirate = math.ceil(total_wells/num_transfer_groups)
                logger.info(f'EVENLY SPLIT TRANSFERS, max_wells_per_aspirate = {max_wells_per_aspirate}')
            
            # while there are still destination wells left to pipette into
            while len(working_destination_list) > 0:
            
                # move wells from working_destination_list to temp_transfer_list
                temp_transfer_list = []
                for i in range(max_wells_per_aspirate):
                    if len(working_destination_list) > 0:
                        temp_transfer_list.append(working_destination_list.pop())
                
                logger.debug(f'Distributing to {len(temp_transfer_list)} of {len(destination_list)} wells')
    
                # calculate volume to aspirate
                aspirate_volume = (transfer_volume * len(temp_transfer_list)) + (pipette.max_volume * carryover)
                pipette.aspirate(aspirate_volume - pipette.current_volume, source) # account for volume already in pipette
                logger.info(f'Aspirating {aspirate_volume:.1f} ul from {source}')
                if do_touch_tip:
                    pipette.touch_tip()
                
                for well in temp_transfer_list:
                    # perform transfer
                    logger.debug(f'Performing DISTRIBUTE transfer from well {source} to well {well}')
                    offset_dispense(pipette, 
                                    well, 
                                    transfer_volume, 
                                    dispense_height = dispense_height,
                                    offset_direction = offset_direction, 
                                    do_blowout = False, # never do a blowout step for a custom distribute 
                                    do_touch_tip = do_touch_tip,
                                    dispense_rate = dispense_rate,
                                    inspect_tips = inspect_tips,
                                   )
                    logger.debug(f'Pipette {pipette} current volume is {pipette.current_volume:.2f} after dispense')
    
                    # adjust well volumes
                    # note that this doesn't account for the disposal volume correctly
                    adjust_well_volume(transfer_volume, source, well, pipette.channels)
    
            if drop_tip:
                logger.debug('Dropping tips')
                pipette.drop_tip()
       
    def offset_dispense(pipette, 
                        well, 
                        transfer_volume,
                        offset_direction = 'center', 
                        custom_offset_distance = None,
                        dispense_height = 1,
                        do_blowout = False,
                        dest_blowout_height = 8, # height from bottom of well
                        do_touch_tip = False,
                        dispense_rate = 1.0, # dispense rate, as a multiple of the default dispense rate 
                        inspect_tips = False, # move tips up and pause after dispense to allow inspection of tips for troubleshooting
                       ):
        """
        Dispense liquid into a well at a position offset from the center. Pipette should already
        have aspirated liquid.
        """
    
        assert isinstance(well, opentrons.protocol_api.labware.Well), f'Well {well} must be a well object'
      
        # calculate offset for pipetting against edge of well
        if custom_offset_distance:
            offset_distance = custom_offset_distance # override the automated well offset calculaton
        else:
            offset_distance = get_well_edge_offset(pipette = pipette, well = well, height_above_bottom = dispense_height)
            logger.debug(f'offset is: {offset_distance:.2f}')
    
        # calculate the x and y offsets for a given direction
        if offset_direction == 'center':
            x_offset = 0
            y_offset = 0
        elif offset_direction == 'right':
            x_offset = offset_distance
            y_offset = 0
        elif offset_direction == 'left':
            x_offset = -offset_distance
            y_offset = 0
        elif offset_direction == 'top':
            x_offset = 0
            y_offset = offset_distance
        elif offset_direction == 'bottom':
            x_offset = 0
            y_offset = -offset_distance
        else:
            raise RuntimeError('Unknown well destination location, must be "center", "left", "right", "top", or "bottom"')    
        
        logger.debug(f'Performing offset dispense into well {well}, offset {offset_direction}, distance {offset_distance:.2f}')
        
        # perform dispense operation
        if (inspect_tips):
            logger.info(f'Moving tips up to allow inspection for debugging')
            pipette.move_to(well.top(z=10)) # move to 10 mm above the well top
            protocol.delay(seconds=3) # pause 3 seconds to visually inspect pipette tips
        
        pipette.move_to(well.bottom(z=dispense_height)) # move to the center of the well
        pipette.move_to(well.bottom(z=dispense_height).move(Point(x_offset, y_offset, 0)), force_direct = True) # move over to the side of the well
        logger.info(f'Dispense height is {dispense_height:.2f}')
        pipette.dispense(transfer_volume, rate = dispense_rate) # no transfer location specified, dispense in current location
        pipette.move_to(well.bottom(z=dispense_height), force_direct = True) # move back to the center of the well to avoid flicking drops of liquid from the edge of the well
        if (do_blowout):
            logger.debug('Performing blowout')                  
            # do blowout at the dispense location offset by the blowout height
            # first, move to the center of the well at the specified blowout height below the top of the well
            pipette.move_to(well.bottom(z = dest_blowout_height))
            # then move over to the edge, based on the offset_direction parameter
            pipette.move_to(well.bottom(z=dispense_height).move(Point(x = x_offset,y = y_offset)), force_direct = True)
            pipette.blow_out() # blow out at current location
            pipette.move_to(well.bottom(z = dest_blowout_height), force_direct = True) # move pipette back to center
        if (do_touch_tip):
            pipette.touch_tip(speed = 20) # 1/3 of default speed of 60 mm/min      
               
    def adjust_well_volume(transfer_volume, source, destination, num_channels = 1):
        """Keep track of the volume removed from the source and added to the destination
        
        Assume that all_volumes is available as a global variable
        
        Parameters:
            transfer_volume: (float) volume of liquid transferred, in ul
            source:          (Opentrons well object) source well
            destination:     (Opentrons well object) destination well
            num_channels:    (int) number of pipette channels, must be 1 or 8
        """
        # get well if source object is a point location
        if isinstance(source, opentrons.types.Location):
            source = source.labware.as_well() # return well object
            
        # get well if source object is a point location
        if isinstance(destination, opentrons.types.Location):
            destination = destination.labware.as_well() # return well object
            
        # collect source well information
        source_plate = all_volumes[str(source.parent)][0] # opentrons labware object
        source_plate_rows = len(source_plate.rows()) # integer number of rows in the plate
        source_start_row = source.well_name[0]
        
        # collect destination well information
        not_trash = True
        if 'Fixed Trash' in str(destination.parent):
            not_trash = False
            destination_plate_rows = 0
        else:
            destination_plate = all_volumes[str(destination.parent)][0] # opentrons labware object
            destination_plate_rows = len(destination_plate.rows()) # integer number of rows in the plate
            destination_start_row = destination.well_name[0]
        
        # determine the number of pipette channels withdrawing liquid from the source well
        if num_channels == 1: # source and destination are both single wells
            uptake_channels = 1
            dispense_channels = 1
            source_well_list = [source]
            destination_well_list = [destination]
            
        elif num_channels == 8:
            if source_plate_rows == 1: # one column per channel
                uptake_channels = 8 # all 8 channels drawing from the same well
                source_well_list = [source]
            if source_plate_rows == 8:
                assert source_start_row == 'A', f'Multiwell transfers from an {source_plate_rows} row plate must start at row A, not row {source_start_row}'
                uptake_channels = 1
                source_well_list = source_plate.columns_by_name()[source.well_name[1:]] # all of the wells in the column
            if source_plate_rows == 16:
                assert source_start_row in ['A', 'B'], f'Multiwell transfers from a {source_plate_rows} row plate must start at row A or B, not row {source_row}'
                uptake_channels = 1
                start_dict = {'A':0, 'B':1}
                start_idx = start_dict[source_start_row]
                source_well_list = source_plate.columns_by_name()[source.well_name[1:]][start_idx::2] # every other wells in the column         
                
            if destination_plate_rows == 1:
                dispense_channels = 8 # all 8 channels dispensing into the same well
                destination_well_list = [destination]
            if destination_plate_rows == 8:
                assert destination_start_row == 'A', f'Multiwell transfers to an {destination_plate_rows} row plate must start at row A, not row {destination_start_row}'
                dispense_channels = 1 
                destination_well_list = destination_plate.columns_by_name()[destination.well_name[1:]] # all of the wells in the column
            if destination_plate_rows == 16:
                assert destination_start_row in ['A', 'B'], f'Multiwell transfers to a {destination_plate_rows} plate must start at row A or B, not row {destination_row}'
                dispense_channels = 1 
                start_dict = {'A':0, 'B':1}
                start_idx = start_dict[destination_start_row]
                destination_well_list = destination_plate.columns_by_name()[destination.well_name[1:]][start_idx::2] # every other wells in the column   
                    
        else:
            raise ValueError('Wrong number of channels. num_channels should be 1 or 8')
        
        # subtract volume from source well(s)
        for s in source_well_list:
            logger.debug(f'Subtracting {transfer_volume * uptake_channels} ul from {s}')
            all_volumes[str(source.parent)][1][s.well_name] -= transfer_volume * uptake_channels
            
        # add volume to destination well(s)
        if not_trash:
            for d in destination_well_list:
                logger.debug(f'Adding {transfer_volume * dispense_channels} ul to {d}')
                all_volumes[str(destination.parent)][1][d.well_name] += transfer_volume * dispense_channels
        

    ## set up the machine
    protocol.set_rail_lights(True) # turn the lights on so we can see what's going on
    protocol.comment('  ### SETTING UP THE INSTRUMENT ###') # adds a user-readable comment string
    temp_module = protocol.load_module('temperature module gen2', '8') # temperature module loaded into position 7 on the deck
    assay_plate = temp_module.load_labware('corning_384_wellplate_112ul_flat') # 384 well plated loaded onto the temperature module
    deepwell_plate = protocol.load_labware('usascientific_96_wellplate_2.4ml_deep', '5')
    reservoir_plate = protocol.load_labware('nest_12_reservoir_15ml', '2') # 12 well reservoir loaded onto position 5
    tip300_1 = protocol.load_labware('opentrons_96_tiprack_300ul', '6')
    tip300_2 = protocol.load_labware('opentrons_96_tiprack_300ul', '3')
    tip20_1 = protocol.load_labware('opentrons_96_tiprack_20ul', '9')
    p300m = protocol.load_instrument('p300_multi_gen2', 'right', tip_racks=[tip300_1, tip300_2])
    p20 = protocol.load_instrument('p20_single_gen2', 'left', tip_racks=[tip20_1])
    
    ## Set up default liquids and locations
    # individual named wells
    trash = protocol.fixed_trash['A1']
    
    # ## Keep track of the volume in each well
    # * The main purpose of this is to allow for adjusting the pipette height to pipette at the top of the liquid volume to reduce cross-contamination
    # * It is also useful to make sure the initial well volume is high enough for the expected pipetting
    # * For each plate, there is a dictionary with well names as keys and well volume as values
    ## keep track of well volumes
    assay_plate_volumes = dict(zip(assay_plate.wells_by_name(), itertools.repeat(0))) # dictionary to hold well volumes
    deepwell_plate_volumes = dict(zip(deepwell_plate.wells_by_name(), itertools.repeat(0))) # dictionary to hold well volumes
    reservoir_plate_volumes = dict(zip(reservoir_plate.wells_by_name(), itertools.repeat(0))) # dictionary to hold well volumes
    ## make dictionary to hold volume information for all plates
    all_volumes = {'NEST 12 Well Reservoir 15 mL on 2': [reservoir_plate, reservoir_plate_volumes],
                   'USA Scientific 96 Deep Well Plate 2.4 mL on 5': [deepwell_plate, deepwell_plate_volumes],
                   'Corning 384 Well Plate 112 µL Flat on Temperature Module GEN2 on 8': [assay_plate, assay_plate_volumes]
                  }
    
    
    ############ PROCESS EXCEL FILE ################
    ### Process 'Assay plate wells' sheet to generatre well_info dataframe ###
    # import the level information from the Excel file
    # eventually I want to end up with a dataframe that has one row for each component being added to each well
    
    # to get around Excel sharing violations, read the excel file to an object and manually close it
    # see more details at github https://github.com/pandas-dev/pandas/issues/29803
    with open(experiment_setup_excel_file, "rb") as f:
        file_io_obj = io.BytesIO(f.read())
    
    well_raw = pd.read_excel(file_io_obj, 
                       sheet_name = 'Assay plate wells', 
                       skiprows = 2,
                       header = [0,1],
                       #index_col = 1,
                       engine = 'openpyxl'
                      )
    
    well_raw.dropna(subset = [('Metadata', 'Well')], inplace = True) # get rid of rows with no well ID (i.e. subtotal rows on the bottom)
    
    # make a dataframe for the metadata
    df_met = well_raw['Metadata'].copy()
    df_met['Col'] = df_met['Col'].astype('int') # make sure column values are stored as integers
    df_met.Well = df_met.Well.astype('int')
    
    # Turn the Compound, Level, and Volume columns into a long-form dataframe
    # Choose the columns of interest
    wr2 = well_raw.loc[:, ['Compound', 'Level', 'Volume']]
    wr2 = wr2.droplevel(level=0, axis=1).drop('Total', axis=1)
    wr2.columns = pd.MultiIndex.from_tuples(wr2.columns.str.split('_').tolist())
    
    # Convert to long-form using stack()
    wr3 = wr2.stack(0)
    wr3 = wr3.reset_index(level=1)
    
    #wr3 = wr3.droplevel(1)
    wr3.columns = ['Category', 'Level', 'Compound', 'Volume']
    wr3['Level'].fillna(1, inplace = True) # fill na values from the water column
    
    # combine the dataframes
    well_info = df_met.merge(wr3, left_index = True, right_index = True) # merge metadata with well data
    
    # fix level column
    well_info.Level = well_info.Level.astype('int') # set to integer, since it will be used for indexing later
    
    # get rid of rows we don't need
    well_info = well_info[~(well_info.Volume == 0)] # get rid of rows where the transfer volume is 0
    well_info.dropna(subset = ['Volume'], inplace = True)
    
    # add a well name column
    well_info.insert(3, 'Name', well_info.Row + well_info.Col.astype(str))
    
    ### Process 'Parameters to test' sheet to make source_well dataframe ###
    # import the level information from the Excel file
    source_well = pd.read_excel(file_io_obj, 
                       sheet_name = 'Parameters to test', 
                       skiprows = 8,
                       header = [0],
                       #index_col = 0,
                       engine = 'openpyxl'
                      )
    
    # find the first row with no value in the "Component" column
    last_row = source_well.Component[source_well.Component.isnull()].index[0]
    source_well = source_well.iloc[0:last_row].dropna(axis=1, how = 'all') 
    
    ### Add source well and pipette information to the well_info dataframe ###
    # Set source wells, this takes a few seconds to run
    source_well_list = [] # empty list to hold wells
    pipette_list = [] # empty list to hold pipette info
    
    # loop through the well_info dataframe and add the source well and pipette information
    for index, row in well_info.iterrows():
        component_ser = source_well.loc[source_well.Component == row.Compound, :]
        plate_str = component_ser['Plate'].values[0] # source plate string
        well_str = component_ser['Starting well'].values[0] # source well level 1
        offset_well_str = offset_well(well_str,row.Level-1,0) # offset the well based on the desired dilution level
        
        # add the source well to the list
        source_well_list.append(eval(plate_str).wells_by_name()[offset_well_str]) # the Opentrons Well object for the desired well 
        
        # add the pipette to the list
        pipette_list.append(eval(component_ser['Pipette'].values[0]))
        
    well_info['Source'] = source_well_list # add source info to dataframe
    well_info['Pipette'] = pipette_list # add pipette info to dataframe  
    # Note that the pipette assignments are not correct for the Std curves. This doesn't matter because the pipetting protocol ignores the "Source" and "Pipette" columns of the well_info dataframe.

    
    ########### PREPARE CFE DILUTIONS ################
    ## Reset tips and well volumes ##
    if p300m.has_tip: p300m.drop_tip() # drop tips if possible
    tip300_1.reset()
    tip300_2.reset()
    tip20_1.reset()
    reset_well_volumes()
    ot_logger.setLevel('WARNING') # for troubleshooting this step
    logger.setLevel('INFO')
    
    number_serial_dilutions = 8
    protocol.home()
    transfer_volume = 300
    aspirates_per_well = 0
    max_aspirates_per_well = 10000/(transfer_volume * number_serial_dilutions)
    source_wells = iter(reservoir_plate.wells()[10:12])
    serial_dilution_wells = deepwell_plate.rows_by_name()['A'][0:number_serial_dilutions] # all wells in serial dilution series
    source = next(source_wells) # starting well
    
    # add buffer to wells of deepwell plate, skip first well
    for destination_well in serial_dilution_wells[1:]:
        if aspirates_per_well >= max_aspirates_per_well-1: # if volume is too low, go to the next well
            source = next(source_wells)
            aspirates_per_well = 0
        logger.info(f'Adding buffer from well {source} to well {destination_well}')
        custom_transfer(p300m, 
                        transfer_volume, 
                        source, 
                        [destination_well], 
                        drop_tip = False, 
                        prewet_tip = False,
                        do_blowout = False,
                        evenly_split_transfers = True
                       )
        aspirates_per_well += 1
    
    # perform serial dilution
    serial_dilution(pipette = p300m, 
                    transfer_volume = 300, 
                    mix_volume = 200, 
                    wells = serial_dilution_wells, 
                    mix_steps = 10, 
                    dispense_height = 12, 
                    blowout_height = 20, # this doesn't matter, since do_blowout is False 
                    mix_before = True,
                    do_blowout = False,
                    always_get_new_tip = True
                   ) 
        
    ########### SET UP STANDARD CURVE #####################
    #### Automatically choose standard curves based on Excel file
    # * Choose wells of type "Std curve"
    # * Standard curves will always start in row A
    # * Standard curves will always be made in some kind of assay buffer, so we can ignore water
    # * The main component in a standard curve will always be added by the p20 pipette
    # * The buffer in the standard curve will always be added by the p300 pipette
    # * Replicates will always be in columns that are adjacent to each other
    do_pipetting = True # set this flag to false to troubleshoot without actually piptetting, 
                         # note that well volumes are not updated if the do_pipetting flag is false
    ## Choose standard curve wells from well_info dataframe
    if do_pipetting: protocol.home()
    # Note: backticks are needed for names with spaces
    # There's no good way to do multi-line queries with comments
    query_string = """
                   `Reaction type` == "Std curve" 
                   """
    std_wells = well_info.query(query_string).copy()
    std_wells[:7]
    
    # assume that all standard curves are performed column-wise, and have the same number of dilutions
    dilutions = std_wells.query('Category == "C1"').groupby('Col')['Well'].count().values[0]
    blanks = well_info.query('`Reaction type` == "Blank"').groupby('Col')['Well'].count().values[0]
    
    # assume that the starting wells of std curves have either "Buffer_NADH" or "Buffer_NADPH" in the standard curves
    for compound in std_wells.query('Category == "C3"').Compound.unique():
        compound_source = std_wells.query('Compound == @compound')['Source'].iloc[0]
        start_well_names = std_wells.query('Compound == @compound')['Name'].values
        start_wells = [assay_plate[i] for i in start_well_names]
        buffer_source = std_wells.query('Name == @start_well_names[0] & Category == "C1"')['Source'].values[0] # buffer associated with first well in start_wells
        
        logger.info(f'Preparing standard curve of component {compound} in wells {start_well_names} with {dilutions} dilutions and {blanks} blanks')
        if do_pipetting:
            p20_std_curve(compound_source, buffer_source, start_wells, dilutions = dilutions, blanks = blanks)
        
        
    ################ ADD REAGENTS WITH P300 MULTICHANNEL PIPETTE #############
    ### Prepare wells for pipetting
    #logger.setLevel('DEBUG')
    reaction_types_to_ignore = ['Std curve', 'Blank']
    assay_start_components = source_well.query('`Start component` == "yes"').Component.tolist() # wait until the end to pipette these
    
    # destination wells grouped by compound, for all subsequent pipetting steps
    grouped = (well_info
        .query('`Reaction type` not in @reaction_types_to_ignore')
        .loc[:, ['Compound', 'Source', 'Level', 'Volume', 'Name', 'Row', 'Col']]
        .groupby('Compound')
    )
      
    ### Pause and wait for the user to be ready to start the assays
    blink_rail_lights(3)
    protocol.pause(msg = "The next step will start adding reagets with the p300 multichannel pipette. Press the Run button to resume.")
    
    protocol.home()
    do_pipetting = True # set this flag to false to troubleshoot without actually piptetting
    # make a list of components to be added, sorted by the "order of addition" column
    component_df = (source_well
                    .query('Pipette == "p300m" & Component not in @assay_start_components')
                    .loc[:, ['Component', 'Order of addition']]
                    .sort_values('Order of addition')
                   )
    p300m_component_list = component_df['Component'].tolist()
    
    for component in p300m_component_list:
        logger.info(f'Loading component {component}')
        if component in grouped.groups.keys(): # make sure the component is actually used 
            level_groups = grouped.get_group(component).groupby('Level')
            levels = sorted(level_groups.groups.keys(), reverse = True) # make a sorted list of levels, from highest to lowest
    
            for lev in levels:
                ab_rows = ['A', 'B']
                one_level = level_groups.get_group(lev).query('Row in @ab_rows').sort_values('Col')
                if len(one_level) == 0:
                        # no level groups in this level, skip to the next level
                        logger.debug(f'   Level -- {lev} is empty')
                        continue
                
                logger.info(f'   Level -- {lev}')
                p300m.pick_up_tip()
    
                assert len(one_level) > 0, 'Error, no volume groups in this level.'
                assert len(one_level.Source.unique()) == 1, 'Error, multiple source locations detected for one component and level'
                source = one_level.Source.unique()[0] # this should be unique for a given component and level
                
                # thoroughly mix source well prior to pipetting
                # 1-31-2023, don't bother mixing. wells should be pre-mixed to save time and ensure uniformity
                # 2-2-2024, re-enable this to pre-wet tips
                # 2-12-2024, I disabled this to prevent bubble formation. I don't think it helped accuracy
                #logger.info(f'Mixing source well {source}')
                #p300m.mix(8, 200, source)
                
                # group by destination volumes
                vol_groups = one_level.groupby('Volume')
                vols = sorted(vol_groups.groups.keys(), reverse = True) # pipette larger volumes first, since they're usually control wells
                logger.info(f'      Volume groups: {vols}')
                for v in vols:
                    one_vol = vol_groups.get_group(v)
                    transfer_volume = v
                    destination_locations = [assay_plate.wells_by_name()[well] for well in one_vol.Name]
                    
                    # calculate dispense height based on well volume
                    plate_key = list(all_volumes.keys())[2] # assume that the destination is in the assay_plate, which is the 3rd position in the all_volumes dictionary
                    first_well_name = destination_locations[0].well_name 
                    first_well_volume = all_volumes[plate_key][1][first_well_name] 
                    dispense_height = get_384_well_dispense_height(first_well_volume, 0) # 2-7-2024 changed transfer volume to 0 to dispense at the bottom of the liquid volume
                    logger.info(f'Well {first_well_name} has vol. {first_well_volume:.2f}, adding {transfer_volume}, dispense height: {dispense_height:.2f}')
                    
                    # perform the pipetting
                    if do_pipetting:
                        custom_transfer(pipette = p300m, 
                                        transfer_volume = transfer_volume, 
                                        source = source, 
                                        destination_list = destination_locations,
                                        dispense_height = dispense_height,
                                        offset_direction = 'left',
                                        prewet_tip = False,  
                                        drop_tip = False,
                                        do_blowout = False,
                                        do_touch_tip = False,
                                        transfer_type = 'distribute',
                                        dispense_rate = 0.5, # 2-2-2024, dispense rate reduced from 1 to 0.5 to improve accuracy
                                        inspect_tips = False, # 2-9-2024, for troubleshooting
                                        evenly_split_transfers = True,
                                       )
                if p300m.has_tip:
                    logger.debug('Dropping tips')
                    p300m.drop_tip() # drop tip between levels, not between different volumes of the same level
    
        
    ########## ALERT THE USER TO START THE ASSAY ################
    blink_rail_lights(3) # blink the rail lights 3 times
    protocol.pause("The next step will heat the plate and add the start reagent. Press the Run button to resume.")
       
    ########### START REACTION ################
    temp_module.set_temperature(50) # set temperature, wait for block to reach temperature, then proceed with the next step
    
    
    protocol.home()
    do_pipetting = True # set this flag to false to troubleshoot without actually piptetting
    # make a list of components to be added, sorted by the "order of addition" column
    component_df = (source_well
                    .query('Pipette == "p300m" & Component in @assay_start_components')
                    .loc[:, ['Component', 'Order of addition']]
                    .sort_values('Order of addition')
                   )
    p300m_component_list = component_df['Component'].tolist()
    
    for component in p300m_component_list:
        logger.info(f'Loading component {component}')
        if component in grouped.groups.keys(): # make sure the component is actually used 
            level_groups = grouped.get_group(component).groupby('Level')
            levels = sorted(level_groups.groups.keys(), reverse = True) # make a sorted list of levels, from highest to lowest
    
            for lev in levels:
                ab_rows = ['A', 'B']
                one_level = level_groups.get_group(lev).query('Row in @ab_rows').sort_values('Col')
                if len(one_level) == 0:
                        # no level groups in this level, skip to the next level
                        logger.debug(f'   Level -- {lev} is empty')
                        continue
                
                logger.info(f'   Level -- {lev}')
                p300m.pick_up_tip()
    
                assert len(one_level) > 0, 'Error, no volume groups in this level.'
                assert len(one_level.Source.unique()) == 1, 'Error, multiple source locations detected for one component and level'
                source = one_level.Source.unique()[0] # this should be unique for a given component and level
                              
                # group by destination volumes
                vol_groups = one_level.groupby('Volume')
                vols = sorted(vol_groups.groups.keys(), reverse = True) # pipette larger volumes first, since they're usually control wells
                logger.info(f'      Volume groups: {vols}')
                for v in vols:
                    one_vol = vol_groups.get_group(v)
                    transfer_volume = v
                    destination_locations = [assay_plate.wells_by_name()[well] for well in one_vol.Name]
                    
                    # calculate dispense height based on well volume
                    plate_key = list(all_volumes.keys())[2] # assume that the destination is in the assay_plate, which is the 3rd position in the all_volumes dictionary
                    first_well_name = destination_locations[0].well_name 
                    first_well_volume = all_volumes[plate_key][1][first_well_name] 
                    dispense_height = get_384_well_dispense_height(first_well_volume, 0) # 2-14-2024 dispense below top of liquid height to reduce bubble formation
                    logger.info(f'Well {first_well_name} has vol. {first_well_volume:.2f}, adding {transfer_volume}, dispense height: {dispense_height:.2f}')
                    
                    # perform the pipetting
                    if do_pipetting:
                        custom_transfer(pipette = p300m, 
                                        transfer_volume = transfer_volume, 
                                        source = source, 
                                        destination_list = destination_locations,
                                        dispense_height = dispense_height,
                                        offset_direction = 'left',
                                        prewet_tip = False, 
                                        drop_tip = False,
                                        do_blowout = False,
                                        do_touch_tip = False,
                                        transfer_type = 'distribute',
                                        dispense_rate = 0.5, # 2-2-2024, dispense rate reduced from 1 to 0.5 to improve accuracy
                                        evenly_split_transfers = True,
                                       )
                if p300m.has_tip:
                    logger.debug('Dropping tips')
                    p300m.drop_tip() # drop tip between levels, not between different volumes of the same level

    blink_rail_lights(1) # indicate that the plate is ready to be removed
    protocol.comment('Remove plate, apply sealing film, centrifuge, load into plate reader.')
        
    # Turn off the temperature module at the end
    temp_module.deactivate()   
    protocol.home()
    protocol.set_rail_lights(False) # turn off the rail lights when the protocol is done
    protocol.comment('The protocol has completed')
    
        
    
    
