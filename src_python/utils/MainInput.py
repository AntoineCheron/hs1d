# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 13:28:02 2017

@author: Quentin Courtois
"""
import lxml.etree as ET


class MainInput(object):
    """
        Class managing XML files to read and construct inputs for Boussinesq
        Simulmation.
    """
    
    def __init__(self):
        
        #Write needed input XML
        self.ask_input()
        
    def read_xml(self):
        
    
        
    def ask_input(self):
        """
            Method writing a XML file that contains all needed inputs for
            the simulation
        """
        #Creating the XML object
        root = ET.Element("root")
        
        #######################################################################
        #Writting Morpho
        
        ctg = "morphologic"
        names = ["nx", "x", "discretization_type", "x_custom", "angle", "z_custom", "w"]
        types = ["Integer", "Float", "String", "Float", "Float", "Float", "Float"]
        signed = ["True", "True", "", "False", "False", "True", "False" ]
        min_value = ["10", "", "", "-1", "0", "-10", "1"]
        max_value = ["1000", "", "", "", "0.1", "10000", ""]
        description = [ "Number of cells to define over the Hillslope if x_custom is empty", \
                        "Coordinates along the hillslope", \
                        "Type of discretization along the Hillslope", \
                        "A vector of coordinates along the hillslope or a value : -1 if not known and computed based on xmin, xmax and temp",\
                        "A vector of slopes along the Hillslope for each x coordinate or a single value corresponding to a constant slope" , \
                        "Altitude of the top of the hillslope at each x coordinate or -1 if not known", \
                        "Width of the top of the hillslope at each x coordinate"]
        range_val = ["", "", "linear, logarithmic, square, custom", "", "", "", ""]
        adjust = ["False", "False", "False", "False", "False","False", "False"]
        calib  = ["False", "False", "False", "False", "False","False"," False"]
        default = ["100", "0:100", "linear", "-1", "0.05", "-1", "500"]
        
        self.write_category(root, ctg, names, types, signed, min_value, max_value, range_val, description, adjust, calib, default)
        
        #######################################################################
        #Writting Geol
        ctg = "geologic"
        names = ["k", "f", "soil_depth"]
        types = ["Float", "Float", "Float"]
        signed = ["False", "False", "False"]
        min_value = ["10e-8","10e-6","1"]
        max_value = ["10e-1","0.3","100"]
        description = ["Hydraulic conductivity along the hillslope", \
                       "Effective Porosity along the hillslope", \
                       "A vector of soil's depth along the hillslope or a value if constant"]
        range_val = ["", "", ""]
        adjust = ["False", "False", "True"]
        calib = ["True", "True", "False"]
        default = ["1", "0.3", "2"]
        
        self.write_category(root, ctg, names, types, signed, min_value, max_value, range_val, description, adjust, calib, default)
        
        #######################################################################
        #Writting Hydro
        ctg = "hydrlogic"
        names = ["t", "Nt", "unit", "recharge_values", "time_custom", "recharge_period", "recharge_type", "perc_loaded"]
        types = ["Float", "Integer", "String", "Float", "Float", "Float", "String", "Float"]
        signed = ["True", "False", "", "False", "True", "True", "", "False"]
        min_value = ["", "10",  "", "0", "", "", "", "0"]
        max_value = ["", "1000", "", "2", "", "", "", "1"]
        range_val = ["", "", "days, hour, min, sec, year", "", "", "", "periodical, square, steady, random or databased", ""]
        description = ["time value of the recharge", \
                       "Number of time values to describe the recharge", \
                       "Time unit of the time serie", \
                       "A vector containing recharge values in m/unit or a single value if recomputed", \
                       "Time vector corresponding to the recharge chronicle or -1 if computed from tmin, tmax, Nt", \
                       "Period used to compute a periodical recharge or -1 if recharge databased, 0 if random, inf if steady. Possible values : -1, 0, 'inf' or a -inf to inf", \
                       "Type of used recharged in the simulation", \
                       "Percentage of the Hillslope Stock initially filled with water"]
        adjust = ["False", "False", "False", "False", "False", "False", "False", "True"]
        calib = ["False", "False", "False", "False", "False","False","False","False"]
        default = ["0:35", "8400", "days", "30", "-1", "None", "squared", "0"]
        self.write_category(root, ctg, names, types, signed, min_value, max_value, range_val, description, adjust, calib, default)
        
        #######################################################################
        #Writting XML file
        tree = ET.ElementTree(root)
        tree.write("inputs.xml")        
            
    def write_category(self, root, ctg, names, types, signed, min_value, max_value, range_val, description, adjust, calib, default):
        
        #Creating category sub-element
        category = ET.SubElement(root,'category')
        category.set("name",ctg)
        
        for i in range(len(names)):
            #Create input category
            inp = ET.SubElement(category, "input")
            #Name field
            name = ET.SubElement(inp, "name")
            name.text = names[i]
            #Type field
            type_in = ET.SubElement(inp, "type")
            type_in.text = types[i]
            #Signed field
            if signed[i] != "" and types[i] != "String":
                sig = ET.SubElement(inp, "signed")
                sig.text = signed[i]
            #Min value
            if min_value[i] != "":
                min_val = ET.SubElement(inp, "min-value")
                min_val.text = min_value[i]
            #Max value
            if max_value[i] != "":
                max_val = ET.SubElement(inp, "max-value")
                max_val.text = max_value[i]
            #Range value
            if range_val[i] != "":
                ran = ET.SubElement(inp, "value-range")
                ran.text = range_val[i]
            #Description
            if description[i] != "":
                des = ET.SubElement(inp, "description")
                des.text = description[i]
            #Adjustable
            if adjust[i].lower() != "false" and adjust[i] != "":
                ad = ET.SubElement(inp, "adjustable")
                ad.text = adjust[i]
            #Calibration
            if calib[i].lower() != "false" and calib[i] != "":
                cal = ET.SubElement(inp, "ToCalibrate")
                cal.text = calib[i]
            #Default value
            if default[i] != "":
                defa = ET.SubElement(inp, "Default-value")
                defa.text = default[i]
        