# Chap2and3_standalone.py
# PURPOSE : To generate least cost paths as polylines and identify major flow paths
#           in each assesment unit
# Inputs  : BF_Combined_0325.shp
#           AU_Combined_0325.shp
#           SL_Combined_0325.shp
#           CAN_DEM.tif
# OUTPUTS : CostPath_assess_*.shp
# LAST EDIT: 14 Aug 2021 by Sara

import arcpy
import time
import os
import shutil
import csv
import numpy
import pandas as pd
from arcpy import env
from arcpy.sa import *

def flow_paths(inBF, inAU, inSL, inRaster):

    # Work directory
    inPath = os.path.dirname(os.path.realpath(inBF))
    newPath = inPath.replace(os.sep, '/')
    gdbPath = newPath + "/results/flowpaths.gdb"
    mergePath = newPath + "/results/flowpathsMerge.gdb"
    # Environment Settings
    env.workspace = inPath
    env.overwriteOutput=True
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference("North America Albers Equal Area Conic")
    arcpy.env.geographicTransformations = "WGS_1984_To_NAD_1983"
    arcpy.CheckOutExtension("Spatial")
    arcpy.CheckOutExtension("3D")

    # Part 1: Select data needed for the analysis and export to workspace, create basin DEM rasters
    # SELECT LAYER BY ATTRIBUTE AND COPY FEATURES (basin)
    inBFfl = arcpy.management.MakeFeatureLayer(inBF, "inBFfl")
    selectOutput_BF = gdbPath + "/BF_" + GRID # output basin feature class
    arcpy.management.SelectLayerByAttribute(inBFfl, "NEW_SELECTION", "\"gridcode\" = {}".format(GRID))
    arcpy.CopyFeatures_management(inBFfl, selectOutput_BF)

    # SELECT LAYER BY ATTRIBUTE AND COPY FEATURES (assessment unit)
    inAUfl = arcpy.management.MakeFeatureLayer(inAU, "inAUfl")
    selectOutput_AU = gdbPath + "/AU_" + GRID
    arcpy.management.SelectLayerByAttribute(inAUfl, "NEW_SELECTION", "\"basin_id\" = {}".format(GRID))
    arcpy.CopyFeatures_management(inAUfl, selectOutput_AU)

    # SELECT LAYER BY ATTRIBUTE AND COPY FEATURES (stream link)
    inSLfl = arcpy.management.MakeFeatureLayer(inSL, "inSLfl")
    selectOutput_SL = gdbPath + "/SL_" + GRID
    arcpy.management.SelectLayerByAttribute(inSLfl, "NEW_SELECTION", "\"basin_id\" = {}".format(GRID))
    arcpy.CopyFeatures_management(inSLfl, selectOutput_SL)

    # Clear Selection (ensure all selections are cleared)
    arcpy.SelectLayerByAttribute_management(inBFfl, "CLEAR_SELECTION", "\"GRID\" = \"\"")
    arcpy.SelectLayerByAttribute_management(inAUfl, "CLEAR_SELECTION", "\"basin_id\" = \"\"")
    arcpy.SelectLayerByAttribute_management(inSLfl, "CLEAR_SELECTION", "\"basin_id\" = \"\"")

    # ADD FIELD and CALCULATE FIELD (assessment units)
    fieldName_assess = "AU_char"
    fieldType_assess = "TEXT"
    fieldValue = "!assess_id!"
    arcpy.AddField_management(selectOutput_AU, fieldName_assess, fieldType_assess, "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
    arcpy.CalculateField_management(selectOutput_AU, fieldName_assess, fieldValue, "PYTHON_9.3")

    # SPLIT (assessment units)
    outWorkspace = gdbPath
    arcpy.Split_analysis(selectOutput_AU, selectOutput_AU, fieldName_assess, outWorkspace, "")
    arcpy.management.Delete(selectOutput_AU)

    # ADD GEOMETRY ATTRIBUTES: Add extent to basin attribute table to reduce processing time.
    basinExtent = selectOutput_BF
    arcpy.AddGeometryAttributes_management(basinExtent, "EXTENT", "", "", "")
    print("Basin extent: " + basinExtent)

    # SEARCH CURSOR: X Min
    Xmin = "EXT_MIN_X"
    num = [row[0] for row in arcpy.da.SearchCursor(basinExtent, Xmin)]
    XminExtent = str(num)[1:-1]

    # SEARCH CURSOR: Y Min
    Xmax = "EXT_MAX_X"
    num = [row[0] for row in arcpy.da.SearchCursor(basinExtent, Xmax)]
    XmaxExtent = str(num)[1:-1]

    # SEARCH CURSOR: X Max
    Ymin = "EXT_MIN_Y"
    num = [row[0] for row in arcpy.da.SearchCursor(basinExtent, Ymin)]
    YminExtent = str(num)[1:-1]

    # SEARCH CURSOR: Y Max
    Ymax = "EXT_MAX_Y"
    num = [row[0] for row in arcpy.da.SearchCursor(basinExtent, Ymax)]
    YmaxExtent = str(num)[1:-1]

    env.workspace = gdbPath

    # CLIP: Clip raster to extent of basin +/- 5 meter.
    clipRaster = gdbPath + "/ClipDEM_" + GRID
    Xmin = float(XminExtent) - 5
    Ymin = float(YminExtent) - 5
    Xmax = float(XmaxExtent) + 5
    Ymax = float(YmaxExtent) + 5
    rectangle = str(Xmin) + " " + str(Ymin) + " " + str(Xmax) + " " + str(Ymax)
    arcpy.management.Clip(inRaster, rectangle, clipRaster, basinExtent, "-9999", "#", "MAINTAIN_EXTENT")

    # Part 2: Create slope rasters, reclassify them, and create destination point layers for each basin
    arcpy.env.mask = selectOutput_BF + ".shp"

    # gather basin ID associated with each input basin
    gridcodeField = "gridcode"
    gridcodeNum = [row[0] for row in arcpy.da.SearchCursor(selectOutput_BF, gridcodeField)]
    gridcode = str(gridcodeNum)[1:-1]

    # SLOPE
    basinDEM = clipRaster
    outSlope = gdbPath + "/Slope" + GRID
    arcpy.gp.Slope_sa(basinDEM, outSlope, "DEGREE", "1")

    # Reclass Variables: Equal Interval classification scheme
    minSlopeValue = arcpy.GetRasterProperties_management(outSlope, "MINIMUM")
    minSlope = float(minSlopeValue.getOutput(0))
    maxSlopeValue = arcpy.GetRasterProperties_management(outSlope, "MAXIMUM")
    maxSlope = float(maxSlopeValue.getOutput(0))
    classSize = (maxSlope - minSlope) / 10
    c1 = minSlope + classSize # upper class limit of first class
    c2 = c1 + classSize
    c3 = c2 + classSize
    c4 = c3 + classSize
    c5 = c4 + classSize
    c6 = c5 + classSize
    c7 = c6 + classSize
    c8 = c7 + classSize
    c9 = c8 + classSize
    c10 = c9 + classSize # upper class limit of 10th class
    print('Upper class break values:')
    print(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10) # new upper class limits for each slope raster
    reclass = RemapValue([[minSlope, c1, 1],[c1, c2, 2],[c2, c3, 3],[c3, c4, 4],[c4, c5, 5],[c5, c6, 6],[c6, c7, 7],[c7, c8, 8],[c8, c9, 9], [c9, c10, 10]])

    # RECLASSIFY
    outReclass = gdbPath + "/Reclass_" + GRID
    reclassField = "Value"
    slopeReclass = Reclassify(outSlope, reclassField, reclass)
    slopeReclass.save(outReclass)

    # RASTER TO POINT
    inRaster = clipRaster
    outSourcePoint = gdbPath + "/SourcePoint_" + GRID
    arcpy.RasterToPoint_conversion(inRaster, outSourcePoint, "Value")

    # Part 3: Create least cost path polylines for each assessment unit, identify major flow paths
    env.workspace = gdbPath

    # Local Variables
    auName =  "T*"
    assessUnits = arcpy.ListFeatureClasses(auName)

    for assessUnit in assessUnits:
        arcpy.env.mask = assessUnit

        # gather assess ID associated with each input assessment unit
        assessField = "assess_id"
        assessNum = [row[0] for row in arcpy.da.SearchCursor(assessUnit, assessField)]
        assessID = str(assessNum)[1:-1] # used to identify inputs for Path Distance tool

        # gather basin ID associated with each input assessment unit
        basinField = "basin_id"
        basinNum = [row[0] for row in arcpy.da.SearchCursor(assessUnit, basinField)]
        basinID = str(basinNum)[1:-1] # used to identify inputs for Path Distance tool

        # PATH DISTANCE
        inStreamlink = selectOutput_SL
        inSlope = gdbPath + "/Reclass_" + GRID
        inDEM = clipRaster
        outPathDist = resultsPath + "/PathDist_" + assessID + ".tif"
        outBackLink = resultsPath + "/BackLink_" + assessID + ".tif"
        arcpy.gp.PathDistance_sa(inStreamlink, outPathDist, inSlope, inDEM, "", "TABLE ", inDEM, "LINEAR 2 -90 90 -0.022222", "", outBackLink, "", "", "", "", "")

        # COST PATH AS POLYLINE
        outCostPath = resultsPath + "/CostPath_" + assessID
        # sourcePoint = gdbPath + "/SourcePoint_" + GRID
        arcpy.gp.CostPathAsPolyline_sa(outSourcePoint, outPathDist, outBackLink, outCostPath, "EACH_CELL", "pointid")

        # ADD SURFACE INFORMATION: calculates surface distance of flow paths
        inCostPath = outCostPath + ".shp"
        arcpy.AddSurfaceInformation_3d(inCostPath, clipRaster, "SURFACE_LENGTH", "BILINEAR", "", "1", "0", "NO_FILTER")

        # ADD GEOMETRY ATTRIBUTES
        fieldNameX = "END_X_char"
        fieldNameY = "END_Y_char"
        fieldName = "JoinField"
        fieldTypeSTR = "TEXT"
        fieldTypeINT = "LONG"
        fieldValue = "!END_X_char!" + "!END_Y_char!"
        arcpy.AddGeometryAttributes_management(inCostPath, "LINE_START_MID_END", "", "", "")
        arcpy.AddField_management(inCostPath, fieldNameX, fieldTypeSTR, "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(inCostPath, fieldNameY, fieldTypeSTR, "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.AddField_management(inCostPath, fieldName, fieldTypeSTR, "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(inCostPath, fieldNameX, "!END_X!", "PYTHON_9.3")
        arcpy.CalculateField_management(inCostPath, fieldNameY, "!END_Y!", "PYTHON_9.3")
        arcpy.CalculateField_management(inCostPath, fieldName, fieldValue, "PYTHON_9.3")
        # arcpy.CalculateField_management(inCostPath, "DestID", "!FID!", "PYTHON_9.3")

        # DISSOlVE
        outMajorFlowPathGroups = resultsPath + "/MajorFlowPathGroups_assess_" + assessID
        arcpy.Dissolve_management(inCostPath, outMajorFlowPathGroups, fieldName, "", "MULTI_PART", "DISSOLVE_LINES")

        # ADD FIELD
        inMajorFlowPathGroups = outMajorFlowPathGroups + ".shp"
        fieldNameID = "MFPGID" # i.e., MajorFlowPathGroupID
        fieldValueID = "!FID!"
        arcpy.AddField_management(inMajorFlowPathGroups, fieldNameID, fieldTypeINT, "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(inMajorFlowPathGroups, fieldNameID, fieldValueID, "PYTHON_9.3")

        # JOIN FIELD
        inCostPath = outCostPath + ".shp"
        arcpy.JoinField_management(inCostPath, fieldName, inMajorFlowPathGroups, fieldName, fieldNameID)

        # Delete unecessary fields
        deleteFields = ["START_X", "START_Y", "MID_X", "MID_Y", "END_X", "END_Y", "END_X_char", "END_Y_char", "JoinField"]
        arcpy.DeleteField_management(inCostPath, deleteFields)
        arcpy.DeleteField_management(inMajorFlowPathGroups, "JoinField")

        # get number of polylines in each cost path feature class
        # resultCP = arcpy.GetCount_management(inCostPath)
        # countCP = int(resultCP.getOutput(0))
        # print(countCP)

        cpValues = [row[0] for row in arcpy.da.SearchCursor(inCostPath, "DestID")]
        sliceobject = slice(1)
        cpValues[sliceobject]
        minCP = min(cpValues)
        print minCP

        # calculate mean and std dev Slope and Elevation values for each cost path polyline
        # for num in range(1,countCP + 1):
        for cp in cpValues[:10]:
            inCPfl = arcpy.management.MakeFeatureLayer(inCostPath, "inCPfl")
            selectOutput_CP = mergePath + "/CP_" + str(cp)
            arcpy.management.SelectLayerByAttribute(inCPfl, "NEW_SELECTION", "\"DestID\" = {}".format(cp))
            arcpy.CopyFeatures_management(inCPfl, selectOutput_CP)
            outTableSlope = mergePath + "/SlopeZS" + str(cp)
            arcpy.gp.ZonalStatisticsAsTable(selectOutput_CP, "DestID", outSlope, outTableSlope, "DATA", "MEAN_STD")
            arcpy.AlterField_management(outTableSlope, "MEAN", "SLOPE_M")
            arcpy.AlterField_management(outTableSlope, "STD", "SLOPE_SD")
            outTableDEM = mergePath + "/ElevationZS" + str(cp)
            arcpy.gp.ZonalStatisticsAsTable(selectOutput_CP, "DestID", inDEM, outTableDEM, "DATA", "MEAN_STD")
            arcpy.AlterField_management(outTableDEM, "MEAN", "DEM_M")
            arcpy.AlterField_management(outTableDEM, "STD", "DEM_SD")
            arcpy.SelectLayerByAttribute_management(inCPfl, "CLEAR_SELECTION", "\"DestID\" = \"\"")
            arcpy.JoinField_management(outTableSlope, "DestID", outTableDEM, "DestID", ["DEM_M", "DEM_SD"])
            outTableZS = mergePath + "/ZS" + str(cp)
            arcpy.Rename_management(outTableSlope, outTableZS)
            arcpy.management.Delete(outTableDEM)
            arcpy.management.Delete(outTableDEM) # for stubborn files
            arcpy.management.Delete(selectOutput_CP)
            arcpy.management.Delete(selectOutput_CP)

        # Local Variables
        zsName =  "ZS*"
        env.workspace = mergePath
        zonalStats = arcpy.ListTables(zsName)
        minZS = "ZS" + str(minCP)
        zonalStats.remove(minZS) # remove ZS1 from list and use it as the target table
        print(zonalStats)
        arcpy.Append_management(zonalStats, minZS, "NO_TEST") # append all tables to ZS1

        # join ZS1 to cost path shapefile
        zsTableCP = "ZonalStats_CostPath_" + assessID
        arcpy.Rename_management(minZS, zsTableCP)
        arcpy.JoinField_management(inCostPath, "DestID", zsTableCP, "DestID")
        #arcpy.JoinField_management(inCostPath, "DestID", zsTableCP, "DestID", ["SLOPE_M", "SLOPE_SD", "DEM_M", "DEM_SD"])

        # Save to csv
        nparrCP = arcpy.da.FeatureClassToNumPyArray(inCostPath, ["FID", "PathCost", "DestID", "SLength", "MFPGID", "SLOPE_M", "SLOPE_SD", "DEM_M", "DEM_SD"])
        pd.DataFrame(nparrCP).to_csv(resultsPath + '/CostPath_ZonalStats.csv', index=None)

        for zonalStat in zonalStats: # delete all tables except ZS1
            arcpy.management.Delete(zonalStat)
            arcpy.management.Delete(zonalStat)

        # Major Flow Path Groups
        arcpy.JoinField_management(inMajorFlowPathGroups, "MFPGID", inCostPath, "MFPGID", "DestID")
        mfpgValues = [row[0] for row in arcpy.da.SearchCursor(inMajorFlowPathGroups, "DestID")]
        sliceobject = slice(1)
        mfpgValues[sliceobject]
        minMFPG = min(mfpgValues)
        print minMFPG

        # calculate mean and std dev Slope and Elevation values for each major flow path group
        # for num in range(1,countMFPG + 1):
        for mfpg in mfpgValues[:5]:
            inMFPGfl = arcpy.management.MakeFeatureLayer(inMajorFlowPathGroups, "inMFPGfl")
            selectOutput_MFPG = mergePath + "/MFPG_" + str(mfpg)
            arcpy.management.SelectLayerByAttribute(inMFPGfl, "NEW_SELECTION", "\"MFPGID\" = {}".format(mfpg))
            arcpy.CopyFeatures_management(inMFPGfl, selectOutput_MFPG)
            outTableSlope = mergePath + "/SlopeZS" + str(mfpg)
            inSlope = gdbPath + "/Slope" + str(GRID)
            print(inSlope)
            arcpy.gp.ZonalStatisticsAsTable(selectOutput_MFPG, "MFPGID", inSlope, outTableSlope, "DATA", "MEAN_STD")
            arcpy.AlterField_management(outTableSlope, "MEAN", "SLOPE_M")
            arcpy.AlterField_management(outTableSlope, "STD", "SLOPE_SD")
            outTableDEM = mergePath + "/ElevationZS" + str(mfpg)
            arcpy.gp.ZonalStatisticsAsTable(selectOutput_MFPG, "MFPGID", inDEM, outTableDEM, "DATA", "MEAN_STD")
            arcpy.AlterField_management(outTableDEM, "MEAN", "DEM_M")
            arcpy.AlterField_management(outTableDEM, "STD", "DEM_SD")
            arcpy.SelectLayerByAttribute_management(inMFPGfl, "CLEAR_SELECTION", "\"MFPGID\" = \"\"")
            arcpy.JoinField_management(outTableSlope, "MFPGID", outTableDEM, "MFPGID", ["DEM_M", "DEM_SD"])
            outTableZS = mergePath + "/ZS" + str(mfpg)
            arcpy.Rename_management(outTableSlope, outTableZS)
            arcpy.management.Delete(outTableDEM)
            arcpy.management.Delete(outTableDEM) # for stubborn files
            arcpy.management.Delete(selectOutput_MFPG)
            arcpy.management.Delete(selectOutput_MFPG)

        # Local Variables
        zsName =  "ZS*"
        env.workspace = mergePath
        zonalStats = arcpy.ListTables(zsName)
        zonalStats.remove(minMFPG) # remove ZS1 from list and use it as the target table
        print(zonalStats)
        arcpy.Append_management(zonalStats, minMFPG, "NO_TEST") # append all tables to ZS1

        # join ZS1 to cost path shapefile
        zsTableMFPG = "ZonalStats_MFPG_" + assessID
        arcpy.Rename_management(minMFPG, zsTableMFPG)
        arcpy.JoinField_management(inMajorFlowPathGroups, "MFPGID", zsTableMFPG)
        #arcpy.JoinField_management(inMajorFlowPathGroups, "MFPGID", zsTableMFPG, "MFPGID", ["SLOPE_M", "SLOPE_SD", "DEM_M", "DEM_SD"])

        # join average shape length value to major flow path groups attribute table
        cpArray = arcpy.da.FeatureClassToNumPyArray(inCostPath, ["FID", "PathCost", "DestID", "SLength", "MFPGID", "SLOPE_M", "SLOPE_SD", "DEM_M", "DEM_SD"])
        cpDataFrame = pd.DataFrame(cpArray) # create pandas data frame
        shapeLengthMean = cpDataFrame.groupby('MFPGID').mean() # group data frame by MFPGID
        shapeLengthFilter = shapeLengthMean.filter(items=['MFPGID', 'SLength']) # filter data frame to only include MGFPID and mean shape length
        outTable = resultsPath + "/outTable.csv"
        shapeLengthFilter.to_csv(outTable) # export dataframe to csv
        arcpy.TableToTable_conversion(outTable, mergePath, "avgSlength")
        os.remove(outTable)
        arcpy.JoinField_management(inMajorFlowPathGroups, "MFPGID", mergePath + "/avgSlength", "MFPGID", ["SLength"])

        # export table to csv
        nparrMFPG = arcpy.da.FeatureClassToNumPyArray(inMajorFlowPathGroups, ["FID", "MFPGID", "SLOPE_M", "SLOPE_SD", "DEM_M", "DEM_SD", "SLength"])
        pd.DataFrame(nparrMFPG).to_csv(resultsPath + '/MajorFlowPathGroup_ZonalStats.csv', index=None)

        for zonalStat in zonalStats: # delete all tables except ZS1
            arcpy.management.Delete(zonalStat)
            arcpy.management.Delete(zonalStat)

    assessName = "T*" # delete assessment unit shapefiles
    assessUnits = arcpy.ListFeatureClasses(assessName)
    for assessUnit in assessUnits:
        arcpy.management.Delete(assessUnit)
        arcpy.management.Delete(assessUnit)

    arcpy.management.Delete(selectOutput_BF)
    arcpy.management.Delete(selectOutput_BF)

    arcpy.management.Delete(selectOutput_SL)
    arcpy.management.Delete(selectOutput_SL)

    arcpy.management.Delete(outSourcePoint)
    arcpy.management.Delete(outSourcePoint)

    arcpy.management.Delete(clipRaster)
    arcpy.management.Delete(clipRaster)

    arcpy.management.Delete(outSlope)
    arcpy.management.Delete(outSlope)

    arcpy.management.Delete(outReclass)
    arcpy.management.Delete(outReclass)

print("Gathering inputs: {}".format(time.ctime()))

GRIDS = ["1399", "1018"]
inBF = "D:/Sara/Chapter1_Export_210521/BF_Combined_0521.shp"
inAU = "D:/Sara/Chapter1_Export_210521/AU_Combined_0521.shp"
inSL = "D:/Sara/Chapter1_Export_210521/SL_Combined_0521.shp"
inRaster = "D:/Sara/DEM/CAN_DEM.tif"

"""
while True:
    try:
        GRIDcodes = input("Please enter the GRIDS (e.g., [\"1399\", \"1018\"]): ")
        GRIDS = GRIDcodes
    except:
        print("Invalid path")
        continue
    else:
        break
while True:
    try:
        print("Example: \"D:/Sara/Chapter1_Export_210521/BF_Combined_0521.shp\"")
        inBFpath = input("Please enter the path of your basin input: ")
        inBF = inBFpath
    except:
        print("Invalid path")
        continue
    else:
        break
while True:
    try:
        print("Example: \"D:/Sara/Chapter1_Export_210521/AU_Combined_0521.shp\"")
        inAUpath = input("Please enter the path of your assessment unit input: ")
        inAU = inAUpath
    except:
        print("Invalid path")
        continue
    else:
        break
while True:
    try:
        print("Example: \"D:/Sara/Chapter1_Export_210521/SL_Combined_0521.shp\"")
        inSLpath = input("Please enter the path of your stream link input: ")
        inSL = inSLpath
    except:
        print("Invalid path")
        continue
    else:
        break
while True:
    try:
        print("Example: \"D:/Sara/DEM/CAN_DEM.tif\"")
        inRasterPath = input("Please enter the path of your DEM raster: ")
        inRaster = inRasterPath
    except:
        print("Invalid path")
        continue
    else:
        break
"""

inPath = os.path.dirname(os.path.realpath(inBF))
newPath = inPath.replace(os.sep, '/')
resultsPath = newPath + "/results" # "/results"

if os.path.exists(resultsPath): # create results folder
    resultsPath = newPath + "/results"
else:
    #shutil.rmtree(resultsPath)
    os.mkdir(resultsPath)

gdbPath = resultsPath + "/flowpaths.gdb" # create flow path gdb
gdbName = "flowpaths.gdb"

if arcpy.Exists(gdbPath):
    gdbPath = resultsPath + "/flowpaths.gdb"
else:
    #shutil.rmtree(gdbPath)
    arcpy.CreateFileGDB_management(resultsPath, gdbName)

mergePath = resultsPath + "/flowpathsMerge.gdb" # create flow path gdb
mergeName = "flowpathsMerge.gdb"

if arcpy.Exists(mergePath):
    mergePath = resultsPath + "/flowpathsMerge.gdb"
else:
    #shutil.rmtree(mergePath)
    arcpy.CreateFileGDB_management(resultsPath, mergeName)

print("Generating Flow Paths and Major Flow Paths: {}".format(time.ctime()))

for GRID in GRIDS:
    flow_paths(inBF, inAU, inSL, inRaster)


#for GRID in GRIDS:
#    try:
#        flow_paths(inBF, inAU, inSL, inRaster)
#    except arcpy.ExecuteError:
#        print(arcpy.GetMessages())
#    finally:
#        print("Code completed")
