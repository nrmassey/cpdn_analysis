###############################################################################
# File         : cpdn_analysis_test.py
# Author       : Neil Massey
# Created      : 18/05/12
# Purpose      : test suite for cpdn_analysis
# Changes      : 
###############################################################################

from cpdn_box import *
from cpdn_means import *
from cpdn_region import *
from cpdn_ensemble import *
from cpdn_regrid import *
from cpdn_remap_box import *
from cpdn_smooth import *
from map_plot import map_plot
from os.path import expanduser as EU
import matplotlib.cm as cm

###############################################################################

### tests ###

if __name__ == "__main__":
    test_number = 14

    if test_number == 0:
        path = EU("/Volumes/Macintosh HD2/shared/sas_test2/")
        file = "nhqhma.pcg7oct.nc"
        box = cpdn_box()
        box.load(path+file, "field16")
        print box.get_dimension("Y").get_values()

#       # file with no bounds data so will guess bounds
#       path = "/home/orwell/cpdn/massey/HadAM3P_Output/1960/hadam3p_m001_1960_2_006101835_0/"
#       file = "m001ma.pag0dec.nc"
#       box = cpdn_box()
#       box.load(path+file, "field16")
#       x = box["1960-12-01":"1961-01-02", 0, :, :]
#       print x.get_values()

    elif test_number == 1:

        # file with bounds data so will read them in
        path = "/home/ares/mad/rye/data/cmip5/CMIP5/output/MOHC/HadGEM2-ES/historical/mon/atmos/Amon/r1i1p1/v20110329/tas/"
        file = "tas_Amon_HadGEM2-ES_historical_r1i1p1_185912-188411.nc"
        box = cpdn_box()
        box.load(path+file, "tas")

        x = box["1860-12-10"]
        d = datetime(1865, 10, 10)
        e = datetime(1866, 10, 10)
        x = box[d:e:2, 90.0, 12.0]
        x = box[30.0]
        x = box[12.2]

        d = datetime(1998, 12, 10, 23, 12, 00)
        x = box[d]
        print x.get_values()
    
    elif test_number == 2:

        # get dimensions test
        path = "/home/orwell/cpdn/massey/ClimateData/HadISST_NV/"
        file = "HadISST_sst.nc"
        box = cpdn_box()
        box.load(path+file, "sst")
        for d in range(0, len(box.get_dimensions())):
            box.get_dimension(d).get_name()
        print box.get_dimension("T").get_values()

    elif test_number == 3:

        # region test
        path = "/home/ares/mad/rye/data/cmip5/CMIP5/output/MOHC/HadGEM2-ES/historical/mon/atmos/Amon/r1i1p1/v20110329/tas/"
        file = "tas_Amon_HadGEM2-ES_historical_r1i1p1_185912-188411.nc"
        
        box = cpdn_box()
        box.load(path+file, "tas")

        r_box = get_region_from_box(box, "EAf")     
        print r_box.get_dimension("X").get_values(), r_box.get_dimension("Y").get_values()

    elif test_number == 4:

        # meaning test
        path = "/home/ares/mad/rye/data/cmip5/CMIP5/output/MOHC/HadGEM2-ES/historical/mon/atmos/Amon/r1i1p1/v20110329/tas/"
        file = "tas_Amon_HadGEM2-ES_historical_r1i1p1_185912-188411.nc"
        var = "tas"
        
        box = cpdn_box()
        box.load(path+file, var)
        g_means = cpdn_global_mean(box)
        g_means.save("g_mean_box2")
#       z_means = cpdn_zonal_mean(box)
#       z_means.save("z_mean_box2")
#       m_means = cpdn_meridional_mean(box)
#       m_means.save("m_mean_box2")
#       t_means = cpdn_temporal_mean(box)
#       t_means.save("t_mean_box2")
#       m_means = cpdn_climatological_month_mean(box)
#       m_means.save("m_mean_box")
#       s_means = cpdn_climatological_seasonal_mean(box)
#       s_means.save("s_mean_box")
        
    elif test_number == 5:
        # ensemble test
        path = "/home/orwell/cpdn/massey/HadAM3P_Output/1960/"
        ens = cpdn_ensemble("field16")
        ens.add_file(path + "*.pb*", limit=5*12, recursive_level=1)
        dt1 = daytime(1960, 12, 1)
        dt2 = daytime(1961, 12, 15)
        mems = ens.subset_by_date(slice(dt1, dt2, None))
        mems2 = ens.subset_by_month("Jan")
        for m in mems2.get_members():
            d = m.get_t_dim().get_values("daytime")
            print m.get_filename(), m[dt1:dt2].get_values().shape, d[0].isoformat(), d[-1].isoformat()
        dt3 = "1960-12-01"
        mems = ens.subset_by_date(dt3)
        print
        for m in mems:
            print m.get_filename(), m[dt1:dt2].get_values().shape
        print
        mems = ens.get_members_month("Jan")
        for m in mems:
            print m.get_filename()
        print
        mems = ens.get_members_year(1961)
        for m in mems:
            print m.get_filename()
        print "DJF"
        mems = ens.get_members_season("DJF")
        for m in mems:
            print m.get_filename()

        ens.save("test_ensemble.pik")

    elif test_number == 6:
        # load ensemble test
        ens2 = cpdn_ensemble()
        ens2.load("test_ensemble.pik")

    elif test_number == 7:
        # regrid test
#       path = "/home/orwell/cpdn/massey/HadAM3P_Output/1960/"
#       file = "ens_mean_field16_1960.nc"
        path = "./"
        file = "delta_2.nc"
        box = cpdn_box()
        box.load(path+file, "tos_sub_tos")
        flip_box = cpdn_remap_latitude(box)
        target_X = numpy.arange(0, 360, 1.875)
        target_Y = numpy.arange(90, -91.25, -1.25)
        tgt_box = cpdn_regrid(flip_box, target_X, target_Y)
        tgt_box.save("delta_3.nc")

#       target_X = numpy.arange(0, 360, 0.9375)
#       target_Y = numpy.arange(90, -90.625, -0.625)
#       tgt_box2 = cpdn_regrid(box, target_X, target_Y, "cubic")
#       tgt_box2.save("regrid_upscale_test.nc")

    elif test_number == 8:
        # plotting test
        path = "/home/orwell/cpdn/massey/HadAM3P_Output/1960/"
        file = "ens_mean_field16_1960.nc"
        box = cpdn_box()
        box.load(path+file, "field16")
        mp = map_plot.map_plot(projection="robinson")
        mp.block_plot(box[0,0])
        mp.draw_colorbar("Temperature in $^\circ$K")
        mp.draw_continents()
        mp.save("block_plot_test.png")

    elif test_number == 9:
        # ensemble mean test
        path = "/home/orwell/cpdn/massey/HadAM3P_Output/1960/"
        file = "ens_mean_field16_1960.nc"
        box1 = cpdn_box()
        box1.load(path+file, "field16")

        box2 = cpdn_box()
        box2.load(path+file, "field16")

        box3 = box1 + box2
        box3.save("test_box3.nc")
        box4 = box1 - box2
        box4.save("test_box4.nc")
        box5 = box1 / box2
        box5.save("test_box5.nc")
        box6 = box1 * box2
        box6.save("test_box6.nc")

    elif test_number == 10:
        path = "/home/orwell/cpdn/massey/HadAM3P_Output/1960/"
        file = "ens_mean_field16_1960.nc"
        box1 = cpdn_box()
        box1.load(path+file, "field16_ens_mean")

        new_box = cpdn_reshape(box1, [25, 145, 192])
        new_box.save("reshape1.nc")
        new_box = cpdn_reshape(box1, [25, 1, 1,145, 192])
        new_box.save("reshape2.nc")

    elif test_number == 11:
#       path = "/home/orwell/cpdn/massey/HadAM3P_Output/1960/"
#       file = "ens_mean_field16_1960.nc"
        path = "./"
        file = "delta_07.nc"
        box1 = cpdn_box()
        box1.load(path+file, "tos")
        box2 = cpdn_spatial_smooth(box1[0], 2, 1, "gauss")
        box2.save("gauss_box.nc")
        box3 = cpdn_spatial_smooth(box1[0], 2, 1, "flat")
        box3.save("flat_box.nc")

    elif test_number == 12:
        path = EU("~ssparrow/cpdn/gen_natural_ssts/cmip5/")
        file = "tos_nat1.nc"
        box = cpdn_box()
        box.load(path+file, "tos")

    elif test_number == 13:
        path = EU("/home/ares/mad/ssparrow/cmip5/NCC/NorESM1-M/")
        file = "tos_Omon_NorESM1-M_historicalNat_r1i1p1_200601-201212.nc"
        box = cpdn_box()
        box.load(path+file, "tos")
        x = box.get_dimension("T").get_values()
        
    elif test_number == 14:     # rotated grid data
        path = EU("/Volumes/MacintoshHD2/shared/sas_test2/")
        file = "nhqhga.pdg7oct.nc"
        box = cpdn_box()
        box.load(path+file, "field16")
        mp = map_plot.map_plot(projection="cylindrical")
        mp.block_plot(box[0,0], colmap=cm.RdYlBu_r)
        mp.draw_continents()
        mp.draw_colorbar(ptype='integer')
        mp.save("rotate_grid.png")
#        x.save("test_rotated.nc")

    elif test_number == 15:     # not rotated sas grid data
        path = EU("/Volumes/MacintoshHD2/shared/sas_test2/")
        file = "nhqhma.pcg7oct.nc"
        box = cpdn_box()
        box.load(path+file, "field16")
        x = box[0]
        print x.get_dimension_lengths()
        x = box[0:1, 0]
        print x.get_dimension_lengths()
        x = box[0,0,90.0:0.0, 45.0:90.0]        
        print x.get_dimension_lengths()
#        mp = map_plot.map_plot(projection="robinson")
#        mp.block_plot(box[0,0])
#        mp.draw_colorbar("Temperature in $^\circ$K")
#        mp.draw_continents()
#        mp.save("block_plot_test.png")