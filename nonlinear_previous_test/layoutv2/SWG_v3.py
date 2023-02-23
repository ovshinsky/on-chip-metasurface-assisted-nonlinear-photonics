import math

from splayout import *

##define cell
cell = Cell("PBS_v2-2")
##define layer
wg_layer = Layer(1, 0)
grating_layer = Layer(2, 0)
perturabation_layer = Layer(6, 1)
heater_layer = Layer(4, 0)
pad_layer = Layer(5, 0)
inv_layer = Layer(6, 1)

# text_start_point = Point(0,100)
# text = Text(text_start_point,"CELL3")
# text.draw(cell,wg_layer)

##################perturbations: 函数
def one_single_etch(amplitude1, period, offsetx, offsety):
    # points for the polygon

    pointlist = [Point(-period / 4 + offsetx, amplitude1/2 + offsety),
                 Point(period / 4 + offsetx, amplitude1 / 2 + offsety),
                 Point(period / 4 + offsetx, -amplitude1 / 2 + offsety),
                 Point(-period / 4 + offsetx, -amplitude1 / 2 + offsety),
                 # Point(period / 8 + offsetx, -0.707 * amplitude1 / 2 + offsety),
                 # Point(0 + offsetx, -amplitude1 / 2 + offsety),
                 # Point(-period / 8 + offsetx, -0.707 * amplitude1 / 2 + offsety)
                 ]
    # make the polygon
    polygon = Polygon(pointlist)
    # draw the polygon on the layout
    polygon.draw(cell, perturabation_layer)

# 加载预先定义的光栅
# get a AEMD grating definition
TE_grating = MAKE_COMPONENT("FCGCTE_heyuv2.gds")  # port width = 0.5 um
TM_grating = MAKE_COMPONENT("FCGCTM_heyuv2.gds")  # port width = 0.45 um
grating_port_width_TE = 0.45
grating_port_width_TM = 0.5



# parameters
grating_gap = 160
radius = 40
taper_wg_length = 200


ring_width = 0.45
taper_for_grating_length = 50


# PBS parameters
bus_wg_length = 4000
bus_wg_width = 0.5
amplitude1 = 0.1

period1 =518e-3

edge1 = 0  # 微扰和中心的偏移量





def pbs(start_point_x, start_point_y, GRATING, grating_port_width):
    # <editor-fold desc = "bus_wg">
    # ##########bus_wg
    wg_start_point = Point(start_point_x, start_point_y)
    wg_end_point = wg_start_point + (bus_wg_length, 0)
    wg_bus = Waveguide(wg_start_point, wg_end_point, width=bus_wg_width)
    wg_bus.draw(cell, wg_layer)
    # </editor-fold>
    # <editor-fold desc = "Connect taper">
    grating_point = wg_start_point
    left_top_grating = GRATING(grating_point, RIGHT)
    left_top_grating.draw(cell)

    grating_point = wg_end_point
    left_top_grating = GRATING(grating_point, LEFT)
    left_top_grating.draw(cell)
    # # </editor-fold>
    # <editor-fold desc = "Perturbations">
    N1 = int(bus_wg_length / period1)


    start_x = wg_bus.get_start_point().x
    start_y = wg_bus.get_start_point().y
    # TE0 perturbations
    for i in range(N1):
        offSetx1 = i * period1 + start_x + 0.5 * period1
        offSety1 = -edge1 + start_y
        one_single_etch(amplitude1, period1, offSetx1, offSety1)
    for i in range(N1):
        offSetx1 = i * period1 + start_x + 0.5 * period1
        offSety1 = edge1 + start_y
        one_single_etch(amplitude1, period1, offSetx1, offSety1)
    # TE1 perturbations
    # for i in range(N2):
    #     offSetx2 = i * period2 + start_x + 0.5 * period2 + TE1_x
    #     offSety2 = -edge2 + start_y
    #     one_single_etch(amplitude2, period2, offSetx2, offSety2)
    #     offSety2 = edge2 + start_y
    #     one_single_etch(amplitude2, period2, offSetx2, offSety2)
    # # TE2 perturbations
    # for i in range(N3):
    #     offSetx3 = i * period3 + start_x + 0.5 * period3
    #     offSety3 = -edge3 + start_y
    #     one_single_etch(amplitude3, period3, offSetx3, offSety3)
    #     offSety3 = edge3 + start_y
    #     one_single_etch(amplitude3, period3, offSetx3, offSety3)
    # # </editor-fold>



amplitude1_initial = amplitude1

# col
for j in range(3):
    period1_initial = [358e-3, 378e-3,398e-3]
    # period1_initial = period1
    period1 = period1_initial[j]
#     amplitude1 = amplitude1_initial*amplitude_scale_ratio[j]
#     amplitude2 = amplitude2_initial*amplitude_scale_ratio[j]
#     amplitude3 = amplitude3_initial*amplitude_scale_ratio[j]
    # row
    # for i in range(10):
        # waveguide_scale_ratio = [0.96 ,0.97, 0.98, 0.99, 1, 1.01,1.02, 1.03,1.04,1.05]
        # bus_wg_width = waveguide_scale_ratio[i] * bus_wg_width_intial

        # grating_ref(480*(i), 0-1200*(j),TE_grating,grating_port_width_TE)
    pbs(      0, -600-1200*(j), TE_grating,grating_port_width_TE)
        # grating_ref(480*(i), -600-1200*(j),TM_grating,grating_port_width_TM)
        # pbs(        480*(i), -800-1200*(j),TM_grating,grating_port_width_TM)
        # text_start_point = Point(480*(i), 0-1200*(j)+10);
        # text = Text(text_start_point, str(10*(i+1)+j+1));
        # text.draw(cell, wg_layer);

make_gdsii_file("PBS_v2-2.gds", inv_source_layer=wg_layer,inv_target_layer=inv_layer)
# make_gdsii_file("PBS_v1.gds")
