import math
from scipy.io import loadmat
from splayout import *
file = 'x'
#mat_dtype = True
data = loadmat(file,matdtype = True)
# I0 = data['']


##define cell
cell = Cell("ADC_SiN_200nm_TE1_1550nm")
##define layer
wg_layer = Layer(1, 0)
grating_layer = Layer(2, 0)
perturabation_layer = Layer(3, 0)
heater_layer = Layer(4, 0)
pad_layer = Layer(5, 0)
inv_layer = Layer(6, 0)

##define parameters:
#ADC parameters:
leftWidthBusWaveguide= 3
rightWidthBusWaveguide = 2
leftWidthAccessWaveguide = 0
rightWidthAccessWaveguide = 1.462
gapBetweenBusAndAccessWaveguide = 0.463
couplingLength = 200

#other parameters:
radiusOfConnectedBend = 100
widthOfEdgeCouplerOutput = 1.52


text_start_point = Point(0,20)
text = Text(text_start_point,"ADC_SiN_200nm_TE1_1550nm")
text.draw(cell,wg_layer)


#ADC部分：
pointlist = [   Point(0,0),
                Point(0,leftWidthBusWaveguide),
                Point(couplingLength,rightWidthBusWaveguide),
                Point(couplingLength,0),
                ]
# make the polygon
polygon = Polygon(pointlist)
# draw the polygon on the layout
polygon.draw(cell, wg_layer)

pointlist = [   Point(0,-gapBetweenBusAndAccessWaveguide),
                Point(0,-gapBetweenBusAndAccessWaveguide-leftWidthAccessWaveguide),
                Point(couplingLength,-gapBetweenBusAndAccessWaveguide-rightWidthAccessWaveguide),
                Point(couplingLength,-gapBetweenBusAndAccessWaveguide),
                ]
# make the polygon
polygon = Polygon(pointlist)
# draw the polygon on the layout
polygon.draw(cell, wg_layer)


#ADC右下部分的连接：
double_connect_start_point = Point(couplingLength,-gapBetweenBusAndAccessWaveguide-rightWidthAccessWaveguide/2)
double_connect_end_point = double_connect_start_point + (2*radiusOfConnectedBend+0.1,-2*radiusOfConnectedBend-0.1)
connector = DoubleBendConnector(double_connect_start_point,double_connect_end_point, width=rightWidthAccessWaveguide, radius=radiusOfConnectedBend)
connector.draw(cell,wg_layer)

tp_start_point = connector.get_end_point()
tp_end_point = tp_start_point + (100,0)
# make a taper
tp = Taper(tp_start_point,tp_end_point,start_width=rightWidthAccessWaveguide,end_width=widthOfEdgeCouplerOutput)
# draw the taper on the layout
tp.draw(cell,wg_layer)
connector.get_end_point()

#ADC右上部分的连接：
tp_start_point = Point(couplingLength,rightWidthBusWaveguide/2)
tp_end_point = tp_start_point + (300.2,0)
# make a taper
tp = Taper(tp_start_point,tp_end_point,start_width=rightWidthBusWaveguide,end_width=widthOfEdgeCouplerOutput)
# draw the taper on the layout
tp.draw(cell,wg_layer)
connector.get_end_point()

make_gdsii_file("ADC_SiN_200nm_TE1_1550nm.gds")

#
# ##################perturbations: 函数
# def one_single_etch(amplitude1, period, offsetx, offsety):
#     # points for the polygon
#
#     pointlist = [Point(-period / 4 + offsetx, amplitude1/2 + offsety),
#                  Point(period / 4 + offsetx, amplitude1 / 2 + offsety),
#                  Point(period / 4 + offsetx, -amplitude1 / 2 + offsety),
#                  Point(-period / 4 + offsetx, -amplitude1 / 2 + offsety),
#                  # Point(period / 8 + offsetx, -0.707 * amplitude1 / 2 + offsety),
#                  # Point(0 + offsetx, -amplitude1 / 2 + offsety),
#                  # Point(-period / 8 + offsetx, -0.707 * amplitude1 / 2 + offsety)
#                  ]
#
#     # make the polygon
#     polygon = Polygon(pointlist)
#     # draw the polygon on the layout
#     polygon.draw(cell, perturabation_layer)
#
# # 加载预先定义的光栅
# # get a AEMD grating definition
# TE_grating = MAKE_COMPONENT("FCGCTE_heyuv2.gds")  # port width = 0.5 um
# TM_grating = MAKE_COMPONENT("FCGCTM_heyuv2.gds")  # port width = 0.45 um
# grating_port_width_TE = 0.45
# grating_port_width_TM = 0.5
#
#
#
# # parameters
# grating_gap = 160
# radius = 40
# taper_wg_length = 200
#
#
# ring_width = 0.45
# taper_for_grating_length = 50
#
#
# # PBS parameters
# bus_wg_length = 4000
# bus_wg_width = 0.45
# amplitude1 = 0.45
#
# period1 =518e-3
#
# edge1 = 0  # 微扰和中心的偏移量
#
#
#
#
#
# def pbs(start_point_x, start_point_y, GRATING, grating_port_width):
#     # <editor-fold desc = "bus_wg">
#     # ##########bus_wg
#     wg_start_point = Point(start_point_x, start_point_y)
#     wg_end_point = wg_start_point + (bus_wg_length, 0)
#     wg_bus = Waveguide(wg_start_point, wg_end_point, width=bus_wg_width)
#     wg_bus.draw(cell, wg_layer)
#     # </editor-fold>
#     # <editor-fold desc = "Connect taper">
#     grating_point = wg_start_point
#     left_top_grating = GRATING(grating_point, RIGHT)
#     left_top_grating.draw(cell)
#
#     grating_point = wg_end_point
#     left_top_grating = GRATING(grating_point, LEFT)
#     left_top_grating.draw(cell)
#     # # </editor-fold>
#     # <editor-fold desc = "Perturbations">
#     N1 = int(bus_wg_length / period1)
#
#
#     start_x = wg_bus.get_start_point().x
#     start_y = wg_bus.get_start_point().y
#     # TE0 perturbations
#     for i in range(N1):
#         offSetx1 = i * period1 + start_x + 0.5 * period1
#         offSety1 = -edge1 + start_y
#         one_single_etch(amplitude1, period1, offSetx1, offSety1)
#     for i in range(N1):
#         offSetx1 = i * period1 + start_x + 0.5 * period1
#         offSety1 = edge1 + start_y
#         one_single_etch(amplitude1, period1, offSetx1, offSety1)
#     # TE1 perturbations
#     # for i in range(N2):
#     #     offSetx2 = i * period2 + start_x + 0.5 * period2 + TE1_x
#     #     offSety2 = -edge2 + start_y
#     #     one_single_etch(amplitude2, period2, offSetx2, offSety2)
#     #     offSety2 = edge2 + start_y
#     #     one_single_etch(amplitude2, period2, offSetx2, offSety2)
#     # # TE2 perturbations
#     # for i in range(N3):
#     #     offSetx3 = i * period3 + start_x + 0.5 * period3
#     #     offSety3 = -edge3 + start_y
#     #     one_single_etch(amplitude3, period3, offSetx3, offSety3)
#     #     offSety3 = edge3 + start_y
#     #     one_single_etch(amplitude3, period3, offSetx3, offSety3)
#     # # </editor-fold>
#
#
#
# amplitude1_initial = amplitude1
#
# # col
# for j in range(3):
#     period1_initial = [498e-3, 518e-3,538e-3]
#     # period1_initial = period1
#     period1 = period1_initial[j]
# #     amplitude1 = amplitude1_initial*amplitude_scale_ratio[j]
# #     amplitude2 = amplitude2_initial*amplitude_scale_ratio[j]
# #     amplitude3 = amplitude3_initial*amplitude_scale_ratio[j]
#     # row
#     # for i in range(10):
#         # waveguide_scale_ratio = [0.96 ,0.97, 0.98, 0.99, 1, 1.01,1.02, 1.03,1.04,1.05]
#         # bus_wg_width = waveguide_scale_ratio[i] * bus_wg_width_intial
#
#         # grating_ref(480*(i), 0-1200*(j),TE_grating,grating_port_width_TE)
#     pbs(      0, -600-1200*(j), TE_grating,grating_port_width_TE)
#         # grating_ref(480*(i), -600-1200*(j),TM_grating,grating_port_width_TM)
#         # pbs(        480*(i), -800-1200*(j),TM_grating,grating_port_width_TM)
#         # text_start_point = Point(480*(i), 0-1200*(j)+10);
#         # text = Text(text_start_point, str(10*(i+1)+j+1));
#         # text.draw(cell, wg_layer);
#
# # make_gdsii_file("PBS_v2-1.gds", inv_source_layer=wg_layer,inv_target_layer=inv_layer)
