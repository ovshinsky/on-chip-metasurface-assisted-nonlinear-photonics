import math

from splayout import *

##define cell
cell = Cell("SWG_TE0_REFLECTION_MWBG_300nm-filter")
##define layer
wg_layer = Layer(1, 0)
grating_layer = Layer(2, 0)
perturbation_layer = Layer(3, 0)
heater_layer = Layer(4, 0)
pad_layer = Layer(5, 0)
inv_layer = Layer(6, 0)

##define parameters:
#MWBG parameters:
widthBusWaveguide= 3
lengthBusWaveguide = 4010
lengthGrating = 4010
heightPerturbationSquare = 0.2
periodPerturbationSquare = 0.4657
gapBetweenBusWaveguideAndPerturbation = 1
apodizationIndex = 0

#other parameters:
widthOfEdgeCouplerOutput = 0.15
widthOfADCCouplerOutput = widthOfEdgeCouplerOutput

#caculated parameters:
numberPerturbation = int(lengthGrating/periodPerturbationSquare)

text_start_point = Point(0,20)
text = Text(text_start_point,"SWG_TE0_REFLECTION_MWBG_300nm-filter")
text.draw(cell,wg_layer)

for k in range (11):
    periodPerturbationSquare = 0.4657 + (k-5)*0.005
    numberPerturbation = int(lengthGrating / periodPerturbationSquare)
    y_shift = k*50

    #MWBG部分：
    ##bus waveguide
    pointlist = [   Point(0,-widthBusWaveguide/2+y_shift),
                    Point(0,widthBusWaveguide/2+y_shift),
                    Point(lengthBusWaveguide,widthBusWaveguide/2+y_shift),
                    Point(lengthBusWaveguide,-widthBusWaveguide/2+y_shift),
                    ]
    # make the polygon
    polygon = Polygon(pointlist)
    # draw the polygon on the layout
    polygon.draw(cell, wg_layer)


    #MWBG右侧的连接：

    tp_start_point = Point(lengthGrating,y_shift)
    tp_end_point = tp_start_point - (100,0)
    #
    # make a taper
    tp = Taper(tp_start_point,tp_end_point,start_width=widthBusWaveguide,end_width=0)
    # draw the taper on the layout
    tp.draw(cell,wg_layer)

    tp_start_point = Point(lengthGrating,y_shift)
    tp_end_point = tp_start_point + (600,0)
    tp = Taper(tp_start_point,tp_end_point,start_width=widthBusWaveguide,end_width=widthOfEdgeCouplerOutput)
    tp.draw(cell,wg_layer)


    #MWBG左侧的连接：
    tp_start_point = Point(0,y_shift)
    tp_end_point = tp_start_point + (100,0)
    # make a taper
    tp = Taper(tp_start_point,tp_end_point,start_width=widthBusWaveguide,end_width=0)
    # draw the taper on the layout
    tp.draw(cell,wg_layer)

    tp_start_point = Point(0,y_shift)
    tp_end_point = tp_start_point - (600,0)
    tp = Taper(tp_start_point,tp_end_point,start_width=widthBusWaveguide,end_width=widthOfADCCouplerOutput)
    tp.draw(cell,wg_layer)


    #perturbation gratings:
    for i in range (numberPerturbation):
        deltas=periodPerturbationSquare/4*math.exp(-apodizationIndex*((i-numberPerturbation/2)/numberPerturbation)**2);
        squareStartPointX = periodPerturbationSquare*i + deltas
        squareStartPointY = widthBusWaveguide/2 + gapBetweenBusWaveguideAndPerturbation + heightPerturbationSquare/2 +y_shift
        wg_start_point = Point(squareStartPointX, squareStartPointY)
        wg_end_point = wg_start_point + (periodPerturbationSquare/2, 0)
        wg = Waveguide(wg_start_point,wg_end_point,width=heightPerturbationSquare)
        wg.draw(cell, wg_layer)

        squareStartPointX = periodPerturbationSquare*i + deltas
        squareStartPointY = -widthBusWaveguide/2 - gapBetweenBusWaveguideAndPerturbation - heightPerturbationSquare/2+y_shift
        wg_start_point = Point(squareStartPointX, squareStartPointY)
        wg_end_point = wg_start_point + (periodPerturbationSquare/2, 0)
        wg = Waveguide(wg_start_point,wg_end_point,width=heightPerturbationSquare)
        wg.draw(cell, wg_layer)


make_gdsii_file("SWG_TE0_REFLECTION_MWBG_300nm-filter.gds")
