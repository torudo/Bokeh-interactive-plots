"""This programm creates an .html with an
   interactive Bokeh plot to observe
   the behaviour of the function R(lambda)"""

import numpy as np
from bokeh.layouts import column, row
from bokeh.models import CustomJS, Slider
from bokeh.plotting import ColumnDataSource, figure, output_file, show

def R(x, f = -1, k = -0.2, M = 1):
    return (M*(-f*x**2 + 2*k*x**4) )

x = np.linspace(0.51, 5, 500)
x = 2*np.pi/x
y = R(2*np.pi/x)
lambda_0 = np.array([np.sqrt(8*0.2*np.pi**2 / (1) )])
y_0 = R(2*np.pi/lambda_0)

lambda_max = np.array([2*np.pi * np.sqrt(4*(0.2) / (1*1))])
y_p = R(2*np.pi/lambda_max)

source = ColumnDataSource(data=dict(x=x, y=y))
source_point = ColumnDataSource(data=dict(x=lambda_max, y=y_p))
source_0 = ColumnDataSource(data=dict(x=lambda_0, y=y_0))

plot = figure(x_range=(0, 10),
              y_range=(-0.7, 3),
              plot_width=500,
              plot_height=500,
              x_axis_label=r'lambda',
              y_axis_label = 'R(lambda)')

plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)
plot.circle('x', 'y', source=source_point, line_width=3,
                        line_alpha=0.6,
                        legend_label=r"lambda_max",
                        color="green", size=8)
plot.circle('x', 'y', source=source_0, line_width=3,
                        line_alpha=0.6,
                        legend_label=r"lambda_0",
                        color="red", size=8)

amp_slider = Slider(start=-1, end=2, value=1, step=.01, title="M")
freq_slider = Slider(start=-1, end=0.8, value=-0.2, step=.0005, title="k")
phase_slider = Slider(start=-5, end=2, value=-1, step=.0005, title="f")
#offset_slider = Slider(start=-5, end=5, value=0, step=.05, title="k_max")

callback = CustomJS(args=dict(source=source, source_p=source_point,source_0=source_0, amp=amp_slider, freq=freq_slider, phase=phase_slider),
                    code="""
    function f_add(x, A,k,phi) {
        return ( A*(-phi*Math.pow(Math.PI*2/x, 2)+2*k*Math.pow(Math.PI*2/x, 4) ) );   // The function returns the product of p1 and p2
    }

    const data = source.data;
    const A = amp.value;
    const k = freq.value;
    const phi = phase.value;
    const x = data['x']
    const y = data['y']

    var data_p = source_p.data;
    var x_p = data_p['x']
    var y_p = data_p['y']

    var data_0 = source_0.data;
    var x_0 = data_0['x']
    var y_0 = data_0['y']

    for (var i = 0; i < x.length; i++) {
        y[i] = f_add(x[i],A,k,phi) ;
        x_p[i] = 2*Math.PI*Math.sqrt(4*k/(A*phi)) ;
        y_p[i] =  f_add(x_p[i],A,k,phi) ;

        x_0[i] = Math.sqrt( 8*k*Math.pow(Math.PI,2) / (A*phi) ) ;
        y_0[i] =  f_add(x_0[i],A,k,phi) ;

    }

    source.change.emit();
    source_p.change.emit();
    source_0.change.emit();
""")


amp_slider.js_on_change('value', callback)
freq_slider.js_on_change('value', callback)
phase_slider.js_on_change('value', callback)
#offset_slider.js_on_change('value', callback)

layout = row(
    plot,
    column(amp_slider, freq_slider, phase_slider),
)

output_file("R_slider.html", title="Rslider.py example")

show(layout)
