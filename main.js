/* <!--
    ########################################
    ## @author Benjamin Thomas Schwertfeger (October 2021)
    ## copyright by Benjamin Thomas Schwertfeger (October 2021)
    ## https://b-schwertfeger.de
    ############
--> */

/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
// --> IMPORTS
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

const chi = require("chi-squared");
const {
    fft,
    ifft,
    dft,
    idft
} = require("fft-js");


/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
// --> HELPER FUNCTIONS
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

const convolve = (vec1, vec2) => {
    // taken from https://gist.github.com/jdpigeon/1de0b43eed7ae7e4080818cad53be200
    if (vec1.length === 0 || vec2.length === 0) {
        throw new Error("Vectors can not be empty!");
    }
    const volume = vec1;
    const kernel = vec2;
    let displacement = 0;
    const convVec = [];

    for (let i = 0; i < volume.length; i++) {
        for (let j = 0; j < kernel.length; j++) {
            if (displacement + j !== convVec.length) {
                convVec[displacement + j] =
                    convVec[displacement + j] + volume[i] * kernel[j];
            } else {
                convVec.push(volume[i] * kernel[j]);
            }
        }
        displacement++;
    }

    return convVec;
};

function rep(arr, n) {
    // repeat array n times
    let output = new Array(n * arr.length);
    for (let i = 0; i < output.length; i++) {
        output[i] = arr[i % arr.length];
    }
    return output
}

function getAvg(grades) {
    const total = grades.reduce((acc, c) => acc + c, 0);
    return total / grades.length;
}

function spec_pgram() {
    // https://github.com/telmo-correa/time-series-analysis/blob/master/Python/spectrum.py
    // AWI-Workspace -> my notebook
}
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
// --> CALCULATION
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

window.TIMESTEPS = 5000;
window.orbital_global = {
    kyear0: new Array(),
    ecc: new Array(),
    epsilon: new Array(),
    omega: new Array()
};

$(document).ready(function () {
    // === Load orbital parameters (given each kyr for 0-5Mya) ===
    // Load the matrix contains data from Berger and Loutre (1991),
    // downloaded as ORBIT91 from ncdc.noaa.gov

    $.ajax({
        type: "GET",
        url: "orbital_param.csv",
        dataType: "text",
        success: function (data) {
            let kyear = [...new Array(window.TIMESTEPS)].map((elem, index) => index / 10);
            processData(data, kyear);
        }
    });
});

// const Spline = require("cubic-spline");

function processData(allText, kyear) {
    let allTextLines = allText.split(/\r\n|\n/);
    let headers = allTextLines[0].split(' ');

    let kyear0 = new Array(0),
        ecc0 = new Array(0),
        epsilon0 = new Array(0),
        omega0 = new Array(0);

    for (let i = 1; i < allTextLines.length; i++) {
        let data = allTextLines[i].split(' ');
        for (let j = 0; j < headers.length; j++) {
            switch (j) {
                case 0:
                    kyear0.push(-parseInt(data[j]));
                    break;
                case 1:
                    ecc0.push(parseFloat(data[j])); // eccentricity
                    break;
                case 2:
                    // longitude of perihelion (precession angle)
                    // add 180 degrees to omega (see lambda definition, Berger 1978 Appendix)
                    omega0.push(parseFloat(data[j]) + 180);
                    break;
                case 3:
                    epsilon0.push(parseFloat(data[j]));
            }
        }
    }
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ---------- ----- ----- */
    // --> BEGIN TESTING
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ---------- ----- ----- */

    // var ky0TEST = [1, 2, 3, 4, 5, 6];
    // var ecc0TEST = [1, 5, 4, 5, 1, 6];
    // let kyearTEST = [...new Array(60)].map((elem, index) => (1 + index) / 10);
    // const spline = new Spline(ky0TEST, ecc0TEST);

    // let res = [...new Array(kyearTEST.length)].map((elem, index) => 0)
    // for (let i = 0; i < kyearTEST.length; i++) {
    //     res[i] = spline.at(parseFloat(kyear[i]));
    // }
    // console.log(kyearTEST)
    // console.log(res)

    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ---------- ----- ----- */
    // --> END TESTING
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ---------- ----- ----- */

    // remove discontinuities(360 degree jumps)
    omega0 = unwrap(omega0.map((elem) => elem * Math.PI / 180)).map((elem) => elem * 180 / Math.PI);
    // Interpolate to requested dates
    // const eccSpline = new Spline(kyear0, ecc0),
    //     omegaSpline = new Spline(kyear0, omega0),
    //     epsilonSpline = new Spline(kyear0, epsilon0);

    // for (let i = 0; i < kyear.length; i++) {
    //     window.orbital_global.ecc.push(eccSpline.at(parseFloat(kyear[i])));
    //     window.orbital_global.omega.push(omegaSpline.at(parseFloat(kyear[i])) * Math.PI / 180);
    //     window.orbital_global.epsilon.push(epsilonSpline.at(parseFloat(kyear[i])) * Math.PI / 180);
    // }

    for (let i = 0; i < window.TIMESTEPS; i++) {
        window.orbital_global.ecc.push(ecc0[i]);
        window.orbital_global.omega.push(omega0[i] * Math.PI / 180);
        window.orbital_global.epsilon.push(epsilon0[i] * Math.PI / 180);
    }
    // console.log(window.orbital_global.ecc)
    // console.log(window.orbital_global.omega)
    // console.log(window.orbital_global.epsilon)

    plotALL()

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
}


function unwrap(p) {
    //  Q = unwrap(P) corrects the radian phase angles in array P by adding multiples of ±2pi when absolute jumps between consecutive array elements are greater than pi radians.
    //  based on http://ccrma.stanford.edu/~jos/sasp/Matlab_listing_unwrap_m.html

    let N = p.length;
    let up = [...new Array(N)].map(() => 0);
    let pm1 = p[0];
    up[0] = pm1;
    let po = 0,
        thr = Math.PI,
        pi2 = 2 * Math.PI;

    for (let i = 1; i < N; i++) {
        let cp = p[i] + po;
        let dp = cp - pm1;
        pm1 = cp;
        if (dp >= thr) {
            while (dp >= thr) {
                po = po - pi2;
                dp = dp - pi2;
            }
        }
        if (dp <= ((-1) * thr)) {
            while (dp <= thr) {
                po = po + pi2;
                dp = dp + pi2;
            }
        }
        cp = p[i] + po
        pm1 = cp
        up[i] = cp
    }
    return up
}

function tlag(data, ilag) {
    let temp = rep(data, 3);
    for (let i = 366; i < 731; i++) {
        temp[i] -= ilag;
    }
    return temp
}

// Cover functions to return the daily insolation, either using a classical calendar (aligned with March21) or using an alignment with Dec21 summer solstice
function insolMarch21(kyear, LAT) {
    let res = new Array(0);
    for (let day = 1; day < 365 + 1; day++) {
        res.push(daily_insolation(kyear, LAT, day).Fsw);
    }
    return res;
}

function insolDec21(kyear, LAT) {
    let res = {
        Fsw: new Array(0),
        lambda: new Array(0),
    };
    for (let day = 1; day < 365 + 1; day++) { //Eigentlich passt 356 besser
        let tmp = daily_insolation(kyear, LAT, day);
        res.Fsw.push(tmp.Fsw);
        res.lambda.push(tmp.lambda);
    }
    let shift = 355 - Math.min.apply(Math, res.lambda.map((elem) => Math.abs(elem - 270))); //dann entprechend hinschieben
    return tlag(res.Fsw, shift);
}

function insolDec21_param(ecc, obliquity, long_perh, LAT) {
    let res = {
        Fsw: new Array(0),
        lambda: new Array(0),
    };
    for (let day = 1; day < 365 + 1; day++) { //Eigentlich passt 356 besser
        let tmp = daily_insolation_param(LAT, day, ecc, obliquity, long_perh);
        res.Fsw.push(tmp.Fsw);
        res.lambda.push(tmp.lambda);
    }

    let shift = 355 - Math.min.apply(Math, res.lambda.map((elem) => Math.abs(elem - 270))); // which.min(abs(r.lambda - 270)) //dann entprechend hinschieben1
    return tlag(res.Fsw, shift)
}

function daily_insolation_param(lat, day, ecc, obliquity, long_perh, day_type = 1) {
    // Insolation, converted and adapted from Huybers Code, based on Berger 1991

    // === Get orbital parameters ===
    let epsilon = obliquity * Math.PI / 180;
    let omega = long_perh * Math.PI / 180;

    // === Calculate insolation ===
    lat = lat * Math.PI / 180 // latitude


    // lambda (or solar longitude) is the angular distance along Earth's orbit measured from spring equinox (21 March)
    let lambda = null;
    if (day_type === 1) { // calendar days
        // estimate lambda from calendar day using an approximation from Berger 1978 section 3
        let delta_lambda_m = (day[i] - 80) * 2 * Math.PI / 365.2422; // lambda bei gleich langen Tagen
        let beta = Math.pow((1 - Math.pow(ecc, 2)), 1 / 2);
        let lambda_m0 = (-2) * ((1 / 2 * ecc + 1 / 8 * Math.pow(ecc, 3)) * (1 + beta) * Math.sin(-omega) - 1 / 4 * Math.pow(ecc, 2) * (1 / 2 + beta) * Math.sin(-2 * omega) + 1 / 8 * Math.pow(ecc, 3) * (1 / 3 + beta) * (Math.sin(-3 * omega)));
        let lambda_m = lambda_m0 + delta_lambda_m;

        lambda = lambda_m + (2 * ecc - 1 / 4 * Math.pow(ecc, 3)) * Math.sin(lambda_m - omega) + (5 / 4) * Math.pow(ecc, 2) * Math.sin(2 * (lambda_m - omega)) + (13 / 12) * Math.pow(ecc, 3) * Math.sin(3 * (lambda_m - omega));
    } else if (day_type === 2) { // solar longitude (1-360)
        lambda = day[i] * 2 * Math.PI / 360; // lambda=0 for spring equinox
    } else {
        console.log("was geschieht hier?");
    }

    let So = 1365; // solar constant (W/m^2)
    let delta = Math.asin(Math.sin(epsilon) * Math.sin(lambda)); // declination of the sun
    let Ho = Math.acos(-Math.tan(lat) * Math.tan(delta)); // hour angle at sunrise/sunset

    // no sunrise or no sunset: Berger 1978 eqn(8), (9)
    if ((Math.abs(lat) >= Math.PI / 2 - Math.abs(delta)) && (lat * delta > 0)) {
        Ho = Math.PI;
    } else {
        Ho = 0;
    }

    // Insolation: Berger 1978 eq(10)
    let Fsw = So / Math.PI * Math.pow((1 + ecc * Math.cos(lambda - omega)), 2) / Math.pow((1 - Math.pow(ecc, 2)), 2) * (Ho * Math.sin(lat) * Math.sin(delta) + Math.cos(lat) * Math.cos(delta) * Math.sin(Ho));

    return {
        Fsw: Fsw,
        ecc: ecc,
        obliquity: obliquity,
        long_perh: long_perh,
        lambda: lambda / 2 / Math.PI * 360
    };
}


function daily_insolation(kyear, lat, day, day_type = 1, fast = true) {
    // CALCULATE DAILY INSOLATION

    // === Get orbital parameters ===
    let temp = {};
    if (fast) {
        temp = orbital_parameters_fast(kyear)
        // console.log(temp)
    } else {
        temp.ecc = window.orbital_global.ecc[kyear]
        temp.epsilon = window.orbital_global.epsilon[kyear]
        temp.omega = window.orbital_global.omega[kyear]
    }
    let ecc = temp.ecc,
        epsilon = temp.epsilon,
        omega = temp.omega;

    // For output of orbital parameters
    let obliquity = epsilon * 180 / Math.PI,
        long_perh = omega * 180 / Math.PI;

    // console.log(ecc, epsilon, omega, obliquity, long_perh)
    // === Calculate insolation ===
    lat = lat * Math.PI / 180 // latitude

    // lambda(or solar longitude) is the angular distance along Earth 's orbit measured from spring equinox (21 March)
    let lambda = null;
    if (day_type == 1) { //calendar days 
        // estimate lambda from calendar day using an approximation from Berger 1978 section 3
        const delta_lambda_m = (day - 80) * 2 * Math.PI / 365.2422;
        const beta = Math.pow((1 - Math.pow(ecc, 2)), (1 / 2))
        const lambda_m0 = (-2) * ((1 / 2 * ecc + 1 / 8 * Math.pow(ecc, 3)) * (1 + beta) * Math.sin(-omega) - 1 / 4 * Math.pow(ecc, 2) * (1 / 2 + beta) * Math.sin(-2 * omega) + 1 / 8 * Math.pow(ecc, 3) * (1 / 3 + beta) * (Math.sin(-3 * omega)))
        const lambda_m = lambda_m0 + delta_lambda_m
        lambda = lambda_m + (2 * ecc - 1 / 4 * Math.pow(ecc, 3)) * Math.sin(lambda_m - omega) + (5 / 4) * Math.pow(ecc, 2) * Math.sin(2 * (lambda_m - omega)) + (13 / 12) * Math.pow(ecc, 3) * Math.sin(3 * (lambda_m - omega))
    } else if (day_type == 2) { // solar longitude(1 - 360) {
        lambda = day * 2 * Math.PI / 360 // lambda = 0  for spring equinox
    } else {
        console.log("was geschieht hier?");
    }

    let So = 1365; // solar constant(W / m ^ 2)
    let delta = Math.asin(Math.sin(epsilon) * Math.sin(lambda)); // declination of the sun
    let Ho = Math.acos(-Math.tan(lat) * Math.tan(delta)); // hour angle at sunrise / sunset

    // no sunrise or no sunset: Berger 1978 eqn(8), (9)
    if (Math.abs(lat) >= (Math.PI / 2 - Math.abs(delta))) {
        if (lat * delta > 0) {
            Ho = Math.PI;
        } else {
            Ho = 0;
        }
    }

    // Insolation: Berger 1978 eq(10)
    //Fsw=So/pi*(1+ecc*cos(lambda-omega))^2 /(1-ecc^2)^2 * ( Ho*sin(lat)*sin(delta) + cos(lat)*cos(delta)*sin(Ho))
    let Fsw = So / Math.PI * Math.pow(1 + ecc * Math.cos(lambda - omega), 2) / Math.pow(1 - Math.pow(ecc, 2), 2) * (Ho * Math.sin(lat) * Math.sin(delta) + Math.cos(lat) * Math.cos(delta) * Math.sin(Ho));

    return {
        Fsw: Fsw,
        ecc: ecc,
        obliquity: obliquity,
        long_perh: long_perh,
        lambda: lambda / 2 / Math.PI * 360
    };
}

function annual_insolation(kyear, lat) {
    // CALCULATIOE ANNUAL INSOLATION
    let result = new Array(kyear.length);
    for (let year = 0; year < result.length; year++) {
        let daysInYearInsolation = new Array(365);
        for (let day = 0; day < 365; day++) {
            daysInYearInsolation[day] = daily_insolation(kyear[year], lat, day).Fsw
        }
        result[year] = getAvg(daysInYearInsolation);
    }
    return result;
}

function orbital_parameters_fast(kyear) {
    return {
        ecc: window.orbital_global.ecc[parseInt(kyear * 10)],
        epsilon: window.orbital_global.epsilon[parseInt(kyear * 10)],
        omega: window.orbital_global.omega[parseInt(kyear * 10)]
    }
}

/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
// --> PLOTTING
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

function plotALL(input = null) {
    let day = 172,
        lat = 65;
    if (input !== null) {
        day = input.day, lat = input.lat
    }

    console.log(day, lat)
    // ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    const time = [...new Array(window.TIMESTEPS)].map((elem, index) => -index)

    const values = new Array();
    let dailyInsolatoinResult = {
        Fsw: new Array(),
        ecc: new Array(),
        obliquity: new Array(),
        lambda: new Array(),
        long_perh: new Array()
    };

    for (let year = 0; year < 5000; year++) {
        const res = daily_insolation(year, lat, day, 1, false) // false or true for fast and not fast
        dailyInsolatoinResult.Fsw.push(res.Fsw);
        dailyInsolatoinResult.ecc.push(res.ecc);
        dailyInsolatoinResult.obliquity.push(res.obliquity);
        dailyInsolatoinResult.long_perh.push(res.long_perh);
        dailyInsolatoinResult.lambda.push(res.lambda);
    }

    let default_config = {
        type: 'line',
        data: {
            labels: time,
            datasets: [],
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                title: {
                    display: true,
                    text: "",
                    font: {
                        Family: window.font_famliy,
                        size: 18,
                    },
                },
                legend: {
                    position: 'top',
                    display: false,
                },
                tooltip: {
                    usePointStyle: true,
                    callbacks: {
                        labelPointStyle: function (context) {
                            return {
                                pointStyle: 'rectRot',
                                rotation: 0,
                            };
                        },
                    },
                },
            },
            scales: {
                x: {
                    display: true,
                    title: {
                        display: true,
                        text: '',
                        font: {
                            family: window.font_famliy,
                            size: 16,
                        },
                    },
                    // reverse: true
                },
                y: {
                    display: true,
                    title: {
                        display: true,
                        text: '',
                        font: {
                            family: window.font_famliy,
                            size: 16,
                        },
                    },
                },
            },
            animations: {
                radius: {
                    duration: 400,
                    easing: 'linear',
                    loop: (ctx) => ctx.activate,
                },
            },
            hoverRadius: 8,
            hoverBackgroundColor: 'yellow',
            interaction: {
                mode: 'nearest',
                intersect: false,
                axis: 'x',
            },
        },
    };

    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
    // 1. PLOT
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

    const insol5000Data = {
        label: "Insolation for 5.000 kys",
        data: dailyInsolatoinResult.Fsw,
        fill: false,
        borderColor: 'rgb(255, 0, 0)',
        pointRadius: 0,
        tension: 0.1,
        borderWidth: 2
    };

    document.getElementById('orbital_line_plot_1').remove();
    document.getElementById('orbital_line_plot_1_container').innerHTML = '<canvas id="orbital_line_plot_1"></canvas>';
    let ctx1 = document.getElementById('orbital_line_plot_1');

    let config1 = {
        ...default_config
    };

    config1.data.datasets = [insol5000Data];
    config1.options.plugins.title.text = "Insolation for 5.000 kys";
    config1.options.plugins.legend.display = false;
    config1.options.scales.x.title.text = "-(1:5000)";
    config1.options.scales.y.title.text = "june.65N.new";

    window.orbital_line_plot_1 = new Chart(ctx1, config1);

    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
    // 2. PLOT
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

    const insol5000max510 = new Array();
    for (let i = 0; i < dailyInsolatoinResult.Fsw.length; i++) {
        if (dailyInsolatoinResult.Fsw[i] > 510) {
            insol5000max510.push(510);
        } else {
            insol5000max510.push(dailyInsolatoinResult.Fsw[i]);
        }
    }

    const insol5000max510Data = {
        label: "Insolation for 5.000 kys y>510 = 510",
        data: insol5000max510,
        fill: false,
        borderColor: 'rgb(255, 0, 0)',
        pointRadius: 0,
        tension: 0.1,
        borderWidth: 2
    };

    const meanOfInsol = getAvg(insol5000max510);

    let meanInsol5000max510Data = {
        label: "Mean",
        data: [...new Array(window.TIMESTEPS)].map(() => meanOfInsol),
        borderColor: "black",
        pointRadius: 0,
        borderDash: [10, 5],
        fill: false,
        borderWidth: 1
    };

    document.getElementById('orbital_line_plot_2').remove();
    document.getElementById('orbital_line_plot_2_container').innerHTML = '<canvas id="orbital_line_plot_2"></canvas>';
    let ctx2 = document.getElementById('orbital_line_plot_2');
    insol5000Data.label = "Insolation for 5.000 kys "

    let config2 = {
        ...default_config
    };

    config2.data.datasets = [meanInsol5000max510Data, insol5000max510Data];
    config2.options.plugins.title.text = "overflowed, non-linear wave";
    config2.options.plugins.legend.display = true;
    config2.options.scales.x.title.text = "-(1:5000)";
    config2.options.scales.y.title.text = "Wave";

    window.orbital_line_plot_2 = new Chart(ctx2, config2);

    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
    // 3. PLOT
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

    // plot(time, ecc.new, col = "red", type = 'l')

    document.getElementById('orbital_line_plot_3').remove();
    document.getElementById('orbital_line_plot_3_container').innerHTML = '<canvas id="orbital_line_plot_3"></canvas>';
    const ctx3 = document.getElementById('orbital_line_plot_3');

    const dailyInsol_ecc = {
        label: "Eccentricity",
        data: dailyInsolatoinResult.ecc,
        fill: false,
        borderColor: 'rgb(255, 0, 0)',
        pointRadius: 0,
        tension: 0.1,
        borderWidth: 2
    };

    let config3 = {
        ...default_config
    };

    config3.data.datasets = [dailyInsol_ecc];
    config3.options.plugins.title.text = "Plot of orbital parameters";
    config3.options.plugins.legend.display = false;
    config3.options.scales.x.title.text = "-(1:5000)";
    config3.options.scales.y.title.text = "Eccentricity";

    window.orbital_line_plot_3 = new Chart(ctx3, config3);

    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
    // 4. PLOT
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

    document.getElementById('orbital_line_plot_4').remove();
    document.getElementById('orbital_line_plot_4_container').innerHTML = '<canvas id="orbital_line_plot_4"></canvas>';
    const ctx4 = document.getElementById('orbital_line_plot_4');

    const dailyInsol_obliquity = {
        label: "Obliquity",
        data: dailyInsolatoinResult.obliquity,
        fill: false,
        borderColor: 'blue',
        pointRadius: 0,
        tension: 0.1,
        borderWidth: 2
    };

    let config4 = {
        ...default_config
    };

    config4.data.datasets = [dailyInsol_obliquity];
    config4.options.plugins.title.text = "Plot of orbital parameters";
    config4.options.plugins.legend.display = false;
    config4.options.scales.x.title.text = "-(1:5000)";
    config4.options.scales.y.title.text = "Obliquity";

    window.orbital_line_plot_4 = new Chart(ctx4, config4);

    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
    // 5. PLOT
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

    document.getElementById('orbital_line_plot_5').remove();
    document.getElementById('orbital_line_plot_5_container').innerHTML = '<canvas id="orbital_line_plot_5"></canvas>';
    const ctx5 = document.getElementById('orbital_line_plot_5');

    const dailyInsol_lambda = {
        label: "Lambda",
        data: dailyInsolatoinResult.lambda,
        fill: false,
        borderColor: 'black',
        pointRadius: 0,
        tension: 0.1,
        borderWidth: 2
    };

    let config5 = {
        ...default_config
    }

    config5.data.datasets = [dailyInsol_lambda];
    config5.options.plugins.title.text = "Plot of orbital parameters";
    config5.options.plugins.legend.display = false;
    config5.options.scales.x.title.text = "-(1:5000)";
    config5.options.scales.y.title.text = "Lambda";

    window.orbital_line_plot_5 = new Chart(ctx5, config5);

    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
    // 6. PLOT
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

    // xx <-spectrum(obliquity.new,spans=5,main="spans=5"); abline(v=0:5/100,h=0:1000/1000,lty=3);


    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
    // 7. PLOT
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

    // xxx<-spec.ar(obliquity.new, order = 1)

    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
    // 8. PLOT
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

    // plot(xx)
    // lines(xxx$freq,xxx$spec)

    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
    // 10. PLOT
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

    // xx <-spectrum(wave,spans=5,main="spans=5")

    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
    // 11. PLOT
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

    // xxx<-spec.ar(wave, order = 1)

    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
    // 12. PLOT
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

    // plot(xx,col="green"); abline(v=0:5/100,h=0:1000/1000,lty=3);
    // lines(xxx$freq,xxx$spec)

}


/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
// --> MAKE BUTTONS AND SLIDER DO WHAT THEY SHOULD DO
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

const
    orbital_lat_slide = document.getElementById("orbital_lat_slide"),
    orbital_day_slide = document.getElementById("orbital_day_slide"),
    orbital_slide_value_fields = document.getElementsByName("orbital_slide_value"),
    orbital_slider = document.getElementsByName("orbital_slide"),
    orbital_plot_variables = ["lat", "day"];

window.default_orbital_values = {
    lat: 65,
    day: 172
}
// ----- ----- ----- ----- ----- ----- ----- ----- -----

const orbital_RESET_BTN = document.getElementById('orbital_resetBtn');
orbital_RESET_BTN.onclick = function () {
    plotALL(); // resets the plot

    orbital_lat_slide.value = window.default_orbital_values.lat,
        orbital_day_slide.value = window.default_orbital_values.day;

    orbital_slide_value_fields.forEach((element, index) => { // Reset value fields
        const default_value = window.default_orbital_values[orbital_plot_variables[index]];
        document.getElementById(element.id).innerHTML = default_value;
    });
}

// ----- ----- ----- ----- ----- ----- ----- ----- -----
for (let entry = 0; entry < orbital_slider.length; entry++) {
    orbital_slider[entry].oninput = function () {
        let elem_id = orbital_slider[entry].id;
        elem_id = elem_id.substring(0, elem_id.length - 5)
        document.getElementById(elem_id + "value").innerHTML = document.getElementById(orbital_slider[entry].id).value;
    }
    orbital_slider[entry].onchange = function () {
        plotALL({
            day: parseInt(orbital_day_slide.value),
            lat: parseInt(orbital_lat_slide.value)
        });
    }
}


/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
// EOF
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */