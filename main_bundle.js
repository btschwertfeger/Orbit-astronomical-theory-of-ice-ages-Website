(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){
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
},{"chi-squared":4,"fft-js":5}],2:[function(require,module,exports){
/**
 * Bit twiddling hacks for JavaScript.
 *
 * Author: Mikola Lysenko
 *
 * Ported from Stanford bit twiddling hack library:
 *    http://graphics.stanford.edu/~seander/bithacks.html
 */

"use strict"; "use restrict";

//Number of bits in an integer
var INT_BITS = 32;

//Constants
exports.INT_BITS  = INT_BITS;
exports.INT_MAX   =  0x7fffffff;
exports.INT_MIN   = -1<<(INT_BITS-1);

//Returns -1, 0, +1 depending on sign of x
exports.sign = function(v) {
  return (v > 0) - (v < 0);
}

//Computes absolute value of integer
exports.abs = function(v) {
  var mask = v >> (INT_BITS-1);
  return (v ^ mask) - mask;
}

//Computes minimum of integers x and y
exports.min = function(x, y) {
  return y ^ ((x ^ y) & -(x < y));
}

//Computes maximum of integers x and y
exports.max = function(x, y) {
  return x ^ ((x ^ y) & -(x < y));
}

//Checks if a number is a power of two
exports.isPow2 = function(v) {
  return !(v & (v-1)) && (!!v);
}

//Computes log base 2 of v
exports.log2 = function(v) {
  var r, shift;
  r =     (v > 0xFFFF) << 4; v >>>= r;
  shift = (v > 0xFF  ) << 3; v >>>= shift; r |= shift;
  shift = (v > 0xF   ) << 2; v >>>= shift; r |= shift;
  shift = (v > 0x3   ) << 1; v >>>= shift; r |= shift;
  return r | (v >> 1);
}

//Computes log base 10 of v
exports.log10 = function(v) {
  return  (v >= 1000000000) ? 9 : (v >= 100000000) ? 8 : (v >= 10000000) ? 7 :
          (v >= 1000000) ? 6 : (v >= 100000) ? 5 : (v >= 10000) ? 4 :
          (v >= 1000) ? 3 : (v >= 100) ? 2 : (v >= 10) ? 1 : 0;
}

//Counts number of bits
exports.popCount = function(v) {
  v = v - ((v >>> 1) & 0x55555555);
  v = (v & 0x33333333) + ((v >>> 2) & 0x33333333);
  return ((v + (v >>> 4) & 0xF0F0F0F) * 0x1010101) >>> 24;
}

//Counts number of trailing zeros
function countTrailingZeros(v) {
  var c = 32;
  v &= -v;
  if (v) c--;
  if (v & 0x0000FFFF) c -= 16;
  if (v & 0x00FF00FF) c -= 8;
  if (v & 0x0F0F0F0F) c -= 4;
  if (v & 0x33333333) c -= 2;
  if (v & 0x55555555) c -= 1;
  return c;
}
exports.countTrailingZeros = countTrailingZeros;

//Rounds to next power of 2
exports.nextPow2 = function(v) {
  v += v === 0;
  --v;
  v |= v >>> 1;
  v |= v >>> 2;
  v |= v >>> 4;
  v |= v >>> 8;
  v |= v >>> 16;
  return v + 1;
}

//Rounds down to previous power of 2
exports.prevPow2 = function(v) {
  v |= v >>> 1;
  v |= v >>> 2;
  v |= v >>> 4;
  v |= v >>> 8;
  v |= v >>> 16;
  return v - (v>>>1);
}

//Computes parity of word
exports.parity = function(v) {
  v ^= v >>> 16;
  v ^= v >>> 8;
  v ^= v >>> 4;
  v &= 0xf;
  return (0x6996 >>> v) & 1;
}

var REVERSE_TABLE = new Array(256);

(function(tab) {
  for(var i=0; i<256; ++i) {
    var v = i, r = i, s = 7;
    for (v >>>= 1; v; v >>>= 1) {
      r <<= 1;
      r |= v & 1;
      --s;
    }
    tab[i] = (r << s) & 0xff;
  }
})(REVERSE_TABLE);

//Reverse bits in a 32 bit word
exports.reverse = function(v) {
  return  (REVERSE_TABLE[ v         & 0xff] << 24) |
          (REVERSE_TABLE[(v >>> 8)  & 0xff] << 16) |
          (REVERSE_TABLE[(v >>> 16) & 0xff] << 8)  |
           REVERSE_TABLE[(v >>> 24) & 0xff];
}

//Interleave bits of 2 coordinates with 16 bits.  Useful for fast quadtree codes
exports.interleave2 = function(x, y) {
  x &= 0xFFFF;
  x = (x | (x << 8)) & 0x00FF00FF;
  x = (x | (x << 4)) & 0x0F0F0F0F;
  x = (x | (x << 2)) & 0x33333333;
  x = (x | (x << 1)) & 0x55555555;

  y &= 0xFFFF;
  y = (y | (y << 8)) & 0x00FF00FF;
  y = (y | (y << 4)) & 0x0F0F0F0F;
  y = (y | (y << 2)) & 0x33333333;
  y = (y | (y << 1)) & 0x55555555;

  return x | (y << 1);
}

//Extracts the nth interleaved component
exports.deinterleave2 = function(v, n) {
  v = (v >>> n) & 0x55555555;
  v = (v | (v >>> 1))  & 0x33333333;
  v = (v | (v >>> 2))  & 0x0F0F0F0F;
  v = (v | (v >>> 4))  & 0x00FF00FF;
  v = (v | (v >>> 16)) & 0x000FFFF;
  return (v << 16) >> 16;
}


//Interleave bits of 3 coordinates, each with 10 bits.  Useful for fast octree codes
exports.interleave3 = function(x, y, z) {
  x &= 0x3FF;
  x  = (x | (x<<16)) & 4278190335;
  x  = (x | (x<<8))  & 251719695;
  x  = (x | (x<<4))  & 3272356035;
  x  = (x | (x<<2))  & 1227133513;

  y &= 0x3FF;
  y  = (y | (y<<16)) & 4278190335;
  y  = (y | (y<<8))  & 251719695;
  y  = (y | (y<<4))  & 3272356035;
  y  = (y | (y<<2))  & 1227133513;
  x |= (y << 1);
  
  z &= 0x3FF;
  z  = (z | (z<<16)) & 4278190335;
  z  = (z | (z<<8))  & 251719695;
  z  = (z | (z<<4))  & 3272356035;
  z  = (z | (z<<2))  & 1227133513;
  
  return x | (z << 2);
}

//Extracts nth interleaved component of a 3-tuple
exports.deinterleave3 = function(v, n) {
  v = (v >>> n)       & 1227133513;
  v = (v | (v>>>2))   & 3272356035;
  v = (v | (v>>>4))   & 251719695;
  v = (v | (v>>>8))   & 4278190335;
  v = (v | (v>>>16))  & 0x3FF;
  return (v<<22)>>22;
}

//Computes next combination in colexicographic order (this is mistakenly called nextPermutation on the bit twiddling hacks page)
exports.nextCombination = function(v) {
  var t = v | (v - 1);
  return (t + 1) | (((~t & -~t) - 1) >>> (countTrailingZeros(v) + 1));
}


},{}],3:[function(require,module,exports){
var LogGamma = require('gamma').log

// The following code liberated from
// http://www.math.ucla.edu/~tom/distributions/chisq.html

function Gcf(X, A) { // Good for X>A+1
  var A0 = 0;
  var B0 = 1;
  var A1 = 1;
  var B1 = X;
  var AOLD = 0;
  var N = 0;
  while (Math.abs((A1 - AOLD) / A1) > .00001) {
    AOLD = A1;
    N = N + 1;
    A0 = A1 + (N - A) * A0;
    B0 = B1 + (N - A) * B0;
    A1 = X * A0 + N * A1;
    B1 = X * B0 + N * B1;
    A0 = A0 / B1;
    B0 = B0 / B1;
    A1 = A1 / B1;
    B1 = 1;
  }
  var Prob = Math.exp(A * Math.log(X) - X - Math.LogGamma(A)) * A1;

  return 1 - Prob
}

function Gser(X, A) { // Good for X<A+1.
  var T9 = 1 / A;
  var G = T9;
  var I = 1;
  while (T9 > G * .00001) {
    T9 = T9 * X / (A + I);
    G = G + T9;
    I = I + 1;
  }
  G = G * Math.exp(A * Math.log(X) - X - Math.LogGamma(A));

  return G
}

function Gammacdf(x, a) {
  var GI;
  if (x <= 0) {
    GI = 0
  } else if (x < a + 1) {
    GI = Gser(x, a)
  } else {
    GI = Gcf(x, a)
  }
  return GI
}

module.exports = function (Z, DF) {
  if (DF <= 0) {
    throw new Error("Degrees of freedom must be positive")
  }
  return Gammacdf(Z / 2, DF / 2)
}
},{"gamma":12}],4:[function(require,module,exports){
var gamma = require('gamma');

exports.pdf = function (x, k_) {
    if (x < 0) return 0;
    var k = k_ / 2;
    return 1 / (Math.pow(2, k) * gamma(k))
        * Math.pow(x, k - 1)
        * Math.exp(-x / 2)
    ;
};

exports.cdf = require('./cdf')

},{"./cdf":3,"gamma":12}],5:[function(require,module,exports){
/*===========================================================================*\
 * Fast Fourier Transform (Cooley-Tukey Method)
 *
 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
 *
 * This code is not designed to be highly optimized but as an educational
 * tool to understand the Fast Fourier Transform.
\*===========================================================================*/
module.exports = {
    fft: require('./src/fft').fft,
    ifft: require('./src/ifft').ifft,
    fftInPlace: require('./src/fft').fftInPlace,
    util: require('./src/fftutil'),
    dft: require('./src/dft'),
    idft: require('./src/idft')
};

},{"./src/dft":7,"./src/fft":8,"./src/fftutil":9,"./src/idft":10,"./src/ifft":11}],6:[function(require,module,exports){
//-------------------------------------------------
// Add two complex numbers
//-------------------------------------------------
var complexAdd = function (a, b)
{
    return [a[0] + b[0], a[1] + b[1]];
};

//-------------------------------------------------
// Subtract two complex numbers
//-------------------------------------------------
var complexSubtract = function (a, b)
{
    return [a[0] - b[0], a[1] - b[1]];
};

//-------------------------------------------------
// Multiply two complex numbers
//
// (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
//-------------------------------------------------
var complexMultiply = function (a, b) 
{
    return [(a[0] * b[0] - a[1] * b[1]), 
            (a[0] * b[1] + a[1] * b[0])];
};

//-------------------------------------------------
// Calculate |a + bi|
//
// sqrt(a*a + b*b)
//-------------------------------------------------
var complexMagnitude = function (c) 
{
    return Math.sqrt(c[0]*c[0] + c[1]*c[1]); 
};

//-------------------------------------------------
// Exports
//-------------------------------------------------
module.exports = {
    add: complexAdd,
    subtract: complexSubtract,
    multiply: complexMultiply,
    magnitude: complexMagnitude
};

},{}],7:[function(require,module,exports){
/*===========================================================================*\
 * Discrete Fourier Transform (O(n^2) brute-force method)
 *
 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
 *
 * This code is not designed to be highly optimized but as an educational
 * tool to understand the Fast Fourier Transform.
\*===========================================================================*/

//------------------------------------------------
// Note: this code is not optimized and is
// primarily designed as an educational and testing
// tool.
//------------------------------------------------
var complex = require('./complex');
var fftUtil = require('./fftutil');

//-------------------------------------------------
// Calculate brute-force O(n^2) DFT for vector.
//-------------------------------------------------
var dft = function(vector) {
  var X = [],
      N = vector.length;

  for (var k = 0; k < N; k++) {
    X[k] = [0, 0]; //Initialize to a 0-valued complex number.

    for (var i = 0; i < N; i++) {
      var exp = fftUtil.exponent(k * i, N);
      var term;
      if (Array.isArray(vector[i]))
        term = complex.multiply(vector[i], exp)//If input vector contains complex numbers
      else
        term = complex.multiply([vector[i], 0], exp);//Complex mult of the signal with the exponential term.  
      X[k] = complex.add(X[k], term); //Complex summation of X[k] and exponential
    }
  }

  return X;
};

module.exports = dft;
},{"./complex":6,"./fftutil":9}],8:[function(require,module,exports){
/*===========================================================================*\
 * Fast Fourier Transform (Cooley-Tukey Method)
 *
 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
 *
 * This code is not designed to be highly optimized but as an educational
 * tool to understand the Fast Fourier Transform.
\*===========================================================================*/

//------------------------------------------------
// Note: Some of this code is not optimized and is
// primarily designed as an educational and testing
// tool.
// To get high performace would require transforming
// the recursive calls into a loop and then loop
// unrolling. All of this is best accomplished
// in C or assembly.
//-------------------------------------------------

//-------------------------------------------------
// The following code assumes a complex number is
// an array: [real, imaginary]
//-------------------------------------------------
var complex = require('./complex'),
    fftUtil = require('./fftutil'),
    twiddle = require('bit-twiddle');

module.exports = {
  //-------------------------------------------------
  // Calculate FFT for vector where vector.length
  // is assumed to be a power of 2.
  //-------------------------------------------------
  fft: function fft(vector) {
    var X = [],
        N = vector.length;

    // Base case is X = x + 0i since our input is assumed to be real only.
    if (N == 1) {
      if (Array.isArray(vector[0])) //If input vector contains complex numbers
        return [[vector[0][0], vector[0][1]]];      
      else
        return [[vector[0], 0]];
    }

    // Recurse: all even samples
    var X_evens = fft(vector.filter(even)),

        // Recurse: all odd samples
        X_odds  = fft(vector.filter(odd));

    // Now, perform N/2 operations!
    for (var k = 0; k < N / 2; k++) {
      // t is a complex number!
      var t = X_evens[k],
          e = complex.multiply(fftUtil.exponent(k, N), X_odds[k]);

      X[k] = complex.add(t, e);
      X[k + (N / 2)] = complex.subtract(t, e);
    }

    function even(__, ix) {
      return ix % 2 == 0;
    }

    function odd(__, ix) {
      return ix % 2 == 1;
    }

    return X;
  },
  //-------------------------------------------------
  // Calculate FFT for vector where vector.length
  // is assumed to be a power of 2.  This is the in-
  // place implementation, to avoid the memory
  // footprint used by recursion.
  //-------------------------------------------------
  fftInPlace: function(vector) {
    var N = vector.length;

    var trailingZeros = twiddle.countTrailingZeros(N); //Once reversed, this will be leading zeros

    // Reverse bits
    for (var k = 0; k < N; k++) {
      var p = twiddle.reverse(k) >>> (twiddle.INT_BITS - trailingZeros);
      if (p > k) {
        var complexTemp = [vector[k], 0];
        vector[k] = vector[p];
        vector[p] = complexTemp;
      } else {
        vector[p] = [vector[p], 0];
      }
    }

    //Do the DIT now in-place
    for (var len = 2; len <= N; len += len) {
      for (var i = 0; i < len / 2; i++) {
        var w = fftUtil.exponent(i, len);
        for (var j = 0; j < N / len; j++) {
          var t = complex.multiply(w, vector[j * len + i + len / 2]);
          vector[j * len + i + len / 2] = complex.subtract(vector[j * len + i], t);
          vector[j * len + i] = complex.add(vector[j * len + i], t);
        }
      }
    }
  }
};

},{"./complex":6,"./fftutil":9,"bit-twiddle":2}],9:[function(require,module,exports){
/*===========================================================================*\
 * Fast Fourier Transform Frequency/Magnitude passes
 *
 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
 *
 * This code is not designed to be highly optimized but as an educational
 * tool to understand the Fast Fourier Transform.
\*===========================================================================*/

//-------------------------------------------------
// The following code assumes a complex number is
// an array: [real, imaginary]
//-------------------------------------------------
var complex = require('./complex');


//-------------------------------------------------
// By Eulers Formula:
//
// e^(i*x) = cos(x) + i*sin(x)
//
// and in DFT:
//
// x = -2*PI*(k/N)
//-------------------------------------------------
var mapExponent = {},
    exponent = function (k, N) {
      var x = -2 * Math.PI * (k / N);

      mapExponent[N] = mapExponent[N] || {};
      mapExponent[N][k] = mapExponent[N][k] || [Math.cos(x), Math.sin(x)];// [Real, Imaginary]

      return mapExponent[N][k];
};

//-------------------------------------------------
// Calculate FFT Magnitude for complex numbers.
//-------------------------------------------------
var fftMag = function (fftBins) {
    var ret = fftBins.map(complex.magnitude);
    return ret.slice(0, ret.length / 2);
};

//-------------------------------------------------
// Calculate Frequency Bins
// 
// Returns an array of the frequencies (in hertz) of
// each FFT bin provided, assuming the sampleRate is
// samples taken per second.
//-------------------------------------------------
var fftFreq = function (fftBins, sampleRate) {
    var stepFreq = sampleRate / (fftBins.length);
    var ret = fftBins.slice(0, fftBins.length / 2);

    return ret.map(function (__, ix) {
        return ix * stepFreq;
    });
};

//-------------------------------------------------
// Exports
//-------------------------------------------------
module.exports = {
    fftMag: fftMag,
    fftFreq: fftFreq,
    exponent: exponent
};

},{"./complex":6}],10:[function(require,module,exports){
/*===========================================================================*\
 * Inverse Discrete Fourier Transform (O(n^2) brute-force method)
 *
 * (c) Maximilian Bügler. 2016
 *
 * Based on and using the code by
 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
 *
 * This code is not designed to be highly optimized but as an educational
 * tool to understand the Fast Fourier Transform.
\*===========================================================================*/

//------------------------------------------------
// Note: Some of this code is not optimized and is
// primarily designed as an educational and testing
// tool.
//-------------------------------------------------

//-------------------------------------------------
// The following code assumes a complex number is
// an array: [real, imaginary]
//-------------------------------------------------
var dft = require('./dft');

function idft(signal) {
    //Interchange real and imaginary parts
    var csignal = [];
    for (var i = 0; i < signal.length; i++) {
        csignal[i] = [signal[i][1], signal[i][0]];
    }

    //Apply dft
    var ps = dft(csignal);

    //Interchange real and imaginary parts and normalize
    var res = [];
    for (var j = 0; j < ps.length; j++) {
        res[j] = [ps[j][1] / ps.length, ps[j][0] / ps.length];
    }
    return res;
}

module.exports = idft;
},{"./dft":7}],11:[function(require,module,exports){
/*===========================================================================*\
 * Inverse Fast Fourier Transform (Cooley-Tukey Method)
 *
 * (c) Maximilian Bügler. 2016
 *
 * Based on and using the code by
 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
 *
 * This code is not designed to be highly optimized but as an educational
 * tool to understand the Fast Fourier Transform.
\*===========================================================================*/

//------------------------------------------------
// Note: Some of this code is not optimized and is
// primarily designed as an educational and testing
// tool.
// To get high performace would require transforming
// the recursive calls into a loop and then loop
// unrolling. All of this is best accomplished
// in C or assembly.
//-------------------------------------------------

//-------------------------------------------------
// The following code assumes a complex number is
// an array: [real, imaginary]
//-------------------------------------------------

var fft = require('./fft').fft;


module.exports = {
    ifft: function ifft(signal){
        //Interchange real and imaginary parts
        var csignal=[];
        for(var i=0; i<signal.length; i++){
            csignal[i]=[signal[i][1], signal[i][0]];
        }
    
        //Apply fft
        var ps=fft(csignal);
        
        //Interchange real and imaginary parts and normalize
        var res=[];
        for(var j=0; j<ps.length; j++){
            res[j]=[ps[j][1]/ps.length, ps[j][0]/ps.length];
        }
        return res;
    }
};

},{"./fft":8}],12:[function(require,module,exports){
// transliterated from the python snippet here:
// http://en.wikipedia.org/wiki/Lanczos_approximation

var g = 7;
var p = [
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7
];

var g_ln = 607/128;
var p_ln = [
    0.99999999999999709182,
    57.156235665862923517,
    -59.597960355475491248,
    14.136097974741747174,
    -0.49191381609762019978,
    0.33994649984811888699e-4,
    0.46523628927048575665e-4,
    -0.98374475304879564677e-4,
    0.15808870322491248884e-3,
    -0.21026444172410488319e-3,
    0.21743961811521264320e-3,
    -0.16431810653676389022e-3,
    0.84418223983852743293e-4,
    -0.26190838401581408670e-4,
    0.36899182659531622704e-5
];

// Spouge approximation (suitable for large arguments)
function lngamma(z) {

    if(z < 0) return Number('0/0');
    var x = p_ln[0];
    for(var i = p_ln.length - 1; i > 0; --i) x += p_ln[i] / (z + i);
    var t = z + g_ln + 0.5;
    return .5*Math.log(2*Math.PI)+(z+.5)*Math.log(t)-t+Math.log(x)-Math.log(z);
}

module.exports = function gamma (z) {
    if (z < 0.5) {
        return Math.PI / (Math.sin(Math.PI * z) * gamma(1 - z));
    }
    else if(z > 100) return Math.exp(lngamma(z));
    else {
        z -= 1;
        var x = p[0];
        for (var i = 1; i < g + 2; i++) {
            x += p[i] / (z + i);
        }
        var t = z + g + 0.5;

        return Math.sqrt(2 * Math.PI)
            * Math.pow(t, z + 0.5)
            * Math.exp(-t)
            * x
        ;
    }
};

module.exports.log = lngamma;

},{}]},{},[1]);
