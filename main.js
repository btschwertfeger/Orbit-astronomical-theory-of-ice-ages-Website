/**
 * Website to visualize and manipulate orbital parameter (teaching material)
 * 
 * © Alfred-Wegener-Institute Bremerhaven, Germany (2021)
 * @author Benjamin Thomas Schwertfeger (January 2022)
 * @email development@b-schwertfeger.de
 * @link https://github.com/btschwertfeger-AWI-Workspace/Orbit-astronomical-theory-of-ice-ages
 * 
 * most comments are taken from the original R implementation 
 **/

/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
// --> HELPER FUNCTIONS
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

function dateFromDay(year, day) {
    day -= 1;
    Date.prototype.addDays = function (days) {
        const date = new Date(this.valueOf());
        date.setDate(date.getDate() + days);
        return date;
    };

    let date = new Date(year, 0); // initialize a date in `year-01-01`
    let newdate = date.addDays(day)

    return newdate.toLocaleString('default', {
        month: 'long'
    }) + ', ' + newdate.getDate();
}

const convolve = (vec1, vec2) => {
    // (2021) source: https://gist.github.com/jdpigeon/1de0b43eed7ae7e4080818cad53be200
    if (vec1.length === 0 || vec2.length === 0) throw new Error('Vectors can not be empty!')

    const
        volume = vec1,
        kernel = vec2,
        convVec = [];

    let displacement = 0;

    for (let i = 0; i < volume.length; i++) {
        for (let j = 0; j < kernel.length; j++) {
            if (displacement + j !== convVec.length) convVec[displacement + j] = convVec[displacement + j] + volume[i] * kernel[j];
            else convVec.push(volume[i] * kernel[j]);
        }
        displacement++;
    }
    return convVec;
};

function get2dmax(arr) {
    return Math.max.apply(null, arr.map(function (row) {
        return Math.max.apply(Math, row);
    }));
}

function get2dmin(arr) {
    return Math.min.apply(null, arr.map(function (row) {
        return Math.min.apply(Math, row);
    }));
}

function rep(arr, n) {
    // repeat array n times
    let output = new Array(n * arr.length);
    for (let i = 0; i < output.length; i++) {
        output[i] = arr[i % arr.length];
    }
    return output
}

function avg(grades) {
    const total = grades.reduce((acc, c) => acc + c, 0);
    return total / grades.length;
}

/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
// --> CALCULATION
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
const
    FONT_FAMILY = 'Helvetica',
    TIMESTEPS = 5000;
window.orbital_global = {
    kyear0: new Array(),
    ecc: new Array(),
    epsilon: new Array(),
    omega: new Array()
};

$(document).ready(function () {
    document.getElementById('dayval').innerHTML = dateFromDay(2021, parseInt(document.getElementById('orbital_day_slide').value));
    // === Load orbital parameters (given each kyr for 0-5Mya) ===
    // Load the matrix contains data from Berger and Loutre (1991)

    $.ajax({
        type: 'GET',
        url: 'orbital_param.csv',
        dataType: 'text',
        success: function (data) {
            const kyear = [...new Array(TIMESTEPS)].map((elem, index) => index / 10);
            processData(data, kyear);
        }
    });
});

function processData(allText, kyear) {
    let allTextLines = allText.split(/\r\n|\n/);

    let
        headers = allTextLines[0].split(' '),
        kyear0 = new Array(0),
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

    for (let i = 0; i < TIMESTEPS; i++) {
        window.orbital_global.ecc.push(ecc0[i]);
        window.orbital_global.omega.push(omega0[i] * Math.PI / 180);
        window.orbital_global.epsilon.push(epsilon0[i] * Math.PI / 180);
    }

    plotALL(null, true); // line plots (parameter: data<null>, init<true>)

    plot_contour_1();
    plot_contour_2();
    plot_contour_3();

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
}

function unwrap(p) {
    //  Q = unwrap(P) corrects the radian phase angles in array P by adding multiples of ±2pi when absolute jumps between consecutive array elements are greater than pi radians.
    //  based on http://ccrma.stanford.edu/~jos/sasp/Matlab_listing_unwrap_m.html

    let N = p.length;
    let
        up = [...new Array(N)].map(() => 0),
        pm1 = p[0];
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
    //let shift = 355 - Math.min.apply(Math, res.lambda.map((elem) => Math.abs(elem - 270))); //dann entprechend hinschieben
    const min = Math.min.apply(Math, res.lambda.map((elem) => Math.abs(elem - 270)));
    for (let i = 0; i < res.lambda.length; i++) {
        if (res.lambda[i] == min)
            return tlag(res.Fsw, 355 - i);
    }
    return tlag(res.Fsw, shift);
}

function insolDec21_param(ecc, obliquity, long_perh, LAT) {
    let res = {
        Fsw: new Array(0),
        lambda: new Array(0),
    };
    for (let day = 1; day < 365 + 1; day++) {
        let tmp = daily_insolation_param(LAT, day, ecc, obliquity, long_perh);
        res.Fsw.push(tmp.Fsw);
        res.lambda.push(tmp.lambda);
    }

    let shift = 355 - Math.min.apply(Math, res.lambda.map((elem) => Math.abs(elem - 270))); // which.min(abs(r.lambda - 270)) //dann entprechend hinschieben1
    return tlag(res.Fsw, shift)
}

function daily_insolation_param(lat, day, ecc, obliquity, long_perh, day_type = 1) {
    /* Insolation, converted and adapted from Huybers Code, based on Berger 1991
    
    Description:
    Computes daily average insolation as a function of day and latitude at
    any point during the past 5 million years.

    Inputs:
    kyear:    Thousands of years before present (0 to 5000).
    lat:      Latitude in degrees (-90 to 90).
    day:      Indicator of time of year; calendar day by default.
    day_type: Convention for specifying time of year (+/- 1,2) [optional].
        day_type=1 (default): day input is calendar day (1-365.24), where day 1
        is January first.  The calendar is referenced to the vernal equinox
        which always occurs at day 80.
        day_type=2: day input is solar longitude (0-360 degrees). Solar
        longitude is the angle of the Earth's orbit measured from spring
        equinox (21 March). Note that calendar days and solar longitude are
        not linearly related because, by Kepler's Second Law, Earth's
        angular velocity varies according to its distance from the sun.
    Output:
    Fsw = Daily average solar radiation in W/m^2.
    Can also output orbital parameters.

    This script contains orbital parameter data for the past 50000 years
    from Berger and Loutre (1991).

    Detailed description of calculation:
    Values for eccentricity, obliquity, and longitude of perihelion for the
    past 5 Myr are taken from Berger and Loutre 1991 (data from
    ncdc.noaa.gov). If using calendar days, solar longitude is found using an
    approximate solution to the differential equation representing conservation
    of angular momentum (Kepler's Second Law).  Given the orbital parameters
    and solar longitude, daily average insolation is calculated exactly
    following Berger 1978.

    References:
    Berger A. and Loutre M.F. (1991). Insolation values for the climate of
        the last 10 million years. Quaternary Science Reviews, 10(4), 297-317.
    Berger A. (1978). Long-term variations of daily insolation and
        Quaternary climatic changes. Journal of Atmospheric Science, 35(12),
        2362-2367.

    Authors:
        Ian Eisenman and Peter Huybers, Harvard University, August 2006
        eisenman@fas.harvard.edu
        This file is available online at
        http://deas.harvard.edu/~eisenman/downloads
    Translated into JavaScript by Benjamin Thomas Schwertfeger
    Suggested citation:
        P. Huybers and I. Eisenman, 2006. Integrated summer insolation
        calculations. NOAA/NCDC Paleoclimatology Program Data
        Contribution #2006-079.
    
    */

    // // === Get orbital parameters ===
    // let epsilon = obliquity * Math.PI / 180;
    // let omega = long_perh * Math.PI / 180;

    // // === Calculate insolation ===
    // lat = lat * Math.PI / 180 // latitude

    // // lambda (or solar longitude) is the angular distance along Earth's orbit measured from spring equinox (21 March)
    // let lambda = null;
    // if (day_type === 1) { // calendar days
    //     // estimate lambda from calendar day using an approximation from Berger 1978 section 3
    //     let delta_lambda_m = (day[i] - 80) * 2 * Math.PI / 365.2422; // lambda bei gleich langen Tagen
    //     let beta = Math.pow((1 - Math.pow(ecc, 2)), 1 / 2);
    //     let lambda_m0 = (-2) * ((1 / 2 * ecc + 1 / 8 * Math.pow(ecc, 3)) * (1 + beta) * Math.sin(-omega) - 1 / 4 * Math.pow(ecc, 2) * (1 / 2 + beta) * Math.sin(-2 * omega) + 1 / 8 * Math.pow(ecc, 3) * (1 / 3 + beta) * (Math.sin(-3 * omega)));
    //     let lambda_m = lambda_m0 + delta_lambda_m;

    //     lambda = lambda_m + (2 * ecc - 1 / 4 * Math.pow(ecc, 3)) * Math.sin(lambda_m - omega) + (5 / 4) * Math.pow(ecc, 2) * Math.sin(2 * (lambda_m - omega)) + (13 / 12) * Math.pow(ecc, 3) * Math.sin(3 * (lambda_m - omega));
    // } else if (day_type === 2) { // solar longitude (1-360)
    //     lambda = day[i] * 2 * Math.PI / 360; // lambda=0 for spring equinox
    // } else console.log('was geschieht hier?');


    // let So = 1365; // solar constant (W/m^2)
    // let delta = Math.asin(Math.sin(epsilon) * Math.sin(lambda)); // declination of the sun
    // let Ho = Math.acos(-Math.tan(lat) * Math.tan(delta)); // hour angle at sunrise/sunset

    // // no sunrise or no sunset: Berger 1978 eqn(8), (9)
    // if ((Math.abs(lat) >= Math.PI / 2 - Math.abs(delta)) && (lat * delta > 0)) Ho = Math.PI;
    // else Ho = 0;


    // // Insolation: Berger 1978 eq(10)
    // let Fsw = So / Math.PI * Math.pow((1 + ecc * Math.cos(lambda - omega)), 2) / Math.pow((1 - Math.pow(ecc, 2)), 2) * (Ho * Math.sin(lat) * Math.sin(delta) + Math.cos(lat) * Math.cos(delta) * Math.sin(Ho));

    // return {
    //     Fsw: Fsw,
    //     ecc: ecc,
    //     obliquity: obliquity,
    //     long_perh: long_perh,
    //     lambda: lambda / 2 / Math.PI * 360
    // };
    return compute(ecc, obliquity, long_perh, lat, day, day_type);
}

function daily_insolation(kyear, lat, day, day_type = 1, fast = true) {
    /* CALCULATE DAILY INSOLATION
    Description:
    Computes daily average insolation as a function of day and latitude at
    any point during the past 5 million years.

    Inputs:
    kyear:    Thousands of years before present (0 to 5000).
    lat:      Latitude in degrees (-90 to 90).
    day:      Indicator of time of year; calendar day by default.
    day_type: Convention for specifying time of year (+/- 1,2) [optional].
        day_type=1 (default): day input is calendar day (1-365.24), where day 1
        is January first.  The calendar is referenced to the vernal equinox
        which always occurs at day 80.
        day_type=2: day input is solar longitude (0-360 degrees). Solar
        longitude is the angle of the Earth's orbit measured from spring
        equinox (21 March). Note that calendar days and solar longitude are
        not linearly related because, by Kepler's Second Law, Earth's
        angular velocity varies according to its distance from the sun.
    Output:
    Fsw = Daily average solar radiation in W/m^2.
    Can also output orbital parameters.

    This script contains orbital parameter data for the past 50000 years
    from Berger and Loutre (1991).

    Detailed description of calculation:
    Values for eccentricity, obliquity, and longitude of perihelion for the
    past 5 Myr are taken from Berger and Loutre 1991 (data from
    ncdc.noaa.gov). If using calendar days, solar longitude is found using an
    approximate solution to the differential equation representing conservation
    of angular momentum (Kepler's Second Law).  Given the orbital parameters
    and solar longitude, daily average insolation is calculated exactly
    following Berger 1978.

    References:
    Berger A. and Loutre M.F. (1991). Insolation values for the climate of
        the last 10 million years. Quaternary Science Reviews, 10(4), 297-317.
    Berger A. (1978). Long-term variations of daily insolation and
        Quaternary climatic changes. Journal of Atmospheric Science, 35(12),
        2362-2367.

    Authors:
        Ian Eisenman and Peter Huybers, Harvard University, August 2006
        eisenman@fas.harvard.edu
        This file is available online at
        http://deas.harvard.edu/~eisenman/downloads
        Translated into JavaScript by Benjamin Thomas Schwertfeger
    Suggested citation:
        P. Huybers and I. Eisenman, 2006. Integrated summer insolation
        calculations. NOAA/NCDC Paleoclimatology Program Data
        Contribution #2006-079.
    */

    // === Get orbital parameters ===
    let temp = {};
    if (fast) temp = orbital_parameters_fast(kyear)
    else {
        temp.ecc = window.orbital_global.ecc[kyear]
        temp.epsilon = window.orbital_global.epsilon[kyear]
        temp.omega = window.orbital_global.omega[kyear]
    }

    let
        ecc = temp.ecc,
        epsilon = temp.epsilon,
        omega = temp.omega;

    // For output of orbital parameters
    let
        obliquity = epsilon * 180 / Math.PI,
        long_perh = omega * 180 / Math.PI;

    var x = compute(ecc, obliquity, long_perh, lat, day, day_type)
    return x;
}

function compute(ecc, obliquity, long_perh, lat, day, day_type = 1) {
    let
        epsilon = (obliquity * Math.PI) / 180,
        omega = (long_perh * Math.PI) / 180;

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
    } else if (day_type == 2) lambda = day * 2 * Math.PI / 360 // lambda = 0  for spring equinox
    else console.log('was geschieht hier?');

    let So = 1365; // solar constant(W / m ^ 2)
    let delta = Math.asin(Math.sin(epsilon) * Math.sin(lambda)); // declination of the sun
    let Ho = Math.acos(-Math.tan(lat) * Math.tan(delta)); // hour angle at sunrise / sunset

    // no sunrise or no sunset: Berger 1978 eqn(8), (9)
    if (Math.abs(lat) >= (Math.PI / 2 - Math.abs(delta))) {
        if (lat * delta > 0) Ho = Math.PI;
        else Ho = 0;
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

function get_insolation_of_year(year) {
    let result = new Array(0);
    for (var lat = -90; lat < 90; lat++) {
        let inner = new Array(0);
        for (var day = 0; day < 365; day++) inner.push(daily_insolation(year, lat, day, 1, false).Fsw);
        result.push(inner)
    }
    return result;
}

function anno_insol_by_param(ecc, obliquity, long_perh) {
    // CALCULATIOE INSOLATION by parameter
    let result = new Array(0);
    for (var lat = -90; lat < 90; lat++) {
        var inner = new Array(0);
        for (var day = 0; day < 365; day++) inner.push(daily_insolation_param(lat, day, ecc, obliquity, long_perh).Fsw);
        result.push(inner);
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
// --> LINE PLOTS
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

function plotALL(input = null, init = false) {
    let
        day = 172,
        lat = 65;
    if (input !== null) day = input.day, lat = input.lat;

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    const time = [...new Array(TIMESTEPS)].map((elem, index) => -index)

    let dailyInsolationResult = {
        Fsw: new Array(),
        ecc: new Array(),
        obliquity: new Array(),
        lambda: new Array(),
        long_perh: new Array()
    };

    for (let year = 0; year < 5000; year++) {
        const res = daily_insolation(year, lat, day, 1, false) // false or true for fast and not fast
        dailyInsolationResult.Fsw.push(res.Fsw);
        dailyInsolationResult.ecc.push(res.ecc);
        dailyInsolationResult.obliquity.push(res.obliquity);
        dailyInsolationResult.long_perh.push(res.long_perh);
        dailyInsolationResult.lambda.push(res.lambda);
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
                    text: '',
                    font: {
                        Family: FONT_FAMILY,
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
                            family: FONT_FAMILY,
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
                            family: FONT_FAMILY,
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
        label: 'Insolation for 5.000ky',
        data: dailyInsolationResult.Fsw,
        fill: false,
        borderColor: 'rgb(255, 0, 0)',
        pointRadius: 0,
        tension: 0.1,
        borderWidth: 2
    };

    document.getElementById('orbital_line_plot_1').remove();
    document.getElementById('orbital_line_plot_1_container').innerHTML = '<canvas id=\'orbital_line_plot_1\'></canvas>';
    let ctx1 = document.getElementById('orbital_line_plot_1');

    let config1 = {
        ...default_config
    };

    config1.data.datasets = [insol5000Data];
    config1.options.plugins.title.text = 'Insolation for 5.000ky';
    config1.options.plugins.legend.display = false;
    config1.options.scales.x.title.text = 'ky';
    config1.options.scales.y.title.text = 'Insolation';

    window.orbital_line_plot_1 = new Chart(ctx1, config1);

    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
    // 2. PLOT
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

    const insol5000max510 = new Array();
    for (let i = 0; i < dailyInsolationResult.Fsw.length; i++) {
        if (dailyInsolationResult.Fsw[i] > 510) {
            insol5000max510.push(510);
        } else {
            insol5000max510.push(dailyInsolationResult.Fsw[i]);
        }
    }

    const insol5000max510Data = {
        label: 'Insolation for 5.000 ky (y > 510 => 510)',
        data: insol5000max510,
        fill: false,
        borderColor: 'rgb(255, 0, 0)',
        pointRadius: 0,
        tension: 0.1,
        borderWidth: 2
    };

    const meanOfInsol = avg(insol5000max510);

    let meanInsol5000max510Data = {
        label: 'Mean',
        data: [...new Array(TIMESTEPS)].map(() => meanOfInsol),
        borderColor: 'black',
        pointRadius: 0,
        borderDash: [10, 5],
        fill: false,
        borderWidth: 1
    };

    document.getElementById('orbital_line_plot_2').remove();
    document.getElementById('orbital_line_plot_2_container').innerHTML = '<canvas id=\'orbital_line_plot_2\'></canvas>';
    let ctx2 = document.getElementById('orbital_line_plot_2');
    insol5000Data.label = 'Insolation for 5.000ky'

    let config2 = {
        ...default_config
    };

    config2.data.datasets = [meanInsol5000max510Data, insol5000max510Data];
    config2.options.plugins.title.text = 'overflowed, non-linear wave';
    config2.options.plugins.legend.display = true;
    config2.options.scales.x.title.text = 'ky';
    config2.options.scales.y.title.text = 'Wave';

    window.orbital_line_plot_2 = new Chart(ctx2, config2);


    if (!init) return;
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
    // 3. PLOT
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

    // plot(time, ecc.new, col = 'red', type = 'l')

    document.getElementById('orbital_line_plot_3').remove();
    document.getElementById('orbital_line_plot_3_container').innerHTML = '<canvas id=\'orbital_line_plot_3\'></canvas>';
    const ctx3 = document.getElementById('orbital_line_plot_3');

    console.log(dailyInsolationResult.ecc)
    const dailyInsol_ecc = {
        label: 'Eccentricity',
        data: dailyInsolationResult.ecc,
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
    config3.options.plugins.title.text = 'Eccentricity';
    config3.options.plugins.legend.display = false;
    config3.options.scales.x.title.text = 'ky';
    config3.options.scales.y.title.text = 'Eccentricity';

    window.orbital_line_plot_3 = new Chart(ctx3, config3);

    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
    // 4. PLOT
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

    document.getElementById('orbital_line_plot_4').remove();
    document.getElementById('orbital_line_plot_4_container').innerHTML = '<canvas id=\'orbital_line_plot_4\'></canvas>';
    const ctx4 = document.getElementById('orbital_line_plot_4');

    const dailyInsol_obliquity = {
        label: 'Obliquity',
        data: dailyInsolationResult.obliquity,
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
    config4.options.plugins.title.text = 'Obliquity';
    config4.options.plugins.legend.display = false;
    config4.options.scales.x.title.text = 'ky';
    config4.options.scales.y.title.text = 'Obliquity';

    window.orbital_line_plot_4 = new Chart(ctx4, config4);

    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
    // 5. PLOT
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

    document.getElementById('orbital_line_plot_5').remove();
    document.getElementById('orbital_line_plot_5_container').innerHTML = '<canvas id=\'orbital_line_plot_5\'></canvas>';
    const ctx5 = document.getElementById('orbital_line_plot_5');

    const dailyInsol_lambda = {
        label: 'Lambda',
        data: dailyInsolationResult.lambda,
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
    config5.options.plugins.title.text = 'Lambda';
    config5.options.plugins.legend.display = false;
    config5.options.scales.x.title.text = 'ky';
    config5.options.scales.y.title.text = 'Lambda';

    window.orbital_line_plot_5 = new Chart(ctx5, config5);
}

const
    orbital_lat_slide = document.getElementById('orbital_lat_slide'),
    orbital_day_slide = document.getElementById('orbital_day_slide'),
    orbital_slide_value_fields = document.getElementsByName('orbital_slide_value'),
    orbital_slider = document.getElementsByName('orbital_slide'),
    orbital_RESET_BTN = document.getElementById('orbital_resetBtn'),
    orbital_plot_variables = ['lat', 'day'],
    default_orbital_values = {
        lat: 65,
        day: 172
    };
// ----- ----- ----- ----- ----- ----- ----- ----- -----

orbital_RESET_BTN.onclick = () => {
    plotALL(); // resets the plot

    orbital_lat_slide.value = default_orbital_values.lat,
        orbital_day_slide.value = default_orbital_values.day;

    orbital_slide_value_fields.forEach((element, index) => { // Reset value fields
        const default_value = default_orbital_values[orbital_plot_variables[index]];
        document.getElementById(element.id).innerHTML = default_value;
    });
}

orbital_lat_slide.oninput = function () {
    document.getElementById('orbital_lat_value').innerHTML = orbital_lat_slide.value;
}

// ----- ----- ----- ----- ----- ----- ----- ----- -----
for (let entry = 0; entry < orbital_slider.length; entry++) {
    orbital_slider[entry].onchange = () => {
        plotALL({
            day: parseInt(orbital_day_slide.value),
            lat: parseInt(orbital_lat_slide.value)
        });
    }
}

orbital_day_slide.oninput = function () {
    document.getElementById('dayval').innerHTML = dateFromDay(2021, parseInt(document.getElementById('orbital_day_slide').value));
    document.getElementById('orbital_day_value').innerHTML = orbital_day_slide.value;
}

/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
// --> C O N T O U R P L O T S
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

const
    MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'],
    colorscale_str = 'Jet',
    contour_plot_config = {
        toImageButtonOptions: {
            format: 'svg',
            filename: 'contour_1',
            width: 1920,
            height: 1080,
            scale: 1,
        }
    },
    contour_plot_layout = {
        title: {
            text: '',
            font: {
                family: FONT_FAMILY,
                size: 18,
            },
            xref: 'paper',
            x: 0.05,
        },
        xaxis: {
            title: {
                text: 'time',
                font: {
                    family: FONT_FAMILY,
                    size: 18,
                    color: '#7f7f7f',
                },
            },
            // tickformat: '%m',
            dtick: '30',
            /* Set the values at which ticks on this axis appear */
            tickvals: [15, 45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 345],
            /* Set the text displayed at the ticks position via tickvals */
            ticktext: [
                'January', 'February',
                'March', 'April', 'May',
                'June', 'July', 'August',
                'September', 'October', 'November',
                'December'
            ],
            /* Specifies the maximum number of ticks */
            nticks: 12,
        },
        yaxis: {
            title: {
                text: 'latitude',
                font: {
                    family: FONT_FAMILY,
                    size: 18,
                    color: '#7f7f7f',
                },
            },
        },
    };

function plot_contour(input) {
    Plotly.newPlot(
        input.divId,
        input.data,
        input.layout,
        contour_plot_config,
    );
}

/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
// --> CONTOURPLOT 1
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

function plot_contour_1(input = {
    year: 0
}) {
    input.year = (input.year < 0) ? input.year * -1 : input.year;
    const RESULT = get_insolation_of_year(input.year);
    const max = get2dmax(RESULT)

    const contour_plot_data = [{
        z: RESULT,
        x: [...new Array(RESULT[0].length)].map((elem, index) => index), //dateFromDay(2021, index + 1)), // time
        y: [...new Array(180)].map((elem, index) => index - 90), // latitude
        type: 'contour',
        colorscale: colorscale_str,
        line: {
            smoothing: 1
        },
        autocontour: false,
        colorbar: {
            title: 'Insolation',
            tickfont: {
                color: 'black',
            }
        },
        contours: {
            start: 0,
            end: max,
            size: 25,
        }
    }];

    var layout = Object.assign({}, contour_plot_layout);
    layout.title.text = `Insolation ${input.year}ky ago`;

    plot_contour({
        divId: 'contour_1',
        data: contour_plot_data,
        layout: layout
    });
}

const
    contour_1_slide_value = document.getElementById('x_input_contour_1_kyear'),
    contour_1_slider = document.getElementById('contour_1_kyear_slide');
contour_1_slide_value.onchange = function () {
    plot_contour_1({
        year: contour_1_slide_value.value
    });
}
contour_1_slider.oninput = () => {
    contour_1_slide_value.value = contour_1_slider.value
}
contour_1_slider.onchange = () => {
    plot_contour_1({
        year: contour_1_slide_value.value
    });
}

/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
// --> CONTOUR PLOT 2
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

function plot_contour_2(input = {
    year1: -125,
    year2: 0
}) {
    input.year1 = (input.year1 < 0) ? input.year1 * -1 : input.year1;
    input.year2 = (input.year2 < 0) ? input.year2 * -1 : input.year2;

    const
        res_1 = get_insolation_of_year(input.year1),
        res_2 = get_insolation_of_year(input.year2);

    let RESULT = [...new Array(res_1.length)].map((e, i) => [...new Array(res_1[i].length)].map((e, i) => 0));
    for (var i = 0; i < res_1.length; i++)
        for (var j = 0; j < res_1[0].length; j++) RESULT[i][j] = res_1[i][j] - res_2[i][j];


    const
        max = get2dmax(RESULT),
        min = get2dmin(RESULT);

    const range = (max > 0) ? max + (min < 0) ? min * -1 : min : max * -1 + -1 * min;
    const contour_plot_data = [{
        z: RESULT,
        x: [...new Array(RESULT[0].length)].map((elem, index) => index), //dateFromDay(2021, index + 1)), // time
        y: [...new Array(180)].map((elem, index) => index - 90), // latitude
        type: 'contour',
        colorscale: colorscale_str,
        line: {
            smoothing: 1
        },
        autocontour: false,
        colorbar: {
            title: 'Insolation',
            tickfont: {
                color: 'black',
            }
        },
        contours: {
            start: min,
            end: max,
            size: range / 15
        }
    }];

    var layout = Object.assign({}, contour_plot_layout);
    layout.title.text = `Anomalie of insolation: ${input.year1}ky minus ${input.year2}ky ago`;

    plot_contour({
        divId: 'contour_2',
        data: contour_plot_data,
        layout: layout
    });
}

const
    contour_2_year1 = document.getElementById('x_input_contour_2_kyear1'),
    contour_2_year2 = document.getElementById('x_input_contour_2_kyear2'),
    contour_2_kyear1_slide = document.getElementById('contour_2_kyear1_slide'),
    contour_2_kyear2_slide = document.getElementById('contour_2_kyear2_slide'),
    contour_2_resetBtn = document.getElementById('contour_2_resetBtn');

contour_2_kyear1_slide.oninput = () => {
    contour_2_year1.value = contour_2_kyear1_slide.value
}
contour_2_kyear2_slide.oninput = () => {
    contour_2_year2.value = contour_2_kyear2_slide.value
}

[contour_2_kyear1_slide, contour_2_kyear2_slide].forEach((e) => {
    e.onchange = () => {
        plot_contour_2({
            year1: contour_2_year1.value,
            year2: contour_2_year2.value
        });
    }
});

[contour_2_year1, contour_2_year2].forEach((e => {
    e.onchange = function () {
        plot_contour_2({
            year1: contour_2_year1.value,
            year2: contour_2_year2.value
        });
    }
}));
contour_2_resetBtn.onclick = () => {
    const
        y1 = -125,
        y2 = 0;
    contour_2_year1.value = y1;
    contour_2_year2.value = y2;
    contour_2_kyear1_slide.value = y1;
    contour_2_kyear2_slide.value = y2;
    plot_contour_2({
        year1: y1,
        year2: y2
    });
}

/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
// --> CONTOUR PLOT 3
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

const cplt3_default_values = {
    ecc: 0.0268550999,
    obl: 22.525,
    long_perh: -83619.1209
}

function plot_contour_3(ecc = cplt3_default_values.ecc, obl = cplt3_default_values.obl, long_perh = cplt3_default_values.long_perh) {

    const RESULT = anno_insol_by_param(ecc, obl, long_perh);
    const
        max = get2dmax(RESULT),
        min = get2dmin(RESULT);

    const contour_plot_data = [{
        z: RESULT,
        x: [...new Array(RESULT[0].length)].map((elem, index) => index), //dateFromDay(2021, index + 1)), // time
        y: [...new Array(180)].map((elem, index) => index - 90), // latitude
        type: 'contour',
        colorscale: colorscale_str,
        line: {
            smoothing: 1
        },
        autocontour: false,
        colorbar: {
            title: 'Insolation',
            tickfont: {
                color: 'black',
            }
        },
        contours: {
            start: min,
            end: max,
            size: 25,
        }
    }];

    var layout = Object.assign({}, contour_plot_layout);
    layout.title.text = `Insolation over time and space`;

    plot_contour({
        divId: 'contour_3',
        data: contour_plot_data,
        layout: layout
    });
}
const
    contour_3_ecc_input = document.getElementById('x_input_contour_3_ecc'),
    contour_3_obl_input = document.getElementById('x_input_contour_3_obl'),
    contour_3_long_input = document.getElementById('x_input_contour_3_long'),
    contour_3_ecc_slide = document.getElementById('contour_3_ecc_slide'),
    contour_3_obl_slide = document.getElementById('contour_3_obl_slide'),
    contour_3_long_slide = document.getElementById('contour_3_long_slide'),
    contour_3_resetBtn = document.getElementById('contour_3_resetBtn');

contour_3_ecc_slide.oninput = () => {
    contour_3_ecc_input.value = contour_3_ecc_slide.value
}
contour_3_obl_slide.oninput = () => {
    contour_3_obl_input.value = contour_3_obl_slide.value
}
contour_3_long_slide.oninput = () => {
    contour_3_long_input.value = contour_3_long_slide.value
}

[contour_3_ecc_slide, contour_3_obl_slide, contour_3_long_slide].forEach((e) => {
    e.onchange = () => {
        plot_contour_3(
            contour_3_ecc_input.value,
            contour_3_obl_input.value,
            contour_3_long_input.value
        );
    }
});

[contour_3_ecc_input, contour_3_obl_input, contour_3_long_input].forEach((e => {
    e.onchange = function () {
        plot_contour_3(
            contour_3_ecc_input.value,
            contour_3_obl_input.value,
            contour_3_long_input.value
        );
    }
}));

contour_3_resetBtn.onclick = () => {
    contour_3_ecc_input.value = contour_3_ecc_slide.value = cplt3_default_values.ecc;
    contour_3_obl_input.value = contour_3_obl_slide.value = cplt3_default_values.obl;
    contour_3_long_input.value = contour_3_long_slide.value = cplt3_default_values.long_perh;
    plot_contour_3(
        cplt3_default_values.ecc,
        cplt3_default_values.obl,
        cplt3_default_values.ecc
    );
}

/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
// EOF
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */