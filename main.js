/* <!--
    ########################################
    ## @author Benjamin Thomas Schwertfeger (October 2021)
    ## copyright by Benjamin Thomas Schwertfeger (October 2021)
    ## E-Mail: kontakt@b-schwertfeger.de
    ############
--> */

window.orbital_global = {
    kyears: new Array(0),
    ecc: new Array(0),
    epsilon: new Array(0),
    omega: new Array(0)
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
            processData(data);
        }
    });
});

const Spline = require('cubic-spline');
const test = require("tape");

function processData(allText) {
    let allTextLines = allText.split(/\r\n|\n/);
    let headers = allTextLines[0].split(' ');

    for (let i = 1; i < allTextLines.length; i++) {
        let data = allTextLines[i].split(' ');
        for (let j = 0; j < headers.length; j++) {
            switch (j) {
                case 0:
                    window.orbital_global.kyears.push(-parseInt(data[j]));
                    break;
                case 1:
                    window.orbital_global.ecc.push(parseFloat(data[j])); // eccentricity
                    break;
                case 2:
                    window.orbital_global.epsilon.push(parseFloat(data[j]) + 180);
                    break;
                case 3:
                    window.orbital_global.omega.push(parseFloat(data[j]));
            }
        }
        // if (i == 1) {
        //     break;
        // }
    }
    // const xs = [1, 2, 3, 4, 5];
    // const ys = [9, 3, 6, 2, 4];

    // // new a Spline object
    // const spline = new Spline(xs, ys);

    // // get Y at arbitrary X
    // console.log(spline.at(1.4));

    // // interpolate a line at a higher resolution
    // for (let i = 0; i < 50; i++) {
    //     console.log(spline.at(i * 0.1));
    // }

    //console.log(window.orbital_global.kyears, window.orbital_global.ecc)
    console.log(new Spline([1, 2, 3, 4, 5, 6], [1, 5, 4, 5, 1, 1.6]))
    // console.log(new Spline(window.orbital_global.kyears, window.orbital_global.ecc))
    // console.log(new Spline(window.orbital_global.kyears, window.orbital_global.ecc).at(99))
    // console.log(window.orbital_global.ecc)


    // console.log(window.orbital_global.kyears)
    // console.log(window.orbital_global.ecc)
    // console.log(window.orbital_global.epsilon)
    // console.log(window.orbital_global.omega)


    // function interpolateArray(data, fitCount) {
    //     var linearInterpolate = function (before, after, atPoint) {
    //         return before + (after - before) * atPoint;
    //     };

    //     var newData = new Array();
    //     var springFactor = new Number((data.length - 1) / (fitCount - 1));
    //     newData[0] = data[0]; // for new allocation
    //     for (var i = 1; i < fitCount - 1; i++) {
    //         var tmp = i * springFactor;
    //         var before = new Number(Math.floor(tmp)).toFixed();
    //         var after = new Number(Math.ceil(tmp)).toFixed();
    //         var atPoint = tmp - before;
    //         newData[i] = linearInterpolate(data[before], data[after], atPoint);
    //     }
    //     newData[fitCount - 1] = data[data.length - 1]; // for new allocation
    //     return newData;
    // };

    var ky0 = [1, 2, 3, 4, 5, 6];
    var ecc0 = [1, 5, 4, 5, 1, 1.6];
    // var newArry = interpolateArray(ky0, 18);
    // console.log(new Spline(newArry, [1, 5, 4, 5, 1, 1.6]))
    // console.log(newArry)


}




/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

function getAvg(grades) {
    const total = grades.reduce((acc, c) => acc + c, 0);
    return total / grades.length;
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
        let delta_lambda_m = (day - 80) * 2 * Math.PI / 365.2422; // lambda bei gleich langen Tagen
        let beta = Math.pow((1 - Math.pow(ecc, 2)), 1 / 2);
        let lambda_m0 = (-2) * ((1 / 2 * ecc + 1 / 8 * Math.pow(ecc, 3)) * (1 + beta) * Math.sin(-omega) - 1 / 4 * Math.pow(ecc, 2) * (1 / 2 + beta) * Math.sin(-2 * omega) + 1 / 8 * Math.pow(ecc, 3) * (1 / 3 + beta) * (Math.sin(-3 * omega)));
        let lambda_m = lambda_m0 + delta_lambda_m;

        lambda = lambda_m + (2 * ecc - 1 / 4 * Math.pow(ecc, 3)) * Math.sin(lambda_m - omega) + (5 / 4) * Math.pow(ecc, 2) * Math.sin(2 * (lambda_m - omega)) + (13 / 12) * Math.pow(ecc, 3) * Math.sin(3 * (lambda_m - omega));
    } else if (day_type === 2) { // solar longitude (1-360)
        lambda = day * 2 * Math.PI / 360; // lambda=0 for spring equinox
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
    }
}


function daily_insolation(kyear, lat, day, day_type = 1, fast = true) {
    // CALCULATE DAILY INSOLATION

    // === Get orbital parameters ===
    let temp = null;
    if (fast)
        temp = orbital_parameters_fast(kyear)
    else
        temp = orbital_parameters(kyear)

    let ecc = temp.ecc,
        epsilon = temp.epsilon,
        omega = temp.omega;

    // For output of orbital parameters
    let obliquity = epsilon * 180 / Math.PI,
        long_perh = omega * 180 / Math.PI;

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
    if ((Math.abs(lat) >= Math.PI / 2 - Math.abs(delta)) && (lat * delta > 0)) {
        Ho = Math.PI;
    } else {
        Ho = 0;
    }

    // Insolation: Berger 1978 eq(10)
    let Fsw = So / Math.PI * (1 + ecc * Math.pow(Math.cos(lambda - omega)), 2) / Math.pow((1 - Math.pow(ecc, 2)), 2) * (Ho * Math.sin(lat) * Math.sin(delta) + Math.cos(lat) * Math.cos(delta) * Math.sin(Ho));

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

    let result = new Array(kyear.length());
    for (let year = 0; year < result.length(); year++) {
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
        ecc: window.orbital_global.ecc[kyear * 10],
        epsilon: window.orbital_global.epsilon[kyear * 10],
        omega: window.orbital_global.epsilon[kyear * 10]
    }
}