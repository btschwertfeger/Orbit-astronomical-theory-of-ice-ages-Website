(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){
/* <!--
    ########################################
    ## @author Benjamin Thomas Schwertfeger (October 2021)
    ## copyright by Benjamin Thomas Schwertfeger (October 2021)
    ## E-Mail: kontakt@b-schwertfeger.de
    ############
--> */


/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
// HELPER FUNCTIONS
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */


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


/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */
// CALCULATION
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

const Spline = require("cubic-spline");

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
        // if (i == 5000) {
        //     break;
        // }
    }
    /* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ---------- ----- ----- */
    // BEGIN TESTING
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
    // END TESTING
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

    plotDefault()

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
// PLOTTING
/* ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- */

// Plot of insolation for 5.000 kys
// #Validation
// #plot(june.65N)

function plotDefault() {
    // ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    const labels = [...new Array(window.TIMESTEPS)].map((elem, index) => -index)
    const values = new Array();

    for (let year = 0; year < 5000; year++) {
        values.push(daily_insolation(year, 65, 172, 1, false).Fsw); // false or true for fast and not fast
    }

    const data = {
        labels: labels,
        datasets: [{
            data: values,
            fill: false,
            borderColor: 'rgb(255, 0, 0)',
            tension: 0.1
        }]
    };

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    document.getElementById('orbital_line_plot').remove();
    document.getElementById('orbital_line_plot_container').innerHTML = '<canvas id="orbital_line_plot"></canvas>';
    let ctx = document.getElementById('orbital_line_plot');
    const config1 = {
        type: 'line',
        data: data,
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                title: {
                    display: true,
                    text: "Insolation for 5.000 kys",
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
                        text: 't',
                        font: {
                            family: window.font_famliy,
                            size: 16,
                        },
                    },
                },
                y: {
                    display: true,
                    title: {
                        display: true,
                        text: 'y',
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
    window.orbital_line_plot = new Chart(ctx, config1);

    //     june .65 N.new < -vector()
    //     for (i in 1: 5000) {
    //         june .65 N.new[i] < -daily_insolation(kyear = i, lat = 65, day = 172) $Fsw
    //     }

    //     plot(-(1: 5000), june .65 N.new, col = "red", type = 'l')
}
},{"cubic-spline":2}],2:[function(require,module,exports){
module.exports = class Spline {
  constructor(xs, ys) {
    this.xs = xs;
    this.ys = ys;
    this.ks = this.getNaturalKs(new Float64Array(this.xs.length));
  }

  getNaturalKs(ks) {
    const n = this.xs.length - 1;
    const A = zerosMat(n + 1, n + 2);

    for (
      let i = 1;
      i < n;
      i++ // rows
    ) {
      A[i][i - 1] = 1 / (this.xs[i] - this.xs[i - 1]);
      A[i][i] =
        2 *
        (1 / (this.xs[i] - this.xs[i - 1]) + 1 / (this.xs[i + 1] - this.xs[i]));
      A[i][i + 1] = 1 / (this.xs[i + 1] - this.xs[i]);
      A[i][n + 1] =
        3 *
        ((this.ys[i] - this.ys[i - 1]) /
          ((this.xs[i] - this.xs[i - 1]) * (this.xs[i] - this.xs[i - 1])) +
          (this.ys[i + 1] - this.ys[i]) /
            ((this.xs[i + 1] - this.xs[i]) * (this.xs[i + 1] - this.xs[i])));
    }

    A[0][0] = 2 / (this.xs[1] - this.xs[0]);
    A[0][1] = 1 / (this.xs[1] - this.xs[0]);
    A[0][n + 1] =
      (3 * (this.ys[1] - this.ys[0])) /
      ((this.xs[1] - this.xs[0]) * (this.xs[1] - this.xs[0]));

    A[n][n - 1] = 1 / (this.xs[n] - this.xs[n - 1]);
    A[n][n] = 2 / (this.xs[n] - this.xs[n - 1]);
    A[n][n + 1] =
      (3 * (this.ys[n] - this.ys[n - 1])) /
      ((this.xs[n] - this.xs[n - 1]) * (this.xs[n] - this.xs[n - 1]));

    return solve(A, ks);
  }

  /**
   * inspired by https://stackoverflow.com/a/40850313/4417327
   */
  getIndexBefore(target) {
    let low = 0;
    let high = this.xs.length;
    let mid = 0;
    while (low < high) {
      mid = Math.floor((low + high) / 2);
      if (this.xs[mid] < target && mid !== low) {
        low = mid;
      } else if (this.xs[mid] >= target && mid !== high) {
        high = mid;
      } else {
        high = low;
      }
    }
    return low + 1;
  }

  at(x) {
    let i = this.getIndexBefore(x);
    const t = (x - this.xs[i - 1]) / (this.xs[i] - this.xs[i - 1]);
    const a =
      this.ks[i - 1] * (this.xs[i] - this.xs[i - 1]) -
      (this.ys[i] - this.ys[i - 1]);
    const b =
      -this.ks[i] * (this.xs[i] - this.xs[i - 1]) +
      (this.ys[i] - this.ys[i - 1]);
    const q =
      (1 - t) * this.ys[i - 1] +
      t * this.ys[i] +
      t * (1 - t) * (a * (1 - t) + b * t);
    return q;
  }
};

function solve(A, ks) {
  const m = A.length;
  let h = 0;
  let k = 0;
  while (h < m && k <= m) {
    let i_max = 0;
    let max = -Infinity;
    for (let i = h; i < m; i++) {
      const v = Math.abs(A[i][k]);
      if (v > max) {
        i_max = i;
        max = v;
      }
    }

    if (A[i_max][k] === 0) {
      k++;
    } else {
      swapRows(A, h, i_max);
      for (let i = h + 1; i < m; i++) {
        const f = A[i][k] / A[h][k];
        A[i][k] = 0;
        for (let j = k + 1; j <= m; j++) A[i][j] -= A[h][j] * f;
      }
      h++;
      k++;
    }
  }

  for (
    let i = m - 1;
    i >= 0;
    i-- // rows = columns
  ) {
    var v = 0;
    if (A[i][i]) {
      v = A[i][m] / A[i][i];
    }
    ks[i] = v;
    for (
      let j = i - 1;
      j >= 0;
      j-- // rows
    ) {
      A[j][m] -= A[j][i] * v;
      A[j][i] = 0;
    }
  }
  return ks;
}

function zerosMat(r, c) {
  const A = [];
  for (let i = 0; i < r; i++) A.push(new Float64Array(c));
  return A;
}

function swapRows(m, k, l) {
  let p = m[k];
  m[k] = m[l];
  m[l] = p;
}

},{}]},{},[1]);
