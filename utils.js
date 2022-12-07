/**
 * © Alfred-Wegener-Institute Bremerhaven, Germany (2022)
 * @link https://awi.de
 * 
 * @author Benjamin Thomas Schwertfeger (2022)
 * @email benjamin.schwertfeger@awi.de
 * @email development@b-schwertfeger.de
 * @link https://b-schwertfeger.de
 * 
 **/

export function dateFromDay(year, day) {
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

export const convolve = (vec1, vec2) => {
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

export function get2dmax(arr) {
    return Math.max.apply(null, arr.map(function (row) {
        return Math.max.apply(Math, row);
    }));
}

export function get2dmin(arr) {
    return Math.min.apply(null, arr.map(function (row) {
        return Math.min.apply(Math, row);
    }));
}

export function rep(arr, n) {
    // repeat array n times
    let output = new Array(n * arr.length);
    for (let i = 0; i < output.length; i++) {
        output[i] = arr[i % arr.length];
    }
    return output
}

export function avg(grades) {
    const total = grades.reduce((acc, c) => acc + c, 0);
    return total / grades.length;
}