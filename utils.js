/**
 * © Alfred-Wegener-Institute Bremerhaven, Germany (2022)
 * @link https://awi.de
 * 
 * @author Benjamin Thomas Schwertfeger (2022)
 * @email benjamin.schwertfeger@awi.de
 * @email development@b-schwertfeger.de
 * @link https://b-schwertfeger.de
 * @link https://github.com/btschwertfeger/Orbit-astronomical-theory-of-ice-ages-Website
 **/

module.exports = {
    dateFromDay: (year, day) => {
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
    },
    convolve: (vec1, vec2) => {
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
    },
    get2dmax: (arr) => {
        return Math.max.apply(null, arr.map((row) => {
            return Math.max.apply(Math, row);
        }));
    },
    get2dmin: (arr) => {
        return Math.min.apply(null, arr.map((row) => {
            return Math.min.apply(Math, row);
        }));
    },
    rep: (arr, n) => {
        // repeat array n times
        let output = new Array(n * arr.length);
        for (let i = 0; i < output.length; i++)
            output[i] = arr[i % arr.length];

        return output
    },
    avg: (grades) => {
        const total = grades.reduce((acc, c) => acc + c, 0);
        return total / grades.length;
    },
    arange: (start, end, step = 1) => {
        var
            range = [],
            typeofStart = typeof start,
            typeofEnd = typeof end;

        if (step === 0) throw TypeError('Step cannot be zero.');
        if (typeofStart == 'undefined' || typeofEnd == 'undefined') throw TypeError('Must pass start and end arguments.');
        else if (typeofStart != typeofEnd) throw TypeError('Start and end arguments must be of same type.');

        typeof step == 'undefined' && (step = 1);

        if (end < start) step = -step;
        if (typeofStart == 'number') {
            while (step > 0 ? end >= start : end <= start) {
                range.push(start);
                start += step;
            }
            range.push(start);
        } else throw TypeError('Only string and number types are supported');
        return range;
    },
    divmod: (x, y) => { return [Math.floor(x / y), x % y] },
    argMax: (arr) => {
        return arr.map((x, i) => [x, i]).reduce((r, a) => (a[0] > r[0] ? a : r))[1];
    }
};