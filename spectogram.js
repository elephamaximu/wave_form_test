/**
 * Calculate FFT - Based on https://github.com/corbanbrook/dsp.js
 */
/* eslint-disable complexity, no-redeclare, no-var, one-var */

var  jet = [[0,0,128],[0,0,132],[0,0,136],[0,0,141],[0,0,145],[0,0,149],[0,0,155],[0,0,159],[0,0,164],[0,0,168],[0,0,173],[0,0,177],[0,0,182],[0,0,187],[0,0,191],[0,0,195],[0,0,200],[0,0,204],[0,0,209],[0,0,214],[0,0,219],[0,0,223],[0,0,227],[0,0,232],[0,0,237],[0,0,241],[0,0,245],[0,0,250],[0,0,254],[0,0,255],[0,0,255],[0,0,255],[0,1,255],[0,5,255],[0,8,255],[0,12,255],[0,17,255],[0,20,255],[0,24,255],[0,28,255],[0,32,255],[0,36,255],[0,41,255],[0,45,255],[0,48,255],[0,53,255],[0,57,255],[0,60,255],[0,65,255],[0,69,255],[0,72,255],[0,77,255],[0,80,255],[0,84,255],[0,88,255],[0,93,255],[0,96,255],[0,100,255],[0,105,255],[0,108,255],[0,112,255],[0,117,255],[0,121,255],[0,124,255],[0,129,255],[0,133,255],[0,136,255],[0,141,255],[0,145,255],[0,148,255],[0,152,255],[0,157,255],[0,160,255],[0,164,255],[0,169,255],[0,172,255],[0,176,255],[0,181,255],[0,185,255],[0,188,255],[0,193,255],[0,197,255],[0,200,255],[0,205,255],[0,209,255],[0,212,255],[0,216,255],[0,221,254],[0,224,251],[0,228,248],[2,233,244],[6,236,241],[9,240,238],[12,245,235],[15,249,232],[19,252,228],[22,255,225],[25,255,222],[28,255,219],[31,255,215],[35,255,212],[38,255,209],[41,255,206],[44,255,202],[48,255,199],[51,255,196],[54,255,193],[57,255,189],[60,255,186],[64,255,183],[67,255,180],[70,255,177],[73,255,173],[77,255,170],[80,255,167],[83,255,164],[86,255,160],[90,255,157],[93,255,154],[96,255,151],[99,255,148],[103,255,144],[106,255,141],[109,255,138],[112,255,135],[115,255,131],[119,255,128],[122,255,125],[125,255,122],[128,255,119],[131,255,115],[135,255,112],[138,255,109],[141,255,106],[144,255,103],[148,255,99],[151,255,96],[154,255,93],[157,255,90],[160,255,86],[164,255,83],[167,255,80],[170,255,77],[173,255,73],[177,255,70],[180,255,67],[183,255,64],[186,255,60],[189,255,57],[193,255,54],[196,255,51],[199,255,48],[202,255,44],[206,255,41],[209,255,38],[212,255,35],[215,255,31],[219,255,28],[222,255,25],[225,255,22],[228,255,19],[232,255,15],[235,255,12],[238,255,9],[241,252,6],[244,248,2],[248,245,0],[251,241,0],[254,237,0],[255,234,0],[255,230,0],[255,226,0],[255,222,0],[255,219,0],[255,215,0],[255,211,0],[255,208,0],[255,204,0],[255,200,0],[255,197,0],[255,193,0],[255,189,0],[255,185,0],[255,182,0],[255,178,0],[255,174,0],[255,171,0],[255,167,0],[255,163,0],[255,159,0],[255,156,0],[255,152,0],[255,148,0],[255,145,0],[255,141,0],[255,137,0],[255,134,0],[255,130,0],[255,126,0],[255,122,0],[255,119,0],[255,115,0],[255,111,0],[255,107,0],[255,104,0],[255,100,0],[255,96,0],[255,93,0],[255,89,0],[255,85,0],[255,82,0],[255,78,0],[255,74,0],[255,71,0],[255,67,0],[255,63,0],[255,59,0],[255,56,0],[255,52,0],[255,48,0],[255,45,0],[255,41,0],[255,37,0],[255,33,0],[255,30,0],[255,26,0],[255,22,0],[254,18,0],[250,15,0],[245,11,0],[241,8,0],[237,4,0],[232,0,0],[227,0,0],[223,0,0],[219,0,0],[214,0,0],[209,0,0],[204,0,0],[200,0,0],[195,0,0],[191,0,0],[187,0,0],[182,0,0],[177,0,0],[173,0,0],[168,0,0],[164,0,0],[159,0,0],[155,0,0],[149,0,0],[145,0,0],[141,0,0],[136,0,0],[132,0,0],[128,0,0]];

const FFT = function(bufferSize, sampleRate, windowFunc, alpha) {
    this.bufferSize = bufferSize;
    this.sampleRate = sampleRate;
    this.bandwidth = (2 / bufferSize) * (sampleRate / 2);

    this.sinTable = new Float32Array(bufferSize);
    this.cosTable = new Float32Array(bufferSize);
    this.windowValues = new Float32Array(bufferSize);
    this.reverseTable = new Uint32Array(bufferSize);

    this.peakBand = 0;
    this.peak = 0;

    var i;
    switch (windowFunc) {
        case 'bartlett':
            for (i = 0; i < bufferSize; i++) {
                this.windowValues[i] =
                    (2 / (bufferSize - 1)) *
                    ((bufferSize - 1) / 2 - Math.abs(i - (bufferSize - 1) / 2));
            }
            break;
        case 'bartlettHann':
            for (i = 0; i < bufferSize; i++) {
                this.windowValues[i] =
                    0.62 -
                    0.48 * Math.abs(i / (bufferSize - 1) - 0.5) -
                    0.38 * Math.cos((Math.PI * 2 * i) / (bufferSize - 1));
            }
            break;
        case 'blackman':
            alpha = alpha || 0.16;
            for (i = 0; i < bufferSize; i++) {
                this.windowValues[i] =
                    (1 - alpha) / 2 -
                    0.5 * Math.cos((Math.PI * 2 * i) / (bufferSize - 1)) +
                    (alpha / 2) *
                        Math.cos((4 * Math.PI * i) / (bufferSize - 1));
            }
            break;
        case 'cosine':
            for (i = 0; i < bufferSize; i++) {
                this.windowValues[i] = Math.cos(
                    (Math.PI * i) / (bufferSize - 1) - Math.PI / 2
                );
            }
            break;
        case 'gauss':
            alpha = alpha || 0.25;
            for (i = 0; i < bufferSize; i++) {
                this.windowValues[i] = Math.pow(
                    Math.E,
                    -0.5 *
                        Math.pow(
                            (i - (bufferSize - 1) / 2) /
                                ((alpha * (bufferSize - 1)) / 2),
                            2
                        )
                );
            }
            break;
        case 'hamming':
            for (i = 0; i < bufferSize; i++) {
                this.windowValues[i] =
                    (0.54 - 0.46) *
                    Math.cos((Math.PI * 2 * i) / (bufferSize - 1));
            }
            break;
        case 'hann':
        case undefined:
            for (i = 0; i < bufferSize; i++) {
                this.windowValues[i] =
                    0.5 * (1 - Math.cos((Math.PI * 2 * i) / (bufferSize - 1)));
            }
            break;
        case 'lanczoz':
            for (i = 0; i < bufferSize; i++) {
                this.windowValues[i] =
                    Math.sin(Math.PI * ((2 * i) / (bufferSize - 1) - 1)) /
                    (Math.PI * ((2 * i) / (bufferSize - 1) - 1));
            }
            break;
        case 'rectangular':
            for (i = 0; i < bufferSize; i++) {
                this.windowValues[i] = 1;
            }
            break;
        case 'triangular':
            for (i = 0; i < bufferSize; i++) {
                this.windowValues[i] =
                    (2 / bufferSize) *
                    (bufferSize / 2 - Math.abs(i - (bufferSize - 1) / 2));
            }
            break;
        default:
            throw Error("No such window function '" + windowFunc + "'");
    }

    var limit = 1;
    var bit = bufferSize >> 1;
    var i;

    while (limit < bufferSize) {
        for (i = 0; i < limit; i++) {
            this.reverseTable[i + limit] = this.reverseTable[i] + bit;
        }

        limit = limit << 1;
        bit = bit >> 1;
    }

    for (i = 0; i < bufferSize; i++) {
        this.sinTable[i] = Math.sin(-Math.PI / i);
        this.cosTable[i] = Math.cos(-Math.PI / i);
    }

    this.calculateSpectrum = function(buffer) {
        // Locally scope variables for speed up
        var bufferSize = this.bufferSize,
            cosTable = this.cosTable,
            sinTable = this.sinTable,
            reverseTable = this.reverseTable,
            real = new Float32Array(bufferSize),
            imag = new Float32Array(bufferSize),
            bSi = 2 / this.bufferSize,
            sqrt = Math.sqrt,
            rval,
            ival,
            mag,
            spectrum = new Float32Array(bufferSize / 2);

        var k = Math.floor(Math.log(bufferSize) / Math.LN2);

        if (Math.pow(2, k) !== bufferSize) {
            throw 'Invalid buffer size, must be a power of 2.';
        }
        if (bufferSize !== buffer.length) {
            throw 'Supplied buffer is not the same size as defined FFT. FFT Size: ' +
                bufferSize +
                ' Buffer Size: ' +
                buffer.length;
        }

        var halfSize = 1,
            phaseShiftStepReal,
            phaseShiftStepImag,
            currentPhaseShiftReal,
            currentPhaseShiftImag,
            off,
            tr,
            ti,
            tmpReal;

        for (var i = 0; i < bufferSize; i++) {
            real[i] =
                buffer[reverseTable[i]] * this.windowValues[reverseTable[i]];
            imag[i] = 0;
        }

        while (halfSize < bufferSize) {
            phaseShiftStepReal = cosTable[halfSize];
            phaseShiftStepImag = sinTable[halfSize];

            currentPhaseShiftReal = 1;
            currentPhaseShiftImag = 0;

            for (var fftStep = 0; fftStep < halfSize; fftStep++) {
                var i = fftStep;

                while (i < bufferSize) {
                    off = i + halfSize;
                    tr =
                        currentPhaseShiftReal * real[off] -
                        currentPhaseShiftImag * imag[off];
                    ti =
                        currentPhaseShiftReal * imag[off] +
                        currentPhaseShiftImag * real[off];

                    real[off] = real[i] - tr;
                    imag[off] = imag[i] - ti;
                    real[i] += tr;
                    imag[i] += ti;

                    i += halfSize << 1;
                }

                tmpReal = currentPhaseShiftReal;
                currentPhaseShiftReal =
                    tmpReal * phaseShiftStepReal -
                    currentPhaseShiftImag * phaseShiftStepImag;
                currentPhaseShiftImag =
                    tmpReal * phaseShiftStepImag +
                    currentPhaseShiftImag * phaseShiftStepReal;
            }

            halfSize = halfSize << 1;
        }

        for (var i = 0, N = bufferSize / 2; i < N; i++) {
            rval = real[i];
            ival = imag[i];
            mag = bSi * sqrt(rval * rval + ival * ival);

            if (mag > this.peak) {
                this.peakBand = i;
                this.peak = mag;
            }
            spectrum[i] = mag;
        }
        return spectrum;
    };
};
/* eslint-enable complexity, no-redeclare, no-var, one-var */

/**
 * @typedef {Object} SpectrogramPluginParams
 * @property {string|HTMLElement} container Selector of element or element in
 * which to render
 * @property {number} fftSamples=512 Number of samples to fetch to FFT. Must be
 * a power of 2.
 * @property {boolean} labels Set to true to display frequency labels.
 * @property {number} noverlap Size of the overlapping window. Must be <
 * fftSamples. Auto deduced from canvas size by default.
 * @property {string} windowFunc='hann' The window function to be used. One of
 * these: `'bartlett'`, `'bartlettHann'`, `'blackman'`, `'cosine'`, `'gauss'`,
 * `'hamming'`, `'hann'`, `'lanczoz'`, `'rectangular'`, `'triangular'`
 * @property {?number} alpha Some window functions have this extra value.
 * (Between 0 and 1)
 * @property {number} pixelRatio=wavesurfer.params.pixelRatio to control the
 * size of the spectrogram in relation with its canvas. 1 = Draw on the whole
 * canvas. 2 = Draw on a quarter (1/2 the length and 1/2 the width)
 * @property {?boolean} deferInit Set to true to manually call
 * `initPlugin('spectrogram')`
 */

/**
 * Render a spectrogram visualisation of the audio.
 *
 * @implements {PluginClass}
 * @extends {Observer}
 * @example
 * // es6
 * import SpectrogramPlugin from 'wavesurfer.spectrogram.js';
 *
 * // commonjs
 * var SpectrogramPlugin = require('wavesurfer.spectrogram.js');
 *
 * // if you are using <script> tags
 * var SpectrogramPlugin = window.WaveSurfer.spectrogram;
 *
 * // ... initialising wavesurfer with the plugin
 * var wavesurfer = WaveSurfer.create({
 *   // wavesurfer options ...
 *   plugins: [
 *     SpectrogramPlugin.create({
 *       // plugin options ...
 *     })
 *   ]
 * });
 */
export default class SpectrogramPlugin {
    /**
     * Spectrogram plugin definition factory
     *
     * This function must be used to create a plugin definition which can be
     * used by wavesurfer to correctly instantiate the plugin.
     *
     * @param  {SpectrogramPluginParams} params Parameters used to initialise the plugin
     * @return {PluginDefinition} An object representing the plugin.
     */
    static create(params) {
        return {
            name: 'spectrogram',
            deferInit: params && params.deferInit ? params.deferInit : false,
            params: params,
            staticProps: {
                FFT: FFT
            },
            instance: SpectrogramPlugin
        };
    }

    constructor(params, ws) {
        this.params = params;
        this.wavesurfer = ws;
        this.util = ws.util;

        this.frequenciesDataUrl = params.frequenciesDataUrl;
        this._onScroll = e => {
            this.updateScroll(e);
        };
        this._onRender = () => {
            this.render();
        };
        this._onWrapperClick = e => {
            this._wrapperClickHandler(e);
        };
        this._onReady = () => {
            const drawer = (this.drawer = ws.drawer);

            this.container =
                'string' == typeof params.container
                    ? document.querySelector(params.container)
                    : params.container;

            if (!this.container) {
                throw Error('No container for WaveSurfer spectrogram');
            }

            this.width = drawer.width;
            this.pixelRatio = this.params.pixelRatio || ws.params.pixelRatio;
            this.fftSamples =
                this.params.fftSamples || ws.params.fftSamples || 512;
            this.height = this.fftSamples / 2;
            this.noverlap = params.noverlap;
            this.windowFunc = params.windowFunc;
            this.alpha = params.alpha;

            this.createWrapper();
            this.createCanvas();
            this.render();

            drawer.wrapper.addEventListener('scroll', this._onScroll);
            ws.on('redraw', this._onRender);
        };
    }

    init() {
        // Check if wavesurfer is ready
        if (this.wavesurfer.isReady) {
            this._onReady();
        } else {
            this.wavesurfer.once('ready', this._onReady);
        }
    }

    destroy() {
        this.unAll();
        this.wavesurfer.un('ready', this._onReady);
        this.wavesurfer.un('redraw', this._onRender);
        this.drawer.wrapper.removeEventListener('scroll', this._onScroll);
        this.wavesurfer = null;
        this.util = null;
        this.params = null;
        if (this.wrapper) {
            this.wrapper.removeEventListener('click', this._onWrapperClick);
            this.wrapper.parentNode.removeChild(this.wrapper);
            this.wrapper = null;
        }
    }

    createWrapper() {
        const prevSpectrogram = this.container.querySelector('spectrogram');
        if (prevSpectrogram) {
            this.container.removeChild(prevSpectrogram);
        }
        const wsParams = this.wavesurfer.params;
        this.wrapper = document.createElement('spectrogram');
        // if labels are active
        if (this.params.labels) {
            const labelsEl = (this.labelsEl = document.createElement('canvas'));
            labelsEl.classList.add('spec-labels');
            this.drawer.style(labelsEl, {
                left: 0,
                position: 'absolute',
                zIndex: 9,
                height: `${this.height / this.pixelRatio}px`,
                width: `${55 / this.pixelRatio}px`
            });
            this.wrapper.appendChild(labelsEl);
            this.loadLabels(
                'rgba(68,68,68,0.5)',
                '12px',
                '10px',
                '',
                '#fff',
                '#f7f7f7',
                'center',
                '#specLabels'
            );
        }

        this.drawer.style(this.wrapper, {
            display: 'block',
            position: 'relative',
            userSelect: 'none',
            webkitUserSelect: 'none',
            height: `${this.height / this.pixelRatio}px`
        });

        if (wsParams.fillParent || wsParams.scrollParent) {
            this.drawer.style(this.wrapper, {
                width: '100%',
                overflowX: 'hidden',
                overflowY: 'hidden'
            });
        }
        this.container.appendChild(this.wrapper);

        this.wrapper.addEventListener('click', this._onWrapperClick);
    }

    _wrapperClickHandler(event) {
        event.preventDefault();

        const relX = 'offsetX' in event ? event.offsetX : event.layerX;
        this.fireEvent('click', relX / this.scrollWidth || 0);
    }

    createCanvas() {
        const canvas = (this.canvas = this.wrapper.appendChild(
            document.createElement('canvas')
        ));

        this.spectrCc = canvas.getContext('2d');

        this.util.style(canvas, {
            position: 'absolute',
            zIndex: 4
        });
    }

    render() {
        this.updateCanvasStyle();

        if (this.frequenciesDataUrl) {
            this.loadFrequenciesData(this.frequenciesDataUrl);
        } else {
            this.getFrequencies(this.drawSpectrogram);
        }
    }

    updateCanvasStyle() {
        const width = Math.round(this.width / this.pixelRatio) + 'px';
        this.canvas.width = this.width;
        this.canvas.height = this.height;
        this.canvas.style.width = width;
    }

    drawSpectrogram(frequenciesData, my) {
        const spectrCc = my.spectrCc;
        const length = my.wavesurfer.backend.getDuration();
        const height = my.height;
        const pixels = my.resample(frequenciesData);
        const heightFactor = my.buffer ? 2 / my.buffer.numberOfChannels : 1;
        let i;
        let j;

        for (i = 0; i < pixels.length; i++) {
            for (j = 0; j < pixels[i].length; j++) {
                var colorValue = 255 - pixels[i][j],
                            rgb = jet[pixels[i][j]];
                        my.spectrCc.fillStyle = `rgb(${rgb[0]}, ${rgb[1]}, ${rgb[2]})`;


                my.spectrCc.fillRect(
                    i,
                    height - j * heightFactor,
                    1,
                    heightFactor
                );
            }
        }
    }

    getFrequencies(callback) {
        const fftSamples = this.fftSamples;
        const buffer = (this.buffer = this.wavesurfer.backend.buffer);
        const channelOne = buffer.getChannelData(0);
        const bufferLength = buffer.length;
        const sampleRate = buffer.sampleRate;
        const frequencies = [];

        if (!buffer) {
            this.fireEvent('error', 'Web Audio buffer is not available');
            return;
        }

        let noverlap = this.noverlap;
        if (!noverlap) {
            const uniqueSamplesPerPx = buffer.length / this.canvas.width;
            noverlap = Math.max(0, Math.round(fftSamples - uniqueSamplesPerPx));
        }

        const fft = new FFT(
            fftSamples,
            sampleRate,
            this.windowFunc,
            this.alpha
        );
        const maxSlicesCount = Math.floor(
            bufferLength / (fftSamples - noverlap)
        );
        let currentOffset = 0;

        while (currentOffset + fftSamples < channelOne.length) {
            const segment = channelOne.slice(
                currentOffset,
                currentOffset + fftSamples
            );
            const spectrum = fft.calculateSpectrum(segment);
            const array = new Uint8Array(fftSamples / 2);
            let j;
            for (j = 0; j < fftSamples / 2; j++) {
                array[j] = Math.max(-255, Math.log10(spectrum[j]) * 45);
            }
            frequencies.push(array);
            currentOffset += fftSamples - noverlap;
        }
        callback(frequencies, this);
    }

    loadFrequenciesData(url) {
        const ajax = this.util.ajax({ url: url });

        ajax.on('success', data =>
            this.drawSpectrogram(JSON.parse(data), this)
        );
        ajax.on('error', e =>
            this.fireEvent('error', 'XHR error: ' + e.target.statusText)
        );

        return ajax;
    }

    freqType(freq) {
        return freq >= 1000 ? (freq / 1000).toFixed(1) : Math.round(freq);
    }

    unitType(freq) {
        return freq >= 1000 ? 'KHz' : 'Hz';
    }

    loadLabels(
        bgFill,
        fontSizeFreq,
        fontSizeUnit,
        fontType,
        textColorFreq,
        textColorUnit,
        textAlign,
        container
    ) {
        const frequenciesHeight = this.height;
        bgFill = bgFill || 'rgba(68,68,68,0)';
        fontSizeFreq = fontSizeFreq || '12px';
        fontSizeUnit = fontSizeUnit || '10px';
        fontType = fontType || 'Helvetica';
        textColorFreq = textColorFreq || '#fff';
        textColorUnit = textColorUnit || '#fff';
        textAlign = textAlign || 'center';
        container = container || '#specLabels';
        const getMaxY = frequenciesHeight || 512;
        const labelIndex = 5 * (getMaxY / 256);
        const freqStart = 0;
        const step =
            (this.wavesurfer.backend.ac.sampleRate / 2 - freqStart) /
            labelIndex;

        const ctx = this.labelsEl.getContext('2d');
        this.labelsEl.height = this.height;
        this.labelsEl.width = 55;

        ctx.fillStyle = bgFill;
        ctx.fillRect(0, 0, 55, getMaxY);
        ctx.fill();
        let i;

        for (i = 0; i <= labelIndex; i++) {
            ctx.textAlign = textAlign;
            ctx.textBaseline = 'middle';

            const freq = freqStart + step * i;
            const index = Math.round(
                (freq / (this.sampleRate / 2)) * this.fftSamples
            );
            const label = this.freqType(freq);
            const units = this.unitType(freq);
            const x = 16;
            const yLabelOffset = 2;

            if (i == 0) {
                ctx.fillStyle = textColorUnit;
                ctx.font = fontSizeUnit + ' ' + fontType;
                ctx.fillText(units, x + 24, getMaxY + i - 10);
                ctx.fillStyle = textColorFreq;
                ctx.font = fontSizeFreq + ' ' + fontType;
                ctx.fillText(label, x, getMaxY + i - 10);
            } else {
                ctx.fillStyle = textColorUnit;
                ctx.font = fontSizeUnit + ' ' + fontType;
                ctx.fillText(units, x + 24, getMaxY - i * 50 + yLabelOffset);
                ctx.fillStyle = textColorFreq;
                ctx.font = fontSizeFreq + ' ' + fontType;
                ctx.fillText(label, x, getMaxY - i * 50 + yLabelOffset);
            }
        }
    }

    updateScroll(e) {
        if (this.wrapper) {
            this.wrapper.scrollLeft = e.target.scrollLeft;
        }
    }

    resample(oldMatrix) {
        const columnsNumber = this.width;
        const newMatrix = [];

        const oldPiece = 1 / oldMatrix.length;
        const newPiece = 1 / columnsNumber;
        let i;

        for (i = 0; i < columnsNumber; i++) {
            const column = new Array(oldMatrix[0].length);
            let j;

            for (j = 0; j < oldMatrix.length; j++) {
                const oldStart = j * oldPiece;
                const oldEnd = oldStart + oldPiece;
                const newStart = i * newPiece;
                const newEnd = newStart + newPiece;

                const overlap =
                    oldEnd <= newStart || newEnd <= oldStart
                        ? 0
                        : Math.min(
                              Math.max(oldEnd, newStart),
                              Math.max(newEnd, oldStart)
                          ) -
                          Math.max(
                              Math.min(oldEnd, newStart),
                              Math.min(newEnd, oldStart)
                          );
                let k;
                /* eslint-disable max-depth */
                if (overlap > 0) {
                    for (k = 0; k < oldMatrix[0].length; k++) {
                        if (column[k] == null) {
                            column[k] = 0;
                        }
                        column[k] += (overlap / newPiece) * oldMatrix[j][k];
                    }
                }
                /* eslint-enable max-depth */
            }

            const intColumn = new Uint8Array(oldMatrix[0].length);
            let m;

            for (m = 0; m < oldMatrix[0].length; m++) {
                intColumn[m] = column[m];
            }

            newMatrix.push(intColumn);
        }

        return newMatrix;
    }
}