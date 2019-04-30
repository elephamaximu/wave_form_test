var Spectrum = WaveSurfer.create({
    container: '#audio-spectrum',
    progressColor: 'lightseagreen',
    //barWidth: 3.8,
    barHeight: 1.7,
    hideScrollbar: true,
    cursorWidth:0,
});    


window.addEventListener("resize", function(){
    var currentProgress = Spectrum.getCurrentTime() / Spectrum.getDuration();

    Spectrum.empty();
    Spectrum.drawBuffer();

    Spectrum.seekTo(currentProgress);
});
        
Spectrum.load('step-13000-audio.wav');

Spectrum.on('ready', function () {
    var spectrogram = Object.create(WaveSurfer.Spectrogram);
    spectrogram.init({
      wavesurfer: Spectrum,
      container: "#wave-spectrogram",
      fftSamples: 512,
      labels: true
    });
  });

Spectrum.on('audioprocess', function() {
    if(Spectrum.isPlaying()) {
        var totalTime = Spectrum.getDuration(),
        currentTime = Spectrum.getCurrentTime(),
        remainingTime = totalTime - currentTime;
        console.log(Math.round(remainingTime));
        if(Math.round(remainingTime)===0.0) {
            
            Spectrum.stop();
            document.getElementById('btn-pause').value = 'Play';
            document.getElementById('btn-pause').id = 'btn-play';
            document.getElementById('icon').className = 'fa fa-play'
           
        };
        document.getElementById('audio-time').innerText = Math.round(remainingTime).toFixed(1);
    }; 
});


jQuery(function($) {
  $('#btn-play').on('click', function() {
    var $el = $(this)
    
    $el.toggleClass('toggler');
    console.log($el.attr('id'));
    
    $el.attr('id', 'btn-pause')
    $el.attr('value', 'Pause')
    $el.find('i').toggleClass('fa fa-play fa fa-pause');
    });
});

jQuery(function($) {
    $('#btn-pause').on('click', function() {
      var $el = $(this)
      console.log($el.attr('id'));
      
    $el.attr('id', 'btn-play')
    $el.attr('value', 'Play')
    $el.find('i').toggleClass('fa fa-play fa fa-play');
    });
  });