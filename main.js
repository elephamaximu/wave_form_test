jQuery(function($) {
    $('#btn-play').on('click', function() {
      var $el = $(this)
      $el.find('i').toggleClass('fa fa-play fa fa-pause');
      $el.toggleClass('toggler');
      console.log($el.attr('id'));
      console.log($el.attr('value'));
      if($el.attr('id')==='btn-play' & $el.attr('value')==='Play') {
          $el.attr('id', 'btn-pause')
          $el.attr('value', 'Pause')
      } else {
          $el.attr('id', 'btn-play')
          $el.attr('value', 'Play')
      }
    });
  });