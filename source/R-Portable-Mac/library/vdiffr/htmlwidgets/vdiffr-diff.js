
HTMLWidgets.widget({

  name: 'vdiffr-diff',
  type: 'output',

  initialize: function(el, width, height) {
    return {
      width: width
    }
  },

  renderValue: function(el, x, instance) {
    var div = document.getElementById('vdiffr-diff-img');

    if (div) {
      div.innerHTML = "";
    } else {
      div = document.createElement('div');
      div.id = 'vdiffr-diff-img';
      el.appendChild(div);
    }

    var nloaded = 0;
    var images = [];


    function img_onload(event) {
      nloaded++;

      if (nloaded == 2) {
        var diff = imagediff.diff(images[0], images[1]);
        var canvas = imagediff.createCanvas(instance.width, height);

        var context = canvas.getContext('2d');
        context.putImageData(diff, 0, 0);

        div.appendChild(canvas);
      }
    }

    var height;

    for (var source in x.sources) {
      var img = document.createElement('img');
      img.onload = img_onload;
      img.src = x.sources[source];

      // Rescale height to keep aspect ratio in the Canvas
      if (height == null) {
        var ratio = img.naturalHeight / img.naturalWidth;
        height = Math.floor(instance.width * ratio);
      }

      img.width = instance.width;
      img.height = height;
      images.push(img);
    }
  },

  resize: function(el, width, height, instance) { }

});
