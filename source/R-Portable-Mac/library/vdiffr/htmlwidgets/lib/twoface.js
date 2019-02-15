// TwoFace – Canvas based before and after image comparison
//
// Copyright (c) 2012-2016 David Hong
// Distributed under the MIT License
//
// http://jsfiddle.net/davidhong/WkM5z/
//
// Similar to jQuery Before & After plugin, except it has no external
// dependency. It is drawn onto a Canvas.

function slide(id, width, height) {

  var canvas = document.createElement('canvas');
  var container = document.getElementById(id);
  var divide = 0.5;

  this.canvas = canvas;
  this.container = container;
  this.ctx = canvas.getContext('2d');
  this.images = [];

  // Event handlers
  canvas.addEventListener('mousemove', handler, false);
  canvas.addEventListener('mousedown', handler, false);
  canvas.addEventListener('mouseup', handler, false);

  var self = this;

  function handler(ev) {
    if (ev.layerX || ev.layerX == 0) { // Firefox
      ev._x = ev.layerX;
      ev._y = ev.layerY;
    } else if (ev.offsetX || ev.offsetX == 0) { // Opera
      ev._x = ev.offsetX;
      ev._y = ev.offsetY;
    }

    var eventHandler = self[ev.type];
    if (typeof eventHandler == 'function') {
      eventHandler.call(self, ev);
    }
  }

  Object.defineProperty(this, 'ready', {
    get: function() {
      return this.images.length >= 2;
    }
  });

  Object.defineProperty(this, 'width', {
    get: function() {
      return width;
    }
  });

  Object.defineProperty(this, 'height', {
    get: function() {
      return height;
    },
    set: function(value) {
      height = value;
    }
  });

  Object.defineProperty(this, 'divide', {
    get: function() {
      return divide;
    },
    set: function(value) {
      if (value > 1) {
        value = (value / 100);
      }

      divide = value;
      this.draw();
    }
  });
}


slide.prototype = {
  add: function(src) {
    function onload(event) {
      this.images.push(img);

      if (this.ready) {
        // Rescale height of Canvas to keep aspect ratio
        var before = this.images[0];
        var ratio = before.naturalHeight / before.naturalWidth;
        this.height = Math.floor(this.width * ratio);

        // Draw canvas into its container
        this.canvas.setAttribute('width', this.width);
        this.canvas.setAttribute('height', this.height);
        this.container.appendChild(this.canvas);

        this.draw();
      }
    }

    var img = createImage(src, onload.bind(this));
  },

  draw: function() {
    if (!this.ready) {
      return;
    }

    var lastIndex = this.images.length - 1,
        before = this.images[lastIndex - 1],
        after = this.images[lastIndex];

    this.drawImages(this.ctx, before, after);
    this.drawHandle(this.ctx);
  },

  drawImages: function(ctx, before, after) {
    ctx.drawImage(after,
      0, 0, after.naturalWidth, after.naturalHeight,
      0, 0, this.width, this.height
    );
    ctx.drawImage(before,
      0, 0, this.divide * before.naturalWidth, before.naturalHeight,
      0, 0, this.divide * this.width, this.height
    );
  },

  drawHandle: function(ctx) {
    var split = this.divide * this.width;

    ctx.fillStyle = "rgb(220, 50, 50)";
    ctx.fillRect(split - 1, 0, 2, this.height);
  },

  mousedown: function(event) {
    var divide = event._x / this.width;
    this.divide = divide;

    this.dragstart = true;
  },

  mousemove: function(event) {
    if (this.dragstart === true) {
      var divide = event._x / this.width;
      this.divide = divide;
    }
  },

  mouseup: function(event) {
    var divide = event._x / this.width;
    this.divide = divide;

    this.dragstart = false;
  }
};


function createImage(src, onload) {
  var img = document.createElement('img');
  img.src = src;

  if (typeof onload == 'function') {
    img.addEventListener('load', onload);
  }

  return img;
}
