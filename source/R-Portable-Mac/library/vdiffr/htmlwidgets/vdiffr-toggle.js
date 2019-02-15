
HTMLWidgets.widget({

  name: 'vdiffr-toggle',
  type: 'output',

  initialize: function(el, width, height) { },

  renderValue: function(el, x, instance) {
    var switched = false;
    var img = document.getElementById('vdiffr-toggle-img');

    if (!img) {
      img = document.createElement('img');
      img.id = 'vdiffr-toggle-img';
    }

    img.src = x.files['before'];
    img.className = 'vdiffr-before'

    // Just easier to use a document-wide event
    var toggle_event = document.createEvent('Event');
    toggle_event.initEvent('toggled');
    toggle_event.now_active = 'before';
    document.dispatchEvent(toggle_event);

    img.onclick = function () {
      if (switched) {
        img.src = x.files['before'];
        img.className = 'vdiffr-before';
        toggle_event.now_active = 'before';
        document.dispatchEvent(toggle_event);
      } else {
        img.src = x.files['after'];
        img.className = 'vdiffr-after';
        toggle_event.now_active = 'after';
        document.dispatchEvent(toggle_event);
      }

      switched = !switched;
    };

    el.appendChild(img);
  },

  resize: function(el, width, height, instance) { }

});
