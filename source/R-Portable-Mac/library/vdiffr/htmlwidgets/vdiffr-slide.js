
HTMLWidgets.widget({

  name: "vdiffr-slide",
  type: "output",

  initialize: function(el, width, height) {
    var slide_el = document.getElementById("vdiffr-slide");

    if (!slide_el) {
      var slide_div = document.createElement('div');
      slide_div.id = "vdiffr-slide";

      slide_fig = document.createElement('figure');
      slide_fig.appendChild(slide_div);

      el.appendChild(slide_fig);
    }

    var slider = new slide('vdiffr-slide', width, height);

    return {
      slider: slider
    };
  },

  renderValue: function(el, x, instance) {
    var slider = instance.slider;

    for (var name in x.sources)
      slider.add(x.sources[name]);
  },

  resize: function(el, width, height, instance) { }

});
