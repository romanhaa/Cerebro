
document.addEventListener('toggled', function(event) {
  var status = document.getElementById('shiny-toggle-status');
  if (event.now_active === 'before') {
    status.innerHTML = 'Before';
  } else {
    status.innerHTML = 'After';
  }
});
