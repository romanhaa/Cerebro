// Modules to control application life and create native browser window
const { app, BrowserWindow, Menu, MenuItem } = require('electron')

const path = require('path')
const url = require('url')
const port = '8080'
const child = require('child_process')
const MACOS = 'darwin'
const WINDOWS = 'win32'

var killStr = ''
var appPath = path.join(app.getAppPath(), 'app.R' )
var execPath = 'RScript'

// log file
var log = require('electron-log')
log.transports.console.level = true;
log.transports.file.level = true;
log.transports.file.level = 'info'

//
if(process.platform == WINDOWS){
  appPath = appPath.replace(/\\/g, '\\\\')
  execPath = path.join(app.getAppPath(), 'R-Portable-Win', 'bin', 'RScript.exe' )
} else if(process.platform == MACOS){
  var macAbsolutePath = path.join(app.getAppPath(), 'R-Portable-Mac')
  var env_path = macAbsolutePath+((process.env.PATH)?':'+process.env.PATH:'');
  var env_libs_site = macAbsolutePath+'/library'+((process.env.R_LIBS_SITE)?':'+process.env.R_LIBS_SITE:'');
  process.env.PATH = env_path
  process.env.R_LIBS_SITE = env_libs_site
  process.env.NODE_R_HOME = macAbsolutePath
  execPath = path.join(app.getAppPath(), 'R-Portable-Mac', 'bin', 'R' )
} else {
  log.error('Not on Windows or macOS?')
  throw new Error('Not on Windows or macOS?')
}

//
log.info('Environment:\n', process.env)

//
const childProcess = child.spawn(execPath, ['--vanilla', '-e', "shiny::runApp(file.path('"+appPath+"'), port="+port+")"])
childProcess.stdout.on('data', (data) => {
  log.info('stdout:\n', `${data}`)
})
childProcess.stderr.on('data', (data) => {
  log.info('stderr:\n', `${data}`)
})

// Keep a global reference of the window object, if you don't, the window will
// be closed automatically when the JavaScript object is garbage collected.
let mainWindow

// Create the browser window.
function createWindow () {
  log.info('Working directory: '+process.cwd())
  log.info('create-window')

  let loading = new BrowserWindow({show: false, frame: false})

  log.info('showing loading');
  loading.loadURL('data:text/html;charset=utf-8;base64,PGh0bWw+DQo8c3R5bGU+DQpib2R5ew0KICBwYWRkaW5nOiAxZW07DQogIGNvbG9yOiAjNzc3Ow0KICB0ZXh0LWFsaWduOiBjZW50ZXI7DQogIGZvbnQtZmFtaWx5OiAiR2lsbCBzYW5zIiwgc2Fucy1zZXJpZjsNCiAgd2lkdGg6IDgwJTsNCiAgbWFyZ2luOiAwIGF1dG87DQp9DQpoMXsNCiAgbWFyZ2luOiAxZW0gMDsNCiAgYm9yZGVyLWJvdHRvbTogMXB4IGRhc2hlZDsNCiAgcGFkZGluZy1ib3R0b206IDFlbTsNCiAgZm9udC13ZWlnaHQ6IGxpZ2h0ZXI7DQp9DQpwew0KICBmb250LXN0eWxlOiBpdGFsaWM7DQp9DQoubG9hZGVyew0KICBtYXJnaW46IDAgMCAyZW07DQogIGhlaWdodDogMTAwcHg7DQogIHdpZHRoOiAyMCU7DQogIHRleHQtYWxpZ246IGNlbnRlcjsNCiAgcGFkZGluZzogMWVtOw0KICBtYXJnaW46IDAgYXV0byAxZW07DQogIGRpc3BsYXk6IGlubGluZS1ibG9jazsNCiAgdmVydGljYWwtYWxpZ246IHRvcDsNCn0NCg0KLyoNCiAgU2V0IHRoZSBjb2xvciBvZiB0aGUgaWNvbg0KKi8NCnN2ZyBwYXRoLA0Kc3ZnIHJlY3R7DQogIGZpbGw6ICNGRjY3MDA7DQp9DQo8L3N0eWxlPg0KPGJvZHk+PCEtLSAzICAtLT4NCjxkaXYgY2xhc3M9ImxvYWRlciBsb2FkZXItLXN0eWxlMyIgdGl0bGU9IjIiPg0KICA8c3ZnIHZlcnNpb249IjEuMSIgaWQ9ImxvYWRlci0xIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHhtbG5zOnhsaW5rPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hsaW5rIiB4PSIwcHgiIHk9IjBweCINCiAgICAgd2lkdGg9IjgwcHgiIGhlaWdodD0iODBweCIgdmlld0JveD0iMCAwIDUwIDUwIiBzdHlsZT0iZW5hYmxlLWJhY2tncm91bmQ6bmV3IDAgMCA1MCA1MDsiIHhtbDpzcGFjZT0icHJlc2VydmUiPg0KICA8cGF0aCBmaWxsPSIjMDAwIiBkPSJNNDMuOTM1LDI1LjE0NWMwLTEwLjMxOC04LjM2NC0xOC42ODMtMTguNjgzLTE4LjY4M2MtMTAuMzE4LDAtMTguNjgzLDguMzY1LTE4LjY4MywxOC42ODNoNC4wNjhjMC04LjA3MSw2LjU0My0xNC42MTUsMTQuNjE1LTE0LjYxNWM4LjA3MiwwLDE0LjYxNSw2LjU0MywxNC42MTUsMTQuNjE1SDQzLjkzNXoiPg0KICAgIDxhbmltYXRlVHJhbnNmb3JtIGF0dHJpYnV0ZVR5cGU9InhtbCINCiAgICAgIGF0dHJpYnV0ZU5hbWU9InRyYW5zZm9ybSINCiAgICAgIHR5cGU9InJvdGF0ZSINCiAgICAgIGZyb209IjAgMjUgMjUiDQogICAgICB0bz0iMzYwIDI1IDI1Ig0KICAgICAgZHVyPSIwLjZzIg0KICAgICAgcmVwZWF0Q291bnQ9ImluZGVmaW5pdGUiLz4NCiAgICA8L3BhdGg+DQogIDwvc3ZnPg0KPC9kaXY+DQo8L2JvZHk+DQo8L2h0bWw+');

  loading.once('show', () => {
    log.info('show loading')
    mainWindow = new BrowserWindow({webPreferences:{nodeIntegration:false}, show:false, width: 1200, height: 800, title:'Cerebro'})

    mainWindow.webContents.once('dom-ready', () => {
      log.info('mainWindow loaded')
      setTimeout( () => {
        mainWindow.show()
        if(process.platform=MACOS){
          mainWindow.reload()
        }
        loading.hide()
        loading.close()
      }, 2000)
    })

    log.info('Using port '+port)

    // long loading html
    mainWindow.loadURL('http://127.0.0.1:'+port)
    
    mainWindow.webContents.on('did-finish-load', function() {
      log.info('did-finish-load')
    })

    mainWindow.webContents.on('did-start-load', function() {
      log.info('did-start-load')
    })

    mainWindow.webContents.on('did-stop-load', function() {
      log.info('did-stop-load')
    })

    mainWindow.webContents.on('dom-ready', function() {
      log.info('dom-ready')
    })

    // Emitted when the window is closed.
    mainWindow.on('closed', function () {
      log.info('mainWindow.closed()')
      cleanUpApplication()
    })

    // prevent CSS code from appearing in application window title
    mainWindow.on('page-title-updated', function(e) {
      e.preventDefault()
    });

  })

  loading.show()
}


//
function cleanUpApplication(){
  app.quit()
  if(childProcess){
    childProcess.kill();
    if(killStr != '')
      child.execSync(killStr)      
  }
}

//
function setMainMenu() {
  const template = [{
      label: 'Cerebro',
      submenu: [
        { label: 'About', role: 'about' },
        { type: 'separator' },
        { label: 'Quit', accelerator: 'Command+Q', click: function() { app.quit(); }}
      ]}, {
      label: 'Edit',
      submenu: [
        { label: 'Undo', accelerator: 'CmdOrCtrl+Z', role: 'undo' },
        { label: 'Redo', accelerator: 'Shift+CmdOrCtrl+Z', role: 'redo' },
        { type: 'separator' },
        { label: 'Cut', accelerator: 'CmdOrCtrl+X', role: 'cut' },
        { label: 'Copy', accelerator: 'CmdOrCtrl+C', role: 'copy' },
        { label: 'Paste', accelerator: 'CmdOrCtrl+V', role: 'paste' },
        { label: 'Select All', accelerator: 'CmdOrCtrl+A', role: 'selectAll' }
      ]}, {
      label: 'View',
      submenu: [
        { label: 'Refresh', accelerator: 'CmdOrCtrl+R', role: 'reload'}
      ]}, {
      label: 'Window',
      submenu: [
        { label: 'Minimize', accelerator: 'CmdOrCtrl+M', role: 'minimize'},
        { label: 'Close', accelerator: 'CmdOrCtrl+W', role: 'close'}
      ]}
  ];
  Menu.setApplicationMenu(Menu.buildFromTemplate(template));
}



// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.on('ready', function () {
  createWindow();
  setMainMenu();
})

// Quit when all windows are closed.
app.on('window-all-closed', function () {
  log.info('window-all-closed')
  // On OS X it is common for applications and their menu bar
  // to stay active until the user quits explicitly with Cmd + Q
  cleanUpApplication()
})

//
app.on('activate', function () {
  // On OS X it's common to re-create a window in the app when the
  // dock icon is clicked and there are no other windows open.
  if (mainWindow === null) {
    createWindow()
  }
})

// In this file you can include the rest of your app's specific main process
// code. You can also put them in separate files and require them here.
