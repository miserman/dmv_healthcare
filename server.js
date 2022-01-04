'use strict'
const express = require('express'), app = express()
app.use(express.static('docs'))
app.listen(3000, function () {
  console.log('listening on port 3000')
})

