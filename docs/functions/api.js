'use strict'
const DataHandler = require('../data_handler.min.js'),
  data = new DataHandler(require('../settings.json'), void 0, {
    block_group: require('../block_group.json'),
    tract: require('../tract.json'),
    county: require('../county.json'),
  })
module.exports.handler = async function (event) {
  return data.export(event.queryStringParameters)
}
