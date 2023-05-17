'use strict'
const settings = require('../settings.json')
settings.metadata.info = {}
const dp = require('../data/datapackage.json')
if (dp.measure_info) settings.metadata.measure_info = dp.measure_info
dp.resources.forEach(r => (settings.metadata.info[r.name] = r))
const DataHandler = require('../data_handler.v1.min.js'),
  data = new DataHandler(settings, void 0, {
    block_group: require('../block_group.json'),
    tract: require('../tract.json'),
    county: require('../county.json'),
  })
module.exports.handler = async function (event) {
  return data.export(event.queryStringParameters)
}
