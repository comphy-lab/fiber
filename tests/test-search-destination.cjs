const assert = require('node:assert/strict');
const test = require('node:test');

test('search only offers approved destinations', async () => {
  global.window = {
    location: { href: '' },
    searchFuse: {
      search: () => [
        { item: { title: 'Non-lab', url: 'https://example.org/' } },
        { item: { title: 'Non-HTTPS', url: 'http://comphy-lab.org/' } },
        { item: { title: 'Lab', url: 'https://comphy-lab.org/research/' } },
      ],
    },
  };
  require('../.github/assets/js/search-helper.js');
  const results = await window.searchHelper.searchDatabaseForCommandPalette('research');
  assert.deepEqual(results.map(result => result.title), ['Lab']);
  results[0].handler();
  assert.equal(window.location.href, 'https://comphy-lab.org/research/');
});
