// All content is already in the document. Only visibility and URL state change.
(() => {
  const tablist = document.getElementById('homepage-tabs');
  if (!tablist) return;
  const tabs = [...tablist.querySelectorAll('[role="tab"]')];
  const panels = tabs.map(tab => document.getElementById(tab.getAttribute('aria-controls')));

  function selectFromLocation() {
    const hash = location.hash.slice(1);
    // Links to a publication or a heading also reveal its containing panel.
    const target = document.getElementById(hash);
    const panel = target?.closest('.homepage-section') || panels[0];
    document.documentElement.dataset.section = panel.id;
    tabs.forEach((tab, index) => {
      const active = panels[index] === panel;
      tab.setAttribute('aria-selected', String(active));
      tab.tabIndex = active ? 0 : -1;
      tab.parentElement.classList.toggle('active', active);
      panels[index].hidden = !active;
    });
    if (target === panel) window.scrollTo(0, 0);
  }

  tablist.addEventListener('click', event => {
    const tab = event.target.closest('[role="tab"]');
    if (!tab || event.button !== 0 || event.ctrlKey || event.metaKey || event.shiftKey || event.altKey) return;
    event.preventDefault();
    if (location.hash !== tab.hash) history.pushState(null, '', tab.hash);
    selectFromLocation();
  });

  tablist.addEventListener('keydown', event => {
    const index = tabs.indexOf(document.activeElement);
    if (index < 0) return;
    let next;
    if (event.key === 'ArrowRight') next = (index + 1) % tabs.length;
    if (event.key === 'ArrowLeft') next = (index + tabs.length - 1) % tabs.length;
    if (event.key === 'Home') next = 0;
    if (event.key === 'End') next = tabs.length - 1;
    if (next === undefined) return;
    event.preventDefault();
    tabs[next].focus();
    tabs[next].click();
  });

  window.addEventListener('hashchange', selectFromLocation);
  window.addEventListener('popstate', selectFromLocation);
  selectFromLocation();
  // Opening /#publications should retain the page header, just like a tab
  // click. Let genuine links to an individual publication scroll normally.
  window.addEventListener('load', () => {
    if (panels.includes(document.getElementById(location.hash.slice(1)))) {
      window.scrollTo(0, 0);
    }
  });

  const tocToggle = document.getElementById('toc-toggle');
  const toc = document.getElementById('book-toc');
  tocToggle.addEventListener('click', () => {
    const expanded = tocToggle.getAttribute('aria-expanded') !== 'true';
    tocToggle.setAttribute('aria-expanded', String(expanded));
    tocToggle.textContent = expanded ? 'collapse' : 'expand';
    toc.classList.toggle('collapsed', !expanded);
  });
})();
