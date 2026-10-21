// All content is already in the document. Only visibility and URL state change.
(() => {
  const tablist = document.getElementById('homepage-tabs');
  if (!tablist) return;
  const tabs = [...tablist.querySelectorAll('[role="tab"]')];
  const panels = tabs.map(tab => document.getElementById(tab.getAttribute('aria-controls')));
  // Restore the section named by the URL, rather than an old scroll position.
  if ('scrollRestoration' in history) history.scrollRestoration = 'manual';

  function selectFromLocation(scrollBehavior = 'auto') {
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
    else if (target) target.scrollIntoView({ block: 'start', behavior: scrollBehavior });
  }

  tablist.addEventListener('click', event => {
    const tab = event.target.closest('[role="tab"]');
    if (!tab || event.button !== 0 || event.ctrlKey || event.metaKey || event.shiftKey || event.altKey) return;
    event.preventDefault();
    if (location.hash !== tab.hash) history.pushState(null, '', tab.hash);
    selectFromLocation();
  });

  document.querySelector('.publication-sections')?.addEventListener('click', event => {
    const link = event.target.closest('a[href^="#pubs-"]');
    if (!link || event.button !== 0 || event.ctrlKey || event.metaKey || event.shiftKey || event.altKey) return;
    event.preventDefault();
    if (location.hash !== link.hash) history.pushState(null, '', link.hash);
    const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
    selectFromLocation(reducedMotion ? 'auto' : 'smooth');
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

  window.addEventListener('hashchange', () => selectFromLocation());
  window.addEventListener('popstate', () => selectFromLocation());
  selectFromLocation();
  // Reapply after the browser's initial fragment scroll. Main tab URLs retain
  // the page header; category and individual-entry URLs scroll to their target.
  window.addEventListener('load', () => selectFromLocation());

  const tocToggle = document.getElementById('toc-toggle');
  const toc = document.getElementById('book-toc');
  tocToggle.addEventListener('click', () => {
    const expanded = tocToggle.getAttribute('aria-expanded') !== 'true';
    tocToggle.setAttribute('aria-expanded', String(expanded));
    tocToggle.textContent = expanded ? 'collapse' : 'expand';
    toc.classList.toggle('collapsed', !expanded);
  });
})();
