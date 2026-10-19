// Enhance native details so BibTeX slides open/closed without moving its link.
// Native mouse and keyboard toggling still work when JavaScript is disabled.
document.querySelectorAll('.publication details').forEach(details => {
  const summary = details.querySelector('summary');
  const content = details.querySelector('.bibtex-content');
  let expanded = details.open;
  let animation;

  summary.addEventListener('click', event => {
    event.preventDefault();
    const startHeight = details.open ? content.getBoundingClientRect().height : 0;
    expanded = !expanded;
    animation?.cancel();

    if (window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
      details.open = expanded;
      content.classList.remove('is-animating');
      animation = undefined;
      return;
    }

    // Keep details open throughout closing so the browser can show each frame.
    details.open = true;
    const endHeight = expanded ? content.getBoundingClientRect().height : 0;
    content.classList.add('is-animating');
    animation = content.animate(
      [{ height: `${startHeight}px` }, { height: `${endHeight}px` }],
      { duration: 400, easing: 'ease-in-out' }
    );
    animation.onfinish = () => {
      details.open = expanded;
      content.classList.remove('is-animating');
      animation = undefined;
    };
  });
});
