// Apply saved preference changes to tags that are already running. Wait for
// Quarto's consent tool to save its cookies before reloading the page.
window.addEventListener('cc_preferencesCenterSavePressed', () => {
  setTimeout(() => window.location.reload(), 0);
});
// Keep Quarto's native consent controls, with wording specific to this site.
document.addEventListener('DOMContentLoaded', () => {
  const title = document.querySelector('.cc-nb-title');
  const text = document.querySelector('.cc-nb-text');
  if (title) title.textContent = 'Optional analytics cookies';
  if (text) {
    text.textContent = 'Google Analytics helps me understand visits to this website. '
      + 'You can accept or reject analytics cookies.';
  }
});

