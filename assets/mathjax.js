// Quarto includes MathJax when it detects math in Markdown. Publication math
// lives inside generated HTML, so ensure the library also loads in that case.
// Deferred execution lets us inspect all script tags and avoid a second load.
if (!document.querySelector('script[src*="/mathjax@"]') && !window.MathJax?.startup) {
  const script = document.createElement('script');
  script.src = 'https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-chtml.js';
  document.head.appendChild(script);
}
