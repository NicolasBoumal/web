# Publication data format

Edit `publications.yml` for papers and theses. Talks and posters are maintained in
`presentations.yml`, documented in [PRESENTATIONS.md](PRESENTATIONS.md).
The Lua script `scripts/publications.lua` reads both YAML files
at build time and generates `sections/_publications-generated.qmd`; do not edit
that generated file.

This document is repository documentation. It is not included in the explicit
`project.render` or `project.resources` lists in `_quarto.yml`, so Quarto does
not render it or copy it into the website.

## Records and categories

The YAML file is a list of records. Each record requires `id`, `type`, `title`,
and `authors`. Optional fields should be omitted when unused, rather than left
empty.

| Type | Displayed category |
| --- | --- |
| `preprint` | Reports |
| `journal` | Journal papers |
| `conference` | Conference papers |
| `note` | Notes |
| `thesis` | Theses |

Categories appear in this order, followed by the integrated Presentations list.
Within a category, entries retain their order in the YAML file; they are not
automatically sorted by year or ID.

## Fields

| Field | Value and interpretation |
| --- | --- |
| `id` | Unique quoted string, such as `'1041'`. Use letters, digits, underscores, or hyphens. IDs must be unique across both YAML files. Keep existing IDs stable: they identify anchors such as `#publication-1041`. |
| `type` | One of the five types above. |
| `title` | Paper title, supporting inline Markdown and `$math$`. Displayed in italics. |
| `authors` | Nonempty list of display names in author order. Names are plain text. The renderer supplies commas and the final “and”. Authors are retained in thesis data but omitted from the displayed thesis entry. |
| `arxiv` | Optional quoted arXiv ID, such as `'2610.26748'`, without a URL. The current reader accepts modern `YYMM.NNNN` and `YYMM.NNNNN` IDs. |
| `link` | Optional title-link override. Can be an absolute URL or a path relative to the website root, such as `downloads/paper.pdf`. |
| `publication` | Optional publication metadata, described below. |
| `code` | Optional code URL, including a GitHub repository, software website, or relative download path. |
| `video` | Optional video URL or mapping, described below. |
| `comment` | Optional Quarto Markdown displayed below the publication information in smaller text. Supports links, paragraphs, italics, and mathematics. |
| `bibtex` | Optional literal block containing the BibTeX citation. Displayed verbatim behind the `bibtex` control. |
| `tags` | Optional list of strings, listed last in the record. Retained as metadata and not displayed. |

### Publication metadata and year

```yaml
publication:
  year: 2026
  venue: SIAM Journal on Optimization
  link: https://epubs.siam.org/doi/10.1137/24M1672742
  comment: to appear
```

When `publication` is present, `year` and `venue` are required. The year is an
integer. The optional `link` supplies the publisher's link. The optional
`comment` is displayed parenthetically after the venue:

> SIAM Journal on Optimization (to appear), 2026.

Omit `publication.comment` for an ordinary published paper. There is no
`publication.status` field. Write the venue without an introductory “in”; for
conference papers, names such as `Proceedings of COLT` can be used.

For unpublished papers and notes, omit the whole `publication` mapping. The
year is then inferred from `arxiv`: for the supported IDs, `'1506.00575'`
gives 2015. The information line reads `Preprint, 2015.` or `Notes, 2015.`
according to the type. There is no top-level `year` field. Every entry needs
either `publication` with its year and venue, or an arXiv ID.

Theses use `publication.year` and a descriptive `publication.venue`, such as
`Ph.D. thesis` or `Master thesis`. Workshop papers can retain `type: conference`
without a publication mapping; workshop information can be placed in `comment`.

### Title links

The renderer uses the first available destination in this order:

1. Top-level `link`.
2. `https://arxiv.org/abs/` followed by `arxiv`.
3. `publication.link`.
4. An unlinked title.

There is no separate arXiv link in the resource row. An explicit top-level
`link` overrides the arXiv destination even when `arxiv` is present.

### Video information

Use a URL when no annotation is needed:

```yaml
video: https://youtu.be/-VqBzz3Nedc
```

Use a mapping to include an optional comment:

```yaml
video:
  link: https://www.youtube.com/watch?v=N_bR03acjqc
  comment: 17 min
```

This renders as `video (17 min)`. The comment can also identify the speaker or
provide another short annotation. Omit parentheses from the stored comment;
the renderer adds them. The field is named `comment`, not `duration`.

## Display order

Each entry displays authors, title, publication information, the top-level
comment, and finally the resource controls. Available controls appear in this
exact order, with lowercase labels:

> publisher's link, code, video, bibtex

Code and video links have small icons. Missing resources and their separators
are omitted. BibTeX expands and collapses below its control without a triangle.

Awards, interviews, additional versions, corrections, and poster PDF links can
remain Markdown links in the top-level `comment`. A poster link in a paper's
comment is independent of its presentation record in `presentations.yml`.

## Editing conventions

- Quote IDs and arXiv IDs so YAML treats them as strings.
- Quote titles containing a colon followed by a space, for example
  `title: 'Fast convergence: a result for $C^2$ functions'`.
- Use `comment: |` and `bibtex: |` for multiline text. Backslashes in literal
  blocks are preserved, so write `\log`, not `\\log`.
- Indent author lists under `authors` and use spaces for YAML indentation.
  The bundled Lua parser supports a YAML subset; avoid anchors, aliases, and
  custom YAML tags.
- Keep `tags` last. Other field order does not affect rendering.
- Keep BibTeX consistent with the display metadata. The renderer does not
  derive authors, venue, or year from BibTeX, or regenerate the citation from
  those fields.

For example, a comment can contain both a link and mathematics:

```yaml
comment: |
  Recipient of the [paper prize](https://example.org/prize).

  Correction: replace $\log p(X, r | x_k)$ with $\log p(X, r | x)$.
```

Run `quarto run scripts/publications.lua` to regenerate the publication partial,
or `quarto render` to build the complete website.
