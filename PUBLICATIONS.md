# Publication data format

Edit `publications.yml` for papers and theses. Talks and posters are maintained in
`presentations.yml`, documented in [PRESENTATIONS.md](PRESENTATIONS.md).
The Lua script `scripts/publications.lua` reads both YAML files
at build time and generates `sections/_publications-generated.qmd`; do not edit
that generated file.

The old `publications.xml` is a historical reference; the website no longer
reads it. Make additions and corrections in the YAML files.

This document is repository documentation. It is not included in the explicit
`project.render` or `project.resources` lists in `_quarto.yml`, so Quarto does
not render it or copy it into the website.

## Records and categories

The YAML file is a list of records. Each record requires `id`, `type`, `title`,
and `authors`. Optional fields should be omitted when unused, rather than left
empty.

| Type | Displayed category | Section URL |
| --- | --- | --- |
| `preprint` | Reports | `#pubs-reports` |
| `journal` | Journal papers | `#pubs-journals` |
| `conference` | Conference papers | `#pubs-conferences` |
| `note` | Notes | `#pubs-notes` |
| `thesis` | Theses | `#pubs-theses` |

Categories appear in this order, followed by the integrated Presentations list.
Within a category, entries retain their order in the YAML file; they are not
automatically sorted by year or ID.

The section links at the top of Publications open these category headings;
Presentations uses `#pubs-talks`. These URLs activate the Publications tab
before scrolling. `#publications` opens the top of the tab, retaining the
page header. Individual entries remain addressable as `#publication-ID`.

For a new paper, add a record with a new ID and place it where it should
appear within its category. When a report is accepted or published, update
the existing record: keep its ID and arXiv ID, change its `type` as appropriate,
add `publication`, and update its BibTeX. Do not add a second record for the
published version. Reposition it in the YAML if its display order should change.

## Fields

| Field | Value and interpretation |
| --- | --- |
| `id` | Unique quoted string, such as `'1041'`. Use letters, digits, underscores, or hyphens. IDs must be unique across both YAML files. Keep existing IDs stable: they identify anchors such as `#publication-1041`. |
| `type` | One of the five types above. |
| `title` | Paper title, supporting inline Markdown and `$math$`. Displayed in italics. Do not embed a Markdown link in the title; use the link fields described below. |
| `authors` | Nonempty list of display names in author order. Names are plain text. The renderer supplies commas and the final “and”. Authors are retained in thesis data but omitted from the displayed thesis entry. |
| `arxiv` | Optional quoted arXiv ID, such as `'2610.26748'`, without a URL or version suffix. The current reader accepts only modern `YYMM.NNNN` and `YYMM.NNNNN` IDs, not older identifiers such as `math/0601234`. |
| `link` | Optional title-link override. Can be an absolute URL or a path relative to the website root, such as `downloads/paper.pdf`. |
| `publication` | Optional publication metadata, described below. |
| `code` | Optional code URL, including a GitHub repository, software website, or relative download path. |
| `video` | Optional video URL or mapping, described below. |
| `comment` | Optional Quarto Markdown displayed below the publication information in smaller text. Supports links, paragraphs, italics, and mathematics. |
| `bibtex` | Optional literal block containing the BibTeX citation. Displayed verbatim behind the `bibtex` control. |
| `tags` | Optional list of strings, listed last in the record. Retained as metadata and not displayed. |

All link fields accept absolute URLs or paths relative to the website root.
Place new downloadable files in `downloads/`, which is copied to the rendered
site. The renderer does not check whether URLs or files exist.

### Publication metadata and year

```yaml
publication:
  year: 2026
  venue: SIAM Journal on Optimization
  link: https://epubs.siam.org/doi/10.1137/24M1672742
```

When `publication` is present, `year` and `venue` are required. The year is an
integer. The optional `link` supplies the publisher's link. The optional
`comment` supports inline Markdown and is displayed parenthetically after
the venue. For example, adding `comment: to appear` would produce:

> SIAM Journal on Optimization (to appear), 2026.

Omit `publication.comment` for an ordinary published paper. The comment can
also describe a distinction such as an oral presentation. Do not include the
parentheses yourself. Write the venue using inline Markdown if needed, without
an introductory “in”; for conference papers, use names such as
`Proceedings of COLT`.

Omit the whole `publication` mapping when there is no publication metadata,
as for most reports and notes. Without this mapping, the year is calculated
as 2000 plus the first two digits of `arxiv`: `'1506.00575'` gives 2015.
The information line reads `Notes, 2015.` for `type: note`, and
`Preprint, 2015.` for every other type without publication metadata, including
conference workshop papers. The category still follows `type`.

With `publication` present, its year is used instead, even if the paper was
posted to arXiv in a different year. There is no top-level `year` field.
Every entry needs either `publication` with its year and venue, or an arXiv ID;
a title-link override alone is not sufficient to supply the year.

Theses use `publication.year` and a descriptive `publication.venue`, such as
`Ph.D. thesis` or `Master thesis`. Workshop papers use `type: conference`
and can have a `publication` mapping even when there are no formal proceedings.
Use the workshop name as `publication.venue` and its year as `publication.year`;
no publisher link is required. This displays the workshop information instead
of labeling the paper a preprint.

`publication.comment` annotates the venue parenthetically. The top-level
`comment` is a separate block below the publication information; use it for
awards, corrections, or other contextual text.

### Title links

The renderer uses the first available destination in this order:

1. Top-level `link`.
2. `https://arxiv.org/abs/` followed by `arxiv`.
3. `publication.link`.
4. An unlinked title.

There is no separate arXiv link in the resource row. An explicit top-level
`link` overrides the arXiv destination even when `arxiv` is present.
The publisher's link also appears in the resource row whenever
`publication.link` is provided, even if it is also the title's destination.
Title and resource links open in a new tab; links written inside Markdown
comments use normal Markdown link behavior.

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
the renderer adds them. Video comments support inline Markdown and math.

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
- Quote any text containing a colon followed by a space, for example
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

## Examples

A minimal report, taken from an existing record. Assign a new ID when using
this as a template for a different paper; the optional BibTeX and other
resources can be added afterward.

```yaml
- id: '1041'
  type: preprint
  title: Gradient descent reaches arbitrary accuracy in a fixed number of large steps on multiclass logistic regression
  authors:
    - A.-A. Musat
    - N. Boumal
  arxiv: '2610.26748'
```

For example, a comment can contain both a link and mathematics:

```yaml
comment: |
  Recipient of the [paper prize](https://example.org/prize).

  Correction: replace $\log p(X, r | x_k)$ with $\log p(X, r | x)$.
```

Run `quarto run scripts/publications.lua` to regenerate the publication partial,
or `quarto render` to build the complete website.
