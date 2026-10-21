# Presentation data format

Edit `presentations.yml` for talks and posters, and `publications.yml` for
papers and theses (see [PUBLICATIONS.md](PUBLICATIONS.md)). The Lua script
`scripts/publications.lua` reads both YAML files at build time and uses
`scripts/presentations.lua` to validate and render presentations into
`sections/_publications-generated.qmd`. Do not edit that generated file.

The old `publications.xml` is a historical reference, not an editing source;
the website no longer reads it. No Python or XML parser is required to build
the website.

This is repository documentation. It is absent from the explicit
`project.render` and `project.resources` lists in `_quarto.yml`, so Quarto
does not render it or copy it into the website.

## Records and events

The YAML file is a list of presentations. Each record groups related
appearances under one title, with an `events` list describing the individual
appearances. Talks and posters appear together under **Presentations**.
There is no `type` field; a talk is the default.

For a new presentation, add a record with a new ID. For another appearance
of an existing presentation, add an event to that record. Keep existing IDs
stable. Identical titles alone are not a reason to merge records: the content
and associated resources may differ.

The renderer automatically sorts groups by their most recent event, newest
first, and events within each group newest first. This can move an older
presentation upward when it is given again. For sorting only, an unknown day is treated as day 1;
the displayed date remains month-only. Ties preserve source order.

## Record fields

| Field | Value and interpretation |
| --- | --- |
| `id` | Required unique quoted string, such as `'3036'`. Use letters, digits, underscores, or hyphens. IDs must be unique across both YAML files and identify anchors such as `#publication-3036`. |
| `title` | Required title, supporting inline Markdown and `$math$`. Displayed in italics. Do not put links in the title; use the event or resource fields. |
| `coauthors` | Optional nonempty list of plain-text collaborator names, excluding the presenter. Displayed as “with …”; this is not necessarily a paper's complete author list. The renderer supplies commas and the final “and”. |
| `slides` | Optional shared slides URL or relative path. |
| `poster` | Optional shared poster URL or relative path, or `true` if the poster file is unavailable. |
| `video` | Optional shared video URL or mapping, as described below. |
| `notes` | Optional shared lecture-notes URL or relative path. Display label: `notes`. |
| `abstract` | Optional literal block of Quarto Markdown preserving the abstract. Retained as metadata, not currently displayed. |
| `comment` | Optional Quarto Markdown, including links, displayed below the title, its resource links, and any coauthors, before the events, in smaller text. |
| `events` | Required nonempty list of appearances, described below. |
| `meta` | Optional plain-text internal note, at record or event level. Never displayed. Unlike `comment`, it does not appear on the page. |
| `tags` | Optional list of strings, listed last. Retained as metadata, not displayed. |

Use absolute URLs or paths relative to the website root, such as
`downloads/talk.pdf`. HTML slide decks are valid slide resources too.
Place new downloadable files in `downloads/`; that folder is copied to the
rendered site. The renderer does not check whether URLs or files exist.

## Event fields

| Field | Value and interpretation |
| --- | --- |
| `date` | Required quoted ISO date: `'2026-10-21'`, or `'2019-12'` when the day is unknown. Do not invent missing days. |
| `venue` | Required event or seminar name, supporting inline Markdown and math. Put the announcement URL in `link` rather than embedding a link in this text. |
| `location` | Optional plain-text location, such as `Lausanne, Switzerland`. |
| `online` | Optional boolean. `true` appends `(online)` after the venue/location; omit it for an in-person event. A hosting location may still be retained. |
| `link` | Optional event announcement, programme, or tutorial-page URL or relative path. This links the venue name, never the presentation title. |
| `slides`, `poster`, `video`, `notes` | Optional resources specific to this event, shown at the end of its event line. |
| `comment` | Optional Quarto Markdown displayed below this event's line, in smaller text. |
| `abstract` | Optional Quarto Markdown retained as metadata, not displayed, as for a record-level abstract. |
| `meta` | Optional undisplayed internal note, as above. |

Record-level resources appear once immediately after the title on the same
line, before any coauthor line, record comment, or events. Event-level resources appear only beside their
own events. This rule is the same whether a record has one event or several.
There is no resource inheritance or replacement between the two levels:
an event resource does not hide a record-level link. Omit unused fields;
`false` produces no resource at that level and does not affect the other level.

If a resource applies to only some appearances, put it on those events
instead of at record level. Omit unused record-level resource fields.

### Posters and videos

A poster URL produces a lowercase `poster` link. When the file is lost, use
`poster: true`: this produces the plain, unlinked word `poster` in precisely
the position where the link would have appeared. The `poster` resource itself
identifies a poster presentation; no separate classification field is needed.

Video information follows the publication format. Use a URL alone when no
annotation is needed, or a mapping for an optional parenthetical comment:

```yaml
video:
  link: https://youtu.be/-VqBzz3Nedc
  comment: 45 min
```

This displays `video (45 min)`. The comment may be something other than a
duration; omit parentheses from the stored text.

## Display rules

Each group has an unlinked italic title followed by top-level resources,
separated by middle dots, even for a single-event presentation:

> *My title in italics* · slides · notes

Collaborators, when provided, appear on the next line, followed by any
record comment and the event list. Event-specific resources appear at the end of
the corresponding event line. Do not repeat shared links beside every event.

Resource labels are lowercase and ordered
`slides, poster, video, notes`. Missing
resources and their separators are omitted. The plain `poster` label for
a lost file participates in the same ordering.

Format dates as `2026, Oct. 21`; a month-only date becomes `2019, Dec.`.
Use a middle dot (`·`) between the date, venue/location information, and
the resource links. Within an event's resource list, use commas:

> 2026, Oct. 21 · Seminar name, Lausanne, Switzerland · slides, video

The venue name is linked only when `event.link` is present. The words
`slides`, `poster`, `video`, and `notes` link to their respective resources,
except for the plain `poster` label when its file is lost. Resource and
venue links open in a new tab. Online appearances show `(online)` once;
do not repeat “Online” or “(remote)” in the venue or location. Any retained
online location describes the host, not necessarily the presenter's location.

## Examples

These examples are existing records. Assign a new ID when using one as a
template for a new presentation.

Different resources for different appearances of the same presentation:

```yaml
- id: '3036'
  title: What is the role of curvature in complexity of nonconvex optimization on manifolds?
  events:
    - date: '2020-02-21'
      venue: Caltech EE Systems Seminar
      location: Pasadena, CA, USA
      slides: downloads/CalTech_2020_curvature_optimization.pdf
    - date: '2019-10-31'
      venue: CMO-BIRS Workshop on Computational Harmonic Analysis and Data Science
      location: Oaxaca, Mexico
      video: http://www.birs.ca/events/2019/5-day-workshops/19w5061/videos/watch/201910311103-Boumal.html
```

A poster whose file is unavailable:

```yaml
- id: '4004'
  title: The non-convex Burer-Monteiro approach works on smooth semidefinite programs
  coauthors:
    - V. Voroninski
    - A. Bandeira
  events:
    - date: '2016-12-05'
      venue: Annual Conference on Neural Information Processing Systems (NIPS)
      location: Barcelona, Spain
      poster: true
  tags: [manopt, optim, sdplr]
```

## Editing

Quote IDs and dates. Quote any text containing a colon followed by a space,
not just titles. Keep booleans (`true`, `false`) unquoted.
Use literal blocks (`abstract: |`, `comment: |`) for multiline Markdown;
backslashes are preserved. Use spaces for indentation and omit unused fields.
Keep lists indented under their field names, as in the examples. Keep `tags`
last; other field order does not affect display. For compatibility with the
bundled Lua YAML parser, avoid anchors, aliases, and custom YAML tags.

Use `comment` for visible contextual text, `meta` for undisplayed internal
notes, and `abstract` for undisplayed abstract prose. These are text fields;
only `video` accepts a mapping with its own `link` and `comment`. Do not use
the paper-specific fields from `publications.yml` in presentation records.

Run `quarto run scripts/publications.lua` to regenerate the partial, or
`quarto render` to build the website. Presentation titles do not have a
top-level `link` field; the validator rejects it to avoid ambiguous links.

## Notes on existing records

Some migrated presentations use one slide deck for several appearances.
The shared assignments for 3047, 3023, and 3002 have been confirmed; other
shared assignments were retained from the XML. When adding or revising an
appearance, check whether the shared deck is appropriate, and use event-level
resources when different versions are needed.

For 3036, the Caltech slides belong only to the Caltech event. For 3003,
the Imperial PDF is slides, not an abstract. Presentation 3049 uses its
local tutorial page as the event's `link`; 3053 uses `notes` for lecture notes.

The date of 4007 was verified from the
[official NeurIPS schedule](https://neurips.cc/Conferences/2019/Schedule?showEvent=13718).
