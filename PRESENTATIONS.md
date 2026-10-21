# Presentation data format

Edit `presentations.yml` for talks and posters. It contains 64 presentation
groups and 99 appearances migrated from `publications.xml`. The Lua script
`scripts/publications.lua` reads both YAML files at build time and uses
`scripts/presentations.lua` to validate and render presentations into
`sections/_publications-generated.qmd`. Do not edit that generated file.

The XML is retained for comparison during review; it is no longer read by
the website. No Python or XML parser is required to build the website.

This is repository documentation. It is absent from the explicit
`project.render` and `project.resources` lists in `_quarto.yml`, so Quarto
does not render it or copy it into the website.

## Records and events

The YAML file is a list of presentations. Each record groups related
appearances under one title, with an `events` list describing the individual
appearances. Talks and posters will appear together under **Presentations**.
There is no required `type` or `format` field; a talk is the default.

Keep existing IDs and groupings. Identical titles alone are not a reason to
merge records: the content and associated resources may differ.

Order groups by their most recent event, newest first, and events within
each group newest first. This can move an older presentation upward when
it is given again. For sorting only, an unknown day is treated as day 1;
the displayed date remains month-only. Ties preserve source order.

## Record fields

| Field | Value and interpretation |
| --- | --- |
| `id` | Required unique quoted string, such as `'3036'`. Preserve existing IDs; IDs must remain unique across publications and presentations. |
| `title` | Required title, supporting inline Markdown and `$math$`. Displayed in italics, without a link. |
| `coauthors` | Optional list of collaborators, excluding the presenter. Displayed as “with …”; this is not necessarily a paper's complete author list. |
| `slides` | Optional shared slides URL or relative path. |
| `poster` | Optional shared poster URL or relative path, or `true` if the poster file is unavailable. |
| `video` | Optional shared video URL or mapping, as described below. |
| `notes` | Optional shared lecture-notes URL. Display label: `notes`. |
| `webpage` | Optional shared resource-page URL, distinct from an event announcement. Display label: `webpage`. |
| `abstract_link` | Optional shared URL for an abstract document. Display label: `abstract`. |
| `abstract` | Optional literal block of Quarto Markdown preserving the abstract. Retained as metadata, not currently displayed. |
| `comment` | Optional Quarto Markdown, including links, displayed below the title/coauthors and before shared resources and events, in smaller text. |
| `events` | Required nonempty list of appearances, described below. |
| `format` | Optional internal metadata, at record or event level. Never affects display, ordering, grouping, or resource labels. Omit obvious values. |
| `tags` | Optional list of strings, listed last. Retained as metadata, not displayed. |

Use absolute URLs or paths relative to the website root, such as
`downloads/talk.pdf`. HTML slide decks are valid slide resources too.

## Event fields

| Field | Value and interpretation |
| --- | --- |
| `date` | Required quoted ISO date: `'2026-10-21'`, or `'2019-12'` when the day is unknown. Do not invent missing days. |
| `venue` | Required event or seminar name. Displayed as a link when the event has a `link`. |
| `location` | Optional location text, such as `Lausanne, Switzerland`. |
| `online` | Optional boolean. Use `true` for a remote appearance rather than encoding “Online” as a city. A hosting location may still be retained. |
| `link` | Optional event announcement or programme URL. This links the venue name, never the presentation title. |
| `slides`, `poster`, `video`, `notes`, `webpage`, `abstract_link` | Optional resources specific to this event, taking precedence over the corresponding shared resources. |
| `comment` | Optional Quarto Markdown displayed below this event's line, in smaller text. |
| `format` | Optional internal metadata only, as above. |

Shared resources apply to the group unless an event overrides the same
resource field. An omitted event field inherits the shared value. Use
`false` to explicitly suppress a resource for an event. Shared resources
are displayed once at group level for multi-event groups, so event-level
exceptions do not remove that shared line. If a resource only applies to
some appearances, prefer assigning it directly to those events.

### Posters and videos

A poster URL produces a lowercase `poster` link. When the file is lost, use
`poster: true`: this produces the plain, unlinked word `poster` in precisely
the position where the link would have appeared. No `format: poster` is
needed to obtain that label. Likewise, `format: talk` is normally unnecessary.

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

Each group has an unlinked italic title, with collaborators when provided.
For a group with multiple events, shared resources appear once on a line
below the title and before the events. Event-specific resources appear at
the end of the corresponding event line. Do not repeat shared links beside
every event. For a single-event group, put all its resources on the event
line, regardless of whether they are stored at record or event level.

Resource labels are lowercase and ordered
`slides, poster, video, notes, webpage, abstract`. Missing
resources and their separators are omitted. The plain `poster` label for
a lost file participates in the same ordering.

Format dates as `2026, Oct. 21`; a month-only date becomes `2019, Dec.`.
Use a middle dot (`·`) between the date, venue/location information, and
the resource links. Within the resource list, use commas:

> 2026, Oct. 21 · Seminar name, Lausanne, Switzerland · slides, video

The venue name links to the event announcement. The words `slides`, `poster`,
and `video` link to their respective resources. Titles never inherit a link
from any of these fields. Online appearances should be identified without
repeating legacy “Online” and “(remote)” text in the venue and location.

## Examples

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

Quote IDs and dates. Quote titles containing a colon followed by a space.
Use literal blocks (`abstract: |`, `comment: |`) for multiline Markdown;
backslashes are preserved. Use spaces for indentation and omit unused fields.
For compatibility with the bundled Lua YAML parser, avoid anchors, aliases,
and custom YAML tags.

Run `quarto run scripts/publications.lua` to regenerate the partial, or
`quarto render` to build the website. Presentation titles do not have a
top-level `link` field; the validator rejects it to avoid ambiguous links.

## Migration review

The following are choices made during migration, rather than newly verified
historical facts. They do not prevent rendering, but deserve review.

- **Shared slide versions:** XML group-level slide links remain shared for
  3057, 3047, 3045, 3040, 3037, 3028, 3023, 3014, 3010, and 3002.
  In particular, 3047's February 2023 deck is also associated with a December
  2022 appearance; 3023's Gimli deck covers four appearances, including two
  before Gimli; 3002's November 2010 deck covers appearances through 2013.
  Confirm whether these should remain shared or be assigned to specific events.
  3036's Caltech-only assignment was already confirmed and is preserved.
- **Abstract PDF assignments:** 3003's Imperial PDF is assigned only to the
  Imperial event on 2011-06-14; 3010's Benelux abstract PDF only to the
  Benelux event on 2013-03-26. These assignments follow filenames and event
  names, rather than applying the documents to every appearance.
- **3028's whiteboard comment:** kept at group level, applying to both
  London and MIT. It may belong to only one event.
- **Resource fields and placement:** adopted `notes` for 3053, `webpage` for
  3049, and `abstract_link` for the two PDF abstracts. They follow video in
  the resource order above. Decide whether these labels/order are satisfactory.
- **Abstracts and comments:** prose abstracts are preserved but not shown,
  as before. Comments are now visible, including 3019's prize information
  and 3010's Manopt link. Decide whether abstracts should become accessible
  on the page, and whether this comment placement works well.
- **Event URLs:** former title links to announcements were assigned to the
  corresponding event. This includes 3041's Swiss OR Days link (only its
  2023 appearance), and links extracted from abstracts for 3043 and 3039.
  Trailing announcement URLs in 3045 and 3050's prose abstracts were moved
  to the ETHZ and SIGMA-MODE events respectively. 3047's three comment URLs
  retain the assignments already made in the draft. 3023's embedded Gimli
  link now links its venue; “invited speaker” is an event comment.
- **3043's duplicate video:** retained the event URL and dropped the otherwise
  identical title-link URL ending in `&t=0`.
- **Online appearances:** retain hosting cities/countries where available,
  remove literal “Online” cities and “(remote)” venue suffixes, and display
  `(online)` once. `Internet` is not retained as a country. These locations
  describe the host, not necessarily the presenter's physical location.
- **Incomplete dates:** 4007 retains December 2019 without a day. No date
  was invented; month-only sorting uses day 1 as described above.
- **Grouping:** all original IDs/groupings are preserved, even for identical
  titles. Posters are interleaved chronologically with talks. Records with
  recent repeat appearances, such as 3045, move upward accordingly.
