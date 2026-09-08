# raymondxuanhe.github.io

Personal academic website for **Raymond He** (economics, UT Austin), hosted at
<https://raymondxuanhe.github.io> via GitHub Pages. It is a customized fork of the
[academicpages](https://academicpages.github.io/) Jekyll template (derived from Minimal
Mistakes). GitHub Pages builds and deploys automatically on push to `master`.

## Quick start

```bash
bundle install                              # install Ruby/Jekyll deps (delete Gemfile.lock if it errors)
bundle exec jekyll serve -l -H localhost    # serve + live-reload at localhost:4000
docker compose up                           # alternative: serve via the provided Dockerfile
npm run build:js                            # rebuild assets/js/main.min.js after editing assets/js/**
```

`_config.yml` is **not** reloaded by `jekyll serve` — restart the server after editing it.

## Directory structure

```
.
├── _config.yml              # site-wide configuration (Jekyll settings, collections, author info)
├── Gemfile                  # Ruby/Jekyll dependencies
├── package.json             # JS build tooling (npm run build:js)
├── Dockerfile / docker-compose.yaml
├── Project.toml / Manifest.toml   # Julia project for the computational benchmarks
│
├── _pages/                  # standalone site pages (about, research, teaching, publications, CV, benchmarks…)
├── _data/                   # site data: navigation.yml (top nav), authors.yml, ui-text.yml, cv.json
│
│                            # Jekyll collections (configured in _config.yml):
├── _publications/           #   one markdown file per publication
├── _talks/                  #   one markdown file per talk
├── _teaching/               #   one markdown file per teaching entry
├── _portfolio/              #   portfolio items
├── _posts/ , _drafts/       #   blog posts / drafts (mostly unused)
│
│                            # Template internals — change only for structural/theming work:
├── _layouts/                #   page layouts (default, single, cv-layout, talk, splash…)
├── _includes/               #   reusable partials (masthead, footer, head, cv-template…)
├── _sass/                   #   SCSS styles (theme, layout, vendor)
├── assets/                  #   compiled CSS/JS, fonts, webfonts
├── images/                  #   template/site images
│
├── files/                   # all downloadable assets, served at /files/...
│   ├── CV/                  #   CV authored in LaTeX (raymond_he_CV.tex) + compiled PDF & build artifacts
│   ├── Bio/                 #   bio (LaTeX source + PDF)
│   ├── research/            #   papers: jmp/, pay_transparency/, research_statement/
│   ├── research_slides/     #   talk/seminar slide PDFs
│   ├── teaching/            #   jupyter_notebooks/, teaching_evaluations/, teaching_statement/
│   ├── benchmarks/          #   benchmark chart PNGs (cpu, gpu, lightcast)
│   └── julia_files/         #   Julia sources (opt_growth.jl, generate_benchmark_charts.jl)
│
├── markdown_generator/      # optional: generate publication/talk markdown from *.tsv (publications.py, talks.py)
├── scripts/                 # cv_markdown_to_json.py + update_cv_json.sh (inherited template pipeline, not the live CV workflow)
├── talkmap/ , talkmap.py    # generates the map of talk locations
│
├── .github/                 # GitHub Actions / Pages config
└── .devcontainer/           # VS Code dev container config
```

## Content model

Content lives in front-matter markdown/HTML, not code. Add/edit pages in `_pages/` and
collection entries in `_publications/`, `_talks/`, `_teaching/`, `_portfolio/`. Site-wide
config and the top navigation are in `_config.yml` and `_data/navigation.yml`.

## Notable specifics

- **CV** is authored in **LaTeX** at `files/CV/raymond_he_CV.tex`, committed alongside its
  compiled PDF and latexmk build artifacts. The `scripts/` markdown→JSON CV pipeline is
  inherited from the template and is *not* the live workflow.
- **Job-market paper** PDF lives under `files/research/jmp/`.
- **Julia project** (`Project.toml`, `Manifest.toml`, `files/julia_files/`) produces the
  economics computation benchmarks surfaced by `_pages/computational_benchmarks.*`. It is
  separate from the website build — you do not need Julia to serve or deploy the site.

There is no test suite.
