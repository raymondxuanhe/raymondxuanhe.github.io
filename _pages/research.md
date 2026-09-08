---
permalink: /research/
title: "Research"
author_profile: true
redirect_from: 
#  - /research/
  - /research.html
---

<style>
/* Light (default) tokens */
.research-list {
  --rl-heading-border: #e6e6e6;
  --rl-hover-bg:       #f2f3f5;
  --rl-hover-shadow:   0 1px 3px rgba(0,0,0,0.06);
  --rl-meta:           #666;
  --rl-note:           #333;
  --rl-abstract:       #2f2f2f;
  --rl-toggle:         #555;
  --rl-btn-blue-fg:    #2f6fb3;
  --rl-btn-blue-hover: #2f6fb3;
  --rl-btn-red-fg:     #c0392b;
  --rl-btn-red-hover:  #c0392b;
  --rl-btn-hover-fg:   #ffffff;
  --rl-wip:            #777;
  margin: 0;
  padding: 0;
}
/* Dark-mode overrides (theme toggles html[data-theme="dark"]) */
html[data-theme="dark"] .research-list {
  --rl-heading-border: #5c5c5c;
  --rl-hover-bg:       #525252;
  --rl-hover-shadow:   0 1px 3px rgba(0,0,0,0.35);
  --rl-meta:           #c2c7ca;
  --rl-note:           #e6e6e6;
  --rl-abstract:       #ececec;
  --rl-toggle:         #cfcfcf;
  --rl-btn-blue-fg:    #7fb6e6;
  --rl-btn-blue-hover: #3f7fc0;
  --rl-btn-red-fg:     #ef8a80;
  --rl-btn-red-hover:  #c0392b;
  --rl-btn-hover-fg:   #ffffff;
  --rl-wip:            #b8bdc0;
}

.research-list h2.section-heading {
  margin: 1.6em 0 0.2em;
  padding-bottom: 0.3em;
  border-bottom: 1px solid var(--rl-heading-border);
  font-size: 1.35em;
}
.research-list h2.section-heading:first-child { margin-top: 0.4em; }

.paper {
  border-radius: 5px;
  padding: 0.65em 0.85em;
  margin: 0.15em 0;
  transition: background-color 0.25s ease, box-shadow 0.25s ease;
  outline: none;
}
.paper:hover,
.paper:focus-within {
  background-color: var(--rl-hover-bg);
  box-shadow: var(--rl-hover-shadow);
}
/* Works-in-progress entries have nothing to expand: no hover highlight */
.paper--plain:hover,
.paper--plain:focus-within {
  background-color: transparent;
  box-shadow: none;
}

.paper-title {
  font-weight: 700;
  font-size: 1.08em;
  line-height: 1.35;
}

.paper-meta {
  color: var(--rl-meta);
  font-size: 0.93em;
  margin-top: 0.15em;
}
.paper-meta a { color: inherit; text-decoration: underline; }
.paper-note {
  color: var(--rl-note);
  font-size: 0.93em;
  margin-top: 0.15em;
}
.paper-note a { color: inherit; text-decoration: underline; }

/* Abstract: collapsed by default, expands on hover/focus (sits ABOVE the controls) */
.paper-abstract {
  max-height: 0;
  opacity: 0;
  overflow: hidden;
  transition: max-height 0.35s ease, opacity 0.3s ease, margin 0.3s ease;
  margin: 0;
  font-size: 0.93em;
  line-height: 1.55;
  color: var(--rl-abstract);
  text-align: justify;
}
.paper:hover .paper-abstract,
.paper:focus-within .paper-abstract {
  max-height: 40em;
  opacity: 1;
  margin: 0.6em 0 0.2em;
}

/* Control row: "Abstract" toggle + always-visible link buttons */
.paper-controls {
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  gap: 0.5em;
  margin-top: 0.55em;
}
.paper-toggle {
  display: inline-flex;
  align-items: center;
  font-size: 0.92em;
  font-weight: 600;
  color: var(--rl-toggle);
  margin-right: 0.15em;
  user-select: none;
}
.paper-toggle::before {
  content: "\25B6";              /* ▶ */
  font-size: 0.7em;
  margin-right: 0.45em;
  transition: transform 0.25s ease;
  display: inline-block;
}
.paper:hover .paper-toggle::before,
.paper:focus-within .paper-toggle::before {
  transform: rotate(90deg);     /* ▶ -> ▼ */
}

.paper-links a {
  display: inline-block;
  font-size: 0.85em;
  font-weight: 600;
  text-decoration: none;
  color: var(--rl-btn-blue-fg);
  border: 1px solid var(--rl-btn-blue-fg);
  border-radius: 4px;
  padding: 0.12em 0.65em;
  transition: background-color 0.2s ease, color 0.2s ease;
}
.paper-links a:hover {
  background-color: var(--rl-btn-blue-hover);
  color: var(--rl-btn-hover-fg);
  text-decoration: none;
}
.paper-links a.btn-paper {
  color: var(--rl-btn-red-fg);
  border-color: var(--rl-btn-red-fg);
}
.paper-links a.btn-paper:hover {
  background-color: var(--rl-btn-red-hover);
  color: var(--rl-btn-hover-fg);
}

.paper-wip { color: var(--rl-wip); }

/* Small screens have no hover: show abstracts, keep toggle pointing down */
@media (max-width: 600px) {
  .paper-abstract { max-height: none; opacity: 1; margin: 0.6em 0 0.2em; }
  .paper-toggle::before { transform: rotate(90deg); }
}
</style>

<div class="research-list" markdown="0">

  <h2 class="section-heading">Working Papers</h2>

  <div class="paper" tabindex="0">
    <div class="paper-title">How Much Can I Make? Cross-Firm Pay Transparency's Effects On the US Labor Market</div>
    <div class="paper-meta">with <a href="https://upmanyu-suryansh.github.io">Suryansh Upmanyu</a> &middot; Submitted</div>
    <div class="paper-abstract">
      Many major US jurisdictions have implemented pay transparency laws that require firms to advertise wage offers in vacancy postings. Using data on the near universe of job postings and representative survey data in a difference-in-differences framework, we find that these laws increased the fraction of postings with wage information by 24.7 percentage points. This translated to average real wage increases of 3.2%-4.4% in Colorado and 0.6% - 1.3% in California and Washington. We further consistently find significant positive effects on real wages for workers who are male, have a college-degree, or are over forty years old.
    </div>
    <div class="paper-controls">
      <span class="paper-toggle">Abstract</span>
      <span class="paper-links">
        <a class="btn-paper" href="/files/research/pay_transparency/pay_transparency_latest_version.pdf">Paper</a>
        <a href="/files/pay_transparency/pay_transparency_online_appendix.pdf">Online Appendix</a>
        <a href="/files/research/pay_transparency/pay_transparency_slides.pdf">Slides</a>
      </span>
    </div>
  </div>

  <h2 class="section-heading">Works in Progress</h2>

  <div class="paper paper--plain">
    <div class="paper-title">Everything In Moderation: The Currency and Maturity Composition of Sovereign Debt</div>
  </div>

  <div class="paper paper--plain">
    <div class="paper-title">Optimal Portfolio Uniqueness in Sovereign Debt Models</div>
    <div class="paper-meta">with <a href="https://www.zachstangebye.com">Zachary Stangebye</a></div>
  </div>

</div>
