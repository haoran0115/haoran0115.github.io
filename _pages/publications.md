---
layout: archive
title: "Publications"
permalink: /publications/
author_profile: true
---

No publications till now. However, I am currently working on a QM/MM study of reaction mechanism of CYP enzyme-catalyzed reactions.

{% if author.googlescholar %}
  You can also find my articles on <u><a href="{{author.googlescholar}}">my Google Scholar profile</a>.</u>
{% endif %}

{% include base_path %}

{% for post in site.publications reversed %}
  {% include archive-single.html %}
{% endfor %}
