---
layout: archive
title: "Publications"
permalink: /publications/
author_profile: true
---

No publications till now. However, I am currently working on a QM/MM study of the reaction mechanism of CYP enzyme-catalyzed reactions.

<script src="http://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>

<div id="container-01" class="mol-container"></div>

<script>
let element = $('#container-01');
let config = {};
let viewer = $3Dmol.createViewer( element, config );
let pdbUri = '{{ site.baseurl}}/files/pdb/6j83.pdb';
jQuery.ajax( pdbUri, { 
  success: function(data) {
    let v = viewer;
    v.addModel( data, "pdb" );
    v.setStyle({chain: 'A'}, {cartoon: {color: 'spectrum'}});
    v.setStyle({resn: ["HEM", "B9R"]}, {stick: {}});
    v.setStyle({resi: "347"}, {cartoon: {color: 'spectrum'}, stick: {}});
    v.setStyle({resi: "592"}, {stick: {}});
    v.zoomTo();
    v.render();
  },
  error: function(hdr, status, err) {
    console.error( "Failed to load PDB " + pdbUri + ": " + err );
  },
});
</script>

<style>
.mol-container {
  width: 100%;
  height: 400px;
  position: relative;
}
</style>

Structure reference: [6j83](https://www.rcsb.org/structure/6J83).


{% if author.googlescholar %}
  You can also find my articles on <u><a href="{{author.googlescholar}}">my Google Scholar profile</a>.</u>
{% endif %}

{% include base_path %}

{% for post in site.publications reversed %}
  {% include archive-single.html %}
{% endfor %}
