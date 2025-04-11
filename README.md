# My personal website:

[marcopxt.github.io](https://filipe-manoel.github.io/)


# How to modify

1. Change '_data/navigation.yml' for add/delete tabs
2. If you add one, create a new .rd file and replace it.
3. Modify the top header of the .rd file in the field /permalink/ with the name of the tab that you added in Step 1. This will ensure that you will have a connection to the website.
3.1 It should look like this:  

---
layout: archive  
title: "Research topics"   #Title on the page  
permalink: /projects/  #Conect with step 1  
author_profile: true  
---  

4. You should add the following to the .md file, after the box above:

{% include base_path %}
