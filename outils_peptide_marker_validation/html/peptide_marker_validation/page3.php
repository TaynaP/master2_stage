<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
    <link style="text/css" rel="stylesheet" href="/Style/css/bioinfo.css" />
    <link style="text/css" rel="stylesheet" href="/Style/css/page_theme.css" />
    <!--modifier le nom du dossier et du fichier css ci dessous, mettre le même nom que vous lui avez donné -->
    <link href="/peptide_marker_validation/peptide_marker_validation.css" rel="stylesheet" style="text/css"/>
    <script type="text/javascript" src="/libs/jquery-1.11.3.min.js"></script>
    <script type="text/javascript" src="/libs/jquery.history.js"></script>
    <script type="text/javascript" src="/scripts/bioinfo.js"></script>
   <!--modifier le nom du dossier ci dessous -->
   <script type="text/javascript" src="/peptide_marker_validation/js/script.js"></script>   
   <title>Bonsai  :: Bioinformatics Software Server</title>
    <script type="text/javascript">
        var i_am_old_ie = false;
    </script>
    <!--[if LT IE  9]>
    <script type="text/javascript">
            i_am_old_ie = true;
    </script>
    <![endif]-->
  </head>
  <body>
   
   
   <div class="frametitle">
   <h1 id="title">Peptide marker validation</h1>                 
   </div>
   <div id="center_sup">
     <div class="theme-border" style="display:none"></div>
     <div id="link_home" style="display:inline-block"><a href="/" class="text_onglet"><img src="/Style/icon/home_w.png" alt="home_general"/></a></div>
   <div class="tabs" id="menu_central" style="display:inline-block"><?php include("menu_central.txt")?></div>
   </div>
    <div id="main">
 <div id="center">
<br />

<h2> Retrieve result with an ID</h2>

<br />

<div class="formulaire">
<form id="form_id" method="post">
<p>The ID remains valid 24 hours after sequence submission.</p>
<p>

<b>  Enter the ID :</b>  
<input id="run_id" type="text" name="run_id" size="20" />
<input type="hidden" name="command" value="result" />
<input type="submit" value="Go" />
</p>

</form>
</div> <!-- form -->
</div> <!-- center -->
</div> <!-- main -->
<!-- chargement de la librairie php lib.inc -->
   <?php require("../lib.inc")?>
<!-- appel de la fonction footer qui permet d'afficher au bas de la page (nom du logiciel, un lien vers le mail, la date de modif -->
   <?php footer("peptide_marker_validation","peptide_marker_validation", "tayna.pellegri.etu@univ-lille.fr","2022"); ?>

</body>                                        
</html>        


