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
   <div class="tabs" id="menu_central" style="display:inline-block"><?php include("menu_central.txt")?></div> <!--php pas necessaire-->
   </div>
    <div id="main">
     <div id="center">
<form id="formulaire" method="post">
  <div class="formulaire">
    <table class="vide">
      <tr>
	<td class="label">Enter a <b>name</b> for the sequences 
	  <i>(optional) </i> : 
	  <input id="seq_name" type="text" name="seq_name" size="20" />
	</td>
      </tr>
    </table>
  </div>
  <!-- Partie à adapter en fonction de votre programme, essayer d'avoir des noms (name) et id (id) en adéquation avec ce que vous demandez  -->
  <div class="formulaire">
    <table class="vide">
      <tr>
	<td class="label"><b>Enter a list of putative markers</b>
	</td>
      </tr>
      <tr>
	<td>
	  <textarea id="paste_pep" name="peptides" rows="15" cols="40"></textarea>
	</td>
      </tr>
  <tr>
  <td>
    Format:
    <ul>
      <li>NAME of the peptide : =</li>
      <li>OR : [ , ]</li>
      <li>AND : , </li>
    </ul>
  </td>
      </tr>
      <tr>
	<td>
	</td>
      </tr>
      <tr>
	<td><input id="ex1" type="submit" name="example" value="Example" /></td>
      </tr>
    </table>
  </div>

  <div class="formulaire">
    <table class="vide">
      <tr>
	<td><b>Enter target sequences (multiFASTA format)</b>
	</td>
      </tr>
      <tr>
	<td class="label">Protein sequences
	</td>
      </tr>
      <tr>
	<td>
	  <textarea id="paste_prot" name="proteins" rows="15" cols="40"></textarea>
	</td>
      </tr>
      <tr>
	<td class="label"> 
	  <B>upload</B> a file
	  <input type="file" name="file_prot" value="file_prot"></input>
	</td>
      </tr> 

      <tr>
  <td class="label">and/or RNA sequences (messenger RNA)
	</td>
      </tr>
      <tr>
	<td>
	  <textarea id="paste_rna" name="rna" rows="15" cols="40"></textarea>
	</td>
      </tr>
      <tr>
	<td class="label"> 
	  <B>upload</B> a file
	  <input type="file" name="file_rna" value="file_rna"></input>
	</td>
      </tr> 

      <tr>
  <td class="label">and/or DNA sequences (genomic DNA)
	</td>
      </tr>
      <tr>
	<td>
	  <textarea id="paste_dna" name="dna" rows="15" cols="40"></textarea>
	</td>
      </tr>
      <tr>
	<td class="label"> 
	  <B>upload</B> a file
	  <input type="file" name="file_dna" value="file_dna"></input>
	</td>
      </tr> 
    </table>
  </div>

    <div class="formulaire">
    <table class="vide">
    <tr><B>Sequence model</B></tr>
    <tr>
    <td>This model is used to construct the mature chain
    </td>
    </tr>
      <tr>
	<td class="label">
	  <input id="col1A1" type="radio" name="model" value="col1a1"/>
	  col1a1
	</td>
      </tr>
      <tr>
	<td class="label">
	  <input id="col1A2" type="radio" name="model" value="col1a2" />
	  col1a2
	</td>
      </tr>
      <tr>
	<td class="label">
	  <input id="no_model" type="radio" name="model" value="no_model"/>
	  no model
	</td>
      </tr>
      <!--rajouter un bouton other pour donner la possinilité d'ajouter une séquence modèle-->
      <tr>
	<td class="label">
	  <input id="other" type="radio" name="model" value="other"/>
	  other
	</td>
      </tr>
    </table>
  </div>

  <div class="formulaire">
    <table class="vide">
      <tr>
	<td class="label"><b>Upload the taxonomy</b> <i>(optional, Newick format)</i>
	</td>
      </tr>
      <tr>
	<td>
	</td>
      </tr>
      <tr>
	<td class="label"> <B>upload</B> a file
	  <input type="file" name="file_taxonomy" value="file_taxonomy"></input>
	</td>
      </tr> 
    </table>
  </div>

  <div class="formulaire">
    <table class="vide">
     <tr><B>Parameters</B></tr>
     <tr>
    <td>In silico digestion
    </td>
     </tr>
      <tr>
	<td class="label">
	  <input id="digestion_none" type="radio" name="digestion" value="no_digestion" />
	  None
	</td>
    </tr>
    <tr>
    <td class="label">
	  <input id="digestion_enzyme" type="radio" name="digestion" value="digestion" />
    Select enzyme :
    <select id="enzyme" name="enzyme">
        <option value="42">Trypsin</option> <!--rajouter d'autres enzymes à voir sur la doc d'rpg-->
      </select>
	</td>
      </tr>
    </table>
  </div>

  <div class="formulaire">
    <table class="vide">
     <tr><B>Mass compatibility</B></tr>
      <tr>
	<td class="label">
	  <input id="check_mass" type="checkbox" name="check_mass" value="check_mass" />
	  Also check mass for comparison
  </td>
  <td>
    Select resolution :
    <select id="resolution" name="resolution">
        <option value="0.1">0.1</option>
        <option value="0.001">0.001</option>
      </select>
	</td>
    </tr>

  
  <div class="formulaire">
    <table class="vide">
      <tr>
	<td class="label"> 
	  Enter your <b>E-mail</b> address <i>(optional)</i>: 
	  <input id="email" type="text" name="email" size="20" />
	</td>
      </tr>
    </table>
  </div>

  <div class="center">
    <input type="submit" id="reset" name="reset" value="Reset" /> 
    <input type="submit" id="run" name="button" value="Run Example" /> <!--faire un exemple-->
    <input type="hidden" name="command" value="request" /> 
  </div>

</form>

</div><!--bloc -->
</div><!-- main-->

<!-- chargement de la librairie php lib.inc -->
   <?php require("../lib.inc")?>
<!-- appel de la fonction footer qui permet d'afficher au bas de la page (nom du logiciel, un lien vers le mail, la date de modif -->
<!-- A modifier en fonction de votre logiciel -->
   <?php footer("peptide_marker_validation","peptide_marker_validation", "tayna.pellegri.etu@univ-lille.fr","2022"); ?>
   <!-- php à modifer-->
                                                                                


</body>                                        
</html>        

