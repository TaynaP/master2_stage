/*
 * bioinfo.js
 *
 *
 * *** A FAIRE ****************************************************************
 * - le reset des formulaires par rechargement de la page est insatisfaisant
 * ****************************************************************************
 *
 *
 * Les URLs sont construites ainsi:
 * nom_du_serveur/nom_du_soft/onglet_actif.php
 * nom_du_serveur/nom_du_soft/result/id/results.php
 *
 * Le chargement direct d'une page d'un logiciel déclenche un traitement qui :
 * - affiche l'arborescence
 * - modifie le titre de la page pour qu'il corresponde bien à l'onglet
 * - donne la classe active à l'élément (onglet) de la page d'accueil du 
 *   logiciel du menu central
 * - vérifie si la gestion de l'historique peut être réalisée avec les 
 *   fonctionalités HTML5
 *   Si cela est possible alors le site est géré dynamiquement, sinon
 *   statiquement.
 *
 * Le clic sur un des éléments du menu central declenche un traitement de 
 * l'URL cible du lien href qui:
 * - dans le cas où le site est géré dynamiquement
 *   - donne la classe active à l'élément (ou onglet) du menu central
 *   - modifie le titre de la page pour qu'il corresponde bien à l'onglet
 *   - modifie le contenu HTML de l'élément #main de la page pour charger 
 *     celui correspondant à l'élément
 *   - ajoute à l'historique de navigation la page 
 * - dans le cas où le site n'est pas géré dynamiquement
 *   - suit le lien et procède à un chargement classique
 *
 */

var soft_theme = {
    /* theme DNA */
    "yass":"dna", 
    "magnolia":"dna", 
    "procars":"dna", 
    /* theme RNA */
    "carnac":"rna", 
    "gardenia":"rna", 
    "regliss":"rna", 
    "CGseq":"rna", 
    "mirkwood":"rna", 
    /* theme HTS */
    "sortmerna":"hts", 
    "crac":"hts", 
    "vidjil":"hts", 
    "storm":"hts",
    "matam":"hts",
    "dinamo":"hts", 
    /* theme PROTEINS */
    "path":"proteins", 
    "protea":"proteins", 
    "reblosum":"proteins", 
    /* theme TFM */
    "tfm-explorer":"tfm", 
    "tfm-scan":"tfm", 
    "tfm-pvalue":"tfm", 
    "tfm-cuda":"tfm", 
    /* theme NRP */
    "norine":"nrp", 
    "mynorine":"nrp", 
    "smiles2mono":"nrp", 
    "smilescolor":"nrp"
};

// la variable gloable dynamic est mise a jour après que le DOM de la
// page du logiciel soit chargée et permet de définir si le site est
// géré dynamiquement ou non
var dynamic = true;

/******************************************************************************
 * Fonctions de mise en place de la page
 ******************************************************************************/

function is_a_soft_page (url) {
    var url_tab = url.split("/");
    return (url_tab[1] in soft_theme);
}

/*
 * Une fois que le DOM est chargé :
 * - vérifie si le site doit être dynamique ou non
 * - mets à jour le fil d'arianne si on est sur une page de logiciel
 * - active l'onglet
 */
$(document).ready(function(){
    if (! (history.pushstate == undefined)) {
	dynamic = false;
    }
//    alert("DYNAMIC " + dynamic);
    var url=window.location.pathname;
    // si c'est un logiciel, on affiche le fil d'arianne
    if (is_a_soft_page(url)) {
	display_arborescence(url);
    }
    set_appropriate_title(url);
    // active l'onglet
    setTimeout(function(){active_class(url);}, 200); // ??? pourquoi un timeout ???
});

/*
 * Affiche le fil d'arianne sur la page.
 *
 * Parametres
 *  url: URL complète d'une des pages d'un logiciel
 */
function display_arborescence(url){
    var tab_url=url.split("/");
    var name_soft=tab_url[1];
    if(soft_theme.hasOwnProperty(name_soft)){
	var theme = soft_theme[name_soft];
	$("#arborescence").html("<a href='/'>home</a> &gt; <a href='/topic/"+theme+".html'>"+theme.toUpperCase()+"</a> &gt; "+"<a href='/"+name_soft+"/"+name_soft.toLowerCase()+".php'>"+name_soft+"</a>");
    }
}

/*
 * Modifie le titre de la page dynamiquement;
 *
 * Parametres
 *  url: URL complète d'une des pages 
 */
function set_appropriate_title (url) {
    var tab_url = url.split("/");
    var name_soft = tab_url[1];
    if (is_a_soft_page (url)) {
	var onglet = tab_url[2].replace(/.php/,"");
	var title = "Bonsai bioinformatics - " + name_soft
	if (onglet != name_soft) {
	    // on va chercher le texte de l'onglet affiché sur la page
	    var elt_current=".tabs"; 
	    var val = $(elt_current).find('a[href="'+url+'"]').text();
	    title = title + " - " + val;
	}    
	$(document).attr("title", title);
    } else if (name_soft == "topic") {
	var onglet = tab_url[2].replace(/.php/,"");
	onglet = onglet.replace(/.html/,"");
	var title = "Bonsai bioinformatics - " + onglet;
	$(document).attr("title", title);
    }	

}

/******************************************************************************
 * Fonctions de mise en place des evenements
 ******************************************************************************/

/*
 * evt      : click sur un element de classe aLoad
 * concerne : les boutons du menu central
 *
 * Recupere le nom du soft et XXX a partir de l'URL de la page puis
 * lance le chargement effectif.
 */
$(document).on("click", ".aLoad", function(){
    var url=$(this).attr("href");
    var tab=url.split("/");
    loadPageSoft(tab[1],tab[2]);
    return false;
});

/*
 * evt      : clic sur le bouton home
 * concerne : le bouton-image "home" a gauche du menu central
 */
$(document).on("click", ".returnHome", function(){
    window.location="/";
    return false;
});

/*
 * evt      : clic sur un elements de classe normal_link
 * concerne : ???
 */
$(document).on("click", ".normal_link", function(){
    var url=$(this).attr("href");
    window.location=url;
    return false;
});


/*
 * evt      : clic sur un element de classe button_form_fill
 * concerne : ???
 */
$(document).on('click', '.button_form_fill', function(){
    var url = document.location.href;
    tab=url.split("/");
    // récupère l'id et le nom du logiciel à partir de l'url 
    var name_soft = tab[3].substring(1);
    var id = tab[5];
    // appelle loafFormFill qui permet d'afficher le formulaire prè-rempli
    url_script="cgi-bin/"+name_soft+"/read_param_txt.py";
    loadFormFill(name_soft, id, url_script);
    return false;
});


/******************************************************************************
 * Gestion du chargement des pages
 ******************************************************************************/

/* 
 * Afficher le corps de la page correspondant a un onglet
 *
 * Parametres:
 *   name_soft : le nom du soft
 *   onglet    : l'onglet a activer pour le soft
 */
function loadPageSoft(name_soft, onglet){
    var url_soft = "/"+name_soft+"/"+onglet;
//    alert("LOAD PAGE SOFT : " + url_soft);
    if (dynamic) {
	// cas ou l'historique peut être géré manuellement
	set_appropriate_title (url_soft);
	push_in_history (url_soft); 
	// #main est l'id de la div du contenu de la page du soft, de plus
	// l'ajout de '#center' précédé d'un espace permet de ne charger
	// que le code HTML correspondant à l'élément d'id 'center' de la
	// page 'url_soft'
	$("#main").load(url_soft + " #center"); 
    } else {
	// cas ou l'historique ne peut etre géré manuellement
	window.location = url_soft;
    }
    // ??? pourquoi active_class est appelé avec un timeout et pas en callback du load ???
    setTimeout(function(){active_class(url_soft);}, 200);
    // ??? gestion des erreurs dans la fonction de callback également : par exemple si lancement du programme plante ???
}

/* 
 * Activer l'onglet sélectionne
 *
 * Parametres:
 *   url : l'URL nom du soft avec l'onglet a activer
 */
function active_class(url){
    if (url.search("results.php") >= 0) {
	var tab_url = url.split("/");
	url = "/" + tab_url[1] + "/id.php";
    }
//    alert("ACTIVE CLASS: " + url);
    // tabs est la classe de la div qui contient le menu central
    var elt_current=".tabs"; 
    // ajoute la classe active a l'element parent de l'element d'interet
    $(elt_current).find('a[href="'+url+'"]').parent().addClass("active");
    // supprime la classe active a tous les elements freres
    $(elt_current).find('a[href="'+url+'"]').parent().siblings().removeClass("active");
}

/******************************************************************************
 * Gestion de l'historique des pages
 ******************************************************************************/

/*
 * Permet d'utiliser les boutons back et forward du navigateur
 *
 * Parametres:
 *   event: l'evenement declenchant, envoye autamatiquement par le navigateur
 * 
 * NOTE: d'apres le MDN ne fonctionne qu'a partir de HTML 5 et
 * navigateurs recents
 * https://developer.mozilla.org/fr/docs/Web/Guide/DOM/Manipuler_historique_du_navigateur
 */
window.onpopstate = function(event) {
    var url=window.location.pathname;
// JSV    var state_url = event.state;
// JSV   if (state_url != null) { alert("POP STATE " + state_url["url"] + " vs. " + url); }
    loadPageSoftAfterHistoryChange(url);
};


/*
 * Permet d'enregistrer la page accédée dans l'historique.
 *
 * Parametres:
 *  url_soft: l'URL complète de la page à enregister
 */
function push_in_history (url_soft) {
    // construit la chaîne qui sera visible dans l'historique
    var targetTitle = $(this).attr('title');
    // /!\ la manipulation de la pile de l'historique n'est valide que sous HTML5 et navigateurs recents 
    history.pushState({url : ""+url_soft}, targetTitle, url_soft);
}

/*
 * Permet d'afficher le corps de la page correspondant à l'onglet
 * après le bouton back ou forward du navigateur
 *
 * Parametre:
 *   url: l'URL du HTML a charger dans la div #main
 *
 * ??? cette fonction est en fait une sous-fonction de loadPageSoft ???
 */
function loadPageSoftAfterHistoryChange(url){
//    alert("LOAD PAGE AFTER HISTORY CHANGE: " + url);
    $("#main").load(url+" #center");
    set_appropriate_title(url);
    setTimeout(function(){active_class(url);}, 200); // ??? ne fonctionne pas sur els resultats car l'url n'a pas le bon format ???
}

/******************************************************************************
 * Gestion des formulaires
 ******************************************************************************/

/* A modifier si utilisation des url normal ??? je ne comprends pas ???
 *
 * Parametres:
 *   data_form:
 *   url_soft:
 *   name_soft:
 */
function sendForm(data_form, url_soft, name_soft){
    $.ajax({
        url: url_soft,
        type: "POST",            // Type of request to be send, called as method                                                                                                        
	data: data_form,
        dataType:"json", // Data sent to server, a set of key/value pairs (i.e. form fields and values)                                                                      
        contentType: false,       // The content type used when sending data to the server.                                                                                
        cache: false,             // To unable request pages to be cached                                                                                            
        processData:false,        // To send DOMDocument or non processed data file it is set to false                      
        success: function(response){
	    setTimeout(function(){  
		url_result="/"+name_soft+"/result/"+response.run_id+"/results.php";
		push_in_history(url_result);
		$("#main").load(url_result+" #center");
            }, 200);
        },
        error: function(request, status, error){
            $("#main").html(request.responseText);
        }
    });
}


/*
 * Envoie les données du formulaire au wrapper, récupère le résultat
 * au format html, affiche l'id du résultat dans la page wait
 *
 * Parametres
 *   data_form:
 *   name_soft: nom imple du logiciel
 *   wrapper: nom du script permettant de lancer le calcul
 *
 * NOTE: cette fonction est ensuite associée au clic sur le bouton de
 * soumission du formulaire dans le script.js spécifique au
 * logiciel.
 * Typiquement :
 * $(document).on("submit","#formulaire", function(){
 *     var param = loadParaSoft(soft);
 *     var formData=new FormData(this);
 *     loadFormulaire(formData, param[1], param[3]);
 *     return false;
 * });
 */
function loadFormulaire(data_form, name_soft, wrapper){
//    alert("LOAD FORMULAIRE " + name_soft + " " + wrapper);
    var url_wait="/wait.html";
    var url_cgi="/cgi-bin/"+name_soft+"/";
    $('#main').load(url_wait);
    var url_soft = url_cgi+wrapper;
    var url_create_id = url_cgi+"create_id.py";
    $.ajax({
	url:url_create_id,
	async:false,
        datatype : "json",
	cache : false,
	success: function(data) {
	    setTimeout(function(){
		document.getElementById("name_soft_wait").innerHTML = "You request has been successfully submitted to <B>"+name_soft.toUpperCase()+"</B>";
		document.getElementById("id_wait").innerHTML = "Your ID is <B>"+data.run_id+"</B><br/>";
	    }, 100);
	    data_form.append("run_id", data.run_id);
	}
    });
    sendForm(data_form, url_soft, name_soft);
}


/* 
 * Charger les paramêtres nécessaires à partir du fichier "para_software.txt" 
 * Chargee du nom du dossier et des deux wrappers (pour le formulaire et l'id)
 *
 * Parametres:
 *   rep-html:
 *
 * NOTE: cette fonction doit être appelée avant la fonction loadFormulaire 
 * au moment où cette dernière est appelée à la soumission du formulaire.
 */
function loadParaSoft(rep_html){
//    alert("LOAD PARA SOFT");
    url_para=rep_html+"para_software.txt";
    var tab_para = new Array();
    $.ajax({
        type:"GET",
        url:url_para,
	async:false,
        success: function(data){
	    tab_para = data.split(":");
	},
	error: function(){
            alert("Probleme dans la lecture de para_software.txt");
        }
    });
    return tab_para
}


/*
 * Afficher les résultats sauvegardés pour un ID.
 *
 * Parametres
 *   name_soft: nom simple du logiciel
 *   id: numero du job
 *
 * NOTE: cette fonction est ensuite associée au clic sur le bouton de
 * chargement des resultats par ID dans le script.js spécifique au
 * logiciel.
 * Typiquement :
 * $(document).on("submit","#form_id", function(){
 *   var id = document.getElementById("run_id").value;
 *   loadResultByID(soft, id);    
 *   return false;
 * });
 */
function loadResultByID(name_soft, id){
    url_result = name_soft+"result/"+id+"/results.php";
    $("#main").load(url_result+" #center");
//    History.pushState("soft", "bioinfo.lifl.fr-Bonsai bioinformatics", url_result);
}


/******************************************************************************/


/* Fonction utilisée pour les box se trouvant dans la page d'index et pour les pages de thèmes
 * Permet d'aller directement vers la page, lors d'un clique sur un lien se situant dans les box.
 * 
 * Paramètre
 *   url : path pour aller vers la page
 *
 */ 
function loadLink(url){
    window.location=url;
    return false;
}



/*
 * Permet de charger les données prédéfinies d'un example dans un
 * formulaire.
 *
 * Parametres:
 *   name_file: fichier contenant les données de l'exemple
 *   id: id du bouton exemple du formulaire
 *
 * NOTE: cette fonction est ensuite associée au bouton exemple
 * correspondant dans le script.js spécifique au logiciel. 
 * Typiquement pour un bouton d'id ex1
 * $(document).on("click","#ex1", function(){
 *   loadExample(soft+"example1", "#paste_seq");
 *   return false;
 *   });
 *
 */
function loadExample(name_file, id){
    $.ajax({
        type:"GET",
        url:name_file,
        success: function(data){
	    if (i_am_old_ie) {
		document.getElementById(id.substring(1)).value=data;
	    } else {
		$(id).html(data);
	    }
        },
        error: function(){
	    $(id).html('Une erreur est survenue.');
        }
    });
}


/*
 * Permet de modifier les données arrivant du formulaire au format JSON.
 * Ceci se réalise dans le fichier scripts.js de chaque logiciel.
 *
 * Parametres:
 *   les données de base du formulaire
 * 
 * Retour:
 *   o : données du formulaire au format JSON
 *  
 * NOTE: cette fonction est appelée dans le script.js spécifique au logiciel. 
 * 
 */
$.fn.serializeObject = function()
{
    var o = {};
    var a = this.serializeArray();
    $.each(a, function() {
        if (o[this.name] !== undefined) {
            if (!o[this.name].push) {
                o[this.name] = [o[this.name]];
            }
            o[this.name].push(this.value || '');
        } else {
            o[this.name] = this.value || '';
        }
    });
    return o;
};

/******Peut être utile pour plus tard**************/

// Cette fonction devait servir à afficher le formulaire pré-rempli lors d'un retour sur ce dernier à partir des résultats.
// Cette fonctionnalité n'est plus au gout du jour, mais on peut la conserver si un jour l'envie vous dit de la remettre.
function loadFormFill(name_soft, id, url_script){
    url_form=name_soft+"/form.php";
    window.location.href = "#"+url_form;
    $("#main").load("/"+url_form);
    $.post(url_script, { id: id, name_soft: name_soft }, function(data ){
	for(var key in data){
	    if((key=="sequence")||(key=="seq")&&(document.getElementById("paste_seq"))){
		document.getElementById("paste_seq").value = data[key];
	    }
	    else if ((document.getElementsByName(key))&&(data[key]!="None")){
		var x = document.getElementsByName(key);
		var i;
		for (i = 0; i < x.length; i++) {
		    if(x[i].getAttribute("type")=="checkbox"){
			if(data[key]=="False"){
			    x[i].checked = false;
			}
			else{
			    x[i].checked = true;
			}
		    }
		    else if (x[i].getAttribute("type")=="radio"){
			document.getElementById(data[key]).checked = true;
		    }
		    else{
			document.getElementById(key).value = data[key];
		    }
		}
	    }	    
	}
    });
};







