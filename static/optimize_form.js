//jQuery time
var current_fs, next_fs, previous_fs; //fieldsets
var left, opacity, scale; //fieldset properties which we will animate
var animating; //flag to prevent quick multi-click glitches

$("#next_1").click(function(){
	// Check number of optimized organisms
    let valid = $('#num1').val();
    if(valid == '' || valid == 0){
        alert("Please upload at least 1 organism to optimize!");
        return false;
    }
    for (i = 0; i < valid; i++) {
        let seq_file_id = "optimized_seq_file_#".concat((i+1).toString());
        let seq_file = document.getElementById(seq_file_id).value;
        if (seq_file == '') {
            alert("Please upload .gb files for all organisms");
            return false;
        }
    }

	if(animating) return false;
	animating = true;

	current_fs = $(this).parent();
	next_fs = $(this).parent().next();

	//activate next step on progressbar using the index of next_fs
	$("#progressbar li").eq($("fieldset").index(next_fs)).addClass("active");

	//show the next fieldset
	next_fs.show();
	//hide the current fieldset with style
	current_fs.animate({opacity: 0}, {
		step: function(now, mx) {
			//as the opacity of current_fs reduces to 0 - stored in "now"
			//1. scale current_fs down to 80%
			scale = 1 - (1 - now) * 0.2;
			//2. bring next_fs from the right(50%)
			left = (now * 50)+"%";
			//3. increase opacity of next_fs to 1 as it moves in
			opacity = 1 - now;
			current_fs.css({'transform': 'scale('+scale+')'});
			next_fs.css({'left': left, 'opacity': opacity});
		},
		duration: 800,
		complete: function(){
			current_fs.hide();
			animating = false;
		},
		//this comes from the custom easing plugin
		easing: 'easeInOutBack'
	});
})

$("#next_2").click(function(){
	// Check number of deoptimized organisms
    let valid = $('#num2').val();
    if(valid == '' || valid == 0){
        alert("Please upload at least 1 organism to optimize!");
        return false;
    }
    for (i = 0; i < valid; i++) {
        let seq_file_id = "deoptimized_seq_file_#".concat((i+1).toString());
        let seq_file = document.getElementById(seq_file_id).value;
        if (seq_file == '') {
            alert("Please upload .gb files for all organisms");
            return false;
        }
    }

	if(animating) return false;
	animating = true;

	current_fs = $(this).parent();
	next_fs = $(this).parent().next();

	//activate next step on progressbar using the index of next_fs
	$("#progressbar li").eq($("fieldset").index(next_fs)).addClass("active");

	//show the next fieldset
	next_fs.show();
	//hide the current fieldset with style
	current_fs.animate({opacity: 0}, {
		step: function(now, mx) {
			//as the opacity of current_fs reduces to 0 - stored in "now"
			//1. scale current_fs down to 80%
			scale = 1 - (1 - now) * 0.2;
			//2. bring next_fs from the right(50%)
			left = (now * 50)+"%";
			//3. increase opacity of next_fs to 1 as it moves in
			opacity = 1 - now;
			current_fs.css({'transform': 'scale('+scale+')'});
			next_fs.css({'left': left, 'opacity': opacity});
		},
		duration: 800,
		complete: function(){
			current_fs.hide();
			animating = false;
		},
		//this comes from the custom easing plugin
		easing: 'easeInOutBack'
	});
})

$(".next").click(function(){
	if(animating) return false;
	animating = true;

	current_fs = $(this).parent();
	next_fs = $(this).parent().next();

	//activate next step on progressbar using the index of next_fs
	$("#progressbar li").eq($("fieldset").index(next_fs)).addClass("active");

	//show the next fieldset
	next_fs.show();
	//hide the current fieldset with style
	current_fs.animate({opacity: 0}, {
		step: function(now, mx) {
			//as the opacity of current_fs reduces to 0 - stored in "now"
			//1. scale current_fs down to 80%
			scale = 1 - (1 - now) * 0.2;
			//2. bring next_fs from the right(50%)
			left = (now * 50)+"%";
			//3. increase opacity of next_fs to 1 as it moves in
			opacity = 1 - now;
			current_fs.css({'transform': 'scale('+scale+')'});
			next_fs.css({'left': left, 'opacity': opacity});
		},
		duration: 800,
		complete: function(){
			current_fs.hide();
			animating = false;
		},
		//this comes from the custom easing plugin
		easing: 'easeInOutBack'
	});
});

$(".previous").click(function(){
	if(animating) return false;
	animating = true;

	current_fs = $(this).parent();
	previous_fs = $(this).parent().prev();

	//de-activate current step on progressbar
	$("#progressbar li").eq($("fieldset").index(current_fs)).removeClass("active");

	//show the previous fieldset
	previous_fs.show();
	//hide the current fieldset with style
	current_fs.animate({opacity: 0}, {
		step: function(now, mx) {
			//as the opacity of current_fs reduces to 0 - stored in "now"
			//1. scale previous_fs from 80% to 100%
			scale = 0.8 + (1 - now) * 0.2;
			//2. take current_fs to the right(50%) - from 0%
			left = ((1-now) * 50)+"%";
			//3. increase opacity of previous_fs to 1 as it moves in
			opacity = 1 - now;
			current_fs.css({'left': left});
			previous_fs.css({'transform': 'scale('+scale+')', 'opacity': opacity});
		},
		duration: 800,
		complete: function(){
			current_fs.hide();
			animating = false;
		},
		//this comes from the custom easing plugin
		easing: 'easeInOutBack'
	});
});

//$(".submit").click(function(){
//	return false;
//})

function optimize_inputbox() {
  $('#optimize_inputboxes').empty();
  var num1 = document.getElementById("num1").value;
  for (i = 0; i < num1; i++) {
    //Add text-boxes for organism names
    var namebox = document.createElement("input");
    let label_name = document.createElement('h4');
    label_name.textContent = "Optimized Organism #".concat((i+1).toString(),":");
    namebox.setAttribute("type", "text");
    namebox.setAttribute("name","optimized_name_#".concat((i+1).toString()))
    namebox.setAttribute("placeholder", "Enter Organism Name #".concat((i+1).toString()));
    document.getElementById("optimized_inputboxes").appendChild(document.createElement('hr'));
    document.getElementById("optimized_inputboxes").appendChild(document.createElement('br'));
    document.getElementById("optimized_inputboxes").appendChild(label_name);
    document.getElementById("optimized_inputboxes").appendChild(namebox);
    //Add file-boxes for organism sequence
    let inputbox = document.createElement("input");
    let temp_id = "optimized_seq_file_#".concat((i+1).toString());
    let label_sequence = document.createElement('h5');
    let label_express_lvl = document.createElement('h5');
    label_sequence.textContent = 'Upload Genebank File #'.concat((i+1).toString(),' (mandatory):');
    label_express_lvl.textContent = 'Upload CSV File #'.concat((i+1).toString(),' (optional):');

    inputbox.setAttribute("type", "file");
    inputbox.setAttribute("accept", ".gb");
    inputbox.setAttribute("id",temp_id);
    inputbox.setAttribute("name",temp_id);
    inputbox.setAttribute("onchange", `return fileValidation('${temp_id}', '.gb')`);
    document.getElementById("optimized_inputboxes").appendChild(label_sequence);
    document.getElementById("optimized_inputboxes").appendChild(inputbox);
    // Add file-boxes for organism gene-expression
    inputbox = document.createElement("input");
    temp_id = "optimized_express_lvl_file_#".concat((i+1).toString());
    inputbox.setAttribute("type", "file");
    inputbox.setAttribute("accept", ".csv");
    inputbox.setAttribute("id",temp_id);
    inputbox.setAttribute("name",temp_id);
    inputbox.setAttribute("onchange", `return fileValidation('${temp_id}', '.csv')`);
    document.getElementById("optimized_inputboxes").appendChild(label_express_lvl);
    document.getElementById("optimized_inputboxes").appendChild(inputbox);
  }
}

function deoptimize_inputbox() {
  $('#deoptimize_inputboxes').empty();
  var num2 = document.getElementById("num2").value;
  for (i = 0; i < num2; i++) {
    //Add text-boxes for organism names
    var namebox = document.createElement("input");
    let label_name = document.createElement('h4');
    label_name.textContent = "Deoptimized Organism #".concat((i+1).toString(),":");
    namebox.setAttribute("type", "text");
    namebox.setAttribute("name","deoptimized_name_#".concat((i+1).toString()))
    namebox.setAttribute("placeholder", "Enter Organism Name #".concat((i+1).toString()));
    document.getElementById("deoptimized_inputboxes").appendChild(document.createElement('hr'));
    document.getElementById("deoptimized_inputboxes").appendChild(document.createElement('br'));
    document.getElementById("deoptimized_inputboxes").appendChild(label_name);
    document.getElementById("deoptimized_inputboxes").appendChild(namebox);
    //Add file-boxes for organism sequence
    var inputbox = document.createElement("input");
    let temp_id = "deoptimized_seq_file_#".concat((i+1).toString());
    let label_sequence = document.createElement('h5');
    let label_express_lvl = document.createElement('h5');
    label_sequence.textContent = 'Upload Genebank File #'.concat((i+1).toString(),' (mandatory):');
    label_express_lvl.textContent = 'Upload CSV File #'.concat((i+1).toString(),' (optional):');
    inputbox.setAttribute("type", "file");
    inputbox.setAttribute("accept", ".gb");
    inputbox.setAttribute("id",temp_id);
    inputbox.setAttribute("name", temp_id);
    inputbox.setAttribute("onchange", `return fileValidation('${temp_id}', '.gb')`);
    document.getElementById("deoptimized_inputboxes").appendChild(label_sequence);
    document.getElementById("deoptimized_inputboxes").appendChild(inputbox);
    // Add file-boxes for organism gene-expression
    inputbox = document.createElement("input");
    temp_id = "deoptimized_express_lvl_file_#".concat((i+1).toString());
    inputbox.setAttribute("type", "file");
    inputbox.setAttribute("accept", ".csv");
    inputbox.setAttribute("id",temp_id);
    inputbox.setAttribute("name",temp_id);
    inputbox.setAttribute("onchange", `return fileValidation('${temp_id}', '.csv')`);
    document.getElementById("deoptimized_inputboxes").appendChild(label_express_lvl);
    document.getElementById("deoptimized_inputboxes").appendChild(inputbox);
  }
}
function empty_inputbox(id_box, id_num) {
    $(id_box).empty();
    $(id_num).val('');
}
function fileValidation(element_id, allowed_type) {
    var fileInput = document.getElementById(element_id);
    var filePath = fileInput.value;
    // Allowing file type
    switch (allowed_type){
        case '.fa':
            var allowedExtensions = /(.fa)$/i;
            break;
        case '.csv':
            var allowedExtensions = /(.csv)$/i;
            break;
        case '.gb':
            var allowedExtensions = /(.gb)$/i;
            break;
        default:
            console.log('Please add allowed types');
    }
    if (!allowedExtensions.exec(filePath)) {
        alert('Invalid file type');
        fileInput.value = '';
        return false;
    }
}

var slider = document.getElementById("tuning_param");
var output = document.getElementById("demo");
output.innerHTML = slider.value;

slider.oninput = function() {
  output.innerHTML = this.value;
  $("#tuning_param_text").val(this.value)
}
