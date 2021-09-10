//jQuery time
var current_fs, next_fs, previous_fs; //fieldsets
var left, opacity, scale; //fieldset properties which we will animate
var animating; //flag to prevent quick multi-click glitches

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
    namebox.setAttribute("type", "text");
    namebox.setAttribute("name","optimized_name_#".concat((i+1).toString()))
    namebox.setAttribute("placeholder", "Enter Organism Name #".concat((i+1).toString()));
    document.getElementById("optimized_inputboxes").appendChild(namebox);
    //Add file-boxes for organism sequence
    var inputbox = document.createElement("input");
    inputbox.setAttribute("type", "file");
    inputbox.setAttribute("accept", ".gb");
    inputbox.setAttribute("name","optimized_seq_file_#".concat((i+1).toString()))
    document.getElementById("optimized_inputboxes").appendChild(inputbox);
  }
}

function deoptimize_inputbox() {
  $('#deoptimize_inputboxes').empty();
  var num2 = document.getElementById("num2").value;
  for (i = 0; i < num2; i++) {
    //Add text-boxes for organism names
    var namebox = document.createElement("input");
    namebox.setAttribute("type", "text");
    namebox.setAttribute("name","deoptimized_name_#".concat((i+1).toString()))
    namebox.setAttribute("placeholder", "Enter Organism Name #".concat((i+1).toString()));
    document.getElementById("deoptimized_inputboxes").appendChild(namebox);
    //Add file-boxes for organism sequence
    var inputbox = document.createElement("input");
    inputbox.setAttribute("type", "file");
    inputbox.setAttribute("accept", ".gb");
    inputbox.setAttribute("name","deoptimized_seq_file_#".concat((i+1).toString()))
    document.getElementById("deoptimized_inputboxes").appendChild(inputbox);
  }
}
function empty_inputbox(id_box, id_num) {
    $(id_box).empty();
    $(id_num).val('');
}