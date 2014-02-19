delayedRefreshSession = function(sessionKey) {
	setTimeout(function() { 
		var currentURL = window.location.href.split( '?' );
		window.location = currentURL[0] + "?session=" + sessionKey;
		}, 3000); 
}

$(document).ready(function() {
    $(window).load(function() {
        // this code will run after all other $(document).ready() scripts
        // have completely finished, AND all page elements are fully loaded.
		$( ":button[id='restore.session']" ).on( "click", function() { 
			var sessionKey = $( ":input[id='session.key']" ).val();
			var currentURL = window.location.href.split( '?' );
			window.location = currentURL[0] + "?session=" + sessionKey;
		});


		$( ":button[id='add.input']" ).on("click", function() {
			if ($( ":input[id='input.list']").val() == "") {
				$( ":input[id='input.list']").val( $( ":input[id='lookup']" ).val());
			} else {
				$( ":input[id='input.list']").val(  $( ":input[id='input.list']").val() + "\n" +  $( ":input[id='lookup']" ).val());

			}
	//		$( ":input[id='input.list']").val(  $( ":input[id='input.list']").val() + "\n" +  $( ":input[id='lookup']" ).val());
	//		alert( $( ":input[id='lookup']").val() );
		});
    });
});


