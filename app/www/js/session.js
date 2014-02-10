$(document).ready(function() {
    $(window).load(function() {
        // this code will run after all other $(document).ready() scripts
        // have completely finished, AND all page elements are fully loaded.
		$( ":button[id='restore.session']" ).on( "click", function() { 
			var sessionKey = $( ":input[id='session.key']" ).val()
			var currentURL = window.location.href.split( '?' );
			window.location = currentURL[0] + "?session=" + sessionKey
		});
    });
});

