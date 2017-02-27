// Subroutine to calculate the number of photons incident on the flux monitor
//
void
fluxcal( Float_t factor, Float_t factor_err,
	Float_t monitor_counts, Float_t veto_counts, Float_t beam_bunches,
	Float_t *photons, Float_t *photons_err)
{
	Float_t r1 = (beam_bunches - veto_counts);
	Float_t r2 = (beam_bunches - veto_counts - monitor_counts);
	Float_t ratio = r1 / r2;
	Float_t gammas = factor * beam_bunches * log(ratio);

	// estimate error
	Float_t a = gammas / factor;
	Float_t b = factor * beam_bunches * monitor_counts / r1 /r2;
	Float_t c = factor * beam_bunches / r2;
	Float_t sum = pow( (a * factor_err), 2)
			+ pow(b, 2) * veto_counts
			+ pow(c, 2) * monitor_counts;
	Float_t gammas_err = sqrt(sum);

	*photons = gammas / 1.060;// 1.060 +/- 0.006 = attenuation from air to molly, from Rob!
	*photons_err = gammas_err / 1.060;
}

Int_t
test_flux()
{
	Float_t factor, factor_err;
	Float_t monitor_counts, veto_counts, beam_bunches;

	printf("Enter Calibration Factor: ");
	scanf("%f", &factor);
	printf("Enter Calibration Factor error: ");
	scanf("%f", &factor_err);

	printf("Enter 5-Paddle counts during Live Time: ");
	scanf("%f", &monitor_counts);
	printf("Enter Veto Hit Counts during Live Time: ");
	scanf("%f", &veto_counts);
	printf("Enter Number of beam bunches during Live Time: ");
	scanf("%f", &beam_bunches);

	Float_t photons, photons_err;
	(void) fluxcal( factor, factor_err, monitor_counts,
		veto_counts, beam_bunches,
		&photons, &photons_err);

	printf ("\nNumber of photons = %.0f +/- %.0f\n", photons, photons_err);
	return 0;
}

