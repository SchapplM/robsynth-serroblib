% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRR11
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR11_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR11_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR11_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (14->13), mult. (48->27), div. (0->0), fcn. (38->6), ass. (0->12)
	t114 = sin(pkin(6));
	t124 = t114 * (r_i_i_C(3) + qJ(2));
	t116 = cos(pkin(6));
	t117 = sin(qJ(1));
	t122 = t116 * t117;
	t118 = cos(qJ(1));
	t121 = t116 * t118;
	t120 = qJD(1) * t114;
	t119 = t114 * qJD(2);
	t115 = cos(pkin(12));
	t113 = sin(pkin(12));
	t1 = [t118 * t119 + ((t113 * t122 - t115 * t118) * r_i_i_C(1) + (t113 * t118 + t115 * t122) * r_i_i_C(2) - t118 * pkin(1) - t117 * t124) * qJD(1), t118 * t120, 0, 0, 0, 0; t117 * t119 + ((-t113 * t121 - t115 * t117) * r_i_i_C(1) + (t113 * t117 - t115 * t121) * r_i_i_C(2) - t117 * pkin(1) + t118 * t124) * qJD(1), t117 * t120, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:01
	% EndTime: 2019-10-10 01:04:01
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (103->49), mult. (358->92), div. (0->0), fcn. (352->10), ass. (0->47)
	t248 = cos(pkin(7));
	t278 = r_i_i_C(3) + pkin(9);
	t279 = t278 * t248 + qJ(2);
	t250 = sin(qJ(3));
	t277 = t250 * r_i_i_C(2);
	t246 = sin(pkin(6));
	t251 = sin(qJ(1));
	t276 = t246 * t251;
	t253 = cos(qJ(1));
	t275 = t246 * t253;
	t274 = t248 * t250;
	t252 = cos(qJ(3));
	t273 = t248 * t252;
	t244 = sin(pkin(12));
	t272 = t251 * t244;
	t247 = cos(pkin(12));
	t271 = t251 * t247;
	t270 = t253 * t244;
	t269 = t253 * t247;
	t268 = qJD(1) * t251;
	t267 = qJD(1) * t253;
	t245 = sin(pkin(7));
	t266 = qJD(3) * t245;
	t265 = t278 * t245;
	t249 = cos(pkin(6));
	t264 = t249 * t272;
	t263 = t246 * t268;
	t262 = t246 * t267;
	t261 = t266 * t275;
	t260 = r_i_i_C(1) * t250 + r_i_i_C(2) * t252;
	t236 = -qJD(1) * t264 + t247 * t267;
	t237 = -t249 * t269 + t272;
	t259 = qJD(3) * t237 * t248 - t236;
	t258 = t260 * t245;
	t256 = t249 * t271 + t270;
	t257 = t245 * t276 - t248 * t256;
	t238 = t249 * t270 + t271;
	t233 = t237 * qJD(1);
	t255 = t233 * t248 + t245 * t262;
	t235 = t256 * qJD(1);
	t254 = -qJD(3) * t238 - t235 * t248 + t245 * t263;
	t241 = t252 * t261;
	t240 = -t264 + t269;
	t234 = t238 * qJD(1);
	t232 = -t234 * t252 + t255 * t250 + (-t240 * t250 + t257 * t252) * qJD(3);
	t231 = t234 * t250 + t255 * t252 + (-t240 * t252 - t257 * t250) * qJD(3);
	t1 = [-pkin(1) * t267 + t241 * r_i_i_C(1) + (-t252 * r_i_i_C(1) - pkin(2) + t277) * t236 + (t260 * t248 - t265) * t235 + ((t237 * t273 + t238 * t250) * r_i_i_C(1) + (-t237 * t274 + t238 * t252) * r_i_i_C(2)) * qJD(3) + ((-t266 * t277 + qJD(2)) * t253 + (-t258 - t279) * t268) * t246, t262, t231 * r_i_i_C(1) - t232 * r_i_i_C(2), 0, 0, 0; qJD(2) * t276 - t234 * pkin(2) + t232 * r_i_i_C(1) + t231 * r_i_i_C(2) - t233 * t265 + (-t251 * pkin(1) + t279 * t275) * qJD(1), t263, t241 * r_i_i_C(2) + (t254 * r_i_i_C(1) + t259 * r_i_i_C(2)) * t252 + ((t259 + t261) * r_i_i_C(1) - t254 * r_i_i_C(2)) * t250, 0, 0, 0; 0, 0, (-t249 * t258 + ((-t244 * t252 - t247 * t274) * r_i_i_C(1) + (t244 * t250 - t247 * t273) * r_i_i_C(2)) * t246) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:02
	% EndTime: 2019-10-10 01:04:03
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (281->64), mult. (951->111), div. (0->0), fcn. (978->12), ass. (0->56)
	t357 = sin(qJ(3));
	t351 = sin(pkin(7));
	t352 = sin(pkin(6));
	t360 = cos(qJ(1));
	t389 = t352 * t360;
	t377 = t351 * t389;
	t355 = cos(pkin(7));
	t356 = cos(pkin(6));
	t354 = cos(pkin(12));
	t383 = t360 * t354;
	t350 = sin(pkin(12));
	t358 = sin(qJ(1));
	t386 = t358 * t350;
	t367 = t356 * t383 - t386;
	t401 = t367 * t355;
	t370 = -t401 + t377;
	t403 = t370 * t357;
	t378 = t356 * t386;
	t381 = qJD(1) * t360;
	t338 = -qJD(1) * t378 + t354 * t381;
	t402 = -qJD(3) * t377 + t338;
	t359 = cos(qJ(3));
	t388 = t355 * t357;
	t391 = t351 * t356;
	t400 = (t350 * t359 + t354 * t388) * t352 + t357 * t391;
	t342 = -t378 + t383;
	t384 = t360 * t350;
	t385 = t358 * t354;
	t366 = t356 * t385 + t384;
	t390 = t352 * t358;
	t369 = t351 * t390 - t355 * t366;
	t399 = t342 * t359 + t369 * t357;
	t337 = t366 * qJD(1);
	t382 = qJD(1) * t359;
	t375 = t351 * t352 * t382;
	t387 = t355 * t359;
	t340 = t356 * t384 + t385;
	t393 = t340 * t359;
	t398 = (-t367 * t388 - t393) * qJD(3) - t337 * t387 + t358 * t375 - t402 * t357;
	t376 = qJD(1) * t390;
	t397 = (-qJD(3) * t340 - t337 * t355 + t351 * t376) * t357 + (qJD(3) * t401 + t402) * t359;
	t396 = r_i_i_C(3) + qJ(4);
	t395 = t337 * t351;
	t380 = t352 * qJD(2);
	t349 = sin(pkin(13));
	t353 = cos(pkin(13));
	t372 = -t353 * r_i_i_C(1) + t349 * r_i_i_C(2) - pkin(3);
	t362 = -t342 * t357 + t369 * t359;
	t361 = t367 * t351 + t355 * t389;
	t336 = t340 * qJD(1);
	t333 = -t355 * t376 - t395;
	t332 = t361 * qJD(1);
	t331 = t400 * qJD(3);
	t325 = qJD(1) * t403 + t362 * qJD(3) - t336 * t359;
	t324 = t399 * qJD(3) - t336 * t357 - t360 * t375 + t382 * t401;
	t1 = [(t333 * t349 - t353 * t397) * r_i_i_C(1) + (t333 * t353 + t349 * t397) * r_i_i_C(2) - t397 * pkin(3) - (t340 * t357 + t370 * t359) * qJD(4) - t338 * pkin(2) - pkin(9) * t395 + t360 * t380 + t396 * t398 + (-t360 * pkin(1) + (-pkin(9) * t355 - qJ(2)) * t390) * qJD(1), t352 * t381, t399 * qJD(4) + t372 * t324 + t396 * t325, t324, 0, 0; (t325 * t353 + t332 * t349) * r_i_i_C(1) + (-t325 * t349 + t332 * t353) * r_i_i_C(2) + t325 * pkin(3) - t362 * qJD(4) - t336 * pkin(2) + t358 * t380 + t396 * t324 + (-t358 * pkin(1) + t361 * pkin(9) + qJ(2) * t389) * qJD(1), t376, -(-t393 + t403) * qJD(4) + t396 * t397 - t372 * t398, -t398, 0, 0; 0, 0, t400 * qJD(4) - t396 * (-t359 * t391 + (t350 * t357 - t354 * t387) * t352) * qJD(3) + t372 * t331, t331, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:03
	% EndTime: 2019-10-10 01:04:03
	% DurationCPUTime: 0.74s
	% Computational Cost: add. (502->90), mult. (1462->155), div. (0->0), fcn. (1557->14), ass. (0->73)
	t422 = cos(pkin(12));
	t424 = cos(pkin(6));
	t419 = sin(pkin(12));
	t427 = sin(qJ(1));
	t458 = t427 * t419;
	t448 = t424 * t458;
	t429 = cos(qJ(1));
	t453 = qJD(1) * t429;
	t399 = -qJD(1) * t448 + t422 * t453;
	t420 = sin(pkin(7));
	t421 = sin(pkin(6));
	t461 = t421 * t429;
	t447 = t420 * t461;
	t484 = -qJD(3) * t447 + t399;
	t423 = cos(pkin(7));
	t455 = t429 * t422;
	t436 = t424 * t455 - t458;
	t478 = t436 * t423;
	t483 = -t478 + t447;
	t456 = t429 * t419;
	t457 = t427 * t422;
	t435 = t424 * t457 + t456;
	t398 = t435 * qJD(1);
	t402 = t424 * t456 + t457;
	t426 = sin(qJ(3));
	t428 = cos(qJ(3));
	t462 = t421 * t427;
	t445 = qJD(1) * t462;
	t377 = (-qJD(3) * t402 - t398 * t423 + t420 * t445) * t426 + (qJD(3) * t478 + t484) * t428;
	t468 = t398 * t420;
	t389 = t423 * t445 + t468;
	t417 = pkin(13) + qJ(5);
	t415 = sin(t417);
	t416 = cos(t417);
	t482 = t377 * t415 - t389 * t416;
	t481 = -t377 * t416 - t389 * t415;
	t477 = t436 * t420 + t423 * t461;
	t459 = t423 * t428;
	t464 = t420 * t424;
	t476 = (-t419 * t426 + t422 * t459) * t421 + t428 * t464;
	t460 = t423 * t426;
	t439 = -t402 * t428 - t436 * t460;
	t382 = t426 * t447 + t439;
	t475 = t382 * t416 + t415 * t477;
	t474 = -t382 * t415 + t416 * t477;
	t472 = t402 * t426 + t483 * t428;
	t454 = qJD(1) * t428;
	t463 = t421 * t420;
	t444 = t454 * t463;
	t471 = t439 * qJD(3) - t398 * t459 - t484 * t426 + t427 * t444;
	t470 = sin(pkin(13)) * pkin(4);
	t469 = r_i_i_C(3) + pkin(10) + qJ(4);
	t451 = t421 * qJD(2);
	t442 = t415 * r_i_i_C(1) + t416 * r_i_i_C(2);
	t414 = cos(pkin(13)) * pkin(4) + pkin(3);
	t440 = -t416 * r_i_i_C(1) + t415 * r_i_i_C(2) - t414;
	t437 = t420 * t462 - t423 * t435;
	t434 = qJD(5) * t442;
	t404 = -t448 + t455;
	t383 = -t404 * t426 + t437 * t428;
	t384 = t404 * t428 + t437 * t426;
	t391 = t426 * t464 + (t419 * t428 + t422 * t460) * t421;
	t400 = -t422 * t463 + t424 * t423;
	t397 = t402 * qJD(1);
	t394 = t420 * t435 + t423 * t462;
	t387 = t477 * qJD(1);
	t386 = t391 * qJD(3);
	t385 = t476 * qJD(3);
	t375 = t483 * t426 * qJD(1) + t383 * qJD(3) - t397 * t428;
	t374 = t384 * qJD(3) - t397 * t426 - t429 * t444 + t454 * t478;
	t373 = t375 * t416 + t387 * t415 + (-t384 * t415 + t394 * t416) * qJD(5);
	t372 = -t375 * t415 + t387 * t416 + (-t384 * t416 - t394 * t415) * qJD(5);
	t1 = [t481 * r_i_i_C(1) + t482 * r_i_i_C(2) - t377 * t414 - t472 * qJD(4) - t389 * t470 - t399 * pkin(2) - pkin(9) * t468 + t429 * t451 + t469 * t471 + (t474 * r_i_i_C(1) - t475 * r_i_i_C(2)) * qJD(5) + (-t429 * pkin(1) + (-pkin(9) * t423 - qJ(2)) * t462) * qJD(1), t421 * t453, t384 * qJD(4) + t440 * t374 + t469 * t375 - t383 * t434, t374, t372 * r_i_i_C(1) - t373 * r_i_i_C(2), 0; t387 * t470 + t427 * t451 - t397 * pkin(2) + t373 * r_i_i_C(1) + t372 * r_i_i_C(2) - t383 * qJD(4) + t375 * t414 + t469 * t374 + (-t427 * pkin(1) + pkin(9) * t477 + qJ(2) * t461) * qJD(1), t445, -qJD(4) * t382 + t469 * t377 + t472 * t434 - t440 * t471, -t471, -t482 * r_i_i_C(1) + t481 * r_i_i_C(2) + (t475 * r_i_i_C(1) + t474 * r_i_i_C(2)) * qJD(5), 0; 0, 0, t391 * qJD(4) + t469 * t385 + t440 * t386 - t476 * t434, t386, -t442 * t385 + ((-t391 * t416 - t400 * t415) * r_i_i_C(1) + (t391 * t415 - t400 * t416) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:06
	% EndTime: 2019-10-10 01:04:08
	% DurationCPUTime: 1.67s
	% Computational Cost: add. (1421->136), mult. (3937->227), div. (0->0), fcn. (4375->16), ass. (0->89)
	t601 = cos(qJ(1));
	t654 = cos(pkin(12));
	t656 = cos(pkin(6));
	t632 = t656 * t654;
	t652 = sin(pkin(12));
	t659 = sin(qJ(1));
	t613 = t601 * t652 + t659 * t632;
	t571 = t613 * qJD(1);
	t630 = t656 * t652;
	t576 = t601 * t654 - t659 * t630;
	t572 = t576 * qJD(1);
	t599 = sin(qJ(3));
	t611 = -t601 * t630 - t659 * t654;
	t596 = sin(pkin(6));
	t653 = sin(pkin(7));
	t641 = t596 * t653;
	t623 = t659 * t641;
	t655 = cos(pkin(7));
	t660 = cos(qJ(3));
	t612 = t601 * t632 - t659 * t652;
	t636 = t601 * t641;
	t621 = t660 * t636;
	t635 = t655 * t660;
	t666 = t612 * t635 - t621;
	t538 = (qJD(1) * t623 + qJD(3) * t611 - t655 * t571) * t599 + t572 * t660 + t666 * qJD(3);
	t634 = t655 * t659;
	t624 = t596 * t634;
	t642 = t571 * t653;
	t560 = qJD(1) * t624 + t642;
	t594 = pkin(13) + qJ(5);
	t592 = sin(t594);
	t593 = cos(t594);
	t648 = t596 * t601;
	t665 = t612 * t653 + t655 * t648;
	t584 = t599 * t636;
	t667 = t612 * t655;
	t674 = -t599 * t667 + t611 * t660 + t584;
	t664 = -t592 * t674 + t593 * t665;
	t534 = t664 * qJD(5) - t538 * t593 - t560 * t592;
	t617 = t660 * t623;
	t539 = qJD(1) * t617 + t674 * qJD(3) - t571 * t635 - t572 * t599;
	t598 = sin(qJ(6));
	t600 = cos(qJ(6));
	t681 = t534 * t598 - t539 * t600;
	t680 = t534 * t600 + t539 * t598;
	t544 = -t592 * t665 - t593 * t674;
	t550 = -t599 * t611 - t666;
	t679 = t544 * t598 - t550 * t600;
	t678 = t544 * t600 + t550 * t598;
	t677 = -t544 * qJD(5) - t538 * t592 + t560 * t593;
	t661 = r_i_i_C(3) + pkin(11);
	t607 = t613 * t655;
	t555 = t576 * t660 + (-t607 + t623) * t599;
	t565 = t613 * t653 + t624;
	t663 = -t555 * t592 + t565 * t593;
	t622 = qJD(6) * (t598 * r_i_i_C(1) + t600 * r_i_i_C(2));
	t662 = -t576 * t599 - t660 * t607 + t617;
	t629 = t655 * t654;
	t631 = t656 * t653;
	t561 = -t660 * t631 + (t599 * t652 - t629 * t660) * t596;
	t602 = t665 * qJD(1);
	t658 = pkin(4) * sin(pkin(13));
	t647 = qJD(6) * t598;
	t646 = qJD(6) * t600;
	t645 = t596 * qJD(2);
	t644 = qJD(1) * t648;
	t643 = qJD(1) * t659;
	t547 = t555 * t593 + t565 * t592;
	t562 = t599 * t631 + (t599 * t629 + t660 * t652) * t596;
	t573 = -t654 * t641 + t656 * t655;
	t549 = t562 * t593 + t573 * t592;
	t627 = -t562 * t592 + t573 * t593;
	t626 = t600 * r_i_i_C(1) - t598 * r_i_i_C(2) + pkin(5);
	t591 = cos(pkin(13)) * pkin(4) + pkin(3);
	t610 = -t661 * t592 - t626 * t593 - t591;
	t605 = qJD(1) * t667;
	t603 = t593 * t622 + (t626 * t592 - t661 * t593) * qJD(5);
	t597 = -pkin(10) - qJ(4);
	t570 = t611 * qJD(1);
	t558 = t562 * qJD(3);
	t557 = t561 * qJD(3);
	t542 = t627 * qJD(5) - t557 * t593;
	t536 = qJD(1) * t584 + t662 * qJD(3) + t570 * t660 - t599 * t605;
	t535 = -qJD(1) * t621 + t555 * qJD(3) + t570 * t599 + t660 * t605;
	t530 = t663 * qJD(5) + t536 * t593 + t592 * t602;
	t529 = t547 * qJD(5) + t536 * t592 - t593 * t602;
	t528 = t530 * t600 + t535 * t598 + (-t547 * t598 - t600 * t662) * qJD(6);
	t527 = -t530 * t598 + t535 * t600 + (-t547 * t600 + t598 * t662) * qJD(6);
	t1 = [t680 * r_i_i_C(1) - t681 * r_i_i_C(2) + t534 * pkin(5) - t538 * t591 - t539 * t597 - t550 * qJD(4) - t560 * t658 - t572 * pkin(2) - pkin(9) * t642 + t601 * t645 + t661 * t677 + (t679 * r_i_i_C(1) + t678 * r_i_i_C(2)) * qJD(6) + (-t601 * pkin(1) + (-pkin(9) * t634 - t659 * qJ(2)) * t596) * qJD(1), t644, (t536 * t598 + t555 * t646) * r_i_i_C(1) + (t536 * t600 - t555 * t647) * r_i_i_C(2) - t536 * t597 + t555 * qJD(4) + t610 * t535 - t603 * t662, t535, -t626 * t529 + t661 * t530 - t663 * t622, t527 * r_i_i_C(1) - t528 * r_i_i_C(2); -pkin(1) * t643 + t570 * pkin(2) + t530 * pkin(5) + t528 * r_i_i_C(1) + t527 * r_i_i_C(2) + qJ(2) * t644 - t662 * qJD(4) + t661 * t529 - t535 * t597 + t536 * t591 + t659 * t645 + (pkin(9) + t658) * t602, t596 * t643, (t538 * t598 - t646 * t674) * r_i_i_C(1) + (t538 * t600 + t647 * t674) * r_i_i_C(2) - t538 * t597 - t674 * qJD(4) - t610 * t539 + t603 * t550, -t539, -t534 * t661 + t664 * t622 + t626 * t677, t681 * r_i_i_C(1) + t680 * r_i_i_C(2) + (-t678 * r_i_i_C(1) + t679 * r_i_i_C(2)) * qJD(6); 0, 0, (-t557 * t598 + t562 * t646) * r_i_i_C(1) + (-t557 * t600 - t562 * t647) * r_i_i_C(2) + t557 * t597 + t562 * qJD(4) + t610 * t558 + t603 * t561, t558, t661 * t542 - t627 * t622 + t626 * (-t549 * qJD(5) + t557 * t592), (-t542 * t598 + t558 * t600) * r_i_i_C(1) + (-t542 * t600 - t558 * t598) * r_i_i_C(2) + ((-t549 * t600 - t561 * t598) * r_i_i_C(1) + (t549 * t598 - t561 * t600) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end