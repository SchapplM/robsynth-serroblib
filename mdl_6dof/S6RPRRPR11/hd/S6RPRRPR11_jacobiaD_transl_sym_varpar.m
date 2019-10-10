% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR11_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR11_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR11_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
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
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:21
	% DurationCPUTime: 0.11s
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
	% StartTime: 2019-10-10 01:41:21
	% EndTime: 2019-10-10 01:41:21
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
	% StartTime: 2019-10-10 01:41:23
	% EndTime: 2019-10-10 01:41:23
	% DurationCPUTime: 0.65s
	% Computational Cost: add. (393->78), mult. (1301->142), div. (0->0), fcn. (1376->12), ass. (0->64)
	t408 = sin(qJ(3));
	t402 = sin(pkin(7));
	t403 = sin(pkin(6));
	t412 = cos(qJ(1));
	t442 = t403 * t412;
	t431 = t402 * t442;
	t405 = cos(pkin(7));
	t404 = cos(pkin(12));
	t406 = cos(pkin(6));
	t440 = t406 * t412;
	t401 = sin(pkin(12));
	t409 = sin(qJ(1));
	t445 = t401 * t409;
	t421 = t404 * t440 - t445;
	t454 = t421 * t405;
	t424 = -t454 + t431;
	t463 = t424 * t408;
	t432 = t406 * t445;
	t438 = qJD(1) * t412;
	t389 = -qJD(1) * t432 + t404 * t438;
	t439 = t409 * t404;
	t392 = t401 * t440 + t439;
	t411 = cos(qJ(3));
	t420 = t401 * t412 + t406 * t439;
	t388 = t420 * qJD(1);
	t443 = t403 * t409;
	t429 = qJD(1) * t443;
	t418 = -t388 * t405 + t402 * t429;
	t428 = t411 * t431;
	t368 = -qJD(3) * t428 + (-qJD(3) * t392 + t418) * t408 - (-qJD(3) * t454 - t389) * t411;
	t449 = t388 * t402;
	t380 = t405 * t429 + t449;
	t407 = sin(qJ(4));
	t410 = cos(qJ(4));
	t462 = t368 * t407 - t380 * t410;
	t461 = -t368 * t410 - t380 * t407;
	t373 = -t392 * t411 + t463;
	t453 = t421 * t402 + t405 * t442;
	t460 = -t373 * t407 + t453 * t410;
	t459 = t373 * t410 + t453 * t407;
	t441 = t404 * t405;
	t444 = t402 * t406;
	t452 = (-t401 * t408 + t411 * t441) * t403 + t411 * t444;
	t382 = (t401 * t411 + t408 * t441) * t403 + t408 * t444;
	t450 = r_i_i_C(3) + pkin(10);
	t436 = t403 * qJD(2);
	t427 = t407 * r_i_i_C(1) + t410 * r_i_i_C(2);
	t425 = t410 * r_i_i_C(1) - t407 * r_i_i_C(2) + pkin(3);
	t423 = t402 * t443 - t405 * t420;
	t419 = qJD(4) * t427;
	t394 = t404 * t412 - t432;
	t415 = -t394 * t408 + t423 * t411;
	t375 = t394 * t411 + t423 * t408;
	t413 = t373 * qJD(3) - t389 * t408 + t418 * t411;
	t390 = -t402 * t403 * t404 + t405 * t406;
	t387 = t392 * qJD(1);
	t385 = t402 * t420 + t405 * t443;
	t378 = t453 * qJD(1);
	t376 = t452 * qJD(3);
	t366 = qJD(1) * t463 + t415 * qJD(3) - t387 * t411;
	t365 = t375 * qJD(3) - t387 * t408 + (t411 * t454 - t428) * qJD(1);
	t364 = t366 * t410 + t378 * t407 + (-t375 * t407 + t385 * t410) * qJD(4);
	t363 = -t366 * t407 + t378 * t410 + (-t375 * t410 - t385 * t407) * qJD(4);
	t1 = [t461 * r_i_i_C(1) + t462 * r_i_i_C(2) - t368 * pkin(3) - t389 * pkin(2) - pkin(9) * t449 + t412 * t436 + t450 * t413 + (t460 * r_i_i_C(1) - t459 * r_i_i_C(2)) * qJD(4) + (-t412 * pkin(1) + (-pkin(9) * t405 - qJ(2)) * t443) * qJD(1), t403 * t438, -t425 * t365 + t450 * t366 - t415 * t419, r_i_i_C(1) * t363 - r_i_i_C(2) * t364, 0, 0; t409 * t436 - t387 * pkin(2) + t366 * pkin(3) + t364 * r_i_i_C(1) + t363 * r_i_i_C(2) + t450 * t365 + (-t409 * pkin(1) + pkin(9) * t453 + qJ(2) * t442) * qJD(1), t429, t450 * t368 - (-t392 * t408 - t424 * t411) * t419 + t425 * t413, -t462 * r_i_i_C(1) + t461 * r_i_i_C(2) + (t459 * r_i_i_C(1) + t460 * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, -t425 * t382 * qJD(3) + t450 * t376 - t452 * t419, -t427 * t376 + ((-t382 * t410 - t390 * t407) * r_i_i_C(1) + (t382 * t407 - t390 * t410) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:26
	% EndTime: 2019-10-10 01:41:27
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (917->99), mult. (2979->171), div. (0->0), fcn. (3248->14), ass. (0->75)
	t518 = sin(qJ(1));
	t521 = cos(qJ(1));
	t511 = sin(pkin(12));
	t571 = cos(pkin(6));
	t553 = t511 * t571;
	t570 = cos(pkin(12));
	t503 = -t518 * t553 + t521 * t570;
	t498 = t503 * qJD(1);
	t517 = sin(qJ(3));
	t520 = cos(qJ(3));
	t529 = -t518 * t570 - t521 * t553;
	t542 = t571 * t570;
	t530 = t521 * t511 + t518 * t542;
	t497 = t530 * qJD(1);
	t512 = sin(pkin(7));
	t515 = cos(pkin(7));
	t513 = sin(pkin(6));
	t560 = qJD(1) * t518;
	t555 = t513 * t560;
	t532 = -t497 * t515 + t512 * t555;
	t558 = qJD(3) * t520;
	t563 = t513 * t512;
	t548 = t558 * t563;
	t528 = -t518 * t511 + t521 * t542;
	t566 = t528 * t515;
	t475 = (qJD(3) * t529 + t532) * t517 - (-qJD(3) * t566 - t498) * t520 - t521 * t548;
	t516 = sin(qJ(4));
	t519 = cos(qJ(4));
	t567 = t497 * t512;
	t533 = t515 * t555 + t567;
	t561 = t513 * t521;
	t535 = t512 * t561 - t566;
	t482 = t535 * t517 + t520 * t529;
	t492 = t528 * t512 + t515 * t561;
	t581 = t482 * t519 + t492 * t516;
	t586 = t581 * qJD(4) - t475 * t516 + t533 * t519;
	t582 = t482 * t516 - t492 * t519;
	t585 = t582 * qJD(4) + t475 * t519 + t533 * t516;
	t572 = r_i_i_C(3) + qJ(5);
	t526 = t528 * qJD(1);
	t554 = qJD(1) * t561;
	t576 = t512 * t554 - t515 * t526;
	t562 = t513 * t518;
	t564 = t530 * t515;
	t534 = t512 * t562 - t564;
	t484 = t503 * t520 + t534 * t517;
	t494 = t512 * t530 + t515 * t562;
	t575 = -t484 * t516 + t494 * t519;
	t574 = pkin(9) * t515 + qJ(2);
	t550 = t571 * t512;
	t552 = t515 * t570;
	t491 = (t511 * t520 + t517 * t552) * t513 + t517 * t550;
	t510 = sin(pkin(13));
	t514 = cos(pkin(13));
	t537 = t514 * r_i_i_C(1) - t510 * r_i_i_C(2) + pkin(4);
	t524 = t572 * t516 + t537 * t519 + pkin(3);
	t559 = qJD(3) * t517;
	t557 = t513 * qJD(2);
	t545 = t520 * t550;
	t544 = t520 * t552;
	t539 = t484 * t519 + t494 * t516;
	t499 = t571 * t515 - t570 * t563;
	t538 = t491 * t519 + t499 * t516;
	t536 = t510 * r_i_i_C(1) + t514 * r_i_i_C(2) + pkin(10);
	t523 = t516 * qJD(5) + (-t537 * t516 + t572 * t519) * qJD(4);
	t522 = t492 * qJD(1);
	t476 = t482 * qJD(3) - t498 * t517 + t532 * t520;
	t496 = t529 * qJD(1);
	t487 = t513 * t511 * t559 + (-t513 * t544 - t545) * qJD(3);
	t478 = t538 * qJD(4) - t487 * t516;
	t473 = t496 * t520 - t503 * t559 + t576 * t517 + t518 * t548 - t558 * t564;
	t472 = t484 * qJD(3) + t496 * t517 - t576 * t520;
	t467 = t575 * qJD(4) + t473 * t519 + t516 * t522;
	t466 = t539 * qJD(4) + t473 * t516 - t519 * t522;
	t1 = [(t476 * t510 - t514 * t585) * r_i_i_C(1) + (t476 * t514 + t510 * t585) * r_i_i_C(2) - t585 * pkin(4) + t582 * qJD(5) - t475 * pkin(3) + t476 * pkin(10) - t498 * pkin(2) - pkin(9) * t567 + t521 * t557 + t572 * t586 + (-t521 * pkin(1) - t574 * t562) * qJD(1), t554, t536 * t473 + t523 * (-t503 * t517 + t534 * t520) - t524 * t472, t539 * qJD(5) - t537 * t466 + t572 * t467, t466, 0; (t467 * t514 + t472 * t510) * r_i_i_C(1) + (-t467 * t510 + t472 * t514) * r_i_i_C(2) + t467 * pkin(4) - t575 * qJD(5) + t473 * pkin(3) + t472 * pkin(10) + t496 * pkin(2) + t512 * pkin(9) * t526 - pkin(1) * t560 + t518 * t557 + t574 * t554 + t572 * t466, t555, t536 * t475 + t523 * (t517 * t529 - t535 * t520) + t524 * t476, -t581 * qJD(5) + t537 * t586 + t572 * t585, -t586, 0; 0, 0, -t536 * t487 + t523 * (t545 + (-t511 * t517 + t544) * t513) - t524 * t491 * qJD(3), t538 * qJD(5) + t572 * (-t487 * t519 + (-t491 * t516 + t499 * t519) * qJD(4)) - t537 * t478, t478, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:27
	% EndTime: 2019-10-10 01:41:28
	% DurationCPUTime: 1.58s
	% Computational Cost: add. (1421->127), mult. (4224->209), div. (0->0), fcn. (4706->16), ass. (0->89)
	t599 = sin(qJ(3));
	t596 = sin(pkin(6));
	t659 = sin(pkin(7));
	t650 = t596 * t659;
	t668 = cos(qJ(1));
	t628 = t668 * t650;
	t660 = cos(pkin(12));
	t662 = cos(pkin(6));
	t640 = t662 * t660;
	t658 = sin(pkin(12));
	t666 = sin(qJ(1));
	t575 = -t640 * t668 + t666 * t658;
	t661 = cos(pkin(7));
	t649 = t661 * t575;
	t638 = t662 * t658;
	t576 = t638 * t668 + t666 * t660;
	t667 = cos(qJ(3));
	t653 = t576 * t667;
	t555 = t599 * (t628 + t649) - t653;
	t652 = t596 * t668;
	t566 = t575 * t659 - t661 * t652;
	t598 = sin(qJ(4));
	t600 = cos(qJ(4));
	t547 = t555 * t600 - t566 * t598;
	t616 = t667 * t628;
	t674 = -t667 * t649 - t616;
	t552 = t576 * t599 - t674;
	t594 = pkin(13) + qJ(6);
	t592 = sin(t594);
	t593 = cos(t594);
	t683 = -t547 * t592 - t552 * t593;
	t682 = t547 * t593 - t552 * t592;
	t577 = -t666 * t638 + t668 * t660;
	t573 = t577 * qJD(1);
	t626 = t666 * t650;
	t608 = t666 * t640 + t668 * t658;
	t654 = t608 * qJD(1);
	t671 = t654 * t661;
	t542 = -t674 * qJD(3) - t573 * t667 + (-qJD(1) * t626 + t576 * qJD(3) + t671) * t599;
	t646 = t661 * t666;
	t627 = t596 * t646;
	t634 = t654 * t659;
	t607 = qJD(1) * t627 + t634;
	t681 = qJD(4) * t547 + t542 * t598 + t600 * t607;
	t633 = t555 * t598 + t566 * t600;
	t534 = qJD(4) * t633 - t542 * t600 + t598 * t607;
	t615 = t667 * t626;
	t622 = t599 * t628;
	t539 = t667 * t671 - qJD(1) * t615 + t573 * t599 - (t599 * t649 + t622 - t653) * qJD(3);
	t672 = pkin(10) + sin(pkin(13)) * pkin(5);
	t664 = r_i_i_C(3) + pkin(11) + qJ(5);
	t606 = t608 * t661;
	t557 = t577 * t667 + (-t606 + t626) * t599;
	t567 = t608 * t659 + t627;
	t548 = -t557 * t598 + t567 * t600;
	t642 = t592 * r_i_i_C(1) + t593 * r_i_i_C(2);
	t623 = qJD(6) * t642;
	t670 = -t577 * t599 - t667 * t606 + t615;
	t637 = t661 * t660;
	t639 = t662 * t659;
	t563 = -t667 * t639 + (t599 * t658 - t637 * t667) * t596;
	t602 = t566 * qJD(1);
	t651 = qJD(1) * t666;
	t648 = qJD(1) * t652;
	t643 = t593 * r_i_i_C(1) - t592 * r_i_i_C(2);
	t549 = t557 * t600 + t567 * t598;
	t564 = t599 * t639 + (t599 * t637 + t658 * t667) * t596;
	t574 = -t650 * t660 + t661 * t662;
	t551 = t564 * t600 + t574 * t598;
	t631 = -t564 * t598 + t574 * t600;
	t591 = cos(pkin(13)) * pkin(5) + pkin(4);
	t625 = t591 + t643;
	t624 = qJD(6) * t643;
	t617 = t642 + t672;
	t609 = -t598 * t664 - t600 * t625 - pkin(3);
	t604 = qJD(1) * t649;
	t601 = -t598 * qJD(5) + t600 * t623 + (t625 * t598 - t664 * t600) * qJD(4);
	t572 = t576 * qJD(1);
	t561 = t564 * qJD(3);
	t560 = t563 * qJD(3);
	t544 = qJD(4) * t631 - t560 * t600;
	t543 = qJD(4) * t551 - t560 * t598;
	t538 = qJD(1) * t622 + t670 * qJD(3) - t572 * t667 + t599 * t604;
	t537 = -qJD(1) * t616 + qJD(3) * t557 - t572 * t599 - t604 * t667;
	t532 = t548 * qJD(4) + t538 * t600 - t598 * t602;
	t531 = qJD(4) * t549 + t538 * t598 + t600 * t602;
	t530 = t532 * t593 + t537 * t592 + (-t549 * t592 - t593 * t670) * qJD(6);
	t529 = -t532 * t592 + t537 * t593 + (-t549 * t593 + t592 * t670) * qJD(6);
	t1 = [t633 * qJD(5) + t542 * pkin(3) - t573 * pkin(2) - pkin(9) * t634 + qJD(2) * t652 - t625 * t534 - t617 * t539 + t664 * t681 + (t683 * r_i_i_C(1) - t682 * r_i_i_C(2)) * qJD(6) + (-t668 * pkin(1) + (-pkin(9) * t646 - qJ(2) * t666) * t596) * qJD(1), t648, t537 * t609 + t538 * t617 + t557 * t624 - t601 * t670, t549 * qJD(5) - t531 * t625 + t532 * t664 - t548 * t623, t531, t529 * r_i_i_C(1) - t530 * r_i_i_C(2); qJD(2) * t596 * t666 - pkin(1) * t651 - t572 * pkin(2) + t538 * pkin(3) - pkin(9) * t602 + t530 * r_i_i_C(1) + t529 * r_i_i_C(2) + qJ(2) * t648 - t548 * qJD(5) + t664 * t531 + t532 * t591 + t672 * t537, t596 * t651, t539 * t609 - t542 * t617 + t552 * t601 - t555 * t624, -qJD(5) * t547 + t534 * t664 - t623 * t633 + t625 * t681, -t681, (-t534 * t592 + t539 * t593) * r_i_i_C(1) + (-t534 * t593 - t539 * t592) * r_i_i_C(2) + (t682 * r_i_i_C(1) + t683 * r_i_i_C(2)) * qJD(6); 0, 0, -t560 * t617 + t561 * t609 + t563 * t601 + t564 * t624, t551 * qJD(5) - t543 * t625 + t544 * t664 - t623 * t631, t543, (-t544 * t592 + t561 * t593) * r_i_i_C(1) + (-t544 * t593 - t561 * t592) * r_i_i_C(2) + ((-t551 * t593 - t563 * t592) * r_i_i_C(1) + (t551 * t592 - t563 * t593) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end