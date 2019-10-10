% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:12
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPRR12_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR12_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:18
	% EndTime: 2019-10-10 12:12:18
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:18
	% EndTime: 2019-10-10 12:12:18
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:19
	% EndTime: 2019-10-10 12:12:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (27->13), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t158 = sin(qJ(2));
	t159 = sin(qJ(1));
	t171 = t158 * t159;
	t161 = cos(qJ(1));
	t170 = t158 * t161;
	t160 = cos(qJ(2));
	t169 = t159 * t160;
	t168 = t160 * t161;
	t156 = sin(pkin(6));
	t167 = qJD(1) * t156;
	t166 = qJD(2) * t156;
	t157 = cos(pkin(6));
	t165 = -t157 * t168 + t171;
	t164 = t157 * t169 + t170;
	t163 = t157 * t170 + t169;
	t162 = t157 * t171 - t168;
	t155 = t162 * qJD(1) + t165 * qJD(2);
	t154 = t164 * qJD(1) + t163 * qJD(2);
	t153 = t163 * qJD(1) + t164 * qJD(2);
	t152 = t165 * qJD(1) + t162 * qJD(2);
	t1 = [t155, t152, 0, 0, 0, 0; -t153, -t154, 0, 0, 0, 0; 0, -t158 * t166, 0, 0, 0, 0; t154, t153, 0, 0, 0, 0; t152, t155, 0, 0, 0, 0; 0, -t160 * t166, 0, 0, 0, 0; -t159 * t167, 0, 0, 0, 0, 0; t161 * t167, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:20
	% EndTime: 2019-10-10 12:12:20
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (94->35), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->37)
	t284 = cos(pkin(6));
	t287 = sin(qJ(1));
	t286 = sin(qJ(2));
	t307 = t287 * t286;
	t298 = t284 * t307;
	t302 = qJD(2) * t286;
	t289 = cos(qJ(2));
	t290 = cos(qJ(1));
	t304 = t290 * t289;
	t275 = -qJD(1) * t298 - t287 * t302 + (qJD(2) * t284 + qJD(1)) * t304;
	t305 = t290 * t286;
	t306 = t287 * t289;
	t277 = t284 * t305 + t306;
	t285 = sin(qJ(3));
	t288 = cos(qJ(3));
	t283 = sin(pkin(6));
	t303 = qJD(1) * t283;
	t297 = t287 * t303;
	t308 = t283 * t290;
	t311 = (-t277 * t288 + t285 * t308) * qJD(3) - t275 * t285 + t288 * t297;
	t310 = t283 * t285;
	t309 = t283 * t288;
	t301 = qJD(3) * t285;
	t300 = qJD(3) * t288;
	t299 = qJD(3) * t289;
	t296 = t290 * t303;
	t295 = t283 * qJD(2) * t289;
	t276 = t284 * t304 - t307;
	t278 = -t284 * t306 - t305;
	t293 = t298 - t304;
	t291 = -t275 * t288 + t300 * t308 + (qJD(3) * t277 - t297) * t285;
	t274 = t278 * qJD(1) - t277 * qJD(2);
	t273 = -t277 * qJD(1) + t278 * qJD(2);
	t272 = -t276 * qJD(1) + t293 * qJD(2);
	t271 = t285 * t296 + t273 * t288 + (t285 * t293 + t287 * t309) * qJD(3);
	t270 = t288 * t296 - t273 * t285 + (-t287 * t310 + t288 * t293) * qJD(3);
	t1 = [t291, t272 * t288 - t278 * t301, t270, 0, 0, 0; t271, t274 * t288 - t276 * t301, t311, 0, 0, 0; 0, (-t285 * t299 - t288 * t302) * t283, -t285 * t295 + (-t284 * t285 - t286 * t309) * qJD(3), 0, 0, 0; -t311, -t272 * t285 - t278 * t300, -t271, 0, 0, 0; t270, -t274 * t285 - t276 * t300, t291, 0, 0, 0; 0, (t285 * t302 - t288 * t299) * t283, -t288 * t295 + (-t284 * t288 + t286 * t310) * qJD(3), 0, 0, 0; t274, t273, 0, 0, 0, 0; -t272, t275, 0, 0, 0, 0; 0, t295, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:21
	% EndTime: 2019-10-10 12:12:21
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (160->48), mult. (518->103), div. (0->0), fcn. (536->10), ass. (0->46)
	t390 = sin(qJ(1));
	t387 = cos(pkin(6));
	t398 = qJD(2) * t387 + qJD(1);
	t389 = sin(qJ(2));
	t414 = t390 * t389;
	t404 = t387 * t414;
	t409 = qJD(2) * t389;
	t392 = cos(qJ(2));
	t393 = cos(qJ(1));
	t411 = t393 * t392;
	t371 = -qJD(1) * t404 - t390 * t409 + t398 * t411;
	t412 = t393 * t389;
	t413 = t390 * t392;
	t374 = t387 * t412 + t413;
	t388 = sin(qJ(3));
	t391 = cos(qJ(3));
	t385 = sin(pkin(6));
	t410 = qJD(1) * t385;
	t402 = t390 * t410;
	t416 = t385 * t393;
	t367 = (t374 * t388 + t391 * t416) * qJD(3) - t371 * t391 - t388 * t402;
	t417 = t385 * t390;
	t415 = t389 * t391;
	t408 = qJD(2) * t392;
	t407 = qJD(3) * t388;
	t406 = qJD(3) * t391;
	t405 = qJD(3) * t392;
	t403 = t387 * t411;
	t401 = t393 * t410;
	t400 = t385 * t408;
	t399 = t388 * t405;
	t375 = -t387 * t413 - t412;
	t368 = -qJD(1) * t403 - t393 * t408 + t398 * t414;
	t396 = -t368 * t391 + t375 * t407;
	t370 = t375 * qJD(1) - t374 * qJD(2);
	t373 = t403 - t414;
	t395 = -t370 * t391 + t373 * t407;
	t366 = -t371 * t388 - t374 * t406 + t391 * t402 + t407 * t416;
	t386 = cos(pkin(12));
	t384 = sin(pkin(12));
	t376 = -t404 + t411;
	t372 = -t388 * t400 + (-t385 * t415 - t387 * t388) * qJD(3);
	t369 = -t374 * qJD(1) + t375 * qJD(2);
	t365 = t388 * t401 + t369 * t391 + (-t376 * t388 + t391 * t417) * qJD(3);
	t364 = t369 * t388 - t391 * t401 + (t376 * t391 + t388 * t417) * qJD(3);
	t1 = [t367 * t386 + t370 * t384, t369 * t384 - t396 * t386, -t364 * t386, 0, 0, 0; t365 * t386 - t368 * t384, t371 * t384 - t395 * t386, t366 * t386, 0, 0, 0; 0, (-t386 * t399 + (t384 * t392 - t386 * t415) * qJD(2)) * t385, t372 * t386, 0, 0, 0; -t367 * t384 + t370 * t386, t369 * t386 + t396 * t384, t364 * t384, 0, 0, 0; -t365 * t384 - t368 * t386, t371 * t386 + t395 * t384, -t366 * t384, 0, 0, 0; 0, (t384 * t399 + (t384 * t415 + t386 * t392) * qJD(2)) * t385, -t372 * t384, 0, 0, 0; t366, t368 * t388 + t375 * t406, t365, 0, 0, 0; t364, t370 * t388 + t373 * t406, -t367, 0, 0, 0; 0, (-t388 * t409 + t391 * t405) * t385, t391 * t400 + (-t385 * t388 * t389 + t387 * t391) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:22
	% EndTime: 2019-10-10 12:12:22
	% DurationCPUTime: 0.61s
	% Computational Cost: add. (371->74), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->67)
	t473 = sin(qJ(1));
	t470 = cos(pkin(6));
	t485 = qJD(2) * t470 + qJD(1);
	t472 = sin(qJ(2));
	t506 = t473 * t472;
	t493 = t470 * t506;
	t501 = qJD(2) * t472;
	t475 = cos(qJ(2));
	t476 = cos(qJ(1));
	t503 = t476 * t475;
	t444 = -qJD(1) * t493 - t473 * t501 + t485 * t503;
	t471 = sin(qJ(3));
	t474 = cos(qJ(3));
	t504 = t476 * t472;
	t505 = t473 * t475;
	t455 = t470 * t504 + t505;
	t469 = sin(pkin(6));
	t507 = t469 * t476;
	t480 = t455 * t471 + t474 * t507;
	t502 = qJD(1) * t469;
	t490 = t473 * t502;
	t440 = t480 * qJD(3) - t444 * t474 - t471 * t490;
	t456 = t470 * t505 + t504;
	t443 = t456 * qJD(1) + t455 * qJD(2);
	t492 = t471 * t507;
	t449 = -t455 * t474 + t492;
	t491 = t470 * t503;
	t454 = -t491 + t506;
	t468 = pkin(12) + qJ(5);
	t466 = sin(t468);
	t467 = cos(t468);
	t519 = t440 * t467 - t443 * t466 + (-t449 * t466 - t454 * t467) * qJD(5);
	t518 = (t449 * t467 - t454 * t466) * qJD(5) + t440 * t466 + t443 * t467;
	t497 = qJD(3) * t475;
	t515 = (qJD(2) * t474 - qJD(5)) * t472 + t471 * t497;
	t510 = t469 * t471;
	t509 = t469 * t474;
	t508 = t469 * t475;
	t500 = qJD(2) * t475;
	t499 = qJD(3) * t471;
	t498 = qJD(3) * t474;
	t496 = qJD(5) * t466;
	t495 = qJD(5) * t467;
	t494 = qJD(5) * t474;
	t489 = t476 * t502;
	t488 = t469 * t501;
	t487 = t469 * t500;
	t442 = -t455 * qJD(1) - t456 * qJD(2);
	t483 = t456 * t494 + t442;
	t482 = t454 * t494 + t444;
	t481 = (qJD(2) - t494) * t475;
	t457 = -t493 + t503;
	t450 = -t457 * t471 + t473 * t509;
	t451 = t457 * t474 + t473 * t510;
	t453 = t470 * t471 + t472 * t509;
	t452 = t470 * t474 - t472 * t510;
	t438 = qJD(3) * t492 - t444 * t471 - t455 * t498 + t474 * t490;
	t441 = -qJD(1) * t491 - t476 * t500 + t485 * t506;
	t478 = qJD(5) * t457 + t441 * t474 + t456 * t499;
	t477 = qJD(5) * t455 - t443 * t474 + t454 * t499;
	t446 = t452 * qJD(3) + t474 * t487;
	t445 = -t453 * qJD(3) - t471 * t487;
	t437 = t450 * qJD(3) + t442 * t474 + t471 * t489;
	t436 = t451 * qJD(3) + t442 * t471 - t474 * t489;
	t435 = t437 * t467 - t441 * t466 + (-t451 * t466 + t456 * t467) * qJD(5);
	t434 = -t437 * t466 - t441 * t467 + (-t451 * t467 - t456 * t466) * qJD(5);
	t1 = [t519, t483 * t466 + t478 * t467, -t436 * t467 - t450 * t496, 0, t434, 0; t435, t482 * t466 + t477 * t467, t438 * t467 + t480 * t496, 0, t518, 0; 0, (t466 * t481 - t515 * t467) * t469, t445 * t467 - t452 * t496, 0, t467 * t488 - t446 * t466 + (-t453 * t467 + t466 * t508) * qJD(5), 0; -t518, -t478 * t466 + t483 * t467, t436 * t466 - t450 * t495, 0, -t435, 0; t434, -t477 * t466 + t482 * t467, -t438 * t466 + t480 * t495, 0, t519, 0; 0, (t515 * t466 + t467 * t481) * t469, -t445 * t466 - t452 * t495, 0, -t466 * t488 - t446 * t467 + (t453 * t466 + t467 * t508) * qJD(5), 0; t438, t441 * t471 - t456 * t498, t437, 0, 0, 0; t436, -t443 * t471 - t454 * t498, -t440, 0, 0, 0; 0, (-t471 * t501 + t474 * t497) * t469, t446, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:12:22
	% EndTime: 2019-10-10 12:12:23
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (627->69), mult. (1126->128), div. (0->0), fcn. (1210->10), ass. (0->73)
	t514 = pkin(12) + qJ(5) + qJ(6);
	t512 = sin(t514);
	t513 = cos(t514);
	t520 = sin(qJ(1));
	t517 = cos(pkin(6));
	t535 = qJD(2) * t517 + qJD(1);
	t519 = sin(qJ(2));
	t558 = t520 * t519;
	t548 = t517 * t558;
	t553 = qJD(2) * t519;
	t522 = cos(qJ(2));
	t523 = cos(qJ(1));
	t555 = t523 * t522;
	t491 = -qJD(1) * t548 - t520 * t553 + t535 * t555;
	t518 = sin(qJ(3));
	t521 = cos(qJ(3));
	t556 = t523 * t519;
	t557 = t520 * t522;
	t501 = t517 * t556 + t557;
	t516 = sin(pkin(6));
	t559 = t516 * t523;
	t530 = t501 * t518 + t521 * t559;
	t554 = qJD(1) * t516;
	t545 = t520 * t554;
	t485 = -t530 * qJD(3) + t491 * t521 + t518 * t545;
	t546 = t517 * t555;
	t500 = -t546 + t558;
	t515 = qJD(5) + qJD(6);
	t540 = -t500 * t515 - t485;
	t502 = t517 * t557 + t556;
	t490 = t502 * qJD(1) + t501 * qJD(2);
	t547 = t518 * t559;
	t567 = -(-t501 * t521 + t547) * t515 - t490;
	t478 = t540 * t512 - t513 * t567;
	t479 = t512 * t567 + t540 * t513;
	t552 = qJD(2) * t522;
	t488 = -qJD(1) * t546 - t523 * t552 + t535 * t558;
	t503 = -t548 + t555;
	t561 = t516 * t518;
	t529 = t503 * t521 + t520 * t561;
	t568 = -t529 * t515 - t488;
	t549 = qJD(3) * t522;
	t566 = (qJD(2) * t521 - t515) * t519 + t518 * t549;
	t564 = t512 * t515;
	t563 = t513 * t515;
	t562 = t515 * t521;
	t560 = t516 * t521;
	t551 = qJD(3) * t518;
	t550 = qJD(3) * t521;
	t544 = t523 * t554;
	t543 = t516 * t552;
	t489 = -t501 * qJD(1) - t502 * qJD(2);
	t497 = -t503 * t518 + t520 * t560;
	t483 = t497 * qJD(3) + t489 * t521 + t518 * t544;
	t541 = t502 * t515 + t483;
	t534 = t502 * t562 + t489;
	t533 = t500 * t562 + t491;
	t498 = t517 * t521 - t519 * t561;
	t493 = t498 * qJD(3) + t521 * t543;
	t532 = t515 * t516 * t522 - t493;
	t531 = (qJD(2) - t562) * t522;
	t499 = t517 * t518 + t519 * t560;
	t528 = -t499 * t515 + t516 * t553;
	t484 = qJD(3) * t547 - t491 * t518 - t501 * t550 + t521 * t545;
	t525 = t488 * t521 + t502 * t551 + t503 * t515;
	t524 = -t490 * t521 + t500 * t551 + t501 * t515;
	t492 = -t499 * qJD(3) - t518 * t543;
	t482 = t529 * qJD(3) + t489 * t518 - t521 * t544;
	t481 = -t528 * t512 + t532 * t513;
	t480 = t532 * t512 + t528 * t513;
	t477 = t568 * t512 + t541 * t513;
	t476 = -t541 * t512 + t568 * t513;
	t1 = [t479, t534 * t512 + t525 * t513, -t482 * t513 - t497 * t564, 0, t476, t476; t477, t533 * t512 + t524 * t513, t484 * t513 + t530 * t564, 0, t478, t478; 0, (t512 * t531 - t566 * t513) * t516, t492 * t513 - t498 * t564, 0, t480, t480; -t478, -t525 * t512 + t534 * t513, t482 * t512 - t497 * t563, 0, -t477, -t477; t476, -t524 * t512 + t533 * t513, -t484 * t512 + t530 * t563, 0, t479, t479; 0, (t566 * t512 + t513 * t531) * t516, -t492 * t512 - t498 * t563, 0, t481, t481; t484, t488 * t518 - t502 * t550, t483, 0, 0, 0; t482, -t490 * t518 - t500 * t550, t485, 0, 0, 0; 0, (-t518 * t553 + t521 * t549) * t516, t493, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end