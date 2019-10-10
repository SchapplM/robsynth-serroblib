% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRPR9_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR9_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:10
	% EndTime: 2019-10-10 12:44:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:10
	% EndTime: 2019-10-10 12:44:10
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
	% StartTime: 2019-10-10 12:44:10
	% EndTime: 2019-10-10 12:44:10
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
	% StartTime: 2019-10-10 12:44:11
	% EndTime: 2019-10-10 12:44:11
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
	% StartTime: 2019-10-10 12:44:11
	% EndTime: 2019-10-10 12:44:11
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (236->33), mult. (414->59), div. (0->0), fcn. (430->8), ass. (0->43)
	t345 = cos(pkin(6));
	t347 = sin(qJ(1));
	t346 = sin(qJ(2));
	t367 = t347 * t346;
	t358 = t345 * t367;
	t361 = qJD(2) * t346;
	t348 = cos(qJ(2));
	t349 = cos(qJ(1));
	t363 = t349 * t348;
	t332 = -qJD(1) * t358 - t347 * t361 + (qJD(2) * t345 + qJD(1)) * t363;
	t342 = qJD(3) + qJD(4);
	t344 = sin(pkin(6));
	t368 = t342 * t344;
	t372 = t349 * t368 - t332;
	t353 = t358 - t363;
	t362 = qJD(1) * t344;
	t371 = t353 * t342 + t349 * t362;
	t343 = qJ(3) + qJ(4);
	t340 = sin(t343);
	t341 = cos(t343);
	t364 = t349 * t346;
	t366 = t347 * t348;
	t335 = t345 * t364 + t366;
	t352 = -t335 * t342 + t347 * t362;
	t325 = t372 * t340 + t352 * t341;
	t370 = t340 * t342;
	t369 = t341 * t342;
	t365 = t348 * t342;
	t360 = t346 * t368;
	t356 = t344 * qJD(2) * t348;
	t336 = -t345 * t366 - t364;
	t330 = -t335 * qJD(1) + t336 * qJD(2);
	t355 = t347 * t368 + t330;
	t334 = t345 * t363 - t367;
	t351 = -t342 * t345 - t356;
	t326 = -t352 * t340 + t372 * t341;
	t331 = t336 * qJD(1) - t335 * qJD(2);
	t329 = -t334 * qJD(1) + t353 * qJD(2);
	t328 = t340 * t360 + t351 * t341;
	t327 = t351 * t340 - t341 * t360;
	t324 = t371 * t340 + t355 * t341;
	t323 = -t355 * t340 + t371 * t341;
	t1 = [t326, t329 * t341 - t336 * t370, t323, t323, 0, 0; t324, t331 * t341 - t334 * t370, t325, t325, 0, 0; 0, (-t340 * t365 - t341 * t361) * t344, t327, t327, 0, 0; -t325, -t329 * t340 - t336 * t369, -t324, -t324, 0, 0; t323, -t331 * t340 - t334 * t369, t326, t326, 0, 0; 0, (t340 * t361 - t341 * t365) * t344, t328, t328, 0, 0; t331, t330, 0, 0, 0, 0; -t329, t332, 0, 0, 0, 0; 0, t356, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:13
	% EndTime: 2019-10-10 12:44:13
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (376->49), mult. (680->94), div. (0->0), fcn. (704->10), ass. (0->57)
	t469 = sin(qJ(1));
	t467 = cos(pkin(6));
	t478 = qJD(2) * t467 + qJD(1);
	t468 = sin(qJ(2));
	t492 = t469 * t468;
	t482 = t467 * t492;
	t487 = qJD(2) * t468;
	t470 = cos(qJ(2));
	t471 = cos(qJ(1));
	t489 = t471 * t470;
	t448 = -qJD(1) * t482 - t469 * t487 + t478 * t489;
	t462 = qJD(3) + qJD(4);
	t465 = sin(pkin(6));
	t494 = t462 * t465;
	t501 = t471 * t494 - t448;
	t490 = t471 * t468;
	t491 = t469 * t470;
	t450 = t467 * t490 + t491;
	t463 = qJ(3) + qJ(4);
	t460 = sin(t463);
	t461 = cos(t463);
	t488 = qJD(1) * t465;
	t480 = t469 * t488;
	t440 = (-t450 * t462 + t480) * t460 - t501 * t461;
	t452 = -t482 + t489;
	t451 = -t467 * t491 - t490;
	t446 = -t450 * qJD(1) + t451 * qJD(2);
	t477 = t469 * t494 + t446;
	t479 = t471 * t488;
	t496 = t461 * t462;
	t437 = t452 * t496 + t477 * t460 - t461 * t479;
	t464 = sin(pkin(12));
	t500 = t437 * t464;
	t439 = -t450 * t496 + t501 * t460 + t461 * t480;
	t499 = t439 * t464;
	t486 = qJD(2) * t470;
	t472 = t462 * t467 + t465 * t486;
	t484 = t468 * t494;
	t443 = -t472 * t460 - t461 * t484;
	t498 = t443 * t464;
	t497 = t460 * t462;
	t495 = t461 * t468;
	t493 = t462 * t470;
	t485 = t460 * t493;
	t481 = t467 * t489;
	t445 = -qJD(1) * t481 - t471 * t486 + t478 * t492;
	t475 = -t445 * t461 + t451 * t497;
	t447 = t451 * qJD(1) - t450 * qJD(2);
	t449 = t481 - t492;
	t474 = -t447 * t461 + t449 * t497;
	t466 = cos(pkin(12));
	t444 = -t460 * t484 + t472 * t461;
	t442 = t443 * t466;
	t438 = t477 * t461 + (-t452 * t462 + t479) * t460;
	t436 = t439 * t466;
	t435 = t437 * t466;
	t1 = [-t440 * t466 + t447 * t464, t446 * t464 - t475 * t466, -t435, -t435, 0, 0; t438 * t466 - t445 * t464, t448 * t464 - t474 * t466, t436, t436, 0, 0; 0, (-t466 * t485 + (t464 * t470 - t466 * t495) * qJD(2)) * t465, t442, t442, 0, 0; t440 * t464 + t447 * t466, t446 * t466 + t475 * t464, t500, t500, 0, 0; -t438 * t464 - t445 * t466, t448 * t466 + t474 * t464, -t499, -t499, 0, 0; 0, (t464 * t485 + (t464 * t495 + t466 * t470) * qJD(2)) * t465, -t498, -t498, 0, 0; t439, t445 * t460 + t451 * t496, t438, t438, 0, 0; t437, t447 * t460 + t449 * t496, t440, t440, 0, 0; 0, (-t460 * t487 + t461 * t493) * t465, t444, t444, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:14
	% EndTime: 2019-10-10 12:44:15
	% DurationCPUTime: 0.68s
	% Computational Cost: add. (696->78), mult. (1102->143), div. (0->0), fcn. (1184->10), ass. (0->78)
	t560 = sin(qJ(1));
	t558 = cos(pkin(6));
	t573 = qJD(2) * t558 + qJD(1);
	t559 = sin(qJ(2));
	t592 = t560 * t559;
	t579 = t558 * t592;
	t587 = qJD(2) * t559;
	t561 = cos(qJ(2));
	t562 = cos(qJ(1));
	t589 = t562 * t561;
	t530 = -qJD(1) * t579 - t560 * t587 + t573 * t589;
	t590 = t562 * t559;
	t591 = t560 * t561;
	t539 = t558 * t590 + t591;
	t556 = qJ(3) + qJ(4);
	t552 = sin(t556);
	t553 = cos(t556);
	t555 = qJD(3) + qJD(4);
	t557 = sin(pkin(6));
	t588 = qJD(1) * t557;
	t576 = t560 * t588;
	t593 = t557 * t562;
	t521 = (-t539 * t555 + t576) * t552 - (t555 * t593 - t530) * t553;
	t540 = t558 * t591 + t590;
	t529 = t540 * qJD(1) + t539 * qJD(2);
	t578 = t552 * t593;
	t533 = -t539 * t553 + t578;
	t577 = t558 * t589;
	t538 = -t577 + t592;
	t554 = pkin(12) + qJ(6);
	t550 = sin(t554);
	t551 = cos(t554);
	t608 = -t521 * t551 + (-t533 * t550 - t538 * t551) * qJD(6) - t529 * t550;
	t607 = (t533 * t551 - t538 * t550) * qJD(6) - t521 * t550 + t529 * t551;
	t597 = t555 * t561;
	t604 = (qJD(2) * t553 - qJD(6)) * t559 + t552 * t597;
	t599 = t552 * t555;
	t598 = t553 * t555;
	t596 = t557 * t559;
	t595 = t557 * t560;
	t594 = t557 * t561;
	t586 = qJD(2) * t561;
	t585 = qJD(6) * t550;
	t584 = qJD(6) * t551;
	t583 = qJD(6) * t553;
	t581 = t552 * t596;
	t580 = t553 * t596;
	t575 = t562 * t588;
	t574 = t557 * t587;
	t528 = -t539 * qJD(1) - t540 * qJD(2);
	t571 = t555 * t595 + t528;
	t569 = t540 * t583 + t528;
	t568 = t538 * t583 + t530;
	t567 = (qJD(2) - t583) * t561;
	t565 = t555 * t558 + t557 * t586;
	t520 = -t530 * t552 - t539 * t598 + t553 * t576 + t555 * t578;
	t527 = -qJD(1) * t577 - t562 * t586 + t573 * t592;
	t541 = -t579 + t589;
	t564 = qJD(6) * t541 + t527 * t553 + t540 * t599;
	t563 = qJD(6) * t539 - t529 * t553 + t538 * t599;
	t537 = t558 * t552 + t580;
	t536 = t558 * t553 - t581;
	t535 = t541 * t553 + t552 * t595;
	t534 = -t541 * t552 + t553 * t595;
	t531 = -t539 * t552 - t553 * t593;
	t526 = t565 * t553 - t555 * t581;
	t525 = -t565 * t552 - t555 * t580;
	t524 = t525 * t551 - t536 * t585;
	t523 = -t525 * t550 - t536 * t584;
	t519 = t571 * t553 + (-t541 * t555 + t575) * t552;
	t518 = t541 * t598 + t571 * t552 - t553 * t575;
	t517 = t520 * t551 - t531 * t585;
	t516 = -t520 * t550 - t531 * t584;
	t515 = -t518 * t551 - t534 * t585;
	t514 = t518 * t550 - t534 * t584;
	t513 = t519 * t551 - t527 * t550 + (-t535 * t550 + t540 * t551) * qJD(6);
	t512 = -t519 * t550 - t527 * t551 + (-t535 * t551 - t540 * t550) * qJD(6);
	t1 = [t608, t569 * t550 + t564 * t551, t515, t515, 0, t512; t513, t568 * t550 + t563 * t551, t517, t517, 0, t607; 0, (t550 * t567 - t604 * t551) * t557, t524, t524, 0, t551 * t574 - t526 * t550 + (-t537 * t551 + t550 * t594) * qJD(6); -t607, -t564 * t550 + t569 * t551, t514, t514, 0, -t513; t512, -t563 * t550 + t568 * t551, t516, t516, 0, t608; 0, (t604 * t550 + t551 * t567) * t557, t523, t523, 0, -t550 * t574 - t526 * t551 + (t537 * t550 + t551 * t594) * qJD(6); t520, t527 * t552 - t540 * t598, t519, t519, 0, 0; t518, -t529 * t552 - t538 * t598, t521, t521, 0, 0; 0, (-t552 * t587 + t553 * t597) * t557, t526, t526, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end