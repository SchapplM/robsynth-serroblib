% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRPR10_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR10_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:06
	% EndTime: 2019-10-10 12:46:06
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
	% StartTime: 2019-10-10 12:46:07
	% EndTime: 2019-10-10 12:46:07
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
	% StartTime: 2019-10-10 12:46:07
	% EndTime: 2019-10-10 12:46:07
	% DurationCPUTime: 0.22s
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
	% StartTime: 2019-10-10 12:46:08
	% EndTime: 2019-10-10 12:46:08
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (236->34), mult. (414->60), div. (0->0), fcn. (430->8), ass. (0->43)
	t408 = cos(pkin(6));
	t409 = sin(qJ(2));
	t412 = cos(qJ(1));
	t426 = t412 * t409;
	t410 = sin(qJ(1));
	t411 = cos(qJ(2));
	t427 = t410 * t411;
	t397 = t408 * t426 + t427;
	t406 = qJ(3) + qJ(4);
	t403 = sin(t406);
	t404 = cos(t406);
	t405 = qJD(3) + qJD(4);
	t428 = t410 * t409;
	t421 = t408 * t428;
	t423 = qJD(2) * t409;
	t425 = t412 * t411;
	t395 = -qJD(1) * t421 - t410 * t423 + (qJD(2) * t408 + qJD(1)) * t425;
	t407 = sin(pkin(6));
	t430 = t405 * t407;
	t417 = t412 * t430 - t395;
	t424 = qJD(1) * t407;
	t420 = t410 * t424;
	t433 = t417 * t403 + (-t397 * t405 + t420) * t404;
	t432 = t403 * t405;
	t431 = t404 * t405;
	t429 = t405 * t411;
	t422 = t409 * t430;
	t419 = t407 * qJD(2) * t411;
	t398 = -t408 * t427 - t426;
	t393 = -t397 * qJD(1) + t398 * qJD(2);
	t418 = t410 * t430 + t393;
	t396 = t408 * t425 - t428;
	t416 = t421 - t425;
	t414 = t405 * t416 + t412 * t424;
	t413 = t405 * t408 + t419;
	t389 = -t397 * t432 + t403 * t420 - t417 * t404;
	t394 = t398 * qJD(1) - t397 * qJD(2);
	t392 = -t396 * qJD(1) + t416 * qJD(2);
	t391 = -t403 * t422 + t413 * t404;
	t390 = t413 * t403 + t404 * t422;
	t387 = t414 * t403 + t418 * t404;
	t386 = t418 * t403 - t414 * t404;
	t1 = [t394, t393, 0, 0, 0, 0; -t392, t395, 0, 0, 0, 0; 0, t419, 0, 0, 0, 0; t389, -t392 * t404 + t398 * t432, t386, t386, 0, 0; -t387, -t394 * t404 + t396 * t432, -t433, -t433, 0, 0; 0, (t403 * t429 + t404 * t423) * t407, t390, t390, 0, 0; t433, t392 * t403 + t398 * t431, t387, t387, 0, 0; t386, t394 * t403 + t396 * t431, t389, t389, 0, 0; 0, (-t403 * t423 + t404 * t429) * t407, t391, t391, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:09
	% EndTime: 2019-10-10 12:46:10
	% DurationCPUTime: 0.63s
	% Computational Cost: add. (602->84), mult. (1102->143), div. (0->0), fcn. (1184->10), ass. (0->78)
	t554 = sin(qJ(1));
	t551 = cos(pkin(6));
	t568 = qJD(2) * t551 + qJD(1);
	t553 = sin(qJ(2));
	t587 = t554 * t553;
	t574 = t551 * t587;
	t581 = qJD(2) * t553;
	t556 = cos(qJ(2));
	t557 = cos(qJ(1));
	t583 = t557 * t556;
	t526 = -qJD(1) * t574 - t554 * t581 + t568 * t583;
	t584 = t557 * t553;
	t586 = t554 * t556;
	t535 = t551 * t584 + t586;
	t549 = qJ(3) + qJ(4);
	t546 = sin(t549);
	t547 = cos(t549);
	t548 = qJD(3) + qJD(4);
	t550 = sin(pkin(6));
	t582 = qJD(1) * t550;
	t571 = t554 * t582;
	t589 = t550 * t557;
	t573 = t546 * t589;
	t593 = t547 * t548;
	t516 = t526 * t546 + t535 * t593 - t547 * t571 - t548 * t573;
	t536 = t551 * t586 + t584;
	t525 = t536 * qJD(1) + t535 * qJD(2);
	t527 = t535 * t546 + t547 * t589;
	t572 = t551 * t583;
	t534 = -t572 + t587;
	t552 = sin(qJ(6));
	t555 = cos(qJ(6));
	t602 = (t527 * t552 + t534 * t555) * qJD(6) - t516 * t555 + t525 * t552;
	t599 = -t516 * t552 - t525 * t555 + (-t527 * t555 + t534 * t552) * qJD(6);
	t517 = (-t535 * t548 + t571) * t546 - (t548 * t589 - t526) * t547;
	t594 = t546 * t548;
	t592 = t548 * t556;
	t591 = t550 * t553;
	t590 = t550 * t554;
	t588 = t552 * t556;
	t585 = t555 * t556;
	t580 = qJD(2) * t556;
	t579 = qJD(6) * t546;
	t578 = qJD(6) * t552;
	t577 = qJD(6) * t555;
	t576 = t546 * t591;
	t575 = t547 * t591;
	t570 = t557 * t582;
	t569 = t550 * t581;
	t567 = qJD(2) + t579;
	t524 = -t535 * qJD(1) - t536 * qJD(2);
	t566 = t548 * t590 + t524;
	t564 = t536 * t579 - t524;
	t563 = t534 * t579 - t526;
	t561 = t548 * t551 + t550 * t580;
	t523 = -qJD(1) * t572 - t557 * t580 + t568 * t587;
	t537 = -t574 + t583;
	t560 = -qJD(6) * t537 + t523 * t546 - t536 * t593;
	t559 = -qJD(6) * t535 - t525 * t546 - t534 * t593;
	t558 = t547 * t592 + (-qJD(2) * t546 - qJD(6)) * t553;
	t533 = t551 * t546 + t575;
	t532 = -t551 * t547 + t576;
	t531 = t537 * t547 + t546 * t590;
	t530 = t537 * t546 - t547 * t590;
	t528 = t535 * t547 - t573;
	t522 = t561 * t547 - t548 * t576;
	t521 = t561 * t546 + t548 * t575;
	t520 = t522 * t555 - t533 * t578;
	t519 = t522 * t552 + t533 * t577;
	t515 = t566 * t547 + (-t537 * t548 + t570) * t546;
	t514 = t537 * t593 + t566 * t546 - t547 * t570;
	t513 = t517 * t555 - t528 * t578;
	t512 = t517 * t552 + t528 * t577;
	t511 = t515 * t555 - t531 * t578;
	t510 = t515 * t552 + t531 * t577;
	t509 = t514 * t552 - t523 * t555 + (t530 * t555 - t536 * t552) * qJD(6);
	t508 = t514 * t555 + t523 * t552 + (-t530 * t552 - t536 * t555) * qJD(6);
	t1 = [t599, t560 * t552 - t564 * t555, t510, t510, 0, t508; t509, t559 * t552 - t563 * t555, t512, t512, 0, -t602; 0, (t558 * t552 + t567 * t585) * t550, t519, t519, 0, -t552 * t569 + t521 * t555 + (-t532 * t552 + t550 * t585) * qJD(6); t602, t564 * t552 + t560 * t555, t511, t511, 0, -t509; t508, t563 * t552 + t559 * t555, t513, t513, 0, t599; 0, (t558 * t555 - t567 * t588) * t550, t520, t520, 0, -t555 * t569 - t521 * t552 + (-t532 * t555 - t550 * t588) * qJD(6); -t517, t523 * t547 + t536 * t594, -t514, -t514, 0, 0; t515, -t525 * t547 + t534 * t594, -t516, -t516, 0, 0; 0, (-t546 * t592 - t547 * t581) * t550, -t521, -t521, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end