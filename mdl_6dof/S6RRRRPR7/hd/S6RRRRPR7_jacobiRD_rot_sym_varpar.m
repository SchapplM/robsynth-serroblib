% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR7
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
% Datum: 2019-10-10 12:42
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRPR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:14
	% EndTime: 2019-10-10 12:42:14
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
	% StartTime: 2019-10-10 12:42:15
	% EndTime: 2019-10-10 12:42:15
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
	% StartTime: 2019-10-10 12:42:16
	% EndTime: 2019-10-10 12:42:16
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
	% StartTime: 2019-10-10 12:42:16
	% EndTime: 2019-10-10 12:42:16
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
	% StartTime: 2019-10-10 12:42:16
	% EndTime: 2019-10-10 12:42:16
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (308->33), mult. (414->59), div. (0->0), fcn. (430->8), ass. (0->43)
	t362 = cos(pkin(6));
	t364 = sin(qJ(1));
	t363 = sin(qJ(2));
	t384 = t364 * t363;
	t375 = t362 * t384;
	t378 = qJD(2) * t363;
	t365 = cos(qJ(2));
	t366 = cos(qJ(1));
	t380 = t366 * t365;
	t349 = -qJD(1) * t375 - t364 * t378 + (qJD(2) * t362 + qJD(1)) * t380;
	t360 = qJD(3) + qJD(4);
	t361 = sin(pkin(6));
	t385 = t360 * t361;
	t389 = t366 * t385 - t349;
	t370 = t375 - t380;
	t379 = qJD(1) * t361;
	t388 = t370 * t360 + t366 * t379;
	t359 = qJ(3) + qJ(4) + pkin(12);
	t357 = sin(t359);
	t358 = cos(t359);
	t381 = t366 * t363;
	t383 = t364 * t365;
	t353 = t362 * t381 + t383;
	t369 = -t353 * t360 + t364 * t379;
	t342 = t389 * t357 + t369 * t358;
	t387 = t357 * t360;
	t386 = t358 * t360;
	t382 = t365 * t360;
	t377 = t363 * t385;
	t373 = t361 * qJD(2) * t365;
	t354 = -t362 * t383 - t381;
	t347 = -t353 * qJD(1) + t354 * qJD(2);
	t372 = t364 * t385 + t347;
	t352 = t362 * t380 - t384;
	t368 = -t360 * t362 - t373;
	t343 = -t369 * t357 + t389 * t358;
	t348 = t354 * qJD(1) - t353 * qJD(2);
	t346 = -t352 * qJD(1) + t370 * qJD(2);
	t345 = t357 * t377 + t368 * t358;
	t344 = t368 * t357 - t358 * t377;
	t341 = t388 * t357 + t372 * t358;
	t340 = -t372 * t357 + t388 * t358;
	t1 = [t343, t346 * t358 - t354 * t387, t340, t340, 0, 0; t341, t348 * t358 - t352 * t387, t342, t342, 0, 0; 0, (-t357 * t382 - t358 * t378) * t361, t344, t344, 0, 0; -t342, -t346 * t357 - t354 * t386, -t341, -t341, 0, 0; t340, -t348 * t357 - t352 * t386, t343, t343, 0, 0; 0, (t357 * t378 - t358 * t382) * t361, t345, t345, 0, 0; t348, t347, 0, 0, 0, 0; -t346, t349, 0, 0, 0, 0; 0, t373, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:42:19
	% EndTime: 2019-10-10 12:42:19
	% DurationCPUTime: 0.66s
	% Computational Cost: add. (782->76), mult. (1102->143), div. (0->0), fcn. (1184->10), ass. (0->76)
	t560 = cos(pkin(6));
	t563 = sin(qJ(1));
	t565 = cos(qJ(2));
	t593 = t563 * t565;
	t562 = sin(qJ(2));
	t566 = cos(qJ(1));
	t594 = t562 * t566;
	t545 = t560 * t594 + t593;
	t557 = qJ(3) + qJ(4) + pkin(12);
	t555 = sin(t557);
	t556 = cos(t557);
	t558 = qJD(3) + qJD(4);
	t559 = sin(pkin(6));
	t591 = qJD(1) * t559;
	t580 = t563 * t591;
	t577 = qJD(2) * t560 + qJD(1);
	t595 = t562 * t563;
	t582 = t560 * t595;
	t590 = qJD(2) * t562;
	t592 = t565 * t566;
	t535 = -qJD(1) * t582 - t563 * t590 + t577 * t592;
	t596 = t559 * t566;
	t610 = t558 * t596 - t535;
	t526 = (-t545 * t558 + t580) * t555 - t610 * t556;
	t546 = t560 * t593 + t594;
	t534 = t546 * qJD(1) + t545 * qJD(2);
	t538 = -t545 * t556 + t555 * t596;
	t581 = t560 * t592;
	t544 = -t581 + t595;
	t561 = sin(qJ(6));
	t564 = cos(qJ(6));
	t612 = -t526 * t564 + (-t538 * t561 - t544 * t564) * qJD(6) - t534 * t561;
	t611 = (t538 * t564 - t544 * t561) * qJD(6) - t526 * t561 + t534 * t564;
	t600 = t558 * t565;
	t607 = (qJD(2) * t556 - qJD(6)) * t562 + t555 * t600;
	t602 = t555 * t558;
	t601 = t556 * t558;
	t599 = t559 * t562;
	t598 = t559 * t563;
	t597 = t559 * t565;
	t589 = qJD(2) * t565;
	t588 = qJD(6) * t556;
	t587 = qJD(6) * t561;
	t586 = qJD(6) * t564;
	t584 = t558 * t599;
	t579 = t566 * t591;
	t578 = t559 * t590;
	t533 = -t545 * qJD(1) - t546 * qJD(2);
	t575 = t558 * t598 + t533;
	t573 = t546 * t588 + t533;
	t572 = t544 * t588 + t535;
	t571 = (qJD(2) - t588) * t565;
	t569 = t558 * t560 + t559 * t589;
	t525 = -t545 * t601 + t610 * t555 + t556 * t580;
	t532 = -qJD(1) * t581 - t566 * t589 + t577 * t595;
	t547 = -t582 + t592;
	t568 = qJD(6) * t547 + t532 * t556 + t546 * t602;
	t567 = qJD(6) * t545 - t534 * t556 + t544 * t602;
	t542 = t555 * t560 + t556 * t599;
	t541 = -t555 * t599 + t556 * t560;
	t540 = t547 * t556 + t555 * t598;
	t539 = -t547 * t555 + t556 * t598;
	t536 = -t545 * t555 - t556 * t596;
	t531 = -t555 * t584 + t569 * t556;
	t530 = -t569 * t555 - t556 * t584;
	t529 = t530 * t564 - t541 * t587;
	t528 = -t530 * t561 - t541 * t586;
	t524 = t575 * t556 + (-t547 * t558 + t579) * t555;
	t523 = t547 * t601 + t575 * t555 - t556 * t579;
	t522 = t525 * t564 - t536 * t587;
	t521 = -t525 * t561 - t536 * t586;
	t520 = -t523 * t564 - t539 * t587;
	t519 = t523 * t561 - t539 * t586;
	t518 = t524 * t564 - t532 * t561 + (-t540 * t561 + t546 * t564) * qJD(6);
	t517 = -t524 * t561 - t532 * t564 + (-t540 * t564 - t546 * t561) * qJD(6);
	t1 = [t612, t573 * t561 + t568 * t564, t520, t520, 0, t517; t518, t572 * t561 + t567 * t564, t522, t522, 0, t611; 0, (t561 * t571 - t607 * t564) * t559, t529, t529, 0, t564 * t578 - t531 * t561 + (-t542 * t564 + t561 * t597) * qJD(6); -t611, -t568 * t561 + t573 * t564, t519, t519, 0, -t518; t517, -t567 * t561 + t572 * t564, t521, t521, 0, t612; 0, (t607 * t561 + t564 * t571) * t559, t528, t528, 0, -t561 * t578 - t531 * t564 + (t542 * t561 + t564 * t597) * qJD(6); t525, t532 * t555 - t546 * t601, t524, t524, 0, 0; t523, -t534 * t555 - t544 * t601, t526, t526, 0, 0; 0, (-t555 * t590 + t556 * t600) * t559, t531, t531, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end