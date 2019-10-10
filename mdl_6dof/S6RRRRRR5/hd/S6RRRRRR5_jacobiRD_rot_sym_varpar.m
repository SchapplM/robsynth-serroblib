% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:22
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
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
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:27
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
	% StartTime: 2019-10-10 13:22:27
	% EndTime: 2019-10-10 13:22:28
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
	% StartTime: 2019-10-10 13:22:28
	% EndTime: 2019-10-10 13:22:28
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
	% StartTime: 2019-10-10 13:22:28
	% EndTime: 2019-10-10 13:22:28
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (446->34), mult. (518->59), div. (0->0), fcn. (538->8), ass. (0->43)
	t363 = cos(pkin(6));
	t365 = sin(qJ(1));
	t364 = sin(qJ(2));
	t385 = t365 * t364;
	t376 = t363 * t385;
	t379 = qJD(2) * t364;
	t366 = cos(qJ(2));
	t367 = cos(qJ(1));
	t381 = t367 * t366;
	t350 = -qJD(1) * t376 - t365 * t379 + (qJD(2) * t363 + qJD(1)) * t381;
	t360 = qJD(3) + qJD(4) + qJD(5);
	t362 = sin(pkin(6));
	t386 = t360 * t362;
	t390 = t367 * t386 - t350;
	t371 = t376 - t381;
	t380 = qJD(1) * t362;
	t389 = t371 * t360 + t367 * t380;
	t361 = qJ(3) + qJ(4) + qJ(5);
	t358 = sin(t361);
	t359 = cos(t361);
	t382 = t367 * t364;
	t384 = t365 * t366;
	t354 = t363 * t382 + t384;
	t370 = -t354 * t360 + t365 * t380;
	t343 = t390 * t358 + t370 * t359;
	t388 = t358 * t360;
	t387 = t359 * t360;
	t383 = t366 * t360;
	t378 = t364 * t386;
	t374 = t362 * qJD(2) * t366;
	t355 = -t363 * t384 - t382;
	t348 = -t354 * qJD(1) + t355 * qJD(2);
	t373 = t365 * t386 + t348;
	t353 = t363 * t381 - t385;
	t369 = -t360 * t363 - t374;
	t344 = -t370 * t358 + t390 * t359;
	t349 = t355 * qJD(1) - t354 * qJD(2);
	t347 = -t353 * qJD(1) + t371 * qJD(2);
	t346 = t358 * t378 + t369 * t359;
	t345 = t369 * t358 - t359 * t378;
	t342 = t389 * t358 + t373 * t359;
	t341 = -t373 * t358 + t389 * t359;
	t1 = [t344, t347 * t359 - t355 * t388, t341, t341, t341, 0; t342, t349 * t359 - t353 * t388, t343, t343, t343, 0; 0, (-t358 * t383 - t359 * t379) * t362, t345, t345, t345, 0; -t343, -t347 * t358 - t355 * t387, -t342, -t342, -t342, 0; t341, -t349 * t358 - t353 * t387, t344, t344, t344, 0; 0, (t358 * t379 - t359 * t383) * t362, t346, t346, t346, 0; t349, t348, 0, 0, 0, 0; -t347, t350, 0, 0, 0, 0; 0, t374, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:22:31
	% EndTime: 2019-10-10 13:22:31
	% DurationCPUTime: 0.61s
	% Computational Cost: add. (1041->77), mult. (1306->143), div. (0->0), fcn. (1404->10), ass. (0->77)
	t565 = sin(qJ(1));
	t562 = cos(pkin(6));
	t579 = qJD(2) * t562 + qJD(1);
	t564 = sin(qJ(2));
	t598 = t565 * t564;
	t585 = t562 * t598;
	t593 = qJD(2) * t564;
	t567 = cos(qJ(2));
	t568 = cos(qJ(1));
	t595 = t568 * t567;
	t537 = -qJD(1) * t585 - t565 * t593 + t579 * t595;
	t596 = t568 * t564;
	t597 = t565 * t567;
	t547 = t562 * t596 + t597;
	t560 = qJ(3) + qJ(4) + qJ(5);
	t557 = sin(t560);
	t558 = cos(t560);
	t559 = qJD(3) + qJD(4) + qJD(5);
	t561 = sin(pkin(6));
	t594 = qJD(1) * t561;
	t582 = t565 * t594;
	t599 = t561 * t568;
	t528 = (-t547 * t559 + t582) * t557 - (t559 * t599 - t537) * t558;
	t548 = t562 * t597 + t596;
	t536 = t548 * qJD(1) + t547 * qJD(2);
	t584 = t557 * t599;
	t540 = -t547 * t558 + t584;
	t583 = t562 * t595;
	t546 = -t583 + t598;
	t563 = sin(qJ(6));
	t566 = cos(qJ(6));
	t614 = -t528 * t566 + (-t540 * t563 - t546 * t566) * qJD(6) - t536 * t563;
	t613 = (t540 * t566 - t546 * t563) * qJD(6) - t528 * t563 + t536 * t566;
	t603 = t559 * t567;
	t610 = (qJD(2) * t558 - qJD(6)) * t564 + t557 * t603;
	t605 = t557 * t559;
	t604 = t558 * t559;
	t602 = t561 * t564;
	t601 = t561 * t565;
	t600 = t561 * t567;
	t592 = qJD(2) * t567;
	t591 = qJD(6) * t558;
	t590 = qJD(6) * t563;
	t589 = qJD(6) * t566;
	t587 = t557 * t602;
	t586 = t558 * t602;
	t581 = t568 * t594;
	t580 = t561 * t593;
	t535 = -t547 * qJD(1) - t548 * qJD(2);
	t577 = t559 * t601 + t535;
	t575 = t548 * t591 + t535;
	t574 = t546 * t591 + t537;
	t573 = (qJD(2) - t591) * t567;
	t571 = t559 * t562 + t561 * t592;
	t527 = -t537 * t557 - t547 * t604 + t558 * t582 + t559 * t584;
	t534 = -qJD(1) * t583 - t568 * t592 + t579 * t598;
	t549 = -t585 + t595;
	t570 = qJD(6) * t549 + t534 * t558 + t548 * t605;
	t569 = qJD(6) * t547 - t536 * t558 + t546 * t605;
	t544 = t562 * t557 + t586;
	t543 = t562 * t558 - t587;
	t542 = t549 * t558 + t557 * t601;
	t541 = -t549 * t557 + t558 * t601;
	t538 = -t547 * t557 - t558 * t599;
	t533 = t571 * t558 - t559 * t587;
	t532 = -t571 * t557 - t559 * t586;
	t531 = t532 * t566 - t543 * t590;
	t530 = -t532 * t563 - t543 * t589;
	t526 = t577 * t558 + (-t549 * t559 + t581) * t557;
	t525 = t549 * t604 + t577 * t557 - t558 * t581;
	t524 = t527 * t566 - t538 * t590;
	t523 = -t527 * t563 - t538 * t589;
	t522 = -t525 * t566 - t541 * t590;
	t521 = t525 * t563 - t541 * t589;
	t520 = t526 * t566 - t534 * t563 + (-t542 * t563 + t548 * t566) * qJD(6);
	t519 = -t526 * t563 - t534 * t566 + (-t542 * t566 - t548 * t563) * qJD(6);
	t1 = [t614, t575 * t563 + t570 * t566, t522, t522, t522, t519; t520, t574 * t563 + t569 * t566, t524, t524, t524, t613; 0, (t563 * t573 - t610 * t566) * t561, t531, t531, t531, t566 * t580 - t533 * t563 + (-t544 * t566 + t563 * t600) * qJD(6); -t613, -t570 * t563 + t575 * t566, t521, t521, t521, -t520; t519, -t569 * t563 + t574 * t566, t523, t523, t523, t614; 0, (t610 * t563 + t566 * t573) * t561, t530, t530, t530, -t563 * t580 - t533 * t566 + (t544 * t563 + t566 * t600) * qJD(6); t527, t534 * t557 - t548 * t604, t526, t526, t526, 0; t525, -t536 * t557 - t546 * t604, t528, t528, t528, 0; 0, (-t557 * t593 + t558 * t603) * t561, t533, t533, t533, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end