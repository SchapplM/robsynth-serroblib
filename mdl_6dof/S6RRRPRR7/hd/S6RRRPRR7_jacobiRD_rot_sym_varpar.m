% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR7
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
% Datum: 2019-10-10 12:04
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPRR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
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
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
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
	% StartTime: 2019-10-10 12:04:13
	% EndTime: 2019-10-10 12:04:14
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
	% StartTime: 2019-10-10 12:04:14
	% EndTime: 2019-10-10 12:04:14
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (144->36), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->38)
	t310 = cos(pkin(6));
	t312 = sin(qJ(1));
	t311 = sin(qJ(2));
	t331 = t312 * t311;
	t322 = t310 * t331;
	t326 = qJD(2) * t311;
	t313 = cos(qJ(2));
	t314 = cos(qJ(1));
	t328 = t314 * t313;
	t298 = -qJD(1) * t322 - t312 * t326 + (qJD(2) * t310 + qJD(1)) * t328;
	t329 = t314 * t311;
	t330 = t312 * t313;
	t300 = t310 * t329 + t330;
	t308 = qJ(3) + pkin(12);
	t306 = sin(t308);
	t307 = cos(t308);
	t309 = sin(pkin(6));
	t327 = qJD(1) * t309;
	t321 = t312 * t327;
	t332 = t309 * t314;
	t335 = (-t300 * t307 + t306 * t332) * qJD(3) - t298 * t306 + t307 * t321;
	t334 = t309 * t311;
	t333 = t309 * t312;
	t325 = qJD(3) * t306;
	t324 = qJD(3) * t307;
	t323 = qJD(3) * t313;
	t320 = t314 * t327;
	t319 = t309 * qJD(2) * t313;
	t299 = t310 * t328 - t331;
	t301 = -t310 * t330 - t329;
	t317 = t322 - t328;
	t315 = -t298 * t307 + t324 * t332 + (qJD(3) * t300 - t321) * t306;
	t297 = t301 * qJD(1) - t300 * qJD(2);
	t296 = -t300 * qJD(1) + t301 * qJD(2);
	t295 = -t299 * qJD(1) + t317 * qJD(2);
	t294 = t306 * t320 + t296 * t307 + (t306 * t317 + t307 * t333) * qJD(3);
	t293 = t307 * t320 - t296 * t306 + (-t306 * t333 + t307 * t317) * qJD(3);
	t1 = [t315, t295 * t307 - t301 * t325, t293, 0, 0, 0; t294, t297 * t307 - t299 * t325, t335, 0, 0, 0; 0, (-t306 * t323 - t307 * t326) * t309, -t306 * t319 + (-t306 * t310 - t307 * t334) * qJD(3), 0, 0, 0; -t335, -t295 * t306 - t301 * t324, -t294, 0, 0, 0; t293, -t297 * t306 - t299 * t324, t315, 0, 0, 0; 0, (t306 * t326 - t307 * t323) * t309, -t307 * t319 + (t306 * t334 - t307 * t310) * qJD(3), 0, 0, 0; t297, t296, 0, 0, 0, 0; -t295, t298, 0, 0, 0, 0; 0, t319, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:14
	% EndTime: 2019-10-10 12:04:14
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (308->33), mult. (414->59), div. (0->0), fcn. (430->8), ass. (0->43)
	t358 = cos(pkin(6));
	t360 = sin(qJ(1));
	t359 = sin(qJ(2));
	t379 = t360 * t359;
	t371 = t358 * t379;
	t374 = qJD(2) * t359;
	t361 = cos(qJ(2));
	t362 = cos(qJ(1));
	t376 = t362 * t361;
	t345 = -qJD(1) * t371 - t360 * t374 + (qJD(2) * t358 + qJD(1)) * t376;
	t356 = qJD(3) + qJD(5);
	t357 = sin(pkin(6));
	t381 = t356 * t357;
	t385 = t362 * t381 - t345;
	t366 = t371 - t376;
	t375 = qJD(1) * t357;
	t384 = t366 * t356 + t362 * t375;
	t355 = qJ(3) + pkin(12) + qJ(5);
	t353 = sin(t355);
	t354 = cos(t355);
	t377 = t362 * t359;
	t378 = t360 * t361;
	t349 = t358 * t377 + t378;
	t365 = -t349 * t356 + t360 * t375;
	t338 = t385 * t353 + t365 * t354;
	t383 = t353 * t356;
	t382 = t354 * t356;
	t380 = t356 * t361;
	t373 = t359 * t381;
	t369 = t357 * qJD(2) * t361;
	t350 = -t358 * t378 - t377;
	t343 = -t349 * qJD(1) + t350 * qJD(2);
	t368 = t360 * t381 + t343;
	t348 = t358 * t376 - t379;
	t364 = -t356 * t358 - t369;
	t339 = -t365 * t353 + t385 * t354;
	t344 = t350 * qJD(1) - t349 * qJD(2);
	t342 = -t348 * qJD(1) + t366 * qJD(2);
	t341 = t353 * t373 + t364 * t354;
	t340 = t364 * t353 - t354 * t373;
	t337 = t384 * t353 + t368 * t354;
	t336 = -t368 * t353 + t384 * t354;
	t1 = [t339, t342 * t354 - t350 * t383, t336, 0, t336, 0; t337, t344 * t354 - t348 * t383, t338, 0, t338, 0; 0, (-t353 * t380 - t354 * t374) * t357, t340, 0, t340, 0; -t338, -t342 * t353 - t350 * t382, -t337, 0, -t337, 0; t336, -t344 * t353 - t348 * t382, t339, 0, t339, 0; 0, (t353 * t374 - t354 * t380) * t357, t341, 0, t341, 0; t344, t343, 0, 0, 0, 0; -t342, t345, 0, 0, 0, 0; 0, t369, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:16
	% EndTime: 2019-10-10 12:04:17
	% DurationCPUTime: 0.66s
	% Computational Cost: add. (782->77), mult. (1102->143), div. (0->0), fcn. (1184->10), ass. (0->77)
	t558 = sin(qJ(1));
	t555 = cos(pkin(6));
	t572 = qJD(2) * t555 + qJD(1);
	t557 = sin(qJ(2));
	t591 = t558 * t557;
	t578 = t555 * t591;
	t586 = qJD(2) * t557;
	t560 = cos(qJ(2));
	t561 = cos(qJ(1));
	t588 = t561 * t560;
	t530 = -qJD(1) * t578 - t558 * t586 + t572 * t588;
	t589 = t561 * t557;
	t590 = t558 * t560;
	t540 = t555 * t589 + t590;
	t552 = qJ(3) + pkin(12) + qJ(5);
	t550 = sin(t552);
	t551 = cos(t552);
	t553 = qJD(3) + qJD(5);
	t554 = sin(pkin(6));
	t587 = qJD(1) * t554;
	t575 = t558 * t587;
	t592 = t554 * t561;
	t521 = (-t540 * t553 + t575) * t550 - (t553 * t592 - t530) * t551;
	t541 = t555 * t590 + t589;
	t529 = t541 * qJD(1) + t540 * qJD(2);
	t577 = t550 * t592;
	t533 = -t540 * t551 + t577;
	t576 = t555 * t588;
	t539 = -t576 + t591;
	t556 = sin(qJ(6));
	t559 = cos(qJ(6));
	t607 = -t521 * t559 + (-t533 * t556 - t539 * t559) * qJD(6) - t529 * t556;
	t606 = (t533 * t559 - t539 * t556) * qJD(6) - t521 * t556 + t529 * t559;
	t596 = t553 * t560;
	t603 = (qJD(2) * t551 - qJD(6)) * t557 + t550 * t596;
	t598 = t550 * t553;
	t597 = t551 * t553;
	t595 = t554 * t557;
	t594 = t554 * t558;
	t593 = t554 * t560;
	t585 = qJD(2) * t560;
	t584 = qJD(6) * t551;
	t583 = qJD(6) * t556;
	t582 = qJD(6) * t559;
	t580 = t550 * t595;
	t579 = t551 * t595;
	t574 = t561 * t587;
	t573 = t554 * t586;
	t528 = -t540 * qJD(1) - t541 * qJD(2);
	t570 = t553 * t594 + t528;
	t568 = t541 * t584 + t528;
	t567 = t539 * t584 + t530;
	t566 = (qJD(2) - t584) * t560;
	t564 = t553 * t555 + t554 * t585;
	t520 = -t530 * t550 - t540 * t597 + t551 * t575 + t553 * t577;
	t527 = -qJD(1) * t576 - t561 * t585 + t572 * t591;
	t542 = -t578 + t588;
	t563 = qJD(6) * t542 + t527 * t551 + t541 * t598;
	t562 = qJD(6) * t540 - t529 * t551 + t539 * t598;
	t537 = t555 * t550 + t579;
	t536 = t555 * t551 - t580;
	t535 = t542 * t551 + t550 * t594;
	t534 = -t542 * t550 + t551 * t594;
	t531 = -t540 * t550 - t551 * t592;
	t526 = t564 * t551 - t553 * t580;
	t525 = -t564 * t550 - t553 * t579;
	t524 = t525 * t559 - t536 * t583;
	t523 = -t525 * t556 - t536 * t582;
	t519 = t570 * t551 + (-t542 * t553 + t574) * t550;
	t518 = t542 * t597 + t570 * t550 - t551 * t574;
	t517 = t520 * t559 - t531 * t583;
	t516 = -t520 * t556 - t531 * t582;
	t515 = -t518 * t559 - t534 * t583;
	t514 = t518 * t556 - t534 * t582;
	t513 = t519 * t559 - t527 * t556 + (-t535 * t556 + t541 * t559) * qJD(6);
	t512 = -t519 * t556 - t527 * t559 + (-t535 * t559 - t541 * t556) * qJD(6);
	t1 = [t607, t568 * t556 + t563 * t559, t515, 0, t515, t512; t513, t567 * t556 + t562 * t559, t517, 0, t517, t606; 0, (t556 * t566 - t603 * t559) * t554, t524, 0, t524, t559 * t573 - t526 * t556 + (-t537 * t559 + t556 * t593) * qJD(6); -t606, -t563 * t556 + t568 * t559, t514, 0, t514, -t513; t512, -t562 * t556 + t567 * t559, t516, 0, t516, t607; 0, (t603 * t556 + t559 * t566) * t554, t523, 0, t523, -t556 * t573 - t526 * t559 + (t537 * t556 + t559 * t593) * qJD(6); t520, t527 * t550 - t541 * t597, t519, 0, t519, 0; t518, -t529 * t550 - t539 * t597, t521, 0, t521, 0; 0, (-t550 * t586 + t551 * t596) * t554, t526, 0, t526, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end