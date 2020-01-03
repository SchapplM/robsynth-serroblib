% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRR10_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR10_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:37:49
	% EndTime: 2019-12-31 22:37:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:37:49
	% EndTime: 2019-12-31 22:37:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:37:50
	% EndTime: 2019-12-31 22:37:50
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (27->13), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t158 = sin(qJ(2));
	t159 = sin(qJ(1));
	t171 = t158 * t159;
	t161 = cos(qJ(1));
	t170 = t158 * t161;
	t160 = cos(qJ(2));
	t169 = t159 * t160;
	t168 = t160 * t161;
	t156 = sin(pkin(5));
	t167 = qJD(1) * t156;
	t166 = qJD(2) * t156;
	t157 = cos(pkin(5));
	t165 = -t157 * t168 + t171;
	t164 = t157 * t169 + t170;
	t163 = t157 * t170 + t169;
	t162 = t157 * t171 - t168;
	t155 = t162 * qJD(1) + t165 * qJD(2);
	t154 = t164 * qJD(1) + t163 * qJD(2);
	t153 = t163 * qJD(1) + t164 * qJD(2);
	t152 = t165 * qJD(1) + t162 * qJD(2);
	t1 = [t155, t152, 0, 0, 0; -t153, -t154, 0, 0, 0; 0, -t158 * t166, 0, 0, 0; t154, t153, 0, 0, 0; t152, t155, 0, 0, 0; 0, -t160 * t166, 0, 0, 0; -t159 * t167, 0, 0, 0, 0; t161 * t167, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:37:50
	% EndTime: 2019-12-31 22:37:50
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (94->35), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->37)
	t284 = cos(pkin(5));
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
	t283 = sin(pkin(5));
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
	t1 = [t291, t272 * t288 - t278 * t301, t270, 0, 0; t271, t274 * t288 - t276 * t301, t311, 0, 0; 0, (-t285 * t299 - t288 * t302) * t283, -t285 * t295 + (-t284 * t285 - t286 * t309) * qJD(3), 0, 0; -t311, -t272 * t285 - t278 * t300, -t271, 0, 0; t270, -t274 * t285 - t276 * t300, t291, 0, 0; 0, (t285 * t302 - t288 * t299) * t283, -t288 * t295 + (-t284 * t288 + t286 * t310) * qJD(3), 0, 0; t274, t273, 0, 0, 0; -t272, t275, 0, 0, 0; 0, t295, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:37:50
	% EndTime: 2019-12-31 22:37:51
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (236->33), mult. (414->59), div. (0->0), fcn. (430->8), ass. (0->43)
	t345 = cos(pkin(5));
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
	t344 = sin(pkin(5));
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
	t1 = [t326, t329 * t341 - t336 * t370, t323, t323, 0; t324, t331 * t341 - t334 * t370, t325, t325, 0; 0, (-t340 * t365 - t341 * t361) * t344, t327, t327, 0; -t325, -t329 * t340 - t336 * t369, -t324, -t324, 0; t323, -t331 * t340 - t334 * t369, t326, t326, 0; 0, (t340 * t361 - t341 * t365) * t344, t328, t328, 0; t331, t330, 0, 0, 0; -t329, t332, 0, 0, 0; 0, t356, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:37:52
	% EndTime: 2019-12-31 22:37:52
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (602->77), mult. (1102->143), div. (0->0), fcn. (1184->10), ass. (0->77)
	t547 = sin(qJ(1));
	t544 = cos(pkin(5));
	t561 = qJD(2) * t544 + qJD(1);
	t546 = sin(qJ(2));
	t580 = t547 * t546;
	t567 = t544 * t580;
	t575 = qJD(2) * t546;
	t549 = cos(qJ(2));
	t550 = cos(qJ(1));
	t577 = t550 * t549;
	t519 = -qJD(1) * t567 - t547 * t575 + t561 * t577;
	t578 = t550 * t546;
	t579 = t547 * t549;
	t528 = t544 * t578 + t579;
	t542 = qJ(3) + qJ(4);
	t539 = sin(t542);
	t540 = cos(t542);
	t541 = qJD(3) + qJD(4);
	t543 = sin(pkin(5));
	t576 = qJD(1) * t543;
	t564 = t547 * t576;
	t581 = t543 * t550;
	t510 = (-t528 * t541 + t564) * t539 - (t541 * t581 - t519) * t540;
	t529 = t544 * t579 + t578;
	t518 = t529 * qJD(1) + t528 * qJD(2);
	t566 = t539 * t581;
	t522 = -t528 * t540 + t566;
	t565 = t544 * t577;
	t527 = -t565 + t580;
	t545 = sin(qJ(5));
	t548 = cos(qJ(5));
	t596 = -t510 * t548 + (-t522 * t545 - t527 * t548) * qJD(5) - t518 * t545;
	t595 = (t522 * t548 - t527 * t545) * qJD(5) - t510 * t545 + t518 * t548;
	t585 = t541 * t549;
	t592 = (qJD(2) * t540 - qJD(5)) * t546 + t539 * t585;
	t587 = t539 * t541;
	t586 = t540 * t541;
	t584 = t543 * t546;
	t583 = t543 * t547;
	t582 = t543 * t549;
	t574 = qJD(2) * t549;
	t573 = qJD(5) * t540;
	t572 = qJD(5) * t545;
	t571 = qJD(5) * t548;
	t569 = t539 * t584;
	t568 = t540 * t584;
	t563 = t550 * t576;
	t562 = t543 * t575;
	t517 = -t528 * qJD(1) - t529 * qJD(2);
	t559 = t541 * t583 + t517;
	t557 = t529 * t573 + t517;
	t556 = t527 * t573 + t519;
	t555 = (qJD(2) - t573) * t549;
	t553 = t541 * t544 + t543 * t574;
	t509 = -t519 * t539 - t528 * t586 + t540 * t564 + t541 * t566;
	t516 = -qJD(1) * t565 - t550 * t574 + t561 * t580;
	t530 = -t567 + t577;
	t552 = qJD(5) * t530 + t516 * t540 + t529 * t587;
	t551 = qJD(5) * t528 - t518 * t540 + t527 * t587;
	t526 = t544 * t539 + t568;
	t525 = t544 * t540 - t569;
	t524 = t530 * t540 + t539 * t583;
	t523 = -t530 * t539 + t540 * t583;
	t520 = -t528 * t539 - t540 * t581;
	t515 = t553 * t540 - t541 * t569;
	t514 = -t553 * t539 - t541 * t568;
	t513 = t514 * t548 - t525 * t572;
	t512 = -t514 * t545 - t525 * t571;
	t508 = t559 * t540 + (-t530 * t541 + t563) * t539;
	t507 = t530 * t586 + t559 * t539 - t540 * t563;
	t506 = t509 * t548 - t520 * t572;
	t505 = -t509 * t545 - t520 * t571;
	t504 = -t507 * t548 - t523 * t572;
	t503 = t507 * t545 - t523 * t571;
	t502 = t508 * t548 - t516 * t545 + (-t524 * t545 + t529 * t548) * qJD(5);
	t501 = -t508 * t545 - t516 * t548 + (-t524 * t548 - t529 * t545) * qJD(5);
	t1 = [t596, t557 * t545 + t552 * t548, t504, t504, t501; t502, t556 * t545 + t551 * t548, t506, t506, t595; 0, (t545 * t555 - t592 * t548) * t543, t513, t513, t548 * t562 - t515 * t545 + (-t526 * t548 + t545 * t582) * qJD(5); -t595, -t552 * t545 + t557 * t548, t503, t503, -t502; t501, -t551 * t545 + t556 * t548, t505, t505, t596; 0, (t592 * t545 + t548 * t555) * t543, t512, t512, -t545 * t562 - t515 * t548 + (t526 * t545 + t548 * t582) * qJD(5); t509, t516 * t539 - t529 * t586, t508, t508, 0; t507, -t518 * t539 - t527 * t586, t510, t510, 0; 0, (-t539 * t575 + t540 * t585) * t543, t515, t515, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end