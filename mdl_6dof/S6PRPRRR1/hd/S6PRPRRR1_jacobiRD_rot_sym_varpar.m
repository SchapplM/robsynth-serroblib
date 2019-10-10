% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:53
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(11));
	t58 = sin(pkin(11));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (14->7), mult. (60->18), div. (0->0), fcn. (60->8), ass. (0->13)
	t105 = sin(pkin(12));
	t108 = cos(pkin(12));
	t111 = sin(qJ(2));
	t112 = cos(qJ(2));
	t115 = (t105 * t112 + t108 * t111) * qJD(2);
	t103 = (t105 * t111 - t108 * t112) * qJD(2);
	t110 = cos(pkin(6));
	t109 = cos(pkin(11));
	t107 = sin(pkin(6));
	t106 = sin(pkin(11));
	t102 = t110 * t115;
	t101 = t110 * t103;
	t1 = [0, t106 * t102 + t109 * t103, 0, 0, 0, 0; 0, -t109 * t102 + t106 * t103, 0, 0, 0, 0; 0, -t107 * t115, 0, 0, 0, 0; 0, -t106 * t101 + t109 * t115, 0, 0, 0, 0; 0, t109 * t101 + t106 * t115, 0, 0, 0, 0; 0, t107 * t103, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:49
	% EndTime: 2019-10-09 21:53:49
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (81->31), mult. (282->76), div. (0->0), fcn. (310->10), ass. (0->38)
	t294 = sin(pkin(6));
	t298 = sin(qJ(4));
	t310 = t294 * t298;
	t300 = cos(qJ(4));
	t309 = t294 * t300;
	t299 = sin(qJ(2));
	t308 = qJD(2) * t299;
	t301 = cos(qJ(2));
	t307 = qJD(2) * t301;
	t306 = qJD(4) * t298;
	t305 = qJD(4) * t300;
	t292 = sin(pkin(12));
	t295 = cos(pkin(12));
	t297 = cos(pkin(6));
	t279 = (t292 * t308 - t295 * t307) * t297;
	t286 = -t292 * t307 - t295 * t308;
	t293 = sin(pkin(11));
	t296 = cos(pkin(11));
	t270 = -t296 * t279 + t293 * t286;
	t272 = t293 * t279 + t296 * t286;
	t304 = t301 * t292 + t299 * t295;
	t303 = t299 * t292 - t301 * t295;
	t281 = t303 * t294;
	t302 = qJD(2) * t304;
	t285 = t303 * qJD(2);
	t284 = t304 * t297;
	t283 = t303 * t297;
	t282 = t304 * t294;
	t280 = t297 * t302;
	t278 = qJD(2) * t281;
	t277 = t294 * t302;
	t276 = -t293 * t284 - t296 * t303;
	t275 = t293 * t283 - t296 * t304;
	t274 = t296 * t284 - t293 * t303;
	t273 = -t296 * t283 - t293 * t304;
	t271 = t293 * t280 + t296 * t285;
	t269 = -t296 * t280 + t293 * t285;
	t1 = [0, t271 * t300 - t275 * t306, 0, -t272 * t298 + (-t276 * t300 - t293 * t310) * qJD(4), 0, 0; 0, t269 * t300 - t273 * t306, 0, -t270 * t298 + (-t274 * t300 + t296 * t310) * qJD(4), 0, 0; 0, -t277 * t300 + t281 * t306, 0, t278 * t298 + (-t282 * t300 - t297 * t298) * qJD(4), 0, 0; 0, -t271 * t298 - t275 * t305, 0, -t272 * t300 + (t276 * t298 - t293 * t309) * qJD(4), 0, 0; 0, -t269 * t298 - t273 * t305, 0, -t270 * t300 + (t274 * t298 + t296 * t309) * qJD(4), 0, 0; 0, t277 * t298 + t281 * t305, 0, t278 * t300 + (t282 * t298 - t297 * t300) * qJD(4), 0, 0; 0, t272, 0, 0, 0, 0; 0, t270, 0, 0, 0, 0; 0, -t278, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:49
	% EndTime: 2019-10-09 21:53:49
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (199->30), mult. (416->66), div. (0->0), fcn. (458->10), ass. (0->48)
	t350 = qJ(4) + qJ(5);
	t347 = sin(t350);
	t349 = qJD(4) + qJD(5);
	t369 = t347 * t349;
	t348 = cos(t350);
	t368 = t348 * t349;
	t353 = sin(pkin(6));
	t367 = t349 * t353;
	t357 = sin(qJ(2));
	t366 = qJD(2) * t357;
	t358 = cos(qJ(2));
	t365 = qJD(2) * t358;
	t351 = sin(pkin(12));
	t354 = cos(pkin(12));
	t360 = t357 * t351 - t358 * t354;
	t336 = t360 * t353;
	t333 = qJD(2) * t336;
	t356 = cos(pkin(6));
	t364 = -t349 * t356 + t333;
	t334 = (t351 * t366 - t354 * t365) * t356;
	t341 = -t351 * t365 - t354 * t366;
	t352 = sin(pkin(11));
	t355 = cos(pkin(11));
	t325 = -t355 * t334 + t352 * t341;
	t363 = t355 * t367 - t325;
	t327 = t352 * t334 + t355 * t341;
	t362 = -t352 * t367 - t327;
	t361 = t358 * t351 + t357 * t354;
	t359 = qJD(2) * t361;
	t340 = t360 * qJD(2);
	t339 = t361 * t356;
	t338 = t360 * t356;
	t337 = t361 * t353;
	t335 = t356 * t359;
	t332 = t353 * t359;
	t331 = -t352 * t339 - t355 * t360;
	t330 = t352 * t338 - t355 * t361;
	t329 = t355 * t339 - t352 * t360;
	t328 = -t355 * t338 - t352 * t361;
	t326 = t352 * t335 + t355 * t340;
	t324 = -t355 * t335 + t352 * t340;
	t323 = t337 * t369 + t364 * t348;
	t322 = -t337 * t368 + t364 * t347;
	t321 = t331 * t369 + t362 * t348;
	t320 = -t331 * t368 + t362 * t347;
	t319 = t329 * t369 + t363 * t348;
	t318 = -t329 * t368 + t363 * t347;
	t1 = [0, t326 * t348 - t330 * t369, 0, t320, t320, 0; 0, t324 * t348 - t328 * t369, 0, t318, t318, 0; 0, -t332 * t348 + t336 * t369, 0, t322, t322, 0; 0, -t326 * t347 - t330 * t368, 0, t321, t321, 0; 0, -t324 * t347 - t328 * t368, 0, t319, t319, 0; 0, t332 * t347 + t336 * t368, 0, t323, t323, 0; 0, t327, 0, 0, 0, 0; 0, t325, 0, 0, 0, 0; 0, -t333, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:52
	% EndTime: 2019-10-09 21:53:52
	% DurationCPUTime: 0.48s
	% Computational Cost: add. (578->71), mult. (1241->141), div. (0->0), fcn. (1420->12), ass. (0->71)
	t552 = cos(pkin(6));
	t547 = sin(pkin(12));
	t550 = cos(pkin(12));
	t554 = sin(qJ(2));
	t556 = cos(qJ(2));
	t562 = t556 * t547 + t554 * t550;
	t533 = t562 * t552;
	t536 = t554 * t547 - t556 * t550;
	t534 = t536 * qJD(2);
	t546 = qJ(4) + qJ(5);
	t543 = sin(t546);
	t545 = qJD(4) + qJD(5);
	t583 = t543 * t545;
	t544 = cos(t546);
	t582 = t544 * t545;
	t548 = sin(pkin(11));
	t549 = sin(pkin(6));
	t581 = t548 * t549;
	t551 = cos(pkin(11));
	t580 = t549 * t551;
	t577 = qJD(2) * t554;
	t576 = qJD(2) * t556;
	t575 = qJD(6) * t544;
	t553 = sin(qJ(6));
	t574 = qJD(6) * t553;
	t555 = cos(qJ(6));
	t573 = qJD(6) * t555;
	t529 = t549 * t534;
	t572 = t545 * t552 - t529;
	t530 = (t547 * t577 - t550 * t576) * t552;
	t535 = -t547 * t576 - t550 * t577;
	t566 = -t551 * t530 + t548 * t535;
	t571 = t545 * t580 - t566;
	t565 = t548 * t530 + t551 * t535;
	t570 = t545 * t581 + t565;
	t561 = t536 * t552;
	t519 = -t548 * t562 - t551 * t561;
	t569 = -t519 * t575 + t566;
	t522 = t548 * t561 - t551 * t562;
	t568 = -t522 * t575 + t565;
	t531 = t536 * t549;
	t567 = t531 * t575 - t529;
	t564 = t551 * t533 - t548 * t536;
	t563 = -t548 * t533 - t551 * t536;
	t532 = t562 * t549;
	t560 = qJD(2) * t533;
	t512 = t548 * t534 - t551 * t560;
	t559 = -qJD(6) * t564 - t512 * t544 + t519 * t583;
	t515 = t551 * t534 + t548 * t560;
	t558 = -qJD(6) * t563 - t515 * t544 + t522 * t583;
	t528 = qJD(2) * t532;
	t557 = qJD(6) * t532 - t528 * t544 + t531 * t583;
	t525 = t532 * t544 + t552 * t543;
	t524 = -t532 * t543 + t552 * t544;
	t511 = t543 * t581 + t544 * t563;
	t510 = -t543 * t563 + t544 * t581;
	t509 = -t543 * t580 + t544 * t564;
	t508 = -t543 * t564 - t544 * t580;
	t507 = -t532 * t583 + t572 * t544;
	t506 = -t532 * t582 - t572 * t543;
	t505 = t506 * t555 - t524 * t574;
	t504 = -t506 * t553 - t524 * t573;
	t503 = t570 * t544 - t563 * t583;
	t502 = -t570 * t543 - t563 * t582;
	t501 = -t571 * t544 - t564 * t583;
	t500 = t571 * t543 - t564 * t582;
	t499 = t502 * t555 - t510 * t574;
	t498 = -t502 * t553 - t510 * t573;
	t497 = t500 * t555 - t508 * t574;
	t496 = -t500 * t553 - t508 * t573;
	t1 = [0, t568 * t553 - t558 * t555, 0, t499, t499, -t503 * t553 - t515 * t555 + (-t511 * t555 + t522 * t553) * qJD(6); 0, t569 * t553 - t559 * t555, 0, t497, t497, -t501 * t553 - t512 * t555 + (-t509 * t555 + t519 * t553) * qJD(6); 0, t567 * t553 + t557 * t555, 0, t505, t505, -t507 * t553 + t528 * t555 + (-t525 * t555 - t531 * t553) * qJD(6); 0, t558 * t553 + t568 * t555, 0, t498, t498, -t503 * t555 + t515 * t553 + (t511 * t553 + t522 * t555) * qJD(6); 0, t559 * t553 + t569 * t555, 0, t496, t496, -t501 * t555 + t512 * t553 + (t509 * t553 + t519 * t555) * qJD(6); 0, -t557 * t553 + t567 * t555, 0, t504, t504, -t507 * t555 - t528 * t553 + (t525 * t553 - t531 * t555) * qJD(6); 0, t515 * t543 + t522 * t582, 0, t503, t503, 0; 0, t512 * t543 + t519 * t582, 0, t501, t501, 0; 0, -t528 * t543 - t531 * t582, 0, t507, t507, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end