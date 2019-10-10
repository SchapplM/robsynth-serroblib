% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR7
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
% Datum: 2019-10-10 13:27
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRRR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:27:07
	% EndTime: 2019-10-10 13:27:07
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:27:07
	% EndTime: 2019-10-10 13:27:07
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
	% StartTime: 2019-10-10 13:27:07
	% EndTime: 2019-10-10 13:27:07
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
	% StartTime: 2019-10-10 13:27:08
	% EndTime: 2019-10-10 13:27:08
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
	% StartTime: 2019-10-10 13:27:10
	% EndTime: 2019-10-10 13:27:11
	% DurationCPUTime: 0.52s
	% Computational Cost: add. (289->73), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->66)
	t464 = sin(qJ(1));
	t460 = cos(pkin(6));
	t477 = qJD(2) * t460 + qJD(1);
	t463 = sin(qJ(2));
	t498 = t464 * t463;
	t485 = t460 * t498;
	t493 = qJD(2) * t463;
	t467 = cos(qJ(2));
	t468 = cos(qJ(1));
	t495 = t468 * t467;
	t437 = -qJD(1) * t485 - t464 * t493 + t477 * t495;
	t462 = sin(qJ(3));
	t466 = cos(qJ(3));
	t496 = t468 * t463;
	t497 = t464 * t467;
	t448 = t460 * t496 + t497;
	t459 = sin(pkin(6));
	t499 = t459 * t468;
	t472 = t448 * t462 + t466 * t499;
	t494 = qJD(1) * t459;
	t482 = t464 * t494;
	t433 = t472 * qJD(3) - t437 * t466 - t462 * t482;
	t449 = t460 * t497 + t496;
	t436 = t449 * qJD(1) + t448 * qJD(2);
	t484 = t462 * t499;
	t442 = -t448 * t466 + t484;
	t483 = t460 * t495;
	t447 = -t483 + t498;
	t461 = sin(qJ(4));
	t465 = cos(qJ(4));
	t511 = t433 * t465 - t436 * t461 + (-t442 * t461 - t447 * t465) * qJD(4);
	t510 = (t442 * t465 - t447 * t461) * qJD(4) + t433 * t461 + t436 * t465;
	t489 = qJD(3) * t467;
	t507 = (qJD(2) * t466 - qJD(4)) * t463 + t462 * t489;
	t502 = t459 * t462;
	t501 = t459 * t466;
	t500 = t459 * t467;
	t492 = qJD(2) * t467;
	t491 = qJD(3) * t462;
	t490 = qJD(3) * t466;
	t488 = qJD(4) * t461;
	t487 = qJD(4) * t465;
	t486 = qJD(4) * t466;
	t481 = t468 * t494;
	t480 = t459 * t493;
	t479 = t459 * t492;
	t435 = -t448 * qJD(1) - t449 * qJD(2);
	t475 = t449 * t486 + t435;
	t474 = t447 * t486 + t437;
	t473 = (qJD(2) - t486) * t467;
	t450 = -t485 + t495;
	t443 = -t450 * t462 + t464 * t501;
	t444 = t450 * t466 + t464 * t502;
	t446 = t460 * t462 + t463 * t501;
	t445 = t460 * t466 - t463 * t502;
	t431 = qJD(3) * t484 - t437 * t462 - t448 * t490 + t466 * t482;
	t434 = -qJD(1) * t483 - t468 * t492 + t477 * t498;
	t470 = qJD(4) * t450 + t434 * t466 + t449 * t491;
	t469 = qJD(4) * t448 - t436 * t466 + t447 * t491;
	t439 = t445 * qJD(3) + t466 * t479;
	t438 = -t446 * qJD(3) - t462 * t479;
	t430 = t443 * qJD(3) + t435 * t466 + t462 * t481;
	t429 = t444 * qJD(3) + t435 * t462 - t466 * t481;
	t428 = t430 * t465 - t434 * t461 + (-t444 * t461 + t449 * t465) * qJD(4);
	t427 = -t430 * t461 - t434 * t465 + (-t444 * t465 - t449 * t461) * qJD(4);
	t1 = [t511, t475 * t461 + t470 * t465, -t429 * t465 - t443 * t488, t427, 0, 0; t428, t474 * t461 + t469 * t465, t431 * t465 + t472 * t488, t510, 0, 0; 0, (t461 * t473 - t507 * t465) * t459, t438 * t465 - t445 * t488, t465 * t480 - t439 * t461 + (-t446 * t465 + t461 * t500) * qJD(4), 0, 0; -t510, -t470 * t461 + t475 * t465, t429 * t461 - t443 * t487, -t428, 0, 0; t427, -t469 * t461 + t474 * t465, -t431 * t461 + t472 * t487, t511, 0, 0; 0, (t507 * t461 + t465 * t473) * t459, -t438 * t461 - t445 * t487, -t461 * t480 - t439 * t465 + (t446 * t461 + t465 * t500) * qJD(4), 0, 0; t431, t434 * t462 - t449 * t490, t430, 0, 0, 0; t429, -t436 * t462 - t447 * t490, -t433, 0, 0, 0; 0, (-t462 * t493 + t466 * t489) * t459, t439, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:27:11
	% EndTime: 2019-10-10 13:27:11
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (521->69), mult. (1126->128), div. (0->0), fcn. (1210->10), ass. (0->73)
	t514 = qJ(4) + qJ(5);
	t511 = sin(t514);
	t512 = cos(t514);
	t519 = sin(qJ(1));
	t516 = cos(pkin(6));
	t534 = qJD(2) * t516 + qJD(1);
	t518 = sin(qJ(2));
	t557 = t519 * t518;
	t547 = t516 * t557;
	t552 = qJD(2) * t518;
	t521 = cos(qJ(2));
	t522 = cos(qJ(1));
	t554 = t522 * t521;
	t490 = -qJD(1) * t547 - t519 * t552 + t534 * t554;
	t517 = sin(qJ(3));
	t520 = cos(qJ(3));
	t555 = t522 * t518;
	t556 = t519 * t521;
	t500 = t516 * t555 + t556;
	t515 = sin(pkin(6));
	t558 = t515 * t522;
	t529 = t500 * t517 + t520 * t558;
	t553 = qJD(1) * t515;
	t544 = t519 * t553;
	t484 = -t529 * qJD(3) + t490 * t520 + t517 * t544;
	t545 = t516 * t554;
	t499 = -t545 + t557;
	t513 = qJD(4) + qJD(5);
	t539 = -t499 * t513 - t484;
	t501 = t516 * t556 + t555;
	t489 = t501 * qJD(1) + t500 * qJD(2);
	t546 = t517 * t558;
	t566 = -(-t500 * t520 + t546) * t513 - t489;
	t477 = t539 * t511 - t512 * t566;
	t478 = t511 * t566 + t539 * t512;
	t551 = qJD(2) * t521;
	t487 = -qJD(1) * t545 - t522 * t551 + t534 * t557;
	t502 = -t547 + t554;
	t560 = t515 * t517;
	t528 = t502 * t520 + t519 * t560;
	t567 = -t528 * t513 - t487;
	t548 = qJD(3) * t521;
	t565 = (qJD(2) * t520 - t513) * t518 + t517 * t548;
	t563 = t511 * t513;
	t562 = t512 * t513;
	t561 = t513 * t520;
	t559 = t515 * t520;
	t550 = qJD(3) * t517;
	t549 = qJD(3) * t520;
	t543 = t522 * t553;
	t542 = t515 * t551;
	t488 = -t500 * qJD(1) - t501 * qJD(2);
	t496 = -t502 * t517 + t519 * t559;
	t482 = t496 * qJD(3) + t488 * t520 + t517 * t543;
	t540 = t501 * t513 + t482;
	t533 = t501 * t561 + t488;
	t532 = t499 * t561 + t490;
	t497 = t516 * t520 - t518 * t560;
	t492 = t497 * qJD(3) + t520 * t542;
	t531 = t513 * t515 * t521 - t492;
	t530 = (qJD(2) - t561) * t521;
	t498 = t516 * t517 + t518 * t559;
	t527 = -t498 * t513 + t515 * t552;
	t483 = qJD(3) * t546 - t490 * t517 - t500 * t549 + t520 * t544;
	t524 = t487 * t520 + t501 * t550 + t502 * t513;
	t523 = -t489 * t520 + t499 * t550 + t500 * t513;
	t491 = -t498 * qJD(3) - t517 * t542;
	t481 = t528 * qJD(3) + t488 * t517 - t520 * t543;
	t480 = -t527 * t511 + t531 * t512;
	t479 = t531 * t511 + t527 * t512;
	t476 = t567 * t511 + t540 * t512;
	t475 = -t540 * t511 + t567 * t512;
	t1 = [t478, t533 * t511 + t524 * t512, -t481 * t512 - t496 * t563, t475, t475, 0; t476, t532 * t511 + t523 * t512, t483 * t512 + t529 * t563, t477, t477, 0; 0, (t511 * t530 - t565 * t512) * t515, t491 * t512 - t497 * t563, t479, t479, 0; -t477, -t524 * t511 + t533 * t512, t481 * t511 - t496 * t562, -t476, -t476, 0; t475, -t523 * t511 + t532 * t512, -t483 * t511 + t529 * t562, t478, t478, 0; 0, (t565 * t511 + t512 * t530) * t515, -t491 * t511 - t497 * t562, t480, t480, 0; t483, t487 * t517 - t501 * t549, t482, 0, 0, 0; t481, -t489 * t517 - t499 * t549, t484, 0, 0, 0; 0, (-t517 * t552 + t520 * t548) * t515, t492, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:27:11
	% EndTime: 2019-10-10 13:27:11
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (825->70), mult. (1354->128), div. (0->0), fcn. (1456->10), ass. (0->73)
	t520 = qJ(4) + qJ(5) + qJ(6);
	t517 = sin(t520);
	t518 = cos(t520);
	t525 = sin(qJ(1));
	t522 = cos(pkin(6));
	t540 = qJD(2) * t522 + qJD(1);
	t524 = sin(qJ(2));
	t563 = t525 * t524;
	t553 = t522 * t563;
	t558 = qJD(2) * t524;
	t527 = cos(qJ(2));
	t528 = cos(qJ(1));
	t560 = t528 * t527;
	t496 = -qJD(1) * t553 - t525 * t558 + t540 * t560;
	t523 = sin(qJ(3));
	t526 = cos(qJ(3));
	t561 = t528 * t524;
	t562 = t525 * t527;
	t506 = t522 * t561 + t562;
	t521 = sin(pkin(6));
	t564 = t521 * t528;
	t535 = t506 * t523 + t526 * t564;
	t559 = qJD(1) * t521;
	t550 = t525 * t559;
	t490 = -t535 * qJD(3) + t496 * t526 + t523 * t550;
	t551 = t522 * t560;
	t505 = -t551 + t563;
	t519 = qJD(4) + qJD(5) + qJD(6);
	t545 = -t505 * t519 - t490;
	t507 = t522 * t562 + t561;
	t495 = t507 * qJD(1) + t506 * qJD(2);
	t552 = t523 * t564;
	t572 = -(-t506 * t526 + t552) * t519 - t495;
	t483 = t545 * t517 - t518 * t572;
	t484 = t517 * t572 + t545 * t518;
	t557 = qJD(2) * t527;
	t493 = -qJD(1) * t551 - t528 * t557 + t540 * t563;
	t508 = -t553 + t560;
	t566 = t521 * t523;
	t534 = t508 * t526 + t525 * t566;
	t573 = -t534 * t519 - t493;
	t554 = qJD(3) * t527;
	t571 = (qJD(2) * t526 - t519) * t524 + t523 * t554;
	t569 = t517 * t519;
	t568 = t518 * t519;
	t567 = t519 * t526;
	t565 = t521 * t526;
	t556 = qJD(3) * t523;
	t555 = qJD(3) * t526;
	t549 = t528 * t559;
	t548 = t521 * t557;
	t494 = -t506 * qJD(1) - t507 * qJD(2);
	t502 = -t508 * t523 + t525 * t565;
	t488 = t502 * qJD(3) + t494 * t526 + t523 * t549;
	t546 = t507 * t519 + t488;
	t539 = t507 * t567 + t494;
	t538 = t505 * t567 + t496;
	t503 = t522 * t526 - t524 * t566;
	t498 = t503 * qJD(3) + t526 * t548;
	t537 = t519 * t521 * t527 - t498;
	t536 = (qJD(2) - t567) * t527;
	t504 = t522 * t523 + t524 * t565;
	t533 = -t504 * t519 + t521 * t558;
	t489 = qJD(3) * t552 - t496 * t523 - t506 * t555 + t526 * t550;
	t530 = t493 * t526 + t507 * t556 + t508 * t519;
	t529 = -t495 * t526 + t505 * t556 + t506 * t519;
	t497 = -t504 * qJD(3) - t523 * t548;
	t487 = t534 * qJD(3) + t494 * t523 - t526 * t549;
	t486 = -t533 * t517 + t537 * t518;
	t485 = t537 * t517 + t533 * t518;
	t482 = t573 * t517 + t546 * t518;
	t481 = -t546 * t517 + t573 * t518;
	t1 = [t484, t539 * t517 + t530 * t518, -t487 * t518 - t502 * t569, t481, t481, t481; t482, t538 * t517 + t529 * t518, t489 * t518 + t535 * t569, t483, t483, t483; 0, (t517 * t536 - t571 * t518) * t521, t497 * t518 - t503 * t569, t485, t485, t485; -t483, -t530 * t517 + t539 * t518, t487 * t517 - t502 * t568, -t482, -t482, -t482; t481, -t529 * t517 + t538 * t518, -t489 * t517 + t535 * t568, t484, t484, t484; 0, (t571 * t517 + t518 * t536) * t521, -t497 * t517 - t503 * t568, t486, t486, t486; t489, t493 * t523 - t507 * t555, t488, 0, 0, 0; t487, -t495 * t523 - t505 * t555, t490, 0, 0, 0; 0, (-t523 * t558 + t526 * t554) * t521, t498, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end