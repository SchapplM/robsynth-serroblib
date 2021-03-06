% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:33
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRPP9_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP9_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:10
	% EndTime: 2019-10-10 12:33:11
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:11
	% EndTime: 2019-10-10 12:33:11
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
	% StartTime: 2019-10-10 12:33:11
	% EndTime: 2019-10-10 12:33:11
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
	% StartTime: 2019-10-10 12:33:12
	% EndTime: 2019-10-10 12:33:12
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
	% StartTime: 2019-10-10 12:33:14
	% EndTime: 2019-10-10 12:33:14
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
	% StartTime: 2019-10-10 12:33:16
	% EndTime: 2019-10-10 12:33:17
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (289->73), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->66)
	t543 = sin(qJ(1));
	t539 = cos(pkin(6));
	t556 = qJD(2) * t539 + qJD(1);
	t542 = sin(qJ(2));
	t577 = t543 * t542;
	t564 = t539 * t577;
	t572 = qJD(2) * t542;
	t546 = cos(qJ(2));
	t547 = cos(qJ(1));
	t574 = t547 * t546;
	t516 = -qJD(1) * t564 - t543 * t572 + t556 * t574;
	t541 = sin(qJ(3));
	t545 = cos(qJ(3));
	t575 = t547 * t542;
	t576 = t543 * t546;
	t527 = t539 * t575 + t576;
	t538 = sin(pkin(6));
	t578 = t538 * t547;
	t551 = t527 * t541 + t545 * t578;
	t573 = qJD(1) * t538;
	t561 = t543 * t573;
	t512 = t551 * qJD(3) - t516 * t545 - t541 * t561;
	t528 = t539 * t576 + t575;
	t515 = t528 * qJD(1) + t527 * qJD(2);
	t563 = t541 * t578;
	t521 = -t527 * t545 + t563;
	t562 = t539 * t574;
	t526 = -t562 + t577;
	t540 = sin(qJ(4));
	t544 = cos(qJ(4));
	t590 = -t512 * t544 + t515 * t540 + (t521 * t540 + t526 * t544) * qJD(4);
	t589 = t512 * t540 + t515 * t544 + (t521 * t544 - t526 * t540) * qJD(4);
	t568 = qJD(3) * t546;
	t586 = (qJD(2) * t545 - qJD(4)) * t542 + t541 * t568;
	t581 = t538 * t541;
	t580 = t538 * t545;
	t579 = t538 * t546;
	t571 = qJD(2) * t546;
	t570 = qJD(3) * t541;
	t569 = qJD(3) * t545;
	t567 = qJD(4) * t540;
	t566 = qJD(4) * t544;
	t565 = qJD(4) * t545;
	t560 = t547 * t573;
	t559 = t538 * t572;
	t558 = t538 * t571;
	t514 = -t527 * qJD(1) - t528 * qJD(2);
	t554 = -t528 * t565 - t514;
	t553 = -t526 * t565 - t516;
	t552 = (-qJD(2) + t565) * t546;
	t529 = -t564 + t574;
	t522 = -t529 * t541 + t543 * t580;
	t523 = t529 * t545 + t543 * t581;
	t525 = t539 * t541 + t542 * t580;
	t524 = t539 * t545 - t542 * t581;
	t510 = qJD(3) * t563 - t516 * t541 - t527 * t569 + t545 * t561;
	t513 = -qJD(1) * t562 - t547 * t571 + t556 * t577;
	t549 = qJD(4) * t529 + t513 * t545 + t528 * t570;
	t548 = qJD(4) * t527 - t515 * t545 + t526 * t570;
	t518 = t524 * qJD(3) + t545 * t558;
	t517 = -t525 * qJD(3) - t541 * t558;
	t509 = t522 * qJD(3) + t514 * t545 + t541 * t560;
	t508 = t523 * qJD(3) + t514 * t541 - t545 * t560;
	t507 = t509 * t544 - t513 * t540 + (-t523 * t540 + t528 * t544) * qJD(4);
	t506 = t509 * t540 + t513 * t544 + (t523 * t544 + t528 * t540) * qJD(4);
	t1 = [t510, t513 * t541 - t528 * t569, t509, 0, 0, 0; t508, -t515 * t541 - t526 * t569, -t512, 0, 0, 0; 0, (-t541 * t572 + t545 * t568) * t538, t518, 0, 0, 0; t590, t554 * t540 - t549 * t544, t508 * t544 + t522 * t567, t506, 0, 0; -t507, t553 * t540 - t548 * t544, -t510 * t544 - t551 * t567, -t589, 0, 0; 0, (t540 * t552 + t586 * t544) * t538, -t517 * t544 + t524 * t567, -t544 * t559 + t518 * t540 + (t525 * t544 - t540 * t579) * qJD(4), 0, 0; t589, t549 * t540 + t554 * t544, -t508 * t540 + t522 * t566, t507, 0, 0; t506, t548 * t540 + t553 * t544, t510 * t540 - t551 * t566, t590, 0, 0; 0, (-t586 * t540 + t544 * t552) * t538, t517 * t540 + t524 * t566, t540 * t559 + t518 * t544 + (-t525 * t540 - t544 * t579) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:16
	% EndTime: 2019-10-10 12:33:17
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (289->73), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->67)
	t549 = sin(qJ(1));
	t545 = cos(pkin(6));
	t562 = qJD(2) * t545 + qJD(1);
	t548 = sin(qJ(2));
	t583 = t549 * t548;
	t569 = t545 * t583;
	t577 = qJD(2) * t548;
	t552 = cos(qJ(2));
	t553 = cos(qJ(1));
	t579 = t553 * t552;
	t522 = -qJD(1) * t569 - t549 * t577 + t562 * t579;
	t547 = sin(qJ(3));
	t551 = cos(qJ(3));
	t580 = t553 * t548;
	t582 = t549 * t552;
	t533 = t545 * t580 + t582;
	t544 = sin(pkin(6));
	t585 = t544 * t553;
	t558 = t533 * t547 + t551 * t585;
	t578 = qJD(1) * t544;
	t566 = t549 * t578;
	t518 = t558 * qJD(3) - t522 * t551 - t547 * t566;
	t534 = t545 * t582 + t580;
	t521 = t534 * qJD(1) + t533 * qJD(2);
	t568 = t547 * t585;
	t527 = -t533 * t551 + t568;
	t567 = t545 * t579;
	t532 = -t567 + t583;
	t546 = sin(qJ(4));
	t550 = cos(qJ(4));
	t595 = t518 * t546 + t521 * t550 + (t527 * t550 - t532 * t546) * qJD(4);
	t594 = (t527 * t546 + t532 * t550) * qJD(4) - t518 * t550 + t521 * t546;
	t587 = t544 * t547;
	t586 = t544 * t551;
	t584 = t546 * t552;
	t581 = t550 * t552;
	t576 = qJD(2) * t552;
	t575 = qJD(3) * t547;
	t574 = qJD(3) * t551;
	t573 = qJD(3) * t552;
	t572 = qJD(4) * t546;
	t571 = qJD(4) * t550;
	t570 = qJD(4) * t551;
	t565 = t553 * t578;
	t564 = t544 * t577;
	t563 = t544 * t576;
	t561 = -qJD(2) + t570;
	t520 = -t533 * qJD(1) - t534 * qJD(2);
	t560 = t534 * t570 + t520;
	t559 = t532 * t570 + t522;
	t535 = -t569 + t579;
	t528 = -t535 * t547 + t549 * t586;
	t529 = t535 * t551 + t549 * t587;
	t531 = t545 * t547 + t548 * t586;
	t530 = t545 * t551 - t548 * t587;
	t516 = qJD(3) * t568 - t522 * t547 - t533 * t574 + t551 * t566;
	t519 = -qJD(1) * t567 - t553 * t576 + t562 * t583;
	t556 = qJD(4) * t535 + t519 * t551 + t534 * t575;
	t555 = qJD(4) * t533 - t521 * t551 + t532 * t575;
	t554 = -t547 * t573 + (-qJD(2) * t551 + qJD(4)) * t548;
	t524 = t530 * qJD(3) + t551 * t563;
	t523 = -t531 * qJD(3) - t547 * t563;
	t515 = t528 * qJD(3) + t520 * t551 + t547 * t565;
	t514 = t529 * qJD(3) + t520 * t547 - t551 * t565;
	t513 = t515 * t550 - t519 * t546 + (-t529 * t546 + t534 * t550) * qJD(4);
	t512 = t515 * t546 + t519 * t550 + (t529 * t550 + t534 * t546) * qJD(4);
	t1 = [t516, t519 * t547 - t534 * t574, t515, 0, 0, 0; t514, -t521 * t547 - t532 * t574, -t518, 0, 0, 0; 0, (-t547 * t577 + t551 * t573) * t544, t524, 0, 0, 0; t595, t556 * t546 - t560 * t550, -t514 * t546 + t528 * t571, t513, 0, 0; t512, t555 * t546 - t559 * t550, t516 * t546 - t558 * t571, t594, 0, 0; 0, (t554 * t546 + t561 * t581) * t544, t523 * t546 + t530 * t571, t546 * t564 + t524 * t550 + (-t531 * t546 - t544 * t581) * qJD(4), 0, 0; -t594, t560 * t546 + t556 * t550, -t514 * t550 - t528 * t572, -t512, 0, 0; t513, t559 * t546 + t555 * t550, t516 * t550 + t558 * t572, t595, 0, 0; 0, (t554 * t550 - t561 * t584) * t544, t523 * t550 - t530 * t572, t550 * t564 - t524 * t546 + (-t531 * t550 + t544 * t584) * qJD(4), 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end