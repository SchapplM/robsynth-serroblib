% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP8
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
% Datum: 2019-10-10 12:31
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRPP8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:15
	% EndTime: 2019-10-10 12:31:16
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:16
	% EndTime: 2019-10-10 12:31:16
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
	% StartTime: 2019-10-10 12:31:16
	% EndTime: 2019-10-10 12:31:16
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
	% StartTime: 2019-10-10 12:31:17
	% EndTime: 2019-10-10 12:31:17
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
	% StartTime: 2019-10-10 12:31:19
	% EndTime: 2019-10-10 12:31:19
	% DurationCPUTime: 0.55s
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
	% StartTime: 2019-10-10 12:31:22
	% EndTime: 2019-10-10 12:31:22
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (289->73), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->67)
	t564 = sin(qJ(1));
	t560 = cos(pkin(6));
	t577 = qJD(2) * t560 + qJD(1);
	t563 = sin(qJ(2));
	t598 = t564 * t563;
	t584 = t560 * t598;
	t592 = qJD(2) * t563;
	t567 = cos(qJ(2));
	t568 = cos(qJ(1));
	t594 = t568 * t567;
	t537 = -qJD(1) * t584 - t564 * t592 + t577 * t594;
	t562 = sin(qJ(3));
	t566 = cos(qJ(3));
	t595 = t568 * t563;
	t597 = t564 * t567;
	t548 = t560 * t595 + t597;
	t559 = sin(pkin(6));
	t600 = t559 * t568;
	t573 = t548 * t562 + t566 * t600;
	t593 = qJD(1) * t559;
	t581 = t564 * t593;
	t533 = t573 * qJD(3) - t537 * t566 - t562 * t581;
	t549 = t560 * t597 + t595;
	t536 = t549 * qJD(1) + t548 * qJD(2);
	t583 = t562 * t600;
	t542 = -t548 * t566 + t583;
	t582 = t560 * t594;
	t547 = -t582 + t598;
	t561 = sin(qJ(4));
	t565 = cos(qJ(4));
	t610 = t533 * t561 + t536 * t565 + (t542 * t565 - t547 * t561) * qJD(4);
	t609 = (t542 * t561 + t547 * t565) * qJD(4) - t533 * t565 + t536 * t561;
	t602 = t559 * t562;
	t601 = t559 * t566;
	t599 = t561 * t567;
	t596 = t565 * t567;
	t591 = qJD(2) * t567;
	t590 = qJD(3) * t562;
	t589 = qJD(3) * t566;
	t588 = qJD(3) * t567;
	t587 = qJD(4) * t561;
	t586 = qJD(4) * t565;
	t585 = qJD(4) * t566;
	t580 = t568 * t593;
	t579 = t559 * t592;
	t578 = t559 * t591;
	t576 = -qJD(2) + t585;
	t535 = -t548 * qJD(1) - t549 * qJD(2);
	t575 = t549 * t585 + t535;
	t574 = t547 * t585 + t537;
	t550 = -t584 + t594;
	t543 = -t550 * t562 + t564 * t601;
	t544 = t550 * t566 + t564 * t602;
	t546 = t560 * t562 + t563 * t601;
	t545 = t560 * t566 - t563 * t602;
	t531 = qJD(3) * t583 - t537 * t562 - t548 * t589 + t566 * t581;
	t534 = -qJD(1) * t582 - t568 * t591 + t577 * t598;
	t571 = qJD(4) * t550 + t534 * t566 + t549 * t590;
	t570 = qJD(4) * t548 - t536 * t566 + t547 * t590;
	t569 = -t562 * t588 + (-qJD(2) * t566 + qJD(4)) * t563;
	t539 = t545 * qJD(3) + t566 * t578;
	t538 = -t546 * qJD(3) - t562 * t578;
	t530 = t543 * qJD(3) + t535 * t566 + t562 * t580;
	t529 = t544 * qJD(3) + t535 * t562 - t566 * t580;
	t528 = t530 * t565 - t534 * t561 + (-t544 * t561 + t549 * t565) * qJD(4);
	t527 = t530 * t561 + t534 * t565 + (t544 * t565 + t549 * t561) * qJD(4);
	t1 = [-t609, t575 * t561 + t571 * t565, -t529 * t565 - t543 * t587, -t527, 0, 0; t528, t574 * t561 + t570 * t565, t531 * t565 + t573 * t587, t610, 0, 0; 0, (t569 * t565 - t576 * t599) * t559, t538 * t565 - t545 * t587, t565 * t579 - t539 * t561 + (-t546 * t565 + t559 * t599) * qJD(4), 0, 0; t531, t534 * t562 - t549 * t589, t530, 0, 0, 0; t529, -t536 * t562 - t547 * t589, -t533, 0, 0, 0; 0, (-t562 * t592 + t566 * t588) * t559, t539, 0, 0, 0; t610, t571 * t561 - t575 * t565, -t529 * t561 + t543 * t586, t528, 0, 0; t527, t570 * t561 - t574 * t565, t531 * t561 - t573 * t586, t609, 0, 0; 0, (t569 * t561 + t576 * t596) * t559, t538 * t561 + t545 * t586, t561 * t579 + t539 * t565 + (-t546 * t561 - t559 * t596) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:19
	% EndTime: 2019-10-10 12:31:20
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (289->76), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->67)
	t520 = sin(qJ(1));
	t516 = cos(pkin(6));
	t533 = qJD(2) * t516 + qJD(1);
	t519 = sin(qJ(2));
	t555 = t520 * t519;
	t541 = t516 * t555;
	t549 = qJD(2) * t519;
	t523 = cos(qJ(2));
	t524 = cos(qJ(1));
	t551 = t524 * t523;
	t495 = -qJD(1) * t541 - t520 * t549 + t533 * t551;
	t552 = t524 * t519;
	t554 = t520 * t523;
	t506 = t516 * t552 + t554;
	t518 = sin(qJ(3));
	t522 = cos(qJ(3));
	t515 = sin(pkin(6));
	t550 = qJD(1) * t515;
	t538 = t520 * t550;
	t557 = t515 * t524;
	t540 = t522 * t557;
	t490 = (-qJD(3) * t506 + t538) * t518 - qJD(3) * t540 + t495 * t522;
	t507 = t516 * t554 + t552;
	t494 = t507 * qJD(1) + t506 * qJD(2);
	t500 = -t506 * t522 + t518 * t557;
	t539 = t516 * t551;
	t505 = -t539 + t555;
	t517 = sin(qJ(4));
	t521 = cos(qJ(4));
	t567 = -t490 * t517 + t494 * t521 + (t500 * t521 - t505 * t517) * qJD(4);
	t566 = (t500 * t517 + t505 * t521) * qJD(4) + t490 * t521 + t494 * t517;
	t489 = t500 * qJD(3) - t495 * t518 + t522 * t538;
	t559 = t515 * t518;
	t558 = t515 * t522;
	t556 = t517 * t523;
	t553 = t521 * t523;
	t548 = qJD(2) * t523;
	t547 = qJD(3) * t518;
	t546 = qJD(3) * t522;
	t545 = qJD(3) * t523;
	t544 = qJD(4) * t517;
	t543 = qJD(4) * t521;
	t542 = qJD(4) * t522;
	t537 = t524 * t550;
	t536 = t515 * t549;
	t535 = t515 * t548;
	t532 = -qJD(2) + t542;
	t493 = -t506 * qJD(1) - t507 * qJD(2);
	t531 = t507 * t542 + t493;
	t530 = t505 * t542 + t495;
	t508 = -t541 + t551;
	t501 = -t508 * t518 + t520 * t558;
	t502 = t508 * t522 + t520 * t559;
	t504 = t516 * t518 + t519 * t558;
	t503 = t516 * t522 - t519 * t559;
	t492 = -qJD(1) * t539 - t524 * t548 + t533 * t555;
	t527 = qJD(4) * t508 + t492 * t522 + t507 * t547;
	t526 = qJD(4) * t506 - t494 * t522 + t505 * t547;
	t525 = -t518 * t545 + (-qJD(2) * t522 + qJD(4)) * t519;
	t498 = -t506 * t518 - t540;
	t497 = t503 * qJD(3) + t522 * t535;
	t496 = -t504 * qJD(3) - t518 * t535;
	t488 = t501 * qJD(3) + t493 * t522 + t518 * t537;
	t487 = -t502 * qJD(3) - t493 * t518 + t522 * t537;
	t486 = t488 * t521 - t492 * t517 + (-t502 * t517 + t507 * t521) * qJD(4);
	t485 = t488 * t517 + t492 * t521 + (t502 * t521 + t507 * t517) * qJD(4);
	t1 = [-t566, t531 * t517 + t527 * t521, t487 * t521 - t501 * t544, -t485, 0, 0; t486, t530 * t517 + t526 * t521, t489 * t521 - t498 * t544, t567, 0, 0; 0, (t525 * t521 - t532 * t556) * t515, t496 * t521 - t503 * t544, t521 * t536 - t497 * t517 + (-t504 * t521 + t515 * t556) * qJD(4), 0, 0; t567, t527 * t517 - t531 * t521, t487 * t517 + t501 * t543, t486, 0, 0; t485, t526 * t517 - t530 * t521, t489 * t517 + t498 * t543, t566, 0, 0; 0, (t525 * t517 + t532 * t553) * t515, t496 * t517 + t503 * t543, t517 * t536 + t497 * t521 + (-t504 * t517 - t515 * t553) * qJD(4), 0, 0; -t489, -t492 * t518 + t507 * t546, -t488, 0, 0, 0; t487, t494 * t518 + t505 * t546, -t490, 0, 0, 0; 0, (t518 * t549 - t522 * t545) * t515, -t497, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end