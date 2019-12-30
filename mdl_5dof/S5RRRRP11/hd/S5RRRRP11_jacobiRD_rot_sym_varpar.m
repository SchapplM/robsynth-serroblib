% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRP11
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:47
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRP11_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP11_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:47:10
	% EndTime: 2019-12-29 20:47:10
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:47:10
	% EndTime: 2019-12-29 20:47:10
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:47:11
	% EndTime: 2019-12-29 20:47:11
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-12-29 20:47:13
	% EndTime: 2019-12-29 20:47:13
	% DurationCPUTime: 0.26s
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
	% StartTime: 2019-12-29 20:47:20
	% EndTime: 2019-12-29 20:47:21
	% DurationCPUTime: 0.93s
	% Computational Cost: add. (289->73), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->66)
	t464 = sin(qJ(1));
	t460 = cos(pkin(5));
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
	t459 = sin(pkin(5));
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
	t1 = [t511, t475 * t461 + t470 * t465, -t429 * t465 - t443 * t488, t427, 0; t428, t474 * t461 + t469 * t465, t431 * t465 + t472 * t488, t510, 0; 0, (t461 * t473 - t507 * t465) * t459, t438 * t465 - t445 * t488, t465 * t480 - t439 * t461 + (-t446 * t465 + t461 * t500) * qJD(4), 0; -t510, -t470 * t461 + t475 * t465, t429 * t461 - t443 * t487, -t428, 0; t427, -t469 * t461 + t474 * t465, -t431 * t461 + t472 * t487, t511, 0; 0, (t507 * t461 + t465 * t473) * t459, -t438 * t461 - t445 * t487, -t461 * t480 - t439 * t465 + (t446 * t461 + t465 * t500) * qJD(4), 0; t431, t434 * t462 - t449 * t490, t430, 0, 0; t429, -t436 * t462 - t447 * t490, -t433, 0, 0; 0, (-t462 * t493 + t466 * t489) * t459, t439, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:47:19
	% EndTime: 2019-12-29 20:47:20
	% DurationCPUTime: 0.92s
	% Computational Cost: add. (289->73), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->67)
	t564 = sin(qJ(1));
	t560 = cos(pkin(5));
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
	t559 = sin(pkin(5));
	t600 = t559 * t568;
	t573 = t548 * t562 + t566 * t600;
	t593 = qJD(1) * t559;
	t581 = t564 * t593;
	t533 = qJD(3) * t573 - t537 * t566 - t562 * t581;
	t549 = t560 * t597 + t595;
	t536 = qJD(1) * t549 + qJD(2) * t548;
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
	t535 = -qJD(1) * t548 - qJD(2) * t549;
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
	t539 = qJD(3) * t545 + t566 * t578;
	t538 = -qJD(3) * t546 - t562 * t578;
	t530 = qJD(3) * t543 + t535 * t566 + t562 * t580;
	t529 = qJD(3) * t544 + t535 * t562 - t566 * t580;
	t528 = t530 * t565 - t534 * t561 + (-t544 * t561 + t549 * t565) * qJD(4);
	t527 = t530 * t561 + t534 * t565 + (t544 * t565 + t549 * t561) * qJD(4);
	t1 = [-t609, t561 * t575 + t565 * t571, -t529 * t565 - t543 * t587, -t527, 0; t528, t561 * t574 + t565 * t570, t531 * t565 + t573 * t587, t610, 0; 0, (t565 * t569 - t576 * t599) * t559, t538 * t565 - t545 * t587, t565 * t579 - t539 * t561 + (-t546 * t565 + t559 * t599) * qJD(4), 0; t531, t534 * t562 - t549 * t589, t530, 0, 0; t529, -t536 * t562 - t547 * t589, -t533, 0, 0; 0, (-t562 * t592 + t566 * t588) * t559, t539, 0, 0; t610, t561 * t571 - t565 * t575, -t529 * t561 + t543 * t586, t528, 0; t527, t561 * t570 - t565 * t574, t531 * t561 - t573 * t586, t609, 0; 0, (t561 * t569 + t576 * t596) * t559, t538 * t561 + t545 * t586, t561 * t579 + t539 * t565 + (-t546 * t561 - t559 * t596) * qJD(4), 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end