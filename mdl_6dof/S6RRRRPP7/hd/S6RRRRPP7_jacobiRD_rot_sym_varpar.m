% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRPP7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:29:20
	% EndTime: 2019-10-10 12:29:20
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:29:20
	% EndTime: 2019-10-10 12:29:20
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
	% StartTime: 2019-10-10 12:29:20
	% EndTime: 2019-10-10 12:29:21
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
	% StartTime: 2019-10-10 12:29:21
	% EndTime: 2019-10-10 12:29:22
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
	% StartTime: 2019-10-10 12:29:23
	% EndTime: 2019-10-10 12:29:24
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
	% StartTime: 2019-10-10 12:29:23
	% EndTime: 2019-10-10 12:29:24
	% DurationCPUTime: 0.61s
	% Computational Cost: add. (371->74), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->67)
	t479 = sin(qJ(1));
	t476 = cos(pkin(6));
	t491 = qJD(2) * t476 + qJD(1);
	t478 = sin(qJ(2));
	t512 = t479 * t478;
	t499 = t476 * t512;
	t507 = qJD(2) * t478;
	t481 = cos(qJ(2));
	t482 = cos(qJ(1));
	t509 = t482 * t481;
	t450 = -qJD(1) * t499 - t479 * t507 + t491 * t509;
	t477 = sin(qJ(3));
	t480 = cos(qJ(3));
	t510 = t482 * t478;
	t511 = t479 * t481;
	t461 = t476 * t510 + t511;
	t475 = sin(pkin(6));
	t513 = t475 * t482;
	t486 = t461 * t477 + t480 * t513;
	t508 = qJD(1) * t475;
	t496 = t479 * t508;
	t446 = t486 * qJD(3) - t450 * t480 - t477 * t496;
	t462 = t476 * t511 + t510;
	t449 = t462 * qJD(1) + t461 * qJD(2);
	t498 = t477 * t513;
	t455 = -t461 * t480 + t498;
	t497 = t476 * t509;
	t460 = -t497 + t512;
	t474 = qJ(4) + pkin(11);
	t472 = sin(t474);
	t473 = cos(t474);
	t525 = t446 * t473 - t449 * t472 + (-t455 * t472 - t460 * t473) * qJD(4);
	t524 = (t455 * t473 - t460 * t472) * qJD(4) + t446 * t472 + t449 * t473;
	t503 = qJD(3) * t481;
	t521 = (qJD(2) * t480 - qJD(4)) * t478 + t477 * t503;
	t516 = t475 * t477;
	t515 = t475 * t480;
	t514 = t475 * t481;
	t506 = qJD(2) * t481;
	t505 = qJD(3) * t477;
	t504 = qJD(3) * t480;
	t502 = qJD(4) * t472;
	t501 = qJD(4) * t473;
	t500 = qJD(4) * t480;
	t495 = t482 * t508;
	t494 = t475 * t507;
	t493 = t475 * t506;
	t448 = -t461 * qJD(1) - t462 * qJD(2);
	t489 = t462 * t500 + t448;
	t488 = t460 * t500 + t450;
	t487 = (qJD(2) - t500) * t481;
	t463 = -t499 + t509;
	t456 = -t463 * t477 + t479 * t515;
	t457 = t463 * t480 + t479 * t516;
	t459 = t476 * t477 + t478 * t515;
	t458 = t476 * t480 - t478 * t516;
	t444 = qJD(3) * t498 - t450 * t477 - t461 * t504 + t480 * t496;
	t447 = -qJD(1) * t497 - t482 * t506 + t491 * t512;
	t484 = qJD(4) * t463 + t447 * t480 + t462 * t505;
	t483 = qJD(4) * t461 - t449 * t480 + t460 * t505;
	t452 = t458 * qJD(3) + t480 * t493;
	t451 = -t459 * qJD(3) - t477 * t493;
	t443 = t456 * qJD(3) + t448 * t480 + t477 * t495;
	t442 = t457 * qJD(3) + t448 * t477 - t480 * t495;
	t441 = t443 * t473 - t447 * t472 + (-t457 * t472 + t462 * t473) * qJD(4);
	t440 = -t443 * t472 - t447 * t473 + (-t457 * t473 - t462 * t472) * qJD(4);
	t1 = [t525, t489 * t472 + t484 * t473, -t442 * t473 - t456 * t502, t440, 0, 0; t441, t488 * t472 + t483 * t473, t444 * t473 + t486 * t502, t524, 0, 0; 0, (t472 * t487 - t521 * t473) * t475, t451 * t473 - t458 * t502, t473 * t494 - t452 * t472 + (-t459 * t473 + t472 * t514) * qJD(4), 0, 0; -t524, -t484 * t472 + t489 * t473, t442 * t472 - t456 * t501, -t441, 0, 0; t440, -t483 * t472 + t488 * t473, -t444 * t472 + t486 * t501, t525, 0, 0; 0, (t521 * t472 + t473 * t487) * t475, -t451 * t472 - t458 * t501, -t472 * t494 - t452 * t473 + (t459 * t472 + t473 * t514) * qJD(4), 0, 0; t444, t447 * t477 - t462 * t504, t443, 0, 0, 0; t442, -t449 * t477 - t460 * t504, -t446, 0, 0, 0; 0, (-t477 * t507 + t480 * t503) * t475, t452, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:29:26
	% EndTime: 2019-10-10 12:29:27
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (371->74), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->67)
	t586 = sin(qJ(1));
	t583 = cos(pkin(6));
	t598 = qJD(2) * t583 + qJD(1);
	t585 = sin(qJ(2));
	t618 = t586 * t585;
	t605 = t583 * t618;
	t613 = qJD(2) * t585;
	t588 = cos(qJ(2));
	t589 = cos(qJ(1));
	t615 = t589 * t588;
	t557 = -qJD(1) * t605 - t586 * t613 + t598 * t615;
	t584 = sin(qJ(3));
	t587 = cos(qJ(3));
	t616 = t589 * t585;
	t617 = t586 * t588;
	t568 = t583 * t616 + t617;
	t582 = sin(pkin(6));
	t619 = t582 * t589;
	t594 = t568 * t584 + t587 * t619;
	t614 = qJD(1) * t582;
	t602 = t586 * t614;
	t553 = t594 * qJD(3) - t557 * t587 - t584 * t602;
	t569 = t583 * t617 + t616;
	t556 = t569 * qJD(1) + t568 * qJD(2);
	t604 = t584 * t619;
	t562 = -t568 * t587 + t604;
	t603 = t583 * t615;
	t567 = -t603 + t618;
	t581 = qJ(4) + pkin(11);
	t579 = sin(t581);
	t580 = cos(t581);
	t631 = t553 * t579 + t556 * t580 + (t562 * t580 - t567 * t579) * qJD(4);
	t630 = (t562 * t579 + t567 * t580) * qJD(4) - t553 * t580 + t556 * t579;
	t606 = qJD(4) * t587;
	t627 = t588 * (-qJD(2) + t606);
	t622 = t582 * t584;
	t621 = t582 * t587;
	t620 = t582 * t588;
	t612 = qJD(2) * t588;
	t611 = qJD(3) * t584;
	t610 = qJD(3) * t587;
	t609 = qJD(3) * t588;
	t608 = qJD(4) * t579;
	t607 = qJD(4) * t580;
	t601 = t589 * t614;
	t600 = t582 * t613;
	t599 = t582 * t612;
	t555 = -t568 * qJD(1) - t569 * qJD(2);
	t596 = t569 * t606 + t555;
	t595 = t567 * t606 + t557;
	t570 = -t605 + t615;
	t563 = -t570 * t584 + t586 * t621;
	t564 = t570 * t587 + t586 * t622;
	t566 = t583 * t584 + t585 * t621;
	t565 = t583 * t587 - t585 * t622;
	t551 = qJD(3) * t604 - t557 * t584 - t568 * t610 + t587 * t602;
	t554 = -qJD(1) * t603 - t589 * t612 + t598 * t618;
	t592 = qJD(4) * t570 + t554 * t587 + t569 * t611;
	t591 = qJD(4) * t568 - t556 * t587 + t567 * t611;
	t590 = -t584 * t609 + (-qJD(2) * t587 + qJD(4)) * t585;
	t559 = t565 * qJD(3) + t587 * t599;
	t558 = -t566 * qJD(3) - t584 * t599;
	t550 = t563 * qJD(3) + t555 * t587 + t584 * t601;
	t549 = t564 * qJD(3) + t555 * t584 - t587 * t601;
	t548 = t550 * t580 - t554 * t579 + (-t564 * t579 + t569 * t580) * qJD(4);
	t547 = t550 * t579 + t554 * t580 + (t564 * t580 + t569 * t579) * qJD(4);
	t1 = [-t630, t596 * t579 + t592 * t580, -t549 * t580 - t563 * t608, -t547, 0, 0; t548, t595 * t579 + t591 * t580, t551 * t580 + t594 * t608, t631, 0, 0; 0, (-t579 * t627 + t590 * t580) * t582, t558 * t580 - t565 * t608, t580 * t600 - t559 * t579 + (-t566 * t580 + t579 * t620) * qJD(4), 0, 0; t551, t554 * t584 - t569 * t610, t550, 0, 0, 0; t549, -t556 * t584 - t567 * t610, -t553, 0, 0, 0; 0, (-t584 * t613 + t587 * t609) * t582, t559, 0, 0, 0; t631, t592 * t579 - t596 * t580, -t549 * t579 + t563 * t607, t548, 0, 0; t547, t591 * t579 - t595 * t580, t551 * t579 - t594 * t607, t630, 0, 0; 0, (t590 * t579 + t580 * t627) * t582, t558 * t579 + t565 * t607, t579 * t600 + t559 * t580 + (-t566 * t579 - t580 * t620) * qJD(4), 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end