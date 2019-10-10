% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR2
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
% Datum: 2019-10-09 21:55
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRRR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:55:41
	% EndTime: 2019-10-09 21:55:41
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:55:41
	% EndTime: 2019-10-09 21:55:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:55:41
	% EndTime: 2019-10-09 21:55:41
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
	% StartTime: 2019-10-09 21:55:41
	% EndTime: 2019-10-09 21:55:42
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
	% StartTime: 2019-10-09 21:55:42
	% EndTime: 2019-10-09 21:55:43
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
	% StartTime: 2019-10-09 21:55:44
	% EndTime: 2019-10-09 21:55:45
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (303->66), mult. (966->138), div. (0->0), fcn. (1104->12), ass. (0->60)
	t463 = cos(pkin(6));
	t458 = sin(pkin(12));
	t461 = cos(pkin(12));
	t466 = sin(qJ(2));
	t469 = cos(qJ(2));
	t476 = t469 * t458 + t466 * t461;
	t448 = t476 * t463;
	t451 = t466 * t458 - t469 * t461;
	t449 = t451 * qJD(2);
	t460 = sin(pkin(6));
	t465 = sin(qJ(4));
	t494 = t460 * t465;
	t468 = cos(qJ(4));
	t493 = t460 * t468;
	t490 = qJD(2) * t466;
	t489 = qJD(2) * t469;
	t488 = qJD(4) * t465;
	t487 = qJD(4) * t468;
	t464 = sin(qJ(5));
	t486 = qJD(5) * t464;
	t467 = cos(qJ(5));
	t485 = qJD(5) * t467;
	t484 = qJD(5) * t468;
	t459 = sin(pkin(11));
	t462 = cos(pkin(11));
	t474 = t451 * t463;
	t434 = -t459 * t476 - t462 * t474;
	t445 = (t458 * t490 - t461 * t489) * t463;
	t450 = -t458 * t489 - t461 * t490;
	t480 = -t462 * t445 + t459 * t450;
	t483 = -t434 * t484 + t480;
	t437 = t459 * t474 - t462 * t476;
	t479 = t459 * t445 + t462 * t450;
	t482 = -t437 * t484 + t479;
	t444 = t460 * t449;
	t446 = t451 * t460;
	t481 = t446 * t484 - t444;
	t447 = t476 * t460;
	t440 = t447 * t468 + t463 * t465;
	t439 = -t447 * t465 + t463 * t468;
	t478 = t462 * t448 - t459 * t451;
	t477 = -t459 * t448 - t462 * t451;
	t423 = -t462 * t493 - t465 * t478;
	t475 = t462 * t494 - t468 * t478;
	t425 = t459 * t493 - t465 * t477;
	t426 = t459 * t494 + t468 * t477;
	t473 = qJD(2) * t448;
	t427 = t459 * t449 - t462 * t473;
	t472 = -qJD(5) * t478 - t427 * t468 + t434 * t488;
	t430 = t462 * t449 + t459 * t473;
	t471 = -qJD(5) * t477 - t430 * t468 + t437 * t488;
	t443 = qJD(2) * t447;
	t470 = qJD(5) * t447 - t443 * t468 + t446 * t488;
	t422 = t439 * qJD(4) - t444 * t468;
	t421 = -t440 * qJD(4) + t444 * t465;
	t420 = t425 * qJD(4) + t468 * t479;
	t419 = -t426 * qJD(4) - t465 * t479;
	t418 = t423 * qJD(4) + t468 * t480;
	t417 = t475 * qJD(4) - t465 * t480;
	t1 = [0, t482 * t464 - t471 * t467, 0, t419 * t467 - t425 * t486, -t420 * t464 - t430 * t467 + (-t426 * t467 + t437 * t464) * qJD(5), 0; 0, t483 * t464 - t472 * t467, 0, t417 * t467 - t423 * t486, -t418 * t464 - t427 * t467 + (t434 * t464 + t467 * t475) * qJD(5), 0; 0, t481 * t464 + t470 * t467, 0, t421 * t467 - t439 * t486, -t422 * t464 + t443 * t467 + (-t440 * t467 - t446 * t464) * qJD(5), 0; 0, t471 * t464 + t482 * t467, 0, -t419 * t464 - t425 * t485, -t420 * t467 + t430 * t464 + (t426 * t464 + t437 * t467) * qJD(5), 0; 0, t472 * t464 + t483 * t467, 0, -t417 * t464 - t423 * t485, -t418 * t467 + t427 * t464 + (t434 * t467 - t464 * t475) * qJD(5), 0; 0, -t470 * t464 + t481 * t467, 0, -t421 * t464 - t439 * t485, -t422 * t467 - t443 * t464 + (t440 * t464 - t446 * t467) * qJD(5), 0; 0, t430 * t465 + t437 * t487, 0, t420, 0, 0; 0, t427 * t465 + t434 * t487, 0, t418, 0, 0; 0, -t443 * t465 - t446 * t487, 0, t422, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:55:45
	% EndTime: 2019-10-09 21:55:45
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (541->62), mult. (1292->126), div. (0->0), fcn. (1478->12), ass. (0->74)
	t530 = cos(pkin(6));
	t525 = sin(pkin(12));
	t528 = cos(pkin(12));
	t532 = sin(qJ(2));
	t534 = cos(qJ(2));
	t541 = t534 * t525 + t532 * t528;
	t511 = t541 * t530;
	t514 = t532 * t525 - t534 * t528;
	t512 = t514 * qJD(2);
	t524 = qJ(5) + qJ(6);
	t521 = sin(t524);
	t523 = qJD(5) + qJD(6);
	t565 = t521 * t523;
	t522 = cos(t524);
	t564 = t522 * t523;
	t533 = cos(qJ(4));
	t563 = t523 * t533;
	t527 = sin(pkin(6));
	t531 = sin(qJ(4));
	t562 = t527 * t531;
	t561 = t527 * t533;
	t558 = qJD(2) * t532;
	t557 = qJD(2) * t534;
	t556 = qJD(4) * t531;
	t555 = qJD(4) * t533;
	t529 = cos(pkin(11));
	t526 = sin(pkin(11));
	t543 = t529 * t511 - t526 * t514;
	t486 = -t529 * t561 - t531 * t543;
	t508 = (t525 * t558 - t528 * t557) * t530;
	t513 = -t525 * t557 - t528 * t558;
	t545 = -t529 * t508 + t526 * t513;
	t481 = t486 * qJD(4) + t533 * t545;
	t539 = t514 * t530;
	t497 = -t526 * t541 - t529 * t539;
	t554 = t497 * t523 - t481;
	t542 = -t526 * t511 - t529 * t514;
	t488 = t526 * t561 - t531 * t542;
	t544 = t526 * t508 + t529 * t513;
	t483 = t488 * qJD(4) + t533 * t544;
	t500 = t526 * t539 - t529 * t541;
	t553 = t500 * t523 - t483;
	t510 = t541 * t527;
	t502 = -t510 * t531 + t530 * t533;
	t507 = t527 * t512;
	t485 = t502 * qJD(4) - t507 * t533;
	t509 = t514 * t527;
	t552 = -t509 * t523 - t485;
	t538 = qJD(2) * t511;
	t490 = t526 * t512 - t529 * t538;
	t540 = t529 * t562 - t533 * t543;
	t551 = -t523 * t540 + t490;
	t489 = t526 * t562 + t533 * t542;
	t493 = t529 * t512 + t526 * t538;
	t550 = t489 * t523 + t493;
	t503 = t510 * t533 + t530 * t531;
	t506 = qJD(2) * t510;
	t549 = t503 * t523 - t506;
	t548 = -t497 * t563 + t545;
	t547 = -t500 * t563 + t544;
	t546 = t509 * t563 - t507;
	t537 = -t490 * t533 + t497 * t556 - t523 * t543;
	t536 = -t493 * t533 + t500 * t556 - t523 * t542;
	t535 = -t506 * t533 + t509 * t556 + t510 * t523;
	t484 = -t503 * qJD(4) + t507 * t531;
	t482 = -t489 * qJD(4) - t531 * t544;
	t480 = t540 * qJD(4) - t531 * t545;
	t479 = t549 * t521 + t552 * t522;
	t478 = t552 * t521 - t549 * t522;
	t477 = t550 * t521 + t553 * t522;
	t476 = t553 * t521 - t550 * t522;
	t475 = t551 * t521 + t554 * t522;
	t474 = t554 * t521 - t551 * t522;
	t1 = [0, t547 * t521 - t536 * t522, 0, t482 * t522 - t488 * t565, t476, t476; 0, t548 * t521 - t537 * t522, 0, t480 * t522 - t486 * t565, t474, t474; 0, t546 * t521 + t535 * t522, 0, t484 * t522 - t502 * t565, t478, t478; 0, t536 * t521 + t547 * t522, 0, -t482 * t521 - t488 * t564, t477, t477; 0, t537 * t521 + t548 * t522, 0, -t480 * t521 - t486 * t564, t475, t475; 0, -t535 * t521 + t546 * t522, 0, -t484 * t521 - t502 * t564, t479, t479; 0, t493 * t531 + t500 * t555, 0, t483, 0, 0; 0, t490 * t531 + t497 * t555, 0, t481, 0, 0; 0, -t506 * t531 - t509 * t555, 0, t485, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end