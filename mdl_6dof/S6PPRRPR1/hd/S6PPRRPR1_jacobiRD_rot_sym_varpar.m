% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PPRRPR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PPRRPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (20->14), mult. (80->39), div. (0->0), fcn. (88->10), ass. (0->20)
	t132 = sin(pkin(11));
	t138 = cos(pkin(6));
	t147 = t132 * t138;
	t133 = sin(pkin(7));
	t134 = sin(pkin(6));
	t146 = t133 * t134;
	t145 = t133 * t138;
	t135 = cos(pkin(12));
	t137 = cos(pkin(7));
	t144 = t135 * t137;
	t136 = cos(pkin(11));
	t143 = t136 * t138;
	t131 = sin(pkin(12));
	t142 = -(-t132 * t131 + t135 * t143) * t137 + t136 * t146;
	t141 = -(-t136 * t131 - t135 * t147) * t137 - t132 * t146;
	t140 = cos(qJ(3));
	t139 = sin(qJ(3));
	t130 = -t131 * t147 + t136 * t135;
	t128 = t131 * t143 + t132 * t135;
	t1 = [0, 0, (-t130 * t140 + t141 * t139) * qJD(3), 0, 0, 0; 0, 0, (-t128 * t140 + t142 * t139) * qJD(3), 0, 0, 0; 0, 0, (-t139 * t145 + (-t131 * t140 - t139 * t144) * t134) * qJD(3), 0, 0, 0; 0, 0, (t130 * t139 + t141 * t140) * qJD(3), 0, 0, 0; 0, 0, (t128 * t139 + t142 * t140) * qJD(3), 0, 0, 0; 0, 0, (-t140 * t145 + (t131 * t139 - t140 * t144) * t134) * qJD(3), 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:41
	% EndTime: 2019-10-09 21:10:42
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (118->35), mult. (402->84), div. (0->0), fcn. (472->12), ass. (0->42)
	t336 = sin(pkin(12));
	t339 = sin(pkin(6));
	t345 = sin(qJ(3));
	t347 = cos(qJ(3));
	t340 = cos(pkin(12));
	t342 = cos(pkin(7));
	t355 = t340 * t342;
	t338 = sin(pkin(7));
	t343 = cos(pkin(6));
	t357 = t338 * t343;
	t328 = (t336 * t347 + t345 * t355) * t339 + t345 * t357;
	t341 = cos(pkin(11));
	t337 = sin(pkin(11));
	t359 = t337 * t343;
	t335 = -t336 * t359 + t341 * t340;
	t334 = -t341 * t336 - t340 * t359;
	t358 = t338 * t339;
	t349 = t334 * t342 + t337 * t358;
	t324 = t335 * t347 + t349 * t345;
	t354 = t341 * t343;
	t333 = t336 * t354 + t337 * t340;
	t332 = -t337 * t336 + t340 * t354;
	t350 = -t332 * t342 + t341 * t358;
	t362 = -t333 * t347 + t350 * t345;
	t356 = t339 * t342;
	t344 = sin(qJ(4));
	t353 = qJD(4) * t344;
	t346 = cos(qJ(4));
	t352 = qJD(4) * t346;
	t321 = -t333 * t345 - t350 * t347;
	t323 = -t335 * t345 + t349 * t347;
	t327 = t347 * t357 + (-t336 * t345 + t347 * t355) * t339;
	t331 = -t340 * t358 + t343 * t342;
	t330 = -t334 * t338 + t337 * t356;
	t329 = -t332 * t338 - t341 * t356;
	t326 = t328 * qJD(3);
	t325 = t327 * qJD(3);
	t320 = t324 * qJD(3);
	t319 = t323 * qJD(3);
	t318 = t362 * qJD(3);
	t317 = t321 * qJD(3);
	t1 = [0, 0, -t320 * t346 - t323 * t353, -t319 * t344 + (-t324 * t346 - t330 * t344) * qJD(4), 0, 0; 0, 0, t318 * t346 - t321 * t353, -t317 * t344 + (-t329 * t344 + t346 * t362) * qJD(4), 0, 0; 0, 0, -t326 * t346 - t327 * t353, -t325 * t344 + (-t328 * t346 - t331 * t344) * qJD(4), 0, 0; 0, 0, t320 * t344 - t323 * t352, -t319 * t346 + (t324 * t344 - t330 * t346) * qJD(4), 0, 0; 0, 0, -t318 * t344 - t321 * t352, -t317 * t346 + (-t329 * t346 - t344 * t362) * qJD(4), 0, 0; 0, 0, t326 * t344 - t327 * t352, -t325 * t346 + (t328 * t344 - t331 * t346) * qJD(4), 0, 0; 0, 0, t319, 0, 0, 0; 0, 0, t317, 0, 0, 0; 0, 0, t325, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:43
	% EndTime: 2019-10-09 21:10:43
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (191->44), mult. (647->102), div. (0->0), fcn. (754->14), ass. (0->50)
	t425 = sin(pkin(12));
	t428 = sin(pkin(6));
	t435 = sin(qJ(3));
	t437 = cos(qJ(3));
	t430 = cos(pkin(12));
	t432 = cos(pkin(7));
	t448 = t430 * t432;
	t427 = sin(pkin(7));
	t433 = cos(pkin(6));
	t450 = t427 * t433;
	t416 = (t425 * t437 + t435 * t448) * t428 + t435 * t450;
	t431 = cos(pkin(11));
	t426 = sin(pkin(11));
	t452 = t426 * t433;
	t423 = -t425 * t452 + t431 * t430;
	t422 = -t431 * t425 - t430 * t452;
	t451 = t427 * t428;
	t442 = t422 * t432 + t426 * t451;
	t412 = t423 * t437 + t442 * t435;
	t447 = t431 * t433;
	t421 = t425 * t447 + t426 * t430;
	t420 = -t426 * t425 + t430 * t447;
	t443 = -t420 * t432 + t431 * t451;
	t455 = -t421 * t437 + t443 * t435;
	t449 = t428 * t432;
	t434 = sin(qJ(4));
	t446 = qJD(4) * t434;
	t436 = cos(qJ(4));
	t445 = qJD(4) * t436;
	t406 = t455 * qJD(3);
	t409 = -t421 * t435 - t443 * t437;
	t440 = -t406 * t436 + t409 * t446;
	t408 = t412 * qJD(3);
	t411 = -t423 * t435 + t442 * t437;
	t439 = t408 * t436 + t411 * t446;
	t414 = t416 * qJD(3);
	t415 = t437 * t450 + (-t425 * t435 + t437 * t448) * t428;
	t438 = t414 * t436 + t415 * t446;
	t429 = cos(pkin(13));
	t424 = sin(pkin(13));
	t419 = -t430 * t451 + t433 * t432;
	t418 = -t422 * t427 + t426 * t449;
	t417 = -t420 * t427 - t431 * t449;
	t413 = t415 * qJD(3);
	t407 = t411 * qJD(3);
	t405 = t409 * qJD(3);
	t404 = -t413 * t434 + (-t416 * t436 - t419 * t434) * qJD(4);
	t403 = -t407 * t434 + (-t412 * t436 - t418 * t434) * qJD(4);
	t402 = -t405 * t434 + (-t417 * t434 + t436 * t455) * qJD(4);
	t1 = [0, 0, t407 * t424 - t439 * t429, t403 * t429, 0, 0; 0, 0, t405 * t424 - t440 * t429, t402 * t429, 0, 0; 0, 0, t413 * t424 - t438 * t429, t404 * t429, 0, 0; 0, 0, t407 * t429 + t439 * t424, -t403 * t424, 0, 0; 0, 0, t405 * t429 + t440 * t424, -t402 * t424, 0, 0; 0, 0, t413 * t429 + t438 * t424, -t404 * t424, 0, 0; 0, 0, -t408 * t434 + t411 * t445, t407 * t436 + (-t412 * t434 + t418 * t436) * qJD(4), 0, 0; 0, 0, t406 * t434 + t409 * t445, t405 * t436 + (t417 * t436 + t434 * t455) * qJD(4), 0, 0; 0, 0, -t414 * t434 + t415 * t445, t413 * t436 + (-t416 * t434 + t419 * t436) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:44
	% EndTime: 2019-10-09 21:10:44
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (504->74), mult. (1401->151), div. (0->0), fcn. (1686->14), ass. (0->67)
	t530 = sin(pkin(12));
	t533 = sin(pkin(6));
	t537 = cos(pkin(6));
	t539 = sin(qJ(3));
	t532 = sin(pkin(7));
	t568 = cos(qJ(3));
	t554 = t532 * t568;
	t534 = cos(pkin(12));
	t536 = cos(pkin(7));
	t563 = t534 * t536;
	t570 = (-t530 * t539 + t568 * t563) * t533 + t537 * t554;
	t531 = sin(pkin(11));
	t535 = cos(pkin(11));
	t562 = t535 * t537;
	t521 = t530 * t562 + t531 * t534;
	t547 = -t531 * t530 + t534 * t562;
	t545 = t547 * t536;
	t565 = t533 * t532;
	t507 = t521 * t568 + (-t535 * t565 + t545) * t539;
	t566 = t531 * t537;
	t564 = t533 * t536;
	t538 = sin(qJ(4));
	t561 = qJD(4) * t538;
	t540 = cos(qJ(4));
	t560 = qJD(4) * t540;
	t529 = pkin(13) + qJ(6);
	t527 = sin(t529);
	t559 = qJD(6) * t527;
	t528 = cos(t529);
	t558 = qJD(6) * t528;
	t557 = qJD(6) * t540;
	t553 = t533 * t554;
	t506 = t521 * t539 + t535 * t553 - t568 * t545;
	t502 = t506 * qJD(3);
	t550 = t506 * t557 - t502;
	t522 = -t530 * t566 + t535 * t534;
	t546 = t535 * t530 + t534 * t566;
	t544 = t546 * t536;
	t508 = t522 * t539 - t531 * t553 + t568 * t544;
	t504 = t508 * qJD(3);
	t549 = t508 * t557 - t504;
	t512 = t570 * qJD(3);
	t548 = -t557 * t570 + t512;
	t516 = -t547 * t532 - t535 * t564;
	t499 = t507 * t540 + t516 * t538;
	t498 = -t507 * t538 + t516 * t540;
	t509 = t522 * t568 + (t531 * t565 - t544) * t539;
	t517 = t531 * t564 + t546 * t532;
	t501 = t509 * t540 + t517 * t538;
	t500 = -t509 * t538 + t517 * t540;
	t515 = t537 * t532 * t539 + (t568 * t530 + t539 * t563) * t533;
	t520 = -t534 * t565 + t537 * t536;
	t511 = t515 * t540 + t520 * t538;
	t510 = -t515 * t538 + t520 * t540;
	t503 = t507 * qJD(3);
	t543 = qJD(6) * t507 - t503 * t540 + t506 * t561;
	t505 = t509 * qJD(3);
	t542 = qJD(6) * t509 - t505 * t540 + t508 * t561;
	t513 = t515 * qJD(3);
	t541 = qJD(6) * t515 - t513 * t540 - t561 * t570;
	t497 = t510 * qJD(4) + t512 * t540;
	t496 = -t511 * qJD(4) - t512 * t538;
	t495 = t500 * qJD(4) - t504 * t540;
	t494 = -t501 * qJD(4) + t504 * t538;
	t493 = t498 * qJD(4) - t502 * t540;
	t492 = -t499 * qJD(4) + t502 * t538;
	t1 = [0, 0, t549 * t527 + t542 * t528, t494 * t528 - t500 * t559, 0, -t495 * t527 + t505 * t528 + (-t501 * t528 - t508 * t527) * qJD(6); 0, 0, t550 * t527 + t543 * t528, t492 * t528 - t498 * t559, 0, -t493 * t527 + t503 * t528 + (-t499 * t528 - t506 * t527) * qJD(6); 0, 0, t548 * t527 + t541 * t528, t496 * t528 - t510 * t559, 0, -t497 * t527 + t513 * t528 + (-t511 * t528 + t527 * t570) * qJD(6); 0, 0, -t542 * t527 + t549 * t528, -t494 * t527 - t500 * t558, 0, -t495 * t528 - t505 * t527 + (t501 * t527 - t508 * t528) * qJD(6); 0, 0, -t543 * t527 + t550 * t528, -t492 * t527 - t498 * t558, 0, -t493 * t528 - t503 * t527 + (t499 * t527 - t506 * t528) * qJD(6); 0, 0, -t541 * t527 + t548 * t528, -t496 * t527 - t510 * t558, 0, -t497 * t528 - t513 * t527 + (t511 * t527 + t528 * t570) * qJD(6); 0, 0, -t505 * t538 - t508 * t560, t495, 0, 0; 0, 0, -t503 * t538 - t506 * t560, t493, 0, 0; 0, 0, -t513 * t538 + t560 * t570, t497, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end