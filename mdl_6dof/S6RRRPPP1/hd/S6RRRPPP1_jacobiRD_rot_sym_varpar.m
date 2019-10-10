% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:15
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPPP1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPP1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
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
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t34 = sin(qJ(1));
	t41 = qJD(1) * t34;
	t36 = cos(qJ(1));
	t40 = qJD(1) * t36;
	t33 = sin(qJ(2));
	t39 = qJD(2) * t33;
	t35 = cos(qJ(2));
	t38 = qJD(2) * t35;
	t37 = qJD(2) * t36;
	t32 = t34 * t39 - t35 * t40;
	t31 = t33 * t40 + t34 * t38;
	t30 = t33 * t37 + t35 * t41;
	t29 = t33 * t41 - t35 * t37;
	t1 = [t32, t29, 0, 0, 0, 0; -t30, -t31, 0, 0, 0, 0; 0, -t39, 0, 0, 0, 0; t31, t30, 0, 0, 0, 0; t29, t32, 0, 0, 0, 0; 0, -t38, 0, 0, 0, 0; -t41, 0, 0, 0, 0, 0; t40, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:13
	% EndTime: 2019-10-10 11:15:14
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (48->26), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t232 = cos(qJ(3));
	t234 = cos(qJ(1));
	t256 = t232 * t234;
	t231 = sin(qJ(1));
	t255 = qJD(1) * t231;
	t233 = cos(qJ(2));
	t254 = qJD(1) * t233;
	t253 = qJD(1) * t234;
	t230 = sin(qJ(2));
	t252 = qJD(2) * t230;
	t251 = qJD(2) * t233;
	t250 = qJD(2) * t234;
	t229 = sin(qJ(3));
	t249 = qJD(3) * t229;
	t248 = qJD(3) * t230;
	t247 = qJD(3) * t233;
	t246 = t232 * t252;
	t245 = t232 * t248;
	t244 = t231 * t252;
	t243 = t231 * t251;
	t242 = t230 * t250;
	t241 = t233 * t250;
	t240 = -qJD(1) + t247;
	t239 = -qJD(3) + t254;
	t238 = t240 * t229;
	t237 = t230 * t253 + t243;
	t236 = -t230 * t255 + t241;
	t235 = t239 * t231 + t242;
	t228 = -t239 * t256 + (t238 + t246) * t231;
	t227 = t240 * t232 * t231 + (t239 * t234 - t244) * t229;
	t226 = t235 * t232 + t234 * t238;
	t225 = t235 * t229 - t240 * t256;
	t1 = [t228, -t232 * t241 + (t232 * t255 + t234 * t249) * t230, t225, 0, 0, 0; -t226, -t232 * t243 + (t231 * t249 - t232 * t253) * t230, -t227, 0, 0, 0; 0, -t229 * t247 - t246, -t229 * t251 - t245, 0, 0, 0; t227, t236 * t229 + t234 * t245, t226, 0, 0, 0; t225, t237 * t229 + t231 * t245, t228, 0, 0, 0; 0, t229 * t252 - t232 * t247, t229 * t248 - t232 * t251, 0, 0, 0; -t237, -t231 * t254 - t242, 0, 0, 0, 0; t236, t233 * t253 - t244, 0, 0, 0, 0; 0, t251, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:16
	% EndTime: 2019-10-10 11:15:16
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (155->57), mult. (584->124), div. (0->0), fcn. (584->10), ass. (0->57)
	t401 = sin(qJ(3));
	t402 = sin(qJ(2));
	t405 = cos(qJ(2));
	t406 = cos(qJ(1));
	t437 = qJD(2) * t406;
	t429 = t405 * t437;
	t403 = sin(qJ(1));
	t442 = qJD(1) * t403;
	t416 = t402 * t442 - t429;
	t433 = qJD(3) * t406;
	t404 = cos(qJ(3));
	t444 = t402 * t404;
	t453 = t416 * t401 - t433 * t444;
	t438 = qJD(2) * t405;
	t440 = qJD(1) * t406;
	t417 = t402 * t440 + t403 * t438;
	t435 = qJD(3) * t403;
	t428 = t404 * t435;
	t452 = t417 * t401 + t402 * t428;
	t427 = t401 * t433;
	t414 = t404 * t442 + t427;
	t451 = t414 * t402 - t404 * t429;
	t430 = t404 * t438;
	t450 = (-t401 * t435 + t404 * t440) * t402 + t403 * t430;
	t397 = sin(pkin(10));
	t400 = cos(pkin(6));
	t449 = t397 * t400;
	t399 = cos(pkin(10));
	t448 = t399 * t400;
	t447 = t400 * t401;
	t446 = t400 * t404;
	t445 = t401 * t402;
	t443 = t404 * t406;
	t441 = qJD(1) * t405;
	t439 = qJD(2) * t402;
	t436 = qJD(3) * t402;
	t434 = qJD(3) * t405;
	t432 = t403 * t439;
	t431 = t402 * t437;
	t426 = -qJD(1) + t434;
	t425 = -qJD(3) + t441;
	t398 = sin(pkin(6));
	t420 = t398 * t405 + t400 * t445;
	t419 = t405 * t440 - t432;
	t418 = t403 * t441 + t431;
	t413 = t401 * t440 + t428;
	t412 = (t397 * t401 - t399 * t446) * t405;
	t411 = (-t397 * t446 - t399 * t401) * t405;
	t388 = -t426 * t443 + (t425 * t403 + t431) * t401;
	t410 = t388 * t400 - t416 * t398;
	t390 = -t401 * t432 + t413 * t405 - t414;
	t409 = t390 * t400 - t417 * t398;
	t408 = t419 * t398 + t452 * t400;
	t407 = -t418 * t398 - t453 * t400;
	t391 = -t425 * t443 + (t426 * t401 + t404 * t439) * t403;
	t389 = t418 * t404 + t405 * t427 - t413;
	t1 = [t391 * t399 + t409 * t397, t407 * t397 + t451 * t399, t388 * t399 + t389 * t449, 0, 0, 0; -t389 * t399 + t410 * t397, t408 * t397 - t450 * t399, -t390 * t399 + t391 * t449, 0, 0, 0; 0, qJD(3) * t411 + (t420 * t397 - t399 * t444) * qJD(2), (t397 * t447 - t399 * t404) * t436 + qJD(2) * t411, 0, 0, 0; -t391 * t397 + t409 * t399, -t451 * t397 + t407 * t399, -t388 * t397 + t389 * t448, 0, 0, 0; t389 * t397 + t410 * t399, t450 * t397 + t408 * t399, t390 * t397 + t391 * t448, 0, 0, 0; 0, qJD(3) * t412 + (t397 * t444 + t420 * t399) * qJD(2), (t397 * t404 + t399 * t447) * t436 + qJD(2) * t412, 0, 0, 0; -t390 * t398 - t417 * t400, t453 * t398 - t418 * t400, -t389 * t398, 0, 0, 0; -t388 * t398 - t416 * t400, -t452 * t398 + t419 * t400, -t391 * t398, 0, 0, 0; 0, t404 * t398 * t434 + (-t398 * t445 + t400 * t405) * qJD(2), (-t401 * t436 + t430) * t398, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:17
	% EndTime: 2019-10-10 11:15:18
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (155->57), mult. (584->124), div. (0->0), fcn. (584->10), ass. (0->57)
	t509 = sin(qJ(2));
	t511 = cos(qJ(3));
	t508 = sin(qJ(3));
	t513 = cos(qJ(1));
	t540 = qJD(3) * t513;
	t534 = t508 * t540;
	t510 = sin(qJ(1));
	t549 = qJD(1) * t510;
	t523 = t511 * t549 + t534;
	t512 = cos(qJ(2));
	t544 = qJD(2) * t513;
	t536 = t512 * t544;
	t558 = t523 * t509 - t511 * t536;
	t545 = qJD(2) * t512;
	t537 = t511 * t545;
	t542 = qJD(3) * t510;
	t547 = qJD(1) * t513;
	t557 = (-t508 * t542 + t511 * t547) * t509 + t510 * t537;
	t504 = sin(pkin(10));
	t507 = cos(pkin(6));
	t556 = t504 * t507;
	t506 = cos(pkin(10));
	t555 = t506 * t507;
	t554 = t507 * t508;
	t553 = t507 * t511;
	t552 = t508 * t509;
	t551 = t509 * t511;
	t550 = t511 * t513;
	t548 = qJD(1) * t512;
	t546 = qJD(2) * t509;
	t543 = qJD(3) * t509;
	t541 = qJD(3) * t512;
	t539 = t510 * t546;
	t538 = t509 * t544;
	t535 = t511 * t542;
	t533 = -qJD(1) + t541;
	t532 = -qJD(3) + t548;
	t505 = sin(pkin(6));
	t529 = -t505 * t512 - t507 * t552;
	t528 = t512 * t547 - t539;
	t527 = t510 * t548 + t538;
	t526 = t509 * t547 + t510 * t545;
	t525 = t509 * t549 - t536;
	t522 = t508 * t547 + t535;
	t521 = (-t504 * t508 + t506 * t553) * t512;
	t520 = (t504 * t553 + t506 * t508) * t512;
	t495 = -t533 * t550 + (t532 * t510 + t538) * t508;
	t519 = -t495 * t507 + t525 * t505;
	t497 = -t508 * t539 + t522 * t512 - t523;
	t518 = -t497 * t507 + t526 * t505;
	t517 = -t526 * t508 - t509 * t535;
	t516 = t525 * t508 - t540 * t551;
	t515 = -t528 * t505 + t517 * t507;
	t514 = t527 * t505 + t516 * t507;
	t498 = -t532 * t550 + (t533 * t508 + t511 * t546) * t510;
	t496 = t527 * t511 + t512 * t534 - t522;
	t1 = [-t497 * t505 - t526 * t507, t516 * t505 - t527 * t507, -t496 * t505, 0, 0, 0; -t495 * t505 - t525 * t507, t517 * t505 + t528 * t507, -t498 * t505, 0, 0, 0; 0, t511 * t505 * t541 + (-t505 * t552 + t507 * t512) * qJD(2), (-t508 * t543 + t537) * t505, 0, 0, 0; -t498 * t506 + t518 * t504, t514 * t504 - t558 * t506, -t495 * t506 - t496 * t556, 0, 0, 0; t496 * t506 + t519 * t504, t515 * t504 + t557 * t506, t497 * t506 - t498 * t556, 0, 0, 0; 0, qJD(3) * t520 + (t529 * t504 + t506 * t551) * qJD(2), (-t504 * t554 + t506 * t511) * t543 + qJD(2) * t520, 0, 0, 0; t498 * t504 + t518 * t506, t558 * t504 + t514 * t506, t495 * t504 - t496 * t555, 0, 0, 0; -t496 * t504 + t519 * t506, -t557 * t504 + t515 * t506, -t497 * t504 - t498 * t555, 0, 0, 0; 0, qJD(3) * t521 + (-t504 * t551 + t529 * t506) * qJD(2), (-t504 * t511 - t506 * t554) * t543 + qJD(2) * t521, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:17
	% EndTime: 2019-10-10 11:15:18
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (155->57), mult. (584->124), div. (0->0), fcn. (584->10), ass. (0->57)
	t547 = sin(pkin(6));
	t549 = cos(pkin(6));
	t554 = cos(qJ(2));
	t552 = sin(qJ(1));
	t551 = sin(qJ(2));
	t585 = qJD(2) * t551;
	t578 = t552 * t585;
	t555 = cos(qJ(1));
	t586 = qJD(1) * t555;
	t567 = t554 * t586 - t578;
	t550 = sin(qJ(3));
	t584 = qJD(2) * t554;
	t565 = t551 * t586 + t552 * t584;
	t553 = cos(qJ(3));
	t581 = qJD(3) * t552;
	t574 = t553 * t581;
	t598 = t565 * t550 + t551 * t574;
	t605 = t567 * t547 + t598 * t549;
	t583 = qJD(2) * t555;
	t575 = t554 * t583;
	t588 = qJD(1) * t552;
	t564 = t551 * t588 - t575;
	t579 = qJD(3) * t555;
	t590 = t551 * t553;
	t556 = t564 * t550 - t579 * t590;
	t577 = t551 * t583;
	t587 = qJD(1) * t554;
	t566 = t552 * t587 + t577;
	t604 = t566 * t547 + t556 * t549;
	t571 = -qJD(3) + t587;
	t580 = qJD(3) * t554;
	t572 = -qJD(1) + t580;
	t589 = t553 * t555;
	t537 = -t572 * t589 + (t571 * t552 + t577) * t550;
	t600 = -t537 * t549 + t564 * t547;
	t562 = t550 * t586 + t574;
	t573 = t550 * t579;
	t563 = t553 * t588 + t573;
	t539 = -t550 * t578 + t562 * t554 - t563;
	t599 = -t539 * t549 + t565 * t547;
	t546 = sin(pkin(10));
	t595 = t546 * t549;
	t548 = cos(pkin(10));
	t594 = t548 * t549;
	t593 = t549 * t550;
	t592 = t549 * t553;
	t591 = t550 * t551;
	t582 = qJD(3) * t551;
	t576 = t553 * t584;
	t568 = t547 * t554 + t549 * t591;
	t561 = (-t546 * t550 + t548 * t592) * t554;
	t560 = (-t546 * t592 - t548 * t550) * t554;
	t559 = -t552 * t576 + (t550 * t581 - t553 * t586) * t551;
	t558 = t563 * t551 - t553 * t575;
	t540 = -t571 * t589 + (t572 * t550 + t553 * t585) * t552;
	t538 = t566 * t553 + t554 * t573 - t562;
	t1 = [-t539 * t547 - t565 * t549, t556 * t547 - t566 * t549, -t538 * t547, 0, 0, 0; -t537 * t547 - t564 * t549, -t547 * t598 + t567 * t549, -t540 * t547, 0, 0, 0; 0, t553 * t547 * t580 + (-t547 * t591 + t549 * t554) * qJD(2), (-t550 * t582 + t576) * t547, 0, 0, 0; t540 * t546 + t599 * t548, t558 * t546 + t604 * t548, t537 * t546 - t538 * t594, 0, 0, 0; -t538 * t546 + t600 * t548, t559 * t546 - t605 * t548, -t539 * t546 - t540 * t594, 0, 0, 0; 0, qJD(3) * t561 + (-t546 * t590 - t568 * t548) * qJD(2), (-t546 * t553 - t548 * t593) * t582 + qJD(2) * t561, 0, 0, 0; t540 * t548 - t599 * t546, -t604 * t546 + t558 * t548, t537 * t548 + t538 * t595, 0, 0, 0; -t538 * t548 - t600 * t546, t605 * t546 + t559 * t548, -t539 * t548 + t540 * t595, 0, 0, 0; 0, qJD(3) * t560 + (t568 * t546 - t548 * t590) * qJD(2), (t546 * t593 - t548 * t553) * t582 + qJD(2) * t560, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end