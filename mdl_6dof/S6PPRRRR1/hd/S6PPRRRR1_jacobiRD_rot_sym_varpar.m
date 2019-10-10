% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:18
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PPRRRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:21
	% EndTime: 2019-10-09 21:18:21
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (20->14), mult. (80->39), div. (0->0), fcn. (88->10), ass. (0->20)
	t132 = sin(pkin(12));
	t138 = cos(pkin(6));
	t147 = t132 * t138;
	t133 = sin(pkin(7));
	t134 = sin(pkin(6));
	t146 = t133 * t134;
	t145 = t133 * t138;
	t135 = cos(pkin(13));
	t137 = cos(pkin(7));
	t144 = t135 * t137;
	t136 = cos(pkin(12));
	t143 = t136 * t138;
	t131 = sin(pkin(13));
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
	% StartTime: 2019-10-09 21:18:22
	% EndTime: 2019-10-09 21:18:22
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (118->35), mult. (402->84), div. (0->0), fcn. (472->12), ass. (0->42)
	t336 = sin(pkin(13));
	t339 = sin(pkin(6));
	t345 = sin(qJ(3));
	t347 = cos(qJ(3));
	t340 = cos(pkin(13));
	t342 = cos(pkin(7));
	t355 = t340 * t342;
	t338 = sin(pkin(7));
	t343 = cos(pkin(6));
	t357 = t338 * t343;
	t328 = (t336 * t347 + t345 * t355) * t339 + t345 * t357;
	t341 = cos(pkin(12));
	t337 = sin(pkin(12));
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
	% StartTime: 2019-10-09 21:18:22
	% EndTime: 2019-10-09 21:18:23
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (258->34), mult. (604->75), div. (0->0), fcn. (712->12), ass. (0->50)
	t398 = sin(pkin(13));
	t401 = sin(pkin(6));
	t406 = sin(qJ(3));
	t407 = cos(qJ(3));
	t402 = cos(pkin(13));
	t404 = cos(pkin(7));
	t416 = t402 * t404;
	t400 = sin(pkin(7));
	t405 = cos(pkin(6));
	t418 = t400 * t405;
	t386 = (t398 * t407 + t406 * t416) * t401 + t406 * t418;
	t403 = cos(pkin(12));
	t399 = sin(pkin(12));
	t420 = t399 * t405;
	t393 = -t398 * t420 + t403 * t402;
	t392 = -t403 * t398 - t402 * t420;
	t419 = t400 * t401;
	t409 = t392 * t404 + t399 * t419;
	t382 = t393 * t407 + t409 * t406;
	t415 = t403 * t405;
	t391 = t398 * t415 + t399 * t402;
	t390 = -t399 * t398 + t402 * t415;
	t410 = -t390 * t404 + t403 * t419;
	t425 = -t391 * t407 + t410 * t406;
	t397 = qJ(4) + qJ(5);
	t394 = sin(t397);
	t396 = qJD(4) + qJD(5);
	t422 = t394 * t396;
	t395 = cos(t397);
	t421 = t395 * t396;
	t417 = t401 * t404;
	t379 = -t391 * t406 - t410 * t407;
	t375 = t379 * qJD(3);
	t413 = -(-t390 * t400 - t403 * t417) * t396 - t375;
	t381 = -t393 * t406 + t409 * t407;
	t377 = t381 * qJD(3);
	t412 = -(-t392 * t400 + t399 * t417) * t396 - t377;
	t385 = t407 * t418 + (-t398 * t406 + t407 * t416) * t401;
	t383 = t385 * qJD(3);
	t411 = -(-t402 * t419 + t405 * t404) * t396 - t383;
	t384 = t386 * qJD(3);
	t378 = t382 * qJD(3);
	t376 = t425 * qJD(3);
	t374 = t386 * t422 + t411 * t395;
	t373 = -t386 * t421 + t411 * t394;
	t372 = t382 * t422 + t412 * t395;
	t371 = -t382 * t421 + t412 * t394;
	t370 = t413 * t395 - t422 * t425;
	t369 = t413 * t394 + t421 * t425;
	t1 = [0, 0, -t378 * t395 - t381 * t422, t371, t371, 0; 0, 0, t376 * t395 - t379 * t422, t369, t369, 0; 0, 0, -t384 * t395 - t385 * t422, t373, t373, 0; 0, 0, t378 * t394 - t381 * t421, t372, t372, 0; 0, 0, -t376 * t394 - t379 * t421, t370, t370, 0; 0, 0, t384 * t394 - t385 * t421, t374, t374, 0; 0, 0, t377, 0, 0, 0; 0, 0, t375, 0, 0, 0; 0, 0, t383, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:26
	% EndTime: 2019-10-09 21:18:26
	% DurationCPUTime: 0.68s
	% Computational Cost: add. (762->78), mult. (1826->154), div. (0->0), fcn. (2204->14), ass. (0->77)
	t591 = sin(pkin(13));
	t594 = sin(pkin(6));
	t598 = cos(pkin(6));
	t600 = sin(qJ(3));
	t593 = sin(pkin(7));
	t632 = cos(qJ(3));
	t618 = t593 * t632;
	t595 = cos(pkin(13));
	t597 = cos(pkin(7));
	t625 = t595 * t597;
	t634 = (-t591 * t600 + t632 * t625) * t594 + t598 * t618;
	t592 = sin(pkin(12));
	t596 = cos(pkin(12));
	t624 = t596 * t598;
	t581 = t591 * t624 + t592 * t595;
	t608 = -t591 * t592 + t595 * t624;
	t606 = t608 * t597;
	t627 = t593 * t594;
	t569 = t581 * t632 + (-t596 * t627 + t606) * t600;
	t590 = qJ(4) + qJ(5);
	t587 = sin(t590);
	t589 = qJD(4) + qJD(5);
	t631 = t587 * t589;
	t588 = cos(t590);
	t630 = t588 * t589;
	t628 = t592 * t598;
	t626 = t594 * t597;
	t623 = qJD(6) * t588;
	t599 = sin(qJ(6));
	t622 = qJD(6) * t599;
	t601 = cos(qJ(6));
	t621 = qJD(6) * t601;
	t614 = t594 * t618;
	t568 = t581 * t600 + t596 * t614 - t632 * t606;
	t562 = t568 * qJD(3);
	t576 = -t608 * t593 - t596 * t626;
	t617 = t576 * t589 - t562;
	t582 = -t591 * t628 + t595 * t596;
	t607 = t591 * t596 + t595 * t628;
	t605 = t607 * t597;
	t570 = t582 * t600 - t592 * t614 + t632 * t605;
	t564 = t570 * qJD(3);
	t577 = t592 * t626 + t607 * t593;
	t616 = t577 * t589 - t564;
	t572 = t634 * qJD(3);
	t580 = -t595 * t627 + t597 * t598;
	t615 = t580 * t589 + t572;
	t611 = t568 * t623 - t562;
	t610 = t570 * t623 - t564;
	t609 = -t623 * t634 + t572;
	t563 = t569 * qJD(3);
	t604 = qJD(6) * t569 - t563 * t588 + t568 * t631;
	t571 = t582 * t632 + (t592 * t627 - t605) * t600;
	t565 = t571 * qJD(3);
	t603 = qJD(6) * t571 - t565 * t588 + t570 * t631;
	t575 = t598 * t593 * t600 + (t632 * t591 + t600 * t625) * t594;
	t573 = t575 * qJD(3);
	t602 = qJD(6) * t575 - t573 * t588 - t631 * t634;
	t567 = t575 * t588 + t580 * t587;
	t566 = -t575 * t587 + t580 * t588;
	t561 = t571 * t588 + t577 * t587;
	t560 = -t571 * t587 + t577 * t588;
	t559 = t569 * t588 + t576 * t587;
	t558 = -t569 * t587 + t576 * t588;
	t557 = -t575 * t631 + t615 * t588;
	t556 = -t575 * t630 - t615 * t587;
	t555 = -t571 * t631 + t616 * t588;
	t554 = -t571 * t630 - t616 * t587;
	t553 = -t569 * t631 + t617 * t588;
	t552 = -t569 * t630 - t617 * t587;
	t551 = t556 * t601 - t566 * t622;
	t550 = -t556 * t599 - t566 * t621;
	t549 = t554 * t601 - t560 * t622;
	t548 = -t554 * t599 - t560 * t621;
	t547 = t552 * t601 - t558 * t622;
	t546 = -t552 * t599 - t558 * t621;
	t1 = [0, 0, t610 * t599 + t603 * t601, t549, t549, -t555 * t599 + t565 * t601 + (-t561 * t601 - t570 * t599) * qJD(6); 0, 0, t611 * t599 + t604 * t601, t547, t547, -t553 * t599 + t563 * t601 + (-t559 * t601 - t568 * t599) * qJD(6); 0, 0, t609 * t599 + t602 * t601, t551, t551, -t557 * t599 + t573 * t601 + (-t567 * t601 + t599 * t634) * qJD(6); 0, 0, -t603 * t599 + t610 * t601, t548, t548, -t555 * t601 - t565 * t599 + (t561 * t599 - t570 * t601) * qJD(6); 0, 0, -t604 * t599 + t611 * t601, t546, t546, -t553 * t601 - t563 * t599 + (t559 * t599 - t568 * t601) * qJD(6); 0, 0, -t602 * t599 + t609 * t601, t550, t550, -t557 * t601 - t573 * t599 + (t567 * t599 + t601 * t634) * qJD(6); 0, 0, -t565 * t587 - t570 * t630, t555, t555, 0; 0, 0, -t563 * t587 - t568 * t630, t553, t553, 0; 0, 0, -t573 * t587 + t630 * t634, t557, t557, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end