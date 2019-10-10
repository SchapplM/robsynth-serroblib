% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:31
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPPR9_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR9_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
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
	% StartTime: 2019-10-10 11:31:25
	% EndTime: 2019-10-10 11:31:25
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
	% StartTime: 2019-10-10 11:31:26
	% EndTime: 2019-10-10 11:31:26
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
	% StartTime: 2019-10-10 11:31:27
	% EndTime: 2019-10-10 11:31:27
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (160->48), mult. (518->103), div. (0->0), fcn. (536->10), ass. (0->46)
	t390 = sin(qJ(1));
	t387 = cos(pkin(6));
	t398 = qJD(2) * t387 + qJD(1);
	t389 = sin(qJ(2));
	t414 = t390 * t389;
	t404 = t387 * t414;
	t409 = qJD(2) * t389;
	t392 = cos(qJ(2));
	t393 = cos(qJ(1));
	t411 = t393 * t392;
	t371 = -qJD(1) * t404 - t390 * t409 + t398 * t411;
	t412 = t393 * t389;
	t413 = t390 * t392;
	t374 = t387 * t412 + t413;
	t388 = sin(qJ(3));
	t391 = cos(qJ(3));
	t385 = sin(pkin(6));
	t410 = qJD(1) * t385;
	t402 = t390 * t410;
	t416 = t385 * t393;
	t367 = (t374 * t388 + t391 * t416) * qJD(3) - t371 * t391 - t388 * t402;
	t417 = t385 * t390;
	t415 = t389 * t391;
	t408 = qJD(2) * t392;
	t407 = qJD(3) * t388;
	t406 = qJD(3) * t391;
	t405 = qJD(3) * t392;
	t403 = t387 * t411;
	t401 = t393 * t410;
	t400 = t385 * t408;
	t399 = t388 * t405;
	t375 = -t387 * t413 - t412;
	t368 = -qJD(1) * t403 - t393 * t408 + t398 * t414;
	t396 = -t368 * t391 + t375 * t407;
	t370 = t375 * qJD(1) - t374 * qJD(2);
	t373 = t403 - t414;
	t395 = -t370 * t391 + t373 * t407;
	t366 = -t371 * t388 - t374 * t406 + t391 * t402 + t407 * t416;
	t386 = cos(pkin(11));
	t384 = sin(pkin(11));
	t376 = -t404 + t411;
	t372 = -t388 * t400 + (-t385 * t415 - t387 * t388) * qJD(3);
	t369 = -t374 * qJD(1) + t375 * qJD(2);
	t365 = t388 * t401 + t369 * t391 + (-t376 * t388 + t391 * t417) * qJD(3);
	t364 = t369 * t388 - t391 * t401 + (t376 * t391 + t388 * t417) * qJD(3);
	t1 = [t367 * t386 + t370 * t384, t369 * t384 - t396 * t386, -t364 * t386, 0, 0, 0; t365 * t386 - t368 * t384, t371 * t384 - t395 * t386, t366 * t386, 0, 0, 0; 0, (-t386 * t399 + (t384 * t392 - t386 * t415) * qJD(2)) * t385, t372 * t386, 0, 0, 0; -t367 * t384 + t370 * t386, t369 * t386 + t396 * t384, t364 * t384, 0, 0, 0; -t365 * t384 - t368 * t386, t371 * t386 + t395 * t384, -t366 * t384, 0, 0, 0; 0, (t384 * t399 + (t384 * t415 + t386 * t392) * qJD(2)) * t385, -t372 * t384, 0, 0, 0; t366, t368 * t388 + t375 * t406, t365, 0, 0, 0; t364, t370 * t388 + t373 * t406, -t367, 0, 0, 0; 0, (-t388 * t409 + t391 * t405) * t385, t391 * t400 + (-t385 * t388 * t389 + t387 * t391) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:28
	% EndTime: 2019-10-10 11:31:28
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (157->47), mult. (518->103), div. (0->0), fcn. (536->10), ass. (0->46)
	t447 = sin(qJ(1));
	t444 = cos(pkin(6));
	t455 = qJD(2) * t444 + qJD(1);
	t446 = sin(qJ(2));
	t471 = t447 * t446;
	t461 = t444 * t471;
	t466 = qJD(2) * t446;
	t449 = cos(qJ(2));
	t450 = cos(qJ(1));
	t468 = t450 * t449;
	t428 = -qJD(1) * t461 - t447 * t466 + t455 * t468;
	t469 = t450 * t446;
	t470 = t447 * t449;
	t431 = t444 * t469 + t470;
	t445 = sin(qJ(3));
	t448 = cos(qJ(3));
	t442 = sin(pkin(6));
	t467 = qJD(1) * t442;
	t459 = t447 * t467;
	t473 = t442 * t450;
	t424 = (t431 * t445 + t448 * t473) * qJD(3) - t428 * t448 - t445 * t459;
	t474 = t442 * t447;
	t472 = t446 * t448;
	t465 = qJD(2) * t449;
	t464 = qJD(3) * t445;
	t463 = qJD(3) * t448;
	t462 = qJD(3) * t449;
	t460 = t444 * t468;
	t458 = t450 * t467;
	t457 = t442 * t465;
	t456 = t445 * t462;
	t432 = -t444 * t470 - t469;
	t425 = -qJD(1) * t460 - t450 * t465 + t455 * t471;
	t453 = t425 * t448 - t432 * t464;
	t427 = t432 * qJD(1) - t431 * qJD(2);
	t430 = t460 - t471;
	t452 = t427 * t448 - t430 * t464;
	t423 = -t428 * t445 - t431 * t463 + t448 * t459 + t464 * t473;
	t443 = cos(pkin(11));
	t441 = sin(pkin(11));
	t433 = -t461 + t468;
	t429 = -t445 * t457 + (-t442 * t472 - t444 * t445) * qJD(3);
	t426 = -t431 * qJD(1) + t432 * qJD(2);
	t422 = t445 * t458 + t426 * t448 + (-t433 * t445 + t448 * t474) * qJD(3);
	t421 = t426 * t445 - t448 * t458 + (t433 * t448 + t445 * t474) * qJD(3);
	t1 = [t424 * t443 + t427 * t441, t426 * t441 + t453 * t443, -t421 * t443, 0, 0, 0; t422 * t443 - t425 * t441, t428 * t441 + t452 * t443, t423 * t443, 0, 0, 0; 0, (-t443 * t456 + (t441 * t449 - t443 * t472) * qJD(2)) * t442, t429 * t443, 0, 0, 0; t423, t425 * t445 + t432 * t463, t422, 0, 0, 0; t421, t427 * t445 + t430 * t463, -t424, 0, 0, 0; 0, (-t445 * t466 + t448 * t462) * t442, t448 * t457 + (-t442 * t445 * t446 + t444 * t448) * qJD(3), 0, 0, 0; t424 * t441 - t427 * t443, -t426 * t443 + t453 * t441, -t421 * t441, 0, 0, 0; t422 * t441 + t425 * t443, -t428 * t443 + t452 * t441, t423 * t441, 0, 0, 0; 0, (-t441 * t456 + (-t441 * t472 - t443 * t449) * qJD(2)) * t442, t429 * t441, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:29
	% EndTime: 2019-10-10 11:31:30
	% DurationCPUTime: 1.19s
	% Computational Cost: add. (559->109), mult. (1718->211), div. (0->0), fcn. (1876->12), ass. (0->94)
	t596 = sin(qJ(1));
	t592 = cos(pkin(6));
	t609 = qJD(2) * t592 + qJD(1);
	t595 = sin(qJ(2));
	t629 = t596 * t595;
	t618 = t592 * t629;
	t623 = qJD(2) * t595;
	t599 = cos(qJ(2));
	t600 = cos(qJ(1));
	t625 = t600 * t599;
	t565 = -qJD(1) * t618 - t596 * t623 + t609 * t625;
	t626 = t600 * t595;
	t628 = t596 * t599;
	t580 = t592 * t626 + t628;
	t594 = sin(qJ(3));
	t598 = cos(qJ(3));
	t590 = sin(pkin(6));
	t624 = qJD(1) * t590;
	t615 = t596 * t624;
	t633 = t590 * t600;
	t617 = t598 * t633;
	t546 = (-qJD(3) * t580 + t615) * t594 - qJD(3) * t617 + t565 * t598;
	t581 = t592 * t628 + t626;
	t564 = t581 * qJD(1) + t580 * qJD(2);
	t589 = sin(pkin(11));
	t591 = cos(pkin(11));
	t535 = t546 * t589 - t564 * t591;
	t536 = t546 * t591 + t564 * t589;
	t572 = -t580 * t598 + t594 * t633;
	t616 = t592 * t625;
	t579 = -t616 + t629;
	t550 = t572 * t589 + t579 * t591;
	t551 = t572 * t591 - t579 * t589;
	t593 = sin(qJ(6));
	t597 = cos(qJ(6));
	t644 = -t535 * t593 - t536 * t597 + (t550 * t597 - t551 * t593) * qJD(6);
	t643 = t535 * t597 - t536 * t593 + (t550 * t593 + t551 * t597) * qJD(6);
	t545 = t572 * qJD(3) - t565 * t594 + t598 * t615;
	t636 = t589 * t598;
	t635 = t589 * t599;
	t634 = t590 * t596;
	t632 = t591 * t598;
	t631 = t591 * t599;
	t630 = t595 * t598;
	t627 = t598 * t599;
	t622 = qJD(2) * t599;
	t621 = qJD(3) * t594;
	t620 = qJD(3) * t598;
	t619 = qJD(3) * t599;
	t614 = t600 * t624;
	t613 = t590 * t623;
	t612 = t590 * t622;
	t611 = t594 * t619;
	t608 = t589 * t597 - t591 * t593;
	t607 = t589 * t593 + t591 * t597;
	t582 = -t618 + t625;
	t573 = -t582 * t594 + t598 * t634;
	t574 = t582 * t598 + t594 * t634;
	t578 = t590 * t630 + t592 * t594;
	t577 = -t590 * t595 * t594 + t592 * t598;
	t562 = -qJD(1) * t616 - t600 * t622 + t609 * t629;
	t606 = t562 * t598 + t581 * t621;
	t605 = -t564 * t598 + t579 * t621;
	t603 = qJD(6) * t608;
	t602 = qJD(6) * t607;
	t576 = (t589 * t595 + t591 * t627) * t590;
	t575 = (t589 * t627 - t591 * t595) * t590;
	t570 = -t580 * t594 - t617;
	t569 = t577 * qJD(3) + t598 * t612;
	t568 = -t578 * qJD(3) - t594 * t612;
	t567 = t578 * t591 - t590 * t635;
	t566 = t578 * t589 + t590 * t631;
	t563 = -t580 * qJD(1) - t581 * qJD(2);
	t561 = (-t591 * t611 + (-t591 * t630 + t635) * qJD(2)) * t590;
	t560 = (-t589 * t611 + (-t589 * t630 - t631) * qJD(2)) * t590;
	t559 = -t581 * t632 + t582 * t589;
	t558 = -t581 * t636 - t582 * t591;
	t557 = -t579 * t632 + t580 * t589;
	t556 = -t579 * t636 - t580 * t591;
	t555 = t569 * t591 + t589 * t613;
	t554 = t569 * t589 - t591 * t613;
	t553 = t574 * t591 + t581 * t589;
	t552 = t574 * t589 - t581 * t591;
	t544 = t573 * qJD(3) + t563 * t598 + t594 * t614;
	t543 = -t574 * qJD(3) - t563 * t594 + t598 * t614;
	t542 = t565 * t589 + t605 * t591;
	t541 = -t565 * t591 + t605 * t589;
	t540 = t563 * t589 + t606 * t591;
	t539 = -t563 * t591 + t606 * t589;
	t534 = t544 * t591 - t562 * t589;
	t533 = t544 * t589 + t562 * t591;
	t532 = t533 * t593 + t534 * t597 + (t552 * t597 - t553 * t593) * qJD(6);
	t531 = t533 * t597 - t534 * t593 + (-t552 * t593 - t553 * t597) * qJD(6);
	t1 = [t644, t539 * t593 + t540 * t597 + (t558 * t597 - t559 * t593) * qJD(6), t607 * t543 + t573 * t603, 0, 0, t531; t532, t541 * t593 + t542 * t597 + (t556 * t597 - t557 * t593) * qJD(6), t607 * t545 + t570 * t603, 0, 0, t643; 0, t560 * t593 + t561 * t597 + (t575 * t597 - t576 * t593) * qJD(6), t607 * t568 + t577 * t603, 0, 0, t554 * t597 - t555 * t593 + (-t566 * t593 - t567 * t597) * qJD(6); -t643, t539 * t597 - t540 * t593 + (-t558 * t593 - t559 * t597) * qJD(6), t608 * t543 - t573 * t602, 0, 0, -t532; t531, t541 * t597 - t542 * t593 + (-t556 * t593 - t557 * t597) * qJD(6), t608 * t545 - t570 * t602, 0, 0, t644; 0, t560 * t597 - t561 * t593 + (-t575 * t593 - t576 * t597) * qJD(6), t608 * t568 - t577 * t602, 0, 0, -t554 * t593 - t555 * t597 + (-t566 * t597 + t567 * t593) * qJD(6); -t545, -t562 * t594 + t581 * t620, -t544, 0, 0, 0; t543, t564 * t594 + t579 * t620, -t546, 0, 0, 0; 0, (t594 * t623 - t598 * t619) * t590, -t569, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end