% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:14
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPRR13_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR13_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
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
	% StartTime: 2019-10-10 12:14:18
	% EndTime: 2019-10-10 12:14:18
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
	% StartTime: 2019-10-10 12:14:20
	% EndTime: 2019-10-10 12:14:20
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (178->54), mult. (587->106), div. (0->0), fcn. (611->10), ass. (0->49)
	t364 = cos(pkin(6));
	t366 = sin(qJ(2));
	t367 = sin(qJ(1));
	t392 = t367 * t366;
	t385 = t364 * t392;
	t369 = cos(qJ(2));
	t370 = cos(qJ(1));
	t388 = t370 * t369;
	t351 = -qJD(1) * t385 - qJD(2) * t392 + (qJD(2) * t364 + qJD(1)) * t388;
	t389 = t370 * t366;
	t391 = t367 * t369;
	t353 = t364 * t389 + t391;
	t361 = sin(pkin(7));
	t365 = sin(qJ(3));
	t368 = cos(qJ(3));
	t375 = t364 * t391 + t389;
	t350 = t375 * qJD(1) + t353 * qJD(2);
	t363 = cos(pkin(7));
	t362 = sin(pkin(6));
	t387 = qJD(1) * t362;
	t384 = t367 * t387;
	t372 = -t350 * t363 + t361 * t384;
	t398 = t362 * t370;
	t352 = -t364 * t388 + t392;
	t401 = t352 * t363;
	t405 = ((t361 * t398 + t401) * t365 - t353 * t368) * qJD(3) - t351 * t365 + t372 * t368;
	t399 = t361 * t362;
	t397 = t363 * t365;
	t396 = t363 * t368;
	t395 = t365 * t366;
	t394 = t365 * t369;
	t393 = t366 * t368;
	t390 = t368 * t369;
	t386 = qJD(3) * t361;
	t383 = t370 * t387;
	t382 = t368 * t386;
	t380 = -t363 * t375 + t367 * t399;
	t379 = -t363 * t390 + t395;
	t378 = -t363 * t393 - t394;
	t377 = -t363 * t394 - t393;
	t376 = t363 * t395 - t390;
	t374 = t385 - t388;
	t348 = t352 * qJD(1) + t374 * qJD(2);
	t373 = t348 * t363 + t361 * t383;
	t371 = t382 * t398 + (qJD(3) * t401 - t351) * t368 + (qJD(3) * t353 - t372) * t365;
	t349 = t353 * qJD(1) + t375 * qJD(2);
	t347 = -t349 * t368 + t373 * t365 + (t365 * t374 + t380 * t368) * qJD(3);
	t346 = t349 * t365 + t373 * t368 + (-t380 * t365 + t368 * t374) * qJD(3);
	t1 = [t371, t349 * t397 + t348 * t368 + (t365 * t375 + t374 * t396) * qJD(3), t346, 0, 0, 0; t347, -t351 * t397 - t350 * t368 + (t352 * t365 - t353 * t396) * qJD(3), t405, 0, 0, 0; 0, (t377 * qJD(2) + t378 * qJD(3)) * t362, -t364 * t365 * t386 + (t378 * qJD(2) + t377 * qJD(3)) * t362, 0, 0, 0; -t405, t349 * t396 - t348 * t365 + (t368 * t375 - t374 * t397) * qJD(3), -t347, 0, 0, 0; t346, -t351 * t396 + t350 * t365 + (t352 * t368 + t353 * t397) * qJD(3), t371, 0, 0, 0; 0, (t379 * qJD(2) + t376 * qJD(3)) * t362, -t364 * t382 + (t376 * qJD(2) + t379 * qJD(3)) * t362, 0, 0, 0; -t350 * t361 - t363 * t384, -t349 * t361, 0, 0, 0, 0; -t348 * t361 + t363 * t383, t351 * t361, 0, 0, 0, 0; 0, qJD(2) * t369 * t399, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:22
	% EndTime: 2019-10-10 12:14:22
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (284->67), mult. (933->134), div. (0->0), fcn. (969->12), ass. (0->59)
	t537 = cos(pkin(6));
	t539 = sin(qJ(2));
	t543 = cos(qJ(1));
	t563 = t543 * t539;
	t540 = sin(qJ(1));
	t542 = cos(qJ(2));
	t565 = t540 * t542;
	t522 = t537 * t563 + t565;
	t546 = t537 * t565 + t563;
	t519 = t546 * qJD(1) + t522 * qJD(2);
	t566 = t540 * t539;
	t559 = t537 * t566;
	t562 = t543 * t542;
	t520 = -qJD(1) * t559 - qJD(2) * t566 + (qJD(2) * t537 + qJD(1)) * t562;
	t521 = -t537 * t562 + t566;
	t533 = sin(pkin(7));
	t536 = cos(pkin(7));
	t538 = sin(qJ(3));
	t541 = cos(qJ(3));
	t534 = sin(pkin(6));
	t561 = qJD(1) * t534;
	t558 = t540 * t561;
	t555 = t533 * t558;
	t572 = t534 * t543;
	t509 = ((t521 * t536 + t533 * t572) * t541 + t522 * t538) * qJD(3) - (-t519 * t536 + t555) * t538 - t520 * t541;
	t532 = sin(pkin(13));
	t575 = t532 * t533;
	t574 = t533 * t534;
	t535 = cos(pkin(13));
	t573 = t533 * t535;
	t571 = t536 * t538;
	t570 = t536 * t541;
	t569 = t538 * t539;
	t568 = t538 * t542;
	t567 = t539 * t541;
	t564 = t541 * t542;
	t560 = qJD(3) * t533;
	t557 = t543 * t561;
	t556 = t538 * t560;
	t554 = t533 * t557;
	t553 = qJD(2) * t542 * t574;
	t551 = -t536 * t546 + t540 * t574;
	t550 = t536 * t564 - t569;
	t549 = -t536 * t567 - t568;
	t548 = -t536 * t568 - t567;
	t547 = -t536 * t569 + t564;
	t545 = t559 - t562;
	t508 = -t519 * t570 - t520 * t538 + t556 * t572 + t541 * t555 + (t521 * t571 - t522 * t541) * qJD(3);
	t518 = t522 * qJD(1) + t546 * qJD(2);
	t517 = t521 * qJD(1) + t545 * qJD(2);
	t515 = (t548 * qJD(2) + t549 * qJD(3)) * t534;
	t514 = -t519 * t533 - t536 * t558;
	t513 = -t517 * t533 + t536 * t557;
	t512 = -t537 * t556 + (t549 * qJD(2) + t548 * qJD(3)) * t534;
	t511 = -t520 * t571 - t519 * t541 + (t521 * t538 - t522 * t570) * qJD(3);
	t510 = t518 * t571 + t517 * t541 + (t538 * t546 + t545 * t570) * qJD(3);
	t507 = -t518 * t541 + (t517 * t536 + t554) * t538 + (t538 * t545 + t551 * t541) * qJD(3);
	t506 = -t518 * t538 - t517 * t570 - t541 * t554 + (t551 * t538 - t541 * t545) * qJD(3);
	t1 = [t509 * t535 + t514 * t532, t510 * t535 - t518 * t575, -t506 * t535, 0, 0, 0; t507 * t535 + t513 * t532, t511 * t535 + t520 * t575, t508 * t535, 0, 0, 0; 0, t515 * t535 + t532 * t553, t512 * t535, 0, 0, 0; -t509 * t532 + t514 * t535, -t510 * t532 - t518 * t573, t506 * t532, 0, 0, 0; -t507 * t532 + t513 * t535, -t511 * t532 + t520 * t573, -t508 * t532, 0, 0, 0; 0, -t515 * t532 + t535 * t553, -t512 * t532, 0, 0, 0; t508, -t518 * t570 + t517 * t538 + (-t541 * t546 + t545 * t571) * qJD(3), t507, 0, 0, 0; t506, t520 * t570 - t519 * t538 + (-t521 * t541 - t522 * t571) * qJD(3), -t509, 0, 0, 0; 0, (t550 * qJD(2) + t547 * qJD(3)) * t534, t537 * t541 * t560 + (t547 * qJD(2) + t550 * qJD(3)) * t534, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:23
	% EndTime: 2019-10-10 12:14:24
	% DurationCPUTime: 1.12s
	% Computational Cost: add. (569->105), mult. (1577->199), div. (0->0), fcn. (1713->12), ass. (0->86)
	t640 = sin(qJ(2));
	t641 = sin(qJ(1));
	t643 = cos(qJ(2));
	t644 = cos(qJ(1));
	t687 = cos(pkin(6));
	t661 = t644 * t687;
	t623 = t640 * t661 + t641 * t643;
	t662 = t641 * t687;
	t645 = t644 * t640 + t643 * t662;
	t613 = t645 * qJD(1) + t623 * qJD(2);
	t655 = t640 * t662;
	t671 = t644 * t643;
	t673 = t641 * t640;
	t614 = -qJD(1) * t655 - qJD(2) * t673 + (qJD(2) * t687 + qJD(1)) * t671;
	t638 = cos(pkin(7));
	t639 = sin(qJ(3));
	t642 = cos(qJ(3));
	t636 = sin(pkin(7));
	t637 = sin(pkin(6));
	t670 = qJD(1) * t637;
	t665 = t641 * t670;
	t659 = t636 * t665;
	t622 = -t643 * t661 + t673;
	t679 = t637 * t644;
	t666 = t636 * t679;
	t653 = t622 * t638 + t666;
	t688 = t623 * t639 + t653 * t642;
	t594 = t688 * qJD(3) - (-t613 * t638 + t659) * t639 - t614 * t642;
	t684 = t623 * t642;
	t601 = t653 * t639 - t684;
	t606 = t613 * t636 + t638 * t665;
	t617 = -t622 * t636 + t638 * t679;
	t635 = pkin(13) + qJ(5);
	t633 = sin(t635);
	t634 = cos(t635);
	t698 = t594 * t633 + (t601 * t634 + t617 * t633) * qJD(5) + t606 * t634;
	t697 = t594 * t634 - t606 * t633 + (-t601 * t633 + t617 * t634) * qJD(5);
	t683 = t633 * t636;
	t682 = t634 * t636;
	t681 = t636 * t637;
	t680 = t637 * t641;
	t678 = t638 * t639;
	t677 = t638 * t642;
	t676 = t639 * t640;
	t675 = t639 * t643;
	t674 = t640 * t642;
	t672 = t642 * t643;
	t669 = qJD(5) * t633;
	t668 = qJD(5) * t634;
	t667 = t640 * t681;
	t664 = t644 * t670;
	t663 = qJD(2) * t681;
	t660 = t687 * t636;
	t658 = t636 * t664;
	t657 = t640 * t663;
	t656 = t643 * t663;
	t654 = qJD(3) * t660;
	t609 = -t622 * t642 - t623 * t678;
	t646 = t655 - t671;
	t610 = -t642 * t645 + t646 * t678;
	t652 = t636 * t680 - t638 * t645;
	t651 = t638 * t672 - t676;
	t650 = -t638 * t674 - t675;
	t649 = t638 * t675 + t674;
	t648 = -t638 * t676 + t672;
	t602 = t639 * t646 + t652 * t642;
	t603 = t652 * t639 - t642 * t646;
	t592 = -t613 * t677 - t614 * t639 + t642 * t659 + (t622 * t678 + t639 * t666 - t684) * qJD(3);
	t621 = t687 * t638 - t643 * t681;
	t620 = t648 * t637;
	t619 = t636 * t645 + t638 * t680;
	t616 = t649 * t637 + t639 * t660;
	t615 = t651 * t637 + t642 * t660;
	t612 = t623 * qJD(1) + t645 * qJD(2);
	t611 = t622 * qJD(1) + t646 * qJD(2);
	t607 = (-t649 * qJD(2) + t650 * qJD(3)) * t637;
	t604 = -t611 * t636 + t638 * t664;
	t598 = t642 * t654 + (t648 * qJD(2) + t651 * qJD(3)) * t637;
	t597 = -t639 * t654 + (t650 * qJD(2) - t649 * qJD(3)) * t637;
	t596 = -t614 * t678 - t613 * t642 + (t622 * t639 - t623 * t677) * qJD(3);
	t595 = t612 * t678 + t611 * t642 + (t639 * t645 + t646 * t677) * qJD(3);
	t591 = -t612 * t642 + (t611 * t638 + t658) * t639 + t602 * qJD(3);
	t590 = t603 * qJD(3) - t611 * t677 - t612 * t639 - t642 * t658;
	t589 = t591 * t634 + t604 * t633 + (-t603 * t633 + t619 * t634) * qJD(5);
	t588 = -t591 * t633 + t604 * t634 + (-t603 * t634 - t619 * t633) * qJD(5);
	t1 = [t697, -t612 * t683 + t595 * t634 + (-t610 * t633 - t646 * t682) * qJD(5), -t590 * t634 - t602 * t669, 0, t588, 0; t589, t614 * t683 + t596 * t634 + (-t609 * t633 + t623 * t682) * qJD(5), t592 * t634 + t669 * t688, 0, t698, 0; 0, t633 * t656 + t607 * t634 + (-t620 * t633 + t634 * t667) * qJD(5), t597 * t634 - t615 * t669, 0, t634 * t657 - t598 * t633 + (-t616 * t634 - t621 * t633) * qJD(5), 0; -t698, -t612 * t682 - t595 * t633 + (-t610 * t634 + t646 * t683) * qJD(5), t590 * t633 - t602 * t668, 0, -t589, 0; t588, t614 * t682 - t596 * t633 + (-t609 * t634 - t623 * t683) * qJD(5), -t592 * t633 + t668 * t688, 0, t697, 0; 0, t634 * t656 - t607 * t633 + (-t620 * t634 - t633 * t667) * qJD(5), -t597 * t633 - t615 * t668, 0, -t633 * t657 - t598 * t634 + (t616 * t633 - t621 * t634) * qJD(5), 0; t592, t610 * qJD(3) + t611 * t639 - t612 * t677, t591, 0, 0, 0; t590, t609 * qJD(3) - t613 * t639 + t614 * t677, -t594, 0, 0, 0; 0, (t651 * qJD(2) + t648 * qJD(3)) * t637, t598, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:32
	% EndTime: 2019-10-10 12:14:34
	% DurationCPUTime: 2.64s
	% Computational Cost: add. (1438->166), mult. (3901->300), div. (0->0), fcn. (4371->14), ass. (0->126)
	t878 = sin(qJ(2));
	t881 = cos(qJ(2));
	t882 = cos(qJ(1));
	t939 = cos(pkin(6));
	t911 = t882 * t939;
	t940 = sin(qJ(1));
	t857 = t878 * t911 + t940 * t881;
	t904 = t939 * t940;
	t886 = t882 * t878 + t881 * t904;
	t844 = t886 * qJD(1) + t857 * qJD(2);
	t898 = t878 * t904;
	t914 = t940 * t878;
	t924 = t882 * t881;
	t845 = -qJD(1) * t898 - qJD(2) * t914 + (qJD(2) * t939 + qJD(1)) * t924;
	t875 = cos(pkin(7));
	t877 = sin(qJ(3));
	t880 = cos(qJ(3));
	t873 = sin(pkin(7));
	t874 = sin(pkin(6));
	t915 = t874 * t940;
	t905 = qJD(1) * t915;
	t899 = t873 * t905;
	t856 = -t881 * t911 + t914;
	t931 = t874 * t882;
	t917 = t873 * t931;
	t896 = t856 * t875 + t917;
	t937 = t857 * t877;
	t941 = t896 * t880 + t937;
	t805 = t941 * qJD(3) - (-t844 * t875 + t899) * t877 - t845 * t880;
	t936 = t857 * t880;
	t828 = t896 * t877 - t936;
	t848 = -t856 * t873 + t875 * t931;
	t872 = pkin(13) + qJ(5);
	t870 = sin(t872);
	t871 = cos(t872);
	t813 = t828 * t870 - t848 * t871;
	t832 = t844 * t873 + t875 * t905;
	t797 = -t813 * qJD(5) + t805 * t871 - t832 * t870;
	t929 = t875 * t880;
	t930 = t875 * t877;
	t804 = (t856 * t930 - t936) * qJD(3) - t844 * t929 + t880 * t899 - (-qJD(3) * t917 + t845) * t877;
	t876 = sin(qJ(6));
	t879 = cos(qJ(6));
	t954 = -t797 * t876 + t804 * t879;
	t953 = t797 * t879 + t804 * t876;
	t815 = t828 * t871 + t848 * t870;
	t952 = t815 * t876;
	t951 = t815 * t879;
	t795 = t815 * qJD(5) + t805 * t870 + t832 * t871;
	t925 = t880 * t881;
	t928 = t877 * t878;
	t895 = t875 * t925 - t928;
	t887 = t898 - t924;
	t935 = t887 * t877;
	t934 = t870 * t873;
	t933 = t871 * t873;
	t932 = t873 * t874;
	t927 = t877 * t881;
	t926 = t878 * t880;
	t923 = qJD(5) * t870;
	t922 = qJD(5) * t871;
	t921 = qJD(6) * t871;
	t920 = qJD(6) * t876;
	t919 = qJD(6) * t879;
	t918 = t878 * t932;
	t913 = qJD(1) * t931;
	t912 = qJD(2) * t932;
	t910 = t939 * t873;
	t909 = t873 * t915;
	t908 = t873 * t913;
	t907 = t878 * t912;
	t906 = t881 * t912;
	t903 = qJD(3) * t910;
	t842 = t856 * qJD(1) + t887 * qJD(2);
	t843 = t857 * qJD(1) + t886 * qJD(2);
	t891 = -t875 * t886 + t909;
	t801 = -t843 * t880 + (t842 * t875 + t908) * t877 + (t891 * t880 + t935) * qJD(3);
	t829 = -t880 * t909 + t886 * t929 - t935;
	t902 = t829 * t921 + t801;
	t825 = t856 * t929 + t880 * t917 + t937;
	t901 = t825 * t921 - t805;
	t892 = -t875 * t928 + t925;
	t823 = t880 * t903 + (t892 * qJD(2) + t895 * qJD(3)) * t874;
	t846 = -t895 * t874 - t880 * t910;
	t900 = t846 * t921 + t823;
	t830 = t891 * t877 - t880 * t887;
	t850 = t873 * t886 + t875 * t915;
	t817 = t830 * t871 + t850 * t870;
	t816 = -t830 * t870 + t850 * t871;
	t893 = t875 * t927 + t926;
	t847 = t893 * t874 + t877 * t910;
	t855 = t939 * t875 - t881 * t932;
	t821 = t847 * t871 + t855 * t870;
	t820 = -t847 * t870 + t855 * t871;
	t839 = -t856 * t880 - t857 * t930;
	t818 = t839 * t871 + t857 * t934;
	t841 = -t880 * t886 + t887 * t930;
	t819 = t841 * t871 - t887 * t934;
	t838 = -t856 * t877 + t857 * t929;
	t840 = -t877 * t886 - t887 * t929;
	t894 = t875 * t926 + t927;
	t854 = t892 * t874;
	t837 = t854 * t871 + t870 * t918;
	t890 = -t842 * t873 + t875 * t913;
	t800 = qJD(3) * t830 - t842 * t929 - t843 * t877 - t880 * t908;
	t885 = qJD(6) * t830 - t800 * t871 + t829 * t923;
	t884 = -qJD(6) * t828 + t804 * t871 + t825 * t923;
	t822 = t877 * t903 + (t894 * qJD(2) + t893 * qJD(3)) * t874;
	t883 = qJD(6) * t847 - t822 * t871 + t846 * t923;
	t853 = t894 * t874;
	t834 = (-t893 * qJD(2) - t894 * qJD(3)) * t874;
	t833 = (t895 * qJD(2) + t892 * qJD(3)) * t874;
	t812 = t870 * t906 + t834 * t871 + (-t854 * t870 + t871 * t918) * qJD(5);
	t811 = -t838 * qJD(3) - t844 * t880 - t845 * t930;
	t810 = t839 * qJD(3) - t844 * t877 + t845 * t929;
	t809 = -t840 * qJD(3) + t842 * t880 + t843 * t930;
	t808 = t841 * qJD(3) + t842 * t877 - t843 * t929;
	t807 = t820 * qJD(5) + t823 * t871 + t870 * t907;
	t806 = -t821 * qJD(5) - t823 * t870 + t871 * t907;
	t799 = t845 * t934 + t811 * t871 + (-t839 * t870 + t857 * t933) * qJD(5);
	t798 = -t843 * t934 + t809 * t871 + (-t841 * t870 - t887 * t933) * qJD(5);
	t794 = t816 * qJD(5) + t801 * t871 + t890 * t870;
	t793 = t817 * qJD(5) + t801 * t870 - t890 * t871;
	t792 = t794 * t879 + t800 * t876 + (-t817 * t876 + t829 * t879) * qJD(6);
	t791 = -t794 * t876 + t800 * t879 + (-t817 * t879 - t829 * t876) * qJD(6);
	t1 = [(-t879 * t941 - t952) * qJD(6) + t953, t798 * t879 + t808 * t876 + (-t819 * t876 + t840 * t879) * qJD(6), t902 * t876 + t885 * t879, 0, -t793 * t879 - t816 * t920, t791; t792, t799 * t879 + t810 * t876 + (-t818 * t876 + t838 * t879) * qJD(6), t901 * t876 + t884 * t879, 0, t795 * t879 - t813 * t920, (-t825 * t876 + t951) * qJD(6) - t954; 0, t812 * t879 + t833 * t876 + (-t837 * t876 + t853 * t879) * qJD(6), t900 * t876 + t883 * t879, 0, t806 * t879 - t820 * t920, -t807 * t876 + t822 * t879 + (-t821 * t879 - t846 * t876) * qJD(6); (t876 * t941 - t951) * qJD(6) + t954, -t798 * t876 + t808 * t879 + (-t819 * t879 - t840 * t876) * qJD(6), -t885 * t876 + t902 * t879, 0, t793 * t876 - t816 * t919, -t792; t791, -t799 * t876 + t810 * t879 + (-t818 * t879 - t838 * t876) * qJD(6), -t884 * t876 + t901 * t879, 0, -t795 * t876 - t813 * t919, (-t825 * t879 - t952) * qJD(6) + t953; 0, -t812 * t876 + t833 * t879 + (-t837 * t879 - t853 * t876) * qJD(6), -t883 * t876 + t900 * t879, 0, -t806 * t876 - t820 * t919, -t807 * t879 - t822 * t876 + (t821 * t876 - t846 * t879) * qJD(6); t795, t819 * qJD(5) + t809 * t870 + t843 * t933, -t800 * t870 - t829 * t922, 0, t794, 0; t793, t818 * qJD(5) + t811 * t870 - t845 * t933, t804 * t870 - t825 * t922, 0, -t797, 0; 0, t837 * qJD(5) + t834 * t870 - t871 * t906, -t822 * t870 - t846 * t922, 0, t807, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end