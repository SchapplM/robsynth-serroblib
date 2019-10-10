% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:50
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRPR12_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR12_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
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
	% StartTime: 2019-10-10 12:49:59
	% EndTime: 2019-10-10 12:49:59
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
	% StartTime: 2019-10-10 12:50:00
	% EndTime: 2019-10-10 12:50:01
	% DurationCPUTime: 0.32s
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
	% StartTime: 2019-10-10 12:50:03
	% EndTime: 2019-10-10 12:50:04
	% DurationCPUTime: 1.11s
	% Computational Cost: add. (493->104), mult. (1577->199), div. (0->0), fcn. (1713->12), ass. (0->85)
	t613 = sin(qJ(2));
	t614 = sin(qJ(1));
	t617 = cos(qJ(2));
	t618 = cos(qJ(1));
	t661 = cos(pkin(6));
	t635 = t618 * t661;
	t598 = t613 * t635 + t614 * t617;
	t636 = t614 * t661;
	t619 = t618 * t613 + t617 * t636;
	t588 = t619 * qJD(1) + t598 * qJD(2);
	t629 = t613 * t636;
	t645 = t618 * t617;
	t647 = t614 * t613;
	t589 = -qJD(1) * t629 - qJD(2) * t647 + (qJD(2) * t661 + qJD(1)) * t645;
	t610 = cos(pkin(7));
	t612 = sin(qJ(3));
	t616 = cos(qJ(3));
	t608 = sin(pkin(7));
	t609 = sin(pkin(6));
	t644 = qJD(1) * t609;
	t639 = t614 * t644;
	t633 = t608 * t639;
	t597 = -t617 * t635 + t647;
	t653 = t609 * t618;
	t640 = t608 * t653;
	t627 = t597 * t610 + t640;
	t662 = t598 * t612 + t627 * t616;
	t569 = t662 * qJD(3) - (-t588 * t610 + t633) * t612 - t589 * t616;
	t658 = t598 * t616;
	t576 = t627 * t612 - t658;
	t581 = t588 * t608 + t610 * t639;
	t592 = -t597 * t608 + t610 * t653;
	t611 = sin(qJ(4));
	t615 = cos(qJ(4));
	t672 = t569 * t611 + (t576 * t615 + t592 * t611) * qJD(4) + t581 * t615;
	t671 = t569 * t615 - t581 * t611 + (-t576 * t611 + t592 * t615) * qJD(4);
	t657 = t608 * t609;
	t656 = t608 * t611;
	t655 = t608 * t615;
	t654 = t609 * t614;
	t652 = t610 * t612;
	t651 = t610 * t616;
	t650 = t612 * t613;
	t649 = t612 * t617;
	t648 = t613 * t616;
	t646 = t616 * t617;
	t643 = qJD(4) * t611;
	t642 = qJD(4) * t615;
	t641 = t613 * t657;
	t638 = t618 * t644;
	t637 = qJD(2) * t657;
	t634 = t661 * t608;
	t632 = t608 * t638;
	t631 = t613 * t637;
	t630 = t617 * t637;
	t628 = qJD(3) * t634;
	t584 = -t597 * t616 - t598 * t652;
	t620 = t629 - t645;
	t585 = -t616 * t619 + t620 * t652;
	t626 = t608 * t654 - t610 * t619;
	t625 = t610 * t646 - t650;
	t624 = -t610 * t648 - t649;
	t623 = t610 * t649 + t648;
	t622 = -t610 * t650 + t646;
	t577 = t612 * t620 + t626 * t616;
	t578 = t626 * t612 - t616 * t620;
	t567 = -t588 * t651 - t589 * t612 + t616 * t633 + (t597 * t652 + t612 * t640 - t658) * qJD(3);
	t596 = t661 * t610 - t617 * t657;
	t595 = t622 * t609;
	t594 = t608 * t619 + t610 * t654;
	t591 = t623 * t609 + t612 * t634;
	t590 = t625 * t609 + t616 * t634;
	t587 = t598 * qJD(1) + t619 * qJD(2);
	t586 = t597 * qJD(1) + t620 * qJD(2);
	t582 = (-t623 * qJD(2) + t624 * qJD(3)) * t609;
	t579 = -t586 * t608 + t610 * t638;
	t573 = t616 * t628 + (t622 * qJD(2) + t625 * qJD(3)) * t609;
	t572 = -t612 * t628 + (t624 * qJD(2) - t623 * qJD(3)) * t609;
	t571 = -t589 * t652 - t588 * t616 + (t597 * t612 - t598 * t651) * qJD(3);
	t570 = t587 * t652 + t586 * t616 + (t612 * t619 + t620 * t651) * qJD(3);
	t566 = -t587 * t616 + (t586 * t610 + t632) * t612 + t577 * qJD(3);
	t565 = t578 * qJD(3) - t586 * t651 - t587 * t612 - t616 * t632;
	t564 = t566 * t615 + t579 * t611 + (-t578 * t611 + t594 * t615) * qJD(4);
	t563 = -t566 * t611 + t579 * t615 + (-t578 * t615 - t594 * t611) * qJD(4);
	t1 = [t671, -t587 * t656 + t570 * t615 + (-t585 * t611 - t620 * t655) * qJD(4), -t565 * t615 - t577 * t643, t563, 0, 0; t564, t589 * t656 + t571 * t615 + (-t584 * t611 + t598 * t655) * qJD(4), t567 * t615 + t643 * t662, t672, 0, 0; 0, t611 * t630 + t582 * t615 + (-t595 * t611 + t615 * t641) * qJD(4), t572 * t615 - t590 * t643, t615 * t631 - t573 * t611 + (-t591 * t615 - t596 * t611) * qJD(4), 0, 0; -t672, -t587 * t655 - t570 * t611 + (-t585 * t615 + t620 * t656) * qJD(4), t565 * t611 - t577 * t642, -t564, 0, 0; t563, t589 * t655 - t571 * t611 + (-t584 * t615 - t598 * t656) * qJD(4), -t567 * t611 + t642 * t662, t671, 0, 0; 0, t615 * t630 - t582 * t611 + (-t595 * t615 - t611 * t641) * qJD(4), -t572 * t611 - t590 * t642, -t611 * t631 - t573 * t615 + (t591 * t611 - t596 * t615) * qJD(4), 0, 0; t567, t585 * qJD(3) + t586 * t612 - t587 * t651, t566, 0, 0, 0; t565, t584 * qJD(3) - t588 * t612 + t589 * t651, -t569, 0, 0, 0; 0, (t625 * qJD(2) + t622 * qJD(3)) * t609, t573, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:50:04
	% EndTime: 2019-10-10 12:50:05
	% DurationCPUTime: 1.13s
	% Computational Cost: add. (569->105), mult. (1577->199), div. (0->0), fcn. (1713->12), ass. (0->86)
	t642 = sin(qJ(2));
	t643 = sin(qJ(1));
	t645 = cos(qJ(2));
	t646 = cos(qJ(1));
	t689 = cos(pkin(6));
	t662 = t646 * t689;
	t625 = t642 * t662 + t643 * t645;
	t663 = t643 * t689;
	t647 = t646 * t642 + t645 * t663;
	t615 = t647 * qJD(1) + t625 * qJD(2);
	t657 = t642 * t663;
	t673 = t646 * t645;
	t675 = t643 * t642;
	t616 = -qJD(1) * t657 - qJD(2) * t675 + (qJD(2) * t689 + qJD(1)) * t673;
	t640 = cos(pkin(7));
	t641 = sin(qJ(3));
	t644 = cos(qJ(3));
	t638 = sin(pkin(7));
	t639 = sin(pkin(6));
	t672 = qJD(1) * t639;
	t667 = t643 * t672;
	t661 = t638 * t667;
	t624 = -t645 * t662 + t675;
	t681 = t639 * t646;
	t668 = t638 * t681;
	t655 = t624 * t640 + t668;
	t690 = t625 * t641 + t655 * t644;
	t596 = t690 * qJD(3) - (-t615 * t640 + t661) * t641 - t616 * t644;
	t686 = t625 * t644;
	t603 = t655 * t641 - t686;
	t608 = t615 * t638 + t640 * t667;
	t619 = -t624 * t638 + t640 * t681;
	t637 = qJ(4) + pkin(13);
	t635 = sin(t637);
	t636 = cos(t637);
	t700 = t596 * t635 + (t603 * t636 + t619 * t635) * qJD(4) + t608 * t636;
	t699 = t596 * t636 - t608 * t635 + (-t603 * t635 + t619 * t636) * qJD(4);
	t685 = t635 * t638;
	t684 = t636 * t638;
	t683 = t638 * t639;
	t682 = t639 * t643;
	t680 = t640 * t641;
	t679 = t640 * t644;
	t678 = t641 * t642;
	t677 = t641 * t645;
	t676 = t642 * t644;
	t674 = t644 * t645;
	t671 = qJD(4) * t635;
	t670 = qJD(4) * t636;
	t669 = t642 * t683;
	t666 = t646 * t672;
	t665 = qJD(2) * t683;
	t664 = t638 * t689;
	t660 = t638 * t666;
	t659 = t642 * t665;
	t658 = t645 * t665;
	t656 = qJD(3) * t664;
	t611 = -t624 * t644 - t625 * t680;
	t648 = t657 - t673;
	t612 = -t644 * t647 + t648 * t680;
	t654 = t638 * t682 - t640 * t647;
	t653 = t640 * t674 - t678;
	t652 = -t640 * t676 - t677;
	t651 = t640 * t677 + t676;
	t650 = -t640 * t678 + t674;
	t604 = t641 * t648 + t654 * t644;
	t605 = t654 * t641 - t644 * t648;
	t594 = -t615 * t679 - t616 * t641 + t644 * t661 + (t624 * t680 + t641 * t668 - t686) * qJD(3);
	t623 = t689 * t640 - t645 * t683;
	t622 = t650 * t639;
	t621 = t638 * t647 + t640 * t682;
	t618 = t651 * t639 + t641 * t664;
	t617 = t653 * t639 + t644 * t664;
	t614 = t625 * qJD(1) + t647 * qJD(2);
	t613 = t624 * qJD(1) + t648 * qJD(2);
	t609 = (-t651 * qJD(2) + t652 * qJD(3)) * t639;
	t606 = -t613 * t638 + t640 * t666;
	t600 = t644 * t656 + (t650 * qJD(2) + t653 * qJD(3)) * t639;
	t599 = -t641 * t656 + (t652 * qJD(2) - t651 * qJD(3)) * t639;
	t598 = -t616 * t680 - t615 * t644 + (t624 * t641 - t625 * t679) * qJD(3);
	t597 = t614 * t680 + t613 * t644 + (t641 * t647 + t648 * t679) * qJD(3);
	t593 = -t614 * t644 + (t613 * t640 + t660) * t641 + t604 * qJD(3);
	t592 = t605 * qJD(3) - t613 * t679 - t614 * t641 - t644 * t660;
	t591 = t593 * t636 + t606 * t635 + (-t605 * t635 + t621 * t636) * qJD(4);
	t590 = -t593 * t635 + t606 * t636 + (-t605 * t636 - t621 * t635) * qJD(4);
	t1 = [t699, -t614 * t685 + t597 * t636 + (-t612 * t635 - t648 * t684) * qJD(4), -t592 * t636 - t604 * t671, t590, 0, 0; t591, t616 * t685 + t598 * t636 + (-t611 * t635 + t625 * t684) * qJD(4), t594 * t636 + t671 * t690, t700, 0, 0; 0, t635 * t658 + t609 * t636 + (-t622 * t635 + t636 * t669) * qJD(4), t599 * t636 - t617 * t671, t636 * t659 - t600 * t635 + (-t618 * t636 - t623 * t635) * qJD(4), 0, 0; -t700, -t614 * t684 - t597 * t635 + (-t612 * t636 + t648 * t685) * qJD(4), t592 * t635 - t604 * t670, -t591, 0, 0; t590, t616 * t684 - t598 * t635 + (-t611 * t636 - t625 * t685) * qJD(4), -t594 * t635 + t670 * t690, t699, 0, 0; 0, t636 * t658 - t609 * t635 + (-t622 * t636 - t635 * t669) * qJD(4), -t599 * t635 - t617 * t670, -t635 * t659 - t600 * t636 + (t618 * t635 - t623 * t636) * qJD(4), 0, 0; t594, t612 * qJD(3) + t613 * t641 - t614 * t679, t593, 0, 0, 0; t592, t611 * qJD(3) - t615 * t641 + t616 * t679, -t596, 0, 0, 0; 0, (t653 * qJD(2) + t650 * qJD(3)) * t639, t600, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:50:13
	% EndTime: 2019-10-10 12:50:15
	% DurationCPUTime: 2.68s
	% Computational Cost: add. (1438->166), mult. (3901->300), div. (0->0), fcn. (4371->14), ass. (0->126)
	t886 = sin(qJ(2));
	t889 = cos(qJ(2));
	t890 = cos(qJ(1));
	t947 = cos(pkin(6));
	t919 = t890 * t947;
	t948 = sin(qJ(1));
	t865 = t886 * t919 + t948 * t889;
	t912 = t947 * t948;
	t894 = t890 * t886 + t889 * t912;
	t852 = t894 * qJD(1) + t865 * qJD(2);
	t906 = t886 * t912;
	t922 = t948 * t886;
	t932 = t890 * t889;
	t853 = -qJD(1) * t906 - qJD(2) * t922 + (qJD(2) * t947 + qJD(1)) * t932;
	t883 = cos(pkin(7));
	t885 = sin(qJ(3));
	t888 = cos(qJ(3));
	t881 = sin(pkin(7));
	t882 = sin(pkin(6));
	t923 = t882 * t948;
	t913 = qJD(1) * t923;
	t907 = t881 * t913;
	t864 = -t889 * t919 + t922;
	t939 = t882 * t890;
	t925 = t881 * t939;
	t904 = t864 * t883 + t925;
	t945 = t865 * t885;
	t949 = t904 * t888 + t945;
	t813 = t949 * qJD(3) - (-t852 * t883 + t907) * t885 - t853 * t888;
	t944 = t865 * t888;
	t836 = t904 * t885 - t944;
	t856 = -t864 * t881 + t883 * t939;
	t880 = qJ(4) + pkin(13);
	t878 = sin(t880);
	t879 = cos(t880);
	t821 = t836 * t878 - t856 * t879;
	t840 = t852 * t881 + t883 * t913;
	t805 = -t821 * qJD(4) + t813 * t879 - t840 * t878;
	t937 = t883 * t888;
	t938 = t883 * t885;
	t812 = (t864 * t938 - t944) * qJD(3) - t852 * t937 + t888 * t907 - (-qJD(3) * t925 + t853) * t885;
	t884 = sin(qJ(6));
	t887 = cos(qJ(6));
	t962 = -t805 * t884 + t812 * t887;
	t961 = t805 * t887 + t812 * t884;
	t823 = t836 * t879 + t856 * t878;
	t960 = t823 * t884;
	t959 = t823 * t887;
	t803 = t823 * qJD(4) + t813 * t878 + t840 * t879;
	t933 = t888 * t889;
	t936 = t885 * t886;
	t903 = t883 * t933 - t936;
	t895 = t906 - t932;
	t943 = t895 * t885;
	t942 = t878 * t881;
	t941 = t879 * t881;
	t940 = t881 * t882;
	t935 = t885 * t889;
	t934 = t886 * t888;
	t931 = qJD(4) * t878;
	t930 = qJD(4) * t879;
	t929 = qJD(6) * t879;
	t928 = qJD(6) * t884;
	t927 = qJD(6) * t887;
	t926 = t886 * t940;
	t921 = qJD(1) * t939;
	t920 = qJD(2) * t940;
	t918 = t947 * t881;
	t917 = t881 * t923;
	t916 = t881 * t921;
	t915 = t886 * t920;
	t914 = t889 * t920;
	t911 = qJD(3) * t918;
	t850 = t864 * qJD(1) + t895 * qJD(2);
	t851 = t865 * qJD(1) + t894 * qJD(2);
	t899 = -t883 * t894 + t917;
	t809 = -t851 * t888 + (t850 * t883 + t916) * t885 + (t899 * t888 + t943) * qJD(3);
	t837 = -t888 * t917 + t894 * t937 - t943;
	t910 = t837 * t929 + t809;
	t833 = t864 * t937 + t888 * t925 + t945;
	t909 = t833 * t929 - t813;
	t900 = -t883 * t936 + t933;
	t831 = t888 * t911 + (t900 * qJD(2) + t903 * qJD(3)) * t882;
	t854 = -t903 * t882 - t888 * t918;
	t908 = t854 * t929 + t831;
	t838 = t899 * t885 - t888 * t895;
	t858 = t881 * t894 + t883 * t923;
	t825 = t838 * t879 + t858 * t878;
	t824 = -t838 * t878 + t858 * t879;
	t901 = t883 * t935 + t934;
	t855 = t901 * t882 + t885 * t918;
	t863 = t947 * t883 - t889 * t940;
	t829 = t855 * t879 + t863 * t878;
	t828 = -t855 * t878 + t863 * t879;
	t847 = -t864 * t888 - t865 * t938;
	t826 = t847 * t879 + t865 * t942;
	t849 = -t888 * t894 + t895 * t938;
	t827 = t849 * t879 - t895 * t942;
	t846 = -t864 * t885 + t865 * t937;
	t848 = -t885 * t894 - t895 * t937;
	t902 = t883 * t934 + t935;
	t862 = t900 * t882;
	t845 = t862 * t879 + t878 * t926;
	t898 = -t850 * t881 + t883 * t921;
	t808 = t838 * qJD(3) - t850 * t937 - t851 * t885 - t888 * t916;
	t893 = qJD(6) * t838 - t808 * t879 + t837 * t931;
	t892 = -qJD(6) * t836 + t812 * t879 + t833 * t931;
	t830 = t885 * t911 + (t902 * qJD(2) + t901 * qJD(3)) * t882;
	t891 = qJD(6) * t855 - t830 * t879 + t854 * t931;
	t861 = t902 * t882;
	t842 = (-t901 * qJD(2) - t902 * qJD(3)) * t882;
	t841 = (t903 * qJD(2) + t900 * qJD(3)) * t882;
	t820 = t878 * t914 + t842 * t879 + (-t862 * t878 + t879 * t926) * qJD(4);
	t819 = -t846 * qJD(3) - t852 * t888 - t853 * t938;
	t818 = t847 * qJD(3) - t852 * t885 + t853 * t937;
	t817 = -t848 * qJD(3) + t850 * t888 + t851 * t938;
	t816 = t849 * qJD(3) + t850 * t885 - t851 * t937;
	t815 = t828 * qJD(4) + t831 * t879 + t878 * t915;
	t814 = -t829 * qJD(4) - t831 * t878 + t879 * t915;
	t807 = t853 * t942 + t819 * t879 + (-t847 * t878 + t865 * t941) * qJD(4);
	t806 = -t851 * t942 + t817 * t879 + (-t849 * t878 - t895 * t941) * qJD(4);
	t802 = t824 * qJD(4) + t809 * t879 + t898 * t878;
	t801 = t825 * qJD(4) + t809 * t878 - t898 * t879;
	t800 = t802 * t887 + t808 * t884 + (-t825 * t884 + t837 * t887) * qJD(6);
	t799 = -t802 * t884 + t808 * t887 + (-t825 * t887 - t837 * t884) * qJD(6);
	t1 = [(-t887 * t949 - t960) * qJD(6) + t961, t806 * t887 + t816 * t884 + (-t827 * t884 + t848 * t887) * qJD(6), t884 * t910 + t887 * t893, -t801 * t887 - t824 * t928, 0, t799; t800, t807 * t887 + t818 * t884 + (-t826 * t884 + t846 * t887) * qJD(6), t884 * t909 + t887 * t892, t803 * t887 - t821 * t928, 0, (-t833 * t884 + t959) * qJD(6) - t962; 0, t820 * t887 + t841 * t884 + (-t845 * t884 + t861 * t887) * qJD(6), t884 * t908 + t887 * t891, t814 * t887 - t828 * t928, 0, -t815 * t884 + t830 * t887 + (-t829 * t887 - t854 * t884) * qJD(6); (t884 * t949 - t959) * qJD(6) + t962, -t806 * t884 + t816 * t887 + (-t827 * t887 - t848 * t884) * qJD(6), -t884 * t893 + t887 * t910, t801 * t884 - t824 * t927, 0, -t800; t799, -t807 * t884 + t818 * t887 + (-t826 * t887 - t846 * t884) * qJD(6), -t884 * t892 + t887 * t909, -t803 * t884 - t821 * t927, 0, (-t833 * t887 - t960) * qJD(6) + t961; 0, -t820 * t884 + t841 * t887 + (-t845 * t887 - t861 * t884) * qJD(6), -t884 * t891 + t887 * t908, -t814 * t884 - t828 * t927, 0, -t815 * t887 - t830 * t884 + (t829 * t884 - t854 * t887) * qJD(6); t803, t827 * qJD(4) + t817 * t878 + t851 * t941, -t808 * t878 - t837 * t930, t802, 0, 0; t801, t826 * qJD(4) + t819 * t878 - t853 * t941, t812 * t878 - t833 * t930, -t805, 0, 0; 0, qJD(4) * t845 + t842 * t878 - t879 * t914, -t830 * t878 - t854 * t930, t815, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end