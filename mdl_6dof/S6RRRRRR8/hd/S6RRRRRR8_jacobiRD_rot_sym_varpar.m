% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRRR8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
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
	% StartTime: 2019-10-10 13:29:31
	% EndTime: 2019-10-10 13:29:31
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
	% StartTime: 2019-10-10 13:29:33
	% EndTime: 2019-10-10 13:29:33
	% DurationCPUTime: 0.30s
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
	% StartTime: 2019-10-10 13:29:35
	% EndTime: 2019-10-10 13:29:37
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
	% StartTime: 2019-10-10 13:29:36
	% EndTime: 2019-10-10 13:29:37
	% DurationCPUTime: 0.98s
	% Computational Cost: add. (767->94), mult. (1969->176), div. (0->0), fcn. (2143->12), ass. (0->86)
	t666 = qJ(4) + qJ(5);
	t663 = sin(t666);
	t664 = cos(t666);
	t671 = sin(qJ(2));
	t672 = sin(qJ(1));
	t674 = cos(qJ(2));
	t675 = cos(qJ(1));
	t729 = cos(pkin(6));
	t705 = t675 * t729;
	t653 = t671 * t705 + t672 * t674;
	t706 = t672 * t729;
	t678 = t675 * t671 + t674 * t706;
	t643 = t678 * qJD(1) + t653 * qJD(2);
	t665 = qJD(4) + qJD(5);
	t667 = sin(pkin(7));
	t669 = cos(pkin(7));
	t670 = sin(qJ(3));
	t714 = t672 * t671;
	t652 = -t674 * t705 + t714;
	t668 = sin(pkin(6));
	t720 = t668 * t675;
	t710 = t667 * t720;
	t688 = t652 * t669 + t710;
	t711 = qJD(1) * t668;
	t709 = t672 * t711;
	t673 = cos(qJ(3));
	t726 = t653 * t673;
	t735 = -t643 * t667 - (t688 * t670 - t726) * t665 - t669 * t709;
	t695 = t671 * t706;
	t712 = t675 * t674;
	t644 = -qJD(1) * t695 - qJD(2) * t714 + (qJD(2) * t729 + qJD(1)) * t712;
	t697 = t667 * t709;
	t730 = t653 * t670 + t688 * t673;
	t621 = -t730 * qJD(3) + (-t643 * t669 + t697) * t670 + t644 * t673;
	t736 = -t621 + (-t652 * t667 + t669 * t720) * t665;
	t616 = t663 * t736 - t664 * t735;
	t617 = t663 * t735 + t664 * t736;
	t679 = t695 - t712;
	t641 = t652 * qJD(1) + t679 * qJD(2);
	t721 = t668 * t672;
	t687 = t667 * t721 - t669 * t678;
	t677 = t687 * t670 - t673 * t679;
	t708 = t675 * t711;
	t733 = -t641 * t667 - t677 * t665 + t669 * t708;
	t725 = t663 * t665;
	t724 = t664 * t665;
	t723 = t665 * t667;
	t722 = t667 * t668;
	t719 = t669 * t670;
	t718 = t669 * t673;
	t717 = t670 * t671;
	t716 = t670 * t674;
	t715 = t671 * t673;
	t713 = t673 * t674;
	t707 = qJD(2) * t722;
	t704 = t729 * t667;
	t633 = t670 * t679 + t687 * t673;
	t642 = t653 * qJD(1) + t678 * qJD(2);
	t696 = t667 * t708;
	t619 = -t642 * t673 + (t641 * t669 + t696) * t670 + t633 * qJD(3);
	t703 = (t667 * t678 + t669 * t721) * t665 + t619;
	t683 = -t669 * t717 + t713;
	t686 = t669 * t713 - t717;
	t694 = qJD(3) * t704;
	t629 = t673 * t694 + (t683 * qJD(2) + t686 * qJD(3)) * t668;
	t700 = -(t729 * t669 - t674 * t722) * t665 - t629;
	t693 = -t679 * t723 + t642 * t719 + t641 * t673 + (t670 * t678 + t679 * t718) * qJD(3);
	t692 = t653 * t723 - t644 * t719 - t643 * t673 + (t652 * t670 - t653 * t718) * qJD(3);
	t639 = -t652 * t673 - t653 * t719;
	t691 = -t639 * t665 + t644 * t667;
	t640 = -t673 * t678 + t679 * t719;
	t690 = -t640 * t665 - t642 * t667;
	t684 = t669 * t716 + t715;
	t685 = -t669 * t715 - t716;
	t689 = t665 * t671 * t722 + (-t684 * qJD(2) + t685 * qJD(3)) * t668;
	t681 = -(t684 * t668 + t670 * t704) * t665 + t671 * t707;
	t680 = -t665 * t668 * t683 + t674 * t707;
	t620 = -t643 * t718 - t644 * t670 + t673 * t697 + (t652 * t719 + t670 * t710 - t726) * qJD(3);
	t645 = t686 * t668 + t673 * t704;
	t628 = -t670 * t694 + (t685 * qJD(2) - t684 * qJD(3)) * t668;
	t624 = -t681 * t663 + t700 * t664;
	t623 = t700 * t663 + t681 * t664;
	t618 = t677 * qJD(3) - t641 * t718 - t642 * t670 - t673 * t696;
	t615 = t733 * t663 + t703 * t664;
	t614 = -t703 * t663 + t733 * t664;
	t1 = [t617, t690 * t663 + t693 * t664, -t618 * t664 - t633 * t725, t614, t614, 0; t615, t691 * t663 + t692 * t664, t620 * t664 + t725 * t730, t616, t616, 0; 0, t680 * t663 + t689 * t664, t628 * t664 - t645 * t725, t623, t623, 0; -t616, -t693 * t663 + t690 * t664, t618 * t663 - t633 * t724, -t615, -t615, 0; t614, -t692 * t663 + t691 * t664, -t620 * t663 + t724 * t730, t617, t617, 0; 0, -t689 * t663 + t680 * t664, -t628 * t663 - t645 * t724, t624, t624, 0; t620, t640 * qJD(3) + t641 * t670 - t642 * t718, t619, 0, 0, 0; t618, t639 * qJD(3) - t643 * t670 + t644 * t718, t621, 0, 0, 0; 0, (t686 * qJD(2) + t683 * qJD(3)) * t668, t629, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:46
	% EndTime: 2019-10-10 13:29:49
	% DurationCPUTime: 2.73s
	% Computational Cost: add. (1801->164), mult. (4611->297), div. (0->0), fcn. (5174->14), ass. (0->138)
	t1031 = cos(pkin(6));
	t963 = cos(qJ(1));
	t1003 = t963 * t1031;
	t1032 = sin(qJ(1));
	t959 = sin(qJ(2));
	t962 = cos(qJ(2));
	t937 = t959 * t1003 + t1032 * t962;
	t958 = sin(qJ(3));
	t1029 = t937 * t958;
	t961 = cos(qJ(3));
	t955 = sin(pkin(6));
	t1021 = t955 * t963;
	t954 = sin(pkin(7));
	t1009 = t954 * t1021;
	t1006 = t1032 * t959;
	t936 = -t962 * t1003 + t1006;
	t956 = cos(pkin(7));
	t979 = t936 * t956 + t1009;
	t1034 = t979 * t961 + t1029;
	t992 = t1031 * t1032;
	t967 = t963 * t959 + t962 * t992;
	t924 = t967 * qJD(1) + t937 * qJD(2);
	t1014 = t963 * t962;
	t982 = t959 * t992;
	t925 = -qJD(1) * t982 - qJD(2) * t1006 + (qJD(2) * t1031 + qJD(1)) * t1014;
	t1007 = t955 * t1032;
	t993 = qJD(1) * t1007;
	t983 = t954 * t993;
	t883 = -t1034 * qJD(3) + (-t924 * t956 + t983) * t958 + t925 * t961;
	t928 = t956 * t1021 - t936 * t954;
	t952 = qJD(4) + qJD(5);
	t1000 = -t928 * t952 + t883;
	t1028 = t937 * t961;
	t908 = t979 * t958 - t1028;
	t1037 = t908 * t952 + t924 * t954 + t956 * t993;
	t953 = qJ(4) + qJ(5);
	t950 = sin(t953);
	t951 = cos(t953);
	t874 = t1000 * t951 + t1037 * t950;
	t1019 = t956 * t961;
	t1020 = t956 * t958;
	t884 = (t936 * t1020 - t1028) * qJD(3) - t924 * t1019 + t961 * t983 - (-qJD(3) * t1009 + t925) * t958;
	t957 = sin(qJ(6));
	t960 = cos(qJ(6));
	t1045 = -t874 * t957 - t884 * t960;
	t1044 = -t874 * t960 + t884 * t957;
	t895 = t908 * t951 + t928 * t950;
	t1043 = t895 * t957;
	t1042 = t895 * t960;
	t873 = -t1000 * t950 + t1037 * t951;
	t1015 = t961 * t962;
	t1018 = t958 * t959;
	t978 = t956 * t1015 - t1018;
	t968 = t982 - t1014;
	t1027 = t968 * t958;
	t1026 = t950 * t952;
	t1025 = t950 * t954;
	t1024 = t951 * t952;
	t1023 = t952 * t954;
	t1022 = t954 * t955;
	t1017 = t958 * t962;
	t1016 = t959 * t961;
	t1013 = qJD(6) * t951;
	t1012 = qJD(6) * t957;
	t1011 = qJD(6) * t960;
	t1010 = t959 * t1022;
	t1005 = qJD(1) * t1021;
	t1004 = qJD(2) * t1022;
	t1002 = t1031 * t954;
	t922 = t936 * qJD(1) + t968 * qJD(2);
	t923 = t937 * qJD(1) + t967 * qJD(2);
	t995 = t954 * t1007;
	t974 = -t956 * t967 + t995;
	t994 = t954 * t1005;
	t881 = -t923 * t961 + (t922 * t956 + t994) * t958 + (t974 * t961 + t1027) * qJD(3);
	t930 = t956 * t1007 + t954 * t967;
	t1001 = t930 * t952 + t881;
	t975 = -t956 * t1018 + t1015;
	t991 = qJD(3) * t1002;
	t901 = t961 * t991 + (t975 * qJD(2) + t978 * qJD(3)) * t955;
	t935 = -t962 * t1022 + t1031 * t956;
	t998 = t935 * t952 + t901;
	t919 = -t1019 * t968 - t958 * t967;
	t990 = -t919 * qJD(3) + t923 * t1020 - t1023 * t968 + t922 * t961;
	t917 = t937 * t1019 - t936 * t958;
	t989 = -t917 * qJD(3) - t925 * t1020 + t937 * t1023 - t924 * t961;
	t909 = t1019 * t967 - t961 * t995 - t1027;
	t988 = t909 * t1013 + t881;
	t905 = t961 * t1009 + t936 * t1019 + t1029;
	t987 = t905 * t1013 + t883;
	t926 = -t961 * t1002 - t978 * t955;
	t986 = t926 * t1013 + t901;
	t918 = -t937 * t1020 - t936 * t961;
	t985 = t918 * t952 - t925 * t954;
	t920 = t1020 * t968 - t961 * t967;
	t984 = t920 * t952 + t923 * t954;
	t976 = t956 * t1017 + t1016;
	t977 = t956 * t1016 + t1017;
	t981 = t952 * t1010 + (-t976 * qJD(2) - t977 * qJD(3)) * t955;
	t973 = t956 * t1005 - t922 * t954;
	t927 = t958 * t1002 + t976 * t955;
	t972 = t959 * t1004 - t927 * t952;
	t934 = t975 * t955;
	t971 = t962 * t1004 - t934 * t952;
	t910 = t974 * t958 - t961 * t968;
	t880 = t910 * qJD(3) - t922 * t1019 - t923 * t958 - t961 * t994;
	t966 = qJD(6) * t910 + t909 * t1026 - t880 * t951;
	t965 = -qJD(6) * t908 + t905 * t1026 + t884 * t951;
	t900 = t958 * t991 + (t977 * qJD(2) + t976 * qJD(3)) * t955;
	t964 = qJD(6) * t927 + t926 * t1026 - t900 * t951;
	t933 = t977 * t955;
	t921 = t950 * t1010 + t934 * t951;
	t913 = (t978 * qJD(2) + t975 * qJD(3)) * t955;
	t903 = t927 * t951 + t935 * t950;
	t902 = -t927 * t950 + t935 * t951;
	t899 = -t1025 * t968 + t920 * t951;
	t898 = t937 * t1025 + t918 * t951;
	t897 = t910 * t951 + t930 * t950;
	t896 = -t910 * t950 + t930 * t951;
	t893 = t908 * t950 - t928 * t951;
	t892 = t971 * t950 + t981 * t951;
	t890 = t918 * qJD(3) + t925 * t1019 - t924 * t958;
	t888 = t920 * qJD(3) - t923 * t1019 + t922 * t958;
	t887 = t972 * t950 + t998 * t951;
	t886 = -t998 * t950 + t972 * t951;
	t879 = -t902 * t1012 + t886 * t960;
	t878 = -t902 * t1011 - t886 * t957;
	t877 = -t985 * t950 + t989 * t951;
	t876 = -t984 * t950 + t990 * t951;
	t872 = t1001 * t951 + (-t910 * t952 + t973) * t950;
	t871 = t1001 * t950 + t910 * t1024 - t973 * t951;
	t870 = -t893 * t1012 + t873 * t960;
	t869 = -t893 * t1011 - t873 * t957;
	t868 = -t896 * t1012 - t871 * t960;
	t867 = -t896 * t1011 + t871 * t957;
	t866 = t872 * t960 + t880 * t957 + (-t897 * t957 + t909 * t960) * qJD(6);
	t865 = -t872 * t957 + t880 * t960 + (-t897 * t960 - t909 * t957) * qJD(6);
	t1 = [(-t1034 * t960 - t1043) * qJD(6) + t1044, t876 * t960 + t888 * t957 + (-t899 * t957 + t919 * t960) * qJD(6), t957 * t988 + t960 * t966, t868, t868, t865; t866, t877 * t960 + t890 * t957 + (-t898 * t957 + t917 * t960) * qJD(6), t957 * t987 + t960 * t965, t870, t870, (-t905 * t957 + t1042) * qJD(6) + t1045; 0, t892 * t960 + t913 * t957 + (-t921 * t957 + t933 * t960) * qJD(6), t957 * t986 + t960 * t964, t879, t879, -t887 * t957 + t900 * t960 + (-t903 * t960 - t926 * t957) * qJD(6); (t1034 * t957 - t1042) * qJD(6) - t1045, -t876 * t957 + t888 * t960 + (-t899 * t960 - t919 * t957) * qJD(6), -t957 * t966 + t960 * t988, t867, t867, -t866; t865, -t877 * t957 + t890 * t960 + (-t898 * t960 - t917 * t957) * qJD(6), -t957 * t965 + t960 * t987, t869, t869, (-t905 * t960 - t1043) * qJD(6) + t1044; 0, -t892 * t957 + t913 * t960 + (-t921 * t960 - t933 * t957) * qJD(6), -t957 * t964 + t960 * t986, t878, t878, -t887 * t960 - t900 * t957 + (t903 * t957 - t926 * t960) * qJD(6); t873, t950 * t990 + t951 * t984, -t909 * t1024 - t880 * t950, t872, t872, 0; t871, t950 * t989 + t951 * t985, -t905 * t1024 + t884 * t950, t874, t874, 0; 0, t950 * t981 - t951 * t971, -t926 * t1024 - t900 * t950, t887, t887, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end