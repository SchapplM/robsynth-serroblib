% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR14
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
% Datum: 2019-10-10 12:54
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRPR14_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR14_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
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
	% StartTime: 2019-10-10 12:54:13
	% EndTime: 2019-10-10 12:54:13
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
	% StartTime: 2019-10-10 12:54:15
	% EndTime: 2019-10-10 12:54:15
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
	% StartTime: 2019-10-10 12:54:18
	% EndTime: 2019-10-10 12:54:19
	% DurationCPUTime: 1.12s
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
	% StartTime: 2019-10-10 12:54:24
	% EndTime: 2019-10-10 12:54:26
	% DurationCPUTime: 1.32s
	% Computational Cost: add. (799->121), mult. (2543->237), div. (0->0), fcn. (2753->14), ass. (0->97)
	t758 = sin(qJ(2));
	t759 = sin(qJ(1));
	t762 = cos(qJ(2));
	t763 = cos(qJ(1));
	t811 = cos(pkin(6));
	t785 = t763 * t811;
	t741 = t758 * t785 + t759 * t762;
	t786 = t759 * t811;
	t765 = t763 * t758 + t762 * t786;
	t731 = t765 * qJD(1) + t741 * qJD(2);
	t779 = t758 * t786;
	t795 = t763 * t762;
	t797 = t759 * t758;
	t732 = -qJD(1) * t779 - qJD(2) * t797 + (qJD(2) * t811 + qJD(1)) * t795;
	t755 = cos(pkin(7));
	t757 = sin(qJ(3));
	t761 = cos(qJ(3));
	t752 = sin(pkin(7));
	t753 = sin(pkin(6));
	t794 = qJD(1) * t753;
	t789 = t759 * t794;
	t783 = t752 * t789;
	t740 = -t762 * t785 + t797;
	t803 = t753 * t763;
	t790 = t752 * t803;
	t777 = t740 * t755 + t790;
	t812 = t741 * t757 + t777 * t761;
	t707 = t812 * qJD(3) - (-t731 * t755 + t783) * t757 - t732 * t761;
	t808 = t741 * t761;
	t719 = t777 * t757 - t808;
	t723 = t731 * t752 + t755 * t789;
	t735 = -t740 * t752 + t755 * t803;
	t756 = sin(qJ(4));
	t760 = cos(qJ(4));
	t698 = t707 * t756 + t723 * t760 + (t719 * t760 + t735 * t756) * qJD(4);
	t699 = t707 * t760 - (t719 * t756 - t735 * t760) * qJD(4) - t723 * t756;
	t807 = t752 * t753;
	t806 = t752 * t756;
	t805 = t752 * t760;
	t804 = t753 * t759;
	t802 = t755 * t757;
	t801 = t755 * t761;
	t800 = t757 * t758;
	t799 = t757 * t762;
	t798 = t758 * t761;
	t796 = t761 * t762;
	t793 = qJD(4) * t756;
	t792 = qJD(4) * t760;
	t791 = t758 * t807;
	t788 = t763 * t794;
	t787 = qJD(2) * t807;
	t784 = t811 * t752;
	t782 = t752 * t788;
	t781 = t758 * t787;
	t780 = t762 * t787;
	t778 = qJD(3) * t784;
	t727 = -t740 * t761 - t741 * t802;
	t766 = t779 - t795;
	t728 = -t761 * t765 + t766 * t802;
	t776 = t752 * t804 - t755 * t765;
	t775 = t755 * t796 - t800;
	t774 = -t755 * t798 - t799;
	t773 = t755 * t799 + t798;
	t772 = -t755 * t800 + t796;
	t721 = t776 * t757 - t761 * t766;
	t729 = t740 * qJD(1) + t766 * qJD(2);
	t730 = t741 * qJD(1) + t765 * qJD(2);
	t702 = t721 * qJD(3) - t729 * t801 - t730 * t757 - t761 * t782;
	t720 = t757 * t766 + t776 * t761;
	t771 = t702 * t760 + t720 * t793;
	t764 = -t731 * t801 - t732 * t757 + t761 * t783 + (t740 * t802 + t757 * t790 - t808) * qJD(3);
	t770 = -t760 * t764 - t793 * t812;
	t714 = -t757 * t778 + (t774 * qJD(2) - t773 * qJD(3)) * t753;
	t733 = t775 * t753 + t761 * t784;
	t769 = -t714 * t760 + t733 * t793;
	t768 = -t729 * t752 + t755 * t788;
	t754 = cos(pkin(13));
	t751 = sin(pkin(13));
	t739 = t811 * t755 - t762 * t807;
	t738 = t772 * t753;
	t737 = t752 * t765 + t755 * t804;
	t734 = t773 * t753 + t757 * t784;
	t725 = (-t773 * qJD(2) + t774 * qJD(3)) * t753;
	t724 = (t775 * qJD(2) + t772 * qJD(3)) * t753;
	t715 = t761 * t778 + (t772 * qJD(2) + t775 * qJD(3)) * t753;
	t713 = t756 * t780 + t725 * t760 + (-t738 * t756 + t760 * t791) * qJD(4);
	t712 = -t732 * t802 - t731 * t761 + (t740 * t757 - t741 * t801) * qJD(3);
	t711 = t727 * qJD(3) - t731 * t757 + t732 * t801;
	t710 = t730 * t802 + t729 * t761 + (t757 * t765 + t766 * t801) * qJD(3);
	t709 = t728 * qJD(3) + t729 * t757 - t730 * t801;
	t708 = t760 * t781 - t715 * t756 + (-t734 * t760 - t739 * t756) * qJD(4);
	t703 = -t730 * t761 + (t729 * t755 + t782) * t757 + t720 * qJD(3);
	t701 = t732 * t806 + t712 * t760 + (-t727 * t756 + t741 * t805) * qJD(4);
	t700 = -t730 * t806 + t710 * t760 + (-t728 * t756 - t766 * t805) * qJD(4);
	t697 = t703 * t760 + t768 * t756 + (-t721 * t756 + t737 * t760) * qJD(4);
	t696 = t703 * t756 - t768 * t760 + (t721 * t760 + t737 * t756) * qJD(4);
	t1 = [t699 * t754 + t751 * t764, t700 * t754 + t709 * t751, t703 * t751 - t771 * t754, -t696 * t754, 0, 0; t697 * t754 + t702 * t751, t701 * t754 + t711 * t751, -t707 * t751 - t770 * t754, t698 * t754, 0, 0; 0, t713 * t754 + t724 * t751, t715 * t751 - t769 * t754, t708 * t754, 0, 0; -t699 * t751 + t754 * t764, -t700 * t751 + t709 * t754, t703 * t754 + t771 * t751, t696 * t751, 0, 0; -t697 * t751 + t702 * t754, -t701 * t751 + t711 * t754, -t707 * t754 + t770 * t751, -t698 * t751, 0, 0; 0, -t713 * t751 + t724 * t754, t715 * t754 + t769 * t751, -t708 * t751, 0, 0; t698, t730 * t805 + t710 * t756 + (t728 * t760 - t766 * t806) * qJD(4), -t702 * t756 + t720 * t792, t697, 0, 0; t696, -t732 * t805 + t712 * t756 + (t727 * t760 + t741 * t806) * qJD(4), t756 * t764 - t792 * t812, -t699, 0, 0; 0, -t760 * t780 + t725 * t756 + (t738 * t760 + t756 * t791) * qJD(4), t714 * t756 + t733 * t792, t756 * t781 + t715 * t760 + (-t734 * t756 + t739 * t760) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:26
	% EndTime: 2019-10-10 12:54:29
	% DurationCPUTime: 2.63s
	% Computational Cost: add. (1356->166), mult. (3901->300), div. (0->0), fcn. (4371->14), ass. (0->126)
	t872 = sin(qJ(2));
	t875 = cos(qJ(2));
	t876 = cos(qJ(1));
	t933 = cos(pkin(6));
	t905 = t876 * t933;
	t934 = sin(qJ(1));
	t851 = t872 * t905 + t934 * t875;
	t898 = t933 * t934;
	t880 = t876 * t872 + t875 * t898;
	t838 = t880 * qJD(1) + t851 * qJD(2);
	t892 = t872 * t898;
	t908 = t934 * t872;
	t918 = t876 * t875;
	t839 = -qJD(1) * t892 - qJD(2) * t908 + (qJD(2) * t933 + qJD(1)) * t918;
	t869 = cos(pkin(7));
	t871 = sin(qJ(3));
	t874 = cos(qJ(3));
	t867 = sin(pkin(7));
	t868 = sin(pkin(6));
	t909 = t868 * t934;
	t899 = qJD(1) * t909;
	t893 = t867 * t899;
	t850 = -t875 * t905 + t908;
	t925 = t868 * t876;
	t911 = t867 * t925;
	t890 = t850 * t869 + t911;
	t931 = t851 * t871;
	t935 = t890 * t874 + t931;
	t799 = t935 * qJD(3) - (-t838 * t869 + t893) * t871 - t839 * t874;
	t930 = t851 * t874;
	t822 = t890 * t871 - t930;
	t842 = -t850 * t867 + t869 * t925;
	t870 = sin(qJ(4));
	t873 = cos(qJ(4));
	t807 = t822 * t870 - t842 * t873;
	t826 = t838 * t867 + t869 * t899;
	t791 = -t807 * qJD(4) + t799 * t873 - t826 * t870;
	t923 = t869 * t874;
	t924 = t869 * t871;
	t798 = (t850 * t924 - t930) * qJD(3) - t838 * t923 + t874 * t893 - (-qJD(3) * t911 + t839) * t871;
	t866 = pkin(13) + qJ(6);
	t864 = sin(t866);
	t865 = cos(t866);
	t948 = -t791 * t864 + t798 * t865;
	t947 = t791 * t865 + t798 * t864;
	t809 = t822 * t873 + t842 * t870;
	t946 = t809 * t864;
	t945 = t809 * t865;
	t789 = t809 * qJD(4) + t799 * t870 + t826 * t873;
	t919 = t874 * t875;
	t922 = t871 * t872;
	t889 = t869 * t919 - t922;
	t881 = t892 - t918;
	t929 = t881 * t871;
	t928 = t867 * t868;
	t927 = t867 * t870;
	t926 = t867 * t873;
	t921 = t871 * t875;
	t920 = t872 * t874;
	t917 = qJD(4) * t870;
	t916 = qJD(4) * t873;
	t915 = qJD(6) * t864;
	t914 = qJD(6) * t865;
	t913 = qJD(6) * t873;
	t912 = t872 * t928;
	t907 = qJD(1) * t925;
	t906 = qJD(2) * t928;
	t904 = t933 * t867;
	t903 = t867 * t909;
	t902 = t867 * t907;
	t901 = t872 * t906;
	t900 = t875 * t906;
	t897 = qJD(3) * t904;
	t836 = t850 * qJD(1) + t881 * qJD(2);
	t837 = t851 * qJD(1) + t880 * qJD(2);
	t885 = -t869 * t880 + t903;
	t795 = -t837 * t874 + (t836 * t869 + t902) * t871 + (t885 * t874 + t929) * qJD(3);
	t823 = -t874 * t903 + t880 * t923 - t929;
	t896 = t823 * t913 + t795;
	t819 = t850 * t923 + t874 * t911 + t931;
	t895 = t819 * t913 - t799;
	t886 = -t869 * t922 + t919;
	t815 = t874 * t897 + (t886 * qJD(2) + t889 * qJD(3)) * t868;
	t840 = -t889 * t868 - t874 * t904;
	t894 = t840 * t913 + t815;
	t824 = t885 * t871 - t874 * t881;
	t844 = t867 * t880 + t869 * t909;
	t811 = t824 * t873 + t844 * t870;
	t810 = -t824 * t870 + t844 * t873;
	t887 = t869 * t921 + t920;
	t841 = t887 * t868 + t871 * t904;
	t849 = t933 * t869 - t875 * t928;
	t817 = t841 * t873 + t849 * t870;
	t816 = -t841 * t870 + t849 * t873;
	t832 = -t850 * t874 - t851 * t924;
	t812 = t832 * t873 + t851 * t927;
	t834 = -t874 * t880 + t881 * t924;
	t813 = t834 * t873 - t881 * t927;
	t831 = -t850 * t871 + t851 * t923;
	t833 = -t871 * t880 - t881 * t923;
	t888 = t869 * t920 + t921;
	t848 = t886 * t868;
	t835 = t848 * t873 + t870 * t912;
	t884 = -t836 * t867 + t869 * t907;
	t794 = t824 * qJD(3) - t836 * t923 - t837 * t871 - t874 * t902;
	t879 = qJD(6) * t824 - t794 * t873 + t823 * t917;
	t878 = -qJD(6) * t822 + t798 * t873 + t819 * t917;
	t814 = t871 * t897 + (t888 * qJD(2) + t887 * qJD(3)) * t868;
	t877 = qJD(6) * t841 - t814 * t873 + t840 * t917;
	t847 = t888 * t868;
	t828 = (-t887 * qJD(2) - t888 * qJD(3)) * t868;
	t827 = (t889 * qJD(2) + t886 * qJD(3)) * t868;
	t806 = t870 * t900 + t828 * t873 + (-t848 * t870 + t873 * t912) * qJD(4);
	t805 = -t831 * qJD(3) - t838 * t874 - t839 * t924;
	t804 = t832 * qJD(3) - t838 * t871 + t839 * t923;
	t803 = -t833 * qJD(3) + t836 * t874 + t837 * t924;
	t802 = t834 * qJD(3) + t836 * t871 - t837 * t923;
	t801 = t816 * qJD(4) + t815 * t873 + t870 * t901;
	t800 = -t817 * qJD(4) - t815 * t870 + t873 * t901;
	t793 = t839 * t927 + t805 * t873 + (-t832 * t870 + t851 * t926) * qJD(4);
	t792 = -t837 * t927 + t803 * t873 + (-t834 * t870 - t881 * t926) * qJD(4);
	t788 = t810 * qJD(4) + t795 * t873 + t884 * t870;
	t787 = t811 * qJD(4) + t795 * t870 - t884 * t873;
	t786 = t788 * t865 + t794 * t864 + (-t811 * t864 + t823 * t865) * qJD(6);
	t785 = -t788 * t864 + t794 * t865 + (-t811 * t865 - t823 * t864) * qJD(6);
	t1 = [(-t865 * t935 - t946) * qJD(6) + t947, t792 * t865 + t802 * t864 + (-t813 * t864 + t833 * t865) * qJD(6), t896 * t864 + t879 * t865, -t787 * t865 - t810 * t915, 0, t785; t786, t793 * t865 + t804 * t864 + (-t812 * t864 + t831 * t865) * qJD(6), t895 * t864 + t878 * t865, t789 * t865 - t807 * t915, 0, (-t819 * t864 + t945) * qJD(6) - t948; 0, t806 * t865 + t827 * t864 + (-t835 * t864 + t847 * t865) * qJD(6), t894 * t864 + t877 * t865, t800 * t865 - t816 * t915, 0, -t801 * t864 + t814 * t865 + (-t817 * t865 - t840 * t864) * qJD(6); (t864 * t935 - t945) * qJD(6) + t948, -t792 * t864 + t802 * t865 + (-t813 * t865 - t833 * t864) * qJD(6), -t879 * t864 + t896 * t865, t787 * t864 - t810 * t914, 0, -t786; t785, -t793 * t864 + t804 * t865 + (-t812 * t865 - t831 * t864) * qJD(6), -t878 * t864 + t895 * t865, -t789 * t864 - t807 * t914, 0, (-t819 * t865 - t946) * qJD(6) + t947; 0, -t806 * t864 + t827 * t865 + (-t835 * t865 - t847 * t864) * qJD(6), -t877 * t864 + t894 * t865, -t800 * t864 - t816 * t914, 0, -t801 * t865 - t814 * t864 + (t817 * t864 - t840 * t865) * qJD(6); t789, t813 * qJD(4) + t803 * t870 + t837 * t926, -t794 * t870 - t823 * t916, t788, 0, 0; t787, t812 * qJD(4) + t805 * t870 - t839 * t926, t798 * t870 - t819 * t916, -t791, 0, 0; 0, t835 * qJD(4) + t828 * t870 - t873 * t900, -t814 * t870 - t840 * t916, t801, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end