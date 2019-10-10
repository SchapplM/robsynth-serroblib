% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:56
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRPR15_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR15_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
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
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
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
	% StartTime: 2019-10-10 12:56:28
	% EndTime: 2019-10-10 12:56:28
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
	% StartTime: 2019-10-10 12:56:31
	% EndTime: 2019-10-10 12:56:32
	% DurationCPUTime: 1.10s
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
	% StartTime: 2019-10-10 12:56:36
	% EndTime: 2019-10-10 12:56:37
	% DurationCPUTime: 1.07s
	% Computational Cost: add. (493->104), mult. (1577->199), div. (0->0), fcn. (1713->12), ass. (0->85)
	t691 = sin(qJ(2));
	t692 = sin(qJ(1));
	t695 = cos(qJ(2));
	t696 = cos(qJ(1));
	t739 = cos(pkin(6));
	t713 = t696 * t739;
	t676 = t691 * t713 + t692 * t695;
	t714 = t692 * t739;
	t697 = t691 * t696 + t695 * t714;
	t666 = qJD(1) * t697 + qJD(2) * t676;
	t707 = t691 * t714;
	t723 = t696 * t695;
	t725 = t692 * t691;
	t667 = -qJD(1) * t707 - qJD(2) * t725 + (qJD(2) * t739 + qJD(1)) * t723;
	t688 = cos(pkin(7));
	t690 = sin(qJ(3));
	t694 = cos(qJ(3));
	t686 = sin(pkin(7));
	t687 = sin(pkin(6));
	t722 = qJD(1) * t687;
	t717 = t692 * t722;
	t711 = t686 * t717;
	t675 = -t695 * t713 + t725;
	t731 = t687 * t696;
	t718 = t686 * t731;
	t705 = t675 * t688 + t718;
	t740 = t676 * t690 + t694 * t705;
	t647 = t740 * qJD(3) - (-t666 * t688 + t711) * t690 - t667 * t694;
	t736 = t676 * t694;
	t654 = t690 * t705 - t736;
	t659 = t666 * t686 + t688 * t717;
	t670 = -t675 * t686 + t688 * t731;
	t689 = sin(qJ(4));
	t693 = cos(qJ(4));
	t750 = t647 * t689 + t659 * t693 + (t654 * t693 + t670 * t689) * qJD(4);
	t749 = -t647 * t693 + t659 * t689 + (t654 * t689 - t670 * t693) * qJD(4);
	t735 = t686 * t687;
	t734 = t686 * t689;
	t733 = t686 * t693;
	t732 = t687 * t692;
	t730 = t688 * t690;
	t729 = t688 * t694;
	t728 = t690 * t691;
	t727 = t690 * t695;
	t726 = t691 * t694;
	t724 = t694 * t695;
	t721 = qJD(4) * t689;
	t720 = qJD(4) * t693;
	t719 = t691 * t735;
	t716 = t696 * t722;
	t715 = qJD(2) * t735;
	t712 = t739 * t686;
	t710 = t686 * t716;
	t709 = t691 * t715;
	t708 = t695 * t715;
	t706 = qJD(3) * t712;
	t662 = -t675 * t694 - t676 * t730;
	t698 = t707 - t723;
	t663 = -t694 * t697 + t698 * t730;
	t704 = t686 * t732 - t688 * t697;
	t703 = t688 * t724 - t728;
	t702 = -t688 * t726 - t727;
	t701 = t688 * t727 + t726;
	t700 = -t688 * t728 + t724;
	t655 = t690 * t698 + t694 * t704;
	t656 = t690 * t704 - t694 * t698;
	t645 = -t666 * t729 - t667 * t690 + t694 * t711 + (t675 * t730 + t690 * t718 - t736) * qJD(3);
	t674 = t688 * t739 - t695 * t735;
	t673 = t700 * t687;
	t672 = t686 * t697 + t688 * t732;
	t669 = t687 * t701 + t690 * t712;
	t668 = t687 * t703 + t694 * t712;
	t665 = qJD(1) * t676 + qJD(2) * t697;
	t664 = qJD(1) * t675 + qJD(2) * t698;
	t660 = (-qJD(2) * t701 + qJD(3) * t702) * t687;
	t657 = -t664 * t686 + t688 * t716;
	t651 = t694 * t706 + (qJD(2) * t700 + qJD(3) * t703) * t687;
	t650 = -t690 * t706 + (qJD(2) * t702 - qJD(3) * t701) * t687;
	t649 = -t667 * t730 - t666 * t694 + (t675 * t690 - t676 * t729) * qJD(3);
	t648 = t665 * t730 + t664 * t694 + (t690 * t697 + t698 * t729) * qJD(3);
	t644 = -t665 * t694 + (t664 * t688 + t710) * t690 + t655 * qJD(3);
	t643 = qJD(3) * t656 - t664 * t729 - t665 * t690 - t694 * t710;
	t642 = t644 * t693 + t657 * t689 + (-t656 * t689 + t672 * t693) * qJD(4);
	t641 = t644 * t689 - t657 * t693 + (t656 * t693 + t672 * t689) * qJD(4);
	t1 = [t645, qJD(3) * t663 + t664 * t690 - t665 * t729, t644, 0, 0, 0; t643, qJD(3) * t662 - t666 * t690 + t667 * t729, -t647, 0, 0, 0; 0, (qJD(2) * t703 + qJD(3) * t700) * t687, t651, 0, 0, 0; t749, t665 * t734 - t648 * t693 + (t663 * t689 + t698 * t733) * qJD(4), t643 * t693 + t655 * t721, t641, 0, 0; -t642, -t667 * t734 - t649 * t693 + (t662 * t689 - t676 * t733) * qJD(4), -t645 * t693 - t721 * t740, -t750, 0, 0; 0, -t689 * t708 - t660 * t693 + (t673 * t689 - t693 * t719) * qJD(4), -t650 * t693 + t668 * t721, -t693 * t709 + t651 * t689 + (t669 * t693 + t674 * t689) * qJD(4), 0, 0; t750, t665 * t733 + t648 * t689 + (t663 * t693 - t698 * t734) * qJD(4), -t643 * t689 + t655 * t720, t642, 0, 0; t641, -t667 * t733 + t649 * t689 + (t662 * t693 + t676 * t734) * qJD(4), t645 * t689 - t720 * t740, t749, 0, 0; 0, -t693 * t708 + t660 * t689 + (t673 * t693 + t689 * t719) * qJD(4), t650 * t689 + t668 * t720, t689 * t709 + t651 * t693 + (-t669 * t689 + t674 * t693) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:38
	% EndTime: 2019-10-10 12:56:40
	% DurationCPUTime: 2.55s
	% Computational Cost: add. (1250->167), mult. (3901->300), div. (0->0), fcn. (4371->14), ass. (0->125)
	t854 = sin(qJ(2));
	t858 = cos(qJ(2));
	t859 = cos(qJ(1));
	t918 = cos(pkin(6));
	t890 = t859 * t918;
	t919 = sin(qJ(1));
	t835 = t854 * t890 + t919 * t858;
	t883 = t918 * t919;
	t863 = t859 * t854 + t858 * t883;
	t822 = t863 * qJD(1) + t835 * qJD(2);
	t877 = t854 * t883;
	t893 = t919 * t854;
	t903 = t859 * t858;
	t823 = -qJD(1) * t877 - qJD(2) * t893 + (qJD(2) * t918 + qJD(1)) * t903;
	t850 = cos(pkin(7));
	t853 = sin(qJ(3));
	t857 = cos(qJ(3));
	t848 = sin(pkin(7));
	t849 = sin(pkin(6));
	t894 = t849 * t919;
	t884 = qJD(1) * t894;
	t878 = t848 * t884;
	t834 = -t858 * t890 + t893;
	t910 = t849 * t859;
	t896 = t848 * t910;
	t873 = t834 * t850 + t896;
	t916 = t835 * t853;
	t920 = t873 * t857 + t916;
	t783 = t920 * qJD(3) - (-t822 * t850 + t878) * t853 - t823 * t857;
	t810 = t822 * t848 + t850 * t884;
	t852 = sin(qJ(4));
	t856 = cos(qJ(4));
	t915 = t835 * t857;
	t805 = t873 * t853 - t915;
	t826 = -t834 * t848 + t850 * t910;
	t926 = t805 * t856 + t826 * t852;
	t775 = t926 * qJD(4) + t783 * t852 + t810 * t856;
	t908 = t850 * t857;
	t909 = t850 * t853;
	t782 = (t834 * t909 - t915) * qJD(3) - t822 * t908 + t857 * t878 - (-qJD(3) * t896 + t823) * t853;
	t851 = sin(qJ(6));
	t855 = cos(qJ(6));
	t932 = t775 * t851 + t782 * t855;
	t931 = -t775 * t855 + t782 * t851;
	t793 = t805 * t852 - t826 * t856;
	t930 = t793 * t851;
	t929 = t793 * t855;
	t774 = t793 * qJD(4) - t783 * t856 + t810 * t852;
	t904 = t857 * t858;
	t907 = t853 * t854;
	t872 = t850 * t904 - t907;
	t864 = t877 - t903;
	t914 = t864 * t853;
	t913 = t848 * t849;
	t912 = t848 * t852;
	t911 = t848 * t856;
	t906 = t853 * t858;
	t905 = t854 * t857;
	t902 = qJD(4) * t852;
	t901 = qJD(4) * t856;
	t900 = qJD(6) * t851;
	t899 = qJD(6) * t852;
	t898 = qJD(6) * t855;
	t897 = t854 * t913;
	t892 = qJD(1) * t910;
	t891 = qJD(2) * t913;
	t889 = t918 * t848;
	t888 = t848 * t894;
	t887 = t848 * t892;
	t886 = t854 * t891;
	t885 = t858 * t891;
	t882 = qJD(3) * t889;
	t820 = t834 * qJD(1) + t864 * qJD(2);
	t821 = t835 * qJD(1) + t863 * qJD(2);
	t867 = -t850 * t863 + t888;
	t779 = -t821 * t857 + (t820 * t850 + t887) * t853 + (t867 * t857 + t914) * qJD(3);
	t806 = -t857 * t888 + t863 * t908 - t914;
	t881 = t806 * t899 - t779;
	t802 = t834 * t908 + t857 * t896 + t916;
	t880 = t802 * t899 + t783;
	t869 = -t850 * t907 + t904;
	t799 = t857 * t882 + (t869 * qJD(2) + t872 * qJD(3)) * t849;
	t824 = -t872 * t849 - t857 * t889;
	t879 = t824 * t899 - t799;
	t807 = t867 * t853 - t857 * t864;
	t828 = t848 * t863 + t850 * t894;
	t795 = t807 * t856 + t828 * t852;
	t794 = t807 * t852 - t828 * t856;
	t870 = t850 * t906 + t905;
	t825 = t870 * t849 + t853 * t889;
	t833 = t918 * t850 - t858 * t913;
	t801 = t825 * t856 + t833 * t852;
	t800 = t825 * t852 - t833 * t856;
	t816 = -t834 * t857 - t835 * t909;
	t876 = -t816 * t852 + t835 * t911;
	t818 = -t857 * t863 + t864 * t909;
	t875 = -t818 * t852 - t864 * t911;
	t815 = -t834 * t853 + t835 * t908;
	t817 = -t853 * t863 - t864 * t908;
	t871 = t850 * t905 + t906;
	t832 = t869 * t849;
	t868 = -t832 * t852 + t856 * t897;
	t778 = t807 * qJD(3) - t820 * t908 - t821 * t853 - t857 * t887;
	t862 = -qJD(6) * t807 - t778 * t852 - t806 * t901;
	t861 = qJD(6) * t805 + t782 * t852 - t802 * t901;
	t798 = t853 * t882 + (t871 * qJD(2) + t870 * qJD(3)) * t849;
	t860 = -qJD(6) * t825 - t798 * t852 - t824 * t901;
	t831 = t871 * t849;
	t812 = (-t870 * qJD(2) - t871 * qJD(3)) * t849;
	t811 = (t872 * qJD(2) + t869 * qJD(3)) * t849;
	t808 = -t820 * t848 + t850 * t892;
	t790 = -t856 * t885 + t812 * t852 + (t832 * t856 + t852 * t897) * qJD(4);
	t789 = -t815 * qJD(3) - t822 * t857 - t823 * t909;
	t788 = t816 * qJD(3) - t822 * t853 + t823 * t908;
	t787 = -t817 * qJD(3) + t820 * t857 + t821 * t909;
	t786 = t818 * qJD(3) + t820 * t853 - t821 * t908;
	t785 = -t800 * qJD(4) + t799 * t856 + t852 * t886;
	t784 = t801 * qJD(4) + t799 * t852 - t856 * t886;
	t777 = -t823 * t911 + t789 * t852 + (t816 * t856 + t835 * t912) * qJD(4);
	t776 = t821 * t911 + t787 * t852 + (t818 * t856 - t864 * t912) * qJD(4);
	t772 = -t794 * qJD(4) + t779 * t856 + t808 * t852;
	t771 = t795 * qJD(4) + t779 * t852 - t808 * t856;
	t770 = t771 * t851 + t778 * t855 + (t794 * t855 - t806 * t851) * qJD(6);
	t769 = t771 * t855 - t778 * t851 + (-t794 * t851 - t806 * t855) * qJD(6);
	t1 = [(t851 * t920 + t929) * qJD(6) + t932, t776 * t851 + t786 * t855 + (-t817 * t851 - t855 * t875) * qJD(6), t862 * t851 - t881 * t855, t772 * t851 + t795 * t898, 0, t769; t770, t777 * t851 + t788 * t855 + (-t815 * t851 - t855 * t876) * qJD(6), t861 * t851 - t880 * t855, t774 * t851 - t898 * t926, 0, (-t802 * t855 + t930) * qJD(6) + t931; 0, t790 * t851 + t811 * t855 + (-t831 * t851 - t855 * t868) * qJD(6), t860 * t851 - t879 * t855, t785 * t851 + t801 * t898, 0, t784 * t855 - t798 * t851 + (-t800 * t851 - t824 * t855) * qJD(6); (t855 * t920 - t930) * qJD(6) - t931, t776 * t855 - t786 * t851 + (-t817 * t855 + t851 * t875) * qJD(6), t881 * t851 + t862 * t855, t772 * t855 - t795 * t900, 0, -t770; t769, t777 * t855 - t788 * t851 + (-t815 * t855 + t851 * t876) * qJD(6), t880 * t851 + t861 * t855, t774 * t855 + t900 * t926, 0, (t802 * t851 + t929) * qJD(6) + t932; 0, t790 * t855 - t811 * t851 + (-t831 * t855 + t851 * t868) * qJD(6), t879 * t851 + t860 * t855, t785 * t855 - t801 * t900, 0, -t784 * t851 - t798 * t855 + (-t800 * t855 + t824 * t851) * qJD(6); -t774, t875 * qJD(4) + t787 * t856 - t821 * t912, -t778 * t856 + t806 * t902, -t771, 0, 0; t772, t876 * qJD(4) + t789 * t856 + t823 * t912, t782 * t856 + t802 * t902, t775, 0, 0; 0, t868 * qJD(4) + t812 * t856 + t852 * t885, -t798 * t856 + t824 * t902, -t784, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end