% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:16
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRRP12_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP12_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:15:47
	% EndTime: 2019-10-10 13:15:47
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:15:47
	% EndTime: 2019-10-10 13:15:47
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
	% StartTime: 2019-10-10 13:15:48
	% EndTime: 2019-10-10 13:15:48
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
	% StartTime: 2019-10-10 13:15:49
	% EndTime: 2019-10-10 13:15:50
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
	% StartTime: 2019-10-10 13:15:52
	% EndTime: 2019-10-10 13:15:54
	% DurationCPUTime: 1.14s
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
	% StartTime: 2019-10-10 13:16:01
	% EndTime: 2019-10-10 13:16:03
	% DurationCPUTime: 2.64s
	% Computational Cost: add. (1250->165), mult. (3901->300), div. (0->0), fcn. (4371->14), ass. (0->125)
	t860 = sin(qJ(2));
	t864 = cos(qJ(2));
	t865 = cos(qJ(1));
	t922 = cos(pkin(6));
	t894 = t865 * t922;
	t923 = sin(qJ(1));
	t841 = t860 * t894 + t923 * t864;
	t887 = t922 * t923;
	t869 = t865 * t860 + t864 * t887;
	t828 = t869 * qJD(1) + t841 * qJD(2);
	t881 = t860 * t887;
	t897 = t923 * t860;
	t907 = t865 * t864;
	t829 = -qJD(1) * t881 - qJD(2) * t897 + (qJD(2) * t922 + qJD(1)) * t907;
	t856 = cos(pkin(7));
	t859 = sin(qJ(3));
	t863 = cos(qJ(3));
	t854 = sin(pkin(7));
	t855 = sin(pkin(6));
	t898 = t855 * t923;
	t888 = qJD(1) * t898;
	t882 = t854 * t888;
	t840 = -t864 * t894 + t897;
	t914 = t855 * t865;
	t900 = t854 * t914;
	t879 = t840 * t856 + t900;
	t920 = t841 * t859;
	t924 = t879 * t863 + t920;
	t789 = t924 * qJD(3) - (-t828 * t856 + t882) * t859 - t829 * t863;
	t919 = t841 * t863;
	t812 = t879 * t859 - t919;
	t832 = -t840 * t854 + t856 * t914;
	t858 = sin(qJ(4));
	t862 = cos(qJ(4));
	t797 = t812 * t858 - t832 * t862;
	t816 = t828 * t854 + t856 * t888;
	t781 = -t797 * qJD(4) + t789 * t862 - t816 * t858;
	t912 = t856 * t863;
	t913 = t856 * t859;
	t788 = (t840 * t913 - t919) * qJD(3) - t828 * t912 + t863 * t882 - (-qJD(3) * t900 + t829) * t859;
	t857 = sin(qJ(5));
	t861 = cos(qJ(5));
	t937 = -t781 * t857 + t788 * t861;
	t936 = t781 * t861 + t788 * t857;
	t799 = t812 * t862 + t832 * t858;
	t935 = t799 * t857;
	t934 = t799 * t861;
	t779 = t799 * qJD(4) + t789 * t858 + t816 * t862;
	t908 = t863 * t864;
	t911 = t859 * t860;
	t878 = t856 * t908 - t911;
	t870 = t881 - t907;
	t918 = t870 * t859;
	t917 = t854 * t855;
	t916 = t854 * t858;
	t915 = t854 * t862;
	t910 = t859 * t864;
	t909 = t860 * t863;
	t906 = qJD(4) * t858;
	t905 = qJD(4) * t862;
	t904 = qJD(5) * t857;
	t903 = qJD(5) * t861;
	t902 = qJD(5) * t862;
	t901 = t860 * t917;
	t896 = qJD(1) * t914;
	t895 = qJD(2) * t917;
	t893 = t922 * t854;
	t892 = t854 * t898;
	t891 = t854 * t896;
	t890 = t860 * t895;
	t889 = t864 * t895;
	t886 = qJD(3) * t893;
	t826 = t840 * qJD(1) + t870 * qJD(2);
	t827 = t841 * qJD(1) + t869 * qJD(2);
	t874 = -t856 * t869 + t892;
	t785 = -t827 * t863 + (t826 * t856 + t891) * t859 + (t874 * t863 + t918) * qJD(3);
	t813 = -t863 * t892 + t869 * t912 - t918;
	t885 = t813 * t902 + t785;
	t809 = t840 * t912 + t863 * t900 + t920;
	t884 = t809 * t902 - t789;
	t875 = -t856 * t911 + t908;
	t805 = t863 * t886 + (t875 * qJD(2) + t878 * qJD(3)) * t855;
	t830 = -t878 * t855 - t863 * t893;
	t883 = t830 * t902 + t805;
	t814 = t874 * t859 - t863 * t870;
	t834 = t854 * t869 + t856 * t898;
	t801 = t814 * t862 + t834 * t858;
	t800 = -t814 * t858 + t834 * t862;
	t876 = t856 * t910 + t909;
	t831 = t876 * t855 + t859 * t893;
	t839 = t922 * t856 - t864 * t917;
	t807 = t831 * t862 + t839 * t858;
	t806 = -t831 * t858 + t839 * t862;
	t822 = -t840 * t863 - t841 * t913;
	t802 = t822 * t862 + t841 * t916;
	t824 = -t863 * t869 + t870 * t913;
	t803 = t824 * t862 - t870 * t916;
	t821 = -t840 * t859 + t841 * t912;
	t823 = -t859 * t869 - t870 * t912;
	t877 = t856 * t909 + t910;
	t838 = t875 * t855;
	t825 = t838 * t862 + t858 * t901;
	t873 = -t826 * t854 + t856 * t896;
	t784 = t814 * qJD(3) - t826 * t912 - t827 * t859 - t863 * t891;
	t868 = qJD(5) * t814 - t784 * t862 + t813 * t906;
	t867 = -qJD(5) * t812 + t788 * t862 + t809 * t906;
	t804 = t859 * t886 + (t877 * qJD(2) + t876 * qJD(3)) * t855;
	t866 = qJD(5) * t831 - t804 * t862 + t830 * t906;
	t837 = t877 * t855;
	t818 = (-t876 * qJD(2) - t877 * qJD(3)) * t855;
	t817 = (t878 * qJD(2) + t875 * qJD(3)) * t855;
	t796 = t858 * t889 + t818 * t862 + (-t838 * t858 + t862 * t901) * qJD(4);
	t795 = -t821 * qJD(3) - t828 * t863 - t829 * t913;
	t794 = t822 * qJD(3) - t828 * t859 + t829 * t912;
	t793 = -t823 * qJD(3) + t826 * t863 + t827 * t913;
	t792 = t824 * qJD(3) + t826 * t859 - t827 * t912;
	t791 = t806 * qJD(4) + t805 * t862 + t858 * t890;
	t790 = -t807 * qJD(4) - t805 * t858 + t862 * t890;
	t783 = t829 * t916 + t795 * t862 + (-t822 * t858 + t841 * t915) * qJD(4);
	t782 = -t827 * t916 + t793 * t862 + (-t824 * t858 - t870 * t915) * qJD(4);
	t778 = t800 * qJD(4) + t785 * t862 + t873 * t858;
	t777 = t801 * qJD(4) + t785 * t858 - t873 * t862;
	t776 = t778 * t861 + t784 * t857 + (-t801 * t857 + t813 * t861) * qJD(5);
	t775 = -t778 * t857 + t784 * t861 + (-t801 * t861 - t813 * t857) * qJD(5);
	t1 = [(-t861 * t924 - t935) * qJD(5) + t936, t782 * t861 + t792 * t857 + (-t803 * t857 + t823 * t861) * qJD(5), t885 * t857 + t868 * t861, -t777 * t861 - t800 * t904, t775, 0; t776, t783 * t861 + t794 * t857 + (-t802 * t857 + t821 * t861) * qJD(5), t884 * t857 + t867 * t861, t779 * t861 - t797 * t904, (-t809 * t857 + t934) * qJD(5) - t937, 0; 0, t796 * t861 + t817 * t857 + (-t825 * t857 + t837 * t861) * qJD(5), t883 * t857 + t866 * t861, t790 * t861 - t806 * t904, -t791 * t857 + t804 * t861 + (-t807 * t861 - t830 * t857) * qJD(5), 0; (t857 * t924 - t934) * qJD(5) + t937, -t782 * t857 + t792 * t861 + (-t803 * t861 - t823 * t857) * qJD(5), -t868 * t857 + t885 * t861, t777 * t857 - t800 * t903, -t776, 0; t775, -t783 * t857 + t794 * t861 + (-t802 * t861 - t821 * t857) * qJD(5), -t867 * t857 + t884 * t861, -t779 * t857 - t797 * t903, (-t809 * t861 - t935) * qJD(5) + t936, 0; 0, -t796 * t857 + t817 * t861 + (-t825 * t861 - t837 * t857) * qJD(5), -t866 * t857 + t883 * t861, -t790 * t857 - t806 * t903, -t791 * t861 - t804 * t857 + (t807 * t857 - t830 * t861) * qJD(5), 0; t779, t803 * qJD(4) + t793 * t858 + t827 * t915, -t784 * t858 - t813 * t905, t778, 0, 0; t777, t802 * qJD(4) + t795 * t858 - t829 * t915, t788 * t858 - t809 * t905, -t781, 0, 0; 0, t825 * qJD(4) + t818 * t858 - t862 * t889, -t804 * t858 - t830 * t905, t791, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:16:06
	% EndTime: 2019-10-10 13:16:08
	% DurationCPUTime: 2.49s
	% Computational Cost: add. (1250->165), mult. (3901->300), div. (0->0), fcn. (4371->14), ass. (0->125)
	t1034 = cos(pkin(6));
	t977 = cos(qJ(1));
	t1006 = t977 * t1034;
	t1035 = sin(qJ(1));
	t972 = sin(qJ(2));
	t976 = cos(qJ(2));
	t953 = t972 * t1006 + t1035 * t976;
	t971 = sin(qJ(3));
	t1032 = t953 * t971;
	t975 = cos(qJ(3));
	t967 = sin(pkin(6));
	t1026 = t967 * t977;
	t966 = sin(pkin(7));
	t1012 = t966 * t1026;
	t1009 = t1035 * t972;
	t952 = -t976 * t1006 + t1009;
	t968 = cos(pkin(7));
	t991 = t952 * t968 + t1012;
	t1037 = t991 * t975 + t1032;
	t999 = t1034 * t1035;
	t981 = t977 * t972 + t976 * t999;
	t940 = t981 * qJD(1) + t953 * qJD(2);
	t1019 = t977 * t976;
	t993 = t972 * t999;
	t941 = -qJD(1) * t993 - qJD(2) * t1009 + (qJD(2) * t1034 + qJD(1)) * t1019;
	t1010 = t967 * t1035;
	t1000 = qJD(1) * t1010;
	t994 = t966 * t1000;
	t901 = t1037 * qJD(3) - (-t940 * t968 + t994) * t971 - t941 * t975;
	t1031 = t953 * t975;
	t924 = t991 * t971 - t1031;
	t944 = t968 * t1026 - t952 * t966;
	t970 = sin(qJ(4));
	t974 = cos(qJ(4));
	t909 = t924 * t970 - t944 * t974;
	t928 = t968 * t1000 + t940 * t966;
	t893 = -t909 * qJD(4) + t901 * t974 - t928 * t970;
	t1024 = t968 * t975;
	t1025 = t968 * t971;
	t900 = (t952 * t1025 - t1031) * qJD(3) - t940 * t1024 + t975 * t994 - (-qJD(3) * t1012 + t941) * t971;
	t969 = sin(qJ(5));
	t973 = cos(qJ(5));
	t1049 = t893 * t969 - t900 * t973;
	t1048 = -t893 * t973 - t900 * t969;
	t911 = t924 * t974 + t944 * t970;
	t1047 = t911 * t969;
	t1046 = t911 * t973;
	t891 = t911 * qJD(4) + t901 * t970 + t928 * t974;
	t1020 = t975 * t976;
	t1023 = t971 * t972;
	t990 = t968 * t1020 - t1023;
	t982 = t993 - t1019;
	t1030 = t982 * t971;
	t1029 = t966 * t967;
	t1028 = t966 * t970;
	t1027 = t966 * t974;
	t1022 = t971 * t976;
	t1021 = t972 * t975;
	t1018 = qJD(4) * t970;
	t1017 = qJD(4) * t974;
	t1016 = qJD(5) * t969;
	t1015 = qJD(5) * t973;
	t1014 = qJD(5) * t974;
	t1013 = t972 * t1029;
	t1008 = qJD(1) * t1026;
	t1007 = qJD(2) * t1029;
	t1005 = t1034 * t966;
	t1004 = t966 * t1010;
	t1003 = t966 * t1008;
	t1002 = t972 * t1007;
	t1001 = t976 * t1007;
	t998 = qJD(3) * t1005;
	t938 = t952 * qJD(1) + t982 * qJD(2);
	t939 = t953 * qJD(1) + t981 * qJD(2);
	t986 = -t968 * t981 + t1004;
	t897 = -t939 * t975 + (t938 * t968 + t1003) * t971 + (t986 * t975 + t1030) * qJD(3);
	t925 = -t975 * t1004 + t1024 * t981 - t1030;
	t997 = t925 * t1014 + t897;
	t921 = t975 * t1012 + t952 * t1024 + t1032;
	t996 = t921 * t1014 - t901;
	t987 = -t968 * t1023 + t1020;
	t917 = t975 * t998 + (t987 * qJD(2) + t990 * qJD(3)) * t967;
	t942 = -t975 * t1005 - t990 * t967;
	t995 = t942 * t1014 + t917;
	t926 = t986 * t971 - t975 * t982;
	t946 = t968 * t1010 + t966 * t981;
	t913 = t926 * t974 + t946 * t970;
	t912 = -t926 * t970 + t946 * t974;
	t988 = t968 * t1022 + t1021;
	t943 = t971 * t1005 + t988 * t967;
	t951 = -t976 * t1029 + t1034 * t968;
	t919 = t943 * t974 + t951 * t970;
	t918 = -t943 * t970 + t951 * t974;
	t934 = -t953 * t1025 - t952 * t975;
	t914 = t953 * t1028 + t934 * t974;
	t936 = t1025 * t982 - t975 * t981;
	t915 = -t1028 * t982 + t936 * t974;
	t933 = t953 * t1024 - t952 * t971;
	t935 = -t1024 * t982 - t971 * t981;
	t989 = t968 * t1021 + t1022;
	t950 = t987 * t967;
	t937 = t970 * t1013 + t950 * t974;
	t985 = t968 * t1008 - t938 * t966;
	t896 = t926 * qJD(3) - t975 * t1003 - t938 * t1024 - t939 * t971;
	t980 = qJD(5) * t926 + t925 * t1018 - t896 * t974;
	t979 = -qJD(5) * t924 + t921 * t1018 + t900 * t974;
	t916 = t971 * t998 + (t989 * qJD(2) + t988 * qJD(3)) * t967;
	t978 = qJD(5) * t943 + t942 * t1018 - t916 * t974;
	t949 = t989 * t967;
	t930 = (-t988 * qJD(2) - t989 * qJD(3)) * t967;
	t929 = (t990 * qJD(2) + t987 * qJD(3)) * t967;
	t908 = t970 * t1001 + t930 * t974 + (t974 * t1013 - t950 * t970) * qJD(4);
	t907 = -t933 * qJD(3) - t941 * t1025 - t940 * t975;
	t906 = t934 * qJD(3) + t941 * t1024 - t940 * t971;
	t905 = -t935 * qJD(3) + t939 * t1025 + t938 * t975;
	t904 = t936 * qJD(3) - t939 * t1024 + t938 * t971;
	t903 = t918 * qJD(4) + t970 * t1002 + t917 * t974;
	t902 = -t919 * qJD(4) + t974 * t1002 - t917 * t970;
	t895 = t941 * t1028 + t907 * t974 + (t953 * t1027 - t934 * t970) * qJD(4);
	t894 = -t939 * t1028 + t905 * t974 + (-t1027 * t982 - t936 * t970) * qJD(4);
	t890 = t912 * qJD(4) + t897 * t974 + t985 * t970;
	t889 = t913 * qJD(4) + t897 * t970 - t985 * t974;
	t888 = t890 * t973 + t896 * t969 + (-t913 * t969 + t925 * t973) * qJD(5);
	t887 = t890 * t969 - t896 * t973 + (t913 * t973 + t925 * t969) * qJD(5);
	t1 = [(-t1037 * t973 - t1047) * qJD(5) - t1048, t894 * t973 + t904 * t969 + (-t915 * t969 + t935 * t973) * qJD(5), t997 * t969 + t980 * t973, -t912 * t1016 - t889 * t973, -t887, 0; t888, t895 * t973 + t906 * t969 + (-t914 * t969 + t933 * t973) * qJD(5), t996 * t969 + t979 * t973, -t909 * t1016 + t891 * t973, (-t921 * t969 + t1046) * qJD(5) + t1049, 0; 0, t908 * t973 + t929 * t969 + (-t937 * t969 + t949 * t973) * qJD(5), t995 * t969 + t978 * t973, -t918 * t1016 + t902 * t973, -t903 * t969 + t916 * t973 + (-t919 * t973 - t942 * t969) * qJD(5), 0; t891, t915 * qJD(4) + t939 * t1027 + t905 * t970, -t925 * t1017 - t896 * t970, t890, 0, 0; t889, t914 * qJD(4) - t941 * t1027 + t907 * t970, -t921 * t1017 + t900 * t970, -t893, 0, 0; 0, t937 * qJD(4) - t974 * t1001 + t930 * t970, -t942 * t1017 - t916 * t970, t903, 0, 0; (-t1037 * t969 + t1046) * qJD(5) + t1049, t894 * t969 - t904 * t973 + (t915 * t973 + t935 * t969) * qJD(5), t980 * t969 - t997 * t973, t912 * t1015 - t889 * t969, t888, 0; t887, t895 * t969 - t906 * t973 + (t914 * t973 + t933 * t969) * qJD(5), t979 * t969 - t996 * t973, t909 * t1015 + t891 * t969, (t921 * t973 + t1047) * qJD(5) + t1048, 0; 0, t908 * t969 - t929 * t973 + (t937 * t973 + t949 * t969) * qJD(5), t978 * t969 - t995 * t973, t918 * t1015 + t902 * t969, t903 * t973 + t916 * t969 + (-t919 * t969 + t942 * t973) * qJD(5), 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end