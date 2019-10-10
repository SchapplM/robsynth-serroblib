% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR9
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
% Datum: 2019-10-10 13:32
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRRR9_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR9_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:11
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:11
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
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:11
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
	% StartTime: 2019-10-10 13:32:13
	% EndTime: 2019-10-10 13:32:13
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
	% StartTime: 2019-10-10 13:32:16
	% EndTime: 2019-10-10 13:32:17
	% DurationCPUTime: 1.27s
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
	% StartTime: 2019-10-10 13:32:24
	% EndTime: 2019-10-10 13:32:27
	% DurationCPUTime: 2.65s
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
	t781 = -qJD(4) * t797 + t789 * t862 - t816 * t858;
	t912 = t856 * t863;
	t913 = t856 * t859;
	t788 = qJD(3) * (t840 * t913 - t919) - t828 * t912 + t863 * t882 - (-qJD(3) * t900 + t829) * t859;
	t857 = sin(qJ(5));
	t861 = cos(qJ(5));
	t937 = -t781 * t857 + t788 * t861;
	t936 = t781 * t861 + t788 * t857;
	t799 = t812 * t862 + t832 * t858;
	t935 = t799 * t857;
	t934 = t799 * t861;
	t779 = qJD(4) * t799 + t789 * t858 + t816 * t862;
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
	t785 = -t827 * t863 + (t826 * t856 + t891) * t859 + (t863 * t874 + t918) * qJD(3);
	t813 = -t863 * t892 + t869 * t912 - t918;
	t885 = t813 * t902 + t785;
	t809 = t840 * t912 + t863 * t900 + t920;
	t884 = t809 * t902 - t789;
	t875 = -t856 * t911 + t908;
	t805 = t863 * t886 + (qJD(2) * t875 + qJD(3) * t878) * t855;
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
	t784 = qJD(3) * t814 - t826 * t912 - t827 * t859 - t863 * t891;
	t868 = qJD(5) * t814 - t784 * t862 + t813 * t906;
	t867 = -qJD(5) * t812 + t788 * t862 + t809 * t906;
	t804 = t859 * t886 + (qJD(2) * t877 + qJD(3) * t876) * t855;
	t866 = qJD(5) * t831 - t804 * t862 + t830 * t906;
	t837 = t877 * t855;
	t818 = (-t876 * qJD(2) - t877 * qJD(3)) * t855;
	t817 = (t878 * qJD(2) + t875 * qJD(3)) * t855;
	t796 = t858 * t889 + t818 * t862 + (-t838 * t858 + t862 * t901) * qJD(4);
	t795 = -qJD(3) * t821 - t828 * t863 - t829 * t913;
	t794 = qJD(3) * t822 - t828 * t859 + t829 * t912;
	t793 = -qJD(3) * t823 + t826 * t863 + t827 * t913;
	t792 = qJD(3) * t824 + t826 * t859 - t827 * t912;
	t791 = qJD(4) * t806 + t805 * t862 + t858 * t890;
	t790 = -qJD(4) * t807 - t805 * t858 + t862 * t890;
	t783 = t829 * t916 + t795 * t862 + (-t822 * t858 + t841 * t915) * qJD(4);
	t782 = -t827 * t916 + t793 * t862 + (-t824 * t858 - t870 * t915) * qJD(4);
	t778 = qJD(4) * t800 + t785 * t862 + t858 * t873;
	t777 = qJD(4) * t801 + t785 * t858 - t862 * t873;
	t776 = t778 * t861 + t784 * t857 + (-t801 * t857 + t813 * t861) * qJD(5);
	t775 = -t778 * t857 + t784 * t861 + (-t801 * t861 - t813 * t857) * qJD(5);
	t1 = [(-t861 * t924 - t935) * qJD(5) + t936, t782 * t861 + t792 * t857 + (-t803 * t857 + t823 * t861) * qJD(5), t857 * t885 + t861 * t868, -t777 * t861 - t800 * t904, t775, 0; t776, t783 * t861 + t794 * t857 + (-t802 * t857 + t821 * t861) * qJD(5), t857 * t884 + t861 * t867, t779 * t861 - t797 * t904, (-t809 * t857 + t934) * qJD(5) - t937, 0; 0, t796 * t861 + t817 * t857 + (-t825 * t857 + t837 * t861) * qJD(5), t857 * t883 + t861 * t866, t790 * t861 - t806 * t904, -t791 * t857 + t804 * t861 + (-t807 * t861 - t830 * t857) * qJD(5), 0; (t857 * t924 - t934) * qJD(5) + t937, -t782 * t857 + t792 * t861 + (-t803 * t861 - t823 * t857) * qJD(5), -t857 * t868 + t861 * t885, t777 * t857 - t800 * t903, -t776, 0; t775, -t783 * t857 + t794 * t861 + (-t802 * t861 - t821 * t857) * qJD(5), -t857 * t867 + t861 * t884, -t779 * t857 - t797 * t903, (-t809 * t861 - t935) * qJD(5) + t936, 0; 0, -t796 * t857 + t817 * t861 + (-t825 * t861 - t837 * t857) * qJD(5), -t857 * t866 + t861 * t883, -t790 * t857 - t806 * t903, -t791 * t861 - t804 * t857 + (t807 * t857 - t830 * t861) * qJD(5), 0; t779, qJD(4) * t803 + t793 * t858 + t827 * t915, -t784 * t858 - t813 * t905, t778, 0, 0; t777, qJD(4) * t802 + t795 * t858 - t829 * t915, t788 * t858 - t809 * t905, -t781, 0, 0; 0, qJD(4) * t825 + t818 * t858 - t862 * t889, -t804 * t858 - t830 * t905, t791, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:32:25
	% EndTime: 2019-10-10 13:32:27
	% DurationCPUTime: 2.24s
	% Computational Cost: add. (1700->154), mult. (4697->275), div. (0->0), fcn. (5269->14), ass. (0->135)
	t1005 = sin(qJ(1));
	t928 = sin(qJ(2));
	t931 = cos(qJ(2));
	t1004 = cos(pkin(6));
	t932 = cos(qJ(1));
	t975 = t932 * t1004;
	t906 = t1005 * t931 + t928 * t975;
	t930 = cos(qJ(3));
	t1001 = t906 * t930;
	t956 = t1004 * t1005;
	t936 = t932 * t928 + t931 * t956;
	t893 = t936 * qJD(1) + t906 * qJD(2);
	t949 = t928 * t956;
	t979 = t1005 * t928;
	t986 = t932 * t931;
	t894 = -qJD(1) * t949 - qJD(2) * t979 + (qJD(2) * t1004 + qJD(1)) * t986;
	t905 = -t931 * t975 + t979;
	t927 = sin(qJ(3));
	t923 = sin(pkin(7));
	t924 = sin(pkin(6));
	t980 = t924 * t1005;
	t957 = qJD(1) * t980;
	t950 = t923 * t957;
	t993 = t924 * t932;
	t982 = t923 * t993;
	t925 = cos(pkin(7));
	t991 = t925 * t930;
	t992 = t925 * t927;
	t851 = -(t905 * t992 - t1001) * qJD(3) + t893 * t991 - t930 * t950 + (-qJD(3) * t982 + t894) * t927;
	t947 = t905 * t925 + t982;
	t877 = t947 * t927 - t1001;
	t897 = -t905 * t923 + t925 * t993;
	t926 = sin(qJ(4));
	t929 = cos(qJ(4));
	t865 = t877 * t929 + t897 * t926;
	t921 = qJD(5) + qJD(6);
	t1015 = -t865 * t921 - t851;
	t922 = qJ(5) + qJ(6);
	t919 = sin(t922);
	t1017 = t1015 * t919;
	t920 = cos(t922);
	t1016 = t1015 * t920;
	t1002 = t906 * t927;
	t1007 = t947 * t930 + t1002;
	t854 = t1007 * qJD(3) - (-t893 * t925 + t950) * t927 - t894 * t930;
	t863 = t877 * t926 - t897 * t929;
	t881 = t893 * t923 + t925 * t957;
	t845 = t863 * qJD(4) - t854 * t929 + t881 * t926;
	t844 = t865 * qJD(4) + t854 * t926 + t881 * t929;
	t937 = t949 - t986;
	t961 = t923 * t980;
	t941 = -t925 * t936 + t961;
	t879 = t941 * t927 - t930 * t937;
	t891 = t905 * qJD(1) + t937 * qJD(2);
	t892 = t906 * qJD(1) + t936 * qJD(2);
	t978 = qJD(1) * t993;
	t960 = t923 * t978;
	t849 = t879 * qJD(3) - t891 * t991 - t892 * t927 - t930 * t960;
	t899 = t923 * t936 + t925 * t980;
	t951 = t879 * t929 + t899 * t926;
	t1008 = -t951 * t921 + t849;
	t987 = t930 * t931;
	t990 = t927 * t928;
	t946 = t925 * t987 - t990;
	t1000 = t937 * t927;
	t999 = t919 * t921;
	t998 = t920 * t921;
	t997 = t921 * t929;
	t996 = t923 * t924;
	t995 = t923 * t926;
	t994 = t923 * t929;
	t989 = t927 * t931;
	t988 = t928 * t930;
	t985 = qJD(4) * t926;
	t984 = qJD(4) * t929;
	t983 = t928 * t996;
	t977 = qJD(2) * t996;
	t976 = t923 * t1004;
	t850 = -t892 * t930 + (t891 * t925 + t960) * t927 + (t941 * t930 + t1000) * qJD(3);
	t866 = -t879 * t926 + t899 * t929;
	t940 = -t891 * t923 + t925 * t978;
	t843 = t866 * qJD(4) + t850 * t929 + t940 * t926;
	t878 = -t930 * t961 + t936 * t991 - t1000;
	t974 = t878 * t921 + t843;
	t874 = t905 * t991 + t930 * t982 + t1002;
	t973 = -t874 * t921 - t845;
	t972 = -t1007 * t921 - t845;
	t888 = -t927 * t936 - t937 * t991;
	t858 = -t888 * qJD(3) + t891 * t930 + t892 * t992;
	t889 = -t930 * t936 + t937 * t992;
	t971 = t888 * t921 - t892 * t995 + t858 * t929 + (-t889 * t926 - t937 * t994) * qJD(4);
	t886 = -t905 * t927 + t906 * t991;
	t860 = -t886 * qJD(3) - t893 * t930 - t894 * t992;
	t887 = -t905 * t930 - t906 * t992;
	t970 = t886 * t921 + t894 * t995 + t860 * t929 + (-t887 * t926 + t906 * t994) * qJD(4);
	t943 = -t925 * t990 + t987;
	t955 = qJD(3) * t976;
	t870 = t930 * t955 + (t943 * qJD(2) + t946 * qJD(3)) * t924;
	t944 = t925 * t989 + t988;
	t896 = t944 * t924 + t927 * t976;
	t904 = t1004 * t925 - t931 * t996;
	t871 = -t896 * t926 + t904 * t929;
	t959 = t928 * t977;
	t856 = t871 * qJD(4) + t870 * t929 + t926 * t959;
	t895 = -t946 * t924 - t930 * t976;
	t967 = -t895 * t921 - t856;
	t868 = t889 * t929 - t937 * t995;
	t966 = t889 * qJD(3) - t868 * t921 + t891 * t927 - t892 * t991;
	t867 = t887 * t929 + t906 * t995;
	t965 = t887 * qJD(3) - t867 * t921 - t893 * t927 + t894 * t991;
	t945 = t925 * t988 + t989;
	t883 = (-t944 * qJD(2) - t945 * qJD(3)) * t924;
	t903 = t943 * t924;
	t958 = t931 * t977;
	t964 = t945 * t924 * t921 + t926 * t958 + t883 * t929 + (-t903 * t926 + t929 * t983) * qJD(4);
	t869 = t927 * t955 + (t945 * qJD(2) + t944 * qJD(3)) * t924;
	t872 = t896 * t929 + t904 * t926;
	t963 = t872 * t921 - t869;
	t890 = t903 * t929 + t926 * t983;
	t962 = -t890 * t921 + (t946 * qJD(2) + t943 * qJD(3)) * t924;
	t954 = t878 * t997 + t850;
	t953 = t874 * t997 - t854;
	t952 = t895 * t997 + t870;
	t935 = -t849 * t929 + t878 * t985 + t879 * t921;
	t934 = -t851 * t929 + t874 * t985 - t877 * t921;
	t933 = -t869 * t929 + t895 * t985 + t896 * t921;
	t855 = -t872 * qJD(4) - t870 * t926 + t929 * t959;
	t842 = t951 * qJD(4) + t850 * t926 - t940 * t929;
	t841 = t963 * t919 + t967 * t920;
	t840 = t967 * t919 - t963 * t920;
	t839 = t973 * t920 + t1017;
	t838 = t973 * t919 - t1016;
	t837 = t1008 * t919 + t974 * t920;
	t836 = t1008 * t920 - t974 * t919;
	t1 = [t972 * t920 + t1017, t966 * t919 + t971 * t920, t919 * t954 + t920 * t935, -t842 * t920 - t866 * t999, t836, t836; t837, t965 * t919 + t970 * t920, t919 * t953 + t920 * t934, t844 * t920 - t863 * t999, t838, t838; 0, t962 * t919 + t964 * t920, t919 * t952 + t920 * t933, t855 * t920 - t871 * t999, t840, t840; -t972 * t919 + t1016, -t971 * t919 + t966 * t920, -t919 * t935 + t920 * t954, t842 * t919 - t866 * t998, -t837, -t837; t836, -t970 * t919 + t965 * t920, -t919 * t934 + t920 * t953, -t844 * t919 - t863 * t998, t839, t839; 0, -t964 * t919 + t920 * t962, -t919 * t933 + t920 * t952, -t855 * t919 - t871 * t998, t841, t841; t844, qJD(4) * t868 + t858 * t926 + t892 * t994, -t849 * t926 - t878 * t984, t843, 0, 0; t842, qJD(4) * t867 + t860 * t926 - t894 * t994, -t851 * t926 - t874 * t984, t845, 0, 0; 0, qJD(4) * t890 + t883 * t926 - t929 * t958, -t869 * t926 - t895 * t984, t856, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end