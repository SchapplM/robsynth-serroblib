% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 21:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRR12_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR12_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:07:12
	% EndTime: 2019-12-29 21:07:12
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:07:18
	% EndTime: 2019-12-29 21:07:18
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:07:13
	% EndTime: 2019-12-29 21:07:13
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (27->13), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t158 = sin(qJ(2));
	t159 = sin(qJ(1));
	t171 = t158 * t159;
	t161 = cos(qJ(1));
	t170 = t158 * t161;
	t160 = cos(qJ(2));
	t169 = t159 * t160;
	t168 = t160 * t161;
	t156 = sin(pkin(5));
	t167 = qJD(1) * t156;
	t166 = qJD(2) * t156;
	t157 = cos(pkin(5));
	t165 = -t157 * t168 + t171;
	t164 = t157 * t169 + t170;
	t163 = t157 * t170 + t169;
	t162 = t157 * t171 - t168;
	t155 = t162 * qJD(1) + t165 * qJD(2);
	t154 = t164 * qJD(1) + t163 * qJD(2);
	t153 = t163 * qJD(1) + t164 * qJD(2);
	t152 = t165 * qJD(1) + t162 * qJD(2);
	t1 = [t155, t152, 0, 0, 0; -t153, -t154, 0, 0, 0; 0, -t158 * t166, 0, 0, 0; t154, t153, 0, 0, 0; t152, t155, 0, 0, 0; 0, -t160 * t166, 0, 0, 0; -t159 * t167, 0, 0, 0, 0; t161 * t167, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:07:21
	% EndTime: 2019-12-29 21:07:21
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (178->54), mult. (587->106), div. (0->0), fcn. (611->10), ass. (0->49)
	t364 = cos(pkin(5));
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
	t361 = sin(pkin(6));
	t365 = sin(qJ(3));
	t368 = cos(qJ(3));
	t375 = t364 * t391 + t389;
	t350 = t375 * qJD(1) + t353 * qJD(2);
	t363 = cos(pkin(6));
	t362 = sin(pkin(5));
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
	t1 = [t371, t349 * t397 + t348 * t368 + (t365 * t375 + t374 * t396) * qJD(3), t346, 0, 0; t347, -t351 * t397 - t350 * t368 + (t352 * t365 - t353 * t396) * qJD(3), t405, 0, 0; 0, (t377 * qJD(2) + t378 * qJD(3)) * t362, -t364 * t365 * t386 + (t378 * qJD(2) + t377 * qJD(3)) * t362, 0, 0; -t405, t349 * t396 - t348 * t365 + (t368 * t375 - t374 * t397) * qJD(3), -t347, 0, 0; t346, -t351 * t396 + t350 * t365 + (t352 * t368 + t353 * t397) * qJD(3), t371, 0, 0; 0, (t379 * qJD(2) + t376 * qJD(3)) * t362, -t364 * t382 + (t376 * qJD(2) + t379 * qJD(3)) * t362, 0, 0; -t350 * t361 - t363 * t384, -t349 * t361, 0, 0, 0; -t348 * t361 + t363 * t383, t351 * t361, 0, 0, 0; 0, qJD(2) * t369 * t399, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:07:21
	% EndTime: 2019-12-29 21:07:22
	% DurationCPUTime: 1.67s
	% Computational Cost: add. (493->104), mult. (1577->199), div. (0->0), fcn. (1713->12), ass. (0->85)
	t613 = sin(qJ(2));
	t614 = sin(qJ(1));
	t617 = cos(qJ(2));
	t618 = cos(qJ(1));
	t661 = cos(pkin(5));
	t635 = t618 * t661;
	t598 = t613 * t635 + t614 * t617;
	t636 = t614 * t661;
	t619 = t618 * t613 + t617 * t636;
	t588 = t619 * qJD(1) + t598 * qJD(2);
	t629 = t613 * t636;
	t645 = t618 * t617;
	t647 = t614 * t613;
	t589 = -qJD(1) * t629 - qJD(2) * t647 + (qJD(2) * t661 + qJD(1)) * t645;
	t610 = cos(pkin(6));
	t612 = sin(qJ(3));
	t616 = cos(qJ(3));
	t608 = sin(pkin(6));
	t609 = sin(pkin(5));
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
	t1 = [t671, -t587 * t656 + t570 * t615 + (-t585 * t611 - t620 * t655) * qJD(4), -t565 * t615 - t577 * t643, t563, 0; t564, t589 * t656 + t571 * t615 + (-t584 * t611 + t598 * t655) * qJD(4), t567 * t615 + t643 * t662, t672, 0; 0, t611 * t630 + t582 * t615 + (-t595 * t611 + t615 * t641) * qJD(4), t572 * t615 - t590 * t643, t615 * t631 - t573 * t611 + (-t591 * t615 - t596 * t611) * qJD(4), 0; -t672, -t587 * t655 - t570 * t611 + (-t585 * t615 + t620 * t656) * qJD(4), t565 * t611 - t577 * t642, -t564, 0; t563, t589 * t655 - t571 * t611 + (-t584 * t615 - t598 * t656) * qJD(4), -t567 * t611 + t642 * t662, t671, 0; 0, t615 * t630 - t582 * t611 + (-t595 * t615 - t611 * t641) * qJD(4), -t572 * t611 - t590 * t642, -t611 * t631 - t573 * t615 + (t591 * t611 - t596 * t615) * qJD(4), 0; t567, t585 * qJD(3) + t586 * t612 - t587 * t651, t566, 0, 0; t565, t584 * qJD(3) - t588 * t612 + t589 * t651, -t569, 0, 0; 0, (t625 * qJD(2) + t622 * qJD(3)) * t609, t573, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 21:07:28
	% EndTime: 2019-12-29 21:07:32
	% DurationCPUTime: 3.95s
	% Computational Cost: add. (1250->165), mult. (3901->300), div. (0->0), fcn. (4371->14), ass. (0->125)
	t858 = sin(qJ(2));
	t862 = cos(qJ(2));
	t863 = cos(qJ(1));
	t920 = cos(pkin(5));
	t892 = t863 * t920;
	t921 = sin(qJ(1));
	t839 = t858 * t892 + t921 * t862;
	t885 = t920 * t921;
	t867 = t863 * t858 + t862 * t885;
	t826 = t867 * qJD(1) + t839 * qJD(2);
	t879 = t858 * t885;
	t895 = t921 * t858;
	t905 = t863 * t862;
	t827 = -qJD(1) * t879 - qJD(2) * t895 + (qJD(2) * t920 + qJD(1)) * t905;
	t854 = cos(pkin(6));
	t857 = sin(qJ(3));
	t861 = cos(qJ(3));
	t852 = sin(pkin(6));
	t853 = sin(pkin(5));
	t896 = t853 * t921;
	t886 = qJD(1) * t896;
	t880 = t852 * t886;
	t838 = -t862 * t892 + t895;
	t912 = t853 * t863;
	t898 = t852 * t912;
	t877 = t838 * t854 + t898;
	t918 = t839 * t857;
	t922 = t877 * t861 + t918;
	t787 = t922 * qJD(3) - (-t826 * t854 + t880) * t857 - t827 * t861;
	t917 = t839 * t861;
	t810 = t877 * t857 - t917;
	t830 = -t838 * t852 + t854 * t912;
	t856 = sin(qJ(4));
	t860 = cos(qJ(4));
	t795 = t810 * t856 - t830 * t860;
	t814 = t826 * t852 + t854 * t886;
	t779 = -t795 * qJD(4) + t787 * t860 - t814 * t856;
	t910 = t854 * t861;
	t911 = t854 * t857;
	t786 = (t838 * t911 - t917) * qJD(3) - t826 * t910 + t861 * t880 - (-qJD(3) * t898 + t827) * t857;
	t855 = sin(qJ(5));
	t859 = cos(qJ(5));
	t935 = -t779 * t855 + t786 * t859;
	t934 = t779 * t859 + t786 * t855;
	t797 = t810 * t860 + t830 * t856;
	t933 = t797 * t855;
	t932 = t797 * t859;
	t777 = t797 * qJD(4) + t787 * t856 + t814 * t860;
	t906 = t861 * t862;
	t909 = t857 * t858;
	t876 = t854 * t906 - t909;
	t868 = t879 - t905;
	t916 = t868 * t857;
	t915 = t852 * t853;
	t914 = t852 * t856;
	t913 = t852 * t860;
	t908 = t857 * t862;
	t907 = t858 * t861;
	t904 = qJD(4) * t856;
	t903 = qJD(4) * t860;
	t902 = qJD(5) * t855;
	t901 = qJD(5) * t859;
	t900 = qJD(5) * t860;
	t899 = t858 * t915;
	t894 = qJD(1) * t912;
	t893 = qJD(2) * t915;
	t891 = t920 * t852;
	t890 = t852 * t896;
	t889 = t852 * t894;
	t888 = t858 * t893;
	t887 = t862 * t893;
	t884 = qJD(3) * t891;
	t824 = t838 * qJD(1) + t868 * qJD(2);
	t825 = t839 * qJD(1) + t867 * qJD(2);
	t872 = -t854 * t867 + t890;
	t783 = -t825 * t861 + (t824 * t854 + t889) * t857 + (t872 * t861 + t916) * qJD(3);
	t811 = -t861 * t890 + t867 * t910 - t916;
	t883 = t811 * t900 + t783;
	t807 = t838 * t910 + t861 * t898 + t918;
	t882 = t807 * t900 - t787;
	t873 = -t854 * t909 + t906;
	t803 = t861 * t884 + (t873 * qJD(2) + t876 * qJD(3)) * t853;
	t828 = -t876 * t853 - t861 * t891;
	t881 = t828 * t900 + t803;
	t812 = t872 * t857 - t861 * t868;
	t832 = t852 * t867 + t854 * t896;
	t799 = t812 * t860 + t832 * t856;
	t798 = -t812 * t856 + t832 * t860;
	t874 = t854 * t908 + t907;
	t829 = t874 * t853 + t857 * t891;
	t837 = t920 * t854 - t862 * t915;
	t805 = t829 * t860 + t837 * t856;
	t804 = -t829 * t856 + t837 * t860;
	t820 = -t838 * t861 - t839 * t911;
	t800 = t820 * t860 + t839 * t914;
	t822 = -t861 * t867 + t868 * t911;
	t801 = t822 * t860 - t868 * t914;
	t819 = -t838 * t857 + t839 * t910;
	t821 = -t857 * t867 - t868 * t910;
	t875 = t854 * t907 + t908;
	t836 = t873 * t853;
	t823 = t836 * t860 + t856 * t899;
	t871 = -t824 * t852 + t854 * t894;
	t782 = t812 * qJD(3) - t824 * t910 - t825 * t857 - t861 * t889;
	t866 = qJD(5) * t812 - t782 * t860 + t811 * t904;
	t865 = -qJD(5) * t810 + t786 * t860 + t807 * t904;
	t802 = t857 * t884 + (t875 * qJD(2) + t874 * qJD(3)) * t853;
	t864 = qJD(5) * t829 - t802 * t860 + t828 * t904;
	t835 = t875 * t853;
	t816 = (-t874 * qJD(2) - t875 * qJD(3)) * t853;
	t815 = (t876 * qJD(2) + t873 * qJD(3)) * t853;
	t794 = t856 * t887 + t816 * t860 + (-t836 * t856 + t860 * t899) * qJD(4);
	t793 = -t819 * qJD(3) - t826 * t861 - t827 * t911;
	t792 = t820 * qJD(3) - t826 * t857 + t827 * t910;
	t791 = -t821 * qJD(3) + t824 * t861 + t825 * t911;
	t790 = t822 * qJD(3) + t824 * t857 - t825 * t910;
	t789 = t804 * qJD(4) + t803 * t860 + t856 * t888;
	t788 = -t805 * qJD(4) - t803 * t856 + t860 * t888;
	t781 = t827 * t914 + t793 * t860 + (-t820 * t856 + t839 * t913) * qJD(4);
	t780 = -t825 * t914 + t791 * t860 + (-t822 * t856 - t868 * t913) * qJD(4);
	t776 = t798 * qJD(4) + t783 * t860 + t871 * t856;
	t775 = t799 * qJD(4) + t783 * t856 - t871 * t860;
	t774 = t776 * t859 + t782 * t855 + (-t799 * t855 + t811 * t859) * qJD(5);
	t773 = -t776 * t855 + t782 * t859 + (-t799 * t859 - t811 * t855) * qJD(5);
	t1 = [(-t859 * t922 - t933) * qJD(5) + t934, t780 * t859 + t790 * t855 + (-t801 * t855 + t821 * t859) * qJD(5), t883 * t855 + t866 * t859, -t775 * t859 - t798 * t902, t773; t774, t781 * t859 + t792 * t855 + (-t800 * t855 + t819 * t859) * qJD(5), t882 * t855 + t865 * t859, t777 * t859 - t795 * t902, (-t807 * t855 + t932) * qJD(5) - t935; 0, t794 * t859 + t815 * t855 + (-t823 * t855 + t835 * t859) * qJD(5), t881 * t855 + t864 * t859, t788 * t859 - t804 * t902, -t789 * t855 + t802 * t859 + (-t805 * t859 - t828 * t855) * qJD(5); (t855 * t922 - t932) * qJD(5) + t935, -t780 * t855 + t790 * t859 + (-t801 * t859 - t821 * t855) * qJD(5), -t866 * t855 + t883 * t859, t775 * t855 - t798 * t901, -t774; t773, -t781 * t855 + t792 * t859 + (-t800 * t859 - t819 * t855) * qJD(5), -t865 * t855 + t882 * t859, -t777 * t855 - t795 * t901, (-t807 * t859 - t933) * qJD(5) + t934; 0, -t794 * t855 + t815 * t859 + (-t823 * t859 - t835 * t855) * qJD(5), -t864 * t855 + t881 * t859, -t788 * t855 - t804 * t901, -t789 * t859 - t802 * t855 + (t805 * t855 - t828 * t859) * qJD(5); t777, t801 * qJD(4) + t791 * t856 + t825 * t913, -t782 * t856 - t811 * t903, t776, 0; t775, t800 * qJD(4) + t793 * t856 - t827 * t913, t786 * t856 - t807 * t903, -t779, 0; 0, t823 * qJD(4) + t816 * t856 - t860 * t887, -t802 * t856 - t828 * t903, t789, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end