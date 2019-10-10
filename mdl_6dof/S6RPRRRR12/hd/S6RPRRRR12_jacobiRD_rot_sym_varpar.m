% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR12
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:16
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRRR12_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR12_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobiRD_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
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
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->5), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t128 = cos(pkin(6));
	t129 = sin(qJ(1));
	t133 = t128 * t129;
	t130 = cos(qJ(1));
	t132 = t128 * t130;
	t131 = qJD(1) * sin(pkin(6));
	t127 = cos(pkin(14));
	t125 = sin(pkin(14));
	t1 = [(t125 * t133 - t127 * t130) * qJD(1), 0, 0, 0, 0, 0; (-t125 * t132 - t127 * t129) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; (t125 * t130 + t127 * t133) * qJD(1), 0, 0, 0, 0, 0; (t125 * t129 - t127 * t132) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t129 * t131, 0, 0, 0, 0, 0; t130 * t131, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:47
	% EndTime: 2019-10-10 09:15:47
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (80->29), mult. (294->61), div. (0->0), fcn. (310->10), ass. (0->38)
	t298 = cos(pkin(6));
	t296 = cos(pkin(14));
	t302 = cos(qJ(1));
	t314 = t302 * t296;
	t293 = sin(pkin(14));
	t300 = sin(qJ(1));
	t317 = t300 * t293;
	t286 = -t298 * t314 + t317;
	t294 = sin(pkin(7));
	t295 = sin(pkin(6));
	t297 = cos(pkin(7));
	t327 = t294 * t295 * t302 + t286 * t297;
	t311 = t298 * t317;
	t313 = qJD(1) * t302;
	t285 = -qJD(1) * t311 + t296 * t313;
	t315 = t302 * t293;
	t316 = t300 * t296;
	t287 = t298 * t315 + t316;
	t299 = sin(qJ(3));
	t301 = cos(qJ(3));
	t306 = t298 * t316 + t315;
	t284 = t306 * qJD(1);
	t319 = t295 * t300;
	t310 = qJD(1) * t319;
	t304 = -t284 * t297 + t294 * t310;
	t326 = (-t287 * t301 + t327 * t299) * qJD(3) - t285 * t299 + t304 * t301;
	t320 = t294 * t298;
	t318 = t296 * t297;
	t309 = t295 * t313;
	t307 = t294 * t319 - t297 * t306;
	t282 = t286 * qJD(1);
	t305 = t282 * t297 + t294 * t309;
	t303 = (qJD(3) * t287 - t304) * t299 + (t327 * qJD(3) - t285) * t301;
	t289 = -t311 + t314;
	t283 = t287 * qJD(1);
	t281 = -t283 * t301 + t305 * t299 + (-t289 * t299 + t307 * t301) * qJD(3);
	t280 = t283 * t299 + t305 * t301 + (-t289 * t301 - t307 * t299) * qJD(3);
	t1 = [t303, 0, t280, 0, 0, 0; t281, 0, t326, 0, 0, 0; 0, 0, (-t299 * t320 + (-t293 * t301 - t299 * t318) * t295) * qJD(3), 0, 0, 0; -t326, 0, -t281, 0, 0, 0; t280, 0, t303, 0, 0, 0; 0, 0, (-t301 * t320 + (t293 * t299 - t301 * t318) * t295) * qJD(3), 0, 0, 0; -t284 * t294 - t297 * t310, 0, 0, 0, 0, 0; -t282 * t294 + t297 * t309, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:51
	% EndTime: 2019-10-10 09:15:52
	% DurationCPUTime: 0.90s
	% Computational Cost: add. (491->89), mult. (1655->168), div. (0->0), fcn. (1841->14), ass. (0->75)
	t605 = cos(pkin(6));
	t602 = cos(pkin(14));
	t611 = cos(qJ(1));
	t636 = t611 * t602;
	t598 = sin(pkin(14));
	t608 = sin(qJ(1));
	t639 = t608 * t598;
	t588 = -t605 * t636 + t639;
	t607 = sin(qJ(3));
	t600 = sin(pkin(7));
	t601 = sin(pkin(6));
	t644 = t601 * t611;
	t632 = t600 * t644;
	t604 = cos(pkin(7));
	t640 = t604 * t607;
	t637 = t611 * t598;
	t638 = t608 * t602;
	t589 = t605 * t637 + t638;
	t610 = cos(qJ(3));
	t648 = t589 * t610;
	t568 = t588 * t640 + t607 * t632 - t648;
	t606 = sin(qJ(4));
	t609 = cos(qJ(4));
	t629 = t605 * t639;
	t635 = qJD(1) * t611;
	t586 = -qJD(1) * t629 + t602 * t635;
	t614 = t605 * t638 + t637;
	t585 = t614 * qJD(1);
	t645 = t601 * t608;
	t628 = qJD(1) * t645;
	t625 = t600 * t628;
	t612 = -t585 * t604 + t625;
	t650 = t588 * t604;
	t617 = t632 + t650;
	t561 = -t586 * t607 + t612 * t610 + (t617 * t607 - t648) * qJD(3);
	t575 = t585 * t600 + t604 * t628;
	t599 = sin(pkin(8));
	t603 = cos(pkin(8));
	t623 = t561 * t603 + t575 * t599;
	t652 = t589 * t607 + t617 * t610;
	t659 = t652 * t603 + (-t588 * t600 + t604 * t644) * t599;
	t664 = -t623 * t606 + (-t568 * t606 + t609 * t659) * qJD(4);
	t663 = t623 * t609 + (t568 * t609 + t606 * t659) * qJD(4);
	t646 = t600 * t605;
	t577 = (t598 * t610 + t602 * t640) * t601 + t607 * t646;
	t591 = -t629 + t636;
	t616 = t600 * t645 - t604 * t614;
	t570 = t591 * t610 + t616 * t607;
	t643 = t602 * t604;
	t642 = t603 * t606;
	t641 = t603 * t609;
	t634 = qJD(3) * t607;
	t633 = qJD(3) * t610;
	t630 = t610 * t646;
	t627 = t601 * t635;
	t626 = t601 * t633;
	t584 = t589 * qJD(1);
	t583 = t588 * qJD(1);
	t613 = t583 * t604 + t600 * t627;
	t559 = -t570 * qJD(3) + t584 * t607 + t613 * t610;
	t573 = -t583 * t600 + t604 * t627;
	t624 = t559 * t603 + t573 * t599;
	t569 = -t591 * t607 + t616 * t610;
	t619 = t569 * t603 + (t600 * t614 + t604 * t645) * t599;
	t576 = t630 + (-t598 * t607 + t610 * t643) * t601;
	t618 = -t576 * t603 - (-t601 * t602 * t600 + t605 * t604) * t599;
	t593 = t611 * t600 * t626;
	t572 = t577 * qJD(3);
	t571 = t601 * t598 * t634 - qJD(3) * t630 - t626 * t643;
	t564 = t593 + (qJD(3) * t650 - t586) * t610 + (qJD(3) * t589 - t612) * t607;
	t562 = t607 * t625 + t586 * t610 - t589 * t634 - t593 + (-t585 * t607 - t588 * t633) * t604;
	t560 = -t591 * t634 + t613 * t607 + (t616 * qJD(3) - t584) * t610;
	t558 = t560 * t609 + t624 * t606 + (-t570 * t606 + t619 * t609) * qJD(4);
	t557 = -t560 * t606 + t624 * t609 + (-t570 * t609 - t619 * t606) * qJD(4);
	t1 = [t564 * t609 + t664, 0, -t560 * t642 + t559 * t609 + (-t569 * t606 - t570 * t641) * qJD(4), t557, 0, 0; t558, 0, -t562 * t642 + t561 * t609 + (t568 * t641 + t606 * t652) * qJD(4), -t562 * t606 + t663, 0, 0; 0, 0, t571 * t642 - t572 * t609 + (-t576 * t606 - t577 * t641) * qJD(4), -t572 * t641 + t571 * t606 + (-t577 * t609 + t618 * t606) * qJD(4), 0, 0; -t564 * t606 - t663, 0, -t560 * t641 - t559 * t606 + (-t569 * t609 + t570 * t642) * qJD(4), -t558, 0, 0; t557, 0, -t562 * t641 - t561 * t606 + (-t568 * t642 + t609 * t652) * qJD(4), -t562 * t609 + t664, 0, 0; 0, 0, t571 * t641 + t572 * t606 + (-t576 * t609 + t577 * t642) * qJD(4), t572 * t642 + t571 * t609 + (t577 * t606 + t618 * t609) * qJD(4), 0, 0; t561 * t599 - t575 * t603, 0, t560 * t599, 0, 0, 0; -t559 * t599 + t573 * t603, 0, t562 * t599, 0, 0, 0; 0, 0, -t571 * t599, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:58
	% EndTime: 2019-10-10 09:16:00
	% DurationCPUTime: 2.25s
	% Computational Cost: add. (1328->136), mult. (4334->250), div. (0->0), fcn. (4993->16), ass. (0->111)
	t855 = cos(pkin(14));
	t858 = cos(pkin(6));
	t851 = sin(pkin(14));
	t862 = sin(qJ(1));
	t898 = t862 * t851;
	t886 = t858 * t898;
	t866 = cos(qJ(1));
	t894 = qJD(1) * t866;
	t840 = -qJD(1) * t886 + t855 * t894;
	t861 = sin(qJ(3));
	t865 = cos(qJ(3));
	t896 = t866 * t851;
	t897 = t862 * t855;
	t872 = t858 * t897 + t896;
	t839 = t872 * qJD(1);
	t857 = cos(pkin(7));
	t853 = sin(pkin(7));
	t854 = sin(pkin(6));
	t904 = t854 * t862;
	t885 = qJD(1) * t904;
	t882 = t853 * t885;
	t869 = -t839 * t857 + t882;
	t903 = t854 * t866;
	t889 = t853 * t903;
	t895 = t866 * t855;
	t842 = -t858 * t895 + t898;
	t910 = t842 * t857;
	t875 = t889 + t910;
	t843 = t858 * t896 + t897;
	t908 = t843 * t865;
	t805 = -t840 * t861 + t869 * t865 + (t875 * t861 - t908) * qJD(3);
	t829 = t839 * t853 + t857 * t885;
	t852 = sin(pkin(8));
	t856 = cos(pkin(8));
	t795 = t805 * t852 - t829 * t856;
	t899 = t857 * t861;
	t821 = t842 * t899 + t861 * t889 - t908;
	t860 = sin(qJ(4));
	t864 = cos(qJ(4));
	t820 = t843 * t861 + t875 * t865;
	t832 = -t842 * t853 + t857 * t903;
	t878 = t820 * t856 + t832 * t852;
	t799 = t821 * t864 + t878 * t860;
	t813 = t820 * t852 - t832 * t856;
	t859 = sin(qJ(5));
	t863 = cos(qJ(5));
	t934 = (t799 * t863 - t813 * t859) * qJD(5) - t795 * t863;
	t881 = t805 * t856 + t829 * t852;
	t932 = t799 * qJD(4) + t881 * t864;
	t923 = -t821 * t860 + t878 * t864;
	t931 = -t923 * qJD(4) + t881 * t860;
	t930 = t795 * t859 + (-t799 * t859 - t813 * t863) * qJD(5);
	t905 = t853 * t858;
	t831 = (t851 * t865 + t855 * t899) * t854 + t861 * t905;
	t887 = t865 * t905;
	t902 = t855 * t857;
	t830 = t887 + (-t851 * t861 + t865 * t902) * t854;
	t841 = -t854 * t855 * t853 + t858 * t857;
	t876 = t830 * t856 + t841 * t852;
	t812 = t831 * t864 + t876 * t860;
	t871 = t886 - t895;
	t874 = t853 * t904 - t857 * t872;
	t823 = t874 * t861 - t871 * t865;
	t837 = t842 * qJD(1);
	t884 = t854 * t894;
	t827 = -t837 * t853 + t857 * t884;
	t913 = t827 * t852;
	t907 = t852 * t859;
	t906 = t852 * t863;
	t901 = t856 * t860;
	t900 = t856 * t864;
	t893 = qJD(3) * t861;
	t892 = qJD(3) * t865;
	t891 = qJD(5) * t859;
	t890 = qJD(5) * t863;
	t883 = t854 * t892;
	t868 = t871 * t861;
	t822 = t874 * t865 + t868;
	t834 = t853 * t872 + t857 * t904;
	t877 = t822 * t856 + t834 * t852;
	t809 = -t820 * t864 + t821 * t901;
	t810 = t822 * t864 - t823 * t901;
	t816 = t830 * t864 - t831 * t901;
	t870 = t837 * t857 + t853 * t884;
	t800 = -t823 * t860 + t877 * t864;
	t801 = t823 * t864 + t877 * t860;
	t811 = -t831 * t860 + t876 * t864;
	t846 = t866 * t853 * t883;
	t838 = t843 * qJD(1);
	t826 = t831 * qJD(3);
	t825 = t854 * t851 * t893 - qJD(3) * t887 - t883 * t902;
	t817 = -t830 * t852 + t841 * t856;
	t815 = -t822 * t852 + t834 * t856;
	t808 = t846 + (qJD(3) * t910 - t840) * t865 + (qJD(3) * t843 - t869) * t861;
	t806 = t861 * t882 + t840 * t865 - t843 * t893 - t846 + (-t839 * t861 - t842 * t892) * t857;
	t804 = qJD(3) * t868 + t870 * t861 + (t874 * qJD(3) - t838) * t865;
	t803 = -t823 * qJD(3) + t838 * t861 + t870 * t865;
	t796 = t825 * t901 - t826 * t864 + (-t830 * t860 - t831 * t900) * qJD(4);
	t793 = -t803 * t852 + t827 * t856;
	t792 = t811 * qJD(4) - t825 * t864 - t826 * t901;
	t791 = -t812 * qJD(4) + t825 * t860 - t826 * t900;
	t790 = -t806 * t901 + t805 * t864 + (t820 * t860 + t821 * t900) * qJD(4);
	t789 = -t804 * t901 + t803 * t864 + (-t822 * t860 - t823 * t900) * qJD(4);
	t788 = t808 * t864 - t931;
	t787 = t806 * t864 + t931;
	t786 = -t806 * t860 + t932;
	t785 = t804 * t864 + (t803 * t856 + t913) * t860 + t800 * qJD(4);
	t784 = t801 * qJD(4) - t803 * t900 + t804 * t860 - t864 * t913;
	t783 = t785 * t863 + t793 * t859 + (-t801 * t859 + t815 * t863) * qJD(5);
	t782 = -t785 * t859 + t793 * t863 + (-t801 * t863 - t815 * t859) * qJD(5);
	t1 = [t788 * t863 + t930, 0, t804 * t907 + t789 * t863 + (-t810 * t859 + t823 * t906) * qJD(5), -t784 * t863 - t800 * t891, t782, 0; t783, 0, t806 * t907 + t790 * t863 + (-t809 * t859 - t821 * t906) * qJD(5), t786 * t863 + t891 * t923, -t787 * t859 + t934, 0; 0, 0, -t825 * t907 + t796 * t863 + (-t816 * t859 + t831 * t906) * qJD(5), t791 * t863 - t811 * t891, t826 * t906 - t792 * t859 + (-t812 * t863 - t817 * t859) * qJD(5), 0; -t788 * t859 - t934, 0, t804 * t906 - t789 * t859 + (-t810 * t863 - t823 * t907) * qJD(5), t784 * t859 - t800 * t890, -t783, 0; t782, 0, t806 * t906 - t790 * t859 + (-t809 * t863 + t821 * t907) * qJD(5), -t786 * t859 + t890 * t923, -t787 * t863 + t930, 0; 0, 0, -t825 * t906 - t796 * t859 + (-t816 * t863 - t831 * t907) * qJD(5), -t791 * t859 - t811 * t890, -t826 * t907 - t792 * t863 + (t812 * t859 - t817 * t863) * qJD(5), 0; t808 * t860 + t932, 0, t810 * qJD(4) + t803 * t860 + t804 * t900, t785, 0, 0; t784, 0, t809 * qJD(4) + t805 * t860 + t806 * t900, t787, 0, 0; 0, 0, t816 * qJD(4) - t825 * t900 - t826 * t860, t792, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:16:08
	% EndTime: 2019-10-10 09:16:12
	% DurationCPUTime: 4.43s
	% Computational Cost: add. (3315->201), mult. (10558->351), div. (0->0), fcn. (12469->18), ass. (0->155)
	t1113 = cos(pkin(6));
	t1110 = cos(pkin(14));
	t1123 = cos(qJ(1));
	t1163 = t1123 * t1110;
	t1106 = sin(pkin(14));
	t1118 = sin(qJ(1));
	t1166 = t1118 * t1106;
	t1097 = -t1113 * t1163 + t1166;
	t1117 = sin(qJ(3));
	t1108 = sin(pkin(7));
	t1109 = sin(pkin(6));
	t1170 = t1109 * t1123;
	t1162 = t1108 * t1170;
	t1112 = cos(pkin(7));
	t1167 = t1112 * t1117;
	t1164 = t1123 * t1106;
	t1165 = t1118 * t1110;
	t1098 = t1113 * t1164 + t1165;
	t1122 = cos(qJ(3));
	t1176 = t1098 * t1122;
	t1075 = t1097 * t1167 + t1117 * t1162 - t1176;
	t1116 = sin(qJ(4));
	t1121 = cos(qJ(4));
	t1178 = t1097 * t1112;
	t1148 = t1162 + t1178;
	t1177 = t1098 * t1117;
	t1074 = t1148 * t1122 + t1177;
	t1111 = cos(pkin(8));
	t1155 = t1112 * t1170;
	t1089 = -t1097 * t1108 + t1155;
	t1107 = sin(pkin(8));
	t1181 = t1089 * t1107;
	t1149 = t1074 * t1111 + t1181;
	t1046 = t1075 * t1121 + t1149 * t1116;
	t1064 = t1074 * t1107 - t1089 * t1111;
	t1115 = sin(qJ(5));
	t1120 = cos(qJ(5));
	t1029 = t1046 * t1120 - t1064 * t1115;
	t1114 = sin(qJ(6));
	t1212 = t1029 * t1114;
	t1119 = cos(qJ(6));
	t1211 = t1029 * t1119;
	t1143 = t1113 * t1166 - t1163;
	t1096 = t1143 * qJD(1);
	t1144 = t1113 * t1165 + t1164;
	t1095 = t1144 * qJD(1);
	t1171 = t1109 * t1118;
	t1158 = t1108 * t1171;
	t1154 = qJD(1) * t1158;
	t1138 = -t1095 * t1112 + t1154;
	t1053 = (t1148 * t1117 - t1176) * qJD(3) + t1096 * t1117 + t1138 * t1122;
	t1157 = t1112 * t1171;
	t1139 = qJD(1) * t1157 + t1095 * t1108;
	t1038 = t1053 * t1107 - t1139 * t1111;
	t1210 = t1029 * qJD(5) - t1038 * t1120;
	t1027 = t1046 * t1115 + t1064 * t1120;
	t1209 = t1027 * qJD(5) - t1038 * t1115;
	t1207 = t1046 * qJD(4);
	t1185 = t1075 * t1116;
	t1195 = t1149 * t1121 - t1185;
	t1131 = t1139 * t1107;
	t1200 = t1053 * t1111 + t1131;
	t1204 = -t1195 * qJD(4) + t1200 * t1116;
	t1172 = t1108 * t1113;
	t1087 = (t1110 * t1112 * t1122 - t1106 * t1117) * t1109 + t1122 * t1172;
	t1088 = (t1106 * t1122 + t1110 * t1167) * t1109 + t1117 * t1172;
	t1146 = -t1112 * t1144 + t1158;
	t1077 = t1146 * t1117 - t1143 * t1122;
	t1192 = qJD(3) * t1122;
	t1191 = qJD(5) * t1115;
	t1190 = qJD(6) * t1114;
	t1189 = qJD(6) * t1119;
	t1168 = t1111 * t1121;
	t1043 = t1074 * t1168 + t1121 * t1181 - t1185;
	t1188 = t1043 * t1120;
	t1135 = t1143 * t1117;
	t1076 = t1146 * t1122 + t1135;
	t1147 = t1108 * t1144 + t1157;
	t1136 = t1147 * t1107;
	t1183 = t1077 * t1116;
	t1047 = -t1076 * t1168 - t1121 * t1136 + t1183;
	t1187 = t1047 * t1120;
	t1142 = -t1109 * t1110 * t1108 + t1113 * t1112;
	t1133 = t1142 * t1107;
	t1182 = t1088 * t1116;
	t1062 = -t1087 * t1168 - t1121 * t1133 + t1182;
	t1186 = t1062 * t1120;
	t1174 = t1107 * t1115;
	t1173 = t1107 * t1120;
	t1169 = t1111 * t1116;
	t1094 = t1098 * qJD(1);
	t1093 = t1097 * qJD(1);
	t1140 = qJD(1) * t1162 + t1093 * t1112;
	t1051 = -t1077 * qJD(3) + t1094 * t1117 + t1140 * t1122;
	t1052 = qJD(3) * t1135 + t1140 * t1117 + (t1146 * qJD(3) - t1094) * t1122;
	t1126 = t1076 * t1111 + t1136;
	t1141 = qJD(1) * t1155 - t1093 * t1108;
	t1132 = t1141 * t1107;
	t1015 = t1052 * t1121 + (t1051 * t1111 + t1132) * t1116 + (t1126 * t1121 - t1183) * qJD(4);
	t1153 = qJD(6) * t1187 + t1015;
	t1101 = t1162 * t1192;
	t1054 = t1117 * t1154 - t1096 * t1122 - qJD(3) * t1177 - t1101 + (-t1095 * t1117 - t1097 * t1192) * t1112;
	t1017 = t1054 * t1121 + t1204;
	t1152 = qJD(6) * t1188 + t1017;
	t1084 = t1087 * qJD(3);
	t1085 = t1088 * qJD(3);
	t1125 = t1087 * t1111 + t1133;
	t1035 = -t1085 * t1169 + t1084 * t1121 + (t1125 * t1121 - t1182) * qJD(4);
	t1151 = qJD(6) * t1186 + t1035;
	t1048 = t1077 * t1121 + t1126 * t1116;
	t1066 = -t1076 * t1107 + t1147 * t1111;
	t1031 = t1048 * t1120 + t1066 * t1115;
	t1030 = -t1048 * t1115 + t1066 * t1120;
	t1063 = t1088 * t1121 + t1125 * t1116;
	t1071 = -t1087 * t1107 + t1142 * t1111;
	t1040 = t1063 * t1120 + t1071 * t1115;
	t1039 = -t1063 * t1115 + t1071 * t1120;
	t1058 = -t1074 * t1121 + t1075 * t1169;
	t1032 = t1058 * t1120 - t1075 * t1174;
	t1060 = t1076 * t1121 - t1077 * t1169;
	t1033 = t1060 * t1120 + t1077 * t1174;
	t1068 = t1087 * t1121 - t1088 * t1169;
	t1061 = t1068 * t1120 + t1088 * t1174;
	t1057 = -t1074 * t1116 - t1075 * t1168;
	t1059 = t1076 * t1116 + t1077 * t1168;
	t1067 = t1087 * t1116 + t1088 * t1168;
	t1014 = t1048 * qJD(4) - t1051 * t1168 + t1052 * t1116 - t1121 * t1132;
	t1130 = qJD(6) * t1048 - t1014 * t1120 + t1047 * t1191;
	t1016 = -t1053 * t1168 + t1054 * t1116 - t1121 * t1131 - t1207;
	t1129 = -qJD(6) * t1046 - t1016 * t1120 + t1043 * t1191;
	t1034 = t1063 * qJD(4) + t1084 * t1116 + t1085 * t1168;
	t1128 = qJD(6) * t1063 - t1034 * t1120 + t1062 * t1191;
	t1124 = -t1051 * t1107 + t1141 * t1111;
	t1056 = t1101 + (qJD(3) * t1178 + t1096) * t1122 + (qJD(3) * t1098 - t1138) * t1117;
	t1042 = -t1067 * qJD(4) - t1084 * t1169 - t1085 * t1121;
	t1041 = t1068 * qJD(4) + t1084 * t1168 - t1085 * t1116;
	t1026 = t1084 * t1174 + t1042 * t1120 + (-t1068 * t1115 + t1088 * t1173) * qJD(5);
	t1025 = -t1057 * qJD(4) + t1053 * t1121 - t1054 * t1169;
	t1024 = t1058 * qJD(4) + t1053 * t1116 + t1054 * t1168;
	t1023 = -t1059 * qJD(4) + t1051 * t1121 - t1052 * t1169;
	t1022 = t1060 * qJD(4) + t1051 * t1116 + t1052 * t1168;
	t1021 = t1039 * qJD(5) + t1035 * t1120 + t1085 * t1174;
	t1020 = -t1040 * qJD(5) - t1035 * t1115 + t1085 * t1173;
	t1019 = t1056 * t1121 - t1204;
	t1018 = t1056 * t1116 + t1121 * t1200 + t1207;
	t1013 = t1054 * t1174 + t1025 * t1120 + (-t1058 * t1115 - t1075 * t1173) * qJD(5);
	t1012 = t1052 * t1174 + t1023 * t1120 + (-t1060 * t1115 + t1077 * t1173) * qJD(5);
	t1011 = t1019 * t1120 - t1209;
	t1010 = t1017 * t1120 + t1209;
	t1009 = -t1017 * t1115 + t1210;
	t1008 = t1030 * qJD(5) + t1015 * t1120 + t1124 * t1115;
	t1007 = t1031 * qJD(5) + t1015 * t1115 - t1124 * t1120;
	t1006 = t1008 * t1119 + t1014 * t1114 + (-t1031 * t1114 + t1047 * t1119) * qJD(6);
	t1005 = -t1008 * t1114 + t1014 * t1119 + (-t1031 * t1119 - t1047 * t1114) * qJD(6);
	t1 = [t1011 * t1119 + t1018 * t1114 + (-t1119 * t1195 - t1212) * qJD(6), 0, t1012 * t1119 + t1022 * t1114 + (-t1033 * t1114 + t1059 * t1119) * qJD(6), t1153 * t1114 + t1130 * t1119, -t1007 * t1119 - t1030 * t1190, t1005; t1006, 0, t1013 * t1119 + t1024 * t1114 + (-t1032 * t1114 + t1057 * t1119) * qJD(6), t1152 * t1114 + t1129 * t1119, t1009 * t1119 - t1027 * t1190, -t1010 * t1114 + t1016 * t1119 + (-t1043 * t1114 + t1211) * qJD(6); 0, 0, t1026 * t1119 + t1041 * t1114 + (-t1061 * t1114 + t1067 * t1119) * qJD(6), t1151 * t1114 + t1128 * t1119, t1020 * t1119 - t1039 * t1190, -t1021 * t1114 + t1034 * t1119 + (-t1040 * t1119 - t1062 * t1114) * qJD(6); -t1011 * t1114 + t1018 * t1119 + (t1114 * t1195 - t1211) * qJD(6), 0, -t1012 * t1114 + t1022 * t1119 + (-t1033 * t1119 - t1059 * t1114) * qJD(6), -t1130 * t1114 + t1153 * t1119, t1007 * t1114 - t1030 * t1189, -t1006; t1005, 0, -t1013 * t1114 + t1024 * t1119 + (-t1032 * t1119 - t1057 * t1114) * qJD(6), -t1129 * t1114 + t1152 * t1119, -t1009 * t1114 - t1027 * t1189, -t1010 * t1119 - t1016 * t1114 + (-t1043 * t1119 - t1212) * qJD(6); 0, 0, -t1026 * t1114 + t1041 * t1119 + (-t1061 * t1119 - t1067 * t1114) * qJD(6), -t1128 * t1114 + t1151 * t1119, -t1020 * t1114 - t1039 * t1189, -t1021 * t1119 - t1034 * t1114 + (t1040 * t1114 - t1062 * t1119) * qJD(6); t1019 * t1115 + t1210, 0, t1033 * qJD(5) + t1023 * t1115 - t1052 * t1173, -qJD(5) * t1187 - t1014 * t1115, t1008, 0; t1007, 0, t1032 * qJD(5) + t1025 * t1115 - t1054 * t1173, -qJD(5) * t1188 - t1016 * t1115, t1010, 0; 0, 0, t1061 * qJD(5) + t1042 * t1115 - t1084 * t1173, -qJD(5) * t1186 - t1034 * t1115, t1021, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end