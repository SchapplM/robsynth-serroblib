% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRRP12_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP12_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
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
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->5), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t128 = cos(pkin(6));
	t129 = sin(qJ(1));
	t133 = t128 * t129;
	t130 = cos(qJ(1));
	t132 = t128 * t130;
	t131 = qJD(1) * sin(pkin(6));
	t127 = cos(pkin(12));
	t125 = sin(pkin(12));
	t1 = [(t125 * t133 - t127 * t130) * qJD(1), 0, 0, 0, 0, 0; (-t125 * t132 - t127 * t129) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; (t125 * t130 + t127 * t133) * qJD(1), 0, 0, 0, 0, 0; (t125 * t129 - t127 * t132) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t129 * t131, 0, 0, 0, 0, 0; t130 * t131, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:25
	% EndTime: 2019-10-10 08:58:25
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (80->29), mult. (294->61), div. (0->0), fcn. (310->10), ass. (0->38)
	t298 = cos(pkin(6));
	t296 = cos(pkin(12));
	t302 = cos(qJ(1));
	t314 = t302 * t296;
	t293 = sin(pkin(12));
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
	% StartTime: 2019-10-10 08:58:27
	% EndTime: 2019-10-10 08:58:28
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (278->59), mult. (936->116), div. (0->0), fcn. (1042->12), ass. (0->63)
	t506 = cos(pkin(6));
	t501 = sin(pkin(12));
	t512 = cos(qJ(1));
	t529 = t512 * t501;
	t504 = cos(pkin(12));
	t509 = sin(qJ(1));
	t530 = t509 * t504;
	t514 = t506 * t530 + t529;
	t489 = t514 * qJD(1);
	t531 = t509 * t501;
	t522 = t506 * t531;
	t527 = qJD(1) * t512;
	t490 = -qJD(1) * t522 + t504 * t527;
	t505 = cos(pkin(7));
	t508 = sin(qJ(3));
	t511 = cos(qJ(3));
	t502 = sin(pkin(7));
	t503 = sin(pkin(6));
	t535 = t503 * t509;
	t521 = qJD(1) * t535;
	t519 = t502 * t521;
	t493 = t506 * t529 + t530;
	t528 = t512 * t504;
	t492 = -t506 * t528 + t531;
	t534 = t503 * t512;
	t524 = t502 * t534;
	t517 = t492 * t505 + t524;
	t540 = t493 * t508 + t517 * t511;
	t470 = t540 * qJD(3) - (-t489 * t505 + t519) * t508 - t490 * t511;
	t537 = t493 * t511;
	t473 = t517 * t508 - t537;
	t480 = t489 * t502 + t505 * t521;
	t483 = -t492 * t502 + t505 * t534;
	t507 = sin(qJ(4));
	t510 = cos(qJ(4));
	t550 = t470 * t507 + (t473 * t510 + t483 * t507) * qJD(4) + t480 * t510;
	t549 = t470 * t510 - t480 * t507 + (-t473 * t507 + t483 * t510) * qJD(4);
	t533 = t505 * t508;
	t536 = t502 * t506;
	t482 = (t501 * t511 + t504 * t533) * t503 + t508 * t536;
	t532 = t505 * t511;
	t526 = qJD(4) * t507;
	t525 = qJD(4) * t510;
	t520 = t503 * t527;
	t518 = t502 * t520;
	t516 = t502 * t535 - t505 * t514;
	t495 = -t522 + t528;
	t474 = -t495 * t508 + t516 * t511;
	t475 = t495 * t511 + t516 * t508;
	t481 = t511 * t536 + (-t501 * t508 + t504 * t532) * t503;
	t468 = -t489 * t532 - t490 * t508 + t511 * t519 + (t492 * t533 + t508 * t524 - t537) * qJD(3);
	t491 = -t503 * t504 * t502 + t506 * t505;
	t488 = t493 * qJD(1);
	t487 = t492 * qJD(1);
	t485 = t502 * t514 + t505 * t535;
	t478 = -t487 * t502 + t505 * t520;
	t477 = t482 * qJD(3);
	t476 = t481 * qJD(3);
	t467 = -t488 * t511 + (t487 * t505 + t518) * t508 + t474 * qJD(3);
	t466 = t475 * qJD(3) - t487 * t532 - t488 * t508 - t511 * t518;
	t465 = t467 * t510 + t478 * t507 + (-t475 * t507 + t485 * t510) * qJD(4);
	t464 = -t467 * t507 + t478 * t510 + (-t475 * t510 - t485 * t507) * qJD(4);
	t1 = [t549, 0, -t466 * t510 - t474 * t526, t464, 0, 0; t465, 0, t468 * t510 + t526 * t540, t550, 0, 0; 0, 0, -t477 * t510 - t481 * t526, -t476 * t507 + (-t482 * t510 - t491 * t507) * qJD(4), 0, 0; -t550, 0, t466 * t507 - t474 * t525, -t465, 0, 0; t464, 0, -t468 * t507 + t525 * t540, t549, 0, 0; 0, 0, t477 * t507 - t481 * t525, -t476 * t510 + (t482 * t507 - t491 * t510) * qJD(4), 0, 0; t468, 0, t467, 0, 0, 0; t466, 0, -t470, 0, 0, 0; 0, 0, t476, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:32
	% EndTime: 2019-10-10 08:58:33
	% DurationCPUTime: 1.35s
	% Computational Cost: add. (800->98), mult. (2579->181), div. (0->0), fcn. (2964->14), ass. (0->89)
	t729 = cos(pkin(6));
	t724 = sin(pkin(12));
	t737 = cos(qJ(1));
	t769 = t737 * t724;
	t727 = cos(pkin(12));
	t733 = sin(qJ(1));
	t770 = t733 * t727;
	t747 = t729 * t770 + t769;
	t711 = t747 * qJD(1);
	t768 = t737 * t727;
	t771 = t733 * t724;
	t716 = -t729 * t771 + t768;
	t712 = t716 * qJD(1);
	t725 = sin(pkin(7));
	t728 = cos(pkin(7));
	t732 = sin(qJ(3));
	t736 = cos(qJ(3));
	t726 = sin(pkin(6));
	t775 = t726 * t733;
	t755 = qJD(1) * t775;
	t748 = t729 * t768 - t771;
	t745 = t748 * t728;
	t774 = t726 * t737;
	t757 = t725 * t774;
	t749 = -t745 + t757;
	t715 = t729 * t769 + t770;
	t781 = t715 * t732;
	t783 = t749 * t736 + t781;
	t680 = qJD(3) * t783 - (-t711 * t728 + t725 * t755) * t732 - t712 * t736;
	t780 = t715 * t736;
	t793 = t749 * t732;
	t693 = -t780 + t793;
	t703 = t748 * t725 + t728 * t774;
	t731 = sin(qJ(4));
	t735 = cos(qJ(4));
	t683 = t693 * t731 - t703 * t735;
	t700 = t711 * t725 + t728 * t755;
	t674 = -t683 * qJD(4) + t680 * t735 - t700 * t731;
	t767 = qJD(1) * t736;
	t776 = t726 * t725;
	t754 = t767 * t776;
	t772 = t728 * t736;
	t773 = t728 * t732;
	t679 = (-t748 * t773 - t780) * qJD(3) - t711 * t772 + t733 * t754 - (-qJD(3) * t757 + t712) * t732;
	t730 = sin(qJ(5));
	t734 = cos(qJ(5));
	t800 = -t674 * t730 + t679 * t734;
	t799 = t674 * t734 + t679 * t730;
	t685 = t693 * t735 + t703 * t731;
	t798 = t685 * t730;
	t797 = t685 * t734;
	t672 = t685 * qJD(4) + t680 * t731 + t700 * t735;
	t786 = -t725 * t775 + t747 * t728;
	t694 = t716 * t732 + t736 * t786;
	t777 = t725 * t729;
	t787 = t726 * (-t724 * t732 + t727 * t772) + t736 * t777;
	t766 = qJD(4) * t731;
	t765 = qJD(4) * t735;
	t764 = qJD(5) * t730;
	t763 = qJD(5) * t734;
	t762 = qJD(5) * t735;
	t710 = t715 * qJD(1);
	t676 = qJD(1) * t793 - t694 * qJD(3) - t710 * t736;
	t753 = t694 * t762 + t676;
	t690 = t736 * t757 - t748 * t772 + t781;
	t752 = t690 * t762 - t680;
	t697 = t787 * qJD(3);
	t751 = -t762 * t787 + t697;
	t695 = t716 * t736 - t732 * t786;
	t705 = t747 * t725 + t728 * t775;
	t687 = t695 * t735 + t705 * t731;
	t702 = t732 * t777 + (t724 * t736 + t727 * t773) * t726;
	t713 = -t727 * t776 + t729 * t728;
	t689 = t702 * t735 + t713 * t731;
	t688 = -t702 * t731 + t713 * t735;
	t675 = t695 * qJD(3) - t710 * t732 - t737 * t754 + t745 * t767;
	t742 = qJD(5) * t695 - t675 * t735 + t694 * t766;
	t741 = -qJD(5) * t693 + t679 * t735 + t690 * t766;
	t698 = t702 * qJD(3);
	t740 = qJD(5) * t702 - t698 * t735 - t766 * t787;
	t738 = t703 * qJD(1);
	t686 = -t695 * t731 + t705 * t735;
	t682 = t688 * qJD(4) + t697 * t735;
	t681 = -t689 * qJD(4) - t697 * t731;
	t671 = t676 * t735 - t695 * t766 + t705 * t765 + t731 * t738;
	t670 = t687 * qJD(4) + t676 * t731 - t735 * t738;
	t669 = t671 * t734 + t675 * t730 + (-t687 * t730 + t694 * t734) * qJD(5);
	t668 = -t671 * t730 + t675 * t734 + (-t687 * t734 - t694 * t730) * qJD(5);
	t1 = [(-t734 * t783 - t798) * qJD(5) + t799, 0, t753 * t730 + t742 * t734, -t670 * t734 - t686 * t764, t668, 0; t669, 0, t752 * t730 + t741 * t734, t672 * t734 - t683 * t764, (-t690 * t730 + t797) * qJD(5) - t800, 0; 0, 0, t751 * t730 + t740 * t734, t681 * t734 - t688 * t764, -t682 * t730 + t698 * t734 + (-t689 * t734 + t730 * t787) * qJD(5), 0; (t730 * t783 - t797) * qJD(5) + t800, 0, -t742 * t730 + t753 * t734, t670 * t730 - t686 * t763, -t669, 0; t668, 0, -t741 * t730 + t752 * t734, -t672 * t730 - t683 * t763, (-t690 * t734 - t798) * qJD(5) + t799, 0; 0, 0, -t740 * t730 + t751 * t734, -t681 * t730 - t688 * t763, -t682 * t734 - t698 * t730 + (t689 * t730 + t734 * t787) * qJD(5), 0; t672, 0, -t675 * t731 - t694 * t765, t671, 0, 0; t670, 0, t679 * t731 - t690 * t765, -t674, 0, 0; 0, 0, -t698 * t731 + t765 * t787, t682, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:36
	% EndTime: 2019-10-10 08:58:38
	% DurationCPUTime: 1.28s
	% Computational Cost: add. (800->98), mult. (2579->181), div. (0->0), fcn. (2964->14), ass. (0->89)
	t810 = cos(pkin(6));
	t805 = sin(pkin(12));
	t818 = cos(qJ(1));
	t850 = t818 * t805;
	t808 = cos(pkin(12));
	t814 = sin(qJ(1));
	t851 = t814 * t808;
	t828 = t810 * t851 + t850;
	t792 = t828 * qJD(1);
	t849 = t818 * t808;
	t852 = t814 * t805;
	t797 = -t810 * t852 + t849;
	t793 = t797 * qJD(1);
	t806 = sin(pkin(7));
	t809 = cos(pkin(7));
	t813 = sin(qJ(3));
	t817 = cos(qJ(3));
	t807 = sin(pkin(6));
	t856 = t807 * t814;
	t836 = qJD(1) * t856;
	t829 = t810 * t849 - t852;
	t826 = t829 * t809;
	t855 = t807 * t818;
	t838 = t806 * t855;
	t830 = -t826 + t838;
	t796 = t810 * t850 + t851;
	t862 = t796 * t813;
	t864 = t830 * t817 + t862;
	t761 = t864 * qJD(3) - (-t792 * t809 + t806 * t836) * t813 - t793 * t817;
	t861 = t796 * t817;
	t874 = t830 * t813;
	t774 = -t861 + t874;
	t784 = t829 * t806 + t809 * t855;
	t812 = sin(qJ(4));
	t816 = cos(qJ(4));
	t764 = t774 * t812 - t784 * t816;
	t781 = t792 * t806 + t809 * t836;
	t755 = -t764 * qJD(4) + t761 * t816 - t781 * t812;
	t848 = qJD(1) * t817;
	t857 = t807 * t806;
	t835 = t848 * t857;
	t853 = t809 * t817;
	t854 = t809 * t813;
	t760 = (-t829 * t854 - t861) * qJD(3) - t792 * t853 + t814 * t835 - (-qJD(3) * t838 + t793) * t813;
	t811 = sin(qJ(5));
	t815 = cos(qJ(5));
	t881 = t755 * t811 - t760 * t815;
	t880 = -t755 * t815 - t760 * t811;
	t766 = t774 * t816 + t784 * t812;
	t879 = t766 * t811;
	t878 = t766 * t815;
	t753 = t766 * qJD(4) + t761 * t812 + t781 * t816;
	t867 = -t806 * t856 + t828 * t809;
	t775 = t797 * t813 + t867 * t817;
	t858 = t806 * t810;
	t868 = (-t805 * t813 + t808 * t853) * t807 + t817 * t858;
	t847 = qJD(4) * t812;
	t846 = qJD(4) * t816;
	t845 = qJD(5) * t811;
	t844 = qJD(5) * t815;
	t843 = qJD(5) * t816;
	t791 = t796 * qJD(1);
	t757 = qJD(1) * t874 - t775 * qJD(3) - t791 * t817;
	t834 = t775 * t843 + t757;
	t771 = t817 * t838 - t829 * t853 + t862;
	t833 = t771 * t843 - t761;
	t778 = t868 * qJD(3);
	t832 = -t843 * t868 + t778;
	t776 = t797 * t817 - t813 * t867;
	t786 = t828 * t806 + t809 * t856;
	t768 = t776 * t816 + t786 * t812;
	t783 = t813 * t858 + (t805 * t817 + t808 * t854) * t807;
	t794 = -t808 * t857 + t810 * t809;
	t770 = t783 * t816 + t794 * t812;
	t769 = -t783 * t812 + t794 * t816;
	t756 = t776 * qJD(3) - t791 * t813 - t818 * t835 + t826 * t848;
	t823 = qJD(5) * t776 - t756 * t816 + t775 * t847;
	t822 = -qJD(5) * t774 + t760 * t816 + t771 * t847;
	t779 = t783 * qJD(3);
	t821 = qJD(5) * t783 - t779 * t816 - t847 * t868;
	t819 = t784 * qJD(1);
	t767 = -t776 * t812 + t786 * t816;
	t763 = t769 * qJD(4) + t778 * t816;
	t762 = -t770 * qJD(4) - t778 * t812;
	t752 = t757 * t816 - t776 * t847 + t786 * t846 + t812 * t819;
	t751 = t768 * qJD(4) + t757 * t812 - t816 * t819;
	t750 = t752 * t815 + t756 * t811 + (-t768 * t811 + t775 * t815) * qJD(5);
	t749 = t752 * t811 - t756 * t815 + (t768 * t815 + t775 * t811) * qJD(5);
	t1 = [(-t815 * t864 - t879) * qJD(5) - t880, 0, t834 * t811 + t823 * t815, -t751 * t815 - t767 * t845, -t749, 0; t750, 0, t833 * t811 + t822 * t815, t753 * t815 - t764 * t845, (-t771 * t811 + t878) * qJD(5) + t881, 0; 0, 0, t832 * t811 + t821 * t815, t762 * t815 - t769 * t845, -t763 * t811 + t779 * t815 + (-t770 * t815 + t811 * t868) * qJD(5), 0; t753, 0, -t756 * t812 - t775 * t846, t752, 0, 0; t751, 0, t760 * t812 - t771 * t846, -t755, 0, 0; 0, 0, -t779 * t812 + t846 * t868, t763, 0, 0; (-t811 * t864 + t878) * qJD(5) + t881, 0, t823 * t811 - t834 * t815, -t751 * t811 + t767 * t844, t750, 0; t749, 0, t822 * t811 - t833 * t815, t753 * t811 + t764 * t844, (t771 * t815 + t879) * qJD(5) + t880, 0; 0, 0, t821 * t811 - t832 * t815, t762 * t811 + t769 * t844, t763 * t815 + t779 * t811 + (-t770 * t811 - t815 * t868) * qJD(5), 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end