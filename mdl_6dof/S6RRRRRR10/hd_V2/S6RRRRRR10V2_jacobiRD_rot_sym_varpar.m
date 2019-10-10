% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:38
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRRRR10V2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10V2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobiRD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
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
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t32 = sin(qJ(1));
	t39 = qJD(1) * t32;
	t34 = cos(qJ(1));
	t38 = qJD(1) * t34;
	t31 = sin(qJ(2));
	t37 = qJD(2) * t31;
	t33 = cos(qJ(2));
	t36 = qJD(2) * t33;
	t35 = qJD(2) * t34;
	t30 = t32 * t37 - t33 * t38;
	t29 = t31 * t38 + t32 * t36;
	t28 = t31 * t35 + t33 * t39;
	t27 = t31 * t39 - t33 * t35;
	t1 = [t30, t27, 0, 0, 0, 0; -t28, -t29, 0, 0, 0, 0; 0, -t37, 0, 0, 0, 0; t29, t28, 0, 0, 0, 0; t27, t30, 0, 0, 0, 0; 0, -t36, 0, 0, 0, 0; -t39, 0, 0, 0, 0, 0; t38, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (61->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t67 = qJ(2) + qJ(3);
	t64 = sin(t67);
	t66 = qJD(2) + qJD(3);
	t75 = t66 * t64;
	t65 = cos(t67);
	t74 = t66 * t65;
	t68 = sin(qJ(1));
	t73 = t66 * t68;
	t69 = cos(qJ(1));
	t72 = t66 * t69;
	t71 = qJD(1) * t68;
	t70 = qJD(1) * t69;
	t63 = t64 * t73 - t65 * t70;
	t62 = t64 * t70 + t65 * t73;
	t61 = t64 * t72 + t65 * t71;
	t60 = t64 * t71 - t65 * t72;
	t1 = [t63, t60, t60, 0, 0, 0; -t61, -t62, -t62, 0, 0, 0; 0, -t75, -t75, 0, 0, 0; t62, t61, t61, 0, 0, 0; t60, t63, t63, 0, 0, 0; 0, -t74, -t74, 0, 0, 0; -t71, 0, 0, 0, 0, 0; t70, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:03
	% EndTime: 2019-10-10 13:38:04
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (164->29), mult. (226->57), div. (0->0), fcn. (226->6), ass. (0->41)
	t327 = qJD(2) + qJD(3);
	t329 = sin(qJ(4));
	t351 = t327 * t329;
	t330 = sin(qJ(1));
	t350 = t327 * t330;
	t331 = cos(qJ(4));
	t349 = t327 * t331;
	t332 = cos(qJ(1));
	t348 = t327 * t332;
	t347 = t331 * t332;
	t346 = qJD(1) * t330;
	t345 = qJD(1) * t332;
	t344 = qJD(4) * t329;
	t343 = qJD(4) * t331;
	t342 = qJD(4) * t332;
	t328 = qJ(2) + qJ(3);
	t325 = sin(t328);
	t341 = t325 * t349;
	t340 = t325 * t350;
	t326 = cos(t328);
	t339 = t326 * t350;
	t338 = t325 * t348;
	t337 = t326 * t348;
	t336 = qJD(4) * t326 - qJD(1);
	t335 = qJD(1) * t326 - qJD(4);
	t334 = t336 * t329;
	t333 = t335 * t330 + t338;
	t324 = t327 * t326;
	t323 = t326 * t345 - t340;
	t322 = -t326 * t346 - t338;
	t321 = -t326 * t344 - t341;
	t320 = t325 * t351 - t326 * t343;
	t319 = -t331 * t339 + (t330 * t344 - t331 * t345) * t325;
	t318 = t329 * t339 + (t329 * t345 + t330 * t343) * t325;
	t317 = -t331 * t337 + (t329 * t342 + t331 * t346) * t325;
	t316 = t329 * t337 + (-t329 * t346 + t331 * t342) * t325;
	t315 = -t335 * t347 + (t334 + t341) * t330;
	t314 = t336 * t331 * t330 + (t335 * t332 - t340) * t329;
	t313 = t333 * t331 + t332 * t334;
	t312 = t333 * t329 - t336 * t347;
	t1 = [t315, t317, t317, t312, 0, 0; -t313, t319, t319, -t314, 0, 0; 0, t321, t321, -t325 * t343 - t326 * t351, 0, 0; t314, t316, t316, t313, 0, 0; t312, t318, t318, t315, 0, 0; 0, t320, t320, t325 * t344 - t326 * t349, 0, 0; -t325 * t345 - t339, t322, t322, 0, 0, 0; -t325 * t346 + t337, t323, t323, 0, 0, 0; 0, t324, t324, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:05
	% EndTime: 2019-10-10 13:38:06
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (454->58), mult. (696->109), div. (0->0), fcn. (708->8), ass. (0->72)
	t465 = qJ(2) + qJ(3);
	t462 = sin(t465);
	t466 = sin(qJ(5));
	t469 = cos(qJ(5));
	t464 = qJD(2) + qJD(3);
	t470 = cos(qJ(4));
	t487 = qJD(5) * t470 - t464;
	t467 = sin(qJ(4));
	t501 = qJD(4) * t467;
	t515 = -t466 * t501 + t487 * t469;
	t518 = t515 * t462;
	t517 = t487 * t466 + t469 * t501;
	t468 = sin(qJ(1));
	t463 = cos(t465);
	t485 = qJD(1) * t463 - qJD(4);
	t471 = cos(qJ(1));
	t508 = t464 * t471;
	t516 = t462 * t508 + t485 * t468;
	t502 = qJD(1) * t471;
	t505 = t471 * t467;
	t506 = t468 * t470;
	t510 = t464 * t468;
	t476 = -qJD(5) * (t463 * t506 - t505) + t462 * t502 + t463 * t510;
	t504 = t471 * t470;
	t507 = t468 * t467;
	t457 = t463 * t504 + t507;
	t500 = qJD(4) * t470;
	t489 = t471 * t500;
	t493 = t463 * t507;
	t496 = t462 * t510;
	t450 = t457 * qJD(1) - qJD(4) * t493 - t470 * t496 - t489;
	t499 = qJD(5) * t462;
	t483 = t468 * t499 + t450;
	t514 = t483 * t466 - t476 * t469;
	t513 = t462 * t470;
	t512 = t463 * t466;
	t511 = t464 * t467;
	t509 = t464 * t470;
	t503 = qJD(1) * t468;
	t498 = qJD(5) * t466;
	t497 = qJD(5) * t469;
	t494 = t463 * t511;
	t492 = t463 * t505;
	t488 = -qJD(5) + t509;
	t486 = -qJD(4) * t463 + qJD(1);
	t480 = t486 * t471;
	t448 = t467 * t480 - t516 * t470;
	t484 = t471 * t499 + t448;
	t482 = t488 * t469;
	t481 = t488 * t466;
	t479 = t469 * t513 - t512;
	t478 = t463 * t469 + t466 * t513;
	t475 = -qJD(5) * t457 - t462 * t503 + t463 * t508;
	t474 = t517 * t462 - t463 * t482;
	t473 = t463 * t481 + t518;
	t472 = -t476 * t466 - t483 * t469;
	t456 = -t492 + t506;
	t454 = -t493 - t504;
	t453 = -t462 * t511 + t463 * t500;
	t452 = -t464 * t493 + (-t467 * t502 - t468 * t500) * t462;
	t451 = -t464 * t492 + (t467 * t503 - t489) * t462;
	t449 = t486 * t506 + (-t485 * t471 + t496) * t467;
	t447 = t516 * t467 + t470 * t480;
	t446 = -t462 * t482 - t517 * t463;
	t445 = t462 * t481 - t463 * t515;
	t444 = t474 * t468 - t479 * t502;
	t443 = t473 * t468 + t478 * t502;
	t442 = t474 * t471 + t479 * t503;
	t441 = t473 * t471 - t478 * t503;
	t440 = t475 * t466 + t484 * t469;
	t439 = -t484 * t466 + t475 * t469;
	t1 = [t472, t442, t442, t447 * t469 - t456 * t498, t439, 0; t440, t444, t444, t449 * t469 - t454 * t498, -t514, 0; 0, t446, t446, -t469 * t494 + (t467 * t498 - t469 * t500) * t462, -t488 * t512 - t518, 0; t514, t441, t441, -t447 * t466 - t456 * t497, -t440, 0; t439, t443, t443, -t449 * t466 - t454 * t497, t472, 0; 0, t445, t445, t466 * t494 + (t466 * t500 + t467 * t497) * t462, t474, 0; t449, t451, t451, t448, 0, 0; -t447, t452, t452, t450, 0, 0; 0, t453, t453, -t462 * t501 + t463 * t509, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:08
	% EndTime: 2019-10-10 13:38:09
	% DurationCPUTime: 1.19s
	% Computational Cost: add. (1085->118), mult. (1726->213), div. (0->0), fcn. (1804->10), ass. (0->110)
	t690 = qJ(2) + qJ(3);
	t687 = sin(t690);
	t692 = sin(qJ(5));
	t693 = sin(qJ(4));
	t689 = qJD(2) + qJD(3);
	t697 = cos(qJ(4));
	t721 = qJD(5) * t697 - t689;
	t696 = cos(qJ(5));
	t745 = qJD(4) * t696;
	t769 = t721 * t692 + t693 * t745;
	t776 = t769 * t687;
	t688 = cos(t690);
	t698 = cos(qJ(1));
	t749 = t698 * t697;
	t694 = sin(qJ(1));
	t754 = t694 * t693;
	t673 = t688 * t749 + t754;
	t743 = qJD(4) * t698;
	t723 = t697 * t743;
	t746 = qJD(4) * t693;
	t724 = t694 * t746;
	t757 = t689 * t694;
	t735 = t687 * t757;
	t659 = t673 * qJD(1) - t688 * t724 - t697 * t735 - t723;
	t750 = t698 * t693;
	t753 = t694 * t697;
	t671 = t688 * t753 - t750;
	t747 = qJD(1) * t698;
	t728 = t687 * t747;
	t742 = qJD(5) * t687;
	t646 = (-qJD(5) * t671 + t688 * t757 + t728) * t692 + (t694 * t742 + t659) * t696;
	t744 = qJD(4) * t697;
	t708 = t693 * t747 + t694 * t744;
	t748 = qJD(1) * t694;
	t658 = t708 * t688 - t697 * t748 + (-t735 - t743) * t693;
	t762 = t687 * t694;
	t661 = t671 * t696 + t692 * t762;
	t731 = t688 * t754;
	t670 = t731 + t749;
	t691 = sin(qJ(6));
	t695 = cos(qJ(6));
	t775 = -t646 * t695 - t658 * t691 + (t661 * t691 - t670 * t695) * qJD(6);
	t774 = (t661 * t695 + t670 * t691) * qJD(6) + t646 * t691 - t658 * t695;
	t722 = t689 * t697 - qJD(5);
	t712 = t722 * t696;
	t738 = qJD(6) * t693;
	t725 = t687 * t738;
	t773 = t688 * t712 + t725 - t776;
	t751 = t696 * t697;
	t668 = t687 * t751 - t688 * t692;
	t770 = qJD(6) * t668;
	t768 = -t692 * t746 + t721 * t696;
	t763 = t687 * t692;
	t761 = t687 * t698;
	t760 = t688 * t689;
	t759 = t688 * t696;
	t758 = t689 * t693;
	t756 = t689 * t698;
	t755 = t693 * t696;
	t752 = t695 * t697;
	t741 = qJD(5) * t692;
	t740 = qJD(5) * t696;
	t739 = qJD(6) * t691;
	t737 = qJD(6) * t695;
	t736 = qJD(6) * t696;
	t734 = t687 * t756;
	t733 = t688 * t758;
	t732 = t689 * t759;
	t730 = t688 * t750;
	t729 = t687 * t748;
	t720 = -qJD(4) + t736;
	t719 = -qJD(6) + t745;
	t717 = -t687 * t712 + (-t769 + t738) * t688;
	t655 = t722 * t759 - t776;
	t716 = -t655 - t725;
	t657 = (-qJD(4) * t688 + qJD(1)) * t750 + (-t734 + (-qJD(1) * t688 + qJD(4)) * t694) * t697;
	t672 = t730 - t753;
	t715 = t672 * t736 + t657;
	t714 = t670 * t736 + t659;
	t713 = t722 * t692;
	t711 = -t668 * t748 + t773 * t698;
	t666 = t668 * t698;
	t710 = qJD(1) * t666 + t773 * t694;
	t664 = t673 * t696 + t692 * t761;
	t709 = t697 * t763 + t759;
	t656 = t670 * qJD(1) - t688 * t723 + t693 * t734 - t724;
	t707 = qJD(6) * t673 + t656 * t696 + t672 * t741;
	t706 = qJD(6) * t671 - t658 * t696 + t670 * t741;
	t703 = t687 * t744 + t733 - t770;
	t702 = -qJD(6) * (t688 * t751 + t763) - t687 * t758 + t688 * t744;
	t645 = -t661 * qJD(5) - t659 * t692 + t694 * t732 + t696 * t728;
	t701 = -t708 * t687 - t689 * t731 + t694 * t770;
	t700 = -t689 * t730 + qJD(6) * t666 + (t693 * t748 - t723) * t687;
	t654 = -t768 * t687 - t688 * t713;
	t663 = -t673 * t692 + t696 * t761;
	t660 = -t671 * t692 + t696 * t762;
	t652 = -t687 * t713 + t768 * t688;
	t650 = t654 * t694 - t709 * t747;
	t648 = t654 * t698 + t709 * t748;
	t644 = (t698 * t742 + t657) * t696 + (-qJD(5) * t673 + t688 * t756 - t729) * t692;
	t643 = t664 * qJD(5) + t657 * t692 + t696 * t729 - t698 * t732;
	t642 = t702 * t691 + t717 * t695;
	t641 = -t717 * t691 + t702 * t695;
	t640 = t701 * t691 - t710 * t695;
	t639 = t710 * t691 + t701 * t695;
	t638 = t700 * t691 - t711 * t695;
	t637 = t711 * t691 + t700 * t695;
	t636 = t644 * t695 - t656 * t691 + (-t664 * t691 + t672 * t695) * qJD(6);
	t635 = -t644 * t691 - t656 * t695 + (-t664 * t695 - t672 * t691) * qJD(6);
	t1 = [t775, t638, t638, t715 * t691 + t707 * t695, -t643 * t695 - t663 * t739, t635; t636, t640, t640, t714 * t691 + t706 * t695, t645 * t695 - t660 * t739, -t774; 0, t642, t642, (t691 * t697 - t695 * t755) * t760 + (-t719 * t752 + (t720 * t691 + t695 * t741) * t693) * t687, t654 * t695 + t709 * t739, t716 * t691 + t703 * t695; t774, t637, t637, -t707 * t691 + t715 * t695, t643 * t691 - t663 * t737, -t636; t635, t639, t639, -t706 * t691 + t714 * t695, -t645 * t691 - t660 * t737, t775; 0, t641, t641, (t691 * t755 + t752) * t760 + (t720 * t695 * t693 + (-t693 * t741 + t719 * t697) * t691) * t687, -t654 * t691 + t709 * t737, -t703 * t691 + t716 * t695; t645, t648, t648, t656 * t692 - t672 * t740, t644, 0; t643, t650, t650, -t658 * t692 - t670 * t740, t646, 0; 0, t652, t652, -t692 * t733 + (-t692 * t744 - t693 * t740) * t687, t655, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end