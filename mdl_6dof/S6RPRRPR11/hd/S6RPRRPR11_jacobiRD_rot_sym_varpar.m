% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRPR11_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR11_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
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
	% StartTime: 2019-10-10 01:41:21
	% EndTime: 2019-10-10 01:41:21
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
	% StartTime: 2019-10-10 01:41:21
	% EndTime: 2019-10-10 01:41:21
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
	% StartTime: 2019-10-10 01:41:23
	% EndTime: 2019-10-10 01:41:24
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
	% StartTime: 2019-10-10 01:41:27
	% EndTime: 2019-10-10 01:41:28
	% DurationCPUTime: 0.64s
	% Computational Cost: add. (449->67), mult. (1505->136), div. (0->0), fcn. (1668->14), ass. (0->68)
	t634 = cos(pkin(6));
	t628 = sin(pkin(12));
	t640 = cos(qJ(1));
	t663 = t640 * t628;
	t632 = cos(pkin(12));
	t637 = sin(qJ(1));
	t664 = t637 * t632;
	t648 = t634 * t664 + t663;
	t615 = t648 * qJD(1);
	t662 = t640 * t632;
	t665 = t637 * t628;
	t621 = -t634 * t665 + t662;
	t616 = t621 * qJD(1);
	t629 = sin(pkin(7));
	t633 = cos(pkin(7));
	t636 = sin(qJ(3));
	t639 = cos(qJ(3));
	t630 = sin(pkin(6));
	t669 = t630 * t637;
	t654 = qJD(1) * t669;
	t619 = t634 * t663 + t664;
	t649 = t634 * t662 - t665;
	t643 = t649 * t633;
	t668 = t630 * t640;
	t656 = t629 * t668;
	t652 = -t643 + t656;
	t675 = t619 * t636 + t652 * t639;
	t596 = t675 * qJD(3) - (-t615 * t633 + t629 * t654) * t636 - t616 * t639;
	t672 = t619 * t639;
	t683 = t652 * t636;
	t600 = -t672 + t683;
	t607 = t615 * t629 + t633 * t654;
	t610 = t649 * t629 + t633 * t668;
	t635 = sin(qJ(4));
	t638 = cos(qJ(4));
	t589 = t596 * t635 + t607 * t638 + (t600 * t638 + t610 * t635) * qJD(4);
	t590 = t596 * t638 - (t600 * t635 - t610 * t638) * qJD(4) - t607 * t635;
	t667 = t633 * t636;
	t671 = t629 * t634;
	t609 = (t628 * t639 + t632 * t667) * t630 + t636 * t671;
	t670 = t630 * t629;
	t666 = t633 * t639;
	t661 = qJD(1) * t639;
	t660 = qJD(4) * t635;
	t659 = qJD(4) * t638;
	t653 = t661 * t670;
	t651 = t629 * t669 - t633 * t648;
	t602 = t621 * t639 + t651 * t636;
	t614 = t619 * qJD(1);
	t591 = t602 * qJD(3) - t614 * t636 - t640 * t653 + t643 * t661;
	t601 = -t621 * t636 + t651 * t639;
	t647 = t591 * t638 + t601 * t660;
	t641 = -t615 * t666 - t616 * t636 + t637 * t653 + (t636 * t656 - t649 * t667 - t672) * qJD(3);
	t646 = -t638 * t641 - t660 * t675;
	t605 = t609 * qJD(3);
	t608 = t639 * t671 + (-t628 * t636 + t632 * t666) * t630;
	t645 = t605 * t638 + t608 * t660;
	t642 = t610 * qJD(1);
	t631 = cos(pkin(13));
	t627 = sin(pkin(13));
	t617 = -t632 * t670 + t634 * t633;
	t612 = t629 * t648 + t633 * t669;
	t604 = t608 * qJD(3);
	t597 = -t604 * t635 + (-t609 * t638 - t617 * t635) * qJD(4);
	t592 = qJD(1) * t683 + t601 * qJD(3) - t614 * t639;
	t588 = t592 * t638 - t602 * t660 + t612 * t659 + t635 * t642;
	t587 = t592 * t635 - t638 * t642 + (t602 * t638 + t612 * t635) * qJD(4);
	t1 = [t590 * t631 + t627 * t641, 0, t592 * t627 - t647 * t631, -t587 * t631, 0, 0; t588 * t631 + t591 * t627, 0, -t596 * t627 - t646 * t631, t589 * t631, 0, 0; 0, 0, t604 * t627 - t645 * t631, t597 * t631, 0, 0; -t590 * t627 + t631 * t641, 0, t592 * t631 + t647 * t627, t587 * t627, 0, 0; -t588 * t627 + t591 * t631, 0, -t596 * t631 + t646 * t627, -t589 * t627, 0, 0; 0, 0, t604 * t631 + t645 * t627, -t597 * t627, 0, 0; t589, 0, -t591 * t635 + t601 * t659, t588, 0, 0; t587, 0, t635 * t641 - t659 * t675, -t590, 0, 0; 0, 0, -t605 * t635 + t608 * t659, t604 * t638 + (-t609 * t635 + t617 * t638) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:28
	% EndTime: 2019-10-10 01:41:30
	% DurationCPUTime: 1.30s
	% Computational Cost: add. (882->99), mult. (2579->181), div. (0->0), fcn. (2964->14), ass. (0->90)
	t748 = cos(pkin(6));
	t743 = sin(pkin(12));
	t754 = cos(qJ(1));
	t786 = t754 * t743;
	t746 = cos(pkin(12));
	t751 = sin(qJ(1));
	t787 = t751 * t746;
	t764 = t748 * t787 + t786;
	t727 = t764 * qJD(1);
	t785 = t754 * t746;
	t788 = t751 * t743;
	t732 = -t748 * t788 + t785;
	t728 = t732 * qJD(1);
	t744 = sin(pkin(7));
	t747 = cos(pkin(7));
	t750 = sin(qJ(3));
	t753 = cos(qJ(3));
	t745 = sin(pkin(6));
	t792 = t745 * t751;
	t772 = qJD(1) * t792;
	t765 = t748 * t785 - t788;
	t762 = t765 * t747;
	t791 = t745 * t754;
	t774 = t744 * t791;
	t766 = -t762 + t774;
	t731 = t748 * t786 + t787;
	t798 = t731 * t750;
	t800 = t766 * t753 + t798;
	t696 = t800 * qJD(3) - (-t727 * t747 + t744 * t772) * t750 - t728 * t753;
	t797 = t731 * t753;
	t810 = t766 * t750;
	t709 = -t797 + t810;
	t719 = t765 * t744 + t747 * t791;
	t749 = sin(qJ(4));
	t752 = cos(qJ(4));
	t699 = t709 * t749 - t719 * t752;
	t716 = t727 * t744 + t747 * t772;
	t690 = -t699 * qJD(4) + t696 * t752 - t716 * t749;
	t784 = qJD(1) * t753;
	t793 = t745 * t744;
	t771 = t784 * t793;
	t789 = t747 * t753;
	t790 = t747 * t750;
	t695 = (-t765 * t790 - t797) * qJD(3) - t727 * t789 + t751 * t771 - (-qJD(3) * t774 + t728) * t750;
	t742 = pkin(13) + qJ(6);
	t740 = sin(t742);
	t741 = cos(t742);
	t817 = -t690 * t740 + t695 * t741;
	t816 = t690 * t741 + t695 * t740;
	t701 = t709 * t752 + t719 * t749;
	t815 = t701 * t740;
	t814 = t701 * t741;
	t688 = t701 * qJD(4) + t696 * t749 + t716 * t752;
	t803 = -t744 * t792 + t764 * t747;
	t710 = t732 * t750 + t803 * t753;
	t794 = t744 * t748;
	t804 = (-t743 * t750 + t746 * t789) * t745 + t753 * t794;
	t783 = qJD(4) * t749;
	t782 = qJD(4) * t752;
	t781 = qJD(6) * t740;
	t780 = qJD(6) * t741;
	t779 = qJD(6) * t752;
	t726 = t731 * qJD(1);
	t692 = qJD(1) * t810 - t710 * qJD(3) - t726 * t753;
	t770 = t710 * t779 + t692;
	t706 = t753 * t774 - t765 * t789 + t798;
	t769 = t706 * t779 - t696;
	t713 = t804 * qJD(3);
	t768 = -t779 * t804 + t713;
	t711 = t732 * t753 - t750 * t803;
	t721 = t764 * t744 + t747 * t792;
	t703 = t711 * t752 + t721 * t749;
	t718 = t750 * t794 + (t743 * t753 + t746 * t790) * t745;
	t729 = -t746 * t793 + t748 * t747;
	t705 = t718 * t752 + t729 * t749;
	t704 = -t718 * t749 + t729 * t752;
	t691 = t711 * qJD(3) - t726 * t750 - t754 * t771 + t762 * t784;
	t759 = qJD(6) * t711 - t691 * t752 + t710 * t783;
	t758 = -qJD(6) * t709 + t695 * t752 + t706 * t783;
	t714 = t718 * qJD(3);
	t757 = qJD(6) * t718 - t714 * t752 - t783 * t804;
	t755 = t719 * qJD(1);
	t702 = -t711 * t749 + t721 * t752;
	t698 = t704 * qJD(4) + t713 * t752;
	t697 = -t705 * qJD(4) - t713 * t749;
	t687 = t692 * t752 - t711 * t783 + t721 * t782 + t749 * t755;
	t686 = t703 * qJD(4) + t692 * t749 - t752 * t755;
	t685 = t687 * t741 + t691 * t740 + (-t703 * t740 + t710 * t741) * qJD(6);
	t684 = -t687 * t740 + t691 * t741 + (-t703 * t741 - t710 * t740) * qJD(6);
	t1 = [(-t741 * t800 - t815) * qJD(6) + t816, 0, t770 * t740 + t759 * t741, -t686 * t741 - t702 * t781, 0, t684; t685, 0, t769 * t740 + t758 * t741, t688 * t741 - t699 * t781, 0, (-t706 * t740 + t814) * qJD(6) - t817; 0, 0, t768 * t740 + t757 * t741, t697 * t741 - t704 * t781, 0, -t698 * t740 + t714 * t741 + (-t705 * t741 + t740 * t804) * qJD(6); (t740 * t800 - t814) * qJD(6) + t817, 0, -t759 * t740 + t770 * t741, t686 * t740 - t702 * t780, 0, -t685; t684, 0, -t758 * t740 + t769 * t741, -t688 * t740 - t699 * t780, 0, (-t706 * t741 - t815) * qJD(6) + t816; 0, 0, -t757 * t740 + t768 * t741, -t697 * t740 - t704 * t780, 0, -t698 * t741 - t714 * t740 + (t705 * t740 + t741 * t804) * qJD(6); t688, 0, -t691 * t749 - t710 * t782, t687, 0, 0; t686, 0, t695 * t749 - t706 * t782, -t690, 0, 0; 0, 0, -t714 * t749 + t782 * t804, t698, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end