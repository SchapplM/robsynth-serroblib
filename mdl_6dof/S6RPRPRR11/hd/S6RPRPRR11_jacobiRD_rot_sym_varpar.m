% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPRR11_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR11_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
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
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.06s
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
	% StartTime: 2019-10-10 01:04:01
	% EndTime: 2019-10-10 01:04:01
	% DurationCPUTime: 0.17s
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
	% StartTime: 2019-10-10 01:04:03
	% EndTime: 2019-10-10 01:04:03
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (129->38), mult. (464->79), div. (0->0), fcn. (488->12), ass. (0->44)
	t432 = cos(pkin(6));
	t426 = sin(pkin(12));
	t436 = cos(qJ(1));
	t449 = t436 * t426;
	t430 = cos(pkin(12));
	t434 = sin(qJ(1));
	t450 = t434 * t430;
	t438 = t432 * t450 + t449;
	t414 = t438 * qJD(1);
	t451 = t434 * t426;
	t445 = t432 * t451;
	t447 = qJD(1) * t436;
	t415 = -qJD(1) * t445 + t430 * t447;
	t448 = t436 * t430;
	t416 = -t432 * t448 + t451;
	t417 = t432 * t449 + t450;
	t431 = cos(pkin(7));
	t433 = sin(qJ(3));
	t435 = cos(qJ(3));
	t427 = sin(pkin(7));
	t428 = sin(pkin(6));
	t454 = t428 * t434;
	t444 = qJD(1) * t454;
	t442 = t427 * t444;
	t446 = t427 * t428 * t436;
	t407 = ((t416 * t431 + t446) * t435 + t417 * t433) * qJD(3) - (-t414 * t431 + t442) * t433 - t415 * t435;
	t455 = t427 * t432;
	t453 = t431 * t433;
	t452 = t431 * t435;
	t443 = t428 * t447;
	t441 = t427 * t443;
	t439 = t427 * t454 - t431 * t438;
	t406 = -t414 * t452 - t415 * t433 + t435 * t442 + (t416 * t453 - t417 * t435 + t433 * t446) * qJD(3);
	t429 = cos(pkin(13));
	t425 = sin(pkin(13));
	t419 = -t445 + t448;
	t413 = t417 * qJD(1);
	t412 = t416 * qJD(1);
	t410 = -t414 * t427 - t431 * t444;
	t409 = -t412 * t427 + t431 * t443;
	t408 = (-t433 * t455 + (-t426 * t435 - t430 * t453) * t428) * qJD(3);
	t405 = -t413 * t435 + (t412 * t431 + t441) * t433 + (-t419 * t433 + t439 * t435) * qJD(3);
	t404 = -t413 * t433 - t412 * t452 - t435 * t441 + (t419 * t435 + t439 * t433) * qJD(3);
	t1 = [t407 * t429 + t410 * t425, 0, -t404 * t429, 0, 0, 0; t405 * t429 + t409 * t425, 0, t406 * t429, 0, 0, 0; 0, 0, t408 * t429, 0, 0, 0; -t407 * t425 + t410 * t429, 0, t404 * t425, 0, 0, 0; -t405 * t425 + t409 * t429, 0, -t406 * t425, 0, 0, 0; 0, 0, -t408 * t425, 0, 0, 0; t406, 0, t405, 0, 0, 0; t404, 0, -t407, 0, 0, 0; 0, 0, (t435 * t455 + (-t426 * t433 + t430 * t452) * t428) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:03
	% EndTime: 2019-10-10 01:04:04
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (328->60), mult. (936->116), div. (0->0), fcn. (1042->12), ass. (0->64)
	t530 = cos(pkin(6));
	t525 = sin(pkin(12));
	t534 = cos(qJ(1));
	t551 = t534 * t525;
	t528 = cos(pkin(12));
	t532 = sin(qJ(1));
	t552 = t532 * t528;
	t536 = t530 * t552 + t551;
	t510 = t536 * qJD(1);
	t553 = t532 * t525;
	t544 = t530 * t553;
	t549 = qJD(1) * t534;
	t511 = -qJD(1) * t544 + t528 * t549;
	t529 = cos(pkin(7));
	t531 = sin(qJ(3));
	t533 = cos(qJ(3));
	t526 = sin(pkin(7));
	t527 = sin(pkin(6));
	t557 = t527 * t532;
	t543 = qJD(1) * t557;
	t541 = t526 * t543;
	t514 = t530 * t551 + t552;
	t550 = t534 * t528;
	t513 = -t530 * t550 + t553;
	t556 = t527 * t534;
	t546 = t526 * t556;
	t539 = t513 * t529 + t546;
	t562 = t514 * t531 + t539 * t533;
	t491 = t562 * qJD(3) - (-t510 * t529 + t541) * t531 - t511 * t533;
	t559 = t514 * t533;
	t494 = t539 * t531 - t559;
	t501 = t510 * t526 + t529 * t543;
	t504 = -t513 * t526 + t529 * t556;
	t524 = pkin(13) + qJ(5);
	t522 = sin(t524);
	t523 = cos(t524);
	t572 = t491 * t522 + (t494 * t523 + t504 * t522) * qJD(5) + t501 * t523;
	t571 = t491 * t523 - t501 * t522 + (-t494 * t522 + t504 * t523) * qJD(5);
	t555 = t529 * t531;
	t558 = t526 * t530;
	t503 = (t525 * t533 + t528 * t555) * t527 + t531 * t558;
	t554 = t529 * t533;
	t548 = qJD(5) * t522;
	t547 = qJD(5) * t523;
	t542 = t527 * t549;
	t540 = t526 * t542;
	t538 = t526 * t557 - t529 * t536;
	t516 = -t544 + t550;
	t495 = -t516 * t531 + t538 * t533;
	t496 = t516 * t533 + t538 * t531;
	t502 = t533 * t558 + (-t525 * t531 + t528 * t554) * t527;
	t489 = -t510 * t554 - t511 * t531 + t533 * t541 + (t513 * t555 + t531 * t546 - t559) * qJD(3);
	t512 = -t527 * t528 * t526 + t530 * t529;
	t509 = t514 * qJD(1);
	t508 = t513 * qJD(1);
	t506 = t526 * t536 + t529 * t557;
	t499 = -t508 * t526 + t529 * t542;
	t498 = t503 * qJD(3);
	t497 = t502 * qJD(3);
	t488 = -t509 * t533 + (t508 * t529 + t540) * t531 + t495 * qJD(3);
	t487 = t496 * qJD(3) - t508 * t554 - t509 * t531 - t533 * t540;
	t486 = t488 * t523 + t499 * t522 + (-t496 * t522 + t506 * t523) * qJD(5);
	t485 = -t488 * t522 + t499 * t523 + (-t496 * t523 - t506 * t522) * qJD(5);
	t1 = [t571, 0, -t487 * t523 - t495 * t548, 0, t485, 0; t486, 0, t489 * t523 + t548 * t562, 0, t572, 0; 0, 0, -t498 * t523 - t502 * t548, 0, -t497 * t522 + (-t503 * t523 - t512 * t522) * qJD(5), 0; -t572, 0, t487 * t522 - t495 * t547, 0, -t486, 0; t485, 0, -t489 * t522 + t547 * t562, 0, t571, 0; 0, 0, t498 * t522 - t502 * t547, 0, -t497 * t523 + (t503 * t522 - t512 * t523) * qJD(5), 0; t489, 0, t488, 0, 0, 0; t487, 0, -t491, 0, 0, 0; 0, 0, t497, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:08
	% EndTime: 2019-10-10 01:04:09
	% DurationCPUTime: 1.31s
	% Computational Cost: add. (935->99), mult. (2579->181), div. (0->0), fcn. (2964->14), ass. (0->90)
	t746 = cos(pkin(6));
	t741 = sin(pkin(12));
	t752 = cos(qJ(1));
	t784 = t752 * t741;
	t744 = cos(pkin(12));
	t749 = sin(qJ(1));
	t785 = t749 * t744;
	t762 = t746 * t785 + t784;
	t725 = t762 * qJD(1);
	t783 = t752 * t744;
	t786 = t749 * t741;
	t730 = -t746 * t786 + t783;
	t726 = t730 * qJD(1);
	t742 = sin(pkin(7));
	t745 = cos(pkin(7));
	t748 = sin(qJ(3));
	t751 = cos(qJ(3));
	t743 = sin(pkin(6));
	t790 = t743 * t749;
	t770 = qJD(1) * t790;
	t763 = t746 * t783 - t786;
	t760 = t763 * t745;
	t789 = t743 * t752;
	t772 = t742 * t789;
	t764 = -t760 + t772;
	t729 = t746 * t784 + t785;
	t796 = t729 * t748;
	t798 = t764 * t751 + t796;
	t694 = t798 * qJD(3) - (-t725 * t745 + t742 * t770) * t748 - t726 * t751;
	t795 = t729 * t751;
	t808 = t764 * t748;
	t707 = -t795 + t808;
	t717 = t763 * t742 + t745 * t789;
	t740 = pkin(13) + qJ(5);
	t738 = sin(t740);
	t739 = cos(t740);
	t697 = t707 * t738 - t717 * t739;
	t714 = t725 * t742 + t745 * t770;
	t688 = -t697 * qJD(5) + t694 * t739 - t714 * t738;
	t782 = qJD(1) * t751;
	t791 = t743 * t742;
	t769 = t782 * t791;
	t787 = t745 * t751;
	t788 = t745 * t748;
	t693 = (-t763 * t788 - t795) * qJD(3) - t725 * t787 + t749 * t769 - (-qJD(3) * t772 + t726) * t748;
	t747 = sin(qJ(6));
	t750 = cos(qJ(6));
	t815 = -t688 * t747 + t693 * t750;
	t814 = t688 * t750 + t693 * t747;
	t699 = t707 * t739 + t717 * t738;
	t813 = t699 * t747;
	t812 = t699 * t750;
	t686 = t699 * qJD(5) + t694 * t738 + t714 * t739;
	t801 = -t742 * t790 + t762 * t745;
	t708 = t730 * t748 + t801 * t751;
	t792 = t742 * t746;
	t802 = (-t741 * t748 + t744 * t787) * t743 + t751 * t792;
	t781 = qJD(5) * t738;
	t780 = qJD(5) * t739;
	t779 = qJD(6) * t739;
	t778 = qJD(6) * t747;
	t777 = qJD(6) * t750;
	t724 = t729 * qJD(1);
	t690 = qJD(1) * t808 - t708 * qJD(3) - t724 * t751;
	t768 = t708 * t779 + t690;
	t704 = t751 * t772 - t763 * t787 + t796;
	t767 = t704 * t779 - t694;
	t711 = t802 * qJD(3);
	t766 = -t779 * t802 + t711;
	t709 = t730 * t751 - t748 * t801;
	t719 = t762 * t742 + t745 * t790;
	t701 = t709 * t739 + t719 * t738;
	t716 = t748 * t792 + (t741 * t751 + t744 * t788) * t743;
	t727 = -t744 * t791 + t746 * t745;
	t703 = t716 * t739 + t727 * t738;
	t702 = -t716 * t738 + t727 * t739;
	t689 = t709 * qJD(3) - t724 * t748 - t752 * t769 + t760 * t782;
	t757 = qJD(6) * t709 - t689 * t739 + t708 * t781;
	t756 = -qJD(6) * t707 + t693 * t739 + t704 * t781;
	t712 = t716 * qJD(3);
	t755 = qJD(6) * t716 - t712 * t739 - t781 * t802;
	t753 = t717 * qJD(1);
	t700 = -t709 * t738 + t719 * t739;
	t696 = t702 * qJD(5) + t711 * t739;
	t695 = -t703 * qJD(5) - t711 * t738;
	t685 = t690 * t739 - t709 * t781 + t719 * t780 + t738 * t753;
	t684 = t701 * qJD(5) + t690 * t738 - t739 * t753;
	t683 = t685 * t750 + t689 * t747 + (-t701 * t747 + t708 * t750) * qJD(6);
	t682 = -t685 * t747 + t689 * t750 + (-t701 * t750 - t708 * t747) * qJD(6);
	t1 = [(-t750 * t798 - t813) * qJD(6) + t814, 0, t768 * t747 + t757 * t750, 0, -t684 * t750 - t700 * t778, t682; t683, 0, t767 * t747 + t756 * t750, 0, t686 * t750 - t697 * t778, (-t704 * t747 + t812) * qJD(6) - t815; 0, 0, t766 * t747 + t755 * t750, 0, t695 * t750 - t702 * t778, -t696 * t747 + t712 * t750 + (-t703 * t750 + t747 * t802) * qJD(6); (t747 * t798 - t812) * qJD(6) + t815, 0, -t757 * t747 + t768 * t750, 0, t684 * t747 - t700 * t777, -t683; t682, 0, -t756 * t747 + t767 * t750, 0, -t686 * t747 - t697 * t777, (-t704 * t750 - t813) * qJD(6) + t814; 0, 0, -t755 * t747 + t766 * t750, 0, -t695 * t747 - t702 * t777, -t696 * t750 - t712 * t747 + (t703 * t747 + t750 * t802) * qJD(6); t686, 0, -t689 * t738 - t708 * t780, 0, t685, 0; t684, 0, t693 * t738 - t704 * t780, 0, -t688, 0; 0, 0, -t712 * t738 + t780 * t802, 0, t696, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end