% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPRR13_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR13_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
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
	% StartTime: 2019-10-10 01:07:45
	% EndTime: 2019-10-10 01:07:45
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
	% StartTime: 2019-10-10 01:07:45
	% EndTime: 2019-10-10 01:07:46
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
	% StartTime: 2019-10-10 01:07:46
	% EndTime: 2019-10-10 01:07:47
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (80->30), mult. (294->64), div. (0->0), fcn. (310->10), ass. (0->38)
	t389 = cos(pkin(6));
	t384 = sin(pkin(12));
	t393 = cos(qJ(1));
	t406 = t393 * t384;
	t387 = cos(pkin(12));
	t391 = sin(qJ(1));
	t407 = t391 * t387;
	t397 = t389 * t407 + t406;
	t375 = t397 * qJD(1);
	t408 = t391 * t384;
	t403 = t389 * t408;
	t404 = qJD(1) * t393;
	t376 = -qJD(1) * t403 + t387 * t404;
	t378 = t389 * t406 + t407;
	t388 = cos(pkin(7));
	t390 = sin(qJ(3));
	t392 = cos(qJ(3));
	t405 = t393 * t387;
	t377 = -t389 * t405 + t408;
	t385 = sin(pkin(7));
	t386 = sin(pkin(6));
	t399 = t385 * t386 * t393 + t377 * t388;
	t410 = t386 * t391;
	t402 = qJD(1) * t410;
	t400 = t385 * t402;
	t416 = (-t378 * t392 + t399 * t390) * qJD(3) + (-t375 * t388 + t400) * t392 - t376 * t390;
	t411 = t385 * t389;
	t409 = t388 * t390;
	t401 = t386 * t404;
	t398 = t385 * t410 - t388 * t397;
	t373 = t377 * qJD(1);
	t396 = t373 * t388 + t385 * t401;
	t394 = -t375 * t409 + t376 * t392 + t390 * t400 + (-t378 * t390 - t399 * t392) * qJD(3);
	t380 = -t403 + t405;
	t374 = t378 * qJD(1);
	t372 = -t374 * t392 + t396 * t390 + (-t380 * t390 + t398 * t392) * qJD(3);
	t371 = -t374 * t390 - t396 * t392 + (t380 * t392 + t398 * t390) * qJD(3);
	t1 = [-t375 * t385 - t388 * t402, 0, 0, 0, 0, 0; -t373 * t385 + t388 * t401, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t394, 0, t371, 0, 0, 0; -t372, 0, -t416, 0, 0, 0; 0, 0, (t390 * t411 + (t384 * t392 + t387 * t409) * t386) * qJD(3), 0, 0, 0; t416, 0, t372, 0, 0, 0; t371, 0, t394, 0, 0, 0; 0, 0, (t392 * t411 + (t387 * t388 * t392 - t384 * t390) * t386) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:47
	% EndTime: 2019-10-10 01:07:48
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (278->59), mult. (936->115), div. (0->0), fcn. (1042->12), ass. (0->63)
	t503 = cos(pkin(6));
	t498 = sin(pkin(12));
	t509 = cos(qJ(1));
	t528 = t509 * t498;
	t501 = cos(pkin(12));
	t506 = sin(qJ(1));
	t529 = t506 * t501;
	t513 = t503 * t529 + t528;
	t487 = t513 * qJD(1);
	t530 = t506 * t498;
	t521 = t503 * t530;
	t526 = qJD(1) * t509;
	t488 = -qJD(1) * t521 + t501 * t526;
	t527 = t509 * t501;
	t490 = -t503 * t527 + t530;
	t505 = sin(qJ(3));
	t508 = cos(qJ(3));
	t499 = sin(pkin(7));
	t500 = sin(pkin(6));
	t534 = t500 * t506;
	t520 = qJD(1) * t534;
	t518 = t499 * t520;
	t533 = t500 * t509;
	t523 = t499 * t533;
	t502 = cos(pkin(7));
	t531 = t502 * t508;
	t532 = t502 * t505;
	t491 = t503 * t528 + t529;
	t537 = t491 * t508;
	t469 = (t490 * t532 - t537) * qJD(3) - t487 * t531 + t508 * t518 + (qJD(3) * t523 - t488) * t505;
	t516 = t490 * t502 + t523;
	t470 = t491 * t505 + t516 * t508;
	t479 = t487 * t499 + t502 * t520;
	t482 = -t490 * t499 + t502 * t533;
	t504 = sin(qJ(5));
	t507 = cos(qJ(5));
	t551 = (t470 * t504 - t482 * t507) * qJD(5) + t469 * t507 + t479 * t504;
	t548 = t470 * qJD(3) - (-t487 * t502 + t518) * t505 - t488 * t508;
	t547 = t469 * t504 - t479 * t507 + (-t470 * t507 - t482 * t504) * qJD(5);
	t535 = t499 * t503;
	t542 = (-t498 * t505 + t501 * t531) * t500 + t508 * t535;
	t493 = -t521 + t527;
	t515 = t499 * t534 - t502 * t513;
	t540 = -t493 * t505 + t515 * t508;
	t525 = qJD(5) * t504;
	t524 = qJD(5) * t507;
	t519 = t500 * t526;
	t485 = t490 * qJD(1);
	t512 = t485 * t502 + t499 * t519;
	t474 = t493 * t508 + t515 * t505;
	t481 = t505 * t535 + (t498 * t508 + t501 * t532) * t500;
	t489 = -t500 * t501 * t499 + t503 * t502;
	t486 = t491 * qJD(1);
	t484 = t499 * t513 + t502 * t534;
	t477 = -t485 * t499 + t502 * t519;
	t476 = t481 * qJD(3);
	t475 = t542 * qJD(3);
	t471 = -t516 * t505 + t537;
	t466 = t540 * qJD(3) - t486 * t508 + t512 * t505;
	t465 = t474 * qJD(3) - t486 * t505 - t512 * t508;
	t464 = t465 * t504 + t477 * t507 + (-t484 * t504 - t507 * t540) * qJD(5);
	t463 = t465 * t507 - t477 * t504 + (-t484 * t507 + t504 * t540) * qJD(5);
	t1 = [t547, 0, t466 * t504 + t474 * t524, 0, t463, 0; t464, 0, t471 * t524 - t504 * t548, 0, -t551, 0; 0, 0, t475 * t504 + t481 * t524, 0, t476 * t507 + (-t489 * t507 + t504 * t542) * qJD(5), 0; t551, 0, t466 * t507 - t474 * t525, 0, -t464, 0; t463, 0, -t471 * t525 - t507 * t548, 0, t547, 0; 0, 0, t475 * t507 - t481 * t525, 0, -t476 * t504 + (t489 * t504 + t507 * t542) * qJD(5), 0; t548, 0, -t465, 0, 0, 0; t466, 0, t469, 0, 0, 0; 0, 0, -t476, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:52
	% EndTime: 2019-10-10 01:07:53
	% DurationCPUTime: 1.24s
	% Computational Cost: add. (800->107), mult. (2579->193), div. (0->0), fcn. (2964->14), ass. (0->98)
	t720 = cos(pkin(6));
	t715 = sin(pkin(12));
	t728 = cos(qJ(1));
	t757 = t728 * t715;
	t718 = cos(pkin(12));
	t724 = sin(qJ(1));
	t758 = t724 * t718;
	t734 = t720 * t758 + t757;
	t700 = t734 * qJD(1);
	t759 = t724 * t715;
	t745 = t720 * t759;
	t755 = qJD(1) * t728;
	t701 = -qJD(1) * t745 + t718 * t755;
	t719 = cos(pkin(7));
	t723 = sin(qJ(3));
	t727 = cos(qJ(3));
	t716 = sin(pkin(7));
	t717 = sin(pkin(6));
	t763 = t717 * t724;
	t744 = qJD(1) * t763;
	t742 = t716 * t744;
	t756 = t728 * t718;
	t703 = -t720 * t756 + t759;
	t762 = t717 * t728;
	t748 = t716 * t762;
	t736 = t703 * t719 + t748;
	t704 = t720 * t757 + t758;
	t768 = t704 * t723;
	t770 = t736 * t727 + t768;
	t668 = t770 * qJD(3) - (-t700 * t719 + t742) * t723 - t701 * t727;
	t721 = sin(qJ(6));
	t783 = t668 * t721;
	t725 = cos(qJ(6));
	t782 = t668 * t725;
	t767 = t704 * t727;
	t681 = t736 * t723 - t767;
	t781 = t681 * t721;
	t780 = t681 * t725;
	t760 = t719 * t727;
	t761 = t719 * t723;
	t667 = (t703 * t761 - t767) * qJD(3) - t700 * t760 + t727 * t742 - (-qJD(3) * t748 + t701) * t723;
	t688 = t700 * t716 + t719 * t744;
	t722 = sin(qJ(5));
	t726 = cos(qJ(5));
	t779 = t667 * t722 - t688 * t726;
	t778 = -t667 * t726 - t688 * t722;
	t691 = -t703 * t716 + t719 * t762;
	t775 = t691 * t722;
	t774 = t691 * t726;
	t764 = t716 * t720;
	t773 = (-t715 * t723 + t718 * t760) * t717 + t727 * t764;
	t706 = -t745 + t756;
	t766 = t706 * t723;
	t754 = qJD(5) * t722;
	t753 = qJD(5) * t726;
	t752 = qJD(6) * t721;
	t751 = qJD(6) * t722;
	t750 = qJD(6) * t725;
	t749 = t716 * t763;
	t743 = t717 * t755;
	t741 = t716 * t743;
	t735 = -t719 * t734 + t749;
	t683 = t706 * t727 + t735 * t723;
	t698 = t703 * qJD(1);
	t699 = t704 * qJD(1);
	t663 = t683 * qJD(3) - t698 * t760 - t699 * t723 - t727 * t741;
	t740 = -t683 * t751 - t663;
	t739 = t681 * t751 + t667;
	t690 = t723 * t764 + (t715 * t727 + t718 * t761) * t717;
	t685 = t690 * qJD(3);
	t738 = -t690 * t751 - t685;
	t678 = t703 * t760 + t727 * t748 + t768;
	t672 = t678 * t726 + t775;
	t673 = t678 * t722 - t774;
	t671 = -t722 * t770 + t774;
	t682 = -t727 * t749 + t734 * t760 + t766;
	t693 = t716 * t734 + t719 * t763;
	t674 = t682 * t726 - t693 * t722;
	t675 = t682 * t722 + t693 * t726;
	t702 = -t717 * t718 * t716 + t720 * t719;
	t676 = -t702 * t722 - t726 * t773;
	t677 = t702 * t726 - t722 * t773;
	t664 = -t699 * t727 + (t698 * t719 + t741) * t723 + (t735 * t727 - t766) * qJD(3);
	t731 = -qJD(6) * t682 + t664 * t722 + t683 * t753;
	t730 = -qJD(6) * t678 - t668 * t722 - t681 * t753;
	t684 = t773 * qJD(3);
	t729 = qJD(6) * t773 + t684 * t722 + t690 * t753;
	t686 = -t698 * t716 + t719 * t743;
	t670 = -t677 * qJD(5) + t685 * t726;
	t669 = t676 * qJD(5) + t685 * t722;
	t661 = t672 * qJD(5) - t779;
	t660 = -t673 * qJD(5) + t778;
	t659 = (-t726 * t770 - t775) * qJD(5) + t779;
	t658 = t674 * qJD(5) + t663 * t722 + t686 * t726;
	t657 = t675 * qJD(5) - t663 * t726 + t686 * t722;
	t656 = t658 * t725 + t664 * t721 + (-t675 * t721 + t683 * t725) * qJD(6);
	t655 = -t658 * t721 + t664 * t725 + (-t675 * t725 - t683 * t721) * qJD(6);
	t1 = [t659 * t725 + t783 + (-t671 * t721 + t780) * qJD(6), 0, t740 * t721 + t731 * t725, 0, -t657 * t725 - t674 * t752, t655; t656, 0, t739 * t721 + t730 * t725, 0, t660 * t725 - t672 * t752, -t661 * t721 - t782 + (-t673 * t725 + t781) * qJD(6); 0, 0, t738 * t721 + t729 * t725, 0, t670 * t725 - t676 * t752, -t669 * t721 + t684 * t725 + (-t677 * t725 - t690 * t721) * qJD(6); -t659 * t721 + t782 + (-t671 * t725 - t781) * qJD(6), 0, -t731 * t721 + t740 * t725, 0, t657 * t721 - t674 * t750, -t656; t655, 0, -t730 * t721 + t739 * t725, 0, -t660 * t721 - t672 * t750, -t661 * t725 + t783 + (t673 * t721 + t780) * qJD(6); 0, 0, -t729 * t721 + t738 * t725, 0, -t670 * t721 - t676 * t750, -t669 * t725 - t684 * t721 + (t677 * t721 - t690 * t725) * qJD(6); t671 * qJD(5) + t778, 0, -t664 * t726 + t683 * t754, 0, t658, 0; t657, 0, t668 * t726 - t681 * t754, 0, t661, 0; 0, 0, -t684 * t726 + t690 * t754, 0, t669, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end