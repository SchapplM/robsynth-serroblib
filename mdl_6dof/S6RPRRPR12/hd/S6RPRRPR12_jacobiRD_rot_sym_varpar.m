% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:43
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRPR12_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR12_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
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
	% StartTime: 2019-10-10 01:43:26
	% EndTime: 2019-10-10 01:43:26
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
	% StartTime: 2019-10-10 01:43:26
	% EndTime: 2019-10-10 01:43:27
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
	% StartTime: 2019-10-10 01:43:28
	% EndTime: 2019-10-10 01:43:29
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
	% StartTime: 2019-10-10 01:43:31
	% EndTime: 2019-10-10 01:43:32
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (278->59), mult. (936->117), div. (0->0), fcn. (1042->12), ass. (0->62)
	t577 = sin(pkin(12));
	t582 = cos(pkin(6));
	t588 = cos(qJ(1));
	t580 = cos(pkin(12));
	t585 = sin(qJ(1));
	t604 = t585 * t580;
	t590 = t577 * t588 + t582 * t604;
	t565 = t590 * qJD(1);
	t611 = t577 * t585;
	t598 = t582 * t611;
	t603 = qJD(1) * t588;
	t566 = -qJD(1) * t598 + t580 * t603;
	t581 = cos(pkin(7));
	t584 = sin(qJ(3));
	t587 = cos(qJ(3));
	t578 = sin(pkin(7));
	t579 = sin(pkin(6));
	t609 = t579 * t585;
	t597 = qJD(1) * t609;
	t595 = t578 * t597;
	t605 = t582 * t588;
	t569 = t577 * t605 + t604;
	t568 = -t580 * t605 + t611;
	t608 = t579 * t588;
	t600 = t578 * t608;
	t593 = t568 * t581 + t600;
	t615 = t569 * t584 + t593 * t587;
	t546 = t615 * qJD(3) - (-t565 * t581 + t595) * t584 - t566 * t587;
	t612 = t569 * t587;
	t549 = t593 * t584 - t612;
	t556 = t565 * t578 + t581 * t597;
	t559 = -t568 * t578 + t581 * t608;
	t583 = sin(qJ(4));
	t586 = cos(qJ(4));
	t625 = t546 * t583 + t556 * t586 + (t549 * t586 + t559 * t583) * qJD(4);
	t624 = -t546 * t586 + t556 * t583 + (t549 * t583 - t559 * t586) * qJD(4);
	t607 = t581 * t584;
	t610 = t578 * t582;
	t558 = (t577 * t587 + t580 * t607) * t579 + t584 * t610;
	t606 = t581 * t587;
	t602 = qJD(4) * t583;
	t601 = qJD(4) * t586;
	t596 = t579 * t603;
	t594 = t578 * t596;
	t592 = t578 * t609 - t581 * t590;
	t571 = t580 * t588 - t598;
	t550 = -t571 * t584 + t592 * t587;
	t551 = t571 * t587 + t592 * t584;
	t557 = t587 * t610 + (-t577 * t584 + t580 * t606) * t579;
	t544 = -t565 * t606 - t566 * t584 + t587 * t595 + (t568 * t607 + t584 * t600 - t612) * qJD(3);
	t567 = -t578 * t579 * t580 + t581 * t582;
	t564 = t569 * qJD(1);
	t563 = t568 * qJD(1);
	t561 = t578 * t590 + t581 * t609;
	t554 = -t563 * t578 + t581 * t596;
	t553 = t558 * qJD(3);
	t552 = t557 * qJD(3);
	t543 = -t564 * t587 + (t563 * t581 + t594) * t584 + t550 * qJD(3);
	t542 = t551 * qJD(3) - t563 * t606 - t564 * t584 - t587 * t594;
	t541 = t543 * t586 + t554 * t583 + (-t551 * t583 + t561 * t586) * qJD(4);
	t540 = t543 * t583 - t554 * t586 + (t551 * t586 + t561 * t583) * qJD(4);
	t1 = [t544, 0, t543, 0, 0, 0; t542, 0, -t546, 0, 0, 0; 0, 0, t552, 0, 0, 0; t624, 0, t542 * t586 + t550 * t602, t540, 0, 0; -t541, 0, -t544 * t586 - t602 * t615, -t625, 0, 0; 0, 0, t553 * t586 + t557 * t602, t552 * t583 + (t558 * t586 + t567 * t583) * qJD(4), 0, 0; t625, 0, -t542 * t583 + t550 * t601, t541, 0, 0; t540, 0, t544 * t583 - t601 * t615, t624, 0, 0; 0, 0, -t553 * t583 + t557 * t601, t552 * t586 + (-t558 * t583 + t567 * t586) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:33
	% EndTime: 2019-10-10 01:43:34
	% DurationCPUTime: 1.30s
	% Computational Cost: add. (800->99), mult. (2579->180), div. (0->0), fcn. (2964->14), ass. (0->89)
	t724 = cos(pkin(6));
	t719 = sin(pkin(12));
	t732 = cos(qJ(1));
	t763 = t732 * t719;
	t722 = cos(pkin(12));
	t728 = sin(qJ(1));
	t764 = t728 * t722;
	t741 = t724 * t764 + t763;
	t706 = t741 * qJD(1);
	t762 = t732 * t722;
	t765 = t728 * t719;
	t711 = -t724 * t765 + t762;
	t707 = t711 * qJD(1);
	t720 = sin(pkin(7));
	t723 = cos(pkin(7));
	t727 = sin(qJ(3));
	t731 = cos(qJ(3));
	t721 = sin(pkin(6));
	t769 = t721 * t728;
	t749 = qJD(1) * t769;
	t742 = t724 * t762 - t765;
	t739 = t742 * t723;
	t768 = t721 * t732;
	t751 = t720 * t768;
	t743 = -t739 + t751;
	t710 = t724 * t763 + t764;
	t775 = t710 * t727;
	t777 = t743 * t731 + t775;
	t675 = t777 * qJD(3) - (-t706 * t723 + t720 * t749) * t727 - t707 * t731;
	t695 = t706 * t720 + t723 * t749;
	t726 = sin(qJ(4));
	t730 = cos(qJ(4));
	t774 = t710 * t731;
	t787 = t743 * t727;
	t688 = -t774 + t787;
	t698 = t742 * t720 + t723 * t768;
	t786 = t688 * t730 + t698 * t726;
	t669 = t786 * qJD(4) + t675 * t726 + t695 * t730;
	t761 = qJD(1) * t731;
	t770 = t721 * t720;
	t748 = t761 * t770;
	t766 = t723 * t731;
	t767 = t723 * t727;
	t674 = (-t742 * t767 - t774) * qJD(3) - t706 * t766 + t728 * t748 - (-qJD(3) * t751 + t707) * t727;
	t725 = sin(qJ(6));
	t729 = cos(qJ(6));
	t793 = t669 * t725 + t674 * t729;
	t792 = -t669 * t729 + t674 * t725;
	t680 = t688 * t726 - t698 * t730;
	t791 = t680 * t725;
	t790 = t680 * t729;
	t668 = t680 * qJD(4) - t675 * t730 + t695 * t726;
	t780 = -t720 * t769 + t741 * t723;
	t689 = t711 * t727 + t780 * t731;
	t771 = t720 * t724;
	t781 = (-t719 * t727 + t722 * t766) * t721 + t731 * t771;
	t760 = qJD(4) * t726;
	t759 = qJD(4) * t730;
	t758 = qJD(6) * t725;
	t757 = qJD(6) * t726;
	t756 = qJD(6) * t729;
	t705 = t710 * qJD(1);
	t671 = qJD(1) * t787 - t689 * qJD(3) - t705 * t731;
	t747 = t689 * t757 - t671;
	t685 = t731 * t751 - t742 * t766 + t775;
	t746 = t685 * t757 + t675;
	t691 = t781 * qJD(3);
	t745 = -t757 * t781 - t691;
	t690 = t711 * t731 - t727 * t780;
	t700 = t741 * t720 + t723 * t769;
	t682 = t690 * t730 + t700 * t726;
	t681 = t690 * t726 - t700 * t730;
	t697 = t727 * t771 + (t719 * t731 + t722 * t767) * t721;
	t708 = -t722 * t770 + t724 * t723;
	t684 = t697 * t730 + t708 * t726;
	t683 = t697 * t726 - t708 * t730;
	t670 = t690 * qJD(3) - t705 * t727 - t732 * t748 + t739 * t761;
	t736 = -qJD(6) * t690 - t670 * t726 - t689 * t759;
	t735 = qJD(6) * t688 + t674 * t726 - t685 * t759;
	t692 = t697 * qJD(3);
	t734 = -qJD(6) * t697 - t692 * t726 + t759 * t781;
	t693 = t698 * qJD(1);
	t677 = -t683 * qJD(4) + t691 * t730;
	t676 = t684 * qJD(4) + t691 * t726;
	t666 = -t681 * qJD(4) + t671 * t730 + t693 * t726;
	t665 = t682 * qJD(4) + t671 * t726 - t693 * t730;
	t664 = t665 * t725 + t670 * t729 + (t681 * t729 - t689 * t725) * qJD(6);
	t663 = t665 * t729 - t670 * t725 + (-t681 * t725 - t689 * t729) * qJD(6);
	t1 = [(t725 * t777 + t790) * qJD(6) + t793, 0, t736 * t725 - t747 * t729, t666 * t725 + t682 * t756, 0, t663; t664, 0, t735 * t725 - t746 * t729, t668 * t725 - t756 * t786, 0, (-t685 * t729 + t791) * qJD(6) + t792; 0, 0, t734 * t725 - t745 * t729, t677 * t725 + t684 * t756, 0, t676 * t729 - t692 * t725 + (-t683 * t725 + t729 * t781) * qJD(6); (t729 * t777 - t791) * qJD(6) - t792, 0, t747 * t725 + t736 * t729, t666 * t729 - t682 * t758, 0, -t664; t663, 0, t746 * t725 + t735 * t729, t668 * t729 + t758 * t786, 0, (t685 * t725 + t790) * qJD(6) + t793; 0, 0, t745 * t725 + t734 * t729, t677 * t729 - t684 * t758, 0, -t676 * t725 - t692 * t729 + (-t683 * t729 - t725 * t781) * qJD(6); -t668, 0, -t670 * t730 + t689 * t760, -t665, 0, 0; t666, 0, t674 * t730 + t685 * t760, t669, 0, 0; 0, 0, -t692 * t730 - t760 * t781, -t676, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end