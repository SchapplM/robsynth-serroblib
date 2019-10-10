% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR9
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
% Datum: 2019-10-10 01:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRPR9_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR9_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
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
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:32
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
	% StartTime: 2019-10-10 01:37:33
	% EndTime: 2019-10-10 01:37:33
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
	% StartTime: 2019-10-10 01:37:35
	% EndTime: 2019-10-10 01:37:35
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
	% StartTime: 2019-10-10 01:37:35
	% EndTime: 2019-10-10 01:37:35
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (328->60), mult. (936->116), div. (0->0), fcn. (1042->12), ass. (0->64)
	t534 = cos(pkin(6));
	t529 = sin(pkin(12));
	t538 = cos(qJ(1));
	t555 = t538 * t529;
	t532 = cos(pkin(12));
	t536 = sin(qJ(1));
	t556 = t536 * t532;
	t540 = t534 * t556 + t555;
	t514 = t540 * qJD(1);
	t557 = t536 * t529;
	t548 = t534 * t557;
	t553 = qJD(1) * t538;
	t515 = -qJD(1) * t548 + t532 * t553;
	t533 = cos(pkin(7));
	t535 = sin(qJ(3));
	t537 = cos(qJ(3));
	t530 = sin(pkin(7));
	t531 = sin(pkin(6));
	t561 = t531 * t536;
	t547 = qJD(1) * t561;
	t545 = t530 * t547;
	t518 = t534 * t555 + t556;
	t554 = t538 * t532;
	t517 = -t534 * t554 + t557;
	t560 = t531 * t538;
	t550 = t530 * t560;
	t543 = t517 * t533 + t550;
	t566 = t518 * t535 + t543 * t537;
	t495 = t566 * qJD(3) - (-t514 * t533 + t545) * t535 - t515 * t537;
	t563 = t518 * t537;
	t498 = t543 * t535 - t563;
	t505 = t514 * t530 + t533 * t547;
	t508 = -t517 * t530 + t533 * t560;
	t528 = qJ(4) + pkin(13);
	t526 = sin(t528);
	t527 = cos(t528);
	t576 = t495 * t526 + (t498 * t527 + t508 * t526) * qJD(4) + t505 * t527;
	t575 = t495 * t527 - t505 * t526 + (-t498 * t526 + t508 * t527) * qJD(4);
	t559 = t533 * t535;
	t562 = t530 * t534;
	t507 = (t529 * t537 + t532 * t559) * t531 + t535 * t562;
	t558 = t533 * t537;
	t552 = qJD(4) * t526;
	t551 = qJD(4) * t527;
	t546 = t531 * t553;
	t544 = t530 * t546;
	t542 = t530 * t561 - t533 * t540;
	t520 = -t548 + t554;
	t499 = -t520 * t535 + t542 * t537;
	t500 = t520 * t537 + t542 * t535;
	t506 = t537 * t562 + (-t529 * t535 + t532 * t558) * t531;
	t493 = -t514 * t558 - t515 * t535 + t537 * t545 + (t517 * t559 + t535 * t550 - t563) * qJD(3);
	t516 = -t531 * t532 * t530 + t534 * t533;
	t513 = t518 * qJD(1);
	t512 = t517 * qJD(1);
	t510 = t530 * t540 + t533 * t561;
	t503 = -t512 * t530 + t533 * t546;
	t502 = t507 * qJD(3);
	t501 = t506 * qJD(3);
	t492 = -t513 * t537 + (t512 * t533 + t544) * t535 + t499 * qJD(3);
	t491 = t500 * qJD(3) - t512 * t558 - t513 * t535 - t537 * t544;
	t490 = t492 * t527 + t503 * t526 + (-t500 * t526 + t510 * t527) * qJD(4);
	t489 = -t492 * t526 + t503 * t527 + (-t500 * t527 - t510 * t526) * qJD(4);
	t1 = [t575, 0, -t491 * t527 - t499 * t552, t489, 0, 0; t490, 0, t493 * t527 + t552 * t566, t576, 0, 0; 0, 0, -t502 * t527 - t506 * t552, -t501 * t526 + (-t507 * t527 - t516 * t526) * qJD(4), 0, 0; -t576, 0, t491 * t526 - t499 * t551, -t490, 0, 0; t489, 0, -t493 * t526 + t551 * t566, t575, 0, 0; 0, 0, t502 * t526 - t506 * t551, -t501 * t527 + (t507 * t526 - t516 * t527) * qJD(4), 0, 0; t493, 0, t492, 0, 0, 0; t491, 0, -t495, 0, 0, 0; 0, 0, t501, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:40
	% EndTime: 2019-10-10 01:37:41
	% DurationCPUTime: 1.30s
	% Computational Cost: add. (935->99), mult. (2579->181), div. (0->0), fcn. (2964->14), ass. (0->90)
	t752 = cos(pkin(6));
	t747 = sin(pkin(12));
	t758 = cos(qJ(1));
	t790 = t758 * t747;
	t750 = cos(pkin(12));
	t755 = sin(qJ(1));
	t791 = t755 * t750;
	t768 = t752 * t791 + t790;
	t731 = t768 * qJD(1);
	t789 = t758 * t750;
	t792 = t755 * t747;
	t736 = -t752 * t792 + t789;
	t732 = t736 * qJD(1);
	t748 = sin(pkin(7));
	t751 = cos(pkin(7));
	t754 = sin(qJ(3));
	t757 = cos(qJ(3));
	t749 = sin(pkin(6));
	t796 = t749 * t755;
	t776 = qJD(1) * t796;
	t769 = t752 * t789 - t792;
	t766 = t769 * t751;
	t795 = t749 * t758;
	t778 = t748 * t795;
	t770 = -t766 + t778;
	t735 = t752 * t790 + t791;
	t802 = t735 * t754;
	t804 = t770 * t757 + t802;
	t700 = t804 * qJD(3) - (-t731 * t751 + t748 * t776) * t754 - t732 * t757;
	t801 = t735 * t757;
	t814 = t770 * t754;
	t713 = -t801 + t814;
	t723 = t769 * t748 + t751 * t795;
	t746 = qJ(4) + pkin(13);
	t744 = sin(t746);
	t745 = cos(t746);
	t703 = t713 * t744 - t723 * t745;
	t720 = t731 * t748 + t751 * t776;
	t694 = -t703 * qJD(4) + t700 * t745 - t720 * t744;
	t788 = qJD(1) * t757;
	t797 = t749 * t748;
	t775 = t788 * t797;
	t793 = t751 * t757;
	t794 = t751 * t754;
	t699 = (-t769 * t794 - t801) * qJD(3) - t731 * t793 + t755 * t775 - (-qJD(3) * t778 + t732) * t754;
	t753 = sin(qJ(6));
	t756 = cos(qJ(6));
	t821 = -t694 * t753 + t699 * t756;
	t820 = t694 * t756 + t699 * t753;
	t705 = t713 * t745 + t723 * t744;
	t819 = t705 * t753;
	t818 = t705 * t756;
	t692 = qJD(4) * t705 + t700 * t744 + t720 * t745;
	t807 = -t748 * t796 + t768 * t751;
	t714 = t736 * t754 + t807 * t757;
	t798 = t748 * t752;
	t808 = (-t747 * t754 + t750 * t793) * t749 + t757 * t798;
	t787 = qJD(4) * t744;
	t786 = qJD(4) * t745;
	t785 = qJD(6) * t745;
	t784 = qJD(6) * t753;
	t783 = qJD(6) * t756;
	t730 = t735 * qJD(1);
	t696 = qJD(1) * t814 - t714 * qJD(3) - t730 * t757;
	t774 = t714 * t785 + t696;
	t710 = t757 * t778 - t769 * t793 + t802;
	t773 = t710 * t785 - t700;
	t717 = t808 * qJD(3);
	t772 = -t785 * t808 + t717;
	t715 = t736 * t757 - t754 * t807;
	t725 = t768 * t748 + t751 * t796;
	t707 = t715 * t745 + t725 * t744;
	t722 = t754 * t798 + (t747 * t757 + t750 * t794) * t749;
	t733 = -t750 * t797 + t752 * t751;
	t709 = t722 * t745 + t733 * t744;
	t708 = -t722 * t744 + t733 * t745;
	t695 = t715 * qJD(3) - t730 * t754 - t758 * t775 + t766 * t788;
	t763 = qJD(6) * t715 - t695 * t745 + t714 * t787;
	t762 = -qJD(6) * t713 + t699 * t745 + t710 * t787;
	t718 = t722 * qJD(3);
	t761 = qJD(6) * t722 - t718 * t745 - t787 * t808;
	t759 = t723 * qJD(1);
	t706 = -t715 * t744 + t725 * t745;
	t702 = t708 * qJD(4) + t717 * t745;
	t701 = -t709 * qJD(4) - t717 * t744;
	t691 = t696 * t745 - t715 * t787 + t725 * t786 + t744 * t759;
	t690 = t707 * qJD(4) + t696 * t744 - t745 * t759;
	t689 = t691 * t756 + t695 * t753 + (-t707 * t753 + t714 * t756) * qJD(6);
	t688 = -t691 * t753 + t695 * t756 + (-t707 * t756 - t714 * t753) * qJD(6);
	t1 = [(-t756 * t804 - t819) * qJD(6) + t820, 0, t753 * t774 + t756 * t763, -t690 * t756 - t706 * t784, 0, t688; t689, 0, t753 * t773 + t756 * t762, t692 * t756 - t703 * t784, 0, (-t710 * t753 + t818) * qJD(6) - t821; 0, 0, t753 * t772 + t756 * t761, t701 * t756 - t708 * t784, 0, -t702 * t753 + t718 * t756 + (-t709 * t756 + t753 * t808) * qJD(6); (t753 * t804 - t818) * qJD(6) + t821, 0, -t753 * t763 + t756 * t774, t690 * t753 - t706 * t783, 0, -t689; t688, 0, -t753 * t762 + t756 * t773, -t692 * t753 - t703 * t783, 0, (-t710 * t756 - t819) * qJD(6) + t820; 0, 0, -t753 * t761 + t756 * t772, -t701 * t753 - t708 * t783, 0, -t702 * t756 - t718 * t753 + (t709 * t753 + t756 * t808) * qJD(6); t692, 0, -t695 * t744 - t714 * t786, t691, 0, 0; t690, 0, t699 * t744 - t710 * t786, -t694, 0, 0; 0, 0, -t718 * t744 + t786 * t808, t702, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end