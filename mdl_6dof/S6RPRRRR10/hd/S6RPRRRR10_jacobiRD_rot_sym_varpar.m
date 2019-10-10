% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRRR10_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR10_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
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
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->5), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t128 = cos(pkin(6));
	t129 = sin(qJ(1));
	t133 = t128 * t129;
	t130 = cos(qJ(1));
	t132 = t128 * t130;
	t131 = qJD(1) * sin(pkin(6));
	t127 = cos(pkin(13));
	t125 = sin(pkin(13));
	t1 = [(t125 * t133 - t127 * t130) * qJD(1), 0, 0, 0, 0, 0; (-t125 * t132 - t127 * t129) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; (t125 * t130 + t127 * t133) * qJD(1), 0, 0, 0, 0, 0; (t125 * t129 - t127 * t132) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t129 * t131, 0, 0, 0, 0, 0; t130 * t131, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:31
	% EndTime: 2019-10-10 09:11:31
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (80->29), mult. (294->61), div. (0->0), fcn. (310->10), ass. (0->38)
	t298 = cos(pkin(6));
	t296 = cos(pkin(13));
	t302 = cos(qJ(1));
	t314 = t302 * t296;
	t293 = sin(pkin(13));
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
	% StartTime: 2019-10-10 09:11:33
	% EndTime: 2019-10-10 09:11:33
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (278->59), mult. (936->116), div. (0->0), fcn. (1042->12), ass. (0->63)
	t506 = cos(pkin(6));
	t501 = sin(pkin(13));
	t512 = cos(qJ(1));
	t529 = t512 * t501;
	t504 = cos(pkin(13));
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
	% StartTime: 2019-10-10 09:11:33
	% EndTime: 2019-10-10 09:11:34
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (482->56), mult. (1246->105), div. (0->0), fcn. (1390->12), ass. (0->66)
	t557 = qJ(4) + qJ(5);
	t554 = sin(t557);
	t555 = cos(t557);
	t563 = cos(pkin(6));
	t558 = sin(pkin(13));
	t567 = cos(qJ(1));
	t590 = t567 * t558;
	t561 = cos(pkin(13));
	t565 = sin(qJ(1));
	t591 = t565 * t561;
	t571 = t563 * t591 + t590;
	t542 = t571 * qJD(1);
	t556 = qJD(4) + qJD(5);
	t559 = sin(pkin(7));
	t562 = cos(pkin(7));
	t564 = sin(qJ(3));
	t589 = t567 * t561;
	t592 = t565 * t558;
	t545 = -t563 * t589 + t592;
	t560 = sin(pkin(6));
	t595 = t560 * t567;
	t587 = t559 * t595;
	t574 = t545 * t562 + t587;
	t596 = t560 * t565;
	t584 = qJD(1) * t596;
	t546 = t563 * t590 + t591;
	t566 = cos(qJ(3));
	t600 = t546 * t566;
	t608 = -t542 * t559 - (t574 * t564 - t600) * t556 - t562 * t584;
	t585 = t563 * t592;
	t588 = qJD(1) * t567;
	t543 = -qJD(1) * t585 + t561 * t588;
	t576 = t559 * t584;
	t603 = t546 * t564 + t574 * t566;
	t520 = -t603 * qJD(3) + (-t542 * t562 + t576) * t564 + t543 * t566;
	t609 = -t520 + (-t545 * t559 + t562 * t595) * t556;
	t515 = t554 * t609 - t555 * t608;
	t516 = t554 * t608 + t555 * t609;
	t540 = t545 * qJD(1);
	t548 = -t585 + t589;
	t573 = t559 * t596 - t562 * t571;
	t569 = t548 * t566 + t573 * t564;
	t583 = t560 * t588;
	t606 = -t540 * t559 - t569 * t556 + t562 * t583;
	t594 = t562 * t564;
	t597 = t559 * t563;
	t535 = (t558 * t566 + t561 * t594) * t560 + t564 * t597;
	t599 = t554 * t556;
	t598 = t555 * t556;
	t593 = t562 * t566;
	t528 = -t548 * t564 + t573 * t566;
	t541 = t546 * qJD(1);
	t575 = t559 * t583;
	t518 = -t541 * t566 + (t540 * t562 + t575) * t564 + t528 * qJD(3);
	t582 = (t559 * t571 + t562 * t596) * t556 + t518;
	t534 = t566 * t597 + (-t558 * t564 + t561 * t593) * t560;
	t529 = t534 * qJD(3);
	t577 = -(-t560 * t561 * t559 + t563 * t562) * t556 - t529;
	t519 = -t542 * t593 - t543 * t564 + t566 * t576 + (t545 * t594 + t564 * t587 - t600) * qJD(3);
	t530 = t535 * qJD(3);
	t523 = t535 * t599 + t577 * t555;
	t522 = -t535 * t598 + t577 * t554;
	t517 = t569 * qJD(3) - t540 * t593 - t541 * t564 - t566 * t575;
	t514 = t606 * t554 + t582 * t555;
	t513 = -t582 * t554 + t606 * t555;
	t1 = [t516, 0, -t517 * t555 - t528 * t599, t513, t513, 0; t514, 0, t519 * t555 + t599 * t603, t515, t515, 0; 0, 0, -t530 * t555 - t534 * t599, t522, t522, 0; -t515, 0, t517 * t554 - t528 * t598, -t514, -t514, 0; t513, 0, -t519 * t554 + t598 * t603, t516, t516, 0; 0, 0, t530 * t554 - t534 * t598, t523, t523, 0; t519, 0, t518, 0, 0, 0; t517, 0, t520, 0, 0, 0; 0, 0, t529, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:39
	% EndTime: 2019-10-10 09:11:41
	% DurationCPUTime: 1.46s
	% Computational Cost: add. (1232->101), mult. (3166->183), div. (0->0), fcn. (3644->14), ass. (0->99)
	t805 = qJ(4) + qJ(5);
	t802 = sin(t805);
	t803 = cos(t805);
	t811 = cos(pkin(6));
	t806 = sin(pkin(13));
	t817 = cos(qJ(1));
	t852 = t817 * t806;
	t809 = cos(pkin(13));
	t814 = sin(qJ(1));
	t853 = t814 * t809;
	t827 = t811 * t853 + t852;
	t789 = t827 * qJD(1);
	t851 = t817 * t809;
	t854 = t814 * t806;
	t794 = -t811 * t854 + t851;
	t790 = t794 * qJD(1);
	t807 = sin(pkin(7));
	t810 = cos(pkin(7));
	t813 = sin(qJ(3));
	t816 = cos(qJ(3));
	t808 = sin(pkin(6));
	t858 = t808 * t814;
	t840 = qJD(1) * t858;
	t828 = t811 * t851 - t854;
	t825 = t828 * t810;
	t857 = t808 * t817;
	t842 = t807 * t857;
	t829 = -t825 + t842;
	t793 = t811 * t852 + t853;
	t866 = t793 * t813;
	t868 = t829 * t816 + t866;
	t756 = -t868 * qJD(3) + (-t789 * t810 + t807 * t840) * t813 + t790 * t816;
	t781 = t828 * t807 + t810 * t857;
	t804 = qJD(4) + qJD(5);
	t839 = -t781 * t804 + t756;
	t865 = t793 * t816;
	t876 = t829 * t813;
	t771 = -t865 + t876;
	t875 = t771 * t804 + t789 * t807 + t810 * t840;
	t749 = t802 * t875 + t839 * t803;
	t850 = qJD(1) * t816;
	t859 = t808 * t807;
	t834 = t850 * t859;
	t855 = t810 * t816;
	t856 = t810 * t813;
	t757 = (-t828 * t856 - t865) * qJD(3) - t789 * t855 + t814 * t834 - (-qJD(3) * t842 + t790) * t813;
	t812 = sin(qJ(6));
	t815 = cos(qJ(6));
	t884 = -t749 * t812 - t757 * t815;
	t883 = -t749 * t815 + t757 * t812;
	t763 = t771 * t803 + t781 * t802;
	t882 = t763 * t812;
	t881 = t763 * t815;
	t748 = -t839 * t802 + t803 * t875;
	t871 = -t807 * t858 + t827 * t810;
	t772 = t794 * t813 + t871 * t816;
	t860 = t807 * t811;
	t872 = (-t806 * t813 + t809 * t855) * t808 + t816 * t860;
	t863 = t802 * t804;
	t862 = t803 * t804;
	t849 = qJD(6) * t803;
	t848 = qJD(6) * t812;
	t847 = qJD(6) * t815;
	t775 = t872 * qJD(3);
	t791 = -t809 * t859 + t811 * t810;
	t835 = t791 * t804 + t775;
	t788 = t793 * qJD(1);
	t754 = qJD(1) * t876 - t772 * qJD(3) - t788 * t816;
	t833 = t772 * t849 + t754;
	t768 = t816 * t842 - t828 * t855 + t866;
	t832 = t768 * t849 + t756;
	t831 = -t849 * t872 + t775;
	t773 = t794 * t816 - t813 * t871;
	t753 = t773 * qJD(3) - t788 * t813 - t817 * t834 + t825 * t850;
	t822 = qJD(6) * t773 - t753 * t803 + t772 * t863;
	t821 = -qJD(6) * t771 + t757 * t803 + t768 * t863;
	t780 = t813 * t860 + (t806 * t816 + t809 * t856) * t808;
	t776 = t780 * qJD(3);
	t820 = qJD(6) * t780 - t776 * t803 - t863 * t872;
	t818 = t781 * qJD(1);
	t783 = t827 * t807 + t810 * t858;
	t767 = t780 * t803 + t791 * t802;
	t766 = -t780 * t802 + t791 * t803;
	t765 = t773 * t803 + t783 * t802;
	t764 = -t773 * t802 + t783 * t803;
	t761 = t771 * t802 - t781 * t803;
	t760 = -t780 * t863 + t835 * t803;
	t759 = -t780 * t862 - t835 * t802;
	t752 = t759 * t815 - t766 * t848;
	t751 = -t759 * t812 - t766 * t847;
	t747 = t754 * t803 - t773 * t863 + t783 * t862 + t802 * t818;
	t746 = t773 * t862 - t803 * t818 + (t783 * t804 + t754) * t802;
	t745 = t748 * t815 - t761 * t848;
	t744 = -t748 * t812 - t761 * t847;
	t743 = -t746 * t815 - t764 * t848;
	t742 = t746 * t812 - t764 * t847;
	t741 = t747 * t815 + t753 * t812 + (-t765 * t812 + t772 * t815) * qJD(6);
	t740 = -t747 * t812 + t753 * t815 + (-t765 * t815 - t772 * t812) * qJD(6);
	t1 = [(-t815 * t868 - t882) * qJD(6) + t883, 0, t833 * t812 + t822 * t815, t743, t743, t740; t741, 0, t832 * t812 + t821 * t815, t745, t745, (-t768 * t812 + t881) * qJD(6) + t884; 0, 0, t831 * t812 + t820 * t815, t752, t752, -t760 * t812 + t776 * t815 + (-t767 * t815 + t812 * t872) * qJD(6); (t812 * t868 - t881) * qJD(6) - t884, 0, -t822 * t812 + t833 * t815, t742, t742, -t741; t740, 0, -t821 * t812 + t832 * t815, t744, t744, (-t768 * t815 - t882) * qJD(6) + t883; 0, 0, -t820 * t812 + t831 * t815, t751, t751, -t760 * t815 - t776 * t812 + (t767 * t812 + t815 * t872) * qJD(6); t748, 0, -t753 * t802 - t772 * t862, t747, t747, 0; t746, 0, t757 * t802 - t768 * t862, t749, t749, 0; 0, 0, -t776 * t802 + t862 * t872, t760, t760, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end