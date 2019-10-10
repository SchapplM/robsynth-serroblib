% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR9
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
% Datum: 2019-10-10 01:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPRR9_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR9_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
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
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
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
	% StartTime: 2019-10-10 01:00:11
	% EndTime: 2019-10-10 01:00:11
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
	% StartTime: 2019-10-10 01:00:11
	% EndTime: 2019-10-10 01:00:12
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (134->37), mult. (468->76), div. (0->0), fcn. (508->12), ass. (0->49)
	t369 = sin(pkin(13));
	t373 = cos(pkin(13));
	t379 = cos(qJ(3));
	t393 = qJD(3) * t379;
	t377 = sin(qJ(3));
	t394 = qJD(3) * t377;
	t403 = t369 * t394 - t373 * t393;
	t371 = sin(pkin(7));
	t388 = t379 * t369 + t377 * t373;
	t384 = qJD(3) * t388;
	t343 = t371 * t384;
	t375 = cos(pkin(7));
	t345 = t375 * t384;
	t360 = t377 * t369 - t379 * t373;
	t346 = t360 * t371;
	t348 = t360 * t375;
	t376 = cos(pkin(6));
	t370 = sin(pkin(12));
	t380 = cos(qJ(1));
	t398 = t380 * t370;
	t374 = cos(pkin(12));
	t378 = sin(qJ(1));
	t399 = t378 * t374;
	t387 = t376 * t399 + t398;
	t352 = t387 * qJD(1);
	t400 = t378 * t370;
	t392 = t376 * t400;
	t395 = qJD(1) * t380;
	t353 = -qJD(1) * t392 + t374 * t395;
	t397 = t380 * t374;
	t354 = -t376 * t397 + t400;
	t355 = t376 * t398 + t399;
	t358 = t360 * qJD(3);
	t372 = sin(pkin(6));
	t396 = qJD(1) * t378;
	t402 = t354 * t345 + t352 * t348 - t353 * t388 + t355 * t358 + (t343 * t380 - t346 * t396) * t372;
	t342 = t403 * t371;
	t344 = t403 * t375;
	t347 = t388 * t371;
	t349 = t388 * t375;
	t350 = t354 * qJD(1);
	t351 = t355 * qJD(1);
	t357 = -t392 + t397;
	t359 = -t369 * t393 - t373 * t394;
	t401 = t387 * t344 + t350 * t349 + t351 * t360 + t357 * t359 + (-t342 * t378 + t347 * t395) * t372;
	t391 = qJD(1) * t372 * t375;
	t381 = -t354 * t344 + t352 * t349 + t353 * t360 - t355 * t359 + (-t342 * t380 - t347 * t396) * t372;
	t341 = t387 * t345 - t350 * t348 + t351 * t388 + t357 * t358 + (-t343 * t378 - t346 * t395) * t372;
	t1 = [t381, 0, t341, 0, 0, 0; t401, 0, t402, 0, 0, 0; 0, 0, -t376 * t343 + (-t345 * t374 + t358 * t370) * t372, 0, 0, 0; -t402, 0, -t401, 0, 0, 0; t341, 0, t381, 0, 0, 0; 0, 0, t376 * t342 + (t344 * t374 - t359 * t370) * t372, 0, 0, 0; -t352 * t371 - t378 * t391, 0, 0, 0, 0, 0; -t350 * t371 + t380 * t391, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:14
	% EndTime: 2019-10-10 01:00:14
	% DurationCPUTime: 0.65s
	% Computational Cost: add. (437->73), mult. (1423->140), div. (0->0), fcn. (1625->14), ass. (0->69)
	t565 = sin(pkin(13));
	t573 = sin(qJ(3));
	t602 = cos(pkin(13));
	t603 = cos(qJ(3));
	t579 = -t573 * t565 + t603 * t602;
	t611 = t579 * qJD(3);
	t567 = sin(pkin(7));
	t537 = t611 * t567;
	t570 = cos(pkin(7));
	t539 = t611 * t570;
	t580 = t603 * t565 + t573 * t602;
	t542 = t580 * t567;
	t544 = t580 * t570;
	t571 = cos(pkin(6));
	t566 = sin(pkin(12));
	t576 = cos(qJ(1));
	t596 = t576 * t566;
	t569 = cos(pkin(12));
	t574 = sin(qJ(1));
	t597 = t574 * t569;
	t582 = t571 * t597 + t596;
	t547 = t582 * qJD(1);
	t598 = t574 * t566;
	t590 = t571 * t598;
	t593 = qJD(1) * t576;
	t548 = -qJD(1) * t590 + t569 * t593;
	t595 = t576 * t569;
	t550 = -t571 * t595 + t598;
	t551 = t571 * t596 + t597;
	t568 = sin(pkin(6));
	t594 = qJD(1) * t574;
	t607 = qJD(3) * t580;
	t519 = -t550 * t539 - t547 * t544 + t548 * t579 - t551 * t607 + (-t537 * t576 + t542 * t594) * t568;
	t600 = t568 * t576;
	t523 = t542 * t600 + t550 * t544 - t551 * t579;
	t587 = t568 * t594;
	t532 = t547 * t567 + t570 * t587;
	t589 = t570 * t600;
	t533 = -t550 * t567 + t589;
	t572 = sin(qJ(5));
	t575 = cos(qJ(5));
	t610 = -t519 * t575 - t532 * t572 + (-t523 * t572 + t533 * t575) * qJD(5);
	t609 = -t519 * t572 + t532 * t575 + (t523 * t575 + t533 * t572) * qJD(5);
	t601 = t568 * t574;
	t592 = qJD(5) * t572;
	t591 = qJD(5) * t575;
	t527 = t571 * t537 + (t539 * t569 - t566 * t607) * t568;
	t538 = t567 * t607;
	t540 = t570 * t607;
	t541 = t579 * t567;
	t543 = t579 * t570;
	t518 = t538 * t600 + t550 * t540 + t541 * t587 - t547 * t543 - t548 * t580 - t551 * t611;
	t545 = t550 * qJD(1);
	t546 = t551 * qJD(1);
	t553 = -t590 + t595;
	t517 = -t582 * t539 + t545 * t544 - t546 * t579 - t553 * t607 + (t537 * t574 + t542 * t593) * t568;
	t549 = -t568 * t569 * t567 + t571 * t570;
	t535 = t567 * t582 + t570 * t601;
	t530 = qJD(1) * t589 - t545 * t567;
	t529 = t571 * t542 + (t544 * t569 + t566 * t579) * t568;
	t528 = t571 * t541 + (t543 * t569 - t566 * t580) * t568;
	t526 = -t571 * t538 + (-t540 * t569 - t566 * t611) * t568;
	t525 = t542 * t601 - t544 * t582 + t553 * t579;
	t524 = t541 * t601 - t543 * t582 - t553 * t580;
	t521 = -t541 * t600 - t550 * t543 - t551 * t580;
	t516 = t582 * t540 + t545 * t543 + t546 * t580 - t553 * t611 + (-t538 * t574 + t541 * t593) * t568;
	t515 = t517 * t575 + t530 * t572 + (-t525 * t572 + t535 * t575) * qJD(5);
	t514 = -t517 * t572 + t530 * t575 + (-t525 * t575 - t535 * t572) * qJD(5);
	t1 = [t610, 0, t516 * t575 - t524 * t592, 0, t514, 0; t515, 0, t518 * t575 - t521 * t592, 0, t609, 0; 0, 0, t526 * t575 - t528 * t592, 0, -t527 * t572 + (-t529 * t575 - t549 * t572) * qJD(5), 0; -t609, 0, -t516 * t572 - t524 * t591, 0, -t515, 0; t514, 0, -t518 * t572 - t521 * t591, 0, t610, 0; 0, 0, -t526 * t572 - t528 * t591, 0, -t527 * t575 + (t529 * t572 - t549 * t575) * qJD(5), 0; t518, 0, t517, 0, 0, 0; -t516, 0, t519, 0, 0, 0; 0, 0, t527, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:19
	% EndTime: 2019-10-10 01:00:21
	% DurationCPUTime: 1.32s
	% Computational Cost: add. (1244->108), mult. (3901->203), div. (0->0), fcn. (4592->16), ass. (0->92)
	t781 = sin(pkin(13));
	t790 = sin(qJ(3));
	t836 = cos(pkin(13));
	t837 = cos(qJ(3));
	t773 = -t837 * t781 - t790 * t836;
	t783 = sin(pkin(7));
	t758 = t773 * t783;
	t786 = cos(pkin(7));
	t760 = t773 * t786;
	t787 = cos(pkin(6));
	t785 = cos(pkin(12));
	t794 = cos(qJ(1));
	t829 = t794 * t785;
	t782 = sin(pkin(12));
	t791 = sin(qJ(1));
	t832 = t791 * t782;
	t766 = -t787 * t829 + t832;
	t830 = t794 * t782;
	t831 = t791 * t785;
	t767 = t787 * t830 + t831;
	t806 = -t790 * t781 + t837 * t836;
	t784 = sin(pkin(6));
	t834 = t784 * t794;
	t735 = -t758 * t834 - t766 * t760 - t767 * t806;
	t820 = t786 * t834;
	t748 = -t766 * t783 + t820;
	t789 = sin(qJ(5));
	t793 = cos(qJ(5));
	t724 = t735 * t789 - t748 * t793;
	t810 = t787 * t831 + t830;
	t763 = t810 * qJD(1);
	t828 = qJD(1) * t791;
	t818 = t784 * t828;
	t747 = t763 * t783 + t786 * t818;
	t850 = t806 * qJD(3);
	t753 = t850 * t783;
	t755 = t850 * t786;
	t821 = t787 * t832;
	t827 = qJD(1) * t794;
	t764 = -qJD(1) * t821 + t785 * t827;
	t844 = qJD(3) * t773;
	t796 = -t766 * t755 + t763 * t760 + t764 * t806 + t767 * t844 + (-t753 * t794 - t758 * t828) * t784;
	t712 = t724 * qJD(5) + t747 * t789 + t793 * t796;
	t726 = t735 * t793 + t748 * t789;
	t788 = sin(qJ(6));
	t792 = cos(qJ(6));
	t754 = t783 * t844;
	t756 = t786 * t844;
	t757 = t806 * t783;
	t759 = t806 * t786;
	t797 = -t754 * t834 - t766 * t756 + t757 * t818 - t763 * t759 + t764 * t773 - t767 * t850;
	t809 = -t757 * t834 - t766 * t759 + t767 * t773;
	t852 = -t712 * t792 + t797 * t788 + (-t726 * t788 + t792 * t809) * qJD(6);
	t851 = (t726 * t792 + t788 * t809) * qJD(6) - t712 * t788 - t797 * t792;
	t711 = t726 * qJD(5) + t747 * t793 - t789 * t796;
	t835 = t784 * t791;
	t826 = qJD(5) * t789;
	t825 = qJD(5) * t793;
	t824 = qJD(6) * t788;
	t823 = qJD(6) * t792;
	t822 = qJD(6) * t793;
	t769 = -t821 + t829;
	t737 = t757 * t835 - t759 * t810 + t769 * t773;
	t761 = t766 * qJD(1);
	t762 = t767 * qJD(1);
	t795 = -t810 * t755 - t761 * t760 - t762 * t806 + t769 * t844 + (t753 * t791 - t758 * t827) * t784;
	t814 = -t737 * t822 + t795;
	t813 = -t809 * t822 + t796;
	t743 = t787 * t757 + (t759 * t785 + t773 * t782) * t784;
	t800 = t787 * t753 + (t755 * t785 + t782 * t844) * t784;
	t812 = -t743 * t822 + t800;
	t750 = t783 * t810 + t786 * t835;
	t805 = -t758 * t835 + t760 * t810 + t769 * t806;
	t728 = t750 * t789 + t793 * t805;
	t727 = t750 * t793 - t789 * t805;
	t765 = -t784 * t785 * t783 + t787 * t786;
	t799 = -t787 * t758 + (-t760 * t785 + t782 * t806) * t784;
	t730 = t765 * t789 + t793 * t799;
	t729 = t765 * t793 - t789 * t799;
	t807 = qJD(1) * t820 - t761 * t783;
	t715 = -t810 * t756 + t761 * t759 - t762 * t773 - t769 * t850 + (t754 * t791 + t757 * t827) * t784;
	t804 = -qJD(6) * t805 - t715 * t793 + t737 * t826;
	t803 = qJD(6) * t735 - t793 * t797 + t809 * t826;
	t739 = t787 * t754 + (t756 * t785 - t782 * t850) * t784;
	t802 = -qJD(6) * t799 - t739 * t793 + t743 * t826;
	t723 = t729 * qJD(5) + t793 * t800;
	t722 = -t730 * qJD(5) - t789 * t800;
	t710 = t727 * qJD(5) + t807 * t789 + t793 * t795;
	t709 = t728 * qJD(5) + t789 * t795 - t807 * t793;
	t708 = t710 * t792 - t715 * t788 + (-t728 * t788 - t737 * t792) * qJD(6);
	t707 = -t710 * t788 - t715 * t792 + (-t728 * t792 + t737 * t788) * qJD(6);
	t1 = [t852, 0, t814 * t788 - t804 * t792, 0, -t709 * t792 - t727 * t824, t707; t708, 0, t813 * t788 - t803 * t792, 0, t711 * t792 - t724 * t824, t851; 0, 0, t812 * t788 - t802 * t792, 0, t722 * t792 - t729 * t824, -t723 * t788 - t739 * t792 + (-t730 * t792 + t743 * t788) * qJD(6); -t851, 0, t804 * t788 + t814 * t792, 0, t709 * t788 - t727 * t823, -t708; t707, 0, t803 * t788 + t813 * t792, 0, -t711 * t788 - t724 * t823, t852; 0, 0, t802 * t788 + t812 * t792, 0, -t722 * t788 - t729 * t823, -t723 * t792 + t739 * t788 + (t730 * t788 + t743 * t792) * qJD(6); t711, 0, t715 * t789 + t737 * t825, 0, t710, 0; t709, 0, t789 * t797 + t809 * t825, 0, t712, 0; 0, 0, t739 * t789 + t743 * t825, 0, t723, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end