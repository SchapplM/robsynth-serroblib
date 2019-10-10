% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:22
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PPRRRR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobiRD_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (20->14), mult. (80->39), div. (0->0), fcn. (88->10), ass. (0->20)
	t132 = sin(pkin(13));
	t138 = cos(pkin(6));
	t147 = t132 * t138;
	t133 = sin(pkin(7));
	t134 = sin(pkin(6));
	t146 = t133 * t134;
	t145 = t133 * t138;
	t135 = cos(pkin(14));
	t137 = cos(pkin(7));
	t144 = t135 * t137;
	t136 = cos(pkin(13));
	t143 = t136 * t138;
	t131 = sin(pkin(14));
	t142 = -(-t132 * t131 + t135 * t143) * t137 + t136 * t146;
	t141 = -(-t136 * t131 - t135 * t147) * t137 - t132 * t146;
	t140 = cos(qJ(3));
	t139 = sin(qJ(3));
	t130 = -t131 * t147 + t136 * t135;
	t128 = t131 * t143 + t132 * t135;
	t1 = [0, 0, (-t130 * t140 + t141 * t139) * qJD(3), 0, 0, 0; 0, 0, (-t128 * t140 + t142 * t139) * qJD(3), 0, 0, 0; 0, 0, (-t139 * t145 + (-t131 * t140 - t139 * t144) * t134) * qJD(3), 0, 0, 0; 0, 0, (t130 * t139 + t141 * t140) * qJD(3), 0, 0, 0; 0, 0, (t128 * t139 + t142 * t140) * qJD(3), 0, 0, 0; 0, 0, (-t140 * t145 + (t131 * t139 - t140 * t144) * t134) * qJD(3), 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:23
	% EndTime: 2019-10-09 21:22:23
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (225->67), mult. (755->133), div. (0->0), fcn. (881->14), ass. (0->52)
	t382 = sin(pkin(14));
	t386 = sin(pkin(6));
	t387 = cos(pkin(14));
	t393 = sin(qJ(3));
	t395 = cos(qJ(3));
	t390 = cos(pkin(7));
	t408 = t390 * t393;
	t385 = sin(pkin(7));
	t391 = cos(pkin(6));
	t414 = t385 * t391;
	t368 = (t382 * t395 + t387 * t408) * t386 + t393 * t414;
	t388 = cos(pkin(13));
	t383 = sin(pkin(13));
	t416 = t383 * t391;
	t377 = -t382 * t416 + t388 * t387;
	t376 = -t388 * t382 - t387 * t416;
	t415 = t385 * t386;
	t397 = t376 * t390 + t383 * t415;
	t364 = t377 * t395 + t397 * t393;
	t411 = t388 * t391;
	t375 = t382 * t411 + t383 * t387;
	t419 = t375 * t393;
	t418 = t375 * t395;
	t413 = t386 * t387;
	t412 = t386 * t390;
	t389 = cos(pkin(8));
	t392 = sin(qJ(4));
	t410 = t389 * t392;
	t394 = cos(qJ(4));
	t409 = t389 * t394;
	t407 = qJD(3) * t393;
	t406 = qJD(3) * t395;
	t405 = t388 * t415;
	t403 = t385 * t406;
	t402 = t390 * t406;
	t374 = -t383 * t382 + t387 * t411;
	t398 = -t374 * t390 + t405;
	t361 = -t398 * t395 - t419;
	t384 = sin(pkin(8));
	t401 = -t361 * t389 - (-t374 * t385 - t388 * t412) * t384;
	t363 = -t377 * t393 + t397 * t395;
	t400 = -t363 * t389 - (-t376 * t385 + t383 * t412) * t384;
	t367 = t395 * t414 + (t387 * t390 * t395 - t382 * t393) * t386;
	t399 = -t367 * t389 - (-t385 * t413 + t391 * t390) * t384;
	t366 = t368 * qJD(3);
	t365 = t386 * t382 * t407 - t391 * t403 - t402 * t413;
	t362 = t374 * t408 - t393 * t405 + t418;
	t360 = t364 * qJD(3);
	t359 = -t383 * t386 * t403 - t376 * t402 + t377 * t407;
	t358 = (t398 * t393 - t418) * qJD(3);
	t357 = -t374 * t402 + (t395 * t405 + t419) * qJD(3);
	t1 = [0, 0, t359 * t410 - t360 * t394 + (-t363 * t392 - t364 * t409) * qJD(4), -t360 * t409 + t359 * t392 + (-t364 * t394 + t400 * t392) * qJD(4), 0, 0; 0, 0, t357 * t410 + t358 * t394 + (-t361 * t392 - t362 * t409) * qJD(4), t358 * t409 + t357 * t392 + (-t362 * t394 + t401 * t392) * qJD(4), 0, 0; 0, 0, t365 * t410 - t366 * t394 + (-t367 * t392 - t368 * t409) * qJD(4), -t366 * t409 + t365 * t392 + (-t368 * t394 + t399 * t392) * qJD(4), 0, 0; 0, 0, t359 * t409 + t360 * t392 + (-t363 * t394 + t364 * t410) * qJD(4), t360 * t410 + t359 * t394 + (t364 * t392 + t400 * t394) * qJD(4), 0, 0; 0, 0, t357 * t409 - t358 * t392 + (-t361 * t394 + t362 * t410) * qJD(4), -t358 * t410 + t357 * t394 + (t362 * t392 + t401 * t394) * qJD(4), 0, 0; 0, 0, t365 * t409 + t366 * t392 + (-t367 * t394 + t368 * t410) * qJD(4), t366 * t410 + t365 * t394 + (t368 * t392 + t399 * t394) * qJD(4), 0, 0; 0, 0, -t359 * t384, 0, 0, 0; 0, 0, -t357 * t384, 0, 0, 0; 0, 0, -t365 * t384, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:26
	% EndTime: 2019-10-09 21:22:27
	% DurationCPUTime: 0.80s
	% Computational Cost: add. (740->109), mult. (2392->212), div. (0->0), fcn. (2883->16), ass. (0->82)
	t609 = sin(pkin(14));
	t613 = sin(pkin(6));
	t614 = cos(pkin(14));
	t621 = sin(qJ(3));
	t624 = cos(qJ(3));
	t617 = cos(pkin(7));
	t639 = t617 * t621;
	t612 = sin(pkin(7));
	t618 = cos(pkin(6));
	t645 = t612 * t618;
	t595 = (t609 * t624 + t614 * t639) * t613 + t621 * t645;
	t620 = sin(qJ(4));
	t623 = cos(qJ(4));
	t594 = t624 * t645 + (t614 * t617 * t624 - t609 * t621) * t613;
	t644 = t613 * t614;
	t600 = -t612 * t644 + t617 * t618;
	t611 = sin(pkin(8));
	t616 = cos(pkin(8));
	t628 = t594 * t616 + t600 * t611;
	t579 = t595 * t623 + t628 * t620;
	t615 = cos(pkin(13));
	t610 = sin(pkin(13));
	t649 = t610 * t618;
	t604 = -t609 * t649 + t614 * t615;
	t603 = -t609 * t615 - t614 * t649;
	t646 = t612 * t613;
	t626 = t603 * t617 + t610 * t646;
	t590 = t604 * t624 + t626 * t621;
	t589 = -t604 * t621 + t626 * t624;
	t643 = t613 * t617;
	t597 = -t603 * t612 + t610 * t643;
	t629 = t589 * t616 + t597 * t611;
	t574 = t590 * t623 + t629 * t620;
	t642 = t615 * t618;
	t601 = -t609 * t610 + t614 * t642;
	t634 = t615 * t646;
	t602 = t609 * t642 + t610 * t614;
	t651 = t602 * t624;
	t588 = t601 * t639 - t621 * t634 + t651;
	t627 = -t601 * t617 + t634;
	t652 = t602 * t621;
	t587 = -t627 * t624 - t652;
	t596 = -t601 * t612 - t615 * t643;
	t630 = t587 * t616 + t596 * t611;
	t572 = t588 * t623 + t630 * t620;
	t619 = sin(qJ(5));
	t648 = t611 * t619;
	t622 = cos(qJ(5));
	t647 = t611 * t622;
	t641 = t616 * t620;
	t640 = t616 * t623;
	t638 = qJD(3) * t621;
	t637 = qJD(3) * t624;
	t636 = qJD(5) * t619;
	t635 = qJD(5) * t622;
	t632 = t612 * t637;
	t631 = t617 * t637;
	t576 = t587 * t623 - t588 * t641;
	t577 = t589 * t623 - t590 * t641;
	t582 = t594 * t623 - t595 * t641;
	t571 = -t588 * t620 + t630 * t623;
	t573 = -t590 * t620 + t629 * t623;
	t578 = -t595 * t620 + t628 * t623;
	t593 = t595 * qJD(3);
	t592 = t609 * t613 * t638 - t618 * t632 - t631 * t644;
	t591 = -t594 * t611 + t600 * t616;
	t586 = t590 * qJD(3);
	t585 = -t610 * t613 * t632 - t603 * t631 + t604 * t638;
	t584 = (t627 * t621 - t651) * qJD(3);
	t583 = -t601 * t631 + (t624 * t634 + t652) * qJD(3);
	t581 = -t589 * t611 + t597 * t616;
	t580 = -t587 * t611 + t596 * t616;
	t575 = t592 * t641 - t593 * t623 + (-t594 * t620 - t595 * t640) * qJD(4);
	t570 = t578 * qJD(4) - t592 * t623 - t593 * t641;
	t569 = -t579 * qJD(4) + t592 * t620 - t593 * t640;
	t568 = t585 * t641 - t586 * t623 + (-t589 * t620 - t590 * t640) * qJD(4);
	t567 = t583 * t641 + t584 * t623 + (-t587 * t620 - t588 * t640) * qJD(4);
	t566 = t573 * qJD(4) - t585 * t623 - t586 * t641;
	t565 = -t574 * qJD(4) + t585 * t620 - t586 * t640;
	t564 = t571 * qJD(4) - t583 * t623 + t584 * t641;
	t563 = -t572 * qJD(4) + t583 * t620 + t584 * t640;
	t1 = [0, 0, -t585 * t648 + t568 * t622 + (-t577 * t619 + t590 * t647) * qJD(5), t565 * t622 - t573 * t636, t586 * t647 - t566 * t619 + (-t574 * t622 - t581 * t619) * qJD(5), 0; 0, 0, -t583 * t648 + t567 * t622 + (-t576 * t619 + t588 * t647) * qJD(5), t563 * t622 - t571 * t636, -t584 * t647 - t564 * t619 + (-t572 * t622 - t580 * t619) * qJD(5), 0; 0, 0, -t592 * t648 + t575 * t622 + (-t582 * t619 + t595 * t647) * qJD(5), t569 * t622 - t578 * t636, t593 * t647 - t570 * t619 + (-t579 * t622 - t591 * t619) * qJD(5), 0; 0, 0, -t585 * t647 - t568 * t619 + (-t577 * t622 - t590 * t648) * qJD(5), -t565 * t619 - t573 * t635, -t586 * t648 - t566 * t622 + (t574 * t619 - t581 * t622) * qJD(5), 0; 0, 0, -t583 * t647 - t567 * t619 + (-t576 * t622 - t588 * t648) * qJD(5), -t563 * t619 - t571 * t635, t584 * t648 - t564 * t622 + (t572 * t619 - t580 * t622) * qJD(5), 0; 0, 0, -t592 * t647 - t575 * t619 + (-t582 * t622 - t595 * t648) * qJD(5), -t569 * t619 - t578 * t635, -t593 * t648 - t570 * t622 + (t579 * t619 - t591 * t622) * qJD(5), 0; 0, 0, t577 * qJD(4) - t585 * t640 - t586 * t620, t566, 0, 0; 0, 0, t576 * qJD(4) - t583 * t640 + t584 * t620, t564, 0, 0; 0, 0, t582 * qJD(4) - t592 * t640 - t593 * t620, t570, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:33
	% EndTime: 2019-10-09 21:22:34
	% DurationCPUTime: 1.50s
	% Computational Cost: add. (2099->169), mult. (6590->311), div. (0->0), fcn. (8107->18), ass. (0->123)
	t808 = sin(pkin(14));
	t812 = sin(pkin(6));
	t813 = cos(pkin(14));
	t821 = sin(qJ(3));
	t825 = cos(qJ(3));
	t816 = cos(pkin(7));
	t855 = t816 * t821;
	t811 = sin(pkin(7));
	t817 = cos(pkin(6));
	t861 = t811 * t817;
	t796 = (t808 * t825 + t813 * t855) * t812 + t821 * t861;
	t814 = cos(pkin(13));
	t809 = sin(pkin(13));
	t865 = t809 * t817;
	t803 = -t808 * t865 + t813 * t814;
	t802 = -t808 * t814 - t813 * t865;
	t862 = t811 * t812;
	t837 = t802 * t816 + t809 * t862;
	t787 = t803 * t825 + t837 * t821;
	t858 = t814 * t817;
	t800 = -t808 * t809 + t813 * t858;
	t847 = t814 * t862;
	t801 = t808 * t858 + t809 * t813;
	t867 = t801 * t825;
	t785 = t800 * t855 - t821 * t847 + t867;
	t820 = sin(qJ(4));
	t871 = t785 * t820;
	t870 = t787 * t820;
	t869 = t796 * t820;
	t868 = t801 * t821;
	t810 = sin(pkin(8));
	t819 = sin(qJ(5));
	t864 = t810 * t819;
	t823 = cos(qJ(5));
	t863 = t810 * t823;
	t860 = t812 * t813;
	t859 = t812 * t816;
	t815 = cos(pkin(8));
	t857 = t815 * t820;
	t824 = cos(qJ(4));
	t856 = t815 * t824;
	t854 = qJD(3) * t821;
	t853 = qJD(3) * t825;
	t852 = qJD(5) * t819;
	t851 = qJD(5) * t823;
	t818 = sin(qJ(6));
	t850 = qJD(6) * t818;
	t822 = cos(qJ(6));
	t849 = qJD(6) * t822;
	t848 = qJD(6) * t823;
	t845 = t811 * t853;
	t844 = t816 * t853;
	t780 = -t800 * t844 + (t825 * t847 + t868) * qJD(3);
	t839 = -t800 * t816 + t847;
	t781 = (t839 * t821 - t867) * qJD(3);
	t784 = -t839 * t825 - t868;
	t840 = -t800 * t811 - t814 * t859;
	t834 = t840 * t810;
	t828 = t784 * t815 + t834;
	t742 = t781 * t857 - t780 * t824 + (t828 * t824 - t871) * qJD(4);
	t757 = -t784 * t856 - t824 * t834 + t871;
	t843 = t757 * t848 + t742;
	t782 = -t809 * t812 * t845 - t802 * t844 + t803 * t854;
	t783 = t787 * qJD(3);
	t786 = -t803 * t821 + t837 * t825;
	t838 = -t802 * t811 + t809 * t859;
	t833 = t838 * t810;
	t827 = t786 * t815 + t833;
	t744 = -t783 * t857 - t782 * t824 + (t827 * t824 - t870) * qJD(4);
	t759 = -t786 * t856 - t824 * t833 + t870;
	t842 = t759 * t848 + t744;
	t793 = t808 * t812 * t854 - t817 * t845 - t844 * t860;
	t794 = t796 * qJD(3);
	t795 = t825 * t861 + (t813 * t816 * t825 - t808 * t821) * t812;
	t835 = -t811 * t860 + t816 * t817;
	t832 = t835 * t810;
	t826 = t795 * t815 + t832;
	t756 = -t794 * t857 - t793 * t824 + (t826 * t824 - t869) * qJD(4);
	t770 = -t795 * t856 - t824 * t832 + t869;
	t841 = t770 * t848 + t756;
	t758 = t785 * t824 + t828 * t820;
	t772 = -t784 * t810 + t840 * t815;
	t746 = t758 * t823 + t772 * t819;
	t745 = -t758 * t819 + t772 * t823;
	t760 = t787 * t824 + t827 * t820;
	t773 = -t786 * t810 + t838 * t815;
	t748 = t760 * t823 + t773 * t819;
	t747 = -t760 * t819 + t773 * t823;
	t771 = t796 * t824 + t826 * t820;
	t788 = -t795 * t810 + t835 * t815;
	t762 = t771 * t823 + t788 * t819;
	t761 = -t771 * t819 + t788 * t823;
	t766 = t784 * t824 - t785 * t857;
	t753 = t766 * t823 + t785 * t864;
	t768 = t786 * t824 - t787 * t857;
	t754 = t768 * t823 + t787 * t864;
	t775 = t795 * t824 - t796 * t857;
	t769 = t775 * t823 + t796 * t864;
	t765 = t784 * t820 + t785 * t856;
	t767 = t786 * t820 + t787 * t856;
	t774 = t795 * t820 + t796 * t856;
	t741 = t758 * qJD(4) - t780 * t820 - t781 * t856;
	t831 = qJD(6) * t758 - t741 * t823 + t757 * t852;
	t743 = t760 * qJD(4) - t782 * t820 + t783 * t856;
	t830 = qJD(6) * t760 - t743 * t823 + t759 * t852;
	t755 = t771 * qJD(4) - t793 * t820 + t794 * t856;
	t829 = qJD(6) * t771 - t755 * t823 + t770 * t852;
	t764 = -t774 * qJD(4) + t793 * t857 - t794 * t824;
	t763 = t775 * qJD(4) - t793 * t856 - t794 * t820;
	t752 = -t767 * qJD(4) + t782 * t857 - t783 * t824;
	t751 = t768 * qJD(4) - t782 * t856 - t783 * t820;
	t750 = -t765 * qJD(4) + t780 * t857 + t781 * t824;
	t749 = t766 * qJD(4) - t780 * t856 + t781 * t820;
	t740 = -t793 * t864 + t764 * t823 + (-t775 * t819 + t796 * t863) * qJD(5);
	t739 = t761 * qJD(5) + t756 * t823 + t794 * t864;
	t738 = -t762 * qJD(5) - t756 * t819 + t794 * t863;
	t737 = -t782 * t864 + t752 * t823 + (-t768 * t819 + t787 * t863) * qJD(5);
	t736 = -t780 * t864 + t750 * t823 + (-t766 * t819 + t785 * t863) * qJD(5);
	t735 = t747 * qJD(5) + t744 * t823 + t783 * t864;
	t734 = -t748 * qJD(5) - t744 * t819 + t783 * t863;
	t733 = t745 * qJD(5) + t742 * t823 - t781 * t864;
	t732 = -t746 * qJD(5) - t742 * t819 - t781 * t863;
	t1 = [0, 0, t737 * t822 + t751 * t818 + (-t754 * t818 + t767 * t822) * qJD(6), t842 * t818 + t830 * t822, t734 * t822 - t747 * t850, -t735 * t818 + t743 * t822 + (-t748 * t822 - t759 * t818) * qJD(6); 0, 0, t736 * t822 + t749 * t818 + (-t753 * t818 + t765 * t822) * qJD(6), t843 * t818 + t831 * t822, t732 * t822 - t745 * t850, -t733 * t818 + t741 * t822 + (-t746 * t822 - t757 * t818) * qJD(6); 0, 0, t740 * t822 + t763 * t818 + (-t769 * t818 + t774 * t822) * qJD(6), t841 * t818 + t829 * t822, t738 * t822 - t761 * t850, -t739 * t818 + t755 * t822 + (-t762 * t822 - t770 * t818) * qJD(6); 0, 0, -t737 * t818 + t751 * t822 + (-t754 * t822 - t767 * t818) * qJD(6), -t830 * t818 + t842 * t822, -t734 * t818 - t747 * t849, -t735 * t822 - t743 * t818 + (t748 * t818 - t759 * t822) * qJD(6); 0, 0, -t736 * t818 + t749 * t822 + (-t753 * t822 - t765 * t818) * qJD(6), -t831 * t818 + t843 * t822, -t732 * t818 - t745 * t849, -t733 * t822 - t741 * t818 + (t746 * t818 - t757 * t822) * qJD(6); 0, 0, -t740 * t818 + t763 * t822 + (-t769 * t822 - t774 * t818) * qJD(6), -t829 * t818 + t841 * t822, -t738 * t818 - t761 * t849, -t739 * t822 - t755 * t818 + (t762 * t818 - t770 * t822) * qJD(6); 0, 0, t754 * qJD(5) + t752 * t819 + t782 * t863, -t743 * t819 - t759 * t851, t735, 0; 0, 0, t753 * qJD(5) + t750 * t819 + t780 * t863, -t741 * t819 - t757 * t851, t733, 0; 0, 0, t769 * qJD(5) + t764 * t819 + t793 * t863, -t755 * t819 - t770 * t851, t739, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end