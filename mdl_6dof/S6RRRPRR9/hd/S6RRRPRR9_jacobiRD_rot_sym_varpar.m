% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:08
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPRR9_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR9_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
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
	% StartTime: 2019-10-10 12:08:08
	% EndTime: 2019-10-10 12:08:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (27->13), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t158 = sin(qJ(2));
	t159 = sin(qJ(1));
	t171 = t158 * t159;
	t161 = cos(qJ(1));
	t170 = t158 * t161;
	t160 = cos(qJ(2));
	t169 = t159 * t160;
	t168 = t160 * t161;
	t156 = sin(pkin(6));
	t167 = qJD(1) * t156;
	t166 = qJD(2) * t156;
	t157 = cos(pkin(6));
	t165 = -t157 * t168 + t171;
	t164 = t157 * t169 + t170;
	t163 = t157 * t170 + t169;
	t162 = t157 * t171 - t168;
	t155 = t162 * qJD(1) + t165 * qJD(2);
	t154 = t164 * qJD(1) + t163 * qJD(2);
	t153 = t163 * qJD(1) + t164 * qJD(2);
	t152 = t165 * qJD(1) + t162 * qJD(2);
	t1 = [t155, t152, 0, 0, 0, 0; -t153, -t154, 0, 0, 0, 0; 0, -t158 * t166, 0, 0, 0, 0; t154, t153, 0, 0, 0, 0; t152, t155, 0, 0, 0, 0; 0, -t160 * t166, 0, 0, 0, 0; -t159 * t167, 0, 0, 0, 0, 0; t161 * t167, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:09
	% EndTime: 2019-10-10 12:08:10
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (178->54), mult. (587->106), div. (0->0), fcn. (611->10), ass. (0->49)
	t364 = cos(pkin(6));
	t366 = sin(qJ(2));
	t367 = sin(qJ(1));
	t392 = t367 * t366;
	t385 = t364 * t392;
	t369 = cos(qJ(2));
	t370 = cos(qJ(1));
	t388 = t370 * t369;
	t351 = -qJD(1) * t385 - qJD(2) * t392 + (qJD(2) * t364 + qJD(1)) * t388;
	t389 = t370 * t366;
	t391 = t367 * t369;
	t353 = t364 * t389 + t391;
	t361 = sin(pkin(7));
	t365 = sin(qJ(3));
	t368 = cos(qJ(3));
	t375 = t364 * t391 + t389;
	t350 = t375 * qJD(1) + t353 * qJD(2);
	t363 = cos(pkin(7));
	t362 = sin(pkin(6));
	t387 = qJD(1) * t362;
	t384 = t367 * t387;
	t372 = -t350 * t363 + t361 * t384;
	t398 = t362 * t370;
	t352 = -t364 * t388 + t392;
	t401 = t352 * t363;
	t405 = ((t361 * t398 + t401) * t365 - t353 * t368) * qJD(3) - t351 * t365 + t372 * t368;
	t399 = t361 * t362;
	t397 = t363 * t365;
	t396 = t363 * t368;
	t395 = t365 * t366;
	t394 = t365 * t369;
	t393 = t366 * t368;
	t390 = t368 * t369;
	t386 = qJD(3) * t361;
	t383 = t370 * t387;
	t382 = t368 * t386;
	t380 = -t363 * t375 + t367 * t399;
	t379 = -t363 * t390 + t395;
	t378 = -t363 * t393 - t394;
	t377 = -t363 * t394 - t393;
	t376 = t363 * t395 - t390;
	t374 = t385 - t388;
	t348 = t352 * qJD(1) + t374 * qJD(2);
	t373 = t348 * t363 + t361 * t383;
	t371 = t382 * t398 + (qJD(3) * t401 - t351) * t368 + (qJD(3) * t353 - t372) * t365;
	t349 = t353 * qJD(1) + t375 * qJD(2);
	t347 = -t349 * t368 + t373 * t365 + (t365 * t374 + t380 * t368) * qJD(3);
	t346 = t349 * t365 + t373 * t368 + (-t380 * t365 + t368 * t374) * qJD(3);
	t1 = [t371, t349 * t397 + t348 * t368 + (t365 * t375 + t374 * t396) * qJD(3), t346, 0, 0, 0; t347, -t351 * t397 - t350 * t368 + (t352 * t365 - t353 * t396) * qJD(3), t405, 0, 0, 0; 0, (t377 * qJD(2) + t378 * qJD(3)) * t362, -t364 * t365 * t386 + (t378 * qJD(2) + t377 * qJD(3)) * t362, 0, 0, 0; -t405, t349 * t396 - t348 * t365 + (t368 * t375 - t374 * t397) * qJD(3), -t347, 0, 0, 0; t346, -t351 * t396 + t350 * t365 + (t352 * t368 + t353 * t397) * qJD(3), t371, 0, 0, 0; 0, (t379 * qJD(2) + t376 * qJD(3)) * t362, -t364 * t382 + (t376 * qJD(2) + t379 * qJD(3)) * t362, 0, 0, 0; -t350 * t361 - t363 * t384, -t349 * t361, 0, 0, 0, 0; -t348 * t361 + t363 * t383, t351 * t361, 0, 0, 0, 0; 0, qJD(2) * t369 * t399, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:10
	% EndTime: 2019-10-10 12:08:10
	% DurationCPUTime: 0.48s
	% Computational Cost: add. (260->63), mult. (843->114), div. (0->0), fcn. (907->12), ass. (0->45)
	t424 = sin(pkin(13));
	t429 = sin(qJ(3));
	t458 = cos(pkin(13));
	t459 = cos(qJ(3));
	t438 = -t429 * t424 + t459 * t458;
	t464 = t438 * qJD(3);
	t439 = t424 * t459 + t429 * t458;
	t463 = qJD(3) * t439;
	t428 = cos(pkin(6));
	t430 = sin(qJ(2));
	t433 = cos(qJ(1));
	t454 = t433 * t430;
	t431 = sin(qJ(1));
	t432 = cos(qJ(2));
	t455 = t431 * t432;
	t410 = t428 * t454 + t455;
	t443 = t428 * t455 + t454;
	t397 = qJD(1) * t443 + qJD(2) * t410;
	t456 = t431 * t430;
	t450 = t428 * t456;
	t453 = t433 * t432;
	t398 = -qJD(1) * t450 - qJD(2) * t456 + (qJD(2) * t428 + qJD(1)) * t453;
	t425 = sin(pkin(7));
	t400 = t425 * t463;
	t427 = cos(pkin(7));
	t402 = t427 * t463;
	t403 = t438 * t425;
	t405 = t438 * t427;
	t409 = -t428 * t453 + t456;
	t426 = sin(pkin(6));
	t452 = qJD(1) * t431;
	t461 = -t397 * t405 - t398 * t439 + t409 * t402 - t410 * t464 + (t400 * t433 + t403 * t452) * t426;
	t442 = t450 - t453;
	t395 = qJD(1) * t409 + qJD(2) * t442;
	t396 = qJD(1) * t410 + qJD(2) * t443;
	t399 = t464 * t425;
	t401 = t464 * t427;
	t404 = t439 * t425;
	t406 = t439 * t427;
	t451 = qJD(1) * t433;
	t460 = t395 * t406 - t396 * t438 - t443 * t401 + t442 * t463 + (t399 * t431 + t404 * t451) * t426;
	t448 = qJD(1) * t426 * t427;
	t434 = t397 * t406 - t398 * t438 + t409 * t401 + t410 * t463 + (t399 * t433 - t404 * t452) * t426;
	t394 = t395 * t405 + t396 * t439 + t443 * t402 + t442 * t464 + (-t400 * t431 + t403 * t451) * t426;
	t1 = [t434, t395 * t438 + t396 * t406 + t401 * t442 + t443 * t463, t394, 0, 0, 0; t460, -t397 * t438 - t398 * t406 - t401 * t410 + t409 * t463, t461, 0, 0, 0; 0, (-t401 * t430 - t463 * t432 + (-t406 * t432 - t430 * t438) * qJD(2)) * t426, -t428 * t400 + (-t402 * t432 - t464 * t430 + (-t405 * t430 - t432 * t439) * qJD(2)) * t426, 0, 0, 0; -t461, -t395 * t439 + t396 * t405 - t402 * t442 + t443 * t464, -t460, 0, 0, 0; t394, t397 * t439 - t398 * t405 + t402 * t410 + t409 * t464, t434, 0, 0, 0; 0, (t402 * t430 - t464 * t432 + (-t405 * t432 + t430 * t439) * qJD(2)) * t426, -t428 * t399 + (-t401 * t432 + t463 * t430 + (t406 * t430 - t432 * t438) * qJD(2)) * t426, 0, 0, 0; -t397 * t425 - t431 * t448, -t396 * t425, 0, 0, 0, 0; -t395 * t425 + t433 * t448, t398 * t425, 0, 0, 0, 0; 0, t426 * qJD(2) * t432 * t425, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:14
	% EndTime: 2019-10-10 12:08:16
	% DurationCPUTime: 1.17s
	% Computational Cost: add. (710->120), mult. (2227->220), div. (0->0), fcn. (2499->14), ass. (0->85)
	t659 = sin(pkin(13));
	t664 = sin(qJ(3));
	t706 = cos(pkin(13));
	t708 = cos(qJ(3));
	t672 = -t664 * t659 + t708 * t706;
	t716 = t672 * qJD(3);
	t665 = sin(qJ(2));
	t666 = sin(qJ(1));
	t668 = cos(qJ(2));
	t669 = cos(qJ(1));
	t707 = cos(pkin(6));
	t685 = t669 * t707;
	t645 = t665 * t685 + t666 * t668;
	t686 = t666 * t707;
	t674 = t669 * t665 + t668 * t686;
	t627 = t674 * qJD(1) + t645 * qJD(2);
	t682 = t665 * t686;
	t698 = t669 * t668;
	t699 = t666 * t665;
	t628 = -qJD(1) * t682 - qJD(2) * t699 + (qJD(2) * t707 + qJD(1)) * t698;
	t660 = sin(pkin(7));
	t633 = t716 * t660;
	t662 = cos(pkin(7));
	t635 = t716 * t662;
	t673 = t708 * t659 + t664 * t706;
	t638 = t673 * t660;
	t640 = t673 * t662;
	t644 = -t668 * t685 + t699;
	t661 = sin(pkin(6));
	t697 = qJD(1) * t666;
	t712 = qJD(3) * t673;
	t605 = -t627 * t640 + t628 * t672 - t644 * t635 - t645 * t712 + (-t633 * t669 + t638 * t697) * t661;
	t701 = t661 * t669;
	t614 = t638 * t701 + t644 * t640 - t645 * t672;
	t690 = t661 * t697;
	t623 = t627 * t660 + t662 * t690;
	t692 = t662 * t701;
	t629 = -t644 * t660 + t692;
	t663 = sin(qJ(5));
	t667 = cos(qJ(5));
	t715 = -t605 * t667 - t623 * t663 + (-t614 * t663 + t629 * t667) * qJD(5);
	t714 = -t605 * t663 + t623 * t667 + (t614 * t667 + t629 * t663) * qJD(5);
	t705 = t660 * t661;
	t704 = t660 * t663;
	t703 = t660 * t667;
	t702 = t661 * t666;
	t696 = qJD(1) * t669;
	t695 = qJD(5) * t663;
	t694 = qJD(5) * t667;
	t693 = t665 * t705;
	t689 = qJD(2) * t705;
	t684 = t665 * t689;
	t683 = t668 * t689;
	t639 = t672 * t662;
	t680 = t639 * t668 - t665 * t673;
	t679 = t640 * t668 + t665 * t672;
	t678 = -t640 * t665 + t668 * t672;
	t675 = t682 - t698;
	t634 = t660 * t712;
	t636 = t662 * t712;
	t637 = t672 * t660;
	t604 = -t627 * t639 - t628 * t673 + t634 * t701 + t644 * t636 + t637 * t690 - t645 * t716;
	t610 = t707 * t633 + (t678 * qJD(2) + t635 * t668 - t665 * t712) * t661;
	t625 = t644 * qJD(1) + t675 * qJD(2);
	t626 = t645 * qJD(1) + t674 * qJD(2);
	t603 = t625 * t640 - t626 * t672 - t674 * t635 + t675 * t712 + (t633 * t666 + t638 * t696) * t661;
	t641 = t707 * t662 - t668 * t705;
	t631 = t660 * t674 + t662 * t702;
	t624 = t678 * t661;
	t621 = qJD(1) * t692 - t625 * t660;
	t620 = t640 * t675 - t672 * t674;
	t619 = -t645 * t640 - t644 * t672;
	t618 = t707 * t638 + t679 * t661;
	t617 = t707 * t637 + t680 * t661;
	t616 = t638 * t702 - t640 * t674 - t672 * t675;
	t615 = t637 * t702 - t639 * t674 + t673 * t675;
	t612 = -t637 * t701 - t644 * t639 - t645 * t673;
	t611 = (-t679 * qJD(2) - t635 * t665 - t668 * t712) * t661;
	t609 = -t707 * t634 + (-t636 * t668 - t716 * t665 + (-t639 * t665 - t668 * t673) * qJD(2)) * t661;
	t608 = -t627 * t672 - t628 * t640 - t645 * t635 + t644 * t712;
	t607 = t625 * t672 + t626 * t640 + t635 * t675 + t674 * t712;
	t602 = t625 * t639 + t626 * t673 + t674 * t636 + t675 * t716 + (-t634 * t666 + t637 * t696) * t661;
	t601 = t603 * t667 + t621 * t663 + (-t616 * t663 + t631 * t667) * qJD(5);
	t600 = -t603 * t663 + t621 * t667 + (-t616 * t667 - t631 * t663) * qJD(5);
	t1 = [t715, -t626 * t704 + t607 * t667 + (-t620 * t663 - t675 * t703) * qJD(5), t602 * t667 - t615 * t695, 0, t600, 0; t601, t628 * t704 + t608 * t667 + (-t619 * t663 + t645 * t703) * qJD(5), t604 * t667 - t612 * t695, 0, t714, 0; 0, t663 * t683 + t611 * t667 + (-t624 * t663 + t667 * t693) * qJD(5), t609 * t667 - t617 * t695, 0, t667 * t684 - t610 * t663 + (-t618 * t667 - t641 * t663) * qJD(5), 0; -t714, -t626 * t703 - t607 * t663 + (-t620 * t667 + t675 * t704) * qJD(5), -t602 * t663 - t615 * t694, 0, -t601, 0; t600, t628 * t703 - t608 * t663 + (-t619 * t667 - t645 * t704) * qJD(5), -t604 * t663 - t612 * t694, 0, t715, 0; 0, t667 * t683 - t611 * t663 + (-t624 * t667 - t663 * t693) * qJD(5), -t609 * t663 - t617 * t694, 0, -t663 * t684 - t610 * t667 + (t618 * t663 - t641 * t667) * qJD(5), 0; t604, t625 * t673 - t626 * t639 + t636 * t675 - t674 * t716, t603, 0, 0, 0; -t602, -t627 * t673 + t628 * t639 - t645 * t636 - t644 * t716, t605, 0, 0, 0; 0, (t680 * qJD(2) - t636 * t665 + t668 * t716) * t661, t610, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:22
	% EndTime: 2019-10-10 12:08:25
	% DurationCPUTime: 2.44s
	% Computational Cost: add. (1820->176), mult. (5568->320), div. (0->0), fcn. (6440->16), ass. (0->121)
	t899 = sin(pkin(13));
	t905 = sin(qJ(3));
	t965 = cos(pkin(13));
	t967 = cos(qJ(3));
	t889 = -t967 * t899 - t905 * t965;
	t900 = sin(pkin(7));
	t878 = t889 * t900;
	t902 = cos(pkin(7));
	t880 = t889 * t902;
	t906 = sin(qJ(2));
	t907 = sin(qJ(1));
	t910 = cos(qJ(2));
	t911 = cos(qJ(1));
	t966 = cos(pkin(6));
	t940 = t911 * t966;
	t884 = t907 * t906 - t910 * t940;
	t885 = t906 * t940 + t907 * t910;
	t923 = -t905 * t899 + t967 * t965;
	t901 = sin(pkin(6));
	t960 = t901 * t911;
	t846 = -t878 * t960 - t884 * t880 - t885 * t923;
	t947 = t902 * t960;
	t867 = -t884 * t900 + t947;
	t904 = sin(qJ(5));
	t909 = cos(qJ(5));
	t828 = t846 * t904 - t867 * t909;
	t941 = t907 * t966;
	t924 = t911 * t906 + t910 * t941;
	t865 = t924 * qJD(1) + t885 * qJD(2);
	t957 = qJD(1) * t907;
	t945 = t901 * t957;
	t860 = t865 * t900 + t902 * t945;
	t937 = t906 * t941;
	t955 = qJD(2) * t906;
	t958 = t911 * t910;
	t866 = -qJD(1) * t937 - t907 * t955 + (qJD(2) * t966 + qJD(1)) * t958;
	t980 = t923 * qJD(3);
	t873 = t980 * t900;
	t875 = t980 * t902;
	t974 = qJD(3) * t889;
	t912 = t865 * t880 + t866 * t923 - t884 * t875 + t885 * t974 + (-t873 * t911 - t878 * t957) * t901;
	t809 = t828 * qJD(5) + t860 * t904 + t909 * t912;
	t830 = t846 * t909 + t867 * t904;
	t903 = sin(qJ(6));
	t908 = cos(qJ(6));
	t874 = t900 * t974;
	t876 = t902 * t974;
	t877 = t923 * t900;
	t879 = t923 * t902;
	t915 = -t865 * t879 + t866 * t889 - t874 * t960 - t884 * t876 + t877 * t945 - t885 * t980;
	t928 = -t877 * t960 - t884 * t879 + t885 * t889;
	t982 = -t809 * t908 + t915 * t903 + (-t830 * t903 + t908 * t928) * qJD(6);
	t981 = (t830 * t908 + t903 * t928) * qJD(6) - t809 * t903 - t915 * t908;
	t808 = t830 * qJD(5) + t860 * t909 - t904 * t912;
	t964 = t900 * t901;
	t963 = t900 * t904;
	t962 = t900 * t909;
	t961 = t901 * t907;
	t956 = qJD(1) * t911;
	t954 = qJD(5) * t904;
	t953 = qJD(5) * t909;
	t952 = qJD(6) * t903;
	t951 = qJD(6) * t908;
	t950 = qJD(6) * t909;
	t949 = t906 * t964;
	t948 = t910 * t964;
	t944 = t901 * t955;
	t939 = t900 * t944;
	t938 = qJD(2) * t948;
	t925 = t937 - t958;
	t848 = t877 * t961 - t879 * t924 - t889 * t925;
	t863 = t884 * qJD(1) + t925 * qJD(2);
	t864 = t885 * qJD(1) + t924 * qJD(2);
	t913 = -t863 * t880 - t864 * t923 - t924 * t875 - t925 * t974 + (t873 * t907 - t878 * t956) * t901;
	t935 = -t848 * t950 + t913;
	t934 = -t928 * t950 + t912;
	t932 = t879 * t910 + t889 * t906;
	t851 = t966 * t877 + t932 * t901;
	t930 = t880 * t906 + t910 * t923;
	t914 = t966 * t873 + (t930 * qJD(2) + t875 * t910 + t906 * t974) * t901;
	t933 = -t851 * t950 + t914;
	t869 = t900 * t924 + t902 * t961;
	t922 = -t878 * t961 + t880 * t924 - t923 * t925;
	t832 = t869 * t904 + t909 * t922;
	t831 = t869 * t909 - t904 * t922;
	t881 = t966 * t902 - t948;
	t931 = -t880 * t910 + t906 * t923;
	t916 = -t966 * t878 + t931 * t901;
	t839 = t881 * t904 + t909 * t916;
	t838 = t881 * t909 - t904 * t916;
	t854 = t885 * t880 - t884 * t923;
	t840 = t854 * t909 + t885 * t963;
	t856 = -t880 * t925 - t923 * t924;
	t841 = t856 * t909 - t925 * t963;
	t862 = t930 * t901;
	t857 = t862 * t909 + t904 * t949;
	t926 = qJD(1) * t947 - t863 * t900;
	t816 = t863 * t879 - t864 * t889 - t924 * t876 + t925 * t980 + (t874 * t907 + t877 * t956) * t901;
	t921 = -qJD(6) * t922 - t816 * t909 + t848 * t954;
	t920 = qJD(6) * t846 - t909 * t915 + t928 * t954;
	t834 = t966 * t874 - t879 * t944 + (-t980 * t906 + (qJD(2) * t889 + t876) * t910) * t901;
	t919 = -qJD(6) * t916 - t834 * t909 + t851 * t954;
	t861 = (t879 * t906 - t889 * t910) * t901;
	t855 = -t879 * t925 + t889 * t924;
	t853 = t885 * t879 + t884 * t889;
	t837 = (-t931 * qJD(2) - t875 * t906 + t910 * t974) * t901;
	t836 = (t932 * qJD(2) + t876 * t906 + t910 * t980) * t901;
	t827 = -t865 * t923 + t866 * t880 - t885 * t875 - t884 * t974;
	t826 = t865 * t889 + t866 * t879 + t885 * t876 - t884 * t980;
	t825 = t863 * t923 - t864 * t880 + t875 * t925 - t924 * t974;
	t824 = -t863 * t889 - t864 * t879 - t876 * t925 - t924 * t980;
	t823 = t904 * t938 + t837 * t909 + (-t862 * t904 + t909 * t949) * qJD(5);
	t814 = t838 * qJD(5) + t904 * t939 + t909 * t914;
	t813 = -t839 * qJD(5) - t904 * t914 + t909 * t939;
	t812 = t866 * t963 + t827 * t909 + (-t854 * t904 + t885 * t962) * qJD(5);
	t811 = -t864 * t963 + t825 * t909 + (-t856 * t904 - t925 * t962) * qJD(5);
	t807 = t831 * qJD(5) + t926 * t904 + t909 * t913;
	t806 = t832 * qJD(5) + t904 * t913 - t926 * t909;
	t805 = t807 * t908 - t816 * t903 + (-t832 * t903 - t848 * t908) * qJD(6);
	t804 = -t807 * t903 - t816 * t908 + (-t832 * t908 + t848 * t903) * qJD(6);
	t1 = [t982, t811 * t908 + t824 * t903 + (-t841 * t903 + t855 * t908) * qJD(6), t935 * t903 - t921 * t908, 0, -t806 * t908 - t831 * t952, t804; t805, t812 * t908 + t826 * t903 + (-t840 * t903 + t853 * t908) * qJD(6), t934 * t903 - t920 * t908, 0, t808 * t908 - t828 * t952, t981; 0, t823 * t908 + t836 * t903 + (-t857 * t903 + t861 * t908) * qJD(6), t933 * t903 - t919 * t908, 0, t813 * t908 - t838 * t952, -t814 * t903 - t834 * t908 + (-t839 * t908 + t851 * t903) * qJD(6); -t981, -t811 * t903 + t824 * t908 + (-t841 * t908 - t855 * t903) * qJD(6), t921 * t903 + t935 * t908, 0, t806 * t903 - t831 * t951, -t805; t804, -t812 * t903 + t826 * t908 + (-t840 * t908 - t853 * t903) * qJD(6), t920 * t903 + t934 * t908, 0, -t808 * t903 - t828 * t951, t982; 0, -t823 * t903 + t836 * t908 + (-t857 * t908 - t861 * t903) * qJD(6), t919 * t903 + t933 * t908, 0, -t813 * t903 - t838 * t951, -t814 * t908 + t834 * t903 + (t839 * t903 + t851 * t908) * qJD(6); t808, t841 * qJD(5) + t825 * t904 + t864 * t962, t816 * t904 + t848 * t953, 0, t807, 0; t806, t840 * qJD(5) + t827 * t904 - t866 * t962, t904 * t915 + t928 * t953, 0, t809, 0; 0, t857 * qJD(5) + t837 * t904 - t909 * t938, t834 * t904 + t851 * t953, 0, t814, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end