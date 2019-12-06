% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRRRR10_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR10_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:27:37
	% EndTime: 2019-12-05 17:27:37
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:27:37
	% EndTime: 2019-12-05 17:27:37
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:27:37
	% EndTime: 2019-12-05 17:27:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(5));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(5));
	t60 = cos(pkin(11));
	t58 = sin(pkin(11));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0; 0, -t62 * t64, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0; 0, -t63 * t64, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:27:38
	% EndTime: 2019-12-05 17:27:38
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (78->41), mult. (285->88), div. (0->0), fcn. (301->10), ass. (0->35)
	t259 = sin(pkin(11));
	t262 = cos(pkin(11));
	t266 = sin(qJ(2));
	t264 = cos(pkin(5));
	t268 = cos(qJ(2));
	t283 = t264 * t268;
	t253 = -t259 * t266 + t262 * t283;
	t260 = sin(pkin(6));
	t261 = sin(pkin(5));
	t287 = t260 * t261;
	t263 = cos(pkin(6));
	t265 = sin(qJ(3));
	t286 = t263 * t265;
	t267 = cos(qJ(3));
	t285 = t263 * t267;
	t284 = t264 * t266;
	t282 = t265 * t266;
	t281 = t265 * t268;
	t280 = t266 * t267;
	t279 = t267 * t268;
	t277 = qJD(3) * t260 * t264;
	t276 = -t253 * t263 + t262 * t287;
	t274 = t259 * t283 + t262 * t266;
	t275 = -t259 * t287 + t263 * t274;
	t254 = t259 * t268 + t262 * t284;
	t273 = t259 * t284 - t262 * t268;
	t272 = -t263 * t279 + t282;
	t271 = -t263 * t280 - t281;
	t270 = -t263 * t281 - t280;
	t269 = t263 * t282 - t279;
	t252 = t273 * qJD(2);
	t251 = t274 * qJD(2);
	t250 = t254 * qJD(2);
	t249 = t253 * qJD(2);
	t1 = [0, t251 * t286 + t252 * t267 + (t265 * t274 + t273 * t285) * qJD(3), t252 * t285 + t251 * t265 + (t275 * t265 + t267 * t273) * qJD(3), 0, 0; 0, -t249 * t286 - t250 * t267 + (-t253 * t265 - t254 * t285) * qJD(3), -t250 * t285 - t249 * t265 + (-t254 * t267 + t276 * t265) * qJD(3), 0, 0; 0, (t270 * qJD(2) + t271 * qJD(3)) * t261, -t265 * t277 + (t271 * qJD(2) + t270 * qJD(3)) * t261, 0, 0; 0, t251 * t285 - t252 * t265 + (t267 * t274 - t273 * t286) * qJD(3), -t252 * t286 + t251 * t267 + (-t265 * t273 + t275 * t267) * qJD(3), 0, 0; 0, -t249 * t285 + t250 * t265 + (-t253 * t267 + t254 * t286) * qJD(3), t250 * t286 - t249 * t267 + (t254 * t265 + t276 * t267) * qJD(3), 0, 0; 0, (t272 * qJD(2) + t269 * qJD(3)) * t261, -t267 * t277 + (t269 * qJD(2) + t272 * qJD(3)) * t261, 0, 0; 0, -t251 * t260, 0, 0, 0; 0, t249 * t260, 0, 0, 0; 0, qJD(2) * t268 * t287, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:27:40
	% EndTime: 2019-12-05 17:27:40
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (269->87), mult. (921->177), div. (0->0), fcn. (1021->12), ass. (0->68)
	t472 = sin(pkin(11));
	t475 = cos(pkin(11));
	t480 = sin(qJ(2));
	t477 = cos(pkin(5));
	t483 = cos(qJ(2));
	t504 = t477 * t483;
	t466 = -t472 * t480 + t475 * t504;
	t479 = sin(qJ(3));
	t482 = cos(qJ(3));
	t505 = t477 * t480;
	t488 = t472 * t505 - t475 * t483;
	t476 = cos(pkin(6));
	t489 = t472 * t504 + t475 * t480;
	t473 = sin(pkin(6));
	t474 = sin(pkin(5));
	t512 = t473 * t474;
	t490 = t472 * t512 - t476 * t489;
	t452 = t490 * t479 - t482 * t488;
	t467 = t472 * t483 + t475 * t505;
	t491 = -t466 * t476 + t475 * t512;
	t516 = -t467 * t482 + t491 * t479;
	t511 = t473 * t477;
	t478 = sin(qJ(4));
	t510 = t473 * t478;
	t481 = cos(qJ(4));
	t509 = t473 * t481;
	t508 = t474 * t476;
	t507 = t476 * t479;
	t506 = t476 * t482;
	t503 = t479 * t480;
	t502 = t479 * t483;
	t501 = t480 * t482;
	t500 = t482 * t483;
	t499 = qJD(4) * t478;
	t498 = qJD(4) * t481;
	t497 = t480 * t512;
	t495 = qJD(2) * t512;
	t494 = qJD(3) * t511;
	t493 = t480 * t495;
	t492 = t483 * t495;
	t454 = t466 * t482 - t467 * t507;
	t455 = -t482 * t489 + t488 * t507;
	t487 = t476 * t500 - t503;
	t486 = -t476 * t501 - t502;
	t485 = t476 * t502 + t501;
	t484 = -t476 * t503 + t500;
	t449 = -t467 * t479 - t491 * t482;
	t451 = t479 * t488 + t490 * t482;
	t465 = t477 * t476 - t483 * t512;
	t464 = t488 * qJD(2);
	t463 = t489 * qJD(2);
	t462 = t467 * qJD(2);
	t461 = t466 * qJD(2);
	t460 = t484 * t474;
	t459 = t472 * t508 + t473 * t489;
	t458 = -t466 * t473 - t475 * t508;
	t457 = t485 * t474 + t479 * t511;
	t456 = t487 * t474 + t482 * t511;
	t453 = (-t485 * qJD(2) + t486 * qJD(3)) * t474;
	t448 = t482 * t494 + (t484 * qJD(2) + t487 * qJD(3)) * t474;
	t447 = -t479 * t494 + (t486 * qJD(2) - t485 * qJD(3)) * t474;
	t446 = t463 * t507 + t464 * t482 + (t479 * t489 + t488 * t506) * qJD(3);
	t445 = -t461 * t507 - t462 * t482 + (-t466 * t479 - t467 * t506) * qJD(3);
	t444 = t451 * qJD(3) - t463 * t482 + t464 * t507;
	t443 = -t452 * qJD(3) + t463 * t479 + t464 * t506;
	t442 = t449 * qJD(3) + t461 * t482 - t462 * t507;
	t441 = t516 * qJD(3) - t461 * t479 - t462 * t506;
	t1 = [0, -t463 * t510 + t446 * t481 + (-t455 * t478 - t488 * t509) * qJD(4), t443 * t481 - t451 * t499, -t464 * t509 - t444 * t478 + (-t452 * t481 - t459 * t478) * qJD(4), 0; 0, t461 * t510 + t445 * t481 + (-t454 * t478 + t467 * t509) * qJD(4), t441 * t481 - t449 * t499, t462 * t509 - t442 * t478 + (-t458 * t478 + t481 * t516) * qJD(4), 0; 0, t478 * t492 + t453 * t481 + (-t460 * t478 + t481 * t497) * qJD(4), t447 * t481 - t456 * t499, t481 * t493 - t448 * t478 + (-t457 * t481 - t465 * t478) * qJD(4), 0; 0, -t463 * t509 - t446 * t478 + (-t455 * t481 + t488 * t510) * qJD(4), -t443 * t478 - t451 * t498, t464 * t510 - t444 * t481 + (t452 * t478 - t459 * t481) * qJD(4), 0; 0, t461 * t509 - t445 * t478 + (-t454 * t481 - t467 * t510) * qJD(4), -t441 * t478 - t449 * t498, -t462 * t510 - t442 * t481 + (-t458 * t481 - t478 * t516) * qJD(4), 0; 0, t481 * t492 - t453 * t478 + (-t460 * t481 - t478 * t497) * qJD(4), -t447 * t478 - t456 * t498, -t478 * t493 - t448 * t481 + (t457 * t478 - t465 * t481) * qJD(4), 0; 0, t455 * qJD(3) - t463 * t506 + t464 * t479, t444, 0, 0; 0, t454 * qJD(3) + t461 * t506 - t462 * t479, t442, 0, 0; 0, (t487 * qJD(2) + t484 * qJD(3)) * t474, t448, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:27:42
	% EndTime: 2019-12-05 17:27:42
	% DurationCPUTime: 0.74s
	% Computational Cost: add. (784->147), mult. (2555->280), div. (0->0), fcn. (2925->14), ass. (0->110)
	t677 = cos(pkin(6));
	t685 = cos(qJ(3));
	t686 = cos(qJ(2));
	t715 = t685 * t686;
	t681 = sin(qJ(3));
	t682 = sin(qJ(2));
	t718 = t681 * t682;
	t693 = t677 * t715 - t718;
	t673 = sin(pkin(11));
	t676 = cos(pkin(11));
	t678 = cos(pkin(5));
	t719 = t678 * t686;
	t663 = -t673 * t682 + t676 * t719;
	t720 = t678 * t682;
	t664 = t673 * t686 + t676 * t720;
	t733 = t664 * t681;
	t732 = t664 * t685;
	t694 = t673 * t720 - t676 * t686;
	t731 = t694 * t681;
	t674 = sin(pkin(6));
	t675 = sin(pkin(5));
	t729 = t674 * t675;
	t680 = sin(qJ(4));
	t728 = t674 * t680;
	t727 = t674 * t681;
	t684 = cos(qJ(4));
	t726 = t674 * t684;
	t725 = t674 * t685;
	t724 = t675 * t676;
	t723 = t675 * t677;
	t722 = t677 * t681;
	t721 = t677 * t685;
	t717 = t681 * t686;
	t716 = t682 * t685;
	t714 = qJD(4) * t680;
	t713 = qJD(4) * t684;
	t679 = sin(qJ(5));
	t712 = qJD(5) * t679;
	t683 = cos(qJ(5));
	t711 = qJD(5) * t683;
	t710 = qJD(5) * t684;
	t709 = t682 * t729;
	t708 = t675 * t725;
	t705 = t678 * t725;
	t704 = qJD(2) * t729;
	t703 = qJD(3) * t727;
	t702 = t682 * t704;
	t701 = t686 * t704;
	t658 = t663 * qJD(2);
	t659 = t664 * qJD(2);
	t697 = t663 * t677 - t674 * t724;
	t619 = -t659 * t722 + t658 * t685 + (t697 * t685 - t733) * qJD(3);
	t635 = -t663 * t721 + t676 * t708 + t733;
	t700 = t635 * t710 + t619;
	t695 = t673 * t719 + t676 * t682;
	t660 = t695 * qJD(2);
	t661 = t694 * qJD(2);
	t696 = t673 * t729 - t677 * t695;
	t621 = t661 * t722 - t660 * t685 + (t696 * t685 + t731) * qJD(3);
	t637 = -t673 * t708 + t695 * t721 - t731;
	t699 = t637 * t710 + t621;
	t690 = -t677 * t718 + t715;
	t634 = qJD(3) * t705 + (t690 * qJD(2) + t693 * qJD(3)) * t675;
	t648 = -t693 * t675 - t705;
	t698 = t648 * t710 + t634;
	t636 = t697 * t681 + t732;
	t650 = -t663 * t674 - t676 * t723;
	t628 = t636 * t684 + t650 * t680;
	t627 = -t636 * t680 + t650 * t684;
	t638 = t696 * t681 - t685 * t694;
	t651 = t673 * t723 + t674 * t695;
	t630 = t638 * t684 + t651 * t680;
	t629 = -t638 * t680 + t651 * t684;
	t691 = t677 * t717 + t716;
	t649 = t691 * t675 + t678 * t727;
	t662 = t678 * t677 - t686 * t729;
	t640 = t649 * t684 + t662 * t680;
	t639 = -t649 * t680 + t662 * t684;
	t644 = t663 * t685 - t664 * t722;
	t631 = t644 * t684 + t664 * t728;
	t646 = -t685 * t695 + t694 * t722;
	t632 = t646 * t684 - t694 * t728;
	t643 = t663 * t681 + t664 * t721;
	t645 = -t681 * t695 - t694 * t721;
	t692 = t677 * t716 + t717;
	t657 = t690 * t675;
	t647 = t657 * t684 + t680 * t709;
	t618 = t658 * t681 + t659 * t721 - t703 * t724 + (t663 * t722 + t732) * qJD(3);
	t689 = qJD(5) * t636 - t618 * t684 + t635 * t714;
	t620 = t638 * qJD(3) - t660 * t681 - t661 * t721;
	t688 = qJD(5) * t638 - t620 * t684 + t637 * t714;
	t633 = t678 * t703 + (t692 * qJD(2) + t691 * qJD(3)) * t675;
	t687 = qJD(5) * t649 - t633 * t684 + t648 * t714;
	t656 = t692 * t675;
	t642 = (-t691 * qJD(2) - t692 * qJD(3)) * t675;
	t641 = (t693 * qJD(2) + t690 * qJD(3)) * t675;
	t626 = -t645 * qJD(3) + t660 * t722 + t661 * t685;
	t625 = t646 * qJD(3) - t660 * t721 + t661 * t681;
	t624 = -t643 * qJD(3) - t658 * t722 - t659 * t685;
	t623 = t644 * qJD(3) + t658 * t721 - t659 * t681;
	t622 = t680 * t701 + t642 * t684 + (-t657 * t680 + t684 * t709) * qJD(4);
	t617 = t639 * qJD(4) + t634 * t684 + t680 * t702;
	t616 = -t640 * qJD(4) - t634 * t680 + t684 * t702;
	t615 = -t660 * t728 + t626 * t684 + (-t646 * t680 - t694 * t726) * qJD(4);
	t614 = t658 * t728 + t624 * t684 + (-t644 * t680 + t664 * t726) * qJD(4);
	t613 = t629 * qJD(4) + t621 * t684 - t661 * t728;
	t612 = -t630 * qJD(4) - t621 * t680 - t661 * t726;
	t611 = t627 * qJD(4) + t619 * t684 + t659 * t728;
	t610 = -t628 * qJD(4) - t619 * t680 + t659 * t726;
	t1 = [0, t615 * t683 + t625 * t679 + (-t632 * t679 + t645 * t683) * qJD(5), t699 * t679 + t688 * t683, t612 * t683 - t629 * t712, -t613 * t679 + t620 * t683 + (-t630 * t683 - t637 * t679) * qJD(5); 0, t614 * t683 + t623 * t679 + (-t631 * t679 + t643 * t683) * qJD(5), t700 * t679 + t689 * t683, t610 * t683 - t627 * t712, -t611 * t679 + t618 * t683 + (-t628 * t683 - t635 * t679) * qJD(5); 0, t622 * t683 + t641 * t679 + (-t647 * t679 + t656 * t683) * qJD(5), t698 * t679 + t687 * t683, t616 * t683 - t639 * t712, -t617 * t679 + t633 * t683 + (-t640 * t683 - t648 * t679) * qJD(5); 0, -t615 * t679 + t625 * t683 + (-t632 * t683 - t645 * t679) * qJD(5), -t688 * t679 + t699 * t683, -t612 * t679 - t629 * t711, -t613 * t683 - t620 * t679 + (t630 * t679 - t637 * t683) * qJD(5); 0, -t614 * t679 + t623 * t683 + (-t631 * t683 - t643 * t679) * qJD(5), -t689 * t679 + t700 * t683, -t610 * t679 - t627 * t711, -t611 * t683 - t618 * t679 + (t628 * t679 - t635 * t683) * qJD(5); 0, -t622 * t679 + t641 * t683 + (-t647 * t683 - t656 * t679) * qJD(5), -t687 * t679 + t698 * t683, -t616 * t679 - t639 * t711, -t617 * t683 - t633 * t679 + (t640 * t679 - t648 * t683) * qJD(5); 0, t632 * qJD(4) + t626 * t680 + t660 * t726, -t620 * t680 - t637 * t713, t613, 0; 0, t631 * qJD(4) + t624 * t680 - t658 * t726, -t618 * t680 - t635 * t713, t611, 0; 0, t647 * qJD(4) + t642 * t680 - t684 * t701, -t633 * t680 - t648 * t713, t617, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end