% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRR12
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:16
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR12_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR12_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR12_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobiaD_transl_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (14->13), mult. (48->27), div. (0->0), fcn. (38->6), ass. (0->12)
	t114 = sin(pkin(6));
	t124 = t114 * (r_i_i_C(3) + qJ(2));
	t116 = cos(pkin(6));
	t117 = sin(qJ(1));
	t122 = t116 * t117;
	t118 = cos(qJ(1));
	t121 = t116 * t118;
	t120 = qJD(1) * t114;
	t119 = t114 * qJD(2);
	t115 = cos(pkin(14));
	t113 = sin(pkin(14));
	t1 = [t118 * t119 + ((t113 * t122 - t115 * t118) * r_i_i_C(1) + (t113 * t118 + t115 * t122) * r_i_i_C(2) - t118 * pkin(1) - t117 * t124) * qJD(1), t118 * t120, 0, 0, 0, 0; t117 * t119 + ((-t113 * t121 - t115 * t117) * r_i_i_C(1) + (t113 * t117 - t115 * t121) * r_i_i_C(2) - t117 * pkin(1) + t118 * t124) * qJD(1), t117 * t120, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:47
	% EndTime: 2019-10-10 09:15:47
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (103->49), mult. (358->92), div. (0->0), fcn. (352->10), ass. (0->47)
	t248 = cos(pkin(7));
	t278 = pkin(10) + r_i_i_C(3);
	t279 = t278 * t248 + qJ(2);
	t250 = sin(qJ(3));
	t277 = t250 * r_i_i_C(2);
	t246 = sin(pkin(6));
	t251 = sin(qJ(1));
	t276 = t246 * t251;
	t253 = cos(qJ(1));
	t275 = t246 * t253;
	t274 = t248 * t250;
	t252 = cos(qJ(3));
	t273 = t248 * t252;
	t244 = sin(pkin(14));
	t272 = t251 * t244;
	t247 = cos(pkin(14));
	t271 = t251 * t247;
	t270 = t253 * t244;
	t269 = t253 * t247;
	t268 = qJD(1) * t251;
	t267 = qJD(1) * t253;
	t245 = sin(pkin(7));
	t266 = qJD(3) * t245;
	t265 = t278 * t245;
	t249 = cos(pkin(6));
	t264 = t249 * t272;
	t263 = t246 * t268;
	t262 = t246 * t267;
	t261 = t266 * t275;
	t260 = r_i_i_C(1) * t250 + r_i_i_C(2) * t252;
	t236 = -qJD(1) * t264 + t247 * t267;
	t237 = -t249 * t269 + t272;
	t259 = qJD(3) * t237 * t248 - t236;
	t258 = t260 * t245;
	t256 = t249 * t271 + t270;
	t257 = t245 * t276 - t248 * t256;
	t238 = t249 * t270 + t271;
	t233 = t237 * qJD(1);
	t255 = t233 * t248 + t245 * t262;
	t235 = t256 * qJD(1);
	t254 = -qJD(3) * t238 - t235 * t248 + t245 * t263;
	t241 = t252 * t261;
	t240 = -t264 + t269;
	t234 = t238 * qJD(1);
	t232 = -t234 * t252 + t255 * t250 + (-t240 * t250 + t257 * t252) * qJD(3);
	t231 = t234 * t250 + t255 * t252 + (-t240 * t252 - t257 * t250) * qJD(3);
	t1 = [-pkin(1) * t267 + t241 * r_i_i_C(1) + (-t252 * r_i_i_C(1) - pkin(2) + t277) * t236 + (t260 * t248 - t265) * t235 + ((t237 * t273 + t238 * t250) * r_i_i_C(1) + (-t237 * t274 + t238 * t252) * r_i_i_C(2)) * qJD(3) + ((-t266 * t277 + qJD(2)) * t253 + (-t258 - t279) * t268) * t246, t262, t231 * r_i_i_C(1) - t232 * r_i_i_C(2), 0, 0, 0; qJD(2) * t276 - t234 * pkin(2) + t232 * r_i_i_C(1) + t231 * r_i_i_C(2) - t233 * t265 + (-t251 * pkin(1) + t279 * t275) * qJD(1), t263, t241 * r_i_i_C(2) + (t254 * r_i_i_C(1) + t259 * r_i_i_C(2)) * t252 + ((t259 + t261) * r_i_i_C(1) - t254 * r_i_i_C(2)) * t250, 0, 0, 0; 0, 0, (-t249 * t258 + ((-t244 * t252 - t247 * t274) * r_i_i_C(1) + (t244 * t250 - t247 * t273) * r_i_i_C(2)) * t246) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:50
	% EndTime: 2019-10-10 09:15:51
	% DurationCPUTime: 1.16s
	% Computational Cost: add. (609->114), mult. (2041->203), div. (0->0), fcn. (2200->14), ass. (0->86)
	t488 = cos(pkin(14));
	t491 = cos(pkin(6));
	t484 = sin(pkin(14));
	t494 = sin(qJ(1));
	t530 = t494 * t484;
	t519 = t491 * t530;
	t497 = cos(qJ(1));
	t526 = qJD(1) * t497;
	t473 = -qJD(1) * t519 + t488 * t526;
	t493 = sin(qJ(3));
	t496 = cos(qJ(3));
	t528 = t497 * t484;
	t529 = t494 * t488;
	t504 = t491 * t529 + t528;
	t472 = t504 * qJD(1);
	t490 = cos(pkin(7));
	t486 = sin(pkin(7));
	t487 = sin(pkin(6));
	t536 = t487 * t494;
	t518 = qJD(1) * t536;
	t514 = t486 * t518;
	t501 = -t472 * t490 + t514;
	t535 = t487 * t497;
	t522 = t486 * t535;
	t527 = t497 * t488;
	t475 = -t491 * t527 + t530;
	t542 = t475 * t490;
	t507 = t522 + t542;
	t476 = t491 * t528 + t529;
	t540 = t476 * t496;
	t448 = -t473 * t493 + t501 * t496 + (t507 * t493 - t540) * qJD(3);
	t544 = t472 * t486;
	t462 = t490 * t518 + t544;
	t489 = cos(pkin(8));
	t495 = cos(qJ(4));
	t532 = t489 * t495;
	t485 = sin(pkin(8));
	t538 = t485 * t495;
	t556 = t448 * t532 + t462 * t538;
	t492 = sin(qJ(4));
	t533 = t489 * t492;
	t539 = t485 * t492;
	t555 = -t448 * t533 - t462 * t539;
	t454 = t476 * t493 + t507 * t496;
	t531 = t490 * t493;
	t455 = t475 * t531 + t493 * t522 - t540;
	t465 = -t475 * t486 + t490 * t535;
	t554 = t454 * t532 - t455 * t492 + t465 * t538;
	t553 = t454 * t533 + t455 * t495 + t465 * t539;
	t546 = pkin(11) + r_i_i_C(3);
	t537 = t486 * t491;
	t464 = (t484 * t496 + t488 * t531) * t487 + t493 * t537;
	t503 = t519 - t527;
	t506 = t486 * t536 - t490 * t504;
	t457 = t506 * t493 - t503 * t496;
	t498 = t546 * t485 - (t492 * r_i_i_C(1) + t495 * r_i_i_C(2)) * t489;
	t470 = t475 * qJD(1);
	t545 = t470 * t486;
	t534 = t488 * t490;
	t525 = qJD(3) * t493;
	t524 = qJD(3) * t496;
	t523 = t487 * qJD(2);
	t520 = t496 * t537;
	t517 = t487 * t526;
	t516 = t487 * t524;
	t515 = pkin(10) * t490 + qJ(2);
	t471 = t476 * qJD(1);
	t502 = t470 * t490 + t486 * t517;
	t446 = -t457 * qJD(3) + t471 * t493 + t502 * t496;
	t460 = t490 * t517 - t545;
	t512 = t446 * t489 + t460 * t485;
	t500 = t503 * t493;
	t456 = t506 * t496 + t500;
	t509 = t456 * t489 + (t486 * t504 + t490 * t536) * t485;
	t508 = t495 * r_i_i_C(1) - t492 * r_i_i_C(2) + pkin(3);
	t479 = t497 * t486 * t516;
	t474 = -t487 * t488 * t486 + t491 * t490;
	t463 = t520 + (-t484 * t493 + t496 * t534) * t487;
	t459 = t464 * qJD(3);
	t458 = t487 * t484 * t525 - qJD(3) * t520 - t516 * t534;
	t451 = t479 + (qJD(3) * t542 - t473) * t496 + (qJD(3) * t476 - t501) * t493;
	t449 = t493 * t514 + t473 * t496 - t476 * t525 - t479 + (-t472 * t493 - t475 * t524) * t490;
	t447 = qJD(3) * t500 + t502 * t493 + (t506 * qJD(3) - t471) * t496;
	t445 = t447 * t495 + t512 * t492 + (-t457 * t492 + t509 * t495) * qJD(4);
	t444 = -t447 * t492 + t512 * t495 + (-t457 * t495 - t509 * t492) * qJD(4);
	t1 = [(t451 * t495 + t555) * r_i_i_C(1) + (-t451 * t492 - t556) * r_i_i_C(2) + t451 * pkin(3) - t473 * pkin(2) - pkin(10) * t544 + t497 * t523 + (-t497 * pkin(1) - t515 * t536) * qJD(1) + (t554 * r_i_i_C(1) - t553 * r_i_i_C(2)) * qJD(4) + t546 * (t448 * t485 - t462 * t489), t517, t508 * t446 + t498 * t447 + ((-t456 * t492 - t457 * t532) * r_i_i_C(1) + (-t456 * t495 + t457 * t533) * r_i_i_C(2)) * qJD(4), t444 * r_i_i_C(1) - t445 * r_i_i_C(2), 0, 0; t445 * r_i_i_C(1) + t444 * r_i_i_C(2) + t447 * pkin(3) - t471 * pkin(2) - pkin(10) * t545 + t494 * t523 + (-t494 * pkin(1) + t515 * t535) * qJD(1) + t546 * (-t446 * t485 + t460 * t489), t518, t508 * t448 + t498 * t449 + ((t454 * t492 + t455 * t532) * r_i_i_C(1) + (t454 * t495 - t455 * t533) * r_i_i_C(2)) * qJD(4), (-t449 * t492 + t556) * r_i_i_C(1) + (-t449 * t495 + t555) * r_i_i_C(2) + (t553 * r_i_i_C(1) + t554 * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, -t508 * t459 - t498 * t458 + ((-t463 * t492 - t464 * t532) * r_i_i_C(1) + (-t463 * t495 + t464 * t533) * r_i_i_C(2)) * qJD(4), (t458 * t492 - t459 * t532) * r_i_i_C(1) + (t458 * t495 + t459 * t533) * r_i_i_C(2) + ((-t463 * t533 - t464 * t495 - t474 * t539) * r_i_i_C(1) + (-t463 * t532 + t464 * t492 - t474 * t538) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:55
	% EndTime: 2019-10-10 09:15:58
	% DurationCPUTime: 2.67s
	% Computational Cost: add. (1915->180), mult. (6243->309), div. (0->0), fcn. (7022->16), ass. (0->115)
	t698 = cos(pkin(6));
	t695 = cos(pkin(14));
	t706 = cos(qJ(1));
	t743 = t706 * t695;
	t691 = sin(pkin(14));
	t702 = sin(qJ(1));
	t746 = t702 * t691;
	t682 = -t698 * t743 + t746;
	t701 = sin(qJ(3));
	t693 = sin(pkin(7));
	t694 = sin(pkin(6));
	t751 = t694 * t706;
	t736 = t693 * t751;
	t697 = cos(pkin(7));
	t747 = t697 * t701;
	t744 = t706 * t691;
	t745 = t702 * t695;
	t683 = t698 * t744 + t745;
	t705 = cos(qJ(3));
	t754 = t683 * t705;
	t662 = t682 * t747 + t701 * t736 - t754;
	t700 = sin(qJ(4));
	t704 = cos(qJ(4));
	t756 = t682 * t697;
	t721 = t736 + t756;
	t661 = t683 * t701 + t721 * t705;
	t672 = -t682 * t693 + t697 * t751;
	t692 = sin(pkin(8));
	t696 = cos(pkin(8));
	t725 = t661 * t696 + t672 * t692;
	t640 = t662 * t704 + t725 * t700;
	t654 = t661 * t692 - t672 * t696;
	t699 = sin(qJ(5));
	t703 = cos(qJ(5));
	t781 = -t640 * t699 - t654 * t703;
	t780 = t640 * t703 - t654 * t699;
	t733 = t698 * t746;
	t742 = qJD(1) * t706;
	t680 = -qJD(1) * t733 + t695 * t742;
	t718 = t698 * t745 + t744;
	t679 = t718 * qJD(1);
	t752 = t694 * t702;
	t732 = qJD(1) * t752;
	t729 = t693 * t732;
	t713 = -t679 * t697 + t729;
	t646 = -t680 * t701 + t713 * t705 + (t721 * t701 - t754) * qJD(3);
	t758 = t679 * t693;
	t669 = t697 * t732 + t758;
	t728 = t646 * t696 + t669 * t692;
	t779 = t640 * qJD(4) + t728 * t704;
	t771 = -t662 * t700 + t725 * t704;
	t778 = -t771 * qJD(4) + t728 * t700;
	t635 = t646 * t692 - t669 * t696;
	t777 = t635 * t699;
	t776 = t635 * t703;
	t753 = t693 * t698;
	t671 = (t691 * t705 + t695 * t747) * t694 + t701 * t753;
	t734 = t705 * t753;
	t750 = t695 * t697;
	t670 = t734 + (-t691 * t701 + t705 * t750) * t694;
	t681 = -t694 * t695 * t693 + t698 * t697;
	t723 = t670 * t696 + t681 * t692;
	t653 = t671 * t704 + t723 * t700;
	t717 = t733 - t743;
	t720 = t693 * t752 - t697 * t718;
	t664 = t720 * t701 - t717 * t705;
	t762 = r_i_i_C(3) + pkin(12);
	t667 = t671 * qJD(3);
	t760 = t667 * t692;
	t749 = t696 * t700;
	t748 = t696 * t704;
	t741 = qJD(3) * t701;
	t740 = qJD(3) * t705;
	t739 = qJD(5) * t699;
	t738 = qJD(5) * t703;
	t737 = t694 * qJD(2);
	t731 = t694 * t742;
	t730 = t694 * t740;
	t712 = t717 * t701;
	t663 = t720 * t705 + t712;
	t674 = t693 * t718 + t697 * t752;
	t724 = t663 * t696 + t674 * t692;
	t722 = t703 * r_i_i_C(1) - t699 * r_i_i_C(2) + pkin(4);
	t650 = -t661 * t704 + t662 * t749;
	t651 = t663 * t704 - t664 * t749;
	t657 = t670 * t704 - t671 * t749;
	t716 = qJD(5) * (-t699 * r_i_i_C(1) - t703 * r_i_i_C(2));
	t677 = t682 * qJD(1);
	t715 = -t677 * t693 + t697 * t731;
	t714 = t677 * t697 + t693 * t731;
	t710 = t715 * t692;
	t708 = -t664 * t700 + t724 * t704;
	t642 = t664 * t704 + t724 * t700;
	t707 = -t671 * t700 + t723 * t704;
	t678 = t683 * qJD(1);
	t644 = -t664 * qJD(3) + t678 * t701 + t714 * t705;
	t633 = -t644 * t692 + t715 * t696;
	t686 = t706 * t693 * t730;
	t666 = t694 * t691 * t741 - qJD(3) * t734 - t730 * t750;
	t658 = -t670 * t692 + t681 * t696;
	t656 = -t663 * t692 + t674 * t696;
	t649 = t686 + (qJD(3) * t756 - t680) * t705 + (qJD(3) * t683 - t713) * t701;
	t647 = t701 * t729 + t680 * t705 - t683 * t741 - t686 + (-t679 * t701 - t682 * t740) * t697;
	t645 = qJD(3) * t712 + t714 * t701 + (t720 * qJD(3) - t678) * t705;
	t637 = t666 * t749 - t667 * t704 + (-t670 * t700 - t671 * t748) * qJD(4);
	t632 = t707 * qJD(4) - t666 * t704 - t667 * t749;
	t630 = -t647 * t749 + t646 * t704 + (t661 * t700 + t662 * t748) * qJD(4);
	t628 = -t645 * t749 + t644 * t704 + (-t663 * t700 - t664 * t748) * qJD(4);
	t626 = t649 * t704 - t778;
	t624 = t647 * t704 + t778;
	t622 = t645 * t704 + (t644 * t696 + t710) * t700 + t708 * qJD(4);
	t621 = t642 * qJD(4) - t644 * t748 + t645 * t700 - t704 * t710;
	t620 = t622 * t703 + t633 * t699 + (-t642 * t699 + t656 * t703) * qJD(5);
	t619 = -t622 * t699 + t633 * t703 + (-t642 * t703 - t656 * t699) * qJD(5);
	t1 = [(t626 * t703 + t777) * r_i_i_C(1) + (-t626 * t699 + t776) * r_i_i_C(2) + t626 * pkin(4) + t649 * pkin(3) - t680 * pkin(2) - pkin(10) * t758 + t706 * t737 + t762 * (t649 * t700 + t779) + (t781 * r_i_i_C(1) - t780 * r_i_i_C(2)) * qJD(5) + t635 * pkin(11) + (-t706 * pkin(1) + (-pkin(10) * t697 - qJ(2)) * t752) * qJD(1), t731, (t628 * t703 - t651 * t739) * r_i_i_C(1) + (-t628 * t699 - t651 * t738) * r_i_i_C(2) + t628 * pkin(4) + t644 * pkin(3) + t762 * (t651 * qJD(4) + t644 * t700 + t645 * t748) + ((t645 * t699 + t664 * t738) * r_i_i_C(1) + (t645 * t703 - t664 * t739) * r_i_i_C(2) + t645 * pkin(11)) * t692, -t722 * t621 + t762 * t622 + t708 * t716, t619 * r_i_i_C(1) - t620 * r_i_i_C(2), 0; t702 * t737 - t678 * pkin(2) + t645 * pkin(3) + t622 * pkin(4) + t620 * r_i_i_C(1) + t619 * r_i_i_C(2) + t762 * t621 + (-t702 * pkin(1) + qJ(2) * t751) * qJD(1) + t633 * pkin(11) + t715 * pkin(10), t732, (t630 * t703 - t650 * t739) * r_i_i_C(1) + (-t630 * t699 - t650 * t738) * r_i_i_C(2) + t630 * pkin(4) + t646 * pkin(3) + t762 * (t650 * qJD(4) + t646 * t700 + t647 * t748) + ((t647 * t699 - t662 * t738) * r_i_i_C(1) + (t647 * t703 + t662 * t739) * r_i_i_C(2) + t647 * pkin(11)) * t692, t762 * t624 - t771 * t716 + t722 * (-t647 * t700 + t779), (-t624 * t699 - t776) * r_i_i_C(1) + (-t624 * t703 + t777) * r_i_i_C(2) + (t780 * r_i_i_C(1) + t781 * r_i_i_C(2)) * qJD(5), 0; 0, 0, (t637 * t703 - t657 * t739) * r_i_i_C(1) + (-t637 * t699 - t657 * t738) * r_i_i_C(2) + t637 * pkin(4) - t667 * pkin(3) + t762 * (t657 * qJD(4) - t666 * t748 - t667 * t700) + ((-t666 * t699 + t671 * t738) * r_i_i_C(1) + (-t666 * t703 - t671 * t739) * r_i_i_C(2) - t666 * pkin(11)) * t692, t762 * t632 + t707 * t716 + t722 * (-t653 * qJD(4) + t666 * t700 - t667 * t748), (-t632 * t699 + t703 * t760) * r_i_i_C(1) + (-t632 * t703 - t699 * t760) * r_i_i_C(2) + ((-t653 * t703 - t658 * t699) * r_i_i_C(1) + (t653 * t699 - t658 * t703) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:16:02
	% EndTime: 2019-10-10 09:16:08
	% DurationCPUTime: 5.42s
	% Computational Cost: add. (5036->257), mult. (16095->423), div. (0->0), fcn. (18656->18), ass. (0->155)
	t876 = sin(qJ(3));
	t880 = cos(qJ(3));
	t877 = sin(qJ(1));
	t868 = sin(pkin(14));
	t970 = cos(pkin(6));
	t943 = t868 * t970;
	t969 = cos(pkin(14));
	t974 = cos(qJ(1));
	t910 = -t877 * t969 - t974 * t943;
	t870 = sin(pkin(6));
	t968 = sin(pkin(7));
	t941 = t968 * t870;
	t926 = t974 * t941;
	t931 = t970 * t969;
	t860 = t877 * t868 - t974 * t931;
	t872 = cos(pkin(7));
	t963 = t860 * t872;
	t914 = t926 + t963;
	t839 = -t876 * t910 + t914 * t880;
	t840 = t914 * t876 + t910 * t880;
	t871 = cos(pkin(8));
	t875 = sin(qJ(4));
	t869 = sin(pkin(8));
	t898 = t860 * t968;
	t948 = t870 * t974;
	t915 = -t872 * t948 + t898;
	t908 = t915 * t869;
	t973 = cos(qJ(4));
	t813 = t840 * t973 + (t839 * t871 - t908) * t875;
	t829 = t839 * t869 + t915 * t871;
	t874 = sin(qJ(5));
	t879 = cos(qJ(5));
	t796 = t813 * t879 - t829 * t874;
	t873 = sin(qJ(6));
	t1004 = t796 * t873;
	t878 = cos(qJ(6));
	t1003 = t796 * t878;
	t911 = t974 * t868 + t877 * t931;
	t858 = t911 * qJD(1);
	t944 = t858 * t968;
	t956 = qJD(1) * t877;
	t946 = t870 * t956;
	t916 = t872 * t946 + t944;
	t912 = -t877 * t943 + t974 * t969;
	t859 = t912 * qJD(1);
	t936 = t877 * t941;
	t924 = qJD(1) * t936;
	t913 = -t858 * t872 + t924;
	t985 = t840 * qJD(3) - t859 * t876 + t913 * t880;
	t805 = t869 * t985 - t916 * t871;
	t1002 = t796 * qJD(5) - t805 * t879;
	t929 = t813 * t874 + t829 * t879;
	t1001 = t929 * qJD(5) - t805 * t874;
	t998 = t813 * qJD(4);
	t909 = t916 * t869;
	t994 = (t871 * t985 + t909) * t875;
	t895 = -t911 * t872 + t936;
	t947 = t871 * t973;
	t991 = t839 * t947 - t840 * t875;
	t975 = r_i_i_C(3) + pkin(13);
	t983 = pkin(10) * t872 + qJ(2);
	t978 = t914 * qJD(1) - t912 * qJD(3);
	t904 = qJD(1) * t910;
	t979 = t895 * qJD(3) + t904;
	t885 = t979 * t876 - t978 * t880;
	t841 = -t912 * t876 + t895 * t880;
	t950 = t841 * t973;
	t982 = qJD(4) * t950 - t885 * t875;
	t892 = qJD(1) * t915;
	t881 = t885 * t869 - t871 * t892;
	t842 = t895 * t876 + t912 * t880;
	t958 = t870 * t877;
	t894 = t872 * t958 + t911 * t968;
	t893 = t894 * t869;
	t815 = t842 * t973 + (t841 * t871 + t893) * t875;
	t831 = -t841 * t869 + t894 * t871;
	t980 = -t815 * t874 + t831 * t879;
	t921 = qJD(6) * (t873 * r_i_i_C(1) + t878 * r_i_i_C(2));
	t930 = t970 * t968;
	t942 = t872 * t969;
	t853 = (t868 * t880 + t876 * t942) * t870 + t876 * t930;
	t971 = pkin(11) * t869;
	t960 = t869 * t874;
	t959 = t869 * t879;
	t957 = t871 * t875;
	t955 = qJD(3) * t876;
	t954 = qJD(3) * t880;
	t953 = qJD(4) * t875;
	t952 = qJD(6) * t873;
	t951 = qJD(6) * t878;
	t949 = t869 * t973;
	t945 = t842 * t953;
	t940 = t842 * t947;
	t938 = qJD(1) * t948;
	t934 = t880 * t942;
	t798 = t815 * t879 + t831 * t874;
	t922 = t880 * t930;
	t852 = t922 + (-t868 * t876 + t934) * t870;
	t907 = t970 * t872 - t969 * t941;
	t902 = t907 * t869;
	t828 = t853 * t973 + (t852 * t871 + t902) * t875;
	t836 = -t852 * t869 + t907 * t871;
	t807 = t828 * t879 + t836 * t874;
	t928 = -t828 * t874 + t836 * t879;
	t927 = t878 * r_i_i_C(1) - t873 * r_i_i_C(2) + pkin(5);
	t823 = -t839 * t973 + t840 * t957;
	t799 = t823 * t879 - t840 * t960;
	t825 = -t842 * t957 + t950;
	t800 = t825 * t879 + t842 * t960;
	t833 = t852 * t973 - t853 * t957;
	t826 = t833 * t879 + t853 * t960;
	t919 = t839 * t875 + t840 * t947;
	t918 = -t852 * t875 - t853 * t947;
	t905 = -t975 * t874 - t927 * t879 - pkin(4);
	t901 = -t915 * t949 + t991;
	t810 = -t973 * t908 + t991;
	t891 = t973 * t893;
	t890 = t869 * t892;
	t888 = t879 * t921 + (t927 * t874 - t975 * t879) * qJD(5);
	t827 = -t852 * t947 + t853 * t875 - t973 * t902;
	t886 = t985 * t973;
	t882 = t885 * t973;
	t863 = t926 * t954;
	t850 = t853 * qJD(3);
	t849 = t868 * t870 * t955 + (-t870 * t934 - t922) * qJD(3);
	t824 = t841 * t875 + t940;
	t821 = t863 + (qJD(3) * t963 - t859) * t880 + (-qJD(3) * t910 - t913) * t876;
	t819 = t876 * t924 + t859 * t880 + t910 * t955 - t863 + (-t858 * t876 - t860 * t954) * t872;
	t818 = t978 * t876 + t979 * t880;
	t814 = -t841 * t947 + t842 * t875 - t891;
	t809 = t918 * qJD(4) + t849 * t957 - t850 * t973;
	t808 = t833 * qJD(4) - t849 * t947 - t850 * t875;
	t802 = -t827 * qJD(4) - t849 * t973 - t850 * t957;
	t801 = t828 * qJD(4) - t849 * t875 + t850 * t947;
	t793 = -t849 * t960 + t809 * t879 + (-t833 * t874 + t853 * t959) * qJD(5);
	t791 = t919 * qJD(4) - t819 * t957 + t886;
	t790 = t823 * qJD(4) + t819 * t947 + t875 * t985;
	t789 = -qJD(4) * t940 - t818 * t957 - t841 * t953 - t882;
	t788 = t818 * t947 - t871 * t945 + t982;
	t787 = t928 * qJD(5) + t802 * t879 + t850 * t960;
	t785 = t901 * qJD(4) + t821 * t973 - t994;
	t784 = t821 * t875 + t916 * t949 + t947 * t985 + t998;
	t783 = -t810 * qJD(4) + t819 * t973 + t994;
	t782 = t819 * t875 - t871 * t886 - t973 * t909 - t998;
	t781 = qJD(4) * t891 + t818 * t973 + t982 * t871 - t875 * t890 - t945;
	t780 = t815 * qJD(4) + t818 * t875 + t871 * t882 + t973 * t890;
	t779 = t819 * t960 + t791 * t879 + (-t823 * t874 - t840 * t959) * qJD(5);
	t777 = t818 * t960 + t789 * t879 + (-t825 * t874 + t842 * t959) * qJD(5);
	t775 = t785 * t879 - t1001;
	t773 = t783 * t879 + t1001;
	t771 = t980 * qJD(5) + t781 * t879 + t881 * t874;
	t770 = t798 * qJD(5) + t781 * t874 - t881 * t879;
	t769 = t771 * t878 + t780 * t873 + (-t798 * t873 + t814 * t878) * qJD(6);
	t768 = -t771 * t873 + t780 * t878 + (-t798 * t878 - t814 * t873) * qJD(6);
	t1 = [(t775 * t878 + t784 * t873) * r_i_i_C(1) + (-t775 * t873 + t784 * t878) * r_i_i_C(2) + t775 * pkin(5) + t785 * pkin(4) + t784 * pkin(12) + t821 * pkin(3) - t859 * pkin(2) - pkin(10) * t944 + qJD(2) * t948 + t975 * (t785 * t874 + t1002) + ((-t878 * t901 - t1004) * r_i_i_C(1) + (t873 * t901 - t1003) * r_i_i_C(2)) * qJD(6) + t805 * pkin(11) + (-t974 * pkin(1) - t983 * t958) * qJD(1), t938, (t777 * t878 + t788 * t873 + (-t800 * t873 + t824 * t878) * qJD(6)) * r_i_i_C(1) + (-t777 * t873 + t788 * t878 + (-t800 * t878 - t824 * t873) * qJD(6)) * r_i_i_C(2) + t777 * pkin(5) + t789 * pkin(4) + t788 * pkin(12) - t885 * pkin(3) + t818 * t971 + t975 * (t800 * qJD(5) + t789 * t874 - t818 * t959), (t781 * t873 + t815 * t951) * r_i_i_C(1) + (t781 * t878 - t815 * t952) * r_i_i_C(2) + t781 * pkin(12) + t905 * t780 + t888 * t814, -t927 * t770 + t975 * t771 - t980 * t921, r_i_i_C(1) * t768 - r_i_i_C(2) * t769; -qJD(1) * pkin(10) * t898 - pkin(1) * t956 + pkin(2) * t904 + t818 * pkin(3) + t781 * pkin(4) + t771 * pkin(5) + t881 * pkin(11) + t780 * pkin(12) + t769 * r_i_i_C(1) + t768 * r_i_i_C(2) + qJD(2) * t958 + t975 * t770 + t983 * t938, t946, (t779 * t878 + t790 * t873 + (-t799 * t873 - t878 * t919) * qJD(6)) * r_i_i_C(1) + (-t779 * t873 + t790 * t878 + (-t799 * t878 + t873 * t919) * qJD(6)) * r_i_i_C(2) + t779 * pkin(5) + t791 * pkin(4) + t790 * pkin(12) + t985 * pkin(3) + t819 * t971 + t975 * (t799 * qJD(5) + t791 * t874 - t819 * t959), (t783 * t873 - t813 * t951) * r_i_i_C(1) + (t783 * t878 + t813 * t952) * r_i_i_C(2) + t783 * pkin(12) + t905 * t782 + t888 * t810, t975 * t773 - t929 * t921 + t927 * (-t783 * t874 + t1002), (-t773 * t873 + t782 * t878) * r_i_i_C(1) + (-t773 * t878 - t782 * t873) * r_i_i_C(2) + ((-t810 * t873 + t1003) * r_i_i_C(1) + (-t810 * t878 - t1004) * r_i_i_C(2)) * qJD(6); 0, 0, (t793 * t878 + t808 * t873) * r_i_i_C(1) + (-t793 * t873 + t808 * t878) * r_i_i_C(2) + t793 * pkin(5) + t809 * pkin(4) + t808 * pkin(12) - t850 * pkin(3) - t849 * t971 + t975 * (t826 * qJD(5) + t809 * t874 + t849 * t959) + ((-t826 * t873 - t878 * t918) * r_i_i_C(1) + (-t826 * t878 + t873 * t918) * r_i_i_C(2)) * qJD(6), (t802 * t873 + t828 * t951) * r_i_i_C(1) + (t802 * t878 - t828 * t952) * r_i_i_C(2) + t802 * pkin(12) + t905 * t801 + t888 * t827, t975 * t787 - t928 * t921 + t927 * (-t807 * qJD(5) - t802 * t874 + t850 * t959), (-t787 * t873 + t801 * t878) * r_i_i_C(1) + (-t787 * t878 - t801 * t873) * r_i_i_C(2) + ((-t807 * t878 - t827 * t873) * r_i_i_C(1) + (t807 * t873 - t827 * t878) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end