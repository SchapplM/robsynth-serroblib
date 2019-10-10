% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:23
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobiaD_transl_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:36
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(14));
	t50 = sin(pkin(14));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:36
	% EndTime: 2019-10-09 23:23:37
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (95->48), mult. (333->96), div. (0->0), fcn. (328->10), ass. (0->36)
	t218 = sin(pkin(7));
	t248 = (pkin(10) + r_i_i_C(3)) * t218;
	t217 = sin(pkin(14));
	t220 = cos(pkin(14));
	t224 = sin(qJ(2));
	t222 = cos(pkin(6));
	t226 = cos(qJ(2));
	t241 = t222 * t226;
	t211 = -t217 * t224 + t220 * t241;
	t219 = sin(pkin(6));
	t245 = t218 * t219;
	t221 = cos(pkin(7));
	t223 = sin(qJ(3));
	t244 = t221 * t223;
	t225 = cos(qJ(3));
	t243 = t221 * t225;
	t242 = t222 * t224;
	t240 = t223 * t224;
	t239 = t223 * t226;
	t238 = t224 * t225;
	t237 = t225 * t226;
	t236 = t223 * t245;
	t235 = t225 * t245;
	t233 = r_i_i_C(1) * t223 + r_i_i_C(2) * t225;
	t232 = r_i_i_C(1) * t225 - r_i_i_C(2) * t223 + pkin(2);
	t212 = t217 * t226 + t220 * t242;
	t231 = t217 * t241 + t220 * t224;
	t230 = t217 * t242 - t220 * t226;
	t229 = t233 * t221 - t248;
	t228 = (-t221 * t238 - t239) * r_i_i_C(1) + (t221 * t240 - t237) * r_i_i_C(2);
	t227 = (-t221 * t239 - t238) * r_i_i_C(1) + (-t221 * t237 + t240) * r_i_i_C(2);
	t210 = t230 * qJD(2);
	t209 = t231 * qJD(2);
	t208 = t212 * qJD(2);
	t207 = t211 * qJD(2);
	t1 = [0, t232 * t210 + t229 * t209 + ((t223 * t231 + t230 * t243) * r_i_i_C(1) + (t225 * t231 - t230 * t244) * r_i_i_C(2)) * qJD(3), (t209 * t223 + t210 * t243) * r_i_i_C(1) + (t209 * t225 - t210 * t244) * r_i_i_C(2) + ((-t217 * t236 + t225 * t230 + t231 * t244) * r_i_i_C(1) + (-t217 * t235 - t223 * t230 + t231 * t243) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, -t232 * t208 - t229 * t207 + ((-t211 * t223 - t212 * t243) * r_i_i_C(1) + (-t211 * t225 + t212 * t244) * r_i_i_C(2)) * qJD(3), (-t207 * t223 - t208 * t243) * r_i_i_C(1) + (-t207 * t225 + t208 * t244) * r_i_i_C(2) + ((-t211 * t244 - t212 * t225 + t220 * t236) * r_i_i_C(1) + (-t211 * t243 + t212 * t223 + t220 * t235) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, (t228 * qJD(3) + (-t224 * pkin(2) + t226 * t248 + t227) * qJD(2)) * t219, -t233 * t222 * t218 * qJD(3) + (t228 * qJD(2) + t227 * qJD(3)) * t219, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:39
	% EndTime: 2019-10-09 23:23:40
	% DurationCPUTime: 1.17s
	% Computational Cost: add. (594->150), mult. (2011->271), div. (0->0), fcn. (2160->14), ass. (0->89)
	t508 = pkin(11) + r_i_i_C(3);
	t455 = sin(qJ(3));
	t458 = cos(qJ(3));
	t446 = sin(pkin(14));
	t450 = cos(pkin(14));
	t459 = cos(qJ(2));
	t453 = cos(pkin(6));
	t456 = sin(qJ(2));
	t492 = t453 * t456;
	t465 = t446 * t492 - t450 * t459;
	t452 = cos(pkin(7));
	t491 = t453 * t459;
	t466 = t446 * t491 + t450 * t456;
	t448 = sin(pkin(7));
	t449 = sin(pkin(6));
	t501 = t448 * t449;
	t467 = t446 * t501 - t452 * t466;
	t416 = t467 * t455 - t458 * t465;
	t447 = sin(pkin(8));
	t451 = cos(pkin(8));
	t454 = sin(qJ(4));
	t457 = cos(qJ(4));
	t460 = t508 * t447 - (t454 * r_i_i_C(1) + t457 * r_i_i_C(2)) * t451;
	t507 = t448 * pkin(10);
	t437 = t446 * t459 + t450 * t492;
	t432 = t437 * qJD(2);
	t506 = t432 * t455;
	t505 = t437 * t458;
	t503 = t447 * t454;
	t502 = t447 * t457;
	t500 = t448 * t451;
	t499 = t448 * t453;
	t498 = t449 * t452;
	t497 = t449 * t459;
	t496 = t451 * t454;
	t495 = t451 * t457;
	t494 = t452 * t455;
	t493 = t452 * t458;
	t490 = t455 * t456;
	t489 = t455 * t459;
	t488 = t456 * t458;
	t487 = t458 * t459;
	t486 = qJD(2) * t456;
	t485 = qJD(2) * t459;
	t484 = qJD(3) * t455;
	t483 = qJD(3) * t458;
	t482 = qJD(4) * t456;
	t481 = t448 * t503;
	t480 = t448 * t502;
	t479 = t450 * t501;
	t478 = t450 * t491;
	t477 = t455 * t499;
	t476 = t449 * t486;
	t474 = t447 * t448 * t476;
	t473 = t457 * r_i_i_C(1) - t454 * r_i_i_C(2) + pkin(3);
	t436 = -t446 * t456 + t478;
	t419 = -t436 * t455 - t437 * t493;
	t470 = -t436 * t458 + t437 * t494;
	t469 = -t436 * t452 + t479;
	t421 = t455 * t466 + t465 * t493;
	t468 = t458 * t466 - t465 * t494;
	t464 = t452 * t487 - t490;
	t463 = -t452 * t488 - t489;
	t462 = t452 * t489 + t488;
	t461 = t452 * t490 - t487;
	t435 = -t448 * t497 + t453 * t452;
	t434 = t465 * qJD(2);
	t433 = t466 * qJD(2);
	t431 = -qJD(2) * t478 + t446 * t486;
	t428 = t461 * t449;
	t427 = t463 * t449;
	t426 = t446 * t498 + t448 * t466;
	t425 = -t436 * t448 - t450 * t498;
	t424 = t462 * t449 + t477;
	t423 = t464 * t449 + t458 * t499;
	t415 = t455 * t465 + t467 * t458;
	t414 = t436 * t494 - t455 * t479 + t505;
	t413 = -t437 * t455 - t469 * t458;
	t412 = -t476 * t494 - t449 * t456 * t484 + (t449 * t485 + (t452 * t497 + t499) * qJD(3)) * t458;
	t411 = -qJD(3) * t477 + (t463 * qJD(2) - t462 * qJD(3)) * t449;
	t410 = t421 * qJD(3) + t433 * t494 + t434 * t458;
	t409 = t468 * qJD(3) + t433 * t493 - t434 * t455;
	t408 = t419 * qJD(3) + t431 * t494 - t432 * t458;
	t407 = t470 * qJD(3) + t431 * t493 + t506;
	t406 = t434 * t494 + t465 * t484 + (t467 * qJD(3) - t433) * t458;
	t405 = -t416 * qJD(3) + t433 * t455 + t434 * t493;
	t404 = -t431 * t458 - t437 * t484 - t479 * t483 + (t436 * t483 - t506) * t452;
	t403 = -t432 * t493 + t431 * t455 + (t469 * t455 - t505) * qJD(3);
	t1 = [0, (t409 * t496 + t410 * t457 - t433 * t481) * r_i_i_C(1) + (t409 * t495 - t410 * t454 - t433 * t480) * r_i_i_C(2) + t410 * pkin(3) + t434 * pkin(2) - t433 * t507 + ((t421 * t495 + t454 * t468 - t465 * t480) * r_i_i_C(1) + (-t421 * t496 + t457 * t468 + t465 * t481) * r_i_i_C(2)) * qJD(4) + t508 * (-t409 * t447 - t433 * t500), t473 * t405 + t460 * t406 + ((-t415 * t454 - t416 * t495) * r_i_i_C(1) + (-t415 * t457 + t416 * t496) * r_i_i_C(2)) * qJD(4), (t405 * t495 - t406 * t454 - t434 * t480) * r_i_i_C(1) + (-t405 * t496 - t406 * t457 + t434 * t481) * r_i_i_C(2) + ((-t415 * t496 - t416 * t457 - t426 * t503) * r_i_i_C(1) + (-t415 * t495 + t416 * t454 - t426 * t502) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (t407 * t496 + t408 * t457 - t431 * t481) * r_i_i_C(1) + (t407 * t495 - t408 * t454 - t431 * t480) * r_i_i_C(2) + t408 * pkin(3) - t432 * pkin(2) - t431 * t507 + ((t419 * t495 + t437 * t480 + t454 * t470) * r_i_i_C(1) + (-t419 * t496 - t437 * t481 + t457 * t470) * r_i_i_C(2)) * qJD(4) + t508 * (-t407 * t447 - t431 * t500), t473 * t403 + t460 * t404 + ((-t413 * t454 - t414 * t495) * r_i_i_C(1) + (-t413 * t457 + t414 * t496) * r_i_i_C(2)) * qJD(4), (t403 * t495 - t404 * t454 + t432 * t480) * r_i_i_C(1) + (-t403 * t496 - t404 * t457 - t432 * t481) * r_i_i_C(2) + ((-t413 * t496 - t414 * t457 - t425 * t503) * r_i_i_C(1) + (-t413 * t495 + t414 * t454 - t425 * t502) * r_i_i_C(2)) * qJD(4), 0, 0; 0, ((t427 * t495 + t428 * t454) * r_i_i_C(1) + (-t427 * t496 + t428 * t457) * r_i_i_C(2)) * qJD(4) + (t473 * (-t462 * qJD(2) + t463 * qJD(3)) - t460 * (-t464 * qJD(2) + t461 * qJD(3)) - pkin(2) * t486 + ((t508 * t451 + pkin(10)) * t485 + ((t454 * t485 + t457 * t482) * r_i_i_C(1) + (-t454 * t482 + t457 * t485) * r_i_i_C(2)) * t447) * t448) * t449, t473 * t411 + t460 * t412 + ((-t423 * t454 - t424 * t495) * r_i_i_C(1) + (-t423 * t457 + t424 * t496) * r_i_i_C(2)) * qJD(4), (t411 * t495 - t412 * t454 + t457 * t474) * r_i_i_C(1) + (-t411 * t496 - t412 * t457 - t454 * t474) * r_i_i_C(2) + ((-t423 * t496 - t424 * t457 - t435 * t503) * r_i_i_C(1) + (-t423 * t495 + t424 * t454 - t435 * t502) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:43
	% EndTime: 2019-10-09 23:23:45
	% DurationCPUTime: 2.02s
	% Computational Cost: add. (1882->235), mult. (6202->424), div. (0->0), fcn. (6949->16), ass. (0->141)
	t673 = sin(pkin(6));
	t676 = cos(pkin(7));
	t681 = sin(qJ(2));
	t684 = cos(qJ(3));
	t728 = t681 * t684;
	t680 = sin(qJ(3));
	t685 = cos(qJ(2));
	t729 = t680 * t685;
	t694 = t676 * t729 + t728;
	t672 = sin(pkin(7));
	t677 = cos(pkin(6));
	t740 = t672 * t677;
	t717 = t680 * t740;
	t648 = t694 * t673 + t717;
	t679 = sin(qJ(4));
	t683 = cos(qJ(4));
	t727 = t684 * t685;
	t730 = t680 * t681;
	t696 = t676 * t727 - t730;
	t647 = t696 * t673 + t684 * t740;
	t737 = t673 * t685;
	t659 = -t672 * t737 + t677 * t676;
	t671 = sin(pkin(8));
	t675 = cos(pkin(8));
	t710 = t647 * t675 + t659 * t671;
	t620 = t648 * t683 + t710 * t679;
	t670 = sin(pkin(14));
	t674 = cos(pkin(14));
	t732 = t677 * t681;
	t697 = t670 * t732 - t674 * t685;
	t731 = t677 * t685;
	t698 = t670 * t731 + t674 * t681;
	t742 = t672 * t673;
	t699 = t670 * t742 - t676 * t698;
	t639 = t699 * t680 - t684 * t697;
	t638 = t680 * t697 + t699 * t684;
	t739 = t673 * t676;
	t650 = t670 * t739 + t672 * t698;
	t711 = t638 * t675 + t650 * t671;
	t606 = t639 * t683 + t711 * t679;
	t718 = t674 * t731;
	t660 = -t670 * t681 + t718;
	t720 = t674 * t742;
	t734 = t676 * t680;
	t661 = t670 * t685 + t674 * t732;
	t745 = t661 * t684;
	t637 = t660 * t734 - t680 * t720 + t745;
	t701 = -t660 * t676 + t720;
	t636 = -t661 * t680 - t701 * t684;
	t649 = -t660 * t672 - t674 * t739;
	t712 = t636 * t675 + t649 * t671;
	t604 = t637 * t683 + t712 * t679;
	t752 = r_i_i_C(3) + pkin(12);
	t751 = t672 * pkin(10);
	t693 = t676 * t730 - t727;
	t640 = (-t696 * qJD(2) + t693 * qJD(3)) * t673;
	t748 = t640 * t671;
	t656 = t661 * qJD(2);
	t746 = t656 * t680;
	t743 = t671 * t672;
	t741 = t672 * t675;
	t738 = t673 * t681;
	t736 = t675 * t679;
	t735 = t675 * t683;
	t733 = t676 * t684;
	t726 = qJD(2) * t673;
	t725 = qJD(2) * t681;
	t724 = qJD(3) * t680;
	t723 = qJD(3) * t684;
	t678 = sin(qJ(5));
	t722 = qJD(5) * t678;
	t682 = cos(qJ(5));
	t721 = qJD(5) * t682;
	t719 = t672 * t738;
	t716 = t673 * t725;
	t715 = t685 * t726;
	t714 = t672 * t716;
	t713 = t672 * t715;
	t709 = t682 * r_i_i_C(1) - t678 * r_i_i_C(2) + pkin(4);
	t655 = -qJD(2) * t718 + t670 * t725;
	t615 = -t656 * t733 + t655 * t680 + (t701 * t680 - t745) * qJD(3);
	t708 = t615 * t675 + t656 * t743;
	t657 = t698 * qJD(2);
	t658 = t697 * qJD(2);
	t617 = -t639 * qJD(3) + t657 * t680 + t658 * t733;
	t707 = t617 * t675 - t658 * t743;
	t702 = -t660 * t684 + t661 * t734;
	t621 = t702 * qJD(3) + t655 * t733 + t746;
	t609 = -t621 * t671 - t655 * t741;
	t706 = -t621 * t675 + t655 * t743;
	t700 = t684 * t698 - t697 * t734;
	t623 = t700 * qJD(3) + t657 * t733 - t658 * t680;
	t610 = -t623 * t671 - t657 * t741;
	t705 = -t623 * t675 + t657 * t743;
	t613 = t636 * t683 - t637 * t736;
	t614 = t638 * t683 - t639 * t736;
	t642 = -t660 * t680 - t661 * t733;
	t704 = t642 * t675 + t661 * t743;
	t644 = t680 * t698 + t697 * t733;
	t703 = t644 * t675 - t697 * t743;
	t628 = t647 * t683 - t648 * t736;
	t695 = -t676 * t728 - t729;
	t692 = qJD(5) * (-t678 * r_i_i_C(1) - t682 * r_i_i_C(2));
	t651 = t695 * t673;
	t691 = t651 * t675 + t671 * t719;
	t633 = -qJD(3) * t717 + (t695 * qJD(2) - t694 * qJD(3)) * t673;
	t690 = t633 * t675 + t671 * t714;
	t689 = t640 * t675 + t671 * t713;
	t688 = -t637 * t679 + t712 * t683;
	t687 = -t639 * t679 + t711 * t683;
	t686 = -t648 * t679 + t710 * t683;
	t611 = t704 * t679 - t683 * t702;
	t612 = t703 * t679 - t683 * t700;
	t652 = t693 * t673;
	t631 = -t652 * t683 + t691 * t679;
	t646 = -t651 * t671 + t675 * t719;
	t641 = (-t694 * qJD(2) + t695 * qJD(3)) * t673;
	t635 = -t647 * t671 + t659 * t675;
	t634 = -t716 * t734 - t724 * t738 + (t715 + (t676 * t737 + t740) * qJD(3)) * t684;
	t632 = t675 * t713 - t748;
	t630 = -t644 * t671 - t697 * t741;
	t629 = -t642 * t671 + t661 * t741;
	t627 = -t633 * t671 + t675 * t714;
	t626 = -t638 * t671 + t650 * t675;
	t625 = -t636 * t671 + t649 * t675;
	t624 = t644 * qJD(3) + t657 * t734 + t658 * t684;
	t622 = t642 * qJD(3) + t655 * t734 - t656 * t684;
	t618 = t658 * t734 + t697 * t724 + (t699 * qJD(3) - t657) * t684;
	t616 = -t655 * t684 - t661 * t724 - t720 * t723 + (t660 * t723 - t746) * t676;
	t608 = -t617 * t671 - t658 * t741;
	t607 = -t615 * t671 + t656 * t741;
	t602 = t641 * t683 + t689 * t679 + (t652 * t679 + t691 * t683) * qJD(4);
	t600 = -t634 * t736 + t633 * t683 + (-t647 * t679 - t648 * t735) * qJD(4);
	t598 = t686 * qJD(4) + t634 * t683 + t690 * t679;
	t596 = -t618 * t736 + t617 * t683 + (-t638 * t679 - t639 * t735) * qJD(4);
	t594 = -t616 * t736 + t615 * t683 + (-t636 * t679 - t637 * t735) * qJD(4);
	t592 = t624 * t683 - t705 * t679 + (t679 * t700 + t703 * t683) * qJD(4);
	t590 = t622 * t683 - t706 * t679 + (t679 * t702 + t704 * t683) * qJD(4);
	t588 = t687 * qJD(4) + t618 * t683 + t707 * t679;
	t586 = t688 * qJD(4) + t616 * t683 + t708 * t679;
	t1 = [0, (t592 * t682 + t610 * t678) * r_i_i_C(1) + (-t592 * t678 + t610 * t682) * r_i_i_C(2) + t592 * pkin(4) + t624 * pkin(3) + t658 * pkin(2) - t657 * t751 + t752 * (t612 * qJD(4) + t624 * t679 + t705 * t683) + ((-t612 * t678 + t630 * t682) * r_i_i_C(1) + (-t612 * t682 - t630 * t678) * r_i_i_C(2)) * qJD(5) + t610 * pkin(11), (t596 * t682 - t614 * t722) * r_i_i_C(1) + (-t596 * t678 - t614 * t721) * r_i_i_C(2) + t596 * pkin(4) + t617 * pkin(3) + t752 * (t614 * qJD(4) + t617 * t679 + t618 * t735) + ((t618 * t678 + t639 * t721) * r_i_i_C(1) + (t618 * t682 - t639 * t722) * r_i_i_C(2) + t618 * pkin(11)) * t671, t752 * t588 + t687 * t692 + t709 * (-t606 * qJD(4) - t618 * t679 + t707 * t683), (-t588 * t678 + t608 * t682) * r_i_i_C(1) + (-t588 * t682 - t608 * t678) * r_i_i_C(2) + ((-t606 * t682 - t626 * t678) * r_i_i_C(1) + (t606 * t678 - t626 * t682) * r_i_i_C(2)) * qJD(5), 0; 0, (t590 * t682 + t609 * t678) * r_i_i_C(1) + (-t590 * t678 + t609 * t682) * r_i_i_C(2) + t590 * pkin(4) + t622 * pkin(3) - t656 * pkin(2) - t655 * t751 + t752 * (t611 * qJD(4) + t622 * t679 + t706 * t683) + ((-t611 * t678 + t629 * t682) * r_i_i_C(1) + (-t611 * t682 - t629 * t678) * r_i_i_C(2)) * qJD(5) + t609 * pkin(11), (t594 * t682 - t613 * t722) * r_i_i_C(1) + (-t594 * t678 - t613 * t721) * r_i_i_C(2) + t594 * pkin(4) + t615 * pkin(3) + t752 * (t613 * qJD(4) + t615 * t679 + t616 * t735) + ((t616 * t678 + t637 * t721) * r_i_i_C(1) + (t616 * t682 - t637 * t722) * r_i_i_C(2) + t616 * pkin(11)) * t671, t752 * t586 + t688 * t692 + t709 * (-t604 * qJD(4) - t616 * t679 + t708 * t683), (-t586 * t678 + t607 * t682) * r_i_i_C(1) + (-t586 * t682 - t607 * t678) * r_i_i_C(2) + ((-t604 * t682 - t625 * t678) * r_i_i_C(1) + (t604 * t678 - t625 * t682) * r_i_i_C(2)) * qJD(5), 0; 0, (t602 * t682 + t632 * t678) * r_i_i_C(1) + (-t602 * t678 + t632 * t682) * r_i_i_C(2) + t602 * pkin(4) + t641 * pkin(3) - pkin(11) * t748 + t752 * (t631 * qJD(4) + t641 * t679 - t689 * t683) + ((-t631 * t678 + t646 * t682) * r_i_i_C(1) + (-t631 * t682 - t646 * t678) * r_i_i_C(2)) * qJD(5) + (-pkin(2) * t681 + (pkin(11) * t675 + pkin(10)) * t685 * t672) * t726, (t600 * t682 - t628 * t722) * r_i_i_C(1) + (-t600 * t678 - t628 * t721) * r_i_i_C(2) + t600 * pkin(4) + t633 * pkin(3) + t752 * (t628 * qJD(4) + t633 * t679 + t634 * t735) + ((t634 * t678 + t648 * t721) * r_i_i_C(1) + (t634 * t682 - t648 * t722) * r_i_i_C(2) + t634 * pkin(11)) * t671, t752 * t598 + t686 * t692 + t709 * (-t620 * qJD(4) - t634 * t679 + t690 * t683), (-t598 * t678 + t627 * t682) * r_i_i_C(1) + (-t598 * t682 - t627 * t678) * r_i_i_C(2) + ((-t620 * t682 - t635 * t678) * r_i_i_C(1) + (t620 * t678 - t635 * t682) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:49
	% EndTime: 2019-10-09 23:23:54
	% DurationCPUTime: 4.32s
	% Computational Cost: add. (4964->335), mult. (16032->574), div. (0->0), fcn. (18521->18), ass. (0->192)
	t1000 = r_i_i_C(3) + pkin(13);
	t892 = sin(pkin(6));
	t894 = cos(pkin(7));
	t900 = sin(qJ(2));
	t904 = cos(qJ(3));
	t970 = t900 * t904;
	t899 = sin(qJ(3));
	t905 = cos(qJ(2));
	t971 = t899 * t905;
	t929 = t894 * t971 + t970;
	t918 = t929 * qJD(3);
	t930 = -t894 * t970 - t971;
	t909 = (t930 * qJD(2) - t918) * t892;
	t891 = sin(pkin(7));
	t895 = cos(pkin(6));
	t979 = t891 * t895;
	t960 = t899 * t979;
	t951 = qJD(3) * t960;
	t1006 = t909 - t951;
	t997 = cos(pkin(14));
	t954 = t997 * t905;
	t889 = sin(pkin(14));
	t985 = t889 * t900;
	t878 = t895 * t954 - t985;
	t874 = t878 * qJD(2);
	t955 = t997 * t900;
	t984 = t889 * t905;
	t921 = -t895 * t955 - t984;
	t875 = t921 * qJD(2);
	t956 = t892 * t997;
	t950 = t891 * t956;
	t922 = -t878 * t894 + t950;
	t973 = t894 * t904;
	t988 = t921 * t904;
	t908 = (t922 * t899 + t988) * qJD(3) - t874 * t899 + t875 * t973;
	t924 = t895 * t985 - t954;
	t923 = t895 * t984 + t955;
	t986 = t889 * t892;
	t933 = t891 * t986 - t894 * t923;
	t856 = t933 * t899 - t904 * t924;
	t876 = t923 * qJD(2);
	t877 = t924 * qJD(2);
	t1005 = -t856 * qJD(3) + t876 * t899 + t877 * t973;
	t896 = sin(qJ(6));
	t901 = cos(qJ(6));
	t927 = qJD(6) * (t896 * r_i_i_C(1) + t901 * r_i_i_C(2));
	t853 = t899 * t921 - t922 * t904;
	t974 = t894 * t899;
	t854 = t878 * t974 - t899 * t950 - t988;
	t898 = sin(qJ(4));
	t903 = cos(qJ(4));
	t890 = sin(pkin(8));
	t925 = -t878 * t891 - t894 * t956;
	t916 = t925 * t890;
	t893 = cos(pkin(8));
	t975 = t893 * t903;
	t1004 = t853 * t975 - t854 * t898 + t903 * t916;
	t972 = t899 * t900;
	t961 = t894 * t972;
	t969 = t904 * t905;
	t928 = t961 - t969;
	t871 = t928 * t892;
	t870 = t930 * t892;
	t978 = t891 * t900;
	t959 = t890 * t978;
	t953 = t892 * t959;
	t926 = t870 * t893 + t953;
	t1003 = t871 * t898 + t926 * t903;
	t935 = t904 * t923 - t924 * t974;
	t861 = t899 * t923 + t924 * t973;
	t981 = t891 * t890;
	t939 = t861 * t893 - t924 * t981;
	t1002 = t898 * t935 + t939 * t903;
	t936 = -t878 * t904 - t921 * t974;
	t859 = -t878 * t899 + t921 * t973;
	t940 = t859 * t893 - t921 * t981;
	t1001 = t898 * t936 + t940 * t903;
	t999 = pkin(10) * t891;
	t998 = pkin(11) * t890;
	t995 = t856 * t898;
	t931 = t894 * t969 - t972;
	t857 = (-t931 * qJD(2) + t928 * qJD(3)) * t892;
	t994 = t857 * t890;
	t868 = t929 * t892 + t960;
	t991 = t868 * t898;
	t989 = t875 * t899;
	t897 = sin(qJ(5));
	t983 = t890 * t897;
	t902 = cos(qJ(5));
	t982 = t890 * t902;
	t980 = t891 * t893;
	t977 = t892 * t905;
	t976 = t893 * t898;
	t968 = qJD(2) * t892;
	t967 = qJD(3) * t899;
	t966 = qJD(3) * t904;
	t965 = qJD(6) * t896;
	t964 = qJD(6) * t901;
	t962 = t903 * t981;
	t958 = t893 * t978;
	t957 = t905 * t968;
	t952 = t891 * t957;
	t947 = t898 * t951;
	t811 = t854 * t903 + (t853 * t893 + t916) * t898;
	t838 = -t853 * t890 + t925 * t893;
	t795 = t811 * t902 + t838 * t897;
	t946 = -t811 * t897 + t838 * t902;
	t855 = t899 * t924 + t933 * t904;
	t934 = t891 * t923 + t894 * t986;
	t920 = t934 * t890;
	t914 = t855 * t893 + t920;
	t813 = t856 * t903 + t914 * t898;
	t839 = -t855 * t890 + t934 * t893;
	t797 = t813 * t902 + t839 * t897;
	t945 = -t813 * t897 + t839 * t902;
	t819 = t940 * t898 - t903 * t936;
	t843 = -t859 * t890 - t921 * t980;
	t802 = t819 * t902 + t843 * t897;
	t821 = t939 * t898 - t903 * t935;
	t844 = -t861 * t890 - t924 * t980;
	t803 = t821 * t902 + t844 * t897;
	t867 = t931 * t892 + t904 * t979;
	t932 = -t891 * t977 + t895 * t894;
	t919 = t932 * t890;
	t913 = t867 * t893 + t919;
	t833 = t868 * t903 + t913 * t898;
	t852 = -t867 * t890 + t932 * t893;
	t809 = t833 * t902 + t852 * t897;
	t944 = -t833 * t897 + t852 * t902;
	t846 = -t871 * t903 + t926 * t898;
	t863 = -t870 * t890 + t892 * t958;
	t826 = t846 * t902 + t863 * t897;
	t943 = t901 * r_i_i_C(1) - t896 * r_i_i_C(2) + pkin(5);
	t823 = t853 * t903 - t854 * t976;
	t806 = t823 * t902 + t854 * t983;
	t825 = t855 * t903 - t856 * t976;
	t807 = t825 * t902 + t856 * t983;
	t834 = t936 * qJD(3) - t874 * t973 - t989;
	t816 = -t834 * t890 + t874 * t980;
	t942 = -t834 * t893 - t874 * t981;
	t836 = t935 * qJD(3) + t876 * t973 - t877 * t899;
	t817 = -t836 * t890 - t876 * t980;
	t941 = -t836 * t893 + t876 * t981;
	t842 = t867 * t903 - t868 * t976;
	t829 = t842 * t902 + t868 * t983;
	t822 = t853 * t898 + t854 * t975;
	t824 = t855 * t898 + t856 * t975;
	t841 = t867 * t898 + t868 * t975;
	t917 = t857 * t893 + t890 * t952;
	t912 = -t1000 * t897 - t943 * t902 - pkin(4);
	t907 = t902 * t927 + (-t1000 * t902 + t943 * t897) * qJD(5);
	t906 = t893 * t908;
	t858 = (-t929 * qJD(2) + t930 * qJD(3)) * t892;
	t851 = -t961 * t968 - t892 * t900 * t967 + (t957 + (t894 * t977 + t979) * qJD(3)) * t904;
	t847 = t893 * t952 - t994;
	t840 = t890 * t951 + (t890 * t918 + (-t930 * t890 + t958) * qJD(2)) * t892;
	t837 = t861 * qJD(3) + t876 * t974 + t877 * t904;
	t835 = t859 * qJD(3) - t874 * t974 + t875 * t904;
	t832 = -t867 * t975 - t903 * t919 + t991;
	t831 = t877 * t974 + t924 * t967 + (t933 * qJD(3) - t876) * t904;
	t830 = t874 * t904 + t921 * t967 - t950 * t966 + (t878 * t966 + t989) * t894;
	t815 = -t1005 * t890 - t877 * t980;
	t814 = -t875 * t980 - t908 * t890;
	t812 = -t855 * t975 - t903 * t920 + t995;
	t805 = t1003 * qJD(4) + t858 * t903 + t917 * t898;
	t804 = t846 * qJD(4) + t858 * t898 - t917 * t903;
	t801 = -t841 * qJD(4) + t1006 * t903 - t851 * t976;
	t800 = t842 * qJD(4) + t851 * t975 + t898 * t909 - t947;
	t799 = -t893 * t947 + t851 * t903 + (-t893 * t918 + (t930 * t893 + t959) * qJD(2)) * t898 * t892 + (t913 * t903 - t991) * qJD(4);
	t798 = -qJD(2) * t903 * t953 + t833 * qJD(4) - t1006 * t975 + t851 * t898;
	t793 = -t824 * qJD(4) + t1005 * t903 - t831 * t976;
	t792 = t825 * qJD(4) + t1005 * t898 + t831 * t975;
	t791 = -t822 * qJD(4) - t830 * t976 + t908 * t903;
	t790 = t823 * qJD(4) + t830 * t975 + t908 * t898;
	t789 = t1002 * qJD(4) + t837 * t903 - t941 * t898;
	t788 = t821 * qJD(4) + t837 * t898 + t941 * t903;
	t787 = t1001 * qJD(4) + t835 * t903 - t942 * t898;
	t786 = t819 * qJD(4) + t835 * t898 + t942 * t903;
	t785 = t805 * t902 + t847 * t897 + (-t846 * t897 + t863 * t902) * qJD(5);
	t783 = t851 * t983 + t801 * t902 + (-t842 * t897 + t868 * t982) * qJD(5);
	t781 = t831 * t903 + (t1005 * t893 - t877 * t981) * t898 + (t914 * t903 - t995) * qJD(4);
	t780 = t813 * qJD(4) - t1005 * t975 + t831 * t898 + t877 * t962;
	t779 = t830 * t903 + (-t875 * t981 + t906) * t898 + t1004 * qJD(4);
	t778 = qJD(4) * t811 + t830 * t898 + t875 * t962 - t903 * t906;
	t777 = t944 * qJD(5) + t799 * t902 + t840 * t897;
	t775 = t831 * t983 + t793 * t902 + (-t825 * t897 + t856 * t982) * qJD(5);
	t773 = t830 * t983 + t791 * t902 + (-t823 * t897 + t854 * t982) * qJD(5);
	t771 = t789 * t902 + t817 * t897 + (-t821 * t897 + t844 * t902) * qJD(5);
	t769 = t787 * t902 + t816 * t897 + (-t819 * t897 + t843 * t902) * qJD(5);
	t767 = t945 * qJD(5) + t781 * t902 + t815 * t897;
	t765 = t946 * qJD(5) + t779 * t902 + t814 * t897;
	t1 = [0, (t771 * t901 + t788 * t896) * r_i_i_C(1) + (-t771 * t896 + t788 * t901) * r_i_i_C(2) + t771 * pkin(5) + t789 * pkin(4) + t788 * pkin(12) + t837 * pkin(3) + t877 * pkin(2) - t876 * t999 + t1000 * (t803 * qJD(5) + t789 * t897 - t817 * t902) + ((-t1002 * t901 - t803 * t896) * r_i_i_C(1) + (t1002 * t896 - t803 * t901) * r_i_i_C(2)) * qJD(6) + t817 * pkin(11), (t775 * t901 + t792 * t896) * r_i_i_C(1) + (-t775 * t896 + t792 * t901) * r_i_i_C(2) + t775 * pkin(5) + t793 * pkin(4) + t792 * pkin(12) + t831 * t998 + t1000 * (t807 * qJD(5) + t793 * t897 - t831 * t982) + ((-t807 * t896 + t824 * t901) * r_i_i_C(1) + (-t807 * t901 - t824 * t896) * r_i_i_C(2)) * qJD(6) + t1005 * pkin(3), (t781 * t896 + t813 * t964) * r_i_i_C(1) + (t781 * t901 - t813 * t965) * r_i_i_C(2) + t781 * pkin(12) + t912 * t780 + t907 * t812, t1000 * t767 - t945 * t927 + t943 * (-t797 * qJD(5) - t781 * t897 + t815 * t902), (-t767 * t896 + t780 * t901) * r_i_i_C(1) + (-t767 * t901 - t780 * t896) * r_i_i_C(2) + ((-t797 * t901 - t812 * t896) * r_i_i_C(1) + (t797 * t896 - t812 * t901) * r_i_i_C(2)) * qJD(6); 0, (t769 * t901 + t786 * t896) * r_i_i_C(1) + (-t769 * t896 + t786 * t901) * r_i_i_C(2) + t769 * pkin(5) + t787 * pkin(4) + t786 * pkin(12) + t835 * pkin(3) + t875 * pkin(2) + t874 * t999 + t1000 * (t802 * qJD(5) + t787 * t897 - t816 * t902) + ((-t1001 * t901 - t802 * t896) * r_i_i_C(1) + (t1001 * t896 - t802 * t901) * r_i_i_C(2)) * qJD(6) + t816 * pkin(11), (t773 * t901 + t790 * t896 + (-t806 * t896 + t822 * t901) * qJD(6)) * r_i_i_C(1) + (-t773 * t896 + t790 * t901 + (-t806 * t901 - t822 * t896) * qJD(6)) * r_i_i_C(2) + t773 * pkin(5) + t791 * pkin(4) + t790 * pkin(12) + t908 * pkin(3) + t830 * t998 + t1000 * (t806 * qJD(5) + t791 * t897 - t830 * t982), (t779 * t896 + t811 * t964) * r_i_i_C(1) + (t779 * t901 - t811 * t965) * r_i_i_C(2) + t779 * pkin(12) + t912 * t778 - t907 * t1004, t1000 * t765 - t946 * t927 + t943 * (-t795 * qJD(5) - t779 * t897 + t814 * t902), (-t765 * t896 + t778 * t901) * r_i_i_C(1) + (-t765 * t901 - t778 * t896) * r_i_i_C(2) + ((t1004 * t896 - t795 * t901) * r_i_i_C(1) + (t1004 * t901 + t795 * t896) * r_i_i_C(2)) * qJD(6); 0, (t785 * t901 + t804 * t896) * r_i_i_C(1) + (-t785 * t896 + t804 * t901) * r_i_i_C(2) + t785 * pkin(5) + t805 * pkin(4) + t804 * pkin(12) + t858 * pkin(3) - pkin(11) * t994 + t1000 * (t826 * qJD(5) + t805 * t897 - t847 * t902) + ((-t1003 * t901 - t826 * t896) * r_i_i_C(1) + (t1003 * t896 - t826 * t901) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t900 + (pkin(11) * t893 + pkin(10)) * t905 * t891) * t968, (t783 * t901 + t800 * t896) * r_i_i_C(1) + (-t783 * t896 + t800 * t901) * r_i_i_C(2) + t783 * pkin(5) + t801 * pkin(4) + t800 * pkin(12) + t851 * t998 + t1000 * (t829 * qJD(5) + t801 * t897 - t851 * t982) + ((-t829 * t896 + t841 * t901) * r_i_i_C(1) + (-t829 * t901 - t841 * t896) * r_i_i_C(2)) * qJD(6) + t1006 * pkin(3), (t799 * t896 + t833 * t964) * r_i_i_C(1) + (t799 * t901 - t833 * t965) * r_i_i_C(2) + t799 * pkin(12) + t912 * t798 + t907 * t832, t1000 * t777 - t944 * t927 + t943 * (-t809 * qJD(5) - t799 * t897 + t840 * t902), (-t777 * t896 + t798 * t901) * r_i_i_C(1) + (-t777 * t901 - t798 * t896) * r_i_i_C(2) + ((-t809 * t901 - t832 * t896) * r_i_i_C(1) + (t809 * t896 - t832 * t901) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end