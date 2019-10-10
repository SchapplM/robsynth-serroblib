% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_transl_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
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
	% StartTime: 2019-10-10 11:10:59
	% EndTime: 2019-10-10 11:10:59
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (35->18), mult. (110->35), div. (0->0), fcn. (94->6), ass. (0->20)
	t136 = sin(pkin(6));
	t151 = t136 * (pkin(10) + r_i_i_C(3));
	t138 = sin(qJ(2));
	t139 = sin(qJ(1));
	t149 = t138 * t139;
	t141 = cos(qJ(1));
	t148 = t138 * t141;
	t140 = cos(qJ(2));
	t147 = t139 * t140;
	t146 = t140 * t141;
	t137 = cos(pkin(6));
	t145 = -t137 * t146 + t149;
	t144 = t137 * t147 + t148;
	t143 = t137 * t148 + t147;
	t142 = t137 * t149 - t146;
	t135 = t142 * qJD(1) + t145 * qJD(2);
	t134 = t144 * qJD(1) + t143 * qJD(2);
	t133 = t143 * qJD(1) + t144 * qJD(2);
	t132 = t145 * qJD(1) + t142 * qJD(2);
	t1 = [t135 * r_i_i_C(1) + t134 * r_i_i_C(2) + (-pkin(1) * t141 - t139 * t151) * qJD(1), t132 * r_i_i_C(1) + t133 * r_i_i_C(2), 0, 0, 0, 0; -t133 * r_i_i_C(1) + t132 * r_i_i_C(2) + (-pkin(1) * t139 + t141 * t151) * qJD(1), -t134 * r_i_i_C(1) + t135 * r_i_i_C(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t138 - r_i_i_C(2) * t140) * t136 * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:59
	% EndTime: 2019-10-10 11:11:00
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (145->52), mult. (459->96), div. (0->0), fcn. (432->10), ass. (0->34)
	t237 = sin(pkin(14));
	t238 = sin(pkin(7));
	t240 = cos(pkin(14));
	t265 = (r_i_i_C(1) * t237 + r_i_i_C(2) * t240) * t238 + pkin(10);
	t264 = r_i_i_C(3) + qJ(3);
	t242 = cos(pkin(6));
	t245 = cos(qJ(2));
	t246 = cos(qJ(1));
	t253 = t246 * t245;
	t243 = sin(qJ(2));
	t244 = sin(qJ(1));
	t256 = t244 * t243;
	t247 = t242 * t256 - t253;
	t250 = -t242 * t253 + t256;
	t233 = t250 * qJD(1) + t247 * qJD(2);
	t263 = t233 * t238;
	t241 = cos(pkin(7));
	t262 = t237 * t241;
	t239 = sin(pkin(6));
	t261 = t239 * t244;
	t260 = t239 * t246;
	t259 = t240 * t241;
	t258 = t241 * t245;
	t257 = t243 * t238;
	t255 = t244 * t245;
	t254 = t246 * t243;
	t252 = qJD(1) * t239 * t241;
	t249 = t242 * t255 + t254;
	t248 = t242 * t254 + t255;
	t236 = t247 * qJD(1) + t250 * qJD(2);
	t235 = t249 * qJD(1) + t248 * qJD(2);
	t234 = t248 * qJD(1) + t249 * qJD(2);
	t232 = t246 * t252 - t263;
	t1 = [(t235 * t262 + t236 * t240) * r_i_i_C(1) + (t235 * t259 - t236 * t237) * r_i_i_C(2) + t236 * pkin(2) + t241 * qJD(3) * t260 + (-t250 * qJD(3) - t264 * t235) * t238 + (-t246 * pkin(1) + (-t264 * t241 - t265) * t261) * qJD(1), (t233 * t240 + t234 * t262) * r_i_i_C(1) + (-t233 * t237 + t234 * t259) * r_i_i_C(2) + t233 * pkin(2) + (-t247 * qJD(3) - t264 * t234) * t238, t232, 0, 0, 0; (t233 * t262 - t234 * t240) * r_i_i_C(1) + (t233 * t259 + t234 * t237) * r_i_i_C(2) + t232 * r_i_i_C(3) - t234 * pkin(2) - qJ(3) * t263 + (t249 * t238 + t241 * t261) * qJD(3) + (-t244 * pkin(1) + (qJ(3) * t241 + t265) * t260) * qJD(1), (-t235 * t240 + t236 * t262) * r_i_i_C(1) + (t235 * t237 + t236 * t259) * r_i_i_C(2) - t235 * pkin(2) + (t248 * qJD(3) - t264 * t236) * t238, t235 * t238 + t244 * t252, 0, 0, 0; 0, (qJD(3) * t257 + ((-t237 * t258 - t240 * t243) * r_i_i_C(1) + (t237 * t243 - t240 * t258) * r_i_i_C(2) - t243 * pkin(2) + t264 * t245 * t238) * qJD(2)) * t239, t239 * qJD(2) * t257, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:11:02
	% EndTime: 2019-10-10 11:11:03
	% DurationCPUTime: 1.39s
	% Computational Cost: add. (651->136), mult. (2094->256), div. (0->0), fcn. (2168->14), ass. (0->90)
	t495 = sin(qJ(1));
	t494 = sin(qJ(2));
	t546 = cos(pkin(6));
	t518 = t495 * t546;
	t514 = t494 * t518;
	t525 = qJD(2) * t494;
	t497 = cos(qJ(2));
	t498 = cos(qJ(1));
	t528 = t498 * t497;
	t465 = -qJD(1) * t514 - t495 * t525 + (qJD(2) * t546 + qJD(1)) * t528;
	t486 = sin(pkin(14));
	t490 = cos(pkin(14));
	t517 = t498 * t546;
	t479 = t494 * t517 + t495 * t497;
	t500 = t498 * t494 + t497 * t518;
	t464 = t500 * qJD(1) + t479 * qJD(2);
	t488 = sin(pkin(7));
	t492 = cos(pkin(7));
	t489 = sin(pkin(6));
	t527 = qJD(1) * t489;
	t521 = t495 * t527;
	t502 = -t464 * t492 + t488 * t521;
	t441 = -t465 * t486 + t502 * t490;
	t442 = t465 * t490 + t502 * t486;
	t457 = t464 * t488 + t492 * t521;
	t493 = sin(qJ(4));
	t491 = cos(pkin(8));
	t496 = cos(qJ(4));
	t531 = t491 * t496;
	t487 = sin(pkin(8));
	t539 = t487 * t496;
	t555 = t441 * t531 - t442 * t493 + t457 * t539;
	t532 = t491 * t493;
	t540 = t487 * t493;
	t554 = -t441 * t532 - t442 * t496 - t457 * t540;
	t478 = t495 * t494 - t497 * t517;
	t535 = t489 * t498;
	t506 = t478 * t492 + t488 * t535;
	t451 = t479 * t486 + t506 * t490;
	t452 = -t479 * t490 + t506 * t486;
	t468 = -t478 * t488 + t492 * t535;
	t553 = t451 * t531 - t452 * t493 + t468 * t539;
	t552 = t451 * t532 + t452 * t496 + t468 * t540;
	t547 = pkin(11) + r_i_i_C(3);
	t513 = t547 * t491 + qJ(3);
	t541 = t486 * t492;
	t538 = t488 * qJ(3);
	t537 = t488 * t489;
	t536 = t489 * t495;
	t534 = t490 * t492;
	t533 = t490 * t494;
	t530 = t492 * t494;
	t529 = t492 * t497;
	t524 = qJD(2) * t497;
	t523 = qJD(4) * t493;
	t522 = qJD(4) * t496;
	t520 = t498 * t527;
	t519 = qJ(3) * t492 + pkin(10);
	t516 = t546 * t488;
	t515 = t525 * t537;
	t512 = t487 * t515;
	t463 = t479 * qJD(1) + t500 * qJD(2);
	t501 = t514 - t528;
	t462 = t478 * qJD(1) + t501 * qJD(2);
	t503 = t462 * t492 + t488 * t520;
	t439 = t463 * t486 + t503 * t490;
	t455 = -t462 * t488 + t492 * t520;
	t511 = t439 * t491 + t455 * t487;
	t470 = t488 * t500 + t492 * t536;
	t505 = t488 * t536 - t492 * t500;
	t508 = (t486 * t501 + t505 * t490) * t491 + t470 * t487;
	t507 = t496 * r_i_i_C(1) - t493 * r_i_i_C(2) + pkin(3);
	t504 = -t486 * t494 + t490 * t529;
	t475 = (-t486 * t497 - t490 * t530) * t489;
	t476 = (-t486 * t530 + t490 * t497) * t489;
	t499 = (t493 * r_i_i_C(1) + t496 * r_i_i_C(2)) * t491 - t547 * t487;
	t477 = t546 * t492 - t497 * t537;
	t474 = qJD(2) * t476;
	t473 = qJD(2) * t475;
	t467 = t489 * t533 + (t489 * t529 + t516) * t486;
	t466 = t504 * t489 + t490 * t516;
	t461 = -t490 * t500 + t501 * t541;
	t460 = t486 * t500 + t501 * t534;
	t459 = -t478 * t490 - t479 * t541;
	t458 = t478 * t486 - t479 * t534;
	t454 = t505 * t486 - t490 * t501;
	t440 = -t463 * t490 + t503 * t486;
	t438 = t440 * t496 + t511 * t493 + (-t454 * t493 + t508 * t496) * qJD(4);
	t437 = -t440 * t493 + t511 * t496 + (-t454 * t496 - t508 * t493) * qJD(4);
	t1 = [t554 * r_i_i_C(1) - t555 * r_i_i_C(2) - t442 * pkin(3) - t465 * pkin(2) - t464 * t538 + t468 * qJD(3) + (-t498 * pkin(1) - t519 * t536) * qJD(1) + (t553 * r_i_i_C(1) - t552 * r_i_i_C(2)) * qJD(4) + t547 * (t441 * t487 - t457 * t491), t462 * pkin(2) + t507 * (t462 * t490 + t463 * t541) + t499 * (-t462 * t486 + t463 * t534) + ((t460 * t531 - t461 * t493) * r_i_i_C(1) + (-t460 * t532 - t461 * t496) * r_i_i_C(2)) * qJD(4) + (-t501 * qJD(3) - t513 * t463 + ((-t463 * t493 - t501 * t522) * r_i_i_C(1) + (-t463 * t496 + t501 * t523) * r_i_i_C(2)) * t487) * t488, t455, t437 * r_i_i_C(1) - t438 * r_i_i_C(2), 0, 0; t438 * r_i_i_C(1) + t437 * r_i_i_C(2) + t440 * pkin(3) - t463 * pkin(2) - t462 * t538 + t470 * qJD(3) + (-t495 * pkin(1) + t519 * t535) * qJD(1) + t547 * (-t439 * t487 + t455 * t491), -t464 * pkin(2) + t507 * (-t464 * t490 - t465 * t541) + t499 * (t464 * t486 - t465 * t534) + ((t458 * t531 - t459 * t493) * r_i_i_C(1) + (-t458 * t532 - t459 * t496) * r_i_i_C(2)) * qJD(4) + (t479 * qJD(3) + t513 * t465 + ((t465 * t493 + t479 * t522) * r_i_i_C(1) + (t465 * t496 - t479 * t523) * r_i_i_C(2)) * t487) * t488, t457, t555 * r_i_i_C(1) + t554 * r_i_i_C(2) + (t552 * r_i_i_C(1) + t553 * r_i_i_C(2)) * qJD(4), 0, 0; 0, ((t475 * t531 - t476 * t493) * r_i_i_C(1) + (-t475 * t532 - t476 * t496) * r_i_i_C(2)) * qJD(4) + (-pkin(2) * t525 + (t494 * qJD(3) + t513 * t524 + ((t493 * t524 + t494 * t522) * r_i_i_C(1) + (-t494 * t523 + t496 * t524) * r_i_i_C(2)) * t487) * t488 + (t507 * (-t486 * t529 - t533) - t499 * t504) * qJD(2)) * t489, t515, (t473 * t531 - t474 * t493 + t496 * t512) * r_i_i_C(1) + (-t473 * t532 - t474 * t496 - t493 * t512) * r_i_i_C(2) + ((-t466 * t532 - t467 * t496 - t477 * t540) * r_i_i_C(1) + (-t466 * t531 + t467 * t493 - t477 * t539) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:11:08
	% EndTime: 2019-10-10 11:11:11
	% DurationCPUTime: 3.78s
	% Computational Cost: add. (1958->200), mult. (6206->363), div. (0->0), fcn. (6786->16), ass. (0->140)
	t739 = sin(qJ(2));
	t809 = cos(pkin(6));
	t810 = sin(qJ(1));
	t772 = t809 * t810;
	t764 = t739 * t772;
	t781 = t810 * t739;
	t742 = cos(qJ(2));
	t743 = cos(qJ(1));
	t786 = t743 * t742;
	t709 = -qJD(1) * t764 - qJD(2) * t781 + (qJD(2) * t809 + qJD(1)) * t786;
	t730 = sin(pkin(14));
	t734 = cos(pkin(14));
	t777 = t743 * t809;
	t723 = t739 * t777 + t810 * t742;
	t749 = t743 * t739 + t742 * t772;
	t708 = t749 * qJD(1) + t723 * qJD(2);
	t732 = sin(pkin(7));
	t736 = cos(pkin(7));
	t733 = sin(pkin(6));
	t783 = t733 * t810;
	t773 = qJD(1) * t783;
	t751 = -t708 * t736 + t732 * t773;
	t675 = t709 * t734 + t751 * t730;
	t738 = sin(qJ(4));
	t741 = cos(qJ(4));
	t674 = -t709 * t730 + t751 * t734;
	t803 = t708 * t732;
	t698 = t736 * t773 + t803;
	t731 = sin(pkin(8));
	t735 = cos(pkin(8));
	t771 = t674 * t735 + t698 * t731;
	t722 = -t742 * t777 + t781;
	t791 = t733 * t743;
	t759 = t722 * t736 + t732 * t791;
	t693 = -t723 * t734 + t759 * t730;
	t692 = t723 * t730 + t759 * t734;
	t712 = -t722 * t732 + t736 * t791;
	t768 = t692 * t735 + t712 * t731;
	t822 = -t693 * t738 + t768 * t741;
	t649 = t822 * qJD(4) - t675 * t741 - t771 * t738;
	t658 = t674 * t731 - t698 * t735;
	t737 = sin(qJ(5));
	t740 = cos(qJ(5));
	t832 = t649 * t737 - t658 * t740;
	t831 = t649 * t740 + t658 * t737;
	t661 = t693 * t741 + t768 * t738;
	t678 = t692 * t731 - t712 * t735;
	t830 = -t661 * t737 - t678 * t740;
	t829 = t661 * t740 - t678 * t737;
	t828 = t661 * qJD(4) - t675 * t738 + t771 * t741;
	t780 = qJD(1) * t791;
	t750 = t764 - t786;
	t706 = t722 * qJD(1) + t750 * qJD(2);
	t804 = t706 * t732;
	t696 = t736 * t780 - t804;
	t707 = t723 * qJD(1) + t749 * qJD(2);
	t754 = t706 * t736 + t732 * t780;
	t745 = t707 * t730 + t754 * t734;
	t812 = t696 * t731 + t745 * t735;
	t776 = t809 * t732;
	t787 = t736 * t742;
	t789 = t734 * t739;
	t711 = t733 * t789 + (t733 * t787 + t776) * t730;
	t758 = -t730 * t739 + t734 * t787;
	t710 = t758 * t733 + t734 * t776;
	t792 = t732 * t742;
	t721 = -t733 * t792 + t809 * t736;
	t766 = t710 * t735 + t721 * t731;
	t672 = t711 * t741 + t766 * t738;
	t811 = r_i_i_C(3) + pkin(12);
	t790 = t734 * t736;
	t681 = -t706 * t730 + t707 * t790;
	t808 = t681 * t731;
	t683 = t708 * t730 - t709 * t790;
	t807 = t683 * t731;
	t785 = qJD(2) * t733;
	t715 = t758 * t785;
	t799 = t715 * t731;
	t796 = t730 * t736;
	t795 = t731 * t732;
	t794 = t732 * t735;
	t793 = t732 * t739;
	t788 = t736 * t739;
	t784 = t733 * t793;
	t782 = t736 * t810;
	t779 = t732 * t785;
	t778 = pkin(11) * t735 + qJ(3);
	t775 = t742 * t779;
	t774 = t739 * t779;
	t755 = t732 * t783 - t736 * t749;
	t694 = t730 * t750 + t755 * t734;
	t714 = t732 * t749 + t733 * t782;
	t767 = t694 * t735 + t714 * t731;
	t765 = t740 * r_i_i_C(1) - t737 * r_i_i_C(2) + pkin(4);
	t763 = -t681 * t735 + t707 * t795;
	t762 = t683 * t735 + t709 * t795;
	t701 = t722 * t730 - t723 * t790;
	t761 = t701 * t735 + t723 * t795;
	t703 = t730 * t749 + t750 * t790;
	t760 = t703 * t735 - t750 * t795;
	t757 = qJD(5) * (-t737 * r_i_i_C(1) - t740 * r_i_i_C(2));
	t719 = (-t730 * t742 - t734 * t788) * t733;
	t756 = t719 * t735 + t731 * t784;
	t720 = (-t730 * t788 + t734 * t742) * t733;
	t753 = -t715 * t735 + t731 * t775;
	t717 = qJD(2) * t719;
	t752 = t717 * t735 + t731 * t774;
	t695 = t755 * t730 - t734 * t750;
	t747 = -t695 * t738 + t767 * t741;
	t663 = t695 * t741 + t767 * t738;
	t746 = -t711 * t738 + t766 * t741;
	t702 = -t722 * t734 - t723 * t796;
	t668 = t702 * t741 + t761 * t738;
	t704 = -t734 * t749 + t750 * t796;
	t669 = t704 * t741 + t760 * t738;
	t685 = t720 * t741 + t756 * t738;
	t656 = t696 * t735 - t745 * t731;
	t718 = qJD(2) * t720;
	t716 = (-t730 * t787 - t789) * t785;
	t705 = -t719 * t731 + t735 * t784;
	t700 = -t717 * t731 + t735 * t774;
	t699 = t735 * t775 + t799;
	t689 = -t710 * t731 + t721 * t735;
	t687 = -t703 * t731 - t750 * t794;
	t686 = -t701 * t731 + t723 * t794;
	t684 = -t708 * t734 - t709 * t796;
	t682 = t706 * t734 + t707 * t796;
	t680 = -t694 * t731 + t714 * t735;
	t673 = -t707 * t734 + t754 * t730;
	t667 = t709 * t794 - t807;
	t666 = -t707 * t794 - t808;
	t665 = t716 * t741 + t753 * t738 + (-t720 * t738 + t756 * t741) * qJD(4);
	t655 = t746 * qJD(4) + t718 * t741 + t752 * t738;
	t653 = t684 * t741 + t762 * t738 + (-t702 * t738 + t761 * t741) * qJD(4);
	t651 = t682 * t741 - t763 * t738 + (-t704 * t738 + t760 * t741) * qJD(4);
	t645 = t747 * qJD(4) + t673 * t741 + t812 * t738;
	t644 = t663 * qJD(4) + t673 * t738 - t812 * t741;
	t643 = t645 * t740 + t656 * t737 + (-t663 * t737 + t680 * t740) * qJD(5);
	t642 = -t645 * t737 + t656 * t740 + (-t663 * t740 - t680 * t737) * qJD(5);
	t1 = [t831 * r_i_i_C(1) - t832 * r_i_i_C(2) + t649 * pkin(4) - t675 * pkin(3) - t709 * pkin(2) - qJ(3) * t803 + t811 * t828 + (t830 * r_i_i_C(1) - t829 * r_i_i_C(2)) * qJD(5) + t712 * qJD(3) + t658 * pkin(11) + (-t743 * pkin(1) + (-t810 * pkin(10) - qJ(3) * t782) * t733) * qJD(1), (t651 * t740 + t666 * t737) * r_i_i_C(1) + (-t651 * t737 + t666 * t740) * r_i_i_C(2) + t651 * pkin(4) + t682 * pkin(3) - pkin(11) * t808 + t706 * pkin(2) + t811 * (t669 * qJD(4) + t682 * t738 + t763 * t741) + ((-t669 * t737 + t687 * t740) * r_i_i_C(1) + (-t669 * t740 - t687 * t737) * r_i_i_C(2)) * qJD(5) + (-qJD(3) * t750 - t778 * t707) * t732, t696, -t765 * t644 + t811 * t645 + t747 * t757, t642 * r_i_i_C(1) - t643 * r_i_i_C(2), 0; -qJ(3) * t804 - t707 * pkin(2) + t673 * pkin(3) + t645 * pkin(4) + t643 * r_i_i_C(1) + t642 * r_i_i_C(2) + t811 * t644 + t714 * qJD(3) + (-t810 * pkin(1) + (qJ(3) * t736 + pkin(10)) * t791) * qJD(1) + t656 * pkin(11), (t653 * t740 + t667 * t737) * r_i_i_C(1) + (-t653 * t737 + t667 * t740) * r_i_i_C(2) + t653 * pkin(4) + t684 * pkin(3) - pkin(11) * t807 - t708 * pkin(2) + t811 * (t668 * qJD(4) + t684 * t738 - t762 * t741) + ((-t668 * t737 + t686 * t740) * r_i_i_C(1) + (-t668 * t740 - t686 * t737) * r_i_i_C(2)) * qJD(5) + (t723 * qJD(3) + t778 * t709) * t732, t698, -t649 * t811 - t822 * t757 + t765 * t828, t832 * r_i_i_C(1) + t831 * r_i_i_C(2) + (t829 * r_i_i_C(1) + t830 * r_i_i_C(2)) * qJD(5), 0; 0, (t665 * t740 + t699 * t737) * r_i_i_C(1) + (-t665 * t737 + t699 * t740) * r_i_i_C(2) + t665 * pkin(4) + t716 * pkin(3) + pkin(11) * t799 + t811 * (t685 * qJD(4) + t716 * t738 - t753 * t741) + ((-t685 * t737 + t705 * t740) * r_i_i_C(1) + (-t685 * t740 - t705 * t737) * r_i_i_C(2)) * qJD(5) + (qJD(3) * t793 + (-pkin(2) * t739 + t778 * t792) * qJD(2)) * t733, t774, t811 * t655 + t746 * t757 + t765 * (-t672 * qJD(4) - t718 * t738 + t752 * t741), (-t655 * t737 + t700 * t740) * r_i_i_C(1) + (-t655 * t740 - t700 * t737) * r_i_i_C(2) + ((-t672 * t740 - t689 * t737) * r_i_i_C(1) + (t672 * t737 - t689 * t740) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:11:16
	% EndTime: 2019-10-10 11:11:22
	% DurationCPUTime: 6.35s
	% Computational Cost: add. (5091->293), mult. (15916->498), div. (0->0), fcn. (18069->18), ass. (0->185)
	t1061 = cos(qJ(4));
	t1058 = cos(pkin(6));
	t1059 = sin(qJ(2));
	t1008 = t1058 * t1059;
	t950 = cos(qJ(1));
	t1000 = t950 * t1008;
	t1060 = sin(qJ(1));
	t949 = cos(qJ(2));
	t927 = t1060 * t949 + t1000;
	t937 = sin(pkin(14));
	t941 = cos(pkin(14));
	t940 = sin(pkin(6));
	t1039 = t940 * t950;
	t1010 = t1060 * t1059;
	t1017 = t949 * t1058;
	t926 = -t950 * t1017 + t1010;
	t939 = sin(pkin(7));
	t943 = cos(pkin(7));
	t998 = t939 * t1039 + t926 * t943;
	t1064 = t927 * t937 + t998 * t941;
	t918 = t943 * t1039 - t926 * t939;
	t938 = sin(pkin(8));
	t942 = cos(pkin(8));
	t1075 = t1064 * t942 + t918 * t938;
	t900 = -t927 * t941 + t998 * t937;
	t946 = sin(qJ(4));
	t862 = t900 * t1061 + t1075 * t946;
	t882 = t1064 * t938 - t918 * t942;
	t945 = sin(qJ(5));
	t948 = cos(qJ(5));
	t847 = t862 * t948 - t882 * t945;
	t944 = sin(qJ(6));
	t1088 = t847 * t944;
	t947 = cos(qJ(6));
	t1087 = t847 * t947;
	t1036 = t950 * t949;
	t996 = t1060 * t1008;
	t913 = -qJD(1) * t996 - qJD(2) * t1010 + (qJD(2) * t1058 + qJD(1)) * t1036;
	t1020 = qJD(1) * t1060;
	t1009 = t940 * t1020;
	t983 = t1060 * t1017 + t950 * t1059;
	t912 = t983 * qJD(1) + t927 * qJD(2);
	t987 = t939 * t1009 - t912 * t943;
	t1063 = t913 * t937 - t987 * t941;
	t1052 = t912 * t939;
	t988 = t943 * t1009 + t1052;
	t857 = t1063 * t938 + t988 * t942;
	t1086 = t847 * qJD(5) + t857 * t948;
	t1004 = t862 * t945 + t882 * t948;
	t1085 = t1004 * qJD(5) + t857 * t945;
	t879 = t913 * t941 + t987 * t937;
	t1083 = -t862 * qJD(4) + t879 * t946;
	t1077 = t900 * t946;
	t1016 = t1058 * t939;
	t1023 = t1059 * t937;
	t1037 = t943 * t949;
	t993 = t941 * t1037 - t1023;
	t967 = t941 * t1016 + t993 * t940;
	t1040 = t939 * t949;
	t1029 = t940 * t1040;
	t990 = t1058 * t943 - t1029;
	t1070 = t990 * t938 + t967 * t942;
	t982 = t996 - t1036;
	t992 = t1060 * t939 * t940 - t943 * t983;
	t974 = -t937 * t982 - t992 * t941;
	t1025 = t943 * t1060;
	t1013 = t940 * t1025;
	t994 = t939 * t983 + t1013;
	t1066 = t994 * t938 - t974 * t942;
	t1068 = t1063 * t942 - t988 * t938;
	t1074 = -t879 * t1061 + t1068 * t946;
	t1062 = r_i_i_C(3) + pkin(13);
	t911 = qJD(1) * t1000 + t983 * qJD(2) + t949 * t1020;
	t976 = t982 * qJD(2);
	t957 = t998 * qJD(1) + t943 * t976;
	t954 = t911 * t937 + t957 * t941;
	t958 = t918 * qJD(1) - t939 * t976;
	t951 = -t954 * t938 + t958 * t942;
	t1072 = t938 * t958 + t942 * t954;
	t1071 = t1070 * t1061;
	t1069 = t1066 * t1061;
	t1067 = t1075 * t1061;
	t997 = qJD(6) * (t944 * r_i_i_C(1) + t947 * r_i_i_C(2));
	t1041 = t939 * t942;
	t1038 = t941 * t943;
	t959 = t926 * qJD(1) + t976;
	t885 = t911 * t1038 - t959 * t937;
	t867 = -t911 * t1041 - t885 * t938;
	t901 = t992 * t937 - t941 * t982;
	t864 = t901 * t1061 + t1066 * t946;
	t884 = t974 * t938 + t994 * t942;
	t1065 = -t864 * t945 + t884 * t948;
	t1057 = qJ(3) * t939;
	t887 = -t913 * t1038 + t912 * t937;
	t1053 = t887 * t938;
	t1035 = qJD(2) * t940;
	t919 = t993 * t1035;
	t1047 = t919 * t938;
	t1022 = t1059 * t941;
	t991 = t943 * t1022 + t937 * t949;
	t923 = t991 * t940;
	t1046 = t923 * t942;
	t1043 = t937 * t943;
	t1042 = t939 * t938;
	t1034 = qJD(3) * t939;
	t1033 = qJD(4) * t946;
	t1032 = qJD(6) * t944;
	t1031 = qJD(6) * t947;
	t1028 = t938 * t1061;
	t1027 = t939 * t1059;
	t1026 = t942 * t1061;
	t1019 = qJD(4) * t1061;
	t1018 = pkin(11) * t942 + qJ(3);
	t1015 = t939 * t1028;
	t1014 = t940 * t1027;
	t1012 = t942 * t1027;
	t1011 = qJD(2) * t1029;
	t1006 = t938 * t1014;
	t1005 = t938 * t1011;
	t849 = t864 * t948 + t884 * t945;
	t905 = -t927 * t1038 + t926 * t937;
	t906 = -t927 * t1043 - t926 * t941;
	t870 = t906 * t1061 + (t927 * t1042 + t905 * t942) * t946;
	t891 = t927 * t1041 - t905 * t938;
	t850 = t870 * t948 + t891 * t945;
	t907 = t1038 * t982 + t937 * t983;
	t908 = t1043 * t982 - t941 * t983;
	t872 = t908 * t1061 + (-t1042 * t982 + t907 * t942) * t946;
	t892 = -t1041 * t982 - t907 * t938;
	t851 = t872 * t948 + t892 * t945;
	t916 = t940 * t1022 + (t940 * t1037 + t1016) * t937;
	t877 = t916 * t1061 + t1070 * t946;
	t897 = -t967 * t938 + t990 * t942;
	t855 = t877 * t948 + t897 * t945;
	t1003 = -t877 * t945 + t897 * t948;
	t924 = (-t943 * t1023 + t941 * t949) * t940;
	t890 = t924 * t1061 + (t1006 - t1046) * t946;
	t909 = t940 * t1012 + t923 * t938;
	t873 = t890 * t948 + t909 * t945;
	t1002 = qJD(2) * t1014;
	t1001 = t947 * r_i_i_C(1) - t944 * r_i_i_C(2) + pkin(5);
	t995 = t1061 * t1006;
	t978 = qJD(2) * t1046;
	t975 = -t1001 * t948 - t1062 * t945 - pkin(4);
	t973 = t1026 * t1064 + t918 * t1028 - t1077;
	t971 = t927 * t1015 + t905 * t1026 - t906 * t946;
	t970 = -t1015 * t982 + t907 * t1026 - t908 * t946;
	t964 = -t923 * t1026 - t924 * t946 + t995;
	t960 = t948 * t997 + (t1001 * t945 - t1062 * t948) * qJD(5);
	t921 = qJD(2) * t924;
	t920 = (-t937 * t1037 - t1022) * t1035;
	t904 = (t991 * t938 + t1012) * t1035;
	t903 = t942 * t1011 + t1047;
	t888 = -t913 * t1043 - t912 * t941;
	t886 = t911 * t1043 + t959 * t941;
	t878 = -t911 * t941 + t957 * t937;
	t876 = t916 * t946 - t1071;
	t868 = t913 * t1041 - t1053;
	t866 = t920 * t1061 + (-t919 * t942 + t1005) * t946 + t964 * qJD(4);
	t865 = t890 * qJD(4) - t1061 * t1005 + t919 * t1026 + t920 * t946;
	t863 = t901 * t946 - t1069;
	t859 = t1067 - t1077;
	t853 = -t916 * t1033 + t921 * t1061 + (t1002 * t938 - t978) * t946 + t1071 * qJD(4);
	t852 = -qJD(2) * t995 + t916 * t1019 + t1070 * t1033 + t1061 * t978 + t921 * t946;
	t844 = t866 * t948 + t903 * t945 + (-t890 * t945 + t909 * t948) * qJD(5);
	t842 = t888 * t1061 + (t913 * t1042 + t887 * t942) * t946 + t971 * qJD(4);
	t841 = t870 * qJD(4) - t913 * t1015 - t887 * t1026 + t888 * t946;
	t840 = t886 * t1061 + (-t911 * t1042 + t885 * t942) * t946 + t970 * qJD(4);
	t839 = t872 * qJD(4) + t911 * t1015 - t885 * t1026 + t886 * t946;
	t838 = t1003 * qJD(5) + t853 * t948 + t904 * t945;
	t836 = t973 * qJD(4) + t1074;
	t835 = -t1026 * t1063 + t1028 * t988 - t1083;
	t834 = -t1067 * qJD(4) + t1033 * t900 - t1074;
	t833 = t1068 * t1061 + t1083;
	t832 = t1069 * qJD(4) - t901 * t1033 + t878 * t1061 + t1072 * t946;
	t831 = t901 * t1019 + t1066 * t1033 - t1072 * t1061 + t878 * t946;
	t830 = t842 * t948 + t868 * t945 + (-t870 * t945 + t891 * t948) * qJD(5);
	t828 = t840 * t948 + t867 * t945 + (-t872 * t945 + t892 * t948) * qJD(5);
	t826 = t836 * t948 - t1085;
	t824 = t834 * t948 + t1085;
	t822 = t1065 * qJD(5) + t832 * t948 + t951 * t945;
	t821 = t849 * qJD(5) + t832 * t945 - t951 * t948;
	t820 = t822 * t947 + t831 * t944 + (-t849 * t944 + t863 * t947) * qJD(6);
	t819 = -t822 * t944 + t831 * t947 + (-t849 * t947 - t863 * t944) * qJD(6);
	t1 = [(t826 * t947 + t835 * t944) * r_i_i_C(1) + (-t826 * t944 + t835 * t947) * r_i_i_C(2) + t826 * pkin(5) + t836 * pkin(4) + t835 * pkin(12) - t879 * pkin(3) - t913 * pkin(2) - qJ(3) * t1052 + t1062 * (t836 * t945 + t1086) + ((-t947 * t973 - t1088) * r_i_i_C(1) + (t944 * t973 - t1087) * r_i_i_C(2)) * qJD(6) + t918 * qJD(3) - t857 * pkin(11) + (-t950 * pkin(1) + (-t1060 * pkin(10) - qJ(3) * t1025) * t940) * qJD(1), (t828 * t947 + t839 * t944 + (-t851 * t944 - t947 * t970) * qJD(6)) * r_i_i_C(1) + (-t828 * t944 + t839 * t947 + (-t851 * t947 + t944 * t970) * qJD(6)) * r_i_i_C(2) + t828 * pkin(5) + t840 * pkin(4) + t839 * pkin(12) + t886 * pkin(3) + t959 * pkin(2) - t911 * t1057 - t982 * t1034 + t1062 * (t851 * qJD(5) + t840 * t945 - t867 * t948) + t867 * pkin(11), t958, (t864 * t1031 + t832 * t944) * r_i_i_C(1) + (-t864 * t1032 + t832 * t947) * r_i_i_C(2) + t832 * pkin(12) + t975 * t831 + t960 * t863, -t1001 * t821 + t1062 * t822 - t1065 * t997, t819 * r_i_i_C(1) - t820 * r_i_i_C(2); -pkin(1) * t1020 - t911 * pkin(2) + t878 * pkin(3) + t832 * pkin(4) + t822 * pkin(5) + t831 * pkin(12) + t820 * r_i_i_C(1) + t819 * r_i_i_C(2) + qJD(3) * t1013 + t983 * t1034 - t959 * t1057 + t1062 * t821 + (qJ(3) * t943 + pkin(10)) * qJD(1) * t1039 + t951 * pkin(11), (t830 * t947 + t841 * t944) * r_i_i_C(1) + (-t830 * t944 + t841 * t947) * r_i_i_C(2) + t830 * pkin(5) + t842 * pkin(4) + t841 * pkin(12) + t888 * pkin(3) - pkin(11) * t1053 - t912 * pkin(2) + t1062 * (t850 * qJD(5) + t842 * t945 - t868 * t948) + ((-t850 * t944 - t947 * t971) * r_i_i_C(1) + (-t850 * t947 + t944 * t971) * r_i_i_C(2)) * qJD(6) + (t927 * qJD(3) + t1018 * t913) * t939, t988, (-t1031 * t862 + t834 * t944) * r_i_i_C(1) + (t1032 * t862 + t834 * t947) * r_i_i_C(2) + t834 * pkin(12) + t975 * t833 + t960 * t859, t1062 * t824 - t1004 * t997 + t1001 * (-t834 * t945 + t1086), (-t824 * t944 + t833 * t947) * r_i_i_C(1) + (-t824 * t947 - t833 * t944) * r_i_i_C(2) + ((-t859 * t944 + t1087) * r_i_i_C(1) + (-t859 * t947 - t1088) * r_i_i_C(2)) * qJD(6); 0, (t844 * t947 + t865 * t944) * r_i_i_C(1) + (-t844 * t944 + t865 * t947) * r_i_i_C(2) + t844 * pkin(5) + t866 * pkin(4) + t865 * pkin(12) + t920 * pkin(3) + pkin(11) * t1047 + t1062 * (t873 * qJD(5) + t866 * t945 - t903 * t948) + ((-t873 * t944 - t947 * t964) * r_i_i_C(1) + (-t873 * t947 + t944 * t964) * r_i_i_C(2)) * qJD(6) + (qJD(3) * t1027 + (-t1059 * pkin(2) + t1018 * t1040) * qJD(2)) * t940, t1002, (t877 * t1031 + t853 * t944) * r_i_i_C(1) + (-t877 * t1032 + t853 * t947) * r_i_i_C(2) + t853 * pkin(12) + t975 * t852 + t960 * t876, t1062 * t838 - t1003 * t997 + t1001 * (-t855 * qJD(5) - t853 * t945 + t904 * t948), (-t838 * t944 + t852 * t947) * r_i_i_C(1) + (-t838 * t947 - t852 * t944) * r_i_i_C(2) + ((-t855 * t947 - t876 * t944) * r_i_i_C(1) + (t855 * t944 - t876 * t947) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end