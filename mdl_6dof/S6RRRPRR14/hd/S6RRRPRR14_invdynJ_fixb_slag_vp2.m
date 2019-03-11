% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR14_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR14_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR14_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR14_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR14_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:12:25
% EndTime: 2019-03-09 20:14:02
% DurationCPUTime: 60.64s
% Computational Cost: add. (22367->1113), mult. (52204->1464), div. (0->0), fcn. (41065->14), ass. (0->508)
t401 = cos(qJ(2));
t392 = sin(pkin(6));
t541 = qJD(1) * t392;
t513 = t401 * t541;
t351 = -qJD(3) + t513;
t396 = sin(qJ(2));
t573 = cos(pkin(6));
t490 = t573 * qJD(1);
t475 = pkin(1) * t490;
t301 = pkin(8) * t513 + t396 * t475;
t395 = sin(qJ(3));
t478 = t395 * t513;
t536 = qJD(3) * t395;
t782 = -qJD(4) * t395 - t301 + (-t478 + t536) * pkin(3);
t540 = qJD(1) * t396;
t514 = t392 * t540;
t298 = -pkin(8) * t514 + t401 * t475;
t432 = (pkin(2) * t396 - pkin(9) * t401) * t392;
t299 = qJD(1) * t432;
t400 = cos(qJ(3));
t209 = -t395 * t298 + t299 * t400;
t535 = qJD(3) * t400;
t547 = t400 * t401;
t650 = pkin(4) + pkin(9);
t651 = pkin(3) + pkin(10);
t781 = -(pkin(4) * t547 - t396 * t651) * t541 + t209 + t650 * t535;
t780 = -t782 + t351 * (pkin(10) * t395 - qJ(4) * t400);
t394 = sin(qJ(5));
t399 = cos(qJ(5));
t548 = t399 * t401;
t431 = -t394 * t396 + t395 * t548;
t262 = t431 * t541;
t532 = qJD(5) * t400;
t509 = t394 * t532;
t779 = t262 - t509;
t711 = -t399 * t536 + t779;
t778 = Ifges(4,4) + Ifges(5,6);
t777 = t780 * t394 + t399 * t781;
t571 = qJ(4) * t395;
t496 = -pkin(2) - t571;
t326 = -t400 * t651 + t496;
t354 = t650 * t395;
t533 = qJD(5) * t399;
t534 = qJD(5) * t394;
t721 = -t326 * t534 + t354 * t533 + t394 * t781 - t780 * t399;
t393 = sin(qJ(6));
t530 = qJD(6) * t393;
t398 = cos(qJ(6));
t551 = t398 * t399;
t692 = qJD(5) + qJD(6);
t229 = -t393 * t534 - t394 * t530 + t551 * t692;
t438 = t490 + qJD(2);
t512 = t400 * t540;
t273 = t392 * t512 + t395 * t438;
t440 = t393 * t394 - t551;
t720 = t440 * t273;
t714 = t229 - t720;
t529 = qJD(1) * qJD(2);
t425 = qJDD(1) * t396 + t401 * t529;
t413 = t425 * t392;
t421 = qJD(3) * t438;
t486 = t573 * qJDD(1);
t435 = t486 + qJDD(2);
t560 = t392 * t396;
t523 = t395 * t560;
t476 = qJD(3) * t523;
t179 = qJD(1) * t476 - t395 * t435 + (-t413 - t421) * t400;
t640 = -t179 / 0.2e1;
t180 = t395 * t421 - t400 * t435 + (qJD(3) * t512 + t395 * t425) * t392;
t638 = -t180 / 0.2e1;
t307 = (-qJDD(1) * t401 + t396 * t529) * t392;
t294 = qJDD(3) + t307;
t623 = t294 / 0.2e1;
t736 = Ifges(4,5) - Ifges(5,4);
t735 = Ifges(4,6) - Ifges(5,5);
t734 = Ifges(4,3) + Ifges(5,1);
t552 = t395 * t401;
t430 = t394 * t552 + t396 * t399;
t263 = t430 * t541;
t332 = t394 * t354;
t477 = t400 * t513;
t493 = pkin(11) * t400 - t326;
t555 = t394 * t395;
t776 = -pkin(5) * t477 + pkin(11) * t263 + (pkin(5) * t400 - pkin(11) * t555) * qJD(3) + (t399 * t493 - t332) * qJD(5) + t777;
t775 = -pkin(11) * t711 + t721;
t256 = pkin(9) * t438 + t301;
t264 = (-pkin(2) * t401 - pkin(9) * t396 - pkin(1)) * t541;
t174 = t400 * t256 + t395 * t264;
t272 = t395 * t514 - t400 * t438;
t129 = -pkin(4) * t272 + t174;
t123 = t399 * t129;
t572 = qJ(4) * t272;
t148 = t273 * t651 + t572;
t604 = pkin(11) + t651;
t774 = t604 * t534 + pkin(5) * t272 - t123 - (-pkin(11) * t273 - t148) * t394;
t343 = t604 * t399;
t567 = t273 * t399;
t72 = t394 * t129 + t399 * t148;
t773 = pkin(11) * t567 + qJD(5) * t343 + t72;
t441 = t393 * t399 + t398 * t394;
t230 = t692 * t441;
t426 = t441 * t273;
t713 = -t230 - t426;
t267 = qJD(5) + t273;
t216 = t272 * t394 - t351 * t399;
t173 = t256 * t395 - t400 * t264;
t433 = pkin(4) * t273 + t173;
t758 = qJD(4) + t433;
t104 = t351 * t651 + t758;
t255 = -pkin(2) * t438 - t298;
t405 = -t273 * qJ(4) + t255;
t109 = t272 * t651 + t405;
t53 = t399 * t104 - t109 * t394;
t46 = -pkin(11) * t216 + t53;
t40 = pkin(5) * t267 + t46;
t215 = t272 * t399 + t351 * t394;
t54 = t104 * t394 + t109 * t399;
t47 = pkin(11) * t215 + t54;
t581 = t393 * t47;
t15 = t398 * t40 - t581;
t577 = t398 * t47;
t16 = t393 * t40 + t577;
t761 = -pkin(8) * t392 * t529 + pkin(1) * t486;
t528 = qJDD(1) * t392;
t763 = pkin(8) * t528 + qJD(2) * t475;
t217 = t396 * t761 + t401 * t763;
t201 = pkin(9) * t435 + t217;
t208 = t307 * pkin(2) + (-qJDD(1) * pkin(1) - pkin(9) * t425) * t392;
t70 = -t395 * t201 + t208 * t400 - t256 * t535 - t264 * t536;
t422 = qJDD(4) - t70;
t43 = -pkin(4) * t179 - t294 * t651 + t422;
t218 = -t396 * t763 + t401 * t761;
t202 = -pkin(2) * t435 - t218;
t404 = t179 * qJ(4) - t273 * qJD(4) + t202;
t49 = t180 * t651 + t404;
t12 = -qJD(5) * t54 - t394 * t49 + t399 * t43;
t172 = qJDD(5) - t179;
t91 = qJD(5) * t215 + t180 * t394 + t294 * t399;
t6 = pkin(5) * t172 - pkin(11) * t91 + t12;
t11 = t104 * t533 - t109 * t534 + t394 * t43 + t399 * t49;
t92 = -qJD(5) * t216 + t180 * t399 - t294 * t394;
t7 = pkin(11) * t92 + t11;
t2 = qJD(6) * t15 + t393 * t6 + t398 * t7;
t3 = -qJD(6) * t16 - t393 * t7 + t398 * t6;
t397 = sin(qJ(1));
t618 = cos(qJ(1));
t467 = t573 * t618;
t318 = t396 * t467 + t397 * t401;
t515 = t392 * t618;
t242 = t318 * t395 + t400 * t515;
t492 = t396 * t573;
t320 = -t397 * t492 + t401 * t618;
t558 = t392 * t400;
t246 = t320 * t395 - t397 * t558;
t315 = -t400 * t573 + t523;
t676 = g(1) * t246 + g(2) * t242 + g(3) * t315;
t772 = -t15 * t713 - t16 * t714 - t2 * t441 + t3 * t440 + t676;
t143 = t272 * pkin(3) + t405;
t338 = t351 * qJ(4);
t149 = t338 - t174;
t696 = t149 * mrSges(5,1) - t174 * mrSges(4,3);
t771 = t255 * mrSges(4,1) - t143 * mrSges(5,2) + t696;
t770 = t255 * mrSges(4,2) - t143 * mrSges(5,3);
t694 = -qJD(4) - t173;
t147 = pkin(3) * t351 - t694;
t769 = t147 * mrSges(5,1) + t53 * mrSges(6,1) + t15 * mrSges(7,1) - t54 * mrSges(6,2) - t16 * mrSges(7,2) + t173 * mrSges(4,3);
t35 = Ifges(6,5) * t91 + Ifges(6,6) * t92 + Ifges(6,3) * t172;
t639 = t179 / 0.2e1;
t641 = t172 / 0.2e1;
t161 = qJDD(6) + t172;
t644 = t161 / 0.2e1;
t652 = t92 / 0.2e1;
t653 = t91 / 0.2e1;
t121 = t215 * t393 + t216 * t398;
t32 = -qJD(6) * t121 - t393 * t91 + t398 * t92;
t661 = t32 / 0.2e1;
t488 = t398 * t215 - t216 * t393;
t31 = qJD(6) * t488 + t393 * t92 + t398 * t91;
t662 = t31 / 0.2e1;
t682 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t684 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t766 = -t294 / 0.2e1;
t8 = Ifges(7,5) * t31 + Ifges(7,6) * t32 + Ifges(7,3) * t161;
t768 = t682 + t684 + 0.2e1 * Ifges(4,1) * t640 + Ifges(5,4) * t766 + Ifges(6,5) * t653 + Ifges(7,5) * t662 + Ifges(6,6) * t652 + Ifges(7,6) * t661 + Ifges(6,3) * t641 + Ifges(7,3) * t644 + (-t639 + t640) * Ifges(5,2) + t35 / 0.2e1 + t8 / 0.2e1 + t778 * t638 + (t736 + Ifges(4,5)) * t623;
t594 = Ifges(7,4) * t121;
t258 = qJD(6) + t267;
t631 = -t258 / 0.2e1;
t646 = -t121 / 0.2e1;
t110 = t129 - t338;
t84 = -pkin(5) * t215 + t110;
t767 = t684 + t8 + (Ifges(7,5) * t488 - Ifges(7,6) * t121) * t631 + (t121 * t16 + t15 * t488) * mrSges(7,3) - t84 * (mrSges(7,1) * t121 + mrSges(7,2) * t488) + (Ifges(7,1) * t488 - t594) * t646;
t637 = t180 / 0.2e1;
t765 = mrSges(6,2) * t399;
t458 = mrSges(5,2) * t400 - mrSges(5,3) * t395;
t391 = qJ(5) + qJ(6);
t386 = sin(t391);
t387 = cos(t391);
t460 = t386 * mrSges(7,1) + t387 * mrSges(7,2);
t465 = mrSges(4,1) * t400 - mrSges(4,2) * t395;
t614 = pkin(5) * t394;
t402 = -pkin(11) - pkin(10);
t744 = m(7) * t402;
t764 = -t395 * (m(7) * t614 + t460) + t400 * (-mrSges(7,3) + t744) + t458 - t465;
t462 = mrSges(6,1) * t394 + t765;
t762 = -t460 - t462;
t749 = m(6) * pkin(10);
t679 = mrSges(4,1) - mrSges(5,2) - t744 + t749;
t693 = m(6) + m(7) + m(5);
t759 = pkin(3) * t693 + t679;
t384 = pkin(5) * t399 + pkin(4);
t605 = pkin(9) + t384;
t738 = mrSges(5,1) + mrSges(4,3);
t675 = m(6) * t650 + m(7) * t605 - mrSges(3,2) + t738;
t117 = Ifges(7,4) * t488;
t756 = -Ifges(7,2) * t121 + t117;
t754 = Ifges(4,6) * t766 + 0.2e1 * Ifges(5,3) * t637 + t778 * t639 + (t637 - t638) * Ifges(4,2) + (-t735 + Ifges(5,5)) * t623;
t603 = mrSges(6,3) * t215;
t150 = -mrSges(6,2) * t267 + t603;
t602 = mrSges(6,3) * t216;
t151 = mrSges(6,1) * t267 - t602;
t443 = t150 * t399 - t151 * t394;
t63 = mrSges(6,1) * t172 - mrSges(6,3) * t91;
t64 = -mrSges(6,2) * t172 + mrSges(6,3) * t92;
t753 = t443 * qJD(5) + t394 * t64 + t399 * t63;
t265 = Ifges(5,6) * t272;
t164 = -t351 * Ifges(5,4) - t273 * Ifges(5,2) + t265;
t630 = t258 / 0.2e1;
t645 = t121 / 0.2e1;
t647 = t488 / 0.2e1;
t266 = Ifges(4,4) * t272;
t689 = t273 * Ifges(4,1) - t351 * Ifges(4,5) + t216 * Ifges(6,5) + t121 * Ifges(7,5) + t215 * Ifges(6,6) + Ifges(7,6) * t488 + t267 * Ifges(6,3) + t258 * Ifges(7,3) - t266;
t751 = Ifges(7,5) * t645 + Ifges(7,6) * t647 + Ifges(7,3) * t630 - t164 / 0.2e1 + t689 / 0.2e1 + t769;
t748 = -t307 / 0.2e1;
t747 = t413 / 0.2e1;
t746 = t435 / 0.2e1;
t381 = qJ(4) + t614;
t745 = m(7) * t381;
t737 = mrSges(4,2) - mrSges(5,3);
t333 = t399 * t354;
t211 = pkin(5) * t395 + t394 * t493 + t333;
t228 = t399 * t326 + t332;
t549 = t399 * t400;
t219 = -pkin(11) * t549 + t228;
t115 = t211 * t398 - t219 * t393;
t733 = qJD(6) * t115 + t393 * t776 + t398 * t775;
t116 = t211 * t393 + t219 * t398;
t732 = -qJD(6) * t116 - t393 * t775 + t398 * t776;
t342 = t604 * t394;
t241 = -t342 * t398 - t343 * t393;
t725 = -qJD(6) * t241 + t393 * t773 + t398 * t774;
t240 = t342 * t393 - t343 * t398;
t724 = qJD(6) * t240 + t393 * t774 - t398 * t773;
t664 = m(7) * pkin(5);
t723 = -t664 - mrSges(6,1);
t722 = -qJD(5) * t228 + t777;
t557 = t392 * t401;
t325 = pkin(1) * t492 + pkin(8) * t557;
t290 = pkin(9) * t573 + t325;
t543 = pkin(2) * t557 + pkin(9) * t560;
t617 = pkin(1) * t392;
t291 = -t543 - t617;
t206 = -t395 * t290 + t291 * t400;
t375 = pkin(3) * t557;
t192 = -t206 + t375;
t316 = t395 * t573 + t396 * t558;
t134 = pkin(4) * t316 + pkin(10) * t557 + t192;
t491 = t401 * t573;
t324 = pkin(1) * t491 - pkin(8) * t560;
t289 = -pkin(2) * t573 - t324;
t308 = t315 * pkin(3);
t489 = t316 * qJ(4) - t308;
t190 = t289 - t489;
t607 = t315 * pkin(10);
t144 = t190 + t607;
t74 = t394 * t134 + t399 * t144;
t686 = t400 * t692;
t158 = t440 * t686 + t441 * t536;
t186 = t262 * t393 + t263 * t398;
t719 = t158 - t186;
t159 = -t440 * t536 + t441 * t686;
t185 = t262 * t398 - t263 * t393;
t718 = t159 - t185;
t597 = Ifges(6,4) * t216;
t102 = t215 * Ifges(6,2) + t267 * Ifges(6,6) + t597;
t582 = t273 * Ifges(5,6);
t162 = -t351 * Ifges(5,5) + t272 * Ifges(5,3) - t582;
t717 = t399 * t102 + t162;
t716 = -t272 * t735 + t273 * t736 - t351 * t734;
t127 = -mrSges(6,1) * t215 + mrSges(6,2) * t216;
t222 = mrSges(5,1) * t272 + mrSges(5,3) * t351;
t715 = -t222 + t127;
t210 = t400 * t298 + t395 * t299;
t193 = -qJ(4) * t514 - t210;
t156 = -pkin(4) * t478 - t193;
t712 = pkin(5) * t779 - t536 * t605 - t156;
t710 = -t394 * t536 + t399 * t532 + t263;
t317 = t396 * t397 - t401 * t467;
t562 = t317 * t400;
t709 = -pkin(3) * t562 - t317 * t571;
t319 = t396 * t618 + t397 * t491;
t561 = t319 * t400;
t708 = -pkin(3) * t561 - t319 * t571;
t707 = (t477 - t535) * qJ(4) + t782;
t706 = -t650 * t536 - t156;
t705 = pkin(5) * t533 + t273 * t384 - t694;
t704 = -m(6) * qJ(4) - t745 + t762;
t703 = -t143 * (-mrSges(5,2) * t395 - mrSges(5,3) * t400) - t255 * (mrSges(4,1) * t395 + mrSges(4,2) * t400);
t702 = -t395 * t735 + t400 * t736;
t700 = -t179 * t736 - t180 * t735 + t294 * t734;
t69 = t400 * t201 + t395 * t208 - t256 * t536 + t264 * t535;
t698 = -t395 * t70 + t400 * t69;
t60 = -t294 * qJ(4) + t351 * qJD(4) - t69;
t66 = -pkin(3) * t294 + t422;
t697 = t395 * t66 - t400 * t60;
t695 = -t11 * t394 - t12 * t399;
t601 = Ifges(3,4) * t396;
t667 = t392 ^ 2;
t690 = (pkin(1) * (mrSges(3,1) * t396 + mrSges(3,2) * t401) - t396 * (Ifges(3,1) * t401 - t601) / 0.2e1) * t667;
t463 = t399 * mrSges(6,1) - t394 * mrSges(6,2);
t214 = Ifges(6,4) * t215;
t103 = t216 * Ifges(6,1) + t267 * Ifges(6,5) + t214;
t556 = t394 * t103;
t687 = t110 * t463 - t556 / 0.2e1;
t685 = t737 - t745;
t483 = mrSges(3,3) * t514;
t683 = -m(4) * t255 + mrSges(3,1) * t438 - mrSges(4,1) * t272 - mrSges(4,2) * t273 - t483;
t678 = -m(5) * qJ(4) + t704 + t737;
t673 = -mrSges(6,3) - mrSges(7,3) - t679;
t672 = mrSges(6,1) * t555 + t395 * t765 + mrSges(3,1) - t764;
t671 = -t70 * mrSges(4,1) + t69 * mrSges(4,2) - t66 * mrSges(5,2) + t60 * mrSges(5,3);
t461 = t387 * mrSges(7,1) - t386 * mrSges(7,2);
t670 = -t461 - t463 - t675;
t665 = Ifges(7,4) * t662 + Ifges(7,2) * t661 + Ifges(7,6) * t644;
t663 = Ifges(7,1) * t662 + Ifges(7,4) * t661 + Ifges(7,5) * t644;
t36 = Ifges(6,4) * t91 + Ifges(6,2) * t92 + Ifges(6,6) * t172;
t660 = -t36 / 0.2e1;
t37 = Ifges(6,1) * t91 + Ifges(6,4) * t92 + Ifges(6,5) * t172;
t659 = -t37 / 0.2e1;
t58 = Ifges(7,2) * t488 + t258 * Ifges(7,6) + t594;
t658 = -t58 / 0.2e1;
t657 = t58 / 0.2e1;
t59 = t121 * Ifges(7,1) + t258 * Ifges(7,5) + t117;
t656 = -t59 / 0.2e1;
t655 = t59 / 0.2e1;
t649 = t102 / 0.2e1;
t648 = -t488 / 0.2e1;
t583 = t273 * Ifges(4,4);
t165 = -t272 * Ifges(4,2) - t351 * Ifges(4,6) + t583;
t642 = -t165 / 0.2e1;
t636 = -t215 / 0.2e1;
t635 = t215 / 0.2e1;
t634 = -t216 / 0.2e1;
t633 = t216 / 0.2e1;
t629 = -t267 / 0.2e1;
t628 = t267 / 0.2e1;
t627 = -t272 / 0.2e1;
t626 = t272 / 0.2e1;
t625 = -t273 / 0.2e1;
t624 = t273 / 0.2e1;
t621 = -t351 / 0.2e1;
t620 = t351 / 0.2e1;
t616 = pkin(5) * t216;
t238 = t315 * t399 + t394 * t557;
t615 = pkin(5) * t238;
t612 = pkin(9) * t317;
t611 = pkin(9) * t319;
t388 = t400 * pkin(9);
t600 = Ifges(3,4) * t401;
t599 = Ifges(4,4) * t395;
t598 = Ifges(4,4) * t400;
t596 = Ifges(6,4) * t394;
t595 = Ifges(6,4) * t399;
t593 = Ifges(5,6) * t395;
t592 = Ifges(5,6) * t400;
t566 = t317 * t386;
t565 = t317 * t387;
t564 = t317 * t394;
t563 = t317 * t399;
t559 = t392 * t397;
t554 = t394 * t400;
t546 = (t242 * t387 - t566) * mrSges(7,1) + (-t242 * t386 - t565) * mrSges(7,2);
t187 = t246 * t387 - t319 * t386;
t188 = t246 * t386 + t319 * t387;
t545 = t187 * mrSges(7,1) - t188 * mrSges(7,2);
t544 = (t315 * t387 + t386 * t557) * mrSges(7,1) + (-t315 * t386 + t387 * t557) * mrSges(7,2);
t207 = t400 * t290 + t395 * t291;
t355 = t400 * pkin(4) + t388;
t542 = t618 * pkin(1) + pkin(8) * t559;
t539 = qJD(2) * t401;
t527 = pkin(9) * t536;
t526 = pkin(9) * t535;
t524 = qJ(4) * t557;
t309 = t317 * pkin(2);
t520 = -t309 + t709;
t311 = t319 * pkin(2);
t519 = -t311 + t708;
t518 = Ifges(3,5) * t413 - Ifges(3,6) * t307 + Ifges(3,3) * t435;
t517 = t320 * pkin(2) + t542;
t511 = qJD(2) * t560;
t510 = t392 * t539;
t508 = t560 / 0.2e1;
t504 = -t541 / 0.2e1;
t503 = t541 / 0.2e1;
t498 = -t533 / 0.2e1;
t497 = -pkin(1) * t397 + pkin(8) * t515;
t495 = pkin(9) * t318 - t309;
t494 = pkin(9) * t320 - t311;
t133 = -t179 * mrSges(5,1) + t294 * mrSges(5,2);
t73 = t399 * t134 - t144 * t394;
t243 = t318 * t400 - t395 * t515;
t482 = mrSges(3,3) * t513;
t247 = t320 * t400 + t395 * t559;
t481 = t247 * pkin(3) + t517;
t472 = t401 * t504;
t471 = t401 * t503;
t469 = -t318 * pkin(2) + t497;
t466 = mrSges(4,1) * t315 + mrSges(4,2) * t316;
t239 = -t315 * t394 + t392 * t548;
t464 = mrSges(6,1) * t238 + mrSges(6,2) * t239;
t459 = -t315 * mrSges(5,2) - t316 * mrSges(5,3);
t457 = Ifges(4,1) * t400 - t599;
t456 = Ifges(6,1) * t399 - t596;
t455 = Ifges(6,1) * t394 + t595;
t454 = -Ifges(4,2) * t395 + t598;
t452 = -Ifges(6,2) * t394 + t595;
t451 = Ifges(6,2) * t399 + t596;
t449 = Ifges(6,5) * t399 - Ifges(6,6) * t394;
t448 = Ifges(6,5) * t394 + Ifges(6,6) * t399;
t447 = -Ifges(5,2) * t400 + t593;
t446 = Ifges(5,3) * t395 - t592;
t56 = pkin(5) * t316 + pkin(11) * t239 + t73;
t61 = pkin(11) * t238 + t74;
t24 = -t393 * t61 + t398 * t56;
t25 = t393 * t56 + t398 * t61;
t444 = t394 * t53 - t399 * t54;
t154 = t238 * t398 + t239 * t393;
t155 = t238 * t393 - t239 * t398;
t196 = t246 * t399 - t319 * t394;
t437 = -pkin(3) * t243 + t469;
t191 = t524 - t207;
t300 = qJD(2) * t432;
t302 = t324 * qJD(2);
t112 = -t290 * t535 - t291 * t536 + t300 * t400 - t395 * t302;
t237 = -t476 + (qJD(3) * t573 + t510) * t400;
t78 = pkin(4) * t237 - t511 * t651 - t112;
t236 = qJD(3) * t316 + t395 * t510;
t303 = t325 * qJD(2);
t406 = -t237 * qJ(4) - t316 * qJD(4) + t303;
t88 = t236 * t651 + t406;
t20 = t134 * t533 - t144 * t534 + t394 * t78 + t399 * t88;
t427 = qJ(4) * t246 + t481;
t111 = -t290 * t536 + t291 * t535 + t395 * t300 + t400 * t302;
t145 = -pkin(4) * t315 - t191;
t414 = -qJ(4) * t242 + t437;
t21 = -qJD(5) * t74 - t394 * t88 + t399 * t78;
t410 = Ifges(3,6) * t573 + (Ifges(3,2) * t401 + t601) * t392;
t45 = -pkin(4) * t180 - t60;
t98 = -qJ(4) * t511 + qJD(4) * t557 - t111;
t409 = -qJD(5) * t444 - t695;
t408 = t392 * t438 * (Ifges(3,5) * t401 - Ifges(3,6) * t396);
t79 = -pkin(4) * t236 - t98;
t365 = Ifges(3,4) * t513;
t345 = -pkin(3) * t400 + t496;
t321 = (-mrSges(3,1) * t401 + mrSges(3,2) * t396) * t392;
t314 = pkin(5) * t549 + t355;
t305 = t441 * t400;
t304 = t440 * t400;
t297 = -mrSges(3,2) * t438 + t482;
t252 = Ifges(3,1) * t514 + Ifges(3,5) * t438 + t365;
t251 = Ifges(3,6) * qJD(2) + qJD(1) * t410;
t227 = -t326 * t394 + t333;
t223 = mrSges(5,1) * t273 - mrSges(5,2) * t351;
t221 = -mrSges(4,1) * t351 - mrSges(4,3) * t273;
t220 = mrSges(4,2) * t351 - mrSges(4,3) * t272;
t205 = -mrSges(5,2) * t272 - mrSges(5,3) * t273;
t203 = pkin(3) * t273 + t572;
t197 = t246 * t394 + t319 * t399;
t195 = -pkin(3) * t514 - t209;
t141 = qJD(5) * t238 + t236 * t394 + t399 * t511;
t140 = qJD(5) * t239 + t236 * t399 - t394 * t511;
t132 = mrSges(5,1) * t180 - mrSges(5,3) * t294;
t131 = -mrSges(4,2) * t294 - mrSges(4,3) * t180;
t130 = mrSges(4,1) * t294 + mrSges(4,3) * t179;
t107 = t236 * pkin(3) + t406;
t106 = -pkin(3) * t511 - t112;
t105 = t145 - t615;
t96 = mrSges(7,1) * t258 - mrSges(7,3) * t121;
t95 = -mrSges(7,2) * t258 + mrSges(7,3) * t488;
t94 = mrSges(4,1) * t180 - mrSges(4,2) * t179;
t93 = -mrSges(5,2) * t180 + mrSges(5,3) * t179;
t71 = -t148 * t394 + t123;
t67 = -mrSges(7,1) * t488 + mrSges(7,2) * t121;
t62 = t180 * pkin(3) + t404;
t52 = -pkin(5) * t140 + t79;
t51 = -qJD(6) * t155 + t140 * t398 - t141 * t393;
t50 = qJD(6) * t154 + t140 * t393 + t141 * t398;
t44 = -mrSges(6,1) * t92 + mrSges(6,2) * t91;
t28 = -pkin(5) * t92 + t45;
t23 = -mrSges(7,2) * t161 + mrSges(7,3) * t32;
t22 = mrSges(7,1) * t161 - mrSges(7,3) * t31;
t19 = t398 * t46 - t581;
t18 = -t393 * t46 - t577;
t17 = pkin(11) * t140 + t20;
t14 = pkin(5) * t237 - pkin(11) * t141 + t21;
t13 = -mrSges(7,1) * t32 + mrSges(7,2) * t31;
t5 = -qJD(6) * t25 + t14 * t398 - t17 * t393;
t4 = qJD(6) * t24 + t14 * t393 + t17 * t398;
t1 = [-t45 * t464 + t202 * t466 + t62 * t459 + (t60 * mrSges(5,1) - t69 * mrSges(4,3) - Ifges(4,4) * t640 + Ifges(5,6) * t639 + t754) * t315 + (-t15 * t50 + t154 * t2 - t155 * t3 + t16 * t51) * mrSges(7,3) - t690 * t529 + t573 * t518 / 0.2e1 + (Ifges(7,5) * t50 + Ifges(7,6) * t51) * t630 + (Ifges(7,5) * t155 + Ifges(7,6) * t154) * t644 + (t218 * t573 - t307 * t617 + t324 * t435) * mrSges(3,1) + (Ifges(3,1) * t413 - Ifges(3,4) * t307 + Ifges(3,5) * t435) * t508 + (t11 * t238 + t12 * t239 + t140 * t54 - t141 * t53) * mrSges(6,3) + (Ifges(6,4) * t141 + Ifges(6,2) * t140) * t635 + (-Ifges(6,4) * t239 + Ifges(6,2) * t238) * t652 + (-t217 * t573 - t325 * t435 - t413 * t617) * mrSges(3,2) + (-m(6) * t427 - t197 * mrSges(6,1) - t196 * mrSges(6,2) - m(3) * t542 - t320 * mrSges(3,1) - mrSges(3,3) * t559 - mrSges(2,1) * t618 + t397 * mrSges(2,2) - m(4) * (t517 + t611) - m(5) * (t427 + t611) - m(7) * t481 - t188 * mrSges(7,1) - t187 * mrSges(7,2) + t685 * t246 - t675 * t319 + t673 * t247) * g(2) + (-m(3) * t298 - t683) * t303 + m(4) * (t111 * t174 - t112 * t173 + t202 * t289 + t206 * t70 + t207 * t69) + (-m(6) * t414 + t563 * mrSges(6,1) - t564 * mrSges(6,2) - m(3) * t497 + t318 * mrSges(3,1) - mrSges(3,3) * t515 + t397 * mrSges(2,1) + mrSges(2,2) * t618 - m(5) * (t414 - t612) - m(4) * (t469 - t612) - m(7) * t437 + t565 * mrSges(7,1) - t566 * mrSges(7,2) - (t685 + t762) * t242 + t675 * t317 - t673 * t243) * g(1) + (Ifges(7,4) * t50 + Ifges(7,2) * t51) * t647 + (Ifges(7,4) * t155 + Ifges(7,2) * t154) * t661 + (Ifges(4,5) * t624 + Ifges(5,4) * t625 + Ifges(5,5) * t626 + Ifges(4,6) * t627 - t251 / 0.2e1 - t149 * mrSges(5,3) - t174 * mrSges(4,2) + t147 * mrSges(5,2) - t173 * mrSges(4,1) + t734 * t621) * t511 + (-t218 * t560 - t298 * t510 - t301 * t511 - t307 * t325 - t324 * t413) * mrSges(3,3) + (Ifges(3,4) * t747 - Ifges(5,4) * t639 - Ifges(4,5) * t640 - Ifges(5,5) * t637 + Ifges(3,2) * t748 + Ifges(3,6) * t746 - Ifges(4,6) * t638 - t623 * t734 + t671 + t217 * mrSges(3,3) - t700 / 0.2e1) * t557 + (Ifges(3,3) * t573 + (Ifges(3,5) * t396 + Ifges(3,6) * t401) * t392) * t746 + (Ifges(3,5) * t573 + (t396 * Ifges(3,1) + t600) * t392) * t747 + t410 * t748 + m(5) * (t106 * t147 + t107 * t143 + t149 * t98 + t190 * t62 + t191 * t60 + t192 * t66) + m(6) * (t11 * t74 + t110 * t79 + t12 * t73 + t145 * t45 + t20 * t54 + t21 * t53) + m(7) * (t105 * t28 + t15 * t5 + t16 * t4 + t2 * t25 + t24 * t3 + t52 * t84) + (Ifges(6,5) * t141 + Ifges(6,6) * t140) * t628 + (-Ifges(6,5) * t239 + Ifges(6,6) * t238) * t641 + (Ifges(7,1) * t50 + Ifges(7,4) * t51) * t645 + (Ifges(7,1) * t155 + Ifges(7,4) * t154) * t662 + (t716 * t508 + t408 / 0.2e1) * qJD(2) + (Ifges(6,1) * t141 + Ifges(6,4) * t140) * t633 + (-Ifges(6,1) * t239 + Ifges(6,4) * t238) * t653 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t667 + t217 * t325 + t218 * t324 + t301 * t302) + (t667 * qJD(1) * (-Ifges(3,2) * t396 + t600) + t392 * t252) * t539 / 0.2e1 + Ifges(2,3) * qJDD(1) + t24 * t22 + t25 * t23 - pkin(1) * t321 * t528 + t302 * t297 + t289 * t94 + t238 * t36 / 0.2e1 + t106 * t223 + t111 * t220 + t112 * t221 + t98 * t222 + t107 * t205 + t206 * t130 + t207 * t131 + (t66 * mrSges(5,1) - t70 * mrSges(4,3) + Ifges(4,4) * t638 - Ifges(5,6) * t637 + t768) * t316 + (Ifges(4,1) * t624 + Ifges(4,4) * t627 + Ifges(6,5) * t633 - Ifges(5,2) * t625 - Ifges(5,6) * t626 + Ifges(6,6) * t635 + Ifges(6,3) * t628 + t621 * t736 + t751 + t770) * t237 + (-Ifges(4,4) * t624 + Ifges(5,6) * t625 + Ifges(5,3) * t626 - Ifges(4,2) * t627 + t162 / 0.2e1 + t642 - t735 * t621 + t771) * t236 + t52 * t67 + t73 * t63 + t74 * t64 + t84 * (-mrSges(7,1) * t51 + mrSges(7,2) * t50) + t4 * t95 + t5 * t96 + t105 * t13 + t140 * t649 + t50 * t655 + t51 * t657 + t239 * t659 + t155 * t663 + t154 * t665 + t79 * t127 + t141 * t103 / 0.2e1 + t110 * (-mrSges(6,1) * t140 + mrSges(6,2) * t141) + t145 * t44 + t20 * t150 + t21 * t151 + t28 * (-mrSges(7,1) * t154 + mrSges(7,2) * t155) + t190 * t93 + t191 * t132 + t192 * t133; -t593 * t637 + t599 * t638 + ((t457 / 0.2e1 - t447 / 0.2e1) * t273 + t702 * t621 + (-t454 / 0.2e1 + t446 / 0.2e1) * t272 + (Ifges(6,3) * t400 + t395 * t448) * t628 + (Ifges(6,5) * t400 + t395 * t455) * t633 + (Ifges(6,6) * t400 + t395 * t451) * t635 - t703 + ((t147 * t400 + t149 * t395) * m(5) + (t173 * t400 - t174 * t395) * m(4)) * pkin(9)) * qJD(3) - t202 * t465 + (-t15 * t477 - t28 * t304 - t718 * t84) * mrSges(7,1) + t62 * t458 + (t103 * t498 + t164 * t471 - t448 * t641 + t45 * t463 - t451 * t652 - t455 * t653 - t754) * t400 + (-t408 / 0.2e1 + t690 * qJD(1)) * qJD(1) + (Ifges(7,4) * t158 + Ifges(7,2) * t159) * t647 + (-t526 - t209) * t221 + (t697 * pkin(9) + t707 * t143 - t147 * t195 - t149 * t193 + t345 * t62) * m(5) + (Ifges(7,5) * t158 + Ifges(7,6) * t159) * t630 + (t173 * (mrSges(4,1) * t396 - mrSges(4,3) * t547) - t147 * (mrSges(5,1) * t547 + mrSges(5,2) * t396) - t174 * (-mrSges(4,2) * t396 - mrSges(4,3) * t552) - t149 * (mrSges(5,1) * t552 - mrSges(5,3) * t396)) * t541 + (-m(5) * (t494 + t708) - m(4) * t494 - m(6) * (-pkin(10) * t561 + t519) + mrSges(6,3) * t561 - m(7) * t519 + t670 * t320 + t672 * t319) * g(1) + (-m(5) * (t495 + t709) - m(4) * t495 - m(6) * (-pkin(10) * t562 + t520) + mrSges(6,3) * t562 - m(7) * t520 + t670 * t318 + t672 * t317) * g(2) + (-t11 * t549 + t12 * t554 + t53 * t710 - t54 * t711) * mrSges(6,3) + t712 * t67 + (-Ifges(7,4) * t305 + Ifges(7,2) * t304) * t661 + (-Ifges(7,5) * t305 + Ifges(7,6) * t304) * t644 + (-Ifges(7,1) * t305 + Ifges(7,4) * t304) * t662 + (t16 * t477 - t28 * t305 + t719 * t84) * mrSges(7,2) + (-t195 + t526) * t223 + (-Ifges(3,2) * t514 + t400 * t689 + t252 + t365) * t472 - t592 * t639 + t598 * t640 + (t483 + t683) * t301 + (-t110 * t710 + t477 * t54) * mrSges(6,2) + (Ifges(7,1) * t158 + Ifges(7,4) * t159) * t645 + t721 * t150 + (-t449 * t628 - t452 * t635 - t456 * t633) * t532 + (t273 * (Ifges(4,5) * t396 + t401 * t457) + t272 * (Ifges(5,5) * t396 + t401 * t446) + t716 * t396) * t504 + (-m(4) * t543 + t321 - t693 * (t400 * t375 + t543) + (-t430 * mrSges(6,1) - t431 * mrSges(6,2) + (-mrSges(6,3) - t749) * t547 + t764 * t401 + (-m(6) * pkin(4) - m(7) * t384 - t461 - t738) * t396) * t392) * g(3) + (Ifges(6,5) * t263 + Ifges(6,6) * t262 + Ifges(6,3) * t477) * t629 + (-pkin(2) * t202 + pkin(9) * t698 + t173 * t209 - t174 * t210) * m(4) + t703 * t513 + t706 * t127 + t707 * t205 + (t556 + t717) * t536 / 0.2e1 + (-t15 * t719 + t16 * t718 + t2 * t304 + t3 * t305) * mrSges(7,3) + (t273 * (Ifges(5,4) * t396 + t401 * t447) + t272 * (Ifges(4,6) * t396 + t401 * t454) + t396 * t251 + (t396 * t734 + t401 * t702) * t351) * t503 + (-t193 + t527) * t222 + t732 * t96 + (t115 * t3 + t116 * t2 + t15 * t732 + t16 * t733 + t28 * t314 + t712 * t84) * m(7) + t733 * t95 + (t482 - t297) * t298 + (-t527 - t210) * t220 + (t110 * t711 - t477 * t53) * mrSges(6,1) + t518 + (t642 + t696) * t536 + t697 * mrSges(5,1) + t698 * mrSges(4,3) + t751 * t535 + (t11 * t228 + t110 * t706 + t12 * t227 + t355 * t45 + t53 * t722 + t54 * t721) * m(6) + t722 * t151 + t355 * t44 + t345 * t93 + t314 * t13 - t262 * t102 / 0.2e1 - t263 * t103 / 0.2e1 + t227 * t63 + t228 * t64 - t217 * mrSges(3,2) + t218 * mrSges(3,1) + (-t132 + t131) * t388 + (Ifges(7,5) * t186 + Ifges(7,6) * t185 + Ifges(7,3) * t477) * t631 + (Ifges(6,1) * t263 + Ifges(6,4) * t262 + Ifges(6,5) * t477) * t634 + (Ifges(6,4) * t263 + Ifges(6,2) * t262 + Ifges(6,6) * t477) * t636 + (t165 * t471 + t162 * t472 - t693 * t524 * g(3) + (-t130 + t133) * pkin(9) + t768) * t395 - pkin(2) * t94 + t115 * t22 + t116 * t23 + (Ifges(7,1) * t186 + Ifges(7,4) * t185 + Ifges(7,5) * t477) * t646 + (Ifges(7,4) * t186 + Ifges(7,2) * t185 + Ifges(7,6) * t477) * t648 + t509 * t649 + t158 * t655 + t186 * t656 + t159 * t657 + t185 * t658 + t554 * t659 + t549 * t660 - t305 * t663 + t304 * t665; t102 * t498 + t45 * t462 + t433 * t127 + (Ifges(7,1) * t426 - Ifges(7,4) * t720) * t646 + (Ifges(7,4) * t426 - Ifges(7,2) * t720) * t648 + (Ifges(7,5) * t426 - Ifges(7,6) * t720) * t631 - t441 * t665 + t28 * (mrSges(7,1) * t441 - mrSges(7,2) * t440) + (-Ifges(7,5) * t440 - Ifges(7,6) * t441) * t644 + (-Ifges(7,4) * t440 - Ifges(7,2) * t441) * t661 + (-Ifges(7,1) * t440 - Ifges(7,4) * t441) * t662 - t440 * t663 + (-Ifges(7,5) * t230 - Ifges(7,6) * t229) * t630 + (-Ifges(7,1) * t230 - Ifges(7,4) * t229) * t645 + (-Ifges(7,4) * t230 - Ifges(7,2) * t229) * t647 + (-pkin(3) * t66 - qJ(4) * t60 - t143 * t203 - t147 * t174 + t149 * t694) * m(5) + (t45 * qJ(4) + t758 * t110 - t53 * t71 - t54 * t72) * m(6) + (t165 + t582) * t624 + (mrSges(7,1) * t714 + mrSges(7,2) * t713) * t84 + t714 * t658 + t715 * qJD(4) + t687 * qJD(5) + (t221 - t223) * t174 + (t466 - m(5) * t489 + t459 - m(6) * (-t308 - t607) - m(7) * (t315 * t402 - t308) + t704 * t316) * g(3) + t705 * t67 + (-t583 + t717) * t625 + (t265 + t164) * t627 + (t220 - t222) * t173 + (t44 - t132) * qJ(4) - t671 + (-t266 + t689) * t626 + ((-t533 - t567) * t54 + (t273 * t394 + t534) * t53 + t676 + t695) * mrSges(6,3) + t724 * t95 + t725 * t96 + (t15 * t725 + t16 * t724 + t2 * t241 + t240 * t3 + t28 * t381 + t705 * t84) * m(7) - (t215 * t451 + t216 * t455 + t267 * t448) * qJD(5) / 0.2e1 + t426 * t656 + (-m(6) * t409 - t753) * t651 + (t242 * t759 + t243 * t678) * g(2) + (t246 * t759 + t247 * t678) * g(1) + t399 * t37 / 0.2e1 + t381 * t13 + t240 * t22 + t241 * t23 - t203 * t205 + t700 + (-Ifges(4,1) * t625 - Ifges(6,5) * t634 - Ifges(7,5) * t646 + Ifges(5,2) * t624 - Ifges(6,6) * t636 - Ifges(7,6) * t648 - Ifges(6,3) * t629 - Ifges(7,3) * t631 - t620 * t736 + t769 + t770) * t272 + (-Ifges(4,2) * t626 + Ifges(5,3) * t627 + t448 * t629 + t451 * t636 + t455 * t634 - t620 * t735 + t687 - t771) * t273 + t772 * mrSges(7,3) + t449 * t641 + t452 * t652 + t456 * t653 - t230 * t655 + t394 * t660 - pkin(3) * t133 - t72 * t150 - t71 * t151; -t440 * t22 + t441 * t23 + t713 * t96 + t714 * t95 + (t67 + t715) * t351 + (t205 + t443) * t273 + t133 + (t351 * t84 - t772) * m(7) + (t110 * t351 - t273 * t444 + t409 - t676) * m(6) + (t143 * t273 - t149 * t351 + t66 - t676) * m(5) + t753; t756 * t648 - t121 * t658 + (-m(7) * t615 - t464 - t544) * g(3) + (t602 + t151) * t54 - t67 * t616 - m(7) * (t15 * t18 + t16 * t19 + t616 * t84) + (t23 * t393 - t530 * t96 + (qJD(6) * t95 + t22) * t398) * pkin(5) + (Ifges(6,5) * t215 - Ifges(6,6) * t216) * t629 + t682 + (mrSges(6,2) * t197 + t196 * t723 - t545) * g(1) + (-t546 - (-t242 * t394 - t563) * mrSges(6,2) + t723 * (t242 * t399 - t564)) * g(2) + t35 + t488 * t656 + t767 + (t603 - t150) * t53 - t110 * (mrSges(6,1) * t216 + mrSges(6,2) * t215) + t102 * t633 + (Ifges(6,1) * t215 - t597) * t634 - t19 * t95 - t18 * t96 + (-Ifges(6,2) * t216 + t103 + t214) * t636 + (t2 * t393 + t3 * t398 + (-t15 * t393 + t16 * t398) * qJD(6)) * t664; t58 * t645 - t15 * t95 + t16 * t96 - g(1) * t545 - g(2) * t546 - g(3) * t544 + (t59 + t756) * t648 + t767;];
tau  = t1;
