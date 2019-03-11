% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR13_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR13_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR13_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR13_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:51:08
% EndTime: 2019-03-09 23:52:58
% DurationCPUTime: 65.92s
% Computational Cost: add. (23194->1171), mult. (55221->1533), div. (0->0), fcn. (43946->12), ass. (0->510)
t430 = sin(qJ(3));
t435 = cos(qJ(2));
t427 = sin(pkin(6));
t583 = qJD(1) * t427;
t551 = t435 * t583;
t520 = t430 * t551;
t581 = qJD(3) * t430;
t720 = t520 - t581;
t429 = sin(qJ(4));
t431 = sin(qJ(2));
t433 = cos(qJ(4));
t434 = cos(qJ(3));
t588 = t434 * t435;
t308 = (t429 * t431 + t433 * t588) * t427;
t293 = qJD(1) * t308;
t577 = qJD(4) * t430;
t806 = t429 * t577 + t293;
t770 = Ifges(5,1) + Ifges(6,1);
t768 = Ifges(5,5) + Ifges(6,4);
t802 = Ifges(5,6) - Ifges(6,6);
t766 = Ifges(5,3) + Ifges(6,2);
t607 = cos(pkin(6));
t531 = t607 * qJD(1);
t518 = pkin(1) * t531;
t555 = t431 * t583;
t343 = -pkin(8) * t555 + t435 * t518;
t471 = (pkin(2) * t431 - pkin(9) * t435) * t427;
t344 = qJD(1) * t471;
t232 = t434 * t343 + t430 * t344;
t208 = pkin(10) * t555 + t232;
t558 = pkin(1) * t607;
t420 = t431 * t558;
t510 = pkin(3) * t430 - pkin(10) * t434;
t596 = t427 * t435;
t241 = (t420 + (pkin(8) + t510) * t596) * qJD(1);
t379 = t510 * qJD(3);
t639 = pkin(3) * t434;
t389 = -pkin(10) * t430 - pkin(2) - t639;
t589 = t433 * t434;
t422 = pkin(9) * t589;
t578 = qJD(4) * t429;
t790 = -qJD(4) * t422 + t429 * t208 - t389 * t578 + (-t241 + t379) * t433;
t127 = t433 * t208 + t429 * t241;
t805 = -t720 * qJ(5) - t127;
t573 = qJD(1) * qJD(2);
t352 = (-qJDD(1) * t435 + t431 * t573) * t427;
t331 = qJDD(3) + t352;
t648 = t331 / 0.2e1;
t412 = t531 + qJD(2);
t306 = t412 * t430 + t434 * t555;
t353 = (qJDD(1) * t431 + t435 * t573) * t427;
t528 = t607 * qJDD(1);
t411 = t528 + qJDD(2);
t192 = -qJD(3) * t306 - t353 * t430 + t434 * t411;
t665 = t192 / 0.2e1;
t186 = qJDD(4) - t192;
t667 = t186 / 0.2e1;
t173 = qJDD(6) - t186;
t668 = t173 / 0.2e1;
t305 = t434 * t412 - t430 * t555;
t191 = qJD(3) * t305 + t353 * t434 + t411 * t430;
t484 = -qJD(3) + t551;
t243 = t433 * t306 - t429 * t484;
t100 = qJD(4) * t243 + t429 * t191 - t433 * t331;
t675 = t100 / 0.2e1;
t676 = -t100 / 0.2e1;
t242 = t429 * t306 + t433 * t484;
t579 = qJD(4) * t242;
t99 = t433 * t191 + t429 * t331 - t579;
t680 = t99 / 0.2e1;
t428 = sin(qJ(6));
t432 = cos(qJ(6));
t481 = t242 * t428 + t243 * t432;
t34 = -qJD(6) * t481 + t100 * t432 - t428 * t99;
t688 = t34 / 0.2e1;
t132 = t242 * t432 - t243 * t428;
t33 = qJD(6) * t132 + t100 * t428 + t432 * t99;
t689 = t33 / 0.2e1;
t297 = qJD(4) - t305;
t678 = pkin(4) + pkin(5);
t279 = -t412 * pkin(2) - t343;
t153 = -t305 * pkin(3) - t306 * pkin(10) + t279;
t585 = pkin(8) * t596 + t420;
t346 = t585 * qJD(1);
t280 = t412 * pkin(9) + t346;
t294 = (-pkin(2) * t435 - pkin(9) * t431 - pkin(1)) * t583;
t188 = t434 * t280 + t430 * t294;
t157 = -pkin(10) * t484 + t188;
t81 = t433 * t153 - t429 * t157;
t711 = -t81 + qJD(5);
t797 = -pkin(11) * t243 + t711;
t49 = -t297 * t678 + t797;
t287 = t297 * qJ(5);
t82 = t429 * t153 + t433 * t157;
t60 = pkin(11) * t242 + t82;
t54 = t287 + t60;
t17 = -t428 * t54 + t432 * t49;
t576 = qJD(4) * t433;
t787 = -pkin(8) * t427 * t573 + pkin(1) * t528;
t572 = qJDD(1) * t427;
t788 = pkin(8) * t572 + qJD(2) * t518;
t244 = t431 * t787 + t435 * t788;
t221 = pkin(9) * t411 + t244;
t565 = pkin(1) * t572;
t229 = pkin(2) * t352 - pkin(9) * t353 - t565;
t580 = qJD(3) * t434;
t78 = t434 * t221 + t430 * t229 - t280 * t581 + t294 * t580;
t69 = pkin(10) * t331 + t78;
t245 = -t431 * t788 + t435 * t787;
t222 = -t411 * pkin(2) - t245;
t76 = -t192 * pkin(3) - t191 * pkin(10) + t222;
t15 = -t153 * t578 - t157 * t576 - t429 * t69 + t433 * t76;
t462 = qJDD(5) - t15;
t8 = -pkin(11) * t99 - t186 * t678 + t462;
t14 = t153 * t576 - t157 * t578 + t429 * t76 + t433 * t69;
t12 = t186 * qJ(5) + t297 * qJD(5) + t14;
t9 = pkin(11) * t100 + t12;
t1 = qJD(6) * t17 + t428 * t8 + t432 * t9;
t18 = t428 * t49 + t432 * t54;
t2 = -qJD(6) * t18 - t428 * t9 + t432 * t8;
t703 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t13 = -pkin(4) * t186 + t462;
t778 = t15 * mrSges(5,1) - t13 * mrSges(6,1) - t14 * mrSges(5,2) + t12 * mrSges(6,3);
t804 = -Ifges(7,5) * t689 - Ifges(4,2) * t665 - Ifges(4,6) * t648 + Ifges(5,6) * t676 + Ifges(6,6) * t675 - Ifges(7,6) * t688 - Ifges(7,3) * t668 + t667 * t766 + t680 * t768 - t703 + t778;
t803 = t768 * t667 + t770 * t680;
t769 = -Ifges(5,4) + Ifges(6,5);
t559 = -pkin(9) * t429 - pkin(4);
t801 = t520 * t678 + (-pkin(5) + t559) * t581 - t790 + (-qJD(3) * t589 + t806) * pkin(11);
t562 = t429 * t596;
t524 = t434 * t562;
t292 = qJD(1) * t524 - t433 * t555;
t587 = t429 * t379 + t389 * t576;
t591 = t430 * t433;
t800 = -pkin(11) * t292 + (-pkin(9) * qJD(3) + pkin(11) * qJD(4)) * t591 + (-qJD(5) + (-pkin(9) * qJD(4) + pkin(11) * qJD(3)) * t429) * t434 + t587 + t805;
t603 = t305 * t429;
t677 = pkin(10) - pkin(11);
t187 = -t430 * t280 + t434 * t294;
t226 = pkin(3) * t306 - pkin(10) * t305;
t109 = t433 * t187 + t429 * t226;
t88 = t306 * qJ(5) + t109;
t799 = pkin(11) * t603 + t677 * t578 + t88;
t165 = t429 * t187;
t397 = t677 * t433;
t798 = qJD(4) * t397 - t165 - (-pkin(11) * t305 - t226) * t433 + t678 * t306;
t605 = qJ(5) * t429;
t638 = pkin(4) * t433;
t486 = t605 + t638;
t130 = Ifges(7,4) * t132;
t796 = Ifges(7,2) * t481 - t130;
t690 = m(7) * pkin(5);
t795 = -mrSges(5,1) - mrSges(6,1) - t690;
t71 = -pkin(4) * t297 + t711;
t72 = t287 + t82;
t755 = t484 * Ifges(4,6);
t793 = -t279 * mrSges(4,1) - t81 * mrSges(5,1) + t71 * mrSges(6,1) + t17 * mrSges(7,1) + t82 * mrSges(5,2) - t18 * mrSges(7,2) - t72 * mrSges(6,3) - t755 / 0.2e1;
t5 = Ifges(7,5) * t33 + Ifges(7,6) * t34 + Ifges(7,3) * t173;
t156 = pkin(3) * t484 - t187;
t441 = t243 * qJ(5) - t156;
t58 = -t242 * t678 + t441;
t621 = Ifges(7,4) * t481;
t288 = qJD(6) - t297;
t657 = -t288 / 0.2e1;
t671 = -t481 / 0.2e1;
t792 = (t132 * t17 + t18 * t481) * mrSges(7,3) - t58 * (mrSges(7,1) * t481 + mrSges(7,2) * t132) + (Ifges(7,5) * t132 - Ifges(7,6) * t481) * t657 + (Ifges(7,1) * t132 - t621) * t671 + t5 + t703;
t791 = t769 * t675 + t803;
t626 = Ifges(4,4) * t306;
t751 = Ifges(7,5) * t481 + Ifges(4,2) * t305 + t132 * Ifges(7,6) + t288 * Ifges(7,3) + t626 - t755;
t744 = -t242 * t802 + t243 * t768 + t297 * t766;
t240 = Ifges(5,4) * t242;
t614 = t242 * Ifges(6,5);
t743 = t243 * t770 + t768 * t297 - t240 + t614;
t789 = qJD(5) * t429 + t188;
t696 = -mrSges(7,1) * t432 + mrSges(7,2) * t428 + t795;
t786 = -t99 * Ifges(6,5) / 0.2e1 - t186 * Ifges(6,6) / 0.2e1 + Ifges(5,4) * t680 + Ifges(5,6) * t667 + (Ifges(6,3) + Ifges(5,2)) * t676;
t642 = cos(qJ(1));
t512 = t607 * t642;
t641 = sin(qJ(1));
t363 = t431 * t512 + t435 * t641;
t557 = t427 * t642;
t269 = t363 * t434 - t430 * t557;
t362 = t431 * t641 - t435 * t512;
t209 = t269 * t429 - t362 * t433;
t210 = t269 * t433 + t362 * t429;
t479 = t428 * t429 + t432 * t433;
t480 = t428 * t433 - t429 * t432;
t781 = mrSges(7,1) * t479 - mrSges(7,2) * t480;
t666 = t191 / 0.2e1;
t780 = Ifges(4,1) * t666 + Ifges(4,5) * t648;
t779 = mrSges(6,2) * t13 - mrSges(5,3) * t15 + t791;
t679 = m(7) + m(6);
t775 = t191 * Ifges(4,4) + t192 * Ifges(4,2) + t331 * Ifges(4,6) + t5;
t773 = mrSges(5,2) - mrSges(6,3);
t772 = -mrSges(6,2) - mrSges(5,3);
t771 = mrSges(4,3) - mrSges(3,2);
t765 = -t100 * t802 + t186 * t766 + t768 * t99;
t592 = t429 * t434;
t421 = pkin(9) * t592;
t425 = t434 * pkin(4);
t246 = pkin(5) * t434 + t421 + t425 + (-pkin(11) * t430 - t389) * t433;
t319 = t429 * t389 + t422;
t284 = -qJ(5) * t434 + t319;
t593 = t429 * t430;
t253 = pkin(11) * t593 + t284;
t147 = t246 * t428 + t253 * t432;
t763 = -qJD(6) * t147 - t428 * t800 + t432 * t801;
t146 = t246 * t432 - t253 * t428;
t762 = qJD(6) * t146 + t428 * t801 + t432 * t800;
t758 = Ifges(4,3) * t484;
t756 = t484 * Ifges(4,5);
t231 = -t430 * t343 + t434 * t344;
t463 = pkin(3) * t555 + t231;
t444 = qJ(5) * t293 + t463;
t604 = qJ(5) * t433;
t470 = -t429 * t678 + t604;
t460 = -pkin(9) + t470;
t574 = qJD(5) * t433;
t705 = -t433 * t678 - t605;
t754 = t292 * t678 + (qJD(4) * t705 + t574) * t430 + t460 * t580 - t444;
t396 = t677 * t429;
t282 = t396 * t432 - t397 * t428;
t753 = qJD(6) * t282 + t428 * t798 - t432 * t799;
t283 = t396 * t428 + t397 * t432;
t752 = -qJD(6) * t283 + t428 * t799 + t432 * t798;
t750 = t297 * t470 + t789;
t381 = -t428 * qJ(5) - t432 * t678;
t749 = qJD(6) * t381 - t428 * t60 + t432 * t797;
t382 = t432 * qJ(5) - t428 * t678;
t748 = -qJD(6) * t382 - t428 * t797 - t432 * t60;
t507 = mrSges(4,1) * t434 - mrSges(4,2) * t430;
t747 = t507 + mrSges(3,1);
t519 = m(7) * t677 - mrSges(7,3);
t746 = -t519 + mrSges(4,2);
t349 = t480 * t430;
t745 = t484 * (Ifges(4,5) * t434 - Ifges(4,6) * t430);
t710 = qJD(4) - qJD(6);
t168 = t349 * t710 + t479 * t580;
t196 = t292 * t428 + t293 * t432;
t742 = t168 - t196;
t257 = t710 * t479;
t169 = t257 * t430 - t480 * t580;
t195 = t292 * t432 - t293 * t428;
t741 = t169 - t195;
t239 = Ifges(6,5) * t243;
t115 = t297 * Ifges(6,6) + t242 * Ifges(6,3) + t239;
t296 = Ifges(4,4) * t305;
t180 = Ifges(4,1) * t306 + t296 - t756;
t740 = t429 * t115 + t180;
t224 = (-t433 * t581 - t434 * t578) * pkin(9) + t587;
t739 = -qJD(5) * t434 + t224 + t805;
t193 = t480 * t305;
t258 = t710 * t480;
t738 = t193 - t258;
t194 = t479 * t305;
t737 = t194 - t257;
t736 = pkin(4) * t520 + t559 * t581 - t790;
t485 = pkin(4) * t429 - t604;
t472 = pkin(9) + t485;
t735 = -pkin(4) * t292 + (qJD(4) * t486 - t574) * t430 + t472 * t580 + t444;
t734 = t224 - t127;
t569 = pkin(9) * t581;
t733 = t429 * t569 + t790;
t529 = -t363 * t430 - t434 * t557;
t732 = t486 * t529;
t511 = t607 * t641;
t365 = -t431 * t511 + t435 * t642;
t556 = t427 * t641;
t272 = t365 * t430 - t434 * t556;
t731 = t486 * t272;
t730 = -t429 * t580 - t430 * t576 + t292;
t729 = -t433 * t580 + t806;
t597 = t427 * t431;
t360 = t430 * t597 - t434 * t607;
t351 = t360 * qJ(5);
t728 = -t429 * t351 - t360 * t638;
t727 = t297 * t485 - t789;
t368 = -pkin(8) * t597 + t435 * t558;
t726 = pkin(4) * t679 - t696;
t503 = mrSges(6,1) * t429 - mrSges(6,3) * t433;
t505 = mrSges(5,1) * t429 + mrSges(5,2) * t433;
t80 = t242 * pkin(4) - t441;
t725 = t156 * t505 + t80 * t503;
t724 = t429 * t768 + t433 * t802;
t723 = -t429 * t802 + t433 * t768;
t619 = Ifges(6,5) * t433;
t622 = Ifges(5,4) * t433;
t722 = t429 * t770 - t619 + t622;
t620 = Ifges(6,5) * t429;
t623 = Ifges(5,4) * t429;
t721 = t433 * t770 + t620 - t623;
t602 = t305 * t433;
t716 = t576 - t602;
t715 = -t578 + t603;
t79 = -t430 * t221 + t434 * t229 - t280 * t580 - t294 * t581;
t714 = -t430 * t79 + t434 * t78;
t713 = t14 * t433 - t15 * t429;
t712 = t12 * t433 + t13 * t429;
t628 = Ifges(3,4) * t431;
t694 = t427 ^ 2;
t709 = (t431 * (Ifges(3,1) * t435 - t628) / 0.2e1 - pkin(1) * (mrSges(3,1) * t431 + mrSges(3,2) * t435)) * t694;
t504 = mrSges(6,1) * t433 + mrSges(6,3) * t429;
t506 = mrSges(5,1) * t433 - mrSges(5,2) * t429;
t708 = t433 * t690 + t504 + t506 + t781;
t706 = t746 + t772;
t704 = mrSges(4,1) + t708;
t527 = mrSges(3,3) * t555;
t702 = -m(4) * t279 + mrSges(3,1) * t412 + mrSges(4,1) * t305 - mrSges(4,2) * t306 - t527;
t701 = m(7) * pkin(11) + mrSges(7,3) + t772;
t699 = -t428 * mrSges(7,1) - t432 * mrSges(7,2) + t773;
t698 = -mrSges(6,2) * t12 - mrSges(5,3) * t14 - t786;
t273 = t365 * t434 + t430 * t556;
t361 = t430 * t607 + t434 * t597;
t697 = -g(1) * t273 - g(2) * t269 - g(3) * t361;
t692 = Ifges(7,4) * t689 + Ifges(7,2) * t688 + Ifges(7,6) * t668;
t691 = Ifges(7,1) * t689 + Ifges(7,4) * t688 + Ifges(7,5) * t668;
t56 = Ifges(7,2) * t132 + Ifges(7,6) * t288 + t621;
t685 = -t56 / 0.2e1;
t684 = t56 / 0.2e1;
t57 = Ifges(7,1) * t481 + Ifges(7,5) * t288 + t130;
t683 = -t57 / 0.2e1;
t682 = t57 / 0.2e1;
t681 = Ifges(4,4) * t665 + t780;
t674 = t115 / 0.2e1;
t673 = -t132 / 0.2e1;
t672 = t132 / 0.2e1;
t670 = t481 / 0.2e1;
t664 = -t242 / 0.2e1;
t663 = t242 / 0.2e1;
t662 = -t243 / 0.2e1;
t661 = t243 / 0.2e1;
t656 = t288 / 0.2e1;
t654 = -t297 / 0.2e1;
t653 = t297 / 0.2e1;
t652 = -t305 / 0.2e1;
t651 = t305 / 0.2e1;
t650 = -t306 / 0.2e1;
t649 = t306 / 0.2e1;
t640 = pkin(1) * t427;
t636 = pkin(10) * t272;
t633 = t529 * pkin(10);
t630 = mrSges(5,3) * t242;
t629 = mrSges(5,3) * t243;
t627 = Ifges(3,4) * t435;
t625 = Ifges(4,4) * t430;
t624 = Ifges(4,4) * t434;
t613 = t243 * Ifges(5,4);
t612 = t305 * mrSges(4,3);
t611 = t306 * mrSges(4,3);
t468 = pkin(3) * t331 + t79;
t610 = t430 * t468;
t606 = qJ(5) * t242;
t600 = t362 * t430;
t364 = t431 * t642 + t435 * t511;
t598 = t364 * t430;
t118 = -t242 * Ifges(5,2) + t297 * Ifges(5,6) + t613;
t594 = t429 * t118;
t590 = t430 * t435;
t326 = -pkin(2) * t607 - t368;
t354 = t360 * pkin(3);
t532 = t361 * pkin(10) - t354;
t199 = t326 - t532;
t327 = pkin(9) * t607 + t585;
t586 = pkin(2) * t596 + pkin(9) * t597;
t328 = -t586 - t640;
t228 = t434 * t327 + t430 * t328;
t201 = -pkin(10) * t596 + t228;
t104 = t429 * t199 + t433 * t201;
t227 = -t430 * t327 + t434 * t328;
t584 = t642 * pkin(1) + pkin(8) * t556;
t582 = qJD(2) * t435;
t416 = pkin(3) * t596;
t568 = pkin(9) * t580;
t567 = pkin(10) * t578;
t566 = pkin(10) * t576;
t563 = t427 * t590;
t561 = Ifges(4,5) * t191 + Ifges(4,6) * t192 + Ifges(4,3) * t331;
t91 = t351 + t104;
t200 = t416 - t227;
t560 = Ifges(3,5) * t353 - Ifges(3,6) * t352 + Ifges(3,3) * t411;
t554 = qJD(2) * t597;
t553 = t427 * t582;
t549 = t597 / 0.2e1;
t548 = -t594 / 0.2e1;
t546 = -t583 / 0.2e1;
t545 = t583 / 0.2e1;
t541 = t580 / 0.2e1;
t538 = -t577 / 0.2e1;
t537 = t576 / 0.2e1;
t63 = -t186 * mrSges(6,1) + t99 * mrSges(6,2);
t536 = -t362 * pkin(2) + pkin(9) * t363;
t535 = -t364 * pkin(2) + pkin(9) * t365;
t259 = t529 * pkin(3);
t534 = pkin(10) * t269 + t259;
t261 = t272 * pkin(3);
t533 = pkin(10) * t273 - t261;
t108 = t226 * t433 - t165;
t103 = t199 * t433 - t429 * t201;
t318 = t389 * t433 - t421;
t526 = mrSges(3,3) * t551;
t345 = qJD(2) * t471;
t347 = t368 * qJD(2);
t129 = -t327 * t580 - t328 * t581 + t434 * t345 - t430 * t347;
t523 = pkin(10) * t563 + t434 * t416 + t586;
t515 = t435 * t546;
t513 = -pkin(1) * t641 + pkin(8) * t557;
t508 = mrSges(4,1) * t360 + mrSges(4,2) * t361;
t266 = t361 * t429 + t433 * t596;
t267 = t361 * t433 - t562;
t162 = t266 * t432 - t267 * t428;
t163 = t266 * t428 + t267 * t432;
t502 = mrSges(7,1) * t162 - mrSges(7,2) * t163;
t501 = Ifges(4,1) * t434 - t625;
t496 = -Ifges(4,2) * t430 + t624;
t495 = -Ifges(5,2) * t429 + t622;
t494 = Ifges(5,2) * t433 + t623;
t488 = Ifges(6,3) * t429 + t619;
t487 = -Ifges(6,3) * t433 + t620;
t68 = -pkin(11) * t267 - t360 * t678 - t103;
t75 = pkin(11) * t266 + t91;
t27 = -t428 * t75 + t432 * t68;
t29 = t428 * t68 + t432 * t75;
t106 = -mrSges(7,2) * t288 + mrSges(7,3) * t132;
t107 = mrSges(7,1) * t288 - mrSges(7,3) * t481;
t482 = t106 * t432 - t107 * t428;
t476 = -pkin(10) * t600 - t362 * t639 + t536;
t475 = -pkin(10) * t598 - t364 * t639 + t535;
t474 = t365 * pkin(2) + pkin(9) * t364 + t584;
t128 = -t327 * t581 + t328 * t580 + t430 * t345 + t434 * t347;
t124 = pkin(10) * t554 + t128;
t264 = qJD(3) * t361 + t430 * t553;
t265 = -qJD(3) * t360 + t434 * t553;
t348 = t585 * qJD(2);
t141 = t264 * pkin(3) - t265 * pkin(10) + t348;
t43 = -t429 * t124 + t141 * t433 - t199 * t578 - t201 * t576;
t467 = t273 * pkin(3) + t474;
t465 = t279 * (mrSges(4,1) * t430 + mrSges(4,2) * t434);
t464 = (t435 * Ifges(3,2) + t628) * t427;
t256 = t266 * pkin(4);
t105 = -qJ(5) * t267 + t200 + t256;
t42 = t433 * t124 + t429 * t141 + t199 * t576 - t201 * t578;
t213 = t273 * t429 - t364 * t433;
t459 = -g(1) * t213 - g(2) * t209 - g(3) * t266;
t456 = t412 * t427 * (Ifges(3,5) * t435 - Ifges(3,6) * t431);
t451 = -t363 * pkin(2) - t362 * pkin(9) + t513;
t449 = pkin(3) * t554 + t129;
t447 = -pkin(3) * t269 + t451;
t26 = t264 * qJ(5) + t360 * qJD(5) + t42;
t214 = t273 * t433 + t364 * t429;
t443 = t214 * pkin(4) + qJ(5) * t213 + t467;
t442 = -qJ(5) * t679 + t699;
t440 = qJ(5) * t99 + qJD(5) * t243 + t468;
t439 = -pkin(4) * t210 - qJ(5) * t209 + t447;
t150 = -qJD(4) * t266 + t265 * t433 + t429 * t554;
t438 = qJ(5) * t150 + qJD(5) * t267 + t449;
t408 = Ifges(3,4) * t551;
t385 = -pkin(3) - t486;
t372 = pkin(3) - t705;
t366 = (-mrSges(3,1) * t435 + mrSges(3,2) * t431) * t427;
t350 = t479 * t430;
t341 = -mrSges(3,2) * t412 + t526;
t332 = t472 * t430;
t307 = -t433 * t597 + t524;
t285 = -t318 + t425;
t281 = t460 * t430;
t277 = Ifges(3,1) * t555 + Ifges(3,5) * t412 + t408;
t276 = t412 * Ifges(3,6) + qJD(1) * t464;
t248 = -mrSges(4,1) * t484 - t611;
t247 = mrSges(4,2) * t484 + t612;
t238 = -t364 * t589 + t365 * t429;
t237 = -t364 * t592 - t365 * t433;
t236 = -t362 * t589 + t363 * t429;
t235 = -t362 * t592 - t363 * t433;
t178 = Ifges(4,5) * t306 + Ifges(4,6) * t305 - t758;
t161 = -mrSges(6,1) * t297 + mrSges(6,2) * t243;
t160 = mrSges(5,1) * t297 - t629;
t159 = -mrSges(5,2) * t297 - t630;
t158 = -mrSges(6,2) * t242 + mrSges(6,3) * t297;
t149 = -qJD(4) * t562 + t265 * t429 + t361 * t576 - t433 * t554;
t140 = -mrSges(4,2) * t331 + mrSges(4,3) * t192;
t139 = mrSges(4,1) * t331 - mrSges(4,3) * t191;
t138 = mrSges(5,1) * t242 + mrSges(5,2) * t243;
t137 = mrSges(6,1) * t242 - mrSges(6,3) * t243;
t136 = pkin(4) * t243 + t606;
t113 = t213 * t428 + t214 * t432;
t112 = t213 * t432 - t214 * t428;
t102 = -mrSges(4,1) * t192 + mrSges(4,2) * t191;
t101 = -t243 * t678 - t606;
t92 = -pkin(4) * t360 - t103;
t90 = -pkin(4) * t306 - t108;
t84 = -pkin(5) * t266 - t105;
t73 = -mrSges(7,1) * t132 + mrSges(7,2) * t481;
t64 = -mrSges(5,2) * t186 - mrSges(5,3) * t100;
t62 = mrSges(5,1) * t186 - mrSges(5,3) * t99;
t61 = -mrSges(6,2) * t100 + mrSges(6,3) * t186;
t53 = qJD(6) * t162 + t149 * t428 + t150 * t432;
t52 = -qJD(6) * t163 + t149 * t432 - t150 * t428;
t48 = mrSges(5,1) * t100 + mrSges(5,2) * t99;
t47 = mrSges(6,1) * t100 - mrSges(6,3) * t99;
t46 = pkin(4) * t149 - t438;
t35 = -pkin(4) * t264 - t43;
t25 = -t149 * t678 + t438;
t22 = -mrSges(7,2) * t173 + mrSges(7,3) * t34;
t21 = mrSges(7,1) * t173 - mrSges(7,3) * t33;
t20 = pkin(11) * t149 + t26;
t19 = -pkin(11) * t150 - t264 * t678 - t43;
t16 = pkin(4) * t100 - t440;
t11 = -t100 * t678 + t440;
t10 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t4 = -qJD(6) * t29 + t19 * t432 - t20 * t428;
t3 = qJD(6) * t27 + t19 * t428 + t20 * t432;
t6 = [(-mrSges(5,1) * t468 + mrSges(6,1) * t16 - Ifges(5,2) * t676 + Ifges(6,3) * t675 - t667 * t802 + t680 * t769 + t698) * t266 + (t674 + Ifges(6,3) * t663 - Ifges(5,2) * t664 - t72 * mrSges(6,2) - t82 * mrSges(5,3) + t156 * mrSges(5,1) - t118 / 0.2e1 + t80 * mrSges(6,1) + t769 * t661 - t802 * t653) * t149 + (-Ifges(4,4) * t666 - t78 * mrSges(4,3) - t775 / 0.2e1 + t765 / 0.2e1 + t804) * t360 + (t156 * mrSges(5,2) + Ifges(6,5) * t663 + t768 * t653 + t770 * t661 - t80 * mrSges(6,3) + t743 / 0.2e1 + t71 * mrSges(6,2) - t81 * mrSges(5,3) + Ifges(5,4) * t664) * t150 + (-m(5) * (t447 + t633) - m(6) * (t439 + t633) + mrSges(2,1) * t641 + mrSges(2,2) * t642 - m(7) * t439 - m(4) * t451 + t269 * mrSges(4,1) - m(3) * t513 + t363 * mrSges(3,1) - mrSges(3,3) * t557 + t771 * t362 - t696 * t210 - t699 * t209 + t706 * t529) * g(1) + (-m(4) * t474 - t273 * mrSges(4,1) - mrSges(2,1) * t642 + mrSges(2,2) * t641 - m(7) * t443 - t113 * mrSges(7,1) - t112 * mrSges(7,2) - m(3) * t584 - t365 * mrSges(3,1) - mrSges(3,3) * t556 - m(5) * (t467 + t636) - m(6) * (t443 + t636) - t771 * t364 + t795 * t214 + t773 * t213 + t706 * t272) * g(2) + (-t188 * t554 + t265 * t279 + t596 * t78) * mrSges(4,2) + (t1 * t162 - t163 * t2 - t17 * t53 + t18 * t52) * mrSges(7,3) + (Ifges(5,6) * t664 - Ifges(7,5) * t670 - Ifges(7,6) * t672 - Ifges(4,4) * t649 - Ifges(4,2) * t651 - Ifges(7,3) * t656 + Ifges(6,6) * t663 + t766 * t653 + t768 * t661 - t188 * mrSges(4,3) - t751 / 0.2e1 + t744 / 0.2e1 - t793) * t264 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t694 + t244 * t585 + t245 * t368 + t346 * t347) + (-t244 * t607 - t353 * t640 - t411 * t585) * mrSges(3,2) + m(5) * (t103 * t15 + t104 * t14 - t156 * t449 - t200 * t468 + t42 * t82 + t43 * t81) - t449 * t138 + (t245 * t607 - t352 * t640 + t368 * t411) * mrSges(3,1) + (t244 * t596 - t245 * t597 - t343 * t553 - t346 * t554 - t352 * t585 - t353 * t368) * mrSges(3,3) + (Ifges(3,4) * t353 - Ifges(3,2) * t352 + Ifges(3,6) * t411) * t596 / 0.2e1 + (Ifges(3,1) * t353 - Ifges(3,4) * t352 + Ifges(3,5) * t411) * t549 + (-t468 * mrSges(5,2) - t16 * mrSges(6,3) + Ifges(5,4) * t676 + Ifges(6,5) * t675 + t779 + t803) * t267 - t561 * t596 / 0.2e1 + t79 * (-mrSges(4,1) * t596 - mrSges(4,3) * t361) + (Ifges(7,4) * t163 + Ifges(7,2) * t162) * t688 + (Ifges(4,4) * t361 - Ifges(4,6) * t596) * t665 + (Ifges(4,4) * t265 + Ifges(4,6) * t554) * t651 - t366 * t565 + (-m(3) * t343 - t702) * t348 + (Ifges(7,5) * t163 + Ifges(7,6) * t162) * t668 + (Ifges(7,5) * t53 + Ifges(7,6) * t52) * t656 + t353 * (Ifges(3,5) * t607 + (Ifges(3,1) * t431 + t627) * t427) / 0.2e1 + (t178 * t549 + t456 / 0.2e1) * qJD(2) + t222 * t508 + t709 * t573 - t276 * t554 / 0.2e1 + t187 * (mrSges(4,1) * t554 - mrSges(4,3) * t265) - t352 * (Ifges(3,6) * t607 + t464) / 0.2e1 + t411 * (Ifges(3,3) * t607 + (Ifges(3,5) * t431 + Ifges(3,6) * t435) * t427) / 0.2e1 + t607 * t560 / 0.2e1 + (Ifges(4,5) * t361 - Ifges(4,3) * t596) * t648 - t484 * (Ifges(4,5) * t265 + Ifges(4,3) * t554) / 0.2e1 - t11 * t502 + (t694 * qJD(1) * (-Ifges(3,2) * t431 + t627) + t427 * t277) * t582 / 0.2e1 + Ifges(2,3) * qJDD(1) + t347 * t341 + t326 * t102 + t265 * t180 / 0.2e1 + t129 * t248 + t128 * t247 + t227 * t139 + t228 * t140 + t200 * t48 + t26 * t158 + t42 * t159 + t43 * t160 + t35 * t161 + t46 * t137 + (Ifges(4,1) * t361 - Ifges(4,5) * t596) * t666 + (Ifges(4,1) * t265 + Ifges(4,5) * t554) * t649 + t103 * t62 + t104 * t64 + t105 * t47 + t3 * t106 + t4 * t107 + m(6) * (t105 * t16 + t12 * t91 + t13 * t92 + t26 * t72 + t35 * t71 + t46 * t80) + m(7) * (t1 * t29 + t11 * t84 + t17 * t4 + t18 * t3 + t2 * t27 + t25 * t58) + t91 * t61 + t92 * t63 + t84 * t10 + t361 * t681 + t53 * t682 + t52 * t684 + m(4) * (t128 * t188 + t129 * t187 + t222 * t326 + t227 * t79 + t228 * t78) + (Ifges(7,4) * t53 + Ifges(7,2) * t52) * t672 + t27 * t21 + t29 * t22 + t163 * t691 + t162 * t692 + t58 * (-mrSges(7,1) * t52 + mrSges(7,2) * t53) + t25 * t73 + (Ifges(7,1) * t53 + Ifges(7,4) * t52) * t670 + (Ifges(7,1) * t163 + Ifges(7,4) * t162) * t689; (-t292 * t802 + t293 * t768 + t520 * t766) * t654 + (t115 * t537 + t16 * t503 + t488 * t675 + t495 * t676 + t667 * t723 + t680 * t721 + t681 + t780) * t430 + t779 * t591 + (pkin(9) * t140 - t804) * t434 + (-t456 / 0.2e1 - t709 * qJD(1)) * qJD(1) + t625 * t665 + (t568 + t463) * t138 + (t156 * t463 + t14 * t319 + t15 * t318 + (t156 * t580 - t610) * pkin(9) + t734 * t82 + t733 * t81) * m(5) + (-Ifges(3,2) * t555 + t180 * t434 + t277 + t408) * t515 - t505 * t610 + (Ifges(7,5) * t350 - Ifges(7,6) * t349) * t668 + (Ifges(7,4) * t350 - Ifges(7,2) * t349) * t688 + (Ifges(7,1) * t350 - Ifges(7,4) * t349) * t689 + t548 * t580 + (Ifges(7,1) * t168 + Ifges(7,4) * t169 - Ifges(7,5) * t581) * t670 + (Ifges(7,1) * t196 + Ifges(7,4) * t195 - Ifges(7,5) * t520) * t671 + (Ifges(7,4) * t168 + Ifges(7,2) * t169 - Ifges(7,6) * t581) * t672 + (Ifges(7,4) * t196 + Ifges(7,2) * t195 - Ifges(7,6) * t520) * t673 + (Ifges(5,4) * t293 - Ifges(5,2) * t292 + Ifges(5,6) * t520) * t663 + (Ifges(6,5) * t293 + Ifges(6,6) * t520 + Ifges(6,3) * t292) * t664 + (Ifges(7,5) * t168 + Ifges(7,6) * t169 - Ifges(7,3) * t581) * t656 + (Ifges(7,5) * t196 + Ifges(7,6) * t195 - Ifges(7,3) * t520) * t657 + (t306 * (Ifges(4,5) * t431 + t435 * t501) + t305 * (Ifges(4,6) * t431 + t435 * t496) + t431 * t178) * t546 + (t433 * t538 + t292 / 0.2e1) * t118 + ((Ifges(6,6) * t430 + t434 * t488) * t663 + (Ifges(5,6) * t430 + t434 * t495) * t664 + t465 - t745 / 0.2e1 + (t430 * t768 + t434 * t721) * t661 + (t430 * t766 + t434 * t723) * t653) * qJD(3) + (t292 * t769 + t293 * t770 + t520 * t768) * t662 + (-t487 * t663 - t494 * t664 - t653 * t724 - t661 * t722) * t577 + (-mrSges(5,1) * t720 + mrSges(5,3) * t729) * t81 + (mrSges(6,1) * t720 - mrSges(6,2) * t729) * t71 + t624 * t666 + (t527 + t702) * t346 + t762 * t106 + t763 * t107 + (t1 * t147 + t11 * t281 + t146 * t2 + t17 * t763 + t18 * t762 + t58 * t754) * m(7) - t765 * t434 / 0.2e1 + t435 * t545 * t745 + ((-mrSges(4,1) * t187 + mrSges(4,2) * t188) * t583 + (t276 + t758) * t545) * t431 + (t48 - t139) * pkin(9) * t430 + (mrSges(5,2) * t720 + mrSges(5,3) * t730) * t82 + (mrSges(6,2) * t730 - mrSges(6,3) * t720) * t72 + (-mrSges(6,1) * t730 + mrSges(6,3) * t729) * t80 + (-mrSges(5,1) * t730 - mrSges(5,2) * t729) * t156 + t733 * t160 + t734 * t159 + t735 * t137 + t736 * t161 - t465 * t551 - t222 * t507 + t775 * t434 / 0.2e1 + (-m(4) * t535 - m(5) * t475 - t679 * (t238 * pkin(4) + qJ(5) * t237 + t475) - t771 * t365 + t747 * t364 + t696 * t238 + t699 * t237 - t701 * t598) * g(1) + (-m(4) * t536 - m(5) * t476 - t679 * (t236 * pkin(4) + qJ(5) * t235 + t476) - t771 * t363 + t747 * t362 + t696 * t236 + t699 * t235 - t701 * t600) * g(2) + (-m(5) * t523 - m(4) * t586 - (t431 * mrSges(4,3) + t435 * t507) * t427 + t366 - t679 * (t308 * pkin(4) + qJ(5) * t307 + t523) + t696 * t308 + t699 * t307 + t701 * t563) * g(3) + (t12 * t284 + t13 * t285 + t16 * t332 + t71 * t736 + t72 * t739 + t735 * t80) * m(6) + t739 * t158 + t740 * t541 + (t11 * t349 + t17 * t720 - t58 * t741) * mrSges(7,1) + (t11 * t350 - t18 * t720 + t58 * t742) * mrSges(7,2) + (-t1 * t349 - t17 * t742 + t18 * t741 - t2 * t350) * mrSges(7,3) + t560 + (t305 * t496 + t306 * t501) * qJD(3) / 0.2e1 + t698 * t593 + (t526 - t341) * t343 + t332 * t47 + t318 * t62 + t319 * t64 + (-t568 - t231) * t248 - t292 * t115 / 0.2e1 + t284 * t61 + t285 * t63 + t281 * t10 - t244 * mrSges(3,2) + t245 * mrSges(3,1) + t146 * t21 + t147 * t22 + (-t569 - t232) * t247 - pkin(2) * t102 + t168 * t682 + t196 * t683 + t169 * t684 + t195 * t685 + t743 * (t429 * t538 + t433 * t541 - t293 / 0.2e1) + t744 * (t430 * t515 + t581 / 0.2e1) + t751 * (-t581 / 0.2e1 + t520 / 0.2e1) + t350 * t691 - t349 * t692 + ((t583 * t590 - t581) * t188 + (t583 * t588 - t580) * t187 + t714) * mrSges(4,3) + (-t187 * t231 - t188 * t232 - pkin(2) * t222 + ((-t187 * t434 - t188 * t430) * qJD(3) + t714) * pkin(9)) * m(4) + t754 * t73; t11 * t781 + (-m(6) * (t533 - t731) - m(5) * t533 - m(7) * (-t261 - t731) + t746 * t273 + t704 * t272) * g(1) + (-Ifges(7,5) * t671 - Ifges(4,2) * t652 + Ifges(5,6) * t663 + Ifges(6,6) * t664 - Ifges(7,6) * t673 - Ifges(7,3) * t657 + t766 * t654 + t768 * t662 + t793) * t306 + (-t567 - t109) * t159 + (Ifges(7,4) * t194 - Ifges(7,2) * t193) * t673 + (-Ifges(7,5) * t480 - Ifges(7,6) * t479) * t668 + (-t1 * t479 + t17 * t737 + t18 * t738 + t2 * t480) * mrSges(7,3) - t480 * t691 + (Ifges(7,1) * t257 - Ifges(7,4) * t258) * t670 + (Ifges(7,4) * t257 - Ifges(7,2) * t258) * t672 + (Ifges(7,5) * t257 - Ifges(7,6) * t258) * t656 + (t243 * t721 + t297 * t723) * qJD(4) / 0.2e1 + (pkin(3) * t468 - t108 * t81 - t109 * t82) * m(5) + t468 * t506 - t479 * t692 + (-Ifges(7,4) * t480 - Ifges(7,2) * t479) * t688 + (-Ifges(7,1) * t480 - Ifges(7,4) * t479) * t689 + t786 * t433 + (Ifges(7,5) * t194 - Ifges(7,6) * t193) * t657 + (Ifges(7,1) * t194 - Ifges(7,4) * t193) * t671 + t578 * t674 + t487 * t675 + t494 * t676 + t594 * t651 + (t16 * t385 - t71 * t90 - t72 * t88 + t727 * t80) * m(6) + ((t64 + t61) * t433 + (t63 - t62) * t429 + ((-t429 * t72 + t433 * t71) * qJD(4) + t712) * m(6) + ((-t429 * t82 - t433 * t81) * qJD(4) + t713) * m(5)) * pkin(10) + (-m(5) * t156 - t138 + t248 + t611) * t188 + (-t567 - t88) * t158 + (t566 - t90) * t161 + t724 * t667 + (t548 + t725) * qJD(4) + t727 * t137 + (-m(6) * (t532 + t728) - m(5) * t532 - m(7) * (-t354 + t728) - t519 * t361 + t508 + t708 * t360) * g(3) + (-t566 - t108) * t160 + t429 * t791 + (t495 * t663 + t488 * t664 + Ifges(4,1) * t650 - t279 * mrSges(4,2) + t756 / 0.2e1 + t721 * t662 + t723 * t654 - t725) * t305 + (t612 - t247) * t187 + (-mrSges(7,1) * t738 - mrSges(7,2) * t737) * t58 + t561 + (t296 + t740) * t652 + (-t626 + t744) * t650 - t16 * t504 + (-t495 / 0.2e1 + t488 / 0.2e1) * t579 + t385 * t47 + t372 * t10 + t282 * t21 + t283 * t22 - t78 * mrSges(4,2) + t79 * mrSges(4,1) + t257 * t682 + t194 * t683 - t258 * t684 - t193 * t685 + t743 * (-t602 / 0.2e1 + t537) + t722 * t680 - pkin(3) * t48 + (t71 * t716 + t715 * t72 + t697 + t712) * mrSges(6,2) + (t715 * t82 - t716 * t81 + t697 + t713) * mrSges(5,3) + (-m(6) * (t534 + t732) - m(5) * t534 - m(7) * (t259 + t732) + t746 * t269 - t704 * t529) * g(2) + t750 * t73 + t751 * t649 + t752 * t107 + t753 * t106 + (t1 * t283 + t11 * t372 + t17 * t752 + t18 * t753 + t2 * t282 + t58 * t750) * m(7); (-t242 * t768 - t243 * t802) * t654 + (t256 * t679 - t266 * t696 + t267 * t442) * g(3) + t796 * t673 + t778 - t132 * t683 + t118 * t661 + (Ifges(6,3) * t243 - t614) * t664 - t792 + (-t242 * t770 + t115 + t239 - t613) * t662 + (t209 * t726 + t210 * t442) * g(2) + (t213 * t726 + t214 * t442) * g(1) + t765 + (t242 * t71 + t243 * t72) * mrSges(6,2) + (t160 - t161 + t629) * t82 + (-t158 - t159 - t630) * t81 + (-Ifges(5,2) * t243 - t240 + t743) * t663 + t481 * t685 + t381 * t21 + t382 * t22 - t156 * (mrSges(5,1) * t243 - mrSges(5,2) * t242) - t80 * (mrSges(6,1) * t243 + mrSges(6,3) * t242) + qJD(5) * t158 - t136 * t137 - t101 * t73 + qJ(5) * t61 - pkin(4) * t63 + (-pkin(4) * t13 + qJ(5) * t12 - t136 * t80 - t71 * t82 + t711 * t72) * m(6) + t748 * t107 + t749 * t106 + (t1 * t382 - t101 * t58 + t17 * t748 + t18 * t749 + t2 * t381) * m(7); t432 * t21 + t428 * t22 + (t137 - t73) * t243 + t482 * qJD(6) + (-t158 - t482) * t297 + t63 + (t1 * t428 + t2 * t432 - t243 * t58 + t459 + t288 * (-t17 * t428 + t18 * t432)) * m(7) + (t243 * t80 - t297 * t72 + t13 + t459) * m(6); t56 * t670 - t17 * t106 + t18 * t107 - g(1) * (mrSges(7,1) * t112 - mrSges(7,2) * t113) - g(2) * ((t209 * t432 - t210 * t428) * mrSges(7,1) + (-t209 * t428 - t210 * t432) * mrSges(7,2)) - g(3) * t502 + (t57 - t796) * t673 + t792;];
tau  = t6;
