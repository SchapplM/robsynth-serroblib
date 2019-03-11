% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:06:24
% EndTime: 2019-03-09 22:07:18
% DurationCPUTime: 33.40s
% Computational Cost: add. (29483->971), mult. (63552->1257), div. (0->0), fcn. (47592->18), ass. (0->454)
t498 = sin(qJ(3));
t503 = cos(qJ(3));
t504 = cos(qJ(2));
t588 = qJD(1) * t504;
t499 = sin(qJ(2));
t589 = qJD(1) * t499;
t398 = -t498 * t589 + t503 * t588;
t421 = t498 * t504 + t499 * t503;
t399 = t421 * qJD(1);
t313 = pkin(3) * t399 - pkin(9) * t398;
t573 = pkin(2) * t589;
t293 = t313 + t573;
t506 = -pkin(8) - pkin(7);
t446 = t506 * t504;
t426 = qJD(1) * t446;
t400 = t498 * t426;
t445 = t506 * t499;
t425 = qJD(1) * t445;
t326 = t425 * t503 + t400;
t497 = sin(qJ(4));
t502 = cos(qJ(4));
t223 = t497 * t293 + t502 * t326;
t585 = qJD(3) * t503;
t571 = pkin(2) * t585;
t796 = t502 * t571 - t223;
t622 = t398 * t497;
t795 = -qJ(5) * t622 - t502 * qJD(5);
t222 = t502 * t293 - t326 * t497;
t486 = t502 * qJ(5);
t541 = t399 * pkin(4) - t398 * t486;
t662 = pkin(2) * t498;
t471 = pkin(9) + t662;
t599 = -qJ(5) - t471;
t545 = qJD(4) * t599;
t794 = -t222 - t541 + (-qJD(5) - t571) * t497 + t502 * t545;
t407 = qJD(2) * pkin(2) + t425;
t317 = t407 * t503 + t400;
t225 = t502 * t313 - t317 * t497;
t495 = -qJ(5) - pkin(9);
t553 = qJD(4) * t495;
t793 = -qJD(5) * t497 + t502 * t553 - t225 - t541;
t792 = -t497 * t545 + t795 - t796;
t226 = t497 * t313 + t502 * t317;
t791 = -t497 * t553 + t226 + t795;
t493 = sin(pkin(11));
t494 = cos(pkin(11));
t416 = t493 * t502 + t494 * t497;
t280 = t416 * t398;
t394 = t416 * qJD(4);
t790 = t280 - t394;
t526 = t493 * t497 - t494 * t502;
t281 = t526 * t398;
t395 = t526 * qJD(4);
t789 = t281 - t395;
t491 = qJ(4) + pkin(11);
t482 = qJ(6) + t491;
t465 = sin(t482);
t480 = sin(t491);
t788 = -mrSges(6,2) * t480 - mrSges(7,2) * t465;
t757 = Ifges(6,3) + Ifges(5,3);
t742 = t792 * t493 + t494 * t794;
t741 = t493 * t794 - t792 * t494;
t787 = t790 * pkin(10);
t740 = t493 * t791 + t494 * t793;
t739 = t493 * t793 - t494 * t791;
t466 = cos(t482);
t481 = cos(t491);
t785 = mrSges(6,1) * t481 + mrSges(7,1) * t466;
t784 = -t399 * pkin(5) - pkin(10) * t789;
t783 = -mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t782 = t784 + t742;
t781 = t787 + t741;
t780 = t784 + t740;
t779 = t739 + t787;
t600 = t503 * t426;
t325 = t425 * t498 - t600;
t586 = qJD(3) * t498;
t778 = pkin(2) * t586 - t325;
t496 = sin(qJ(6));
t501 = cos(qJ(6));
t391 = qJD(4) - t398;
t490 = qJD(2) + qJD(3);
t342 = -t399 * t497 + t490 * t502;
t343 = t399 * t502 + t490 * t497;
t247 = t342 * t493 + t343 * t494;
t772 = pkin(10) * t247;
t487 = t504 * pkin(2);
t474 = t487 + pkin(1);
t444 = t474 * qJD(1);
t284 = -pkin(3) * t398 - pkin(9) * t399 - t444;
t318 = t498 * t407 - t600;
t298 = pkin(9) * t490 + t318;
t202 = t502 * t284 - t298 * t497;
t158 = -qJ(5) * t343 + t202;
t140 = pkin(4) * t391 + t158;
t203 = t284 * t497 + t298 * t502;
t159 = qJ(5) * t342 + t203;
t150 = t493 * t159;
t87 = t494 * t140 - t150;
t59 = pkin(5) * t391 - t772 + t87;
t549 = t494 * t342 - t343 * t493;
t758 = pkin(10) * t549;
t612 = t494 * t159;
t88 = t493 * t140 + t612;
t61 = t88 + t758;
t20 = -t496 * t61 + t501 * t59;
t21 = t496 * t59 + t501 * t61;
t752 = t490 * Ifges(4,6);
t777 = t444 * mrSges(4,1) - t202 * mrSges(5,1) - t87 * mrSges(6,1) - t20 * mrSges(7,1) + t203 * mrSges(5,2) + t88 * mrSges(6,2) + t21 * mrSges(7,2) + t752 / 0.2e1;
t753 = t490 * Ifges(4,5);
t776 = t444 * mrSges(4,2) + t317 * mrSges(4,3) - t753 / 0.2e1;
t492 = qJ(2) + qJ(3);
t484 = sin(t492);
t640 = mrSges(5,2) * t497;
t775 = (-t640 + t788) * t484;
t774 = -m(5) * pkin(9) + t783;
t584 = qJD(4) * t497;
t477 = pkin(4) * t584;
t771 = -pkin(5) * t790 + t477;
t628 = t399 * mrSges(4,3);
t736 = mrSges(4,1) * t490 + mrSges(5,1) * t342 - mrSges(5,2) * t343 - t628;
t368 = pkin(4) * t622;
t770 = -t368 + t778;
t582 = qJD(1) * qJD(2);
t430 = qJDD(1) * t504 - t499 * t582;
t420 = t498 * t499 - t503 * t504;
t516 = t420 * qJD(3);
t330 = -qJD(2) * t420 - t516;
t583 = qJD(4) * t502;
t561 = t421 * t583;
t519 = t330 * t497 + t561;
t431 = qJDD(1) * t499 + t504 * t582;
t273 = -qJD(1) * t516 + t430 * t498 + t431 * t503;
t274 = -qJD(3) * t399 + t430 * t503 - t498 * t431;
t624 = qJDD(1) * pkin(1);
t385 = -pkin(2) * t430 - t624;
t160 = -pkin(3) * t274 - pkin(9) * t273 + t385;
t414 = t431 * pkin(7);
t334 = qJDD(2) * pkin(2) - pkin(8) * t431 - t414;
t413 = t430 * pkin(7);
t341 = pkin(8) * t430 + t413;
t181 = t498 * t334 + t503 * t341 + t407 * t585 + t426 * t586;
t488 = qJDD(2) + qJDD(3);
t174 = pkin(9) * t488 + t181;
t64 = t497 * t160 + t502 * t174 + t284 * t583 - t298 * t584;
t65 = -qJD(4) * t203 + t502 * t160 - t174 * t497;
t769 = -t497 * t65 + t502 * t64;
t644 = mrSges(5,1) * t502;
t768 = t640 - t644;
t146 = t247 * t501 + t496 * t549;
t371 = qJD(6) + t391;
t754 = t342 * Ifges(5,6);
t764 = -t247 * t496 + t501 * t549;
t767 = t343 * Ifges(5,5) + t247 * Ifges(6,5) + t146 * Ifges(7,5) + Ifges(6,6) * t549 + Ifges(7,6) * t764 + t371 * Ifges(7,3) + t391 * t757 + t754;
t766 = (t644 + t785) * t484;
t485 = cos(t492);
t763 = -t485 * mrSges(4,1) + (mrSges(4,2) + t783) * t484;
t730 = t485 * pkin(3) + t484 * pkin(9);
t654 = pkin(4) * t502;
t472 = pkin(3) + t654;
t732 = t485 * t472 - t484 * t495;
t435 = pkin(5) * t481 + t654;
t424 = pkin(3) + t435;
t489 = -pkin(10) + t495;
t735 = t485 * t424 - t484 * t489;
t762 = -m(5) * t730 - m(6) * t732 - m(7) * t735;
t710 = m(6) * pkin(4);
t206 = qJD(4) * t342 + t273 * t502 + t488 * t497;
t207 = -qJD(4) * t343 - t273 * t497 + t488 * t502;
t122 = -t206 * t493 + t207 * t494;
t123 = t206 * t494 + t207 * t493;
t37 = qJD(6) * t764 + t122 * t496 + t123 * t501;
t707 = t37 / 0.2e1;
t38 = -qJD(6) * t146 + t122 * t501 - t123 * t496;
t706 = t38 / 0.2e1;
t760 = m(6) + m(7);
t698 = t122 / 0.2e1;
t697 = t123 / 0.2e1;
t688 = t206 / 0.2e1;
t687 = t207 / 0.2e1;
t272 = qJDD(4) - t274;
t267 = qJDD(6) + t272;
t682 = t267 / 0.2e1;
t681 = t272 / 0.2e1;
t759 = t430 / 0.2e1;
t667 = t504 / 0.2e1;
t409 = t599 * t497;
t410 = t471 * t502 + t486;
t309 = t494 * t409 - t410 * t493;
t649 = pkin(10) * t416;
t263 = t309 - t649;
t310 = t493 * t409 + t494 * t410;
t408 = t526 * pkin(10);
t264 = -t408 + t310;
t171 = t263 * t501 - t264 * t496;
t756 = qJD(6) * t171 + t496 * t782 + t501 * t781;
t172 = t263 * t496 + t264 * t501;
t755 = -qJD(6) * t172 - t496 * t781 + t501 * t782;
t751 = t504 * Ifges(3,2);
t441 = t495 * t497;
t443 = pkin(9) * t502 + t486;
t339 = t494 * t441 - t443 * t493;
t288 = t339 - t649;
t340 = t493 * t441 + t494 * t443;
t289 = -t408 + t340;
t199 = t288 * t501 - t289 * t496;
t750 = qJD(6) * t199 + t496 * t780 + t501 * t779;
t200 = t288 * t496 + t289 * t501;
t749 = -qJD(6) * t200 - t496 * t779 + t501 * t780;
t656 = pkin(4) * t494;
t467 = pkin(5) + t656;
t657 = pkin(4) * t493;
t379 = t467 * t501 - t496 * t657;
t91 = -t158 * t493 - t612;
t69 = t91 - t758;
t92 = t494 * t158 - t150;
t70 = t92 - t772;
t748 = t379 * qJD(6) - t496 * t69 - t501 * t70;
t380 = t467 * t496 + t501 * t657;
t747 = -t380 * qJD(6) + t496 * t70 - t501 * t69;
t746 = t710 + mrSges(5,1);
t297 = -pkin(3) * t490 - t317;
t537 = mrSges(5,1) * t497 + mrSges(5,2) * t502;
t745 = t297 * t537;
t738 = t770 + t771;
t316 = pkin(3) * t420 - pkin(9) * t421 - t474;
t349 = t445 * t498 - t446 * t503;
t337 = t502 * t349;
t239 = t497 * t316 + t337;
t258 = t368 + t318;
t737 = -t258 + t771;
t734 = t503 * t445 + t446 * t498;
t733 = t477 + t770;
t729 = t477 - t258;
t727 = Ifges(5,5) * t206 + Ifges(6,5) * t123 + Ifges(5,6) * t207 + Ifges(6,6) * t122 + t272 * t757;
t651 = pkin(7) * t504;
t652 = pkin(7) * t499;
t726 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t589) * t651 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t588) * t652;
t725 = t413 * t504 + t414 * t499;
t500 = sin(qJ(1));
t505 = cos(qJ(1));
t724 = g(1) * t505 + g(2) * t500;
t723 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t722 = m(5) + m(4) + t760;
t527 = -t472 * t484 - t485 * t495;
t529 = -t424 * t484 - t485 * t489;
t659 = pkin(3) * t484;
t661 = pkin(2) * t499;
t721 = -m(7) * (t529 - t661) - m(6) * (t527 - t661) - m(5) * (-t659 - t661) + t766;
t720 = t763 + (t768 - t785 - t788) * t485;
t41 = pkin(4) * t272 - qJ(5) * t206 - qJD(5) * t343 + t65;
t46 = qJ(5) * t207 + qJD(5) * t342 + t64;
t15 = t494 * t41 - t46 * t493;
t7 = pkin(5) * t272 - pkin(10) * t123 + t15;
t16 = t493 * t41 + t494 * t46;
t8 = pkin(10) * t122 + t16;
t3 = qJD(6) * t20 + t496 * t7 + t501 * t8;
t4 = -qJD(6) * t21 - t496 * t8 + t501 * t7;
t719 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t655 = pkin(4) * t497;
t434 = pkin(5) * t480 + t655;
t718 = -m(6) * t655 - m(7) * t434;
t614 = t485 * t500;
t717 = t500 * t775 + t614 * t774;
t613 = t485 * t505;
t716 = t505 * t775 + t613 * t774;
t715 = m(5) * t659 - m(6) * t527 - m(7) * t529 + t766;
t442 = -mrSges(3,1) * t504 + mrSges(3,2) * t499;
t714 = m(3) * pkin(1) + mrSges(2,1) - t442 - t763;
t136 = mrSges(5,1) * t272 - mrSges(5,3) * t206;
t137 = -mrSges(5,2) * t272 + mrSges(5,3) * t207;
t637 = mrSges(5,3) * t342;
t276 = -mrSges(5,2) * t391 + t637;
t636 = mrSges(5,3) * t343;
t277 = mrSges(5,1) * t391 - t636;
t713 = m(5) * ((-t202 * t502 - t203 * t497) * qJD(4) + t769) - t277 * t583 - t276 * t584 + t502 * t137 - t497 * t136;
t712 = t65 * mrSges(5,1) + t15 * mrSges(6,1) - t64 * mrSges(5,2) - t16 * mrSges(6,2);
t709 = Ifges(7,4) * t707 + Ifges(7,2) * t706 + Ifges(7,6) * t682;
t708 = Ifges(7,1) * t707 + Ifges(7,4) * t706 + Ifges(7,5) * t682;
t705 = Ifges(6,4) * t697 + Ifges(6,2) * t698 + Ifges(6,6) * t681;
t704 = Ifges(6,1) * t697 + Ifges(6,4) * t698 + Ifges(6,5) * t681;
t630 = Ifges(7,4) * t146;
t79 = Ifges(7,2) * t764 + Ifges(7,6) * t371 + t630;
t703 = -t79 / 0.2e1;
t702 = t79 / 0.2e1;
t142 = Ifges(7,4) * t764;
t80 = Ifges(7,1) * t146 + Ifges(7,5) * t371 + t142;
t701 = -t80 / 0.2e1;
t700 = t80 / 0.2e1;
t699 = Ifges(5,1) * t688 + Ifges(5,4) * t687 + Ifges(5,5) * t681;
t134 = Ifges(6,4) * t247 + Ifges(6,2) * t549 + Ifges(6,6) * t391;
t696 = -t134 / 0.2e1;
t695 = t134 / 0.2e1;
t135 = Ifges(6,1) * t247 + Ifges(6,4) * t549 + Ifges(6,5) * t391;
t694 = -t135 / 0.2e1;
t693 = t135 / 0.2e1;
t692 = -t764 / 0.2e1;
t691 = t764 / 0.2e1;
t690 = -t146 / 0.2e1;
t689 = t146 / 0.2e1;
t686 = -t549 / 0.2e1;
t685 = t549 / 0.2e1;
t684 = -t247 / 0.2e1;
t683 = t247 / 0.2e1;
t679 = -t342 / 0.2e1;
t678 = -t343 / 0.2e1;
t677 = t343 / 0.2e1;
t676 = -t371 / 0.2e1;
t675 = t371 / 0.2e1;
t674 = -t391 / 0.2e1;
t673 = t391 / 0.2e1;
t671 = t398 / 0.2e1;
t669 = t399 / 0.2e1;
t666 = mrSges(6,3) * t87;
t665 = mrSges(6,3) * t88;
t664 = mrSges(7,3) * t20;
t663 = mrSges(7,3) * t21;
t660 = pkin(2) * t503;
t658 = pkin(4) * t343;
t646 = g(3) * t484;
t331 = t490 * t421;
t525 = -qJ(5) * t330 - qJD(5) * t421;
t587 = qJD(2) * t499;
t572 = pkin(2) * t587;
t234 = pkin(3) * t331 - pkin(9) * t330 + t572;
t562 = qJD(2) * t506;
t428 = t499 * t562;
t429 = t504 * t562;
t250 = qJD(3) * t734 + t428 * t503 + t429 * t498;
t550 = t502 * t234 - t250 * t497;
t68 = pkin(4) * t331 + t525 * t502 + (-t337 + (qJ(5) * t421 - t316) * t497) * qJD(4) + t550;
t563 = t497 * t234 + t502 * t250 + t316 * t583;
t84 = -qJ(5) * t561 + (-qJD(4) * t349 + t525) * t497 + t563;
t25 = t493 * t68 + t494 * t84;
t635 = Ifges(3,4) * t499;
t634 = Ifges(3,4) * t504;
t633 = Ifges(5,4) * t343;
t632 = Ifges(5,4) * t497;
t631 = Ifges(5,4) * t502;
t627 = t399 * Ifges(4,4);
t619 = t421 * t497;
t618 = t421 * t502;
t230 = Ifges(5,2) * t342 + Ifges(5,6) * t391 + t633;
t610 = t497 * t230;
t609 = t497 * t500;
t608 = t497 * t505;
t605 = t500 * t434;
t604 = t500 * t502;
t338 = Ifges(5,4) * t342;
t231 = t343 * Ifges(5,1) + t391 * Ifges(5,5) + t338;
t602 = t502 * t231;
t601 = t502 * t505;
t238 = t502 * t316 - t349 * t497;
t190 = pkin(4) * t420 - t421 * t486 + t238;
t215 = -qJ(5) * t619 + t239;
t116 = t493 * t190 + t494 * t215;
t351 = t465 * t614 + t466 * t505;
t352 = t465 * t505 - t466 * t614;
t598 = -t351 * mrSges(7,1) + t352 * mrSges(7,2);
t353 = -t465 * t613 + t466 * t500;
t354 = t465 * t500 + t466 * t613;
t597 = t353 * mrSges(7,1) - t354 * mrSges(7,2);
t577 = Ifges(7,5) * t37 + Ifges(7,6) * t38 + Ifges(7,3) * t267;
t558 = t602 / 0.2e1;
t13 = -t38 * mrSges(7,1) + t37 * mrSges(7,2);
t555 = -t584 / 0.2e1;
t24 = -t493 * t84 + t494 * t68;
t552 = t582 / 0.2e1;
t54 = -t122 * mrSges(6,1) + t123 * mrSges(6,2);
t115 = t494 * t190 - t215 * t493;
t285 = pkin(4) * t619 - t734;
t540 = mrSges(3,1) * t499 + mrSges(3,2) * t504;
t538 = mrSges(4,1) * t484 + mrSges(4,2) * t485;
t536 = -mrSges(7,1) * t465 - mrSges(7,2) * t466;
t535 = Ifges(5,1) * t502 - t632;
t534 = t635 + t751;
t533 = -Ifges(5,2) * t497 + t631;
t532 = Ifges(3,5) * t504 - Ifges(3,6) * t499;
t531 = Ifges(5,5) * t502 - Ifges(5,6) * t497;
t300 = t416 * t421;
t100 = -pkin(10) * t300 + t116;
t301 = t526 * t421;
t90 = pkin(5) * t420 + pkin(10) * t301 + t115;
t43 = t100 * t501 + t496 * t90;
t42 = -t100 * t496 + t501 * t90;
t213 = -t300 * t501 + t301 * t496;
t214 = -t300 * t496 - t301 * t501;
t319 = -t416 * t496 - t501 * t526;
t320 = t416 * t501 - t496 * t526;
t362 = pkin(5) * t526 - t472;
t182 = t334 * t503 - t498 * t341 - t407 * t586 + t426 * t585;
t521 = t577 + t719;
t520 = pkin(1) * t540;
t383 = -t485 * t608 + t604;
t381 = t485 * t609 + t601;
t518 = -t330 * t502 + t421 * t584;
t517 = t499 * (Ifges(3,1) * t504 - t635);
t175 = -pkin(3) * t488 - t182;
t235 = -pkin(4) * t342 + qJD(5) + t297;
t251 = qJD(3) * t349 + t428 * t498 - t503 * t429;
t113 = -pkin(4) * t207 + qJDD(5) + t175;
t170 = pkin(4) * t519 + t251;
t141 = -pkin(5) * t549 + t235;
t193 = -t280 * t501 + t281 * t496;
t194 = -t280 * t496 - t281 * t501;
t227 = qJD(6) * t319 - t394 * t496 - t395 * t501;
t228 = -qJD(6) * t320 - t394 * t501 + t395 * t496;
t294 = t398 * Ifges(4,2) + t627 + t752;
t386 = Ifges(4,4) * t398;
t295 = t399 * Ifges(4,1) + t386 + t753;
t51 = -pkin(5) * t122 + t113;
t98 = Ifges(5,4) * t206 + Ifges(5,2) * t207 + Ifges(5,6) * t272;
t507 = (-Ifges(6,4) * t683 - Ifges(6,2) * t685 - Ifges(6,6) * t673 - t665 - t695) * t394 + (t531 * t674 + t533 * t679 + t535 * t678 - t745 + t776) * t398 + (Ifges(5,5) * t678 + Ifges(6,5) * t684 + Ifges(7,5) * t690 + Ifges(5,6) * t679 + Ifges(6,6) * t686 + Ifges(7,6) * t692 + Ifges(7,3) * t676 + t757 * t674 + t777) * t399 + (-Ifges(6,5) * t281 - Ifges(6,6) * t280) * t674 + (Ifges(7,5) * t194 + Ifges(7,6) * t193) * t676 + (t342 * t533 + t343 * t535 + t391 * t531) * qJD(4) / 0.2e1 - (-Ifges(4,2) * t399 + t295 + t386 + t602) * t398 / 0.2e1 + (Ifges(5,5) * t497 + Ifges(6,5) * t416 + Ifges(5,6) * t502 - Ifges(6,6) * t526) * t681 - t526 * t705 + (Ifges(6,1) * t416 - Ifges(6,4) * t526) * t697 + (Ifges(6,4) * t416 - Ifges(6,2) * t526) * t698 + t113 * (mrSges(6,1) * t526 + mrSges(6,2) * t416) + (-t15 * t416 - t16 * t526 + t280 * t88 - t281 * t87) * mrSges(6,3) - (Ifges(4,1) * t398 - t627 + t767) * t399 / 0.2e1 + t175 * t768 + ((-t584 + t622) * t203 + (t398 * t502 - t583) * t202 + t769) * mrSges(5,3) + (-Ifges(6,1) * t683 - Ifges(6,4) * t685 - Ifges(6,5) * t673 + t666 - t693) * t395 + (Ifges(7,4) * t689 + Ifges(7,2) * t691 + Ifges(7,6) * t675 + t663 + t702) * t228 + (Ifges(7,1) * t689 + Ifges(7,4) * t691 + Ifges(7,5) * t675 - t664 + t700) * t227 + (Ifges(7,5) * t320 + Ifges(7,6) * t319) * t682 + (Ifges(7,4) * t194 + Ifges(7,2) * t193) * t692 + t497 * t699 + t194 * t701 + t193 * t703 + t416 * t704 + (Ifges(7,4) * t320 + Ifges(7,2) * t319) * t706 + (Ifges(7,1) * t320 + Ifges(7,4) * t319) * t707 + t320 * t708 + t319 * t709 + t502 * t98 / 0.2e1 + (-Ifges(6,1) * t281 - Ifges(6,4) * t280) * t684 + ((-t194 + t227) * mrSges(7,2) + (t193 - t228) * mrSges(7,1)) * t141 + (-Ifges(6,4) * t281 - Ifges(6,2) * t280) * t686 + t230 * t555 + (Ifges(5,2) * t502 + t632) * t687 + (Ifges(5,1) * t497 + t631) * t688 - t281 * t694 - t280 * t696 + (t745 + t558) * qJD(4) + Ifges(4,3) * t488 + (Ifges(7,1) * t194 + Ifges(7,4) * t193) * t690 + t294 * t669 + t610 * t671 + t318 * t628 + t51 * (-mrSges(7,1) * t319 + mrSges(7,2) * t320) + Ifges(4,5) * t273 + Ifges(4,6) * t274 - t181 * mrSges(4,2) + t182 * mrSges(4,1) + (-t193 * t21 + t194 * t20 + t3 * t319 - t320 * t4) * mrSges(7,3) + (-mrSges(6,1) * t790 + mrSges(6,2) * t789) * t235;
t476 = Ifges(3,4) * t588;
t473 = -pkin(3) - t660;
t438 = -t472 - t660;
t397 = Ifges(3,1) * t589 + Ifges(3,5) * qJD(2) + t476;
t396 = Ifges(3,6) * qJD(2) + qJD(1) * t534;
t384 = t485 * t601 + t609;
t382 = -t485 * t604 + t608;
t361 = t480 * t500 + t481 * t613;
t360 = -t480 * t613 + t481 * t500;
t359 = t480 * t505 - t481 * t614;
t358 = t480 * t614 + t481 * t505;
t356 = -mrSges(4,2) * t490 + mrSges(4,3) * t398;
t350 = t362 - t660;
t312 = -mrSges(4,1) * t398 + mrSges(4,2) * t399;
t256 = -mrSges(4,2) * t488 + mrSges(4,3) * t274;
t255 = mrSges(4,1) * t488 - mrSges(4,3) * t273;
t219 = mrSges(6,1) * t391 - mrSges(6,3) * t247;
t218 = -mrSges(6,2) * t391 + mrSges(6,3) * t549;
t217 = pkin(5) * t300 + t285;
t201 = pkin(5) * t247 + t658;
t168 = -t330 * t526 - t394 * t421;
t167 = -t330 * t416 + t395 * t421;
t149 = -mrSges(6,1) * t549 + mrSges(6,2) * t247;
t132 = mrSges(7,1) * t371 - mrSges(7,3) * t146;
t131 = -mrSges(7,2) * t371 + mrSges(7,3) * t764;
t124 = -mrSges(5,1) * t207 + mrSges(5,2) * t206;
t112 = -qJD(4) * t239 + t550;
t111 = -t349 * t584 + t563;
t102 = -pkin(5) * t167 + t170;
t94 = mrSges(6,1) * t272 - mrSges(6,3) * t123;
t93 = -mrSges(6,2) * t272 + mrSges(6,3) * t122;
t89 = -mrSges(7,1) * t764 + mrSges(7,2) * t146;
t63 = -qJD(6) * t214 + t167 * t501 - t168 * t496;
t62 = qJD(6) * t213 + t167 * t496 + t168 * t501;
t29 = -mrSges(7,2) * t267 + mrSges(7,3) * t38;
t28 = mrSges(7,1) * t267 - mrSges(7,3) * t37;
t19 = pkin(10) * t167 + t25;
t18 = pkin(5) * t331 - pkin(10) * t168 + t24;
t6 = -qJD(6) * t43 + t18 * t501 - t19 * t496;
t5 = qJD(6) * t42 + t18 * t496 + t19 * t501;
t1 = [(Ifges(3,4) * t431 + Ifges(3,2) * t430) * t667 + (t385 * mrSges(4,2) - t182 * mrSges(4,3) + Ifges(4,1) * t273 + Ifges(4,4) * t274 + Ifges(4,5) * t488 + t175 * t537 + t231 * t555 + t531 * t681 + t533 * t687 + t535 * t688) * t421 + (-Ifges(5,1) * t518 - Ifges(5,4) * t519) * t677 + t342 * (-Ifges(5,4) * t518 - Ifges(5,2) * t519) / 0.2e1 + (-t609 * t710 - m(7) * t605 - t384 * mrSges(5,1) - t361 * mrSges(6,1) - t354 * mrSges(7,1) - t383 * mrSges(5,2) - t360 * mrSges(6,2) - t353 * mrSges(7,2) - t722 * (t505 * t474 - t500 * t506) + t723 * t500 + (-t714 + t762) * t505) * g(2) + (t558 - t610 / 0.2e1 + Ifges(4,1) * t669 + Ifges(4,4) * t671 + t295 / 0.2e1 - t776) * t330 + (t767 / 0.2e1 + Ifges(7,3) * t675 + Ifges(5,5) * t677 + Ifges(6,5) * t683 + Ifges(6,6) * t685 + Ifges(7,5) * t689 + Ifges(7,6) * t691 - Ifges(4,4) * t669 - Ifges(4,2) * t671 - t294 / 0.2e1 + t754 / 0.2e1 - t318 * mrSges(4,3) + t757 * t673 - t777) * t331 + m(4) * (t181 * t349 + t250 * t318 - t385 * t474 - t444 * t572) + (-t382 * mrSges(5,1) - t359 * mrSges(6,1) - t352 * mrSges(7,1) - t381 * mrSges(5,2) - t358 * mrSges(6,2) - t351 * mrSges(7,2) + (t506 * t722 + t718 + t723) * t505 + (-m(7) * (-t474 - t735) - m(5) * (-t474 - t730) - m(6) * (-t474 - t732) + m(4) * t474 + t714) * t500) * g(1) + (t504 * t634 + t517) * t552 + (Ifges(7,4) * t214 + Ifges(7,2) * t213) * t706 + (Ifges(7,4) * t62 + Ifges(7,2) * t63) * t691 + (Ifges(7,5) * t62 + Ifges(7,6) * t63) * t675 + (Ifges(7,5) * t214 + Ifges(7,6) * t213) * t682 + (-Ifges(6,5) * t301 - Ifges(6,6) * t300) * t681 + (t15 * t301 - t16 * t300 + t167 * t88 - t168 * t87) * mrSges(6,3) + (-Ifges(6,1) * t301 - Ifges(6,4) * t300) * t697 + (-Ifges(6,4) * t301 - Ifges(6,2) * t300) * t698 + t113 * (mrSges(6,1) * t300 - mrSges(6,2) * t301) + (-Ifges(5,5) * t518 + Ifges(6,5) * t168 - Ifges(5,6) * t519 + Ifges(6,6) * t167) * t673 + m(7) * (t102 * t141 + t20 * t6 + t21 * t5 + t217 * t51 + t3 * t43 + t4 * t42) + m(6) * (t113 * t285 + t115 * t15 + t116 * t16 + t170 * t235 + t24 * t87 + t25 * t88) - t520 * t582 + (-t20 * t62 + t21 * t63 + t213 * t3 - t214 * t4) * mrSges(7,3) + t312 * t572 + (Ifges(7,1) * t214 + Ifges(7,4) * t213) * t707 + t618 * t699 + t62 * t700 + t63 * t702 - t301 * t704 - t300 * t705 + t214 * t708 + t213 * t709 - (-m(4) * t182 + m(5) * t175 + t124 - t255) * t734 + m(5) * (t111 * t203 + t112 * t202 + t238 * t65 + t239 * t64) - t396 * t587 / 0.2e1 + t431 * t634 / 0.2e1 + (Ifges(6,1) * t168 + Ifges(6,4) * t167) * t683 - t230 * t561 / 0.2e1 + (Ifges(7,1) * t62 + Ifges(7,4) * t63) * t689 - t442 * t624 + (t202 * t518 - t203 * t519 - t618 * t65 - t619 * t64) * mrSges(5,3) - t98 * t619 / 0.2e1 + t168 * t693 + t167 * t695 + (Ifges(6,4) * t168 + Ifges(6,2) * t167) * t685 - t474 * (-mrSges(4,1) * t274 + mrSges(4,2) * t273) + Ifges(2,3) * qJDD(1) - pkin(1) * (-mrSges(3,1) * t430 + mrSges(3,2) * t431) + t250 * t356 + t349 * t256 + t285 * t54 + t111 * t276 + t112 * t277 + t238 * t136 + t239 * t137 + t235 * (-mrSges(6,1) * t167 + mrSges(6,2) * t168) + t217 * t13 + t25 * t218 + t24 * t219 + t51 * (-mrSges(7,1) * t213 + mrSges(7,2) * t214) + t297 * (mrSges(5,1) * t519 - mrSges(5,2) * t518) + t170 * t149 + t141 * (-mrSges(7,1) * t63 + mrSges(7,2) * t62) + t5 * t131 + t6 * t132 + t115 * t94 + t116 * t93 + t102 * t89 + t534 * t759 + (t430 * t651 + t431 * t652 + t725) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t725) + (t397 * t667 + t532 * qJD(2) / 0.2e1 - t726) * qJD(2) + (t577 + t727) * t420 / 0.2e1 + (-m(4) * t317 + m(5) * t297 - t736) * t251 + (t385 * mrSges(4,1) - t181 * mrSges(4,3) - Ifges(4,4) * t273 + Ifges(5,5) * t688 + Ifges(6,5) * t697 + Ifges(7,5) * t707 - Ifges(4,2) * t274 - Ifges(4,6) * t488 + Ifges(5,6) * t687 + Ifges(6,6) * t698 + Ifges(7,6) * t706 + Ifges(7,3) * t682 + t681 * t757 + t712 + t719) * t420 + t42 * t28 + t43 * t29 + (-mrSges(3,1) * t652 - mrSges(3,2) * t651 + 0.2e1 * Ifges(3,6) * t667) * qJDD(2) + (Ifges(3,1) * t431 + Ifges(3,4) * t759 + Ifges(3,5) * qJDD(2) - t552 * t751) * t499; -(-Ifges(3,2) * t589 + t397 + t476) * t588 / 0.2e1 + (m(4) * t661 + t538 + t540) * t724 + ((t181 * t498 + t182 * t503 + (-t317 * t498 + t318 * t503) * qJD(3)) * pkin(2) + t317 * t325 - t318 * t326 + t444 * t573) * m(4) - t736 * t778 + (-m(6) * (t487 + t732) - m(7) * (t487 + t735) - m(5) * (t487 + t730) + t442 - m(4) * t487 + t720) * g(3) + t713 * t471 - t532 * t582 / 0.2e1 + (t505 * t721 + t716) * g(1) + (t500 * t721 + t717) * g(2) + (-t202 * t222 - t203 * t223 - t297 * t325 + t175 * t473 + (t297 * t498 + (-t202 * t497 + t203 * t502) * t503) * qJD(3) * pkin(2)) * m(5) + t396 * t589 / 0.2e1 + t507 + (t726 + (-t517 / 0.2e1 + t520) * qJD(1)) * qJD(1) + t473 * t124 + t438 * t54 + Ifges(3,6) * t430 + Ifges(3,5) * t431 - t413 * mrSges(3,2) - t414 * mrSges(3,1) + (t571 - t326) * t356 + t255 * t660 + t256 * t662 + t350 * t13 + t309 * t94 + t310 * t93 + Ifges(3,3) * qJDD(2) + t171 * t28 + t172 * t29 + (-t497 * t571 - t222) * t277 + t733 * t149 + t738 * t89 + t741 * t218 + t742 * t219 + (t113 * t438 + t15 * t309 + t16 * t310 + t235 * t733 + t741 * t88 + t742 * t87) * m(6) + t796 * t276 + t755 * t132 + t756 * t131 + (t141 * t738 + t171 * t4 + t172 * t3 + t20 * t755 + t21 * t756 + t350 * t51) * m(7) - t312 * t573; (-pkin(3) * t175 - t202 * t225 - t203 * t226 - t297 * t318) * m(5) + (t720 + t762) * g(3) + (t505 * t715 + t716) * g(1) + (t500 * t715 + t717) * g(2) + t713 * pkin(9) + t507 - t472 * t54 + t362 * t13 - t317 * t356 + t339 * t94 + t340 * t93 - t226 * t276 - t225 * t277 + t199 * t28 + t200 * t29 - pkin(3) * t124 + t724 * t538 + t729 * t149 + t736 * t318 + t737 * t89 + t739 * t218 + t740 * t219 + (-t113 * t472 + t15 * t339 + t16 * t340 + t235 * t729 + t739 * t88 + t740 * t87) * m(6) + t749 * t132 + t750 * t131 + (t141 * t737 + t199 * t4 + t20 * t749 + t200 * t3 + t21 * t750 + t362 * t51) * m(7); (-mrSges(6,2) * t235 + Ifges(6,1) * t684 + Ifges(6,4) * t686 + Ifges(6,5) * t674 + t666 + t694) * t549 - (mrSges(6,1) * t235 + Ifges(6,4) * t684 + Ifges(6,2) * t686 + Ifges(6,6) * t674 - t665 + t696) * t247 + t712 + (mrSges(6,1) * t480 + mrSges(6,2) * t481 - t536 + t537 - t718) * t646 + (Ifges(5,5) * t342 - Ifges(5,6) * t343) * t674 + t230 * t677 + (Ifges(5,1) * t342 - t633) * t678 + t727 + (t15 * t494 + t16 * t493) * t710 + (-mrSges(7,2) * t141 + Ifges(7,1) * t690 + Ifges(7,4) * t692 + Ifges(7,5) * t676 + t664 + t701) * t764 - (mrSges(7,1) * t141 + Ifges(7,4) * t690 + Ifges(7,2) * t692 + Ifges(7,6) * t676 - t663 + t703) * t146 + (-Ifges(5,2) * t343 + t231 + t338) * t679 - t149 * t658 - m(6) * (t235 * t658 + t87 * t91 + t88 * t92) + t521 + t93 * t657 + t94 * t656 + t379 * t28 + t380 * t29 - t297 * (mrSges(5,1) * t343 + mrSges(5,2) * t342) + (t636 + t277) * t203 + (t637 - t276) * t202 - t92 * t218 - t91 * t219 - t201 * t89 + (t358 * mrSges(6,1) - t359 * mrSges(6,2) - m(7) * (-t435 * t505 - t485 * t605) - t598 - mrSges(5,2) * t382 + t746 * t381) * g(2) + (-t360 * mrSges(6,1) + t361 * mrSges(6,2) - m(7) * (-t434 * t613 + t435 * t500) - t597 + mrSges(5,2) * t384 - t746 * t383) * g(1) + t747 * t132 + t748 * t131 + (-t141 * t201 + t747 * t20 + t748 * t21 + t3 * t380 + t379 * t4) * m(7); -t764 * t131 + t146 * t132 - t549 * t218 + t247 * t219 + t13 + t54 + (g(3) * t485 - t484 * t724) * t760 + (t146 * t20 - t21 * t764 + t51) * m(7) + (t247 * t87 - t549 * t88 + t113) * m(6); -t141 * (mrSges(7,1) * t146 + mrSges(7,2) * t764) + (Ifges(7,1) * t764 - t630) * t690 + t79 * t689 + (Ifges(7,5) * t764 - Ifges(7,6) * t146) * t676 - t20 * t131 + t21 * t132 - g(1) * t597 - g(2) * t598 - t536 * t646 + (t146 * t21 + t20 * t764) * mrSges(7,3) + t521 + (-Ifges(7,2) * t146 + t142 + t80) * t692;];
tau  = t1;
