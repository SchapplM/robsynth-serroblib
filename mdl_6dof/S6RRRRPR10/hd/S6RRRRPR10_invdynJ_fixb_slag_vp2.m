% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR10
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
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR10_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR10_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR10_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:01:00
% EndTime: 2019-03-09 23:02:33
% DurationCPUTime: 53.53s
% Computational Cost: add. (24691->1048), mult. (58667->1370), div. (0->0), fcn. (46856->14), ass. (0->507)
t413 = cos(qJ(2));
t409 = sin(qJ(2));
t404 = sin(pkin(6));
t567 = qJD(1) * t404;
t523 = t409 * t567;
t405 = cos(pkin(6));
t566 = qJD(1) * t405;
t551 = pkin(1) * t566;
t325 = -pkin(8) * t523 + t413 * t551;
t465 = (pkin(2) * t409 - pkin(9) * t413) * t404;
t326 = qJD(1) * t465;
t408 = sin(qJ(3));
t412 = cos(qJ(3));
t218 = t412 * t325 + t408 * t326;
t522 = t413 * t567;
t494 = t408 * t522;
t414 = -pkin(10) - pkin(9);
t528 = qJD(3) * t414;
t791 = -pkin(10) * t494 - t408 * t528 + t218;
t644 = cos(qJ(4));
t524 = t644 * t412;
t492 = qJD(3) * t524;
t514 = qJD(4) * t644;
t407 = sin(qJ(4));
t586 = t407 * t408;
t722 = qJD(3) + qJD(4);
t258 = -t412 * t514 + t586 * t722 - t492;
t269 = -t407 * t494 + t522 * t524;
t790 = t258 + t269;
t585 = t407 * t412;
t347 = t408 * t644 + t585;
t259 = t722 * t347;
t268 = t347 * t522;
t783 = t259 - t268;
t377 = pkin(8) * t522;
t328 = t409 * t551 + t377;
t563 = qJD(3) * t408;
t727 = -t328 + (-t494 + t563) * pkin(3);
t386 = qJD(2) + t566;
t294 = t386 * t408 + t412 * t523;
t446 = -t386 * t412 + t408 * t523;
t424 = t294 * t644 - t407 * t446;
t658 = -t424 / 0.2e1;
t789 = Ifges(6,2) * t658;
t217 = -t325 * t408 + t412 * t326;
t578 = t412 * t413;
t449 = (pkin(3) * t409 - pkin(10) * t578) * t404;
t177 = qJD(1) * t449 + t217;
t365 = t414 * t408;
t366 = t414 * t412;
t272 = t407 * t365 - t366 * t644;
t731 = -qJD(4) * t272 - t177 * t644 + t407 * t791 + t414 * t492;
t788 = qJ(5) * t790 - qJD(5) * t347 + t727;
t755 = -mrSges(6,3) + mrSges(5,2);
t406 = sin(qJ(6));
t411 = cos(qJ(6));
t772 = mrSges(7,1) * t406 + mrSges(7,2) * t411;
t777 = t755 - t772;
t558 = qJD(1) * qJD(2);
t332 = (qJDD(1) * t409 + t413 * t558) * t404;
t556 = qJDD(1) * t405;
t383 = qJDD(2) + t556;
t432 = t446 * qJD(3);
t192 = t332 * t412 + t383 * t408 - t432;
t273 = t644 * t446;
t470 = -t408 * t332 + t412 * t383;
t421 = qJD(3) * t294 - t470;
t561 = qJD(4) * t407;
t95 = qJD(4) * t273 - t644 * t192 + t294 * t561 + t407 * t421;
t681 = -t95 / 0.2e1;
t96 = qJD(4) * t424 + t407 * t192 + t421 * t644;
t679 = -t96 / 0.2e1;
t203 = t294 * t407 + t273;
t661 = -t203 / 0.2e1;
t565 = qJD(2) * t409;
t521 = t404 * t565;
t557 = qJDD(1) * t404;
t331 = -qJD(1) * t521 + t413 * t557;
t317 = qJDD(3) - t331;
t308 = qJDD(4) + t317;
t652 = t308 / 0.2e1;
t360 = qJD(3) - t522;
t353 = -qJD(4) - t360;
t649 = -t353 / 0.2e1;
t648 = t353 / 0.2e1;
t657 = t424 / 0.2e1;
t787 = -mrSges(5,1) + mrSges(6,2);
t754 = Ifges(6,4) - Ifges(5,5);
t753 = -Ifges(5,6) + Ifges(6,5);
t752 = Ifges(5,3) + Ifges(6,1);
t594 = t404 * t409;
t677 = pkin(4) + pkin(11);
t500 = t677 * t594;
t786 = -pkin(5) * t790 + qJD(1) * t500 - t731;
t314 = (-pkin(2) * t413 - pkin(9) * t409 - pkin(1)) * t404;
t278 = qJD(1) * t314;
t267 = pkin(9) * t386 + t328;
t584 = t408 * t267;
t188 = t278 * t412 - t584;
t639 = pkin(10) * t294;
t156 = t188 - t639;
t189 = t412 * t267 + t408 * t278;
t157 = -pkin(10) * t446 + t189;
t587 = t407 * t157;
t99 = t156 * t644 - t587;
t785 = pkin(3) * t514 - t99;
t784 = t677 * t783 + t788;
t733 = -t407 * t177 + t365 * t514 + t366 * t561 + t528 * t585 - t644 * t791;
t591 = t404 * t413;
t643 = pkin(1) * t405;
t343 = pkin(8) * t591 + t409 * t643;
t313 = pkin(9) * t405 + t343;
t215 = t412 * t313 + t408 * t314;
t542 = t408 * t594;
t590 = t405 * t412;
t463 = t542 - t590;
t782 = t463 * pkin(10) - t215;
t678 = t96 / 0.2e1;
t781 = t678 - t679;
t680 = t95 / 0.2e1;
t780 = t680 - t681;
t165 = t203 * t411 + t353 * t406;
t55 = qJD(6) * t165 + t308 * t411 + t406 * t96;
t166 = t203 * t406 - t353 * t411;
t56 = -qJD(6) * t166 - t308 * t406 + t411 * t96;
t94 = qJDD(6) - t95;
t10 = Ifges(7,5) * t55 + Ifges(7,6) * t56 + Ifges(7,3) * t94;
t387 = pkin(8) * t594;
t550 = qJD(2) * t643;
t497 = qJD(1) * t550;
t544 = pkin(1) * t556;
t226 = -qJD(2) * t377 - qJDD(1) * t387 - t409 * t497 + t413 * t544;
t212 = -t383 * pkin(2) - t226;
t137 = pkin(3) * t421 + t212;
t431 = pkin(3) * t360 + t156;
t428 = t407 * t431;
t225 = pkin(8) * t331 + t409 * t544 + t413 * t497;
t211 = pkin(9) * t383 + t225;
t545 = pkin(1) * t557;
t216 = -pkin(2) * t331 - pkin(9) * t332 - t545;
t103 = -qJD(3) * t189 - t211 * t408 + t412 * t216;
t61 = pkin(3) * t317 - pkin(10) * t192 + t103;
t562 = qJD(3) * t412;
t535 = t412 * t211 + t408 * t216 + t278 * t562;
t68 = t470 * pkin(10) + (-t584 - t639) * qJD(3) + t535;
t19 = -qJD(4) * t428 - t157 * t514 - t407 * t68 + t61 * t644;
t439 = qJDD(5) - t19;
t16 = -t308 * pkin(4) + t439;
t416 = t95 * qJ(5) - qJD(5) * t424 + t137;
t23 = t96 * pkin(4) + t416;
t682 = t94 / 0.2e1;
t686 = t56 / 0.2e1;
t687 = t55 / 0.2e1;
t14 = t677 * t96 + t416;
t756 = pkin(5) * t424;
t145 = t644 * t431;
t84 = -t145 + t587;
t466 = t84 + t756;
t768 = t466 + qJD(5);
t57 = t353 * t677 + t768;
t266 = -t386 * pkin(2) - t325;
t206 = pkin(3) * t446 + t266;
t419 = -qJ(5) * t424 + t206;
t69 = t203 * t677 + t419;
t25 = -t406 * t69 + t411 * t57;
t5 = -t95 * pkin(5) - t308 * t677 + t439;
t1 = qJD(6) * t25 + t14 * t411 + t406 * t5;
t26 = t406 * t57 + t411 * t69;
t608 = qJD(6) * t26;
t2 = -t14 * t406 + t411 * t5 - t608;
t713 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t762 = -t308 / 0.2e1;
t779 = mrSges(6,1) * t16 + mrSges(5,2) * t137 - mrSges(5,3) * t19 - mrSges(6,3) * t23 + 0.2e1 * Ifges(5,1) * t681 + 0.2e1 * Ifges(5,4) * t679 + Ifges(6,4) * t762 + Ifges(7,5) * t687 - Ifges(6,2) * t780 - Ifges(6,6) * t781 + Ifges(7,6) * t686 + Ifges(7,3) * t682 + t713 + t10 / 0.2e1 + (-t754 + Ifges(5,5)) * t652;
t778 = mrSges(7,3) - t787;
t199 = Ifges(5,4) * t203;
t131 = Ifges(5,1) * t424 - t353 * Ifges(5,5) - t199;
t200 = qJD(6) + t424;
t613 = t200 * Ifges(7,3);
t614 = t166 * Ifges(7,5);
t616 = t165 * Ifges(7,6);
t77 = t613 + t614 + t616;
t737 = -qJD(5) - t84;
t80 = pkin(4) * t353 - t737;
t771 = t80 * mrSges(6,1) + t84 * mrSges(5,3);
t766 = t131 / 0.2e1 + t77 / 0.2e1 + t771;
t651 = t317 / 0.2e1;
t664 = t192 / 0.2e1;
t758 = -t421 / 0.2e1;
t676 = Ifges(4,1) * t664 + Ifges(4,4) * t758 + Ifges(4,5) * t651;
t635 = t203 * pkin(5);
t742 = -qJD(5) - t785;
t610 = qJ(5) * t203;
t732 = qJ(5) * t523 - t733;
t27 = mrSges(7,1) * t94 - mrSges(7,3) * t55;
t28 = -mrSges(7,2) * t94 + mrSges(7,3) * t56;
t474 = t25 * t406 - t26 * t411;
t425 = m(7) * (-qJD(6) * t474 + t1 * t406 + t2 * t411);
t775 = t411 * t27 + t406 * t28 + t425;
t483 = mrSges(7,1) * t411 - mrSges(7,2) * t406;
t153 = t644 * t157;
t85 = t153 + t428;
t81 = t353 * qJ(5) - t85;
t59 = -t81 - t635;
t615 = t166 * Ifges(7,4);
t78 = t165 * Ifges(7,2) + t200 * Ifges(7,6) + t615;
t685 = -t78 / 0.2e1;
t774 = t411 * t685 + t59 * t483;
t645 = cos(qJ(1));
t525 = t645 * t413;
t410 = sin(qJ(1));
t582 = t409 * t410;
t337 = -t405 * t582 + t525;
t592 = t404 * t412;
t261 = -t337 * t408 + t410 * t592;
t773 = m(4) * pkin(9) + m(7) * pkin(5) + mrSges(6,1) + mrSges(4,3) + mrSges(5,3);
t770 = t81 * mrSges(6,1) - t85 * mrSges(5,3);
t769 = qJD(3) * t266;
t105 = t203 * pkin(4) + t419;
t198 = Ifges(6,6) * t203;
t701 = Ifges(6,4) * t649 + t789 + t198 / 0.2e1 - t25 * mrSges(7,1) - t206 * mrSges(5,2) + t26 * mrSges(7,2) + t105 * mrSges(6,3);
t660 = t203 / 0.2e1;
t708 = -t206 * mrSges(5,1) + t105 * mrSges(6,2) + Ifges(6,5) * t648 + Ifges(5,6) * t649 + (Ifges(5,2) + Ifges(6,3)) * t661 + (Ifges(5,4) + Ifges(6,6)) * t657;
t767 = -Ifges(5,4) * t657 - Ifges(5,2) * t661 + Ifges(6,6) * t658 + Ifges(6,3) * t660 + t649 * t753 - t708 + t770;
t485 = mrSges(4,1) * t408 + mrSges(4,2) * t412;
t765 = -t266 * t485 * t591 - t386 * t404 * (Ifges(3,5) * t413 - Ifges(3,6) * t409) / 0.2e1;
t476 = Ifges(7,5) * t406 + Ifges(7,6) * t411;
t618 = Ifges(7,4) * t406;
t478 = Ifges(7,2) * t411 + t618;
t617 = Ifges(7,4) * t411;
t480 = Ifges(7,1) * t406 + t617;
t646 = -t406 / 0.2e1;
t663 = -t200 / 0.2e1;
t667 = -t166 / 0.2e1;
t669 = -t165 / 0.2e1;
t164 = Ifges(7,4) * t165;
t79 = Ifges(7,1) * t166 + Ifges(7,5) * t200 + t164;
t764 = t476 * t663 + t478 * t669 + t480 * t667 + t79 * t646 - t770 + t774;
t763 = -m(7) - m(6);
t761 = t331 / 0.2e1;
t760 = t332 / 0.2e1;
t759 = t405 / 0.2e1;
t757 = pkin(4) * t424;
t346 = -t524 + t586;
t398 = pkin(3) * t412 + pkin(2);
t468 = -qJ(5) * t347 - t398;
t207 = t346 * t677 + t468;
t271 = -t644 * t365 - t366 * t407;
t221 = pkin(5) * t347 + t271;
t144 = t207 * t411 + t221 * t406;
t750 = -qJD(6) * t144 - t406 * t784 + t411 * t786;
t143 = -t207 * t406 + t221 * t411;
t749 = qJD(6) * t143 + t406 * t786 + t411 * t784;
t20 = -mrSges(7,1) * t56 + mrSges(7,2) * t55;
t73 = mrSges(6,1) * t96 - mrSges(6,3) * t308;
t748 = -t73 + t20;
t747 = Ifges(4,6) * t421;
t746 = t225 * mrSges(3,2);
t745 = t226 * mrSges(3,1);
t743 = t756 - t742;
t741 = t424 * t677;
t736 = -pkin(5) * t783 - t732;
t735 = pkin(4) * t783 + t788;
t173 = mrSges(6,1) * t424 - mrSges(6,2) * t353;
t175 = -mrSges(5,1) * t353 - mrSges(5,3) * t424;
t734 = t175 - t173;
t730 = pkin(4) * t523 - t731;
t219 = t268 * t411 - t406 * t523;
t560 = qJD(6) * t406;
t580 = t411 * t259;
t460 = t346 * t560 - t580;
t729 = t219 + t460;
t220 = t268 * t406 + t411 * t523;
t559 = qJD(6) * t411;
t589 = t406 * t259;
t461 = t346 * t559 + t589;
t728 = t220 - t461;
t333 = t405 * t408 + t409 * t592;
t434 = t644 * t463;
t223 = t333 * t407 + t434;
t588 = t406 * t413;
t369 = t404 * t588;
t196 = t223 * t411 + t369;
t642 = pkin(1) * t413;
t342 = t405 * t642 - t387;
t726 = t308 * t752 + t753 * t96 + t754 * t95;
t102 = -t267 * t563 + t535;
t725 = t102 * t412 - t103 * t408;
t623 = Ifges(3,4) * t409;
t693 = t404 ^ 2;
t721 = (t409 * (Ifges(3,1) * t413 - t623) / 0.2e1 - pkin(1) * (mrSges(3,1) * t409 + mrSges(3,2) * t413)) * t693;
t720 = -m(7) * pkin(11) - mrSges(7,3);
t719 = Ifges(4,5) * t294 - Ifges(4,6) * t446 + Ifges(4,3) * t360 + t203 * t753 - t353 * t752 - t424 * t754;
t403 = qJ(3) + qJ(4);
t400 = sin(t403);
t401 = cos(t403);
t486 = -mrSges(4,1) * t412 + mrSges(4,2) * t408;
t718 = -m(4) * pkin(2) + t400 * t755 + t401 * t787 + t486;
t310 = t400 * t594 - t405 * t401;
t311 = t400 * t405 + t401 * t594;
t717 = t310 * t778 + t311 * t777;
t593 = t404 * t410;
t253 = t337 * t400 - t401 * t593;
t254 = t337 * t401 + t400 * t593;
t716 = t253 * t778 + t254 * t777;
t526 = t645 * t409;
t581 = t410 * t413;
t335 = t405 * t526 + t581;
t527 = t404 * t645;
t249 = t335 * t400 + t401 * t527;
t250 = t335 * t401 - t400 * t527;
t715 = t249 * t778 + t250 * t777;
t327 = qJD(2) * t465;
t329 = t342 * qJD(2);
t504 = t412 * t327 - t408 * t329;
t116 = qJD(2) * t449 + qJD(3) * t782 + t504;
t564 = qJD(2) * t413;
t519 = t408 * t564;
t493 = t404 * t519;
t532 = t314 * t562 + t408 * t327 + t412 * t329;
t583 = t408 * t313;
t638 = pkin(10) * t333;
t124 = -pkin(10) * t493 + (-t583 - t638) * qJD(3) + t532;
t214 = t412 * t314 - t583;
t163 = -pkin(3) * t591 + t214 - t638;
t452 = -t407 * t116 - t644 * t124 - t163 * t514 - t561 * t782;
t29 = -t404 * (qJ(5) * t565 - qJD(5) * t413) + t452;
t499 = mrSges(3,3) * t523;
t712 = -m(4) * t266 + mrSges(3,1) * t386 - mrSges(4,1) * t446 - t294 * mrSges(4,2) - t499;
t711 = t720 + t787;
t710 = -t103 * mrSges(4,1) + t102 * mrSges(4,2);
t707 = -mrSges(3,2) + t773;
t706 = t400 * t772 + mrSges(3,1) - t718;
t18 = qJD(4) * t145 - t157 * t561 + t407 * t61 + t644 * t68;
t13 = -qJ(5) * t308 + qJD(5) * t353 - t18;
t705 = -t19 * mrSges(5,1) + t18 * mrSges(5,2) - t16 * mrSges(6,2) + t13 * mrSges(6,3);
t704 = -t84 * mrSges(5,1) - t85 * mrSges(5,2) + t80 * mrSges(6,2) - t81 * mrSges(6,3);
t702 = -t483 - t773;
t700 = mrSges(3,2) + t702;
t699 = -Ifges(5,4) * t658 - Ifges(5,2) * t660 + Ifges(6,6) * t657 + Ifges(6,3) * t661 + t648 * t753 + t708;
t697 = mrSges(5,1) * t137 + mrSges(6,1) * t13 - mrSges(6,2) * t23 - mrSges(5,3) * t18 + Ifges(5,6) * t762 + 0.2e1 * Ifges(6,6) * t680 + 0.2e1 * Ifges(6,3) * t678 + (Ifges(6,5) + t753) * t652 + t781 * Ifges(5,2) + t780 * Ifges(5,4);
t696 = Ifges(5,1) * t658 + Ifges(5,4) * t660 + Ifges(7,5) * t667 - Ifges(6,2) * t657 - Ifges(6,6) * t661 + Ifges(7,6) * t669 + Ifges(7,3) * t663 - t648 * t754 + t701;
t662 = t200 / 0.2e1;
t666 = t166 / 0.2e1;
t668 = t165 / 0.2e1;
t695 = -Ifges(5,1) * t657 - Ifges(5,4) * t661 - Ifges(7,5) * t666 + Ifges(6,6) * t660 - Ifges(7,6) * t668 - Ifges(7,3) * t662 + t649 * t754 + t701 + t789;
t691 = Ifges(7,1) * t687 + Ifges(7,4) * t686 + Ifges(7,5) * t682;
t684 = t78 / 0.2e1;
t683 = t79 / 0.2e1;
t279 = Ifges(4,4) * t446;
t185 = t294 * Ifges(4,1) + t360 * Ifges(4,5) - t279;
t665 = t185 / 0.2e1;
t653 = t294 / 0.2e1;
t647 = t360 / 0.2e1;
t641 = pkin(3) * t294;
t640 = pkin(3) * t407;
t637 = pkin(11) * t253;
t636 = pkin(11) * t310;
t634 = t249 * pkin(11);
t626 = mrSges(4,3) * t294;
t625 = mrSges(7,3) * t406;
t624 = mrSges(7,3) * t411;
t622 = Ifges(3,4) * t413;
t621 = Ifges(4,4) * t294;
t620 = Ifges(4,4) * t408;
t619 = Ifges(4,4) * t412;
t609 = qJ(5) * t400;
t334 = -t405 * t525 + t582;
t604 = t334 * t401;
t603 = t334 * t406;
t602 = t334 * t411;
t336 = t405 * t581 + t526;
t601 = t336 * t401;
t599 = t346 * t406;
t598 = t346 * t411;
t579 = t411 * t413;
t110 = -mrSges(7,1) * t165 + mrSges(7,2) * t166;
t172 = mrSges(6,1) * t203 + mrSges(6,3) * t353;
t577 = t110 - t172;
t109 = t407 * t163 - t644 * t782;
t572 = -t334 * t398 - t335 * t414;
t571 = -t336 * t398 - t337 * t414;
t520 = t404 * t564;
t330 = pkin(8) * t520 + t409 * t550;
t568 = t645 * pkin(1) + pkin(8) * t593;
t389 = pkin(4) * t591;
t552 = t644 * pkin(3);
t547 = mrSges(4,3) * t563;
t546 = mrSges(4,3) * t562;
t543 = qJ(5) * t591;
t541 = t408 * t593;
t539 = t404 * t579;
t536 = Ifges(4,5) * t192 + Ifges(4,3) * t317 - t747;
t530 = Ifges(3,5) * t332 + Ifges(3,6) * t331 + Ifges(3,3) * t383;
t518 = t412 * t564;
t516 = t594 / 0.2e1;
t512 = -t567 / 0.2e1;
t509 = -t560 / 0.2e1;
t74 = -t95 * mrSges(6,1) + t308 * mrSges(6,2);
t508 = -pkin(1) * t410 + pkin(8) * t527;
t507 = -t249 * pkin(4) + t250 * qJ(5);
t506 = -t253 * pkin(4) + qJ(5) * t254;
t505 = -t310 * pkin(4) + qJ(5) * t311;
t98 = t156 * t407 + t153;
t370 = t408 * t527;
t503 = -t335 * t412 + t370;
t502 = pkin(3) * t542;
t498 = mrSges(3,3) * t522;
t496 = -pkin(4) * t604 - t334 * t609 + t572;
t495 = -pkin(4) * t601 - t336 * t609 + t571;
t397 = -t552 - pkin(4);
t490 = t413 * t512;
t488 = t261 * pkin(3);
t108 = t163 * t644 + t407 * t782;
t481 = Ifges(4,1) * t412 - t620;
t479 = -Ifges(4,2) * t408 + t619;
t477 = Ifges(4,5) * t412 - Ifges(4,6) * t408;
t475 = t641 + t610;
t473 = pkin(3) * t541 - t336 * t414 + t337 * t398 + t568;
t224 = t333 * t644 - t407 * t463;
t227 = t502 + t387 + (-t398 - t642) * t405;
t420 = -t224 * qJ(5) + t227;
t101 = t223 * t677 + t420;
t107 = t389 - t108;
t72 = t224 * pkin(5) + pkin(11) * t591 + t107;
t43 = t101 * t411 + t406 * t72;
t42 = -t101 * t406 + t411 * t72;
t120 = -mrSges(7,2) * t200 + mrSges(7,3) * t165;
t121 = mrSges(7,1) * t200 - mrSges(7,3) * t166;
t472 = t120 * t411 - t121 * t406;
t469 = pkin(3) * t590 - t502;
t106 = t543 - t109;
t464 = -t223 * t406 + t539;
t458 = pkin(3) * t370 + t334 * t414 - t335 * t398 + t508;
t457 = (t413 * Ifges(3,2) + t623) * t404;
t455 = t481 * qJD(3);
t454 = t477 * qJD(3);
t453 = -t116 * t644 + t407 * t124 + t163 * t561 - t514 * t782;
t450 = -g(1) * t253 - g(2) * t249 - g(3) * t310;
t447 = t335 * t408 + t412 * t527;
t441 = -t409 * t563 + t518;
t440 = -t409 * t562 - t519;
t438 = t488 + t506;
t437 = t447 * pkin(3);
t435 = t446 * mrSges(4,3);
t433 = t469 + t505;
t429 = -mrSges(4,1) * t463 - t333 * mrSges(4,2);
t426 = -t437 + t507;
t423 = qJD(3) * t333 + t493;
t422 = -qJD(3) * t463 + t404 * t518;
t210 = pkin(3) * t423 + t330;
t135 = qJD(4) * t434 + t333 * t561 + t407 * t423 - t422 * t644;
t418 = t135 * qJ(5) - t224 * qJD(5) + t210;
t11 = Ifges(7,4) * t55 + Ifges(7,2) * t56 + Ifges(7,6) * t94;
t7 = -pkin(5) * t96 - t13;
t417 = t25 * mrSges(7,3) * t560 + (Ifges(7,1) * t411 - t618) * t687 + (-Ifges(7,2) * t406 + t617) * t686 + t7 * t772 + t79 * t509 + t11 * t646 + (Ifges(7,5) * t411 - Ifges(7,6) * t406) * t682 + t411 * t691 - t705 + t726 + t774 * qJD(6) - (t165 * t478 + t166 * t480 + t200 * t476) * qJD(6) / 0.2e1;
t395 = qJ(5) + t640;
t375 = Ifges(3,4) * t522;
t350 = t398 * t591;
t338 = (-mrSges(3,1) * t413 + mrSges(3,2) * t409) * t404;
t324 = -mrSges(3,2) * t386 + t498;
t312 = t387 + (-pkin(2) - t642) * t405;
t265 = Ifges(3,1) * t523 + t386 * Ifges(3,5) + t375;
t264 = t386 * Ifges(3,6) + qJD(1) * t457;
t262 = t337 * t412 + t541;
t234 = pkin(4) * t346 + t468;
t233 = mrSges(4,1) * t360 - t626;
t232 = -t360 * mrSges(4,2) - t435;
t222 = -t346 * pkin(5) + t272;
t194 = t253 * t406 + t336 * t411;
t193 = t253 * t411 - t336 * t406;
t184 = -Ifges(4,2) * t446 + Ifges(4,6) * t360 + t621;
t174 = mrSges(5,2) * t353 - mrSges(5,3) * t203;
t159 = -t317 * mrSges(4,2) - mrSges(4,3) * t421;
t158 = mrSges(4,1) * t317 - mrSges(4,3) * t192;
t147 = -qJD(3) * t215 + t504;
t146 = -t313 * t563 + t532;
t140 = -mrSges(6,2) * t203 - mrSges(6,3) * t424;
t139 = mrSges(5,1) * t203 + mrSges(5,2) * t424;
t138 = t610 + t757;
t136 = qJD(4) * t224 + t407 * t422 + t423 * t644;
t132 = mrSges(4,1) * t421 + t192 * mrSges(4,2);
t125 = t223 * pkin(4) + t420;
t117 = t475 + t757;
t111 = Ifges(4,4) * t192 - t421 * Ifges(4,2) + Ifges(4,6) * t317;
t104 = t610 + t741;
t89 = qJD(6) * t196 + t136 * t406 + t411 * t521;
t88 = qJD(6) * t464 + t136 * t411 - t406 * t521;
t86 = t475 + t741;
t82 = -pkin(5) * t223 - t106;
t76 = -mrSges(5,2) * t308 - mrSges(5,3) * t96;
t75 = mrSges(5,1) * t308 + mrSges(5,3) * t95;
t64 = t98 - t635;
t63 = t85 - t635;
t48 = t136 * pkin(4) + t418;
t45 = mrSges(5,1) * t96 - mrSges(5,2) * t95;
t44 = -mrSges(6,2) * t96 + mrSges(6,3) * t95;
t37 = t104 * t411 + t406 * t63;
t36 = -t104 * t406 + t411 * t63;
t33 = t406 * t64 + t411 * t86;
t32 = -t406 * t86 + t411 * t64;
t31 = -pkin(4) * t521 + t453;
t30 = t136 * t677 + t418;
t22 = -pkin(5) * t136 - t29;
t21 = -t135 * pkin(5) - qJD(2) * t500 + t453;
t4 = -qJD(6) * t43 + t21 * t411 - t30 * t406;
t3 = qJD(6) * t42 + t21 * t406 + t30 * t411;
t6 = [(t719 * t516 - t765) * qJD(2) + (t695 - t766) * t135 + t767 * t136 + (Ifges(7,5) * t89 + Ifges(7,6) * t88) * t662 + t721 * t558 + (Ifges(7,4) * t89 + Ifges(7,2) * t88) * t668 - t452 * t174 + m(5) * (t108 * t19 + t109 * t18 + t137 * t227 + t206 * t210 - t452 * t85 + t453 * t84) - t453 * t175 + t7 * (-mrSges(7,1) * t196 - mrSges(7,2) * t464) - t464 * t691 + (-Ifges(7,5) * t464 + Ifges(7,6) * t196) * t682 + (-Ifges(7,4) * t464 + Ifges(7,2) * t196) * t686 + (t1 * t196 + t2 * t464 - t25 * t89 + t26 * t88) * mrSges(7,3) + (-mrSges(4,2) * t769 - t102 * mrSges(4,3) - Ifges(4,4) * t664 - Ifges(4,6) * t651 - t111 / 0.2e1 - Ifges(4,2) * t758) * t463 + m(4) * (t102 * t215 + t103 * t214 + t146 * t189 + t147 * t188 + t212 * t312) + t405 * t745 + (mrSges(4,1) * t769 - mrSges(4,3) * t103 + 0.2e1 * t676) * t333 + (t225 * t591 - t226 * t594 - t325 * t520 - t328 * t521 + t331 * t343 - t332 * t342) * mrSges(3,3) - t423 * t184 / 0.2e1 + t697 * t223 + t189 * (-t405 * t547 + (-mrSges(4,2) * t565 + mrSges(4,3) * t440) * t404) - t212 * t429 - t338 * t545 + (-m(3) * t325 - t712) * t330 + (-Ifges(7,1) * t464 + Ifges(7,4) * t196) * t687 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t693 + t225 * t343 + t226 * t342 + t328 * t329) + (Ifges(3,1) * t332 + Ifges(3,4) * t331) * t516 + (t602 * mrSges(7,1) - t603 * mrSges(7,2) + t410 * mrSges(2,1) + mrSges(2,2) * t645 - m(5) * t458 - m(4) * (-pkin(2) * t335 + t508) - t503 * mrSges(4,1) - t447 * mrSges(4,2) - m(3) * t508 + t335 * mrSges(3,1) - mrSges(3,3) * t527 + t763 * (-pkin(4) * t250 - qJ(5) * t249 + t458) - t777 * t249 - t711 * t250 + t707 * t334) * g(1) + (Ifges(7,1) * t89 + Ifges(7,4) * t88) * t666 + t89 * t683 + t88 * t684 + t188 * (-t405 * t546 + (mrSges(4,1) * t565 - mrSges(4,3) * t441) * t404) - t446 * (t479 * t405 * qJD(3) + (Ifges(4,4) * t441 + Ifges(4,2) * t440 + Ifges(4,6) * t565) * t404) / 0.2e1 + Ifges(2,3) * qJDD(1) + t146 * t232 + t147 * t233 + t227 * t45 + t214 * t158 + t215 * t159 + t210 * t139 + t196 * t11 / 0.2e1 + t29 * t172 + t31 * t173 + t48 * t140 + t3 * t120 + t4 * t121 + t125 * t44 + t106 * t73 + t107 * t74 + t108 * t75 + t109 * t76 + t22 * t110 + t59 * (-mrSges(7,1) * t88 + mrSges(7,2) * t89) + t82 * t20 + (t404 * t265 + t693 * qJD(1) * (-Ifges(3,2) * t409 + t622)) * t564 / 0.2e1 + t42 * t27 + t43 * t28 + (-t194 * mrSges(7,1) - t193 * mrSges(7,2) - mrSges(2,1) * t645 + t410 * mrSges(2,2) - m(3) * t568 - t337 * mrSges(3,1) - mrSges(3,3) * t593 - m(5) * t473 - m(4) * (pkin(2) * t337 + t568) - t262 * mrSges(4,1) - t261 * mrSges(4,2) + t763 * (t254 * pkin(4) + qJ(5) * t253 + t473) + t755 * t253 + t711 * t254 - t707 * t336) * g(2) + (Ifges(3,4) * t760 + Ifges(3,2) * t761 + Ifges(3,6) * t383 / 0.2e1 + t747 / 0.2e1 - Ifges(6,5) * t678 - Ifges(5,6) * t679 - Ifges(6,4) * t680 - Ifges(5,5) * t681 - Ifges(4,3) * t651 - Ifges(4,5) * t664 - t752 * t652 + t705 + t710) * t591 + m(6) * (t105 * t48 + t106 * t13 + t107 * t16 + t125 * t23 + t29 * t81 + t31 * t80) + m(7) * (t1 * t43 + t2 * t42 + t22 * t59 + t25 * t4 + t26 * t3 + t7 * t82) + t779 * t224 - (t536 + t726) * t591 / 0.2e1 + t312 * t132 + t329 * t324 - t404 * pkin(1) * (-mrSges(3,1) * t331 + mrSges(3,2) * t332) + t530 * t759 + (Ifges(3,5) * t405 + (t409 * Ifges(3,1) + t622) * t404) * t760 + (Ifges(3,6) * t405 + t457) * t761 - t405 * t746 + (-t264 / 0.2e1 + Ifges(5,5) * t657 + Ifges(6,4) * t658 + Ifges(6,5) * t660 + Ifges(5,6) * t661 + t752 * t649 + t704) * t521 + (Ifges(3,3) * t759 + (Ifges(3,5) * t409 + Ifges(3,6) * t413) * t404 / 0.2e1 + t342 * mrSges(3,1) - t343 * mrSges(3,2) + Ifges(3,5) * t516) * t383 + (t405 * t454 + (Ifges(4,5) * t441 + Ifges(4,6) * t440 + Ifges(4,3) * t565) * t404) * t647 + (t405 * t455 + (Ifges(4,1) * t441 + Ifges(4,4) * t440 + Ifges(4,5) * t565) * t404) * t653 + t422 * t665; (-t563 / 0.2e1 + t494 / 0.2e1) * t184 + (t74 - t75) * t271 + (t76 - t73) * t272 + (-qJD(1) * t721 + t765) * qJD(1) + (Ifges(7,1) * t461 - Ifges(7,4) * t460) * t666 + (Ifges(7,1) * t220 + Ifges(7,4) * t219) * t667 + t767 * t259 + (t498 - t324) * t325 + (-pkin(2) * t212 - t188 * t217 - t189 * t218) * m(4) + (t699 - t770) * t268 + (t412 * t185 + t265 + t375) * t490 + (t696 - t771) * t269 + (t695 - t771) * t258 + (t360 * (Ifges(4,3) * t409 + t413 * t477) + t294 * (Ifges(4,5) * t409 + t413 * t481) + t719 * t409) * t512 + t485 * t769 + (t131 + t77) * (-t269 / 0.2e1 - t258 / 0.2e1) + (-t546 - (mrSges(4,1) * t409 - mrSges(4,3) * t578) * t567) * t188 + (t476 * t682 + t478 * t686 + t480 * t687 - t483 * t7 + t509 * t78 + t697) * t346 + t745 - t746 + (-t547 - (-mrSges(4,3) * t408 * t413 - mrSges(4,2) * t409) * t567) * t189 + (t499 + t712) * t328 + (-m(6) * t495 - m(7) * (-pkin(11) * t601 + t495) + mrSges(7,3) * t601 - m(5) * t571 + t700 * t337 + t706 * t336) * g(1) + (-m(7) * (-pkin(11) * t604 + t496) + mrSges(7,3) * t604 - m(6) * t496 - m(5) * t572 + t700 * t335 + t706 * t334) * g(2) - t479 * t432 / 0.2e1 + (Ifges(7,4) * t461 - Ifges(7,2) * t460) * t668 + (Ifges(7,4) * t220 + Ifges(7,2) * t219) * t669 + t580 * t684 + t219 * t685 + t599 * t691 + t408 * t676 + t589 * t683 + (Ifges(7,5) * t461 - Ifges(7,6) * t460) * t662 + (Ifges(7,5) * t220 + Ifges(7,6) * t219) * t663 - t218 * t232 - t217 * t233 + t234 * t44 + t222 * t20 - t220 * t79 / 0.2e1 + t212 * t486 + t143 * t27 + t144 * t28 - pkin(2) * t132 + (-m(5) * t350 + t338 + t763 * (t401 * t389 + t400 * t543 + t350) + ((-mrSges(7,1) * t588 - mrSges(7,2) * t579) * t400 + ((m(5) - t763) * t414 + t702) * t409 + (t401 * t720 + t718) * t413) * t404) * g(3) + (t446 * (Ifges(4,6) * t409 + t413 * t479) + t409 * t264) * t567 / 0.2e1 + (qJD(6) * t79 + t11) * t598 / 0.2e1 + t412 * t111 / 0.2e1 + t530 + t779 * t347 + (m(4) * ((-t188 * t412 - t189 * t408) * qJD(3) + t725) - t233 * t562 - t232 * t563 - t408 * t158 + t412 * t159) * pkin(9) + t725 * mrSges(4,3) + t727 * t139 + (mrSges(7,1) * t729 - mrSges(7,2) * t728) * t59 + (t1 * t598 - t2 * t599 + t25 * t728 - t26 * t729) * mrSges(7,3) + t730 * t173 + t731 * t175 + t732 * t172 + t733 * t174 + (-t137 * t398 + t18 * t272 - t19 * t271 + t727 * t206 - t731 * t84 + t733 * t85) * m(5) + t735 * t140 + (t105 * t735 - t13 * t272 + t16 * t271 + t23 * t234 + t730 * t80 + t732 * t81) * m(6) + t736 * t110 + (Ifges(4,2) * t412 + t620) * t758 - t398 * t45 + t749 * t120 + t750 * t121 + (t1 * t144 + t143 * t2 + t222 * t7 + t750 * t25 + t749 * t26 + t736 * t59) * m(7) + (Ifges(6,4) * t657 + Ifges(5,5) * t658 + Ifges(6,5) * t661 - Ifges(3,2) * t490 + Ifges(5,6) * t660 + t648 * t752 - t704) * t523 + t454 * t647 + (Ifges(4,5) * t408 + Ifges(4,6) * t412) * t651 + t455 * t653 + (Ifges(4,1) * t408 + t619) * t664 + t562 * t665; (-t696 + t766) * t203 + (t233 + t626) * t189 + (-m(7) * (t426 - t634) + mrSges(4,1) * t447 - mrSges(4,2) * t503 + m(5) * t437 - m(6) * t426 + t715) * g(2) + (-m(7) * (t438 - t637) - mrSges(4,1) * t261 + mrSges(4,2) * t262 - m(5) * t488 - m(6) * t438 + t716) * g(1) + (-t429 - m(7) * (t433 - t636) - m(5) * t469 - m(6) * t433 + t717) * g(3) + (t120 * t559 - t121 * t560 + t775) * (-pkin(11) + t397) + (t25 * t625 - t26 * t624 + t699 + t764) * t424 + t446 * (-Ifges(4,2) * t294 - t279) / 0.2e1 + (-t435 - t232) * t188 + (-t206 * t641 - t84 * t98 - t85 * t99 + (t644 * t19 + t18 * t407 + (t407 * t84 + t644 * t85) * qJD(4)) * pkin(3)) * m(5) - t1 * t625 - t139 * t641 - t294 * (-Ifges(4,1) * t446 - t621) / 0.2e1 - t2 * t624 + t536 - t117 * t140 - t33 * t120 - t32 * t121 - t26 * mrSges(7,3) * t559 - t266 * (t294 * mrSges(4,1) - mrSges(4,2) * t446) - t360 * (-Ifges(4,5) * t446 - Ifges(4,6) * t294) / 0.2e1 + (t411 * t121 + t406 * t120 + m(7) * (t25 * t411 + t26 * t406) + m(6) * t80 - t734) * pkin(3) * t561 + t734 * t98 + t417 - t710 + t785 * t174 + t742 * t172 + (-t105 * t117 - t13 * t395 + t16 * t397 + t742 * t81 - t80 * t98) * m(6) + (-t25 * t32 - t26 * t33 + t395 * t7 + t59 * t743) * m(7) + t743 * t110 + t75 * t552 + t397 * t74 + t748 * t395 + t76 * t640 + t184 * t653 + t446 * t665; t577 * qJD(5) + (-t1 * mrSges(7,3) - (-qJD(6) * t121 + t28) * t677) * t406 + t734 * t85 + (-(qJD(6) * t120 + t27) * t677 + (-t2 - t608) * mrSges(7,3)) * t411 + t716 * g(1) + t717 * g(3) - t677 * t425 + t748 * qJ(5) + (t613 / 0.2e1 + t614 / 0.2e1 + t616 / 0.2e1 - t199 / 0.2e1 - t198 / 0.2e1 + (Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1) * t353 - t701 + t766) * t203 + t715 * g(2) + ((Ifges(5,4) / 0.2e1 + Ifges(6,6) / 0.2e1) * t424 + (-Ifges(5,6) / 0.2e1 + Ifges(6,5) / 0.2e1) * t353 + t474 * mrSges(7,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t203 + t708 + t764) * t424 - t138 * t140 - t37 * t120 - t36 * t121 + t466 * t110 - pkin(4) * t74 + t417 + (-t172 + t174) * t84 + (-t25 * t36 - t26 * t37 + t7 * qJ(5) + (-t506 + t637) * g(1) + (-t505 + t636) * g(3) + (-t507 + t634) * g(2) + t768 * t59) * m(7) + (-pkin(4) * t16 - t506 * g(1) - t507 * g(2) - t505 * g(3) - qJ(5) * t13 - t105 * t138 + t737 * t81 - t80 * t85) * m(6); t577 * t353 + t472 * qJD(6) + (t140 + t472) * t424 + t74 + (t353 * t59 - t424 * t474 + t450) * m(7) + (t105 * t424 - t353 * t81 + t16 + t450) * m(6) + t775; -t59 * (mrSges(7,1) * t166 + mrSges(7,2) * t165) + (Ifges(7,1) * t165 - t615) * t667 + t78 * t666 + (Ifges(7,5) * t165 - Ifges(7,6) * t166) * t663 - t25 * t120 + t26 * t121 - g(1) * (mrSges(7,1) * t193 - mrSges(7,2) * t194) - g(2) * ((t249 * t411 - t603) * mrSges(7,1) + (-t249 * t406 - t602) * mrSges(7,2)) - g(3) * ((t310 * t411 + t369) * mrSges(7,1) + (-t310 * t406 + t539) * mrSges(7,2)) + (t165 * t25 + t166 * t26) * mrSges(7,3) + t10 + (-Ifges(7,2) * t166 + t164 + t79) * t669 + t713;];
tau  = t6;
