% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR9
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR9_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR9_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR9_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:47:50
% EndTime: 2019-03-09 22:49:44
% DurationCPUTime: 67.85s
% Computational Cost: add. (38187->1123), mult. (90567->1536), div. (0->0), fcn. (73872->18), ass. (0->494)
t497 = cos(qJ(2));
t492 = sin(qJ(2));
t485 = sin(pkin(6));
t599 = qJD(1) * t485;
t570 = t492 * t599;
t487 = cos(pkin(6));
t598 = qJD(1) * t487;
t587 = pkin(1) * t598;
t401 = -pkin(8) * t570 + t497 * t587;
t522 = (pkin(2) * t492 - pkin(9) * t497) * t485;
t402 = qJD(1) * t522;
t491 = sin(qJ(3));
t496 = cos(qJ(3));
t291 = -t401 * t491 + t496 * t402;
t498 = -pkin(10) - pkin(9);
t574 = qJD(3) * t498;
t615 = t496 * t497;
t830 = -(pkin(3) * t492 - pkin(10) * t615) * t599 - t291 + t496 * t574;
t292 = t496 * t401 + t491 * t402;
t569 = t497 * t599;
t547 = t491 * t569;
t829 = -pkin(10) * t547 - t491 * t574 + t292;
t595 = qJD(3) * t491;
t828 = -t595 + t547;
t463 = qJD(2) + t598;
t372 = t463 * t491 + t496 * t570;
t490 = sin(qJ(4));
t495 = cos(qJ(4));
t510 = -t463 * t496 + t491 * t570;
t276 = t490 * t372 + t495 * t510;
t691 = -t276 / 0.2e1;
t442 = qJD(3) - t569;
t432 = qJD(4) + t442;
t678 = t432 / 0.2e1;
t502 = t495 * t372 - t490 * t510;
t687 = t502 / 0.2e1;
t821 = Ifges(5,4) * t687 + Ifges(5,6) * t678;
t827 = Ifges(5,2) * t691 + t821;
t427 = t490 * t496 + t491 * t495;
t746 = qJD(3) + qJD(4);
t336 = t746 * t427;
t351 = t427 * t569;
t826 = t336 - t351;
t800 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t447 = t498 * t491;
t448 = t498 * t496;
t754 = t495 * t447 + t448 * t490;
t758 = qJD(4) * t754 + t490 * t830 - t829 * t495;
t458 = pkin(8) * t569;
t404 = t492 * t587 + t458;
t751 = -pkin(3) * t828 - t404;
t488 = -pkin(11) - qJ(5);
t512 = -m(6) * qJ(5) + m(7) * t488 + t800;
t482 = pkin(12) + qJ(6);
t476 = sin(t482);
t477 = cos(t482);
t749 = -mrSges(7,1) * t477 + mrSges(7,2) * t476;
t824 = t749 - mrSges(5,1);
t823 = qJ(5) * t570 - t758;
t426 = t490 * t491 - t495 * t496;
t335 = t746 * t426;
t352 = t426 * t569;
t822 = -qJD(5) * t427 + t751 + (t335 - t352) * qJ(5) + t826 * pkin(4);
t484 = sin(pkin(12));
t486 = cos(pkin(12));
t294 = t352 * t484 + t486 * t570;
t804 = -t335 * t484 + t294;
t295 = -t352 * t486 + t484 * t570;
t820 = t335 * t486 + t295;
t538 = -mrSges(6,1) * t486 + mrSges(6,2) * t484;
t489 = sin(qJ(6));
t494 = cos(qJ(6));
t524 = t484 * t489 - t486 * t494;
t409 = t524 * qJD(6);
t790 = t524 * t276;
t819 = t790 + t409;
t425 = t484 * t494 + t486 * t489;
t410 = t425 * qJD(6);
t791 = t425 * t276;
t818 = -t791 - t410;
t350 = pkin(9) * t463 + t404;
t390 = (-pkin(2) * t497 - pkin(9) * t492 - pkin(1)) * t485;
t361 = qJD(1) * t390;
t263 = -t350 * t491 + t496 * t361;
t228 = -pkin(10) * t372 + t263;
t212 = pkin(3) * t442 + t228;
t264 = t496 * t350 + t491 * t361;
t229 = -pkin(10) * t510 + t264;
t227 = t495 * t229;
t138 = t212 * t490 + t227;
t239 = t432 * t484 + t486 * t502;
t129 = qJ(5) * t432 + t138;
t349 = -t463 * pkin(2) - t401;
t279 = pkin(3) * t510 + t349;
t153 = t276 * pkin(4) - qJ(5) * t502 + t279;
t81 = -t129 * t484 + t486 * t153;
t61 = pkin(5) * t276 - pkin(11) * t239 + t81;
t552 = t486 * t432 - t484 * t502;
t82 = t486 * t129 + t484 * t153;
t62 = pkin(11) * t552 + t82;
t19 = -t489 * t62 + t494 * t61;
t20 = t489 * t61 + t494 * t62;
t801 = -t279 * mrSges(5,1) - t81 * mrSges(6,1) - t19 * mrSges(7,1) + t82 * mrSges(6,2) + t20 * mrSges(7,2) + t138 * mrSges(5,3) + t827;
t591 = qJD(1) * qJD(2);
t408 = (qJDD(1) * t492 + t497 * t591) * t485;
t589 = qJDD(1) * t487;
t462 = qJDD(2) + t589;
t504 = t510 * qJD(3);
t267 = t408 * t496 + t462 * t491 - t504;
t268 = -qJD(3) * t372 - t408 * t491 + t462 * t496;
t145 = -qJD(4) * t276 + t495 * t267 + t490 * t268;
t561 = t485 * t591;
t590 = qJDD(1) * t485;
t407 = -t492 * t561 + t497 * t590;
t393 = qJDD(3) - t407;
t385 = qJDD(4) + t393;
t120 = -t145 * t484 + t385 * t486;
t718 = t120 / 0.2e1;
t121 = t145 * t486 + t385 * t484;
t717 = t121 / 0.2e1;
t146 = qJD(4) * t502 + t490 * t267 - t495 * t268;
t713 = t146 / 0.2e1;
t764 = t484 * t823 + t486 * t822;
t763 = t484 * t822 - t486 * t823;
t472 = pkin(5) * t486 + pkin(4);
t816 = -m(6) * pkin(4) - m(7) * t472 + t538;
t815 = -t538 - t824;
t143 = qJDD(6) + t146;
t780 = -t239 * t489 + t494 * t552;
t55 = qJD(6) * t780 + t120 * t489 + t121 * t494;
t163 = t239 * t494 + t489 * t552;
t56 = -qJD(6) * t163 + t120 * t494 - t121 * t489;
t10 = Ifges(7,5) * t55 + Ifges(7,6) * t56 + Ifges(7,3) * t143;
t548 = qJD(2) * t587;
t581 = pkin(1) * t589;
t304 = pkin(8) * t407 + t492 * t581 + t497 * t548;
t283 = pkin(9) * t462 + t304;
t290 = -pkin(1) * t590 - pkin(2) * t407 - pkin(9) * t408;
t150 = -qJD(3) * t264 - t283 * t491 + t496 * t290;
t104 = pkin(3) * t393 - pkin(10) * t267 + t150;
t594 = qJD(3) * t496;
t149 = t496 * t283 + t491 * t290 - t350 * t595 + t361 * t594;
t112 = pkin(10) * t268 + t149;
t592 = qJD(4) * t495;
t593 = qJD(4) * t490;
t39 = t490 * t104 + t495 * t112 + t212 * t592 - t229 * t593;
t35 = qJ(5) * t385 + qJD(5) * t432 + t39;
t623 = t485 * t492;
t464 = pkin(8) * t623;
t305 = -qJD(2) * t458 - qJDD(1) * t464 - t492 * t548 + t497 * t581;
t284 = -pkin(2) * t462 - t305;
t205 = -pkin(3) * t268 + t284;
t60 = pkin(4) * t146 - qJ(5) * t145 - qJD(5) * t502 + t205;
t14 = -t35 * t484 + t486 * t60;
t15 = t486 * t35 + t484 * t60;
t682 = t385 / 0.2e1;
t714 = -t146 / 0.2e1;
t715 = t145 / 0.2e1;
t716 = t143 / 0.2e1;
t725 = t56 / 0.2e1;
t726 = t55 / 0.2e1;
t6 = pkin(5) * t146 - pkin(11) * t121 + t14;
t7 = pkin(11) * t120 + t15;
t2 = qJD(6) * t19 + t489 * t6 + t494 * t7;
t3 = -qJD(6) * t20 - t489 * t7 + t494 * t6;
t739 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t813 = t739 + 0.2e1 * Ifges(6,5) * t717 + Ifges(7,5) * t726 + 0.2e1 * Ifges(6,6) * t718 + Ifges(7,6) * t725 + 0.2e1 * Ifges(6,3) * t713 + Ifges(7,3) * t716 + t14 * mrSges(6,1) + t205 * mrSges(5,1) + t10 / 0.2e1 - mrSges(5,3) * t39 - t15 * mrSges(6,2) + (-t385 / 0.2e1 - t682) * Ifges(5,6) + (t713 - t714) * Ifges(5,2) + (-t145 / 0.2e1 - t715) * Ifges(5,4);
t690 = t276 / 0.2e1;
t273 = qJD(6) + t276;
t693 = t273 / 0.2e1;
t699 = t239 / 0.2e1;
t701 = t552 / 0.2e1;
t709 = t163 / 0.2e1;
t711 = t780 / 0.2e1;
t812 = Ifges(6,5) * t699 + Ifges(7,5) * t709 + Ifges(6,6) * t701 + Ifges(7,6) * t711 + Ifges(6,3) * t690 + Ifges(7,3) * t693 - t801 - t827;
t483 = qJ(3) + qJ(4);
t478 = sin(t483);
t479 = cos(t483);
t540 = -mrSges(4,1) * t496 + mrSges(4,2) * t491;
t811 = -m(4) * pkin(2) + t540 + (t816 + t824) * t479 + t512 * t478;
t126 = t239 * Ifges(6,4) + Ifges(6,2) * t552 + Ifges(6,6) * t276;
t810 = -t126 / 0.2e1;
t792 = t276 * t484;
t809 = pkin(5) * t792;
t808 = pkin(11) * t792;
t807 = pkin(5) * t826 + t820 * pkin(11) + t764;
t805 = pkin(11) * t804 - t763;
t356 = t447 * t490 - t448 * t495;
t757 = -qJD(4) * t356 + t829 * t490 + t495 * t830;
t226 = t490 * t229;
t137 = t212 * t495 - t226;
t128 = -pkin(4) * t432 + qJD(5) - t137;
t664 = mrSges(6,2) * t486;
t537 = t484 * mrSges(6,1) + t664;
t802 = t279 * mrSges(5,2) - t137 * mrSges(5,3) + t128 * t537 + t484 * t810;
t102 = -pkin(5) * t552 + t128;
t40 = t104 * t495 - t490 * t112 - t212 * t593 - t229 * t592;
t37 = -pkin(4) * t385 + qJDD(5) - t40;
t22 = -pkin(5) * t120 + t37;
t577 = Ifges(5,5) * t145 - Ifges(5,6) * t146 + Ifges(5,3) * t385;
t647 = t15 * t486;
t654 = Ifges(6,4) * t486;
t655 = Ifges(6,4) * t484;
t694 = -t273 / 0.2e1;
t710 = -t163 / 0.2e1;
t712 = -t780 / 0.2e1;
t159 = Ifges(7,4) * t780;
t85 = t163 * Ifges(7,1) + t273 * Ifges(7,5) + t159;
t719 = t85 / 0.2e1;
t720 = -t85 / 0.2e1;
t653 = Ifges(7,4) * t163;
t84 = Ifges(7,2) * t780 + t273 * Ifges(7,6) + t653;
t721 = t84 / 0.2e1;
t722 = -t84 / 0.2e1;
t727 = Ifges(6,1) * t717 + Ifges(6,4) * t718 + Ifges(6,5) * t713;
t48 = t121 * Ifges(6,4) + t120 * Ifges(6,2) + t146 * Ifges(6,6);
t728 = t48 / 0.2e1;
t729 = Ifges(7,1) * t726 + Ifges(7,4) * t725 + Ifges(7,5) * t716;
t730 = Ifges(7,4) * t726 + Ifges(7,2) * t725 + Ifges(7,6) * t716;
t796 = (-Ifges(7,4) * t409 - Ifges(7,2) * t410) * t711 + (-Ifges(7,1) * t409 - Ifges(7,4) * t410) * t709 + (-Ifges(7,5) * t409 - Ifges(7,6) * t410) * t693 + (Ifges(6,5) * t484 + Ifges(6,6) * t486) * t713 + t791 * t722 + (Ifges(7,1) * t790 + Ifges(7,4) * t791) * t710 + (Ifges(7,4) * t790 + Ifges(7,2) * t791) * t712 + (Ifges(7,5) * t790 + Ifges(7,6) * t791) * t694 + t790 * t720 + (t19 * t819 - t2 * t524 + t20 * t818 - t3 * t425) * mrSges(7,3) + (-mrSges(7,1) * t818 - mrSges(7,2) * t819) * t102 + t577 + t37 * t538 - t39 * mrSges(5,2) + t40 * mrSges(5,1) - t524 * t730 + (Ifges(7,5) * t425 - Ifges(7,6) * t524) * t716 + t22 * (mrSges(7,1) * t524 + mrSges(7,2) * t425) + (Ifges(7,1) * t425 - Ifges(7,4) * t524) * t726 + (Ifges(7,4) * t425 - Ifges(7,2) * t524) * t725 + (Ifges(6,1) * t484 + t654) * t717 + (Ifges(6,2) * t486 + t655) * t718 - t409 * t719 - t410 * t721 + t484 * t727 + t486 * t728 + t425 * t729 + mrSges(6,3) * t647;
t795 = m(7) * pkin(5);
t411 = t487 * t496 - t491 * t623;
t523 = pkin(3) * t411;
t643 = t239 * Ifges(6,5);
t644 = t552 * Ifges(6,6);
t125 = t276 * Ifges(6,3) + t643 + t644;
t642 = t273 * Ifges(7,3);
t645 = t163 * Ifges(7,5);
t646 = t780 * Ifges(7,6);
t83 = t642 + t645 + t646;
t793 = t125 + t83;
t756 = pkin(4) * t570 - t757;
t775 = Ifges(4,6) * t510;
t789 = Ifges(4,5) * t372 + Ifges(5,5) * t502 - t276 * Ifges(5,6) + Ifges(4,3) * t442 + t432 * Ifges(5,3) - t775;
t531 = -Ifges(6,2) * t484 + t654;
t533 = Ifges(6,1) * t486 - t655;
t127 = t239 * Ifges(6,1) + Ifges(6,4) * t552 + t276 * Ifges(6,5);
t619 = t486 * t127;
t272 = Ifges(5,4) * t276;
t639 = t432 * Ifges(5,5);
t191 = Ifges(5,1) * t502 - t272 + t639;
t704 = t191 / 0.2e1;
t788 = t531 * t701 + t533 * t699 + t619 / 0.2e1 + t704 + t802;
t676 = cos(qJ(1));
t571 = t676 * t497;
t493 = sin(qJ(1));
t617 = t492 * t493;
t416 = -t487 * t617 + t571;
t621 = t485 * t496;
t341 = -t416 * t491 + t493 * t621;
t787 = t305 * mrSges(3,1) - t304 * mrSges(3,2);
t679 = -t432 / 0.2e1;
t688 = -t502 / 0.2e1;
t700 = -t239 / 0.2e1;
t702 = -t552 / 0.2e1;
t786 = -Ifges(5,4) * t688 + Ifges(6,5) * t700 + Ifges(7,5) * t710 - Ifges(5,2) * t690 - Ifges(5,6) * t679 + Ifges(6,6) * t702 + Ifges(7,6) * t712 + Ifges(6,3) * t691 + Ifges(7,3) * t694 + t801;
t572 = t676 * t492;
t616 = t493 * t497;
t414 = t487 * t572 + t616;
t573 = t485 * t676;
t328 = t414 * t478 + t479 * t573;
t329 = t414 * t479 - t478 * t573;
t785 = t328 * t815 + t329 * t800;
t622 = t485 * t493;
t332 = t416 * t478 - t479 * t622;
t333 = t416 * t479 + t478 * t622;
t784 = t332 * t815 + t333 * t800;
t386 = t478 * t623 - t487 * t479;
t387 = t478 * t487 + t479 * t623;
t783 = t386 * t815 + t387 * t800;
t536 = t476 * mrSges(7,1) + t477 * mrSges(7,2);
t774 = -m(4) * pkin(9) - mrSges(4,3) - mrSges(5,3);
t782 = -t484 * t795 - t536 - t537 + t774;
t206 = pkin(4) * t502 + qJ(5) * t276;
t723 = Ifges(5,1) * t715 + Ifges(5,4) * t714 + Ifges(5,5) * t682;
t696 = t267 / 0.2e1;
t695 = t268 / 0.2e1;
t681 = t393 / 0.2e1;
t778 = pkin(5) * t502;
t475 = pkin(3) * t496 + pkin(2);
t317 = pkin(4) * t426 - qJ(5) * t427 - t475;
t243 = t486 * t317 - t356 * t484;
t627 = t427 * t486;
t208 = pkin(5) * t426 - pkin(11) * t627 + t243;
t244 = t484 * t317 + t486 * t356;
t628 = t427 * t484;
t224 = -pkin(11) * t628 + t244;
t131 = t208 * t489 + t224 * t494;
t777 = -qJD(6) * t131 + t489 * t805 + t494 * t807;
t130 = t208 * t494 - t224 * t489;
t776 = qJD(6) * t130 + t489 * t807 - t494 * t805;
t672 = pkin(3) * t490;
t471 = qJ(5) + t672;
t421 = (-pkin(11) - t471) * t484;
t480 = t486 * pkin(11);
t625 = t471 * t486;
t422 = t480 + t625;
t308 = t421 * t494 - t422 * t489;
t584 = pkin(3) * t592;
t461 = qJD(5) + t584;
t632 = t276 * t486;
t148 = t228 * t495 - t226;
t673 = pkin(3) * t372;
t176 = t206 + t673;
t92 = -t148 * t484 + t486 * t176;
t66 = pkin(11) * t632 + t778 + t92;
t93 = t486 * t148 + t484 * t176;
t80 = t93 + t808;
t773 = qJD(6) * t308 - t461 * t524 - t489 * t66 - t494 * t80;
t309 = t421 * t489 + t422 * t494;
t772 = -qJD(6) * t309 - t425 * t461 + t489 * t80 - t494 * t66;
t437 = t488 * t484;
t637 = qJ(5) * t486;
t438 = t480 + t637;
t344 = t437 * t494 - t438 * t489;
t97 = -t137 * t484 + t486 * t206;
t70 = t276 * t480 + t778 + t97;
t98 = t486 * t137 + t484 * t206;
t87 = t98 + t808;
t771 = -qJD(5) * t524 + qJD(6) * t344 - t489 * t70 - t494 * t87;
t345 = t437 * t489 + t438 * t494;
t770 = -qJD(5) * t425 - qJD(6) * t345 + t489 * t87 - t494 * t70;
t193 = t335 * t524 - t410 * t427;
t216 = t294 * t489 + t295 * t494;
t762 = t193 - t216;
t194 = t335 * t425 + t409 * t427;
t215 = t294 * t494 - t295 * t489;
t761 = t194 - t215;
t760 = pkin(5) * t804 + t756;
t164 = -mrSges(6,1) * t552 + mrSges(6,2) * t239;
t249 = mrSges(5,1) * t432 - mrSges(5,3) * t502;
t759 = t249 - t164;
t550 = mrSges(3,3) * t570;
t755 = -mrSges(3,1) * t463 + mrSges(4,1) * t510 + t372 * mrSges(4,2) + t550;
t620 = t485 * t497;
t420 = t487 * t492 * pkin(1) + pkin(8) * t620;
t389 = pkin(9) * t487 + t420;
t287 = t496 * t389 + t491 * t390;
t674 = pkin(1) * t497;
t419 = t487 * t674 - t464;
t750 = t149 * t496 - t150 * t491;
t747 = m(5) + m(7) + m(6);
t740 = -t484 * (mrSges(6,1) + t795) + mrSges(3,2) - t664 + t774;
t738 = mrSges(3,2) + t782;
t737 = mrSges(3,1) - t811;
t732 = t485 ^ 2;
t708 = Ifges(4,4) * t696 + Ifges(4,2) * t695 + Ifges(4,6) * t681;
t707 = Ifges(4,1) * t696 + Ifges(4,4) * t695 + Ifges(4,5) * t681;
t659 = Ifges(4,4) * t372;
t259 = -Ifges(4,2) * t510 + Ifges(4,6) * t442 + t659;
t698 = t259 / 0.2e1;
t363 = Ifges(4,4) * t510;
t260 = t372 * Ifges(4,1) + t442 * Ifges(4,5) - t363;
t697 = t260 / 0.2e1;
t683 = t372 / 0.2e1;
t677 = t487 / 0.2e1;
t671 = pkin(3) * t495;
t668 = qJD(3) / 0.2e1;
t597 = qJD(2) * t492;
t403 = qJD(2) * t522;
t405 = t419 * qJD(2);
t214 = -qJD(3) * t287 + t496 * t403 - t405 * t491;
t596 = qJD(2) * t497;
t567 = t485 * t596;
t339 = qJD(3) * t411 + t496 * t567;
t568 = t485 * t597;
t175 = pkin(3) * t568 - pkin(10) * t339 + t214;
t213 = -t389 * t595 + t390 * t594 + t491 * t403 + t496 * t405;
t412 = t487 * t491 + t492 * t621;
t338 = -qJD(3) * t412 - t491 * t567;
t185 = pkin(10) * t338 + t213;
t286 = -t389 * t491 + t496 * t390;
t235 = -pkin(3) * t620 - pkin(10) * t412 + t286;
t250 = pkin(10) * t411 + t287;
t68 = t490 * t175 + t495 * t185 + t235 * t592 - t250 * t593;
t64 = (qJ(5) * t597 - qJD(5) * t497) * t485 + t68;
t525 = t495 * t411 - t412 * t490;
t203 = qJD(4) * t525 + t338 * t490 + t339 * t495;
t297 = t411 * t490 + t412 * t495;
t204 = qJD(4) * t297 - t495 * t338 + t339 * t490;
t406 = t420 * qJD(2);
t282 = -pkin(3) * t338 + t406;
t91 = pkin(4) * t204 - qJ(5) * t203 - qJD(5) * t297 + t282;
t29 = t484 * t91 + t486 * t64;
t662 = mrSges(4,3) * t372;
t661 = Ifges(3,4) * t492;
t660 = Ifges(3,4) * t497;
t658 = Ifges(4,4) * t491;
t657 = Ifges(4,4) * t496;
t652 = Ifges(6,5) * t486;
t651 = Ifges(6,6) * t484;
t648 = t14 * t484;
t166 = t490 * t235 + t495 * t250;
t157 = -qJ(5) * t620 + t166;
t388 = t464 + (-pkin(2) - t674) * t487;
t310 = t388 - t523;
t188 = -pkin(4) * t525 - qJ(5) * t297 + t310;
t100 = t486 * t157 + t484 * t188;
t614 = -t328 * t472 - t329 * t488;
t613 = -t332 * t472 - t333 * t488;
t606 = -t386 * t472 - t387 * t488;
t600 = t676 * pkin(1) + pkin(8) * t622;
t579 = t491 * t622;
t576 = Ifges(4,5) * t267 + Ifges(4,6) * t268 + Ifges(4,3) * t393;
t575 = Ifges(3,5) * t408 + Ifges(3,6) * t407 + Ifges(3,3) * t462;
t16 = -t56 * mrSges(7,1) + t55 * mrSges(7,2);
t559 = -t599 / 0.2e1;
t558 = t599 / 0.2e1;
t556 = -pkin(1) * t493 + pkin(8) * t573;
t28 = -t484 * t64 + t486 * t91;
t67 = -t120 * mrSges(6,1) + t121 * mrSges(6,2);
t555 = -t328 * pkin(4) + t329 * qJ(5);
t554 = -t332 * pkin(4) + qJ(5) * t333;
t553 = -t386 * pkin(4) + qJ(5) * t387;
t99 = -t157 * t484 + t486 * t188;
t147 = t228 * t490 + t227;
t165 = t235 * t495 - t490 * t250;
t451 = t491 * t573;
t551 = -t414 * t496 + t451;
t549 = mrSges(3,3) * t569;
t543 = t341 * pkin(3);
t158 = pkin(4) * t620 - t165;
t542 = mrSges(3,1) * t492 + mrSges(3,2) * t497;
t541 = mrSges(4,1) * t411 - mrSges(4,2) * t412;
t534 = Ifges(4,1) * t496 - t658;
t530 = Ifges(4,5) * t496 - t491 * Ifges(4,6);
t529 = -t651 + t652;
t528 = t647 - t648;
t527 = -t484 * t81 + t486 * t82;
t271 = t297 * t486 - t484 * t620;
t73 = -pkin(5) * t525 - pkin(11) * t271 + t99;
t270 = -t297 * t484 - t486 * t620;
t88 = pkin(11) * t270 + t100;
t32 = -t489 * t88 + t494 * t73;
t33 = t489 * t73 + t494 * t88;
t199 = t270 * t494 - t271 * t489;
t200 = t270 * t489 + t271 * t494;
t69 = t175 * t495 - t490 * t185 - t235 * t593 - t250 * t592;
t518 = (t497 * Ifges(3,2) + t661) * t485;
t515 = -g(1) * t332 - g(2) * t328 - g(3) * t386;
t511 = t414 * t491 + t496 * t573;
t509 = t463 * t485 * (Ifges(3,5) * t497 - Ifges(3,6) * t492);
t508 = t492 * t732 * (Ifges(3,1) * t497 - t661);
t506 = t511 * pkin(3);
t505 = t510 * mrSges(4,3);
t503 = -mrSges(5,1) + t816;
t65 = -pkin(4) * t568 - t69;
t474 = -pkin(4) - t671;
t456 = Ifges(3,4) * t569;
t436 = -t472 - t671;
t417 = (-mrSges(3,1) * t497 + mrSges(3,2) * t492) * t485;
t415 = t487 * t616 + t572;
t413 = -t487 * t571 + t617;
t400 = -mrSges(3,2) * t463 + t549;
t347 = Ifges(3,1) * t570 + t463 * Ifges(3,5) + t456;
t346 = t463 * Ifges(3,6) + qJD(1) * t518;
t342 = t416 * t496 + t579;
t316 = mrSges(4,1) * t442 - t662;
t315 = -t442 * mrSges(4,2) - t505;
t303 = t524 * t427;
t302 = t425 * t427;
t293 = pkin(5) * t628 - t754;
t254 = t333 * t477 + t415 * t476;
t253 = -t333 * t476 + t415 * t477;
t248 = -mrSges(5,2) * t432 - mrSges(5,3) * t276;
t231 = -mrSges(4,2) * t393 + mrSges(4,3) * t268;
t230 = mrSges(4,1) * t393 - mrSges(4,3) * t267;
t207 = mrSges(5,1) * t276 + mrSges(5,2) * t502;
t192 = -mrSges(4,1) * t268 + mrSges(4,2) * t267;
t187 = t203 * t486 + t484 * t568;
t186 = -t203 * t484 + t486 * t568;
t182 = mrSges(6,1) * t276 - mrSges(6,3) * t239;
t181 = -mrSges(6,2) * t276 + mrSges(6,3) * t552;
t123 = -mrSges(5,2) * t385 - mrSges(5,3) * t146;
t122 = mrSges(5,1) * t385 - mrSges(5,3) * t145;
t119 = mrSges(7,1) * t273 - mrSges(7,3) * t163;
t118 = -mrSges(7,2) * t273 + mrSges(7,3) * t780;
t116 = -pkin(5) * t270 + t158;
t108 = t147 - t809;
t105 = t138 - t809;
t96 = -mrSges(7,1) * t780 + mrSges(7,2) * t163;
t86 = mrSges(5,1) * t146 + mrSges(5,2) * t145;
t79 = -qJD(6) * t200 + t186 * t494 - t187 * t489;
t78 = qJD(6) * t199 + t186 * t489 + t187 * t494;
t72 = mrSges(6,1) * t146 - mrSges(6,3) * t121;
t71 = -mrSges(6,2) * t146 + mrSges(6,3) * t120;
t57 = -pkin(5) * t186 + t65;
t31 = -mrSges(7,2) * t143 + mrSges(7,3) * t56;
t30 = mrSges(7,1) * t143 - mrSges(7,3) * t55;
t23 = pkin(11) * t186 + t29;
t17 = pkin(5) * t204 - pkin(11) * t187 + t28;
t5 = -qJD(6) * t33 + t17 * t494 - t23 * t489;
t4 = qJD(6) * t32 + t17 * t489 + t23 * t494;
t1 = [(Ifges(5,1) * t297 - Ifges(5,5) * t620) * t715 + (Ifges(5,1) * t203 + Ifges(5,5) * t568) * t687 + (t485 * t347 + t732 * qJD(1) * (-Ifges(3,2) * t492 + t660)) * t596 / 0.2e1 + (-t138 * t568 + t203 * t279 + t205 * t297 + t39 * t620) * mrSges(5,2) + (Ifges(6,1) * t271 + Ifges(6,4) * t270) * t717 + (Ifges(6,1) * t187 + Ifges(6,4) * t186) * t699 + (-m(4) * (-pkin(2) * t414 + t556) - t551 * mrSges(4,1) - t511 * mrSges(4,2) + t493 * mrSges(2,1) + t676 * mrSges(2,2) - m(3) * t556 + t414 * mrSges(3,1) - mrSges(3,3) * t573 - (t503 + t749) * t329 + (t536 - t740) * t413 - t512 * t328 + t747 * (-pkin(3) * t451 - t413 * t498 + t414 * t475 - t556)) * g(1) + (Ifges(7,5) * t200 + Ifges(7,6) * t199) * t716 + (Ifges(7,5) * t78 + Ifges(7,6) * t79) * t693 + (Ifges(5,5) * t297 - Ifges(5,3) * t620) * t682 + (Ifges(5,5) * t203 + Ifges(5,3) * t568) * t678 + (-t14 * t271 + t15 * t270 + t186 * t82 - t187 * t81) * mrSges(6,3) + t755 * t406 + (-t19 * t78 + t199 * t2 + t20 * t79 - t200 * t3) * mrSges(7,3) + (-t676 * mrSges(2,1) - m(3) * t600 - t416 * mrSges(3,1) - t254 * mrSges(7,1) - t253 * mrSges(7,2) - m(4) * (pkin(2) * t416 + t600) - t342 * mrSges(4,1) - t341 * mrSges(4,2) + (-mrSges(3,3) * t485 + mrSges(2,2)) * t493 + t503 * t333 + t740 * t415 + t512 * t332 - t747 * (pkin(3) * t579 - t415 * t498 + t416 * t475 + t600)) * g(2) + t203 * t704 + t412 * t707 + t411 * t708 + (Ifges(7,4) * t78 + Ifges(7,2) * t79) * t711 + m(6) * (t100 * t15 + t128 * t65 + t14 * t99 + t158 * t37 + t28 * t81 + t29 * t82) + m(7) * (t102 * t57 + t116 * t22 + t19 * t5 + t2 * t33 + t20 * t4 + t3 * t32) + m(5) * (t137 * t69 + t138 * t68 + t165 * t40 + t166 * t39 + t205 * t310 + t279 * t282) + m(4) * (t149 * t287 + t150 * t286 + t213 * t264 + t214 * t263 + t284 * t388 + t349 * t406) + t286 * t230 + t287 * t231 + (Ifges(6,5) * t271 + Ifges(6,6) * t270) * t713 + (Ifges(6,5) * t187 + Ifges(6,6) * t186) * t690 + (Ifges(6,4) * t271 + Ifges(6,2) * t270) * t718 + (Ifges(6,4) * t187 + Ifges(6,2) * t186) * t701 + t40 * (-mrSges(5,1) * t620 - mrSges(5,3) * t297) + t150 * (-mrSges(4,1) * t620 - mrSges(4,3) * t412) + t149 * (mrSges(4,2) * t620 + mrSges(4,3) * t411) + (mrSges(3,1) * t407 - mrSges(3,2) * t408 - qJDD(1) * t417 - t542 * t561) * t485 * pkin(1) + (Ifges(7,1) * t200 + Ifges(7,4) * t199) * t726 + (Ifges(7,1) * t78 + Ifges(7,4) * t79) * t709 + (t407 * t420 - t408 * t419 + (t304 * t497 - t305 * t492 + (-t401 * t497 - t404 * t492) * qJD(2)) * t485) * mrSges(3,3) + t508 * t591 / 0.2e1 + (Ifges(5,4) * t297 - Ifges(5,6) * t620) * t714 + (Ifges(5,4) * t203 + Ifges(5,6) * t568) * t691 - t346 * t568 / 0.2e1 + t264 * (-mrSges(4,2) * t568 + mrSges(4,3) * t338) + t442 * (Ifges(4,5) * t339 + Ifges(4,6) * t338 + Ifges(4,3) * t568) / 0.2e1 - t510 * (Ifges(4,4) * t339 + Ifges(4,2) * t338 + Ifges(4,6) * t568) / 0.2e1 + t137 * (mrSges(5,1) * t568 - mrSges(5,3) * t203) + t263 * (mrSges(4,1) * t568 - mrSges(4,3) * t339) + (Ifges(3,3) * t677 + t419 * mrSges(3,1) - t420 * mrSges(3,2) + (Ifges(3,5) * t492 + Ifges(3,6) * t497) * t485) * t462 + (Ifges(3,6) * t677 + t518) * t407 + (Ifges(3,5) * t677 + (t492 * Ifges(3,1) + t660) * t485) * t408 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t732 + t304 * t420 + t305 * t419 - t401 * t406 + t404 * t405) + t282 * t207 + (Ifges(7,4) * t200 + Ifges(7,2) * t199) * t725 + t37 * (-mrSges(6,1) * t270 + mrSges(6,2) * t271) + t68 * t248 + t69 * t249 - t813 * t525 + (t793 / 0.2e1 + t812) * t204 + Ifges(2,3) * qJDD(1) + t22 * (-mrSges(7,1) * t199 + mrSges(7,2) * t200) + t186 * t126 / 0.2e1 + t128 * (-mrSges(6,1) * t186 + mrSges(6,2) * t187) + t187 * t127 / 0.2e1 + t29 * t181 + t28 * t182 + t65 * t164 + t165 * t122 + t166 * t123 + t158 * t67 + t5 * t119 + t116 * t16 + t4 * t118 + t99 * t72 + t100 * t71 + t102 * (-mrSges(7,1) * t79 + mrSges(7,2) * t78) + t57 * t96 - t284 * t541 + qJD(2) * t509 / 0.2e1 - (t576 + t577) * t620 / 0.2e1 + t310 * t86 + t787 * t487 + t213 * t315 + t214 * t316 + t789 * t568 / 0.2e1 + t349 * (-mrSges(4,1) * t338 + mrSges(4,2) * t339) + t78 * t719 + t79 * t721 + t297 * t723 + t271 * t727 + t270 * t728 + t200 * t729 + t199 * t730 + t388 * t192 + t405 * t400 + t32 * t30 + t33 * t31 + t575 * t677 + (Ifges(4,5) * t412 + Ifges(4,6) * t411 - Ifges(4,3) * t620) * t681 + (Ifges(4,1) * t339 + Ifges(4,4) * t338 + Ifges(4,5) * t568) * t683 + (Ifges(4,4) * t412 + Ifges(4,2) * t411 - Ifges(4,6) * t620) * t695 + (Ifges(4,1) * t412 + Ifges(4,4) * t411 - Ifges(4,5) * t620) * t696 + t339 * t697 + t338 * t698; (Ifges(7,1) * t193 + Ifges(7,4) * t194) * t709 + (Ifges(7,1) * t216 + Ifges(7,4) * t215) * t710 + (t530 * t668 + (Ifges(4,3) * t492 + t497 * t530) * t559 + t349 * (mrSges(4,1) * t491 + mrSges(4,2) * t496)) * t442 + (-mrSges(7,1) * t761 + mrSges(7,2) * t762) * t102 + (-t19 * t762 - t2 * t302 + t20 * t761 + t3 * t303) * mrSges(7,3) + t763 * t181 + t764 * t182 + t294 * t810 + t776 * t118 + (t102 * t760 + t130 * t3 + t131 * t2 + t19 * t777 + t20 * t776 + t22 * t293) * m(7) + t777 * t119 + t352 * t704 + (Ifges(6,4) * t295 + Ifges(6,2) * t294) * t702 + (m(4) * ((-t263 * t496 - t264 * t491) * qJD(3) + t750) - t316 * t594 - t315 * t595 - t491 * t230 + t496 * t231) * pkin(9) + t751 * t207 - t137 * (mrSges(5,1) * t570 + mrSges(5,3) * t352) + (-Ifges(5,4) * t352 + Ifges(5,6) * t570) * t690 + (t138 * t570 + t279 * t352) * mrSges(5,2) + (-m(4) * t349 + t550 - t755) * t404 + t756 * t164 + t757 * t249 + t758 * t248 + t760 * t96 + (t205 * mrSges(5,2) - mrSges(5,3) * t40 + t37 * t537 + t529 * t713 + t531 * t718 + t533 * t717 + 0.2e1 * t723) * t427 + (-t747 * (-t415 * t475 - t416 * t498) + t738 * t416 + t737 * t415) * g(1) + (-t747 * (-t413 * t475 - t414 * t498) + t738 * t414 + t737 * t413) * g(2) + (t549 - t400) * t401 + (-Ifges(7,1) * t303 - Ifges(7,4) * t302) * t726 + (-Ifges(7,5) * t303 - Ifges(7,6) * t302) * t716 + t22 * (mrSges(7,1) * t302 - mrSges(7,2) * t303) + (-Ifges(7,4) * t303 - Ifges(7,2) * t302) * t725 + (-Ifges(5,1) * t352 + Ifges(5,5) * t570) * t688 + (-Ifges(5,5) * t352 + Ifges(5,3) * t570) * t679 + (Ifges(6,1) * t295 + Ifges(6,4) * t294) * t700 + t491 * t707 + t496 * t708 + t293 * t16 - t128 * (-mrSges(6,1) * t294 + mrSges(6,2) * t295) - t295 * t127 / 0.2e1 + (Ifges(7,5) * t193 + Ifges(7,6) * t194) * t693 + (Ifges(7,5) * t216 + Ifges(7,6) * t215) * t694 + (t128 * t756 + t14 * t243 + t15 * t244 - t37 * t754 + t763 * t82 + t764 * t81) * m(6) + (t137 * t757 + t138 * t758 - t205 * t475 + t279 * t751 + t356 * t39 + t40 * t754) * m(5) - (t67 - t122) * t754 + t372 * t534 * t668 + (-t14 * t627 - t15 * t628 - t804 * t82 + t81 * t820) * mrSges(6,3) + t787 + t575 + (t510 * t497 * t558 - t504 / 0.2e1) * (-Ifges(4,2) * t491 + t657) + (Ifges(6,5) * t295 + Ifges(6,6) * t294) * t691 + (-t509 / 0.2e1 + (-t508 / 0.2e1 + t732 * pkin(1) * t542) * qJD(1)) * qJD(1) + ((-mrSges(4,1) * t263 + mrSges(4,2) * t264) * t599 + (t346 + t775) * t558) * t492 + t243 * t72 + t244 * t71 + t813 * t426 + t812 * t336 - pkin(2) * t192 + (t417 - t747 * t475 * t620 + (t811 * t497 + (t498 * t747 + t782) * t492) * t485) * g(3) + t131 * t31 + t130 * t30 + t284 * t540 - t475 * t86 - t48 * t628 / 0.2e1 + t786 * t351 + (t828 * t264 + (t599 * t615 - t594) * t263 + t750) * mrSges(4,3) - t292 * t315 - t291 * t316 - (Ifges(5,1) * t687 + Ifges(5,4) * t691 + Ifges(5,5) * t678 + t529 * t690 + t788) * t335 + (t372 * (Ifges(4,5) * t492 + t497 * t534) + (-Ifges(3,2) * t570 + t496 * t260 + t347 + t456) * t497 + t789 * t492) * t559 + t356 * t123 + t793 * (-t351 / 0.2e1 + t336 / 0.2e1) + (-pkin(2) * t284 - t263 * t291 - t264 * t292) * m(4) + t193 * t719 + t216 * t720 + t194 * t721 + t215 * t722 + t627 * t727 - t303 * t729 - t302 * t730 + (Ifges(7,4) * t193 + Ifges(7,2) * t194) * t711 + (Ifges(7,4) * t216 + Ifges(7,2) * t215) * t712 - t259 * t595 / 0.2e1 + (Ifges(4,5) * t491 + Ifges(4,6) * t496) * t681 + (Ifges(4,2) * t496 + t658) * t695 + (Ifges(4,1) * t491 + t657) * t696 + t594 * t697 + t547 * t698; t772 * t119 + (-t102 * t108 + t19 * t772 + t2 * t309 + t20 * t773 + t22 * t436 + t3 * t308) * m(7) + t773 * t118 + t796 + (-t461 * t484 - t92) * t182 + (t461 * t486 - t93) * t181 + t510 * (-Ifges(4,2) * t372 - t363) / 0.2e1 + t759 * t147 + (t316 + t662) * t264 + (t137 * t147 - t138 * t148 - t279 * t673 + (t39 * t490 + t40 * t495 + (-t137 * t490 + t138 * t495) * qJD(4)) * pkin(3)) * m(5) + (-t148 + t584) * t248 + (-t632 * t81 - t792 * t82 - t648) * mrSges(6,3) + (-t128 * t147 + t37 * t474 + t461 * t527 + t471 * t528 - t81 * t92 - t82 * t93) * m(6) - t207 * t673 + t576 + (-t315 - t505) * t263 - t372 * (-Ifges(4,1) * t510 - t659) / 0.2e1 + (m(6) * t128 + m(7) * t102 - t759 + t96) * pkin(3) * t593 + (t619 + t191) * t690 - t149 * mrSges(4,2) + t150 * mrSges(4,1) - t108 * t96 - t484 * t471 * t72 - t442 * (-Ifges(4,5) * t510 - Ifges(4,6) * t372) / 0.2e1 - t349 * (t372 * mrSges(4,1) - mrSges(4,2) * t510) + t474 * t67 + (-m(7) * (t523 + t606) - m(6) * (t523 + t553) - m(5) * t523 - t541 + t783) * g(3) + t308 * t30 + t309 * t31 + (-m(5) * t543 - m(7) * (t543 + t613) - m(6) * (t543 + t554) - mrSges(4,1) * t341 + mrSges(4,2) * t342 + t784) * g(1) + (mrSges(4,1) * t511 - mrSges(4,2) * t551 + m(5) * t506 - m(7) * (-t506 + t614) - m(6) * (-t506 + t555) + t785) * g(2) + t786 * t502 + t793 * t688 - (Ifges(5,1) * t688 + Ifges(5,4) * t690 + Ifges(5,5) * t679 + t529 * t691 + t531 * t702 + t533 * t700 - t802) * t276 + t436 * t16 + t71 * t625 + t122 * t671 + t123 * t672 + t259 * t683 + t510 * t697; t770 * t119 + (-t613 * g(1) - t614 * g(2) - t606 * g(3) - t102 * t105 + t19 * t770 + t2 * t345 + t20 * t771 - t22 * t472 + t3 * t344) * m(7) + t771 * t118 + t796 + t759 * t138 + (-t642 / 0.2e1 - t646 / 0.2e1 - t645 / 0.2e1 - t125 / 0.2e1 - t83 / 0.2e1 - t643 / 0.2e1 - t644 / 0.2e1 + t801 + t821) * t502 + (-pkin(4) * t37 - g(1) * t554 - g(2) * t555 - g(3) * t553 + qJ(5) * t528 + qJD(5) * t527 - t128 * t138 - t81 * t97 - t82 * t98) * m(6) + (-t14 * mrSges(6,3) - qJ(5) * t72 - qJD(5) * t182) * t484 + (qJD(5) * t486 - t98) * t181 - t137 * t248 - t97 * t182 - t105 * t96 - pkin(4) * t67 - t472 * t16 + t783 * g(3) + t784 * g(1) + t785 * g(2) + (-t272 / 0.2e1 + t639 / 0.2e1 + (t652 / 0.2e1 - t651 / 0.2e1) * t276 + (-t484 * t82 - t486 * t81) * mrSges(6,3) + (-Ifges(5,2) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(6,3) / 0.2e1) * t502 + t788) * t276 + t344 * t30 + t345 * t31 + t71 * t637; -t780 * t118 + t163 * t119 - t552 * t181 + t239 * t182 + t16 + t67 + (t163 * t19 - t20 * t780 + t22 + t515) * m(7) + (t239 * t81 - t552 * t82 + t37 + t515) * m(6); -t102 * (mrSges(7,1) * t163 + mrSges(7,2) * t780) + (Ifges(7,1) * t780 - t653) * t710 + t84 * t709 + (Ifges(7,5) * t780 - Ifges(7,6) * t163) * t694 - t19 * t118 + t20 * t119 - g(1) * (mrSges(7,1) * t253 - mrSges(7,2) * t254) - g(2) * ((-t329 * t476 + t413 * t477) * mrSges(7,1) + (-t329 * t477 - t413 * t476) * mrSges(7,2)) - g(3) * ((-t387 * t476 - t477 * t620) * mrSges(7,1) + (-t387 * t477 + t476 * t620) * mrSges(7,2)) + (t163 * t20 + t19 * t780) * mrSges(7,3) + t10 + (-Ifges(7,2) * t163 + t159 + t85) * t712 + t739;];
tau  = t1;
