% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:16:00
% EndTime: 2019-03-09 10:16:32
% DurationCPUTime: 16.68s
% Computational Cost: add. (46221->719), mult. (90924->992), div. (0->0), fcn. (107232->10), ass. (0->362)
t456 = sin(qJ(6));
t458 = cos(qJ(6));
t455 = sin(pkin(11));
t457 = sin(qJ(4));
t459 = cos(qJ(4));
t590 = cos(pkin(11));
t486 = t455 * t459 + t457 * t590;
t712 = -t455 * t457 + t590 * t459;
t367 = t456 * t712 + t458 * t486;
t510 = -t456 * t486 + t458 * t712;
t556 = Ifges(7,5) * t510 - Ifges(7,6) * t367;
t589 = sin(pkin(10));
t530 = t589 * pkin(2);
t448 = t530 + pkin(8);
t564 = qJ(5) + t448;
t424 = t564 * t457;
t425 = t564 * t459;
t342 = -t455 * t424 + t590 * t425;
t298 = pkin(9) * t712 + t342;
t714 = -t590 * t424 - t455 * t425;
t729 = -pkin(9) * t486 + t714;
t178 = t298 * t458 + t456 * t729;
t738 = -t298 * t456 + t458 * t729;
t764 = -t178 * mrSges(7,1) - t738 * mrSges(7,2);
t28 = t556 + t764;
t765 = t28 * qJD(6);
t591 = cos(pkin(10));
t639 = sin(qJ(2));
t640 = cos(qJ(2));
t428 = t589 * t639 - t591 * t640;
t430 = -t589 * t640 - t591 * t639;
t327 = t712 * t430;
t329 = t486 * t430;
t511 = t327 * t456 + t458 * t329;
t233 = t327 * t458 - t329 * t456;
t620 = Ifges(7,4) * t233;
t116 = Ifges(7,2) * t511 + t428 * Ifges(7,6) - t620;
t227 = Ifges(7,4) * t511;
t118 = -Ifges(7,1) * t233 + t428 * Ifges(7,5) + t227;
t187 = -mrSges(7,2) * t428 + mrSges(7,3) * t511;
t189 = mrSges(7,1) * t428 + mrSges(7,3) * t233;
t546 = t639 * pkin(7);
t434 = -qJ(3) * t639 - t546;
t548 = t640 * pkin(7);
t436 = qJ(3) * t640 + t548;
t384 = -t591 * t434 + t436 * t589;
t573 = t430 * t457;
t308 = -pkin(4) * t573 + t384;
t241 = -t329 * pkin(5) + t308;
t532 = t591 * pkin(2);
t450 = -t532 - pkin(3);
t433 = -t459 * pkin(4) + t450;
t387 = -pkin(5) * t712 + t433;
t690 = -t178 / 0.2e1;
t728 = Ifges(7,2) * t233 + t227;
t733 = t233 * mrSges(7,1);
t742 = t511 * mrSges(7,2) - t733;
t718 = t510 * mrSges(7,2);
t732 = t367 * mrSges(7,1);
t743 = t732 + t718;
t744 = Ifges(7,1) * t511 + t620;
t762 = t738 / 0.2e1;
t763 = t189 * t690 + (t728 / 0.4e1 + t118 / 0.4e1) * t510 + t187 * t762 + t742 * t387 / 0.2e1 + t743 * t241 / 0.2e1 + (-t116 / 0.4e1 + t744 / 0.4e1) * t367;
t451 = -pkin(2) * t640 - pkin(1);
t362 = t428 * pkin(3) + t430 * pkin(8) + t451;
t713 = t589 * t434 + t591 * t436;
t272 = t459 * t362 - t457 * t713;
t572 = t430 * t459;
t228 = qJ(5) * t572 + t272;
t190 = t428 * pkin(4) + t228;
t273 = t457 * t362 + t459 * t713;
t229 = qJ(5) * t573 + t273;
t203 = t455 * t229;
t111 = t590 * t190 - t203;
t636 = pkin(9) * t327;
t72 = pkin(5) * t428 + t111 + t636;
t516 = t590 * t229;
t112 = t455 * t190 + t516;
t635 = pkin(9) * t329;
t75 = t112 + t635;
t50 = t456 * t72 + t458 * t75;
t120 = -t228 * t455 - t516;
t81 = t120 - t635;
t121 = t590 * t228 - t203;
t82 = t121 + t636;
t59 = -t456 * t82 + t458 * t81;
t761 = t50 + t59;
t597 = t367 * mrSges(7,3);
t758 = t241 * t742;
t756 = t387 * t743;
t755 = t111 - t121;
t547 = t639 * pkin(2);
t363 = -t430 * pkin(3) + t428 * pkin(8) + t547;
t274 = t459 * t363 + t384 * t457;
t574 = t428 * t459;
t191 = -t430 * pkin(4) + qJ(5) * t574 + t274;
t275 = t457 * t363 - t384 * t459;
t575 = t428 * t457;
t239 = qJ(5) * t575 + t275;
t113 = t590 * t191 - t239 * t455;
t114 = t455 * t191 + t590 * t239;
t328 = t486 * t428;
t330 = t712 * t428;
t234 = t328 * t458 + t330 * t456;
t186 = mrSges(7,2) * t430 + mrSges(7,3) * t234;
t237 = t328 * t456 - t330 * t458;
t188 = -mrSges(7,1) * t430 - mrSges(7,3) * t237;
t281 = mrSges(6,2) * t430 + mrSges(6,3) * t328;
t283 = -mrSges(6,1) * t430 + mrSges(6,3) * t330;
t531 = t590 * pkin(4);
t449 = t531 + pkin(5);
t638 = pkin(4) * t455;
t412 = t449 * t458 - t456 * t638;
t413 = t449 * t456 + t458 * t638;
t670 = -t330 / 0.2e1;
t672 = t328 / 0.2e1;
t492 = Ifges(6,5) * t670 + Ifges(6,6) * t672;
t681 = t237 / 0.2e1;
t685 = t234 / 0.2e1;
t491 = Ifges(7,5) * t681 + Ifges(7,6) * t685;
t73 = -pkin(5) * t430 + pkin(9) * t330 + t113;
t76 = pkin(9) * t328 + t114;
t53 = -t456 * t76 + t458 * t73;
t54 = t456 * t73 + t458 * t76;
t646 = -t430 / 0.2e1;
t502 = Ifges(7,3) * t646 - t54 * mrSges(7,2) / 0.2e1 + t53 * mrSges(7,1) / 0.2e1 + t491;
t524 = -t574 / 0.2e1;
t525 = t575 / 0.2e1;
t695 = m(6) * pkin(4);
t551 = t695 / 0.2e1;
t653 = t412 / 0.2e1;
t697 = m(7) / 0.2e1;
t723 = Ifges(5,3) + Ifges(6,3);
t754 = (t412 * t53 + t413 * t54) * t697 + t113 * mrSges(6,1) / 0.2e1 - t114 * mrSges(6,2) / 0.2e1 + t274 * mrSges(5,1) / 0.2e1 - t275 * mrSges(5,2) / 0.2e1 + t188 * t653 + t413 * t186 / 0.2e1 + (t113 * t590 + t114 * t455) * t551 + Ifges(5,5) * t524 + Ifges(5,6) * t525 + t281 * t638 / 0.2e1 + t283 * t531 / 0.2e1 + t492 + t723 * t646 + t502;
t619 = Ifges(7,4) * t367;
t265 = Ifges(7,2) * t510 + t619;
t745 = Ifges(7,1) * t510 - t619;
t752 = t265 / 0.4e1 - t745 / 0.4e1;
t541 = -t732 / 0.2e1;
t687 = t733 / 0.2e1;
t750 = m(4) * t547;
t749 = m(6) * t433;
t715 = t112 + t120;
t282 = -mrSges(6,2) * t428 + t329 * mrSges(6,3);
t284 = mrSges(6,1) * t428 + t327 * mrSges(6,3);
t649 = t486 / 0.2e1;
t662 = t367 / 0.2e1;
t666 = -t510 / 0.2e1;
t475 = t189 * t662 + t187 * t666 - t712 * t282 / 0.2e1 + t284 * t649;
t49 = -t456 * t75 + t458 * t72;
t499 = -t367 * t49 + t50 * t510;
t369 = -mrSges(5,2) * t428 + mrSges(5,3) * t573;
t565 = t459 * t369;
t371 = t428 * mrSges(5,1) + mrSges(5,3) * t572;
t567 = t457 * t371;
t60 = t456 * t81 + t458 * t82;
t698 = m(6) / 0.2e1;
t748 = t475 - (-t486 * t755 + t715 * t712) * t698 - (t367 * t60 + t510 * t59 + t499) * t697 + t567 / 0.2e1 - t565 / 0.2e1;
t747 = t475 - (-t111 * t486 + t112 * t712 + t327 * t714 + t329 * t342) * t698 - (t178 * t511 + t233 * t738 + t499) * t697;
t606 = t237 * mrSges(7,2);
t609 = t234 * mrSges(7,1);
t563 = t609 / 0.2e1 - t606 / 0.2e1;
t602 = t330 * mrSges(6,2);
t603 = t328 * mrSges(6,1);
t746 = t563 + t602 / 0.2e1 + t603 / 0.2e1;
t719 = Ifges(7,5) * t511;
t735 = Ifges(7,6) * t233;
t562 = t719 + t735;
t361 = Ifges(7,4) * t510;
t268 = Ifges(7,1) * t367 + t361;
t727 = -Ifges(7,2) * t367 + t361;
t741 = t727 / 0.4e1 + t268 / 0.4e1;
t520 = -t735 / 0.2e1 - t719 / 0.2e1;
t740 = -t430 * mrSges(4,1) - t428 * mrSges(4,2);
t686 = t233 / 0.2e1;
t737 = -t367 / 0.2e1;
t594 = t486 * mrSges(6,3);
t706 = t457 ^ 2 + t459 ^ 2;
t731 = t428 * t706;
t725 = t510 / 0.2e1;
t724 = -t511 / 0.2e1;
t452 = Ifges(5,5) * t459;
t615 = Ifges(5,6) * t457;
t716 = Ifges(4,4) - t452 / 0.2e1 + t615 / 0.2e1;
t453 = Ifges(5,4) * t459;
t439 = Ifges(5,1) * t457 + t453;
t671 = -t329 / 0.2e1;
t673 = t327 / 0.2e1;
t711 = t342 * t673 + t671 * t714;
t710 = Ifges(6,5) * t712 - Ifges(6,6) * t486 + t556;
t709 = Ifges(6,5) * t329 + Ifges(6,6) * t327 + t562;
t500 = Ifges(5,2) * t457 - t453;
t355 = t439 * t430;
t643 = -t448 / 0.2e1;
t708 = t369 * t643 + t355 / 0.4e1;
t623 = Ifges(5,4) * t457;
t437 = Ifges(5,2) * t459 + t623;
t354 = t430 * t437;
t707 = t371 * t643 + t354 / 0.4e1;
t479 = -t718 + 0.2e1 * t541;
t705 = t479 * qJD(6);
t126 = -mrSges(7,1) * t511 - mrSges(7,2) * t233;
t622 = Ifges(6,4) * t327;
t195 = Ifges(6,2) * t329 + t428 * Ifges(6,6) - t622;
t242 = -t327 * mrSges(6,1) + t329 * mrSges(6,2);
t246 = Ifges(6,1) * t329 + t622;
t262 = -mrSges(7,1) * t510 + mrSges(7,2) * t367;
t549 = pkin(4) * t572;
t280 = -t327 * pkin(5) - t549;
t501 = -t459 * mrSges(5,1) + t457 * mrSges(5,2);
t351 = t430 * t501;
t372 = mrSges(6,1) * t486 + mrSges(6,2) * t712;
t422 = Ifges(6,4) * t712;
t374 = -Ifges(6,2) * t486 + t422;
t621 = Ifges(6,4) * t486;
t375 = Ifges(6,2) * t712 + t621;
t376 = Ifges(6,1) * t712 - t621;
t377 = Ifges(6,1) * t486 + t422;
t637 = pkin(4) * t457;
t393 = pkin(5) * t486 + t637;
t435 = t457 * mrSges(5,1) + t459 * mrSges(5,2);
t647 = t428 / 0.4e1;
t654 = t393 / 0.2e1;
t675 = t280 / 0.2e1;
t694 = -t49 / 0.2e1;
t702 = (t374 / 0.4e1 + t377 / 0.4e1) * t329 + t714 * t282 / 0.2e1 - (-t246 / 0.4e1 + t195 / 0.4e1) * t486 + t262 * t675 - t342 * t284 / 0.2e1 + t741 * t511 + (t375 / 0.4e1 - t376 / 0.4e1) * t327 + t126 * t654 + t752 * t233 + (t241 * t393 + t280 * t387 + t761 * t738 + (-t49 + t60) * t178) * t697 + (t178 * t686 + t510 * t694 + t60 * t725 + t738 * t724 + t737 * t761) * mrSges(7,3) + t450 * t351 / 0.2e1 + t433 * t242 / 0.2e1 + t384 * t435 / 0.2e1 + t308 * t372 / 0.2e1 + t710 * t647 + t763;
t699 = m(5) / 0.2e1;
t696 = m(4) * pkin(2);
t693 = t118 / 0.2e1;
t691 = -t738 / 0.2e1;
t683 = t511 / 0.2e1;
t680 = -t233 / 0.2e1;
t677 = t268 / 0.2e1;
t674 = -t327 / 0.2e1;
t660 = t375 / 0.2e1;
t657 = t377 / 0.2e1;
t652 = -t413 / 0.2e1;
t651 = t712 / 0.2e1;
t648 = t428 / 0.2e1;
t645 = t437 / 0.4e1;
t440 = Ifges(5,1) * t459 - t623;
t644 = -t440 / 0.4e1;
t642 = t457 / 0.2e1;
t641 = t459 / 0.2e1;
t633 = t49 * mrSges(7,2);
t632 = t50 * mrSges(7,1);
t629 = t59 * mrSges(7,1);
t628 = t60 * mrSges(7,2);
t627 = mrSges(4,3) * t430;
t626 = mrSges(7,3) * t412;
t625 = mrSges(7,3) * t413;
t618 = Ifges(5,5) * t428;
t616 = Ifges(5,6) * t428;
t115 = Ifges(7,4) * t237 + Ifges(7,2) * t234 - Ifges(7,6) * t430;
t117 = Ifges(7,1) * t237 + Ifges(7,4) * t234 - Ifges(7,5) * t430;
t125 = t606 - t609;
t194 = -Ifges(6,4) * t330 + Ifges(6,2) * t328 - Ifges(6,6) * t430;
t196 = -Ifges(6,1) * t330 + Ifges(6,4) * t328 - Ifges(6,5) * t430;
t322 = Ifges(6,4) * t329;
t197 = -Ifges(6,1) * t327 + Ifges(6,5) * t428 + t322;
t307 = -pkin(4) * t575 + t713;
t240 = -t328 * pkin(5) + t307;
t243 = -t602 - t603;
t244 = -mrSges(6,1) * t329 - mrSges(6,2) * t327;
t299 = -Ifges(5,6) * t430 + t428 * t500;
t301 = -Ifges(5,5) * t430 - t428 * t440;
t352 = t435 * t428;
t353 = t435 * t430;
t368 = mrSges(5,2) * t430 + mrSges(5,3) * t575;
t370 = -t430 * mrSges(5,1) + mrSges(5,3) * t574;
t302 = -t430 * t440 + t618;
t566 = t459 * t302;
t300 = t430 * t500 + t616;
t568 = t457 * t300;
t3 = -t713 * t353 + (Ifges(3,1) - Ifges(3,2)) * t640 * t639 + t118 * t681 + t115 * t683 + t116 * t685 + t195 * t672 + t196 * t674 + t117 * t680 + m(5) * (t272 * t274 + t273 * t275 + t384 * t713) + (t740 + t750) * t451 - pkin(1) * (mrSges(3,1) * t639 + mrSges(3,2) * t640) + t197 * t670 - t384 * t352 + t273 * t368 + t275 * t369 + t272 * t370 + t274 * t371 + t329 * t194 / 0.2e1 + t307 * t244 + t308 * t243 + t112 * t281 + t114 * t282 + t111 * t283 + t113 * t284 + t240 * t126 + t241 * t125 + t50 * t186 + t54 * t187 + t49 * t188 + t53 * t189 + (-t566 / 0.2e1 + t568 / 0.2e1 + mrSges(4,1) * t547 + t491 + t492 + t716 * t428) * t428 + m(6) * (t111 * t113 + t112 * t114 + t307 * t308) + m(7) * (t240 * t241 + t49 * t53 + t50 * t54) + (-t639 ^ 2 + t640 ^ 2) * Ifges(3,4) + (-t459 * t301 / 0.2e1 + t299 * t642 + Ifges(6,5) * t673 + Ifges(6,6) * t671 + Ifges(7,5) * t686 + Ifges(7,6) * t724 - mrSges(4,2) * t547 - t716 * t430 + (Ifges(4,1) - Ifges(4,2) - Ifges(7,3) - t723) * t428) * t430;
t604 = t3 * qJD(1);
t600 = t510 * mrSges(7,3);
t245 = Ifges(6,2) * t327 + t322;
t4 = t384 * t351 + t272 * t369 - t273 * t371 + t308 * t242 + t280 * t126 + t121 * t282 + t120 * t284 + t758 + t728 * t683 + t511 * t693 + t744 * t680 + t116 * t686 + t60 * t187 + t59 * t189 + (t233 * t50 - t49 * t511) * mrSges(7,3) + m(7) * (t241 * t280 + t49 * t59 + t50 * t60) + m(6) * (t111 * t120 + t112 * t121) + (-t111 * mrSges(6,3) + t245 / 0.2e1 + t197 / 0.2e1) * t329 + (t112 * mrSges(6,3) - t246 / 0.2e1 + t195 / 0.2e1) * t327 + ((t354 / 0.2e1 + t302 / 0.2e1 + t618 / 0.2e1 - t272 * mrSges(5,3)) * t457 + (-t355 / 0.2e1 + t300 / 0.2e1 + t616 / 0.2e1 + t273 * mrSges(5,3) + (-m(6) * t308 - t244) * pkin(4)) * t459) * t430 + t709 * t648;
t596 = t4 * qJD(1);
t595 = t712 * mrSges(6,3);
t593 = t428 * mrSges(4,3);
t7 = t49 * t187 - t50 * t189 + t562 * t648 + t758 - (-t50 * mrSges(7,3) + t744 / 0.2e1 - t116 / 0.2e1) * t233 + (-t49 * mrSges(7,3) + t693 + t728 / 0.2e1) * t511;
t592 = t7 * qJD(1);
t373 = -mrSges(6,1) * t712 + mrSges(6,2) * t486;
t537 = -t594 / 0.2e1;
t538 = t595 / 0.2e1;
t540 = -t597 / 0.2e1;
t542 = t600 / 0.2e1;
t461 = (-t450 * t430 - t448 * t731) * t699 + (t328 * t714 - t330 * t342 - t430 * t433) * t698 + (t178 * t237 + t234 * t738 - t387 * t430) * t697 - t351 / 0.2e1 + (-t428 * t589 + t430 * t591) * t696 / 0.2e1 + t234 * t540 + t237 * t542 + t328 * t537 - t330 * t538 + (t262 + t373) * t646 - mrSges(5,3) * t731 / 0.2e1;
t463 = (t274 * t459 + t457 * t275) * t699 + (t113 * t712 + t114 * t486) * t698 + (t367 * t54 + t510 * t53) * t697 + t188 * t725 + t186 * t662 + t283 * t651 + t281 * t649 + t368 * t642 + t370 * t641 + t750 / 0.2e1;
t12 = -t461 + t463 + t740;
t588 = qJD(1) * t12;
t22 = t511 * t187 + t233 * t189 + t329 * t282 + t327 * t284 + m(7) * (t233 * t49 + t50 * t511) + m(6) * (t111 * t327 + t112 * t329);
t587 = qJD(1) * t22;
t578 = t384 * t430;
t16 = t237 * t187 + t234 * t189 - t330 * t282 + t328 * t284 + (-t565 + t567 + t593) * t428 + (-t126 - t244 + t353 + t627) * t430 + m(7) * (t234 * t49 + t237 * t50 - t241 * t430) + m(6) * (t111 * t328 - t112 * t330 - t308 * t430) + m(5) * (-t578 + (t272 * t457 - t273 * t459) * t428) + m(4) * (-t428 * t713 - t578);
t586 = t16 * qJD(1);
t473 = (-t233 * t737 + t511 * t666) * mrSges(7,3) + t187 * t725 + t189 * t737;
t18 = t473 - t563;
t585 = t18 * qJD(1);
t582 = t308 * t457;
t571 = t448 * t457;
t570 = t448 * t459;
t171 = -mrSges(7,1) * t413 - mrSges(7,2) * t412;
t553 = t171 * qJD(6);
t552 = mrSges(6,3) * t638;
t545 = t637 / 0.2e1;
t543 = mrSges(5,3) * t448 / 0.2e1;
t539 = -t595 / 0.2e1;
t512 = t452 - t615;
t509 = -t549 / 0.2e1;
t508 = mrSges(6,3) * t531;
t505 = t329 * t539 - t594 * t674 + t597 * t686 + t600 * t724;
t504 = t233 * t540 + t327 * t537 + t329 * t538 + t511 * t542;
t27 = t756 + (-t265 / 0.2e1 + t745 / 0.2e1) * t367 + (t727 / 0.2e1 + t677) * t510;
t464 = (mrSges(7,3) * t691 + t741) * t511 - (mrSges(7,3) * t690 - t752) * t233 + t556 * t647 + t763;
t5 = t464 - t502;
t498 = t5 * qJD(1) + t27 * qJD(2);
t470 = (t233 * t412 + t413 * t511) * t697 + (t327 * t590 + t329 * t455) * t551;
t483 = m(6) * t509 + m(7) * t675;
t36 = t742 + t242 - t470 + t483;
t478 = (t455 * t712 - t486 * t590) * t695;
t488 = m(7) * (-t367 * t412 + t413 * t510);
t469 = t488 / 0.2e1 + t478 / 0.2e1;
t484 = t743 + t372;
t487 = m(6) * t545 + m(7) * t654;
t71 = -t469 + t484 + t487;
t497 = qJD(1) * t36 + qJD(2) * t71;
t471 = (t233 * t510 + t367 * t511) * t697 + (t327 * t712 + t329 * t486) * t698;
t65 = 0.2e1 * (m(6) / 0.4e1 + m(7) / 0.4e1) * t430 + t471;
t495 = qJD(1) * t65;
t69 = 0.2e1 * t724 * mrSges(7,2) + 0.2e1 * t687;
t494 = qJD(1) * t69 + qJD(2) * t479;
t476 = -t342 * t755 + t715 * t714;
t2 = t715 * t537 + t707 * t459 + t708 * t457 + t711 * mrSges(6,3) + t702 + ((-t433 * t572 + t582) * pkin(4) + t476) * t698 + t512 * t647 + t572 * t644 + t572 * t645 + t706 * t430 * t543 + (t439 - t500) * t573 / 0.4e1 + (t245 + t197) * t712 / 0.4e1 + t566 / 0.4e1 + t373 * t509 - t568 / 0.4e1 + t121 * t538 + t111 * t539 + t244 * t545 - t754;
t21 = t450 * t435 + t433 * t372 + t756 + t745 * t662 + t265 * t737 + (-t500 / 0.2e1 + t439 / 0.2e1) * t459 - (-t376 / 0.2e1 + t660) * t486 + (t374 / 0.2e1 + t657) * t712 + (pkin(4) * t373 - t437 / 0.2e1 + t440 / 0.2e1) * t457 + t637 * t749 + (t727 + t268) * t725 + (m(7) * t387 + t262) * t393;
t482 = t2 * qJD(1) + t21 * qJD(2);
t465 = (t234 * t412 + t237 * t413) * t697 + (t328 * t590 - t330 * t455) * t551 + mrSges(5,1) * t525 + mrSges(5,2) * t574 / 0.2e1 + t746;
t10 = t465 + t504 + t748;
t481 = t10 * qJD(1);
t474 = t240 * t697 + t307 * t698 - t746;
t14 = t474 + t505 + t747;
t39 = (t367 ^ 2 + t510 ^ 2) * mrSges(7,3) + (t486 ^ 2 + t712 ^ 2) * mrSges(6,3) + m(7) * (t178 * t510 - t367 * t738) + m(6) * (t342 * t712 - t486 * t714);
t480 = -qJD(1) * t14 + qJD(2) * t39;
t29 = (t691 + t762) * mrSges(7,2) + (t690 + t178 / 0.2e1) * mrSges(7,1);
t68 = t541 + t732 / 0.2e1 + (t666 + t725) * mrSges(7,2);
t468 = (-t233 * t652 + t412 * t724) * mrSges(7,3) + t187 * t653 + t189 * t652 - t520;
t8 = (t694 + t60 / 0.2e1) * mrSges(7,2) + (-t50 / 0.2e1 - t59 / 0.2e1) * mrSges(7,1) + t468 + t520;
t477 = t8 * qJD(1) + t29 * qJD(2) + t68 * qJD(3) + t171 * qJD(4);
t122 = t469 + t487;
t70 = -t733 / 0.2e1 + t687;
t66 = t470 + t483;
t64 = t471 + (m(6) + m(7)) * t646;
t17 = t473 + t563;
t15 = t474 + t504 - t747;
t13 = t461 + t463;
t11 = t465 + t505 - t748;
t9 = -t633 / 0.2e1 - t632 / 0.2e1 - t628 / 0.2e1 + t629 / 0.2e1 + t468 - t520;
t6 = t464 + t502;
t1 = (pkin(4) * t244 / 0.2e1 - t616 / 0.4e1 - t300 / 0.4e1 + t708) * t457 + (pkin(4) * t582 + t476) * t698 + (t302 / 0.4e1 + t707) * t459 + (t245 / 0.4e1 + t197 / 0.4e1) * t712 + t452 * t647 + ((t439 / 0.4e1 - t500 / 0.4e1 + t457 * t543) * t457 + (t644 + t645 + t459 * t543 + (-t749 / 0.2e1 - t373 / 0.2e1) * pkin(4)) * t459) * t430 + (-(t120 / 0.2e1 + t112 / 0.2e1) * t486 + (t121 / 0.2e1 - t111 / 0.2e1) * t712 + t711) * mrSges(6,3) + t702 + t754;
t19 = [qJD(2) * t3 + qJD(3) * t16 + qJD(4) * t4 + qJD(5) * t22 + qJD(6) * t7, t13 * qJD(3) + t1 * qJD(4) + t15 * qJD(5) + t6 * qJD(6) + t604 + (m(6) * (t113 * t714 + t114 * t342 + t307 * t433) + t714 * t283 + (m(5) * t450 - t591 * t696 - mrSges(4,1) + t501) * t713 + (-t589 * t696 + mrSges(4,2)) * t384 + t115 * t725 + t265 * t685 + t237 * t677 + (Ifges(5,5) * t457 + Ifges(6,5) * t486 + Ifges(7,5) * t367 + Ifges(5,6) * t459 + Ifges(6,6) * t712 + Ifges(7,6) * t510) * t646 + (m(5) * t448 + mrSges(5,3)) * (-t274 * t457 + t275 * t459) - Ifges(3,6) * t639 + Ifges(3,5) * t640 - t53 * t597 - t370 * t571 - mrSges(3,1) * t548 + t117 * t662 + t194 * t651 - t330 * t657 + t328 * t660 + t196 * t649 + t299 * t641 + t301 * t642 + t530 * t627 - t113 * t594 + t114 * t595 + t54 * t600 + m(7) * (t178 * t54 + t240 * t387 + t53 * t738) + t738 * t188 + t532 * t593 + mrSges(3,2) * t546 + t368 * t570 - t450 * t352 + t433 * t243 - Ifges(4,5) * t428 + Ifges(4,6) * t430 + t387 * t125 + t307 * t373 + t342 * t281 + t240 * t262 + t178 * t186 + t439 * t524 + t437 * t525) * qJD(2), t586 + t13 * qJD(2) + 0.2e1 * ((t328 * t712 - t330 * t486) * t698 + (t234 * t510 + t237 * t367) * t697) * qJD(3) + t11 * qJD(4) + t64 * qJD(5) + t17 * qJD(6), t596 + t1 * qJD(2) + t11 * qJD(3) + (Ifges(5,5) * t573 + Ifges(5,6) * t572 + m(7) * (t412 * t59 + t413 * t60) + t233 * t625 - t511 * t626 - t628 + t629 + (t120 * t590 + t121 * t455) * t695 + t120 * mrSges(6,1) - t121 * mrSges(6,2) - t273 * mrSges(5,1) - t272 * mrSges(5,2) + t327 * t552 - t329 * t508 + t709) * qJD(4) + t66 * qJD(5) + t9 * qJD(6), qJD(2) * t15 + qJD(3) * t64 + qJD(4) * t66 + qJD(6) * t70 + t587, t592 + t6 * qJD(2) + t17 * qJD(3) + t9 * qJD(4) + t70 * qJD(5) + (t562 - t632 - t633) * qJD(6); -qJD(3) * t12 + qJD(4) * t2 - qJD(5) * t14 + qJD(6) * t5 - t604, qJD(4) * t21 + qJD(5) * t39 + qJD(6) * t27, -t588 (m(7) * (-t178 * t412 + t413 * t738) - t714 * mrSges(6,2) - t342 * mrSges(6,1) + (-t342 * t590 + t455 * t714) * t695 - t486 * t552 - t712 * t508 - t367 * t625 - t510 * t626 + mrSges(5,2) * t571 - mrSges(5,1) * t570 + t512 + t710 + t764) * qJD(4) + t122 * qJD(5) + t765 + t482, qJD(4) * t122 + t480, t28 * qJD(4) + t498 + t765; qJD(2) * t12 - qJD(4) * t10 + qJD(5) * t65 + qJD(6) * t18 - t586, t588, 0 (-t435 + t478 - t484 + t488) * qJD(4) + t705 - t481, t495, qJD(4) * t479 - qJD(6) * t743 + t585; -qJD(2) * t2 + qJD(3) * t10 - qJD(5) * t36 + qJD(6) * t8 - t596, -qJD(5) * t71 + qJD(6) * t29 - t482, qJD(6) * t68 + t481, t553, -t497, t477 + t553; qJD(2) * t14 - qJD(3) * t65 + qJD(4) * t36 - qJD(6) * t69 - t587, qJD(4) * t71 - t480 - t705, -t495, t497, 0, -t494; -qJD(2) * t5 - qJD(3) * t18 - qJD(4) * t8 + qJD(5) * t69 - t592, -t29 * qJD(4) + qJD(5) * t479 - t498, -t68 * qJD(4) - t585, -t477, t494, 0;];
Cq  = t19;
