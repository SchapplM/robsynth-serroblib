% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:52:42
% EndTime: 2019-03-09 11:53:10
% DurationCPUTime: 15.24s
% Computational Cost: add. (33980->756), mult. (66444->982), div. (0->0), fcn. (74486->8), ass. (0->379)
t445 = sin(qJ(4));
t446 = cos(qJ(4));
t665 = sin(qJ(5));
t667 = cos(qJ(5));
t414 = t665 * t445 - t667 * t446;
t415 = -t445 * t667 - t446 * t665;
t601 = sin(pkin(10));
t602 = cos(pkin(10));
t666 = sin(qJ(2));
t668 = cos(qJ(2));
t411 = -t601 * t668 - t602 * t666;
t299 = t414 * t411;
t409 = t601 * t666 - t602 * t668;
t231 = mrSges(6,1) * t409 - t299 * mrSges(6,3);
t301 = t415 * t411;
t622 = t301 * mrSges(6,3);
t228 = -mrSges(6,2) * t409 - t622;
t610 = t409 * mrSges(7,3);
t623 = t301 * mrSges(7,2);
t233 = t610 - t623;
t577 = t228 + t233;
t674 = -t415 / 0.2e1;
t675 = t414 / 0.2e1;
t737 = mrSges(6,3) + mrSges(7,2);
t627 = t299 * mrSges(7,2);
t232 = -mrSges(7,1) * t409 + t627;
t755 = t232 / 0.2e1;
t472 = -t231 * t674 - t415 * t755 - t577 * t675 - t737 * (t299 * t674 + t301 * t675);
t588 = t411 * t445;
t335 = -mrSges(5,2) * t409 + mrSges(5,3) * t588;
t579 = t446 * t335;
t587 = t411 * t446;
t337 = t409 * mrSges(5,1) + mrSges(5,3) * t587;
t581 = t445 * t337;
t398 = t409 * qJ(6);
t438 = -pkin(2) * t668 - pkin(1);
t330 = t409 * pkin(3) + t411 * pkin(8) + t438;
t562 = t666 * pkin(7);
t418 = -qJ(3) * t666 - t562;
t565 = t668 * pkin(7);
t420 = qJ(3) * t668 + t565;
t734 = t601 * t418 + t602 * t420;
t213 = t446 * t330 - t445 * t734;
t164 = pkin(9) * t587 + t213;
t130 = t409 * pkin(4) + t164;
t542 = t665 * t130;
t591 = t734 * t446;
t165 = t591 + (pkin(9) * t411 + t330) * t445;
t545 = t667 * t165;
t88 = t545 + t542;
t65 = t398 + t88;
t541 = t665 * t165;
t87 = t130 * t667 - t541;
t66 = -t409 * pkin(5) - t87;
t716 = m(7) / 0.2e1;
t718 = m(6) / 0.2e1;
t96 = t164 * t665 + t545;
t97 = t164 * t667 - t541;
t778 = -((t87 - t97) * t415 + (-t88 + t96) * t414) * t718 - ((-t66 - t97) * t415 + (-t65 + t96) * t414) * t716 + t581 / 0.2e1 - t579 / 0.2e1 - t472;
t777 = Ifges(7,4) + Ifges(6,5);
t776 = Ifges(7,2) + Ifges(6,3);
t543 = t602 * pkin(2);
t433 = -t543 - pkin(3);
t660 = t446 * pkin(4);
t417 = t433 - t660;
t497 = t414 * pkin(5) + t415 * qJ(6);
t312 = t417 + t497;
t344 = mrSges(7,1) * t414 + mrSges(7,3) * t415;
t747 = m(7) * t312 + t344;
t563 = t666 * pkin(2);
t331 = -t411 * pkin(3) + t409 * pkin(8) + t563;
t361 = -t602 * t418 + t420 * t601;
t215 = t446 * t331 + t361 * t445;
t216 = t445 * t331 - t361 * t446;
t300 = t415 * t409;
t227 = mrSges(6,2) * t411 - t300 * mrSges(6,3);
t302 = t409 * t414;
t229 = -mrSges(6,1) * t411 - t302 * mrSges(6,3);
t234 = -t300 * mrSges(7,2) - mrSges(7,3) * t411;
t561 = t665 * pkin(4);
t432 = t561 + qJ(6);
t564 = t667 * pkin(4);
t437 = -t564 - pkin(5);
t510 = t561 / 0.2e1;
t512 = t564 / 0.2e1;
t589 = t409 * t446;
t529 = -t589 / 0.2e1;
t590 = t409 * t445;
t530 = t590 / 0.2e1;
t714 = m(6) * pkin(4);
t566 = t714 / 0.2e1;
t678 = -t411 / 0.2e1;
t131 = -t411 * pkin(4) + pkin(9) * t589 + t215;
t166 = pkin(9) * t590 + t216;
t90 = t665 * t131 + t667 * t166;
t68 = -qJ(6) * t411 + t90;
t89 = t131 * t667 - t166 * t665;
t70 = t411 * pkin(5) - t89;
t609 = t411 * mrSges(7,1);
t620 = t302 * mrSges(7,2);
t230 = t609 + t620;
t701 = t230 / 0.2e1;
t697 = t300 / 0.2e1;
t698 = -t300 / 0.2e1;
t770 = -mrSges(7,3) / 0.2e1;
t771 = mrSges(6,2) / 0.2e1;
t772 = mrSges(7,1) / 0.2e1;
t726 = t90 * t771 - t89 * mrSges(6,1) / 0.2e1 + t70 * t772 + t68 * t770;
t769 = t302 / 0.2e1;
t723 = Ifges(6,6) * t698 + Ifges(7,6) * t697 + t678 * t776 + t769 * t777 - t726;
t775 = (t432 * t68 + t437 * t70) * t716 + Ifges(5,3) * t678 + t215 * mrSges(5,1) / 0.2e1 - t216 * mrSges(5,2) / 0.2e1 + t432 * t234 / 0.2e1 + t437 * t701 + (t665 * t90 + t667 * t89) * t566 + Ifges(5,5) * t529 + Ifges(5,6) * t530 + t229 * t512 + t227 * t510 + t723;
t766 = Ifges(6,6) - Ifges(7,6);
t506 = -t414 * t777 + t415 * t766;
t539 = t601 * pkin(2);
t504 = t539 + pkin(8);
t480 = t445 * (-pkin(9) - t504);
t425 = t446 * t504;
t569 = t446 * pkin(9) + t425;
t728 = t665 * t480 + t569 * t667;
t749 = t728 * mrSges(7,1);
t750 = t728 * mrSges(6,1);
t327 = -t667 * t480 + t569 * t665;
t764 = t327 * mrSges(7,3);
t765 = t327 * mrSges(6,2);
t774 = t506 - t749 - t750 + t765 - t764;
t773 = t749 / 0.2e1 + t750 / 0.2e1 + t764 / 0.2e1 - t765 / 0.2e1;
t692 = -t327 / 0.2e1;
t768 = m(4) * t563;
t753 = m(6) * t417;
t730 = t445 ^ 2 + t446 ^ 2;
t763 = t409 * t730;
t762 = (t233 / 0.2e1 + t228 / 0.2e1) * t327;
t758 = -t411 * mrSges(4,1) - t409 * mrSges(4,2);
t757 = -pkin(5) * t728 - qJ(6) * t327;
t342 = -t415 * mrSges(7,1) + t414 * mrSges(7,3);
t343 = -t415 * mrSges(6,1) - t414 * mrSges(6,2);
t635 = Ifges(7,5) * t414;
t346 = -Ifges(7,3) * t415 - t635;
t405 = Ifges(7,5) * t415;
t347 = Ifges(7,3) * t414 - t405;
t408 = Ifges(6,4) * t414;
t348 = Ifges(6,2) * t415 - t408;
t640 = Ifges(6,4) * t415;
t349 = -Ifges(6,2) * t414 - t640;
t350 = -Ifges(7,1) * t414 - t405;
t351 = -Ifges(7,1) * t415 + t635;
t352 = -Ifges(6,1) * t414 + t640;
t353 = -Ifges(6,1) * t415 - t408;
t756 = (t346 / 0.2e1 - t348 / 0.2e1 - t351 / 0.2e1 - t353 / 0.2e1) * t414 + (-t347 / 0.2e1 + t349 / 0.2e1 - t350 / 0.2e1 - t352 / 0.2e1) * t415 + t312 * t342 + t417 * t343;
t689 = t344 / 0.2e1;
t754 = t728 / 0.2e1;
t751 = t66 + t87;
t748 = t234 + t227;
t440 = m(7) * qJ(6) + mrSges(7,3);
t744 = qJD(5) * t440;
t743 = t440 * qJD(6);
t279 = -t623 / 0.2e1;
t740 = 0.2e1 * t279;
t739 = -t504 / 0.2e1;
t441 = Ifges(5,5) * t446;
t631 = Ifges(5,6) * t445;
t736 = Ifges(4,4) - t441 / 0.2e1 + t631 / 0.2e1;
t576 = -t231 + t232;
t442 = Ifges(5,4) * t446;
t423 = Ifges(5,1) * t445 + t442;
t642 = Ifges(5,4) * t445;
t421 = Ifges(5,2) * t446 + t642;
t324 = t411 * t421;
t735 = t337 * t739 + t324 / 0.4e1;
t733 = -t299 * t754 + t301 * t692;
t507 = -t299 * t766 - t301 * t777;
t500 = Ifges(5,2) * t445 - t442;
t174 = mrSges(6,1) * t301 + mrSges(6,2) * t299;
t325 = t411 * t423;
t729 = t335 * t739 + pkin(4) * t174 / 0.2e1 + t325 / 0.4e1;
t481 = -t347 / 0.4e1 + t349 / 0.4e1 - t350 / 0.4e1 - t352 / 0.4e1;
t483 = t346 / 0.4e1 - t348 / 0.4e1 - t351 / 0.4e1 - t353 / 0.4e1;
t551 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t552 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t725 = -t551 * t300 + t552 * t302;
t719 = m(5) / 0.2e1;
t717 = -m(7) / 0.2e1;
t715 = m(4) * pkin(2);
t713 = m(7) * pkin(4);
t712 = -mrSges(7,2) / 0.2e1;
t711 = t65 / 0.2e1;
t710 = -t66 / 0.2e1;
t709 = t87 / 0.2e1;
t708 = t88 / 0.2e1;
t707 = -t96 / 0.2e1;
t706 = -t97 / 0.2e1;
t705 = m(7) * t96;
t641 = Ifges(6,4) * t299;
t138 = -Ifges(6,2) * t301 + Ifges(6,6) * t409 + t641;
t704 = t138 / 0.4e1;
t173 = mrSges(7,1) * t301 - mrSges(7,3) * t299;
t703 = t173 / 0.2e1;
t636 = Ifges(7,5) * t301;
t175 = Ifges(7,3) * t299 - t636;
t702 = t175 / 0.4e1;
t700 = -t231 / 0.2e1;
t699 = -t234 / 0.2e1;
t679 = t409 / 0.4e1;
t677 = -t414 / 0.2e1;
t672 = t421 / 0.4e1;
t424 = Ifges(5,1) * t446 - t642;
t671 = -t424 / 0.4e1;
t670 = t445 / 0.2e1;
t669 = t446 / 0.2e1;
t664 = m(7) * t728;
t340 = pkin(5) * t415 - qJ(6) * t414;
t663 = m(7) * t340;
t662 = m(7) * t415;
t661 = pkin(4) * t445;
t657 = t87 * mrSges(6,2);
t656 = t87 * mrSges(7,3);
t655 = t88 * mrSges(6,1);
t654 = t88 * mrSges(7,1);
t651 = t96 * mrSges(6,1);
t650 = t96 * mrSges(7,1);
t649 = t97 * mrSges(6,2);
t648 = t97 * mrSges(7,3);
t646 = -t65 + t88;
t645 = mrSges(4,3) * t409;
t644 = mrSges(4,3) * t411;
t638 = Ifges(5,5) * t409;
t632 = Ifges(5,6) * t409;
t260 = -pkin(4) * t590 + t734;
t498 = -pkin(5) * t300 + t302 * qJ(6);
t110 = t260 - t498;
t261 = -pkin(4) * t588 + t361;
t499 = t301 * pkin(5) - t299 * qJ(6);
t111 = t261 + t499;
t135 = Ifges(7,5) * t302 - Ifges(7,6) * t411 + Ifges(7,3) * t300;
t137 = Ifges(6,4) * t302 - Ifges(6,2) * t300 - Ifges(6,6) * t411;
t139 = Ifges(7,1) * t302 - Ifges(7,4) * t411 + Ifges(7,5) * t300;
t141 = Ifges(6,1) * t302 - Ifges(6,4) * t300 - Ifges(6,5) * t411;
t619 = t302 * mrSges(7,3);
t624 = t300 * mrSges(7,1);
t171 = -t619 + t624;
t621 = t302 * mrSges(6,2);
t625 = t300 * mrSges(6,1);
t172 = t621 + t625;
t214 = t445 * t330 + t591;
t246 = -Ifges(5,6) * t411 + t409 * t500;
t248 = -Ifges(5,5) * t411 - t409 * t424;
t419 = t445 * mrSges(5,1) + t446 * mrSges(5,2);
t322 = t419 * t409;
t323 = t419 * t411;
t334 = mrSges(5,2) * t411 + mrSges(5,3) * t590;
t336 = -t411 * mrSges(5,1) + mrSges(5,3) * t589;
t140 = Ifges(7,1) * t299 + t409 * Ifges(7,4) + t636;
t295 = Ifges(6,4) * t301;
t142 = Ifges(6,1) * t299 + t409 * Ifges(6,5) - t295;
t524 = t140 / 0.2e1 + t142 / 0.2e1;
t292 = Ifges(7,5) * t299;
t136 = Ifges(7,6) * t409 + Ifges(7,3) * t301 + t292;
t525 = t136 / 0.2e1 - t138 / 0.2e1;
t249 = -t411 * t424 + t638;
t580 = t446 * t249;
t247 = t411 * t500 + t632;
t582 = t445 * t247;
t3 = (-Ifges(3,2) + Ifges(3,1)) * t668 * t666 + (-t137 / 0.2e1 + t135 / 0.2e1) * t301 - (-t139 / 0.2e1 - t141 / 0.2e1) * t299 + m(6) * (t260 * t261 + t87 * t89 + t88 * t90) + m(7) * (t110 * t111 + t65 * t68 + t66 * t70) + m(5) * (t213 * t215 + t214 * t216 + t361 * t734) - pkin(1) * (mrSges(3,1) * t666 + mrSges(3,2) * t668) - t734 * t323 + t524 * t302 + t525 * t300 + (-t446 * t248 / 0.2e1 + t246 * t670 - mrSges(4,2) * t563 - t736 * t411 + t551 * t301 - t552 * t299 + (-Ifges(5,3) - Ifges(4,2) + Ifges(4,1) - t776) * t409) * t411 + (-t580 / 0.2e1 + t582 / 0.2e1 + mrSges(4,1) * t563 + t736 * t409 + t725) * t409 - t361 * t322 + t214 * t334 + t216 * t335 + t213 * t336 + t215 * t337 + t261 * t172 + t260 * t174 + t68 * t233 + t65 * t234 + t88 * t227 + t90 * t228 + t87 * t229 + t66 * t230 + t89 * t231 + t70 * t232 + t110 * t173 + t111 * t171 + (t758 + t768) * t438 + (-t666 ^ 2 + t668 ^ 2) * Ifges(3,4);
t626 = t3 * qJD(1);
t608 = t414 * mrSges(7,2);
t607 = t414 * mrSges(6,3);
t606 = t415 * mrSges(7,2);
t605 = t415 * mrSges(6,3);
t168 = t299 * pkin(5) + t301 * qJ(6);
t132 = -pkin(4) * t587 + t168;
t501 = -mrSges(5,1) * t446 + t445 * mrSges(5,2);
t321 = t411 * t501;
t169 = mrSges(7,1) * t299 + mrSges(7,3) * t301;
t170 = mrSges(6,1) * t299 - mrSges(6,2) * t301;
t176 = -Ifges(6,2) * t299 - t295;
t177 = -Ifges(7,1) * t301 + t292;
t178 = -Ifges(6,1) * t301 - t641;
t453 = -t65 * t627 + t87 * t622 + t111 * t169 + t261 * t170 + (-t176 / 0.2e1 - t66 * mrSges(7,2) + t175 / 0.2e1 - t524) * t301 - (-t178 / 0.2e1 + t88 * mrSges(6,3) - t177 / 0.2e1 - t525) * t299 + t507 * t409 / 0.2e1;
t6 = t453 + t361 * t321 + t213 * t335 - t214 * t337 + t132 * t173 + t577 * t97 + t576 * t96 + ((t324 / 0.2e1 + t249 / 0.2e1 + t638 / 0.2e1 - t213 * mrSges(5,3)) * t445 + (t247 / 0.2e1 - t325 / 0.2e1 + t632 / 0.2e1 + t214 * mrSges(5,3) + (-m(6) * t261 - t174) * pkin(4)) * t446) * t411 + m(7) * (t111 * t132 + t65 * t97 + t66 * t96) + m(6) * (-t87 * t96 + t88 * t97);
t604 = t6 * qJD(1);
t7 = t453 + (m(7) * t111 + t173) * t168 + (m(7) * t66 + t576) * t88 + (m(7) * t65 + t577) * t87;
t603 = t7 * qJD(1);
t494 = t300 * t327 + t302 * t728;
t345 = mrSges(6,1) * t414 - mrSges(6,2) * t415;
t528 = t345 * t678;
t448 = (-t433 * t411 - t504 * t763) * t719 + (-t411 * t417 + t494) * t718 + (-t312 * t411 + t494) * t716 + t344 * t678 + t528 - t321 / 0.2e1 + (-t409 * t601 + t411 * t602) * t715 / 0.2e1 - mrSges(5,3) * t763 / 0.2e1 + t737 * (t300 * t674 + t302 * t677);
t454 = (t215 * t446 + t445 * t216) * t719 + (-t414 * t89 - t415 * t90) * t718 + (t414 * t70 - t415 * t68) * t716 + t334 * t670 + t336 * t669 + t768 / 0.2e1;
t14 = t229 * t677 + t230 * t675 + t674 * t748 - t448 + t454 + t758;
t600 = qJD(1) * t14;
t32 = m(7) * (-t111 * t299 + t409 * t65) - t299 * t173 + t409 * t233;
t599 = qJD(1) * t32;
t459 = (t414 * t646 - t415 * t751) * t716 + t472;
t487 = m(7) * t498;
t12 = (t770 + t771) * t302 + (t772 + mrSges(6,1) / 0.2e1) * t300 - t487 / 0.2e1 + t459;
t598 = t12 * qJD(1);
t592 = t361 * t411;
t15 = t577 * t302 + t576 * t300 + (-t579 + t581 + t645) * t409 + (-t173 - t174 + t323 + t644) * t411 + m(7) * (-t111 * t411 + t300 * t66 + t302 * t65) + m(6) * (-t261 * t411 - t300 * t87 + t302 * t88) + m(5) * (-t592 + (t213 * t445 - t214 * t446) * t409) + m(4) * (-t409 * t734 - t592);
t597 = t15 * qJD(1);
t596 = t261 * t445;
t578 = qJD(6) * t662;
t526 = t409 * t674;
t191 = (t526 + t698) * m(7);
t567 = t191 * qJD(1);
t560 = t432 * t627;
t559 = t663 / 0.2e1;
t553 = mrSges(7,2) / 0.2e1 + mrSges(6,3) / 0.2e1;
t550 = qJ(6) * t712;
t397 = -t608 / 0.2e1;
t548 = t606 / 0.2e1;
t523 = t755 + t700;
t521 = t327 * t96 + t728 * t97;
t519 = t441 - t631;
t514 = mrSges(6,3) * t561;
t513 = t414 * t564;
t511 = -t561 / 0.2e1;
t509 = 0.2e1 * t397;
t508 = -t624 / 0.2e1 - t625 / 0.2e1 + t619 / 0.2e1 - t621 / 0.2e1;
t503 = t564 * t622;
t502 = t299 * t514;
t20 = -t340 * t747 + t756;
t460 = t111 * t342 / 0.2e1 + t261 * t343 / 0.2e1 + t312 * t169 / 0.2e1 + t417 * t170 / 0.2e1 + t506 * t679;
t485 = -t178 / 0.4e1 - t177 / 0.4e1 - t136 / 0.4e1 + t704;
t486 = -t176 / 0.4e1 - t142 / 0.4e1 - t140 / 0.4e1 + t702;
t447 = -t762 + t486 * t414 + t485 * t415 + ((-t88 / 0.2e1 + t711) * t415 + (t710 - t87 / 0.2e1) * t414 + t733) * mrSges(7,2) + (mrSges(6,3) * t692 + t483) * t301 - (mrSges(6,3) * t754 + t481) * t299 + (-t111 * t340 + t168 * t312 + t327 * t646) * t716 + t168 * t689 - t340 * t703 + t460 + (t716 * t751 + t523) * t728;
t465 = (-pkin(5) * t70 + qJ(6) * t68) * t717 + pkin(5) * t701 + qJ(6) * t699;
t4 = t447 + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t411 + t465 - t725 + t726;
t496 = t4 * qJD(1) + t20 * qJD(2);
t455 = (t300 * t437 + t302 * t432) * t716 + (-t300 * t667 + t302 * t665) * t566 + mrSges(5,1) * t530 + mrSges(5,2) * t589 / 0.2e1 + t508;
t10 = t455 + t778;
t495 = t10 * qJD(1);
t492 = -t432 * t414 - t437 * t415;
t149 = t747 * t415;
t461 = (t111 * t415 - t299 * t312 + t409 * t728) * t717 + t299 * t689 + t173 * t674;
t490 = t70 * t716 + t609 / 0.2e1;
t30 = t461 + t490 + t620;
t491 = qJD(1) * t30 - qJD(2) * t149;
t489 = t445 * t504;
t479 = -t342 - t343;
t478 = mrSges(5,3) * t504 / 0.2e1;
t475 = -t414 * t665 + t415 * t667;
t333 = -t340 + t661;
t452 = (t111 * t333 + t132 * t312 - t327 * t65 + t66 * t728 + t521) * t716 + t132 * t689 + t333 * t703 + t361 * t419 / 0.2e1 + t433 * t321 / 0.2e1 + t460;
t469 = -t327 * t88 - t728 * t87 + t521;
t1 = t65 * t548 + t66 * t397 - (t136 + t178 + t177) * t415 / 0.4e1 - (t176 + t142 + t140) * t414 / 0.4e1 + t729 * t445 + t483 * t301 - t481 * t299 + t580 / 0.4e1 + t735 * t446 + ((-t417 * t587 + t596) * pkin(4) + t469) * t718 + t414 * t702 + t415 * t704 + t605 * t708 + t607 * t709 + t519 * t679 + (t423 - t500) * t588 / 0.4e1 + t232 * t754 + t587 * t671 + t587 * t672 + t528 * t660 + t728 * t700 + t577 * t692 - t775 - t582 / 0.4e1 + t452 + t730 * t411 * t478 + t737 * (t96 * t674 + t97 * t677 + t733);
t17 = t433 * t419 + (-t500 / 0.2e1 + t423 / 0.2e1) * t446 + (pkin(4) * t345 - t421 / 0.2e1 + t424 / 0.2e1) * t445 + t661 * t753 + t756 + t747 * t333;
t474 = t1 * qJD(1) + t17 * qJD(2);
t462 = t610 + (t432 * t409 + t65) * t716;
t34 = -t705 / 0.2e1 + t462;
t416 = m(7) * t432 + mrSges(7,3);
t471 = -qJD(1) * t34 - qJD(4) * t416;
t36 = t610 + 0.2e1 * (t545 / 0.4e1 + t542 / 0.4e1 + t398 / 0.2e1 - t88 / 0.4e1) * m(7);
t470 = qJD(1) * t36 + qJD(4) * t440 + t744;
t463 = m(7) * (-pkin(4) * t475 + t492);
t104 = t559 - t463 / 0.2e1;
t464 = (-mrSges(6,1) - mrSges(7,1)) * t561 + (-mrSges(6,2) + mrSges(7,3)) * t564;
t253 = -(t432 * t667 + t437 * t665) * t713 - t464;
t450 = ((t561 - t432) * t327 + (t437 + t564) * t728) * t716 + t432 * t548 + t437 * t397 + t513 * t712 + t511 * t606 - t773;
t456 = pkin(5) * t397 + t415 * t550 + t757 * t717 + t773;
t29 = t450 + t456;
t449 = (t432 * t87 + t437 * t88 + (t65 * t667 + t66 * t665) * pkin(4)) * t716 - t657 / 0.2e1 + t656 / 0.2e1 - t655 / 0.2e1 - t654 / 0.2e1 - t560 / 0.2e1 + t437 * t279 + t231 * t511 + t232 * t510 + t503 / 0.2e1 - t502 / 0.2e1 + t577 * t512;
t458 = (-pkin(5) * t96 + qJ(6) * t97) * t717 + t651 / 0.2e1 + t650 / 0.2e1 + t649 / 0.2e1 - t648 / 0.2e1 + pkin(5) * t279 - t299 * t550;
t9 = t449 + t458;
t466 = t9 * qJD(1) + t29 * qJD(2) - t104 * qJD(3) - t253 * qJD(4);
t395 = mrSges(7,3) + (0.2e1 * t510 + qJ(6)) * m(7);
t203 = t509 + t664;
t190 = (t526 + t697) * m(7);
t127 = t664 / 0.2e1 + m(7) * t754 + t509;
t71 = t463 / 0.2e1 + t559 + t479;
t35 = (0.2e1 * t398 + t88) * t716 + t610 + t740 + m(7) * t708;
t33 = t740 + t705 / 0.2e1 + t462;
t31 = t409 * t397 + t620 / 0.2e1 - t461 + t490;
t24 = t450 - t456 + t506;
t16 = t448 + (t699 - t227 / 0.2e1) * t415 + (t701 - t229 / 0.2e1) * t414 + t454;
t13 = t459 + t487 / 0.2e1 + t508;
t11 = t455 - t778;
t8 = t449 - t458 + t507;
t5 = t447 - t465 + t723;
t2 = (pkin(4) * t596 + t469) * t718 + t441 * t679 + ((t706 + t709) * mrSges(6,3) + (t710 + t706) * mrSges(7,2) + t486) * t414 + ((t707 + t708) * mrSges(6,3) + (t711 + t707) * mrSges(7,2) + t485) * t415 + (-t632 / 0.4e1 - t247 / 0.4e1 + t729) * t445 + (t249 / 0.4e1 + t735) * t446 + ((t423 / 0.4e1 - t500 / 0.4e1 + t445 * t478) * t445 + (t671 + t672 + t446 * t478 + (-t753 / 0.2e1 - t345 / 0.2e1) * pkin(4)) * t446) * t411 - (t553 * t728 + t481) * t299 + (-t327 * t553 + t483) * t301 + t523 * t728 - t762 + t452 + t775;
t18 = [qJD(2) * t3 + qJD(3) * t15 + qJD(4) * t6 + qJD(5) * t7 + qJD(6) * t32, t16 * qJD(3) + t2 * qJD(4) + t5 * qJD(5) + t31 * qJD(6) + t626 + (t421 * t530 + t423 * t529 + (t139 + t141) * t674 + (-t215 * t445 + t216 * t446) * mrSges(5,3) + (-m(6) * t89 + m(7) * t70 - t229 + t230) * t327 + (m(6) * t90 + m(7) * t68 + t748) * t728 - t336 * t489 + m(5) * (-t215 * t489 + t216 * t425) + t349 * t698 + t347 * t697 + t135 * t675 + t137 * t677 + t747 * t110 + t248 * t670 + t246 * t669 + t539 * t644 + t543 * t645 + Ifges(3,5) * t668 - Ifges(3,6) * t666 - t70 * t606 - t90 * t607 - t68 * t608 - mrSges(3,1) * t565 + t89 * t605 + mrSges(3,2) * t562 + (-t601 * t715 + mrSges(4,2)) * t361 + (m(5) * t433 - t602 * t715 - mrSges(4,1) + t501) * t734 + t334 * t425 + (t353 + t351) * t769 + (Ifges(5,5) * t445 + Ifges(5,6) * t446 - t414 * t766 - t415 * t777) * t678 - t433 * t322 + t417 * t172 - Ifges(4,5) * t409 + Ifges(4,6) * t411 + t312 * t171 + (t345 + t753) * t260) * qJD(2), t16 * qJD(2) + t11 * qJD(4) + t13 * qJD(5) + t190 * qJD(6) + t597 + 0.2e1 * (t716 + t718) * qJD(3) * (t300 * t414 - t415 * t302) t604 + t2 * qJD(2) + t11 * qJD(3) + (Ifges(5,5) * t588 + Ifges(5,6) * t587 + m(7) * (t432 * t97 + t437 * t96) - t649 - t650 - t651 - t437 * t623 + t648 - t560 + (t665 * t97 - t667 * t96) * t714 - t213 * mrSges(5,2) - t214 * mrSges(5,1) - t502 + t503 + t507) * qJD(4) + t8 * qJD(5) + t33 * qJD(6), t603 + t5 * qJD(2) + t13 * qJD(3) + t8 * qJD(4) + (m(7) * (-pkin(5) * t88 + qJ(6) * t87) - t657 - t654 - t655 + t656 + t499 * mrSges(7,2) + t507) * qJD(5) + t35 * qJD(6), qJD(2) * t31 + qJD(3) * t190 + qJD(4) * t33 + qJD(5) * t35 + t599; -qJD(3) * t14 + qJD(4) * t1 + qJD(5) * t4 - qJD(6) * t30 - t626, qJD(4) * t17 + qJD(5) * t20 + qJD(6) * t149, -t600 (mrSges(6,3) * t513 + t415 * t514 - t437 * t608 + t432 * t606 + mrSges(5,2) * t489 - mrSges(5,1) * t425 + m(7) * (-t327 * t432 + t437 * t728) + (-t327 * t665 - t667 * t728) * t714 + t519 + t774) * qJD(4) + t24 * qJD(5) + t127 * qJD(6) + t474, t24 * qJD(4) + (m(7) * t757 + t497 * mrSges(7,2) + t774) * qJD(5) + t203 * qJD(6) + t496, qJD(4) * t127 + qJD(5) * t203 - t491; qJD(2) * t14 - qJD(4) * t10 + qJD(5) * t12 + qJD(6) * t191 - t597, t600, 0 (m(7) * t492 + t475 * t714 - t419 + t479) * qJD(4) + t71 * qJD(5) - t495 - t578, t598 + t71 * qJD(4) + (t479 + t663) * qJD(5) - t578, t567 + 0.2e1 * (-qJD(4) / 0.2e1 - qJD(5) / 0.2e1) * t662; -qJD(2) * t1 + qJD(3) * t10 + qJD(5) * t9 + qJD(6) * t34 - t604, qJD(5) * t29 - t474, -qJD(5) * t104 + t495, -qJD(5) * t253 + qJD(6) * t416 ((-pkin(5) * t665 + qJ(6) * t667) * t713 + t464) * qJD(5) + t395 * qJD(6) + t466, qJD(5) * t395 - t471; -qJD(2) * t4 - qJD(3) * t12 - qJD(4) * t9 + qJD(6) * t36 - t603, -qJD(4) * t29 - t496, qJD(4) * t104 - t598, -t466 + t743, t743, t470; qJD(2) * t30 - qJD(3) * t191 - qJD(4) * t34 - qJD(5) * t36 - t599, t491, -t567, t471 - t744, -t470, 0;];
Cq  = t18;
