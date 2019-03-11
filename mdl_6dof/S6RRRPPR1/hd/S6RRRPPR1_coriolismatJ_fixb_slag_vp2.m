% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:21:23
% EndTime: 2019-03-09 15:21:45
% DurationCPUTime: 13.04s
% Computational Cost: add. (43958->587), mult. (85110->780), div. (0->0), fcn. (102700->10), ass. (0->348)
t398 = sin(qJ(3));
t399 = sin(qJ(2));
t400 = cos(qJ(3));
t401 = cos(qJ(2));
t370 = -t398 * t399 + t400 * t401;
t371 = -t398 * t401 - t400 * t399;
t396 = sin(pkin(10));
t563 = cos(pkin(10));
t317 = t563 * t370 + t371 * t396;
t395 = sin(pkin(11));
t397 = cos(pkin(11));
t607 = sin(qJ(6));
t608 = cos(qJ(6));
t433 = t607 * t395 - t608 * t397;
t233 = t433 * t317;
t582 = t233 * mrSges(7,2);
t369 = -t395 * t608 - t397 * t607;
t231 = t369 * t317;
t584 = t231 * mrSges(7,1);
t510 = t584 / 0.2e1 + t582 / 0.2e1;
t629 = m(7) / 0.2e1;
t631 = m(6) / 0.2e1;
t624 = -pkin(8) - pkin(7);
t377 = t624 * t399;
t378 = t624 * t401;
t338 = t377 * t398 - t378 * t400;
t439 = t370 * qJ(4) + t338;
t647 = t400 * t377 + t398 * t378;
t657 = t371 * qJ(4) + t647;
t665 = t396 * t657 + t563 * t439;
t538 = t317 * t395;
t671 = pkin(5) * t538 + t665;
t699 = t629 * t671 + t631 * t665 - t510;
t434 = t396 * t370 - t371 * t563;
t578 = t434 * mrSges(5,3);
t537 = t317 * t397;
t474 = t537 / 0.2e1;
t570 = t395 * mrSges(6,1);
t482 = t570 / 0.2e1;
t428 = t434 * t433;
t158 = -mrSges(7,1) * t317 + mrSges(7,3) * t428;
t232 = t369 * t434;
t156 = mrSges(7,2) * t317 + mrSges(7,3) * t232;
t621 = t156 / 0.2e1;
t663 = -t369 / 0.2e1;
t438 = t158 * t663 + t433 * t621;
t387 = -pkin(2) * t401 - pkin(1);
t341 = -t370 * pkin(3) + t387;
t203 = -pkin(4) * t317 - qJ(5) * t434 + t341;
t102 = t397 * t203 - t395 * t665;
t103 = t395 * t203 + t397 * t665;
t450 = -t102 * t395 + t103 * t397;
t536 = t434 * t395;
t247 = mrSges(6,2) * t317 - mrSges(6,3) * t536;
t516 = t397 * t247;
t535 = t434 * t397;
t249 = -mrSges(6,1) * t317 - mrSges(6,3) * t535;
t517 = t395 * t249;
t649 = t516 / 0.2e1 - t517 / 0.2e1;
t571 = t369 * mrSges(7,3);
t660 = t571 / 0.2e1;
t611 = -t433 / 0.2e1;
t661 = t232 * t611;
t667 = mrSges(7,3) * t661 + t428 * t660 + t450 * t631 - t438 + t649;
t698 = mrSges(6,2) * t474 + t317 * t482 - t667 + t699;
t569 = t397 * mrSges(6,2);
t697 = (t569 / 0.2e1 + t482) * t317 + t667 + t699;
t590 = Ifges(6,2) * t395;
t593 = Ifges(6,4) * t397;
t168 = t434 * Ifges(6,6) + (-t590 + t593) * t317;
t594 = Ifges(6,4) * t395;
t169 = t434 * Ifges(6,5) + (Ifges(6,1) * t397 - t594) * t317;
t591 = Ifges(7,4) * t369;
t324 = -Ifges(7,2) * t433 - t591;
t362 = Ifges(7,4) * t433;
t326 = -Ifges(7,1) * t369 - t362;
t609 = t397 / 0.2e1;
t615 = t434 / 0.2e1;
t619 = -t233 / 0.2e1;
t620 = t231 / 0.2e1;
t662 = t395 / 0.2e1;
t89 = -Ifges(7,4) * t233 + Ifges(7,2) * t231 + Ifges(7,6) * t434;
t91 = -Ifges(7,1) * t233 + Ifges(7,4) * t231 + Ifges(7,5) * t434;
t422 = t89 * t611 + t169 * t662 + t168 * t609 + t324 * t620 + t326 * t619 - (Ifges(6,2) * t397 + t594) * t538 / 0.2e1 + (Ifges(6,1) * t395 + t593) * t474 - Ifges(5,6) * t434 + Ifges(5,5) * t317 + Ifges(4,6) * t371 + Ifges(4,5) * t370 + t91 * t663 + (Ifges(6,5) * t395 - Ifges(7,5) * t369 + Ifges(6,6) * t397 - Ifges(7,6) * t433) * t615;
t658 = t647 * mrSges(4,2);
t659 = t338 * mrSges(4,1);
t375 = -mrSges(6,1) * t397 + mrSges(6,2) * t395;
t675 = t665 * t375;
t683 = t665 * mrSges(5,1);
t664 = -t396 * t439 + t563 * t657;
t684 = t664 * mrSges(5,2);
t322 = mrSges(7,1) * t433 - mrSges(7,2) * t369;
t690 = t671 * t322;
t696 = t422 - t658 - t659 + t675 - t683 - t684 + t690;
t695 = -t658 / 0.2e1 - t659 / 0.2e1 - t684 / 0.2e1 + t690 / 0.2e1;
t694 = m(5) * t341;
t150 = pkin(5) * t536 - t664;
t693 = t150 * t671;
t604 = pkin(2) * t398;
t379 = t396 * t604;
t603 = pkin(2) * t400;
t386 = pkin(3) + t603;
t349 = t386 * t563 - t379;
t348 = -pkin(4) - t349;
t600 = t397 * pkin(5);
t340 = t348 - t600;
t692 = t340 * t671;
t478 = t563 * pkin(3);
t383 = -t478 - pkin(4);
t374 = t383 - t600;
t691 = t374 * t671;
t688 = t675 / 0.2e1 - t683 / 0.2e1;
t119 = -t582 - t584;
t120 = -mrSges(7,1) * t232 - mrSges(7,2) * t428;
t155 = -mrSges(7,2) * t434 + mrSges(7,3) * t231;
t157 = mrSges(7,1) * t434 + mrSges(7,3) * t233;
t455 = t569 + t570;
t242 = t455 * t317;
t243 = t455 * t434;
t246 = -mrSges(6,2) * t434 - mrSges(6,3) * t538;
t248 = mrSges(6,1) * t434 - mrSges(6,3) * t537;
t311 = t317 * mrSges(5,2);
t443 = Ifges(7,5) * t619 + Ifges(7,6) * t620;
t454 = Ifges(6,5) * t397 - Ifges(6,6) * t395;
t75 = -pkin(5) * t317 - pkin(9) * t535 + t102;
t80 = -pkin(9) * t536 + t103;
t46 = -t607 * t80 + t608 * t75;
t47 = t607 * t75 + t608 * t80;
t577 = t317 * mrSges(5,3);
t579 = t434 * mrSges(5,1);
t394 = t397 ^ 2;
t652 = t394 / 0.2e1;
t592 = Ifges(7,4) * t428;
t90 = Ifges(7,2) * t232 - Ifges(7,6) * t317 - t592;
t229 = Ifges(7,4) * t232;
t92 = -Ifges(7,1) * t428 - Ifges(7,5) * t317 + t229;
t687 = t671 * t120 + t102 * t248 + t103 * t246 + t150 * t119 + t341 * (t311 + t579) + t387 * (-mrSges(4,1) * t371 + mrSges(4,2) * t370) + t46 * t157 + t47 * t155 + t90 * t620 + t232 * t89 / 0.2e1 + t92 * t619 - t428 * t91 / 0.2e1 + (-Ifges(7,5) * t428 + Ifges(7,6) * t232) * t615 - (t577 + t242) * t664 + ((Ifges(5,1) + Ifges(6,1) * t652 + (-t593 + t590 / 0.2e1) * t395) * t317 - Ifges(5,4) * t434 + t454 * t615 - t395 * t168 / 0.2e1 + t169 * t609) * t434 - (-mrSges(5,3) * t664 + (Ifges(6,3) + Ifges(5,2) + Ifges(7,3)) * t434 + (-Ifges(5,4) + t454) * t317 + t443) * t317 + t665 * t243;
t393 = t395 ^ 2;
t504 = t393 + t394;
t646 = mrSges(6,3) * t504;
t682 = -mrSges(5,2) + t646;
t681 = t348 * t665;
t680 = t383 * t665;
t679 = t395 * t664;
t678 = t397 * t664;
t464 = t563 * t398;
t350 = pkin(2) * t464 + t396 * t386;
t676 = t664 * t350;
t674 = t665 * t664;
t602 = pkin(3) * t371;
t215 = pkin(4) * t434 - qJ(5) * t317 - t602;
t110 = t397 * t215 - t679;
t458 = pkin(5) * t434 - pkin(9) * t537;
t77 = t110 + t458;
t111 = t395 * t215 + t678;
t500 = pkin(9) * t538;
t82 = -t500 + t111;
t52 = -t607 * t82 + t608 * t77;
t53 = t607 * t77 + t608 * t82;
t626 = -mrSges(7,2) / 0.2e1;
t627 = mrSges(7,1) / 0.2e1;
t639 = Ifges(7,3) * t615 + t443;
t673 = -t52 * t627 - t53 * t626 - t639;
t392 = t399 * pkin(2);
t204 = t215 + t392;
t104 = t397 * t204 - t679;
t76 = t104 + t458;
t105 = t395 * t204 + t678;
t81 = -t500 + t105;
t48 = -t607 * t81 + t608 * t76;
t49 = t607 * t76 + t608 * t81;
t672 = -t48 * t627 - t49 * t626 - t639;
t670 = t396 * t664 - t563 * t665;
t625 = -mrSges(7,3) / 0.2e1;
t656 = t369 ^ 2 + t433 ^ 2;
t357 = t563 * t603 - t379;
t293 = t369 * t357;
t294 = t433 * t357;
t654 = (t293 * t369 + t294 * t433) * mrSges(7,3) + (-mrSges(4,1) * t398 - mrSges(4,2) * t400) * pkin(2);
t653 = t629 + t631;
t601 = pkin(3) * t396;
t380 = qJ(5) + t601;
t462 = t504 * t380;
t651 = m(6) * t462;
t650 = t120 + t243;
t347 = qJ(5) + t350;
t332 = (-pkin(9) - t347) * t395;
t390 = t397 * pkin(9);
t514 = t397 * t347;
t333 = t390 + t514;
t257 = t332 * t608 - t333 * t607;
t258 = t332 * t607 + t333 * t608;
t507 = t257 * t369 - t258 * t433;
t645 = -t347 * t248 / 0.2e1 - t110 * mrSges(6,3) / 0.2e1;
t644 = t347 * t246 / 0.2e1 + t111 * mrSges(6,3) / 0.2e1;
t448 = -t110 * t395 + t111 * t397;
t640 = -mrSges(5,1) + t322 + t375;
t638 = (-t233 * t611 + t369 * t620) * mrSges(7,3);
t321 = -t369 * mrSges(7,1) - mrSges(7,2) * t433;
t583 = t232 * mrSges(7,2);
t415 = t428 * t627 - t583 / 0.2e1;
t581 = t428 * mrSges(7,1);
t442 = t583 / 0.2e1 - t581 / 0.2e1;
t64 = -t415 + t442;
t637 = qJD(1) * t64 + (qJD(2) + qJD(3)) * t321;
t323 = Ifges(7,2) * t369 - t362;
t325 = -Ifges(7,1) * t433 + t591;
t459 = t369 * t324 / 0.2e1 + t325 * t663 + (t323 + t326) * t611;
t460 = t155 * t663 + t157 * t611 + t246 * t662 + t248 * t609;
t356 = (t396 * t400 + t464) * pkin(2);
t530 = t350 * t434;
t531 = t349 * t317;
t532 = t348 * t242;
t534 = t340 * t119;
t543 = t258 * t155;
t544 = t257 * t157;
t547 = t664 * t356;
t566 = t53 * t433;
t567 = t52 * t369;
t612 = -t357 / 0.2e1;
t613 = t356 / 0.2e1;
t617 = t293 / 0.2e1;
t633 = m(5) / 0.2e1;
t635 = (t676 - t547 + (-t349 + t357) * t665) * t633 + (t347 * t448 + t357 * t450 - t547 + t681) * t631 + (t150 * t356 + t257 * t52 + t258 * t53 + t293 * t46 - t294 * t47 + t692) * t629 + t544 / 0.2e1 + t543 / 0.2e1 + t158 * t617 - t294 * t621 + t534 / 0.2e1 + t532 / 0.2e1 + (t567 / 0.2e1 - t566 / 0.2e1) * mrSges(7,3) + (-t317 * t612 + t434 * t613 - t530 / 0.2e1 - t531 / 0.2e1) * mrSges(5,3) + t695;
t632 = -m(6) / 0.2e1;
t630 = -m(7) / 0.2e1;
t628 = m(5) * pkin(3);
t118 = -t581 + t583;
t622 = t118 / 0.2e1;
t618 = -t258 / 0.2e1;
t354 = (-pkin(9) - t380) * t395;
t513 = t397 * t380;
t355 = t390 + t513;
t306 = t354 * t607 + t355 * t608;
t616 = -t306 / 0.2e1;
t614 = t322 / 0.2e1;
t610 = t375 / 0.2e1;
t606 = m(6) * t434;
t605 = m(7) * t434;
t599 = m(7) * qJD(4);
t596 = mrSges(6,3) * t317;
t595 = Ifges(4,4) * t371;
t251 = -mrSges(5,1) * t317 + mrSges(5,2) * t434;
t346 = t392 - t602;
t456 = Ifges(4,4) * t370 + (-Ifges(4,1) + Ifges(4,2)) * t371;
t3 = t105 * t247 + t104 * t249 + t48 * t158 + t49 * t156 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t399) * t399 + (-mrSges(4,1) * t392 + t456) * t370 + (-mrSges(4,2) * t392 - t595) * t371 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t401 + (Ifges(3,1) - Ifges(3,2)) * t399) * t401 + m(4) * t387 * t392 + m(6) * (t102 * t104 + t103 * t105 - t674) + m(7) * (t46 * t48 + t47 * t49 + t693) + t687 + (t251 + t694) * t346;
t580 = t3 * qJD(1);
t572 = t433 * mrSges(7,3);
t4 = t111 * t247 + t110 * t249 + t52 * t158 + t53 * t156 + t456 * t370 + (-pkin(3) * t251 - t595) * t371 - t602 * t694 + m(6) * (t102 * t110 + t103 * t111 - t674) + m(7) * (t46 * t52 + t47 * t53 + t693) + t687;
t568 = t4 * qJD(1);
t121 = Ifges(7,2) * t428 + t229;
t122 = Ifges(7,1) * t232 + t592;
t508 = Ifges(7,5) * t232 + Ifges(7,6) * t428;
t9 = -t317 * t508 / 0.2e1 + t150 * t118 + t46 * t156 - t47 * t158 - (t122 / 0.2e1 - t90 / 0.2e1 - t47 * mrSges(7,3)) * t428 + (t92 / 0.2e1 + t121 / 0.2e1 - t46 * mrSges(7,3)) * t232;
t565 = t9 * qJD(1);
t463 = t504 * t347;
t409 = t638 + (t652 + t393 / 0.2e1) * t596 + (t317 * t350 - t349 * t434) * t633 + (t317 * t463 + t348 * t434) * t631 + (t231 * t257 - t233 * t258 + t340 * t434) * t629;
t419 = t346 * t633 + (t104 * t397 + t105 * t395) * t631 + (-t369 * t49 - t433 * t48) * t629;
t423 = -t311 - t460;
t467 = t614 + t610;
t11 = (-mrSges(5,1) + t467) * t434 + t409 - t419 + t423;
t562 = qJD(1) * t11;
t23 = t156 * t232 + t158 * t428 + (t369 * t47 + t433 * t46) * t605 + (-t102 * t397 - t103 * t395) * t606 - t247 * t536 - t249 * t535;
t561 = qJD(1) * t23;
t548 = t664 * t434;
t10 = -t233 * t156 + t231 * t158 + (t516 - t517 + t577) * t317 + (t650 + t578) * t434 + m(7) * (t150 * t434 + t231 * t46 - t233 * t47) + m(6) * (t317 * t450 - t548) + m(5) * (t317 * t665 - t548);
t559 = t10 * qJD(1);
t558 = t104 * t395;
t557 = t105 * t397;
t545 = t428 * t369;
t305 = t354 * t608 - t355 * t607;
t540 = t305 * t157;
t539 = t306 * t155;
t533 = t340 * t321;
t521 = t374 * t119;
t520 = t383 * t242;
t518 = t395 * t248;
t411 = -t504 * t606 / 0.2e1 - t656 * t605 / 0.2e1;
t62 = -t411 + (m(6) + m(7)) * t615;
t512 = t62 * qJD(1);
t505 = -Ifges(7,5) * t433 + Ifges(7,6) * t369;
t319 = t321 * qJD(6);
t502 = -t628 / 0.2e1;
t499 = mrSges(6,3) * t558;
t498 = mrSges(6,3) * t557;
t161 = -t293 * t433 + t294 * t369;
t493 = t161 * t629;
t492 = qJD(3) * t493;
t491 = t380 * t518;
t490 = t246 * t513;
t484 = t572 / 0.2e1;
t483 = -t571 / 0.2e1;
t466 = t323 / 0.4e1 + t326 / 0.4e1;
t465 = -t324 / 0.4e1 + t325 / 0.4e1;
t461 = t578 * t601;
t457 = t356 * t653;
t453 = t369 * t46 - t433 * t47;
t452 = t478 * t577;
t413 = m(7) * (t232 * t258 + t257 * t428 + t453);
t13 = -t413 / 0.2e1 + t698;
t446 = mrSges(7,3) * t656 + t646;
t74 = m(6) * t463 + m(7) * t507 + t446;
t451 = -qJD(1) * t13 + qJD(2) * t74;
t449 = t557 - t558;
t447 = t305 * t369 - t306 * t433;
t441 = mrSges(7,1) * t617 - t294 * t626;
t19 = (t545 / 0.2e1 + t661) * mrSges(7,3) + t438 + t510;
t436 = t19 * qJD(1);
t405 = (t317 * t462 + t383 * t434) * t631 + (t231 * t305 - t233 * t306 + t374 * t434) * t629 + t434 * t614 + t434 * t610 + (t317 * t396 - t434 * t563) * t628 / 0.2e1 + t504 * t596 / 0.2e1 + t638;
t417 = (t110 * t397 + t111 * t395) * t631 + (-t369 * t53 - t433 * t52) * t629 + t371 * t502;
t15 = t405 - t417 + t423 - t579;
t435 = -t15 * qJD(1) + qJD(2) * t493;
t295 = t374 * t321;
t421 = t295 / 0.2e1 + t533 / 0.2e1 + t459 + (t660 + t483) * t306;
t29 = t421 - t441;
t45 = -t295 - t459;
t418 = -(t121 / 0.4e1 + t92 / 0.4e1) * t433 + (t90 / 0.4e1 - t122 / 0.4e1) * t369 + t150 * t321 / 0.2e1 - t317 * t505 / 0.4e1;
t407 = -(mrSges(7,3) * t616 + t465) * t428 + (t305 * t625 + t466) * t232 + t305 * t621 + t158 * t616 + t374 * t622 + t418;
t8 = t407 + t673;
t432 = t8 * qJD(1) + t29 * qJD(2) - t45 * qJD(3);
t412 = m(7) * (t232 * t306 + t305 * t428 + t453);
t16 = -t412 / 0.2e1 + t698;
t414 = (t462 + t463) * t632 + (t447 + t507) * t630 - t446;
t55 = t457 + t414;
t99 = -m(7) * t447 - t446 - t651;
t430 = -qJD(1) * t16 - qJD(2) * t55 - qJD(3) * t99;
t35 = t459 + t533;
t408 = -(mrSges(7,3) * t618 + t465) * t428 + (t257 * t625 + t466) * t232 + t257 * t621 + t158 * t618 + t340 * t622 + t418;
t6 = t408 + t672;
t429 = -t6 * qJD(1) - t35 * qJD(2);
t402 = (t380 * t449 + t680) * t632 + (t305 * t48 + t306 * t49 + t691) * t630 - t540 / 0.2e1 - t539 / 0.2e1 - t521 / 0.2e1 - t520 / 0.2e1 + t670 * t502 + t499 / 0.2e1 - t498 / 0.2e1 + t491 / 0.2e1 - t490 / 0.2e1 + t48 * t483 + t49 * t484 + t461 / 0.2e1 + t452 / 0.2e1 - t688 - t695;
t2 = t649 * t357 + t645 * t395 + t644 * t397 + t650 * t613 + t402 + t635 + t688;
t34 = t682 * t357 + t640 * t356 + m(7) * (t257 * t293 - t258 * t294 + t340 * t356) + m(5) * (-t349 * t356 + t350 * t357) + m(6) * (t348 * t356 + t357 * t463) + t654;
t427 = -t2 * qJD(1) - t34 * qJD(2) - t161 * t599 / 0.2e1;
t320 = t321 * qJD(5);
t65 = t415 + t442;
t63 = t434 * t653 + t411;
t56 = t457 - t414;
t30 = t421 + t441;
t20 = t232 * t484 + t545 * t625 - t438 + t510;
t18 = t405 + t417 + t460;
t17 = t412 / 0.2e1 + t697;
t14 = t413 / 0.2e1 + t697;
t12 = t434 * t467 + t409 + t419 + t460;
t7 = t407 - t673;
t5 = t408 - t672;
t1 = -t402 + t422 + (t243 / 0.2e1 + t120 / 0.2e1) * t356 + (t610 - mrSges(5,1) / 0.2e1) * t665 + (t249 * t612 + t645) * t395 + (t357 * t247 / 0.2e1 + t644) * t397 + t635;
t21 = [qJD(2) * t3 + qJD(3) * t4 + qJD(4) * t10 + qJD(5) * t23 + qJD(6) * t9, t580 + (t534 - Ifges(3,6) * t399 + Ifges(3,5) * t401 - t49 * t572 + t48 * t571 + (-mrSges(3,1) * t401 + mrSges(3,2) * t399) * pkin(7) - t499 + t246 * t514 - t347 * t518 + (-t531 - t530) * mrSges(5,3) + (-t370 * t603 + t371 * t604) * mrSges(4,3) + t543 + t544 + t498 + t532 + t696) * qJD(2) + t1 * qJD(3) + t12 * qJD(4) + t14 * qJD(5) + t5 * qJD(6) + 0.2e1 * ((t347 * t449 + t681) * t631 + (-t349 * t665 + t676) * t633 + (t257 * t48 + t258 * t49 + t692) * t629 + m(4) * (-t338 * t400 + t398 * t647) * pkin(2) / 0.2e1) * qJD(2), t568 + t1 * qJD(2) + (t670 * t628 + t520 + t521 + t539 + t540 + t490 - t491 + m(7) * (t305 * t52 + t306 * t53 + t691) + m(6) * (t380 * t448 + t680) - t461 - t452 + (t567 - t566) * mrSges(7,3) + t448 * mrSges(6,3) + t696) * qJD(3) + t18 * qJD(4) + t17 * qJD(5) + t7 * qJD(6), t559 + t12 * qJD(2) + t18 * qJD(3) + (-t231 * t433 + t233 * t369) * t599 + t63 * qJD(5) + t20 * qJD(6), qJD(2) * t14 + qJD(3) * t17 + qJD(4) * t63 + qJD(6) * t65 + t561, t565 + t5 * qJD(2) + t7 * qJD(3) + t20 * qJD(4) + t65 * qJD(5) + (-mrSges(7,1) * t47 - mrSges(7,2) * t46 + t508) * qJD(6); qJD(3) * t2 + qJD(4) * t11 - qJD(5) * t13 + qJD(6) * t6 - t580, qJD(3) * t34 + qJD(5) * t74 + qJD(6) * t35, t56 * qJD(5) + t30 * qJD(6) - t427 + (m(7) * (t293 * t305 - t294 * t306) + t654 + (t396 * t628 + t651 + t682) * t357 + (m(6) * t383 + m(7) * t374 - t563 * t628 + t640) * t356) * qJD(3), t492 + t562, qJD(3) * t56 + t451, t30 * qJD(3) + (-mrSges(7,1) * t258 - mrSges(7,2) * t257 + t505) * qJD(6) - t429; -qJD(2) * t2 + qJD(4) * t15 - qJD(5) * t16 + qJD(6) * t8 - t568, -t55 * qJD(5) + t29 * qJD(6) + t427, -qJD(5) * t99 - qJD(6) * t45, -t435, t430 (-mrSges(7,1) * t306 - mrSges(7,2) * t305 + t505) * qJD(6) + t432; -qJD(2) * t11 - qJD(3) * t15 - qJD(5) * t62 - qJD(6) * t19 - t559, t492 - t562, t435, 0, -t512, -t436 - t319; qJD(2) * t13 + qJD(3) * t16 + qJD(4) * t62 + qJD(6) * t64 - t561, qJD(3) * t55 + t319 - t451, -t430 + t319, t512, 0, t637; -qJD(2) * t6 - qJD(3) * t8 + qJD(4) * t19 - qJD(5) * t64 - t565, -t29 * qJD(3) - t320 + t429, -t320 - t432, t436, -t637, 0;];
Cq  = t21;
