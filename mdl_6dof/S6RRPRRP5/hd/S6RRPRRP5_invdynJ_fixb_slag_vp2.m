% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:57:51
% EndTime: 2019-03-09 11:59:14
% DurationCPUTime: 51.67s
% Computational Cost: add. (19761->940), mult. (54732->1257), div. (0->0), fcn. (44777->12), ass. (0->419)
t337 = sin(pkin(6));
t336 = sin(pkin(11));
t342 = sin(qJ(2));
t346 = cos(qJ(2));
t512 = cos(pkin(11));
t363 = -t342 * t336 + t346 * t512;
t353 = t337 * t363;
t263 = qJD(1) * t353;
t601 = -t263 + qJD(4);
t639 = Ifges(6,4) + Ifges(7,4);
t464 = qJD(1) * qJD(2);
t281 = (qJDD(1) * t346 - t342 * t464) * t337;
t282 = (qJDD(1) * t342 + t346 * t464) * t337;
t215 = t336 * t281 + t282 * t512;
t419 = t342 * t512;
t477 = qJD(1) * t337;
t444 = t346 * t477;
t264 = -t336 * t444 - t419 * t477;
t338 = cos(pkin(6));
t325 = qJD(1) * t338 + qJD(2);
t341 = sin(qJ(4));
t345 = cos(qJ(4));
t234 = t264 * t341 + t325 * t345;
t462 = qJDD(1) * t338;
t324 = qJDD(2) + t462;
t133 = qJD(4) * t234 + t215 * t345 + t324 * t341;
t578 = t133 / 0.2e1;
t671 = Ifges(5,4) * t578;
t370 = t264 * t345 - t325 * t341;
t134 = qJD(4) * t370 - t215 * t341 + t324 * t345;
t132 = qJDD(5) - t134;
t259 = -t363 * t477 + qJD(4);
t340 = sin(qJ(5));
t344 = cos(qJ(5));
t164 = t259 * t340 - t344 * t370;
t214 = t281 * t512 - t336 * t282;
t213 = qJDD(4) - t214;
t463 = qJDD(1) * t337;
t250 = -pkin(1) * t463 - pkin(2) * t281 + qJDD(3);
t111 = -pkin(3) * t214 - pkin(9) * t215 + t250;
t552 = pkin(1) * t338;
t329 = t346 * t552;
t322 = qJD(1) * t329;
t539 = pkin(8) + qJ(3);
t431 = t539 * t342;
t408 = t337 * t431;
t254 = -qJD(1) * t408 + t322;
t236 = pkin(2) * t325 + t254;
t497 = t338 * t342;
t328 = pkin(1) * t497;
t499 = t337 * t346;
t255 = (t499 * t539 + t328) * qJD(1);
t417 = t512 * t255;
t158 = t336 * t236 + t417;
t149 = pkin(9) * t325 + t158;
t335 = pkin(2) * t346 + pkin(1);
t287 = -t335 * t477 + qJD(3);
t173 = -pkin(3) * t263 + pkin(9) * t264 + t287;
t470 = qJD(4) * t345;
t472 = qJD(4) * t341;
t478 = pkin(8) * t499 + t328;
t278 = t478 * qJD(2);
t456 = pkin(1) * t462;
t320 = t346 * t456;
t501 = t337 * t342;
t442 = qJD(3) * t501;
t455 = pkin(8) * t463;
t154 = -t342 * t455 + pkin(2) * t324 - qJ(3) * t282 + t320 + (-t278 - t442) * qJD(1);
t461 = qJD(2) * t552;
t409 = qJD(1) * t461;
t447 = t342 * t456 + (t409 + t455) * t346;
t474 = qJD(3) * t346;
t475 = qJD(2) * t342;
t166 = qJ(3) * t281 + (-pkin(8) * t475 + t474) * t477 + t447;
t88 = t336 * t154 + t512 * t166;
t82 = pkin(9) * t324 + t88;
t23 = t341 * t111 - t149 * t472 + t173 * t470 + t345 * t82;
t17 = pkin(10) * t213 + t23;
t87 = t154 * t512 - t336 * t166;
t81 = -t324 * pkin(3) - t87;
t30 = -t134 * pkin(4) - t133 * pkin(10) + t81;
t90 = t149 * t345 + t173 * t341;
t77 = pkin(10) * t259 + t90;
t244 = t336 * t255;
t157 = t236 * t512 - t244;
t148 = -t325 * pkin(3) - t157;
t84 = -t234 * pkin(4) + pkin(10) * t370 + t148;
t33 = t340 * t84 + t344 * t77;
t4 = -qJD(5) * t33 - t17 * t340 + t344 * t30;
t163 = t259 * t344 + t340 * t370;
t61 = qJD(5) * t163 + t133 * t344 + t213 * t340;
t1 = pkin(5) * t132 - qJ(6) * t61 - qJD(6) * t164 + t4;
t567 = t213 / 0.2e1;
t577 = t134 / 0.2e1;
t579 = t132 / 0.2e1;
t62 = -qJD(5) * t164 - t133 * t340 + t213 * t344;
t583 = t62 / 0.2e1;
t584 = t61 / 0.2e1;
t635 = -Ifges(7,3) - Ifges(6,3);
t636 = Ifges(6,6) + Ifges(7,6);
t638 = -Ifges(7,5) - Ifges(6,5);
t467 = qJD(5) * t344;
t469 = qJD(5) * t340;
t3 = t344 * t17 + t340 * t30 + t84 * t467 - t469 * t77;
t2 = qJ(6) * t62 + qJD(6) * t163 + t3;
t663 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t670 = t1 * mrSges(7,1) - 0.2e1 * Ifges(5,2) * t577 - 0.2e1 * Ifges(5,6) * t567 - t635 * t579 + t636 * t583 - t638 * t584 + t663 - t671;
t640 = Ifges(6,1) + Ifges(7,1);
t660 = -t579 * t638 + t583 * t639 + t584 * t640;
t334 = pkin(5) * t344 + pkin(4);
t669 = m(7) * t334;
t637 = Ifges(6,2) + Ifges(7,2);
t178 = t254 * t336 + t417;
t668 = -t178 + t601 * (pkin(4) * t341 - pkin(10) * t345);
t339 = -qJ(6) - pkin(10);
t667 = -m(7) * t339 + mrSges(6,3) + mrSges(7,3);
t633 = t132 * t636 + t61 * t639 + t62 * t637;
t666 = t633 / 0.2e1;
t665 = t639 * t163;
t495 = t340 * t345;
t198 = -t263 * t495 - t264 * t344;
t439 = t340 * t470;
t616 = t341 * t467 + t198 + t439;
t664 = t639 * t164;
t648 = Ifges(5,1) * t578 + Ifges(5,5) * t567;
t585 = Ifges(5,4) * t577 + t648;
t232 = qJD(5) - t234;
t630 = t163 * t636 - t164 * t638 - t232 * t635;
t661 = t630 / 0.2e1;
t634 = -t132 * t635 - t61 * t638 + t62 * t636;
t629 = t163 * t637 + t232 * t636 + t664;
t628 = t164 * t640 - t232 * t638 + t665;
t551 = pkin(2) * t336;
t332 = pkin(9) + t551;
t440 = t332 * t470;
t179 = t254 * t512 - t244;
t445 = t342 * t477;
t413 = pkin(2) * t445;
t191 = -pkin(3) * t264 - pkin(9) * t263 + t413;
t105 = -t341 * t179 + t191 * t345;
t91 = pkin(4) * t264 - t105;
t659 = t440 - t91;
t441 = t332 * t472;
t106 = t345 * t179 + t341 * t191;
t92 = -pkin(10) * t264 + t106;
t658 = t668 * t344 + (t441 + t92) * t340;
t446 = t512 * pkin(2);
t333 = -t446 - pkin(3);
t399 = t345 * pkin(4) + t341 * pkin(10);
t297 = -t399 + t333;
t657 = t297 * t467 + t340 * t668 - t344 * t92;
t300 = -t346 * t336 - t419;
t343 = sin(qJ(1));
t347 = cos(qJ(1));
t352 = t338 * t363;
t225 = t300 * t347 - t343 * t352;
t479 = t300 * t338;
t226 = t343 * t479 + t347 * t363;
t489 = t343 * t346;
t492 = t342 * t347;
t290 = -t338 * t489 - t492;
t357 = t290 * pkin(2);
t656 = t225 * pkin(3) + pkin(9) * t226 + t357;
t32 = -t340 * t77 + t344 * t84;
t397 = t3 * t344 - t340 * t4;
t655 = -t32 * t467 - t33 * t469 + t397;
t486 = t346 * t347;
t490 = t343 * t342;
t654 = t338 * t486 - t490;
t653 = t639 * t344;
t652 = t639 * t340;
t393 = mrSges(5,1) * t345 - mrSges(5,2) * t341;
t651 = m(6) * t399 + t341 * t667 + t345 * t669 + mrSges(4,1) + t393;
t221 = -t343 * t363 + t347 * t479;
t498 = t337 * t347;
t201 = -t221 * t345 - t341 * t498;
t222 = t343 * t300 + t347 * t352;
t650 = t201 * t340 + t222 * t344;
t649 = -t201 * t344 + t222 * t340;
t566 = -t232 / 0.2e1;
t574 = -t164 / 0.2e1;
t576 = -t163 / 0.2e1;
t647 = -t566 * t635 - t574 * t638 + t576 * t636;
t586 = m(7) * pkin(5);
t646 = -m(6) - m(5);
t89 = -t341 * t149 + t173 * t345;
t76 = -pkin(4) * t259 - t89;
t645 = m(6) * t76;
t644 = t89 * mrSges(5,1);
t643 = t90 * mrSges(5,2);
t541 = -mrSges(6,1) - mrSges(7,1);
t642 = mrSges(6,2) + mrSges(7,2);
t641 = -mrSges(5,3) + mrSges(4,2);
t21 = -mrSges(6,1) * t62 + mrSges(6,2) * t61;
t94 = mrSges(5,1) * t213 - mrSges(5,3) * t133;
t631 = t21 - t94;
t626 = pkin(5) * t616 + t659;
t420 = qJD(5) * t339;
t466 = qJD(6) * t344;
t507 = t234 * t340;
t153 = -pkin(4) * t370 - pkin(10) * t234;
t55 = t340 * t153 + t344 * t89;
t625 = qJ(6) * t507 + t340 * t420 + t466 - t55;
t506 = t234 * t344;
t54 = t344 * t153 - t340 * t89;
t624 = pkin(5) * t370 + qJ(6) * t506 - qJD(6) * t340 + t344 * t420 - t54;
t488 = t344 * t345;
t199 = t263 * t488 - t264 * t340;
t302 = t332 * t488;
t505 = t263 * t341;
t623 = -pkin(5) * t505 + qJ(6) * t199 - t341 * t466 + (pkin(5) * t341 - qJ(6) * t488) * qJD(4) + (-t302 + (qJ(6) * t341 - t297) * t340) * qJD(5) + t658;
t494 = t341 * t344;
t622 = -qJ(6) * t198 + (-qJ(6) * qJD(5) - qJD(4) * t332) * t494 + (-qJD(6) * t341 + (-qJ(6) * qJD(4) - qJD(5) * t332) * t345) * t340 + t657;
t621 = -t90 + (t469 - t507) * pkin(5);
t243 = t340 * t297 + t302;
t620 = -qJD(5) * t243 + t658;
t471 = qJD(4) * t344;
t619 = (-t341 * t471 - t345 * t469) * t332 + t657;
t618 = t586 + mrSges(7,1);
t251 = pkin(2) * t338 + t329 - t408;
t260 = qJ(3) * t499 + t478;
t186 = t251 * t512 - t336 * t260;
t174 = -t338 * pkin(3) - t186;
t271 = t300 * t337;
t247 = -t271 * t341 - t338 * t345;
t248 = -t271 * t345 + t338 * t341;
t107 = t247 * pkin(4) - t248 * pkin(10) + t174;
t187 = t336 * t251 + t512 * t260;
t175 = pkin(9) * t338 + t187;
t326 = pkin(2) * t499;
t401 = pkin(3) * t353 - pkin(9) * t271 + t326;
t553 = pkin(1) * t337;
t197 = -t401 - t553;
t109 = t345 * t175 + t341 * t197;
t99 = -pkin(10) * t353 + t109;
t45 = t340 * t107 + t344 * t99;
t101 = -mrSges(6,1) * t163 + mrSges(6,2) * t164;
t535 = mrSges(5,3) * t370;
t172 = mrSges(5,1) * t259 + t535;
t617 = t101 - t172;
t468 = qJD(5) * t341;
t615 = t340 * t468 - t344 * t470 + t199;
t537 = mrSges(4,3) * t264;
t614 = mrSges(4,1) * t325 + mrSges(5,1) * t234 + mrSges(5,2) * t370 + t537;
t388 = -mrSges(7,1) * t344 + mrSges(7,2) * t340;
t391 = -mrSges(6,1) * t344 + mrSges(6,2) * t340;
t613 = m(6) * pkin(4) - t388 - t391 + t669;
t387 = mrSges(7,1) * t340 + mrSges(7,2) * t344;
t390 = mrSges(6,1) * t340 + mrSges(6,2) * t344;
t51 = -pkin(5) * t163 + qJD(6) + t76;
t612 = t51 * t387 + t76 * t390;
t611 = -t340 * t638 + t344 * t636;
t610 = -t340 * t636 - t344 * t638;
t609 = t344 * t637 + t652;
t608 = -t340 * t637 + t653;
t607 = t340 * t640 + t653;
t606 = t344 * t640 - t652;
t605 = -m(6) * pkin(10) - t667;
t603 = t472 - t505;
t24 = t111 * t345 - t149 * t470 - t173 * t472 - t341 * t82;
t602 = t23 * t345 - t24 * t341;
t465 = m(7) - t646;
t530 = Ifges(3,4) * t342;
t600 = pkin(1) * (mrSges(3,1) * t342 + mrSges(3,2) * t346) - t342 * (Ifges(3,1) * t346 - t530) / 0.2e1;
t599 = -mrSges(6,1) - t618;
t598 = mrSges(5,1) + t613;
t359 = mrSges(5,2) + t605;
t548 = pkin(5) * t340;
t595 = m(7) * t548 + pkin(9) * t465 - t641;
t593 = t24 * mrSges(5,1) - t23 * mrSges(5,2) + Ifges(5,5) * t133 + Ifges(5,6) * t134 + Ifges(5,3) * t213;
t25 = -qJ(6) * t164 + t32;
t19 = pkin(5) * t232 + t25;
t26 = qJ(6) * t163 + t33;
t592 = -t1 * t340 - t19 * t467 + t2 * t344 - t26 * t469;
t443 = t337 * t475;
t412 = pkin(8) * t443;
t229 = -qJD(1) * t412 + t447;
t230 = -pkin(8) * t282 - t342 * t409 + t320;
t591 = t230 * mrSges(3,1) + t87 * mrSges(4,1) - t229 * mrSges(3,2) - t88 * mrSges(4,2) + Ifges(3,5) * t282 + Ifges(4,5) * t215 + Ifges(3,6) * t281 + Ifges(4,6) * t214;
t589 = m(5) / 0.2e1;
t588 = m(6) / 0.2e1;
t587 = m(7) / 0.2e1;
t529 = Ifges(5,4) * t370;
t120 = t234 * Ifges(5,2) + t259 * Ifges(5,6) - t529;
t580 = -t120 / 0.2e1;
t575 = t163 / 0.2e1;
t573 = t164 / 0.2e1;
t565 = t232 / 0.2e1;
t564 = -t234 / 0.2e1;
t563 = t370 / 0.2e1;
t562 = -t370 / 0.2e1;
t560 = -t259 / 0.2e1;
t559 = -t263 / 0.2e1;
t558 = -t264 / 0.2e1;
t550 = pkin(5) * t164;
t193 = -t248 * t340 - t344 * t353;
t549 = pkin(5) * t193;
t538 = mrSges(4,3) * t263;
t536 = mrSges(5,3) * t234;
t534 = mrSges(6,3) * t163;
t533 = mrSges(6,3) * t164;
t532 = mrSges(7,3) * t163;
t531 = mrSges(7,3) * t164;
t528 = Ifges(5,4) * t341;
t527 = Ifges(5,4) * t345;
t18 = -pkin(4) * t213 - t24;
t520 = t18 * t341;
t517 = t264 * Ifges(4,4);
t516 = t325 * Ifges(3,5);
t515 = t325 * Ifges(3,6);
t511 = t221 * t340;
t508 = t226 * t340;
t504 = t263 * t345;
t503 = t271 * t340;
t500 = t337 * t343;
t496 = t340 * t341;
t115 = -mrSges(7,2) * t232 + t532;
t116 = -mrSges(6,2) * t232 + t534;
t484 = t115 + t116;
t117 = mrSges(7,1) * t232 - t531;
t118 = mrSges(6,1) * t232 - t533;
t483 = -t117 - t118;
t473 = qJD(4) * t340;
t459 = m(4) + t465;
t458 = -m(3) * pkin(1) - mrSges(2,1);
t20 = -t62 * mrSges(7,1) + t61 * mrSges(7,2);
t432 = t539 * t337;
t428 = t470 / 0.2e1;
t423 = -t468 / 0.2e1;
t416 = -t214 * mrSges(4,1) + t215 * mrSges(4,2);
t44 = t344 * t107 - t340 * t99;
t108 = -t341 * t175 + t197 * t345;
t323 = t346 * t461;
t237 = t323 + (-qJD(2) * t431 + t474) * t337;
t238 = -t442 + (-t346 * t432 - t328) * qJD(2);
t151 = t237 * t336 - t512 * t238;
t283 = pkin(2) * t497 - t432;
t415 = -t283 * t343 + t347 * t335;
t411 = mrSges(3,3) * t445;
t410 = mrSges(3,3) * t444;
t402 = t654 * pkin(2);
t396 = t226 * pkin(3) + t415;
t394 = mrSges(5,1) * t247 + mrSges(5,2) * t248;
t194 = t248 * t344 - t340 * t353;
t392 = mrSges(6,1) * t193 - mrSges(6,2) * t194;
t389 = -t193 * mrSges(7,1) + t194 * mrSges(7,2);
t386 = Ifges(5,1) * t345 - t528;
t381 = -Ifges(5,2) * t341 + t527;
t376 = Ifges(5,5) * t345 - Ifges(5,6) * t341;
t205 = t226 * t345 + t341 * t500;
t128 = -t205 * t340 - t225 * t344;
t152 = t237 * t512 + t336 * t238;
t265 = qJD(2) * t271;
t266 = qJD(2) * t353;
t192 = pkin(2) * t443 - pkin(3) * t265 - pkin(9) * t266;
t50 = -t341 * t152 - t175 * t470 + t192 * t345 - t197 * t472;
t98 = pkin(4) * t353 - t108;
t200 = -t221 * t341 + t345 * t498;
t49 = t345 * t152 - t175 * t472 + t341 * t192 + t197 * t470;
t42 = -pkin(10) * t265 + t49;
t184 = qJD(4) * t248 + t266 * t341;
t185 = -qJD(4) * t247 + t266 * t345;
t68 = pkin(4) * t184 - pkin(10) * t185 + t151;
t7 = t107 * t467 + t340 * t68 + t344 * t42 - t469 * t99;
t43 = pkin(4) * t265 - t50;
t8 = -qJD(5) * t45 - t340 * t42 + t344 * t68;
t321 = Ifges(3,4) * t444;
t314 = Ifges(3,3) * t324;
t313 = Ifges(4,3) * t324;
t309 = t339 * t344;
t308 = t339 * t340;
t301 = -t326 - t553;
t295 = -pkin(8) * t501 + t329;
t292 = (-mrSges(3,1) * t346 + mrSges(3,2) * t342) * t337;
t291 = -t338 * t490 + t486;
t289 = -t338 * t492 - t489;
t286 = (t332 + t548) * t341;
t280 = t344 * t297;
t277 = t323 - t412;
t276 = t478 * qJD(1);
t275 = -pkin(8) * t445 + t322;
t274 = -mrSges(3,2) * t325 + t410;
t273 = mrSges(3,1) * t325 - t411;
t258 = Ifges(4,4) * t263;
t253 = Ifges(3,1) * t445 + t321 + t516;
t252 = t515 + (t346 * Ifges(3,2) + t530) * t477;
t242 = -t332 * t495 + t280;
t239 = -mrSges(4,2) * t325 + t538;
t233 = -qJ(6) * t496 + t243;
t228 = Ifges(5,4) * t234;
t218 = t221 * pkin(3);
t216 = -qJ(6) * t494 + t280 + (-t332 * t340 - pkin(5)) * t345;
t208 = -mrSges(4,1) * t263 - mrSges(4,2) * t264;
t204 = t226 * t341 - t345 * t500;
t196 = mrSges(4,1) * t324 - mrSges(4,3) * t215;
t195 = -mrSges(4,2) * t324 + mrSges(4,3) * t214;
t190 = -t264 * Ifges(4,1) + t325 * Ifges(4,5) + t258;
t189 = t263 * Ifges(4,2) + t325 * Ifges(4,6) - t517;
t171 = -mrSges(5,2) * t259 + t536;
t129 = t205 * t344 - t225 * t340;
t121 = -Ifges(5,1) * t370 + t259 * Ifges(5,5) + t228;
t119 = -Ifges(5,5) * t370 + t234 * Ifges(5,6) + t259 * Ifges(5,3);
t100 = -mrSges(7,1) * t163 + mrSges(7,2) * t164;
t97 = qJD(5) * t193 + t185 * t344 - t265 * t340;
t96 = -qJD(5) * t194 - t185 * t340 - t265 * t344;
t95 = -mrSges(5,2) * t213 + mrSges(5,3) * t134;
t67 = -mrSges(5,1) * t134 + mrSges(5,2) * t133;
t64 = t98 - t549;
t39 = -mrSges(6,2) * t132 + mrSges(6,3) * t62;
t38 = -mrSges(7,2) * t132 + mrSges(7,3) * t62;
t37 = mrSges(6,1) * t132 - mrSges(6,3) * t61;
t36 = mrSges(7,1) * t132 - mrSges(7,3) * t61;
t34 = qJ(6) * t193 + t45;
t27 = pkin(5) * t247 - qJ(6) * t194 + t44;
t22 = -pkin(5) * t96 + t43;
t9 = -pkin(5) * t62 + qJDD(6) + t18;
t6 = qJ(6) * t96 + qJD(6) * t193 + t7;
t5 = pkin(5) * t184 - qJ(6) * t97 - qJD(6) * t194 + t8;
t10 = [(-t289 * mrSges(3,1) + t654 * mrSges(3,2) - m(5) * t218 + t201 * mrSges(5,1) - m(6) * (-pkin(4) * t201 + t218) - m(7) * (-t201 * t334 + t218) - t221 * mrSges(4,1) + (t283 * t459 + mrSges(2,2)) * t347 + (t335 * t459 - t458) * t343 + t541 * t649 - t642 * t650 - t595 * t222 - t359 * t200) * g(1) + (-t671 - t23 * mrSges(5,3) + t634 / 0.2e1 + t670) * t247 + (-mrSges(6,3) * t4 - mrSges(7,3) * t1 + 0.2e1 * t660) * t194 + m(3) * (t229 * t478 + t230 * t295 - t275 * t278 + t276 * t277) + (t313 / 0.2e1 + t314 / 0.2e1 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t324 + t591) * t338 - t265 * t644 + ((mrSges(3,1) * t281 - mrSges(3,2) * t282 + (m(3) * t553 - t292) * qJDD(1)) * pkin(1) + (mrSges(3,3) * t229 + Ifges(3,4) * t282 + Ifges(3,2) * t281 + Ifges(3,6) * t324) * t346 + (-mrSges(3,3) * t230 + Ifges(3,1) * t282 + Ifges(3,4) * t281 + Ifges(3,5) * t324) * t342 + ((t253 / 0.2e1 + t516 / 0.2e1 - t275 * mrSges(3,3)) * t346 + (-t276 * mrSges(3,3) - t515 / 0.2e1 - t252 / 0.2e1 + (m(4) * t287 + t208) * pkin(2)) * t342 + (t346 * (Ifges(3,4) * t346 - Ifges(3,2) * t342) / 0.2e1 - t600) * t477) * qJD(2) + (g(1) * t347 + g(2) * t343) * (-m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3))) * t337 + (-t184 * t90 - t185 * t89) * mrSges(5,3) + t184 * t580 + m(5) * (t108 * t24 + t109 * t23 + t148 * t151 + t174 * t81 + t49 * t90 + t50 * t89) + m(7) * (t1 * t27 + t19 * t5 + t2 * t34 + t22 * t51 + t26 * t6 + t64 * t9) + m(6) * (t18 * t98 + t3 * t45 + t32 * t8 + t33 * t7 + t4 * t44 + t43 * t76) + (-t157 * t266 + t158 * t265) * mrSges(4,3) + t325 * (Ifges(4,5) * t266 + Ifges(4,6) * t265) / 0.2e1 + t287 * (-mrSges(4,1) * t265 + mrSges(4,2) * t266) + t263 * (Ifges(4,4) * t266 + Ifges(4,2) * t265) / 0.2e1 + t259 * (Ifges(5,5) * t185 - Ifges(5,6) * t184 - Ifges(5,3) * t265) / 0.2e1 + t234 * (Ifges(5,4) * t185 - Ifges(5,2) * t184 - Ifges(5,6) * t265) / 0.2e1 + (Ifges(4,1) * t266 + Ifges(4,4) * t265) * t558 + (Ifges(5,1) * t185 - Ifges(5,4) * t184 - Ifges(5,5) * t265) * t562 + t628 * t97 / 0.2e1 + t629 * t96 / 0.2e1 + t301 * t416 + (-t184 * t635 + t636 * t96 - t638 * t97) * t565 + (t184 * t636 + t637 * t96 + t639 * t97) * t575 + (-t184 * t638 + t639 * t96 + t640 * t97) * t573 + t9 * t389 - t18 * t392 + t81 * t394 - (mrSges(4,2) * t250 - mrSges(4,3) * t87 + Ifges(4,1) * t215 + Ifges(4,4) * t214 + Ifges(4,5) * t324) * t271 + t184 * t661 + t478 * (-mrSges(3,2) * t324 + mrSges(3,3) * t281) + t45 * t39 + t44 * t37 + t34 * t38 + t27 * t36 + Ifges(2,3) * qJDD(1) + t295 * (mrSges(3,1) * t324 - mrSges(3,3) * t282) + t277 * t274 - t278 * t273 + t266 * t190 / 0.2e1 - t265 * t119 / 0.2e1 + t265 * t189 / 0.2e1 + t152 * t239 + (t3 * mrSges(6,3) + t2 * mrSges(7,3) + t636 * t579 + t637 * t583 + t639 * t584 + t666) * t193 + t187 * t195 + t186 * t196 + t33 * (-mrSges(6,2) * t184 + mrSges(6,3) * t96) + t26 * (-mrSges(7,2) * t184 + mrSges(7,3) * t96) + t32 * (mrSges(6,1) * t184 - mrSges(6,3) * t97) + t19 * (mrSges(7,1) * t184 - mrSges(7,3) * t97) + t148 * (mrSges(5,1) * t184 + mrSges(5,2) * t185) + t185 * t121 / 0.2e1 - (t250 * mrSges(4,1) - t88 * mrSges(4,3) - Ifges(4,4) * t215 - Ifges(4,2) * t214 - Ifges(4,6) * t324 + t593) * t353 - t614 * t151 + (-t291 * mrSges(3,1) - t290 * mrSges(3,2) + t343 * mrSges(2,2) - m(6) * (pkin(4) * t205 + t396) - m(5) * t396 - t205 * mrSges(5,1) - m(7) * (t205 * t334 + t396) - m(4) * t415 - t226 * mrSges(4,1) + t458 * t347 + t541 * t129 - t642 * t128 + t595 * t225 + t359 * t204) * g(2) + t64 * t20 + t265 * t643 + m(4) * (-t151 * t157 + t152 * t158 + t186 * t87 + t187 * t88 + t250 * t301) + (-mrSges(5,3) * t24 + 0.2e1 * t585) * t248 + t51 * (-mrSges(7,1) * t96 + mrSges(7,2) * t97) + t76 * (-mrSges(6,1) * t96 + mrSges(6,2) * t97) + t98 * t21 + t22 * t100 + t43 * t101 + t108 * t94 + t109 * t95 + t6 * t115 + t7 * t116 + t5 * t117 + t8 * t118 + t49 * t171 + t50 * t172 + t174 * t67; (t332 * t95 - t670) * t345 + (t511 * t586 - m(4) * t402 - mrSges(3,1) * t654 - mrSges(3,2) * t289 - t465 * (t222 * pkin(3) - pkin(9) * t221 + t402) - t641 * t221 + t541 * (t222 * t488 - t511) - t642 * (-t221 * t344 - t222 * t495) - t651 * t222) * g(2) + (-m(7) * (pkin(5) * t508 + t656) - m(4) * t357 - mrSges(3,1) * t290 + mrSges(3,2) * t291 + t646 * t656 + t641 * t226 + t541 * (t225 * t488 + t508) - t642 * (-t225 * t495 + t226 * t344) - t651 * t225) * g(1) + (t411 + t273) * t276 + (-mrSges(6,2) * t603 - mrSges(6,3) * t616) * t33 + (-mrSges(7,2) * t603 - mrSges(7,3) * t616) * t26 + t614 * t178 + (mrSges(6,1) * t603 + mrSges(6,3) * t615) * t32 + (mrSges(7,1) * t603 + mrSges(7,3) * t615) * t19 + (mrSges(7,1) * t616 - mrSges(7,2) * t615) * t51 + (mrSges(6,1) * t616 - mrSges(6,2) * t615) * t76 + t264 * t644 + (-t105 * t89 - t106 * t90 - t148 * t178 + t333 * t81 + ((-t341 * t90 - t345 * t89) * qJD(4) + t602) * t332) * m(5) + (-t603 * t90 + (-t470 + t504) * t89 + t602) * mrSges(5,3) - ((-Ifges(3,2) * t445 + t253 + t321) * t346 + t325 * (Ifges(3,5) * t346 - Ifges(3,6) * t342)) * t477 / 0.2e1 + (t234 * t381 + t259 * t376 - t370 * t386) * qJD(4) / 0.2e1 + t591 + (-t1 * t494 - t2 * t496) * mrSges(7,3) + t619 * t116 + t620 * t118 + (-t76 * t91 + t242 * t4 + t243 * t3 + (t470 * t76 + t520) * t332 + t619 * t33 + t620 * t32) * m(6) + t622 * t115 + t623 * t117 + t626 * t100 + (t1 * t216 + t19 * t623 + t2 * t233 + t26 * t622 + t286 * t9 + t51 * t626) * m(7) + (-t3 * t496 - t4 * t494) * mrSges(6,3) + (t410 - t274) * t275 + (t428 - t504 / 0.2e1) * t121 + t196 * t446 + t528 * t577 + t631 * t332 * t341 + t628 * (t340 * t423 + t344 * t428 - t199 / 0.2e1) + t629 * (t344 * t423 - t439 / 0.2e1 - t198 / 0.2e1) + (-t441 - t106) * t171 + t527 * t578 + t600 * qJD(1) ^ 2 * t337 ^ 2 + t601 * t148 * (mrSges(5,1) * t341 + mrSges(5,2) * t345) + (t198 * t636 - t199 * t638) * t566 + (t198 * t637 + t199 * t639) * t576 + (t198 * t639 + t199 * t640) * t574 + (-t440 - t105) * t172 + ((t336 * t88 + t512 * t87) * pkin(2) + t157 * t178 - t158 * t179 - t287 * t413) * m(4) - t158 * t537 + t314 + t313 - t81 * t393 + t494 * t660 + t252 * t445 / 0.2e1 + (Ifges(4,2) * t264 + t190 + t258) * t559 - t208 * t413 + (Ifges(4,1) * t263 + t119 + t517) * t264 / 0.2e1 + (t580 + t661) * t472 + t659 * t101 + (t387 * t9 + t579 * t610 + t583 * t608 + t584 * t606 + t585 + t648) * t341 + (t120 / 0.2e1 - t630 / 0.2e1 + t647) * t505 + t333 * t67 - t325 * (Ifges(4,5) * t263 + Ifges(4,6) * t264) / 0.2e1 + t286 * t20 - t287 * (-mrSges(4,1) * t264 + mrSges(4,2) * t263) + t242 * t37 + t243 * t39 - t179 * t239 + t233 * t38 + t216 * t36 + (-Ifges(5,6) * t264 + t263 * t381) * t564 + t189 * t558 + (-Ifges(5,3) * t264 + t263 * t376) * t560 + (-Ifges(5,5) * t264 + t263 * t386) * t563 + t195 * t551 + t157 * t538 + (t503 * t586 - m(4) * t326 + t292 - t465 * t401 - t641 * t271 + t541 * (t353 * t488 - t503) - t642 * (-t271 * t344 - t353 * t495) - t651 * t353) * g(3) + t390 * t520 - t633 * t496 / 0.2e1 - t634 * t345 / 0.2e1 + (-t611 * t468 + (-t341 * t635 + t345 * t610) * qJD(4)) * t565 + (-t609 * t468 + (t341 * t636 + t345 * t608) * qJD(4)) * t575 + (-t607 * t468 + (-t341 * t638 + t345 * t606) * qJD(4)) * t573 - t264 * t643; -t263 * t239 - t484 * t199 + t483 * t198 + (-t263 * t171 - t20 + (t340 * t483 + t344 * t484 + t171) * qJD(4) - t631) * t345 + (t95 + (t38 + t39) * t344 + (-t36 - t37) * t340 + (-t340 * t484 + t344 * t483) * qJD(5) + t601 * (t100 + t617)) * t341 - m(7) * (t19 * t198 + t199 * t26) - m(6) * (t198 * t32 + t199 * t33) + 0.2e1 * ((-t19 * t473 + t26 * t471 - t9) * t587 + (-t32 * t473 + t33 * t471 - t18) * t588 + m(5) * t90 * t559 + (qJD(4) * t90 + t24) * t589) * t345 + 0.2e1 * ((qJD(4) * t51 + t592) * t587 + (qJD(4) * t76 + t655) * t588 + (-qJD(4) * t89 + t23) * t589 + (t89 * t589 - m(7) * t51 / 0.2e1 - t645 / 0.2e1) * t263) * t341 + t416 + (-t338 * g(3) + (-g(1) * t343 + g(2) * t347) * t337) * t459 + (m(5) * t148 - t614) * t264 + (-t157 * t264 - t158 * t263 + t250) * m(4); (t19 * t506 + t26 * t507 + t592) * mrSges(7,3) + t609 * t583 + t611 * t579 + (-t148 * mrSges(5,2) + Ifges(5,1) * t563 + Ifges(5,5) * t560 + t566 * t610 + t574 * t606 + t576 * t608 - t612) * t234 + t612 * qJD(5) + (t247 * t613 + t248 * t605 + t394) * g(3) + (t200 * t598 + t201 * t359) * g(2) + (t204 * t598 + t205 * t359) * g(1) + (t163 * t608 + t164 * t606 + t232 * t610) * qJD(5) / 0.2e1 + t593 + t621 * t100 + t624 * t117 + (t1 * t308 + t19 * t624 - t2 * t309 + t26 * t625 - t334 * t9 + t51 * t621) * m(7) + t625 * t115 + (-t469 / 0.2e1 + t507 / 0.2e1) * t629 + t344 * t666 + (t467 / 0.2e1 - t506 / 0.2e1) * t628 + t9 * t388 + t18 * t391 + (m(6) * ((-t32 * t344 - t33 * t340) * qJD(5) + t397) - t118 * t467 - t116 * t469 + t344 * t39 - t340 * t37) * pkin(10) + t340 * t660 + (t228 + t121) * t564 + (t32 * t506 + t33 * t507 + t655) * mrSges(6,3) - pkin(4) * t21 - (-t148 * mrSges(5,1) - t32 * mrSges(6,1) - t19 * mrSges(7,1) + t33 * mrSges(6,2) + t26 * mrSges(7,2) - Ifges(5,2) * t564 - Ifges(5,6) * t560 + t647) * t370 + (t536 - t171) * t89 - t334 * t20 + t308 * t36 - t309 * t38 + t120 * t562 + (-pkin(4) * t18 - t32 * t54 - t33 * t55) * m(6) + (t529 + t630) * t563 + (-t535 - t617 - t645) * t90 - t55 * t116 - t54 * t118 + t607 * t584; (-t163 * t638 - t164 * t636) * t566 + (t163 * t640 - t664) * t574 + t629 * t573 + (-m(7) * t550 - mrSges(7,1) * t164 - mrSges(7,2) * t163) * t51 + (t533 + t118) * t33 + (t534 - t116) * t32 + (-m(7) * (-t19 + t25) + t531 + t117) * t26 + t618 * t1 + (-m(7) * t549 + t389 - t392) * g(3) + t663 + (-t599 * t650 - t642 * t649) * g(2) + (t128 * t599 + t129 * t642) * g(1) - t100 * t550 + (-t164 * t637 + t628 + t665) * t576 + pkin(5) * t36 + t19 * t532 - t25 * t115 - t76 * (mrSges(6,1) * t164 + mrSges(6,2) * t163) + t634; -t163 * t115 + t164 * t117 + (-g(1) * t204 - g(2) * t200 - g(3) * t247 - t26 * t163 + t19 * t164 + t9) * m(7) + t20;];
tau  = t10;
