% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP10
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
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP10_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP10_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP10_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:36:35
% EndTime: 2019-03-09 12:37:57
% DurationCPUTime: 53.73s
% Computational Cost: add. (22056->970), mult. (53717->1281), div. (0->0), fcn. (44039->14), ass. (0->418)
t377 = sin(pkin(11));
t379 = cos(pkin(11));
t383 = sin(qJ(4));
t387 = cos(qJ(4));
t328 = t377 * t387 + t379 * t383;
t378 = sin(pkin(6));
t388 = cos(qJ(2));
t506 = t378 * t388;
t391 = t328 * t506;
t255 = qJD(1) * t391;
t318 = t328 * qJD(4);
t625 = t255 - t318;
t327 = t377 * t383 - t387 * t379;
t390 = t327 * t506;
t256 = qJD(1) * t390;
t317 = t327 * qJD(4);
t624 = -t256 + t317;
t655 = Ifges(6,1) + Ifges(7,1);
t653 = Ifges(7,4) + Ifges(6,5);
t384 = sin(qJ(2));
t421 = pkin(2) * t384 - qJ(3) * t388;
t494 = qJD(1) * t378;
t307 = t421 * t494;
t468 = t384 * t494;
t380 = cos(pkin(6));
t493 = qJD(1) * t380;
t482 = pkin(1) * t493;
t308 = -pkin(8) * t468 + t388 * t482;
t209 = t379 * t307 - t377 * t308;
t504 = t379 * t388;
t393 = (pkin(3) * t384 - pkin(9) * t504) * t378;
t168 = qJD(1) * t393 + t209;
t210 = t377 * t307 + t379 * t308;
t465 = t388 * t494;
t442 = t377 * t465;
t189 = -pkin(9) * t442 + t210;
t547 = pkin(9) + qJ(3);
t342 = t547 * t377;
t343 = t547 * t379;
t623 = -t387 * t342 - t343 * t383;
t630 = -t327 * qJD(3) + qJD(4) * t623 - t383 * t168 - t387 * t189;
t490 = qJD(2) * t388;
t313 = (qJD(1) * t490 + qJDD(1) * t384) * t378;
t483 = qJDD(1) * t380;
t363 = qJDD(2) + t483;
t236 = -t313 * t377 + t363 * t379;
t237 = t313 * t379 + t363 * t377;
t364 = qJD(2) + t493;
t276 = t364 * t379 - t377 * t468;
t277 = t364 * t377 + t379 * t468;
t451 = t387 * t276 - t277 * t383;
t111 = qJD(4) * t451 + t236 * t383 + t237 * t387;
t491 = qJD(2) * t378;
t467 = t384 * t491;
t667 = -qJD(1) * t467 + qJDD(1) * t506;
t295 = qJDD(4) - t667;
t382 = sin(qJ(5));
t386 = cos(qJ(5));
t407 = t276 * t383 + t387 * t277;
t418 = -qJD(4) + t465;
t161 = t382 * t407 + t386 * t418;
t486 = qJD(5) * t161;
t62 = t386 * t111 + t382 * t295 - t486;
t599 = t62 / 0.2e1;
t162 = -t382 * t418 + t386 * t407;
t63 = qJD(5) * t162 + t382 * t111 - t386 * t295;
t598 = -t63 / 0.2e1;
t112 = -qJD(4) * t407 + t236 * t387 - t237 * t383;
t110 = qJDD(5) - t112;
t591 = t110 / 0.2e1;
t654 = -Ifges(6,4) + Ifges(7,5);
t652 = Ifges(6,6) - Ifges(7,6);
t651 = Ifges(6,3) + Ifges(7,2);
t678 = pkin(10) * t468 - t630;
t359 = pkin(8) * t465;
t309 = t384 * t482 + t359;
t248 = pkin(3) * t442 + t309;
t677 = -pkin(4) * t625 + t624 * pkin(10) - t248;
t590 = t111 / 0.2e1;
t589 = t112 / 0.2e1;
t569 = t295 / 0.2e1;
t649 = t110 * t653 + t62 * t655 + t63 * t654;
t676 = t649 / 0.2e1;
t156 = Ifges(6,4) * t161;
t196 = qJD(5) - t451;
t534 = Ifges(7,5) * t161;
t640 = t162 * t655 + t653 * t196 - t156 + t534;
t254 = -t342 * t383 + t343 * t387;
t629 = -qJD(3) * t328 - qJD(4) * t254 - t168 * t387 + t383 * t189;
t675 = -t62 * Ifges(7,5) / 0.2e1 - t110 * Ifges(7,6) / 0.2e1 + Ifges(6,4) * t599 + Ifges(6,6) * t591 + (Ifges(7,3) + Ifges(6,2)) * t598;
t597 = t63 / 0.2e1;
t674 = Ifges(6,2) * t598 - Ifges(7,3) * t597 + t652 * t591;
t673 = t591 * t653 + t599 * t655;
t575 = t407 / 0.2e1;
t577 = t451 / 0.2e1;
t195 = Ifges(5,4) * t451;
t639 = Ifges(5,5) * t418;
t126 = Ifges(5,1) * t407 + t195 - t639;
t587 = t126 / 0.2e1;
t657 = -t418 / 0.2e1;
t672 = -Ifges(5,1) * t575 - Ifges(5,4) * t577 - Ifges(5,5) * t657 - t587;
t508 = t378 * t384;
t365 = pkin(8) * t508;
t559 = pkin(1) * t380;
t481 = qJD(2) * t559;
t447 = qJD(1) * t481;
t477 = pkin(1) * t483;
t223 = -qJD(2) * t359 - qJDD(1) * t365 - t384 * t447 + t388 * t477;
t203 = -t363 * pkin(2) + qJDD(3) - t223;
t141 = -t236 * pkin(3) + t203;
t259 = qJ(3) * t364 + t309;
t290 = (-pkin(2) * t388 - qJ(3) * t384 - pkin(1)) * t378;
t264 = qJD(1) * t290;
t171 = -t377 * t259 + t379 * t264;
t132 = -pkin(3) * t465 - t277 * pkin(9) + t171;
t172 = t379 * t259 + t377 * t264;
t142 = pkin(9) * t276 + t172;
t487 = qJD(4) * t387;
t488 = qJD(4) * t383;
t222 = pkin(8) * t667 + t384 * t477 + t388 * t447;
t180 = qJ(3) * t363 + qJD(3) * t364 + t222;
t489 = qJD(3) * t384;
t192 = -pkin(2) * t667 - qJ(3) * t313 + (-pkin(1) * qJDD(1) - qJD(1) * t489) * t378;
t122 = -t180 * t377 + t379 * t192;
t80 = -pkin(3) * t667 - pkin(9) * t237 + t122;
t123 = t379 * t180 + t377 * t192;
t89 = pkin(9) * t236 + t123;
t20 = t132 * t487 - t142 * t488 + t383 * t80 + t387 * t89;
t18 = pkin(10) * t295 + t20;
t35 = -t112 * pkin(4) - t111 * pkin(10) + t141;
t484 = qJD(5) * t386;
t485 = qJD(5) * t382;
t76 = t383 * t132 + t387 * t142;
t73 = -pkin(10) * t418 + t76;
t250 = -t364 * pkin(2) + qJD(3) - t308;
t197 = -t276 * pkin(3) + t250;
t85 = -pkin(4) * t451 - pkin(10) * t407 + t197;
t3 = t386 * t18 + t382 * t35 + t85 * t484 - t485 * t73;
t1 = qJ(6) * t110 + qJD(6) * t196 + t3;
t27 = t382 * t85 + t386 * t73;
t4 = -qJD(5) * t27 - t18 * t382 + t35 * t386;
t2 = -pkin(5) * t110 + qJDD(6) - t4;
t605 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t650 = t110 * t651 + t62 * t653 - t63 * t652;
t671 = t591 * t651 + t599 * t653 + t605 + t141 * mrSges(5,1) - t20 * mrSges(5,3) - 0.2e1 * Ifges(5,4) * t590 - 0.2e1 * Ifges(5,2) * t589 - 0.2e1 * Ifges(5,6) * t569 + Ifges(6,6) * t598 + Ifges(7,6) * t597 + t650 / 0.2e1;
t540 = Ifges(5,4) * t407;
t638 = Ifges(5,6) * t418;
t125 = Ifges(5,2) * t451 + t540 - t638;
t579 = t196 / 0.2e1;
t582 = t162 / 0.2e1;
t584 = t161 / 0.2e1;
t585 = -t161 / 0.2e1;
t641 = -t161 * t652 + t162 * t653 + t196 * t651;
t670 = -t125 / 0.2e1 + Ifges(6,6) * t585 + Ifges(7,6) * t584 - Ifges(5,4) * t575 - Ifges(5,2) * t577 + t651 * t579 + t653 * t582 + t641 / 0.2e1 - Ifges(5,6) * t657;
t372 = pkin(3) * t379 + pkin(2);
t226 = pkin(4) * t327 - pkin(10) * t328 - t372;
t626 = t382 * t226 + t386 * t254;
t644 = -qJD(5) * t626 + t382 * t678 + t386 * t677;
t643 = t226 * t484 - t254 * t485 + t382 * t677 - t386 * t678;
t669 = t197 * mrSges(5,2);
t668 = t222 * mrSges(3,2);
t628 = pkin(4) * t468 - t629;
t21 = -t132 * t488 - t142 * t487 - t383 * t89 + t387 * t80;
t666 = -t21 * mrSges(5,1) + t20 * mrSges(5,2);
t560 = cos(qJ(1));
t470 = t560 * t384;
t385 = sin(qJ(1));
t502 = t385 * t388;
t320 = t380 * t470 + t502;
t376 = pkin(11) + qJ(4);
t373 = sin(t376);
t374 = cos(t376);
t471 = t378 * t560;
t239 = t320 * t374 - t373 * t471;
t469 = t560 * t388;
t503 = t384 * t385;
t319 = -t380 * t469 + t503;
t183 = t239 * t382 - t319 * t386;
t665 = t239 * t386 + t319 * t382;
t663 = -mrSges(7,2) * t1 - mrSges(6,3) * t3 - t675;
t25 = qJ(6) * t196 + t27;
t75 = t387 * t132 - t383 * t142;
t72 = pkin(4) * t418 - t75;
t36 = t161 * pkin(5) - t162 * qJ(6) + t72;
t662 = t72 * mrSges(6,1) + t36 * mrSges(7,1) - t25 * mrSges(7,2) - t27 * mrSges(6,3);
t26 = -t382 * t73 + t386 * t85;
t632 = qJD(6) - t26;
t24 = -pkin(5) * t196 + t632;
t661 = -t197 * mrSges(5,1) - t26 * mrSges(6,1) + t24 * mrSges(7,1) + t27 * mrSges(6,2) - t25 * mrSges(7,3);
t660 = t141 * mrSges(5,2) - t21 * mrSges(5,3) + 0.2e1 * Ifges(5,1) * t590 + 0.2e1 * Ifges(5,4) * t589 + 0.2e1 * Ifges(5,5) * t569;
t658 = m(6) + m(7);
t656 = mrSges(6,3) + mrSges(7,2);
t28 = -mrSges(7,2) * t63 + mrSges(7,3) * t110;
t31 = -mrSges(6,2) * t110 - mrSges(6,3) * t63;
t648 = t28 + t31;
t29 = mrSges(6,1) * t110 - mrSges(6,3) * t62;
t30 = -t110 * mrSges(7,1) + t62 * mrSges(7,2);
t647 = -t29 + t30;
t646 = -qJ(6) * t625 + qJD(6) * t327 + t643;
t645 = pkin(5) * t625 - t644;
t215 = -t256 * t382 - t386 * t468;
t216 = -t256 * t386 + t382 * t468;
t419 = pkin(5) * t382 - qJ(6) * t386;
t420 = pkin(5) * t386 + qJ(6) * t382;
t642 = -pkin(5) * t215 + qJ(6) * t216 - t419 * t317 + (qJD(5) * t420 - qJD(6) * t386) * t328 + t628;
t637 = t418 * Ifges(5,3);
t636 = m(4) * qJ(3) + mrSges(4,3) + mrSges(5,3);
t526 = t407 * mrSges(5,3);
t167 = -mrSges(5,1) * t418 - t526;
t99 = mrSges(6,1) * t161 + mrSges(6,2) * t162;
t635 = t167 - t99;
t634 = -qJD(6) * t382 + t196 * t419 - t76;
t315 = -t377 * t508 + t379 * t380;
t316 = t377 * t380 + t379 * t508;
t218 = t315 * t383 + t316 * t387;
t558 = pkin(1) * t388;
t292 = t365 + (-pkin(2) - t558) * t380;
t224 = -t315 * pkin(3) + t292;
t406 = t387 * t315 - t316 * t383;
t121 = -pkin(4) * t406 - t218 * pkin(10) + t224;
t325 = pkin(8) * t506 + t384 * t559;
t289 = qJ(3) * t380 + t325;
t205 = -t377 * t289 + t379 * t290;
t148 = -pkin(3) * t506 - t316 * pkin(9) + t205;
t206 = t379 * t289 + t377 * t290;
t165 = pkin(9) * t315 + t206;
t95 = t383 * t148 + t387 * t165;
t91 = -pkin(10) * t506 + t95;
t633 = t382 * t121 + t386 * t91;
t449 = mrSges(3,3) * t468;
t631 = -mrSges(3,1) * t364 - mrSges(4,1) * t276 + mrSges(4,2) * t277 + t449;
t501 = t386 * t317;
t399 = t328 * t485 + t501;
t627 = t216 + t399;
t430 = mrSges(7,1) * t382 - mrSges(7,3) * t386;
t432 = mrSges(6,1) * t382 + mrSges(6,2) * t386;
t622 = t36 * t430 + t72 * t432;
t621 = -t382 * t652 + t386 * t653;
t533 = Ifges(7,5) * t382;
t538 = Ifges(6,4) * t382;
t620 = t386 * t655 + t533 - t538;
t436 = -mrSges(4,1) * t379 + mrSges(4,2) * t377;
t394 = m(4) * pkin(2) - t436;
t619 = -t374 * mrSges(5,1) + mrSges(5,2) * t373 - t394;
t400 = -t317 * t382 + t328 * t484;
t521 = t451 * t386;
t618 = t484 - t521;
t522 = t451 * t382;
t617 = -t485 + t522;
t616 = -t122 * t377 + t123 * t379;
t615 = t3 * t386 - t382 * t4;
t614 = t1 * t386 + t2 * t382;
t155 = Ifges(7,5) * t162;
t65 = t196 * Ifges(7,6) + t161 * Ifges(7,3) + t155;
t539 = Ifges(6,4) * t162;
t68 = -t161 * Ifges(6,2) + t196 * Ifges(6,6) + t539;
t613 = t68 / 0.2e1 - t65 / 0.2e1;
t544 = Ifges(3,4) * t384;
t612 = pkin(1) * (mrSges(3,1) * t384 + mrSges(3,2) * t388) - t384 * (Ifges(3,1) * t388 - t544) / 0.2e1;
t611 = mrSges(3,2) - t636;
t431 = -t386 * mrSges(7,1) - t382 * mrSges(7,3);
t433 = mrSges(6,1) * t386 - mrSges(6,2) * t382;
t610 = -m(7) * t420 - mrSges(5,1) + t431 - t433;
t609 = mrSges(3,1) - t619;
t261 = (qJD(2) * t421 - t489) * t378;
t310 = -pkin(8) * t467 + t388 * t481;
t270 = qJD(3) * t380 + t310;
t181 = t379 * t261 - t377 * t270;
t146 = qJD(2) * t393 + t181;
t182 = t377 * t261 + t379 * t270;
t441 = t377 * t378 * t490;
t163 = -pkin(9) * t441 + t182;
t42 = t383 * t146 + t148 * t487 + t387 * t163 - t165 * t488;
t40 = pkin(10) * t467 + t42;
t157 = -qJD(2) * t390 + qJD(4) * t406;
t158 = qJD(2) * t391 + qJD(4) * t218;
t311 = t325 * qJD(2);
t249 = pkin(3) * t441 + t311;
t81 = t158 * pkin(4) - t157 * pkin(10) + t249;
t9 = -qJD(5) * t633 - t382 * t40 + t386 * t81;
t450 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t444 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t531 = Ifges(3,6) * t364;
t607 = t76 * mrSges(5,2) + t531 / 0.2e1 + (t388 * Ifges(3,2) + t544) * t494 / 0.2e1 - t75 * mrSges(5,1);
t322 = -t380 * t503 + t469;
t507 = t378 * t385;
t243 = t322 * t374 + t373 * t507;
t288 = t373 * t380 + t374 * t508;
t606 = -g(1) * t243 - g(2) * t239 - g(3) * t288;
t604 = t378 ^ 2;
t595 = -t68 / 0.2e1;
t583 = -t162 / 0.2e1;
t580 = -t196 / 0.2e1;
t578 = -t451 / 0.2e1;
t576 = -t407 / 0.2e1;
t572 = t236 / 0.2e1;
t571 = t237 / 0.2e1;
t568 = t315 / 0.2e1;
t567 = t316 / 0.2e1;
t564 = t380 / 0.2e1;
t563 = t382 / 0.2e1;
t561 = t388 / 0.2e1;
t557 = pkin(4) * t374;
t546 = mrSges(6,3) * t161;
t545 = mrSges(6,3) * t162;
t543 = Ifges(3,4) * t388;
t542 = Ifges(4,4) * t377;
t541 = Ifges(4,4) * t379;
t537 = Ifges(6,4) * t386;
t536 = Ifges(4,5) * t277;
t535 = Ifges(5,5) * t407;
t532 = Ifges(7,5) * t386;
t530 = Ifges(4,6) * t276;
t529 = Ifges(5,6) * t451;
t528 = Ifges(4,3) * t384;
t527 = t451 * mrSges(5,3);
t525 = t364 * Ifges(3,5);
t129 = pkin(4) * t407 - pkin(10) * t451;
t45 = t382 * t129 + t386 * t75;
t518 = t319 * t373;
t321 = t380 * t502 + t470;
t515 = t321 * t373;
t513 = t328 * t386;
t511 = t374 * t382;
t510 = t374 * t386;
t509 = t377 * t388;
t505 = t379 * (Ifges(4,1) * t277 + Ifges(4,4) * t276 - Ifges(4,5) * t465);
t500 = t386 * t388;
t116 = -mrSges(7,2) * t161 + mrSges(7,3) * t196;
t117 = -mrSges(6,2) * t196 - t546;
t499 = t116 + t117;
t118 = mrSges(6,1) * t196 - t545;
t119 = -mrSges(7,1) * t196 + mrSges(7,2) * t162;
t498 = t118 - t119;
t497 = -t319 * t372 + t320 * t547;
t496 = -t321 * t372 + t322 * t547;
t495 = t560 * pkin(1) + pkin(8) * t507;
t367 = pkin(4) * t506;
t480 = pkin(10) * t485;
t479 = pkin(10) * t484;
t476 = t373 * t506;
t475 = t378 * t500;
t352 = t382 * t506;
t474 = t65 * t563;
t473 = Ifges(5,5) * t111 + Ifges(5,6) * t112 + Ifges(5,3) * t295;
t472 = Ifges(3,5) * t313 + Ifges(3,6) * t667 + Ifges(3,3) * t363;
t457 = -t485 / 0.2e1;
t456 = t484 / 0.2e1;
t455 = -pkin(1) * t385 + pkin(8) * t471;
t154 = -t236 * mrSges(4,1) + t237 * mrSges(4,2);
t52 = -t112 * mrSges(5,1) + t111 * mrSges(5,2);
t94 = t148 * t387 - t383 * t165;
t238 = -t320 * t373 - t374 * t471;
t448 = mrSges(3,3) * t465;
t443 = t377 * t471;
t90 = t367 - t94;
t435 = mrSges(4,1) * t377 + mrSges(4,2) * t379;
t429 = Ifges(4,1) * t379 - t542;
t426 = -Ifges(4,2) * t377 + t541;
t425 = -Ifges(6,2) * t382 + t537;
t422 = Ifges(7,3) * t382 + t532;
t412 = t377 * pkin(3) * t507 + t321 * t547 + t322 * t372 + t495;
t47 = t121 * t386 - t382 * t91;
t44 = t129 * t386 - t382 * t75;
t151 = t226 * t386 - t254 * t382;
t43 = t146 * t387 - t148 * t488 - t383 * t163 - t165 * t487;
t405 = -pkin(10) * t658 + mrSges(5,2) - t656;
t193 = t382 * t218 + t475;
t8 = t121 * t484 + t382 * t81 + t386 * t40 - t485 * t91;
t398 = pkin(3) * t443 - t319 * t547 - t320 * t372 + t455;
t397 = t250 * t435;
t19 = -pkin(4) * t295 - t21;
t41 = -pkin(4) * t467 - t43;
t357 = Ifges(3,4) * t465;
t338 = -pkin(4) - t420;
t330 = t372 * t506;
t324 = t380 * t558 - t365;
t323 = (-mrSges(3,1) * t388 + mrSges(3,2) * t384) * t378;
t306 = -t364 * mrSges(3,2) + t448;
t287 = -t373 * t508 + t374 * t380;
t258 = Ifges(3,1) * t468 + t357 + t525;
t242 = t322 * t373 - t374 * t507;
t234 = t288 * t382 + t475;
t233 = -mrSges(4,1) * t465 - t277 * mrSges(4,3);
t232 = mrSges(4,2) * t465 + t276 * mrSges(4,3);
t194 = t218 * t386 - t352;
t191 = -mrSges(4,1) * t667 - mrSges(4,3) * t237;
t190 = mrSges(4,2) * t667 + mrSges(4,3) * t236;
t188 = t243 * t386 + t321 * t382;
t187 = t243 * t382 - t321 * t386;
t178 = Ifges(4,4) * t277 + Ifges(4,2) * t276 - Ifges(4,6) * t465;
t177 = -Ifges(4,3) * t465 + t530 + t536;
t169 = t328 * t419 - t623;
t166 = mrSges(5,2) * t418 + t527;
t140 = -pkin(5) * t327 - t151;
t139 = qJ(6) * t327 + t626;
t135 = t237 * Ifges(4,1) + t236 * Ifges(4,4) - Ifges(4,5) * t667;
t134 = t237 * Ifges(4,4) + t236 * Ifges(4,2) - Ifges(4,6) * t667;
t128 = -mrSges(5,1) * t451 + mrSges(5,2) * t407;
t124 = t529 + t535 - t637;
t101 = -qJD(5) * t352 + t157 * t382 + t218 * t484 - t386 * t467;
t100 = -qJD(5) * t193 + t386 * t157 + t382 * t467;
t98 = mrSges(7,1) * t161 - mrSges(7,3) * t162;
t97 = pkin(5) * t162 + qJ(6) * t161;
t93 = -mrSges(5,2) * t295 + mrSges(5,3) * t112;
t92 = mrSges(5,1) * t295 - mrSges(5,3) * t111;
t51 = pkin(5) * t193 - qJ(6) * t194 + t90;
t38 = pkin(5) * t406 - t47;
t37 = -qJ(6) * t406 + t633;
t33 = -pkin(5) * t407 - t44;
t32 = qJ(6) * t407 + t45;
t23 = mrSges(6,1) * t63 + mrSges(6,2) * t62;
t22 = mrSges(7,1) * t63 - mrSges(7,3) * t62;
t10 = pkin(5) * t101 - qJ(6) * t100 - qJD(6) * t194 + t41;
t7 = -pkin(5) * t158 - t9;
t6 = qJ(6) * t158 - qJD(6) * t406 + t8;
t5 = pkin(5) * t63 - qJ(6) * t62 - qJD(6) * t162 + t19;
t11 = [-t122 * t316 * mrSges(4,3) - (Ifges(4,5) * t567 + Ifges(4,6) * t568 - Ifges(3,6) * t380 / 0.2e1 - t325 * mrSges(3,3) + (-pkin(1) * mrSges(3,1) - t544 + (-Ifges(3,2) - Ifges(4,3)) * t388) * t378) * t667 + t472 * t564 + t135 * t567 + t134 * t568 + (Ifges(4,1) * t316 + Ifges(4,4) * t315) * t571 + (Ifges(4,4) * t316 + Ifges(4,2) * t315) * t572 - t380 * t668 + (t640 / 0.2e1 + mrSges(6,2) * t72 + mrSges(7,2) * t24 - mrSges(6,3) * t26 - mrSges(7,3) * t36 + Ifges(6,4) * t585 + Ifges(7,5) * t584 + t653 * t579 + t655 * t582) * t100 + t660 * t218 + t123 * t315 * mrSges(4,3) + (-Ifges(5,6) * t589 - Ifges(5,5) * t590 - Ifges(5,3) * t569 - t237 * Ifges(4,5) - Ifges(4,6) * t236 - t473 / 0.2e1 - t122 * mrSges(4,1) + t222 * mrSges(3,3) + t123 * mrSges(4,2) + t666) * t506 + m(6) * (t19 * t90 + t26 * t9 + t27 * t8 + t3 * t633 + t4 * t47 + t41 * t72) + t633 * t31 + (mrSges(6,1) * t19 + mrSges(7,1) * t5 + t599 * t654 + t663 - t674) * t193 + (mrSges(6,2) * t19 + mrSges(7,2) * t2 - mrSges(6,3) * t4 - mrSges(7,3) * t5 + Ifges(6,4) * t598 + Ifges(7,5) * t597 + t673 + t676) * t194 + (-pkin(1) * t323 * t378 + Ifges(2,3)) * qJDD(1) + m(7) * (t1 * t37 + t10 * t36 + t2 * t38 + t24 * t7 + t25 * t6 + t5 * t51) + m(5) * (t141 * t224 + t197 * t249 + t20 * t95 + t21 * t94 + t42 * t76 + t43 * t75) + m(4) * (t122 * t205 + t123 * t206 + t171 * t181 + t172 * t182 + t203 * t292 + t250 * t311) - t671 * t406 + (-mrSges(5,3) * t76 - t661 + t670) * t158 + (Ifges(7,3) * t584 - Ifges(6,2) * t585 + t595 + t65 / 0.2e1 + t654 * t582 - t652 * t579 + t662) * t101 + t47 * t29 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t604 + t222 * t325 + t223 * t324 - t308 * t311 + t309 * t310) + (Ifges(3,3) * t564 + t324 * mrSges(3,1) - t325 * mrSges(3,2) + (Ifges(3,5) * t384 + Ifges(3,6) * t388) * t378) * t363 + (Ifges(3,5) * t564 - t324 * mrSges(3,3) + (-pkin(1) * mrSges(3,2) + t384 * Ifges(3,1) + t543) * t378) * t313 + t223 * (mrSges(3,1) * t380 - mrSges(3,3) * t508) + t203 * (-mrSges(4,1) * t315 + mrSges(4,2) * t316) + t310 * t306 + t292 * t154 + t38 * t30 + t37 * t28 + t249 * t128 + t182 * t232 + t181 * t233 + t224 * t52 + t205 * t191 + t206 * t190 + (-m(4) * (-pkin(2) * t320 + t455) - (-t320 * t379 + t443) * mrSges(4,1) - (t320 * t377 + t379 * t471) * mrSges(4,2) - m(3) * t455 + t320 * mrSges(3,1) - mrSges(3,3) * t471 + t385 * mrSges(2,1) + t560 * mrSges(2,2) - m(5) * t398 + t239 * mrSges(5,1) + t450 * t665 - t444 * t183 - t611 * t319 + t405 * t238 + t658 * (pkin(4) * t239 - t398)) * g(1) + (-t560 * mrSges(2,1) - m(5) * t412 - t243 * mrSges(5,1) + (-mrSges(3,1) - t394) * t322 + (mrSges(2,2) + (-mrSges(3,3) - t435) * t378) * t385 - t450 * t188 + t444 * t187 + t611 * t321 + t405 * t242 + (-m(4) - m(3)) * t495 - t658 * (t243 * pkin(4) + t412)) * g(2) + t51 * t22 + t90 * t23 + ((t505 / 0.2e1 - t377 * t178 / 0.2e1 - t308 * mrSges(3,3) + t525 / 0.2e1 + t277 * t429 / 0.2e1 + t276 * t426 / 0.2e1 + t397 + t258 / 0.2e1 + (-t171 * t379 - t172 * t377) * mrSges(4,3)) * t388 + (-t309 * mrSges(3,3) + qJD(4) * Ifges(5,3) / 0.2e1 - t172 * mrSges(4,2) - t531 / 0.2e1 + t171 * mrSges(4,1) + t536 / 0.2e1 + t530 / 0.2e1 + t535 / 0.2e1 + t529 / 0.2e1 + t177 / 0.2e1 + t124 / 0.2e1 - t607) * t384 + ((-Ifges(3,2) * t384 + t543) * t561 - (Ifges(4,5) * t504 - Ifges(4,6) * t509 + Ifges(5,3) * t384 + t528) * t388 / 0.2e1 - t612) * t494) * t491 + t94 * t92 + t95 * t93 + t10 * t98 + t41 * t99 + t6 * t116 + t8 * t117 + t631 * t311 + t9 * t118 + t7 * t119 + t42 * t166 + t43 * t167 + (-mrSges(5,3) * t75 + t669 - t672) * t157; (-Ifges(6,2) * t584 + Ifges(7,3) * t585 - t580 * t652 + t583 * t654 + t613 - t662) * t215 - t667 * (Ifges(4,5) * t377 + Ifges(4,6) * t379) / 0.2e1 + (Ifges(4,1) * t377 + t541) * t571 + (Ifges(4,2) * t379 + t542) * t572 + t256 * t587 + (-t306 + t448) * t308 + (-Ifges(6,4) * t399 + Ifges(7,5) * t216 - Ifges(6,2) * t400 + Ifges(7,6) * t255) * t585 + t626 * t31 + (t19 * t432 + t382 * t663 + t422 * t597 + t425 * t598 + t430 * t5 + t456 * t65 + t591 * t621 + t599 * t620 + t660) * t328 + ((t528 + (Ifges(4,5) * t379 - Ifges(4,6) * t377) * t388) * t561 + t612) * qJD(1) ^ 2 * t604 + t671 * t327 + (t190 * t379 - t191 * t377) * qJ(3) + (-t172 * (-mrSges(4,2) * t384 - mrSges(4,3) * t509) - t171 * (mrSges(4,1) * t384 - mrSges(4,3) * t504)) * t494 + (t232 * t379 - t233 * t377) * qJD(3) + t640 * (t328 * t457 - t501 / 0.2e1 - t216 / 0.2e1) + t203 * t436 + t670 * t318 + t513 * t676 + (Ifges(6,4) * t216 - Ifges(7,5) * t399 + Ifges(6,6) * t255 + Ifges(7,3) * t400) * t584 - t668 + (-m(5) * t330 + t323 - t656 * t476 - t658 * (pkin(10) * t476 + t374 * t367 + t508 * t547 + t330) + t444 * (t352 * t374 - t386 * t508) + (t619 * t388 + (-m(5) * t547 - t636) * t384 - t450 * (t374 * t500 + t382 * t384)) * t378) * g(3) + t418 * (-Ifges(5,5) * t256 - Ifges(5,6) * t255) / 0.2e1 + (-Ifges(5,4) * t256 - Ifges(5,2) * t255) * t578 + (-Ifges(5,1) * t256 - Ifges(5,4) * t255) * t576 + (-t399 * t653 - t400 * t652) * t579 + (-t399 * t655 + t400 * t654) * t582 - t397 * t465 + t379 * t134 / 0.2e1 + t377 * t135 / 0.2e1 - t372 * t52 + t254 * t93 + t255 * t125 / 0.2e1 - t248 * t128 - t210 * t232 - t209 * t233 + t223 * mrSges(3,1) + (-m(5) * t496 + t656 * t515 - t658 * (-pkin(10) * t515 - t321 * t557 + t496) - t450 * (-t321 * t510 + t322 * t382) + t444 * (-t321 * t511 - t322 * t386) + t611 * t322 + t609 * t321) * g(1) + (-m(5) * t497 + t656 * t518 - t658 * (-pkin(10) * t518 - t319 * t557 + t497) - t450 * (-t319 * t510 + t320 * t382) + t444 * (-t319 * t511 - t320 * t386) + t611 * t320 + t609 * t319) * g(2) + (t216 * t653 + t255 * t651) * t580 + (t216 * t655 + t255 * t653) * t583 + t642 * t98 + t643 * t117 + t644 * t118 + t645 * t119 + t646 * t116 + (t1 * t139 + t140 * t2 + t169 * t5 + t24 * t645 + t25 * t646 + t36 * t642) * m(7) - t641 * t255 / 0.2e1 - ((-Ifges(3,2) * t468 + t258 + t357 + t505) * t388 + t277 * (Ifges(4,5) * t384 + t388 * t429) + t276 * (Ifges(4,6) * t384 + t388 * t426) + t364 * (Ifges(3,5) * t388 - Ifges(3,6) * t384) + (t177 + t124) * t384) * t494 / 0.2e1 - (t23 - t92) * t623 + (t151 * t4 - t19 * t623 + t26 * t644 + t27 * t643 + t3 * t626 + t628 * t72) * m(6) + (-t141 * t372 - t197 * t248 + t20 * t254 + t21 * t623 + t629 * t75 + t630 * t76) * m(5) + t178 * t442 / 0.2e1 + t472 + t616 * mrSges(4,3) + (-t171 * t209 - t172 * t210 - pkin(2) * t203 + (-t171 * t377 + t172 * t379) * qJD(3) + t616 * qJ(3)) * m(4) + t400 * t595 + (t24 * t625 + t36 * t400) * mrSges(7,1) + (t624 * t75 + t625 * t76) * mrSges(5,3) + (-t26 * t625 + t400 * t72) * mrSges(6,1) + (-mrSges(5,1) * t625 - mrSges(5,2) * t624) * t197 + (t26 * t627 - t27 * t400 - t4 * t513) * mrSges(6,3) + (-t25 * t625 + t36 * t627) * mrSges(7,3) + (t2 * t513 - t24 * t627 - t25 * t400) * mrSges(7,2) + (t27 * t625 - t627 * t72) * mrSges(6,2) + t628 * t99 + t629 * t167 + t630 * t166 + (-m(4) * t250 + t449 - t631) * t309 + t139 * t28 + t140 * t30 + t151 * t29 - pkin(2) * t154 + (Ifges(5,5) * t576 + Ifges(5,6) * t578 + t637 / 0.2e1 + t607) * t468 + t169 * t22 + (-t474 + t672) * t317; -t451 * t166 - t276 * t232 + t277 * t233 + (-t98 + t635) * t407 + (t196 * t499 - t647) * t386 + (-t196 * t498 + t648) * t382 + t52 + t154 + (t1 * t382 - t407 * t36 - t2 * t386 + t196 * (t24 * t382 + t25 * t386)) * m(7) + (-t407 * t72 + t3 * t382 + t4 * t386 + t196 * (-t26 * t382 + t27 * t386)) * m(6) + (t407 * t75 - t451 * t76 + t141) * m(5) + (t171 * t277 - t172 * t276 + t203) * m(4) + (-g(1) * t321 - g(2) * t319 + g(3) * t506) * (m(4) + m(5) + t658); (-t480 - t32) * t116 + (-t480 - t45) * t117 + t5 * t431 - t19 * t433 + t533 * t597 + t538 * t598 + t125 * t575 + (-t425 / 0.2e1 + t422 / 0.2e1) * t486 + (t425 * t584 + t422 * t585 + Ifges(5,1) * t576 + t639 / 0.2e1 - t669 + t620 * t583 + t621 * t580 - t622) * t451 + (Ifges(6,6) * t584 + Ifges(7,6) * t585 - Ifges(5,2) * t578 - t638 / 0.2e1 + t653 * t583 + t651 * t580 + t661) * t407 + (t479 - t33) * t119 + t673 * t382 + t68 * t457 + (t162 * t620 + t196 * t621) * qJD(5) / 0.2e1 + (t527 - t166) * t75 + (t674 + t675) * t386 + t647 * pkin(10) * t382 + t648 * pkin(10) * t386 + (t456 - t521 / 0.2e1) * t640 + (-t532 + t537) * t599 - t666 + (-t479 - t44) * t118 + (t126 + t195) * t578 + t338 * t22 - pkin(4) * t23 + (mrSges(5,2) * t243 - t658 * (-t242 * pkin(4) + pkin(10) * t243) - t610 * t242) * g(1) + (mrSges(5,2) * t288 - t658 * (t287 * pkin(4) + pkin(10) * t288) + t610 * t287) * g(3) + (mrSges(5,2) * t239 - t658 * (t238 * pkin(4) + pkin(10) * t239) + t610 * t238) * g(2) + t649 * t563 + (-t540 + t641) * t576 + t473 + t613 * t522 + (-pkin(4) * t19 + ((-t26 * t386 - t27 * t382) * qJD(5) + t615) * pkin(10) - t26 * t44 - t27 * t45) * m(6) + (t24 * t618 + t25 * t617 + t606 + t614) * mrSges(7,2) + (-t26 * t618 + t27 * t617 + t606 + t615) * mrSges(6,3) + (t474 + t622) * qJD(5) + t634 * t98 + (-t24 * t33 - t25 * t32 + t338 * t5 + ((t24 * t386 - t25 * t382) * qJD(5) + t614) * pkin(10) + t634 * t36) * m(7) + (-m(6) * t72 + t526 + t635) * t76; (Ifges(7,3) * t162 - t534) * t585 + t68 * t582 + t605 + (t498 + t545) * t27 + (-t499 - t546) * t26 + (-t161 * t655 + t155 - t539 + t65) * t583 + (-pkin(5) * t2 + qJ(6) * t1 - t24 * t27 + t25 * t632 - t36 * t97) * m(7) + (-t161 * t653 - t162 * t652) * t580 + (t450 * t183 + t444 * t665) * g(2) + (-Ifges(6,2) * t162 - t156 + t640) * t584 + (t187 * t450 + t188 * t444) * g(1) + (t444 * (t288 * t386 - t352) + t450 * t234) * g(3) + (t161 * t24 + t162 * t25) * mrSges(7,2) + qJ(6) * t28 - pkin(5) * t30 - t97 * t98 + qJD(6) * t116 - t36 * (mrSges(7,1) * t162 + mrSges(7,3) * t161) - t72 * (mrSges(6,1) * t162 - mrSges(6,2) * t161) + t650; -t196 * t116 + t162 * t98 + (-g(1) * t187 - g(2) * t183 - g(3) * t234 + t162 * t36 - t196 * t25 + t2) * m(7) + t30;];
tau  = t11;
