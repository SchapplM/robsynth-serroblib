% Calculate vector of inverse dynamics joint torques for
% S6RRRRRP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:55:40
% EndTime: 2019-03-10 00:56:28
% DurationCPUTime: 27.71s
% Computational Cost: add. (22615->780), mult. (51739->986), div. (0->0), fcn. (38335->14), ass. (0->384)
t688 = Ifges(6,4) + Ifges(7,4);
t381 = cos(qJ(5));
t651 = -mrSges(7,1) - mrSges(6,1);
t693 = t381 * t651;
t649 = Ifges(6,1) + Ifges(7,1);
t647 = Ifges(6,5) + Ifges(7,5);
t678 = Ifges(6,2) + Ifges(7,2);
t645 = Ifges(6,6) + Ifges(7,6);
t372 = qJD(2) + qJD(3);
t364 = qJD(4) + t372;
t376 = sin(qJ(5));
t378 = sin(qJ(3));
t379 = sin(qJ(2));
t383 = cos(qJ(3));
t384 = cos(qJ(2));
t307 = -t378 * t379 + t383 * t384;
t286 = t307 * qJD(1);
t308 = t378 * t384 + t379 * t383;
t287 = t308 * qJD(1);
t377 = sin(qJ(4));
t382 = cos(qJ(4));
t416 = t286 * t377 + t382 * t287;
t210 = t364 * t381 - t376 * t416;
t692 = t688 * t210;
t488 = qJD(5) * t376;
t445 = t382 * t286 - t287 * t377;
t631 = t445 * t376;
t691 = t488 - t631;
t211 = t364 * t376 + t381 * t416;
t690 = t688 * t211;
t650 = mrSges(6,2) + mrSges(7,2);
t689 = -mrSges(6,3) - mrSges(7,3);
t644 = Ifges(6,3) + Ifges(7,3);
t230 = qJD(5) - t445;
t670 = t210 * t678 + t230 * t645 + t690;
t669 = t211 * t649 + t230 * t647 + t692;
t374 = qJ(2) + qJ(3);
t369 = qJ(4) + t374;
t350 = sin(t369);
t687 = t350 * t693;
t487 = qJD(5) * t381;
t686 = t445 * t381 - t487;
t684 = t350 * t650;
t666 = t691 * pkin(5);
t665 = qJ(6) * t631 + t381 * qJD(6);
t351 = cos(t369);
t565 = pkin(5) * t381;
t353 = pkin(4) + t565;
t375 = -qJ(6) - pkin(10);
t413 = -t350 * t353 - t351 * t375;
t683 = -m(7) * t413 - t687;
t682 = t688 * t381;
t681 = t688 * t376;
t386 = -pkin(8) - pkin(7);
t338 = t386 * t384;
t315 = qJD(1) * t338;
t288 = t378 * t315;
t337 = t386 * t379;
t314 = qJD(1) * t337;
t295 = qJD(2) * pkin(2) + t314;
t240 = t383 * t295 + t288;
t280 = t287 * pkin(9);
t208 = t240 - t280;
t196 = pkin(3) * t372 + t208;
t291 = t383 * t315;
t241 = t295 * t378 - t291;
t562 = pkin(9) * t286;
t209 = t241 + t562;
t198 = t377 * t209;
t129 = t196 * t382 - t198;
t124 = -pkin(4) * t364 - t129;
t628 = mrSges(5,1) * t364 + mrSges(6,1) * t210 - mrSges(6,2) * t211 - mrSges(5,3) * t416;
t680 = -m(5) * t129 + m(6) * t124 - t628;
t486 = qJD(1) * qJD(2);
t319 = qJDD(1) * t384 - t379 * t486;
t320 = qJDD(1) * t379 + t384 * t486;
t400 = t307 * qJD(3);
t214 = qJD(1) * t400 + t319 * t378 + t320 * t383;
t401 = t308 * qJD(3);
t215 = -qJD(1) * t401 + t319 * t383 - t320 * t378;
t119 = -qJD(4) * t416 - t214 * t377 + t215 * t382;
t117 = qJDD(5) - t119;
t302 = t320 * pkin(7);
t254 = qJDD(2) * pkin(2) - pkin(8) * t320 - t302;
t301 = t319 * pkin(7);
t255 = pkin(8) * t319 + t301;
t491 = qJD(3) * t383;
t492 = qJD(3) * t378;
t150 = t378 * t254 + t383 * t255 + t295 * t491 + t315 * t492;
t105 = pkin(9) * t215 + t150;
t489 = qJD(4) * t382;
t490 = qJD(4) * t377;
t151 = -qJD(3) * t241 + t383 * t254 - t255 * t378;
t371 = qJDD(2) + qJDD(3);
t94 = pkin(3) * t371 - pkin(9) * t214 + t151;
t26 = t382 * t105 + t196 * t489 - t209 * t490 + t377 * t94;
t363 = qJDD(4) + t371;
t23 = pkin(10) * t363 + t26;
t118 = qJD(4) * t445 + t214 * t382 + t215 * t377;
t532 = qJDD(1) * pkin(1);
t276 = -pkin(2) * t319 - t532;
t183 = -pkin(3) * t215 + t276;
t37 = -pkin(4) * t119 - pkin(10) * t118 + t183;
t199 = t382 * t209;
t130 = t196 * t377 + t199;
t125 = pkin(10) * t364 + t130;
t370 = t384 * pkin(2);
t356 = t370 + pkin(1);
t336 = t356 * qJD(1);
t256 = -pkin(3) * t286 - t336;
t139 = -pkin(4) * t445 - pkin(10) * t416 + t256;
t55 = t125 * t381 + t139 * t376;
t6 = -qJD(5) * t55 - t23 * t376 + t381 * t37;
t70 = qJD(5) * t210 + t118 * t381 + t363 * t376;
t1 = pkin(5) * t117 - qJ(6) * t70 - qJD(6) * t211 + t6;
t27 = -t377 * t105 - t196 * t490 - t209 * t489 + t382 * t94;
t24 = -pkin(4) * t363 - t27;
t71 = -qJD(5) * t211 - t118 * t376 + t363 * t381;
t12 = -pkin(5) * t71 + qJDD(6) + t24;
t544 = Ifges(5,4) * t416;
t640 = t445 * Ifges(5,2);
t642 = t364 * Ifges(5,6);
t166 = t544 + t640 + t642;
t226 = Ifges(5,4) * t445;
t643 = t364 * Ifges(5,5);
t167 = Ifges(5,1) * t416 + t226 + t643;
t5 = -t125 * t488 + t139 * t487 + t381 * t23 + t376 * t37;
t3 = qJ(6) * t71 + qJD(6) * t210 + t5;
t54 = -t125 * t376 + t381 * t139;
t41 = -qJ(6) * t211 + t54;
t38 = pkin(5) * t230 + t41;
t42 = qJ(6) * t210 + t55;
t452 = -t488 / 0.2e1;
t577 = t376 / 0.2e1;
t578 = -t364 / 0.2e1;
t582 = t416 / 0.2e1;
t584 = -t445 / 0.2e1;
t587 = -t230 / 0.2e1;
t589 = -t211 / 0.2e1;
t591 = -t210 / 0.2e1;
t593 = t117 / 0.2e1;
t594 = t71 / 0.2e1;
t595 = t70 / 0.2e1;
t598 = -t256 * mrSges(5,1) - t54 * mrSges(6,1) - t38 * mrSges(7,1) + t55 * mrSges(6,2) + t42 * mrSges(7,2) + t130 * mrSges(5,3);
t605 = t256 * mrSges(5,2) - t129 * mrSges(5,3);
t616 = t381 * t649 - t681;
t617 = -t376 * t678 + t682;
t618 = -t376 * t645 + t381 * t647;
t629 = t210 * t645 + t211 * t647 + t230 * t644;
t652 = -t416 / 0.2e1;
t660 = -t376 * t6 + t381 * t5;
t553 = mrSges(7,2) * t381;
t427 = mrSges(7,1) * t376 + t553;
t428 = mrSges(6,1) * t376 + mrSges(6,2) * t381;
t87 = -pkin(5) * t210 + qJD(6) + t124;
t664 = t124 * t428 + t87 * t427;
t676 = t117 * t647 + t649 * t70 + t688 * t71;
t677 = t117 * t645 + t678 * t71 + t688 * t70;
t679 = t27 * mrSges(5,1) - t26 * mrSges(5,2) + t24 * (-mrSges(6,1) * t381 + mrSges(6,2) * t376) + t12 * (-mrSges(7,1) * t381 + mrSges(7,2) * t376) + Ifges(5,3) * t363 + Ifges(5,5) * t118 + Ifges(5,6) * t119 + (t376 * t649 + t682) * t595 + (t381 * t678 + t681) * t594 + (t376 * t647 + t381 * t645) * t593 + t676 * t577 + t677 * t381 / 0.2e1 + t670 * t452 + t669 * t487 / 0.2e1 + t664 * qJD(5) + (t210 * t617 + t211 * t616 + t230 * t618) * qJD(5) / 0.2e1 + t166 * t582 + (t226 + t167) * t584 + (-t1 * t376 + t3 * t381 + t38 * t686 - t42 * t691) * mrSges(7,3) + (t54 * t686 - t55 * t691 + t660) * mrSges(6,3) + (-t544 + t629) * t652 + (-t605 + Ifges(5,1) * t652 + Ifges(5,5) * t578 - t664 + t670 * t577 - t669 * t381 / 0.2e1 + t591 * t617 + t589 * t616 + t587 * t618) * t445 + (-Ifges(5,2) * t584 - Ifges(5,6) * t578 + t587 * t644 + t589 * t647 + t591 * t645 + t598) * t416;
t478 = pkin(3) * t489;
t134 = t208 * t382 - t198;
t180 = pkin(4) * t416 - pkin(10) * t445;
t571 = pkin(3) * t287;
t152 = t180 + t571;
t62 = t381 * t134 + t376 * t152;
t675 = t381 * t478 - t62;
t572 = pkin(2) * t383;
t355 = pkin(3) + t572;
t507 = t377 * t378;
t236 = t355 * t489 + (-t378 * t490 + (t382 * t383 - t507) * qJD(3)) * pkin(2);
t247 = -t314 * t378 + t291;
t216 = t247 - t562;
t248 = t383 * t314 + t288;
t217 = -t280 + t248;
t143 = t216 * t377 + t217 * t382;
t495 = qJD(1) * t379;
t359 = pkin(2) * t495;
t144 = t152 + t359;
t64 = t381 * t143 + t376 * t144;
t674 = t236 * t381 - t64;
t133 = t208 * t377 + t199;
t668 = -pkin(3) * t490 + t133;
t506 = t378 * t382;
t667 = t382 * t216 - t217 * t377 + t355 * t490 + (t378 * t489 + (t377 * t383 + t506) * qJD(3)) * pkin(2);
t366 = sin(t374);
t367 = cos(t374);
t429 = mrSges(5,1) * t350 + mrSges(5,2) * t351;
t663 = mrSges(4,1) * t366 + mrSges(4,2) * t367 + t429;
t385 = cos(qJ(1));
t508 = t376 * t385;
t511 = t351 * t385;
t662 = -t508 * t684 + t511 * t689;
t380 = sin(qJ(1));
t509 = t376 * t380;
t513 = t351 * t380;
t661 = -t509 * t684 + t513 * t689;
t610 = g(1) * t385 + g(2) * t380;
t659 = -t351 * mrSges(5,1) + (mrSges(5,2) + t689) * t350;
t620 = t351 * pkin(4) + t350 * pkin(10);
t623 = -t350 * t375 + t351 * t353;
t657 = -m(6) * t620 - m(7) * t623;
t596 = m(7) * pkin(5);
t653 = t319 / 0.2e1;
t574 = t384 / 0.2e1;
t641 = t384 * Ifges(3,2);
t282 = pkin(2) * t506 + t377 * t355;
t275 = pkin(10) + t282;
t503 = -qJ(6) - t275;
t442 = qJD(5) * t503;
t639 = t376 * t442 + t665 + t674;
t368 = t381 * qJ(6);
t434 = pkin(5) * t416 - t368 * t445;
t63 = -t143 * t376 + t381 * t144;
t638 = (-qJD(6) - t236) * t376 + t381 * t442 - t434 - t63;
t569 = pkin(3) * t377;
t352 = pkin(10) + t569;
t502 = -qJ(6) - t352;
t441 = qJD(5) * t502;
t637 = t376 * t441 + t665 + t675;
t61 = -t134 * t376 + t381 * t152;
t636 = (-qJD(6) - t478) * t376 + t381 * t441 - t434 - t61;
t449 = qJD(5) * t375;
t73 = t381 * t129 + t376 * t180;
t635 = t376 * t449 + t665 - t73;
t72 = -t129 * t376 + t381 * t180;
t634 = -qJD(6) * t376 + t381 * t449 - t434 - t72;
t633 = -t130 + t666;
t632 = t596 + mrSges(7,1);
t259 = t383 * t337 + t338 * t378;
t227 = -pkin(9) * t308 + t259;
t260 = t378 * t337 - t383 * t338;
t228 = pkin(9) * t307 + t260;
t169 = t227 * t377 + t228 * t382;
t163 = t381 * t169;
t246 = t307 * t377 + t308 * t382;
t265 = -pkin(3) * t307 - t356;
t415 = t382 * t307 - t308 * t377;
t164 = -pkin(4) * t415 - pkin(10) * t246 + t265;
t82 = t376 * t164 + t163;
t627 = t382 * t227 - t228 * t377;
t626 = t666 + t667;
t625 = t236 - t143;
t622 = t666 - t668;
t447 = t367 * mrSges(4,1) - mrSges(4,2) * t366;
t567 = pkin(4) * t350;
t570 = pkin(3) * t366;
t619 = -m(7) * (t413 - t570) - m(6) * (-t567 - t570) - t687;
t613 = t117 * t644 + t645 * t71 + t647 * t70;
t494 = qJD(1) * t384;
t563 = pkin(7) * t384;
t564 = pkin(7) * t379;
t612 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t495) * t563 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t494) * t564;
t611 = t301 * t384 + t302 * t379;
t609 = m(7) + m(6) + m(5);
t608 = mrSges(6,1) + t632;
t606 = t659 + (t376 * t650 + t693) * t351;
t604 = -m(3) * pkin(7) + m(4) * t386 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t603 = -t447 + t606;
t551 = mrSges(6,3) * t210;
t154 = -mrSges(6,2) * t230 + t551;
t550 = mrSges(6,3) * t211;
t156 = mrSges(6,1) * t230 - t550;
t34 = mrSges(6,1) * t117 - mrSges(6,3) * t70;
t602 = m(6) * ((-t376 * t55 - t381 * t54) * qJD(5) + t660) - t156 * t487 - t154 * t488 - t34 * t376;
t334 = -mrSges(3,1) * t384 + mrSges(3,2) * t379;
t601 = m(3) * pkin(1) + m(4) * t356 + mrSges(2,1) - t334 + t447 - t659;
t600 = t6 * mrSges(6,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2);
t599 = (t385 * t683 + t662) * g(1) + (t380 * t683 + t661) * g(2);
t590 = t210 / 0.2e1;
t588 = t211 / 0.2e1;
t586 = t230 / 0.2e1;
t579 = t287 / 0.2e1;
t573 = pkin(2) * t379;
t349 = pkin(3) * t367;
t568 = pkin(3) * t382;
t561 = pkin(10) * t381;
t552 = mrSges(4,3) * t287;
t549 = mrSges(7,3) * t210;
t548 = mrSges(7,3) * t211;
t547 = Ifges(3,4) * t379;
t546 = Ifges(3,4) * t384;
t545 = Ifges(4,4) * t287;
t535 = t240 * mrSges(4,3);
t251 = qJD(2) * t307 + t400;
t252 = -qJD(2) * t308 - t401;
t145 = qJD(4) * t415 + t251 * t382 + t252 * t377;
t531 = t145 * t376;
t530 = t145 * t381;
t529 = t156 * t376;
t522 = t246 * t376;
t521 = t246 * t381;
t520 = t275 * t381;
t510 = t352 * t381;
t505 = t380 * t381;
t504 = t381 * t385;
t496 = t349 + t370;
t493 = qJD(2) * t379;
t470 = qJD(2) * t386;
t317 = t379 * t470;
t318 = t384 * t470;
t190 = t383 * t317 + t378 * t318 + t337 * t491 + t338 * t492;
t161 = pkin(9) * t252 + t190;
t191 = -qJD(3) * t260 - t317 * t378 + t383 * t318;
t162 = -pkin(9) * t251 + t191;
t45 = qJD(4) * t627 + t161 * t382 + t162 * t377;
t146 = qJD(4) * t246 + t251 * t377 - t382 * t252;
t361 = pkin(2) * t493;
t235 = -pkin(3) * t252 + t361;
t52 = pkin(4) * t146 - pkin(10) * t145 + t235;
t482 = t164 * t487 + t376 * t52 + t381 * t45;
t471 = t349 + t620;
t467 = t246 * t487;
t28 = -t71 * mrSges(7,1) + t70 * mrSges(7,2);
t450 = -t376 * t45 + t381 * t52;
t448 = t486 / 0.2e1;
t81 = t381 * t164 - t169 * t376;
t281 = -pkin(2) * t507 + t355 * t382;
t329 = pkin(10) * t513;
t437 = -t380 * t567 + t329;
t330 = pkin(10) * t511;
t436 = -t385 * t567 + t330;
t274 = -pkin(4) - t281;
t433 = t349 + t623;
t432 = mrSges(3,1) * t379 + mrSges(3,2) * t384;
t424 = t547 + t641;
t421 = Ifges(3,5) * t384 - Ifges(3,6) * t379;
t418 = -t376 * t54 + t381 * t55;
t412 = -qJ(6) * t145 - qJD(6) * t246;
t409 = pkin(1) * t432;
t271 = -t351 * t508 + t505;
t269 = t351 * t509 + t504;
t407 = t467 + t531;
t406 = t246 * t488 - t530;
t404 = t379 * (Ifges(3,1) * t384 - t547);
t46 = qJD(4) * t169 + t161 * t377 - t382 * t162;
t224 = Ifges(4,2) * t286 + Ifges(4,6) * t372 + t545;
t277 = Ifges(4,4) * t286;
t225 = t287 * Ifges(4,1) + t372 * Ifges(4,5) + t277;
t387 = -(-Ifges(4,2) * t287 + t225 + t277) * t286 / 0.2e1 + t336 * (mrSges(4,1) * t287 + mrSges(4,2) * t286) - t287 * (Ifges(4,1) * t286 - t545) / 0.2e1 + Ifges(4,3) * t371 - t372 * (Ifges(4,5) * t286 - Ifges(4,6) * t287) / 0.2e1 + Ifges(4,6) * t215 + Ifges(4,5) * t214 - t150 * mrSges(4,2) + t151 * mrSges(4,1) + t286 * t535 + t241 * t552 + t224 * t579 + t679;
t373 = -pkin(9) + t386;
t358 = Ifges(3,4) * t494;
t354 = -pkin(4) - t568;
t335 = t368 + t561;
t333 = t375 * t376;
t328 = -t353 - t568;
t321 = -t570 - t573;
t313 = pkin(1) + t496;
t299 = t385 * t321;
t298 = t380 * t321;
t297 = t368 + t510;
t296 = t502 * t376;
t285 = Ifges(3,1) * t495 + Ifges(3,5) * qJD(2) + t358;
t284 = Ifges(3,6) * qJD(2) + qJD(1) * t424;
t272 = t351 * t504 + t509;
t270 = -t351 * t505 + t508;
t264 = t274 - t565;
t263 = mrSges(4,1) * t372 - t552;
t262 = -mrSges(4,2) * t372 + mrSges(4,3) * t286;
t261 = t359 + t571;
t258 = t368 + t520;
t257 = t503 * t376;
t239 = -mrSges(4,1) * t286 + mrSges(4,2) * t287;
t218 = -mrSges(5,2) * t364 + mrSges(5,3) * t445;
t195 = -mrSges(4,2) * t371 + mrSges(4,3) * t215;
t194 = mrSges(4,1) * t371 - mrSges(4,3) * t214;
t179 = -mrSges(5,1) * t445 + mrSges(5,2) * t416;
t155 = mrSges(7,1) * t230 - t548;
t153 = -mrSges(7,2) * t230 + t549;
t140 = -mrSges(7,1) * t210 + mrSges(7,2) * t211;
t123 = pkin(5) * t522 - t627;
t97 = -mrSges(5,2) * t363 + mrSges(5,3) * t119;
t96 = mrSges(5,1) * t363 - mrSges(5,3) * t118;
t60 = -qJ(6) * t522 + t82;
t49 = -pkin(5) * t415 - t246 * t368 + t81;
t36 = -mrSges(6,2) * t117 + mrSges(6,3) * t71;
t35 = -mrSges(7,2) * t117 + mrSges(7,3) * t71;
t33 = mrSges(7,1) * t117 - mrSges(7,3) * t70;
t30 = pkin(5) * t407 + t46;
t29 = -mrSges(6,1) * t71 + mrSges(6,2) * t70;
t10 = -qJD(5) * t82 + t450;
t9 = -t169 * t488 + t482;
t8 = -qJ(6) * t467 + (-qJD(5) * t169 + t412) * t376 + t482;
t7 = pkin(5) * t146 + t412 * t381 + (-t163 + (qJ(6) * t246 - t164) * t376) * qJD(5) + t450;
t2 = [(-mrSges(4,1) * t276 + mrSges(4,3) * t150 + Ifges(4,4) * t214 + Ifges(4,2) * t215 + Ifges(4,6) * t371) * t307 + (-t406 * t688 - t407 * t678) * t590 + (-t406 * t649 - t407 * t688) * t588 + (t643 / 0.2e1 + t226 / 0.2e1 + t167 / 0.2e1 + Ifges(5,1) * t582 + t605) * t145 + (-t613 / 0.2e1 - mrSges(5,1) * t183 - t1 * mrSges(7,1) + mrSges(5,3) * t26 + Ifges(5,4) * t118 + Ifges(5,2) * t119 + Ifges(5,6) * t363 - t593 * t644 - t594 * t645 - t595 * t647 - t600) * t415 + (t651 * t270 - t650 * t269 + (t373 * t609 - t376 * t596 + t604) * t385 + (-m(7) * (-t313 - t623) - m(6) * (-t313 - t620) + m(5) * t313 + t601) * t380) * g(1) + m(4) * (t150 * t260 + t151 * t259 + t190 * t241 + t191 * t240 - t276 * t356 - t336 * t361) + (t629 / 0.2e1 - t642 / 0.2e1 - t640 / 0.2e1 - t166 / 0.2e1 - Ifges(5,4) * t582 + t645 * t590 + t647 * t588 + t644 * t586 - t598) * t146 + (Ifges(3,4) * t320 + Ifges(3,2) * t319) * t574 + m(6) * (t10 * t54 + t5 * t82 + t55 * t9 + t6 * t81) + m(5) * (t130 * t45 + t169 * t26 + t183 * t265 + t235 * t256) + (t384 * t546 + t404) * t448 + t320 * t546 / 0.2e1 + (t183 * mrSges(5,2) - t27 * mrSges(5,3) + Ifges(5,1) * t118 + Ifges(5,4) * t119 + Ifges(5,5) * t363 + t12 * t427 + t24 * t428 + t593 * t618 + t594 * t617 + t595 * t616) * t246 - t409 * t486 + (-t1 * t521 - t3 * t522 + t38 * t406 - t407 * t42) * mrSges(7,3) + t680 * t46 - (-m(5) * t27 + m(6) * t24 + t29 - t96) * t627 + (mrSges(4,2) * t276 - mrSges(4,3) * t151 + Ifges(4,1) * t214 + Ifges(4,4) * t215 + Ifges(4,5) * t371) * t308 + m(7) * (t1 * t49 + t12 * t123 + t3 * t60 + t30 * t87 + t38 * t7 + t42 * t8) + (t406 * t54 - t407 * t55 - t5 * t522 - t521 * t6) * mrSges(6,3) - t251 * t535 + (t319 * t563 + t320 * t564 + t611) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t611) + (t285 * t574 + t421 * qJD(2) / 0.2e1 - t612) * qJD(2) - t284 * t493 / 0.2e1 + t676 * t521 / 0.2e1 - t677 * t522 / 0.2e1 + Ifges(2,3) * qJDD(1) + t669 * (t246 * t452 + t530 / 0.2e1) + t670 * (-t467 / 0.2e1 - t531 / 0.2e1) + t87 * (mrSges(7,1) * t407 - mrSges(7,2) * t406) + t124 * (mrSges(6,1) * t407 - mrSges(6,2) * t406) + (-t509 * t596 - t609 * (t385 * t313 - t373 * t380) + t651 * t272 - t650 * t271 + t604 * t380 + (-t601 + t657) * t385) * g(2) + t372 * (Ifges(4,5) * t251 + Ifges(4,6) * t252) / 0.2e1 - t356 * (-mrSges(4,1) * t215 + mrSges(4,2) * t214) - t336 * (-mrSges(4,1) * t252 + mrSges(4,2) * t251) - pkin(1) * (-mrSges(3,1) * t319 + mrSges(3,2) * t320) + t286 * (Ifges(4,4) * t251 + Ifges(4,2) * t252) / 0.2e1 + t265 * (-mrSges(5,1) * t119 + mrSges(5,2) * t118) + t259 * t194 + t260 * t195 + t190 * t262 + t191 * t263 + t251 * t225 / 0.2e1 + t252 * t224 / 0.2e1 + t235 * t179 + t45 * t218 + t169 * t97 + t8 * t153 + t9 * t154 + t7 * t155 + t10 * t156 + (-t406 * t647 - t407 * t645) * t586 + t49 * t33 - t334 * t532 + t239 * t361 + t60 * t35 + (-mrSges(3,1) * t564 - mrSges(3,2) * t563 + 0.2e1 * Ifges(3,6) * t574) * qJDD(2) + (Ifges(3,1) * t320 + Ifges(3,4) * t653 + Ifges(3,5) * qJDD(2) - t448 * t641) * t379 + t81 * t34 + t82 * t36 + (Ifges(4,1) * t251 + Ifges(4,4) * t252) * t579 + t123 * t28 + t30 * t140 + t424 * t653 + t241 * t252 * mrSges(4,3); t599 - m(4) * (t240 * t247 + t241 * t248 - t336 * t359) + (-m(6) * (t370 + t471) - m(4) * t370 - m(7) * (t370 + t433) + t334 - m(5) * t496 + t603) * g(3) - (-Ifges(3,2) * t495 + t285 + t358) * t494 / 0.2e1 + t638 * t155 + (-g(1) * t299 - g(2) * t298 + t1 * t257 + t12 * t264 + t258 * t3 + t38 * t638 + t42 * t639 + t626 * t87) * m(7) + t639 * t153 + t602 * t275 - t421 * t486 / 0.2e1 + t680 * t667 - t236 * t529 + (t236 * t418 + t24 * t274 - t54 * t63 - t55 * t64 - g(1) * (t299 + t436) - g(2) * (t298 + t437)) * m(6) + (t130 * t625 - t256 * t261 + t26 * t282 + t27 * t281) * m(5) + t387 + t625 * t218 + t626 * t140 - t239 * t359 + t674 * t154 + t610 * (m(4) * t573 - m(5) * t321 + t432 + t663) + Ifges(3,3) * qJDD(2) + Ifges(3,6) * t319 + Ifges(3,5) * t320 - t301 * mrSges(3,2) - t302 * mrSges(3,1) + t281 * t96 + t282 * t97 + t274 * t29 + t257 * t33 + t258 * t35 - t261 * t179 - t248 * t262 - t247 * t263 + t264 * t28 + (t612 + (t409 - t404 / 0.2e1) * qJD(1)) * qJD(1) - t63 * t156 + (m(4) * (t150 * t378 + t151 * t383 + (-t240 * t378 + t241 * t383) * qJD(3)) + t262 * t491 + t378 * t195 - t263 * t492) * pkin(2) + t36 * t520 + t194 * t572 + t284 * t495 / 0.2e1; (t24 * t354 + (t124 * t377 + t382 * t418) * qJD(4) * pkin(3) - t124 * t133 - t54 * t61 - t55 * t62 - g(1) * t330 - g(2) * t329) * m(6) + ((t26 * t377 + t27 * t382 + (-t129 * t377 + t130 * t382) * qJD(4)) * pkin(3) + t129 * t133 - t130 * t134 - t256 * t571) * m(5) + (-t134 + t478) * t218 + (-m(5) * t349 - m(6) * t471 - m(7) * t433 + t603) * g(3) + (t1 * t296 + t12 * t328 + t297 * t3 + t38 * t636 + t42 * t637 + t622 * t87) * m(7) + t602 * t352 - t179 * t571 + t622 * t140 - t478 * t529 + t636 * t155 + t637 * t153 + t387 + t675 * t154 + t668 * t628 + (t380 * t619 + t661) * g(2) + (t385 * t619 + t662) * g(1) + (m(5) * t570 + t663) * t610 + t354 * t29 + t328 * t28 + t296 * t33 + t297 * t35 - t240 * t262 + t241 * t263 - t61 * t156 + t36 * t510 + t96 * t568 + t97 * t569; t599 + (-pkin(4) * t24 - g(1) * t436 - g(2) * t437 - t124 * t130 - t54 * t72 - t55 * t73) * m(6) - pkin(4) * t29 + t602 * pkin(10) + t633 * t140 + t634 * t155 + (t1 * t333 - t12 * t353 + t3 * t335 + t38 * t634 + t42 * t635 + t633 * t87) * m(7) + t635 * t153 + t679 + t628 * t130 + t610 * t429 + (t606 + t657) * g(3) - t353 * t28 + t333 * t33 + t335 * t35 - t129 * t218 - t73 * t154 - t72 * t156 + t36 * t561; (-t271 * t608 + t272 * t650) * g(1) + (t269 * t608 - t270 * t650) * g(2) + (t210 * t647 - t211 * t645) * t587 + (t210 * t649 - t690) * t589 + (-t211 * t678 + t669 + t692) * t591 + (t156 + t550) * t55 + t600 + (-t154 + t551) * t54 + (t376 * t632 + t428 + t553) * g(3) * t350 + t632 * t1 + (-m(7) * (-t38 + t41) + t155 + t548) * t42 - t87 * (mrSges(7,1) * t211 + mrSges(7,2) * t210) - t124 * (mrSges(6,1) * t211 + mrSges(6,2) * t210) - t41 * t153 + t38 * t549 + t670 * t588 + t613 + (t33 + (-m(7) * t87 - t140) * t211) * pkin(5); -t210 * t153 + t211 * t155 + (g(3) * t351 - t42 * t210 + t38 * t211 - t350 * t610 + t12) * m(7) + t28;];
tau  = t2;
