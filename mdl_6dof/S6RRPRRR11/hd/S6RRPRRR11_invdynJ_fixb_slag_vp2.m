% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR11_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR11_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR11_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR11_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:29:17
% EndTime: 2019-03-09 14:30:17
% DurationCPUTime: 37.73s
% Computational Cost: add. (17392->935), mult. (35374->1248), div. (0->0), fcn. (23309->14), ass. (0->423)
t352 = sin(qJ(2));
t331 = t352 * qJD(1);
t325 = pkin(2) * t331;
t357 = cos(qJ(2));
t496 = qJ(3) * t357;
t393 = pkin(8) * t352 - t496;
t222 = qJD(1) * t393 + t325;
t461 = qJD(1) * t357;
t326 = pkin(7) * t461;
t278 = pkin(3) * t461 + t326;
t351 = sin(qJ(4));
t356 = cos(qJ(4));
t159 = -t222 * t351 + t356 * t278;
t480 = t351 * t352;
t385 = pkin(4) * t357 - pkin(9) * t480;
t457 = qJD(4) * t351;
t557 = pkin(2) + pkin(8);
t524 = pkin(9) + t557;
t649 = -qJD(1) * t385 + t524 * t457 - t159;
t160 = t356 * t222 + t351 * t278;
t287 = t524 * t356;
t434 = t356 * t331;
t648 = pkin(9) * t434 + qJD(4) * t287 + t160;
t350 = sin(qJ(5));
t444 = qJD(4) + qJD(5);
t453 = qJD(5) * t350;
t355 = cos(qJ(5));
t475 = t355 * t356;
t182 = -t350 * t457 - t351 * t453 + t444 * t475;
t485 = t350 * t351;
t210 = -t331 * t485 + t355 * t434;
t606 = t210 + t182;
t387 = t350 * t356 + t355 * t351;
t183 = t444 * t387;
t376 = t387 * t352;
t211 = qJD(1) * t376;
t605 = -t211 - t183;
t647 = -mrSges(6,3) - mrSges(7,3);
t286 = t524 * t351;
t452 = qJD(5) * t355;
t612 = t286 * t453 - t287 * t452 + t350 * t649 - t648 * t355;
t187 = -t355 * t286 - t350 * t287;
t611 = -qJD(5) * t187 + t648 * t350 + t355 * t649;
t354 = cos(qJ(6));
t301 = t331 + t444;
t334 = t352 * qJ(3);
t424 = -pkin(1) - t334;
t205 = (-t357 * t557 + t424) * qJD(1);
t324 = pkin(7) * t331;
t277 = -pkin(3) * t331 - t324;
t587 = qJD(3) - t277;
t215 = -qJD(2) * t557 + t587;
t135 = -t205 * t351 + t356 * t215;
t459 = qJD(2) * t356;
t377 = t351 * t461 - t459;
t118 = pkin(9) * t377 + t135;
t310 = t331 + qJD(4);
t110 = pkin(4) * t310 + t118;
t136 = t205 * t356 + t215 * t351;
t268 = -qJD(2) * t351 - t356 * t461;
t119 = pkin(9) * t268 + t136;
t114 = t350 * t119;
t58 = t355 * t110 - t114;
t175 = t268 * t350 - t355 * t377;
t634 = pkin(10) * t175;
t46 = t58 - t634;
t43 = pkin(5) * t301 + t46;
t349 = sin(qJ(6));
t116 = t355 * t119;
t59 = t110 * t350 + t116;
t417 = t355 * t268 + t350 * t377;
t623 = pkin(10) * t417;
t47 = t59 + t623;
t504 = t349 * t47;
t14 = t354 * t43 - t504;
t449 = qJD(1) * qJD(2);
t281 = -t357 * qJDD(1) + t352 * t449;
t165 = qJD(4) * t268 + qJDD(2) * t356 + t281 * t351;
t282 = qJDD(1) * t352 + t357 * t449;
t267 = qJDD(4) + t282;
t495 = qJDD(1) * pkin(1);
t370 = -qJ(3) * t282 - qJD(3) * t331 - t495;
t128 = t281 * t557 + t370;
t264 = t282 * pkin(7);
t425 = qJDD(3) + t264;
t168 = pkin(3) * t282 - qJDD(2) * t557 + t425;
t57 = -qJD(4) * t136 - t128 * t351 + t356 * t168;
t41 = pkin(4) * t267 - pkin(9) * t165 + t57;
t166 = qJD(4) * t377 - qJDD(2) * t351 + t281 * t356;
t456 = qJD(4) * t356;
t56 = t356 * t128 + t351 * t168 - t205 * t457 + t215 * t456;
t45 = pkin(9) * t166 + t56;
t13 = -qJD(5) * t59 - t350 * t45 + t355 * t41;
t255 = qJDD(5) + t267;
t68 = qJD(5) * t417 + t165 * t355 + t166 * t350;
t6 = pkin(5) * t255 - pkin(10) * t68 + t13;
t12 = t110 * t452 - t119 * t453 + t350 * t41 + t355 * t45;
t69 = -qJD(5) * t175 - t165 * t350 + t166 * t355;
t9 = pkin(10) * t69 + t12;
t2 = qJD(6) * t14 + t349 * t6 + t354 * t9;
t646 = t2 * mrSges(7,2);
t502 = t354 * t47;
t15 = t349 * t43 + t502;
t3 = -qJD(6) * t15 - t349 * t9 + t354 * t6;
t645 = t3 * mrSges(7,1);
t644 = t12 * mrSges(6,2);
t643 = t13 * mrSges(6,1);
t621 = -Ifges(4,4) + Ifges(3,5);
t620 = Ifges(4,5) - Ifges(3,6);
t642 = -pkin(5) * t461 - pkin(10) * t605 + t611;
t641 = pkin(10) * t606 - t612;
t509 = t135 * mrSges(5,3);
t508 = t136 * mrSges(5,3);
t130 = t210 * t354 - t211 * t349;
t386 = -t475 + t485;
t603 = -t349 * t387 - t354 * t386;
t74 = qJD(6) * t603 + t182 * t354 - t183 * t349;
t632 = t130 + t74;
t405 = t357 * mrSges(4,2) - t352 * mrSges(4,3);
t409 = mrSges(3,1) * t357 - mrSges(3,2) * t352;
t640 = t405 - t409;
t639 = qJD(3) + t324;
t348 = qJ(4) + qJ(5);
t332 = sin(t348);
t532 = pkin(4) * t351;
t283 = pkin(5) * t332 + t532;
t335 = qJ(6) + t348;
t317 = sin(t335);
t318 = cos(t335);
t333 = cos(t348);
t406 = mrSges(5,1) * t351 + mrSges(5,2) * t356;
t638 = -m(6) * t532 - m(7) * t283 - t332 * mrSges(6,1) - t317 * mrSges(7,1) - t333 * mrSges(6,2) - t318 * mrSges(7,2) - t406;
t359 = -pkin(9) - pkin(8);
t345 = -pkin(10) + t359;
t637 = -m(6) * (-pkin(2) + t359) + m(5) * t557 + mrSges(5,3) - m(7) * (-pkin(2) + t345) - t647;
t108 = t175 * t354 + t349 * t417;
t347 = qJD(2) * qJ(3);
t244 = t347 + t278;
t185 = -pkin(4) * t268 + t244;
t117 = -pkin(5) * t417 + t185;
t235 = qJDD(6) + t255;
t420 = -t175 * t349 + t354 * t417;
t26 = qJD(6) * t420 + t349 * t69 + t354 * t68;
t27 = -qJD(6) * t108 - t349 * t68 + t354 * t69;
t443 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t235;
t383 = t443 + t645 - t646;
t442 = Ifges(6,5) * t68 + Ifges(6,6) * t69 + Ifges(6,3) * t255;
t515 = Ifges(6,4) * t175;
t534 = mrSges(7,3) * t15;
t535 = mrSges(7,3) * t14;
t538 = -t301 / 0.2e1;
t291 = qJD(6) + t301;
t540 = -t291 / 0.2e1;
t547 = -t175 / 0.2e1;
t553 = -t108 / 0.2e1;
t555 = -t420 / 0.2e1;
t101 = Ifges(7,4) * t420;
t54 = Ifges(7,1) * t108 + Ifges(7,5) * t291 + t101;
t565 = -t54 / 0.2e1;
t514 = Ifges(7,4) * t108;
t53 = Ifges(7,2) * t420 + Ifges(7,6) * t291 + t514;
t567 = -t53 / 0.2e1;
t636 = t643 - t644 + t383 + t442 + (Ifges(6,5) * t417 - Ifges(6,6) * t175) * t538 - t185 * (mrSges(6,1) * t175 + mrSges(6,2) * t417) + (-mrSges(7,2) * t117 + Ifges(7,1) * t553 + Ifges(7,4) * t555 + Ifges(7,5) * t540 + t535 + t565) * t420 + (Ifges(6,1) * t417 - t515) * t547 - (mrSges(7,1) * t117 + Ifges(7,4) * t553 + Ifges(7,2) * t555 + Ifges(7,6) * t540 - t534 + t567) * t108;
t512 = Ifges(4,6) * t357;
t395 = -t352 * Ifges(4,2) - t512;
t635 = t15 * mrSges(7,2) + t59 * mrSges(6,2) + Ifges(4,4) * qJD(2) / 0.2e1 + qJD(1) * t395 / 0.2e1 - t14 * mrSges(7,1) - t58 * mrSges(6,1);
t530 = pkin(5) * t175;
t353 = sin(qJ(1));
t633 = g(2) * t353;
t337 = t356 * pkin(4);
t320 = t337 + pkin(3);
t604 = pkin(4) * t456 + t320 * t331 + t639;
t629 = -t12 * t387 + t13 * t386 - t58 * t605 - t59 * t606;
t628 = t357 * t647 + t640;
t169 = Ifges(6,4) * t417;
t627 = -Ifges(6,2) * t175 + t169;
t125 = mrSges(5,1) * t267 - mrSges(5,3) * t165;
t126 = -mrSges(5,2) * t267 + mrSges(5,3) * t166;
t192 = -mrSges(5,2) * t310 + mrSges(5,3) * t268;
t193 = mrSges(5,1) * t310 + mrSges(5,3) * t377;
t390 = t356 * t192 - t351 * t193;
t626 = t390 * qJD(4) + t356 * t125 + t351 * t126;
t177 = t349 * t386 - t354 * t387;
t419 = t210 * t349 + t354 * t211;
t625 = -t14 * t419 - t177 * t2 + t3 * t603;
t571 = t26 / 0.2e1;
t570 = t27 / 0.2e1;
t563 = t68 / 0.2e1;
t562 = t69 / 0.2e1;
t545 = t235 / 0.2e1;
t544 = t255 / 0.2e1;
t186 = t286 * t350 - t355 * t287;
t145 = pkin(10) * t386 + t186;
t146 = -pkin(10) * t387 + t187;
t85 = t145 * t349 + t146 * t354;
t619 = -qJD(6) * t85 + t349 * t641 + t354 * t642;
t84 = t145 * t354 - t146 * t349;
t618 = qJD(6) * t84 + t349 * t642 - t354 * t641;
t573 = m(6) * pkin(4);
t613 = -mrSges(5,1) - t573;
t531 = pkin(4) * t355;
t319 = pkin(5) + t531;
t450 = qJD(6) * t354;
t451 = qJD(6) * t349;
t486 = t349 * t350;
t63 = -t118 * t350 - t116;
t50 = t63 - t623;
t64 = t355 * t118 - t114;
t51 = t64 - t634;
t610 = -t349 * t50 - t354 * t51 + t319 * t450 + (-t350 * t451 + (t354 * t355 - t486) * qJD(5)) * pkin(4);
t484 = t350 * t354;
t609 = t349 * t51 - t354 * t50 - t319 * t451 + (-t350 * t450 + (-t349 * t355 - t484) * qJD(5)) * pkin(4);
t607 = pkin(5) * t606 + t604;
t556 = pkin(3) + pkin(7);
t298 = t556 * t352;
t273 = t356 * t298;
t340 = t357 * pkin(2);
t463 = t340 + t334;
t412 = pkin(8) * t357 + t463;
t257 = -pkin(1) - t412;
t423 = pkin(9) * t357 - t257;
t148 = pkin(4) * t352 + t351 * t423 + t273;
t272 = t351 * t298;
t181 = t356 * t257 + t272;
t471 = t356 * t357;
t156 = -pkin(9) * t471 + t181;
t91 = t350 * t148 + t355 * t156;
t289 = t357 * t317 * mrSges(7,2);
t602 = -t357 * t332 * mrSges(6,2) - t289;
t184 = -mrSges(5,1) * t268 - mrSges(5,2) * t377;
t437 = mrSges(4,1) * t461;
t295 = -qJD(2) * mrSges(4,3) - t437;
t601 = -t295 + t184;
t600 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t461 - t295;
t438 = mrSges(4,1) * t331;
t599 = -mrSges(3,3) * t331 - t438 + (mrSges(3,1) - mrSges(4,2)) * qJD(2);
t513 = Ifges(4,6) * t352;
t598 = t352 * (-Ifges(4,2) * t357 + t513) + t357 * (Ifges(4,3) * t352 - t512);
t597 = t352 * t620 + t357 * t621;
t366 = qJD(6) * t177 - t182 * t349 - t354 * t183;
t596 = -t419 + t366;
t358 = cos(qJ(1));
t478 = t352 * t353;
t218 = t332 * t358 + t333 * t478;
t219 = -t332 * t478 + t333 * t358;
t202 = t317 * t358 + t318 * t478;
t203 = -t317 * t478 + t318 * t358;
t466 = t202 * mrSges(7,1) + t203 * mrSges(7,2);
t595 = -t218 * mrSges(6,1) - t219 * mrSges(6,2) - t466;
t477 = t352 * t358;
t216 = -t332 * t353 + t333 * t477;
t217 = t332 * t477 + t333 * t353;
t200 = -t317 * t353 + t318 * t477;
t201 = t317 * t477 + t318 * t353;
t467 = t200 * mrSges(7,1) - t201 * mrSges(7,2);
t594 = -t216 * mrSges(6,1) + t217 * mrSges(6,2) - t467;
t263 = t281 * pkin(7);
t592 = -t263 * t357 + t264 * t352;
t212 = -qJDD(2) * qJ(3) - qJD(2) * qJD(3) + t263;
t220 = -qJDD(2) * pkin(2) + t425;
t591 = -t212 * t357 + t220 * t352;
t590 = -t351 * t56 - t356 * t57;
t522 = mrSges(7,1) * t318;
t589 = mrSges(6,1) * t333 + t522;
t588 = g(1) * t358 + t633;
t439 = m(4) + m(5) + m(6) + m(7);
t505 = t377 * Ifges(5,4);
t152 = t268 * Ifges(5,2) + t310 * Ifges(5,6) - t505;
t258 = Ifges(5,4) * t268;
t153 = -Ifges(5,1) * t377 + t310 * Ifges(5,5) + t258;
t394 = -t357 * Ifges(4,3) - t513;
t585 = Ifges(4,5) * qJD(2) + qJD(1) * t394 + t356 * t152 + t351 * t153;
t323 = Ifges(3,4) * t461;
t584 = Ifges(3,1) * t331 + Ifges(3,5) * qJD(2) - Ifges(5,5) * t377 + t175 * Ifges(6,5) + t108 * Ifges(7,5) + t268 * Ifges(5,6) + Ifges(6,6) * t417 + Ifges(7,6) * t420 + t310 * Ifges(5,3) + t301 * Ifges(6,3) + t291 * Ifges(7,3) + t323;
t583 = t357 * t444;
t500 = t357 * mrSges(4,3);
t582 = -t500 + t638 * t357 + (m(4) * pkin(2) - mrSges(4,2) + t637) * t352;
t314 = pkin(5) * t333;
t284 = t314 + t337;
t579 = -m(5) * pkin(3) - m(6) * t320 - m(7) * (pkin(3) + t284) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t575 = Ifges(7,4) * t571 + Ifges(7,2) * t570 + Ifges(7,6) * t545;
t574 = Ifges(7,1) * t571 + Ifges(7,4) * t570 + Ifges(7,5) * t545;
t572 = m(7) * pkin(5);
t569 = Ifges(6,4) * t563 + Ifges(6,2) * t562 + Ifges(6,6) * t544;
t568 = Ifges(6,1) * t563 + Ifges(6,4) * t562 + Ifges(6,5) * t544;
t566 = t53 / 0.2e1;
t564 = t54 / 0.2e1;
t561 = -t165 * Ifges(5,4) / 0.2e1 - t166 * Ifges(5,2) / 0.2e1 - t267 * Ifges(5,6) / 0.2e1;
t95 = Ifges(6,2) * t417 + Ifges(6,6) * t301 + t515;
t560 = -t95 / 0.2e1;
t96 = Ifges(6,1) * t175 + Ifges(6,5) * t301 + t169;
t559 = -t96 / 0.2e1;
t558 = t96 / 0.2e1;
t554 = t420 / 0.2e1;
t552 = t108 / 0.2e1;
t551 = t165 / 0.2e1;
t550 = t166 / 0.2e1;
t549 = -t417 / 0.2e1;
t548 = t417 / 0.2e1;
t546 = t175 / 0.2e1;
t543 = t267 / 0.2e1;
t541 = -t377 / 0.2e1;
t539 = t291 / 0.2e1;
t537 = t301 / 0.2e1;
t533 = pkin(4) * t377;
t529 = pkin(7) * t352;
t338 = t357 * pkin(7);
t521 = mrSges(6,3) * t417;
t520 = mrSges(6,3) * t175;
t519 = Ifges(3,4) * t352;
t518 = Ifges(3,4) * t357;
t517 = Ifges(5,4) * t351;
t516 = Ifges(5,4) * t356;
t487 = t345 * t357;
t479 = t351 * t357;
t476 = t353 * t356;
t470 = t356 * t358;
t469 = t357 * t358;
t468 = t357 * t359;
t299 = t357 * pkin(3) + t338;
t462 = t358 * pkin(1) + t353 * pkin(7);
t460 = qJD(2) * t352;
t458 = qJD(2) * t357;
t455 = qJD(4) * t357;
t436 = Ifges(5,5) * t165 + Ifges(5,6) * t166 + Ifges(5,3) * t267;
t237 = pkin(4) * t471 + t299;
t433 = t351 * t455;
t427 = -t456 / 0.2e1;
t313 = qJ(3) + t532;
t422 = -t449 / 0.2e1;
t227 = t282 * mrSges(4,1) + qJDD(2) * mrSges(4,2);
t90 = t355 * t148 - t156 * t350;
t415 = pkin(2) * t460 - qJD(3) * t352;
t198 = qJD(2) * t393 + t415;
t280 = t556 * t458;
t418 = -t198 * t351 + t356 * t280;
t408 = mrSges(3,1) * t352 + mrSges(3,2) * t357;
t407 = mrSges(5,1) * t356 - mrSges(5,2) * t351;
t404 = Ifges(5,1) * t356 - t517;
t403 = Ifges(5,1) * t351 + t516;
t402 = t357 * Ifges(3,2) + t519;
t400 = -Ifges(5,2) * t351 + t516;
t399 = Ifges(5,2) * t356 + t517;
t397 = Ifges(5,5) * t356 - Ifges(5,6) * t351;
t396 = Ifges(5,5) * t351 + Ifges(5,6) * t356;
t224 = t387 * t357;
t70 = pkin(5) * t352 + pkin(10) * t224 + t90;
t223 = t386 * t357;
t71 = pkin(10) * t223 + t91;
t37 = -t349 * t71 + t354 * t70;
t38 = t349 * t70 + t354 * t71;
t392 = t135 * t351 - t136 * t356;
t143 = t223 * t354 + t224 * t349;
t144 = t223 * t349 - t224 * t354;
t285 = -qJD(2) * pkin(2) + t639;
t293 = -t326 - t347;
t388 = t285 * t357 + t293 * t352;
t384 = t424 - t340;
t382 = pkin(1) * t408;
t246 = -t351 * t353 + t352 * t470;
t248 = t351 * t358 + t352 * t476;
t245 = t384 * qJD(1);
t381 = t245 * (-mrSges(4,2) * t352 - t500);
t380 = t352 * (Ifges(3,1) * t357 - t519);
t80 = t385 * qJD(2) + (t356 * t423 - t272) * qJD(4) + t418;
t373 = t352 * t459 + t433;
t98 = t356 * t198 - t257 * t457 + t351 * t280 + t298 * t456;
t86 = pkin(9) * t373 + t98;
t30 = t148 * t452 - t156 * t453 + t350 * t80 + t355 * t86;
t372 = t351 * t460 - t356 * t455;
t369 = Ifges(5,5) * t357 + t352 * t403;
t368 = Ifges(5,6) * t357 + t352 * t399;
t367 = Ifges(5,3) * t357 + t352 * t396;
t170 = -pkin(3) * t281 - t212;
t31 = -qJD(5) * t91 - t350 * t86 + t355 * t80;
t188 = -pkin(4) * t433 + (-pkin(7) - t320) * t460;
t109 = -pkin(4) * t166 + t170;
t364 = -qJD(4) * t392 - t590;
t288 = -pkin(1) - t463;
t279 = t556 * t460;
t275 = -qJ(3) * t461 + t325;
t274 = t405 * qJD(1);
t256 = t407 * t357;
t249 = -t351 * t478 + t470;
t247 = t351 * t477 + t476;
t239 = Ifges(3,6) * qJD(2) + qJD(1) * t402;
t233 = pkin(4) * t484 + t319 * t349;
t232 = -pkin(4) * t486 + t319 * t354;
t228 = -qJ(3) * t458 + t415;
t226 = mrSges(4,1) * t281 - qJDD(2) * mrSges(4,3);
t213 = pkin(5) * t387 + t313;
t180 = -t257 * t351 + t273;
t167 = -pkin(5) * t223 + t237;
t161 = pkin(2) * t281 + t370;
t140 = mrSges(6,1) * t301 - t520;
t139 = -mrSges(6,2) * t301 + t521;
t129 = -t533 + t530;
t121 = -t386 * t460 + t387 * t583;
t120 = qJD(2) * t376 + t386 * t583;
t111 = -mrSges(6,1) * t417 + mrSges(6,2) * t175;
t99 = -qJD(4) * t181 + t418;
t97 = -mrSges(5,1) * t166 + mrSges(5,2) * t165;
t93 = mrSges(7,1) * t291 - mrSges(7,3) * t108;
t92 = -mrSges(7,2) * t291 + mrSges(7,3) * t420;
t89 = t165 * Ifges(5,1) + t166 * Ifges(5,4) + t267 * Ifges(5,5);
t87 = -pkin(5) * t121 + t188;
t61 = -mrSges(6,2) * t255 + mrSges(6,3) * t69;
t60 = mrSges(6,1) * t255 - mrSges(6,3) * t68;
t55 = -mrSges(7,1) * t420 + mrSges(7,2) * t108;
t49 = -qJD(6) * t144 - t120 * t349 + t121 * t354;
t48 = qJD(6) * t143 + t120 * t354 + t121 * t349;
t42 = -pkin(5) * t69 + t109;
t36 = -mrSges(6,1) * t69 + mrSges(6,2) * t68;
t23 = -mrSges(7,2) * t235 + mrSges(7,3) * t27;
t22 = mrSges(7,1) * t235 - mrSges(7,3) * t26;
t19 = pkin(10) * t121 + t30;
t18 = pkin(5) * t458 - pkin(10) * t120 + t31;
t17 = t354 * t46 - t504;
t16 = -t349 * t46 - t502;
t10 = -mrSges(7,1) * t27 + mrSges(7,2) * t26;
t5 = -qJD(6) * t38 + t18 * t354 - t19 * t349;
t4 = qJD(6) * t37 + t18 * t349 + t19 * t354;
t1 = [(t585 / 0.2e1 + t293 * mrSges(4,1) - t600 * pkin(7) - t239 / 0.2e1) * t460 + t409 * t495 + m(4) * (t161 * t288 + t228 * t245 + (qJD(2) * t388 + t591) * pkin(7)) + (-t281 * t338 + t282 * t529 + t592) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t592) + t597 * qJD(2) ^ 2 / 0.2e1 + t598 * t422 + t373 * t508 + t352 * t643 + t352 * t645 + (-t249 * mrSges(5,1) - t219 * mrSges(6,1) - t203 * mrSges(7,1) + t248 * mrSges(5,2) + t218 * mrSges(6,2) + t202 * mrSges(7,2) + (-m(6) * (-t313 * t352 - pkin(1)) - m(5) * t424 - m(7) * (-pkin(1) + (-qJ(3) - t283) * t352) + m(3) * pkin(1) - m(4) * t384 + mrSges(2,1) + t637 * t357 - t640) * t353 + ((-m(3) - t439) * pkin(7) + t579) * t358) * g(1) + (Ifges(6,5) * t120 + Ifges(6,6) * t121) * t537 + (-t14 * t48 + t143 * t2 - t144 * t3 + t15 * t49) * mrSges(7,3) + (Ifges(7,4) * t48 + Ifges(7,2) * t49) * t554 + (Ifges(6,1) * t120 + Ifges(6,4) * t121) * t546 + (Ifges(7,5) * t48 + Ifges(7,6) * t49) * t539 + t144 * t574 + t143 * t575 + t471 * t561 + t48 * t564 + t49 * t566 - t224 * t568 + t223 * t569 + (Ifges(7,4) * t144 + Ifges(7,2) * t143 + Ifges(7,6) * t352) * t570 + (Ifges(7,1) * t144 + Ifges(7,4) * t143 + Ifges(7,5) * t352) * t571 + (Ifges(5,6) * t352 - t357 * t399) * t550 + (Ifges(5,5) * t352 - t357 * t403) * t551 + t120 * t558 + t591 * mrSges(4,1) + t56 * (-mrSges(5,2) * t352 - mrSges(5,3) * t471) + (Ifges(6,4) * t120 + Ifges(6,2) * t121) * t548 + (t285 * mrSges(4,1) - t599 * pkin(7) - t136 * mrSges(5,2) + t135 * mrSges(5,1) + Ifges(7,5) * t552 + Ifges(7,6) * t554 + t584 / 0.2e1 + Ifges(6,5) * t546 + Ifges(6,6) * t548 + Ifges(6,3) * t537 + Ifges(7,3) * t539 - t635) * t458 + t268 * (qJD(2) * t368 - t400 * t455) / 0.2e1 + t310 * (qJD(2) * t367 - t397 * t455) / 0.2e1 - t382 * t449 + t282 * t518 / 0.2e1 + (-Ifges(3,4) * t281 + Ifges(3,5) * qJDD(2) + t436 + t442 + t443) * t352 / 0.2e1 + t167 * t10 + (Ifges(7,1) * t48 + Ifges(7,4) * t49) * t552 + (-qJDD(2) * mrSges(3,2) - t226) * t338 + (-qJDD(2) * mrSges(3,1) + t227) * t529 + t244 * (-mrSges(5,1) * t373 + mrSges(5,2) * t372) + (-Ifges(6,4) * t224 + Ifges(6,2) * t223 + Ifges(6,6) * t352) * t562 + (-Ifges(6,1) * t224 + Ifges(6,4) * t223 + Ifges(6,5) * t352) * t563 + (-Ifges(6,5) * t224 + Ifges(6,6) * t223 + Ifges(6,3) * t352) * t544 + t109 * (-mrSges(6,1) * t223 - mrSges(6,2) * t224) + t42 * (-mrSges(7,1) * t143 + mrSges(7,2) * t144) + (t357 * (-Ifges(3,2) * t352 + t518) + t380) * t449 / 0.2e1 + (-m(3) * t462 - t247 * mrSges(5,1) - t217 * mrSges(6,1) - t201 * mrSges(7,1) - t246 * mrSges(5,2) - t216 * mrSges(6,2) - t200 * mrSges(7,2) + (-m(5) * pkin(8) - mrSges(5,3)) * t469 - t439 * (pkin(2) * t469 + qJ(3) * t477 + t462) + t579 * t353 + (-m(6) * (pkin(4) * t480 - t468) - m(7) * (t283 * t352 - t487) - mrSges(2,1) + t628) * t358) * g(2) + t30 * t139 + t31 * t140 + t121 * t95 / 0.2e1 + t117 * (-mrSges(7,1) * t49 + mrSges(7,2) * t48) + t91 * t61 + t4 * t92 + t5 * t93 + Ifges(2,3) * qJDD(1) + t87 * t55 + t90 * t60 + t37 * t22 + t38 * t23 + (qJD(2) * t369 - t404 * t455) * t541 + (Ifges(5,3) * t352 - t357 * t396) * t543 + (Ifges(7,5) * t144 + Ifges(7,6) * t143 + Ifges(7,3) * t352) * t545 + t282 * t352 * Ifges(3,1) + t237 * t36 + t170 * t256 + (t12 * t223 - t120 * t58 + t121 * t59 + t13 * t224) * mrSges(6,3) + t161 * t405 - t281 * t402 / 0.2e1 + t281 * t394 / 0.2e1 - t282 * t395 / 0.2e1 - t89 * t479 / 0.2e1 + t57 * (mrSges(5,1) * t352 + mrSges(5,3) * t479) + t228 * t274 - t279 * t184 - pkin(1) * (mrSges(3,1) * t281 + mrSges(3,2) * t282) + (t352 * t621 - t357 * t620) * qJDD(2) / 0.2e1 + m(6) * (t109 * t237 + t12 * t91 + t13 * t90 + t185 * t188 + t30 * t59 + t31 * t58) + m(7) * (t117 * t87 + t14 * t5 + t15 * t4 + t167 * t42 + t2 * t38 + t3 * t37) + m(5) * (t135 * t99 + t136 * t98 + t170 * t299 + t180 * t57 + t181 * t56 - t244 * t279) + t357 * t153 * t427 + t299 * t97 + t152 * t433 / 0.2e1 + t288 * (-mrSges(4,2) * t281 - mrSges(4,3) * t282) - t372 * t509 - t352 * (Ifges(4,4) * qJDD(2) - Ifges(4,2) * t282 + Ifges(4,6) * t281) / 0.2e1 - t352 * t644 + qJD(2) * t381 + t357 * (Ifges(3,4) * t282 - Ifges(3,2) * t281 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t357 * (Ifges(4,5) * qJDD(2) - Ifges(4,6) * t282 + Ifges(4,3) * t281) / 0.2e1 - t352 * t646 + t180 * t125 + t181 * t126 + t185 * (-mrSges(6,1) * t121 + mrSges(6,2) * t120) + t188 * t111 + t98 * t192 + t99 * t193; (-qJ(3) * t439 * t469 + t358 * t582) * g(1) + t611 * t140 + (t109 * t313 + t12 * t187 + t13 * t186 + t185 * t604 + t58 * t611 + t59 * t612) * m(6) + t612 * t139 + t600 * t324 + t601 * qJD(3) + t604 * t111 + (mrSges(6,1) * t606 + mrSges(6,2) * t605) * t185 + t606 * t560 + t607 * t55 - t585 * t331 / 0.2e1 + t310 * t407 * t244 + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + (Ifges(6,4) * t211 + Ifges(6,2) * t210) * t549 + (t382 - t380 / 0.2e1 + t598 / 0.2e1) * qJD(1) ^ 2 - (t268 * t399 + t310 * t396 - t377 * t403) * qJD(4) / 0.2e1 - (t268 * t368 + t310 * t367 - t369 * t377) * qJD(1) / 0.2e1 - t387 * t569 + (-Ifges(6,4) * t386 - Ifges(6,2) * t387) * t562 + (-Ifges(6,1) * t386 - Ifges(6,4) * t387) * t563 + (-Ifges(6,5) * t386 - Ifges(6,6) * t387) * t544 + t109 * (mrSges(6,1) * t387 - mrSges(6,2) * t386) - t386 * t568 + (Ifges(7,5) * t603 + Ifges(7,6) * t177) * t545 + t42 * (-mrSges(7,1) * t177 + mrSges(7,2) * t603) + (Ifges(7,4) * t603 + Ifges(7,2) * t177) * t570 + (Ifges(7,1) * t603 + Ifges(7,4) * t177) * t571 + t603 * t574 + t588 * t408 + t590 * mrSges(5,3) + t597 * t422 + t599 * t326 + (-Ifges(6,4) * t183 - Ifges(6,2) * t182) * t548 + (-Ifges(6,1) * t183 - Ifges(6,4) * t182) * t546 + (-Ifges(6,5) * t183 - Ifges(6,6) * t182) * t537 + (Ifges(6,5) * t211 + Ifges(6,6) * t210) * t538 + (-t135 * (mrSges(5,1) * t357 - mrSges(5,3) * t480) - t136 * (mrSges(5,3) * t352 * t356 - mrSges(5,2) * t357) - m(4) * t388 * pkin(7) - t381) * qJD(1) + (t170 * qJ(3) - t135 * t159 - t136 * t160 + t244 * t587) * m(5) + t177 * t575 + t351 * t561 + t130 * t567 + t400 * t550 + t404 * t551 - t183 * t558 + t211 * t559 + (Ifges(6,5) * t547 + Ifges(7,5) * t553 + Ifges(6,6) * t549 + Ifges(7,6) * t555 + Ifges(6,3) * t538 + Ifges(7,3) * t540 + t635) * t461 + (-pkin(2) * t220 - qJ(3) * t212 - qJD(3) * t293 - t245 * t275) * m(4) - (Ifges(7,4) * t552 + Ifges(7,2) * t554 + Ifges(7,6) * t539 + t534 + t566) * t74 + (-t130 * t15 - t625) * mrSges(7,3) + (-m(5) * t364 - t626) * t557 - t285 * t437 + t239 * t331 / 0.2e1 + (Ifges(6,1) * t211 + Ifges(6,4) * t210) * t547 + (t97 - t226) * qJ(3) + (-t153 / 0.2e1 + t509) * t457 - t456 * t508 + (Ifges(7,1) * t552 + Ifges(7,4) * t554 + Ifges(7,5) * t539 - t535 + t564) * t366 + (Ifges(7,5) * t419 + Ifges(7,6) * t130) * t540 + t419 * t565 + (Ifges(7,4) * t419 + Ifges(7,2) * t130) * t555 + (Ifges(7,1) * t419 + Ifges(7,4) * t130) * t553 - t293 * t438 + (-m(6) * (t463 - t468) - m(4) * t463 - m(7) * (t463 - t487) - m(5) * t412 - t357 * mrSges(5,3) + t638 * t352 + t628) * g(3) + t629 * mrSges(6,3) - t212 * mrSges(4,3) + t213 * t10 - (-Ifges(3,2) * t331 + t323 + t584) * t461 / 0.2e1 + t84 * t22 + t85 * t23 + t397 * t543 + t220 * mrSges(4,2) - pkin(2) * t227 + t263 * mrSges(3,2) - t264 * mrSges(3,1) + t170 * t406 - t275 * t274 - t277 * t184 + t313 * t36 + t618 * t92 + t619 * t93 + (t607 * t117 + t619 * t14 + t618 * t15 + t2 * t85 + t213 * t42 + t3 * t84) * m(7) + (mrSges(7,1) * t632 + t596 * mrSges(7,2)) * t117 + t620 * t281 + t621 * t282 + t356 * t89 / 0.2e1 + t152 * t427 + (-t439 * t496 + t582) * t633 + t186 * t60 + t187 * t61 - t160 * t192 - t159 * t193; ((t274 + t390) * qJD(1) - t588 * t439) * t352 - t177 * t23 + t605 * t140 + t606 * t139 + t439 * t357 * g(3) + t632 * t92 + t596 * t93 + t603 * t22 + t227 + (-t111 - t55 - t601) * qJD(2) + t387 * t61 - t386 * t60 + (-qJD(2) * t117 + t14 * t366 + t15 * t632 + t625) * m(7) + (-qJD(2) * t185 - t629) * m(6) + (-qJD(2) * t244 - t331 * t392 + t364) * m(5) + (qJD(2) * t293 + t245 * t331 + t220) * m(4) + t626; t609 * t93 + (-t117 * t129 + t14 * t609 + t15 * t610 + t2 * t233 + t232 * t3) * m(7) + t610 * t92 + (t256 + (m(6) * t337 + m(7) * t284 + t589) * t357 + t602) * g(3) + (t139 * t452 - t140 * t453 + t350 * t61) * pkin(4) - t175 * t560 + (t175 * t59 + t417 * t58) * mrSges(6,3) - (Ifges(5,2) * t377 + t153 + t258) * t268 / 0.2e1 - t244 * (-mrSges(5,1) * t377 + mrSges(5,2) * t268) - t310 * (Ifges(5,5) * t268 + Ifges(5,6) * t377) / 0.2e1 - t377 * t508 + t377 * (Ifges(5,1) * t268 + t505) / 0.2e1 + (-m(7) * (-t283 * t353 + t284 * t477) + mrSges(5,2) * t247 + t613 * t246 + t594) * g(1) + (-mrSges(5,2) * t249 - m(7) * (t283 * t358 + t284 * t478) + t613 * t248 + t595) * g(2) + t636 + t436 + (t12 * t350 + t13 * t355 + (-t350 * t58 + t355 * t59) * qJD(5)) * t573 + t627 * t549 + t111 * t533 + t417 * t559 - t64 * t139 - t63 * t140 - t129 * t55 - t56 * mrSges(5,2) + t57 * mrSges(5,1) + t152 * t541 + t60 * t531 + t268 * t509 + t232 * t22 + t233 * t23 - m(6) * (-t185 * t533 + t58 * t63 + t59 * t64) - t135 * t192 + t136 * t193; (t2 * t349 + t3 * t354 + (-t14 * t349 + t15 * t354) * qJD(6)) * t572 - t55 * t530 - m(7) * (t117 * t530 + t14 * t16 + t15 * t17) - t17 * t92 - t16 * t93 + t95 * t546 + (t140 + t520) * t59 + (-t139 + t521) * t58 + (t627 + t96) * t549 + ((m(7) * t314 + t589) * t357 + t602) * g(3) + (-t218 * t572 + t595) * g(2) + (-t216 * t572 + t594) * g(1) + (t22 * t354 + t23 * t349 + t450 * t92 - t451 * t93) * pkin(5) + t636; -t117 * (mrSges(7,1) * t108 + mrSges(7,2) * t420) + (Ifges(7,1) * t420 - t514) * t553 + t53 * t552 + (Ifges(7,5) * t420 - Ifges(7,6) * t108) * t540 - t14 * t92 + t15 * t93 - g(1) * t467 - g(2) * t466 - g(3) * (-t357 * t522 + t289) + (t108 * t15 + t14 * t420) * mrSges(7,3) + t383 + (-Ifges(7,2) * t108 + t101 + t54) * t555;];
tau  = t1;
