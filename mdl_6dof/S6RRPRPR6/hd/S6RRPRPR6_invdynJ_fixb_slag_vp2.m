% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:38:12
% EndTime: 2019-03-09 10:39:17
% DurationCPUTime: 39.91s
% Computational Cost: add. (15772->946), mult. (43580->1207), div. (0->0), fcn. (35398->12), ass. (0->448)
t314 = sin(qJ(4));
t443 = qJD(4) * t314;
t315 = sin(qJ(2));
t319 = cos(qJ(2));
t311 = sin(pkin(6));
t449 = qJD(1) * t311;
t482 = sin(pkin(11));
t386 = t482 * t449;
t483 = cos(pkin(11));
t387 = t483 * t449;
t241 = -t315 * t386 + t319 * t387;
t476 = t241 * t314;
t582 = t443 - t476;
t242 = -t315 * t387 - t319 * t386;
t312 = cos(pkin(6));
t435 = qJDD(1) * t312;
t300 = qJDD(2) + t435;
t318 = cos(qJ(4));
t437 = qJD(1) * qJD(2);
t257 = (qJDD(1) * t319 - t315 * t437) * t311;
t258 = (qJDD(1) * t315 + t319 * t437) * t311;
t332 = t257 * t482 + t258 * t483;
t395 = qJD(1) * t312 + qJD(2);
t352 = t318 * t395;
t106 = -qJD(4) * t352 - t242 * t443 - t314 * t300 - t318 * t332;
t548 = -t106 / 0.2e1;
t632 = Ifges(5,1) * t548;
t210 = -t318 * t242 + t314 * t395;
t444 = qJD(4) * t210;
t107 = -t318 * t300 + t314 * t332 + t444;
t546 = -t107 / 0.2e1;
t631 = Ifges(5,4) * t546;
t235 = qJD(4) - t241;
t189 = t257 * t483 - t258 * t482;
t184 = qJDD(4) - t189;
t539 = t184 / 0.2e1;
t601 = Ifges(6,4) - Ifges(5,5);
t600 = Ifges(5,6) - Ifges(6,5);
t618 = -Ifges(5,3) - Ifges(6,1);
t521 = pkin(1) * t312;
t304 = t319 * t521;
t298 = qJD(1) * t304;
t514 = pkin(8) + qJ(3);
t412 = t514 * t315;
t391 = t311 * t412;
t227 = -qJD(1) * t391 + t298;
t467 = t312 * t315;
t303 = pkin(1) * t467;
t469 = t311 * t319;
t228 = (t469 * t514 + t303) * qJD(1);
t402 = t483 * t228;
t144 = t227 * t482 + t402;
t630 = t582 * pkin(4) - qJD(5) * t314 - t144;
t209 = -t242 * t314 - t352;
t313 = sin(qJ(6));
t317 = cos(qJ(6));
t129 = t209 * t317 - t235 * t313;
t130 = t209 * t313 + t235 * t317;
t204 = Ifges(5,4) * t209;
t208 = qJD(6) + t210;
t598 = t210 * Ifges(5,1) + t235 * Ifges(5,5) + t130 * Ifges(7,5) + t129 * Ifges(7,6) + t208 * Ifges(7,3) - t204;
t203 = Ifges(6,6) * t209;
t94 = t235 * Ifges(6,4) - t210 * Ifges(6,2) + t203;
t629 = t598 / 0.2e1 - t94 / 0.2e1;
t547 = t106 / 0.2e1;
t545 = t107 / 0.2e1;
t628 = -t184 / 0.2e1;
t105 = qJDD(6) - t106;
t42 = qJD(6) * t129 + t107 * t313 + t184 * t317;
t43 = -qJD(6) * t130 + t107 * t317 - t184 * t313;
t7 = Ifges(7,5) * t42 + Ifges(7,6) * t43 + Ifges(7,3) * t105;
t627 = t632 + t631 + Ifges(5,5) * t539 + t7 / 0.2e1;
t626 = t630 + t235 * (pkin(10) * t314 - qJ(5) * t318);
t218 = t482 * t228;
t145 = t227 * t483 - t218;
t136 = t314 * t145;
t421 = t315 * t449;
t397 = pkin(2) * t421;
t157 = -pkin(3) * t242 - pkin(9) * t241 + t397;
t422 = t482 * pkin(2);
t307 = t422 + pkin(9);
t515 = pkin(5) + t307;
t273 = t515 * t318;
t551 = pkin(4) + pkin(10);
t625 = qJD(4) * t273 - t136 - (pkin(5) * t241 - t157) * t318 - t551 * t242;
t442 = qJD(4) * t318;
t475 = t241 * t318;
t583 = t442 - t475;
t339 = t312 * pkin(2) - t391;
t211 = qJD(2) * pkin(2) + qJD(1) * t339 + t298;
t125 = t482 * t211 + t402;
t115 = pkin(9) * t395 + t125;
t309 = pkin(2) * t319 + pkin(1);
t263 = -t309 * t449 + qJD(3);
t325 = -pkin(3) * t241 + pkin(9) * t242 + t263;
t67 = t115 * t314 - t318 * t325;
t350 = pkin(5) * t210 + t67;
t610 = t350 + qJD(5);
t38 = -t235 * t551 + t610;
t124 = t211 * t483 - t218;
t114 = -pkin(3) * t395 - t124;
t323 = -t210 * qJ(5) + t114;
t49 = t209 * t551 + t323;
t14 = -t313 * t49 + t317 * t38;
t15 = t313 * t38 + t317 * t49;
t569 = t14 * mrSges(7,1) - t15 * mrSges(7,2);
t624 = t569 + t629;
t623 = Ifges(5,6) * t628 + 0.2e1 * Ifges(6,3) * t545 + (Ifges(5,4) + Ifges(6,6)) * t547 + (t545 - t546) * Ifges(5,2) + (-t600 + Ifges(6,5)) * t539;
t25 = mrSges(7,1) * t105 - mrSges(7,3) * t42;
t26 = -mrSges(7,2) * t105 + mrSges(7,3) * t43;
t90 = -mrSges(7,2) * t208 + mrSges(7,3) * t129;
t91 = mrSges(7,1) * t208 - mrSges(7,3) * t130;
t580 = -t313 * t91 + t317 * t90;
t622 = -qJD(6) * t580 - t317 * t25 - t313 * t26;
t549 = t105 / 0.2e1;
t557 = t43 / 0.2e1;
t558 = t42 / 0.2e1;
t450 = pkin(8) * t469 + t303;
t256 = t450 * qJD(2);
t430 = pkin(1) * t435;
t296 = t319 * t430;
t471 = t311 * t315;
t418 = qJD(3) * t471;
t436 = qJDD(1) * t311;
t429 = pkin(8) * t436;
t121 = -t315 * t429 + pkin(2) * t300 - qJ(3) * t258 + t296 + (-t256 - t418) * qJD(1);
t434 = qJD(2) * t521;
t394 = qJD(1) * t434;
t424 = t315 * t430 + (t394 + t429) * t319;
t446 = qJD(3) * t319;
t447 = qJD(2) * t315;
t133 = qJ(3) * t257 + (-pkin(8) * t447 + t446) * t449 + t424;
t65 = t121 * t483 - t482 * t133;
t61 = -t300 * pkin(3) - t65;
t324 = t106 * qJ(5) - t210 * qJD(5) + t61;
t11 = t107 * t551 + t324;
t66 = t482 * t121 + t483 * t133;
t614 = pkin(9) * t300 + qJD(4) * t325 + t66;
t223 = -pkin(1) * t436 - t257 * pkin(2) + qJDD(3);
t87 = -t189 * pkin(3) - pkin(9) * t332 + t223;
t17 = -t115 * t442 - t314 * t614 + t318 * t87;
t343 = qJDD(5) - t17;
t5 = -pkin(5) * t106 - t184 * t551 + t343;
t1 = qJD(6) * t14 + t11 * t317 + t313 * t5;
t2 = -qJD(6) * t15 - t11 * t313 + t317 * t5;
t570 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t621 = -t539 * t601 + t570 + t632 + Ifges(6,4) * t628 + Ifges(7,5) * t558 + Ifges(6,6) * t546 + Ifges(7,6) * t557 + Ifges(7,3) * t549 + (-t547 + t548) * Ifges(6,2);
t617 = -t209 * t600 - t210 * t601 - t235 * t618;
t234 = Ifges(4,4) * t241;
t616 = Ifges(4,2) * t241;
t615 = t263 * mrSges(4,2);
t274 = t315 * t482 - t319 * t483;
t248 = t274 * t311;
t333 = t315 * t483 + t319 * t482;
t613 = t106 * t601 - t107 * t600 - t184 * t618;
t320 = cos(qJ(1));
t457 = t319 * t320;
t316 = sin(qJ(1));
t462 = t315 * t316;
t612 = t312 * t457 - t462;
t572 = m(7) * pkin(10) + mrSges(5,1) - mrSges(6,2);
t603 = m(7) + m(6);
t611 = pkin(4) * t603 + t572;
t64 = t209 * pkin(4) + t323;
t609 = t114 * mrSges(5,1) - t64 * mrSges(6,2);
t608 = t114 * mrSges(5,2) - t64 * mrSges(6,3);
t16 = -t115 * t443 + t314 * t87 + t318 * t614;
t10 = -qJ(5) * t184 - qJD(5) * t235 - t16;
t12 = -pkin(4) * t184 + t343;
t606 = t17 * mrSges(5,1) - t16 * mrSges(5,2) + t12 * mrSges(6,2) - t10 * mrSges(6,3);
t55 = -pkin(4) * t235 + qJD(5) + t67;
t68 = t318 * t115 + t314 * t325;
t58 = -t235 * qJ(5) - t68;
t605 = t263 * mrSges(4,1) - t67 * mrSges(5,1) - t68 * mrSges(5,2) + t55 * mrSges(6,2) - t58 * mrSges(6,3);
t602 = mrSges(5,2) - mrSges(6,3);
t72 = -mrSges(5,2) * t184 - mrSges(5,3) * t107;
t73 = mrSges(6,1) * t107 - mrSges(6,3) * t184;
t599 = -t73 + t72;
t597 = t395 * Ifges(4,5);
t596 = t395 * Ifges(4,6);
t423 = t483 * pkin(2);
t308 = -t423 - pkin(3);
t464 = t314 * qJ(5);
t340 = t308 - t464;
t250 = -t318 * t551 + t340;
t272 = t515 * t314;
t194 = -t250 * t313 + t272 * t317;
t595 = qJD(6) * t194 + t313 * t625 + t317 * t626;
t195 = t250 * t317 + t272 * t313;
t594 = -qJD(6) * t195 - t313 * t626 + t317 * t625;
t82 = t318 * t145 + t314 * t157;
t69 = qJ(5) * t242 - t82;
t593 = pkin(5) * t476 - t515 * t443 + t69;
t512 = mrSges(6,1) * t209;
t138 = -mrSges(6,3) * t235 + t512;
t80 = -mrSges(7,1) * t129 + mrSges(7,2) * t130;
t592 = t80 - t138;
t591 = -qJ(5) * t583 + t630;
t490 = t242 * mrSges(4,3);
t455 = -mrSges(4,1) * t395 + mrSges(5,1) * t209 + mrSges(5,2) * t210 - t490;
t511 = mrSges(6,1) * t210;
t139 = mrSges(6,2) * t235 + t511;
t509 = mrSges(5,3) * t210;
t141 = mrSges(5,1) * t235 - t509;
t454 = -t139 + t141;
t463 = t314 * t317;
t168 = t241 * t463 + t242 * t313;
t440 = qJD(6) * t318;
t415 = t313 * t440;
t590 = -t317 * t443 + t168 - t415;
t466 = t313 * t314;
t169 = t241 * t466 - t242 * t317;
t589 = -t313 * t443 + t317 * t440 + t169;
t329 = t312 * t274;
t197 = -t316 * t333 - t320 * t329;
t479 = t197 * t318;
t588 = pkin(4) * t479 + t197 * t464;
t200 = t316 * t329 - t320 * t333;
t478 = t200 * t318;
t587 = pkin(4) * t478 + t200 * t464;
t474 = t248 * t318;
t586 = -pkin(4) * t474 - t248 * t464;
t585 = t114 * (mrSges(5,1) * t314 + mrSges(5,2) * t318) + t64 * (-mrSges(6,2) * t314 - mrSges(6,3) * t318);
t584 = -t314 * t600 - t318 * t601;
t579 = t16 * t318 - t17 * t314;
t578 = -t10 * t318 + t12 * t314;
t507 = Ifges(3,4) * t315;
t577 = -t315 * (Ifges(3,1) * t319 - t507) / 0.2e1 + pkin(1) * (mrSges(3,1) * t315 + mrSges(3,2) * t319);
t576 = mrSges(6,1) + mrSges(5,3) - mrSges(4,2);
t438 = m(5) + t603;
t574 = m(7) * pkin(5) + pkin(9) * t438 + t576;
t378 = t317 * mrSges(7,1) - t313 * mrSges(7,2);
t519 = t209 * pkin(5);
t44 = -t58 - t519;
t126 = Ifges(7,4) * t129;
t54 = t130 * Ifges(7,1) + t208 * Ifges(7,5) + t126;
t488 = t313 * t54;
t573 = t44 * t378 - t488 / 0.2e1;
t452 = t333 * t312;
t201 = -t320 * t274 - t316 * t452;
t196 = t316 * t274 - t320 * t452;
t496 = t14 * t313;
t359 = -t15 * t317 + t496;
t520 = t2 * t317;
t326 = -qJD(6) * t359 + t1 * t313 + t520;
t396 = pkin(8) * t421;
t205 = -qJD(2) * t396 + t424;
t206 = -pkin(8) * t258 - t315 * t394 + t296;
t571 = t206 * mrSges(3,1) + t65 * mrSges(4,1) - t205 * mrSges(3,2) - t66 * mrSges(4,2) + Ifges(3,5) * t258 + Ifges(4,5) * t332 + Ifges(3,6) * t257;
t74 = -t106 * mrSges(6,1) + t184 * mrSges(6,2);
t568 = t622 - t74;
t377 = mrSges(7,1) * t313 + mrSges(7,2) * t317;
t338 = m(7) * qJ(5) + t377;
t567 = -m(6) * qJ(5) - t338 + t602;
t468 = t311 * t320;
t170 = -t196 * t314 + t318 * t468;
t470 = t311 * t316;
t174 = t201 * t314 - t318 * t470;
t249 = t333 * t311;
t221 = t249 * t314 - t312 * t318;
t566 = g(1) * t174 + g(2) * t170 + g(3) * t221;
t375 = mrSges(6,2) * t318 - mrSges(6,3) * t314;
t380 = mrSges(5,1) * t318 - mrSges(5,2) * t314;
t565 = mrSges(7,1) * t466 + mrSges(7,2) * t463 - t375 + t380;
t564 = -mrSges(4,1) - t565;
t563 = m(7) * (pkin(5) + pkin(9)) + t378 + t576;
t538 = -t208 / 0.2e1;
t542 = -t130 / 0.2e1;
t544 = -t129 / 0.2e1;
t562 = Ifges(7,5) * t542 + Ifges(7,6) * t544 + Ifges(7,3) * t538 - t569;
t8 = Ifges(7,4) * t42 + Ifges(7,2) * t43 + Ifges(7,6) * t105;
t561 = -t8 / 0.2e1;
t9 = Ifges(7,1) * t42 + Ifges(7,4) * t43 + Ifges(7,5) * t105;
t560 = t9 / 0.2e1;
t503 = Ifges(7,4) * t130;
t53 = Ifges(7,2) * t129 + Ifges(7,6) * t208 + t503;
t556 = -t53 / 0.2e1;
t555 = t53 / 0.2e1;
t493 = t210 * Ifges(5,4);
t95 = -t209 * Ifges(5,2) + t235 * Ifges(5,6) + t493;
t553 = -t95 / 0.2e1;
t543 = t129 / 0.2e1;
t541 = t130 / 0.2e1;
t537 = t208 / 0.2e1;
t536 = -t209 / 0.2e1;
t535 = t209 / 0.2e1;
t534 = -t210 / 0.2e1;
t533 = t210 / 0.2e1;
t531 = -t235 / 0.2e1;
t530 = t235 / 0.2e1;
t529 = -t241 / 0.2e1;
t528 = -t242 / 0.2e1;
t527 = t242 / 0.2e1;
t524 = t312 / 0.2e1;
t522 = pkin(1) * t311;
t518 = t221 * pkin(10);
t517 = qJD(2) / 0.2e1;
t513 = mrSges(4,1) * t248;
t510 = mrSges(5,3) * t209;
t508 = mrSges(7,3) * t313;
t506 = Ifges(4,4) * t242;
t505 = Ifges(5,4) * t314;
t504 = Ifges(5,4) * t318;
t502 = Ifges(7,4) * t313;
t501 = Ifges(7,4) * t317;
t500 = Ifges(6,6) * t314;
t499 = Ifges(6,6) * t318;
t492 = t210 * Ifges(6,6);
t491 = t241 * mrSges(4,3);
t481 = Ifges(3,6) * qJD(2);
t480 = qJ(5) * t209;
t477 = t210 * t317;
t465 = t313 * t318;
t461 = t315 * t320;
t459 = t316 * t319;
t458 = t317 * t318;
t224 = t304 + t339;
t238 = qJ(3) * t469 + t450;
t153 = t482 * t224 + t483 * t238;
t143 = pkin(9) * t312 + t153;
t301 = pkin(2) * t469;
t451 = -t248 * pkin(3) + t301;
t388 = pkin(9) * t249 + t451;
t167 = -t388 - t522;
t85 = t318 * t143 + t314 * t167;
t441 = qJD(6) * t317;
t432 = m(4) + t438;
t431 = -m(3) * pkin(1) - mrSges(2,1);
t420 = t319 * t449;
t419 = t311 * t447;
t417 = t307 * t443;
t416 = t307 * t442;
t413 = t514 * t311;
t407 = -t441 / 0.2e1;
t216 = t221 * pkin(4);
t222 = t249 * t318 + t312 * t314;
t401 = t222 * qJ(5) - t216;
t84 = -t314 * t143 + t167 * t318;
t81 = t157 * t318 - t136;
t171 = -t196 * t318 - t314 * t468;
t259 = pkin(2) * t467 - t413;
t400 = -t259 * t316 + t320 * t309;
t393 = mrSges(3,3) * t421;
t392 = mrSges(3,3) * t420;
t389 = t612 * pkin(2);
t77 = -qJ(5) * t248 - t85;
t384 = t201 * pkin(3) + t400;
t383 = -mrSges(7,3) - t572;
t381 = mrSges(5,1) * t221 + mrSges(5,2) * t222;
t159 = t221 * t317 - t248 * t313;
t160 = t221 * t313 + t248 * t317;
t379 = mrSges(7,1) * t159 - mrSges(7,2) * t160;
t376 = -t221 * mrSges(6,2) - t222 * mrSges(6,3);
t374 = Ifges(5,1) * t318 - t505;
t373 = Ifges(7,1) * t317 - t502;
t372 = Ifges(7,1) * t313 + t501;
t371 = -Ifges(5,2) * t314 + t504;
t369 = -Ifges(7,2) * t313 + t501;
t368 = Ifges(7,2) * t317 + t502;
t367 = Ifges(3,5) * t319 - Ifges(3,6) * t315;
t243 = qJD(2) * t249;
t244 = qJD(2) * t248;
t366 = -Ifges(4,5) * t244 - Ifges(4,6) * t243;
t364 = Ifges(7,5) * t317 - Ifges(7,6) * t313;
t363 = Ifges(7,5) * t313 + Ifges(7,6) * t317;
t362 = -Ifges(6,2) * t318 + t500;
t361 = Ifges(6,3) * t314 - t499;
t48 = pkin(5) * t222 - t248 * t551 - t84;
t152 = t224 * t483 - t482 * t238;
t142 = -t312 * pkin(3) - t152;
t83 = t142 - t401;
t63 = t83 + t518;
t21 = -t313 * t63 + t317 * t48;
t22 = t313 * t48 + t317 * t63;
t299 = t319 * t434;
t212 = t299 + (-qJD(2) * t412 + t446) * t311;
t213 = -t418 + (-t319 * t413 - t303) * qJD(2);
t119 = t212 * t482 - t483 * t213;
t354 = -qJ(5) * t603 + t602;
t353 = t197 * pkin(3) + t389;
t120 = t212 * t483 + t213 * t482;
t158 = pkin(2) * t419 + pkin(3) * t243 + pkin(9) * t244;
t33 = -t314 * t120 - t143 * t442 + t158 * t318 - t167 * t443;
t266 = -t312 * t459 - t461;
t32 = t318 * t120 - t143 * t443 + t314 * t158 + t167 * t442;
t337 = t266 * pkin(2);
t334 = -pkin(9) * t196 + t353;
t331 = t200 * pkin(3) + t337;
t151 = -t249 * t443 + (qJD(4) * t312 - t244) * t318;
t328 = -t151 * qJ(5) - t222 * qJD(5) + t119;
t327 = pkin(9) * t201 + t331;
t27 = -qJ(5) * t243 - qJD(5) * t248 - t32;
t297 = Ifges(3,4) * t420;
t290 = Ifges(3,3) * t300;
t289 = Ifges(4,3) * t300;
t276 = -t301 - t522;
t271 = -t318 * pkin(4) + t340;
t269 = -pkin(8) * t471 + t304;
t268 = (-mrSges(3,1) * t319 + mrSges(3,2) * t315) * t311;
t267 = -t312 * t462 + t457;
t265 = -t312 * t461 - t459;
t255 = -pkin(8) * t419 + t299;
t254 = t450 * qJD(1);
t253 = t298 - t396;
t252 = -mrSges(3,2) * t395 + t392;
t251 = mrSges(3,1) * t395 - t393;
t226 = Ifges(3,1) * t421 + Ifges(3,5) * t395 + t297;
t225 = t481 + (Ifges(3,6) * t312 + (t319 * Ifges(3,2) + t507) * t311) * qJD(1);
t214 = -mrSges(4,2) * t395 + t491;
t191 = t196 * pkin(3);
t182 = Ifges(4,6) * t189;
t181 = t332 * mrSges(4,2);
t176 = -mrSges(4,1) * t241 - mrSges(4,2) * t242;
t175 = t201 * t318 + t314 * t470;
t162 = t300 * mrSges(4,1) - mrSges(4,3) * t332;
t161 = -mrSges(4,2) * t300 + mrSges(4,3) * t189;
t156 = -Ifges(4,1) * t242 + t234 + t597;
t155 = -t506 + t596 + t616;
t150 = qJD(4) * t222 - t244 * t314;
t140 = -mrSges(5,2) * t235 - t510;
t118 = -mrSges(6,2) * t209 - mrSges(6,3) * t210;
t116 = pkin(4) * t210 + t480;
t100 = t174 * t313 - t200 * t317;
t99 = t174 * t317 + t200 * t313;
t92 = t235 * Ifges(6,5) + t209 * Ifges(6,3) - t492;
t88 = t210 * t551 + t480;
t78 = -pkin(4) * t248 - t84;
t76 = qJD(6) * t159 + t150 * t313 + t243 * t317;
t75 = -qJD(6) * t160 + t150 * t317 - t243 * t313;
t71 = mrSges(5,1) * t184 + mrSges(5,3) * t106;
t70 = pkin(4) * t242 - t81;
t56 = -pkin(5) * t221 - t77;
t51 = t68 - t519;
t47 = mrSges(5,1) * t107 - mrSges(5,2) * t106;
t46 = -mrSges(6,2) * t107 + mrSges(6,3) * t106;
t39 = t150 * pkin(4) + t328;
t31 = t150 * t551 + t328;
t30 = -pkin(4) * t243 - t33;
t29 = t313 * t51 + t317 * t88;
t28 = -t313 * t88 + t317 * t51;
t20 = -pkin(5) * t150 - t27;
t19 = pkin(5) * t151 - t243 * t551 - t33;
t18 = t107 * pkin(4) + t324;
t13 = -mrSges(7,1) * t43 + mrSges(7,2) * t42;
t6 = -pkin(5) * t107 - t10;
t4 = -qJD(6) * t22 + t19 * t317 - t31 * t313;
t3 = qJD(6) * t21 + t19 * t313 + t31 * t317;
t23 = [(t10 * mrSges(6,1) - t16 * mrSges(5,3) - Ifges(5,4) * t548 + Ifges(6,6) * t547 + t623) * t221 + (-t265 * mrSges(3,1) + t612 * mrSges(3,2) - t196 * mrSges(4,1) - m(5) * t191 + (t259 * t432 + mrSges(2,2)) * t320 + (t309 * t432 - t431) * t316 - t383 * t171 - (t354 - t377) * t170 + (-t378 - t574) * t197 - t603 * (-pkin(4) * t171 + t191)) * g(1) + (Ifges(7,5) * t76 + Ifges(7,6) * t75) * t537 + (Ifges(7,5) * t160 + Ifges(7,6) * t159) * t549 + (Ifges(7,4) * t76 + Ifges(7,2) * t75) * t543 + (Ifges(7,4) * t160 + Ifges(7,2) * t159) * t557 + (-t276 * mrSges(4,1) + Ifges(4,6) * t524) * t189 + (-mrSges(4,3) * t66 - Ifges(4,4) * t332 + Ifges(6,4) * t547 + Ifges(5,5) * t548 + Ifges(6,5) * t545 - Ifges(4,2) * t189 - Ifges(4,6) * t300 + Ifges(5,6) * t546 - t539 * t618 + t606 + t613 / 0.2e1) * t248 + (t269 * mrSges(3,1) - t450 * mrSges(3,2) + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t312) * t300 + (t257 * t450 - t258 * t269) * mrSges(3,3) + m(3) * (t205 * t450 + t206 * t269 - t253 * t256 + t254 * t255) + (mrSges(6,1) * t12 - mrSges(5,3) * t17 - Ifges(6,6) * t545 + t621 + t627 + t631) * t222 + (t55 * mrSges(6,1) + t67 * mrSges(5,3) + Ifges(5,1) * t533 + Ifges(5,4) * t536 + Ifges(7,5) * t541 - Ifges(6,2) * t534 - Ifges(6,6) * t535 + Ifges(7,6) * t543 + Ifges(7,3) * t537 - t530 * t601 + t608 + t624) * t151 + (t617 / 0.2e1 - t125 * mrSges(4,3) - Ifges(4,4) * t528 + Ifges(5,5) * t533 + Ifges(6,4) * t534 + Ifges(6,5) * t535 + Ifges(5,6) * t536 - t155 / 0.2e1 - t616 / 0.2e1 - t618 * t530 + t605) * t243 + m(4) * (-t119 * t124 + t120 * t125 + t152 * t65 + t153 * t66 + t223 * t276) + t455 * t119 + (t92 / 0.2e1 - Ifges(5,4) * t533 + Ifges(6,6) * t534 + Ifges(6,3) * t535 - Ifges(5,2) * t536 + t553 + t58 * mrSges(6,1) - t68 * mrSges(5,3) - t600 * t530 + t609) * t150 + t276 * t181 + (Ifges(7,1) * t76 + Ifges(7,4) * t75) * t541 + (Ifges(7,1) * t160 + Ifges(7,4) * t159) * t558 + t39 * t118 + t44 * (-mrSges(7,1) * t75 + mrSges(7,2) * t76) + t76 * t54 / 0.2e1 + t77 * t73 + (t124 * mrSges(4,3) - t615 - Ifges(4,1) * t528 - t234 / 0.2e1 - t156 / 0.2e1) * t244 + t3 * t90 + t4 * t91 + t84 * t71 + t85 * t72 + t78 * t74 + t20 * t80 + t83 * t46 + m(7) * (t1 * t22 + t14 * t4 + t15 * t3 + t2 * t21 + t20 * t44 + t56 * t6) + m(6) * (t10 * t77 + t12 * t78 + t18 * t83 + t27 * t58 + t30 * t55 + t39 * t64) + m(5) * (t114 * t119 + t142 * t61 + t16 * t85 + t17 * t84 + t32 * t68 - t33 * t67) + t56 * t13 + t22 * t26 + t21 * t25 + Ifges(2,3) * qJDD(1) - t6 * t379 + t61 * t381 + t18 * t376 + (-t267 * mrSges(3,1) - t266 * mrSges(3,2) - m(4) * t400 - t201 * mrSges(4,1) + mrSges(2,2) * t316 - m(5) * t384 - t100 * mrSges(7,1) - t99 * mrSges(7,2) + t431 * t320 + t354 * t174 + t383 * t175 + t574 * t200 + t603 * (-t175 * pkin(4) - t384)) * g(2) - t256 * t251 + t255 * t252 + t27 * t138 + t30 * t139 + t32 * t140 + t33 * t141 + t142 * t47 + (t223 * mrSges(4,2) - t65 * mrSges(4,3) + Ifges(4,1) * t332 + Ifges(4,4) * t189 + Ifges(4,5) * t300) * t249 + (t1 * t159 - t14 * t76 + t15 * t75 - t160 * t2) * mrSges(7,3) + t159 * t8 / 0.2e1 + t153 * t161 + t152 * t162 + (t182 / 0.2e1 + t289 / 0.2e1 + t290 / 0.2e1 + qJD(1) * t366 / 0.2e1 + t571) * t312 + t223 * t513 + t366 * t517 + ((mrSges(3,1) * t257 - mrSges(3,2) * t258 + (m(3) * t522 - t268) * qJDD(1)) * pkin(1) + (mrSges(3,3) * t205 + Ifges(3,4) * t258 + Ifges(3,2) * t257 + Ifges(3,6) * t300) * t319 + (-mrSges(3,3) * t206 + Ifges(3,1) * t258 + Ifges(3,4) * t257 + Ifges(3,5) * t300) * t315 + ((t226 / 0.2e1 - t253 * mrSges(3,3) + Ifges(3,5) * t517) * t319 + (-t254 * mrSges(3,3) - t225 / 0.2e1 - t481 / 0.2e1 + (m(4) * t263 + t176) * pkin(2)) * t315 + (t367 * t524 + (t319 * (Ifges(3,4) * t319 - Ifges(3,2) * t315) / 0.2e1 - t577) * t311) * qJD(1)) * qJD(2) + (g(1) * t320 + g(2) * t316) * (-m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3))) * t311 + t120 * t214 + t75 * t555 + t160 * t560; (t124 * t144 - t125 * t145 - t263 * t397 + (t482 * t66 + t483 * t65) * pkin(2)) * m(4) + (-m(7) * (pkin(10) * t478 + t331 + t587) - mrSges(7,3) * t478 - m(5) * t327 - m(6) * (t327 + t587) - mrSges(3,1) * t266 + mrSges(3,2) * t267 - m(4) * t337 + t564 * t200 - t563 * t201) * g(1) + (t374 / 0.2e1 - t362 / 0.2e1) * t444 + (-t363 * t549 - t368 * t557 - t372 * t558 + t6 * t378 + t54 * t407 - t623) * t318 + (-m(7) * (pkin(10) * t479 + t353 + t588) - mrSges(7,3) * t479 - m(4) * t389 - m(5) * t334 - m(6) * (t334 + t588) - mrSges(3,1) * t612 - mrSges(3,2) * t265 + t564 * t197 + t563 * t196) * g(2) + (t562 - t629) * t475 + (Ifges(7,5) * t169 + Ifges(7,6) * t168) * t538 + (-t417 - t82) * t140 + (t95 / 0.2e1 - t92 / 0.2e1) * t476 - t499 * t547 + t504 * t548 + (t393 + t251) * t254 + (Ifges(7,1) * t169 + Ifges(7,4) * t168) * t542 + (t392 - t252) * t253 - ((-Ifges(3,2) * t421 + t226 + t297) * t319 + t395 * t367) * t449 / 0.2e1 + ((-t371 / 0.2e1 + t361 / 0.2e1) * t209 + t584 * t530 + (Ifges(7,3) * t318 + t314 * t363) * t537 + (Ifges(7,5) * t318 + t314 * t372) * t541 + (Ifges(7,6) * t318 + t314 * t368) * t543 + t585) * qJD(4) + t621 * t314 + (t417 - t69) * t138 + (t416 - t70) * t139 + (-t364 * t537 - t369 * t543 - t373 * t541) * t440 + (-t114 * t144 + t308 * t61 + t67 * t81 - t68 * t82) * m(5) + (t18 * t271 - t55 * t70 - t58 * t69 + t591 * t64) * m(6) + ((-t71 + t74) * t314 + t599 * t318 + ((-t314 * t68 + t318 * t67) * qJD(4) + t579) * m(5) + ((t314 * t58 + t318 * t55) * qJD(4) + t578) * m(6)) * t307 - t500 * t545 + t505 * t546 + (t234 + t156) * t529 + (t317 * t53 + t488 + t92) * t443 / 0.2e1 + t624 * t442 - t9 * t465 / 0.2e1 + (Ifges(7,4) * t169 + Ifges(7,2) * t168) * t544 + t225 * t421 / 0.2e1 + (-t416 - t81) * t141 + (t506 + t617) * t527 + (-t596 / 0.2e1 + Ifges(4,2) * t529 - Ifges(6,4) * t533 - Ifges(5,5) * t534 - Ifges(5,6) * t535 - Ifges(6,5) * t536 + t618 * t531 + t605) * t242 - t125 * t490 + t571 + t182 + (-t597 / 0.2e1 - t615 + Ifges(4,1) * t527 + t362 * t533 + t374 * t534 + t371 * t535 + t361 * t536 + t584 * t531 - t585) * t241 - t176 * t397 + t290 + t289 - t61 * t380 + t18 * t375 + t308 * t47 + t271 * t46 + t273 * t13 + t314 * t627 + t161 * t422 + t162 * t423 + t577 * qJD(1) ^ 2 * t311 ^ 2 - t169 * t54 / 0.2e1 + t194 * t25 + t124 * t491 + t155 * t528 + t195 * t26 - t145 * t214 + (t55 * t583 + t58 * t582 + t578) * mrSges(6,1) + (-t582 * t68 + t583 * t67 + t579) * mrSges(5,3) + (-m(4) * t301 + t513 - m(7) * (-pkin(10) * t474 + t451 + t586) + mrSges(7,3) * t474 - m(6) * (t388 + t586) - m(5) * t388 + t268 + t565 * t248 - t563 * t249) * g(3) + (mrSges(7,1) * t590 - mrSges(7,2) * t589) * t44 + (-t1 * t458 + t14 * t589 - t15 * t590 + t2 * t465) * mrSges(7,3) - t455 * t144 + t443 * t553 + t415 * t555 + t168 * t556 + t458 * t561 + t591 * t118 + t593 * t80 + t594 * t91 + t595 * t90 + (t1 * t195 + t14 * t594 + t15 * t595 + t194 * t2 + t273 * t6 + t44 * t593) * m(7); -t189 * mrSges(4,1) - t168 * t91 - t169 * t90 - t241 * t214 + t181 + (t118 + t455) * t242 + (t13 + t454 * t241 + (t313 * t90 + t317 * t91 - t454) * qJD(4) + t599) * t314 + (t71 - t235 * (-t140 - t592) + t568) * t318 + ((t6 + (t14 * t317 + t15 * t313) * qJD(4)) * t314 + (qJD(4) * t44 - t326) * t318 - t14 * t168 - t15 * t169 - t44 * t475) * m(7) + (-t10 * t314 - t12 * t318 + t242 * t64 + t235 * (t314 * t55 - t318 * t58)) * m(6) + (t114 * t242 + t16 * t314 + t17 * t318 + t235 * (t314 * t67 + t318 * t68)) * m(5) + (-t124 * t242 - t125 * t241 + t223) * m(4) + (-t312 * g(3) + (-g(1) * t316 + g(2) * t320) * t311) * t432; (-m(7) * t326 + t622) * t551 + (qJD(6) * t496 - t520 + (-t441 - t477) * t15 + t566) * mrSges(7,3) + (t203 + t94) * t536 + (t492 + t95) * t533 + (-t493 + t92) * t534 + (-m(7) * (-t216 - t518) - t338 * t222 - m(6) * t401 + t376 + t381) * g(3) + t350 * t80 + (-pkin(4) * t12 - qJ(5) * t10 - qJD(5) * t58 - t116 * t64) * m(6) + (-m(6) * t58 - t138 + t140 + t510) * t67 + (-t73 + t13) * qJ(5) - (t129 * t368 + t130 * t372 + t208 * t363) * qJD(6) / 0.2e1 - t1 * t508 - t116 * t118 - pkin(4) * t74 - t29 * t90 - t28 * t91 - t58 * t511 + (t174 * t611 + t175 * t567) * g(1) + (t170 * t611 + t171 * t567) * g(2) + (t6 * qJ(5) - t14 * t28 - t15 * t29 + t44 * t610) * m(7) + (-Ifges(5,1) * t534 + Ifges(6,2) * t533 + t531 * t601 - t562 + t608) * t209 + (Ifges(6,3) * t536 + t14 * t508 + t363 * t538 + t368 * t544 + t372 * t542 - t531 * t600 + t573 - t609) * t210 + t606 + t6 * t377 + t613 + (-Ifges(5,2) * t210 - t204 + t598) * t535 + t53 * t407 + t55 * t512 + t573 * qJD(6) + (-m(6) * t55 + t454 + t509) * t68 + t364 * t549 + t477 * t556 + t369 * t557 + t373 * t558 + t317 * t560 + t313 * t561 + t592 * qJD(5); -t592 * t235 + (t118 + t580) * t210 + (-t210 * t359 - t235 * t44 + t326 - t566) * m(7) + (t210 * t64 + t235 * t58 + t12 - t566) * m(6) - t568; -t44 * (mrSges(7,1) * t130 + mrSges(7,2) * t129) + (Ifges(7,1) * t129 - t503) * t542 + t53 * t541 + (Ifges(7,5) * t129 - Ifges(7,6) * t130) * t538 - t14 * t90 + t15 * t91 - g(1) * (mrSges(7,1) * t99 - mrSges(7,2) * t100) - g(2) * ((t170 * t317 + t197 * t313) * mrSges(7,1) + (-t170 * t313 + t197 * t317) * mrSges(7,2)) - g(3) * t379 + (t129 * t14 + t130 * t15) * mrSges(7,3) + t7 + (-Ifges(7,2) * t130 + t126 + t54) * t544 + t570;];
tau  = t23;
