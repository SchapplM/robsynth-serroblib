% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:21:21
% EndTime: 2019-03-09 10:22:31
% DurationCPUTime: 42.48s
% Computational Cost: add. (24789->966), mult. (68543->1313), div. (0->0), fcn. (56641->16), ass. (0->439)
t359 = cos(pkin(6));
t363 = sin(qJ(2));
t486 = t359 * t363;
t342 = pkin(1) * t486;
t357 = sin(pkin(6));
t367 = cos(qJ(2));
t488 = t357 * t367;
t525 = pkin(8) + qJ(3);
t266 = (t488 * t525 + t342) * qJD(1);
t502 = sin(pkin(11));
t256 = t502 * t266;
t535 = pkin(1) * t359;
t344 = t367 * t535;
t336 = qJD(1) * t344;
t435 = t525 * t363;
t415 = t357 * t435;
t265 = -qJD(1) * t415 + t336;
t503 = cos(pkin(11));
t179 = t265 * t503 - t256;
t310 = t363 * t502 - t367 * t503;
t606 = t357 * t310;
t280 = qJD(1) * t606;
t375 = t363 * t503 + t367 * t502;
t472 = qJD(1) * t357;
t281 = t375 * t472;
t441 = t363 * t472;
t420 = pkin(2) * t441;
t202 = pkin(3) * t281 + pkin(9) * t280 + t420;
t362 = sin(qJ(4));
t366 = cos(qJ(4));
t117 = -t179 * t362 + t366 * t202;
t442 = t502 * pkin(2);
t348 = t442 + pkin(9);
t477 = qJ(5) + t348;
t423 = qJD(4) * t477;
t495 = t280 * t366;
t643 = -pkin(4) * t281 - qJ(5) * t495 - t362 * qJD(5) - t366 * t423 - t117;
t118 = t366 * t179 + t362 * t202;
t496 = t280 * t362;
t642 = qJ(5) * t496 - qJD(5) * t366 + t362 * t423 + t118;
t468 = qJD(4) * t362;
t641 = t468 + t496;
t356 = sin(pkin(12));
t358 = cos(pkin(12));
t309 = t356 * t362 - t358 * t366;
t191 = t309 * t280;
t300 = t309 * qJD(4);
t640 = t191 + t300;
t311 = t356 * t366 + t358 * t362;
t190 = t311 * t280;
t299 = t311 * qJD(4);
t601 = t299 + t190;
t273 = t310 * t472 + qJD(4);
t545 = -t273 / 0.2e1;
t339 = qJD(1) * t359 + qJD(2);
t245 = -t281 * t362 + t339 * t366;
t393 = -t281 * t366 - t339 * t362;
t394 = t245 * t356 - t358 * t393;
t551 = -t394 / 0.2e1;
t425 = t358 * t245 + t356 * t393;
t553 = -t425 / 0.2e1;
t639 = -Ifges(6,4) * t551 - Ifges(6,2) * t553 - Ifges(6,6) * t545;
t149 = qJD(6) - t425;
t554 = t149 / 0.2e1;
t361 = sin(qJ(6));
t365 = cos(qJ(6));
t129 = t273 * t361 + t365 * t394;
t560 = t129 / 0.2e1;
t128 = t273 * t365 - t361 * t394;
t562 = t128 / 0.2e1;
t638 = Ifges(7,5) * t560 + Ifges(7,6) * t562 + Ifges(7,3) * t554;
t624 = Ifges(5,3) + Ifges(6,3);
t406 = -mrSges(7,1) * t365 + mrSges(7,2) * t361;
t611 = m(7) * pkin(5) + mrSges(6,1) - t406;
t421 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t608 = t356 * t643 - t642 * t358;
t427 = t503 * t266;
t178 = t265 * t502 + t427;
t605 = pkin(4) * t641 - t178;
t247 = pkin(2) * t339 + t265;
t164 = t247 * t503 - t256;
t156 = -t339 * pkin(3) - t164;
t125 = -t245 * pkin(4) + qJD(5) + t156;
t165 = t502 * t247 + t427;
t157 = pkin(9) * t339 + t165;
t352 = pkin(2) * t367 + pkin(1);
t298 = -t352 * t472 + qJD(3);
t175 = pkin(3) * t280 - pkin(9) * t281 + t298;
t112 = -t157 * t362 + t366 * t175;
t94 = qJ(5) * t393 + t112;
t76 = pkin(4) * t273 + t94;
t113 = t157 * t366 + t175 * t362;
t95 = qJ(5) * t245 + t113;
t89 = t358 * t95;
t39 = t356 * t76 + t89;
t37 = pkin(10) * t273 + t39;
t63 = -pkin(5) * t425 - pkin(10) * t394 + t125;
t17 = -t361 * t37 + t365 * t63;
t18 = t361 * t63 + t365 * t37;
t637 = mrSges(6,1) * t125 + mrSges(7,1) * t17 - mrSges(7,2) * t18 - mrSges(6,3) * t39 + t638 - t639;
t464 = qJD(1) * qJD(2);
t295 = (qJDD(1) * t367 - t363 * t464) * t357;
t296 = (qJDD(1) * t363 + t367 * t464) * t357;
t222 = t295 * t503 - t296 * t502;
t221 = qJDD(4) - t222;
t549 = t221 / 0.2e1;
t223 = t295 * t502 + t296 * t503;
t462 = qJDD(1) * t359;
t338 = qJDD(2) + t462;
t139 = qJD(4) * t245 + t223 * t366 + t338 * t362;
t140 = qJD(4) * t393 - t223 * t362 + t338 * t366;
t85 = t139 * t358 + t140 * t356;
t568 = t85 / 0.2e1;
t84 = -t139 * t356 + t140 * t358;
t569 = t84 / 0.2e1;
t473 = pkin(8) * t488 + t342;
t294 = t473 * qJD(2);
t455 = pkin(1) * t462;
t334 = t367 * t455;
t490 = t357 * t363;
t438 = qJD(3) * t490;
t463 = qJDD(1) * t357;
t454 = pkin(8) * t463;
t161 = -t363 * t454 + pkin(2) * t338 - qJ(3) * t296 + t334 + (-t294 - t438) * qJD(1);
t460 = qJD(2) * t535;
t416 = qJD(1) * t460;
t444 = t363 * t455 + (t416 + t454) * t367;
t469 = qJD(3) * t367;
t470 = qJD(2) * t363;
t170 = qJ(3) * t295 + (-pkin(8) * t470 + t469) * t472 + t444;
t110 = t161 * t503 - t502 * t170;
t104 = -t338 * pkin(3) - t110;
t66 = -t140 * pkin(4) + qJDD(5) + t104;
t111 = t502 * t161 + t503 * t170;
t105 = pkin(9) * t338 + t111;
t261 = -pkin(1) * t463 - pkin(2) * t295 + qJDD(3);
t126 = -pkin(3) * t222 - pkin(9) * t223 + t261;
t35 = -qJD(4) * t113 - t105 * t362 + t366 * t126;
t21 = pkin(4) * t221 - qJ(5) * t139 + qJD(5) * t393 + t35;
t467 = qJD(4) * t366;
t34 = t366 * t105 + t362 * t126 - t157 * t468 + t175 * t467;
t29 = qJ(5) * t140 + qJD(5) * t245 + t34;
t7 = t21 * t358 - t29 * t356;
t636 = t66 * mrSges(6,2) - t7 * mrSges(6,3) + 0.2e1 * Ifges(6,1) * t568 + 0.2e1 * Ifges(6,4) * t569 + 0.2e1 * Ifges(6,5) * t549;
t634 = t112 * mrSges(5,1);
t633 = t113 * mrSges(5,2);
t632 = -pkin(10) * t281 + t608;
t631 = t601 * pkin(5) + pkin(10) * t640 + t605;
t83 = qJDD(6) - t84;
t570 = t83 / 0.2e1;
t49 = -qJD(6) * t129 + t221 * t365 - t361 * t85;
t578 = t49 / 0.2e1;
t48 = qJD(6) * t128 + t221 * t361 + t365 * t85;
t579 = t48 / 0.2e1;
t19 = -t84 * pkin(5) - t85 * pkin(10) + t66;
t8 = t356 * t21 + t358 * t29;
t6 = pkin(10) * t221 + t8;
t1 = qJD(6) * t17 + t19 * t361 + t365 * t6;
t2 = -qJD(6) * t18 + t19 * t365 - t361 * t6;
t588 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t9 = Ifges(7,5) * t48 + Ifges(7,6) * t49 + Ifges(7,3) * t83;
t629 = t588 + mrSges(6,1) * t66 - mrSges(6,3) * t8 + Ifges(7,5) * t579 + Ifges(7,6) * t578 + Ifges(7,3) * t570 + t9 / 0.2e1 + (-t549 - t221 / 0.2e1) * Ifges(6,6) + (-t569 - t84 / 0.2e1) * Ifges(6,2) + (-t568 - t85 / 0.2e1) * Ifges(6,4);
t615 = m(7) + m(6);
t616 = -m(5) - m(4);
t458 = t615 - t616;
t544 = t273 / 0.2e1;
t550 = t394 / 0.2e1;
t552 = t425 / 0.2e1;
t627 = -Ifges(6,4) * t550 - Ifges(6,2) * t552 - Ifges(6,6) * t544 + t637 + t638;
t625 = t35 * mrSges(5,1) + t7 * mrSges(6,1) - t34 * mrSges(5,2) - t8 * mrSges(6,2);
t239 = Ifges(5,4) * t245;
t136 = -Ifges(5,1) * t393 + t273 * Ifges(5,5) + t239;
t558 = t136 / 0.2e1;
t287 = t375 * t357;
t259 = -t287 * t362 + t359 * t366;
t411 = t259 * pkin(4);
t623 = -Ifges(5,5) * t393 + Ifges(6,5) * t394 + t245 * Ifges(5,6) + Ifges(6,6) * t425 + t273 * t624;
t610 = t642 * t356 + t358 * t643;
t355 = qJ(4) + pkin(12);
t353 = sin(t355);
t354 = cos(t355);
t408 = -mrSges(5,1) * t366 + mrSges(5,2) * t362;
t380 = m(5) * pkin(3) - t408;
t622 = -t353 * t421 + t354 * t611 + t380;
t368 = cos(qJ(1));
t479 = t367 * t368;
t364 = sin(qJ(1));
t481 = t364 * t363;
t621 = t359 * t479 - t481;
t474 = t375 * t359;
t237 = -t368 * t310 - t364 * t474;
t489 = t357 * t364;
t207 = -t237 * t362 + t366 * t489;
t620 = Ifges(5,5) * t139 + Ifges(6,5) * t85 + Ifges(5,6) * t140 + Ifges(6,6) * t84 + t221 * t624;
t555 = -t149 / 0.2e1;
t561 = -t129 / 0.2e1;
t563 = -t128 / 0.2e1;
t619 = Ifges(7,5) * t561 + Ifges(7,6) * t563 + Ifges(7,3) * t555 - t637 + t639;
t557 = t139 / 0.2e1;
t556 = t140 / 0.2e1;
t16 = -mrSges(7,1) * t49 + mrSges(7,2) * t48;
t68 = mrSges(6,1) * t221 - mrSges(6,3) * t85;
t614 = -t68 + t16;
t443 = t503 * pkin(2);
t350 = -t443 - pkin(3);
t529 = t366 * pkin(4);
t319 = t350 - t529;
t231 = t309 * pkin(5) - t311 * pkin(10) + t319;
t308 = t477 * t366;
t426 = t362 * t477;
t244 = t358 * t308 - t356 * t426;
t147 = t231 * t365 - t244 * t361;
t613 = qJD(6) * t147 + t361 * t631 + t365 * t632;
t148 = t231 * t361 + t244 * t365;
t612 = -qJD(6) * t148 - t361 * t632 + t365 * t631;
t609 = pkin(5) * t281 - t610;
t131 = mrSges(6,1) * t273 - mrSges(6,3) * t394;
t74 = -mrSges(7,1) * t128 + mrSges(7,2) * t129;
t607 = t74 - t131;
t145 = -t191 * t361 + t281 * t365;
t465 = qJD(6) * t365;
t388 = -t361 * t300 + t311 * t465;
t604 = t145 + t388;
t146 = t191 * t365 + t281 * t361;
t466 = qJD(6) * t361;
t387 = t365 * t300 + t311 * t466;
t603 = t146 + t387;
t262 = pkin(2) * t359 + t344 - t415;
t275 = qJ(3) * t488 + t473;
t189 = t502 * t262 + t503 * t275;
t177 = pkin(9) * t359 + t189;
t340 = pkin(2) * t488;
t432 = pkin(9) * t287 + t340;
t536 = pkin(1) * t357;
t206 = pkin(3) * t606 - t432 - t536;
t120 = t366 * t177 + t362 * t206;
t521 = mrSges(4,3) * t281;
t602 = mrSges(4,1) * t339 + mrSges(5,1) * t245 + mrSges(5,2) * t393 - t521;
t260 = t287 * t366 + t359 * t362;
t168 = t259 * t356 + t260 * t358;
t278 = t606 * t365;
t141 = -t168 * t361 + t278;
t405 = mrSges(7,1) * t361 + mrSges(7,2) * t365;
t592 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t410 = m(5) * pkin(9) - t592;
t599 = -t405 - t410;
t480 = t364 * t367;
t483 = t363 * t368;
t303 = -t359 * t480 - t483;
t232 = t364 * t310 - t368 * t474;
t487 = t357 * t368;
t598 = t232 * t362 - t366 * t487;
t27 = mrSges(7,1) * t83 - mrSges(7,3) * t48;
t28 = -mrSges(7,2) * t83 + mrSges(7,3) * t49;
t597 = -t361 * t27 + t365 * t28;
t596 = t34 * t366 - t35 * t362;
t595 = t1 * t365 - t2 * t361;
t594 = t280 + qJD(4);
t519 = Ifges(3,4) * t363;
t593 = -t363 * (Ifges(3,1) * t367 - t519) / 0.2e1 + pkin(1) * (mrSges(3,1) * t363 + mrSges(3,2) * t367);
t506 = t356 * t95;
t38 = t358 * t76 - t506;
t591 = t125 * mrSges(6,2) - t38 * mrSges(6,3);
t590 = -mrSges(4,1) - t622;
t439 = t357 * t470;
t419 = pkin(8) * t439;
t240 = -qJD(1) * t419 + t444;
t241 = -pkin(8) * t296 - t363 * t416 + t334;
t589 = t241 * mrSges(3,1) + t110 * mrSges(4,1) - t240 * mrSges(3,2) - t111 * mrSges(4,2) + Ifges(3,5) * t296 + Ifges(4,5) * t223 + Ifges(3,6) * t295;
t583 = m(6) * pkin(4);
t582 = Ifges(7,1) * t579 + Ifges(7,4) * t578 + Ifges(7,5) * t570;
t516 = Ifges(7,4) * t129;
t61 = t128 * Ifges(7,2) + t149 * Ifges(7,6) + t516;
t574 = t61 / 0.2e1;
t127 = Ifges(7,4) * t128;
t62 = t129 * Ifges(7,1) + t149 * Ifges(7,5) + t127;
t573 = -t62 / 0.2e1;
t572 = Ifges(5,4) * t557 + Ifges(5,2) * t556 + Ifges(5,6) * t549;
t571 = Ifges(5,1) * t557 + Ifges(5,4) * t556 + Ifges(5,5) * t549;
t93 = Ifges(6,1) * t394 + Ifges(6,4) * t425 + t273 * Ifges(6,5);
t565 = -t93 / 0.2e1;
t564 = t93 / 0.2e1;
t512 = t393 * Ifges(5,4);
t135 = t245 * Ifges(5,2) + t273 * Ifges(5,6) - t512;
t559 = t135 / 0.2e1;
t548 = -t245 / 0.2e1;
t547 = t393 / 0.2e1;
t546 = -t393 / 0.2e1;
t542 = t281 / 0.2e1;
t538 = t365 / 0.2e1;
t534 = pkin(4) * t393;
t533 = pkin(4) * t356;
t532 = pkin(4) * t358;
t283 = qJD(2) * t606;
t187 = qJD(4) * t259 - t283 * t366;
t282 = qJD(2) * t287;
t337 = t367 * t460;
t248 = t337 + (-qJD(2) * t435 + t469) * t357;
t436 = t525 * t357;
t249 = -t438 + (-t367 * t436 - t342) * qJD(2);
t160 = t248 * t503 + t249 * t502;
t203 = pkin(2) * t439 + pkin(3) * t282 + pkin(9) * t283;
t65 = -qJD(4) * t120 - t160 * t362 + t366 * t203;
t45 = pkin(4) * t282 - qJ(5) * t187 - qJD(5) * t260 + t65;
t186 = -qJD(4) * t260 + t283 * t362;
t64 = t366 * t160 - t177 * t468 + t362 * t203 + t206 * t467;
t55 = qJ(5) * t186 + qJD(5) * t259 + t64;
t15 = t356 * t45 + t358 * t55;
t524 = mrSges(4,1) * t606;
t523 = mrSges(5,2) * t366;
t522 = mrSges(4,3) * t280;
t520 = mrSges(5,3) * t393;
t518 = Ifges(5,4) * t362;
t517 = Ifges(5,4) * t366;
t515 = Ifges(7,4) * t361;
t514 = Ifges(7,4) * t365;
t513 = t112 * mrSges(5,3);
t511 = t281 * Ifges(4,4);
t510 = t339 * Ifges(3,5);
t509 = t339 * Ifges(3,6);
t119 = -t177 * t362 + t366 * t206;
t100 = pkin(4) * t606 - qJ(5) * t260 + t119;
t108 = qJ(5) * t259 + t120;
t57 = t356 * t100 + t358 * t108;
t501 = t425 * t361;
t500 = t425 * t365;
t494 = t606 * t361;
t492 = t311 * t361;
t491 = t311 * t365;
t457 = -m(3) * pkin(1) - mrSges(2,1);
t453 = t362 * t489;
t331 = t362 * t487;
t448 = t62 * t538;
t440 = t367 * t472;
t40 = -t84 * mrSges(6,1) + t85 * mrSges(6,2);
t433 = -t466 / 0.2e1;
t193 = -t232 * t354 - t353 * t487;
t192 = t232 * t353 - t354 * t487;
t297 = pkin(2) * t486 - t436;
t424 = -t297 * t364 + t368 * t352;
t422 = -m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3);
t418 = mrSges(3,3) * t441;
t417 = mrSges(3,3) * t440;
t413 = t621 * pkin(2);
t412 = t207 * pkin(4);
t409 = mrSges(5,1) * t259 - mrSges(5,2) * t260;
t404 = Ifges(5,1) * t366 - t518;
t403 = Ifges(7,1) * t365 - t515;
t402 = -Ifges(5,2) * t362 + t517;
t401 = -Ifges(7,2) * t361 + t514;
t400 = Ifges(5,5) * t366 - Ifges(5,6) * t362;
t399 = Ifges(7,5) * t365 - Ifges(7,6) * t361;
t398 = -t17 * t361 + t18 * t365;
t14 = -t356 * t55 + t358 * t45;
t52 = pkin(10) * t606 + t57;
t188 = t262 * t503 - t502 * t275;
t176 = -t359 * pkin(3) - t188;
t143 = t176 - t411;
t167 = -t358 * t259 + t260 * t356;
t73 = t167 * pkin(5) - t168 * pkin(10) + t143;
t23 = t361 * t73 + t365 * t52;
t22 = -t361 * t52 + t365 * t73;
t86 = -mrSges(7,2) * t149 + mrSges(7,3) * t128;
t87 = mrSges(7,1) * t149 - mrSges(7,3) * t129;
t397 = -t361 * t87 + t365 * t86;
t159 = t248 * t502 - t503 * t249;
t56 = t100 * t358 - t108 * t356;
t142 = t168 * t365 + t494;
t173 = -mrSges(5,2) * t273 + mrSges(5,3) * t245;
t174 = mrSges(5,1) * t273 + t520;
t395 = t173 * t366 - t174 * t362;
t130 = -mrSges(6,2) * t273 + mrSges(6,3) * t425;
t392 = -t130 - t397;
t36 = -pkin(5) * t273 - t38;
t389 = t36 * t405;
t373 = t359 * t310;
t236 = t364 * t373 - t368 * t375;
t351 = pkin(3) + t529;
t360 = -qJ(5) - pkin(9);
t382 = pkin(4) * t453 + t236 * t360 + t237 * t351 + t424;
t233 = -t364 * t375 - t368 * t373;
t378 = g(1) * t236 + g(2) * t233 - g(3) * t606;
t116 = -t186 * pkin(4) + t159;
t370 = (-t17 * t365 - t18 * t361) * qJD(6) + t595;
t349 = -pkin(5) - t532;
t335 = Ifges(3,4) * t440;
t329 = Ifges(3,3) * t338;
t328 = Ifges(4,3) * t338;
t313 = -t340 - t536;
t306 = -pkin(8) * t490 + t344;
t305 = (-mrSges(3,1) * t367 + mrSges(3,2) * t363) * t357;
t304 = -t359 * t481 + t479;
t302 = -t359 * t483 - t480;
t293 = t337 - t419;
t292 = t473 * qJD(1);
t291 = -pkin(8) * t441 + t336;
t289 = -mrSges(3,2) * t339 + t417;
t288 = mrSges(3,1) * t339 - t418;
t272 = Ifges(4,4) * t280;
t264 = Ifges(3,1) * t441 + t335 + t510;
t263 = t509 + (t367 * Ifges(3,2) + t519) * t472;
t252 = -mrSges(4,2) * t339 - t522;
t251 = t287 * t354 + t353 * t359;
t243 = t308 * t356 + t358 * t426;
t219 = Ifges(4,6) * t222;
t218 = t223 * mrSges(4,2);
t211 = mrSges(4,1) * t280 + mrSges(4,2) * t281;
t208 = t237 * t366 + t453;
t205 = mrSges(4,1) * t338 - mrSges(4,3) * t223;
t204 = -mrSges(4,2) * t338 + mrSges(4,3) * t222;
t201 = t281 * Ifges(4,1) + t339 * Ifges(4,5) - t272;
t200 = -Ifges(4,2) * t280 + t339 * Ifges(4,6) + t511;
t197 = t237 * t354 + t353 * t489;
t196 = t237 * t353 - t354 * t489;
t133 = t197 * t365 - t236 * t361;
t132 = -t197 * t361 - t236 * t365;
t122 = t186 * t356 + t187 * t358;
t121 = -t358 * t186 + t187 * t356;
t115 = -mrSges(5,2) * t221 + mrSges(5,3) * t140;
t114 = mrSges(5,1) * t221 - mrSges(5,3) * t139;
t102 = -mrSges(6,1) * t425 + mrSges(6,2) * t394;
t88 = -mrSges(5,1) * t140 + mrSges(5,2) * t139;
t79 = pkin(5) * t394 - pkin(10) * t425 - t534;
t70 = -qJD(6) * t142 - t122 * t361 + t282 * t365;
t69 = qJD(6) * t141 + t122 * t365 + t282 * t361;
t67 = -mrSges(6,2) * t221 + mrSges(6,3) * t84;
t51 = -pkin(5) * t606 - t56;
t43 = t358 * t94 - t506;
t42 = t356 * t94 + t89;
t41 = t121 * pkin(5) - t122 * pkin(10) + t116;
t25 = t361 * t79 + t365 * t43;
t24 = -t361 * t43 + t365 * t79;
t13 = pkin(10) * t282 + t15;
t12 = -pkin(5) * t282 - t14;
t10 = t48 * Ifges(7,4) + t49 * Ifges(7,2) + t83 * Ifges(7,6);
t5 = -pkin(5) * t221 - t7;
t4 = -qJD(6) * t23 - t13 * t361 + t365 * t41;
t3 = qJD(6) * t22 + t13 * t365 + t361 * t41;
t11 = [(t295 * t473 - t296 * t306) * mrSges(3,3) + (Ifges(5,1) * t260 + Ifges(5,4) * t259) * t557 + t282 * t634 + m(4) * (t110 * t188 + t111 * t189 - t159 * t164 + t160 * t165 + t261 * t313) - t602 * t159 + (t306 * mrSges(3,1) - t473 * mrSges(3,2) + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t359) * t338 + (t122 * t125 - t282 * t39) * mrSges(6,2) + t339 * (-Ifges(4,5) * t283 - Ifges(4,6) * t282) / 0.2e1 + t298 * (mrSges(4,1) * t282 - mrSges(4,2) * t283) - t280 * (-Ifges(4,4) * t283 - Ifges(4,2) * t282) / 0.2e1 + t313 * t218 + t636 * t168 + (mrSges(4,2) * t261 - mrSges(4,3) * t110 + Ifges(4,1) * t223 + Ifges(4,4) * t222 + Ifges(4,5) * t338) * t287 + m(5) * (t104 * t176 + t112 * t65 + t113 * t64 + t119 * t35 + t120 * t34 + t156 * t159) + m(6) * (t116 * t125 + t14 * t38 + t143 * t66 + t15 * t39 + t56 * t7 + t57 * t8) + m(7) * (t1 * t23 + t12 * t36 + t17 * t4 + t18 * t3 + t2 * t22 + t5 * t51) + (Ifges(6,1) * t122 + Ifges(6,5) * t282) * t550 + (Ifges(4,6) * t359 / 0.2e1 - t313 * mrSges(4,1)) * t222 + (Ifges(7,4) * t142 + Ifges(7,2) * t141) * t578 + (Ifges(7,4) * t69 + Ifges(7,2) * t70) * t562 - t104 * t409 + (Ifges(5,4) * t260 + Ifges(5,2) * t259) * t556 + (Ifges(5,5) * t260 + Ifges(5,6) * t259) * t549 + (t1 * t141 - t142 * t2 - t17 * t69 + t18 * t70) * mrSges(7,3) + (-t111 * mrSges(4,3) - t223 * Ifges(4,4) - Ifges(4,6) * t338 + Ifges(6,5) * t568 + Ifges(6,6) * t569 - Ifges(4,2) * t222 + Ifges(5,6) * t556 + Ifges(5,5) * t557 + t624 * t549 + t620 / 0.2e1 + t625) * t606 + m(3) * (t240 * t473 + t241 * t306 - t291 * t294 + t292 * t293) + t259 * t572 + t70 * t574 + t142 * t582 + t187 * t558 + t186 * t559 + t122 * t564 + t260 * t571 + (Ifges(7,5) * t142 + Ifges(7,6) * t141) * t570 + (Ifges(7,5) * t69 + Ifges(7,6) * t70) * t554 + (-t112 * t187 + t113 * t186 + t259 * t34 - t260 * t35) * mrSges(5,3) + (Ifges(5,5) * t187 + Ifges(6,5) * t122 + Ifges(5,6) * t186 + t282 * t624) * t544 + (-t331 * mrSges(5,1) - t302 * mrSges(3,1) + t621 * mrSges(3,2) + (t297 * t458 + mrSges(2,2)) * t368 + (t352 * t458 - t457) * t364 - (mrSges(4,1) + t380) * t232 + t421 * t192 + t611 * t193 + t599 * t233 + t615 * (-pkin(4) * t331 - t232 * t351 + t233 * t360)) * g(1) + (t164 * t283 - t165 * t282) * mrSges(4,3) + t623 * t282 / 0.2e1 + (-Ifges(4,1) * t283 - Ifges(4,4) * t282) * t542 + Ifges(2,3) * qJDD(1) + (Ifges(6,4) * t122 + Ifges(6,6) * t282) * t552 + t293 * t289 - t294 * t288 - t283 * t201 / 0.2e1 + t245 * (Ifges(5,4) * t187 + Ifges(5,2) * t186 + Ifges(5,6) * t282) / 0.2e1 + t38 * (mrSges(6,1) * t282 - mrSges(6,3) * t122) - t282 * t200 / 0.2e1 + t160 * t252 + (-m(4) * t424 - t237 * mrSges(4,1) - t304 * mrSges(3,1) - t303 * mrSges(3,2) + t364 * mrSges(2,2) - m(5) * (pkin(3) * t237 + t424) - t208 * mrSges(5,1) - t207 * mrSges(5,2) - m(6) * t382 - t197 * mrSges(6,1) - m(7) * (pkin(5) * t197 + t382) - t133 * mrSges(7,1) - t132 * mrSges(7,2) + t457 * t368 + t421 * t196 + t410 * t236) * g(2) + t189 * t204 + t188 * t205 + t156 * (-mrSges(5,1) * t186 + mrSges(5,2) * t187) + t64 * t173 + t65 * t174 + t176 * t88 + t143 * t40 + t141 * t10 / 0.2e1 + t5 * (-mrSges(7,1) * t141 + mrSges(7,2) * t142) + t15 * t130 + t14 * t131 + t119 * t114 + t120 * t115 + t116 * t102 + t4 * t87 + t3 * t86 + t12 * t74 + t36 * (-mrSges(7,1) * t70 + mrSges(7,2) * t69) + t57 * t67 + t56 * t68 + t69 * t62 / 0.2e1 + t51 * t16 + t627 * t121 + t629 * t167 + (t219 / 0.2e1 + t328 / 0.2e1 + t329 / 0.2e1 + t589) * t359 + (Ifges(7,1) * t142 + Ifges(7,4) * t141) * t579 + (Ifges(7,1) * t69 + Ifges(7,4) * t70) * t560 + (t422 * g(2) * t364 + (mrSges(3,1) * t295 - mrSges(3,2) * t296 + (m(3) * t536 - t305) * qJDD(1)) * pkin(1) + (t422 - t523) * g(1) * t368 + (mrSges(3,3) * t240 + Ifges(3,4) * t296 + Ifges(3,2) * t295 + Ifges(3,6) * t338) * t367 + (-mrSges(3,3) * t241 + Ifges(3,1) * t296 + Ifges(3,4) * t295 + Ifges(3,5) * t338) * t363 + ((t510 / 0.2e1 - t291 * mrSges(3,3) + t264 / 0.2e1) * t367 + (-t509 / 0.2e1 - t292 * mrSges(3,3) - t263 / 0.2e1 + (m(4) * t298 + t211) * pkin(2)) * t363 + (t367 * (Ifges(3,4) * t367 - Ifges(3,2) * t363) / 0.2e1 - t593) * t472) * qJD(2)) * t357 - t282 * t633 + t22 * t27 + t23 * t28 + t261 * t524 + (Ifges(5,1) * t187 + Ifges(5,4) * t186 + Ifges(5,5) * t282) * t546; (t417 - t289) * t291 + t281 * t633 + (-t125 * t191 + t281 * t39) * mrSges(6,2) + (-Ifges(7,4) * t387 - Ifges(7,2) * t388) * t562 + (Ifges(7,4) * t146 + Ifges(7,2) * t145) * t563 - t339 * (-Ifges(4,5) * t280 - Ifges(4,6) * t281) / 0.2e1 - t298 * (mrSges(4,1) * t281 - mrSges(4,2) * t280) + (Ifges(5,5) * t281 - t280 * t404) * t547 + (Ifges(5,6) * t281 - t280 * t402) * t548 + (-Ifges(7,1) * t387 - Ifges(7,4) * t388) * t560 + (Ifges(7,1) * t146 + Ifges(7,4) * t145) * t561 + (t399 * t570 + t401 * t578 + t403 * t579 + t405 * t5 + t433 * t62 + t636) * t311 + (-t112 * t495 - t113 * t641 + t596) * mrSges(5,3) + t594 * t156 * (mrSges(5,1) * t362 + t523) + t593 * qJD(1) ^ 2 * t357 ^ 2 + t104 * t408 + (t245 * t402 + t273 * t400 - t393 * t404) * qJD(4) / 0.2e1 + ((t110 * t503 + t111 * t502) * pkin(2) + t164 * t178 - t165 * t179 - t298 * t420) * m(4) + t204 * t442 + t205 * t443 + (t104 * t350 - t112 * t117 - t113 * t118 - t156 * t178) * m(5) + t328 + t329 + (t558 - t513) * t467 + t366 * t572 + t146 * t573 + t491 * t582 + (Ifges(5,1) * t362 + t517) * t557 - t496 * t559 + t191 * t565 + t362 * t571 + (-Ifges(7,5) * t387 - Ifges(7,6) * t388) * t554 + (Ifges(7,5) * t146 + Ifges(7,6) * t145) * t555 - ((-Ifges(3,2) * t441 + t264 + t335) * t367 + t339 * (Ifges(3,5) * t367 - Ifges(3,6) * t363)) * t472 / 0.2e1 + (-mrSges(3,1) * t621 - mrSges(3,2) * t302 + t616 * t413 - t615 * (t232 * t360 + t233 * t351 + t413) + t590 * t233 - t599 * t232) * g(2) - t619 * t190 + (-m(4) * t340 - m(5) * t432 + t305 + t524 - t615 * (-t287 * t360 - t351 * t606 + t340) + t622 * t606 + (-t405 + t592) * t287) * g(3) - (-Ifges(4,1) * t280 - t511 + t623) * t281 / 0.2e1 + t219 + (t418 + t288) * t292 + t350 * t88 + t165 * t521 + t319 * t40 - t10 * t492 / 0.2e1 - t135 * t468 / 0.2e1 - t179 * t252 + t244 * t67 - t118 * t173 - t117 * t174 + t147 * t27 + t148 * t28 + t263 * t441 / 0.2e1 + t613 * t86 + (t1 * t148 + t147 * t2 + t17 * t612 + t18 * t613 + t243 * t5 + t36 * t609) * m(7) + t614 * t243 + t612 * t87 + (Ifges(6,1) * t191 + Ifges(6,5) * t281) * t551 + (Ifges(6,4) * t191 + Ifges(6,6) * t281) * t553 - t38 * (mrSges(6,1) * t281 - mrSges(6,3) * t191) + (Ifges(6,5) * t191 - t280 * t400 + t281 * t624) * t545 + (-Ifges(4,2) * t281 + t201 - t272) * t280 / 0.2e1 + t627 * t299 + (mrSges(3,2) * t304 - t615 * (t236 * t351 - t237 * t360) + t590 * t236 + t599 * t237 + (-pkin(2) * t458 - mrSges(3,1)) * t303) * g(1) + t629 * t309 + t589 - (Ifges(6,1) * t550 + Ifges(6,4) * t552 + Ifges(6,5) * t544 + t448 + t564 + t591) * t300 - t281 * t634 + (-t174 * t467 - t173 * t468 + m(5) * ((-t112 * t366 - t113 * t362) * qJD(4) + t596) - t362 * t114 + t366 * t115) * t348 + t602 * t178 - t604 * t61 / 0.2e1 + (mrSges(7,1) * t604 - mrSges(7,2) * t603) * t36 + (-t1 * t492 + t17 * t603 - t18 * t604 - t2 * t491) * mrSges(7,3) + t605 * t102 - t211 * t420 + t608 * t130 + t609 * t74 + t610 * t131 + (t125 * t605 - t243 * t7 + t244 * t8 + t319 * t66 + t38 * t610 + t39 * t608) * m(6) + t495 * t558 - t164 * t522 + t200 * t542 + (Ifges(5,5) * t362 + Ifges(5,6) * t366) * t549 + (Ifges(5,2) * t366 + t518) * t556; t392 * t300 - (-t252 - t395) * t280 + t395 * qJD(4) + t218 + t614 * t309 + t366 * t114 + t362 * t115 + (t67 + (-t361 * t86 - t365 * t87) * qJD(6) + t597) * t311 - (t102 - t602) * t281 - t222 * mrSges(4,1) - t191 * t130 - t145 * t87 - t146 * t86 + t607 * t601 + (-t145 * t17 - t146 * t18 - t300 * t398 + t309 * t5 + t311 * t370 + t36 * t601) * m(7) + (-t125 * t281 - t309 * t7 + t311 * t8 - t601 * t38 - t39 * t640) * m(6) + (-t156 * t281 + t34 * t362 + t35 * t366 + t594 * (-t112 * t362 + t113 * t366)) * m(5) + (t164 * t281 + t165 * t280 + t261) * m(4) + (-t359 * g(3) + (-g(1) * t364 + g(2) * t368) * t357) * t458; (t128 * t401 + t129 * t403 + t149 * t399) * qJD(6) / 0.2e1 + (-(t232 * t366 + t331) * mrSges(5,2) + t193 * mrSges(6,2) - m(7) * (pkin(4) * t598 + t193 * pkin(10)) + (-mrSges(5,1) - t583) * t598 - t611 * t192) * g(2) + (Ifges(5,1) * t245 + t512) * t547 + (t174 - t520) * t113 - m(6) * (-t125 * t534 + t39 * t43) + (t389 + t448) * qJD(6) + (-t17 * t24 - t18 * t25 + t349 * t5 - t36 * t42) * m(7) + t620 + t61 * t433 + (Ifges(5,2) * t393 + t136 + t239) * t548 - t156 * (-mrSges(5,1) * t393 + mrSges(5,2) * t245) + (Ifges(5,5) * t245 + Ifges(5,6) * t393) * t545 + t5 * t406 + t500 * t573 + t501 * t574 + (Ifges(7,2) * t365 + t515) * t578 + (Ifges(7,1) * t361 + t514) * t579 + t361 * t582 + (t356 * t8 + t358 * t7) * t583 + (Ifges(7,5) * t361 + Ifges(7,6) * t365) * t570 + t625 + t619 * t394 + t349 * t16 + t102 * t534 - t112 * t173 - t43 * t130 - t25 * t86 - t24 * t87 + (Ifges(6,1) * t551 + Ifges(6,4) * t553 + Ifges(6,5) * t545 + t399 * t555 + t401 * t563 + t403 * t561 - t389 + t565 - t591) * t425 + (-g(1) * t197 - g(2) * t193 - g(3) * t251 + (-t466 + t501) * t18 + (-t465 + t500) * t17 + t595) * mrSges(7,3) + (m(7) * t370 - t465 * t87 - t466 * t86 + t597) * (pkin(10) + t533) + (m(6) * t38 - t607) * t42 + (-m(6) * t412 + t197 * mrSges(6,2) - mrSges(5,1) * t207 + mrSges(5,2) * t208 - m(7) * (pkin(10) * t197 + t412) + t611 * t196) * g(1) + (-m(6) * t411 + t251 * mrSges(6,2) - t409 - m(7) * (pkin(10) * t251 + t411) - t611 * (-t287 * t353 + t354 * t359)) * g(3) + t245 * t513 + t68 * t532 + t67 * t533 + t10 * t538 + t135 * t546; t365 * t27 + t361 * t28 - t607 * t394 + t397 * qJD(6) + t392 * t425 + t40 + (t1 * t361 + t149 * t398 + t2 * t365 - t394 * t36 + t378) * m(7) + (t38 * t394 - t39 * t425 + t378 + t66) * m(6); -t36 * (mrSges(7,1) * t129 + mrSges(7,2) * t128) + (Ifges(7,1) * t128 - t516) * t561 + t61 * t560 + (Ifges(7,5) * t128 - Ifges(7,6) * t129) * t555 - t17 * t86 + t18 * t87 - g(1) * (mrSges(7,1) * t132 - mrSges(7,2) * t133) - g(2) * ((-t193 * t361 - t233 * t365) * mrSges(7,1) + (-t193 * t365 + t233 * t361) * mrSges(7,2)) - g(3) * ((-t251 * t361 + t278) * mrSges(7,1) + (-t251 * t365 - t494) * mrSges(7,2)) + (t128 * t17 + t129 * t18) * mrSges(7,3) + t9 + (-Ifges(7,2) * t129 + t127 + t62) * t563 + t588;];
tau  = t11;
