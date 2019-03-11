% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:54:40
% EndTime: 2019-03-09 08:55:34
% DurationCPUTime: 33.24s
% Computational Cost: add. (23503->903), mult. (65284->1223), div. (0->0), fcn. (54518->16), ass. (0->426)
t346 = sin(pkin(11));
t347 = sin(pkin(6));
t357 = cos(qJ(2));
t485 = cos(pkin(11));
t418 = t357 * t485;
t402 = t347 * t418;
t353 = sin(qJ(2));
t455 = qJD(1) * t347;
t430 = t353 * t455;
t272 = -qJD(1) * t402 + t346 * t430;
t345 = sin(pkin(12));
t348 = cos(pkin(12));
t352 = sin(qJ(5));
t356 = cos(qJ(5));
t303 = t345 * t352 - t356 * t348;
t188 = t303 * t272;
t292 = t303 * qJD(5);
t573 = t188 + t292;
t305 = t345 * t356 + t348 * t352;
t187 = t305 * t272;
t293 = t305 * qJD(5);
t571 = t293 + t187;
t349 = cos(pkin(6));
t331 = qJD(1) * t349 + qJD(2);
t610 = -t331 / 0.2e1;
t469 = t349 * t353;
t334 = pkin(1) * t469;
t471 = t347 * t357;
t513 = pkin(8) + qJ(3);
t259 = (t471 * t513 + t334) * qJD(1);
t249 = t346 * t259;
t519 = pkin(1) * t349;
t335 = t357 * t519;
t328 = qJD(1) * t335;
t422 = t513 * t353;
t406 = t347 * t422;
t258 = -qJD(1) * t406 + t328;
t177 = t258 * t485 - t249;
t306 = -t357 * t346 - t353 * t485;
t280 = t306 * t347;
t275 = qJD(1) * t280;
t411 = pkin(2) * t430;
t197 = -pkin(3) * t275 + qJ(4) * t272 + t411;
t115 = t348 * t177 + t345 * t197;
t479 = t272 * t345;
t104 = pkin(9) * t479 + t115;
t518 = pkin(2) * t346;
t336 = qJ(4) + t518;
t512 = pkin(9) + t336;
t299 = t512 * t345;
t300 = t512 * t348;
t381 = -t356 * t299 - t300 * t352;
t114 = -t177 * t345 + t348 * t197;
t478 = t272 * t348;
t91 = -pkin(4) * t275 + pkin(9) * t478 + t114;
t582 = -qJD(4) * t303 + qJD(5) * t381 - t356 * t104 - t352 * t91;
t524 = t272 / 0.2e1;
t351 = sin(qJ(6));
t355 = cos(qJ(6));
t398 = -mrSges(7,1) * t355 + mrSges(7,2) * t351;
t584 = m(7) * pkin(5) + mrSges(6,1) - t398;
t412 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t609 = pkin(10) * t275 + t582;
t417 = t485 * t259;
t176 = t258 * t346 + t417;
t136 = -pkin(4) * t479 + t176;
t608 = t571 * pkin(5) + pkin(10) * t573 - t136;
t266 = qJD(5) + t272;
t525 = t266 / 0.2e1;
t237 = -t275 * t348 + t331 * t345;
t414 = t275 * t345 + t348 * t331;
t598 = t237 * t356 + t352 * t414;
t531 = t598 / 0.2e1;
t599 = -t237 * t352 + t356 * t414;
t533 = t599 / 0.2e1;
t143 = qJD(6) - t599;
t535 = t143 / 0.2e1;
t122 = t266 * t351 + t355 * t598;
t537 = t122 / 0.2e1;
t121 = t266 * t355 - t351 * t598;
t539 = t121 / 0.2e1;
t240 = pkin(2) * t331 + t258;
t162 = t240 * t485 - t249;
t152 = -t331 * pkin(3) + qJD(4) - t162;
t119 = -pkin(4) * t414 + t152;
t163 = t346 * t240 + t417;
t154 = qJ(4) * t331 + t163;
t341 = pkin(2) * t357 + pkin(1);
t291 = -t341 * t455 + qJD(3);
t171 = pkin(3) * t272 + qJ(4) * t275 + t291;
t105 = -t154 * t345 + t348 * t171;
t68 = pkin(4) * t272 - pkin(9) * t237 + t105;
t106 = t348 * t154 + t345 * t171;
t84 = pkin(9) * t414 + t106;
t36 = t352 * t68 + t356 * t84;
t31 = pkin(10) * t266 + t36;
t57 = -pkin(5) * t599 - pkin(10) * t598 + t119;
t17 = -t31 * t351 + t355 * t57;
t18 = t31 * t355 + t351 * t57;
t558 = -mrSges(6,1) * t119 - mrSges(7,1) * t17 + mrSges(7,2) * t18;
t54 = t122 * Ifges(7,5) + t121 * Ifges(7,6) + t143 * Ifges(7,3);
t506 = Ifges(6,4) * t598;
t82 = Ifges(6,2) * t599 + t266 * Ifges(6,6) + t506;
t603 = -t54 / 0.2e1 + t82 / 0.2e1;
t607 = -Ifges(6,4) * t531 + Ifges(7,5) * t537 - Ifges(6,2) * t533 - Ifges(6,6) * t525 + Ifges(7,6) * t539 + Ifges(7,3) * t535 - t558 - t603;
t447 = qJD(1) * qJD(2);
t288 = (qJDD(1) * t357 - t353 * t447) * t347;
t289 = (qJDD(1) * t353 + t357 * t447) * t347;
t213 = -t485 * t288 + t289 * t346;
t212 = qJDD(5) + t213;
t528 = t212 / 0.2e1;
t214 = t346 * t288 + t289 * t485;
t445 = qJDD(1) * t349;
t330 = qJDD(2) + t445;
t181 = -t214 * t345 + t330 * t348;
t182 = t214 * t348 + t330 * t345;
t75 = -qJD(5) * t598 + t181 * t356 - t182 * t352;
t544 = t75 / 0.2e1;
t74 = qJD(5) * t599 + t181 * t352 + t182 * t356;
t545 = t74 / 0.2e1;
t73 = qJDD(6) - t75;
t546 = t73 / 0.2e1;
t42 = -qJD(6) * t122 + t212 * t355 - t351 * t74;
t551 = t42 / 0.2e1;
t41 = qJD(6) * t121 + t212 * t351 + t355 * t74;
t552 = t41 / 0.2e1;
t456 = pkin(8) * t471 + t334;
t287 = t456 * qJD(2);
t439 = pkin(1) * t445;
t326 = t357 * t439;
t473 = t347 * t353;
t427 = qJD(3) * t473;
t446 = qJDD(1) * t347;
t438 = pkin(8) * t446;
t159 = -t353 * t438 + pkin(2) * t330 - qJ(3) * t289 + t326 + (-t287 - t427) * qJD(1);
t443 = qJD(2) * t519;
t407 = qJD(1) * t443;
t432 = t353 * t439 + (t407 + t438) * t357;
t452 = qJD(3) * t357;
t453 = qJD(2) * t353;
t168 = qJ(3) * t288 + (-pkin(8) * t453 + t452) * t455 + t432;
t107 = t159 * t485 - t346 * t168;
t100 = -t330 * pkin(3) + qJDD(4) - t107;
t64 = -t181 * pkin(4) + t100;
t19 = -t75 * pkin(5) - t74 * pkin(10) + t64;
t254 = -pkin(1) * t446 - pkin(2) * t288 + qJDD(3);
t111 = pkin(3) * t213 - qJ(4) * t214 + qJD(4) * t275 + t254;
t108 = t346 * t159 + t485 * t168;
t93 = qJ(4) * t330 + qJD(4) * t331 + t108;
t52 = t348 * t111 - t345 * t93;
t33 = pkin(4) * t213 - pkin(9) * t182 + t52;
t53 = t345 * t111 + t348 * t93;
t38 = pkin(9) * t181 + t53;
t450 = qJD(5) * t356;
t451 = qJD(5) * t352;
t7 = t352 * t33 + t356 * t38 + t68 * t450 - t451 * t84;
t5 = pkin(10) * t212 + t7;
t1 = qJD(6) * t17 + t19 * t351 + t355 * t5;
t2 = -qJD(6) * t18 + t19 * t355 - t351 * t5;
t561 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t9 = Ifges(7,5) * t41 + Ifges(7,6) * t42 + Ifges(7,3) * t73;
t606 = t561 + mrSges(6,1) * t64 + Ifges(7,5) * t552 + Ifges(7,6) * t551 + Ifges(7,3) * t546 + t9 / 0.2e1 + (-t528 - t212 / 0.2e1) * Ifges(6,6) + (-t544 - t75 / 0.2e1) * Ifges(6,2) + (-t545 - t74 / 0.2e1) * Ifges(6,4);
t585 = t119 * mrSges(6,2);
t142 = Ifges(6,4) * t599;
t83 = Ifges(6,1) * t598 + t266 * Ifges(6,5) + t142;
t605 = t585 + Ifges(6,1) * t531 + Ifges(6,4) * t533 + Ifges(6,5) * t525 + t83 / 0.2e1;
t592 = m(6) + m(7);
t593 = -m(5) - m(4);
t604 = t593 - t592;
t494 = t275 * Ifges(4,4);
t602 = Ifges(4,6) * t610 + t494 / 0.2e1;
t601 = Ifges(4,2) * t524 + t602;
t231 = -t299 * t352 + t300 * t356;
t581 = -qJD(4) * t305 - qJD(5) * t231 + t104 * t352 - t356 * t91;
t358 = cos(qJ(1));
t461 = t357 * t358;
t354 = sin(qJ(1));
t464 = t354 * t353;
t600 = t349 * t461 - t464;
t596 = t64 * mrSges(6,2) + 0.2e1 * Ifges(6,1) * t545 + 0.2e1 * Ifges(6,4) * t544 + 0.2e1 * Ifges(6,5) * t528;
t530 = t181 / 0.2e1;
t529 = t182 / 0.2e1;
t527 = t213 / 0.2e1;
t35 = -t352 * t84 + t356 * t68;
t591 = t35 * mrSges(6,1);
t16 = -mrSges(7,1) * t42 + mrSges(7,2) * t41;
t58 = mrSges(6,1) * t212 - mrSges(6,3) * t74;
t590 = t16 - t58;
t431 = t485 * pkin(2);
t340 = -t431 - pkin(3);
t515 = t348 * pkin(4);
t312 = t340 - t515;
t220 = t303 * pkin(5) - t305 * pkin(10) + t312;
t140 = t220 * t351 + t231 * t355;
t589 = -qJD(6) * t140 - t351 * t609 + t355 * t608;
t139 = t220 * t355 - t231 * t351;
t588 = qJD(6) * t139 + t351 * t608 + t355 * t609;
t587 = t105 * mrSges(5,1);
t586 = t106 * mrSges(5,2);
t401 = -mrSges(5,1) * t348 + mrSges(5,2) * t345;
t365 = m(5) * pkin(3) + mrSges(4,1) - t401;
t255 = pkin(2) * t349 + t335 - t406;
t268 = qJ(3) * t471 + t456;
t186 = t346 * t255 + t485 * t268;
t174 = qJ(4) * t349 + t186;
t467 = t353 * t346;
t279 = t347 * t467 - t402;
t332 = pkin(2) * t471;
t415 = -qJ(4) * t280 + t332;
t520 = pkin(1) * t347;
t202 = pkin(3) * t279 - t415 - t520;
t117 = t348 * t174 + t345 * t202;
t252 = t280 * t345 + t348 * t349;
t102 = pkin(9) * t252 + t117;
t116 = -t174 * t345 + t348 * t202;
t253 = -t280 * t348 + t345 * t349;
t90 = pkin(4) * t279 - pkin(9) * t253 + t116;
t583 = t356 * t102 + t352 * t90;
t580 = -pkin(5) * t275 - t581;
t510 = mrSges(6,3) * t598;
t124 = mrSges(6,1) * t266 - t510;
t63 = -mrSges(7,1) * t121 + mrSges(7,2) * t122;
t579 = t63 - t124;
t400 = mrSges(5,1) * t345 + mrSges(5,2) * t348;
t578 = t152 * t400;
t137 = -t188 * t351 - t275 * t355;
t448 = qJD(6) * t355;
t377 = -t351 * t292 + t305 * t448;
t575 = t137 + t377;
t138 = t188 * t355 - t275 * t351;
t449 = qJD(6) * t351;
t376 = t355 * t292 + t305 * t449;
t574 = t138 + t376;
t572 = mrSges(4,1) * t331 + mrSges(5,1) * t414 - mrSges(5,2) * t237 + mrSges(4,3) * t275;
t165 = t252 * t352 + t253 * t356;
t271 = t279 * t355;
t133 = -t165 * t351 + t271;
t564 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t386 = m(5) * qJ(4) - t564;
t397 = mrSges(7,1) * t351 + mrSges(7,2) * t355;
t570 = -t386 - t397;
t463 = t354 * t357;
t466 = t353 * t358;
t296 = -t349 * t463 - t466;
t172 = -mrSges(5,2) * t272 + mrSges(5,3) * t414;
t173 = mrSges(5,1) * t272 - mrSges(5,3) * t237;
t569 = t172 * t348 - t173 * t345;
t22 = mrSges(7,1) * t73 - mrSges(7,3) * t41;
t23 = -mrSges(7,2) * t73 + mrSges(7,3) * t42;
t568 = -t351 * t22 + t355 * t23;
t567 = -t345 * t52 + t348 * t53;
t566 = t1 * t355 - t2 * t351;
t509 = Ifges(3,4) * t353;
t565 = -t353 * (Ifges(3,1) * t357 - t509) / 0.2e1 + pkin(1) * (mrSges(3,1) * t353 + mrSges(3,2) * t357);
t8 = -qJD(5) * t36 + t33 * t356 - t352 * t38;
t373 = t418 - t467;
t457 = t306 * t349;
t226 = t354 * t457 + t358 * t373;
t221 = -t354 * t373 + t358 * t457;
t563 = t8 * mrSges(6,1) - t7 * mrSges(6,2) + Ifges(6,5) * t74 + Ifges(6,6) * t75 + Ifges(6,3) * t212;
t273 = qJD(2) * t280;
t274 = t373 * t347 * qJD(2);
t329 = t357 * t443;
t241 = t329 + (-qJD(2) * t422 + t452) * t347;
t423 = t513 * t347;
t242 = -t427 + (-t357 * t423 - t334) * qJD(2);
t158 = t485 * t241 + t346 * t242;
t145 = qJD(4) * t349 + t158;
t428 = t347 * t453;
t156 = pkin(2) * t428 - pkin(3) * t273 - qJ(4) * t274 + qJD(4) * t280;
t94 = -t145 * t345 + t348 * t156;
t67 = -pkin(9) * t274 * t348 - pkin(4) * t273 + t94;
t477 = t274 * t345;
t95 = t348 * t145 + t345 * t156;
t80 = -pkin(9) * t477 + t95;
t15 = -qJD(5) * t583 - t352 * t80 + t356 * t67;
t344 = pkin(12) + qJ(5);
t342 = sin(t344);
t343 = cos(t344);
t562 = t342 * t412 - t343 * t584 - t365;
t410 = pkin(8) * t428;
t228 = -qJD(1) * t410 + t432;
t229 = -pkin(8) * t289 - t353 * t407 + t326;
t560 = t229 * mrSges(3,1) + t107 * mrSges(4,1) - t228 * mrSges(3,2) - t108 * mrSges(4,2) + Ifges(3,5) * t289 + Ifges(4,5) * t214 + Ifges(3,6) * t288 - Ifges(4,6) * t213;
t526 = -t266 / 0.2e1;
t534 = -t599 / 0.2e1;
t536 = -t143 / 0.2e1;
t538 = -t122 / 0.2e1;
t540 = -t121 / 0.2e1;
t557 = Ifges(7,5) * t538 - Ifges(6,2) * t534 - Ifges(6,6) * t526 + Ifges(7,6) * t540 + Ifges(7,3) * t536 + t558;
t555 = Ifges(7,1) * t552 + Ifges(7,4) * t551 + Ifges(7,5) * t546;
t505 = Ifges(7,4) * t122;
t55 = Ifges(7,2) * t121 + Ifges(7,6) * t143 + t505;
t548 = t55 / 0.2e1;
t120 = Ifges(7,4) * t121;
t56 = Ifges(7,1) * t122 + Ifges(7,5) * t143 + t120;
t547 = -t56 / 0.2e1;
t541 = Ifges(5,1) * t529 + Ifges(5,4) * t530 + Ifges(5,5) * t527;
t532 = -t598 / 0.2e1;
t522 = t348 / 0.2e1;
t521 = t355 / 0.2e1;
t511 = mrSges(6,3) * t599;
t265 = Ifges(4,4) * t272;
t508 = Ifges(5,4) * t345;
t507 = Ifges(5,4) * t348;
t504 = Ifges(7,4) * t351;
t503 = Ifges(7,4) * t355;
t502 = t599 * Ifges(6,6);
t501 = t598 * Ifges(6,5);
t500 = t162 * mrSges(4,3);
t499 = t163 * mrSges(4,3);
t498 = t414 * Ifges(5,6);
t497 = t237 * Ifges(5,5);
t496 = t266 * Ifges(6,3);
t495 = t275 * Ifges(4,1);
t493 = t331 * Ifges(3,5);
t492 = t331 * Ifges(4,5);
t491 = t331 * Ifges(3,6);
t484 = t599 * t351;
t483 = t599 * t355;
t476 = t279 * t351;
t475 = t305 * t351;
t474 = t305 * t355;
t472 = t347 * t354;
t470 = t347 * t358;
t444 = pkin(4) * t345 * t347;
t441 = -m(3) * pkin(1) - mrSges(2,1);
t435 = t56 * t521;
t429 = t357 * t455;
t425 = -t345 * (t237 * Ifges(5,4) + Ifges(5,2) * t414 + Ifges(5,6) * t272) / 0.2e1;
t424 = (t237 * Ifges(5,1) + Ifges(5,4) * t414 + Ifges(5,5) * t272) * t522;
t34 = -t75 * mrSges(6,1) + t74 * mrSges(6,2);
t420 = -t449 / 0.2e1;
t416 = t213 * mrSges(4,1) + t214 * mrSges(4,2);
t118 = -t181 * mrSges(5,1) + t182 * mrSges(5,2);
t190 = -t221 * t343 - t342 * t470;
t189 = t221 * t342 - t343 * t470;
t157 = t241 * t346 - t485 * t242;
t290 = pkin(2) * t469 - t423;
t413 = t290 * t354 - t358 * t341;
t409 = mrSges(3,3) * t430;
t408 = mrSges(3,3) * t429;
t403 = t600 * pkin(2);
t129 = pkin(4) * t477 + t157;
t396 = Ifges(5,1) * t348 - t508;
t395 = Ifges(7,1) * t355 - t504;
t394 = -Ifges(5,2) * t345 + t507;
t393 = -Ifges(7,2) * t351 + t503;
t392 = Ifges(5,5) * t348 - Ifges(5,6) * t345;
t391 = Ifges(7,5) * t355 - Ifges(7,6) * t351;
t390 = -t17 * t351 + t18 * t355;
t44 = pkin(10) * t279 + t583;
t185 = t255 * t485 - t346 * t268;
t175 = -t349 * pkin(3) - t185;
t135 = -t252 * pkin(4) + t175;
t382 = t356 * t252 - t253 * t352;
t62 = -pkin(5) * t382 - t165 * pkin(10) + t135;
t21 = t351 * t62 + t355 * t44;
t20 = -t351 * t44 + t355 * t62;
t77 = -mrSges(7,2) * t143 + mrSges(7,3) * t121;
t78 = mrSges(7,1) * t143 - mrSges(7,3) * t122;
t388 = -t351 * t78 + t355 * t77;
t48 = -t102 * t352 + t356 * t90;
t384 = t105 * t345 - t106 * t348;
t134 = t165 * t355 + t476;
t123 = -mrSges(6,2) * t266 + t511;
t380 = -t123 - t388;
t30 = -pkin(5) * t266 - t35;
t378 = t30 * t397;
t14 = -t102 * t451 + t352 * t67 + t356 * t80 + t90 * t450;
t364 = t349 * t373;
t225 = t306 * t358 - t354 * t364;
t339 = pkin(3) + t515;
t350 = -pkin(9) - qJ(4);
t371 = t225 * t350 + t226 * t339 + t354 * t444 - t413;
t222 = t354 * t306 + t358 * t364;
t367 = g(1) * t225 + g(2) * t222 - g(3) * t279;
t360 = (-t17 * t355 - t18 * t351) * qJD(6) + t566;
t327 = Ifges(3,4) * t429;
t322 = Ifges(3,3) * t330;
t321 = Ifges(4,3) * t330;
t307 = -t332 - t520;
t301 = -pkin(8) * t473 + t335;
t298 = (-mrSges(3,1) * t357 + mrSges(3,2) * t353) * t347;
t297 = -t349 * t464 + t461;
t295 = -t349 * t466 - t463;
t286 = t329 - t410;
t285 = t456 * qJD(1);
t284 = -pkin(8) * t430 + t328;
t283 = -mrSges(3,2) * t331 + t408;
t282 = mrSges(3,1) * t331 - t409;
t257 = Ifges(3,1) * t430 + t327 + t493;
t256 = t491 + (t357 * Ifges(3,2) + t509) * t455;
t246 = -mrSges(4,2) * t331 - mrSges(4,3) * t272;
t244 = -t280 * t343 + t342 * t349;
t203 = mrSges(4,1) * t272 - mrSges(4,2) * t275;
t201 = mrSges(4,1) * t330 - mrSges(4,3) * t214;
t200 = -mrSges(4,2) * t330 - mrSges(4,3) * t213;
t196 = -t265 + t492 - t495;
t194 = t226 * t343 + t342 * t472;
t193 = t226 * t342 - t343 * t472;
t130 = t272 * Ifges(5,3) + t497 + t498;
t128 = t194 * t355 - t225 * t351;
t127 = -t194 * t351 - t225 * t355;
t126 = mrSges(5,1) * t213 - mrSges(5,3) * t182;
t125 = -mrSges(5,2) * t213 + mrSges(5,3) * t181;
t113 = qJD(5) * t165 + t274 * t305;
t112 = qJD(5) * t382 - t274 * t303;
t98 = t182 * Ifges(5,4) + t181 * Ifges(5,2) + t213 * Ifges(5,6);
t97 = pkin(5) * t598 - pkin(10) * t599;
t96 = -mrSges(6,1) * t599 + mrSges(6,2) * t598;
t81 = t496 + t501 + t502;
t61 = -qJD(6) * t134 - t112 * t351 - t273 * t355;
t60 = qJD(6) * t133 + t112 * t355 - t273 * t351;
t59 = -mrSges(6,2) * t212 + mrSges(6,3) * t75;
t45 = pkin(5) * t113 - pkin(10) * t112 + t129;
t43 = -pkin(5) * t279 - t48;
t25 = t35 * t355 + t351 * t97;
t24 = -t35 * t351 + t355 * t97;
t13 = pkin(5) * t273 - t15;
t12 = -pkin(10) * t273 + t14;
t10 = t41 * Ifges(7,4) + t42 * Ifges(7,2) + t73 * Ifges(7,6);
t6 = -pkin(5) * t212 - t8;
t4 = -qJD(6) * t21 - t12 * t351 + t355 * t45;
t3 = qJD(6) * t20 + t12 * t355 + t351 * t45;
t11 = [m(4) * (t107 * t185 + t108 * t186 - t157 * t162 + t158 * t163 + t254 * t307) + (Ifges(7,4) * t60 + Ifges(7,2) * t61) * t539 + (Ifges(7,4) * t134 + Ifges(7,2) * t133) * t551 - (mrSges(4,2) * t254 - mrSges(4,3) * t107 + Ifges(4,1) * t214 - Ifges(4,4) * t213 + Ifges(4,5) * t330) * t280 + (t252 * t53 - t253 * t52) * mrSges(5,3) + m(3) * (t228 * t456 + t229 * t301 - t284 * t287 + t285 * t286) + t307 * t416 + m(5) * (t100 * t175 + t105 * t94 + t106 * t95 + t116 * t52 + t117 * t53 + t152 * t157) + m(7) * (t1 * t21 + t13 * t30 + t17 * t4 + t18 * t3 + t2 * t20 + t43 * t6) + (Ifges(7,5) * t60 + Ifges(7,6) * t61) * t535 + (Ifges(7,5) * t134 + Ifges(7,6) * t133) * t546 + m(6) * (t119 * t129 + t135 * t64 + t14 * t36 + t15 * t35 + t48 * t8 + t583 * t7) + t583 * t59 - (-t499 + t498 / 0.2e1 + t497 / 0.2e1 + t587 - t586 + t291 * mrSges(4,1) + t496 / 0.2e1 + t81 / 0.2e1 + t591 - t36 * mrSges(6,2) + t130 / 0.2e1 + t502 / 0.2e1 + t501 / 0.2e1 + t601 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t272 + t602) * t273 + t456 * (-mrSges(3,2) * t330 + mrSges(3,3) * t288) + (-t35 * mrSges(6,3) + t605) * t112 + (mrSges(6,3) * t7 - t606) * t382 + (-mrSges(6,3) * t36 + t607) * t113 + (t321 / 0.2e1 + t322 / 0.2e1 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t330 + t560) * t349 + (-Ifges(4,6) * t330 - t108 * mrSges(4,3) - Ifges(4,4) * t214 + t254 * mrSges(4,1) + t52 * mrSges(5,1) + Ifges(5,6) * t181 + Ifges(5,5) * t182 - t53 * mrSges(5,2) + (Ifges(4,2) + Ifges(5,3)) * t213 + t563) * t279 + (-t295 * mrSges(3,1) + t600 * mrSges(3,2) + (-t290 * t604 + mrSges(2,2)) * t358 + (-t341 * t604 - t441) * t354 + t412 * t189 - t365 * t221 + t584 * t190 + t570 * t222 + t592 * (-t221 * t339 + t222 * t350 - t358 * t444)) * g(1) + t30 * (-mrSges(7,1) * t61 + mrSges(7,2) * t60) + (-t8 * mrSges(6,3) + t596) * t165 + (Ifges(7,1) * t60 + Ifges(7,4) * t61) * t537 + (Ifges(7,1) * t134 + Ifges(7,4) * t133) * t552 + t48 * t58 + t60 * t56 / 0.2e1 + Ifges(2,3) * qJDD(1) + t43 * t16 + t20 * t22 + t21 * t23 + (t1 * t133 - t134 * t2 - t17 * t60 + t18 * t61) * mrSges(7,3) - t572 * t157 + t301 * (mrSges(3,1) * t330 - mrSges(3,3) * t289) + t286 * t283 - t287 * t282 + t100 * (-mrSges(5,1) * t252 + mrSges(5,2) * t253) + t252 * t98 / 0.2e1 + t158 * t246 + t186 * t200 + t185 * t201 + t13 * t63 + (t425 + t424 - t500 + t392 * t524 + t414 * t394 / 0.2e1 + t237 * t396 / 0.2e1 + t578 + t492 / 0.2e1 + t291 * mrSges(4,2) - t495 / 0.2e1 - t265 / 0.2e1 + t196 / 0.2e1 + (-t105 * t348 - t106 * t345) * mrSges(5,3)) * t274 + ((mrSges(3,1) * t288 - mrSges(3,2) * t289 + (m(3) * t520 - t298) * qJDD(1)) * pkin(1) + (mrSges(3,3) * t228 + Ifges(3,4) * t289 + Ifges(3,2) * t288 + Ifges(3,6) * t330) * t357 + (-mrSges(3,3) * t229 + Ifges(3,1) * t289 + Ifges(3,4) * t288 + Ifges(3,5) * t330) * t353 + ((t493 / 0.2e1 - t284 * mrSges(3,3) + t257 / 0.2e1) * t357 + (-t491 / 0.2e1 - t285 * mrSges(3,3) - t256 / 0.2e1 + (m(4) * t291 + t203) * pkin(2)) * t353 + (t357 * (Ifges(3,4) * t357 - Ifges(3,2) * t353) / 0.2e1 - t565) * t455) * qJD(2) + (g(1) * t358 + g(2) * t354) * (-m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3) - t400)) * t347 + t3 * t77 + t4 * t78 + t14 * t123 + (Ifges(5,5) * t253 + Ifges(5,6) * t252) * t527 + (Ifges(5,1) * t253 + Ifges(5,4) * t252) * t529 + (Ifges(5,4) * t253 + Ifges(5,2) * t252) * t530 + t253 * t541 + t61 * t548 + t134 * t555 + t15 * t124 + t117 * t125 + t116 * t126 + t129 * t96 + t133 * t10 / 0.2e1 + t6 * (-mrSges(7,1) * t133 + mrSges(7,2) * t134) + t135 * t34 + (-t297 * mrSges(3,1) - t296 * mrSges(3,2) + t354 * mrSges(2,2) - m(6) * t371 - t194 * mrSges(6,1) - m(7) * (pkin(5) * t194 + t371) - t128 * mrSges(7,1) - t127 * mrSges(7,2) + t441 * t358 + t412 * t193 - t365 * t226 + t386 * t225 - t593 * t413) * g(2) + t95 * t172 + t94 * t173 + t175 * t118; (-Ifges(7,4) * t376 - Ifges(7,2) * t377) * t539 + (Ifges(7,4) * t138 + Ifges(7,2) * t137) * t540 + (-Ifges(7,1) * t376 - Ifges(7,4) * t377) * t537 + (Ifges(7,1) * t138 + Ifges(7,4) * t137) * t538 + (t408 - t283) * t284 - t10 * t475 / 0.2e1 + (-Ifges(4,5) * t272 + Ifges(4,6) * t275) * t610 + (t409 + t282) * t285 + t201 * t431 + t256 * t430 / 0.2e1 + t272 * t578 - t275 * t586 - ((-Ifges(3,2) * t430 + t257 + t327) * t357 + t331 * (Ifges(3,5) * t357 - Ifges(3,6) * t353)) * t455 / 0.2e1 + (-mrSges(3,1) * t600 - mrSges(3,2) * t295 + t593 * t403 - t592 * (t221 * t350 + t222 * t339 + t403) + t562 * t222 - t570 * t221) * g(2) + t565 * qJD(1) ^ 2 * t347 ^ 2 + (Ifges(6,5) * t188 - Ifges(6,3) * t275) * t526 + (-t119 * t188 - t275 * t36) * mrSges(6,2) - t272 * (-Ifges(5,3) * t275 - t272 * t392) / 0.2e1 + (Ifges(4,2) * t275 + t196 - t265) * t524 + (Ifges(6,1) * t188 - Ifges(6,5) * t275) * t532 - t414 * (-Ifges(5,6) * t275 - t272 * t394) / 0.2e1 + (Ifges(6,4) * t188 - Ifges(6,6) * t275) * t534 - t237 * (-Ifges(5,5) * t275 - t272 * t396) / 0.2e1 - t203 * t411 + (t348 * t125 - t345 * t126) * t336 + (-t119 * t136 + t231 * t7 + t312 * t64 + t35 * t581 + t36 * t582 + t381 * t8) * m(6) + (t1 * t140 + t139 * t2 + t17 * t589 + t18 * t588 + t30 * t580 - t381 * t6) * m(7) - t590 * t381 - (-Ifges(6,4) * t532 + t557 + t603) * t187 + (t162 * t176 - t163 * t177 - t291 * t411 + (t107 * t485 + t108 * t346) * pkin(2)) * m(4) - t272 * t500 + (mrSges(3,2) * t297 - t592 * (t225 * t339 - t226 * t350) + t562 * t225 + t570 * t226 + (pkin(2) * t604 - mrSges(3,1)) * t296) * g(1) - (t435 + t605) * t292 + t606 * t303 + t607 * t293 + (-Ifges(7,5) * t376 - Ifges(7,6) * t377) * t535 + (Ifges(7,5) * t138 + Ifges(7,6) * t137) * t536 + t100 * t401 + t560 + (t391 * t546 + t393 * t551 + t395 * t552 + t397 * t6 + t420 * t56 + t596) * t305 + t340 * t118 + t312 * t34 - t177 * t246 + t231 * t59 - t291 * (-mrSges(4,1) * t275 - mrSges(4,2) * t272) + (-Ifges(4,1) * t272 + t130 + t494 + t81) * t275 / 0.2e1 + (-t105 * t478 - t106 * t479 + t567) * mrSges(5,3) + (-t384 * qJD(4) + t100 * t340 - t105 * t114 - t106 * t115 - t152 * t176 + t336 * t567) * m(5) + t569 * qJD(4) + t572 * t176 + (-t303 * t7 - t305 * t8 + t35 * t573 - t36 * t571) * mrSges(6,3) - t575 * t55 / 0.2e1 + (mrSges(7,1) * t575 - mrSges(7,2) * t574) * t30 + (-t1 * t475 + t17 * t574 - t18 * t575 - t2 * t474) * mrSges(7,3) + t321 + t322 + t580 * t63 + t581 * t124 + t582 * t123 + t272 * t424 + t272 * t425 - t275 * t499 + t275 * t601 + t200 * t518 + t98 * t522 + (Ifges(5,5) * t345 + Ifges(5,6) * t348) * t527 + (Ifges(5,1) * t345 + t507) * t529 + (Ifges(5,2) * t348 + t508) * t530 + t345 * t541 + t138 * t547 + t474 * t555 + t275 * t587 - t136 * t96 + t139 * t22 + t140 * t23 + t588 * t77 + t589 * t78 + t275 * t591 + (-m(4) * t332 - m(5) * t415 + t298 - t592 * (-t279 * t339 + t280 * t350 + t332) - (-t397 + t564) * t280 - t562 * t279) * g(3) - t115 * t172 - t114 * t173 - t188 * t83 / 0.2e1; (t59 + (-t351 * t77 - t355 * t78) * qJD(6) + t568) * t305 + t380 * t292 + t590 * t303 - (-t96 + t572) * t275 + (t246 + t569) * t272 + t416 + t348 * t126 + t345 * t125 - t137 * t78 - t138 * t77 - t188 * t123 + t579 * t571 + (-t137 * t17 - t138 * t18 - t390 * t292 + t30 * t571 + t303 * t6 + t360 * t305) * m(7) + (t119 * t275 - t303 * t8 + t305 * t7 - t35 * t571 - t36 * t573) * m(6) + (t152 * t275 - t272 * t384 + t345 * t53 + t348 * t52) * m(5) + (-t162 * t275 + t163 * t272 + t254) * m(4) - ((-g(1) * t354 + g(2) * t358) * t347 - t349 * g(3)) * t604; -t414 * t172 + t237 * t173 + t355 * t22 + t351 * t23 - t579 * t598 + t388 * qJD(6) + t380 * t599 + t118 + t34 + (t1 * t351 + t143 * t390 + t2 * t355 - t598 * t30 + t367) * m(7) + (t35 * t598 - t36 * t599 + t367 + t64) * m(6) + (t105 * t237 - t106 * t414 + t100 + t367) * m(5); (t142 + t83) * t534 + (t121 * t393 + t122 * t395 + t143 * t391) * qJD(6) / 0.2e1 + (t378 + t435) * qJD(6) + (-t506 + t54) * t532 + (-t123 + t511) * t35 + (-pkin(5) * t6 - t17 * t24 - t18 * t25) * m(7) + t6 * t398 + t563 - pkin(5) * t16 + (Ifges(6,1) * t532 + Ifges(6,5) * t526 + t391 * t536 + t393 * t540 + t395 * t538 - t378 - t585) * t599 + t557 * t598 + ((-t449 + t484) * t18 + (-t448 + t483) * t17 + t566) * mrSges(7,3) + (m(7) * t360 - t448 * t78 - t449 * t77 + t568) * pkin(10) - t25 * t77 - t24 * t78 + (-m(7) * t30 + t510 - t579) * t36 + (t193 * t584 + t194 * t412) * g(1) + (t412 * t244 - t584 * (t280 * t342 + t343 * t349)) * g(3) + (-t189 * t584 + t190 * t412) * g(2) + t55 * t420 + t10 * t521 + t82 * t531 + (Ifges(7,5) * t351 + Ifges(7,6) * t355) * t546 + t483 * t547 + t484 * t548 + (Ifges(7,2) * t355 + t504) * t551 + (Ifges(7,1) * t351 + t503) * t552 + t351 * t555; -t30 * (mrSges(7,1) * t122 + mrSges(7,2) * t121) + (Ifges(7,1) * t121 - t505) * t538 + t55 * t537 + (Ifges(7,5) * t121 - Ifges(7,6) * t122) * t536 - t17 * t77 + t18 * t78 - g(1) * (mrSges(7,1) * t127 - mrSges(7,2) * t128) - g(2) * ((-t190 * t351 - t222 * t355) * mrSges(7,1) + (-t190 * t355 + t222 * t351) * mrSges(7,2)) - g(3) * ((-t244 * t351 + t271) * mrSges(7,1) + (-t244 * t355 - t476) * mrSges(7,2)) + (t121 * t17 + t122 * t18) * mrSges(7,3) + t9 + (-Ifges(7,2) * t122 + t120 + t56) * t540 + t561;];
tau  = t11;
