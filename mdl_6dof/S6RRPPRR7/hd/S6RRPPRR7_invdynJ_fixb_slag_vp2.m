% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:11
% EndTime: 2019-03-09 09:17:56
% DurationCPUTime: 28.33s
% Computational Cost: add. (9485->872), mult. (22085->1110), div. (0->0), fcn. (15874->10), ass. (0->394)
t274 = sin(qJ(2));
t269 = sin(pkin(6));
t412 = qJD(1) * t269;
t379 = t274 * t412;
t538 = -t379 - qJD(5);
t478 = -pkin(3) - pkin(9);
t268 = -pkin(2) + t478;
t273 = sin(qJ(5));
t407 = qJD(5) * t273;
t278 = cos(qJ(2));
t378 = t278 * t412;
t239 = pkin(8) * t378;
t270 = cos(pkin(6));
t402 = t270 * qJD(1);
t396 = pkin(1) * t402;
t182 = t274 * t396 + t239;
t143 = -qJ(4) * t378 + t182;
t277 = cos(qJ(5));
t225 = qJ(3) * t378;
t289 = t269 * (pkin(4) * t278 + t268 * t274);
t97 = qJD(1) * t289 + t225;
t61 = t277 * t143 + t273 * t97;
t537 = -t268 * t407 - t61;
t272 = sin(qJ(6));
t276 = cos(qJ(6));
t327 = t272 * mrSges(7,1) + t276 * mrSges(7,2);
t458 = -mrSges(4,1) + mrSges(5,2);
t479 = -m(7) - m(6);
t496 = pkin(9) * t479 - mrSges(3,1) - mrSges(6,3) + t458;
t490 = t327 - t496;
t249 = qJD(2) + t402;
t155 = -t249 * t277 + t273 * t378;
t409 = qJD(2) * t274;
t377 = t269 * t409;
t343 = qJD(1) * t377;
t400 = qJDD(1) * t269;
t369 = t278 * t400;
t185 = t343 - t369;
t399 = qJDD(1) * t270;
t246 = qJDD(2) + t399;
t76 = qJD(5) * t155 + t185 * t277 - t246 * t273;
t535 = -t76 / 0.2e1;
t156 = t273 * t249 + t277 * t378;
t77 = qJD(5) * t156 - t185 * t273 - t246 * t277;
t534 = -t77 / 0.2e1;
t223 = qJ(4) * t379;
t415 = pkin(8) * t379 - t278 * t396;
t141 = -t223 + t415;
t532 = qJD(3) + t141 + t538 * (pkin(5) * t273 - pkin(10) * t277);
t531 = pkin(10) * t378 - t537;
t102 = t156 * t272 - t276 * t538;
t153 = qJD(6) - t155;
t311 = t156 * t276 + t272 * t538;
t41 = -Ifges(7,5) * t311 + t102 * Ifges(7,6) + t153 * Ifges(7,3);
t370 = qJD(3) + t415;
t344 = -t223 + t370;
t83 = t249 * t268 + t344;
t147 = -pkin(1) * t412 - pkin(2) * t378 - qJ(3) * t379;
t113 = pkin(3) * t378 + qJD(4) - t147;
t87 = (pkin(4) * t274 + pkin(9) * t278) * t412 + t113;
t40 = t273 * t87 + t277 * t83;
t33 = -pkin(10) * t538 + t40;
t222 = t249 * qJ(3);
t108 = -t143 - t222;
t96 = pkin(4) * t249 - t108;
t50 = -pkin(5) * t155 + pkin(10) * t156 + t96;
t14 = t272 * t50 + t276 * t33;
t516 = t14 * mrSges(7,2);
t13 = -t272 * t33 + t276 * t50;
t517 = t13 * mrSges(7,1);
t491 = t516 - t517;
t451 = Ifges(6,4) * t156;
t526 = t538 * Ifges(6,6);
t528 = t155 * Ifges(6,2);
t70 = -t451 - t526 + t528;
t530 = t70 / 0.2e1 - t41 / 0.2e1 + t491;
t408 = qJD(2) * t278;
t186 = (qJD(1) * t408 + qJDD(1) * t274) * t269;
t170 = qJDD(5) + t186;
t480 = t77 / 0.2e1;
t74 = qJDD(6) - t77;
t482 = t74 / 0.2e1;
t35 = qJD(6) * t311 + t170 * t276 - t272 * t76;
t485 = t35 / 0.2e1;
t34 = qJD(6) * t102 + t170 * t272 + t276 * t76;
t486 = t34 / 0.2e1;
t467 = pkin(1) * t270;
t395 = qJD(2) * t467;
t348 = qJD(1) * t395;
t390 = pkin(1) * t399;
t383 = pkin(8) * t369 + t274 * t390 + t278 * t348;
t503 = t246 * qJ(3) + t249 * qJD(3);
t507 = (pkin(8) * t409 + qJD(4) * t278) * t269;
t54 = -t185 * qJ(4) + qJD(1) * t507 - t383 - t503;
t47 = pkin(4) * t246 - t54;
t15 = -pkin(5) * t77 - pkin(10) * t76 + t47;
t429 = t269 * t274;
t244 = qJD(3) * t429;
t80 = -pkin(1) * t400 + t185 * pkin(2) - t186 * qJ(3) - qJD(1) * t244;
t309 = qJDD(4) - t80;
t38 = pkin(4) * t186 + t185 * t478 + t309;
t406 = qJD(5) * t277;
t250 = pkin(8) * t429;
t105 = -qJD(2) * t239 - qJDD(1) * t250 - t274 * t348 + t278 * t390;
t287 = qJDD(3) - t105;
t375 = qJD(4) * t429;
t283 = -qJ(4) * t186 - qJD(1) * t375 + t287;
t46 = t246 * t268 + t283;
t7 = t273 * t38 + t277 * t46 + t87 * t406 - t407 * t83;
t5 = pkin(10) * t170 + t7;
t1 = qJD(6) * t13 + t15 * t272 + t276 * t5;
t2 = -qJD(6) * t14 + t15 * t276 - t272 * t5;
t493 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t9 = Ifges(7,5) * t34 + Ifges(7,6) * t35 + Ifges(7,3) * t74;
t529 = t493 + Ifges(7,5) * t486 + Ifges(7,6) * t485 + Ifges(7,3) * t482 + Ifges(6,4) * t535 + t9 / 0.2e1 - t170 * Ifges(6,6) + (t534 - t480) * Ifges(6,2);
t152 = Ifges(6,4) * t155;
t527 = t538 * Ifges(6,5);
t328 = -mrSges(7,1) * t276 + mrSges(7,2) * t272;
t296 = m(7) * pkin(5) - t328;
t331 = mrSges(6,1) * t277 - mrSges(6,2) * t273;
t525 = -t296 * t277 - t331;
t334 = t1 * t276 - t2 * t272;
t403 = qJD(6) * t276;
t405 = qJD(6) * t272;
t524 = -t13 * t403 - t14 * t405 + t334;
t345 = t277 * t379;
t150 = t272 * t345 - t276 * t378;
t523 = t272 * t406 + t150;
t242 = t278 * t395;
t262 = t270 * qJD(3);
t100 = -qJ(4) * t377 - t242 - t262 + t507;
t481 = t76 / 0.2e1;
t520 = Ifges(6,4) * t534 - t170 * Ifges(6,5) + (-t481 + t535) * Ifges(6,1);
t39 = -t273 * t83 + t277 * t87;
t32 = pkin(5) * t538 - t39;
t518 = m(7) * t32;
t514 = Ifges(4,4) + Ifges(3,5);
t513 = -Ifges(3,6) + Ifges(4,6);
t12 = -mrSges(7,1) * t35 + mrSges(7,2) * t34;
t58 = mrSges(6,1) * t170 - mrSges(6,3) * t76;
t457 = t12 - t58;
t362 = pkin(10) * t273 + qJ(3);
t209 = pkin(5) * t277 + pkin(4) + t362;
t418 = t276 * t277;
t146 = t209 * t272 + t268 * t418;
t512 = -qJD(6) * t146 + t272 * t531 + t276 * t532;
t430 = t268 * t277;
t145 = t209 * t276 - t272 * t430;
t511 = qJD(6) * t145 + t272 * t532 - t276 * t531;
t510 = t96 * (mrSges(6,1) * t273 + mrSges(6,2) * t277);
t509 = mrSges(6,1) + t296;
t391 = m(7) * pkin(10) + mrSges(7,3);
t357 = mrSges(6,2) - t391;
t335 = -qJ(4) * t429 + t250;
t466 = pkin(1) * t278;
t380 = -pkin(2) - t466;
t356 = -pkin(3) + t380;
t107 = (-pkin(9) + t356) * t270 + t335;
t265 = t269 * pkin(1);
t427 = t269 * t278;
t254 = pkin(3) * t427;
t414 = pkin(2) * t427 + qJ(3) * t429;
t381 = t254 + t414;
t314 = pkin(4) * t429 + pkin(9) * t427 + t381;
t98 = t265 + t314;
t508 = t277 * t107 + t273 * t98;
t455 = mrSges(6,3) * t156;
t110 = -mrSges(6,1) * t538 + t455;
t57 = -mrSges(7,1) * t102 - mrSges(7,2) * t311;
t437 = -t110 + t57;
t349 = mrSges(5,3) * t378;
t436 = -mrSges(5,1) * t249 + mrSges(6,1) * t155 + mrSges(6,2) * t156 + t349;
t506 = t273 * t403 + t523;
t151 = (t272 * t278 + t274 * t418) * t412;
t404 = qJD(6) * t273;
t505 = -t272 * t404 + t276 * t406 + t151;
t352 = mrSges(3,3) * t379;
t354 = mrSges(4,2) * t379;
t504 = t352 + t354 + (-mrSges(3,1) - mrSges(4,1)) * t249;
t434 = qJD(5) * t40;
t8 = -t273 * t46 + t277 * t38 - t434;
t18 = mrSges(7,1) * t74 - mrSges(7,3) * t34;
t19 = -mrSges(7,2) * t74 + mrSges(7,3) * t35;
t501 = -t272 * t18 + t276 * t19;
t500 = -t13 * t272 + t14 * t276;
t499 = -t273 * t8 + t277 * t7;
t454 = Ifges(3,4) * t274;
t498 = -t274 * (Ifges(3,1) * t278 - t454) / 0.2e1 + pkin(1) * (mrSges(3,1) * t274 + mrSges(3,2) * t278);
t398 = mrSges(3,2) - mrSges(5,1) - mrSges(4,3);
t140 = t222 + t182;
t353 = mrSges(4,2) * t378;
t177 = mrSges(4,3) * t249 + t353;
t497 = -m(4) * t140 - t177;
t401 = m(5) - t479;
t392 = m(4) + t401;
t495 = qJ(3) * t392 - t398;
t494 = t8 * mrSges(6,1) - t7 * mrSges(6,2) + Ifges(6,5) * t76 + Ifges(6,6) * t77 + Ifges(6,3) * t170;
t258 = t274 * t467;
t112 = -t375 + (t258 + (pkin(8) - qJ(4)) * t427) * qJD(2);
t376 = t269 * t408;
t416 = qJ(3) * t376 + t244;
t84 = qJD(2) * t289 + t416;
t21 = -qJD(5) * t508 - t112 * t273 + t277 * t84;
t492 = -m(6) * qJ(3) - m(7) * t362 - t273 * mrSges(7,3) + t398 + t525;
t10 = t34 * Ifges(7,4) + t35 * Ifges(7,2) + t74 * Ifges(7,6);
t489 = t10 / 0.2e1;
t11 = t34 * Ifges(7,1) + t35 * Ifges(7,4) + t74 * Ifges(7,5);
t488 = t11 / 0.2e1;
t448 = Ifges(7,4) * t311;
t42 = Ifges(7,2) * t102 + Ifges(7,6) * t153 - t448;
t484 = t42 / 0.2e1;
t99 = Ifges(7,4) * t102;
t43 = -Ifges(7,1) * t311 + Ifges(7,5) * t153 + t99;
t483 = -t43 / 0.2e1;
t280 = -pkin(2) - pkin(3);
t477 = -t102 / 0.2e1;
t476 = t102 / 0.2e1;
t475 = t311 / 0.2e1;
t474 = -t311 / 0.2e1;
t473 = -t153 / 0.2e1;
t472 = t153 / 0.2e1;
t470 = -t156 / 0.2e1;
t469 = t156 / 0.2e1;
t6 = -pkin(5) * t170 - t8;
t463 = t273 * t6;
t456 = mrSges(6,3) * t155;
t453 = Ifges(5,4) * t274;
t452 = Ifges(5,4) * t278;
t450 = Ifges(6,4) * t273;
t449 = Ifges(6,4) * t277;
t447 = Ifges(7,4) * t272;
t446 = Ifges(7,4) * t276;
t445 = Ifges(4,5) * t278;
t442 = t155 * Ifges(6,6);
t441 = t156 * Ifges(6,5);
t440 = t538 * Ifges(6,3);
t435 = qJD(5) * t32;
t433 = t155 * t272;
t432 = t155 * t276;
t275 = sin(qJ(1));
t428 = t269 * t275;
t279 = cos(qJ(1));
t426 = t269 * t279;
t425 = t272 * t273;
t424 = t273 * t274;
t423 = t273 * t276;
t422 = t274 * t275;
t421 = t274 * t277;
t420 = t274 * t279;
t419 = t275 * t278;
t417 = t278 * t279;
t207 = pkin(8) * t427 + t258;
t413 = t279 * pkin(1) + pkin(8) * t428;
t410 = qJD(1) ^ 2 * t269 ^ 2;
t397 = Ifges(4,5) - Ifges(3,4) - Ifges(5,4);
t394 = Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t393 = -Ifges(3,6) / 0.2e1 + Ifges(4,6) / 0.2e1;
t387 = t277 * t427;
t204 = -t270 * t422 + t417;
t384 = t204 * pkin(2) + t413;
t158 = t270 * qJ(3) + t207;
t159 = -t265 - t414;
t373 = t268 * t406;
t367 = t412 / 0.2e1;
t364 = t403 / 0.2e1;
t363 = -pkin(1) * t275 + pkin(8) * t426;
t360 = t186 * mrSges(5,1) + t185 * mrSges(5,2);
t201 = -t270 * t417 + t422;
t189 = t201 * pkin(2);
t202 = t270 * t420 + t419;
t359 = qJ(3) * t202 - t189;
t203 = t270 * t419 + t420;
t195 = t203 * pkin(2);
t358 = qJ(3) * t204 - t195;
t355 = t280 * t429;
t351 = mrSges(3,3) * t378;
t350 = mrSges(5,3) * t379;
t347 = t204 * pkin(3) + t384;
t346 = t273 * t379;
t340 = t274 * t367;
t337 = -t202 * pkin(2) + t363;
t333 = mrSges(4,1) * t278 + mrSges(4,3) * t274;
t199 = -t270 * t277 + t273 * t427;
t200 = t270 * t273 + t387;
t332 = mrSges(6,1) * t199 + mrSges(6,2) * t200;
t123 = t200 * t272 + t276 * t429;
t124 = -t200 * t276 + t272 * t429;
t329 = mrSges(7,1) * t123 - mrSges(7,2) * t124;
t326 = Ifges(6,1) * t277 - t450;
t325 = Ifges(7,1) * t276 - t447;
t324 = Ifges(7,1) * t272 + t446;
t323 = -Ifges(6,2) * t273 + t449;
t322 = -Ifges(7,2) * t272 + t446;
t321 = Ifges(7,2) * t276 + t447;
t320 = Ifges(6,5) * t277 - Ifges(6,6) * t273;
t319 = Ifges(7,5) * t276 - Ifges(7,6) * t272;
t318 = Ifges(7,5) * t272 + Ifges(7,6) * t276;
t49 = pkin(10) * t429 + t508;
t132 = qJ(4) * t427 - t158;
t114 = t270 * pkin(4) - t132;
t65 = -pkin(5) * t199 + pkin(10) * t200 + t114;
t25 = t272 * t65 + t276 * t49;
t24 = -t272 * t49 + t276 * t65;
t316 = t273 * t40 + t277 * t39;
t315 = t203 * pkin(4) + t347;
t55 = -t107 * t273 + t277 * t98;
t60 = -t143 * t273 + t277 * t97;
t308 = -t202 * pkin(3) + t337;
t183 = -pkin(8) * t377 + t242;
t109 = mrSges(6,2) * t538 + t456;
t66 = -mrSges(7,2) * t153 + mrSges(7,3) * t102;
t67 = mrSges(7,1) * t153 + mrSges(7,3) * t311;
t307 = t272 * t67 - t276 * t66 - t109;
t125 = -t201 * t273 + t277 * t426;
t126 = t201 * t277 + t273 * t426;
t300 = t274 * (-Ifges(5,2) * t278 + t453);
t299 = t278 * (Ifges(5,1) * t274 - t452);
t298 = t278 * (Ifges(4,3) * t274 + t445);
t20 = -t107 * t407 + t277 * t112 + t273 * t84 + t98 * t406;
t104 = -pkin(8) * t343 + t383;
t286 = (-t13 * t276 - t14 * t272) * qJD(6) + t334;
t285 = -qJD(5) * t316 + t499;
t51 = t246 * t280 + t283;
t75 = t104 + t503;
t89 = -pkin(2) * t246 + t287;
t284 = t105 * mrSges(3,1) - t89 * mrSges(4,1) - t54 * mrSges(5,1) - t104 * mrSges(3,2) + t51 * mrSges(5,2) + t75 * mrSges(4,3);
t59 = -mrSges(6,2) * t170 + mrSges(6,3) * t77;
t281 = t59 + (-t272 * t66 - t276 * t67) * qJD(6) + t437 * qJD(5) + t501;
t271 = qJ(3) + pkin(4);
t248 = mrSges(5,1) * t429;
t235 = Ifges(3,4) * t378;
t234 = Ifges(4,5) * t379;
t220 = Ifges(4,2) * t246;
t219 = Ifges(3,3) * t246;
t206 = t270 * t466 - t250;
t205 = (-mrSges(3,1) * t278 + mrSges(3,2) * t274) * t269;
t194 = t203 * pkin(3);
t188 = t201 * pkin(3);
t184 = t207 * qJD(2);
t180 = (mrSges(5,1) * t274 - t278 * mrSges(5,2)) * t412;
t179 = pkin(2) * t379 - t225;
t178 = t333 * t412;
t176 = -mrSges(3,2) * t249 + t351;
t172 = mrSges(5,2) * t249 - t350;
t169 = Ifges(4,4) * t186;
t168 = Ifges(3,5) * t186;
t167 = Ifges(3,6) * t185;
t166 = Ifges(4,6) * t185;
t165 = t186 * mrSges(5,3);
t164 = t186 * mrSges(4,2);
t160 = t270 * t380 + t250;
t154 = t183 + t262;
t149 = t249 * t272 + t276 * t346;
t148 = t249 * t276 - t272 * t346;
t144 = pkin(2) * t377 - t416;
t142 = qJD(1) * t355 + t225;
t139 = Ifges(3,1) * t379 + t249 * Ifges(3,5) + t235;
t138 = t249 * Ifges(4,4) + (t274 * Ifges(4,1) - t445) * t412;
t137 = -t249 * Ifges(5,5) + (-t278 * Ifges(5,1) - t453) * t412;
t136 = t249 * Ifges(3,6) + (t278 * Ifges(3,2) + t454) * t412;
t135 = -t249 * Ifges(5,6) + (-t274 * Ifges(5,2) - t452) * t412;
t134 = t249 * Ifges(4,6) - Ifges(4,3) * t378 + t234;
t133 = t254 - t159;
t131 = -pkin(2) * t249 + t370;
t130 = t203 * t277 - t273 * t428;
t129 = t203 * t273 + t277 * t428;
t122 = qJD(5) * t199 + t277 * t377;
t121 = -qJD(5) * t387 - t270 * t407 + t273 * t377;
t120 = -mrSges(4,2) * t185 + mrSges(4,3) * t246;
t119 = t246 * mrSges(5,2) - t165;
t118 = -t246 * mrSges(4,1) + t164;
t117 = -mrSges(5,1) * t246 - mrSges(5,3) * t185;
t115 = t270 * t356 + t335;
t111 = qJD(2) * t355 + t416;
t94 = t249 * t280 + t344;
t92 = -pkin(5) * t156 - pkin(10) * t155;
t86 = t130 * t276 + t204 * t272;
t85 = -t130 * t272 + t204 * t276;
t71 = -t156 * Ifges(6,1) + t152 - t527;
t69 = -t440 - t441 + t442;
t64 = qJD(6) * t123 + t122 * t276 + t272 * t376;
t63 = -qJD(6) * t124 - t122 * t272 + t276 * t376;
t62 = -pkin(3) * t185 + t309;
t52 = -pkin(5) * t378 - t60;
t48 = -pkin(5) * t429 - t55;
t44 = pkin(5) * t121 - pkin(10) * t122 - t100;
t37 = -mrSges(6,1) * t77 + mrSges(6,2) * t76;
t23 = t272 * t92 + t276 * t39;
t22 = -t272 * t39 + t276 * t92;
t17 = -pkin(5) * t376 - t21;
t16 = pkin(10) * t376 + t20;
t4 = -qJD(6) * t25 - t16 * t272 + t276 * t44;
t3 = qJD(6) * t24 + t16 * t276 + t272 * t44;
t26 = [(Ifges(7,1) * t64 + Ifges(7,4) * t63) * t474 + (Ifges(7,1) * t124 + Ifges(7,4) * t123) * t486 + m(6) * (-t100 * t96 + t114 * t47 + t20 * t40 + t21 * t39 + t508 * t7 + t55 * t8) + t508 * t59 + m(3) * (t104 * t207 + t105 * t206 + t182 * t183 + t184 * t415) + (t1 * t123 - t124 * t2 - t13 * t64 + t14 * t63) * mrSges(7,3) + t48 * t12 + t25 * t19 + t24 * t18 - t6 * t329 - t47 * t332 + (mrSges(6,3) * t8 - Ifges(6,4) * t480 + t520) * t200 + t436 * t100 + t62 * t248 + (t168 / 0.2e1 - t167 / 0.2e1 + t219 / 0.2e1 + t169 / 0.2e1 + t220 / 0.2e1 + t166 / 0.2e1 + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(5,3)) * t246 + (Ifges(5,6) + t394) * t186 + (-Ifges(5,5) + t393) * t185 + t284) * t270 + (Ifges(7,4) * t64 + Ifges(7,2) * t63) * t476 + (Ifges(7,4) * t124 + Ifges(7,2) * t123) * t485 + m(7) * (t1 * t25 + t13 * t4 + t14 * t3 + t17 * t32 + t2 * t24 + t48 * t6) + m(5) * (t100 * t108 + t111 * t113 + t112 * t94 + t115 * t51 + t132 * t54 + t133 * t62) + m(4) * (t131 * t184 + t140 * t154 + t144 * t147 + t158 * t75 + t159 * t80 + t160 * t89) + t63 * t484 + (t7 * mrSges(6,3) + Ifges(6,4) * t481 - t529) * t199 + Ifges(2,3) * qJDD(1) + t133 * t360 + t206 * (mrSges(3,1) * t246 - mrSges(3,3) * t186) + t207 * (-mrSges(3,2) * t246 - mrSges(3,3) * t185) + t124 * t488 + t123 * t489 + (Ifges(7,5) * t64 + Ifges(7,6) * t63) * t472 + (Ifges(7,5) * t124 + Ifges(7,6) * t123) * t482 + (-m(3) * t363 - m(4) * t337 - m(5) * t308 + mrSges(2,1) * t275 + mrSges(2,2) * t279 + t357 * t125 + t495 * t201 + t509 * t126 + t490 * t202 - t479 * (t201 * pkin(4) - t308)) * g(1) + t17 * t57 + t55 * t58 + t64 * t43 / 0.2e1 + t32 * (-mrSges(7,1) * t63 + mrSges(7,2) * t64) + t3 * t66 + t4 * t67 + (-m(3) * t413 - m(5) * t347 - m(4) * t384 - m(7) * (pkin(5) * t130 + t315) - t86 * mrSges(7,1) - t85 * mrSges(7,2) - m(6) * t315 - t130 * mrSges(6,1) - mrSges(2,1) * t279 + mrSges(2,2) * t275 + t357 * t129 - t495 * t203 + t496 * t204) * g(2) + t504 * t184 + t20 * t109 + t21 * t110 + t114 * t37 + t115 * t119 + t132 * t117 + ((-mrSges(3,1) * t185 - mrSges(3,2) * t186 + (m(3) * t265 - t205) * qJDD(1)) * pkin(1) + (-t80 * mrSges(4,1) + t75 * mrSges(4,2) - t62 * mrSges(5,2) + t104 * mrSges(3,3) + t54 * mrSges(5,3) + (Ifges(5,5) - t513) * t246 - t397 * t186 + (-Ifges(5,1) - Ifges(3,2) - Ifges(4,3)) * t185) * t278 + (-t80 * mrSges(4,3) - t51 * mrSges(5,3) + t89 * mrSges(4,2) - t105 * mrSges(3,3) + (Ifges(5,6) + t514) * t246 + (Ifges(4,1) + Ifges(3,1) + Ifges(5,2)) * t186 + t397 * t185 + t494) * t274 + ((t113 * mrSges(5,2) + t147 * mrSges(4,1) + t134 / 0.2e1 - t136 / 0.2e1 + t137 / 0.2e1 - t140 * mrSges(4,2) - t182 * mrSges(3,3) - t108 * mrSges(5,3) + (-Ifges(5,5) / 0.2e1 + t393) * t249) * t274 + (t113 * mrSges(5,1) - t147 * mrSges(4,3) + t415 * mrSges(3,3) + t131 * mrSges(4,2) - t94 * mrSges(5,3) - t440 / 0.2e1 - t40 * mrSges(6,2) + t39 * mrSges(6,1) + t442 / 0.2e1 - t441 / 0.2e1 - t135 / 0.2e1 + t138 / 0.2e1 + t139 / 0.2e1 + t69 / 0.2e1 + (Ifges(5,6) / 0.2e1 + t394) * t249) * t278 + (t274 * (Ifges(4,1) * t278 + Ifges(4,5) * t274) / 0.2e1 - t298 / 0.2e1 + t278 * (Ifges(3,4) * t278 - Ifges(3,2) * t274) / 0.2e1 - t299 / 0.2e1 - t300 / 0.2e1 - t498) * t412) * qJD(2) + (g(1) * t279 + g(2) * t275) * (qJ(4) * t401 - mrSges(4,2) - mrSges(3,3) + mrSges(5,3))) * t269 + (-t40 * mrSges(6,3) - Ifges(6,4) * t470 + Ifges(7,3) * t472 + Ifges(7,5) * t474 + Ifges(7,6) * t476 - t528 / 0.2e1 + t526 / 0.2e1 + t96 * mrSges(6,1) - t530) * t121 + (-t39 * mrSges(6,3) + Ifges(6,1) * t470 + t152 / 0.2e1 - t527 / 0.2e1 + t71 / 0.2e1 + t96 * mrSges(6,2)) * t122 + t158 * t120 + t160 * t118 + t112 * t172 + t154 * t177 - t144 * t178 + t111 * t180 + t183 * t176 + t159 * (mrSges(4,1) * t185 - mrSges(4,3) * t186); t523 * t484 + (-m(4) * t414 - m(5) * t381 + mrSges(5,2) * t427 + t205 - t248 + t479 * t314 + (-t333 + (-t327 - mrSges(6,3)) * t278 + (-t273 * t391 + t525) * t274) * t269) * g(3) - t449 * t481 - (-t176 + t351 + t497) * t415 + (Ifges(7,5) * t151 - Ifges(7,6) * t150) * t473 + (-pkin(2) * t89 + qJ(3) * t75 - t147 * t179) * m(4) + (Ifges(7,1) * t151 - Ifges(7,4) * t150) * t475 + (t318 * t472 + t321 * t476 + t324 * t474) * t404 - t167 + t168 + t169 + t166 - (-(Ifges(6,3) * t278 + t274 * t320) * t538 + (Ifges(4,1) * t378 + t273 * t41 + t277 * t71 + t134 + t137 + t234) * t274 + (-Ifges(3,2) * t379 + t138 + t139 + t235 + t69) * t278 + (Ifges(6,6) * t278 + t274 * t323) * t155 + (t274 * t513 + t278 * t514) * t249) * t412 / 0.2e1 - (t155 * t323 - t320 * t538) * qJD(5) / 0.2e1 + (Ifges(7,4) * t151 - Ifges(7,2) * t150) * t477 + (Ifges(7,5) * t475 + Ifges(7,6) * t477 + Ifges(7,3) * t473 + t491) * t346 + (t156 * (Ifges(6,5) * t278 + t274 * t326) + t249 * (Ifges(5,5) * t274 - Ifges(5,6) * t278) + t278 * t135) * t367 - t327 * t463 + t47 * t331 - t131 * t353 - (t276 * t43 + t71) * t406 / 0.2e1 + (t300 + t299 + t298) * t410 / 0.2e1 + (qJD(6) * t43 + t10) * t425 / 0.2e1 + (-t319 * t482 - t322 * t485 - t325 * t486 + t340 * t70 + t364 * t42 + t520) * t273 + (-t117 + t120) * qJ(3) + t219 + t220 - t11 * t423 / 0.2e1 + (t141 * t96 + t268 * t285 + t271 * t47 - t39 * t60 - t40 * t61) * m(6) + (-t54 * qJ(3) - t108 * t141 - t113 * t142 - t143 * t94 + t280 * t51) * m(5) + t284 + (t373 - t52) * t57 - t450 * t480 + t151 * t483 + t59 * t430 + t529 * t277 + t530 * t407 + t94 * t349 + t108 * t350 + t140 * t354 + t280 * t119 + (-t373 - t60) * t110 + t271 * t37 + Ifges(5,3) * t246 + (-t147 * (mrSges(4,1) * t274 - mrSges(4,3) * t278) - t113 * (mrSges(5,1) * t278 + mrSges(5,2) * t274) - t39 * (mrSges(6,1) * t278 - mrSges(6,3) * t421) - t40 * (-mrSges(6,2) * t278 - mrSges(6,3) * t424)) * t412 + t457 * t268 * t273 + t537 * t109 + t498 * t410 + (t39 * t406 + t40 * t407 - t499) * mrSges(6,3) + (-m(4) * t131 + t352 - t504) * t182 + (t1 * t425 + t13 * t505 + t14 * t506 + t2 * t423) * mrSges(7,3) + (-mrSges(7,1) * t506 - mrSges(7,2) * t505) * t32 - pkin(2) * t118 + (-m(5) * t108 + m(6) * t96 - t436 - t497) * qJD(3) - t436 * t141 - t379 * t510 + ((-Ifges(7,3) * t273 - t277 * t319) * t472 + (-Ifges(7,5) * t273 - t277 * t325) * t474 + (-Ifges(7,6) * t273 - t277 * t322) * t476 - t510 + t326 * t469) * qJD(5) + t145 * t18 + t146 * t19 + t511 * t66 + t512 * t67 + (t1 * t146 + t145 * t2 + (t32 * t406 + t463) * t268 - t32 * t52 + t511 * t14 + t512 * t13) * m(7) - t143 * t172 + t179 * t178 - t142 * t180 - Ifges(5,5) * t185 + Ifges(5,6) * t186 + (-m(5) * (-t194 + t358) - m(4) * t358 + t479 * (t204 * pkin(4) - t194 - t195) + t492 * t204 + t490 * t203) * g(1) + (-m(5) * (-t188 + t359) - m(4) * t359 + t479 * (t202 * pkin(4) - t188 - t189) + t492 * t202 + t490 * t201) * g(2) + t136 * t340; -t148 * t67 - t149 * t66 + t164 - t165 + t458 * t246 + (-t177 + t436) * t249 + (t307 * qJD(5) + t457) * t273 + t281 * t277 + (-t273 * t109 + t277 * t437 - t178 - t180) * t379 + (-g(1) * t203 - g(2) * t201 + g(3) * t427) * t392 + ((-qJD(5) * t500 + t6) * t273 + (t286 + t435) * t277 - t13 * t148 - t14 * t149 + t32 * t345) * m(7) + (-t249 * t96 - t316 * t379 + t285) * m(6) + (t108 * t249 - t113 * t379 + t51) * m(5) + (-t140 * t249 + t147 * t379 + t89) * m(4); -t150 * t67 + t151 * t66 - m(7) * (t13 * t150 - t14 * t151) + m(5) * t62 + t401 * t270 * g(3) + (-m(7) * t6 + m(6) * (t8 + t434) - t457 + (m(7) * t500 - t307) * qJD(5)) * t277 + (m(7) * (t435 + t524) + m(6) * (-qJD(5) * t39 + t7) + t281) * t273 + ((-t436 * t278 + (t277 * t109 + t273 * t437 + t172) * t274 - m(6) * (-t278 * t96 + t39 * t424 - t40 * t421) - m(5) * (t108 * t278 - t274 * t94) + t424 * t518) * qJD(1) + (g(1) * t275 - g(2) * t279) * t401) * t269 + t360; (t13 * t432 + t14 * t433 + t524) * mrSges(7,3) + (t102 * t322 + t153 * t319 - t311 * t325) * qJD(6) / 0.2e1 + (-pkin(5) * t6 - t13 * t22 - t14 * t23) * m(7) - pkin(5) * t12 + t538 * (Ifges(6,5) * t155 + Ifges(6,6) * t156) / 0.2e1 + t6 * t328 - (Ifges(6,2) * t156 + t152 + t71) * t155 / 0.2e1 + t70 * t470 + (-Ifges(7,3) * t156 + t155 * t319) * t473 + (-Ifges(7,5) * t156 + t155 * t325) * t475 + (-Ifges(7,6) * t156 + t155 * t322) * t477 + t43 * t364 + t318 * t482 + t432 * t483 + t433 * t484 + t321 * t485 + t324 * t486 + t494 + (-t199 * t296 + t200 * t391 - t332) * g(3) + (Ifges(6,1) * t155 + t41 + t451) * t469 + (t456 - t109) * t39 + t156 * t517 - t42 * t405 / 0.2e1 + t272 * t488 + t276 * t489 + t153 * t32 * t327 - t23 * t66 - t22 * t67 + (m(7) * t286 - t403 * t67 - t405 * t66 + t501) * pkin(10) + (t129 * t509 + t130 * t357) * g(1) + (-t125 * t509 + t126 * t357) * g(2) - t96 * (-mrSges(6,1) * t156 + mrSges(6,2) * t155) - t156 * t516 + (-t437 - t455 - t518) * t40; -t32 * (-mrSges(7,1) * t311 + mrSges(7,2) * t102) + (Ifges(7,1) * t102 + t448) * t475 + t42 * t474 + (Ifges(7,5) * t102 + Ifges(7,6) * t311) * t473 - t13 * t66 + t14 * t67 - g(1) * (mrSges(7,1) * t85 - mrSges(7,2) * t86) - g(2) * ((-t126 * t272 + t202 * t276) * mrSges(7,1) + (-t126 * t276 - t202 * t272) * mrSges(7,2)) - g(3) * t329 + (t102 * t13 - t14 * t311) * mrSges(7,3) + t9 + (Ifges(7,2) * t311 + t43 + t99) * t477 + t493;];
tau  = t26;
