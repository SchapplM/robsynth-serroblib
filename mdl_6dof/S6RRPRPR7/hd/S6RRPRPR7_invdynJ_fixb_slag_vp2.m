% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:45:20
% EndTime: 2019-03-09 10:46:02
% DurationCPUTime: 26.69s
% Computational Cost: add. (12066->758), mult. (25989->971), div. (0->0), fcn. (17905->12), ass. (0->356)
t291 = sin(qJ(4));
t295 = cos(qJ(4));
t296 = cos(qJ(2));
t398 = qJD(1) * t296;
t292 = sin(qJ(2));
t399 = qJD(1) * t292;
t176 = -t291 * t399 - t295 * t398;
t177 = -t291 * t398 + t295 * t399;
t286 = sin(pkin(10));
t287 = cos(pkin(10));
t501 = t176 * t286 + t287 * t177;
t469 = -t501 / 0.2e1;
t107 = -t287 * t176 + t286 * t177;
t470 = t107 / 0.2e1;
t538 = Ifges(6,4) * t469 + Ifges(6,2) * t470;
t525 = Ifges(4,4) + Ifges(3,5);
t524 = Ifges(4,6) - Ifges(3,6);
t281 = -qJD(2) + qJD(4);
t462 = t281 / 0.2e1;
t537 = Ifges(6,5) * t462 + t501 * Ifges(6,1) / 0.2e1 - t107 * Ifges(6,4) / 0.2e1;
t536 = -Ifges(6,6) * t462 + t538;
t97 = qJD(6) + t107;
t473 = -t97 / 0.2e1;
t290 = sin(qJ(6));
t294 = cos(qJ(6));
t86 = t281 * t290 + t294 * t501;
t475 = -t86 / 0.2e1;
t85 = t281 * t294 - t290 * t501;
t477 = -t85 / 0.2e1;
t535 = Ifges(7,5) * t475 + Ifges(7,6) * t477 + Ifges(7,3) * t473;
t466 = t177 / 0.2e1;
t534 = mrSges(6,2) - mrSges(7,3);
t229 = -t294 * mrSges(7,1) + t290 * mrSges(7,2);
t533 = m(7) * pkin(5) + mrSges(6,1) - t229;
t530 = -pkin(5) * t501 - pkin(9) * t107;
t391 = qJD(1) * qJD(2);
t213 = -t296 * qJDD(1) + t292 * t391;
t214 = qJDD(1) * t292 + t296 * t391;
t405 = t295 * t296;
t410 = t291 * t292;
t318 = t405 + t410;
t306 = t318 * qJD(4);
t88 = -qJD(1) * t306 + t213 * t291 + t214 * t295;
t408 = t292 * t295;
t319 = t291 * t296 - t408;
t307 = t319 * qJD(4);
t89 = qJD(1) * t307 + t213 * t295 - t214 * t291;
t56 = -t286 * t88 + t287 * t89;
t57 = t286 * t89 + t287 * t88;
t273 = t292 * qJD(3);
t285 = qJDD(1) * pkin(1);
t111 = t213 * pkin(2) - t214 * qJ(3) - qJD(1) * t273 - t285;
t84 = -pkin(3) * t213 - t111;
t59 = -pkin(4) * t89 + qJDD(5) + t84;
t12 = -pkin(5) * t56 - pkin(9) * t57 + t59;
t418 = qJ(5) * t177;
t268 = pkin(7) * t399;
t209 = pkin(8) * t399 - t268;
t298 = -pkin(2) - pkin(3);
t371 = t298 * qJD(2);
t146 = qJD(3) + t371 - t209;
t269 = pkin(7) * t398;
t210 = -pkin(8) * t398 + t269;
t284 = qJD(2) * qJ(3);
t178 = t210 + t284;
t95 = t295 * t146 - t178 * t291;
t75 = t95 - t418;
t72 = pkin(4) * t281 + t75;
t419 = qJ(5) * t176;
t96 = t146 * t291 + t178 * t295;
t76 = t96 + t419;
t73 = t287 * t76;
t40 = t286 * t72 + t73;
t38 = pkin(9) * t281 + t40;
t179 = -qJD(1) * pkin(1) - pkin(2) * t398 - qJ(3) * t399;
t144 = pkin(3) * t398 - t179;
t101 = -pkin(4) * t176 + qJD(5) + t144;
t43 = pkin(5) * t107 - pkin(9) * t501 + t101;
t15 = -t290 * t38 + t294 * t43;
t280 = -qJDD(2) + qJDD(4);
t192 = t214 * pkin(7);
t361 = qJDD(3) + t192;
t114 = -pkin(8) * t214 + qJDD(2) * t298 + t361;
t191 = t213 * pkin(7);
t145 = qJDD(2) * qJ(3) + qJD(2) * qJD(3) - t191;
t115 = pkin(8) * t213 + t145;
t45 = -qJD(4) * t96 + t295 * t114 - t115 * t291;
t27 = pkin(4) * t280 - qJ(5) * t88 - qJD(5) * t177 + t45;
t394 = qJD(4) * t295;
t395 = qJD(4) * t291;
t44 = t291 * t114 + t295 * t115 + t146 * t394 - t178 * t395;
t29 = qJ(5) * t89 + qJD(5) * t176 + t44;
t10 = t286 * t27 + t287 * t29;
t8 = pkin(9) * t280 + t10;
t1 = qJD(6) * t15 + t12 * t290 + t294 * t8;
t16 = t290 * t43 + t294 * t38;
t163 = Ifges(5,4) * t176;
t164 = Ifges(5,4) * t177;
t2 = -qJD(6) * t16 + t12 * t294 - t290 * t8;
t339 = mrSges(7,1) * t290 + mrSges(7,2) * t294;
t429 = t286 * t76;
t39 = t287 * t72 - t429;
t37 = -pkin(5) * t281 - t39;
t311 = t37 * t339;
t331 = Ifges(7,5) * t294 - Ifges(7,6) * t290;
t439 = Ifges(7,4) * t294;
t334 = -Ifges(7,2) * t290 + t439;
t440 = Ifges(7,4) * t290;
t337 = Ifges(7,1) * t294 - t440;
t448 = t86 * Ifges(7,4);
t35 = t85 * Ifges(7,2) + t97 * Ifges(7,6) + t448;
t393 = qJD(6) * t290;
t362 = -t393 / 0.2e1;
t83 = Ifges(7,4) * t85;
t36 = Ifges(7,1) * t86 + Ifges(7,5) * t97 + t83;
t424 = t294 * t36;
t373 = t424 / 0.2e1;
t392 = qJD(6) * t294;
t431 = t177 * mrSges(5,3);
t432 = t176 * mrSges(5,3);
t444 = mrSges(7,3) * t294;
t445 = mrSges(7,3) * t290;
t461 = t290 / 0.2e1;
t463 = -t281 / 0.2e1;
t472 = t97 / 0.2e1;
t474 = t86 / 0.2e1;
t476 = t85 / 0.2e1;
t52 = qJDD(6) - t56;
t480 = t52 / 0.2e1;
t33 = -qJD(6) * t86 + t280 * t294 - t290 * t57;
t483 = t33 / 0.2e1;
t32 = qJD(6) * t85 + t280 * t290 + t294 * t57;
t484 = t32 / 0.2e1;
t5 = t32 * Ifges(7,4) + t33 * Ifges(7,2) + t52 * Ifges(7,6);
t503 = -mrSges(6,1) * t101 - mrSges(7,1) * t15 + mrSges(7,2) * t16 + t40 * mrSges(6,3) + t535 - t536;
t506 = t101 * mrSges(6,2) - t39 * mrSges(6,3) + t537;
t6 = Ifges(7,1) * t32 + Ifges(7,4) * t33 + Ifges(7,5) * t52;
t9 = t27 * t287 - t286 * t29;
t7 = -pkin(5) * t280 - t9;
t98 = t176 * Ifges(5,2) + t281 * Ifges(5,6) + t164;
t99 = t177 * Ifges(5,1) + t281 * Ifges(5,5) + t163;
t529 = (-Ifges(5,1) * t176 + t164 + t98) * t466 + (Ifges(7,5) * t290 + Ifges(7,6) * t294) * t480 + (Ifges(7,2) * t294 + t440) * t483 + (Ifges(7,1) * t290 + t439) * t484 - (-Ifges(5,2) * t177 + t163 + t99) * t176 / 0.2e1 - (Ifges(6,1) * t469 + Ifges(6,4) * t470 + Ifges(6,5) * t463 + t15 * t444 + t16 * t445 + t331 * t473 + t334 * t477 + t337 * t475 - t311 - t424 / 0.2e1 + t35 * t461 - t506) * t107 + (Ifges(6,3) + Ifges(5,3)) * t280 + (t331 * t472 + t334 * t476 + t337 * t474 + t311 + t373) * qJD(6) + (-Ifges(6,6) * t463 + t503 + t535 - t538) * t501 + (-t15 * t392 - t16 * t393) * mrSges(7,3) + t6 * t461 + t1 * t444 + t96 * t431 + t95 * t432 + t7 * t229 - t144 * (mrSges(5,1) * t177 + mrSges(5,2) * t176) + Ifges(5,6) * t89 + Ifges(5,5) * t88 + Ifges(6,6) * t56 + Ifges(6,5) * t57 - t44 * mrSges(5,2) + t45 * mrSges(5,1) + t9 * mrSges(6,1) - t10 * mrSges(6,2) - t2 * t445 + t294 * t5 / 0.2e1 + t35 * t362;
t293 = sin(qJ(1));
t406 = t293 * t296;
t247 = qJ(3) * t406;
t262 = pkin(4) * t295 + pkin(3);
t447 = -pkin(2) - t262;
t368 = t292 * t447;
t528 = t293 * t368 + t247;
t527 = -m(4) - m(5);
t11 = -mrSges(7,1) * t33 + mrSges(7,2) * t32;
t48 = mrSges(6,1) * t280 - mrSges(6,3) * t57;
t523 = t11 - t48;
t446 = -mrSges(6,1) * t281 - mrSges(7,1) * t85 + mrSges(7,2) * t86 + mrSges(6,3) * t501;
t216 = -qJ(3) * t291 + t295 * t298;
t151 = t295 * qJD(3) + qJD(4) * t216;
t217 = t295 * qJ(3) + t291 * t298;
t152 = -t291 * qJD(3) - qJD(4) * t217;
t123 = -t209 * t291 + t295 * t210;
t313 = t123 + t419;
t124 = t295 * t209 + t291 * t210;
t87 = t124 + t418;
t502 = (-t152 + t313) * t287 + (t151 - t87) * t286;
t195 = t286 * t295 + t287 * t291;
t518 = t281 * t195;
t266 = Ifges(3,4) * t398;
t437 = Ifges(4,5) * t296;
t338 = t292 * Ifges(4,1) - t437;
t517 = Ifges(3,1) * t399 + qJD(1) * t338 + qJD(2) * t525 + t266;
t377 = mrSges(4,2) * t399;
t516 = mrSges(3,3) * t399 + t377 + (-mrSges(3,1) - mrSges(4,1)) * qJD(2);
t515 = t292 * t524 + t296 * t525;
t346 = t296 * mrSges(4,1) + t292 * mrSges(4,3);
t348 = mrSges(3,1) * t296 - mrSges(3,2) * t292;
t514 = t346 + t348;
t387 = qJ(4) + pkin(10);
t353 = sin(t387);
t354 = cos(t387);
t161 = t292 * t353 + t296 * t354;
t162 = t292 * t354 - t296 * t353;
t355 = mrSges(5,1) * t318 - mrSges(5,2) * t319;
t513 = t161 * t533 + t162 * t534 + t355;
t118 = -t286 * t318 - t287 * t319;
t127 = -qJD(2) * t319 + t307;
t128 = qJD(2) * t318 - t306;
t70 = t127 * t286 + t128 * t287;
t315 = t118 * t392 + t290 * t70;
t512 = -t191 * t296 + t192 * t292;
t153 = -qJDD(2) * pkin(2) + t361;
t511 = t145 * t296 + t153 * t292;
t509 = -mrSges(2,1) - t514;
t508 = m(7) * pkin(9) - t534;
t507 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t505 = m(5) * pkin(8) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + mrSges(6,3);
t62 = -mrSges(7,2) * t97 + mrSges(7,3) * t85;
t63 = mrSges(7,1) * t97 - mrSges(7,3) * t86;
t329 = -t290 * t63 + t294 * t62;
t90 = -mrSges(6,2) * t281 - mrSges(6,3) * t107;
t317 = -t329 - t90;
t500 = t151 - t124;
t499 = t152 - t123;
t471 = pkin(7) - pkin(8);
t239 = t471 * t292;
t240 = t471 * t296;
t133 = t291 * t239 + t295 * t240;
t297 = cos(qJ(1));
t278 = t297 * pkin(7);
t289 = -qJ(5) - pkin(8);
t456 = pkin(4) * t291;
t497 = t293 * (-pkin(1) + t447 * t296 + (-qJ(3) - t456) * t292) + t297 * t289 + t278;
t13 = mrSges(7,1) * t52 - mrSges(7,3) * t32;
t14 = -mrSges(7,2) * t52 + mrSges(7,3) * t33;
t496 = -t290 * t13 + t294 * t14;
t494 = g(1) * t297 + g(2) * t293;
t493 = qJD(3) + t268;
t139 = t162 * t293;
t140 = t161 * t293;
t155 = t319 * t293;
t156 = t318 * t293;
t490 = t155 * mrSges(5,1) + t156 * mrSges(5,2) - t139 * t533 + t140 * t534;
t324 = t297 * t353;
t325 = t297 * t354;
t141 = -t292 * t324 - t296 * t325;
t142 = -t292 * t325 + t296 * t324;
t407 = t292 * t297;
t376 = t295 * t407;
t404 = t296 * t297;
t378 = t291 * t404;
t157 = -t376 + t378;
t158 = t318 * t297;
t489 = t157 * mrSges(5,1) + t158 * mrSges(5,2) - t141 * t534 + t142 * t533;
t301 = t1 * t294 - t2 * t290 + (-t15 * t294 - t16 * t290) * qJD(6);
t488 = m(7) * t301 - t63 * t392 - t62 * t393 + t496;
t460 = t292 / 0.2e1;
t459 = pkin(4) * t177;
t458 = pkin(4) * t286;
t457 = pkin(4) * t287;
t455 = pkin(7) * t292;
t454 = pkin(7) * t296;
t453 = pkin(9) * t162;
t277 = t296 * pkin(2);
t442 = Ifges(3,4) * t292;
t441 = Ifges(3,4) * t296;
t438 = Ifges(4,5) * t292;
t422 = t296 * mrSges(4,3);
t417 = t118 * t290;
t416 = t118 * t294;
t274 = t292 * qJ(3);
t208 = -pkin(4) + t216;
t126 = t286 * t208 + t287 * t217;
t396 = qJD(2) * t296;
t402 = qJ(3) * t396 + t273;
t401 = t277 + t274;
t400 = t297 * pkin(1) + t293 * pkin(7);
t397 = qJD(2) * t292;
t386 = Ifges(7,5) * t32 + Ifges(7,6) * t33 + Ifges(7,3) * t52;
t251 = pkin(4) * t410;
t383 = m(6) + m(7) - t527;
t382 = t298 * t292;
t381 = mrSges(4,2) * t398;
t372 = t296 * pkin(3) + t401;
t369 = -t56 * mrSges(6,1) + t57 * mrSges(6,2);
t360 = -pkin(1) - t274;
t235 = t406 * t456;
t359 = pkin(9) * t140 - t235;
t238 = pkin(4) * t378;
t358 = -pkin(9) * t141 - t238;
t357 = -t391 / 0.2e1;
t132 = t295 * t239 - t240 * t291;
t187 = pkin(1) + t372;
t159 = -qJDD(2) * mrSges(4,1) + t214 * mrSges(4,2);
t352 = t296 * t262 + t251 + t401;
t351 = pkin(2) * t404 + qJ(3) * t407 + t400;
t350 = -pkin(4) * t405 - t251;
t349 = g(1) * t293 - g(2) * t297;
t347 = mrSges(3,1) * t292 + mrSges(3,2) * t296;
t336 = t296 * Ifges(3,2) + t442;
t332 = Ifges(5,5) * t176 - Ifges(5,6) * t177;
t330 = -t15 * t290 + t16 * t294;
t117 = -t286 * t319 + t287 * t318;
t129 = pkin(4) * t318 + t187;
t61 = pkin(5) * t117 - pkin(9) * t118 + t129;
t100 = -qJ(5) * t318 + t133;
t312 = qJ(5) * t319 + t132;
t65 = t287 * t100 + t286 * t312;
t25 = -t290 * t65 + t294 * t61;
t26 = t290 * t61 + t294 * t65;
t135 = -mrSges(5,2) * t281 + t432;
t136 = mrSges(5,1) * t281 - t431;
t323 = t135 * t295 - t136 * t291;
t322 = -t140 * t294 - t290 * t297;
t321 = t140 * t290 - t294 * t297;
t125 = t208 * t287 - t217 * t286;
t222 = -qJD(2) * pkin(2) + t493;
t226 = t269 + t284;
t320 = t222 * t296 - t226 * t292;
t194 = t286 * t291 - t287 * t295;
t260 = qJ(3) * t398;
t154 = qJD(1) * t382 + t260;
t316 = pkin(1) * t347;
t314 = t118 * t393 - t294 * t70;
t310 = t179 * (t292 * mrSges(4,1) - t422);
t309 = t292 * (Ifges(3,1) * t296 - t442);
t308 = t296 * (Ifges(4,3) * t292 + t437);
t137 = t292 * t371 + t402;
t211 = t471 * t397;
t212 = qJD(2) * t240;
t77 = -t295 * t211 + t291 * t212 + t239 * t394 - t240 * t395;
t304 = t297 * t251 + t262 * t404 + t293 * t289 + t351;
t303 = t422 + (-m(4) * pkin(2) - mrSges(4,1)) * t292;
t113 = t154 - t459;
t81 = -pkin(4) * t127 + t137;
t78 = -qJD(4) * t133 + t211 * t291 + t295 * t212;
t300 = -qJ(5) * t128 + qJD(5) * t319 + t78;
t265 = Ifges(4,5) * t399;
t255 = -pkin(5) - t457;
t248 = qJ(3) * t404;
t237 = pkin(4) * t376;
t234 = t293 * pkin(4) * t408;
t228 = qJD(2) * mrSges(4,3) + t381;
t227 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t398;
t223 = -pkin(1) - t401;
t203 = pkin(2) * t399 - t260;
t202 = t346 * qJD(1);
t173 = Ifges(3,6) * qJD(2) + qJD(1) * t336;
t172 = Ifges(4,6) * qJD(2) - Ifges(4,3) * t398 + t265;
t170 = t194 * qJD(4);
t169 = t194 * qJD(2);
t165 = pkin(2) * t397 - t402;
t160 = -mrSges(4,2) * t213 + qJDD(2) * mrSges(4,3);
t131 = -t169 * t294 + t290 * t399;
t130 = t169 * t290 + t294 * t399;
t121 = pkin(5) - t125;
t120 = -t141 * t294 - t290 * t293;
t119 = t141 * t290 - t293 * t294;
t112 = -mrSges(5,1) * t176 + mrSges(5,2) * t177;
t80 = -mrSges(5,2) * t280 + mrSges(5,3) * t89;
t79 = mrSges(5,1) * t280 - mrSges(5,3) * t88;
t69 = -t287 * t127 + t128 * t286;
t68 = mrSges(6,1) * t107 + mrSges(6,2) * t501;
t60 = t459 - t530;
t55 = t286 * t313 + t287 * t87;
t53 = qJ(5) * t127 - qJD(5) * t318 + t77;
t47 = -mrSges(6,2) * t280 + mrSges(6,3) * t56;
t46 = t113 + t530;
t42 = t287 * t75 - t429;
t41 = t286 * t75 + t73;
t24 = pkin(5) * t69 - pkin(9) * t70 + t81;
t22 = t286 * t300 + t287 * t53;
t20 = t290 * t46 + t294 * t55;
t19 = -t290 * t55 + t294 * t46;
t18 = t290 * t60 + t294 * t42;
t17 = -t290 * t42 + t294 * t60;
t4 = -qJD(6) * t26 - t22 * t290 + t24 * t294;
t3 = qJD(6) * t25 + t22 * t294 + t24 * t290;
t21 = [(t222 * t396 - t226 * t397 + t511) * mrSges(4,2) + (-t173 / 0.2e1 + t172 / 0.2e1) * t397 - t296 * (Ifges(4,5) * t214 + Ifges(4,6) * qJDD(2)) / 0.2e1 + (-Ifges(7,5) * t314 - Ifges(7,6) * t315) * t472 + (t292 * Ifges(3,1) + t338 + t441) * t214 / 0.2e1 + t348 * t285 + (-m(6) * t39 + m(7) * t37 + t446) * (t286 * t53 - t287 * t300) + (-m(6) * t9 + m(7) * t7 + t523) * (t100 * t286 - t287 * t312) + (-t1 * t417 + t15 * t314 - t16 * t315 - t2 * t416) * mrSges(7,3) + (t292 * (Ifges(4,1) * t296 + t438) + t296 * (-Ifges(3,2) * t292 + t441) + t309) * t391 / 0.2e1 + (t127 * t96 - t128 * t95) * mrSges(5,3) + (t59 * mrSges(6,2) - t9 * mrSges(6,3) + Ifges(6,1) * t57 + Ifges(6,4) * t56 + Ifges(6,5) * t280 + t331 * t480 + t334 * t483 + t337 * t484 + t7 * t339 + t36 * t362) * t118 + (-t44 * mrSges(5,3) - Ifges(5,4) * t88 - Ifges(5,2) * t89 - Ifges(5,6) * t280) * t318 + (t45 * mrSges(5,3) - Ifges(5,1) * t88 - Ifges(5,4) * t89 - Ifges(5,5) * t280) * t319 + (Ifges(7,3) * t480 + Ifges(7,6) * t483 + Ifges(7,5) * t484 - Ifges(6,4) * t57 - Ifges(6,2) * t56 + t59 * mrSges(6,1) - t10 * mrSges(6,3) + t386 / 0.2e1 - Ifges(6,6) * t280 + t507) * t117 + t37 * (mrSges(7,1) * t315 - mrSges(7,2) * t314) + t296 * (Ifges(3,4) * t214 + Ifges(3,6) * qJDD(2)) / 0.2e1 + (-pkin(1) * t214 - qJDD(2) * t454) * mrSges(3,2) + (t525 * t292 - t524 * t296) * qJDD(2) / 0.2e1 + ((Ifges(3,1) + Ifges(4,1)) * t214 + t525 * qJDD(2)) * t460 + (-t322 * mrSges(7,1) - t321 * mrSges(7,2) - (-pkin(5) * t140 + t497) * m(7) - t497 * m(6) + t140 * mrSges(6,1) + t156 * mrSges(5,1) - t155 * mrSges(5,2) - t508 * t139 + (-m(3) + t527) * t278 + (m(3) * pkin(1) - m(5) * (t296 * t298 + t360) - m(4) * (t360 - t277) - t509) * t293 + t505 * t297) * g(1) + (m(4) * (qJD(2) * t320 + t511) + (-t228 - t227) * t397 + t516 * t396 + m(3) * t512) * pkin(7) + (t310 + t515 * qJD(2) / 0.2e1) * qJD(2) + (Ifges(5,1) * t128 + Ifges(5,4) * t127) * t466 - t223 * mrSges(4,3) * t214 + (Ifges(7,5) * t474 + Ifges(7,6) * t476 + Ifges(7,3) * t472 - t503 + t536) * t69 + (t373 + t506 + t537) * t70 + m(5) * (t132 * t45 + t133 * t44 + t137 * t144 + t187 * t84 + t77 * t96 + t78 * t95) + m(4) * (t111 * t223 + t165 * t179) - t111 * t346 + t84 * t355 + (-m(4) * t351 - m(7) * (-pkin(5) * t141 + t304) - t120 * mrSges(7,1) - t119 * mrSges(7,2) - m(6) * t304 + t141 * mrSges(6,1) - m(3) * t400 - m(5) * (pkin(3) * t404 + t351) - t158 * mrSges(5,1) + t157 * mrSges(5,2) - t508 * t142 + t509 * t297 + t505 * t293) * g(2) + t160 * t454 + (Ifges(5,5) * t128 + Ifges(5,6) * t127) * t462 + (-qJDD(2) * mrSges(3,1) + t159) * t455 + (-Ifges(7,4) * t314 - Ifges(7,2) * t315) * t476 + (t214 * t455 + t512) * mrSges(3,3) - t315 * t35 / 0.2e1 + t517 * t396 / 0.2e1 - t165 * t202 + t187 * (-mrSges(5,1) * t89 + mrSges(5,2) * t88) + t176 * (Ifges(5,4) * t128 + Ifges(5,2) * t127) / 0.2e1 + t144 * (-mrSges(5,1) * t127 + mrSges(5,2) * t128) + t77 * t135 + t78 * t136 + t137 * t112 + t128 * t99 / 0.2e1 + t132 * t79 + t133 * t80 + t127 * t98 / 0.2e1 + t22 * t90 + t81 * t68 + t65 * t47 + t3 * t62 + t4 * t63 + t25 * t13 + t26 * t14 + m(7) * (t1 * t26 + t15 * t4 + t16 * t3 + t2 * t25) + m(6) * (t10 * t65 + t101 * t81 + t129 * t59 + t22 * t40) + (-t336 / 0.2e1 - mrSges(3,3) * t454 - pkin(1) * mrSges(3,1) + t223 * mrSges(4,1) + t438 / 0.2e1 + (Ifges(4,5) - Ifges(3,4)) * t460 + (-Ifges(4,3) - Ifges(3,2) / 0.2e1) * t296) * t213 + (m(3) * pkin(1) ^ 2 + Ifges(2,3)) * qJDD(1) + t6 * t416 / 0.2e1 - t5 * t417 / 0.2e1 - t316 * t391 + (-Ifges(7,1) * t314 - Ifges(7,4) * t315) * t474 + t129 * t369 + t308 * t357; -(Ifges(4,1) * t398 + t172 + t265) * t399 / 0.2e1 + t494 * t347 + t502 * t446 + t524 * t213 + t525 * t214 + (-m(4) * t247 - t293 * t303 - m(7) * (-t359 + t528) - m(6) * (t235 + t528) - m(5) * (t293 * t382 + t247) - t490) * g(2) + (t121 * t7 - t15 * t19 - t16 * t20 + t37 * t502) * m(7) + (t10 * t126 - t101 * t113 + t125 * t9 - t39 * t502 - t40 * t55) * m(6) + (-pkin(2) * t153 + qJ(3) * t145 + qJD(3) * t226 - t179 * t203) * m(4) + (-pkin(7) * t320 * m(4) - t310 + (t308 / 0.2e1 - t309 / 0.2e1 + t316) * qJD(1)) * qJD(1) + t332 * t462 + t226 * t377 + t227 * t268 + t493 * t228 - t529 + t488 * (-pkin(9) + t126) + t499 * t136 + t500 * t135 + (-t144 * t154 + t216 * t45 + t217 * t44 + t499 * t95 + t500 * t96) * m(5) + (-m(4) * t248 - t297 * t303 - m(7) * (t407 * t447 + t248 - t358) - m(6) * (t297 * t368 + t238 + t248) - m(5) * (t297 * t382 + t248) - t489) * g(1) + (-m(6) * t352 - m(7) * (t352 - t453) - m(4) * t401 - m(5) * t372 - t513 - t514) * g(3) + t515 * t357 - t516 * t269 - (-Ifges(3,2) * t399 + t266 + t517) * t398 / 0.2e1 + t216 * t79 + t217 * t80 + t203 * t202 + t191 * mrSges(3,2) - t192 * mrSges(3,1) - pkin(2) * t159 + qJ(3) * t160 - t153 * mrSges(4,1) - t154 * t112 + t145 * mrSges(4,3) + t121 * t11 + t125 * t48 + t126 * t47 - t113 * t68 - t55 * t90 - t20 * t62 - t19 * t63 + (Ifges(4,2) + Ifges(3,3)) * qJDD(2) + t173 * t399 / 0.2e1 + (m(6) * t40 + m(7) * t330 - t317) * (t151 * t287 + t152 * t286) - t222 * t381; t317 * t170 + t159 + (-t228 - t323) * qJD(2) + t323 * qJD(4) + t523 * t194 + t383 * t296 * g(3) + t169 * t90 - t130 * t63 - t131 * t62 + ((-t112 - t202 - t68) * qJD(1) - t494 * t383) * t292 + t291 * t80 + t295 * t79 + (t47 + (-t290 * t62 - t294 * t63) * qJD(6) + t496) * t195 + t518 * t446 + (-t130 * t15 - t131 * t16 - t170 * t330 + t194 * t7 + t195 * t301 + t37 * t518) * m(7) + (t10 * t195 - t101 * t399 - t194 * t9 + (t169 - t170) * t40 - t518 * t39) * m(6) + (-t144 * t399 + t291 * t44 + t295 * t45 + t281 * (-t291 * t95 + t295 * t96)) * m(5) + (-qJD(2) * t226 + t179 * t399 + t153) * m(4); (-m(6) * (-t238 + t237) - m(7) * (t237 + t358) + t489) * g(1) + (-m(6) * (-t235 + t234) - m(7) * (t234 + t359) + t490) * g(2) - t446 * t41 + (-t15 * t17 - t16 * t18 + t255 * t7 - t37 * t41) * m(7) + ((t10 * t286 + t287 * t9) * pkin(4) - t101 * t459 + t39 * t41 - t40 * t42) * m(6) + t332 * t463 + t48 * t457 + t47 * t458 + t488 * (pkin(9) + t458) + (-m(6) * t350 - m(7) * (t350 + t453) + t513) * g(3) + t255 * t11 - t95 * t135 + t96 * t136 - t42 * t90 - t18 * t62 - t17 * t63 - t68 * t459 + t529; t329 * qJD(6) - t446 * t501 - t317 * t107 + t294 * t13 + t290 * t14 + t369 + (t1 * t290 + t2 * t294 + t330 * t97 - t501 * t37 + t349) * m(7) + (t107 * t40 + t39 * t501 + t349 + t59) * m(6); -t37 * (mrSges(7,1) * t86 + mrSges(7,2) * t85) + (Ifges(7,1) * t85 - t448) * t475 + t35 * t474 + (Ifges(7,5) * t85 - Ifges(7,6) * t86) * t473 - t15 * t62 + t16 * t63 - g(1) * (mrSges(7,1) * t119 - mrSges(7,2) * t120) - g(2) * (-mrSges(7,1) * t321 + mrSges(7,2) * t322) + g(3) * t339 * t162 + (t15 * t85 + t16 * t86) * mrSges(7,3) + t386 + (-Ifges(7,2) * t86 + t36 + t83) * t477 + t507;];
tau  = t21;
