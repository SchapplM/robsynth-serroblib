% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:13:06
% EndTime: 2019-03-09 09:13:45
% DurationCPUTime: 25.68s
% Computational Cost: add. (11440->752), mult. (24661->972), div. (0->0), fcn. (17410->12), ass. (0->345)
t490 = Ifges(3,5) + Ifges(4,4);
t489 = Ifges(4,6) - Ifges(3,6);
t269 = -qJD(2) + qJD(5);
t278 = sin(qJ(6));
t282 = cos(qJ(6));
t225 = -mrSges(7,1) * t282 + mrSges(7,2) * t278;
t496 = m(7) * pkin(5) + mrSges(6,1) - t225;
t279 = sin(qJ(5));
t283 = cos(qJ(5));
t274 = sin(pkin(10));
t275 = cos(pkin(10));
t284 = cos(qJ(2));
t380 = qJD(1) * t284;
t280 = sin(qJ(2));
t381 = qJD(1) * t280;
t165 = t274 * t380 - t275 * t381;
t433 = pkin(8) * t165;
t256 = pkin(7) * t381;
t198 = qJ(4) * t381 - t256;
t286 = -pkin(2) - pkin(3);
t352 = t286 * qJD(2);
t143 = qJD(3) + t352 - t198;
t257 = pkin(7) * t380;
t201 = -qJ(4) * t380 + t257;
t272 = qJD(2) * qJ(3);
t176 = t201 + t272;
t86 = t275 * t143 - t176 * t274;
t68 = -qJD(2) * pkin(4) + t433 + t86;
t164 = -t274 * t381 - t275 * t380;
t434 = pkin(8) * t164;
t87 = t274 * t143 + t275 * t176;
t70 = t87 + t434;
t35 = -t279 * t70 + t283 * t68;
t33 = -pkin(5) * t269 - t35;
t307 = t164 * t279 - t283 * t165;
t421 = mrSges(6,3) * t307;
t75 = t269 * t282 - t278 * t307;
t76 = t269 * t278 + t282 * t307;
t425 = -mrSges(6,1) * t269 - mrSges(7,1) * t75 + mrSges(7,2) * t76 + t421;
t495 = -m(7) * t33 - t425;
t468 = t283 * t164 + t279 * t165;
t64 = pkin(5) * t307 - pkin(9) * t468;
t281 = sin(qJ(1));
t238 = t284 * t281 * qJ(3);
t244 = pkin(4) * t275 + pkin(3);
t427 = -pkin(2) - t244;
t349 = t280 * t427;
t392 = t274 * t284;
t366 = pkin(4) * t392;
t494 = t238 + (t349 + t366) * t281;
t181 = -qJD(1) * pkin(1) - pkin(2) * t380 - qJ(3) * t381;
t136 = pkin(3) * t380 + qJD(4) - t181;
t97 = -pkin(4) * t164 + t136;
t491 = t97 * mrSges(6,2);
t268 = -qJDD(2) + qJDD(5);
t372 = qJD(1) * qJD(2);
t348 = t280 * t372;
t370 = t284 * qJDD(1);
t207 = t348 - t370;
t208 = qJDD(1) * t280 + t284 * t372;
t122 = t207 * t275 - t208 * t274;
t123 = t207 * t274 + t208 * t275;
t48 = qJD(5) * t468 + t122 * t279 + t123 * t283;
t28 = qJD(6) * t75 + t268 * t278 + t282 * t48;
t29 = -qJD(6) * t76 + t268 * t282 - t278 * t48;
t11 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t42 = mrSges(6,1) * t268 - mrSges(6,3) * t48;
t488 = t11 - t42;
t209 = -qJ(3) * t274 + t275 * t286;
t197 = -pkin(4) + t209;
t210 = t275 * qJ(3) + t274 * t286;
t121 = t279 * t197 + t283 * t210;
t192 = t274 * t283 + t275 * t279;
t116 = -t198 * t274 + t275 * t201;
t300 = t116 + t434;
t117 = t275 * t198 + t274 * t201;
t77 = t117 - t433;
t469 = t192 * qJD(3) + qJD(5) * t121 - t279 * t77 + t283 * t300;
t92 = Ifges(6,4) * t468;
t487 = t269 * Ifges(6,5);
t486 = t269 * Ifges(6,6);
t485 = t468 * Ifges(6,2);
t484 = t269 * t192;
t255 = Ifges(3,4) * t380;
t412 = Ifges(4,5) * t284;
t325 = Ifges(4,1) * t280 - t412;
t483 = Ifges(3,1) * t381 + qJD(1) * t325 + t490 * qJD(2) + t255;
t360 = mrSges(4,2) * t381;
t482 = mrSges(3,3) * t381 + t360 + (-mrSges(3,1) - mrSges(4,1)) * qJD(2);
t359 = mrSges(4,2) * t380;
t224 = qJD(2) * mrSges(4,3) + t359;
t481 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t380 + t224;
t480 = t489 * t280 + t490 * t284;
t330 = t284 * mrSges(4,1) + t280 * mrSges(4,3);
t332 = t284 * mrSges(3,1) - t280 * mrSges(3,2);
t479 = t330 + t332;
t393 = t274 * t280;
t303 = t275 * t284 + t393;
t304 = -t275 * t280 + t392;
t112 = -t279 * t303 - t283 * t304;
t374 = qJD(6) * t282;
t170 = t304 * qJD(2);
t171 = t303 * qJD(2);
t306 = t279 * t304 - t283 * t303;
t65 = qJD(5) * t306 - t170 * t279 + t171 * t283;
t298 = t112 * t374 + t278 * t65;
t253 = pkin(7) * t370;
t186 = -pkin(7) * t348 + t253;
t187 = t208 * pkin(7);
t478 = t186 * t284 + t187 * t280;
t467 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t142 = t186 + t467;
t342 = qJDD(3) + t187;
t153 = -qJDD(2) * pkin(2) + t342;
t477 = t142 * t284 + t153 * t280;
t94 = qJD(6) - t468;
t61 = Ifges(6,1) * t307 + t487 + t92;
t476 = t35 * mrSges(6,3) - t61 / 0.2e1;
t475 = -mrSges(2,1) - t479;
t362 = m(7) * pkin(9) + mrSges(7,3);
t474 = -mrSges(6,2) + t362;
t30 = t76 * Ifges(7,5) + t75 * Ifges(7,6) + t94 * Ifges(7,3);
t36 = t279 * t68 + t283 * t70;
t416 = Ifges(6,4) * t307;
t60 = t416 + t485 + t486;
t473 = t36 * mrSges(6,3) - t30 / 0.2e1 + t60 / 0.2e1;
t49 = -qJD(5) * t307 + t122 * t283 - t123 * t279;
t261 = t280 * qJD(3);
t273 = qJDD(1) * pkin(1);
t108 = t207 * pkin(2) - t208 * qJ(3) - qJD(1) * t261 - t273;
t74 = -pkin(3) * t207 + qJDD(4) - t108;
t62 = -pkin(4) * t122 + t74;
t12 = -pkin(5) * t49 - pkin(9) * t48 + t62;
t34 = pkin(9) * t269 + t36;
t41 = -pkin(5) * t468 - pkin(9) * t307 + t97;
t15 = -t278 * t34 + t282 * t41;
t377 = qJD(4) * t280;
t93 = -qJ(4) * t208 - qJD(1) * t377 + qJDD(2) * t286 + t342;
t379 = qJD(2) * t280;
t365 = pkin(7) * t379;
t376 = qJD(4) * t284;
t95 = qJ(4) * t207 + t253 + (-t365 - t376) * qJD(1) + t467;
t58 = -t274 * t95 + t275 * t93;
t39 = -qJDD(2) * pkin(4) - pkin(8) * t123 + t58;
t59 = t274 * t93 + t275 * t95;
t40 = pkin(8) * t122 + t59;
t9 = qJD(5) * t35 + t279 * t39 + t283 * t40;
t7 = pkin(9) * t268 + t9;
t1 = qJD(6) * t15 + t12 * t278 + t282 * t7;
t16 = t278 * t41 + t282 * t34;
t2 = -qJD(6) * t16 + t12 * t282 - t278 * t7;
t472 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t471 = m(5) * qJ(4) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + mrSges(6,3);
t470 = t97 * mrSges(6,1) + t15 * mrSges(7,1) - t16 * mrSges(7,2);
t326 = mrSges(7,1) * t278 + mrSges(7,2) * t282;
t296 = t33 * t326;
t54 = -mrSges(7,2) * t94 + mrSges(7,3) * t75;
t55 = mrSges(7,1) * t94 - mrSges(7,3) * t76;
t316 = -t278 * t55 + t282 * t54;
t422 = mrSges(6,3) * t468;
t80 = -mrSges(6,2) * t269 + t422;
t302 = -t316 - t80;
t285 = cos(qJ(1));
t266 = t285 * pkin(7);
t277 = -pkin(8) - qJ(4);
t466 = t281 * (-pkin(1) + t427 * t284 + (-pkin(4) * t274 - qJ(3)) * t280) + t285 * t277 + t266;
t368 = pkin(10) + qJ(5);
t335 = sin(t368);
t336 = cos(t368);
t159 = t280 * t336 - t284 * t335;
t47 = qJDD(6) - t49;
t13 = mrSges(7,1) * t47 - mrSges(7,3) * t28;
t14 = -mrSges(7,2) * t47 + mrSges(7,3) * t29;
t465 = -t278 * t13 + t282 * t14;
t464 = g(1) * t285 + g(2) * t281;
t424 = mrSges(5,3) * t164;
t137 = qJD(2) * mrSges(5,2) + t424;
t423 = mrSges(5,3) * t165;
t138 = -qJD(2) * mrSges(5,1) + t423;
t463 = -t275 * t137 + t274 * t138 - t224;
t10 = -qJD(5) * t36 - t279 * t40 + t283 * t39;
t132 = t159 * t281;
t158 = t280 * t335 + t284 * t336;
t133 = t158 * t281;
t462 = -t133 * mrSges(6,2) + t132 * t496;
t311 = t285 * t335;
t312 = t285 * t336;
t134 = -t280 * t311 - t284 * t312;
t135 = -t280 * t312 + t284 * t311;
t461 = t134 * mrSges(6,2) - t135 * t496;
t460 = -t159 * mrSges(6,2) - t158 * t496;
t288 = t1 * t282 - t2 * t278 + (-t15 * t282 - t16 * t278) * qJD(6);
t375 = qJD(6) * t278;
t459 = m(7) * t288 - t55 * t374 - t54 * t375 + t465;
t439 = -t269 / 0.2e1;
t445 = -t468 / 0.2e1;
t447 = -t94 / 0.2e1;
t449 = -t76 / 0.2e1;
t451 = -t75 / 0.2e1;
t458 = Ifges(7,5) * t449 - Ifges(6,2) * t445 - Ifges(6,6) * t439 + Ifges(7,6) * t451 + Ifges(7,3) * t447 - t470;
t428 = t76 * Ifges(7,4);
t31 = t75 * Ifges(7,2) + t94 * Ifges(7,6) + t428;
t319 = Ifges(7,5) * t282 - Ifges(7,6) * t278;
t414 = Ifges(7,4) * t282;
t321 = -Ifges(7,2) * t278 + t414;
t415 = Ifges(7,4) * t278;
t324 = Ifges(7,1) * t282 - t415;
t73 = Ifges(7,4) * t75;
t32 = t76 * Ifges(7,1) + t94 * Ifges(7,5) + t73;
t405 = t282 * t32;
t419 = mrSges(7,3) * t282;
t420 = mrSges(7,3) * t278;
t438 = t278 / 0.2e1;
t444 = -t307 / 0.2e1;
t457 = -t491 + Ifges(6,1) * t444 + Ifges(6,5) * t439 + t15 * t419 + t16 * t420 + t319 * t447 + t321 * t451 + t324 * t449 - t296 - t405 / 0.2e1 + t31 * t438;
t455 = t28 / 0.2e1;
t454 = t29 / 0.2e1;
t452 = t47 / 0.2e1;
t450 = t75 / 0.2e1;
t448 = t76 / 0.2e1;
t446 = t94 / 0.2e1;
t443 = t307 / 0.2e1;
t442 = t164 / 0.2e1;
t437 = t280 / 0.2e1;
t436 = pkin(7) * t280;
t435 = pkin(7) * t284;
t265 = t284 * pkin(2);
t426 = pkin(7) - qJ(4);
t418 = Ifges(3,4) * t280;
t417 = Ifges(3,4) * t284;
t413 = Ifges(4,5) * t280;
t411 = t164 * Ifges(5,4);
t402 = t284 * mrSges(4,3);
t400 = t112 * t278;
t399 = t112 * t282;
t262 = t280 * qJ(3);
t389 = t280 * t285;
t388 = t284 * t285;
t160 = -t379 * t426 - t376;
t228 = t426 * t284;
t162 = qJD(2) * t228 - t377;
t89 = t275 * t160 + t274 * t162;
t226 = t426 * t280;
t128 = t274 * t226 + t275 * t228;
t239 = qJ(3) * t388;
t386 = t285 * t366 + t239;
t378 = qJD(2) * t284;
t384 = qJ(3) * t378 + t261;
t383 = t265 + t262;
t382 = t285 * pkin(1) + t281 * pkin(7);
t373 = m(5) + m(6) + m(7);
t367 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t47;
t237 = pkin(4) * t393;
t363 = m(4) + t373;
t361 = t280 * t286;
t354 = t405 / 0.2e1;
t353 = t284 * pkin(3) + t383;
t350 = -t49 * mrSges(6,1) + t48 * mrSges(6,2);
t343 = -t375 / 0.2e1;
t341 = -pkin(1) - t262;
t340 = -t372 / 0.2e1;
t338 = -t122 * mrSges(5,1) + t123 * mrSges(5,2);
t337 = mrSges(5,1) * t303 - mrSges(5,2) * t304;
t88 = -t160 * t274 + t275 * t162;
t127 = t275 * t226 - t228 * t274;
t184 = pkin(1) + t353;
t156 = -qJDD(2) * mrSges(4,1) + t208 * mrSges(4,2);
t334 = t284 * t244 + t237 + t383;
t333 = pkin(2) * t388 + qJ(3) * t389 + t382;
t331 = mrSges(3,1) * t280 + mrSges(3,2) * t284;
t323 = Ifges(3,2) * t284 + t418;
t318 = -t15 * t278 + t16 * t282;
t317 = -t274 * t86 + t275 * t87;
t124 = pkin(4) * t303 + t184;
t53 = -pkin(5) * t306 - pkin(9) * t112 + t124;
t90 = pkin(8) * t304 + t127;
t91 = -pkin(8) * t303 + t128;
t57 = t279 * t90 + t283 * t91;
t23 = -t278 * t57 + t282 * t53;
t24 = t278 * t53 + t282 * t57;
t56 = t279 * t91 - t283 * t90;
t309 = -t133 * t282 - t278 * t285;
t308 = t133 * t278 - t282 * t285;
t120 = t197 * t283 - t210 * t279;
t215 = -qJD(2) * pkin(2) + qJD(3) + t256;
t222 = t257 + t272;
t305 = t215 * t284 - t222 * t280;
t189 = t274 * t279 - t283 * t275;
t250 = qJ(3) * t380;
t154 = qJD(1) * t361 + t250;
t301 = -pkin(8) * t171 + t88;
t299 = pkin(1) * t331;
t297 = t112 * t375 - t282 * t65;
t295 = t181 * (t280 * mrSges(4,1) - t402);
t294 = t280 * (Ifges(3,1) * t284 - t418);
t293 = t284 * (Ifges(4,3) * t280 + t412);
t131 = t280 * t352 + t384;
t291 = t285 * t237 + t244 * t388 + t281 * t277 + t333;
t290 = t402 + (-m(4) * pkin(2) - mrSges(4,1)) * t280;
t110 = pkin(4) * t165 + t154;
t96 = pkin(4) * t170 + t131;
t5 = Ifges(7,4) * t28 + Ifges(7,2) * t29 + Ifges(7,6) * t47;
t6 = Ifges(7,1) * t28 + Ifges(7,4) * t29 + Ifges(7,5) * t47;
t8 = -pkin(5) * t268 - t10;
t287 = Ifges(6,5) * t48 + Ifges(6,6) * t49 + Ifges(6,3) * t268 + t10 * mrSges(6,1) + t8 * t225 + t1 * t419 + t6 * t438 + (Ifges(7,1) * t278 + t414) * t455 + t282 * t5 / 0.2e1 + (Ifges(7,2) * t282 + t415) * t454 + (Ifges(7,5) * t278 + Ifges(7,6) * t282) * t452 - t9 * mrSges(6,2) + t31 * t343 - t2 * t420 + (-t15 * t374 - t16 * t375) * mrSges(7,3) + (t319 * t446 + t321 * t450 + t324 * t448 + t296 + t354) * qJD(6);
t254 = Ifges(4,5) * t381;
t216 = -pkin(1) - t383;
t200 = pkin(2) * t381 - t250;
t199 = t330 * qJD(1);
t178 = Ifges(3,6) * qJD(2) + qJD(1) * t323;
t177 = Ifges(4,6) * qJD(2) - Ifges(4,3) * t380 + t254;
t168 = t189 * qJD(5);
t167 = t189 * qJD(2);
t161 = pkin(2) * t379 - t384;
t157 = -mrSges(4,2) * t207 + qJDD(2) * mrSges(4,3);
t155 = Ifges(5,4) * t165;
t152 = t303 * t285;
t151 = t304 * t285;
t150 = t303 * t281;
t149 = t304 * t281;
t126 = -t167 * t282 + t278 * t381;
t125 = t167 * t278 + t282 * t381;
t119 = -t134 * t282 - t278 * t281;
t118 = t134 * t278 - t281 * t282;
t114 = pkin(5) - t120;
t109 = -mrSges(5,1) * t164 - mrSges(5,2) * t165;
t107 = -qJDD(2) * mrSges(5,1) - mrSges(5,3) * t123;
t106 = qJDD(2) * mrSges(5,2) + mrSges(5,3) * t122;
t99 = -t165 * Ifges(5,1) - Ifges(5,5) * qJD(2) + t411;
t98 = t164 * Ifges(5,2) - Ifges(5,6) * qJD(2) - t155;
t71 = -pkin(8) * t170 + t89;
t66 = qJD(5) * t112 + t283 * t170 + t171 * t279;
t63 = -mrSges(6,1) * t468 + mrSges(6,2) * t307;
t51 = t279 * t300 + t283 * t77;
t44 = t110 - t64;
t43 = -mrSges(6,2) * t268 + mrSges(6,3) * t49;
t25 = pkin(5) * t66 - pkin(9) * t65 + t96;
t21 = -qJD(5) * t56 + t279 * t301 + t283 * t71;
t20 = t278 * t64 + t282 * t35;
t19 = -t278 * t35 + t282 * t64;
t18 = t278 * t44 + t282 * t51;
t17 = -t278 * t51 + t282 * t44;
t4 = -qJD(6) * t24 - t21 * t278 + t25 * t282;
t3 = qJD(6) * t23 + t21 * t282 + t25 * t278;
t22 = [(t150 * mrSges(5,1) - t149 * mrSges(5,2) - t309 * mrSges(7,1) - t308 * mrSges(7,2) - (-pkin(5) * t133 + t466) * m(7) - m(6) * t466 + t133 * mrSges(6,1) - t474 * t132 + (-m(5) - m(4) - m(3)) * t266 + (-m(5) * (t284 * t286 + t341) - m(4) * (t341 - t265) + m(3) * pkin(1) - t475) * t281 + t471 * t285) * g(1) + (-m(3) * t382 - m(4) * t333 - m(7) * (-pkin(5) * t134 + t291) - t119 * mrSges(7,1) - t118 * mrSges(7,2) - m(6) * t291 + t134 * mrSges(6,1) - m(5) * (pkin(3) * t388 + t333) - t152 * mrSges(5,1) + t151 * mrSges(5,2) - t474 * t135 + t475 * t285 + t471 * t281) * g(2) + m(4) * (t108 * t216 + t161 * t181 + (qJD(2) * t305 + t477) * pkin(7)) + (t208 * t436 + t478) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t478) - t298 * t31 / 0.2e1 + t480 * qJD(2) ^ 2 / 0.2e1 - t481 * t365 + (t62 * mrSges(6,2) - t10 * mrSges(6,3) + Ifges(6,1) * t48 + Ifges(6,4) * t49 + Ifges(6,5) * t268 + t319 * t452 + t32 * t343 + t321 * t454 + t324 * t455 + t326 * t8) * t112 + (-mrSges(3,3) * t435 + t413 / 0.2e1 + t216 * mrSges(4,1) - pkin(1) * mrSges(3,1) - t323 / 0.2e1 + (Ifges(4,5) - Ifges(3,4)) * t437 + (-Ifges(4,3) - Ifges(3,2) / 0.2e1) * t284) * t207 + (-Ifges(7,3) * t452 - Ifges(7,6) * t454 - Ifges(7,5) * t455 + Ifges(6,6) * t268 + Ifges(6,2) * t49 - t62 * mrSges(6,1) + Ifges(6,4) * t48 - t367 / 0.2e1 + t9 * mrSges(6,3) + t472) * t306 + (-qJDD(2) * mrSges(3,1) + t156) * t436 + qJD(2) * t295 + (-Ifges(7,5) * t297 - Ifges(7,6) * t298) * t446 + (t280 * (Ifges(4,1) * t284 + t413) + t284 * (-Ifges(3,2) * t280 + t417) + t294) * t372 / 0.2e1 + (Ifges(3,1) * t280 + t325 + t417) * t208 / 0.2e1 + (-t1 * t400 + t15 * t297 - t16 * t298 - t2 * t399) * mrSges(7,3) + (-t178 / 0.2e1 + t177 / 0.2e1) * t379 + (-t222 * t379 + t477) * mrSges(4,2) + (t215 * mrSges(4,2) + t483 / 0.2e1 + t482 * pkin(7)) * t378 + t157 * t435 + (-Ifges(7,4) * t297 - Ifges(7,2) * t298) * t450 + m(5) * (t127 * t58 + t128 * t59 + t131 * t136 + t184 * t74 + t86 * t88 + t87 * t89) - t216 * mrSges(4,3) * t208 + (-mrSges(5,3) * t59 - Ifges(5,4) * t123 - Ifges(5,2) * t122 + Ifges(5,6) * qJDD(2)) * t303 + (mrSges(5,3) * t58 - Ifges(5,1) * t123 - Ifges(5,4) * t122 + Ifges(5,5) * qJDD(2)) * t304 + t332 * t273 + (-Ifges(6,4) * t443 - t486 / 0.2e1 - t485 / 0.2e1 + Ifges(7,3) * t446 + Ifges(7,5) * t448 + Ifges(7,6) * t450 + t470 - t473) * t66 + (-t170 * t87 - t171 * t86) * mrSges(5,3) - t284 * (Ifges(4,5) * t208 + Ifges(4,6) * qJDD(2)) / 0.2e1 - t299 * t372 + t124 * t350 + t74 * t337 + t184 * t338 + m(7) * (t1 * t24 + t15 * t4 + t16 * t3 + t2 * t23) + m(6) * (t124 * t62 + t21 * t36 + t57 * t9 + t96 * t97) + (-pkin(1) * t208 - qJDD(2) * t435) * mrSges(3,2) + t6 * t399 / 0.2e1 - t5 * t400 / 0.2e1 + Ifges(2,3) * qJDD(1) - t161 * t199 - t170 * t98 / 0.2e1 + t171 * t99 / 0.2e1 + t89 * t137 + t88 * t138 + t131 * t109 + t127 * t107 + t128 * t106 + t96 * t63 + t21 * t80 + t4 * t55 + t57 * t43 + t3 * t54 + t23 * t13 + t24 * t14 - t108 * t330 + t284 * (Ifges(3,4) * t208 + Ifges(3,6) * qJDD(2)) / 0.2e1 + (Ifges(5,4) * t171 - Ifges(5,2) * t170) * t442 - t165 * (Ifges(5,1) * t171 - Ifges(5,4) * t170) / 0.2e1 + t136 * (mrSges(5,1) * t170 + mrSges(5,2) * t171) - qJD(2) * (Ifges(5,5) * t171 - Ifges(5,6) * t170) / 0.2e1 + (-m(6) * t10 + m(7) * t8 + t488) * t56 + ((Ifges(3,1) + Ifges(4,1)) * t208 + t490 * qJDD(2)) * t437 + (t490 * t280 - t489 * t284) * qJDD(2) / 0.2e1 + (Ifges(6,1) * t443 + t487 / 0.2e1 + t92 / 0.2e1 + t491 + t354 - t476) * t65 + (-m(6) * t35 - t495) * (qJD(5) * t57 + t279 * t71 - t283 * t301) + t293 * t340 + (-Ifges(7,1) * t297 - Ifges(7,4) * t298) * t448 + t33 * (mrSges(7,1) * t298 - mrSges(7,2) * t297); (t114 * t8 - t15 * t17 - t16 * t18 + t33 * t469) * m(7) + (t10 * t120 - t110 * t97 + t121 * t9 - t35 * t469 - t36 * t51) * m(6) + (-m(5) * (t285 * t361 + t239) - t151 * mrSges(5,1) - t152 * mrSges(5,2) - m(7) * (pkin(9) * t134 + t389 * t427 + t386) - t134 * mrSges(7,3) - m(4) * t239 - t290 * t285 - m(6) * (t285 * t349 + t386) + t461) * g(1) + (-m(4) * t383 - m(5) * t353 - t337 - m(7) * (-pkin(9) * t159 + t334) + t159 * mrSges(7,3) - m(6) * t334 + t460 - t479) * g(3) + t480 * t340 + t481 * t256 + (m(6) * t36 + m(7) * t318 - t302) * (-qJD(3) * t189 + qJD(5) * t120) - (-Ifges(6,4) * t444 + t458 + t473) * t307 - t287 + t459 * (-pkin(9) + t121) - t482 * t257 - (-Ifges(3,2) * t381 + t255 + t483) * t380 / 0.2e1 + t99 * t442 + t87 * t423 + t464 * t331 + t222 * t360 + (Ifges(4,2) + Ifges(3,3) + Ifges(5,3)) * qJDD(2) + (-m(5) * (t281 * t361 + t238) - t149 * mrSges(5,1) - t150 * mrSges(5,2) - m(7) * (-pkin(9) * t133 + t494) + t133 * mrSges(7,3) - m(4) * t238 - t290 * t281 - m(6) * t494 + t462) * g(2) + (-pkin(2) * t153 + qJ(3) * t142 - t181 * t200) * m(4) + (-m(4) * t305 * pkin(7) - t295 + (t299 + t293 / 0.2e1 - t294 / 0.2e1) * qJD(1)) * qJD(1) + (-t116 * t86 - t117 * t87 - t136 * t154 + t209 * t58 + t210 * t59) * m(5) - t86 * t424 - (Ifges(4,1) * t380 + t177 + t254) * t381 / 0.2e1 + (-Ifges(5,1) * t164 - t155 + t98) * t165 / 0.2e1 + (m(4) * t222 + m(5) * t317 - t463) * qJD(3) - t164 * (-Ifges(5,2) * t165 - t411) / 0.2e1 - (Ifges(6,4) * t445 + t457 + t476) * t468 + t178 * t381 / 0.2e1 - t215 * t359 + t210 * t106 + t209 * t107 + t200 * t199 - t186 * mrSges(3,2) - t187 * mrSges(3,1) - t136 * (mrSges(5,1) * t165 - mrSges(5,2) * t164) + qJD(2) * (-Ifges(5,5) * t164 - Ifges(5,6) * t165) / 0.2e1 - t153 * mrSges(4,1) - t154 * t109 - pkin(2) * t156 + qJ(3) * t157 - t117 * t137 - t116 * t138 + t142 * mrSges(4,3) + t120 * t42 + t121 * t43 - Ifges(5,6) * t122 - Ifges(5,5) * t123 + t114 * t11 - t110 * t63 - t51 * t80 - t17 * t55 - t58 * mrSges(5,1) + t59 * mrSges(5,2) - t18 * t54 + t469 * t425 + t489 * t207 + t490 * t208; ((-t109 - t199 - t63) * qJD(1) - t464 * t363) * t280 + (t43 + (-t278 * t54 - t282 * t55) * qJD(6) + t465) * t192 + t488 * t189 + t463 * qJD(2) + t156 + t274 * t106 + t275 * t107 + t363 * t284 * g(3) + t167 * t80 - t125 * t55 - t126 * t54 + t302 * t168 + t484 * t425 + (-t125 * t15 - t126 * t16 - t168 * t318 + t189 * t8 + t192 * t288 + t33 * t484) * m(7) + (-t10 * t189 + t192 * t9 - t381 * t97 + (t167 - t168) * t36 - t484 * t35) * m(6) + (-qJD(2) * t317 - t136 * t381 + t274 * t59 + t275 * t58) * m(5) + (-qJD(2) * t222 + t181 * t381 + t153) * m(4); t316 * qJD(6) - t425 * t307 + t302 * t468 + t282 * t13 - t164 * t137 - t165 * t138 + t278 * t14 + t338 + t350 + (g(1) * t281 - g(2) * t285) * t373 + (t1 * t278 + t2 * t282 - t307 * t33 + t318 * t94) * m(7) + (t307 * t35 - t36 * t468 + t62) * m(6) + (-t164 * t87 - t165 * t86 + t74) * m(5); -pkin(5) * t11 - t19 * t55 - t20 * t54 + t60 * t443 + t287 + (t92 + t61) * t445 + (-t416 + t30) * t444 + (t422 - t80) * t35 + (-t159 * t362 - t460) * g(3) + (-t133 * t362 - t462) * g(2) + (t134 * t362 - t461) * g(1) + (-pkin(5) * t8 - t15 * t19 - t16 * t20) * m(7) + (t421 + t495) * t36 + t458 * t307 + t457 * t468 + t459 * pkin(9); -t33 * (mrSges(7,1) * t76 + mrSges(7,2) * t75) + (Ifges(7,1) * t75 - t428) * t449 + t31 * t448 + (Ifges(7,5) * t75 - Ifges(7,6) * t76) * t447 - t15 * t54 + t16 * t55 - g(1) * (mrSges(7,1) * t118 - mrSges(7,2) * t119) - g(2) * (-mrSges(7,1) * t308 + mrSges(7,2) * t309) + g(3) * t326 * t159 + (t15 * t75 + t16 * t76) * mrSges(7,3) + t367 + (-Ifges(7,2) * t76 + t32 + t73) * t451 - t472;];
tau  = t22;
