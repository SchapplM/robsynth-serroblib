% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR10_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR10_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR10_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:34:29
% EndTime: 2019-03-09 05:35:12
% DurationCPUTime: 28.30s
% Computational Cost: add. (7103->771), mult. (13666->988), div. (0->0), fcn. (8294->8), ass. (0->357)
t544 = Ifges(5,1) + Ifges(6,1);
t494 = Ifges(6,4) + Ifges(5,5);
t543 = Ifges(5,6) - Ifges(6,6);
t242 = sin(qJ(3));
t377 = qJD(1) * t242;
t220 = qJD(4) + t377;
t216 = qJD(6) - t220;
t442 = t216 / 0.2e1;
t241 = sin(qJ(4));
t245 = cos(qJ(4));
t364 = t245 * qJD(3);
t246 = cos(qJ(3));
t376 = qJD(1) * t246;
t180 = t241 * t376 - t364;
t375 = qJD(3) * t241;
t181 = t245 * t376 + t375;
t240 = sin(qJ(6));
t244 = cos(qJ(6));
t280 = t180 * t240 + t181 * t244;
t453 = t280 / 0.2e1;
t97 = t180 * t244 - t181 * t240;
t455 = t97 / 0.2e1;
t547 = Ifges(7,5) * t453 + Ifges(7,6) * t455 + Ifges(7,3) * t442;
t362 = qJD(1) * qJD(3);
t192 = qJDD(1) * t246 - t242 * t362;
t371 = qJD(4) * t180;
t90 = qJDD(3) * t241 + t192 * t245 - t371;
t459 = t90 / 0.2e1;
t91 = qJD(4) * t181 - t245 * qJDD(3) + t192 * t241;
t457 = t91 / 0.2e1;
t193 = -qJDD(1) * t242 - t246 * t362;
t178 = qJDD(4) - t193;
t448 = t178 / 0.2e1;
t446 = t180 / 0.2e1;
t546 = -t181 / 0.2e1;
t545 = -t220 / 0.2e1;
t537 = Ifges(6,2) + Ifges(5,3);
t347 = t241 * t377;
t370 = qJD(4) * t241;
t450 = pkin(8) - pkin(9);
t312 = pkin(3) * t246 + pkin(8) * t242;
t187 = t312 * qJD(1);
t249 = -pkin(1) - pkin(7);
t217 = qJD(1) * t249 + qJD(2);
t386 = t245 * t246;
t105 = t241 * t187 + t217 * t386;
t89 = qJ(5) * t376 + t105;
t542 = -pkin(9) * t347 + t450 * t370 + t89;
t394 = t241 * t246;
t104 = t187 * t245 - t217 * t394;
t208 = t450 * t245;
t451 = pkin(4) + pkin(5);
t353 = t451 * t246;
t392 = t242 * t245;
t269 = pkin(9) * t392 - t353;
t541 = -qJD(1) * t269 + qJD(4) * t208 + t104;
t237 = t246 * pkin(8);
t197 = pkin(3) * t242 + qJ(2) - t237;
t157 = t197 * qJD(1);
t194 = t242 * t217;
t164 = qJD(3) * pkin(8) + t194;
t80 = t245 * t157 - t241 * t164;
t540 = qJD(5) - t80;
t516 = -t241 * t543 + t245 * t494;
t416 = Ifges(6,5) * t241;
t418 = Ifges(5,4) * t241;
t513 = t245 * t544 + t416 - t418;
t458 = -t91 / 0.2e1;
t539 = (-Ifges(5,4) + Ifges(6,5)) * t457 + t544 * t459 + t494 * t448;
t499 = -m(7) - m(6);
t538 = -m(5) + t499;
t535 = -pkin(9) * t181 + t540;
t534 = qJD(5) * t241 + t194;
t211 = qJDD(1) * t249 + qJDD(2);
t373 = qJD(3) * t246;
t118 = t242 * t211 + t217 * t373;
t111 = qJDD(3) * pkin(8) + t118;
t369 = qJD(4) * t245;
t363 = qJD(1) * qJD(2);
t218 = qJDD(1) * qJ(2) + t363;
t94 = -pkin(3) * t193 - pkin(8) * t192 + t218;
t24 = t245 * t111 + t157 * t369 - t164 * t370 + t241 * t94;
t25 = -t241 * t111 - t157 * t370 - t164 * t369 + t245 * t94;
t285 = t24 * t245 - t241 * t25;
t81 = t241 * t157 + t245 * t164;
t533 = -t80 * t369 - t81 * t370 + t285;
t13 = t178 * qJ(5) + t220 * qJD(5) + t24;
t270 = qJDD(5) - t25;
t16 = -pkin(4) * t178 + t270;
t286 = t13 * t245 + t16 * t241;
t64 = -pkin(4) * t220 + t540;
t210 = t220 * qJ(5);
t65 = t210 + t81;
t532 = t64 * t369 - t65 * t370 + t286;
t306 = t245 * mrSges(6,1) + t241 * mrSges(6,3);
t308 = mrSges(5,1) * t245 - mrSges(5,2) * t241;
t531 = -t306 - t308;
t243 = sin(qJ(1));
t388 = t243 * t245;
t247 = cos(qJ(1));
t391 = t242 * t247;
t160 = t241 * t391 + t388;
t383 = t247 * t245;
t389 = t243 * t241;
t161 = t242 * t383 - t389;
t281 = t160 * t240 + t161 * t244;
t480 = -t160 * t244 + t161 * t240;
t530 = mrSges(7,1) * t480 + t281 * mrSges(7,2);
t41 = -t220 * t451 + t535;
t61 = pkin(9) * t180 + t81;
t45 = t210 + t61;
t11 = -t240 * t45 + t244 * t41;
t12 = t240 * t41 + t244 * t45;
t419 = Ifges(4,4) * t246;
t299 = -t242 * Ifges(4,2) + t419;
t529 = t11 * mrSges(7,1) - t12 * mrSges(7,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t299 / 0.2e1 + t537 * t545 + t494 * t546 + t543 * t446 + t547;
t496 = mrSges(5,3) + mrSges(6,2);
t207 = t450 * t241;
t116 = t207 * t240 + t208 * t244;
t528 = -qJD(6) * t116 + t240 * t542 + t244 * t541;
t115 = t207 * t244 - t208 * t240;
t527 = qJD(6) * t115 + t240 * t541 - t244 * t542;
t174 = Ifges(5,4) * t180;
t411 = t180 * Ifges(6,5);
t492 = t181 * t544 + t494 * t220 - t174 + t411;
t354 = t451 * t241;
t402 = qJ(5) * t245;
t276 = t354 - t402;
t523 = -t220 * t276 + t534;
t393 = t242 * t243;
t395 = t241 * t244;
t279 = t240 * t245 - t395;
t478 = qJD(4) - qJD(6);
t251 = t478 * t279;
t433 = pkin(4) * t241;
t288 = -t402 + t433;
t522 = t220 * t288 - t534;
t420 = Ifges(4,4) * t242;
t304 = t246 * Ifges(4,1) - t420;
t173 = Ifges(6,5) * t181;
t71 = t220 * Ifges(6,6) + t180 * Ifges(6,3) + t173;
t521 = Ifges(4,5) * qJD(3) + qJD(1) * t304 + t241 * t71;
t390 = t242 * t249;
t520 = qJD(3) * t312 - qJD(4) * t390 + qJD(2);
t482 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t180 + mrSges(5,2) * t181 + mrSges(4,3) * t376;
t397 = t241 * qJ(5);
t289 = pkin(4) * t245 + t397;
t198 = -pkin(3) - t289;
t519 = m(5) * pkin(3) - m(6) * t198 - t531;
t518 = -t242 * t516 + t246 * t537;
t517 = -t242 * t513 + t246 * t494;
t515 = t241 * t494 + t245 * t543;
t415 = Ifges(6,5) * t245;
t417 = Ifges(5,4) * t245;
t514 = t241 * t544 - t415 + t417;
t310 = mrSges(4,1) * t246 - mrSges(4,2) * t242;
t509 = t246 * (-Ifges(4,1) * t242 - t419) / 0.2e1 + qJ(2) * t310;
t508 = g(1) * t393;
t507 = -t90 * Ifges(6,5) / 0.2e1 - t178 * Ifges(6,6) / 0.2e1 + Ifges(5,4) * t459 + Ifges(5,6) * t448 + (Ifges(6,3) + Ifges(5,2)) * t458;
t54 = mrSges(5,1) * t178 - mrSges(5,3) * t90;
t55 = -t178 * mrSges(6,1) + t90 * mrSges(6,2);
t506 = (-t54 + t55) * t241;
t56 = -mrSges(5,2) * t178 - mrSges(5,3) * t91;
t57 = -mrSges(6,2) * t91 + mrSges(6,3) * t178;
t505 = (t56 + t57) * t245;
t398 = t217 * t246;
t165 = -qJD(3) * pkin(3) - t398;
t264 = qJ(5) * t181 - t165;
t52 = -t180 * t451 + t264;
t504 = -mrSges(7,1) * t52 + mrSges(7,3) * t12;
t503 = mrSges(7,2) * t52 - mrSges(7,3) * t11;
t278 = t240 * t241 + t244 * t245;
t501 = t278 * mrSges(7,1) - t279 * mrSges(7,2);
t170 = qJDD(6) - t178;
t20 = qJD(6) * t97 + t240 * t91 + t244 * t90;
t21 = -qJD(6) * t280 - t240 * t90 + t244 * t91;
t500 = -t20 * Ifges(7,4) / 0.2e1 - t21 * Ifges(7,2) / 0.2e1 - t170 * Ifges(7,6) / 0.2e1;
t466 = t20 / 0.2e1;
t465 = t21 / 0.2e1;
t449 = t170 / 0.2e1;
t498 = t192 / 0.2e1;
t497 = t193 / 0.2e1;
t495 = mrSges(6,3) - mrSges(5,2);
t275 = -t245 * t451 - t397;
t493 = -pkin(3) + t275;
t196 = qJ(5) * t244 - t240 * t451;
t489 = -qJD(6) * t196 - t240 * t535 - t244 * t61;
t195 = -qJ(5) * t240 - t244 * t451;
t488 = qJD(6) * t195 - t240 * t61 + t244 * t535;
t40 = mrSges(5,1) * t91 + mrSges(5,2) * t90;
t487 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t192 - t40;
t143 = t278 * t246;
t486 = t279 * qJD(1) + qJD(3) * t143 + t242 * t251;
t102 = t478 * t278;
t155 = t278 * qJD(1);
t485 = t102 * t242 - t279 * t373 + t155;
t422 = mrSges(5,3) * t180;
t119 = -mrSges(5,2) * t220 - t422;
t424 = mrSges(6,2) * t180;
t122 = mrSges(6,3) * t220 - t424;
t382 = t119 + t122;
t421 = mrSges(5,3) * t181;
t120 = mrSges(5,1) * t220 - t421;
t423 = mrSges(6,2) * t181;
t121 = -mrSges(6,1) * t220 + t423;
t381 = -t120 + t121;
t142 = t240 * t386 - t244 * t394;
t325 = t142 * mrSges(7,1) + t143 * mrSges(7,2);
t483 = t161 * pkin(4) + qJ(5) * t160;
t481 = t178 * t537 + t494 * t90 - t543 * t91;
t374 = qJD(3) * t242;
t117 = t246 * t211 - t217 * t374;
t479 = t117 * t246 + t118 * t242;
t477 = mrSges(3,2) - mrSges(4,3) - mrSges(2,1);
t309 = mrSges(4,1) * t242 + mrSges(4,2) * t246;
t317 = -m(7) * t450 + mrSges(7,3);
t476 = -t246 * t317 + mrSges(2,2) - mrSges(3,3) - t309;
t475 = m(7) * pkin(5) + mrSges(5,1) + mrSges(6,1);
t8 = -pkin(9) * t90 - t178 * t451 + t270;
t9 = pkin(9) * t91 + t13;
t1 = qJD(6) * t11 + t240 * t8 + t244 * t9;
t2 = -qJD(6) * t12 - t240 * t9 + t244 * t8;
t474 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t472 = t25 * mrSges(5,1) - t16 * mrSges(6,1) - t24 * mrSges(5,2) + t13 * mrSges(6,3);
t471 = qJD(1) ^ 2;
t470 = m(4) / 0.2e1;
t469 = m(5) / 0.2e1;
t468 = m(6) / 0.2e1;
t467 = Ifges(7,1) * t466 + Ifges(7,4) * t465 + Ifges(7,5) * t449;
t434 = Ifges(7,4) * t280;
t33 = Ifges(7,2) * t97 + Ifges(7,6) * t216 + t434;
t464 = -t33 / 0.2e1;
t463 = t33 / 0.2e1;
t93 = Ifges(7,4) * t97;
t34 = Ifges(7,1) * t280 + Ifges(7,5) * t216 + t93;
t462 = -t34 / 0.2e1;
t461 = t34 / 0.2e1;
t410 = t181 * Ifges(5,4);
t74 = -t180 * Ifges(5,2) + t220 * Ifges(5,6) + t410;
t460 = -t74 / 0.2e1;
t456 = -t97 / 0.2e1;
t454 = -t280 / 0.2e1;
t452 = -m(3) - m(4);
t447 = -t180 / 0.2e1;
t444 = t181 / 0.2e1;
t443 = -t216 / 0.2e1;
t432 = pkin(5) * t245;
t429 = pkin(9) * t246;
t428 = g(1) * t243;
t427 = g(2) * t247;
t426 = -qJD(1) / 0.2e1;
t425 = qJD(4) / 0.2e1;
t406 = t245 * mrSges(6,3);
t107 = mrSges(6,1) * t180 - mrSges(6,3) * t181;
t42 = -mrSges(7,1) * t97 + mrSges(7,2) * t280;
t405 = t107 - t42;
t403 = qJ(5) * t180;
t396 = t241 * t242;
t387 = t243 * t246;
t385 = t246 * t247;
t124 = t241 * t197 + t245 * t390;
t380 = pkin(3) * t387 + pkin(8) * t393;
t379 = t247 * pkin(1) + t243 * qJ(2);
t372 = qJD(3) * t249;
t368 = qJD(4) * t246;
t366 = qJD(5) * t245;
t365 = qJDD(1) * mrSges(3,2);
t360 = pkin(8) * t387;
t359 = Ifges(7,5) * t20 + Ifges(7,6) * t21 + Ifges(7,3) * t170;
t356 = pkin(8) * t370;
t355 = pkin(8) * t369;
t113 = t242 * qJ(5) + t124;
t348 = t247 * pkin(7) + t379;
t346 = t241 * t374;
t344 = t246 * t372;
t343 = t241 * t368;
t338 = t377 / 0.2e1;
t330 = t369 / 0.2e1;
t329 = -t368 / 0.2e1;
t236 = t247 * qJ(2);
t328 = -pkin(1) * t243 + t236;
t326 = -t362 / 0.2e1;
t323 = (t218 + t363) * qJ(2);
t214 = t241 * t390;
t123 = t197 * t245 - t214;
t318 = pkin(3) * t393 + t348;
t158 = t242 * t389 - t383;
t159 = t241 * t247 + t242 * t388;
t78 = t158 * t244 - t159 * t240;
t79 = t158 * t240 + t159 * t244;
t313 = mrSges(7,1) * t78 - mrSges(7,2) * t79;
t311 = qJDD(3) * pkin(3) + t117;
t307 = mrSges(5,1) * t241 + mrSges(5,2) * t245;
t305 = t241 * mrSges(6,1) - t406;
t298 = -Ifges(5,2) * t241 + t417;
t297 = Ifges(5,2) * t245 + t418;
t294 = -Ifges(4,5) * t242 - Ifges(4,6) * t246;
t291 = Ifges(6,3) * t241 + t415;
t290 = -Ifges(6,3) * t245 + t416;
t66 = -mrSges(7,2) * t216 + mrSges(7,3) * t97;
t67 = mrSges(7,1) * t216 - mrSges(7,3) * t280;
t284 = -t240 * t67 + t244 * t66;
t77 = t214 + (-t197 - t429) * t245 - t451 * t242;
t95 = pkin(9) * t394 + t113;
t36 = -t240 * t95 + t244 * t77;
t37 = t240 * t77 + t244 * t95;
t283 = t241 * t65 - t245 * t64;
t282 = t241 * t81 + t245 * t80;
t59 = -t197 * t370 - t241 * t344 + t245 * t520;
t273 = t242 * (-Ifges(4,2) * t246 - t420);
t271 = t279 * t242;
t267 = t359 - t474;
t266 = t242 * t364 + t343;
t265 = -t245 * t368 + t346;
t263 = t159 * pkin(4) + qJ(5) * t158 + t318;
t58 = t197 * t369 + t241 * t520 + t245 * t344;
t262 = -g(1) * t158 + g(2) * t160 - g(3) * t394;
t261 = -t241 * t382 + t245 * t381;
t224 = pkin(3) * t391;
t260 = -pkin(8) * t385 + t243 * t249 + t224 + t236;
t256 = Ifges(5,6) * t246 - t242 * t298;
t255 = Ifges(6,6) * t246 - t242 * t291;
t253 = qJ(5) * t90 + qJD(5) * t181 + t311;
t44 = qJ(5) * t373 + t242 * qJD(5) + t58;
t229 = -qJDD(1) * pkin(1) + qJDD(2);
t219 = qJ(5) * t386;
t203 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t377;
t186 = t309 * qJD(1);
t171 = t307 * t246;
t145 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t193;
t141 = t278 * t242;
t140 = t240 * t392 - t242 * t395;
t128 = -t219 + (-t249 + t433) * t246;
t127 = t242 * t155;
t126 = qJD(1) * t271;
t114 = -pkin(4) * t242 - t123;
t112 = t219 + (t249 - t354) * t246;
t106 = pkin(4) * t181 + t403;
t92 = -pkin(4) * t376 - t104;
t70 = -t181 * t451 - t403;
t69 = pkin(4) * t180 - t264;
t63 = (qJD(4) * t289 - t366) * t246 + (t249 - t288) * t374;
t53 = -pkin(4) * t373 - t59;
t49 = t246 * t251 - t278 * t374;
t48 = qJD(3) * t271 + t102 * t246;
t43 = (qJD(4) * t275 + t366) * t246 + (-t249 + t276) * t374;
t39 = mrSges(6,1) * t91 - mrSges(6,3) * t90;
t38 = -pkin(9) * t265 + t44;
t35 = pkin(9) * t343 + qJD(3) * t269 - t59;
t17 = pkin(4) * t91 - t253;
t15 = -mrSges(7,2) * t170 + mrSges(7,3) * t21;
t14 = mrSges(7,1) * t170 - mrSges(7,3) * t20;
t10 = -t451 * t91 + t253;
t7 = -mrSges(7,1) * t21 + mrSges(7,2) * t20;
t4 = -qJD(6) * t37 - t240 * t38 + t244 * t35;
t3 = qJD(6) * t36 + t240 * t35 + t244 * t38;
t5 = [m(4) * (t249 * t479 + t323) - t479 * mrSges(4,3) + (t481 / 0.2e1 + t482 * t372 - Ifges(7,5) * t466 + Ifges(5,6) * t458 + Ifges(6,6) * t457 - Ifges(7,6) * t465 - Ifges(7,3) * t449 + t494 * t459 + t537 * t448 + t472 + t474) * t242 + (t80 * mrSges(5,1) - t64 * mrSges(6,1) - t81 * mrSges(5,2) + t65 * mrSges(6,3) - t529 - t547) * t373 + m(5) * (t165 * t249 * t374 + t123 * t25 + t124 * t24 + t58 * t81 + t59 * t80) + (-m(5) * t260 - m(6) * (t260 + t483) - m(7) * (t224 + t328 + t483) - t281 * mrSges(7,1) + t480 * mrSges(7,2) - m(4) * t236 - m(3) * t328 + t496 * t385 - t475 * t161 - t495 * t160 + (-m(4) * t249 + m(7) * pkin(7) - t477) * t243 + t476 * t247) * g(1) + (t16 * t386 + t265 * t65 - t266 * t64) * mrSges(6,2) - pkin(1) * t365 + (t346 / 0.2e1 + t245 * t329) * t74 + t299 * t497 + t304 * t498 + t142 * t500 + t492 * t241 * t329 - (Ifges(4,4) * t192 + Ifges(4,2) * t193 + 0.2e1 * Ifges(4,6) * qJDD(3) + t359) * t242 / 0.2e1 - t311 * t171 + (t309 + 0.2e1 * mrSges(3,3)) * t218 + t165 * (-mrSges(5,1) * t265 - mrSges(5,2) * t266) + t69 * (-mrSges(6,1) * t265 + mrSges(6,3) * t266) + qJD(3) ^ 2 * t294 / 0.2e1 + (Ifges(7,4) * t49 + Ifges(7,2) * t48) * t455 + (Ifges(7,1) * t49 + Ifges(7,4) * t48) * t453 + (Ifges(7,1) * t143 - Ifges(7,4) * t142) * t466 + (Ifges(7,5) * t49 + Ifges(7,6) * t48) * t442 + (Ifges(7,5) * t143 - Ifges(7,6) * t142) * t449 - (t245 * t492 + t521) * t374 / 0.2e1 + (Ifges(4,5) * qJDD(3) + t17 * t305 + Ifges(4,1) * t498 + Ifges(4,4) * t497 + t71 * t330 + t291 * t457 + t298 * t458 + t513 * t459 + t516 * t448 + (m(5) * t311 + t487) * t249) * t246 + (qJD(3) * t517 - t368 * t514) * t444 + (qJD(3) * t518 - t368 * t515) * t220 / 0.2e1 + (-t13 * mrSges(6,2) - t24 * mrSges(5,3) - t507) * t394 + t509 * t362 + (Ifges(7,4) * t143 - Ifges(7,2) * t142) * t465 + m(7) * (t1 * t37 + t10 * t112 + t11 * t4 + t12 * t3 + t2 * t36 + t43 * t52) + m(6) * (t113 * t13 + t114 * t16 + t128 * t17 + t44 * t65 + t53 * t64 + t63 * t69) + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-t25 * t386 + t265 * t81 + t266 * t80) * mrSges(5,3) + (-t1 * t142 - t11 * t49 + t12 * t48 - t143 * t2) * mrSges(7,3) + t10 * t325 + m(3) * (-pkin(1) * t229 + t323) + t229 * mrSges(3,2) + qJ(2) * (-mrSges(4,1) * t193 + mrSges(4,2) * t192) + qJD(2) * t186 + t145 * t390 + (-m(5) * (t318 - t360) - m(6) * (t263 - t360) - m(3) * t379 - m(7) * t263 - t79 * mrSges(7,1) - t78 * mrSges(7,2) - m(4) * t348 + t496 * t387 - t475 * t159 - t495 * t158 + t477 * t247 + t476 * t243) * g(2) + t203 * t344 + t273 * t326 + t386 * t539 + (qJD(3) * t255 - t290 * t368) * t446 + (qJD(3) * t256 - t297 * t368) * t447 + t49 * t461 + t48 * t463 + t143 * t467 + t36 * t14 + t37 * t15 + t43 * t42 + t52 * (-mrSges(7,1) * t48 + mrSges(7,2) * t49) + t3 * t66 + t4 * t67 + t63 * t107 + t112 * t7 + t113 * t57 + t114 * t55 + t58 * t119 + t59 * t120 + t53 * t121 + t44 * t122 + t123 * t54 + t124 * t56 + t128 * t39; t365 - t140 * t14 + t141 * t15 + t485 * t67 + t486 * t66 + (qJ(2) * t452 - mrSges(3,3)) * t471 + (-t39 + t7 + (t241 * t381 + t245 * t382 + t203) * qJD(3) + t487) * t246 + (t145 + t505 + t506 + t261 * qJD(4) + (t405 + t482) * qJD(3)) * t242 + m(3) * t229 + 0.2e1 * (t117 * t470 + (t364 * t65 + t375 * t64 - t17) * t468 + (t364 * t81 - t375 * t80 + t311) * t469) * t246 + 0.2e1 * (t118 * t470 + (qJD(3) * t69 + t532) * t468 + (qJD(3) * t165 + t533) * t469) * t242 + (t427 - t428) * (-t452 - t538) + (t1 * t141 + t10 * t246 + t11 * t485 + t12 * t486 - t140 * t2 - t374 * t52) * m(7) + (-m(5) * t282 - m(6) * t283 - t186 + t261) * qJD(1); (t355 - t92) * t121 + pkin(8) * t505 + pkin(8) * t506 + (t429 * m(7) + t141 * mrSges(7,1) - t140 * mrSges(7,2) + t309 + (-m(7) * t493 + t519) * t242 + t538 * t237 + (mrSges(7,3) - t496) * t246) * g(3) + (t310 + (-m(7) * (t198 - t432) + t501 + t519) * t246 + (-t317 + (m(5) + m(6)) * pkin(8) + t496) * t242) * t427 + (-m(5) * t380 + t499 * (t243 * pkin(4) * t386 + t387 * t397 + t380) + (-(-m(7) * pkin(9) - mrSges(7,3)) * t242 + (-m(7) * t432 - t501 + t531) * t246) * t243) * g(1) + (t338 * t492 + t507) * t245 + (-t81 * (-mrSges(5,2) * t246 + mrSges(5,3) * t396) - t65 * (mrSges(6,2) * t396 + mrSges(6,3) * t246) - t64 * (-mrSges(6,1) * t246 - mrSges(6,2) * t392) - t80 * (mrSges(5,1) * t246 + mrSges(5,3) * t392) + (t256 / 0.2e1 - t255 / 0.2e1) * t180) * qJD(1) + (-t355 - t104) * t120 + (t425 * t513 + t426 * t517) * t181 + (-Ifges(7,5) * t454 - Ifges(7,6) * t456 - Ifges(7,3) * t443 + t529) * t376 + (t1 * t116 - t10 * t493 + t11 * t528 + t115 * t2 + t12 * t527 + t52 * t523) * m(7) - t493 * t7 + (t71 / 0.2e1 + t460) * t370 + (-t104 * t80 - t105 * t81 - t165 * t194 + pkin(3) * t311 + (-qJD(4) * t282 + t285) * pkin(8)) * m(5) + t311 * t308 + t278 * t500 + (-Ifges(7,5) * t279 - Ifges(7,6) * t278) * t449 + (-Ifges(7,4) * t279 - Ifges(7,2) * t278) * t465 + (-Ifges(7,1) * t279 - Ifges(7,4) * t278) * t466 - t279 * t467 + (-Ifges(7,1) * t127 + Ifges(7,4) * t126) * t454 + (-t1 * t278 - t11 * t127 - t12 * t126 + t2 * t279) * mrSges(7,3) + (-Ifges(7,4) * t127 + Ifges(7,2) * t126) * t456 + (-Ifges(7,5) * t127 + Ifges(7,6) * t126) * t443 - t52 * (-mrSges(7,1) * t126 - mrSges(7,2) * t127) + (-t298 / 0.2e1 + t291 / 0.2e1) * t371 + (-t508 + t532) * mrSges(6,2) + (-t508 + t533) * mrSges(5,3) + (-t356 - t105) * t119 + t492 * t330 + t527 * t66 + t528 * t67 + t10 * t501 - t310 * t428 + t523 * t42 + t521 * t338 + t522 * t107 + (t17 * t198 + (-qJD(4) * t283 + t286) * pkin(8) - t64 * t92 - t65 * t89 + t522 * t69) * m(6) - (Ifges(7,4) * t453 + Ifges(7,2) * t455 + Ifges(7,6) * t442 + t463 + t504) * t251 + t514 * t459 + t515 * t448 - t482 * t194 + (t273 / 0.2e1 - t509) * t471 + (t165 * t307 + t305 * t69 + t425 * t516 + t426 * t518) * t220 + (-t356 - t89) * t122 + (Ifges(7,1) * t453 + Ifges(7,4) * t455 + Ifges(7,5) * t442 + t461 + t503) * t102 + Ifges(4,3) * qJDD(3) + Ifges(4,5) * t192 + Ifges(4,6) * t193 + t198 * t39 - t17 * t306 + t294 * t326 + t241 * t539 + t290 * t457 + t297 * t458 + t347 * t460 - t127 * t462 + t126 * t464 - t203 * t398 - pkin(3) * t40 + t115 * t14 + t116 * t15 + t117 * mrSges(4,1) - t118 * mrSges(4,2); (-t180 * t544 + t173 - t410 + t71) * t546 + t472 + (t499 * (t160 * pkin(4) - qJ(5) * t161) + t495 * t161 - t475 * t160 + t530) * g(2) + (-t494 * t180 - t181 * t543) * t545 + t481 + (-m(7) * (-t241 * t353 + t219) - t325 - m(6) * t219 - (t406 + (-m(6) * pkin(4) - mrSges(6,1)) * t241) * t246 + t171) * g(3) + t488 * t66 + (t1 * t196 + t11 * t489 + t12 * t488 + t195 * t2 - t52 * t70) * m(7) + t489 * t67 - t267 + (-m(6) * t64 - t381 + t421) * t81 + (-m(6) * t65 - t382 - t422) * t80 - (Ifges(7,1) * t454 + Ifges(7,4) * t456 + Ifges(7,5) * t443 + t462 - t503) * t97 + (Ifges(7,4) * t454 + Ifges(7,2) * t456 + Ifges(7,6) * t443 + t464 - t504) * t280 + (-pkin(4) * t16 + qJ(5) * t13 + qJD(5) * t65 - t106 * t69) * m(6) + t195 * t14 + t196 * t15 - t165 * (mrSges(5,1) * t181 - mrSges(5,2) * t180) - t69 * (mrSges(6,1) * t181 + mrSges(6,3) * t180) + (-Ifges(5,2) * t181 - t174 + t492) * t446 + (t313 + t499 * (-t158 * pkin(4) + qJ(5) * t159) - t495 * t159 + t475 * t158) * g(1) + t65 * t423 + t64 * t424 + t74 * t444 + (Ifges(6,3) * t181 - t411) * t447 - pkin(4) * t55 + qJ(5) * t57 - t70 * t42 - t106 * t107 + qJD(5) * t122; t244 * t14 + t240 * t15 + t405 * t181 + t284 * qJD(6) + (-t122 - t284) * t220 + t55 + (t1 * t240 - t181 * t52 + t2 * t244 + t262 + t216 * (-t11 * t240 + t12 * t244)) * m(7) + (t181 * t69 - t220 * t65 + t16 + t262) * m(6); -t52 * (mrSges(7,1) * t280 + mrSges(7,2) * t97) + (Ifges(7,1) * t97 - t434) * t454 + t33 * t453 + (Ifges(7,5) * t97 - Ifges(7,6) * t280) * t443 - t11 * t66 + t12 * t67 - g(1) * t313 - g(2) * t530 + g(3) * t325 + (t11 * t97 + t12 * t280) * mrSges(7,3) + t267 + (-Ifges(7,2) * t280 + t34 + t93) * t456;];
tau  = t5;
