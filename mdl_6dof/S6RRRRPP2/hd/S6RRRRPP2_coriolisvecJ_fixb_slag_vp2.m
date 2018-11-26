% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:05
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:04:51
% EndTime: 2018-11-23 18:05:02
% DurationCPUTime: 11.20s
% Computational Cost: add. (9201->647), mult. (22791->799), div. (0->0), fcn. (15397->6), ass. (0->301)
t504 = Ifges(7,4) + Ifges(6,5);
t505 = Ifges(6,4) + Ifges(5,5);
t503 = Ifges(7,2) + Ifges(6,3);
t502 = Ifges(5,6) - Ifges(6,6);
t497 = -Ifges(7,5) + t505;
t498 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t501 = -Ifges(7,6) + Ifges(6,6);
t300 = cos(qJ(4));
t507 = t504 * t300;
t297 = sin(qJ(4));
t506 = t504 * t297;
t438 = t297 / 0.2e1;
t298 = sin(qJ(3));
t299 = sin(qJ(2));
t301 = cos(qJ(3));
t302 = cos(qJ(2));
t269 = t298 * t302 + t299 * t301;
t256 = t269 * qJD(1);
t376 = qJD(2) + qJD(3);
t232 = t256 * t300 + t297 * t376;
t500 = t504 * t232;
t485 = t297 * t503 + t507;
t483 = -t297 * t502 + t300 * t505;
t380 = qJD(4) * t297;
t382 = qJD(1) * t302;
t383 = qJD(1) * t299;
t255 = -t298 * t383 + t301 * t382;
t395 = t255 * t297;
t499 = t380 - t395;
t422 = Ifges(5,4) * t297;
t478 = t300 * t498 - t422 + t506;
t231 = t256 * t297 - t300 * t376;
t252 = qJD(4) - t255;
t465 = t231 * t503 + t252 * t501 + t500;
t496 = t504 * t231;
t290 = -pkin(2) * t302 - pkin(1);
t280 = qJD(1) * t290;
t353 = Ifges(4,5) * t376;
t495 = -t353 / 0.2e1 - t280 * mrSges(4,2);
t303 = -pkin(4) - pkin(5);
t168 = -pkin(3) * t255 - pkin(9) * t256 + t280;
t456 = -pkin(8) - pkin(7);
t282 = t456 * t302;
t272 = qJD(1) * t282;
t258 = t301 * t272;
t281 = t456 * t299;
t271 = qJD(1) * t281;
t260 = qJD(2) * pkin(2) + t271;
t213 = t260 * t298 - t258;
t185 = pkin(9) * t376 + t213;
t76 = t168 * t300 - t185 * t297;
t44 = qJ(6) * t232 + t76;
t467 = qJD(5) - t44;
t34 = t252 * t303 + t467;
t352 = Ifges(4,6) * t376;
t244 = t252 * qJ(5);
t77 = t168 * t297 + t185 * t300;
t45 = qJ(6) * t231 + t77;
t40 = t244 + t45;
t466 = qJD(5) - t76;
t50 = -pkin(4) * t252 + t466;
t51 = t244 + t77;
t494 = -t280 * mrSges(4,1) - t76 * mrSges(5,1) + t50 * mrSges(6,1) + t34 * mrSges(7,1) + t77 * mrSges(5,2) - t40 * mrSges(7,2) - t51 * mrSges(6,3) + t352 / 0.2e1;
t423 = Ifges(5,4) * t232;
t117 = -Ifges(5,2) * t231 + Ifges(5,6) * t252 + t423;
t257 = t298 * t272;
t212 = t260 * t301 + t257;
t184 = -pkin(3) * t376 - t212;
t344 = -mrSges(7,1) * t297 + mrSges(7,2) * t300;
t346 = mrSges(6,1) * t297 - mrSges(6,3) * t300;
t348 = mrSges(5,1) * t297 + mrSges(5,2) * t300;
t308 = -qJ(5) * t232 + t184;
t43 = t231 * t303 + qJD(6) - t308;
t68 = pkin(4) * t231 + t308;
t493 = t117 * t438 - t184 * t348 - t344 * t43 - t346 * t68;
t268 = t298 * t299 - t301 * t302;
t223 = t376 * t268;
t204 = t223 * qJD(1);
t134 = -qJD(4) * t231 - t204 * t300;
t135 = qJD(4) * t232 - t204 * t297;
t224 = t376 * t269;
t205 = t224 * qJD(1);
t491 = t134 * t504 + t135 * t503 + t205 * t501;
t363 = qJD(2) * t456;
t351 = qJD(1) * t363;
t261 = t299 * t351;
t324 = t302 * t351;
t124 = qJD(3) * t213 + t261 * t298 - t301 * t324;
t489 = -qJ(5) * t134 - qJD(5) * t232 + t124;
t379 = qJD(4) * t300;
t250 = pkin(4) * t380 - qJ(5) * t379 - qJD(5) * t297;
t488 = -pkin(5) * t499 - t250;
t285 = qJ(6) * t380;
t487 = -qJ(6) * t395 + t285;
t486 = -t300 * t503 + t506;
t484 = t297 * t505 + t300 * t502;
t481 = t497 * t205 + (-Ifges(5,4) + t504) * t135 + t498 * t134;
t230 = Ifges(5,4) * t231;
t458 = t232 * t498 + t252 * t497 - t230 + t496;
t421 = Ifges(5,4) * t300;
t479 = t297 * t498 + t421 - t507;
t123 = qJD(3) * t212 + t301 * t261 + t298 * t324;
t381 = qJD(2) * t299;
t373 = pkin(2) * t381;
t87 = pkin(3) * t205 + pkin(9) * t204 + qJD(1) * t373;
t13 = -t123 * t297 - t168 * t380 - t185 * t379 + t300 * t87;
t10 = -pkin(4) * t205 - t13;
t12 = t123 * t300 + t168 * t379 - t185 * t380 + t297 * t87;
t7 = qJ(5) * t205 + qJD(5) * t252 + t12;
t477 = t10 * t297 + t300 * t7 + t379 * t50;
t327 = Ifges(7,5) * t300 + Ifges(7,6) * t297;
t337 = -Ifges(5,2) * t297 + t421;
t443 = t252 / 0.2e1;
t444 = -t252 / 0.2e1;
t445 = t232 / 0.2e1;
t447 = t231 / 0.2e1;
t448 = -t231 / 0.2e1;
t476 = t327 * t444 + t337 * t448 + t443 * t483 + t445 * t478 + t447 * t485 - t493;
t38 = mrSges(5,1) * t135 + mrSges(5,2) * t134;
t474 = m(5) * t124 + t38;
t208 = pkin(3) * t256 - pkin(9) * t255;
t178 = pkin(2) * t383 + t208;
t216 = t271 * t301 + t257;
t209 = t297 * t216;
t288 = pkin(2) * t298 + pkin(9);
t389 = -qJ(6) + t288;
t267 = t389 * t300;
t413 = pkin(2) * qJD(3);
t371 = t301 * t413;
t350 = -qJD(6) + t371;
t367 = t303 * t256;
t400 = qJ(6) * t255;
t473 = -t209 - (-t178 - t400) * t300 - t367 + qJD(4) * t267 + t297 * t350;
t197 = t297 * t212;
t431 = pkin(9) - qJ(6);
t279 = t431 * t300;
t472 = -t197 - (-t208 - t400) * t300 - t367 + qJD(4) * t279 - qJD(6) * t297;
t103 = t178 * t297 + t216 * t300;
t247 = t256 * qJ(5);
t60 = t247 + t103;
t471 = -t288 * t380 + t300 * t350 + t487 - t60;
t111 = t208 * t297 + t212 * t300;
t69 = t247 + t111;
t470 = -pkin(9) * t380 - qJD(6) * t300 + t487 - t69;
t394 = t255 * t300;
t384 = pkin(4) * t395 - qJ(5) * t394;
t122 = t213 + t384;
t469 = t122 + t488;
t215 = t271 * t298 - t258;
t126 = t215 + t384;
t372 = t298 * t413;
t468 = t126 - t372 + t488;
t464 = -t122 + t250;
t463 = qJ(5) * t224 + qJD(5) * t268;
t462 = t281 * t301 + t282 * t298;
t294 = t297 * qJ(5);
t461 = pkin(4) * t300 + t294;
t410 = t13 * t297;
t411 = t12 * t300;
t460 = m(5) * (-t379 * t76 - t380 * t77 - t410 + t411) + m(6) * (-t380 * t51 + t477);
t234 = t281 * t298 - t282 * t301;
t273 = t299 * t363;
t274 = t302 * t363;
t142 = qJD(3) * t234 + t273 * t298 - t274 * t301;
t457 = m(4) / 0.2e1;
t455 = pkin(1) * mrSges(3,1);
t454 = pkin(1) * mrSges(3,2);
t453 = t134 / 0.2e1;
t452 = -t135 / 0.2e1;
t451 = t135 / 0.2e1;
t450 = -t205 / 0.2e1;
t449 = t205 / 0.2e1;
t446 = -t232 / 0.2e1;
t442 = -t255 / 0.2e1;
t441 = -t256 / 0.2e1;
t439 = -t297 / 0.2e1;
t437 = -t300 / 0.2e1;
t436 = t300 / 0.2e1;
t435 = m(4) * t280;
t433 = pkin(4) * t256;
t53 = -mrSges(6,2) * t135 + mrSges(6,3) * t205;
t58 = -mrSges(5,2) * t205 - mrSges(5,3) * t135;
t430 = t53 + t58;
t55 = mrSges(5,1) * t205 - mrSges(5,3) * t134;
t56 = -t205 * mrSges(6,1) + mrSges(6,2) * t134;
t429 = -t55 + t56;
t428 = mrSges(5,3) * t231;
t427 = mrSges(5,3) * t232;
t426 = Ifges(4,1) * t256;
t425 = Ifges(3,4) * t299;
t248 = Ifges(4,4) * t255;
t424 = Ifges(4,4) * t256;
t414 = Ifges(4,2) * t255;
t409 = t255 * mrSges(4,3);
t408 = t256 * mrSges(4,3);
t407 = t297 * t51;
t406 = t300 * t34;
t405 = Ifges(3,5) * qJD(2);
t404 = Ifges(3,6) * qJD(2);
t402 = qJ(5) * t231;
t401 = qJ(5) * t300;
t399 = qJ(6) * t269;
t398 = qJD(2) * mrSges(3,1);
t397 = qJD(2) * mrSges(3,2);
t396 = t124 * t462;
t391 = t297 * t301;
t390 = t300 * t301;
t388 = -mrSges(4,1) * t376 + mrSges(5,1) * t231 + mrSges(5,2) * t232 + t408;
t155 = -mrSges(6,2) * t231 + mrSges(6,3) * t252;
t156 = mrSges(7,2) * t252 + mrSges(7,3) * t231;
t387 = t155 + t156;
t157 = -mrSges(5,2) * t252 - t428;
t386 = t155 + t157;
t159 = mrSges(5,1) * t252 - t427;
t160 = -mrSges(6,1) * t252 + mrSges(6,2) * t232;
t385 = -t159 + t160;
t211 = pkin(3) * t268 - pkin(9) * t269 + t290;
t137 = t211 * t297 + t234 * t300;
t377 = qJD(5) * t300;
t374 = pkin(3) + t461;
t370 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t369 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t368 = -Ifges(5,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t133 = pkin(3) * t224 + pkin(9) * t223 + t373;
t141 = qJD(3) * t462 + t273 * t301 + t274 * t298;
t365 = t133 * t297 + t141 * t300 + t211 * t379;
t364 = t141 * t297 + t211 * t380 + t234 * t379;
t88 = qJ(5) * t268 + t137;
t289 = -pkin(2) * t301 - pkin(3);
t362 = t405 / 0.2e1;
t361 = -t404 / 0.2e1;
t37 = -mrSges(7,1) * t135 + mrSges(7,2) * t134;
t54 = -mrSges(7,1) * t205 - t134 * mrSges(7,3);
t102 = t178 * t300 - t209;
t110 = t208 * t300 - t197;
t225 = t297 * t234;
t136 = t211 * t300 - t225;
t349 = mrSges(5,1) * t300 - mrSges(5,2) * t297;
t347 = mrSges(6,1) * t300 + mrSges(6,3) * t297;
t345 = mrSges(7,1) * t300 + mrSges(7,2) * t297;
t336 = Ifges(5,2) * t300 + t422;
t326 = Ifges(7,5) * t297 - Ifges(7,6) * t300;
t325 = pkin(4) * t297 - t401;
t263 = t289 - t461;
t323 = qJ(6) * t223 - qJD(6) * t269;
t20 = t133 * t300 - t364;
t322 = t297 * t303 + t401;
t319 = (-t297 * t77 - t300 * t76) * mrSges(5,3);
t19 = -t234 * t380 + t365;
t1 = -qJ(6) * t134 - qJD(6) * t232 + t205 * t303 - t13;
t2 = qJ(6) * t135 + qJD(6) * t231 + t7;
t305 = t13 * mrSges(5,1) - t10 * mrSges(6,1) - t1 * mrSges(7,1) - t12 * mrSges(5,2) + t2 * mrSges(7,2) + t7 * mrSges(6,3);
t112 = t232 * Ifges(7,5) + t231 * Ifges(7,6) - t252 * Ifges(7,3);
t114 = Ifges(5,5) * t232 - Ifges(5,6) * t231 + Ifges(5,3) * t252;
t116 = Ifges(6,4) * t232 + Ifges(6,2) * t252 + Ifges(6,6) * t231;
t16 = pkin(4) * t135 + t489;
t179 = t352 + t414 + t424;
t180 = t248 + t353 + t426;
t30 = Ifges(5,4) * t134 - Ifges(5,2) * t135 + Ifges(5,6) * t205;
t8 = t135 * t303 - t489;
t304 = (t394 * t76 + t395 * t77 + t411) * mrSges(5,3) + t465 * (t380 / 0.2e1 - t395 / 0.2e1) + (t179 + t112) * t256 / 0.2e1 + t476 * qJD(4) + t491 * t437 + (-Ifges(4,2) * t442 - Ifges(7,3) * t443 + Ifges(5,6) * t447 + t501 * t448 + (Ifges(5,3) + Ifges(6,2)) * t444 + t497 * t446 + t494) * t256 + (-t349 - mrSges(4,1)) * t124 + (Ifges(4,1) * t441 + t327 * t443 + t337 * t447 + t483 * t444 + t478 * t446 + t485 * t448 + t493 + t495) * t255 + (t248 + t180) * t442 + (t34 * t394 + t40 * t499) * mrSges(7,3) + t8 * t345 - t16 * t347 + (-t424 + t116 + t114) * t441 + t484 * t449 + t486 * t451 + (t379 / 0.2e1 - t394 / 0.2e1) * t458 + t212 * t409 - t123 * mrSges(4,2) + (-t394 * t50 + t395 * t51 + t477) * mrSges(6,2) + t479 * t453 + t481 * t438 - Ifges(4,5) * t204 - Ifges(4,6) * t205 + t30 * t436 + t326 * t450 + t336 * t452;
t295 = t300 * pkin(5);
t291 = Ifges(3,4) * t382;
t278 = t431 * t297;
t277 = mrSges(3,3) * t382 - t397;
t276 = -mrSges(3,3) * t383 + t398;
t266 = t389 * t297;
t264 = t295 + t374;
t254 = Ifges(3,1) * t383 + t291 + t405;
t253 = t404 + (Ifges(3,2) * t302 + t425) * qJD(1);
t243 = -t263 + t295;
t238 = t250 + t372;
t235 = -mrSges(4,2) * t376 + t409;
t207 = -mrSges(4,1) * t255 + mrSges(4,2) * t256;
t202 = Ifges(6,2) * t205;
t200 = Ifges(5,3) * t205;
t158 = -mrSges(7,1) * t252 - mrSges(7,3) * t232;
t147 = -mrSges(7,1) * t231 + mrSges(7,2) * t232;
t146 = mrSges(6,1) * t231 - mrSges(6,3) * t232;
t145 = pkin(4) * t232 + t402;
t140 = t269 * t325 - t462;
t132 = Ifges(6,4) * t134;
t131 = Ifges(5,5) * t134;
t130 = Ifges(5,6) * t135;
t129 = Ifges(6,6) * t135;
t104 = t269 * t322 + t462;
t95 = -pkin(4) * t268 - t136;
t85 = t232 * t303 - t402;
t71 = -t110 - t433;
t61 = -t102 - t433;
t59 = t297 * t399 + t88;
t57 = mrSges(7,2) * t205 + mrSges(7,3) * t135;
t49 = t225 + (-t211 - t399) * t300 + t303 * t268;
t36 = mrSges(6,1) * t135 - mrSges(6,3) * t134;
t27 = -t325 * t223 + (qJD(4) * t461 - t377) * t269 + t142;
t18 = -t322 * t223 + (t377 + (t300 * t303 - t294) * qJD(4)) * t269 - t142;
t17 = -pkin(4) * t224 - t20;
t14 = t19 + t463;
t4 = t379 * t399 + (-qJD(4) * t234 - t323) * t297 + t365 + t463;
t3 = t269 * t285 + t303 * t224 + (-t133 + t323) * t300 + t364;
t5 = [m(6) * (t10 * t95 + t14 * t51 + t140 * t16 + t17 * t50 + t27 * t68 + t7 * t88) + m(7) * (t1 * t49 + t104 * t8 + t18 * t43 + t2 * t59 + t3 * t34 + t4 * t40) + m(4) * (t123 * t234 + t141 * t213 - t142 * t212 - t396) + m(5) * (t12 * t137 + t13 * t136 + t142 * t184 + t19 * t77 + t20 * t76 - t396) - t462 * t38 + (-t179 / 0.2e1 - t424 / 0.2e1 + t116 / 0.2e1 + t114 / 0.2e1 - t112 / 0.2e1 + (-Ifges(7,6) / 0.2e1 + t368) * t231 + (Ifges(7,3) / 0.2e1 + t369) * t252 + (-Ifges(7,5) / 0.2e1 + t370) * t232 - t213 * mrSges(4,3) - t414 / 0.2e1 - t494) * t224 - (t458 * t436 + t476 + t180 / 0.2e1 + t319 + t426 / 0.2e1 + (t297 * t40 - t406) * mrSges(7,3) + (t300 * t50 - t407) * mrSges(6,2) - t212 * mrSges(4,3) + t465 * t438 + t248 / 0.2e1 - t495) * t223 + t388 * t142 + (-pkin(7) * t277 - t253 / 0.2e1 + t361 + (-0.2e1 * t455 - 0.3e1 / 0.2e1 * t425 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t302) * qJD(1) + (t207 + 0.2e1 * t435 + qJD(1) * (mrSges(4,1) * t268 + mrSges(4,2) * t269)) * pkin(2)) * t381 + (-pkin(7) * t276 + t254 / 0.2e1 + t362 + (-0.2e1 * t454 + 0.3e1 / 0.2e1 * Ifges(3,4) * t302) * qJD(1)) * t302 * qJD(2) + (t204 * t462 - t205 * t234) * mrSges(4,3) + (t337 * t452 + t327 * t450 + t30 * t439 + t16 * t346 + t8 * t344 - Ifges(4,1) * t204 - Ifges(4,4) * t205 + (mrSges(4,3) + t348) * t124 + (-t1 * t300 + t2 * t297) * mrSges(7,3) + (-t12 * t297 - t13 * t300) * mrSges(5,3) + (t10 * t300 - t297 * t7) * mrSges(6,2) + (t336 * t447 + t326 * t443 - t43 * t345 + t68 * t347 + t184 * t349 + t117 * t437 + (t297 * t34 + t300 * t40) * mrSges(7,3) + (t297 * t76 - t300 * t77) * mrSges(5,3) + (-t297 * t50 - t300 * t51) * mrSges(6,2) + t486 * t448 + t484 * t444 + t479 * t446 + t458 * t439) * qJD(4) + t485 * t451 + t483 * t449 + t491 * t438 + t478 * t453 + (qJD(4) * t465 + t481) * t436) * t269 + (t202 / 0.2e1 + t200 / 0.2e1 - t123 * mrSges(4,3) + t132 / 0.2e1 + t131 / 0.2e1 - t130 / 0.2e1 + t129 / 0.2e1 + Ifges(4,4) * t204 + (-Ifges(7,6) + t368) * t135 + (-Ifges(7,5) + t370) * t134 + (Ifges(4,2) + Ifges(7,3) + t369) * t205 + t305) * t268 + t290 * (mrSges(4,1) * t205 - mrSges(4,2) * t204) + t49 * t54 + t59 * t57 + t88 * t53 + t95 * t56 + t104 * t37 + t136 * t55 + t137 * t58 + t140 * t36 + t27 * t146 + t18 * t147 + t14 * t155 + t4 * t156 + t19 * t157 + t3 * t158 + t20 * t159 + t17 * t160 + t141 * t235; 0.2e1 * ((t123 * t298 - t124 * t301) * t457 + ((-t212 * t298 + t213 * t301) * t457 + m(6) * (t390 * t51 + t391 * t50) / 0.2e1 + m(5) * (t184 * t298 + t390 * t77 - t391 * t76) / 0.2e1) * qJD(3)) * pkin(2) + (-t13 * mrSges(5,3) - t1 * mrSges(7,3) + t385 * t371 + (-t51 * mrSges(6,2) - t77 * mrSges(5,3)) * qJD(4)) * t297 + (-t2 * mrSges(7,3) + t386 * t371 + (-t76 * mrSges(5,3) - t34 * mrSges(7,3)) * qJD(4)) * t300 - m(5) * (t102 * t76 + t103 * t77 + t184 * t215) - m(6) * (t126 * t68 + t50 * t61 + t51 * t60) + (-t126 + t238) * t146 + t304 - m(4) * (-t212 * t215 + t213 * t216) + m(6) * (t16 * t263 + t238 * t68) - t388 * t215 + ((t204 * t301 - t205 * t298) * mrSges(4,3) + (t235 * t301 + t298 * t388) * qJD(3)) * pkin(2) + t471 * t156 + t473 * t158 + (t1 * t266 + t2 * t267 + t243 * t8 + t34 * t473 + t40 * t471 + t43 * t468) * m(7) + t474 * t289 + t468 * t147 + ((-qJD(4) * t386 + t429) * t297 + (qJD(4) * t385 + t430) * t300 + t460) * t288 + ((t362 - t254 / 0.2e1 - t291 / 0.2e1 + qJD(1) * t454 + (t276 - t398) * pkin(7)) * t302 + (t361 + t253 / 0.2e1 + (t455 + t425 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t302) * qJD(1) + (t277 + t397) * pkin(7) + (-t207 - t435) * pkin(2)) * t299) * qJD(1) + t213 * t408 - t60 * t155 - t103 * t157 - t102 * t159 - t61 * t160 - t216 * t235 + t243 * t37 + t263 * t36 + t266 * t54 + t267 * t57; (t430 * t300 + t429 * t297 + (-t297 * t386 + t300 * t385) * qJD(4) + t460) * pkin(9) + t469 * t147 + t464 * t146 + (-t388 + t408) * t213 - mrSges(5,3) * t410 + t304 + (-mrSges(6,2) * t407 - mrSges(7,3) * t406 + t319) * qJD(4) + t472 * t158 + t470 * t156 - m(5) * (t110 * t76 + t111 * t77 + t184 * t213) + (-t1 * t297 - t2 * t300) * mrSges(7,3) - t69 * t155 - t111 * t157 - t110 * t159 - t71 * t160 - t212 * t235 + t264 * t37 - t374 * t36 + t278 * t54 + t279 * t57 - t474 * pkin(3) + (t1 * t278 + t2 * t279 + t264 * t8 + t34 * t472 + t40 * t470 + t43 * t469) * m(7) + (-t16 * t374 + t464 * t68 - t50 * t71 - t51 * t69) * m(6); t202 + t200 + (-Ifges(5,2) * t232 - t230 + t458) * t447 + (-t231 * t505 - t232 * t502) * t444 + (t232 * t503 - t496) * t448 + (-t231 * t34 - t232 * t40) * mrSges(7,3) + t132 + t131 - t130 + t129 + (-t385 + t427) * t77 + (-t386 - t428) * t76 + (t231 * t50 + t232 * t51) * mrSges(6,2) + t387 * qJD(5) + (-t231 * t498 - t423 + t465 + t500) * t446 + (-pkin(4) * t10 + qJ(5) * t7 - t145 * t68 + t466 * t51 - t50 * t77) * m(6) + (qJ(5) * t2 + t1 * t303 - t34 * t45 + t40 * t467 - t43 * t85) * m(7) + (t53 + t57) * qJ(5) + t305 - pkin(4) * t56 - Ifges(7,5) * t134 - Ifges(7,6) * t135 - t145 * t146 - t85 * t147 - t44 * t156 - t45 * t158 + Ifges(7,3) * t205 - t43 * (-mrSges(7,1) * t232 - mrSges(7,2) * t231) - t68 * (mrSges(6,1) * t232 + mrSges(6,3) * t231) - t184 * (mrSges(5,1) * t232 - mrSges(5,2) * t231) + t303 * t54 + (-Ifges(7,5) * t231 + Ifges(7,6) * t232) * t443 + t117 * t445; -t387 * t252 + (t146 - t147) * t232 + t54 + t56 + (-t232 * t43 - t252 * t40 + t1) * m(7) + (t232 * t68 - t252 * t51 + t10) * m(6); -t231 * t156 + t232 * t158 + 0.2e1 * (t8 / 0.2e1 + t40 * t448 + t34 * t445) * m(7) + t37;];
tauc  = t5(:);
