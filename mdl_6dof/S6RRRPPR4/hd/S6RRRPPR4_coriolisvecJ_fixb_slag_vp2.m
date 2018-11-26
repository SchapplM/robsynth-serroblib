% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2018-11-23 17:35
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:34:38
% EndTime: 2018-11-23 17:34:57
% DurationCPUTime: 19.71s
% Computational Cost: add. (10485->723), mult. (26100->967), div. (0->0), fcn. (17907->8), ass. (0->324)
t290 = sin(qJ(2));
t293 = cos(qJ(2));
t325 = pkin(2) * t290 - pkin(8) * t293;
t245 = t325 * qJD(1);
t292 = cos(qJ(3));
t289 = sin(qJ(3));
t362 = qJD(1) * t290;
t344 = t289 * t362;
t199 = pkin(7) * t344 + t292 * t245;
t367 = t292 * t293;
t306 = pkin(3) * t290 - qJ(4) * t367;
t394 = -qJ(4) - pkin(8);
t333 = qJD(3) * t394;
t502 = -qJD(1) * t306 - qJD(4) * t289 + t292 * t333 - t199;
t227 = t289 * t245;
t354 = qJD(4) * t292;
t368 = t290 * t292;
t369 = t289 * t293;
t501 = t227 + (-pkin(7) * t368 - qJ(4) * t369) * qJD(1) - t289 * t333 - t354;
t287 = sin(pkin(10));
t375 = cos(pkin(10));
t474 = t287 * t502 - t501 * t375;
t331 = t375 * t289;
t237 = t287 * t292 + t331;
t361 = qJD(1) * t293;
t204 = t237 * t361;
t223 = t237 * qJD(3);
t366 = t204 - t223;
t330 = t375 * t292;
t303 = -t287 * t289 + t330;
t301 = t293 * t303;
t205 = qJD(1) * t301;
t224 = t303 * qJD(3);
t365 = t205 - t224;
t475 = -qJ(5) * t362 + t474;
t359 = qJD(2) * t292;
t243 = -t344 + t359;
t343 = t292 * t362;
t244 = qJD(2) * t289 + t343;
t177 = -t375 * t243 + t244 * t287;
t288 = sin(qJ(6));
t291 = cos(qJ(6));
t304 = t287 * t243 + t244 * t375;
t468 = t177 * t288 + t291 * t304;
t99 = t177 * t291 - t288 * t304;
t92 = Ifges(7,4) * t99;
t500 = Ifges(7,2) * t468 - t92;
t472 = t501 * t287 + t375 * t502;
t360 = qJD(2) * t290;
t334 = qJD(1) * t360;
t353 = qJD(2) * qJD(3);
t356 = qJD(3) * t290;
t358 = qJD(2) * t293;
t452 = t289 * t356 - t292 * t358;
t197 = -qJD(1) * t452 + t292 * t353;
t355 = qJD(3) * t292;
t471 = t289 * t358 + t290 * t355;
t198 = -qJD(1) * t471 - t289 * t353;
t132 = t197 * t287 - t198 * t375;
t133 = t197 * t375 + t287 * t198;
t29 = qJD(6) * t99 + t132 * t288 + t133 * t291;
t30 = -qJD(6) * t468 + t132 * t291 - t133 * t288;
t393 = Ifges(7,5) * t29 + Ifges(7,6) * t30;
t272 = qJD(3) - t361;
t261 = qJD(6) - t272;
t411 = -t261 / 0.2e1;
t251 = -pkin(2) * t293 - t290 * pkin(8) - pkin(1);
t231 = t251 * qJD(1);
t283 = pkin(7) * t361;
t258 = qJD(2) * pkin(8) + t283;
t185 = t231 * t289 + t258 * t292;
t246 = t325 * qJD(2);
t232 = qJD(1) * t246;
t329 = pkin(7) * t334;
t115 = -qJD(3) * t185 + t292 * t232 + t289 * t329;
t62 = pkin(3) * t334 - qJ(4) * t197 - qJD(4) * t244 + t115;
t357 = qJD(3) * t289;
t114 = t231 * t355 + t289 * t232 - t258 * t357 - t292 * t329;
t68 = qJ(4) * t198 + qJD(4) * t243 + t114;
t20 = t287 * t62 + t375 * t68;
t15 = qJ(5) * t334 + t272 * qJD(5) + t20;
t10 = pkin(9) * t132 + t15;
t321 = t287 * t68 - t375 * t62;
t434 = pkin(4) + pkin(5);
t349 = t290 * t434;
t328 = qJD(2) * t349;
t11 = -t133 * pkin(9) - qJD(1) * t328 + t321;
t184 = t292 * t231 - t258 * t289;
t147 = -qJ(4) * t244 + t184;
t134 = pkin(3) * t272 + t147;
t148 = qJ(4) * t243 + t185;
t371 = t287 * t148;
t63 = t134 * t375 - t371;
t305 = qJD(5) - t63;
t465 = pkin(9) * t304;
t40 = -t272 * t434 + t305 - t465;
t482 = pkin(9) * t177;
t332 = t375 * t148;
t64 = t287 * t134 + t332;
t57 = t272 * qJ(5) + t64;
t44 = t57 + t482;
t8 = -t288 * t44 + t291 * t40;
t1 = qJD(6) * t8 + t10 * t291 + t11 * t288;
t9 = t288 * t40 + t291 * t44;
t2 = -qJD(6) * t9 - t10 * t288 + t11 * t291;
t451 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t257 = -qJD(2) * pkin(2) + pkin(7) * t362;
t302 = pkin(3) * t243 - qJD(4) - t257;
t298 = qJ(5) * t304 + t302;
t49 = -t177 * t434 + t298;
t499 = (t468 * t9 + t8 * t99) * mrSges(7,3) - Ifges(7,3) * t334 + (Ifges(7,5) * t99 - Ifges(7,6) * t468) * t411 - t49 * (mrSges(7,1) * t468 + mrSges(7,2) * t99) + t393 + t451;
t427 = -t132 / 0.2e1;
t425 = t133 / 0.2e1;
t464 = Ifges(5,1) + Ifges(6,1);
t463 = Ifges(5,4) - Ifges(6,5);
t462 = Ifges(6,4) + Ifges(5,5);
t498 = pkin(9) * t365 + qJD(1) * t349 - t472;
t497 = -pkin(9) * t366 + t475;
t70 = t375 * t147 - t371;
t455 = qJD(5) - t70;
t402 = pkin(3) * t289;
t234 = t361 * t402 + t283;
t496 = -pkin(3) * t357 + t234;
t387 = Ifges(7,4) * t468;
t493 = Ifges(7,1) * t99 - t387;
t38 = Ifges(7,1) * t468 + t261 * Ifges(7,5) + t92;
t491 = t38 / 0.2e1;
t338 = Ifges(3,5) * qJD(2) / 0.2e1;
t429 = -t468 / 0.2e1;
t490 = t462 * t334 / 0.2e1 + t464 * t425 + t463 * t427;
t461 = Ifges(5,6) - Ifges(6,6);
t485 = -t465 + t455;
t484 = qJ(5) * t365 - qJD(5) * t237 - t496;
t483 = Ifges(4,3) + Ifges(5,3) + Ifges(6,2);
t481 = -qJD(2) / 0.2e1;
t255 = t394 * t289;
t256 = t394 * t292;
t191 = -t375 * t255 - t256 * t287;
t160 = -pkin(9) * t237 + t191;
t192 = t287 * t255 - t375 * t256;
t161 = -pkin(9) * t303 + t192;
t75 = t160 * t288 + t161 * t291;
t480 = -qJD(6) * t75 - t288 * t497 + t291 * t498;
t74 = t160 * t291 - t161 * t288;
t479 = qJD(6) * t74 + t288 * t498 + t291 * t497;
t477 = t366 * t434 - t484;
t460 = -t177 * t463 + t272 * t462 + t304 * t464;
t476 = pkin(4) * t362 - t472;
t473 = -pkin(4) * t366 + t484;
t281 = Ifges(3,4) * t361;
t380 = t244 * Ifges(4,4);
t164 = t243 * Ifges(4,2) + t272 * Ifges(4,6) + t380;
t235 = Ifges(4,4) * t243;
t165 = t244 * Ifges(4,1) + t272 * Ifges(4,5) + t235;
t308 = t184 * t292 + t185 * t289;
t388 = Ifges(4,4) * t292;
t315 = -Ifges(4,2) * t289 + t388;
t389 = Ifges(4,4) * t289;
t317 = Ifges(4,1) * t292 - t389;
t318 = mrSges(4,1) * t289 + mrSges(4,2) * t292;
t385 = Ifges(4,6) * t289;
t386 = Ifges(4,5) * t292;
t406 = t292 / 0.2e1;
t407 = -t289 / 0.2e1;
t408 = t272 / 0.2e1;
t412 = t244 / 0.2e1;
t294 = -t308 * mrSges(4,3) + t257 * t318 + t243 * t315 / 0.2e1 + t317 * t412 + (-t385 + t386) * t408 + t164 * t407 + t165 * t406;
t470 = t294 + Ifges(3,1) * t362 / 0.2e1 + t281 / 0.2e1 + t338;
t56 = -t272 * pkin(4) + t305;
t72 = pkin(4) * t177 - t298;
t469 = -mrSges(5,2) * t302 + t56 * mrSges(6,2) - t63 * mrSges(5,3) - t72 * mrSges(6,3) + t460 / 0.2e1;
t467 = Ifges(6,6) / 0.2e1;
t442 = t29 / 0.2e1;
t441 = t30 / 0.2e1;
t426 = t132 / 0.2e1;
t423 = -t177 / 0.2e1;
t422 = t177 / 0.2e1;
t419 = t304 / 0.2e1;
t466 = -t334 / 0.2e1;
t337 = Ifges(3,6) * t481;
t345 = t375 * pkin(3);
t277 = -t345 - pkin(4);
t271 = -pkin(5) + t277;
t403 = pkin(3) * t287;
t275 = qJ(5) + t403;
t208 = t271 * t291 - t275 * t288;
t69 = t147 * t287 + t332;
t47 = t69 + t482;
t459 = qJD(6) * t208 - t288 * t47 + t291 * t485;
t209 = t271 * t288 + t275 * t291;
t458 = -qJD(6) * t209 - t288 * t485 - t291 * t47;
t150 = mrSges(5,1) * t272 - mrSges(5,3) * t304;
t151 = -mrSges(6,1) * t272 + mrSges(6,2) * t304;
t454 = t150 - t151;
t453 = Ifges(4,5) * t197 + Ifges(4,6) * t198;
t274 = pkin(7) * t367;
t207 = t289 * t251 + t274;
t450 = qJD(1) * pkin(1) * mrSges(3,2);
t449 = -t132 * t461 + t133 * t462 + t334 * t483 + t453;
t18 = -pkin(4) * t334 + t321;
t448 = -t115 * mrSges(4,1) + mrSges(5,1) * t321 + t18 * mrSges(6,1) + t114 * mrSges(4,2) + t20 * mrSges(5,2) - t15 * mrSges(6,3);
t327 = Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(4,3) / 0.2e1;
t350 = t467 - Ifges(5,6) / 0.2e1;
t351 = Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1;
t390 = Ifges(3,4) * t290;
t447 = t350 * t177 + t327 * t272 + t351 * t304 + t184 * mrSges(4,1) + t57 * mrSges(6,3) + t63 * mrSges(5,1) + t9 * mrSges(7,2) + t244 * Ifges(4,5) + t243 * Ifges(4,6) + t337 - (t293 * Ifges(3,2) + t390) * qJD(1) / 0.2e1 - t261 * Ifges(7,3) - t468 * Ifges(7,5) - t99 * Ifges(7,6) + Ifges(5,6) * t423 + Ifges(6,6) * t422 - t185 * mrSges(4,2) - t56 * mrSges(6,1) - t64 * mrSges(5,2) - t8 * mrSges(7,1) + t462 * t419 + t483 * t408;
t83 = Ifges(6,5) * t304 + t272 * Ifges(6,6) + t177 * Ifges(6,3);
t86 = Ifges(5,4) * t304 - t177 * Ifges(5,2) + t272 * Ifges(5,6);
t445 = mrSges(5,1) * t302 - t72 * mrSges(6,1) + mrSges(6,2) * t57 + mrSges(5,3) * t64 + t86 / 0.2e1 - t83 / 0.2e1;
t444 = Ifges(7,4) * t442 + Ifges(7,2) * t441 + Ifges(7,6) * t466;
t443 = Ifges(7,1) * t442 + Ifges(7,4) * t441 + Ifges(7,5) * t466;
t440 = Ifges(6,5) * t425 + Ifges(6,3) * t426 + t334 * t467;
t439 = -t133 * Ifges(5,4) / 0.2e1 + Ifges(5,2) * t426 + Ifges(5,6) * t466;
t436 = -t99 / 0.2e1;
t435 = t99 / 0.2e1;
t433 = pkin(1) * mrSges(3,1);
t428 = t468 / 0.2e1;
t420 = -t304 / 0.2e1;
t418 = t197 / 0.2e1;
t417 = t198 / 0.2e1;
t414 = -t243 / 0.2e1;
t413 = -t244 / 0.2e1;
t410 = t261 / 0.2e1;
t409 = -t272 / 0.2e1;
t404 = pkin(3) * t244;
t401 = pkin(7) * t289;
t364 = t289 * t246 + t251 * t355;
t105 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t368 + (-qJD(4) * t290 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t293) * t289 + t364;
t363 = t292 * t246 + t360 * t401;
t91 = -t290 * t354 + t306 * qJD(2) + (-t274 + (qJ(4) * t290 - t251) * t289) * qJD(3) + t363;
t42 = t375 * t105 + t287 * t91;
t392 = mrSges(4,3) * t243;
t391 = mrSges(4,3) * t244;
t138 = t204 * t291 - t205 * t288;
t172 = t237 * t291 - t288 * t303;
t94 = -qJD(6) * t172 + t223 * t291 - t224 * t288;
t377 = t138 - t94;
t139 = t204 * t288 + t205 * t291;
t171 = -t237 * t288 - t291 * t303;
t93 = qJD(6) * t171 + t223 * t288 + t224 * t291;
t376 = t139 - t93;
t372 = qJD(2) * mrSges(3,2);
t370 = t289 * t290;
t239 = t292 * t251;
t181 = -qJ(4) * t368 + t239 + (-pkin(3) - t401) * t293;
t189 = -qJ(4) * t370 + t207;
t109 = t287 * t181 + t375 * t189;
t247 = pkin(3) * t370 + t290 * pkin(7);
t203 = pkin(3) * t471 + pkin(7) * t358;
t280 = -pkin(3) * t292 - pkin(2);
t59 = t132 * mrSges(5,1) + t133 * mrSges(5,2);
t58 = t132 * mrSges(6,1) - t133 * mrSges(6,3);
t326 = m(4) * t257 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t243 + mrSges(4,2) * t244 + mrSges(3,3) * t362;
t7 = -t30 * mrSges(7,1) + t29 * mrSges(7,2);
t102 = -qJ(5) * t293 + t109;
t215 = t303 * t290;
t322 = qJ(5) * t215 - t247;
t320 = t287 * t105 - t375 * t91;
t319 = mrSges(4,1) * t292 - mrSges(4,2) * t289;
t316 = Ifges(4,1) * t289 + t388;
t314 = Ifges(4,2) * t292 + t389;
t313 = Ifges(4,5) * t289 + Ifges(4,6) * t292;
t312 = -qJ(5) * t177 - t404;
t108 = t181 * t375 - t287 * t189;
t104 = t293 * pkin(4) - t108;
t71 = t293 * pkin(5) - t215 * pkin(9) + t104;
t214 = t237 * t290;
t73 = t214 * pkin(9) + t102;
t31 = -t288 * t73 + t291 * t71;
t32 = t288 * t71 + t291 * t73;
t76 = -mrSges(7,2) * t261 + mrSges(7,3) * t99;
t77 = mrSges(7,1) * t261 - mrSges(7,3) * t468;
t311 = -t288 * t77 + t291 * t76;
t309 = t114 * t292 - t115 * t289;
t153 = t214 * t291 - t215 * t288;
t154 = t214 * t288 + t215 * t291;
t35 = qJ(5) * t360 - qJD(5) * t293 + t42;
t175 = -t198 * pkin(3) + qJD(2) * t283;
t307 = qJ(5) * t237 - t280;
t113 = -mrSges(6,1) * t334 + t133 * mrSges(6,2);
t156 = qJD(2) * t301 - t237 * t356;
t300 = qJ(5) * t156 + qJD(5) * t215 - t203;
t33 = t132 * pkin(4) - t133 * qJ(5) - qJD(5) * t304 + t175;
t253 = mrSges(3,3) * t361 - t372;
t206 = -pkin(7) * t369 + t239;
t202 = mrSges(4,1) * t272 - t391;
t201 = -mrSges(4,2) * t272 + t392;
t200 = -pkin(7) * t343 + t227;
t174 = -mrSges(4,2) * t334 + mrSges(4,3) * t198;
t173 = mrSges(4,1) * t334 - mrSges(4,3) * t197;
t168 = -pkin(4) * t303 - t307;
t155 = t287 * t452 - t330 * t356 - t331 * t358;
t152 = -mrSges(6,2) * t177 + mrSges(6,3) * t272;
t149 = -mrSges(5,2) * t272 - mrSges(5,3) * t177;
t146 = -qJD(3) * t207 + t363;
t145 = (-t290 * t359 - t293 * t357) * pkin(7) + t364;
t137 = t303 * t434 + t307;
t136 = pkin(4) * t214 - t322;
t135 = -mrSges(4,1) * t198 + mrSges(4,2) * t197;
t118 = t197 * Ifges(4,1) + t198 * Ifges(4,4) + Ifges(4,5) * t334;
t117 = t197 * Ifges(4,4) + t198 * Ifges(4,2) + Ifges(4,6) * t334;
t112 = mrSges(5,1) * t334 - mrSges(5,3) * t133;
t111 = -mrSges(5,2) * t334 - mrSges(5,3) * t132;
t110 = -mrSges(6,2) * t132 + mrSges(6,3) * t334;
t107 = mrSges(5,1) * t177 + mrSges(5,2) * t304;
t106 = mrSges(6,1) * t177 - mrSges(6,3) * t304;
t103 = -t214 * t434 + t322;
t78 = pkin(4) * t304 - t312;
t55 = -t304 * t434 + t312;
t50 = -pkin(4) * t155 - t300;
t46 = -qJD(6) * t154 - t155 * t291 - t156 * t288;
t45 = qJD(6) * t153 - t155 * t288 + t156 * t291;
t43 = -mrSges(7,1) * t99 + mrSges(7,2) * t468;
t39 = -pkin(4) * t360 + t320;
t37 = t99 * Ifges(7,2) + t261 * Ifges(7,6) + t387;
t34 = t155 * t434 + t300;
t26 = mrSges(7,2) * t334 + mrSges(7,3) * t30;
t25 = -mrSges(7,1) * t334 - mrSges(7,3) * t29;
t24 = -t155 * pkin(9) + t35;
t23 = -t156 * pkin(9) + t320 - t328;
t14 = t132 * pkin(5) + t33;
t4 = -qJD(6) * t32 + t23 * t291 - t24 * t288;
t3 = qJD(6) * t31 + t23 * t288 + t24 * t291;
t5 = [(-t20 * t214 + t215 * t321) * mrSges(5,3) - t320 * t150 + m(5) * (-t108 * t321 + t109 * t20 + t175 * t247 - t203 * t302 - t320 * t63 + t42 * t64) + (Ifges(6,5) * t215 + Ifges(6,3) * t214) * t426 - (t449 + t453) * t293 / 0.2e1 + (pkin(7) * t135 + t315 * t417 + t317 * t418 + t118 * t406 + t117 * t407 + (-t114 * t289 - t115 * t292) * mrSges(4,3) + (-t292 * t164 / 0.2e1 + t165 * t407 + t313 * t409 + t257 * t319 + t316 * t413 + t314 * t414 + (t184 * t289 - t185 * t292) * mrSges(4,3)) * qJD(3) + (-pkin(7) * t253 + t337 + (-Ifges(7,5) * t154 / 0.2e1 - Ifges(7,6) * t153 / 0.2e1 + (-0.3e1 / 0.2e1 * Ifges(3,4) + t386 / 0.2e1 - t385 / 0.2e1) * t290 - 0.2e1 * t433 + t351 * t215 + t350 * t214) * qJD(1) + t447) * qJD(2)) * t290 + (Ifges(7,4) * t154 + Ifges(7,2) * t153) * t441 + (Ifges(5,4) * t423 + Ifges(6,5) * t422 + t462 * t408 + t464 * t419 + t469) * t156 + (t326 * pkin(7) + t338 - 0.2e1 * t450 + t470) * t358 + m(4) * (t114 * t207 + t115 * t206 + t185 * t145 + t184 * t146) + (Ifges(7,1) * t154 + Ifges(7,4) * t153) * t442 + (Ifges(5,4) * t215 - Ifges(5,2) * t214) * t427 + (t1 * t153 - t154 * t2 - t45 * t8 + t46 * t9) * mrSges(7,3) + (-t15 * t214 + t18 * t215) * mrSges(6,2) + t154 * t443 + t153 * t444 + (Ifges(7,1) * t45 + Ifges(7,4) * t46) * t428 + (Ifges(7,4) * t45 + Ifges(7,2) * t46) * t435 + t214 * t439 + t214 * t440 + t247 * t59 + t33 * (mrSges(6,1) * t214 - mrSges(6,3) * t215) + t175 * (mrSges(5,1) * t214 + mrSges(5,2) * t215) + t206 * t173 + t207 * t174 + t145 * t201 + t146 * t202 + t203 * t107 + t42 * t149 + t39 * t151 + t35 * t152 - t14 * (-mrSges(7,1) * t153 + mrSges(7,2) * t154) + t136 * t58 + m(7) * (t1 * t32 - t103 * t14 + t2 * t31 + t3 * t9 + t34 * t49 + t4 * t8) + m(6) * (t102 * t15 + t104 * t18 + t136 * t33 + t35 * t57 + t39 * t56 + t50 * t72) + (t393 / 0.2e1 + Ifges(7,5) * t442 - Ifges(5,6) * t427 + Ifges(7,6) * t441 - Ifges(6,6) * t426 - t462 * t425 + (0.3e1 / 0.2e1 * Ifges(3,4) * t358 + (-Ifges(7,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + (m(4) * pkin(7) + t318) * pkin(7) - t327) * t360) * qJD(1) + t448 + t451) * t293 + (Ifges(5,2) * t423 - Ifges(6,3) * t422 + t408 * t461 + t419 * t463 + t445) * t155 + (-t214 * t463 + t215 * t464) * t425 + t215 * t490 + t45 * t491 + (Ifges(7,5) * t45 + Ifges(7,6) * t46) * t410 + t31 * t25 + t32 * t26 + t34 * t43 + t46 * t37 / 0.2e1 + t49 * (-mrSges(7,1) * t46 + mrSges(7,2) * t45) + t3 * t76 + t4 * t77 + t103 * t7 + t50 * t106 + t102 * t110 + t109 * t111 + t108 * t112 + t104 * t113; (t83 - t86) * (t223 / 0.2e1 - t204 / 0.2e1) - (-mrSges(5,1) * t366 - mrSges(5,2) * t365) * t302 + (t20 * t303 + t237 * t321 + t365 * t63 + t366 * t64) * mrSges(5,3) + (Ifges(6,5) * t237 - Ifges(6,3) * t303) * t426 + (Ifges(5,4) * t237 + Ifges(5,2) * t303) * t427 + (t15 * t303 + t18 * t237 - t365 * t56 + t366 * t57) * mrSges(6,2) + t33 * (-mrSges(6,1) * t303 - mrSges(6,3) * t237) + t175 * (-mrSges(5,1) * t303 + mrSges(5,2) * t237) - t303 * t439 - t303 * t440 + (t175 * t280 + t321 * t191 + t192 * t20 + t302 * t496 + t472 * t63 + t474 * t64) * m(5) + (t113 - t112) * t191 + (t110 + t111) * t192 + (-t138 / 0.2e1 + t94 / 0.2e1) * t37 + (-t139 / 0.2e1 + t93 / 0.2e1) * t38 + t460 * (t224 / 0.2e1 - t205 / 0.2e1) + t477 * t43 + (Ifges(5,4) * t205 + Ifges(6,5) * t224 - Ifges(5,2) * t204 + Ifges(6,3) * t223) * t422 - m(4) * (t184 * t199 + t185 * t200) + t479 * t76 + (t1 * t75 - t137 * t14 + t2 * t74 + t477 * t49 + t479 * t9 + t480 * t8) * m(7) + t480 * t77 + (-t223 * t461 + t224 * t462) * t408 + (-t204 * t461 + t205 * t462) * t409 + t474 * t149 + t475 * t152 + t476 * t151 + (t15 * t192 + t168 * t33 + t18 * t191 + t473 * t72 + t475 * t57 + t476 * t56) * m(6) + t289 * t118 / 0.2e1 + t280 * t59 + t472 * t150 + t473 * t106 + (Ifges(7,4) * t172 + Ifges(7,2) * t171) * t441 + (Ifges(7,1) * t172 + Ifges(7,4) * t171) * t442 + t172 * t443 + t171 * t444 + (Ifges(7,1) * t93 + Ifges(7,4) * t94) * t428 + (Ifges(7,1) * t139 + Ifges(7,4) * t138) * t429 + (Ifges(7,4) * t93 + Ifges(7,2) * t94) * t435 + (Ifges(7,4) * t139 + Ifges(7,2) * t138) * t436 + t316 * t418 + (mrSges(7,1) * t377 - mrSges(7,2) * t376) * t49 + (t1 * t171 - t172 * t2 + t376 * t8 - t377 * t9) * mrSges(7,3) + (-mrSges(6,1) * t366 + mrSges(6,3) * t365) * t72 - t234 * t107 - t200 * t201 - t199 * t202 - t14 * (-mrSges(7,1) * t171 + mrSges(7,2) * t172) + t168 * t58 - pkin(2) * t135 + t137 * t7 + (t237 * t464 + t303 * t463) * t425 + (-t223 * t463 + t224 * t464) * t419 + (-t204 * t463 + t205 * t464) * t420 + ((-t281 / 0.2e1 + t338 + t450 + ((-m(4) * pkin(2) - mrSges(3,1) - t319) * qJD(2) - t326) * pkin(7) - t470) * t293 + ((Ifges(7,5) * t172 + Ifges(7,6) * t171) * t481 + (t253 + t372) * pkin(7) + t337 + (t433 + t390 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t293) * qJD(1) - t447) * t290 + (t237 * t462 + t303 * t461 + t313) * t360 / 0.2e1) * qJD(1) + t237 * t490 + t117 * t406 + (Ifges(7,5) * t93 + Ifges(7,6) * t94) * t410 + (Ifges(7,5) * t139 + Ifges(7,6) * t138) * t411 + t314 * t417 + (t107 * t402 + t294) * qJD(3) + (m(4) * (-qJD(3) * t308 + t309) - t173 * t289 + t174 * t292 + (-t201 * t289 - t202 * t292) * qJD(3)) * pkin(8) + t309 * mrSges(4,3) + (Ifges(5,4) * t224 + Ifges(6,5) * t205 - Ifges(5,2) * t223 + Ifges(6,3) * t204) * t423 + t74 * t25 + t75 * t26; (t302 * t404 + t63 * t69 - t64 * t70 + (t20 * t287 - t321 * t375) * pkin(3)) * m(5) + t99 * t491 + (t37 - t493) * t429 - t499 + (-t201 + t392) * t184 + t454 * t69 + t455 * t152 + (t15 * t275 + t18 * t277 + t455 * t57 - t56 * t69 - t72 * t78) * m(6) + (-Ifges(4,2) * t244 + t165 + t235) * t414 + (-Ifges(5,4) * t422 - Ifges(6,5) * t423 - t462 * t409 - t464 * t420 + t469) * t177 + (t202 + t391) * t185 + t449 + t458 * t77 + (t1 * t209 + t2 * t208 + t458 * t8 + t459 * t9 - t49 * t55) * m(7) + t459 * t76 - t107 * t404 + t277 * t113 + t275 * t110 - t257 * (mrSges(4,1) * t244 + mrSges(4,2) * t243) + t208 * t25 + t209 * t26 - t70 * t149 - t448 + (-Ifges(5,2) * t422 + Ifges(6,3) * t423 - t409 * t461 - t420 * t463 + t445) * t304 + t500 * t436 + t112 * t345 + t111 * t403 + (Ifges(4,5) * t243 - Ifges(4,6) * t244) * t409 + t164 * t412 + (Ifges(4,1) * t243 - t380) * t413 - t55 * t43 - t78 * t106; t99 * t76 - t468 * t77 - (-t149 - t152) * t177 + t454 * t304 - t7 + t58 + t59 + (-t468 * t8 + t9 * t99 + t14) * m(7) + (t177 * t57 - t304 * t56 + t33) * m(6) + (t177 * t64 + t304 * t63 + t175) * m(5); t291 * t25 + t288 * t26 + (t106 - t43) * t304 + t311 * qJD(6) + (-t152 - t311) * t272 + t113 + (t1 * t288 - t304 * t49 + t2 * t291 + t261 * (-t288 * t8 + t291 * t9)) * m(7) + (-t272 * t57 + t304 * t72 + t18) * m(6); t493 * t429 + t37 * t428 - t8 * t76 + t9 * t77 + (t38 - t500) * t436 + t499;];
tauc  = t5(:);
