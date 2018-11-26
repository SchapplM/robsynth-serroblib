% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2018-11-23 17:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:58:21
% EndTime: 2018-11-23 17:58:38
% DurationCPUTime: 17.86s
% Computational Cost: add. (13515->775), mult. (31240->1033), div. (0->0), fcn. (20896->8), ass. (0->363)
t314 = cos(qJ(6));
t310 = sin(qJ(6));
t317 = cos(qJ(2));
t304 = t317 * qJD(1);
t292 = -t304 + qJD(3);
t318 = -pkin(3) - pkin(4);
t312 = sin(qJ(3));
t316 = cos(qJ(3));
t313 = sin(qJ(2));
t384 = qJD(1) * t313;
t366 = t316 * t384;
t255 = qJD(2) * t312 + t366;
t273 = -pkin(2) * t317 - pkin(8) * t313 - pkin(1);
t239 = t273 * qJD(1);
t301 = pkin(7) * t304;
t278 = qJD(2) * pkin(8) + t301;
t178 = t316 * t239 - t312 * t278;
t477 = qJD(4) - t178;
t544 = -pkin(9) * t255 + t477;
t108 = t292 * t318 + t544;
t179 = t312 * t239 + t316 * t278;
t368 = t312 * t384;
t376 = t316 * qJD(2);
t254 = t368 - t376;
t139 = pkin(9) * t254 + t179;
t282 = t292 * qJ(4);
t121 = t139 + t282;
t311 = sin(qJ(5));
t315 = cos(qJ(5));
t51 = t108 * t311 + t121 * t315;
t172 = t254 * t315 - t255 * t311;
t519 = pkin(10) * t172;
t46 = t51 + t519;
t405 = t310 * t46;
t285 = qJD(5) - t292;
t332 = t254 * t311 + t255 * t315;
t495 = pkin(10) * t332;
t50 = t315 * t108 - t121 * t311;
t45 = t50 - t495;
t41 = pkin(5) * t285 + t45;
t13 = t314 * t41 - t405;
t404 = t314 * t46;
t14 = t310 * t41 + t404;
t360 = t314 * t172 - t332 * t310;
t500 = t172 * t310 + t314 * t332;
t419 = Ifges(7,4) * t500;
t274 = qJD(6) + t285;
t43 = Ifges(7,2) * t360 + Ifges(7,6) * t274 + t419;
t94 = Ifges(7,4) * t360;
t44 = Ifges(7,1) * t500 + Ifges(7,5) * t274 + t94;
t441 = -t274 / 0.2e1;
t453 = t500 / 0.2e1;
t454 = -t500 / 0.2e1;
t456 = -t360 / 0.2e1;
t277 = -qJD(2) * pkin(2) + pkin(7) * t384;
t156 = t254 * pkin(3) - t255 * qJ(4) + t277;
t127 = -pkin(4) * t254 - t156;
t77 = -pkin(5) * t172 + t127;
t548 = (t13 * t360 + t14 * t500) * mrSges(7,3) + (Ifges(7,5) * t360 - Ifges(7,6) * t500) * t441 + (Ifges(7,1) * t360 - t419) * t454 - t77 * (mrSges(7,1) * t500 + mrSges(7,2) * t360) + t43 * t453 + (-Ifges(7,2) * t500 + t44 + t94) * t456;
t460 = pkin(8) - pkin(9);
t280 = t460 * t316;
t431 = pkin(7) * t312;
t369 = -pkin(3) - t431;
t329 = (-pkin(4) + t369) * t313;
t393 = t316 * t317;
t323 = -pkin(9) * t393 + t329;
t355 = pkin(2) * t313 - pkin(8) * t317;
t266 = t355 * qJD(1);
t396 = t266 * t316;
t547 = -qJD(1) * t323 + qJD(3) * t280 + t396;
t235 = t312 * t266;
t295 = qJ(4) * t384;
t381 = qJD(3) * t312;
t394 = t313 * t316;
t395 = t312 * t317;
t546 = t235 + t295 + (-pkin(7) * t394 + pkin(9) * t395) * qJD(1) + t460 * t381;
t330 = t311 * t312 + t315 * t316;
t476 = qJD(3) - qJD(5);
t185 = t476 * t330;
t327 = t317 * t330;
t220 = qJD(1) * t327;
t392 = t185 - t220;
t331 = t311 * t316 - t312 * t315;
t186 = t476 * t331;
t326 = t331 * t317;
t219 = qJD(1) * t326;
t537 = t186 - t219;
t375 = qJD(1) * qJD(2);
t361 = t313 * t375;
t427 = Ifges(6,3) + Ifges(7,3);
t439 = -t285 / 0.2e1;
t380 = qJD(3) * t313;
t364 = t312 * t380;
t365 = t317 * t376;
t374 = qJD(2) * qJD(3);
t207 = t316 * t374 + (-t364 + t365) * qJD(1);
t379 = qJD(3) * t316;
t382 = qJD(2) * t317;
t324 = t312 * t382 + t313 * t379;
t208 = qJD(1) * t324 + t312 * t374;
t71 = qJD(5) * t172 + t207 * t315 + t208 * t311;
t72 = -qJD(5) * t332 - t207 * t311 + t208 * t315;
t27 = qJD(6) * t360 + t310 * t72 + t314 * t71;
t28 = -qJD(6) * t500 - t310 * t71 + t314 * t72;
t377 = qJD(5) * t315;
t378 = qJD(5) * t311;
t269 = t355 * qJD(2);
t241 = qJD(1) * t269;
t359 = pkin(7) * t361;
t110 = t239 * t379 + t312 * t241 - t278 * t381 - t316 * t359;
t86 = qJ(4) * t361 + t292 * qJD(4) + t110;
t58 = pkin(9) * t208 + t86;
t354 = t239 * t381 - t241 * t316 + t278 * t379;
t61 = -pkin(9) * t207 + t329 * t375 + t354;
t11 = t108 * t377 - t121 * t378 + t311 * t61 + t315 * t58;
t10 = pkin(10) * t72 + t11;
t12 = -qJD(5) * t51 - t311 * t58 + t315 * t61;
t9 = -pkin(5) * t361 - pkin(10) * t71 + t12;
t3 = -qJD(6) * t14 - t10 * t310 + t314 * t9;
t474 = t3 * mrSges(7,1) + Ifges(7,5) * t27 + Ifges(7,6) * t28;
t473 = Ifges(6,5) * t71 + Ifges(6,6) * t72 + t474;
t543 = (t172 * t50 + t332 * t51) * mrSges(6,3) - t127 * (mrSges(6,1) * t332 + mrSges(6,2) * t172) + (Ifges(6,5) * t172 - Ifges(6,6) * t332) * t439 - t427 * t361 + t473 + t548;
t279 = t460 * t312;
t514 = t279 * t377 - t280 * t378 + t311 * t547 - t546 * t315;
t204 = t311 * t279 + t315 * t280;
t513 = -qJD(5) * t204 + t546 * t311 + t315 * t547;
t541 = qJD(2) / 0.2e1;
t530 = Ifges(4,1) + Ifges(5,1);
t517 = Ifges(5,4) + Ifges(4,5);
t540 = Ifges(4,6) - Ifges(5,6);
t539 = pkin(5) * t384 - pkin(10) * t392 + t513;
t538 = pkin(10) * t537 - t514;
t270 = -qJ(4) * t311 + t315 * t318;
t484 = qJD(5) * t270 - t311 * t139 + t315 * t544;
t271 = t315 * qJ(4) + t311 * t318;
t483 = -qJD(5) * t271 - t315 * t139 - t311 * t544;
t536 = -t312 * qJD(4) - t301 + (t304 * t316 - t379) * qJ(4);
t363 = Ifges(3,5) * t541;
t525 = -t495 + t484;
t524 = -t519 + t483;
t367 = t312 * t304;
t370 = t318 * t312;
t510 = qJD(3) * t370 - t318 * t367 - t536;
t166 = Ifges(6,4) * t172;
t522 = Ifges(6,2) * t332 - t166;
t518 = -qJD(2) / 0.2e1;
t203 = t315 * t279 - t280 * t311;
t151 = pkin(10) * t331 + t203;
t152 = -pkin(10) * t330 + t204;
t81 = t151 * t310 + t152 * t314;
t516 = -qJD(6) * t81 + t310 * t538 + t314 * t539;
t80 = t151 * t314 - t152 * t310;
t515 = qJD(6) * t80 + t310 * t539 - t314 * t538;
t512 = t517 * t361 + (-Ifges(4,4) + Ifges(5,5)) * t208 + t530 * t207;
t511 = pkin(5) * t537 + t510;
t509 = (-t367 + t381) * pkin(3) + t536;
t508 = t312 * t517 + t316 * t540;
t416 = Ifges(5,5) * t316;
t422 = Ifges(4,4) * t316;
t507 = t312 * t530 - t416 + t422;
t298 = Ifges(3,4) * t304;
t247 = Ifges(5,5) * t255;
t157 = t292 * Ifges(5,6) + t254 * Ifges(5,3) + t247;
t408 = t255 * Ifges(4,4);
t160 = -t254 * Ifges(4,2) + t292 * Ifges(4,6) + t408;
t333 = t178 * t316 + t179 * t312;
t148 = -pkin(3) * t292 + t477;
t150 = t282 + t179;
t334 = t148 * t316 - t150 * t312;
t340 = Ifges(5,3) * t312 + t416;
t344 = -Ifges(4,2) * t312 + t422;
t349 = mrSges(5,1) * t312 - mrSges(5,3) * t316;
t351 = mrSges(4,1) * t312 + mrSges(4,2) * t316;
t414 = Ifges(5,6) * t312;
t415 = Ifges(4,6) * t312;
t418 = Ifges(4,5) * t316;
t421 = Ifges(5,4) * t316;
t432 = t316 / 0.2e1;
t434 = t312 / 0.2e1;
t435 = -t312 / 0.2e1;
t442 = t255 / 0.2e1;
t444 = t254 / 0.2e1;
t445 = -t254 / 0.2e1;
t417 = Ifges(5,5) * t312;
t423 = Ifges(4,4) * t312;
t478 = t316 * t530 + t417 - t423;
t248 = Ifges(4,4) * t254;
t409 = t254 * Ifges(5,5);
t481 = t255 * t530 + t292 * t517 - t248 + t409;
t497 = t292 / 0.2e1;
t319 = t334 * mrSges(5,2) - t333 * mrSges(4,3) + t156 * t349 + t157 * t434 + t160 * t435 + t277 * t351 + t340 * t444 + t344 * t445 + t478 * t442 + (-t415 + t418 + t414 + t421) * t497 + t481 * t432;
t506 = t319 + Ifges(3,1) * t384 / 0.2e1 + t298 / 0.2e1 + t363;
t420 = Ifges(6,4) * t332;
t503 = Ifges(6,1) * t172 - t420;
t468 = t27 / 0.2e1;
t467 = t28 / 0.2e1;
t462 = t71 / 0.2e1;
t461 = t72 / 0.2e1;
t84 = Ifges(6,1) * t332 + t285 * Ifges(6,5) + t166;
t498 = t84 / 0.2e1;
t450 = -t332 / 0.2e1;
t496 = -t361 / 0.2e1;
t362 = Ifges(3,6) * t518;
t265 = -pkin(5) + t270;
t181 = t265 * t314 - t271 * t310;
t494 = qJD(6) * t181 + t310 * t524 + t314 * t525;
t182 = t265 * t310 + t271 * t314;
t493 = -qJD(6) * t182 - t310 * t525 + t314 * t524;
t229 = t331 * t313;
t293 = pkin(7) * t395;
t308 = t317 * pkin(3);
t430 = pkin(9) * t313;
t164 = pkin(4) * t317 + t293 + t308 + (-t273 - t430) * t316;
t294 = pkin(7) * t393;
t218 = t312 * t273 + t294;
t205 = -qJ(4) * t317 + t218;
t177 = t312 * t430 + t205;
t93 = t311 * t164 + t315 * t177;
t256 = -t310 * t311 + t314 * t315;
t480 = t274 * t256;
t258 = t310 * t315 + t311 * t314;
t479 = t274 * t258;
t371 = Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1;
t372 = Ifges(4,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t373 = Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t424 = Ifges(3,4) * t313;
t471 = -t372 * t254 + t373 * t255 + t371 * t292 + t14 * mrSges(7,2) + t150 * mrSges(5,3) + t178 * mrSges(4,1) + t51 * mrSges(6,2) + Ifges(4,6) * t445 + Ifges(5,6) * t444 + t362 - (t317 * Ifges(3,2) + t424) * qJD(1) / 0.2e1 - t274 * Ifges(7,3) - t500 * Ifges(7,5) - t360 * Ifges(7,6) - t285 * Ifges(6,3) - t332 * Ifges(6,5) - t172 * Ifges(6,6) - t13 * mrSges(7,1) - t148 * mrSges(5,1) - t179 * mrSges(4,2) - t50 * mrSges(6,1) + (Ifges(4,3) + Ifges(5,2)) * t497 + t517 * t442;
t470 = Ifges(7,4) * t468 + Ifges(7,2) * t467 + Ifges(7,6) * t496;
t469 = Ifges(7,1) * t468 + Ifges(7,4) * t467 + Ifges(7,5) * t496;
t466 = Ifges(6,4) * t462 + Ifges(6,2) * t461 + Ifges(6,6) * t496;
t465 = Ifges(6,1) * t462 + Ifges(6,4) * t461 + Ifges(6,5) * t496;
t459 = pkin(1) * mrSges(3,1);
t458 = pkin(1) * mrSges(3,2);
t2 = qJD(6) * t13 + t10 * t314 + t310 * t9;
t457 = t2 * mrSges(7,2);
t455 = t360 / 0.2e1;
t452 = -t172 / 0.2e1;
t451 = t172 / 0.2e1;
t449 = t332 / 0.2e1;
t448 = t207 / 0.2e1;
t447 = -t208 / 0.2e1;
t446 = t208 / 0.2e1;
t443 = -t255 / 0.2e1;
t440 = t274 / 0.2e1;
t438 = t285 / 0.2e1;
t437 = -t292 / 0.2e1;
t433 = -t316 / 0.2e1;
t426 = mrSges(4,3) * t254;
t425 = mrSges(4,3) * t255;
t131 = -t219 * t314 - t220 * t310;
t176 = -t310 * t330 - t314 * t331;
t60 = -qJD(6) * t176 - t185 * t310 - t186 * t314;
t403 = t131 - t60;
t132 = -t219 * t310 + t220 * t314;
t175 = t310 * t331 - t314 * t330;
t59 = qJD(6) * t175 + t185 * t314 - t186 * t310;
t402 = t132 - t59;
t398 = qJD(2) * mrSges(3,2);
t305 = t312 * qJ(4);
t211 = -mrSges(4,2) * t292 - t426;
t214 = -mrSges(5,2) * t254 + mrSges(5,3) * t292;
t390 = -t211 - t214;
t212 = mrSges(4,1) * t292 - t425;
t213 = -mrSges(5,1) * t292 + mrSges(5,2) * t255;
t389 = -t212 + t213;
t388 = t312 * t269 + t273 * t379;
t189 = t255 * pkin(3) + t254 * qJ(4);
t386 = qJ(4) * t365 + qJD(4) * t394;
t383 = qJD(2) * t313;
t272 = -t316 * pkin(3) - pkin(2) - t305;
t92 = t315 * t164 - t177 * t311;
t217 = t273 * t316 - t293;
t246 = t316 * pkin(4) - t272;
t358 = -pkin(7) + t370;
t357 = m(4) * t277 - qJD(2) * mrSges(3,1) + mrSges(4,1) * t254 + mrSges(4,2) * t255 + mrSges(3,3) * t384;
t143 = -pkin(4) * t255 - t189;
t356 = t313 * t369;
t353 = qJD(3) * t294 - t269 * t316 + t273 * t381;
t352 = mrSges(4,1) * t316 - mrSges(4,2) * t312;
t350 = mrSges(5,1) * t316 + mrSges(5,3) * t312;
t343 = Ifges(4,2) * t316 + t423;
t339 = -Ifges(5,3) * t316 + t417;
t230 = t330 * t313;
t65 = pkin(5) * t317 - pkin(10) * t230 + t92;
t68 = -pkin(10) * t229 + t93;
t33 = -t310 * t68 + t314 * t65;
t34 = t310 * t65 + t314 * t68;
t328 = qJD(2) * t356;
t105 = qJD(1) * t328 + t354;
t337 = t105 * t312 + t316 * t86;
t111 = t312 * t359 - t354;
t336 = t110 * t316 - t111 * t312;
t140 = -mrSges(6,2) * t285 + mrSges(6,3) * t172;
t141 = mrSges(6,1) * t285 - mrSges(6,3) * t332;
t335 = t140 * t315 - t141 * t311;
t145 = -t229 * t314 - t230 * t310;
t146 = -t229 * t310 + t230 * t314;
t210 = -pkin(7) * t366 + t235;
t169 = -mrSges(5,1) * t361 + t207 * mrSges(5,2);
t95 = pkin(9) * t364 + qJD(2) * t323 + t353;
t297 = qJ(4) * t383;
t97 = t297 + (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t394 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t312) * t317 + t388;
t31 = t164 * t377 - t177 * t378 + t311 * t95 + t315 * t97;
t325 = t12 * mrSges(6,1) - t11 * mrSges(6,2) - t457;
t291 = qJ(4) * t394;
t202 = t313 * t358 + t291;
t96 = t208 * pkin(3) - t207 * qJ(4) + qJD(2) * t301 - t255 * qJD(4);
t62 = -pkin(4) * t208 - t96;
t32 = -qJD(5) * t93 - t311 * t97 + t315 * t95;
t134 = (-t313 * t376 - t317 * t381) * pkin(7) + t388;
t322 = -t111 * mrSges(4,1) + t105 * mrSges(5,1) + t110 * mrSges(4,2) - t86 * mrSges(5,3) + t325;
t109 = (t316 * t318 - t305) * t380 + t358 * t382 + t386;
t290 = Ifges(5,2) * t361;
t289 = Ifges(4,3) * t361;
t276 = mrSges(3,3) * t304 - t398;
t223 = -t291 + (pkin(3) * t312 + pkin(7)) * t313;
t209 = pkin(7) * t368 + t396;
t206 = -t217 + t308;
t200 = Ifges(5,4) * t207;
t199 = Ifges(4,5) * t207;
t198 = Ifges(4,6) * t208;
t197 = Ifges(5,6) * t208;
t192 = pkin(5) * t330 + t246;
t190 = mrSges(5,1) * t254 - mrSges(5,3) * t255;
t188 = qJD(1) * t356 - t396;
t187 = t210 + t295;
t170 = -mrSges(4,2) * t361 - mrSges(4,3) * t208;
t168 = mrSges(4,1) * t361 - mrSges(4,3) * t207;
t167 = -mrSges(5,2) * t208 + mrSges(5,3) * t361;
t142 = pkin(5) * t229 + t202;
t135 = t383 * t431 - t353;
t133 = pkin(3) * t324 + pkin(7) * t382 + qJ(4) * t364 - t386;
t126 = t328 + t353;
t125 = mrSges(4,1) * t208 + mrSges(4,2) * t207;
t124 = mrSges(5,1) * t208 - mrSges(5,3) * t207;
t120 = -qJD(4) * t317 + t134 + t297;
t115 = t207 * Ifges(4,4) - t208 * Ifges(4,2) + Ifges(4,6) * t361;
t114 = t207 * Ifges(5,5) + Ifges(5,6) * t361 + t208 * Ifges(5,3);
t113 = qJD(2) * t327 + t229 * t476;
t112 = -qJD(2) * t326 + t185 * t313;
t104 = -mrSges(6,1) * t172 + mrSges(6,2) * t332;
t85 = -pkin(5) * t332 + t143;
t83 = t172 * Ifges(6,2) + t285 * Ifges(6,6) + t420;
t79 = mrSges(7,1) * t274 - mrSges(7,3) * t500;
t78 = -mrSges(7,2) * t274 + mrSges(7,3) * t360;
t64 = mrSges(6,2) * t361 + mrSges(6,3) * t72;
t63 = -mrSges(6,1) * t361 - mrSges(6,3) * t71;
t52 = -pkin(5) * t112 + t109;
t47 = -mrSges(7,1) * t360 + mrSges(7,2) * t500;
t40 = -qJD(6) * t146 + t112 * t314 - t113 * t310;
t39 = qJD(6) * t145 + t112 * t310 + t113 * t314;
t38 = -mrSges(6,1) * t72 + mrSges(6,2) * t71;
t37 = -pkin(5) * t72 + t62;
t24 = mrSges(7,2) * t361 + mrSges(7,3) * t28;
t23 = -mrSges(7,1) * t361 - mrSges(7,3) * t27;
t18 = pkin(10) * t112 + t31;
t17 = -pkin(5) * t383 - pkin(10) * t113 + t32;
t16 = t314 * t45 - t405;
t15 = -t310 * t45 - t404;
t8 = -mrSges(7,1) * t28 + mrSges(7,2) * t27;
t5 = -qJD(6) * t34 + t17 * t314 - t18 * t310;
t4 = qJD(6) * t33 + t17 * t310 + t18 * t314;
t1 = [(-t13 * t39 + t14 * t40 + t145 * t2 - t146 * t3) * mrSges(7,3) + m(4) * (t110 * t218 + t111 * t217 + t134 * t179 + t135 * t178) + t113 * t498 + m(7) * (t13 * t5 + t14 * t4 + t142 * t37 + t2 * t34 + t3 * t33 + t52 * t77) + m(6) * (t109 * t127 + t11 * t93 + t12 * t92 + t202 * t62 + t31 * t51 + t32 * t50) + m(5) * (t105 * t206 + t120 * t150 + t126 * t148 + t133 * t156 + t205 * t86 + t223 * t96) + t223 * t124 + t217 * t168 + t218 * t170 + t134 * t211 + t135 * t212 + t126 * t213 + t120 * t214 + t205 * t167 + t206 * t169 + t202 * t38 + t133 * t190 + t37 * (-mrSges(7,1) * t145 + mrSges(7,2) * t146) + t31 * t140 + t32 * t141 + t142 * t8 + t127 * (-mrSges(6,1) * t112 + mrSges(6,2) * t113) + t109 * t104 + t112 * t83 / 0.2e1 + t92 * t63 + t93 * t64 + t77 * (-mrSges(7,1) * t40 + mrSges(7,2) * t39) + t4 * t78 + t5 * t79 + t52 * t47 + t40 * t43 / 0.2e1 + t39 * t44 / 0.2e1 + t33 * t23 + t34 * t24 + t230 * t465 - t229 * t466 + (Ifges(7,4) * t146 + Ifges(7,2) * t145) * t467 + (Ifges(7,1) * t146 + Ifges(7,4) * t145) * t468 + t146 * t469 + t145 * t470 + (Ifges(6,1) * t113 + Ifges(6,4) * t112) * t449 + (Ifges(6,4) * t113 + Ifges(6,2) * t112) * t451 + (Ifges(7,1) * t39 + Ifges(7,4) * t40) * t453 + (Ifges(7,4) * t39 + Ifges(7,2) * t40) * t455 + (Ifges(6,5) * t113 + Ifges(6,6) * t112) * t438 + (Ifges(7,5) * t39 + Ifges(7,6) * t40) * t440 + (-t11 * t229 + t112 * t51 - t113 * t50 - t12 * t230) * mrSges(6,3) + t62 * (mrSges(6,1) * t229 + mrSges(6,2) * t230) + (Ifges(6,4) * t230 - Ifges(6,2) * t229) * t461 + (Ifges(6,1) * t230 - Ifges(6,4) * t229) * t462 + (-t373 * t207 + t372 * t208 - t289 / 0.2e1 - t290 / 0.2e1 - t199 / 0.2e1 - t200 / 0.2e1 - t197 / 0.2e1 + t198 / 0.2e1 + t322 + ((0.3e1 / 0.2e1 * Ifges(3,4) * t317 - 0.2e1 * t458) * qJD(1) + t363 + t357 * pkin(7) + t506) * qJD(2) + t473) * t317 + (t114 * t434 + pkin(7) * t125 + t96 * t349 + t340 * t446 + t344 * t447 + (-t110 * t312 - t111 * t316) * mrSges(4,3) + (t105 * t316 - t312 * t86) * mrSges(5,2) + (t156 * t350 + t277 * t352 + t343 * t444 + t339 * t445 + t160 * t433 + (t178 * t312 - t179 * t316) * mrSges(4,3) + (-t148 * t312 - t150 * t316) * mrSges(5,2) + t507 * t443 + t508 * t437) * qJD(3) + (-pkin(7) * t276 + t362 + (-Ifges(7,5) * t146 / 0.2e1 - Ifges(7,6) * t145 / 0.2e1 - Ifges(6,5) * t230 / 0.2e1 + Ifges(6,6) * t229 / 0.2e1 + (-0.3e1 / 0.2e1 * Ifges(3,4) + t418 / 0.2e1 - t415 / 0.2e1 + t421 / 0.2e1 + t414 / 0.2e1) * t313 - 0.2e1 * t459 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + (m(4) * pkin(7) + t351) * pkin(7) - t371 - t427) * t317) * qJD(1) + t471) * qJD(2) + t478 * t448 + (qJD(3) * t481 + t115) * t435 + (qJD(3) * t157 + t512) * t432) * t313; (-t132 / 0.2e1 + t59 / 0.2e1) * t44 + (mrSges(6,1) * t537 + mrSges(6,2) * t392) * t127 + (-t11 * t330 + t12 * t331 - t392 * t50 - t51 * t537) * mrSges(6,3) + ((-m(4) * t333 + m(5) * t334 + t312 * t390 + t316 * t389) * qJD(3) + m(5) * t337 + m(4) * t336 + (t167 + t170) * t316 + (-t168 + t169) * t312) * pkin(8) - t330 * t466 + (-t220 / 0.2e1 + t185 / 0.2e1) * t84 + t319 * qJD(3) + (-t131 / 0.2e1 + t60 / 0.2e1) * t43 - m(4) * (t178 * t209 + t179 * t210) + (((t276 + t398) * pkin(7) + (t459 + t424 / 0.2e1) * qJD(1) + t362 + (-Ifges(6,5) * t331 + Ifges(7,5) * t176 - Ifges(6,6) * t330 + Ifges(7,6) * t175) * t518 + t508 * t541 - t471) * t313 + ((t458 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t313) * qJD(1) - t298 / 0.2e1 + t363 + ((-m(4) * pkin(2) - mrSges(3,1) - t352) * qJD(2) - t357) * pkin(7) - t506) * t317) * qJD(1) + t272 * t124 + t246 * t38 - t210 * t211 - t209 * t212 - t188 * t213 - t187 * t214 + t204 * t64 + t203 * t63 + t192 * t8 + t37 * (-mrSges(7,1) * t175 + mrSges(7,2) * t176) + (mrSges(7,1) * t403 - mrSges(7,2) * t402) * t77 + (t13 * t402 - t14 * t403 + t175 * t2 - t176 * t3) * mrSges(7,3) - pkin(2) * t125 + t80 * t23 + t81 * t24 - t96 * t350 + t337 * mrSges(5,2) + t336 * mrSges(4,3) + (Ifges(7,4) * t176 + Ifges(7,2) * t175) * t467 + (Ifges(7,1) * t176 + Ifges(7,4) * t175) * t468 + t176 * t469 + t175 * t470 + t339 * t446 + t343 * t447 + (Ifges(7,1) * t59 + Ifges(7,4) * t60) * t453 + (Ifges(7,1) * t132 + Ifges(7,4) * t131) * t454 + (Ifges(7,4) * t59 + Ifges(7,2) * t60) * t455 + (Ifges(7,4) * t132 + Ifges(7,2) * t131) * t456 + (Ifges(7,5) * t59 + Ifges(7,6) * t60) * t440 + (Ifges(7,5) * t132 + Ifges(7,6) * t131) * t441 + t115 * t432 + t114 * t433 + t62 * (mrSges(6,1) * t330 - mrSges(6,2) * t331) + (-Ifges(6,4) * t331 - Ifges(6,2) * t330) * t461 + (-Ifges(6,1) * t331 - Ifges(6,4) * t330) * t462 - t331 * t465 + (Ifges(6,1) * t185 - Ifges(6,4) * t186) * t449 + (Ifges(6,4) * t185 - Ifges(6,2) * t186) * t451 + (Ifges(6,5) * t185 - Ifges(6,6) * t186) * t438 + (Ifges(6,1) * t220 - Ifges(6,4) * t219) * t450 + (Ifges(6,4) * t220 - Ifges(6,2) * t219) * t452 + (Ifges(6,5) * t220 - Ifges(6,6) * t219) * t439 + (t219 / 0.2e1 - t186 / 0.2e1) * t83 + t507 * t448 + t509 * t190 + (-t148 * t188 - t150 * t187 + t509 * t156 + t272 * t96) * m(5) + t510 * t104 + t511 * t47 + t512 * t434 + t513 * t141 + t514 * t140 + (t11 * t204 + t12 * t203 + t127 * t510 + t246 * t62 + t50 * t513 + t51 * t514) * m(6) + t515 * t78 + t516 * t79 + (t13 * t516 + t14 * t515 + t192 * t37 + t2 * t81 + t3 * t80 + t511 * t77) * m(7); t172 * t498 + t493 * t79 + (t13 * t493 + t14 * t494 + t181 * t3 + t182 * t2 - t77 * t85) * m(7) + t494 * t78 + (-Ifges(4,2) * t255 - t248 + t481) * t444 + (-t254 * t530 + t157 + t247 - t408) * t443 + (t148 * t254 + t150 * t255) * mrSges(5,2) + (-pkin(3) * t105 + qJ(4) * t86 - t148 * t179 + t150 * t477 - t156 * t189) * m(5) + t483 * t141 + t484 * t140 + (t11 * t271 + t12 * t270 - t127 * t143 + t483 * t50 + t484 * t51) * m(6) + t289 + t290 + t522 * t452 + (-t517 * t254 - t255 * t540) * t437 - t277 * (mrSges(4,1) * t255 - mrSges(4,2) * t254) + t270 * t63 + t271 * t64 - t156 * (mrSges(5,1) * t255 + mrSges(5,3) * t254) + qJD(4) * t214 - t189 * t190 + (-t389 + t425) * t179 + (t390 - t426) * t178 + t181 * t23 + t182 * t24 + qJ(4) * t167 - pkin(3) * t169 - t143 * t104 - t85 * t47 + t160 * t442 + (Ifges(5,3) * t255 - t409) * t445 + t199 + t200 + t197 - t198 - t322 + (-t503 + t83) * t450 - t543; t256 * t23 + t258 * t24 + t311 * t64 + t315 * t63 - t479 * t79 + t480 * t78 + t335 * qJD(5) + (-t214 - t335) * t292 + (-t104 + t190 - t47) * t255 + t169 + (-t13 * t479 + t14 * t480 + t2 * t258 - t255 * t77 + t256 * t3) * m(7) + (t11 * t311 + t12 * t315 - t127 * t255 + t285 * (-t311 * t50 + t315 * t51)) * m(6) + (-t150 * t292 + t156 * t255 + t105) * m(5); t325 + (-t332 * t47 + t314 * t23 + t310 * t24 + (-t310 * t79 + t314 * t78) * qJD(6) + (-t332 * t77 + t2 * t310 + t3 * t314 + (-t13 * t310 + t14 * t314) * qJD(6)) * m(7)) * pkin(5) + (t84 - t522) * t452 - t50 * t140 + t51 * t141 - t16 * t78 - t15 * t79 + t83 * t449 - m(7) * (t13 * t15 + t14 * t16) + t503 * t450 + t543; -Ifges(7,3) * t361 - t13 * t78 + t14 * t79 - t457 + t474 + t548;];
tauc  = t1(:);
