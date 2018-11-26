% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 17:22
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:21:46
% EndTime: 2018-11-23 17:22:06
% DurationCPUTime: 20.97s
% Computational Cost: add. (22517->753), mult. (56452->1038), div. (0->0), fcn. (42777->10), ass. (0->352)
t293 = sin(qJ(2));
t297 = cos(qJ(2));
t367 = sin(pkin(11));
t322 = qJD(1) * t367;
t368 = cos(pkin(11));
t323 = qJD(1) * t368;
t257 = -t293 * t323 - t297 * t322;
t292 = sin(qJ(4));
t296 = cos(qJ(4));
t220 = qJD(2) * t296 + t257 * t292;
t221 = qJD(2) * t292 - t257 * t296;
t291 = sin(qJ(5));
t295 = cos(qJ(5));
t158 = t220 * t291 + t221 * t295;
t290 = sin(qJ(6));
t294 = cos(qJ(6));
t319 = t295 * t220 - t221 * t291;
t108 = t158 * t294 + t290 * t319;
t477 = pkin(10) * t319;
t256 = -t293 * t322 + t297 * t323;
t288 = -pkin(2) * t297 - pkin(1);
t351 = qJD(1) * t288;
t275 = qJD(3) + t351;
t177 = -t256 * pkin(3) + t257 * pkin(8) + t275;
t393 = -qJ(3) - pkin(7);
t279 = t393 * t293;
t273 = qJD(1) * t279;
t265 = qJD(2) * pkin(2) + t273;
t280 = t393 * t297;
t274 = qJD(1) * t280;
t325 = t368 * t274;
t208 = t367 * t265 - t325;
t198 = qJD(2) * pkin(8) + t208;
t131 = t177 * t292 + t198 * t296;
t111 = pkin(9) * t220 + t131;
t102 = t295 * t111;
t130 = t296 * t177 - t198 * t292;
t110 = -pkin(9) * t221 + t130;
t248 = qJD(4) - t256;
t91 = pkin(4) * t248 + t110;
t56 = t291 * t91 + t102;
t43 = t56 + t477;
t372 = t290 * t43;
t239 = qJD(5) + t248;
t489 = pkin(10) * t158;
t100 = t291 * t111;
t55 = t295 * t91 - t100;
t42 = t55 - t489;
t38 = pkin(5) * t239 + t42;
t16 = t294 * t38 - t372;
t370 = t294 * t43;
t17 = t290 * t38 + t370;
t269 = t293 * t368 + t297 * t367;
t258 = t269 * qJD(2);
t240 = qJD(1) * t258;
t482 = -t158 * t290 + t294 * t319;
t268 = t293 * t367 - t297 * t368;
t259 = t268 * qJD(2);
t241 = qJD(1) * t259;
t171 = qJD(4) * t220 - t241 * t296;
t172 = -qJD(4) * t221 + t241 * t292;
t80 = qJD(5) * t319 + t171 * t295 + t172 * t291;
t81 = -qJD(5) * t158 - t171 * t291 + t172 * t295;
t32 = qJD(6) * t482 + t290 * t81 + t294 * t80;
t33 = -qJD(6) * t108 - t290 * t80 + t294 * t81;
t340 = Ifges(7,5) * t32 + Ifges(7,6) * t33 + Ifges(7,3) * t240;
t386 = Ifges(7,4) * t108;
t235 = qJD(6) + t239;
t411 = -t235 / 0.2e1;
t423 = -t108 / 0.2e1;
t327 = qJD(2) * t393;
t252 = qJD(3) * t297 + t293 * t327;
t226 = t252 * qJD(1);
t253 = -t293 * qJD(3) + t297 * t327;
t299 = qJD(1) * t253;
t168 = t226 * t368 + t299 * t367;
t341 = qJD(1) * qJD(2);
t330 = t293 * t341;
t318 = pkin(2) * t330;
t174 = pkin(3) * t240 + pkin(8) * t241 + t318;
t73 = -qJD(4) * t131 - t168 * t292 + t296 * t174;
t50 = pkin(4) * t240 - pkin(9) * t171 + t73;
t346 = qJD(4) * t296;
t347 = qJD(4) * t292;
t72 = t296 * t168 + t292 * t174 + t177 * t346 - t198 * t347;
t59 = pkin(9) * t172 + t72;
t13 = -qJD(5) * t56 - t291 * t59 + t295 * t50;
t6 = pkin(5) * t240 - pkin(10) * t80 + t13;
t344 = qJD(5) * t295;
t345 = qJD(5) * t291;
t12 = -t111 * t345 + t291 * t50 + t295 * t59 + t91 * t344;
t7 = pkin(10) * t81 + t12;
t2 = qJD(6) * t16 + t290 * t6 + t294 * t7;
t3 = -qJD(6) * t17 - t290 * t7 + t294 * t6;
t476 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t260 = t367 * t274;
t207 = t265 * t368 + t260;
t197 = -qJD(2) * pkin(3) - t207;
t152 = -t220 * pkin(4) + t197;
t96 = -pkin(5) * t319 + t152;
t508 = t340 + t476 + (Ifges(7,1) * t482 - t386) * t423 + (Ifges(7,5) * t482 - Ifges(7,6) * t108) * t411 + (t108 * t17 + t16 * t482) * mrSges(7,3) - t96 * (mrSges(7,1) * t108 + mrSges(7,2) * t482);
t350 = qJD(1) * t293;
t338 = pkin(2) * t350;
t188 = -pkin(3) * t257 - pkin(8) * t256 + t338;
t211 = t273 * t368 + t260;
t138 = t296 * t188 - t211 * t292;
t334 = t367 * pkin(2);
t285 = t334 + pkin(8);
t394 = pkin(9) + t285;
t326 = qJD(4) * t394;
t359 = t256 * t296;
t507 = pkin(4) * t257 + pkin(9) * t359 - t296 * t326 - t138;
t139 = t292 * t188 + t296 * t211;
t360 = t256 * t292;
t506 = -pkin(9) * t360 + t292 * t326 + t139;
t339 = Ifges(6,5) * t80 + Ifges(6,6) * t81 + Ifges(6,3) * t240;
t387 = Ifges(6,4) * t158;
t409 = -t239 / 0.2e1;
t419 = -t158 / 0.2e1;
t425 = -t482 / 0.2e1;
t97 = Ifges(7,4) * t482;
t54 = Ifges(7,1) * t108 + Ifges(7,5) * t235 + t97;
t434 = -t54 / 0.2e1;
t53 = Ifges(7,2) * t482 + Ifges(7,6) * t235 + t386;
t436 = -t53 / 0.2e1;
t464 = t13 * mrSges(6,1) - t12 * mrSges(6,2);
t484 = -Ifges(7,2) * t108 + t97;
t505 = -t108 * t436 + t339 + t464 + (Ifges(6,5) * t319 - Ifges(6,6) * t158) * t409 + (t158 * t56 + t319 * t55) * mrSges(6,3) - t152 * (mrSges(6,1) * t158 + mrSges(6,2) * t319) + t482 * t434 + (Ifges(6,1) * t319 - t387) * t419 + t484 * t425 + t508;
t305 = t291 * t292 - t295 * t296;
t184 = t305 * t256;
t450 = qJD(4) + qJD(5);
t214 = t450 * t305;
t504 = t214 - t184;
t272 = t291 * t296 + t292 * t295;
t183 = t272 * t256;
t215 = t450 * t272;
t492 = t215 - t183;
t266 = t394 * t292;
t267 = t394 * t296;
t202 = -t291 * t266 + t295 * t267;
t463 = -qJD(5) * t202 + t506 * t291 + t295 * t507;
t462 = -t266 * t344 - t267 * t345 + t291 * t507 - t506 * t295;
t501 = t347 - t360;
t402 = -t292 / 0.2e1;
t405 = -t256 / 0.2e1;
t496 = Ifges(4,2) * t405;
t494 = -pkin(10) * t492 + t462;
t493 = pkin(5) * t257 + pkin(10) * t504 + t463;
t287 = pkin(4) * t295 + pkin(5);
t342 = qJD(6) * t294;
t343 = qJD(6) * t290;
t356 = t290 * t291;
t61 = -t110 * t291 - t102;
t44 = t61 - t477;
t62 = t295 * t110 - t100;
t45 = t62 - t489;
t488 = -t290 * t44 - t294 * t45 + t287 * t342 + (-t291 * t343 + (t294 * t295 - t356) * qJD(5)) * pkin(4);
t355 = t291 * t294;
t487 = t290 * t45 - t294 * t44 - t287 * t343 + (-t291 * t342 + (-t290 * t295 - t355) * qJD(5)) * pkin(4);
t209 = t273 * t367 - t325;
t457 = pkin(4) * t501 - t209;
t486 = -t259 * t292 + t269 * t346;
t380 = t221 * Ifges(5,4);
t146 = t220 * Ifges(5,2) + t248 * Ifges(5,6) + t380;
t316 = mrSges(5,1) * t292 + mrSges(5,2) * t296;
t485 = t146 * t402 + t197 * t316;
t154 = Ifges(6,4) * t319;
t483 = -Ifges(6,2) * t158 + t154;
t479 = -Ifges(4,6) / 0.2e1;
t440 = t32 / 0.2e1;
t439 = t33 / 0.2e1;
t432 = t80 / 0.2e1;
t431 = t81 / 0.2e1;
t407 = t240 / 0.2e1;
t348 = qJD(2) * t293;
t478 = pkin(2) * t348;
t201 = -t295 * t266 - t267 * t291;
t169 = -pkin(10) * t272 + t201;
t170 = -pkin(10) * t305 + t202;
t115 = t169 * t294 - t170 * t290;
t475 = qJD(6) * t115 + t290 * t493 + t294 * t494;
t116 = t169 * t290 + t170 * t294;
t474 = -qJD(6) * t116 - t290 * t494 + t294 * t493;
t466 = t275 * mrSges(4,2);
t461 = Ifges(4,5) * qJD(2);
t194 = t305 * t269;
t212 = t272 * t294 - t290 * t305;
t123 = -qJD(6) * t212 + t214 * t290 - t215 * t294;
t128 = -t183 * t294 + t184 * t290;
t459 = -t128 + t123;
t210 = -t272 * t290 - t294 * t305;
t122 = qJD(6) * t210 - t214 * t294 - t215 * t290;
t129 = -t183 * t290 - t184 * t294;
t458 = -t129 + t122;
t206 = t268 * pkin(3) - t269 * pkin(8) + t288;
t218 = t279 * t367 - t280 * t368;
t150 = t296 * t206 - t218 * t292;
t398 = pkin(9) * t296;
t124 = pkin(4) * t268 - t269 * t398 + t150;
t213 = t296 * t218;
t151 = t292 * t206 + t213;
t357 = t269 * t292;
t133 = -pkin(9) * t357 + t151;
t75 = t291 * t124 + t295 * t133;
t456 = Ifges(5,5) * t171 + Ifges(5,6) * t172;
t455 = pkin(5) * t492 + t457;
t391 = mrSges(4,3) * t257;
t352 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t220 + mrSges(5,2) * t221 - t391;
t452 = -t292 * t73 + t296 * t72;
t451 = t73 * mrSges(5,1) - t72 * mrSges(5,2);
t349 = qJD(1) * t297;
t365 = Ifges(3,6) * qJD(2);
t390 = Ifges(3,4) * t293;
t449 = t365 / 0.2e1 + (t297 * Ifges(3,2) + t390) * qJD(1) / 0.2e1 + pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t349);
t448 = qJD(2) * t479;
t373 = t257 * Ifges(4,4);
t447 = t373 / 0.2e1 + t496;
t311 = Ifges(5,5) * t296 - Ifges(5,6) * t292;
t388 = Ifges(5,4) * t296;
t313 = -Ifges(5,2) * t292 + t388;
t389 = Ifges(5,4) * t292;
t315 = Ifges(5,1) * t296 - t389;
t219 = Ifges(5,4) * t220;
t147 = t221 * Ifges(5,1) + t248 * Ifges(5,5) + t219;
t353 = t296 * t147;
t412 = t221 / 0.2e1;
t446 = t248 * t311 / 0.2e1 + t315 * t412 + t220 * t313 / 0.2e1 + t353 / 0.2e1 + t485;
t445 = t275 * mrSges(4,1) + t130 * mrSges(5,1) + t55 * mrSges(6,1) + t16 * mrSges(7,1) - t131 * mrSges(5,2) - t56 * mrSges(6,2) - t17 * mrSges(7,2) + t447 + t448;
t444 = -0.2e1 * pkin(1);
t442 = Ifges(7,4) * t440 + Ifges(7,2) * t439 + Ifges(7,6) * t407;
t441 = Ifges(7,1) * t440 + Ifges(7,4) * t439 + Ifges(7,5) * t407;
t438 = Ifges(6,4) * t432 + Ifges(6,2) * t431 + Ifges(6,6) * t407;
t437 = Ifges(6,1) * t432 + Ifges(6,4) * t431 + Ifges(6,5) * t407;
t435 = t53 / 0.2e1;
t433 = t54 / 0.2e1;
t89 = Ifges(6,2) * t319 + Ifges(6,6) * t239 + t387;
t430 = -t89 / 0.2e1;
t429 = t89 / 0.2e1;
t90 = Ifges(6,1) * t158 + Ifges(6,5) * t239 + t154;
t428 = -t90 / 0.2e1;
t427 = t90 / 0.2e1;
t424 = t482 / 0.2e1;
t422 = t108 / 0.2e1;
t421 = -t319 / 0.2e1;
t420 = t319 / 0.2e1;
t418 = t158 / 0.2e1;
t417 = t171 / 0.2e1;
t416 = t172 / 0.2e1;
t414 = -t220 / 0.2e1;
t413 = -t221 / 0.2e1;
t410 = t235 / 0.2e1;
t408 = t239 / 0.2e1;
t406 = -t248 / 0.2e1;
t404 = t257 / 0.2e1;
t401 = t296 / 0.2e1;
t400 = pkin(7) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t350);
t392 = mrSges(4,3) * t256;
t244 = Ifges(4,4) * t256;
t385 = t482 * Ifges(7,6);
t384 = t108 * Ifges(7,5);
t383 = t319 * Ifges(6,6);
t382 = t158 * Ifges(6,5);
t381 = t220 * Ifges(5,6);
t379 = t221 * Ifges(5,5);
t378 = t235 * Ifges(7,3);
t377 = t239 * Ifges(6,3);
t376 = t248 * Ifges(5,3);
t374 = t257 * Ifges(4,1);
t366 = Ifges(3,5) * qJD(2);
t167 = t226 * t367 - t368 * t299;
t217 = -t368 * t279 - t280 * t367;
t363 = t167 * t217;
t336 = Ifges(5,3) * t240 + t456;
t335 = t368 * pkin(2);
t329 = t297 * t341;
t324 = t240 * mrSges(4,1) - t241 * mrSges(4,2);
t74 = t295 * t124 - t133 * t291;
t187 = t252 * t368 + t253 * t367;
t189 = pkin(3) * t258 + pkin(8) * t259 + t478;
t320 = -t187 * t292 + t296 * t189;
t286 = -t335 - pkin(3);
t317 = mrSges(5,1) * t296 - mrSges(5,2) * t292;
t314 = Ifges(5,1) * t292 + t388;
t312 = Ifges(5,2) * t296 + t389;
t310 = Ifges(5,5) * t292 + Ifges(5,6) * t296;
t58 = pkin(5) * t268 + pkin(10) * t194 + t74;
t193 = t272 * t269;
t63 = -pkin(10) * t193 + t75;
t26 = -t290 * t63 + t294 * t58;
t27 = t290 * t58 + t294 * t63;
t309 = -t292 * t72 - t296 * t73;
t186 = t252 * t367 - t368 * t253;
t308 = -t130 * t296 - t131 * t292;
t307 = t130 * t292 - t131 * t296;
t175 = -mrSges(5,2) * t248 + mrSges(5,3) * t220;
t176 = mrSges(5,1) * t248 - mrSges(5,3) * t221;
t306 = t175 * t296 - t176 * t292;
t134 = -t193 * t294 + t194 * t290;
t135 = -t193 * t290 - t194 * t294;
t185 = pkin(4) * t357 + t217;
t68 = t259 * t398 + pkin(4) * t258 + (-t213 + (pkin(9) * t269 - t206) * t292) * qJD(4) + t320;
t82 = t296 * t187 + t292 * t189 + t206 * t346 - t218 * t347;
t71 = -pkin(9) * t486 + t82;
t20 = t124 * t344 - t133 * t345 + t291 * t68 + t295 * t71;
t276 = -t296 * pkin(4) + t286;
t140 = pkin(4) * t486 + t186;
t125 = -t172 * pkin(4) + t167;
t21 = -qJD(5) * t75 - t291 * t71 + t295 * t68;
t289 = Ifges(3,4) * t349;
t264 = Ifges(3,1) * t350 + t289 + t366;
t255 = pkin(4) * t355 + t287 * t290;
t254 = -pkin(4) * t356 + t287 * t294;
t224 = -qJD(2) * mrSges(4,2) + t392;
t223 = pkin(5) * t305 + t276;
t200 = -mrSges(4,1) * t256 - mrSges(4,2) * t257;
t192 = t244 - t374 + t461;
t149 = -mrSges(5,2) * t240 + mrSges(5,3) * t172;
t148 = mrSges(5,1) * t240 - mrSges(5,3) * t171;
t145 = t376 + t379 + t381;
t142 = mrSges(6,1) * t239 - mrSges(6,3) * t158;
t141 = -mrSges(6,2) * t239 + mrSges(6,3) * t319;
t137 = t193 * pkin(5) + t185;
t136 = pkin(4) * t221 + pkin(5) * t158;
t121 = -mrSges(5,1) * t172 + mrSges(5,2) * t171;
t112 = -mrSges(6,1) * t319 + mrSges(6,2) * t158;
t95 = t171 * Ifges(5,1) + t172 * Ifges(5,4) + t240 * Ifges(5,5);
t94 = t171 * Ifges(5,4) + t172 * Ifges(5,2) + t240 * Ifges(5,6);
t93 = t194 * t450 + t272 * t259;
t92 = -t215 * t269 + t259 * t305;
t88 = t377 + t382 + t383;
t85 = mrSges(7,1) * t235 - mrSges(7,3) * t108;
t84 = -mrSges(7,2) * t235 + mrSges(7,3) * t482;
t83 = -qJD(4) * t151 + t320;
t77 = -mrSges(6,2) * t240 + mrSges(6,3) * t81;
t76 = mrSges(6,1) * t240 - mrSges(6,3) * t80;
t69 = -t93 * pkin(5) + t140;
t60 = -mrSges(7,1) * t482 + mrSges(7,2) * t108;
t52 = t378 + t384 + t385;
t51 = -t81 * pkin(5) + t125;
t41 = -mrSges(6,1) * t81 + mrSges(6,2) * t80;
t37 = -qJD(6) * t135 - t290 * t92 + t294 * t93;
t36 = qJD(6) * t134 + t290 * t93 + t294 * t92;
t29 = -mrSges(7,2) * t240 + mrSges(7,3) * t33;
t28 = mrSges(7,1) * t240 - mrSges(7,3) * t32;
t19 = t294 * t42 - t372;
t18 = -t290 * t42 - t370;
t15 = pkin(10) * t93 + t20;
t14 = pkin(5) * t258 - pkin(10) * t92 + t21;
t10 = -mrSges(7,1) * t33 + mrSges(7,2) * t32;
t5 = -qJD(6) * t27 + t14 * t294 - t15 * t290;
t4 = qJD(6) * t26 + t14 * t290 + t15 * t294;
t1 = [(t376 / 0.2e1 + t52 / 0.2e1 + t384 / 0.2e1 + t445 + t88 / 0.2e1 + t145 / 0.2e1 + t377 / 0.2e1 + t378 / 0.2e1 + t379 / 0.2e1 + t385 / 0.2e1 + t382 / 0.2e1 + t383 / 0.2e1 + t381 / 0.2e1 + t447) * t258 + (-t365 / 0.2e1 + (mrSges(3,1) * t444 - 0.3e1 / 0.2e1 * t390 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t297) * qJD(1) + (m(4) * (t275 + t351) + t200 + qJD(1) * mrSges(4,2) * t269) * pkin(2) - t449) * t348 + (-mrSges(4,3) * t241 + t121) * t217 + (t95 * t401 + t94 * t402 - Ifges(4,1) * t241 - Ifges(4,4) * t240 + t313 * t416 + t315 * t417 + (mrSges(4,3) + t316) * t167 + t309 * mrSges(5,3) + (-t296 * t146 / 0.2e1 + t147 * t402 + t197 * t317 + t312 * t414 + t314 * t413 + t310 * t406 + t307 * mrSges(5,3)) * qJD(4)) * t269 + (mrSges(4,1) * qJD(1) * t478 - mrSges(4,3) * t168 + Ifges(4,4) * t241 + Ifges(6,5) * t432 + Ifges(7,5) * t440 + Ifges(6,6) * t431 + Ifges(7,6) * t439 + (Ifges(6,3) + Ifges(7,3)) * t407 + t451 + t464 + t476) * t268 + (-Ifges(6,4) * t194 - Ifges(6,2) * t193) * t431 + t125 * (mrSges(6,1) * t193 - mrSges(6,2) * t194) + (t207 * t259 - t208 * t258 - t218 * t240) * mrSges(4,3) + t187 * t224 + t37 * t435 + t36 * t433 + (Ifges(5,3) / 0.2e1 + Ifges(4,2)) * t240 * t268 + (-Ifges(6,1) * t194 - Ifges(6,4) * t193) * t432 + (-Ifges(6,5) * t194 + Ifges(7,5) * t135 - Ifges(6,6) * t193 + Ifges(7,6) * t134 + t269 * t311) * t407 + (-t12 * t193 + t13 * t194 - t55 * t92 + t56 * t93) * mrSges(6,3) + (Ifges(7,4) * t135 + Ifges(7,2) * t134) * t439 + (Ifges(7,1) * t135 + Ifges(7,4) * t134) * t440 + t185 * t41 + t82 * t175 + t83 * t176 + t150 * t148 + t151 * t149 + t152 * (-mrSges(6,1) * t93 + mrSges(6,2) * t92) + t20 * t141 + t21 * t142 + t51 * (-mrSges(7,1) * t134 + mrSges(7,2) * t135) + t137 * t10 + t140 * t112 + t96 * (-mrSges(7,1) * t37 + mrSges(7,2) * t36) + t4 * t84 + t5 * t85 + t74 * t76 + t75 * t77 + t69 * t60 + t27 * t29 + t26 * t28 + m(6) * (t12 * t75 + t125 * t185 + t13 * t74 + t140 * t152 + t20 * t56 + t21 * t55) + m(7) * (t137 * t51 + t16 * t5 + t17 * t4 + t2 * t27 + t26 * t3 + t69 * t96) + m(4) * (t168 * t218 - t186 * t207 + t187 * t208 + t363) + m(5) * (t130 * t83 + t131 * t82 + t150 * t73 + t151 * t72 + t186 * t197 + t363) + t352 * t186 + (Ifges(7,4) * t36 + Ifges(7,2) * t37) * t424 + t92 * t427 + t93 * t429 + (Ifges(6,1) * t92 + Ifges(6,4) * t93) * t418 + (Ifges(6,4) * t92 + Ifges(6,2) * t93) * t420 + (Ifges(7,1) * t36 + Ifges(7,4) * t37) * t422 + (Ifges(6,5) * t92 + Ifges(6,6) * t93) * t408 + (Ifges(7,5) * t36 + Ifges(7,6) * t37) * t410 + (t134 * t2 - t135 * t3 - t16 * t36 + t17 * t37) * mrSges(7,3) - t194 * t437 - t193 * t438 + t135 * t441 + t134 * t442 + (t340 + t339 + t336 + t456) * t268 / 0.2e1 - (t192 / 0.2e1 + t244 / 0.2e1 - t374 / 0.2e1 + t466 + t308 * mrSges(5,3) + t446) * t259 + (-Ifges(4,5) * t259 / 0.2e1 + t258 * t479 + (t264 / 0.2e1 - t400 + t366 / 0.2e1 + (mrSges(3,2) * t444 + 0.3e1 / 0.2e1 * Ifges(3,4) * t297) * qJD(1)) * t297) * qJD(2) + t288 * t324; t51 * (-mrSges(7,1) * t210 + mrSges(7,2) * t212) + (Ifges(7,1) * t129 + Ifges(7,4) * t128) * t423 + (-mrSges(3,1) * t329 + mrSges(3,2) * t330) * pkin(7) + (mrSges(6,1) * t492 - mrSges(6,2) * t504) * t152 + (-t12 * t305 - t13 * t272 - t492 * t56 + t504 * t55) * mrSges(6,3) + t449 * t350 + (-t240 * t334 + t241 * t335) * mrSges(4,3) + (-t501 * t131 + (-t346 + t359) * t130 + t452) * mrSges(5,3) + (Ifges(7,4) * t129 + Ifges(7,2) * t128) * t425 + ((-t167 * t368 + t168 * t367) * pkin(2) + t207 * t209 - t208 * t211 - t275 * t338) * m(4) + (t244 + t192 + t353) * t405 + (Ifges(6,5) * t272 + Ifges(7,5) * t212 - Ifges(6,6) * t305 + Ifges(7,6) * t210 + t310) * t407 + (Ifges(6,1) * t272 - Ifges(6,4) * t305) * t432 + (Ifges(6,4) * t272 - Ifges(6,2) * t305) * t431 + t125 * (mrSges(6,1) * t305 + mrSges(6,2) * t272) - t305 * t438 + (-Ifges(6,1) * t214 - Ifges(6,4) * t215) * t418 + (-Ifges(6,4) * t214 - Ifges(6,2) * t215) * t420 + (-Ifges(6,5) * t214 - Ifges(6,6) * t215) * t408 + (Ifges(7,5) * t129 + Ifges(7,6) * t128) * t411 + (m(5) * t285 * t308 + t446) * qJD(4) + t223 * t10 - t211 * t224 + t128 * t436 + t272 * t437 + t122 * t433 + t129 * t434 + t123 * t435 - m(5) * (t130 * t138 + t131 * t139) + (-Ifges(6,5) * t184 - Ifges(6,6) * t183) * t409 + (-Ifges(6,4) * t184 - Ifges(6,2) * t183) * t421 + (-Ifges(6,1) * t184 - Ifges(6,4) * t183) * t419 + (t373 + t145 + t88 + t52) * t404 + t201 * t76 + t202 * t77 - (-Ifges(3,2) * t350 + t264 + t289) * t349 / 0.2e1 + (-t293 * (Ifges(3,1) * t297 - t390) / 0.2e1 + pkin(1) * (mrSges(3,1) * t293 + mrSges(3,2) * t297)) * qJD(1) ^ 2 - t139 * t175 - t138 * t176 - t168 * mrSges(4,2) + t115 * t28 + t116 * t29 - t208 * t391 + (Ifges(7,4) * t122 + Ifges(7,2) * t123) * t424 - t214 * t427 - t184 * t428 - t215 * t429 - t183 * t430 + t312 * t416 + t314 * t417 + (Ifges(7,1) * t122 + Ifges(7,4) * t123) * t422 + (Ifges(7,5) * t122 + Ifges(7,6) * t123) * t410 + t349 * t400 + t94 * t401 + t207 * t392 - (Ifges(3,5) * t297 - Ifges(3,6) * t293) * t341 / 0.2e1 - t200 * t338 + (Ifges(7,4) * t212 + Ifges(7,2) * t210) * t439 + (Ifges(7,1) * t212 + Ifges(7,4) * t210) * t440 + t212 * t441 + t210 * t442 - Ifges(4,6) * t240 - Ifges(4,5) * t241 + (m(5) * t452 - t292 * t148 + t296 * t149 - t175 * t347 - t176 * t346) * t285 + (-m(5) * t197 - t352) * t209 + t455 * t60 + t457 * t112 + (-t16 * t458 + t17 * t459 + t2 * t210 - t212 * t3) * mrSges(7,3) + (-mrSges(7,1) * t459 + mrSges(7,2) * t458) * t96 + (-Ifges(5,5) * t413 - Ifges(6,5) * t419 - Ifges(7,5) * t423 - Ifges(5,6) * t414 - Ifges(6,6) * t421 - Ifges(7,6) * t425 - Ifges(5,3) * t406 - Ifges(6,3) * t409 - Ifges(7,3) * t411 + t445 + t448 + t496) * t257 + Ifges(3,5) * t329 + (t311 * t406 + t315 * t413 + t313 * t414 + Ifges(4,1) * t404 - t461 / 0.2e1 - t466 - t485) * t256 + t462 * t141 + t463 * t142 + (t12 * t202 + t125 * t276 + t13 * t201 + t152 * t457 + t462 * t56 + t463 * t55) * m(6) + t276 * t41 + t286 * t121 + t474 * t85 + t475 * t84 + (t115 * t3 + t116 * t2 + t16 * t474 + t17 * t475 + t223 * t51 + t455 * t96) * m(7) + t292 * t95 / 0.2e1 + (m(5) * t286 - mrSges(4,1) - t317) * t167 - Ifges(3,6) * t330; t210 * t28 + t212 * t29 + t324 + (-t224 - t306) * t256 + t306 * qJD(4) + t459 * t85 + t458 * t84 + (t112 + t60 + t352) * t257 - t492 * t142 - t504 * t141 - t305 * t76 + t272 * t77 + t292 * t149 + t296 * t148 + (t16 * t459 + t17 * t458 + t2 * t212 + t210 * t3 + t257 * t96) * m(7) + (t12 * t272 - t13 * t305 + t152 * t257 - t492 * t55 - t504 * t56) * m(6) + (t197 * t257 - t248 * t307 - t309) * m(5) + (-t207 * t257 - t208 * t256 + t318) * m(4); t451 - t158 * t430 + t336 + (t130 * t220 + t131 * t221) * mrSges(5,3) - t197 * (mrSges(5,1) * t221 + mrSges(5,2) * t220) + (-Ifges(5,2) * t221 + t147 + t219) * t414 - m(6) * (t55 * t61 + t56 * t62) + t319 * t428 - t130 * t175 + t131 * t176 - t62 * t141 - t61 * t142 - t136 * t60 + (-t221 * t112 + t291 * t77 + t295 * t76 + (t141 * t295 - t142 * t291) * qJD(5) + (t12 * t291 + t13 * t295 - t152 * t221 + t344 * t56 - t345 * t55) * m(6)) * pkin(4) + (Ifges(5,5) * t220 - Ifges(5,6) * t221) * t406 + t146 * t412 + (Ifges(5,1) * t220 - t380) * t413 + t254 * t28 + t255 * t29 + t483 * t421 + t505 + t487 * t85 + t488 * t84 + (-t136 * t96 + t16 * t487 + t17 * t488 + t2 * t255 + t254 * t3) * m(7); -m(7) * (t16 * t18 + t17 * t19) - t55 * t141 + t56 * t142 - t19 * t84 - t18 * t85 + (-t158 * t60 + t294 * t28 + t290 * t29 + (-t290 * t85 + t294 * t84) * qJD(6) + (-t158 * t96 - t16 * t343 + t17 * t342 + t2 * t290 + t294 * t3) * m(7)) * pkin(5) + (t483 + t90) * t421 + t89 * t418 + t505; t53 * t422 - t16 * t84 + t17 * t85 + (t484 + t54) * t425 + t508;];
tauc  = t1(:);
