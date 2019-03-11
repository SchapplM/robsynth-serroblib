% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:42:37
% EndTime: 2019-03-09 00:43:02
% DurationCPUTime: 12.65s
% Computational Cost: add. (13426->636), mult. (32680->882), div. (0->0), fcn. (24562->12), ass. (0->317)
t295 = sin(qJ(4));
t296 = sin(qJ(3));
t300 = cos(qJ(4));
t301 = cos(qJ(3));
t258 = t295 * t296 - t300 * t301;
t290 = qJD(3) + qJD(4);
t216 = t290 * t258;
t260 = t295 * t301 + t296 * t300;
t217 = t290 * t260;
t297 = sin(qJ(2));
t291 = sin(pkin(6));
t364 = qJD(1) * t291;
t347 = t297 * t364;
t358 = qJD(3) * t296;
t353 = pkin(3) * t358;
t487 = pkin(4) * t217 + pkin(10) * t216 - t347 + t353;
t431 = -pkin(9) - pkin(8);
t349 = qJD(3) * t431;
t266 = t296 * t349;
t267 = t301 * t349;
t302 = cos(qJ(2));
t346 = t302 * t364;
t275 = t431 * t296;
t277 = t431 * t301;
t446 = t300 * t275 + t277 * t295;
t462 = qJD(4) * t446 + t258 * t346 + t266 * t300 + t267 * t295;
t287 = -pkin(3) * t301 - pkin(2);
t208 = pkin(4) * t258 - pkin(10) * t260 + t287;
t232 = t275 * t295 - t277 * t300;
t294 = sin(qJ(5));
t299 = cos(qJ(5));
t356 = qJD(5) * t299;
t357 = qJD(5) * t294;
t468 = t208 * t356 - t232 * t357 + t294 * t487 + t462 * t299;
t486 = -t462 * t294 + t299 * t487;
t412 = -t294 / 0.2e1;
t221 = t299 * t232;
t289 = t299 * pkin(11);
t485 = t216 * t289 + pkin(5) * t217 + (-t221 + (pkin(11) * t260 - t208) * t294) * qJD(5) + t486;
t314 = -t216 * t294 + t260 * t356;
t484 = pkin(11) * t314 - t468;
t284 = pkin(3) * t295 + pkin(10);
t408 = -pkin(11) - t284;
t337 = qJD(5) * t408;
t401 = pkin(3) * qJD(4);
t352 = t300 * t401;
t359 = qJD(2) * t301;
t361 = qJD(2) * t296;
t251 = -t295 * t361 + t300 * t359;
t376 = t251 * t294;
t354 = pkin(11) * t376;
t269 = qJD(2) * pkin(8) + t347;
t336 = pkin(9) * qJD(2) + t269;
t292 = cos(pkin(6));
t363 = qJD(1) * t296;
t343 = t292 * t363;
t223 = t301 * t336 + t343;
t210 = t295 * t223;
t369 = t292 * t301;
t280 = qJD(1) * t369;
t321 = t336 * t296;
t222 = t280 - t321;
t142 = t222 * t300 - t210;
t252 = t260 * qJD(2);
t207 = pkin(4) * t252 - pkin(10) * t251;
t180 = pkin(3) * t361 + t207;
t83 = t299 * t142 + t294 * t180;
t483 = t294 * t337 + t299 * t352 + t354 - t83;
t375 = t251 * t299;
t332 = t252 * pkin(5) - pkin(11) * t375;
t82 = -t142 * t294 + t299 * t180;
t482 = -t294 * t352 + t299 * t337 - t332 - t82;
t481 = t290 * Ifges(5,6) / 0.2e1;
t430 = -pkin(11) - pkin(10);
t348 = qJD(5) * t430;
t212 = qJD(3) * pkin(3) + t222;
t136 = t212 * t300 - t210;
t88 = t299 * t136 + t294 * t207;
t480 = t294 * t348 + t354 - t88;
t87 = -t136 * t294 + t299 * t207;
t479 = t299 * t348 - t332 - t87;
t227 = -t252 * t294 + t290 * t299;
t245 = qJD(5) - t251;
t228 = t252 * t299 + t290 * t294;
t394 = t228 * Ifges(6,4);
t123 = t227 * Ifges(6,2) + t245 * Ifges(6,6) + t394;
t127 = -pkin(4) * t290 - t136;
t330 = mrSges(6,1) * t294 + mrSges(6,2) * t299;
t478 = t123 * t412 + t127 * t330;
t243 = qJD(2) * t287 - t346;
t386 = t290 * Ifges(5,5);
t477 = t243 * mrSges(5,2) + t386 / 0.2e1;
t387 = t252 * Ifges(5,4);
t476 = t481 + t387 / 0.2e1 + t251 * Ifges(5,2) / 0.2e1;
t293 = sin(qJ(6));
t298 = cos(qJ(6));
t150 = t227 * t293 + t228 * t298;
t202 = t217 * qJD(2);
t197 = Ifges(7,3) * t202;
t211 = t300 * t223;
t137 = t212 * t295 + t211;
t128 = pkin(10) * t290 + t137;
t164 = -pkin(4) * t251 - pkin(10) * t252 + t243;
t76 = t128 * t299 + t164 * t294;
t65 = pkin(11) * t227 + t76;
t384 = t293 * t65;
t75 = -t128 * t294 + t299 * t164;
t64 = -pkin(11) * t228 + t75;
t50 = pkin(5) * t245 + t64;
t20 = t298 * t50 - t384;
t382 = t298 * t65;
t21 = t293 * t50 + t382;
t334 = t298 * t227 - t228 * t293;
t402 = Ifges(7,4) * t150;
t241 = qJD(6) + t245;
t417 = -t241 / 0.2e1;
t423 = -t150 / 0.2e1;
t93 = -pkin(5) * t227 + t127;
t475 = t197 + (Ifges(7,5) * t334 - Ifges(7,6) * t150) * t417 + (t150 * t21 + t20 * t334) * mrSges(7,3) - t93 * (mrSges(7,1) * t150 + mrSges(7,2) * t334) + (Ifges(7,1) * t334 - t402) * t423;
t140 = t294 * t208 + t221;
t374 = t260 * t294;
t113 = -pkin(11) * t374 + t140;
t139 = t299 * t208 - t232 * t294;
t94 = pkin(5) * t258 - t260 * t289 + t139;
t45 = t113 * t298 + t293 * t94;
t474 = -qJD(6) * t45 + t484 * t293 + t485 * t298;
t44 = -t113 * t293 + t298 * t94;
t473 = qJD(6) * t44 + t485 * t293 - t484 * t298;
t255 = t408 * t294;
t256 = t284 * t299 + t289;
t203 = t255 * t298 - t256 * t293;
t472 = qJD(6) * t203 + t482 * t293 + t483 * t298;
t204 = t255 * t293 + t256 * t298;
t471 = -qJD(6) * t204 - t483 * t293 + t482 * t298;
t469 = -qJD(5) * t140 + t486;
t362 = qJD(2) * t291;
t338 = qJD(1) * t362;
t333 = t302 * t338;
t365 = qJD(3) * t280 + t301 * t333;
t169 = -qJD(3) * t321 + t365;
t445 = -t223 * qJD(3) - t296 * t333;
t61 = qJD(4) * t137 + t169 * t295 - t300 * t445;
t201 = t216 * qJD(2);
t133 = qJD(5) * t227 - t201 * t299;
t134 = -qJD(5) * t228 + t201 * t294;
t73 = -mrSges(6,1) * t134 + mrSges(6,2) * t133;
t467 = -m(6) * t61 - t73;
t154 = qJD(4) * t232 + t266 * t295 - t300 * t267;
t225 = t260 * t346;
t466 = t154 - t225;
t259 = t293 * t299 + t294 * t298;
t172 = t259 * t251;
t443 = qJD(5) + qJD(6);
t215 = t443 * t259;
t465 = t172 - t215;
t318 = t293 * t294 - t298 * t299;
t173 = t318 * t251;
t214 = t443 * t318;
t464 = t173 - t214;
t206 = -mrSges(5,1) * t251 + mrSges(5,2) * t252;
t463 = -m(5) * t243 - t206;
t141 = t222 * t295 + t211;
t239 = pkin(5) * t376;
t351 = pkin(5) * t357;
t461 = t295 * t401 - t141 - t239 + t351;
t246 = qJD(2) * t353 + t297 * t338;
t102 = pkin(4) * t202 + pkin(10) * t201 + t246;
t60 = qJD(4) * t136 + t300 * t169 + t295 * t445;
t17 = t294 * t102 - t128 * t357 + t164 * t356 + t299 * t60;
t18 = -qJD(5) * t76 + t299 * t102 - t294 * t60;
t460 = t17 * t299 - t18 * t294;
t10 = pkin(11) * t134 + t17;
t7 = pkin(5) * t202 - pkin(11) * t133 + t18;
t3 = qJD(6) * t20 + t10 * t298 + t293 * t7;
t38 = qJD(6) * t334 + t133 * t298 + t134 * t293;
t39 = -qJD(6) * t150 - t133 * t293 + t134 * t298;
t4 = -qJD(6) * t21 - t10 * t293 + t298 * t7;
t459 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + Ifges(7,5) * t38 + Ifges(7,6) * t39;
t146 = Ifges(7,4) * t334;
t458 = -Ifges(7,2) * t150 + t146;
t325 = Ifges(6,5) * t299 - Ifges(6,6) * t294;
t403 = Ifges(6,4) * t299;
t327 = -Ifges(6,2) * t294 + t403;
t404 = Ifges(6,4) * t294;
t329 = Ifges(6,1) * t299 - t404;
t224 = Ifges(6,4) * t227;
t124 = t228 * Ifges(6,1) + t245 * Ifges(6,5) + t224;
t367 = t299 * t124;
t418 = t228 / 0.2e1;
t457 = t367 / 0.2e1 + t245 * t325 / 0.2e1 + t329 * t418 + t227 * t327 / 0.2e1 + t478;
t456 = -t243 * mrSges(5,1) - t75 * mrSges(6,1) - t20 * mrSges(7,1) + t76 * mrSges(6,2) + t21 * mrSges(7,2) + t476;
t454 = -Ifges(4,1) / 0.2e1;
t437 = t38 / 0.2e1;
t436 = t39 / 0.2e1;
t421 = t202 / 0.2e1;
t453 = -Ifges(4,4) * t359 / 0.2e1;
t274 = t430 * t294;
t276 = pkin(10) * t299 + t289;
t229 = t274 * t298 - t276 * t293;
t448 = qJD(6) * t229 + t293 * t479 + t298 * t480;
t231 = t274 * t293 + t276 * t298;
t447 = -qJD(6) * t231 - t293 * t480 + t298 * t479;
t190 = t318 * t260;
t444 = -t294 * t75 + t299 * t76;
t441 = t18 * mrSges(6,1) - t17 * mrSges(6,2) + Ifges(6,5) * t133 + Ifges(6,6) * t134 + t459;
t439 = Ifges(7,4) * t437 + Ifges(7,2) * t436 + Ifges(7,6) * t421;
t438 = Ifges(7,1) * t437 + Ifges(7,4) * t436 + Ifges(7,5) * t421;
t71 = Ifges(7,2) * t334 + Ifges(7,6) * t241 + t402;
t435 = -t71 / 0.2e1;
t434 = t71 / 0.2e1;
t72 = Ifges(7,1) * t150 + Ifges(7,5) * t241 + t146;
t433 = -t72 / 0.2e1;
t432 = t72 / 0.2e1;
t427 = t133 / 0.2e1;
t426 = t134 / 0.2e1;
t425 = -t334 / 0.2e1;
t424 = t334 / 0.2e1;
t422 = t150 / 0.2e1;
t420 = -t227 / 0.2e1;
t419 = -t228 / 0.2e1;
t416 = t241 / 0.2e1;
t415 = -t245 / 0.2e1;
t411 = t299 / 0.2e1;
t409 = pkin(3) * t300;
t407 = mrSges(5,3) * t251;
t406 = mrSges(5,3) * t252;
t405 = Ifges(4,4) * t296;
t242 = Ifges(5,4) * t251;
t400 = t334 * Ifges(7,6);
t399 = t150 * Ifges(7,5);
t371 = t291 * t297;
t247 = -t296 * t371 + t369;
t248 = t292 * t296 + t301 * t371;
t319 = t300 * t247 - t248 * t295;
t396 = t319 * t61;
t395 = t227 * Ifges(6,6);
t393 = t228 * Ifges(6,5);
t392 = t446 * t61;
t391 = t241 * Ifges(7,3);
t390 = t245 * Ifges(6,3);
t388 = t252 * Ifges(5,1);
t380 = Ifges(4,5) * qJD(3);
t379 = Ifges(4,6) * qJD(3);
t373 = t269 * t301;
t370 = t291 * t302;
t366 = mrSges(5,1) * t290 + mrSges(6,1) * t227 - mrSges(6,2) * t228 - t406;
t360 = qJD(2) * t297;
t355 = qJD(2) * qJD(3);
t81 = -mrSges(7,1) * t334 + mrSges(7,2) * t150;
t350 = -t81 + t366;
t286 = -pkin(5) * t299 - pkin(4);
t345 = t291 * t360;
t344 = t302 * t362;
t342 = t380 / 0.2e1;
t341 = -t379 / 0.2e1;
t331 = mrSges(6,1) * t299 - mrSges(6,2) * t294;
t328 = Ifges(6,1) * t294 + t403;
t326 = Ifges(6,2) * t299 + t404;
t324 = Ifges(6,5) * t294 + Ifges(6,6) * t299;
t323 = -t294 * t76 - t299 * t75;
t89 = mrSges(6,1) * t202 - mrSges(6,3) * t133;
t90 = -mrSges(6,2) * t202 + mrSges(6,3) * t134;
t322 = -t294 * t89 + t299 * t90;
t185 = t247 * t295 + t248 * t300;
t162 = -t185 * t294 - t299 * t370;
t316 = -t185 * t299 + t294 * t370;
t84 = t162 * t298 + t293 * t316;
t85 = t162 * t293 - t298 * t316;
t191 = -t269 * t358 + t365;
t192 = -qJD(3) * t373 + (-qJD(3) * t292 - t344) * t363;
t320 = t191 * t301 - t192 * t296;
t315 = t323 * mrSges(6,3);
t235 = -t269 * t296 + t280;
t270 = -qJD(2) * pkin(2) - t346;
t312 = t235 * mrSges(4,3) + t361 * t454 + t453 - t380 / 0.2e1 - t270 * mrSges(4,2);
t236 = t343 + t373;
t311 = t236 * mrSges(4,3) + t379 / 0.2e1 + (t301 * Ifges(4,2) + t405) * qJD(2) / 0.2e1 - t270 * mrSges(4,1);
t305 = m(6) * (qJD(5) * t323 + t460);
t122 = t390 + t393 + t395;
t182 = t242 + t386 + t388;
t35 = -pkin(5) * t134 + t61;
t48 = t133 * Ifges(6,4) + t134 * Ifges(6,2) + t202 * Ifges(6,6);
t49 = t133 * Ifges(6,1) + t134 * Ifges(6,4) + t202 * Ifges(6,5);
t70 = t391 + t399 + t400;
t304 = (-t20 * t464 + t21 * t465 - t259 * t4 - t318 * t3) * mrSges(7,3) + (-mrSges(7,1) * t465 + mrSges(7,2) * t464) * t93 + (t375 * t75 + t376 * t76 + t460) * mrSges(6,3) + t48 * t411 + t457 * qJD(5) + (Ifges(6,5) * t419 + Ifges(7,5) * t423 + Ifges(6,6) * t420 + Ifges(7,6) * t425 + Ifges(6,3) * t415 + Ifges(7,3) * t417 + t456 + t481) * t252 + (t325 * t415 + t327 * t420 + t329 * t419 - t477 - t478) * t251 + (-Ifges(7,4) * t214 - Ifges(7,2) * t215) * t424 - (Ifges(5,1) * t251 + t122 - t387 + t70) * t252 / 0.2e1 - (-Ifges(5,2) * t252 + t182 + t242 + t367) * t251 / 0.2e1 + (Ifges(7,4) * t259 - Ifges(7,2) * t318) * t436 + (Ifges(7,1) * t259 - Ifges(7,4) * t318) * t437 + t35 * (mrSges(7,1) * t318 + mrSges(7,2) * t259) - t318 * t439 + (Ifges(7,5) * t259 - Ifges(7,6) * t318 + t324) * t421 + (-Ifges(7,5) * t214 - Ifges(7,6) * t215) * t416 + (-Ifges(7,1) * t214 - Ifges(7,4) * t215) * t422 + t326 * t426 + t328 * t427 - t214 * t432 - t173 * t433 - t215 * t434 - t172 * t435 + t259 * t438 + t294 * t49 / 0.2e1 + (-mrSges(5,1) - t331) * t61 + t136 * t407 - Ifges(5,5) * t201 - Ifges(5,6) * t202 - t60 * mrSges(5,2) + (-Ifges(7,5) * t173 - Ifges(7,6) * t172) * t417 + (-Ifges(7,1) * t173 - Ifges(7,4) * t172) * t423 + (-Ifges(7,4) * t173 - Ifges(7,2) * t172) * t425;
t303 = qJD(2) ^ 2;
t273 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t359;
t272 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t361;
t271 = t286 - t409;
t254 = (mrSges(4,1) * t296 + mrSges(4,2) * t301) * t355;
t237 = -mrSges(5,2) * t290 + t407;
t220 = qJD(3) * t247 + t301 * t344;
t219 = -qJD(3) * t248 - t296 * t344;
t198 = Ifges(6,3) * t202;
t189 = t259 * t260;
t174 = pkin(5) * t374 - t446;
t167 = mrSges(6,1) * t245 - mrSges(6,3) * t228;
t166 = -mrSges(6,2) * t245 + mrSges(6,3) * t227;
t121 = mrSges(5,1) * t202 - mrSges(5,2) * t201;
t115 = mrSges(7,1) * t241 - mrSges(7,3) * t150;
t114 = -mrSges(7,2) * t241 + mrSges(7,3) * t334;
t106 = t137 + t239;
t92 = qJD(4) * t185 - t300 * t219 + t220 * t295;
t91 = qJD(4) * t319 + t219 * t295 + t220 * t300;
t86 = pkin(5) * t314 + t154;
t67 = t190 * t443 + t259 * t216;
t66 = -t215 * t260 + t216 * t318;
t55 = qJD(5) * t316 - t294 * t91 + t299 * t345;
t54 = qJD(5) * t162 + t294 * t345 + t299 * t91;
t32 = -mrSges(7,2) * t202 + mrSges(7,3) * t39;
t31 = mrSges(7,1) * t202 - mrSges(7,3) * t38;
t23 = t298 * t64 - t384;
t22 = -t293 * t64 - t382;
t15 = -mrSges(7,1) * t39 + mrSges(7,2) * t38;
t9 = -qJD(6) * t85 - t293 * t54 + t298 * t55;
t8 = qJD(6) * t84 + t293 * t55 + t298 * t54;
t1 = [-t185 * t202 * mrSges(5,3) + t8 * t114 + t9 * t115 + t162 * t89 - t316 * t90 + t54 * t166 + t55 * t167 + t219 * t272 + t220 * t273 + t91 * t237 + t84 * t31 + t85 * t32 + (-t247 * t301 - t248 * t296) * mrSges(4,3) * t355 - t350 * t92 - (-mrSges(5,3) * t201 + t15 + t73) * t319 + ((-mrSges(3,2) * t303 - t121 - t254) * t302 + (-mrSges(3,1) * t303 + (t206 + qJD(2) * (-mrSges(4,1) * t301 + mrSges(4,2) * t296)) * qJD(2)) * t297) * t291 + m(7) * (t20 * t9 + t21 * t8 + t3 * t85 - t319 * t35 + t4 * t84 + t92 * t93) + m(6) * (t127 * t92 + t162 * t18 - t17 * t316 + t54 * t76 + t55 * t75 - t396) + m(5) * (-t136 * t92 + t137 * t91 - t396 + t185 * t60 + (t243 * t360 - t246 * t302) * t291) + m(4) * (t191 * t248 + t192 * t247 + t219 * t235 + t220 * t236 + (t270 - t346) * t345); t462 * t237 + (-(pkin(2) * t360 + (-t235 * t296 + t236 * t301) * t302) * m(4) + (t272 * t296 - t273 * t301) * t302 + (-t270 * m(4) + t463) * t297) * t364 + (-t136 * t466 + t137 * t462 + t232 * t60 + t243 * t353 + t246 * t287 - t392) * m(5) + (t390 / 0.2e1 + t393 / 0.2e1 + t395 / 0.2e1 + t391 / 0.2e1 + t400 / 0.2e1 + t399 / 0.2e1 + t122 / 0.2e1 + t70 / 0.2e1 - t456 - t476) * t217 + m(4) * ((-t235 * t301 - t236 * t296) * qJD(3) + t320) * pkin(8) + (Ifges(7,5) * t66 + Ifges(7,6) * t67) * t416 + (-0.3e1 / 0.2e1 * t296 ^ 2 + 0.3e1 / 0.2e1 * t301 ^ 2) * Ifges(4,4) * t355 + (-Ifges(7,5) * t190 - Ifges(7,6) * t189) * t421 + (-t189 * t3 + t190 * t4 - t20 * t66 + t21 * t67) * mrSges(7,3) + (-Ifges(7,4) * t190 - Ifges(7,2) * t189) * t436 + (-Ifges(7,1) * t190 - Ifges(7,4) * t189) * t437 + t35 * (mrSges(7,1) * t189 - mrSges(7,2) * t190) + (t136 * t216 - t137 * t217 + t201 * t446 - t202 * t232) * mrSges(5,3) + (t325 * t421 + t329 * t427 + t327 * t426 + t246 * mrSges(5,2) - Ifges(5,4) * t202 - Ifges(5,1) * t201 + t48 * t412 + t49 * t411 + (mrSges(5,3) + t330) * t61 + (-t17 * t294 - t18 * t299) * mrSges(6,3) + (t124 * t412 + t324 * t415 + t328 * t419 + t326 * t420 + t127 * t331 - t299 * t123 / 0.2e1 - t444 * mrSges(6,3)) * qJD(5)) * t260 - (t242 / 0.2e1 + t388 / 0.2e1 + t182 / 0.2e1 + t315 + t457 + t477) * t216 + t468 * t166 + t469 * t167 + (t127 * t466 + t139 * t18 + t140 * t17 + t468 * t76 + t469 * t75 - t392) * m(6) - t446 * t73 + (Ifges(5,4) * t201 + t246 * mrSges(5,1) + t198 / 0.2e1 + t197 / 0.2e1 - t60 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) + Ifges(7,3) / 0.2e1) * t202 + t441) * t258 + (Ifges(7,1) * t66 + Ifges(7,4) * t67) * t422 + (Ifges(7,4) * t66 + Ifges(7,2) * t67) * t424 + t66 * t432 + t67 * t434 - t190 * t438 - t189 * t439 + t473 * t114 + t474 * t115 + (t174 * t35 + t3 * t45 + t4 * t44 + (-t225 + t86) * t93 + t473 * t21 + t474 * t20) * m(7) + t287 * t121 - pkin(2) * t254 + t174 * t15 + t139 * t89 + t140 * t90 + t93 * (-mrSges(7,1) * t67 + mrSges(7,2) * t66) + t86 * t81 + t44 * t31 + t45 * t32 + t320 * mrSges(4,3) + t350 * t225 + ((-pkin(8) * t272 - t312 + t342) * t301 + (-pkin(8) * t273 + pkin(3) * t206 + t341 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t359 - t311) * t296) * qJD(3) - t366 * t154; t461 * t81 + t304 + t471 * t115 + t472 * t114 + t271 * t15 + t236 * t272 - t235 * t273 - t142 * t237 + t137 * t406 + t203 * t31 + t204 * t32 - t191 * mrSges(4,2) + t192 * mrSges(4,1) - t83 * t166 - t82 * t167 - m(5) * (-t136 * t141 + t137 * t142) - m(6) * (t127 * t141 + t75 * t82 + t76 * t83) + t315 * qJD(5) + (m(5) * (t295 * t60 - t300 * t61) + (t201 * t300 - t202 * t295) * mrSges(5,3) + ((-m(5) * t136 + m(6) * t127 - t366) * t295 + (m(5) * t137 + m(6) * t444 + t166 * t299 - t167 * t294 + t237) * t300) * qJD(4)) * pkin(3) + ((t342 + t453 + t312) * t301 + (t341 + (t405 / 0.2e1 + (t454 + Ifges(4,2) / 0.2e1) * t301) * qJD(2) + t463 * pkin(3) + t311) * t296) * qJD(2) + t366 * t141 - t467 * (-pkin(4) - t409) + ((-t166 * t294 - t167 * t299) * qJD(5) + t322 + t305) * t284 + (t20 * t471 + t203 * t4 + t204 * t3 + t21 * t472 + t271 * t35 + t461 * t93) * m(7); t304 + ((-t75 * mrSges(6,3) - pkin(10) * t167) * t299 + (-t76 * mrSges(6,3) + pkin(5) * t81 - pkin(10) * t166) * t294) * qJD(5) + t322 * pkin(10) + t286 * t15 + (t366 + t406) * t137 - t136 * t237 + t229 * t31 + t231 * t32 - t88 * t166 - t87 * t167 - m(6) * (t127 * t137 + t75 * t87 + t76 * t88) - t106 * t81 + pkin(10) * t305 + t448 * t114 + t447 * t115 + t467 * pkin(4) + (t229 * t4 + t231 * t3 + t286 * t35 + (-t106 + t351) * t93 + t448 * t21 + t447 * t20) * m(7); t334 * t433 + (-Ifges(6,2) * t228 + t124 + t224) * t420 + (Ifges(6,5) * t227 - Ifges(6,6) * t228) * t415 + t123 * t418 + (Ifges(6,1) * t227 - t394) * t419 - m(7) * (t20 * t22 + t21 * t23) - t150 * t435 + t198 + t458 * t425 - t127 * (mrSges(6,1) * t228 + mrSges(6,2) * t227) - t75 * t166 + t76 * t167 - t23 * t114 - t22 * t115 + (-t228 * t81 + t293 * t32 + t298 * t31 + (t114 * t298 - t115 * t293) * qJD(6) + (-t228 * t93 + t293 * t3 + t298 * t4 + (-t20 * t293 + t21 * t298) * qJD(6)) * m(7)) * pkin(5) + t441 + (t227 * t75 + t228 * t76) * mrSges(6,3) + t475; -t20 * t114 + t21 * t115 + t71 * t422 + (t458 + t72) * t425 + t459 + t475;];
tauc  = t1(:);
