% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:34
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:34:25
% EndTime: 2018-11-23 15:34:44
% DurationCPUTime: 19.75s
% Computational Cost: add. (18716->784), mult. (49304->1120), div. (0->0), fcn. (40242->14), ass. (0->372)
t290 = sin(qJ(5));
t291 = sin(qJ(4));
t295 = cos(qJ(5));
t296 = cos(qJ(4));
t259 = t290 * t291 - t295 * t296;
t474 = qJD(4) + qJD(5);
t210 = t474 * t259;
t297 = cos(qJ(3));
t285 = sin(pkin(7));
t373 = qJD(2) * t285;
t355 = t297 * t373;
t217 = t259 * t355;
t501 = t210 - t217;
t260 = t290 * t296 + t291 * t295;
t211 = t474 * t260;
t216 = t260 * t355;
t377 = t211 - t216;
t461 = -pkin(11) - pkin(10);
t359 = qJD(4) * t461;
t262 = t291 * t359;
t273 = t461 * t291;
t274 = t461 * t296;
t319 = t295 * t273 + t274 * t290;
t346 = t296 * t359;
t155 = qJD(5) * t319 + t295 * t262 + t290 * t346;
t293 = sin(qJ(2));
t286 = sin(pkin(6));
t375 = qJD(1) * t286;
t357 = t293 * t375;
t258 = pkin(9) * t373 + t357;
t292 = sin(qJ(3));
t250 = t292 * t258;
t298 = cos(qJ(2));
t264 = qJD(2) * pkin(2) + t298 * t375;
t287 = cos(pkin(7));
t288 = cos(pkin(6));
t374 = qJD(1) * t288;
t358 = t285 * t374;
t181 = t297 * (t264 * t287 + t358) - t250;
t337 = pkin(3) * t292 - pkin(10) * t297;
t245 = t337 * t373;
t141 = -t181 * t291 + t296 * t245;
t125 = (-pkin(11) * t296 * t297 + pkin(4) * t292) * t373 + t141;
t142 = t296 * t181 + t291 * t245;
t343 = t291 * t355;
t133 = -pkin(11) * t343 + t142;
t73 = t290 * t125 + t295 * t133;
t500 = t155 - t73;
t387 = t287 * t292;
t389 = t285 * t297;
t257 = pkin(2) * t387 + pkin(9) * t389;
t239 = pkin(10) * t287 + t257;
t338 = -pkin(3) * t297 - pkin(10) * t292;
t240 = (-pkin(2) + t338) * t285;
t313 = t337 * qJD(3);
t246 = t285 * t313;
t390 = t285 * t292;
t278 = pkin(9) * t390;
t386 = t287 * t297;
t256 = pkin(2) * t386 - t278;
t247 = t256 * qJD(3);
t369 = qJD(4) * t296;
t370 = qJD(4) * t291;
t123 = -t239 * t370 + t240 * t369 + t291 * t246 + t296 * t247;
t254 = t287 * t291 + t296 * t390;
t371 = qJD(3) * t297;
t354 = t291 * t371;
t213 = -qJD(4) * t254 - t285 * t354;
t105 = pkin(11) * t213 + t123;
t183 = -t239 * t291 + t296 * t240;
t143 = -pkin(4) * t389 - pkin(11) * t254 + t183;
t184 = t296 * t239 + t291 * t240;
t253 = t287 * t296 - t291 * t390;
t154 = pkin(11) * t253 + t184;
t379 = t297 * t298;
t384 = t292 * t293;
t315 = -t287 * t384 + t379;
t223 = t315 * t375;
t345 = t285 * t357;
t188 = -t223 * t291 + t296 * t345;
t189 = t223 * t296 + t291 * t345;
t367 = qJD(5) * t295;
t368 = qJD(5) * t290;
t124 = -qJD(4) * t184 + t296 * t246 - t247 * t291;
t353 = t296 * t371;
t212 = qJD(4) * t253 + t285 * t353;
t372 = qJD(3) * t285;
t352 = t292 * t372;
t97 = pkin(4) * t352 - pkin(11) * t212 + t124;
t483 = t143 * t367 - t154 * t368 + (t105 - t189) * t295 + (-t188 + t97) * t290;
t356 = t292 * t373;
t499 = -pkin(12) * t356 + t500;
t479 = t297 * t258 + t264 * t387;
t182 = t292 * t358 + t479;
t163 = pkin(4) * t343 + t182;
t498 = pkin(4) * t370 + t377 * pkin(5) + t501 * pkin(12) - t163;
t248 = t257 * qJD(3);
t179 = -pkin(4) * t213 + t248;
t382 = t293 * t297;
t383 = t292 * t298;
t317 = t287 * t382 + t383;
t222 = t317 * t375;
t497 = t179 - t222;
t496 = -pkin(12) * t352 - t483;
t320 = t295 * t253 - t254 * t290;
t116 = qJD(5) * t320 + t212 * t295 + t213 * t290;
t195 = t253 * t290 + t254 * t295;
t117 = qJD(5) * t195 + t212 * t290 - t295 * t213;
t495 = pkin(5) * t117 - pkin(12) * t116 + t497;
t271 = qJD(4) - t355;
t265 = qJD(5) + t271;
t289 = sin(qJ(6));
t294 = cos(qJ(6));
t277 = qJD(2) * t287 + qJD(3);
t230 = t277 * t296 - t291 * t356;
t231 = t277 * t291 + t296 * t356;
t321 = t230 * t290 + t295 * t231;
t146 = t265 * t294 - t289 * t321;
t350 = qJD(2) * t372;
t340 = t292 * t350;
t199 = t277 * t369 + (-t292 * t370 + t353) * t373;
t200 = -t277 * t370 + (-t292 * t369 - t354) * t373;
t347 = t295 * t230 - t231 * t290;
t95 = qJD(5) * t347 + t199 * t295 + t200 * t290;
t60 = qJD(6) * t146 + t289 * t340 + t294 * t95;
t147 = t265 * t289 + t294 * t321;
t61 = -qJD(6) * t147 - t289 * t95 + t294 * t340;
t21 = -mrSges(7,1) * t61 + mrSges(7,2) * t60;
t165 = pkin(10) * t277 + t182;
t276 = t287 * t374;
t190 = t276 + (qJD(2) * t338 - t264) * t285;
t121 = t165 * t296 + t190 * t291;
t101 = pkin(11) * t230 + t121;
t380 = t295 * t101;
t120 = -t165 * t291 + t296 * t190;
t100 = -pkin(11) * t231 + t120;
t87 = pkin(4) * t271 + t100;
t46 = t290 * t87 + t380;
t307 = t315 * qJD(2);
t362 = t288 * t389;
t341 = qJD(3) * t362;
t136 = (t264 * t386 - t250) * qJD(3) + (t286 * t307 + t341) * qJD(1);
t206 = (t313 + t357) * t373;
t57 = -qJD(4) * t121 - t136 * t291 + t296 * t206;
t47 = pkin(4) * t340 - pkin(11) * t199 + t57;
t56 = t296 * t136 - t165 * t370 + t190 * t369 + t291 * t206;
t53 = pkin(11) * t200 + t56;
t11 = -qJD(5) * t46 - t290 * t53 + t295 * t47;
t8 = -pkin(5) * t340 - t11;
t494 = -m(7) * t8 - t21;
t219 = t273 * t290 - t274 * t295;
t156 = qJD(5) * t219 + t262 * t290 - t295 * t346;
t72 = t125 * t295 - t133 * t290;
t493 = -t72 - t156;
t492 = t291 / 0.2e1;
t284 = -pkin(4) * t296 - pkin(3);
t205 = pkin(5) * t259 - pkin(12) * t260 + t284;
t149 = t205 * t294 - t219 * t289;
t491 = qJD(6) * t149 + t289 * t498 + t294 * t499;
t150 = t205 * t289 + t219 * t294;
t490 = -qJD(6) * t150 - t289 * t499 + t294 * t498;
t489 = pkin(5) * t356 - t493;
t39 = pkin(12) * t265 + t46;
t164 = -pkin(3) * t277 - t181;
t135 = -pkin(4) * t230 + t164;
t71 = -pkin(5) * t347 - pkin(12) * t321 + t135;
t19 = -t289 * t39 + t294 * t71;
t308 = t317 * qJD(2);
t342 = t288 * t352;
t137 = t479 * qJD(3) + (t286 * t308 + t342) * qJD(1);
t114 = -pkin(4) * t200 + t137;
t96 = qJD(5) * t321 + t199 * t290 - t295 * t200;
t32 = pkin(5) * t96 - pkin(12) * t95 + t114;
t10 = -t101 * t368 + t290 * t47 + t295 * t53 + t87 * t367;
t7 = pkin(12) * t340 + t10;
t2 = qJD(6) * t19 + t289 * t32 + t294 * t7;
t488 = t2 * mrSges(7,2);
t20 = t289 * t71 + t294 * t39;
t3 = -qJD(6) * t20 - t289 * t7 + t294 * t32;
t487 = t3 * mrSges(7,1);
t486 = -t277 * Ifges(4,6) / 0.2e1;
t238 = t278 + (-pkin(2) * t297 - pkin(3)) * t287;
t198 = -pkin(4) * t253 + t238;
t108 = -pkin(5) * t320 - pkin(12) * t195 + t198;
t482 = t290 * t143 + t295 * t154;
t79 = -pkin(12) * t389 + t482;
t43 = t108 * t289 + t294 * t79;
t485 = -qJD(6) * t43 + t289 * t496 + t495 * t294;
t42 = t108 * t294 - t289 * t79;
t484 = qJD(6) * t42 + t495 * t289 - t294 * t496;
t481 = -t188 + t124;
t480 = -t189 + t123;
t170 = Ifges(6,4) * t347;
t171 = qJD(6) - t347;
t329 = Ifges(7,5) * t294 - Ifges(7,6) * t289;
t309 = t171 * t329;
t417 = Ifges(7,4) * t289;
t333 = Ifges(7,1) * t294 - t417;
t310 = t147 * t333;
t416 = Ifges(7,4) * t294;
t331 = -Ifges(7,2) * t289 + t416;
t311 = t146 * t331;
t334 = mrSges(7,1) * t289 + mrSges(7,2) * t294;
t393 = t101 * t290;
t45 = t295 * t87 - t393;
t38 = -pkin(5) * t265 - t45;
t314 = t38 * t334;
t415 = Ifges(6,5) * t265;
t437 = -t294 / 0.2e1;
t438 = t289 / 0.2e1;
t423 = Ifges(6,1) * t321;
t111 = t170 + t415 + t423;
t470 = t135 * mrSges(6,2) + t111 / 0.2e1 - t45 * mrSges(6,3);
t418 = Ifges(7,4) * t147;
t67 = Ifges(7,2) * t146 + Ifges(7,6) * t171 + t418;
t144 = Ifges(7,4) * t146;
t68 = Ifges(7,1) * t147 + Ifges(7,5) * t171 + t144;
t478 = -t415 / 0.2e1 - t311 / 0.2e1 - t310 / 0.2e1 - t309 / 0.2e1 + t67 * t438 + t68 * t437 - t314 - t470 - t170 / 0.2e1;
t477 = t287 * t379 - t384;
t476 = -t291 * t57 + t296 * t56;
t475 = -t19 * t289 + t20 * t294;
t473 = -t11 * mrSges(6,1) + t10 * mrSges(6,2) - Ifges(6,5) * t95 + Ifges(6,6) * t96;
t402 = t231 * Ifges(5,4);
t158 = t230 * Ifges(5,2) + t271 * Ifges(5,6) + t402;
t226 = Ifges(5,4) * t230;
t159 = t231 * Ifges(5,1) + t271 * Ifges(5,5) + t226;
t324 = t120 * t296 + t121 * t291;
t420 = Ifges(5,4) * t296;
t421 = Ifges(5,4) * t291;
t435 = t296 / 0.2e1;
t440 = t271 / 0.2e1;
t444 = t231 / 0.2e1;
t446 = t230 / 0.2e1;
t472 = t324 * mrSges(5,3) - t164 * (mrSges(5,1) * t291 + mrSges(5,2) * t296) - (Ifges(5,5) * t296 - Ifges(5,6) * t291) * t440 - (-Ifges(5,2) * t291 + t420) * t446 - (Ifges(5,1) * t296 - t421) * t444 - t159 * t435 + t158 * t492;
t25 = -qJD(5) * t482 - t105 * t290 + t295 * t97;
t471 = -t57 * mrSges(5,1) + t56 * mrSges(5,2) - Ifges(5,5) * t199 - Ifges(5,6) * t200;
t119 = pkin(5) * t321 - pkin(12) * t347;
t412 = Ifges(6,6) * t265;
t413 = Ifges(6,2) * t347;
t419 = Ifges(6,4) * t321;
t110 = t412 + t413 + t419;
t410 = Ifges(7,3) * t171;
t411 = Ifges(7,6) * t146;
t414 = Ifges(7,5) * t147;
t66 = t410 + t411 + t414;
t469 = t20 * mrSges(7,2) + t46 * mrSges(6,3) + t110 / 0.2e1 - t66 / 0.2e1 - t135 * mrSges(6,1) - t19 * mrSges(7,1);
t467 = Ifges(6,2) / 0.2e1;
t58 = Ifges(7,6) * t61;
t59 = Ifges(7,5) * t60;
t16 = Ifges(7,3) * t96 + t58 + t59;
t466 = t16 / 0.2e1;
t465 = t60 / 0.2e1;
t464 = t61 / 0.2e1;
t463 = t68 / 0.2e1;
t462 = t96 / 0.2e1;
t460 = -t146 / 0.2e1;
t459 = t146 / 0.2e1;
t458 = -t147 / 0.2e1;
t457 = t147 / 0.2e1;
t455 = -t171 / 0.2e1;
t454 = t171 / 0.2e1;
t453 = t347 / 0.2e1;
t452 = t321 / 0.2e1;
t451 = t320 / 0.2e1;
t450 = t195 / 0.2e1;
t449 = t199 / 0.2e1;
t448 = t200 / 0.2e1;
t445 = -t231 / 0.2e1;
t443 = t253 / 0.2e1;
t442 = t254 / 0.2e1;
t441 = t265 / 0.2e1;
t439 = -t289 / 0.2e1;
t436 = t294 / 0.2e1;
t433 = t2 * t294;
t432 = t3 * t289;
t429 = t95 * Ifges(6,1);
t428 = t95 * Ifges(6,4);
t427 = t96 * Ifges(6,4);
t80 = mrSges(6,1) * t340 - mrSges(6,3) * t95;
t424 = t21 - t80;
t422 = Ifges(4,4) * t292;
t409 = t136 * mrSges(4,2);
t405 = t19 * t294;
t153 = mrSges(6,1) * t265 - mrSges(6,3) * t321;
t86 = -mrSges(7,1) * t146 + mrSges(7,2) * t147;
t394 = t153 - t86;
t208 = -t286 * t477 - t362;
t392 = t137 * t208;
t388 = t286 * t293;
t385 = t289 * t210;
t381 = t294 * t210;
t376 = mrSges(4,1) * t277 + mrSges(5,1) * t230 - mrSges(5,2) * t231 - mrSges(4,3) * t356;
t366 = qJD(6) * t289;
t365 = qJD(6) * t294;
t118 = -mrSges(6,1) * t347 + mrSges(6,2) * t321;
t360 = -t118 + t376;
t191 = t217 * t289 + t294 * t356;
t349 = t191 - t385;
t192 = -t217 * t294 + t289 * t356;
t348 = -t192 - t381;
t344 = t373 * t388;
t339 = t487 - t488;
t336 = -mrSges(4,1) * t297 + mrSges(4,2) * t292;
t335 = mrSges(7,1) * t294 - mrSges(7,2) * t289;
t332 = Ifges(7,1) * t289 + t416;
t330 = Ifges(7,2) * t294 + t417;
t328 = Ifges(7,5) * t289 + Ifges(7,6) * t294;
t327 = t20 * t289 + t405;
t316 = t287 * t383 + t382;
t209 = t286 * t316 + t288 * t390;
t252 = -t285 * t286 * t298 + t287 * t288;
t166 = -t209 * t291 + t252 * t296;
t167 = t209 * t296 + t252 * t291;
t113 = t166 * t290 + t167 * t295;
t77 = t113 * t294 + t208 * t289;
t76 = -t113 * t289 + t208 * t294;
t84 = t143 * t295 - t154 * t290;
t322 = t295 * t166 - t167 * t290;
t168 = -t195 * t289 - t294 * t389;
t318 = -t195 * t294 + t289 * t389;
t225 = -t264 * t285 + t276;
t275 = Ifges(4,4) * t355;
t306 = -t181 * mrSges(4,3) + Ifges(4,1) * t356 / 0.2e1 + t275 / 0.2e1 + t277 * Ifges(4,5) + t225 * mrSges(4,2);
t17 = t60 * Ifges(7,4) + t61 * Ifges(7,2) + t96 * Ifges(7,6);
t18 = t60 * Ifges(7,1) + t61 * Ifges(7,4) + t96 * Ifges(7,5);
t268 = Ifges(6,3) * t340;
t305 = mrSges(7,3) * t433 + qJD(6) * t314 + t17 * t436 + t18 * t438 + t330 * t464 + t332 * t465 - t8 * t335 + t268 - t67 * t366 / 0.2e1 + t365 * t463 + t328 * t462 - t473 + (t311 + t310 + t309) * qJD(6) / 0.2e1;
t103 = -mrSges(7,2) * t171 + mrSges(7,3) * t146;
t104 = mrSges(7,1) * t171 - mrSges(7,3) * t147;
t33 = mrSges(7,1) * t96 - mrSges(7,3) * t60;
t34 = -mrSges(7,2) * t96 + mrSges(7,3) * t61;
t304 = t294 * t34 + m(7) * (-t19 * t365 - t20 * t366 - t432 + t433) - t289 * t33 - t104 * t365 - t103 * t366;
t303 = -t410 / 0.2e1 - t411 / 0.2e1 - t414 / 0.2e1 + t412 / 0.2e1 + t419 / 0.2e1 + t469;
t301 = t120 * mrSges(5,1) + t225 * mrSges(4,1) + t45 * mrSges(6,1) + t265 * Ifges(6,3) + t321 * Ifges(6,5) + t347 * Ifges(6,6) + t271 * Ifges(5,3) + t231 * Ifges(5,5) + t230 * Ifges(5,6) + t486 - (Ifges(4,2) * t297 + t422) * t373 / 0.2e1 - t121 * mrSges(5,2) - t182 * mrSges(4,3) - t46 * mrSges(6,2);
t270 = Ifges(4,5) * t297 * t350;
t269 = Ifges(5,3) * t340;
t244 = t336 * t373;
t243 = -mrSges(4,2) * t277 + mrSges(4,3) * t355;
t237 = (mrSges(4,1) * t292 + mrSges(4,2) * t297) * t350;
t204 = mrSges(5,1) * t271 - mrSges(5,3) * t231;
t203 = -mrSges(5,2) * t271 + mrSges(5,3) * t230;
t178 = -mrSges(5,2) * t340 + mrSges(5,3) * t200;
t177 = mrSges(5,1) * t340 - mrSges(5,3) * t199;
t161 = t342 + (qJD(3) * t316 + t308) * t286;
t160 = t341 + (qJD(3) * t477 + t307) * t286;
t152 = -mrSges(6,2) * t265 + mrSges(6,3) * t347;
t138 = -mrSges(5,1) * t200 + mrSges(5,2) * t199;
t129 = t199 * Ifges(5,1) + t200 * Ifges(5,4) + Ifges(5,5) * t340;
t128 = t199 * Ifges(5,4) + t200 * Ifges(5,2) + Ifges(5,6) * t340;
t126 = -t295 * t188 + t189 * t290;
t99 = pkin(4) * t231 + t119;
t94 = -qJD(4) * t167 - t160 * t291 + t296 * t344;
t93 = qJD(4) * t166 + t160 * t296 + t291 * t344;
t81 = -mrSges(6,2) * t340 - mrSges(6,3) * t96;
t78 = pkin(5) * t389 - t84;
t75 = qJD(6) * t318 - t116 * t289 + t294 * t352;
t74 = qJD(6) * t168 + t116 * t294 + t289 * t352;
t52 = t100 * t295 - t393;
t51 = t100 * t290 + t380;
t48 = mrSges(6,1) * t96 + mrSges(6,2) * t95;
t41 = Ifges(6,5) * t340 - t427 + t429;
t40 = -t96 * Ifges(6,2) + Ifges(6,6) * t340 + t428;
t31 = t119 * t289 + t294 * t45;
t30 = t119 * t294 - t289 * t45;
t29 = qJD(5) * t113 + t290 * t93 - t295 * t94;
t28 = qJD(5) * t322 + t290 * t94 + t295 * t93;
t27 = t289 * t99 + t294 * t52;
t26 = -t289 * t52 + t294 * t99;
t23 = -pkin(5) * t352 - t25;
t13 = -qJD(6) * t77 + t161 * t294 - t28 * t289;
t12 = qJD(6) * t76 + t161 * t289 + t28 * t294;
t1 = [t12 * t103 + t13 * t104 + t113 * t81 + t28 * t152 + t160 * t243 + t166 * t177 + t167 * t178 + t93 * t203 + t94 * t204 + t252 * t237 + t76 * t33 + t77 * t34 - t394 * t29 + (-mrSges(3,1) * t293 - mrSges(3,2) * t298) * qJD(2) ^ 2 * t286 + (t138 + t48) * t208 - t424 * t322 - t360 * t161 + (t244 * t388 + (t208 * t297 - t209 * t292) * qJD(3) * mrSges(4,3)) * t373 + m(4) * (t136 * t209 + t392 + t160 * t182 - t161 * t181 + (qJD(1) * t252 + t225) * t344) + m(5) * (t120 * t94 + t121 * t93 + t161 * t164 + t166 * t57 + t167 * t56 + t392) + m(6) * (t10 * t113 + t11 * t322 + t114 * t208 + t135 * t161 + t28 * t46 - t29 * t45) + m(7) * (t12 * t20 + t13 * t19 + t2 * t77 + t29 * t38 + t3 * t76 - t322 * t8); (Ifges(7,5) * t74 + Ifges(7,6) * t75) * t454 + (Ifges(7,4) * t74 + Ifges(7,2) * t75) * t459 + (Ifges(7,1) * t74 + Ifges(7,4) * t75) * t457 + t482 * t81 - t318 * t18 / 0.2e1 + (t168 * t2 - t19 * t74 + t20 * t75 + t3 * t318) * mrSges(7,3) + t8 * (-mrSges(7,1) * t168 - mrSges(7,2) * t318) + (-Ifges(7,5) * t318 + Ifges(7,6) * t168 - Ifges(7,3) * t320) * t462 + (-Ifges(7,4) * t318 + Ifges(7,2) * t168 - Ifges(7,6) * t320) * t464 + (-Ifges(7,1) * t318 + Ifges(7,4) * t168 - Ifges(7,5) * t320) * t465 - t320 * t487 + (t10 * t320 - t11 * t195) * mrSges(6,3) - t96 * (Ifges(6,4) * t195 + Ifges(6,2) * t320) / 0.2e1 + t95 * (Ifges(6,1) * t195 + Ifges(6,4) * t320) / 0.2e1 + t114 * (-mrSges(6,1) * t320 + mrSges(6,2) * t195) - t320 * t466 + t320 * t488 + m(4) * (t136 * t257 - t137 * t256 - t181 * t248 + t182 * t247) + t137 * (-mrSges(5,1) * t253 + mrSges(5,2) * t254) + t238 * t138 + t212 * t159 / 0.2e1 + t164 * (-mrSges(5,1) * t213 + mrSges(5,2) * t212) + t213 * t158 / 0.2e1 + t198 * t48 + t183 * t177 + t184 * t178 + t179 * t118 + t168 * t17 / 0.2e1 + t25 * t153 + ((-m(4) * pkin(2) + t336) * t345 + ((Ifges(4,5) * t287 / 0.2e1 - t256 * mrSges(4,3) + 0.3e1 / 0.2e1 * Ifges(4,4) * t389) * t297 + (-Ifges(4,6) * t287 + Ifges(5,5) * t442 + Ifges(5,6) * t443 + Ifges(6,5) * t450 + Ifges(6,6) * t451 - t257 * mrSges(4,3) - 0.3e1 / 0.2e1 * Ifges(4,4) * t390 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) - Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1) * t389) * t292) * qJD(3)) * t373 + t84 * t80 + t23 * t86 + t78 * t21 + t38 * (-mrSges(7,1) * t75 + mrSges(7,2) * t74) + t75 * t67 / 0.2e1 + (t306 * t297 + (t301 + t486) * t292) * t372 + (t270 / 0.2e1 - t137 * mrSges(4,1) - t409) * t287 + t42 * t33 + t43 * t34 + t484 * t103 + t485 * t104 + (t2 * t43 + t3 * t42 + t78 * t8 + (-t126 + t23) * t38 + t484 * t20 + t485 * t19) * m(7) + t394 * t126 + t483 * t152 - t376 * t248 + t480 * t203 + t481 * t204 + (t137 * t238 + t183 * t57 + t184 * t56 + (-t222 + t248) * t164 + t480 * t121 + t481 * t120) * m(5) + (t10 * t482 + t11 * t84 + t114 * t198 + t483 * t46 + (t126 + t25) * t45 + t497 * t135) * m(6) + (-t244 * t357 + t137 * mrSges(4,3) * t292 - pkin(2) * t237 + (-t269 / 0.2e1 - t268 / 0.2e1 + t136 * mrSges(4,3) + t471 + t473) * t297) * t285 + t360 * t222 + (Ifges(6,1) * t452 + Ifges(6,4) * t453 + Ifges(6,5) * t441 + t470) * t116 + (-Ifges(6,4) * t452 + Ifges(7,5) * t457 - Ifges(6,2) * t453 - Ifges(6,6) * t441 + Ifges(7,6) * t459 + Ifges(7,3) * t454 - t469) * t117 + (-t120 * t212 + t121 * t213 + t253 * t56 - t254 * t57) * mrSges(5,3) - m(4) * (-t181 * t222 + t182 * t223 + t225 * t345) + (t247 - t223) * t243 + (Ifges(5,5) * t212 + Ifges(5,6) * t213) * t440 + t129 * t442 + t128 * t443 + (Ifges(5,1) * t212 + Ifges(5,4) * t213) * t444 + (Ifges(5,4) * t212 + Ifges(5,2) * t213) * t446 + (Ifges(5,4) * t254 + Ifges(5,2) * t253) * t448 + (Ifges(5,1) * t254 + Ifges(5,4) * t253) * t449 + t41 * t450 + t40 * t451 + t74 * t463; (-t66 + t110) * (t216 / 0.2e1 - t211 / 0.2e1) - t424 * t319 + m(6) * (t10 * t219 + t11 * t319 + t114 * t284 + t155 * t46 - t156 * t45) + (t217 / 0.2e1 - t210 / 0.2e1) * t111 + (-Ifges(6,5) * t210 - Ifges(6,6) * t211) * t441 + (-Ifges(6,1) * t210 - Ifges(6,4) * t211) * t452 + (-Ifges(6,4) * t210 - Ifges(6,2) * t211) * t453 + (-t192 / 0.2e1 - t381 / 0.2e1) * t68 + (-Ifges(7,5) * t381 + Ifges(7,6) * t385 + Ifges(7,3) * t211) * t454 + (-Ifges(7,1) * t381 + Ifges(7,4) * t385 + Ifges(7,5) * t211) * t457 + (-Ifges(7,4) * t381 + Ifges(7,2) * t385 + Ifges(7,6) * t211) * t459 + (-t191 / 0.2e1 + t385 / 0.2e1) * t67 + (-t177 * t291 + t178 * t296) * pkin(10) - t409 - m(5) * (t120 * t141 + t121 * t142 + t164 * t182) - m(6) * (t135 * t163 + t45 * t72 + t46 * t73) + t270 + (mrSges(6,1) * t377 - mrSges(6,2) * t501) * t135 + (-t377 * t46 + t45 * t501) * mrSges(6,3) + (-mrSges(5,1) * t296 + mrSges(5,2) * t291 - mrSges(4,1)) * t137 + t493 * t153 + t284 * t48 + t129 * t492 - t181 * t243 + t219 * t81 - t142 * t203 - t141 * t204 - t163 * t118 + t149 * t33 + t150 * t34 - pkin(3) * t138 + (t58 / 0.2e1 + t466 - t40 / 0.2e1 + t114 * mrSges(6,1) + t59 / 0.2e1 - t10 * mrSges(6,3) - t428 / 0.2e1 + (Ifges(7,3) / 0.2e1 + t467) * t96 + t339) * t259 + (mrSges(7,1) * t377 - mrSges(7,3) * t348) * t19 + (-mrSges(7,2) * t377 - mrSges(7,3) * t349) * t20 + t376 * t182 + (t329 * t462 + t333 * t465 + t331 * t464 + t8 * t334 + t114 * mrSges(6,2) + t429 / 0.2e1 - t427 / 0.2e1 + t41 / 0.2e1 + t18 * t436 + t17 * t439 - t11 * mrSges(6,3) + (-t2 * t289 - t294 * t3) * mrSges(7,3) + (-mrSges(7,3) * t475 + t328 * t455 + t330 * t460 + t332 * t458 + t38 * t335 + t67 * t437 + t68 * t439) * qJD(6)) * t260 + m(5) * (-pkin(3) * t137 + pkin(10) * t476) + t476 * mrSges(5,3) + ((m(6) * t135 + t118) * t291 * pkin(4) + (-m(5) * t324 - t291 * t203 - t296 * t204) * pkin(10) - t472) * qJD(4) + ((-t275 / 0.2e1 - t306 + t472) * t297 + (-t301 + (-qJD(3) + t277 / 0.2e1) * Ifges(4,6) + (t422 / 0.2e1 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t297) * t373 + (Ifges(5,5) * t291 + Ifges(6,5) * t260 + Ifges(5,6) * t296 - Ifges(6,6) * t259) * qJD(3) / 0.2e1) * t292) * t373 + (mrSges(7,1) * t349 + mrSges(7,2) * t348) * t38 + t500 * t152 - t265 * (-Ifges(6,5) * t217 - Ifges(6,6) * t216) / 0.2e1 - t347 * (-Ifges(6,4) * t217 - Ifges(6,2) * t216) / 0.2e1 - t321 * (-Ifges(6,1) * t217 - Ifges(6,4) * t216) / 0.2e1 + t128 * t435 + (Ifges(5,2) * t296 + t421) * t448 + (Ifges(5,1) * t291 + t420) * t449 + t489 * t86 + t490 * t104 + t491 * t103 + (t149 * t3 + t150 * t2 + t19 * t490 + t20 * t491 - t319 * t8 + t38 * t489) * m(7) + (Ifges(7,5) * t192 + Ifges(7,6) * t191 + Ifges(7,3) * t216) * t455 + (Ifges(7,1) * t192 + Ifges(7,4) * t191 + Ifges(7,5) * t216) * t458 + (Ifges(7,4) * t192 + Ifges(7,2) * t191 + Ifges(7,6) * t216) * t460; (-t171 * t405 + (-t171 * t20 - t3) * t289) * mrSges(7,3) - m(6) * (-t45 * t51 + t46 * t52) + t304 * (pkin(4) * t290 + pkin(12)) + t269 + t305 - t471 + (t120 * t230 + t121 * t231) * mrSges(5,3) - t271 * (Ifges(5,5) * t230 - Ifges(5,6) * t231) / 0.2e1 - t164 * (mrSges(5,1) * t231 + mrSges(5,2) * t230) + (-t231 * t118 + t290 * t81 + t295 * t80 + ((m(7) * t38 - t394) * t290 + (m(7) * t475 + t103 * t294 - t104 * t289 + t152) * t295) * qJD(5) + (t10 * t290 + t11 * t295 + 0.2e1 * t135 * t445 + (-t290 * t45 + t295 * t46) * qJD(5)) * m(6)) * pkin(4) - t120 * t203 + t121 * t204 - t52 * t152 - t27 * t103 - t26 * t104 + t394 * t51 + (-t423 / 0.2e1 + t478) * t347 + (t413 / 0.2e1 + t303) * t321 + t158 * t444 + (Ifges(5,1) * t230 - t402) * t445 - m(7) * (t19 * t26 + t20 * t27 + t38 * t51) - (-Ifges(5,2) * t231 + t159 + t226) * t230 / 0.2e1 - t494 * (-pkin(4) * t295 - pkin(5)); (-qJD(6) * t327 - t432) * mrSges(7,3) + t394 * t46 + ((-Ifges(6,1) / 0.2e1 + t467) * t321 + t327 * mrSges(7,3) + t478) * t347 - m(7) * (t19 * t30 + t20 * t31 + t38 * t46) + t303 * t321 + t304 * pkin(12) + t305 - t45 * t152 - t31 * t103 - t30 * t104 + t494 * pkin(5); t67 * t457 + (Ifges(7,5) * t146 - Ifges(7,6) * t147) * t455 - t19 * t103 + t20 * t104 - t38 * (mrSges(7,1) * t147 + mrSges(7,2) * t146) + (Ifges(7,1) * t146 - t418) * t458 + (t146 * t19 + t147 * t20) * mrSges(7,3) + t339 + t16 + (-Ifges(7,2) * t147 + t144 + t68) * t460;];
tauc  = t1(:);
