% Calculate vector of inverse dynamics joint torques for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR8_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:05
% EndTime: 2019-12-31 16:55:14
% DurationCPUTime: 7.13s
% Computational Cost: add. (6074->547), mult. (9611->716), div. (0->0), fcn. (7493->6), ass. (0->305)
t248 = qJ(3) + qJ(4);
t231 = sin(t248);
t232 = cos(t248);
t250 = sin(qJ(1));
t375 = t232 * t250;
t199 = Icges(5,4) * t375;
t376 = t231 * t250;
t252 = cos(qJ(1));
t380 = Icges(5,5) * t252;
t114 = Icges(5,1) * t376 + t199 + t380;
t247 = qJD(3) + qJD(4);
t182 = t247 * t250;
t183 = t247 * t252;
t383 = Icges(5,4) * t231;
t165 = Icges(5,1) * t232 - t383;
t287 = Icges(5,2) * t232 + t383;
t449 = -t287 + t165;
t382 = Icges(5,4) * t232;
t289 = Icges(5,1) * t231 + t382;
t115 = -Icges(5,5) * t250 + t252 * t289;
t163 = -Icges(5,2) * t231 + t382;
t451 = t163 * t252 + t115;
t423 = qJD(1) * t449 - t182 * t451 + t183 * (-Icges(5,2) * t376 + t114 + t199);
t450 = t163 + t289;
t113 = -Icges(5,6) * t250 + t252 * t287;
t452 = -t165 * t252 + t113;
t112 = Icges(5,6) * t252 + t250 * t287;
t453 = -t165 * t250 + t112;
t424 = qJD(1) * t450 - t182 * t452 + t183 * t453;
t455 = t231 * t424 - t232 * t423;
t249 = sin(qJ(3));
t454 = rSges(4,2) * t249;
t297 = rSges(5,1) * t231 + rSges(5,2) * t232;
t251 = cos(qJ(3));
t371 = t250 * t251;
t214 = Icges(4,4) * t371;
t374 = t249 * t250;
t381 = Icges(4,5) * t252;
t132 = Icges(4,1) * t374 + t214 + t381;
t384 = Icges(4,4) * t251;
t290 = Icges(4,1) * t249 + t384;
t133 = -Icges(4,5) * t250 + t252 * t290;
t187 = -Icges(4,2) * t249 + t384;
t154 = t187 * t252;
t261 = t250 * (t133 + t154) - t252 * (-Icges(4,2) * t374 + t132 + t214);
t385 = Icges(4,4) * t249;
t288 = Icges(4,2) * t251 + t385;
t130 = Icges(4,6) * t252 + t250 * t288;
t131 = -Icges(4,6) * t250 + t252 * t288;
t189 = Icges(4,1) * t251 - t385;
t155 = t189 * t250;
t156 = t189 * t252;
t262 = t250 * (t131 - t156) - t252 * (t130 - t155);
t448 = -t262 * t249 + t261 * t251;
t355 = t187 + t290;
t356 = -t288 + t189;
t447 = (t249 * t355 - t251 * t356) * qJD(1);
t196 = -rSges(3,2) * t252 + t250 * rSges(3,3);
t195 = t252 * pkin(1) + t250 * qJ(2);
t335 = qJD(1) * qJD(2);
t341 = qJD(1) * t250;
t229 = qJD(2) * t250;
t340 = qJD(1) * t252;
t350 = qJ(2) * t340 + t229;
t329 = qJD(1) * (-pkin(1) * t341 + t350) + qJDD(1) * t195 + t250 * t335;
t333 = qJDD(2) * t252;
t349 = rSges(3,2) * t341 + rSges(3,3) * t340;
t446 = qJD(1) * t349 + qJDD(1) * t196 - g(2) + t329 - t333;
t445 = g(1) * t250;
t444 = rSges(4,2) * t251;
t202 = rSges(5,1) * t375;
t330 = rSges(5,2) * t376;
t143 = t202 - t330;
t406 = pkin(5) * t252;
t311 = t195 + t406;
t336 = qJD(3) * t252;
t320 = t251 * t336;
t209 = pkin(3) * t320;
t230 = qJD(2) * t252;
t352 = -t209 - t230;
t244 = t252 * rSges(5,3);
t116 = rSges(5,1) * t376 + rSges(5,2) * t375 + t244;
t217 = pkin(3) * t374;
t253 = -pkin(6) - pkin(5);
t401 = pkin(5) + t253;
t148 = -t401 * t252 + t217;
t364 = t116 + t148;
t396 = rSges(5,2) * t231;
t398 = rSges(5,1) * t232;
t167 = -t396 + t398;
t378 = t167 * t183;
t43 = -t378 + (t311 + t364) * qJD(1) + t352;
t443 = t43 * (qJD(1) * t143 + t183 * t297);
t285 = Icges(5,5) * t231 + Icges(5,6) * t232;
t111 = -Icges(5,3) * t250 + t252 * t285;
t45 = -t252 * t111 - t113 * t375 - t115 * t376;
t278 = t163 * t232 + t165 * t231;
t161 = Icges(5,5) * t232 - Icges(5,6) * t231;
t379 = t161 * t252;
t57 = t250 * t278 + t379;
t442 = t57 * qJD(1) + t182 * t45;
t203 = t252 * t396;
t144 = -t252 * t398 + t203;
t70 = -t378 + (t250 * t297 + t244) * qJD(1);
t441 = t143 * t182 - t183 * t144 + t252 * t70;
t283 = t113 * t232 + t115 * t231;
t440 = t252 * t283;
t280 = t131 * t251 + t133 * t249;
t437 = t280 * t252;
t245 = t252 * rSges(4,3);
t134 = rSges(4,1) * t374 + rSges(4,2) * t371 + t245;
t436 = t134 + t311;
t149 = qJD(1) * t195 - t230;
t298 = rSges(4,1) * t249 + t444;
t173 = t298 * qJD(3);
t334 = qJD(1) * qJD(3);
t176 = qJDD(3) * t250 + t252 * t334;
t194 = rSges(4,1) * t251 - t454;
t351 = qJDD(2) * t250 + t252 * t335;
t405 = pkin(5) * qJD(1) ^ 2;
t275 = -t252 * t405 + t351;
t234 = t252 * qJ(2);
t191 = pkin(1) * t250 - t234;
t135 = -t250 * rSges(4,3) + t252 * t298;
t407 = pkin(5) * t250;
t312 = t135 - t407;
t300 = -t191 + t312;
t338 = qJD(3) * t250;
t158 = t194 * t252;
t86 = -qJD(3) * t158 + (t250 * t298 + t245) * qJD(1);
t32 = -t173 * t338 + t176 * t194 + (-t149 - t86) * qJD(1) + t300 * qJDD(1) + t275;
t227 = qJDD(3) * t252;
t177 = -t250 * t334 + t227;
t271 = qJDD(1) * t406 - t250 * t405 + t329;
t337 = qJD(3) * t251;
t321 = t250 * t337;
t325 = t249 * t340;
t327 = t340 * t444 + (t321 + t325) * rSges(4,1);
t339 = qJD(3) * t249;
t87 = (-rSges(4,2) * t339 - rSges(4,3) * qJD(1)) * t250 + t327;
t33 = qJD(1) * t87 + qJDD(1) * t134 - t177 * t194 + (qJD(3) * t173 - qJDD(2)) * t252 + t271;
t435 = t32 * t250 - t33 * t252;
t121 = qJD(4) * t340 + qJDD(4) * t250 + t176;
t127 = t297 * t247;
t241 = t250 * rSges(5,3);
t117 = t252 * t297 - t241;
t409 = pkin(3) * t249;
t147 = t250 * t401 + t252 * t409;
t363 = t117 + t147;
t301 = t363 - t407;
t274 = -t191 + t301;
t373 = t249 * qJD(3) ^ 2;
t97 = qJD(1) * t148 - t209;
t14 = t121 * t167 - t127 * t182 + (t176 * t251 - t250 * t373) * pkin(3) + (-t149 - t70 - t97) * qJD(1) + t274 * qJDD(1) + t275;
t122 = qJDD(4) * t252 - t247 * t341 + t227;
t328 = t247 * t202 + t297 * t340;
t71 = (-rSges(5,3) * qJD(1) - t247 * t396) * t250 + t328;
t207 = pkin(3) * t321;
t326 = pkin(3) * t325 + t253 * t341 + t207;
t96 = pkin(5) * t341 + t326;
t400 = t71 + t96;
t15 = -t333 - t122 * t167 + t127 * t183 + t364 * qJDD(1) + t400 * qJD(1) + (-t177 * t251 + t252 * t373) * pkin(3) + t271;
t434 = t14 * t250 - t15 * t252;
t433 = qJD(1) * t167;
t286 = Icges(4,5) * t249 + Icges(4,6) * t251;
t129 = -Icges(4,3) * t250 + t252 * t286;
t344 = qJD(1) * t129;
t75 = t131 * t249 - t133 * t251;
t82 = qJD(1) * t130 - qJD(3) * t154;
t84 = -qJD(3) * t156 + (t250 * t290 + t381) * qJD(1);
t432 = qJD(3) * t75 + t249 * t84 + t251 * t82 + t344;
t170 = t288 * qJD(3);
t171 = t290 * qJD(3);
t185 = Icges(4,5) * t251 - Icges(4,6) * t249;
t276 = t187 * t249 - t189 * t251;
t431 = qJD(1) * t185 + qJD(3) * t276 + t170 * t251 + t171 * t249;
t281 = t130 * t249 - t132 * t251;
t128 = Icges(4,3) * t252 + t250 * t286;
t345 = qJD(1) * t128;
t83 = qJD(1) * t131 + t187 * t338;
t85 = qJD(1) * t133 + qJD(3) * t155;
t430 = qJD(3) * t281 - t249 * t85 - t251 * t83 + t345;
t429 = g(2) * t252 - t445;
t304 = t450 * t247;
t305 = t449 * t247;
t427 = qJD(1) * t161 + t231 * t304 - t232 * t305;
t307 = -qJD(1) * t112 + t247 * t451;
t309 = (t250 * t289 + t380) * qJD(1) + t452 * t247;
t346 = qJD(1) * t111;
t426 = t231 * t309 - t232 * t307 + t346;
t308 = qJD(1) * t113 + t114 * t247 + t163 * t182;
t310 = -qJD(1) * t115 + t247 * t453;
t110 = Icges(5,3) * t252 + t250 * t285;
t347 = qJD(1) * t110;
t425 = t231 * t310 - t232 * t308 + t347;
t422 = t250 ^ 2;
t421 = -pkin(1) - pkin(5);
t420 = t121 / 0.2e1;
t419 = t122 / 0.2e1;
t418 = t176 / 0.2e1;
t417 = t177 / 0.2e1;
t416 = -t182 / 0.2e1;
t415 = t182 / 0.2e1;
t414 = -t183 / 0.2e1;
t413 = t183 / 0.2e1;
t412 = t250 / 0.2e1;
t411 = t252 / 0.2e1;
t410 = rSges(3,2) - pkin(1);
t408 = pkin(3) * t251;
t404 = -qJD(1) / 0.2e1;
t403 = qJD(1) / 0.2e1;
t402 = -pkin(1) + t253;
t394 = rSges(3,3) * t252;
t299 = t182 * t167 + t207 + t229;
t42 = qJD(1) * t274 + t299;
t390 = t252 * t42;
t159 = t194 * t338;
t61 = qJD(1) * t300 + t159 + t229;
t389 = t252 * t61;
t386 = qJDD(1) / 0.2e1;
t377 = t185 * t252;
t372 = t250 * t111;
t137 = t250 * t161;
t151 = t250 * t185;
t58 = t252 * t278 - t137;
t370 = t58 * qJD(1);
t277 = t251 * t187 + t249 * t189;
t77 = t252 * t277 - t151;
t369 = t77 * qJD(1);
t192 = rSges(3,2) * t250 + t394;
t354 = -t191 + t192;
t146 = t195 + t196;
t179 = qJD(1) * t191;
t348 = t229 - t179;
t342 = qJD(1) * t286;
t332 = -rSges(4,3) + t421;
t218 = pkin(3) * t371;
t44 = t252 * t110 + t112 * t375 + t114 * t376;
t48 = t252 * t128 + t130 * t371 + t132 * t374;
t49 = -t252 * t129 - t131 * t371 - t133 * t374;
t323 = t249 * t338;
t318 = -t341 / 0.2e1;
t317 = t340 / 0.2e1;
t316 = -t338 / 0.2e1;
t315 = t338 / 0.2e1;
t314 = -t336 / 0.2e1;
t313 = t336 / 0.2e1;
t302 = -qJD(1) * t144 - t182 * t297;
t197 = rSges(2,1) * t252 - rSges(2,2) * t250;
t193 = rSges(2,1) * t250 + rSges(2,2) * t252;
t282 = t130 * t251 + t132 * t249;
t263 = qJD(1) * t282 + qJD(3) * t151 + t344;
t264 = -qJD(1) * t280 - qJD(3) * t377 + t345;
t296 = (t263 * t250 + t252 * t430) * t252 + (t264 * t250 - t252 * t432) * t250;
t295 = (-t250 * t430 + t263 * t252) * t252 + (t250 * t432 + t264 * t252) * t250;
t294 = t250 * t49 + t252 * t48;
t118 = t250 * t128;
t50 = -t282 * t252 + t118;
t51 = -t250 * t129 + t437;
t293 = t250 * t51 + t252 * t50;
t62 = qJD(1) * t436 - t194 * t336 - t230;
t292 = t250 * t61 - t252 * t62;
t291 = -t250 * t87 + t252 * t86;
t284 = t112 * t232 + t114 * t231;
t60 = t113 * t231 - t115 * t232;
t279 = -t134 * t250 - t135 * t252;
t273 = t44 + t372;
t272 = t297 + t409;
t268 = -qJD(1) * t285 + t137 * t183 - t182 * t379;
t101 = t250 * t110;
t46 = -t252 * t284 + t101;
t266 = -qJD(1) * t283 - t247 * t379 + t347;
t265 = qJD(1) * t284 + t137 * t247 + t346;
t260 = qJD(1) * t278 - t285 * t247;
t259 = t277 * qJD(1) - t286 * qJD(3);
t10 = t266 * t250 - t252 * t426;
t11 = -t250 * t425 + t265 * t252;
t12 = t250 * t426 + t266 * t252;
t22 = t183 * t44 + t442;
t47 = -t372 + t440;
t23 = t182 * t47 + t183 * t46 - t370;
t28 = t260 * t250 + t252 * t427;
t29 = -t250 * t427 + t260 * t252;
t30 = -t231 * t308 - t232 * t310;
t31 = t231 * t307 + t232 * t309;
t59 = -t112 * t231 + t114 * t232;
t9 = t265 * t250 + t252 * t425;
t258 = (qJD(1) * t28 - qJDD(1) * t58 + t10 * t182 + t121 * t47 + t122 * t46 + t183 * t9) * t412 + (-t231 * t423 - t232 * t424) * t404 + t22 * t318 + t23 * t317 + (qJD(1) * t29 + qJDD(1) * t57 + t11 * t183 + t12 * t182 + t121 * t45 + t122 * t44) * t411 + (t250 * t47 + t252 * t46) * t420 + (t250 * t45 + t252 * t44) * t419 + (t10 * t250 + t252 * t9 + (-t250 * t46 + t252 * t47) * qJD(1)) * t415 + (t250 * t60 + t252 * t59) * t386 + (t11 * t252 + t12 * t250 + (-t250 * t44 + t252 * t45) * qJD(1)) * t413 + (t250 * t31 + t252 * t30 + (-t250 * t59 + t252 * t60) * qJD(1)) * t403 + (t268 * t250 + t455 * t252) * t416 + (-t455 * t250 + t268 * t252) * t414;
t157 = t194 * t250;
t104 = t252 * t117;
t100 = qJD(1) * t146 - t230;
t99 = qJD(1) * t354 + t229;
t76 = t250 * t277 + t377;
t73 = t76 * qJD(1);
t72 = t279 * qJD(3);
t52 = t354 * qJDD(1) + (-qJD(1) * t196 - t149) * qJD(1) + t351;
t39 = -t116 * t182 - t117 * t183 + (-t147 * t252 - t148 * t250) * qJD(3);
t37 = -t250 * t431 + t259 * t252;
t36 = t259 * t250 + t252 * t431;
t35 = t280 * qJD(3) - t249 * t82 + t251 * t84;
t34 = -qJD(3) * t282 - t249 * t83 + t251 * t85;
t27 = qJD(3) * t293 - t369;
t26 = qJD(3) * t294 + t73;
t8 = -t116 * t121 - t117 * t122 - t147 * t177 - t148 * t176 - t182 * t71 + t183 * t70 + (-t250 * t96 + t252 * t97) * qJD(3);
t1 = [((t273 + t47 - t440) * t183 + t442) * t416 + t75 * t418 + t60 * t420 + (t73 + ((-t50 + t118 + t49) * t250 + (t51 - t437 + (t129 - t282) * t250 + t48) * t252) * qJD(3)) * t316 - t121 * t58 / 0.2e1 - t176 * t77 / 0.2e1 - m(2) * (-g(1) * t193 + g(2) * t197) + (t59 + t57) * t419 + (-t281 + t76) * t417 + (t370 + (t283 * t250 - t101 + t45) * t183 + (t273 - t44) * t182 + ((t111 + t284) * t183 - t283 * t182) * t252 + t23) * t414 + (t30 + t29) * t413 + (t369 + (t422 * t129 + (-t118 + t49 + (t129 + t282) * t252) * t252) * qJD(3) + t27) * t314 + (t34 + t37) * t313 + (-qJD(3) * t277 + t170 * t249 - t171 * t251 - t231 * t305 - t232 * t304) * qJD(1) + (t31 + t28 + t22) * t415 + (t35 + t36 + t26) * t315 + (-(qJD(1) * t301 - t179 + t299 - t42) * t43 + t42 * (-t144 * t247 - t352) + t43 * (-t247 * t330 + t326 + t328 + t350) + ((-rSges(5,3) + t402) * t390 + (t42 * (-qJ(2) - t272) + t43 * (-rSges(5,3) - pkin(1))) * t250) * qJD(1) + (t15 - g(2)) * (-t252 * t253 + t116 + t195 + t217) + (t14 - g(1)) * (t250 * t402 + t252 * t272 + t234 - t241)) * m(5) + (t61 * (rSges(4,1) * t320 - t336 * t454 + t230) + t62 * (-rSges(4,2) * t323 + t327 + t350) + (t332 * t389 + (t61 * (-qJ(2) - t298) + t62 * t332) * t250) * qJD(1) - (qJD(1) * t312 + t159 + t348 - t61) * t62 + (t33 - g(2)) * t436 + (t32 - g(1)) * (t250 * t421 + t135 + t234)) * m(4) + (t99 * t230 + t100 * (t349 + t350) + (t99 * t410 * t252 + (t99 * (-rSges(3,3) - qJ(2)) - t100 * pkin(1)) * t250) * qJD(1) - (qJD(1) * t192 + t348 - t99) * t100 + t446 * t146 + (t52 - g(1)) * (t250 * t410 + t234 + t394)) * m(3) + (m(2) * (t193 ^ 2 + t197 ^ 2) - t276 - t163 * t231 + t165 * t232 + Icges(2,3) + Icges(3,1)) * qJDD(1); (-t252 * t446 + 0.2e1 * t52 * t412 - t445) * m(3) + (t429 + t434) * m(5) + (t429 + t435) * m(4); (t35 * t250 + t34 * t252 + (t250 * t281 + t75 * t252) * qJD(1)) * t403 + (qJD(1) * t37 + qJD(3) * t295 + qJDD(1) * t76 + t176 * t49 + t177 * t48) * t411 + (qJD(1) * t36 + qJD(3) * t296 - qJDD(1) * t77 + t176 * t51 + t177 * t50) * t412 + ((t151 * t336 - t342) * t252 + (-t447 + (-t252 * t377 - t448) * qJD(3)) * t250) * t314 + ((-t249 * t356 - t251 * t355) * qJD(1) + (t249 * t261 + t251 * t262) * qJD(3)) * t404 + ((-t338 * t377 - t342) * t250 + (t447 + (t250 * t151 + t448) * qJD(3)) * t252) * t316 + ((-t50 * t250 + t51 * t252) * qJD(1) + t296) * t315 + ((-t48 * t250 + t49 * t252) * qJD(1) + t295) * t313 + t258 + t26 * t318 + t27 * t317 + (t75 * t250 - t252 * t281) * t386 + t293 * t418 + t294 * t417 + (t14 * t218 - t8 * t104 + (t15 * (-t167 - t408) + t43 * t127 - t8 * t147 + t42 * t433) * t252 + (t14 * t167 + t42 * (-pkin(3) * t339 - t127) - t8 * t364 + t43 * t433) * t250 - t42 * (-pkin(3) * t323 + t302) - t443 - g(1) * (t143 + t218) - g(2) * (t203 + (-t398 - t408) * t252) + g(3) * t272 + ((-qJD(1) * t364 + t97) * t252 + (qJD(1) * t363 - t400) * t250 - (-t252 ^ 2 - t422) * pkin(3) * t337 + t441) * t39) * m(5) + (-(t157 * t62 + t158 * t61) * qJD(1) - (t72 * (-t157 * t250 - t158 * t252) - t292 * t298) * qJD(3) + (qJD(3) * t291 - t134 * t176 - t135 * t177) * t279 + t72 * ((-t134 * t252 + t135 * t250) * qJD(1) + t291) - t292 * t173 + ((t250 * t62 + t389) * qJD(1) + t435) * t194 - g(1) * t157 + g(2) * t158 + g(3) * t298) * m(4); t258 + (t8 * (-t116 * t250 - t104) - (t250 * t42 - t252 * t43) * t127 + ((t250 * t43 + t390) * qJD(1) + t434) * t167 - t42 * t302 - t443 - g(1) * t143 - g(2) * t144 + g(3) * t297 + (-t250 * t71 + (-t116 * t252 + t117 * t250) * qJD(1) + t441) * t39) * m(5);];
tau = t1;
