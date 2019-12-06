% Calculate vector of inverse dynamics joint torques for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:41:37
% EndTime: 2019-12-05 17:41:56
% DurationCPUTime: 11.72s
% Computational Cost: add. (14925->629), mult. (10673->782), div. (0->0), fcn. (8267->10), ass. (0->354)
t263 = qJ(1) + pkin(8);
t257 = cos(t263);
t242 = qJDD(4) * t257;
t262 = qJD(4) + qJD(5);
t255 = sin(t263);
t384 = qJD(1) * t255;
t128 = qJDD(5) * t257 - t262 * t384 + t242;
t261 = pkin(9) + qJ(4);
t258 = qJ(5) + t261;
t252 = cos(t258);
t240 = t252 * rSges(6,1);
t251 = sin(t258);
t456 = t251 * rSges(6,2);
t495 = t240 - t456;
t161 = t495 * t262;
t178 = rSges(6,1) * t251 + rSges(6,2) * t252;
t378 = qJD(1) * qJD(4);
t181 = -t255 * t378 + t242;
t189 = t257 * t262;
t268 = qJD(4) ^ 2;
t246 = qJD(3) * t257;
t324 = -t255 * pkin(2) + t257 * qJ(3);
t418 = t255 * qJ(3);
t199 = t257 * pkin(2) + t418;
t187 = qJD(1) * t199;
t398 = t187 - t246;
t336 = qJDD(1) * t324 + qJDD(3) * t255 + (t246 - t398) * qJD(1);
t267 = cos(qJ(1));
t269 = qJD(1) ^ 2;
t413 = t267 * t269;
t375 = pkin(1) * t413;
t256 = cos(t261);
t416 = t256 * t257;
t376 = pkin(4) * t416;
t266 = sin(qJ(1));
t377 = qJDD(1) * t266;
t254 = sin(t261);
t468 = pkin(4) * t254;
t462 = -pkin(6) - qJ(3);
t348 = qJD(1) * t462;
t229 = t255 * t348;
t265 = cos(pkin(9));
t414 = t257 * t265;
t298 = -pkin(3) * t414 + t418;
t469 = pkin(3) * t265;
t302 = t257 * pkin(6) - t255 * t469;
t497 = qJD(1) * (qJD(1) * t298 + t229) + qJDD(1) * t302;
t260 = -pkin(7) + t462;
t417 = t255 * t256;
t288 = (-t260 + t462) * t257 - pkin(4) * t417;
t249 = t257 * rSges(6,3);
t424 = t251 * t255;
t392 = rSges(6,2) * t424 + t249;
t422 = t252 * t255;
t304 = rSges(6,1) * t422 - t392;
t515 = t288 - t304;
t421 = t252 * t257;
t371 = rSges(6,1) * t421;
t299 = -t255 * rSges(6,3) - t371;
t146 = rSges(6,1) * t424 + rSges(6,2) * t422;
t423 = t251 * t257;
t217 = rSges(6,2) * t423;
t369 = qJD(1) * t217 + t146 * t262;
t76 = qJD(1) * t299 + t369;
t382 = qJD(4) * t255;
t363 = t254 * t382;
t210 = pkin(4) * t363;
t228 = t260 * t384;
t78 = -qJD(1) * t376 + t210 + t228 - t229;
t9 = -pkin(1) * t377 - t128 * t178 - t189 * t161 - t181 * t468 - t268 * t376 + t336 - t375 + t515 * qJDD(1) + (t76 + t78) * qJD(1) + t497;
t516 = g(3) - t9;
t208 = Icges(6,4) * t424;
t440 = Icges(6,5) * t257;
t123 = -Icges(6,1) * t422 + t208 + t440;
t209 = Icges(6,4) * t423;
t441 = Icges(6,5) * t255;
t124 = Icges(6,1) * t421 - t209 + t441;
t188 = t255 * t262;
t238 = Icges(6,4) * t252;
t315 = -Icges(6,2) * t251 + t238;
t496 = Icges(6,1) * t251 + t238;
t505 = t496 + t315;
t271 = qJD(1) * t505 + t188 * (-Icges(6,2) * t421 + t124 - t209) + t189 * (Icges(6,2) * t422 + t123 + t208);
t444 = Icges(6,4) * t251;
t174 = Icges(6,2) * t252 + t444;
t177 = Icges(6,1) * t252 - t444;
t509 = t174 - t177;
t436 = Icges(6,6) * t255;
t122 = Icges(6,4) * t421 - Icges(6,2) * t423 + t436;
t510 = t496 * t257 + t122;
t121 = Icges(6,6) * t257 - t255 * t315;
t511 = -t496 * t255 + t121;
t480 = qJD(1) * t509 + t188 * t510 + t189 * t511;
t514 = t271 * t251 + t252 * t480;
t420 = t254 * t255;
t222 = Icges(5,4) * t420;
t442 = Icges(5,5) * t257;
t134 = -Icges(5,1) * t417 + t222 + t442;
t419 = t254 * t257;
t223 = Icges(5,4) * t419;
t443 = Icges(5,5) * t255;
t135 = Icges(5,1) * t416 - t223 + t443;
t282 = t255 * (-Icges(5,2) * t416 + t135 - t223) + t257 * (Icges(5,2) * t417 + t134 + t222);
t247 = Icges(5,4) * t256;
t316 = -Icges(5,2) * t254 + t247;
t132 = Icges(5,6) * t257 - t255 * t316;
t437 = Icges(5,6) * t255;
t133 = Icges(5,4) * t416 - Icges(5,2) * t419 + t437;
t494 = Icges(5,1) * t254 + t247;
t153 = t494 * t255;
t154 = t494 * t257;
t484 = t255 * (t133 + t154) + t257 * (t132 - t153);
t513 = -t282 * t254 - t256 * t484;
t125 = -t217 - t299;
t367 = pkin(2) + t469;
t467 = pkin(4) * t256;
t308 = t367 + t467;
t167 = t257 * t308;
t236 = t255 * t462;
t333 = t257 * t367;
t95 = t255 * t260 - t167 - t236 + t333;
t502 = -t125 + t95;
t508 = qJD(1) * t146 - t178 * t384 - t189 * t495;
t427 = t178 * t257;
t75 = -t262 * t427 + (-t255 * t495 + t249) * qJD(1);
t507 = t146 * t188 + t189 * t427 - t255 * t76 + t257 * t75;
t250 = t257 * rSges(5,3);
t391 = rSges(5,2) * t420 + t250;
t466 = t266 * pkin(1);
t506 = t466 - t391;
t396 = t494 + t316;
t445 = Icges(5,4) * t254;
t192 = Icges(5,2) * t256 + t445;
t195 = Icges(5,1) * t256 - t445;
t397 = t192 - t195;
t504 = (t254 * t396 + t256 * t397) * qJD(1);
t383 = qJD(1) * t257;
t118 = t236 + t298;
t465 = t267 * pkin(1);
t352 = -t199 - t465;
t335 = t118 + t352;
t297 = t335 + t502;
t394 = t210 + t246;
t332 = t188 * t178 + t394;
t42 = qJD(1) * t297 + t332;
t503 = (-qJD(1) * t427 + t161 * t255 + t178 * t383 - t188 * t495) * t42;
t434 = Icges(6,3) * t255;
t120 = Icges(6,5) * t421 - Icges(6,6) * t423 + t434;
t491 = -t257 * t120 + t124 * t422;
t47 = t122 * t424 - t491;
t311 = t174 * t251 - t252 * t496;
t172 = Icges(6,5) * t251 + Icges(6,6) * t252;
t429 = t172 * t257;
t60 = t255 * t311 + t429;
t501 = t60 * qJD(1) + t188 * t47;
t312 = t133 * t254 - t135 * t256;
t500 = t257 * t312;
t461 = rSges(3,1) * t257;
t331 = -t461 - t465;
t166 = t255 * rSges(3,2) + t331;
t493 = t515 * qJD(1);
t245 = qJD(3) * t255;
t395 = t257 * t348 + t367 * t384;
t490 = t395 - t245;
t197 = rSges(5,1) * t254 + rSges(5,2) * t256;
t373 = qJD(1) * t466;
t303 = -qJD(1) * t324 - t245 + t373;
t296 = -qJD(1) * t302 + t303;
t372 = rSges(5,1) * t417;
t305 = t372 - t391;
t380 = qJD(4) * t257;
t55 = qJD(1) * t305 + t197 * t380 + t296;
t164 = t197 * t382;
t225 = rSges(5,2) * t419;
t455 = t255 * rSges(5,3);
t300 = -rSges(5,1) * t416 - t455;
t136 = -t225 - t300;
t307 = -t136 + t335;
t56 = qJD(1) * t307 + t164 + t246;
t489 = t255 * t56 + t257 * t55;
t435 = Icges(5,3) * t255;
t131 = Icges(5,5) * t416 - Icges(5,6) * t419 + t435;
t65 = t133 * t256 + t135 * t254;
t84 = qJD(1) * t132 - t192 * t380;
t86 = -qJD(4) * t154 + (-t195 * t255 + t442) * qJD(1);
t488 = -qJD(1) * t131 + qJD(4) * t65 + t254 * t84 - t256 * t86;
t191 = Icges(5,5) * t256 - Icges(5,6) * t254;
t130 = Icges(5,3) * t257 - t191 * t255;
t387 = qJD(1) * t130;
t64 = t132 * t256 + t134 * t254;
t85 = t192 * t382 + (-t257 * t316 - t437) * qJD(1);
t87 = qJD(4) * t153 + (-t195 * t257 - t443) * qJD(1);
t487 = qJD(4) * t64 + t254 * t85 - t256 * t87 - t387;
t169 = t316 * qJD(4);
t170 = t195 * qJD(4);
t190 = Icges(5,5) * t254 + Icges(5,6) * t256;
t310 = t192 * t256 + t254 * t494;
t486 = -qJD(1) * t190 + qJD(4) * t310 + t169 * t254 - t170 * t256;
t339 = t509 * t262;
t340 = t505 * t262;
t483 = -qJD(1) * t172 + t251 * t340 + t252 * t339;
t344 = qJD(1) * t121 + t124 * t262 - t174 * t189;
t346 = -(-t177 * t255 + t440) * qJD(1) + t510 * t262;
t482 = -qJD(1) * t120 + t251 * t344 + t252 * t346;
t345 = t123 * t262 + t174 * t188 + (-t257 * t315 - t436) * qJD(1);
t347 = -(-t177 * t257 - t441) * qJD(1) + t511 * t262;
t173 = Icges(6,5) * t252 - Icges(6,6) * t251;
t119 = Icges(6,3) * t257 - t173 * t255;
t389 = qJD(1) * t119;
t481 = t251 * t345 + t252 * t347 - t389;
t180 = qJDD(4) * t255 + t257 * t378;
t127 = qJD(5) * t383 + qJDD(5) * t255 + t180;
t479 = t127 / 0.2e1;
t478 = t128 / 0.2e1;
t477 = t180 / 0.2e1;
t476 = t181 / 0.2e1;
t475 = -t188 / 0.2e1;
t474 = t188 / 0.2e1;
t473 = -t189 / 0.2e1;
t472 = t189 / 0.2e1;
t471 = t255 / 0.2e1;
t470 = t257 / 0.2e1;
t464 = -qJD(1) / 0.2e1;
t463 = qJD(1) / 0.2e1;
t460 = rSges(4,1) * t265;
t459 = rSges(5,1) * t256;
t264 = sin(pkin(9));
t458 = rSges(4,2) * t264;
t452 = qJDD(1) / 0.2e1;
t451 = rSges(4,3) + qJ(3);
t450 = -t255 * t119 - t123 * t421;
t449 = t257 * t119 + t121 * t424;
t448 = t260 - rSges(6,3);
t433 = t120 * t255;
t432 = t132 * t254;
t431 = t134 * t256;
t140 = t172 * t255;
t428 = t178 * t255;
t149 = t190 * t255;
t426 = t190 * t257;
t156 = t197 * t255;
t425 = t197 * t257;
t415 = t257 * t260;
t61 = -t257 * t311 + t140;
t412 = t61 * qJD(1);
t309 = t254 * t192 - t256 * t494;
t74 = -t257 * t309 + t149;
t411 = t74 * qJD(1);
t410 = t257 * t130 + t132 * t420;
t409 = t255 * t130 + t134 * t416;
t253 = t269 * t466;
t390 = qJDD(3) * t257 + t253;
t385 = qJD(1) * t191;
t381 = qJD(4) * t256;
t379 = -m(4) - m(5) - m(6);
t226 = pkin(4) * t420;
t374 = qJD(4) * t468;
t231 = t257 * t458;
t368 = rSges(5,1) * t363 + (t254 * t383 + t255 * t381) * rSges(5,2);
t362 = t254 * t380;
t360 = t256 * t380;
t359 = -pkin(2) - t460;
t358 = -t384 / 0.2e1;
t357 = t383 / 0.2e1;
t356 = -t382 / 0.2e1;
t355 = t382 / 0.2e1;
t354 = -t380 / 0.2e1;
t353 = t380 / 0.2e1;
t342 = -t131 - t431;
t337 = pkin(4) * t362;
t301 = -rSges(4,1) * t414 - rSges(4,3) * t255;
t139 = -t231 - t301;
t334 = -t139 + t352;
t239 = pkin(2) * t384;
t330 = qJ(3) * t383 - t239;
t155 = t245 + t330;
t329 = t330 - t155 + t490;
t233 = rSges(2,1) * t267 - t266 * rSges(2,2);
t328 = rSges(2,1) * t266 + rSges(2,2) * t267;
t327 = rSges(3,1) * t255 + rSges(3,2) * t257;
t326 = -t458 + t460;
t198 = -rSges(5,2) * t254 + t459;
t313 = -t431 + t432;
t284 = qJD(4) * t149 + (-t191 * t257 + t313 - t435) * qJD(1);
t285 = qJD(1) * t312 - qJD(4) * t426 + t387;
t323 = (t284 * t255 - t257 * t487) * t257 + (t285 * t255 - t257 * t488) * t255;
t322 = (t255 * t487 + t284 * t257) * t257 + (t255 * t488 + t285 * t257) * t255;
t51 = -t134 * t417 + t410;
t52 = t257 * t131 + t133 * t420 - t135 * t417;
t321 = t255 * t52 + t257 * t51;
t53 = -t132 * t419 + t409;
t54 = t131 * t255 - t500;
t320 = t255 * t54 + t257 * t53;
t88 = -qJD(4) * t425 + (-t198 * t255 + t250) * qJD(1);
t89 = qJD(1) * t300 + t368;
t319 = -t255 * t89 + t257 * t88;
t59 = t122 * t252 + t124 * t251;
t314 = t122 * t251 - t124 * t252;
t306 = -t367 - t459;
t48 = -t121 * t423 - t450;
t294 = t327 + t466;
t291 = -t308 - t240;
t290 = qJD(1) * t173 + t140 * t189 - t188 * t429;
t287 = qJD(1) * t314 - t262 * t429 + t389;
t286 = t262 * t140 + (t121 * t251 - t123 * t252 - t173 * t257 - t434) * qJD(1);
t281 = qJD(1) * t311 + t173 * t262;
t280 = t309 * qJD(1) + t191 * qJD(4);
t279 = t257 * t125 + t255 * t304;
t278 = rSges(4,3) * t257 - t326 * t255;
t277 = t255 * t288;
t275 = t257 * t136 + t255 * t305;
t11 = t286 * t255 - t257 * t481;
t12 = t287 * t255 - t257 * t482;
t13 = t255 * t481 + t286 * t257;
t14 = t255 * t482 + t287 * t257;
t46 = -t123 * t422 + t449;
t22 = t189 * t46 + t501;
t49 = -t257 * t314 + t433;
t23 = t188 * t49 + t189 * t48 + t412;
t30 = t281 * t255 - t257 * t483;
t31 = t255 * t483 + t281 * t257;
t32 = -t251 * t347 + t252 * t345;
t33 = -t251 * t346 + t252 * t344;
t58 = t121 * t252 + t123 * t251;
t274 = (qJD(1) * t30 + qJDD(1) * t61 + t11 * t189 + t12 * t188 + t127 * t49 + t128 * t48) * t471 + (-t251 * t480 + t252 * t271) * t464 + (qJD(1) * t31 + qJDD(1) * t60 + t127 * t47 + t128 * t46 + t13 * t189 + t14 * t188) * t470 + t22 * t358 + t23 * t357 + (t255 * t49 + t257 * t48) * t479 + (t255 * t47 + t257 * t46) * t478 + (t11 * t257 + t12 * t255 + (-t255 * t48 + t257 * t49) * qJD(1)) * t474 + (t13 * t257 + t14 * t255 + (-t255 * t46 + t257 * t47) * qJD(1)) * t472 + (t255 * t59 + t257 * t58) * t452 + (t255 * t33 + t257 * t32 + (-t255 * t58 + t257 * t59) * qJD(1)) * t463 + (t290 * t255 - t257 * t514) * t475 + (t255 * t514 + t290 * t257) * t473;
t273 = (-t377 - t413) * pkin(1) + t336;
t270 = t122 * t423 - t433 + (-t123 * t255 - t124 * t257) * t252 + t449;
t237 = rSges(3,2) * t384;
t227 = qJD(1) * t231;
t171 = t198 * qJD(4);
t162 = t327 * qJD(1) + t373;
t116 = qJD(1) * t118;
t91 = qJD(1) * t334 + t246;
t90 = -qJD(1) * t278 + t303;
t77 = -t337 + (-t255 * t308 - t415) * qJD(1) + t395;
t73 = t255 * t309 + t426;
t66 = t73 * qJD(1);
t62 = qJD(4) * t275 + qJD(2);
t45 = t334 * qJDD(1) + (-rSges(4,3) * t383 - t155 + (qJD(1) * t326 - qJD(3)) * t255) * qJD(1) + t390;
t44 = qJDD(1) * t278 + (qJD(1) * t301 + t227) * qJD(1) + t273;
t41 = t178 * t189 + t296 + t337 - t493;
t39 = t255 * t486 + t280 * t257;
t38 = t280 * t255 - t257 * t486;
t37 = -qJD(4) * t277 + t189 * t125 + t188 * t304 - t380 * t95 + qJD(2);
t36 = -t312 * qJD(4) + t254 * t86 + t256 * t84;
t35 = -qJD(4) * t313 + t254 * t87 + t256 * t85;
t34 = qJD(4) * t319 + t181 * t136 + t180 * t305 + qJDD(2);
t29 = t171 * t382 + t180 * t197 + t307 * qJDD(1) + (t329 - t88) * qJD(1) + t390;
t28 = qJD(1) * t89 - qJDD(1) * t305 - t171 * t380 - t181 * t197 + t273 + t497;
t25 = qJD(4) * t320 + t411;
t24 = qJD(4) * t321 + t66;
t10 = t127 * t178 + t161 * t188 + (t180 * t254 + t268 * t417) * pkin(4) + t297 * qJDD(1) + (t329 - t75 - t77) * qJD(1) + t390;
t8 = t128 * t125 + t127 * t304 - t180 * t288 - t181 * t95 - t188 * t76 + t189 * t75 + t380 * t77 - t382 * t78 + qJDD(2);
t1 = [((t49 + t270) * t189 + t501) * t475 + (t66 + ((t410 + t54 + t500) * t257 + (-t53 + (t342 - t432) * t257 + t52 + t409) * t255) * qJD(4)) * t356 + (t59 + t61) * t479 + (t58 + t60) * t478 + (t65 + t74) * t477 + (t64 + t73) * t476 + (t23 - t412 + (t47 + (t121 * t257 - t122 * t255) * t251 + t450 + t491) * t189 + (-t46 + t270) * t188) * t473 + (t31 + t32) * t472 + (t25 - t411 + ((t52 + (-t131 + t432) * t257 - t409) * t257 + (t255 * t342 + t410 - t51) * t255) * qJD(4)) * t354 + (t35 + t39) * t353 + (-qJD(4) * t309 + t169 * t256 + t170 * t254 - t251 * t339 + t252 * t340) * qJD(1) + (-t162 * t237 + (qJD(1) * t166 * t294 - t162 * t331) * qJD(1) + (t269 * t327 - g(2) + t253) * t166 + (g(3) + t375 + (t461 * qJD(1) - t237) * qJD(1)) * t294) * m(3) + (g(2) * t233 + g(3) * t328) * m(2) + (t33 + t30 + t22) * t474 + (t36 + t38 + t24) * t355 + (-t42 * t245 - t41 * (t228 + t369 + t394) + t42 * (t178 * t262 + t374) * t257 + ((t42 * t266 + t41 * t267) * pkin(1) + (-t291 * t41 + t42 * t448) * t257 + (t42 * (-t291 - t456) + t41 * rSges(6,3)) * t255) * qJD(1) - (-t116 + t187 + t42 + (t465 - t502) * qJD(1) - t332) * t41 + (-g(2) + t10) * (t255 * t448 - t167 + t217 - t371 - t465) - t516 * (t255 * t291 + t392 - t415 - t466)) * m(6) + ((t29 - g(2)) * (t257 * t306 + t225 + t236 - t455 - t465) + (t28 - g(3)) * (t255 * t306 - t257 * t462 - t506) + (rSges(5,1) * t362 + rSges(5,2) * t360 + t490 + (t372 + t506) * qJD(1)) * t56 + (-t229 - t246 - t368 + t116 + t164 - t398 - t56 + (-t300 + t333 - t136) * qJD(1)) * t55) * m(5) + (-(t91 + (t139 + t465) * qJD(1) + t398) * t90 + t91 * (t239 - t245) - t90 * (t227 + t246) + ((t91 * t266 + t90 * t267) * pkin(1) + (-t359 * t90 - t451 * t91) * t257 + (t326 * t91 + t451 * t90) * t255) * qJD(1) + (t44 - g(3)) * (-t466 + t451 * t257 + (-pkin(2) - t326) * t255) + (t45 - g(2)) * (-t451 * t255 + t359 * t257 + t231 - t465)) * m(4) + (t310 + Icges(4,2) * t265 ^ 2 + (Icges(4,1) * t264 + 0.2e1 * Icges(4,4) * t265) * t264 + t174 * t252 + t496 * t251 + m(2) * (t233 ^ 2 + t328 ^ 2) + m(3) * (t166 ^ 2 + t294 ^ 2) + Icges(2,3) + Icges(3,3)) * qJDD(1); (m(3) + m(4)) * qJDD(2) + m(5) * t34 + m(6) * t8 + (-m(3) + t379) * g(1); t379 * (g(2) * t257 + g(3) * t255) + m(4) * (t255 * t44 + t257 * t45) + m(5) * (t255 * t28 + t257 * t29) + m(6) * (t10 * t257 + t255 * t9); ((-t51 * t255 + t52 * t257) * qJD(1) + t322) * t353 + (qJD(1) * t39 + qJD(4) * t322 + qJDD(1) * t73 + t180 * t52 + t181 * t51) * t470 + (qJD(1) * t38 + qJD(4) * t323 + qJDD(1) * t74 + t180 * t54 + t181 * t53) * t471 + ((-t382 * t426 + t385) * t255 + (-t504 + (t149 * t255 + t513) * qJD(4)) * t257) * t356 + ((t149 * t380 + t385) * t257 + (t504 + (-t257 * t426 - t513) * qJD(4)) * t255) * t354 + ((-t254 * t397 + t256 * t396) * qJD(1) + (-t254 * t484 + t256 * t282) * qJD(4)) * t464 + t320 * t477 + t321 * t476 + (t255 * t65 + t257 * t64) * t452 + ((-t53 * t255 + t54 * t257) * qJD(1) + t323) * t355 + t274 + (t255 * t36 + t257 * t35 + (-t255 * t64 + t65 * t257) * qJD(1)) * t463 + t24 * t358 + t25 * t357 + (-g(1) * (t495 + t467) - g(2) * (t226 + t146) + t8 * (-t277 + t279) + t10 * (t226 + t428) + (-t8 * t95 - t516 * (-t178 - t468)) * t257 + (-(-t255 ^ 2 - t257 ^ 2) * t374 - t255 * t78 + (t77 - t493) * t257 + t502 * t384 + t507) * t37 + t503 + (-pkin(4) * t360 + (pkin(4) * t381 + t161) * t257 + t508) * t41) * m(6) + (-g(1) * t198 - g(2) * t156 + g(3) * t425 + t34 * t275 + t62 * ((-t255 * t136 + t257 * t305) * qJD(1) + t319) + t29 * t156 - t28 * t425 + (t383 * t56 - t384 * t55) * t197 + t489 * t171 - (-t156 * t55 + t425 * t56) * qJD(1) - (t62 * (-t156 * t255 - t257 * t425) + t489 * t198) * qJD(4)) * m(5); t274 + (-g(1) * t495 - g(2) * t146 + t10 * t428 + t8 * t279 + ((-t255 * t125 + t257 * t304) * qJD(1) + t507) * t37 + t503 + (t161 * t257 + t508) * t41 + t516 * t427) * m(6);];
tau = t1;
