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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:27:28
% EndTime: 2020-01-03 11:28:04
% DurationCPUTime: 12.88s
% Computational Cost: add. (14925->645), mult. (10673->801), div. (0->0), fcn. (8267->10), ass. (0->345)
t264 = qJ(1) + pkin(8);
t253 = sin(t264);
t255 = cos(t264);
t194 = rSges(3,1) * t253 + rSges(3,2) * t255;
t268 = sin(qJ(1));
t259 = t268 * pkin(1);
t163 = t259 + t194;
t262 = pkin(9) + qJ(4);
t256 = qJ(5) + t262;
t245 = sin(t256);
t246 = cos(t256);
t402 = t245 * t255;
t203 = Icges(6,4) * t402;
t400 = t246 * t255;
t420 = Icges(6,5) * t253;
t119 = Icges(6,1) * t400 - t203 + t420;
t263 = qJD(4) + qJD(5);
t183 = t253 * t263;
t394 = t255 * t263;
t235 = Icges(6,4) * t246;
t304 = -Icges(6,2) * t245 + t235;
t472 = Icges(6,1) * t245 + t235;
t478 = t472 + t304;
t422 = Icges(6,4) * t245;
t176 = Icges(6,1) * t246 - t422;
t118 = -Icges(6,5) * t255 + t176 * t253;
t173 = Icges(6,2) * t246 + t422;
t481 = -t173 * t253 + t118;
t272 = qJD(1) * t478 + t183 * (-Icges(6,2) * t400 + t119 - t203) - t394 * t481;
t480 = t173 - t176;
t416 = Icges(6,6) * t253;
t117 = Icges(6,4) * t400 - Icges(6,2) * t402 + t416;
t482 = t255 * t472 + t117;
t116 = -Icges(6,6) * t255 + t253 * t304;
t483 = t253 * t472 + t116;
t459 = qJD(1) * t480 + t183 * t482 - t394 * t483;
t486 = t272 * t245 + t246 * t459;
t252 = sin(t262);
t254 = cos(t262);
t423 = Icges(5,4) * t252;
t190 = Icges(5,1) * t254 - t423;
t129 = -Icges(5,5) * t255 + t190 * t253;
t398 = t252 * t255;
t214 = Icges(5,4) * t398;
t396 = t254 * t255;
t421 = Icges(5,5) * t253;
t130 = Icges(5,1) * t396 - t214 + t421;
t187 = Icges(5,2) * t254 + t423;
t149 = t187 * t253;
t280 = t253 * (-Icges(5,2) * t396 + t130 - t214) - t255 * (t129 - t149);
t239 = Icges(5,4) * t254;
t305 = -Icges(5,2) * t252 + t239;
t127 = -Icges(5,6) * t255 + t253 * t305;
t417 = Icges(5,6) * t253;
t128 = Icges(5,4) * t396 - Icges(5,2) * t398 + t417;
t470 = Icges(5,1) * t252 + t239;
t151 = t470 * t253;
t152 = t470 * t255;
t463 = t253 * (t128 + t152) - t255 * (t127 + t151);
t485 = -t280 * t252 - t254 * t463;
t437 = rSges(6,2) * t246;
t441 = rSges(6,1) * t245;
t177 = t437 + t441;
t484 = t177 * t394;
t265 = sin(pkin(9));
t440 = rSges(4,2) * t265;
t266 = cos(pkin(9));
t443 = rSges(4,1) * t266;
t479 = (t440 - t443) * t253;
t378 = t470 + t305;
t379 = t187 - t190;
t477 = (t252 * t378 + t254 * t379) * qJD(1);
t448 = pkin(4) * t252;
t172 = Icges(6,5) * t246 - Icges(6,6) * t245;
t114 = -Icges(6,3) * t255 + t172 * t253;
t410 = t117 * t245;
t303 = -t119 * t246 + t410;
t294 = -t114 + t303;
t476 = t394 * t294;
t267 = -pkin(6) - qJ(3);
t227 = t253 * t267;
t412 = qJ(3) * t253;
t247 = pkin(3) * t266 + pkin(2);
t445 = pkin(2) - t247;
t113 = t255 * t445 + t227 + t412;
t196 = pkin(2) * t255 + t412;
t182 = qJD(1) * t196;
t473 = -qJD(1) * t113 + t182;
t438 = rSges(6,2) * t245;
t471 = rSges(6,1) * t246 - t438;
t436 = rSges(4,3) * t255;
t469 = t436 + t479;
t143 = t177 * t253;
t144 = rSges(6,1) * t402 + rSges(6,2) * t400;
t401 = t246 * t253;
t403 = t245 * t253;
t434 = rSges(6,3) * t255;
t120 = rSges(6,1) * t401 - rSges(6,2) * t403 - t434;
t350 = rSges(6,1) * t400;
t318 = -rSges(6,2) * t402 + t350;
t121 = rSges(6,3) * t253 + t318;
t228 = t255 * t267;
t373 = t247 * t253 + t228;
t244 = pkin(4) * t254;
t209 = t244 + t247;
t261 = -pkin(7) + t267;
t382 = t209 * t253 + t255 * t261;
t94 = -t373 + t382;
t166 = t255 * t209;
t95 = t247 * t255 + t253 * t261 - t166 - t227;
t37 = t120 * t183 + t121 * t394 + qJD(2) + (t253 * t94 - t255 * t95) * qJD(4);
t238 = qJD(3) * t253;
t243 = t253 * pkin(2);
t193 = -qJ(3) * t255 + t243;
t112 = -t193 + t373;
t335 = t193 + t259;
t319 = t112 + t335;
t357 = qJD(4) * t255;
t322 = t357 * t448;
t427 = t120 + t94;
t41 = t322 + t484 - t238 + (t319 + t427) * qJD(1);
t269 = cos(qJ(1));
t260 = t269 * pkin(1);
t251 = qJD(1) * t260;
t360 = qJD(3) * t255;
t323 = t251 - t360;
t359 = qJD(4) * t253;
t344 = t252 * t359;
t275 = -pkin(4) * t344 - t177 * t183 + t323;
t376 = t196 - t113;
t426 = t121 - t95;
t321 = t376 + t426;
t42 = qJD(1) * t321 + t275;
t468 = -t41 * (-qJD(1) * t143 + t394 * t471) - t37 * (-t143 * t183 - t144 * t394) - t42 * (-qJD(1) * t144 - t183 * t471);
t186 = Icges(5,5) * t254 - Icges(5,6) * t252;
t125 = -Icges(5,3) * t255 + t186 * t253;
t365 = qJD(1) * t125;
t64 = t127 * t254 + t129 * t252;
t85 = -qJD(4) * t149 + (t255 * t305 + t417) * qJD(1);
t87 = -qJD(4) * t151 + (t190 * t255 + t421) * qJD(1);
t467 = qJD(4) * t64 + t252 * t85 - t254 * t87 - t365;
t168 = t305 * qJD(4);
t169 = t190 * qJD(4);
t185 = Icges(5,5) * t252 + Icges(5,6) * t254;
t297 = t187 * t254 + t252 * t470;
t466 = -qJD(1) * t185 + qJD(4) * t297 + t168 * t252 - t169 * t254;
t415 = Icges(5,3) * t253;
t126 = Icges(5,5) * t396 - Icges(5,6) * t398 + t415;
t301 = t128 * t254 + t130 * t252;
t84 = qJD(1) * t127 + t187 * t357;
t86 = qJD(1) * t129 + qJD(4) * t152;
t465 = -qJD(1) * t126 + qJD(4) * t301 - t252 * t84 + t254 * t86;
t171 = Icges(6,5) * t245 + Icges(6,6) * t246;
t324 = t480 * t263;
t325 = t478 * t263;
t462 = -qJD(1) * t171 + t245 * t325 + t246 * t324;
t414 = Icges(6,3) * t253;
t115 = Icges(6,5) * t400 - Icges(6,6) * t402 + t414;
t329 = -qJD(1) * t116 + t119 * t263 - t173 * t394;
t331 = qJD(1) * t118 + t263 * t482;
t461 = -qJD(1) * t115 + t245 * t329 + t246 * t331;
t330 = (t255 * t304 + t416) * qJD(1) + t481 * t263;
t332 = -(t176 * t255 + t420) * qJD(1) + t483 * t263;
t367 = qJD(1) * t114;
t460 = t245 * t330 + t246 * t332 - t367;
t271 = qJD(1) ^ 2;
t354 = -qJDD(4) - qJDD(5);
t122 = -qJD(1) * t394 + t253 * t354;
t458 = t122 / 0.2e1;
t355 = qJD(1) * qJD(4);
t225 = t253 * t355;
t362 = qJD(1) * t253;
t123 = qJD(5) * t362 + t255 * t354 + t225;
t457 = t123 / 0.2e1;
t179 = -qJDD(4) * t253 - t255 * t355;
t456 = t179 / 0.2e1;
t180 = -qJDD(4) * t255 + t225;
t455 = t180 / 0.2e1;
t454 = -t183 / 0.2e1;
t453 = t183 / 0.2e1;
t452 = t394 / 0.2e1;
t451 = -t394 / 0.2e1;
t450 = -t253 / 0.2e1;
t449 = -t255 / 0.2e1;
t447 = -qJD(1) / 0.2e1;
t446 = qJD(1) / 0.2e1;
t75 = t484 + (t253 * t471 - t434) * qJD(1);
t368 = t261 - t267;
t77 = t322 + (t368 * t255 + (t209 - t247) * t253) * qJD(1);
t444 = -t75 - t77;
t442 = rSges(5,1) * t254;
t439 = rSges(5,2) * t252;
t435 = rSges(5,3) * t255;
t192 = rSges(5,1) * t252 + rSges(5,2) * t254;
t154 = t192 * t253;
t397 = t253 * t254;
t399 = t252 * t253;
t131 = rSges(5,1) * t397 - rSges(5,2) * t399 - t435;
t55 = t192 * t357 - t238 + (t131 + t319) * qJD(1);
t433 = t154 * t55;
t432 = t253 * t41;
t291 = -t192 * t359 + t323;
t333 = -rSges(5,2) * t398 + rSges(5,3) * t253;
t352 = rSges(5,1) * t396;
t132 = t333 + t352;
t345 = t132 + t376;
t56 = qJD(1) * t345 + t291;
t431 = t255 * t56;
t430 = qJDD(1) / 0.2e1;
t429 = rSges(4,3) + qJ(3);
t428 = rSges(6,3) - t261;
t413 = pkin(1) * qJDD(1);
t411 = t116 * t245;
t409 = t118 * t246;
t408 = t127 * t252;
t407 = t128 * t252;
t406 = t129 * t254;
t405 = t171 * t253;
t138 = t171 * t255;
t404 = t185 * t253;
t148 = t185 * t255;
t395 = t254 * qJD(4) ^ 2;
t61 = -t173 * t402 + t400 * t472 + t405;
t393 = t61 * qJD(1);
t74 = -t187 * t398 + t396 * t470 + t404;
t392 = t74 * qJD(1);
t361 = qJD(1) * t255;
t230 = qJ(3) * t361;
t370 = t230 + t238;
t153 = pkin(2) * t362 - t370;
t391 = -t230 - (-t253 * t445 + t228) * qJD(1) - t153;
t377 = rSges(6,3) * t362 + qJD(1) * t350;
t221 = t255 * t440;
t353 = t255 * t443;
t136 = rSges(4,3) * t253 - t221 + t353;
t375 = t196 + t136;
t374 = rSges(5,3) * t362 + qJD(1) * t352;
t372 = rSges(4,3) * t362 + qJD(1) * t353;
t371 = pkin(2) * t361 + qJ(3) * t362;
t369 = t260 * t271 + t268 * t413;
t363 = qJD(1) * t186;
t358 = qJD(4) * t254;
t356 = -m(4) - m(5) - m(6);
t351 = t263 * t441;
t348 = qJD(4) * t448;
t76 = -t253 * t351 + (-t183 * t246 - t245 * t361) * rSges(6,2) + t377;
t346 = t120 * t361 - t121 * t362 + t253 * t76;
t343 = t253 * t358;
t342 = pkin(2) + t443;
t341 = t362 / 0.2e1;
t340 = -t361 / 0.2e1;
t339 = -t359 / 0.2e1;
t338 = t359 / 0.2e1;
t337 = -t357 / 0.2e1;
t336 = t357 / 0.2e1;
t334 = -t177 - t448;
t197 = rSges(3,1) * t255 - rSges(3,2) * t253;
t328 = -t115 - t411;
t327 = -t115 + t409;
t320 = -t259 * t271 + t269 * t413;
t316 = -t360 + t371;
t223 = rSges(2,1) * t269 - rSges(2,2) * t268;
t222 = rSges(2,1) * t268 + rSges(2,2) * t269;
t195 = -t439 + t442;
t302 = t406 - t408;
t282 = qJD(4) * t404 + (-t186 * t255 + t302 - t415) * qJD(1);
t300 = -t130 * t254 + t407;
t283 = qJD(1) * t300 - qJD(4) * t148 - t365;
t313 = -(t253 * t282 + t255 * t467) * t255 - (t253 * t283 - t255 * t465) * t253;
t312 = -(-t253 * t467 + t255 * t282) * t255 - (t253 * t465 + t255 * t283) * t253;
t103 = t129 * t397;
t51 = -t125 * t255 - t127 * t399 + t103;
t104 = t130 * t397;
t52 = t126 * t255 + t128 * t399 - t104;
t311 = -t253 * t52 - t255 * t51;
t105 = t127 * t398;
t53 = -t125 * t253 - t129 * t396 + t105;
t54 = t126 * t253 - t300 * t255;
t310 = -t253 * t54 - t255 * t53;
t309 = -t253 * t56 + t255 * t55;
t155 = t192 * t255;
t88 = qJD(4) * t155 + (t195 * t253 - t435) * qJD(1);
t89 = -rSges(5,1) * t344 + (-t252 * t361 - t343) * rSges(5,2) + t374;
t308 = t253 * t89 - t255 * t88;
t59 = -t117 * t246 - t119 * t245;
t299 = t131 * t253 + t132 * t255;
t298 = -t173 * t245 + t246 * t472;
t296 = -t187 * t252 + t254 * t470;
t293 = (-t255 * t42 - t432) * qJD(1);
t292 = qJD(1) * t316 + qJDD(1) * t193 - qJDD(3) * t253 + t369;
t288 = -qJD(1) * t172 + t138 * t183 - t394 * t405;
t286 = qJD(1) * t238 - qJDD(3) * t255 + t320;
t285 = qJD(1) * t303 - t138 * t263 - t367;
t284 = t263 * t405 + (-t172 * t255 + t409 - t411 - t414) * qJD(1);
t279 = qJD(1) * t298 - t172 * t263;
t278 = t296 * qJD(1) - qJD(4) * t186;
t199 = t247 * t361;
t277 = qJD(1) * (-t267 * t362 + t199 - t371) + qJDD(1) * t112 + t292;
t11 = t253 * t284 + t255 * t460;
t12 = t253 * t285 - t255 * t461;
t13 = -t253 * t460 + t255 * t284;
t14 = t253 * t461 + t255 * t285;
t96 = t118 * t401;
t46 = -t114 * t255 - t116 * t403 + t96;
t97 = t119 * t401;
t47 = t115 * t255 + t117 * t403 - t97;
t60 = t253 * t298 - t138;
t57 = t60 * qJD(1);
t22 = -t183 * t47 - t394 * t46 + t57;
t98 = t116 * t402;
t48 = -t114 * t253 - t118 * t400 + t98;
t49 = t115 * t253 - t255 * t303;
t23 = -t183 * t49 - t394 * t48 - t393;
t30 = t253 * t279 + t255 * t462;
t31 = -t253 * t462 + t255 * t279;
t32 = -t245 * t332 + t246 * t330;
t33 = t245 * t331 - t246 * t329;
t58 = t116 * t246 + t118 * t245;
t276 = (qJD(1) * t30 - qJDD(1) * t61 - t11 * t394 - t12 * t183 + t122 * t49 + t123 * t48) * t450 + (-t245 * t459 + t246 * t272) * t447 + (qJD(1) * t31 + qJDD(1) * t60 + t122 * t47 + t123 * t46 - t13 * t394 - t14 * t183) * t449 + t22 * t341 + t23 * t340 + (-t253 * t49 - t255 * t48) * t458 + (-t253 * t47 - t255 * t46) * t457 + (-t11 * t255 - t12 * t253 + (t253 * t48 - t255 * t49) * qJD(1)) * t454 + (-t13 * t255 - t14 * t253 + (t253 * t46 - t255 * t47) * qJD(1)) * t451 + (-t253 * t59 - t255 * t58) * t430 + (-t253 * t33 - t255 * t32 + (t253 * t58 - t255 * t59) * qJD(1)) * t446 + (t288 * t253 + t255 * t486) * t453 + (-t253 * t486 + t288 * t255) * t452;
t218 = pkin(4) * t398;
t170 = t195 * qJD(4);
t164 = t197 + t260;
t162 = t209 * t361;
t159 = t471 * t263;
t109 = t253 * t120;
t91 = qJD(1) * t375 + t323;
t78 = t162 - t199 + (-qJD(1) * t368 - t348) * t253;
t73 = t253 * t296 - t148;
t66 = t73 * qJD(1);
t62 = qJD(4) * t299 + qJD(2);
t45 = t375 * qJDD(1) + (qJD(1) * t469 - t153) * qJD(1) + t286;
t44 = -qJDD(1) * t469 + ((-qJD(1) * t440 - qJD(3)) * t255 + t372) * qJD(1) + t292;
t39 = -t253 * t466 + t255 * t278;
t38 = t253 * t278 + t255 * t466;
t36 = qJD(4) * t300 + t252 * t86 + t254 * t84;
t35 = qJD(4) * t302 + t252 * t87 + t254 * t85;
t34 = qJD(4) * t308 - t131 * t179 - t132 * t180 + qJDD(2);
t29 = -t170 * t359 + t179 * t192 + t345 * qJDD(1) + (-t88 + t391) * qJD(1) + t286;
t28 = t170 * t357 + qJDD(1) * t131 - t180 * t192 + (t89 - t360) * qJD(1) + t277;
t25 = qJD(4) * t310 - t392;
t24 = qJD(4) * t311 + t66;
t10 = t122 * t177 - t159 * t183 + (t179 * t252 - t253 * t395) * pkin(4) + t321 * qJDD(1) + (t391 + t444) * qJD(1) + t286;
t9 = -t123 * t177 + t159 * t394 + t427 * qJDD(1) + (-t180 * t252 + t255 * t395) * pkin(4) + (t76 + t78 - t360) * qJD(1) + t277;
t8 = -t120 * t122 - t121 * t123 - t179 * t94 + t180 * t95 + t183 * t76 - t394 * t75 + qJDD(2) + (t253 * t78 - t255 * t77) * qJD(4);
t1 = [-t301 * t456 + (t66 + ((t53 + t104 - t105 + (t125 - t407) * t253) * t253 + (-t103 - t54 + (t125 - t300) * t255 + (t406 + t408) * t253) * t255) * qJD(4)) * t338 + t59 * t458 + (t57 - (t253 * t328 + t49 + t96) * t394 + (t48 + t97 - t98 + (t114 - t410) * t253) * t183 + (t183 * t327 - t476) * t255) * t453 - m(2) * (g(2) * t223 + g(3) * t222) - t122 * t61 / 0.2e1 - t179 * t74 / 0.2e1 + (t58 + t60) * t457 + (t64 + t73) * t455 + (t23 + t393 - (t47 - t98) * t394 + (t46 - t96) * t183 + (-t183 * t294 - t327 * t394) * t255 + (-t183 * t328 + t476) * t253) * t452 + (t32 + t31) * t451 + (t35 + t39) * t337 + (t392 + ((t105 - t52 + (t126 - t406) * t255) * t255 + (-t103 + t51 + (t126 + t408) * t253) * t253) * qJD(4) + t25) * t336 + (qJD(4) * t296 + t168 * t254 + t169 * t252 - t245 * t324 + t246 * t325) * qJD(1) + ((-t194 * t271 - g(2) + t320) * t164 + (t369 - g(3) + (0.2e1 * rSges(3,1) * t361 - 0.2e1 * rSges(3,2) * t362 - qJD(1) * t197) * qJD(1)) * t163) * m(3) + (t33 + t30 + t22) * t454 + (t36 + t38 + t24) * t339 + (t42 * t238 + t41 * (t162 + t251 + t377) + (t42 * (-t263 * t437 - t348 - t351) - t41 * qJD(3)) * t255 + (-t177 * t263 - t348) * t432 + (-t42 * t259 + (-t41 * t438 + t42 * t428) * t255 + (t42 * (-t209 - t471) - t41 * t261) * t253) * qJD(1) - (qJD(1) * t426 + t275 - t42 + t473) * t41 + (t10 - g(2)) * (t253 * t428 + t166 + t260 + t318) + (-g(3) + t9) * (t120 + t259 + t382)) * m(6) + (t56 * t238 + t55 * (t199 + t323 + t374) + (-t192 * t431 - t433) * qJD(4) + (-t56 * t259 + (t56 * (rSges(5,3) - t267) - t55 * t439) * t255 + (t56 * (-t195 - t247) - t55 * t267) * t253) * qJD(1) - (qJD(1) * t132 + t291 + t473 - t56) * t55 + (t29 - g(2)) * (-t227 + t260 + (t247 + t442) * t255 + t333) + (t28 - g(3)) * (t131 + t259 + t373)) * m(5) + ((t44 - g(3)) * (-t255 * t429 + t243 + t259 - t479) + (t370 + (t436 - t259 + (-t342 + t440) * t253) * qJD(1)) * t91 + (-t182 + t251 + t316 - t323 + t372 + t91 + (-t136 - t221) * qJD(1)) * (-t238 + (-t469 + t335) * qJD(1)) + (t45 - g(2)) * (t253 * t429 + t255 * t342 - t221 + t260)) * m(4) + (Icges(4,2) * t266 ^ 2 + (Icges(4,1) * t265 + 0.2e1 * Icges(4,4) * t266) * t265 + t173 * t246 + t245 * t472 + m(2) * (t222 ^ 2 + t223 ^ 2) + m(3) * (t163 * t194 + t164 * t197) + t297 + Icges(2,3) + Icges(3,3)) * qJDD(1); (m(3) + m(4)) * qJDD(2) + m(5) * t34 + m(6) * t8 + (-m(3) + t356) * g(1); t356 * (-g(2) * t255 - g(3) * t253) + m(4) * (-t253 * t44 - t255 * t45) + m(5) * (-t253 * t28 - t255 * t29) + m(6) * (-t10 * t255 - t253 * t9); (t253 * t301 - t255 * t64) * t430 + t310 * t456 + t311 * t455 + (qJD(1) * t39 + qJD(4) * t312 + qJDD(1) * t73 + t179 * t52 + t180 * t51) * t449 + (qJD(1) * t38 + qJD(4) * t313 - qJDD(1) * t74 + t179 * t54 + t180 * t53) * t450 + (-t253 * t36 - t255 * t35 + (t253 * t64 + t255 * t301) * qJD(1)) * t446 + t276 + ((t253 * t51 - t255 * t52) * qJD(1) + t312) * t337 + ((t53 * t253 - t255 * t54) * qJD(1) + t313) * t339 + ((-t357 * t404 - t363) * t255 + (-t477 + (t148 * t255 + t485) * qJD(4)) * t253) * t336 + ((t148 * t359 - t363) * t253 + (t477 + (-t253 * t404 - t485) * qJD(4)) * t255) * t338 + t24 * t341 + t25 * t340 + ((-t252 * t379 + t254 * t378) * qJD(1) + (-t252 * t463 + t254 * t280) * qJD(4)) * t447 + (-g(1) * (t471 + t244) - g(3) * (t218 + t144) - g(2) * t334 * t253 - (-t42 * t343 + (t37 * (-t253 ^ 2 - t255 ^ 2) * qJD(4) + t293) * t252) * pkin(4) + t8 * t109 + t37 * t346 + t9 * t218 + (t8 * t426 + t37 * t444 + t9 * t177 + t41 * t159 + (t334 * t42 + t37 * t94) * qJD(1)) * t255 + (t8 * t94 + t37 * t78 + t10 * t334 + t42 * (-pkin(4) * t358 - t159) + (t334 * t41 + t37 * t95) * qJD(1)) * t253 + t468) * m(6) + (t34 * t299 + t62 * ((t131 * t255 - t132 * t253) * qJD(1) + t308) + t309 * t170 + (-t29 * t253 + t28 * t255 + (-t253 * t55 - t431) * qJD(1)) * t192 - (-t155 * t56 - t433) * qJD(1) - (t62 * (-t154 * t253 - t155 * t255) + t309 * t195) * qJD(4) - g(1) * t195 + g(2) * t154 - g(3) * t155) * m(5); t276 + (t8 * (t121 * t255 + t109) + t37 * (-t255 * t75 + t346) + (-t253 * t42 + t255 * t41) * t159 + (-t10 * t253 + t255 * t9 + t293) * t177 - g(1) * t471 + g(2) * t143 - g(3) * t144 + t468) * m(6);];
tau = t1;
