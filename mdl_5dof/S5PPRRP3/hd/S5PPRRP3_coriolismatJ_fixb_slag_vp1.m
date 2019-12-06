% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRRP3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:36
% EndTime: 2019-12-05 15:10:50
% DurationCPUTime: 11.48s
% Computational Cost: add. (23054->716), mult. (60968->1071), div. (0->0), fcn. (74279->8), ass. (0->343)
t315 = sin(pkin(7));
t317 = cos(pkin(7));
t319 = sin(qJ(3));
t388 = t317 * t319;
t316 = cos(pkin(8));
t320 = cos(qJ(3));
t389 = t316 * t320;
t297 = t315 * t389 - t388;
t318 = sin(qJ(4));
t314 = sin(pkin(8));
t411 = cos(qJ(4));
t349 = t314 * t411;
t273 = t297 * t318 - t315 * t349;
t395 = t314 * t318;
t274 = t297 * t411 + t315 * t395;
t387 = t320 * t317;
t391 = t315 * t319;
t296 = t316 * t391 + t387;
t167 = Icges(5,5) * t274 - Icges(5,6) * t273 + Icges(5,3) * t296;
t169 = Icges(6,4) * t274 + Icges(6,2) * t296 + Icges(6,6) * t273;
t468 = t167 + t169;
t299 = t316 * t387 + t391;
t275 = t299 * t318 - t317 * t349;
t276 = t299 * t411 + t317 * t395;
t298 = -t315 * t320 + t316 * t388;
t168 = Icges(5,5) * t276 - Icges(5,6) * t275 + Icges(5,3) * t298;
t170 = Icges(6,4) * t276 + Icges(6,2) * t298 + Icges(6,6) * t275;
t467 = t168 + t170;
t462 = -2 * Icges(4,4);
t461 = 2 * Icges(4,4);
t188 = -Icges(5,5) * t273 - Icges(5,6) * t274;
t190 = -Icges(6,4) * t273 + Icges(6,6) * t274;
t300 = t316 * t411 + t320 * t395;
t301 = -t316 * t318 + t320 * t349;
t271 = Icges(5,4) * t273;
t175 = Icges(5,1) * t274 + Icges(5,5) * t296 - t271;
t380 = Icges(5,2) * t274 - t175 + t271;
t400 = Icges(6,5) * t273;
t173 = Icges(6,1) * t274 + Icges(6,4) * t296 + t400;
t382 = Icges(6,3) * t274 - t173 - t400;
t403 = Icges(5,4) * t274;
t171 = -Icges(5,2) * t273 + Icges(5,6) * t296 + t403;
t384 = -Icges(5,1) * t273 - t171 - t403;
t269 = Icges(6,5) * t274;
t165 = Icges(6,6) * t296 + Icges(6,3) * t273 + t269;
t386 = -Icges(6,1) * t273 + t165 + t269;
t394 = t314 * t319;
t466 = (t188 + t190) * t394 + (t384 + t386) * t301 + (t380 + t382) * t300;
t189 = -Icges(5,5) * t275 - Icges(5,6) * t276;
t191 = -Icges(6,4) * t275 + Icges(6,6) * t276;
t272 = Icges(5,4) * t275;
t176 = Icges(5,1) * t276 + Icges(5,5) * t298 - t272;
t379 = Icges(5,2) * t276 - t176 + t272;
t399 = Icges(6,5) * t275;
t174 = Icges(6,1) * t276 + Icges(6,4) * t298 + t399;
t381 = Icges(6,3) * t276 - t174 - t399;
t402 = Icges(5,4) * t276;
t172 = -Icges(5,2) * t275 + Icges(5,6) * t298 + t402;
t383 = -Icges(5,1) * t275 - t172 - t402;
t270 = Icges(6,5) * t276;
t166 = Icges(6,6) * t298 + Icges(6,3) * t275 + t270;
t385 = -Icges(6,1) * t275 + t166 + t270;
t465 = (t189 + t191) * t394 + (t383 + t385) * t301 + (t379 + t381) * t300;
t464 = -Icges(4,1) + Icges(4,2);
t260 = -Icges(5,5) * t300 - Icges(5,6) * t301;
t261 = -Icges(6,4) * t300 + Icges(6,6) * t301;
t295 = Icges(5,4) * t300;
t238 = Icges(5,1) * t301 + Icges(5,5) * t394 - t295;
t367 = Icges(5,2) * t301 - t238 + t295;
t398 = Icges(6,5) * t300;
t237 = Icges(6,1) * t301 + Icges(6,4) * t394 + t398;
t368 = Icges(6,3) * t301 - t237 - t398;
t401 = Icges(5,4) * t301;
t236 = -Icges(5,2) * t300 + Icges(5,6) * t394 + t401;
t369 = -Icges(5,1) * t300 - t236 - t401;
t294 = Icges(6,5) * t301;
t233 = Icges(6,6) * t394 + Icges(6,3) * t300 + t294;
t370 = -Icges(6,1) * t300 + t233 + t294;
t463 = (t260 + t261) * t394 + (t369 + t370) * t301 + (t367 + t368) * t300;
t460 = rSges(6,1) + pkin(4);
t333 = -Icges(6,5) * t411 - Icges(6,3) * t318;
t208 = Icges(6,6) * t297 + t296 * t333;
t337 = -Icges(6,1) * t411 - Icges(6,5) * t318;
t216 = Icges(6,4) * t297 + t296 * t337;
t335 = -Icges(6,4) * t411 - Icges(6,6) * t318;
t330 = Icges(6,2) * t297 - t165 * t318 - t173 * t411 + t296 * t335;
t55 = t297 * t169 + t273 * t208 + t274 * t216 + t296 * t330;
t209 = Icges(6,6) * t299 + t298 * t333;
t217 = Icges(6,4) * t299 + t298 * t337;
t329 = Icges(6,2) * t299 - t166 * t318 - t174 * t411 + t298 * t335;
t56 = t297 * t170 + t273 * t209 + t274 * t217 + t296 * t329;
t336 = -Icges(5,4) * t411 + Icges(5,2) * t318;
t214 = Icges(5,6) * t297 + t296 * t336;
t338 = -Icges(5,1) * t411 + Icges(5,4) * t318;
t218 = Icges(5,5) * t297 + t296 * t338;
t334 = -Icges(5,5) * t411 + Icges(5,6) * t318;
t328 = Icges(5,3) * t297 + t171 * t318 - t175 * t411 + t296 * t334;
t57 = t297 * t167 - t273 * t214 + t274 * t218 + t296 * t328;
t215 = Icges(5,6) * t299 + t298 * t336;
t219 = Icges(5,5) * t299 + t298 * t338;
t327 = Icges(5,3) * t299 + t172 * t318 - t176 * t411 + t298 * t334;
t58 = t297 * t168 - t273 * t215 + t274 * t219 + t296 * t327;
t235 = Icges(6,4) * t301 + Icges(6,2) * t394 + Icges(6,6) * t300;
t277 = (Icges(6,6) * t320 + t319 * t333) * t314;
t281 = (Icges(6,4) * t320 + t319 * t337) * t314;
t326 = -t233 * t318 - t237 * t411 + (Icges(6,2) * t320 + t319 * t335) * t314;
t86 = t297 * t235 + t273 * t277 + t274 * t281 + t296 * t326;
t234 = Icges(5,5) * t301 - Icges(5,6) * t300 + Icges(5,3) * t394;
t280 = (Icges(5,6) * t320 + t319 * t336) * t314;
t282 = (Icges(5,5) * t320 + t319 * t338) * t314;
t325 = t236 * t318 - t238 * t411 + (Icges(5,3) * t320 + t319 * t334) * t314;
t87 = t297 * t234 - t273 * t280 + t274 * t282 + t296 * t325;
t459 = (-t86 - t87) * t316 + ((t56 + t58) * t317 + (t55 + t57) * t315) * t314;
t59 = t299 * t169 + t275 * t208 + t276 * t216 + t298 * t330;
t60 = t299 * t170 + t275 * t209 + t276 * t217 + t298 * t329;
t61 = t299 * t167 - t275 * t214 + t276 * t218 + t298 * t328;
t62 = t299 * t168 - t275 * t215 + t276 * t219 + t298 * t327;
t88 = t299 * t235 + t275 * t277 + t276 * t281 + t298 * t326;
t89 = t299 * t234 - t275 * t280 + t276 * t282 + t298 * t325;
t458 = (-t89 - t88) * t316 + ((t60 + t62) * t317 + (t59 + t61) * t315) * t314;
t457 = (t468 * t320 + (t328 + t330) * t319) * t314 + (t216 + t218) * t301 + (t208 - t214) * t300;
t456 = (t467 * t320 + (t327 + t329) * t319) * t314 + (t217 + t219) * t301 + (t209 - t215) * t300;
t455 = rSges(6,3) + qJ(5);
t454 = t468 * t394 + (t173 + t175) * t301 + (t165 - t171) * t300;
t453 = t467 * t394 + (t174 + t176) * t301 + (t166 - t172) * t300;
t452 = ((t234 + t235) * t320 + (t325 + t326) * t319) * t314 + (t282 + t281) * t301 + (-t280 + t277) * t300;
t141 = t233 * t300 + t235 * t394 + t237 * t301;
t142 = t234 * t394 - t236 * t300 + t238 * t301;
t451 = t141 + t142;
t339 = -pkin(4) * t411 - qJ(5) * t318;
t340 = -rSges(6,1) * t411 - rSges(6,3) * t318;
t450 = t339 + t340;
t423 = t296 / 0.2e1;
t421 = t297 / 0.2e1;
t420 = t298 / 0.2e1;
t418 = t299 / 0.2e1;
t448 = t314 / 0.2e1;
t415 = -t316 / 0.2e1;
t446 = t394 / 0.2e1;
t445 = t296 * t466 + t298 * t465 + t394 * t463;
t444 = (t453 * t420 + t454 * t423 + t451 * t446) * t320;
t443 = Icges(4,6) * t316 + (t319 * t464 + t320 * t462) * t314;
t442 = Icges(4,5) * t316 + (t319 * t461 + t320 * t464) * t314;
t396 = t314 * t317;
t397 = t314 * t315;
t441 = (-Icges(4,6) * t396 + t298 * t464 + t299 * t462) * t317 + (-Icges(4,6) * t397 + t296 * t464 + t297 * t462) * t315;
t440 = (-Icges(4,5) * t396 + t298 * t461 + t299 * t464) * t317 + (-Icges(4,5) * t397 + t296 * t461 + t297 * t464) * t315;
t439 = 2 * qJD(3);
t438 = 4 * qJD(3);
t437 = 2 * qJD(4);
t436 = 4 * qJD(4);
t435 = m(4) / 0.2e1;
t434 = m(5) / 0.2e1;
t433 = m(5) / 0.4e1;
t432 = m(6) / 0.2e1;
t431 = m(6) / 0.4e1;
t366 = rSges(6,2) * t394 + t455 * t300 + t460 * t301;
t378 = rSges(6,2) * t296 + t455 * t273 + t460 * t274;
t115 = t296 * t366 - t378 * t394;
t377 = rSges(6,2) * t298 + t455 * t275 + t460 * t276;
t116 = -t298 * t366 + t377 * t394;
t372 = t299 * rSges(6,2) + t450 * t298;
t373 = t297 * rSges(6,2) + t450 * t296;
t54 = -t296 * t372 - t297 * t377 + t298 * t373 + t299 * t378;
t362 = (rSges(6,2) * t320 + t319 * t340) * t314 + t339 * t394;
t76 = t366 * t297 + t362 * t296 + (-t319 * t373 - t320 * t378) * t314;
t77 = -t366 * t299 - t362 * t298 + (t319 * t372 + t320 * t377) * t314;
t84 = -t296 * t377 + t298 * t378;
t430 = m(6) * (t273 * t77 + t275 * t76 + t300 * t54 + (-t115 * t298 - t116 * t296 - t394 * t84) * t318);
t429 = m(6) * (t115 * t76 + t116 * t77 + t54 * t84);
t178 = rSges(5,1) * t274 - rSges(5,2) * t273 + rSges(5,3) * t296;
t341 = -rSges(5,1) * t411 + rSges(5,2) * t318;
t221 = t297 * rSges(5,3) + t296 * t341;
t240 = rSges(5,1) * t301 - rSges(5,2) * t300 + rSges(5,3) * t394;
t284 = (rSges(5,3) * t320 + t319 * t341) * t314;
t120 = t297 * t240 + t296 * t284 + (-t178 * t320 - t221 * t319) * t314;
t180 = rSges(5,1) * t276 - rSges(5,2) * t275 + rSges(5,3) * t298;
t223 = t299 * rSges(5,3) + t298 * t341;
t121 = -t299 * t240 - t298 * t284 + (t180 * t320 + t223 * t319) * t314;
t138 = t178 * t298 - t180 * t296;
t148 = -t178 * t394 + t240 * t296;
t149 = t180 * t394 - t240 * t298;
t90 = t178 * t299 - t180 * t297 + t221 * t298 - t223 * t296;
t428 = m(5) * (t120 * t148 + t121 * t149 + t138 * t90);
t256 = pkin(3) * t297 + pkin(6) * t296;
t242 = t256 * t396;
t258 = pkin(3) * t299 + pkin(6) * t298;
t376 = -t180 - t258;
t119 = t242 + (t178 * t317 + t315 * t376) * t314;
t307 = (pkin(3) * t320 + pkin(6) * t319) * t314;
t364 = t316 * t256 + t307 * t397;
t136 = t178 * t316 + t240 * t397 + t364;
t137 = t376 * t316 + (-t240 - t307) * t396;
t200 = -rSges(5,1) * t273 - rSges(5,2) * t274;
t204 = -rSges(5,1) * t275 - rSges(5,2) * t276;
t151 = (t200 * t317 - t204 * t315) * t314;
t267 = -rSges(5,1) * t300 - rSges(5,2) * t301;
t158 = t200 * t316 + t267 * t397;
t159 = -t204 * t316 - t267 * t396;
t427 = m(5) * (t119 * t151 + t136 * t158 + t137 * t159);
t105 = t316 * t378 + t366 * t397 + t364;
t351 = -t258 - t377;
t106 = t351 * t316 + (-t307 - t366) * t396;
t161 = t273 * t298 - t275 * t296;
t182 = (t273 * t317 - t275 * t315) * t314;
t184 = -t273 * t394 + t296 * t300;
t185 = t275 * t394 - t298 * t300;
t224 = t273 * t316 + t300 * t397;
t225 = -t275 * t316 - t300 * t396;
t83 = t242 + (t315 * t351 + t317 * t378) * t314;
t426 = m(6) * (t105 * t184 + t106 * t185 + t115 * t224 + t116 * t225 + t161 * t83 + t182 * t84);
t374 = -t460 * t275 + t455 * t276;
t375 = -t460 * t273 + t455 * t274;
t114 = (-t315 * t374 + t317 * t375) * t314;
t363 = -t460 * t300 + t455 * t301;
t129 = t316 * t375 + t363 * t397;
t130 = -t316 * t374 - t363 * t396;
t425 = m(6) * (t105 * t276 + t106 * t274 + t114 * t300 + t129 * t275 + t130 * t273 + t301 * t83);
t424 = m(6) * (t105 * t129 + t106 * t130 + t114 * t83);
t416 = t315 / 0.2e1;
t413 = t317 / 0.2e1;
t410 = m(6) * qJD(3);
t409 = m(6) * qJD(4);
t408 = m(6) * qJD(5);
t393 = t315 * t316;
t392 = t315 * t317;
t390 = t316 * t317;
t257 = -pkin(3) * t298 + pkin(6) * t299;
t371 = -t223 - t257;
t255 = -pkin(3) * t296 + pkin(6) * t297;
t306 = (-pkin(3) * t319 + pkin(6) * t320) * t314;
t365 = t316 * t255 + t306 * t397;
t361 = qJD(3) * t314;
t360 = qJD(3) * t316;
t133 = t234 * t296 - t236 * t273 + t238 * t274;
t97 = t167 * t296 - t171 * t273 + t175 * t274;
t98 = t168 * t296 - t172 * t273 + t176 * t274;
t10 = t57 * t296 + t97 * t297 + t58 * t298 + t98 * t299 + (t133 * t320 + t319 * t87) * t314;
t132 = t233 * t273 + t235 * t296 + t237 * t274;
t95 = t165 * t273 + t169 * t296 + t173 * t274;
t96 = t166 * t273 + t170 * t296 + t174 * t274;
t9 = t55 * t296 + t95 * t297 + t56 * t298 + t96 * t299 + (t132 * t320 + t319 * t86) * t314;
t359 = t9 / 0.2e1 + t10 / 0.2e1;
t100 = t166 * t275 + t170 * t298 + t174 * t276;
t134 = t233 * t275 + t235 * t298 + t237 * t276;
t99 = t165 * t275 + t169 * t298 + t173 * t276;
t11 = t100 * t299 + t59 * t296 + t99 * t297 + t60 * t298 + (t134 * t320 + t319 * t88) * t314;
t101 = t167 * t298 - t171 * t275 + t175 * t276;
t102 = t168 * t298 - t172 * t275 + t176 * t276;
t135 = t234 * t298 - t236 * t275 + t238 * t276;
t12 = t101 * t297 + t102 * t299 + t61 * t296 + t62 * t298 + (t135 * t320 + t319 * t89) * t314;
t358 = t11 / 0.2e1 + t12 / 0.2e1;
t357 = (t452 * t319 + t451 * t320) * t448 + t457 * t423 + t454 * t421 + t456 * t420 + t453 * t418;
t64 = t190 * t296 + t273 * t382 + t274 * t386;
t65 = t191 * t296 + t273 * t381 + t274 * t385;
t91 = t261 * t296 + t273 * t368 + t274 * t370;
t33 = -t316 * t91 + (t315 * t64 + t317 * t65) * t314;
t66 = t188 * t296 + t273 * t380 + t274 * t384;
t67 = t189 * t296 + t273 * t379 + t274 * t383;
t92 = t260 * t296 + t273 * t367 + t274 * t369;
t34 = -t316 * t92 + (t315 * t66 + t317 * t67) * t314;
t356 = t33 / 0.2e1 + t34 / 0.2e1;
t68 = t190 * t298 + t275 * t382 + t276 * t386;
t69 = t191 * t298 + t275 * t381 + t276 * t385;
t93 = t261 * t298 + t275 * t368 + t276 * t370;
t35 = -t316 * t93 + (t315 * t68 + t317 * t69) * t314;
t70 = t188 * t298 + t275 * t380 + t276 * t384;
t71 = t189 * t298 + t275 * t379 + t276 * t383;
t94 = t260 * t298 + t275 * t367 + t276 * t369;
t36 = -t316 * t94 + (t315 * t70 + t317 * t71) * t314;
t355 = t36 / 0.2e1 + t35 / 0.2e1;
t354 = (t457 * t315 + t456 * t317) * t448 + t452 * t415;
t353 = -(t315 * t466 + t317 * t465) * t314 / 0.2e1 + t463 * t316 / 0.2e1;
t350 = -t257 - t372;
t347 = -t318 * t394 / 0.2e1;
t321 = -m(5) * (t158 * t315 - t159 * t317) / 0.2e1 - m(6) * (t129 * t315 - t130 * t317) / 0.2e1;
t322 = (t120 * t315 - t121 * t317) * t434 + (t315 * t76 - t317 * t77) * t432;
t15 = t321 + t322;
t26 = 0.2e1 * (t54 / 0.4e1 - t114 / 0.4e1) * m(6) + 0.2e1 * (t90 / 0.4e1 - t151 / 0.4e1) * m(5);
t346 = t26 * qJD(1) + t15 * qJD(2);
t253 = -rSges(4,1) * t296 - rSges(4,2) * t297;
t305 = (-rSges(4,1) * t319 - rSges(4,2) * t320) * t314;
t206 = t253 * t316 + t305 * t397;
t254 = -rSges(4,1) * t298 - rSges(4,2) * t299;
t207 = -t254 * t316 - t305 * t396;
t345 = t315 * t206 - t317 * t207;
t324 = (t184 * t315 - t185 * t317) * t432;
t331 = m(6) * (-t274 * t317 + t276 * t315);
t128 = t324 - t331 / 0.2e1;
t156 = 0.2e1 * (t301 / 0.4e1 - t161 / 0.4e1) * m(6);
t344 = -t156 * qJD(1) + t128 * qJD(2);
t323 = (t296 * t317 - t298 * t315) * t318 * t432;
t332 = m(6) * (t224 * t315 - t225 * t317);
t139 = -t332 / 0.2e1 + t323;
t162 = (t347 - t182 / 0.2e1) * m(6);
t343 = -t162 * qJD(1) - t139 * qJD(2);
t342 = -t353 - t357;
t302 = (-Icges(4,5) * t319 - Icges(4,6) * t320) * t314;
t248 = -Icges(4,5) * t298 - Icges(4,6) * t299;
t247 = -Icges(4,5) * t296 - Icges(4,6) * t297;
t241 = t255 * t396;
t232 = rSges(4,1) * t299 - rSges(4,2) * t298 + rSges(4,3) * t396;
t231 = rSges(4,1) * t297 - rSges(4,2) * t296 + rSges(4,3) * t397;
t164 = (t253 * t317 - t254 * t315) * t314;
t163 = m(6) * t347 + t182 * t432;
t157 = (t161 + t301) * t432;
t155 = t204 * t394 - t267 * t298;
t154 = -t200 * t394 + t267 * t296;
t153 = (-t273 * t296 - t275 * t298 - t300 * t394) * t318;
t146 = t273 * t274 + t275 * t276 + t300 * t301;
t145 = t371 * t316 + (-t284 - t306) * t396;
t144 = t221 * t316 + t284 * t397 + t365;
t143 = t200 * t298 - t204 * t296;
t140 = t332 / 0.2e1 + t323;
t131 = t241 + (t221 * t317 + t315 * t371) * t314;
t127 = t324 + t331 / 0.2e1;
t125 = t350 * t316 + (-t306 - t362) * t396;
t124 = t316 * t373 + t362 * t397 + t365;
t123 = -t298 * t363 + t374 * t394;
t122 = t296 * t363 - t375 * t394;
t104 = t241 + (t315 * t350 + t317 * t373) * t314;
t103 = -t296 * t374 + t298 * t375;
t50 = t101 * t296 + t102 * t298 + t135 * t394;
t49 = t100 * t298 + t134 * t394 + t296 * t99;
t48 = t133 * t394 + t296 * t97 + t298 * t98;
t47 = t132 * t394 + t296 * t95 + t298 * t96;
t44 = t115 * t184 + t116 * t185 + t161 * t84;
t39 = t105 * t224 + t106 * t225 + t182 * t83;
t32 = t296 * t70 + t298 * t71 + t394 * t94;
t31 = t296 * t68 + t298 * t69 + t394 * t93;
t30 = t296 * t66 + t298 * t67 + t394 * t92;
t29 = t296 * t64 + t298 * t65 + t394 * t91;
t27 = (t151 + t90) * t434 + (t114 + t54) * t432;
t19 = t425 / 0.2e1;
t16 = -t321 + t322;
t13 = t426 / 0.2e1;
t6 = t430 / 0.2e1;
t5 = t19 + t6 - t426 / 0.2e1;
t4 = t13 + t19 - t430 / 0.2e1;
t3 = t13 + t6 - t425 / 0.2e1;
t2 = t353 * t316 + t427 + t424 + (t315 * t356 + t317 * t355) * t314;
t1 = (t49 / 0.2e1 + t50 / 0.2e1) * t299 + t358 * t298 + (t47 / 0.2e1 + t48 / 0.2e1) * t297 + t359 * t296 + t428 + t429 + (t319 * t357 + t444) * t314;
t7 = [0, 0, t27 * qJD(4) + t163 * qJD(5) + (t104 * t432 + t131 * t434 + t164 * t435) * t439, t27 * qJD(3) + (t103 * t432 + t143 * t434) * t437 + t157 * qJD(5), qJD(3) * t163 + qJD(4) * t157; 0, 0, t16 * qJD(4) + t140 * qJD(5) + (t345 * t435 + (t144 * t315 - t145 * t317) * t434 + (t124 * t315 - t125 * t317) * t432) * t439, t16 * qJD(3) + ((t154 * t315 - t155 * t317) * t434 + (t122 * t315 - t123 * t317) * t432) * t437 + t127 * qJD(5), qJD(3) * t140 + qJD(4) * t127; -qJD(4) * t26 - qJD(5) * t162, -qJD(4) * t15 - qJD(5) * t139, t2 * qJD(4) + (-t316 ^ 2 * t302 / 0.2e1 - t354) * t360 + (-(t298 * t442 + t299 * t443) * t390 / 0.2e1 - (t296 * t442 + t297 * t443) * t393 / 0.2e1 + (-t248 * t390 - t247 * t393 - (t319 * t442 + t320 * t443) * t316 + (t319 * t440 + t320 * t441) * t314) * t415 + ((t296 * t440 + t297 * t441 - t302 * t393 + (t247 * t315 ^ 2 + t248 * t392) * t314) * t314 + t459) * t416 + ((t298 * t440 + t299 * t441 - t302 * t390 + (t248 * t317 ^ 2 + t247 * t392) * t314) * t314 + t458) * t413) * t361 + ((t104 * t83 + t105 * t124 + t106 * t125) * t431 + (t119 * t131 + t136 * t144 + t137 * t145) * t433 + m(4) * ((t206 * t231 - t207 * t232) * t316 + ((t231 * t317 - t232 * t315) * t164 + t345 * (-t316 * rSges(4,3) + (rSges(4,1) * t320 - rSges(4,2) * t319) * t314)) * t314) / 0.4e1) * t438 + t39 * t408, t2 * qJD(3) + t4 * qJD(5) + (-t428 / 0.4e1 - t429 / 0.4e1) * t436 + ((t119 * t143 + t136 * t154 + t137 * t155 + t138 * t151 + t148 * t158 + t149 * t159) * t434 + (t103 * t83 + t105 * t122 + t106 * t123 + t114 * t84 + t115 * t129 + t116 * t130) * t432) * t437 - t346 + ((t355 - t358) * t298 + (t356 - t359) * t296 + (-t444 + (t31 / 0.2e1 + t32 / 0.2e1) * t317 + (t29 / 0.2e1 + t30 / 0.2e1) * t315 + t342 * t319) * t314 - (t47 + t48) * t297 / 0.2e1 - (t49 + t50) * t299 / 0.2e1 + t445 * t415) * qJD(4), t39 * t410 + t4 * qJD(4) + (t182 * t300 + t224 * t275 + t225 * t273 - t153) * t408 + t343; qJD(3) * t26 - qJD(5) * t156, qJD(3) * t15 + qJD(5) * t128, (t458 * t420 + t459 * t423) * qJD(3) + t1 * qJD(4) + t3 * qJD(5) + (-t427 / 0.4e1 - t424 / 0.4e1) * t438 + ((t119 * t90 + t120 * t136 + t121 * t137 + t131 * t138 + t144 * t148 + t145 * t149) * t434 + (t104 * t84 + t105 * t76 + t106 * t77 + t115 * t124 + t116 * t125 + t54 * t83) * t432) * t439 + ((-t134 / 0.2e1 - t135 / 0.2e1) * t299 + (-t132 / 0.2e1 - t133 / 0.2e1) * t297 + t342) * t360 + ((t454 * t315 + t453 * t317) * t320 * t448 + t354 * t319 + (-t142 / 0.2e1 - t141 / 0.2e1) * t389 + ((t96 + t98) * t317 + (t95 + t97) * t315) * t421 + ((t100 + t102) * t317 + (t101 + t99) * t315) * t418 - (t33 + t34) * t315 / 0.2e1 + (t9 + t10) * t416 - (t36 + t35) * t317 / 0.2e1 + (t11 + t12) * t413) * t361 + t346, t1 * qJD(3) + ((t29 + t30) * t423 + (t31 + t32) * t420 + t445 * t446) * qJD(4) + ((t84 * t103 + t115 * t122 + t116 * t123) * t431 + (t138 * t143 + t148 * t154 + t149 * t155) * t433) * t436 + t44 * t408, t3 * qJD(3) + t44 * t409 + (t161 * t300 + t184 * t275 + t185 * t273 - t146) * t408 + t344; qJD(3) * t162 + qJD(4) * t156, qJD(3) * t139 - qJD(4) * t128, (t104 * t300 + t124 * t275 + t125 * t273 + (-t105 * t298 - t106 * t296 - t394 * t83) * t318 - t39) * t410 + t5 * qJD(4) + t153 * t408 - t343, t5 * qJD(3) + (t103 * t300 + t115 * t276 + t116 * t274 + t122 * t275 + t123 * t273 + t301 * t84 - t44) * t409 + t146 * t408 - t344, 0.4e1 * (t153 * qJD(3) / 0.4e1 + t146 * qJD(4) / 0.4e1) * m(6);];
Cq = t7;
