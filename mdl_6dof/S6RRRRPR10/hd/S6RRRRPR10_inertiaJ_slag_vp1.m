% Calculate joint inertia matrix for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR10_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR10_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:01:01
% EndTime: 2019-03-09 23:01:22
% DurationCPUTime: 9.48s
% Computational Cost: add. (29181->651), mult. (49464->868), div. (0->0), fcn. (61905->12), ass. (0->318)
t334 = cos(pkin(6));
t337 = sin(qJ(2));
t333 = sin(pkin(6));
t405 = qJ(3) + qJ(4);
t366 = sin(t405);
t357 = t333 * t366;
t367 = cos(t405);
t302 = -t334 * t367 + t337 * t357;
t358 = t333 * t367;
t303 = t334 * t366 + t337 * t358;
t341 = cos(qJ(2));
t412 = t333 * t341;
t248 = -Icges(6,1) * t412 - Icges(6,4) * t303 + Icges(6,5) * t302;
t249 = Icges(5,5) * t303 - Icges(5,6) * t302 - Icges(5,3) * t412;
t456 = -t248 - t249;
t338 = sin(qJ(1));
t407 = t338 * t341;
t342 = cos(qJ(1));
t408 = t337 * t342;
t316 = t334 * t408 + t407;
t283 = t316 * t366 + t342 * t358;
t284 = t316 * t367 - t342 * t357;
t406 = t341 * t342;
t409 = t337 * t338;
t315 = -t334 * t406 + t409;
t188 = Icges(6,5) * t315 - Icges(6,6) * t284 + Icges(6,3) * t283;
t194 = Icges(5,4) * t284 - Icges(5,2) * t283 + Icges(5,6) * t315;
t455 = t188 - t194;
t318 = -t334 * t409 + t406;
t285 = t318 * t366 - t338 * t358;
t286 = t318 * t367 + t338 * t357;
t317 = t334 * t407 + t408;
t189 = Icges(6,5) * t317 - Icges(6,6) * t286 + Icges(6,3) * t285;
t195 = Icges(5,4) * t286 - Icges(5,2) * t285 + Icges(5,6) * t317;
t454 = t189 - t195;
t190 = Icges(5,5) * t284 - Icges(5,6) * t283 + Icges(5,3) * t315;
t196 = Icges(6,1) * t315 - Icges(6,4) * t284 + Icges(6,5) * t283;
t453 = t190 + t196;
t191 = Icges(5,5) * t286 - Icges(5,6) * t285 + Icges(5,3) * t317;
t197 = Icges(6,1) * t317 - Icges(6,4) * t286 + Icges(6,5) * t285;
t452 = t191 + t197;
t192 = Icges(6,4) * t315 - Icges(6,2) * t284 + Icges(6,6) * t283;
t198 = Icges(5,1) * t284 - Icges(5,4) * t283 + Icges(5,5) * t315;
t451 = -t192 + t198;
t193 = Icges(6,4) * t317 - Icges(6,2) * t286 + Icges(6,6) * t285;
t199 = Icges(5,1) * t286 - Icges(5,4) * t285 + Icges(5,5) * t317;
t450 = -t193 + t199;
t246 = -Icges(6,5) * t412 - Icges(6,6) * t303 + Icges(6,3) * t302;
t247 = -Icges(6,4) * t412 - Icges(6,2) * t303 + Icges(6,6) * t302;
t250 = Icges(5,4) * t303 - Icges(5,2) * t302 - Icges(5,6) * t412;
t251 = Icges(5,1) * t303 - Icges(5,4) * t302 - Icges(5,5) * t412;
t449 = (-t247 + t251) * t303 + (t246 - t250) * t302;
t448 = t412 * t456 + t449;
t447 = t283 * t455 + t451 * t284 + t453 * t315;
t446 = t283 * t454 + t284 * t450 + t315 * t452;
t445 = t285 * t455 + t451 * t286 + t453 * t317;
t444 = t285 * t454 + t286 * t450 + t317 * t452;
t129 = t246 * t283 - t247 * t284 + t248 * t315;
t131 = t249 * t315 - t250 * t283 + t251 * t284;
t443 = t129 + t131;
t130 = t246 * t285 - t247 * t286 + t248 * t317;
t132 = t249 * t317 - t250 * t285 + t251 * t286;
t442 = t130 + t132;
t440 = t448 * t334;
t335 = sin(qJ(6));
t339 = cos(qJ(6));
t238 = t283 * t339 - t315 * t335;
t239 = t283 * t335 + t315 * t339;
t354 = -t239 * rSges(7,1) - t238 * rSges(7,2);
t157 = t284 * rSges(7,3) - t354;
t404 = t315 * pkin(5) + t284 * pkin(11) + t157;
t281 = t302 * t339 + t335 * t412;
t282 = t302 * t335 - t339 * t412;
t185 = rSges(7,1) * t282 + rSges(7,2) * t281 + rSges(7,3) * t303;
t439 = -pkin(5) * t412 + pkin(11) * t303 + t185;
t151 = Icges(7,5) * t239 + Icges(7,6) * t238 + Icges(7,3) * t284;
t153 = Icges(7,4) * t239 + Icges(7,2) * t238 + Icges(7,6) * t284;
t155 = Icges(7,1) * t239 + Icges(7,4) * t238 + Icges(7,5) * t284;
t66 = t151 * t284 + t153 * t238 + t155 * t239;
t240 = t285 * t339 - t317 * t335;
t241 = t285 * t335 + t317 * t339;
t152 = Icges(7,5) * t241 + Icges(7,6) * t240 + Icges(7,3) * t286;
t154 = Icges(7,4) * t241 + Icges(7,2) * t240 + Icges(7,6) * t286;
t156 = Icges(7,1) * t241 + Icges(7,4) * t240 + Icges(7,5) * t286;
t67 = t152 * t284 + t154 * t238 + t156 * t239;
t182 = Icges(7,5) * t282 + Icges(7,6) * t281 + Icges(7,3) * t303;
t183 = Icges(7,4) * t282 + Icges(7,2) * t281 + Icges(7,6) * t303;
t184 = Icges(7,1) * t282 + Icges(7,4) * t281 + Icges(7,5) * t303;
t84 = t182 * t284 + t183 * t238 + t184 * t239;
t11 = t315 * t66 + t317 * t67 - t412 * t84;
t438 = t315 * t447 + t446 * t317 - t443 * t412 + t11;
t68 = t151 * t286 + t153 * t240 + t155 * t241;
t69 = t152 * t286 + t154 * t240 + t156 * t241;
t85 = t182 * t286 + t183 * t240 + t184 * t241;
t12 = t315 * t68 + t317 * t69 - t412 * t85;
t437 = t315 * t445 + t317 * t444 - t412 * t442 + t12;
t15 = t84 * t334 + (t338 * t67 - t342 * t66) * t333;
t436 = t15 + t443 * t334 + (t446 * t338 - t342 * t447) * t333;
t16 = t85 * t334 + (t338 * t69 - t342 * t68) * t333;
t435 = t16 + t442 * t334 + (t338 * t444 - t342 * t445) * t333;
t74 = t152 * t303 + t154 * t281 + t156 * t282;
t421 = t74 * t317;
t73 = t151 * t303 + t153 * t281 + t155 * t282;
t422 = t73 * t315;
t91 = t303 * t182 + t281 * t183 + t282 * t184;
t21 = -t412 * t91 + t421 + t422;
t117 = -t191 * t412 - t195 * t302 + t199 * t303;
t417 = t117 * t317;
t116 = -t190 * t412 - t194 * t302 + t198 * t303;
t418 = t116 * t315;
t115 = t189 * t302 - t193 * t303 - t197 * t412;
t419 = t115 * t317;
t114 = t188 * t302 - t192 * t303 - t196 * t412;
t420 = t114 * t315;
t434 = -t412 * t448 + t21 + t417 + t418 + t419 + t420;
t89 = t91 * t334;
t23 = t89 + (t74 * t338 - t73 * t342) * t333;
t433 = t23 + t440 + ((-t114 - t116) * t342 + (t115 + t117) * t338) * t333;
t411 = t333 * t342;
t432 = t284 / 0.2e1;
t431 = t286 / 0.2e1;
t430 = t303 / 0.2e1;
t429 = t315 / 0.2e1;
t428 = t317 / 0.2e1;
t427 = t334 / 0.2e1;
t426 = t338 / 0.2e1;
t425 = -t342 / 0.2e1;
t340 = cos(qJ(3));
t331 = pkin(3) * t340 + pkin(2);
t424 = -pkin(2) + t331;
t423 = t283 * rSges(6,3);
t260 = Icges(3,5) * t316 - Icges(3,6) * t315 - Icges(3,3) * t411;
t416 = t260 * t342;
t415 = t333 * t337;
t414 = t333 * t338;
t413 = t333 * t340;
t336 = sin(qJ(3));
t410 = t334 * t336;
t158 = t241 * rSges(7,1) + t240 * rSges(7,2) + t286 * rSges(7,3);
t403 = t317 * pkin(5) + pkin(11) * t286 + t158;
t200 = t315 * rSges(6,1) - t284 * rSges(6,2) + t423;
t274 = t283 * qJ(5);
t226 = t284 * pkin(4) + t274;
t206 = t317 * t226;
t402 = t317 * t200 + t206;
t310 = t315 * pkin(9);
t382 = t336 * t411;
t323 = pkin(3) * t382;
t343 = -pkin(10) - pkin(9);
t204 = -t315 * t343 + t316 * t424 - t310 - t323;
t266 = pkin(3) * t410 + ((pkin(9) + t343) * t341 + t424 * t337) * t333;
t401 = t204 * t412 + t315 * t266;
t273 = t318 * pkin(2) + pkin(9) * t317;
t383 = t336 * t414;
t373 = pkin(3) * t383 - t317 * t343 + t318 * t331;
t205 = -t273 + t373;
t271 = t334 * t273;
t400 = t334 * t205 + t271;
t201 = t317 * rSges(6,1) - t286 * rSges(6,2) + t285 * rSges(6,3);
t227 = t286 * pkin(4) + qJ(5) * t285;
t399 = -t201 - t227;
t203 = t286 * rSges(5,1) - t285 * rSges(5,2) + t317 * rSges(5,3);
t398 = -t203 - t205;
t272 = t316 * pkin(2) + t310;
t397 = -t204 - t272;
t290 = -t316 * t336 - t340 * t411;
t291 = t316 * t340 - t382;
t214 = rSges(4,1) * t291 + rSges(4,2) * t290 + rSges(4,3) * t315;
t396 = -t214 - t272;
t258 = pkin(4) * t303 + qJ(5) * t302;
t395 = t226 * t412 + t315 * t258;
t355 = -t284 * rSges(5,1) + t283 * rSges(5,2);
t202 = t315 * rSges(5,3) - t355;
t253 = rSges(5,1) * t303 - rSges(5,2) * t302 - rSges(5,3) * t412;
t148 = t202 * t412 + t315 * t253;
t313 = t334 * t340 - t336 * t415;
t314 = t337 * t413 + t410;
t256 = Icges(4,4) * t314 + Icges(4,2) * t313 - Icges(4,6) * t412;
t257 = Icges(4,1) * t314 + Icges(4,4) * t313 - Icges(4,5) * t412;
t392 = t313 * t256 + t314 * t257;
t252 = -rSges(6,1) * t412 - rSges(6,2) * t303 + rSges(6,3) * t302;
t391 = -t252 - t258;
t390 = -t253 - t266;
t389 = t272 * t414 + t273 * t411;
t388 = t342 * pkin(1) + pkin(8) * t414;
t386 = t73 / 0.2e1 + t84 / 0.2e1;
t385 = t85 / 0.2e1 + t74 / 0.2e1;
t384 = -t91 - t448;
t381 = t317 * t404 + t206;
t380 = -t227 - t403;
t379 = -t258 - t439;
t378 = t334 * t227 + t400;
t377 = -t205 + t399;
t376 = -t226 + t397;
t375 = -t266 + t391;
t292 = -t318 * t336 + t338 * t413;
t293 = t318 * t340 + t383;
t215 = t293 * rSges(4,1) + t292 * rSges(4,2) + t317 * rSges(4,3);
t298 = Icges(3,3) * t334 + (Icges(3,5) * t337 + Icges(3,6) * t341) * t333;
t299 = Icges(3,6) * t334 + (Icges(3,4) * t337 + Icges(3,2) * t341) * t333;
t300 = Icges(3,5) * t334 + (Icges(3,1) * t337 + Icges(3,4) * t341) * t333;
t374 = t334 * t298 + t299 * t412 + t300 * t415;
t268 = t318 * rSges(3,1) - t317 * rSges(3,2) + rSges(3,3) * t414;
t371 = -t412 / 0.2e1;
t208 = Icges(4,5) * t291 + Icges(4,6) * t290 + Icges(4,3) * t315;
t210 = Icges(4,4) * t291 + Icges(4,2) * t290 + Icges(4,6) * t315;
t212 = Icges(4,1) * t291 + Icges(4,4) * t290 + Icges(4,5) * t315;
t122 = -t208 * t412 + t210 * t313 + t212 * t314;
t255 = Icges(4,5) * t314 + Icges(4,6) * t313 - Icges(4,3) * t412;
t135 = t255 * t315 + t256 * t290 + t257 * t291;
t369 = t135 / 0.2e1 + t122 / 0.2e1;
t209 = Icges(4,5) * t293 + Icges(4,6) * t292 + Icges(4,3) * t317;
t211 = Icges(4,4) * t293 + Icges(4,2) * t292 + Icges(4,6) * t317;
t213 = Icges(4,1) * t293 + Icges(4,4) * t292 + Icges(4,5) * t317;
t123 = -t209 * t412 + t211 * t313 + t213 * t314;
t136 = t255 * t317 + t256 * t292 + t257 * t293;
t368 = t136 / 0.2e1 + t123 / 0.2e1;
t365 = -t338 * pkin(1) + pkin(8) * t411;
t259 = rSges(4,1) * t314 + rSges(4,2) * t313 - rSges(4,3) * t412;
t319 = (pkin(2) * t337 - pkin(9) * t341) * t333;
t364 = t333 * (-t259 - t319);
t363 = -t205 + t380;
t362 = t204 * t414 + t205 * t411 + t389;
t361 = -t266 + t379;
t118 = t200 * t412 + t315 * t252 + t395;
t360 = t333 * (-t319 + t390);
t87 = t91 * t303;
t18 = t73 * t284 + t74 * t286 + t87;
t3 = t284 * t66 + t286 * t67 + t303 * t84;
t4 = t284 * t68 + t286 * t69 + t303 * t85;
t359 = t11 * t432 + t12 * t431 + t18 * t371 + t21 * t430 + t3 * t429 + t4 * t428;
t356 = t315 * t438 + t317 * t437;
t353 = t373 + t388;
t352 = t333 * (-t319 + t375);
t351 = t226 * t414 + t227 * t411 + t362;
t75 = t439 * t315 + t404 * t412 + t395;
t350 = t333 * (-t319 + t361);
t349 = -t316 * t331 + t323 + t365;
t348 = -t274 + t349;
t267 = rSges(3,1) * t316 - rSges(3,2) * t315 - rSges(3,3) * t411;
t347 = t227 + t353;
t346 = -t412 * t434 + t356;
t345 = t420 / 0.2e1 + t419 / 0.2e1 + t418 / 0.2e1 + t417 / 0.2e1 + t422 / 0.2e1 + t421 / 0.2e1 + (t84 + t443) * t429 + (t85 + t442) * t428;
t344 = t436 * t429 + t435 * t428 + t434 * t427 + t437 * t414 / 0.2e1 + t433 * t371 - t438 * t411 / 0.2e1;
t325 = rSges(2,1) * t342 - rSges(2,2) * t338;
t324 = -rSges(2,1) * t338 - rSges(2,2) * t342;
t301 = t334 * rSges(3,3) + (rSges(3,1) * t337 + rSges(3,2) * t341) * t333;
t265 = Icges(3,1) * t318 - Icges(3,4) * t317 + Icges(3,5) * t414;
t264 = Icges(3,1) * t316 - Icges(3,4) * t315 - Icges(3,5) * t411;
t263 = Icges(3,4) * t318 - Icges(3,2) * t317 + Icges(3,6) * t414;
t262 = Icges(3,4) * t316 - Icges(3,2) * t315 - Icges(3,6) * t411;
t261 = Icges(3,5) * t318 - Icges(3,6) * t317 + Icges(3,3) * t414;
t243 = t268 + t388;
t242 = -t267 + t365;
t220 = -t267 * t334 - t301 * t411;
t219 = t268 * t334 - t301 * t414;
t207 = t374 * t334;
t176 = (t267 * t338 + t268 * t342) * t333;
t175 = t298 * t414 - t299 * t317 + t300 * t318;
t174 = -t298 * t411 - t299 * t315 + t300 * t316;
t173 = t317 * t204;
t172 = t317 * t202;
t168 = t273 + t215 + t388;
t167 = t365 + t396;
t164 = -t215 * t412 - t259 * t317;
t163 = t214 * t412 + t259 * t315;
t162 = t334 * t261 + (t263 * t341 + t265 * t337) * t333;
t161 = t334 * t260 + (t262 * t341 + t264 * t337) * t333;
t160 = t353 + t203;
t159 = (-rSges(5,3) + t343) * t315 + t349 + t355;
t149 = -t203 * t412 - t317 * t253;
t146 = -t255 * t412 + t392;
t145 = t146 * t334;
t144 = t214 * t317 - t215 * t315;
t143 = t334 * t396 + t342 * t364;
t142 = t334 * t215 + t338 * t364 + t271;
t139 = -t203 * t315 + t172;
t134 = t347 + t201;
t133 = -t423 + (-rSges(6,1) + t343) * t315 + (rSges(6,2) - pkin(4)) * t284 + t348;
t128 = (t214 * t338 + t215 * t342) * t333 + t389;
t121 = t158 * t303 - t185 * t286;
t120 = -t157 * t303 + t185 * t284;
t119 = t317 * t391 + t399 * t412;
t113 = t317 * t390 + t398 * t412;
t112 = t148 + t401;
t111 = t209 * t317 + t211 * t292 + t213 * t293;
t110 = t208 * t317 + t210 * t292 + t212 * t293;
t109 = t209 * t315 + t211 * t290 + t213 * t291;
t108 = t208 * t315 + t210 * t290 + t212 * t291;
t103 = t347 + t403;
t102 = (-pkin(5) + t343) * t315 + (-rSges(7,3) - pkin(4) - pkin(11)) * t284 + t348 + t354;
t93 = (-t202 + t397) * t334 + t342 * t360;
t92 = t334 * t203 + t338 * t360 + t400;
t90 = t157 * t286 - t158 * t284;
t88 = t315 * t399 + t402;
t86 = t315 * t398 + t172 + t173;
t81 = t317 * t375 + t377 * t412;
t80 = t118 + t401;
t79 = (t202 * t338 + t203 * t342) * t333 + t362;
t78 = (-t200 + t376) * t334 + t342 * t352;
t77 = t334 * t201 + t338 * t352 + t378;
t76 = t317 * t379 + t380 * t412;
t70 = t315 * t377 + t173 + t402;
t65 = (t200 * t338 + t201 * t342) * t333 + t351;
t64 = t315 * t380 + t381;
t63 = t317 * t361 + t363 * t412;
t62 = t75 + t401;
t61 = (t376 - t404) * t334 + t342 * t350;
t60 = t334 * t403 + t338 * t350 + t378;
t59 = t145 + (-t122 * t342 + t123 * t338) * t333;
t58 = t122 * t315 + t123 * t317 - t146 * t412;
t55 = t315 * t363 + t173 + t381;
t52 = (t338 * t404 + t342 * t403) * t333 + t351;
t51 = t136 * t334 + (-t110 * t342 + t111 * t338) * t333;
t50 = t135 * t334 + (-t108 * t342 + t109 * t338) * t333;
t45 = t110 * t315 + t111 * t317 - t136 * t412;
t44 = t108 * t315 + t109 * t317 - t135 * t412;
t1 = [Icges(2,3) + (-t255 + t456) * t412 + m(7) * (t102 ^ 2 + t103 ^ 2) + m(5) * (t159 ^ 2 + t160 ^ 2) + m(6) * (t133 ^ 2 + t134 ^ 2) + m(4) * (t167 ^ 2 + t168 ^ 2) + m(3) * (t242 ^ 2 + t243 ^ 2) + m(2) * (t324 ^ 2 + t325 ^ 2) + t374 + t91 + t392 + t449; t89 + t207 + t145 + m(7) * (t102 * t61 + t103 * t60) + m(6) * (t133 * t78 + t134 * t77) + m(5) * (t159 * t93 + t160 * t92) + m(4) * (t142 * t168 + t143 * t167) + m(3) * (t219 * t243 + t220 * t242) + ((-t161 / 0.2e1 - t129 / 0.2e1 - t131 / 0.2e1 - t174 / 0.2e1 - t116 / 0.2e1 - t114 / 0.2e1 - t369 - t386) * t342 + (t162 / 0.2e1 + t130 / 0.2e1 + t132 / 0.2e1 + t175 / 0.2e1 + t117 / 0.2e1 + t115 / 0.2e1 + t368 + t385) * t338) * t333 + t440; (t59 + t207 + t433) * t334 + m(7) * (t52 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(6) * (t65 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(5) * (t79 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(4) * (t128 ^ 2 + t142 ^ 2 + t143 ^ 2) + m(3) * (t176 ^ 2 + t219 ^ 2 + t220 ^ 2) + ((-t50 + (-t262 * t315 + t264 * t316 - t333 * t416) * t411 - t436) * t342 + (t51 + ((-t263 * t317 + t265 * t318 + (t261 * t338 - t416) * t333) * t338 + (t261 * t411 + t262 * t317 + t263 * t315 - t264 * t318 - t265 * t316) * t342) * t333 + t435) * t338 + ((-t161 - t174) * t342 + (t162 + t175) * t338) * t334) * t333; t345 + t368 * t317 + t369 * t315 + m(7) * (t102 * t62 + t103 * t63) + m(5) * (t112 * t159 + t113 * t160) + m(6) * (t133 * t80 + t134 * t81) + m(4) * (t163 * t167 + t164 * t168) + (-t146 + t384) * t412; t344 + m(7) * (t52 * t55 + t60 * t63 + t61 * t62) + m(6) * (t65 * t70 + t77 * t81 + t78 * t80) + m(5) * (t112 * t93 + t113 * t92 + t79 * t86) + m(4) * (t128 * t144 + t142 * t164 + t143 * t163) + (t44 * t425 - t341 * t59 / 0.2e1 + t45 * t426) * t333 + t50 * t429 + t51 * t428 + t58 * t427; t315 * t44 + t317 * t45 + (-t58 - t434) * t412 + m(7) * (t55 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(6) * (t70 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(5) * (t112 ^ 2 + t113 ^ 2 + t86 ^ 2) + m(4) * (t144 ^ 2 + t163 ^ 2 + t164 ^ 2) + t356; t345 + t384 * t412 + m(7) * (t102 * t75 + t103 * t76) + m(5) * (t148 * t159 + t149 * t160) + m(6) * (t118 * t133 + t119 * t134); t344 + m(7) * (t64 * t52 + t60 * t76 + t61 * t75) + m(6) * (t118 * t78 + t119 * t77 + t65 * t88) + m(5) * (t139 * t79 + t148 * t93 + t149 * t92); m(7) * (t64 * t55 + t62 * t75 + t63 * t76) + m(6) * (t118 * t80 + t119 * t81 + t70 * t88) + m(5) * (t112 * t148 + t113 * t149 + t139 * t86) + t346; m(7) * (t64 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(6) * (t118 ^ 2 + t119 ^ 2 + t88 ^ 2) + m(5) * (t139 ^ 2 + t148 ^ 2 + t149 ^ 2) + t346; m(7) * (t102 * t285 + t103 * t283) + m(6) * (t133 * t285 + t134 * t283); m(7) * (t283 * t60 + t285 * t61 + t302 * t52) + m(6) * (t283 * t77 + t285 * t78 + t302 * t65); m(7) * (t283 * t63 + t285 * t62 + t302 * t55) + m(6) * (t283 * t81 + t285 * t80 + t302 * t70); m(7) * (t283 * t76 + t285 * t75 + t302 * t64) + m(6) * (t118 * t285 + t119 * t283 + t302 * t88); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t283 ^ 2 + t285 ^ 2 + t302 ^ 2); t87 + m(7) * (t102 * t120 + t103 * t121) + t385 * t286 + t386 * t284; m(7) * (t120 * t61 + t121 * t60 + t52 * t90) + t18 * t427 + t23 * t430 + t15 * t432 + t16 * t431 + (t3 * t425 + t4 * t426) * t333; m(7) * (t120 * t62 + t121 * t63 + t55 * t90) + t359; m(7) * (t120 * t75 + t121 * t76 + t64 * t90) + t359; m(7) * (t120 * t285 + t121 * t283 + t302 * t90); t286 * t4 + t284 * t3 + t303 * t18 + m(7) * (t120 ^ 2 + t121 ^ 2 + t90 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
