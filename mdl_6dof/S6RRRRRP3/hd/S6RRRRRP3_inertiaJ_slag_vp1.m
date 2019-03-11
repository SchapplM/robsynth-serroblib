% Calculate joint inertia matrix for
% S6RRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:05:53
% EndTime: 2019-03-10 01:06:10
% DurationCPUTime: 7.41s
% Computational Cost: add. (18406->517), mult. (17589->706), div. (0->0), fcn. (19032->10), ass. (0->263)
t415 = Icges(6,4) + Icges(7,4);
t414 = -Icges(7,5) - Icges(6,5);
t413 = -Icges(7,6) - Icges(6,6);
t280 = qJ(4) + qJ(5);
t269 = sin(t280);
t271 = cos(t280);
t287 = cos(qJ(1));
t281 = qJ(2) + qJ(3);
t272 = cos(t281);
t284 = sin(qJ(1));
t355 = t272 * t284;
t221 = -t269 * t355 - t271 * t287;
t222 = -t269 * t287 + t271 * t355;
t270 = sin(t281);
t358 = t270 * t284;
t141 = Icges(7,5) * t222 + Icges(7,6) * t221 + Icges(7,3) * t358;
t143 = Icges(6,5) * t222 + Icges(6,6) * t221 + Icges(6,3) * t358;
t412 = t141 + t143;
t354 = t272 * t287;
t223 = -t269 * t354 + t271 * t284;
t224 = t269 * t284 + t271 * t354;
t357 = t270 * t287;
t142 = Icges(7,5) * t224 + Icges(7,6) * t223 + Icges(7,3) * t357;
t144 = Icges(6,5) * t224 + Icges(6,6) * t223 + Icges(6,3) * t357;
t411 = t142 + t144;
t145 = Icges(7,4) * t222 + Icges(7,2) * t221 + Icges(7,6) * t358;
t147 = Icges(6,4) * t222 + Icges(6,2) * t221 + Icges(6,6) * t358;
t410 = t145 + t147;
t146 = Icges(7,4) * t224 + Icges(7,2) * t223 + Icges(7,6) * t357;
t148 = Icges(6,4) * t224 + Icges(6,2) * t223 + Icges(6,6) * t357;
t409 = t146 + t148;
t149 = Icges(7,1) * t222 + Icges(7,4) * t221 + Icges(7,5) * t358;
t151 = Icges(6,1) * t222 + Icges(6,4) * t221 + Icges(6,5) * t358;
t408 = t149 + t151;
t150 = Icges(7,1) * t224 + Icges(7,4) * t223 + Icges(7,5) * t357;
t152 = Icges(6,1) * t224 + Icges(6,4) * t223 + Icges(6,5) * t357;
t407 = t150 + t152;
t406 = (-Icges(7,3) - Icges(6,3)) * t272 + (t269 * t413 - t271 * t414) * t270;
t405 = t413 * t272 + (t415 * t271 + (-Icges(6,2) - Icges(7,2)) * t269) * t270;
t404 = t414 * t272 + ((Icges(6,1) + Icges(7,1)) * t271 - t415 * t269) * t270;
t403 = t410 * t221 + t408 * t222 + t358 * t412;
t402 = t221 * t409 + t222 * t407 + t358 * t411;
t401 = t410 * t223 + t408 * t224 + t357 * t412;
t400 = t223 * t409 + t224 * t407 + t357 * t411;
t388 = t221 * t405 + t222 * t404 + t358 * t406;
t387 = t223 * t405 + t224 * t404 + t357 * t406;
t395 = t404 * t270 * t271;
t397 = t405 * t269;
t386 = -t270 * t397 - t272 * t406 + t395;
t399 = t386 * t272;
t285 = cos(qJ(4));
t266 = t285 * pkin(4) + pkin(3);
t246 = pkin(5) * t271 + t266;
t282 = sin(qJ(4));
t247 = pkin(4) * t282 + pkin(5) * t269;
t396 = t224 * rSges(7,1) + t223 * rSges(7,2) + rSges(7,3) * t357 + t246 * t354 + t284 * t247;
t394 = -t388 * t272 + (t284 * t403 + t402 * t287) * t270;
t393 = -t387 * t272 + (t284 * t401 + t287 * t400) * t270;
t392 = t402 * t284 - t287 * t403;
t391 = t284 * t400 - t287 * t401;
t65 = -t141 * t272 + (-t145 * t269 + t149 * t271) * t270;
t67 = -t143 * t272 + (-t147 * t269 + t151 * t271) * t270;
t390 = -t65 - t67;
t66 = -t142 * t272 + (-t146 * t269 + t150 * t271) * t270;
t68 = -t144 * t272 + (-t148 * t269 + t152 * t271) * t270;
t389 = -t66 - t68;
t288 = -pkin(10) - pkin(9);
t277 = -qJ(6) + t288;
t311 = -t222 * rSges(7,1) - t221 * rSges(7,2);
t351 = t282 * t287;
t356 = t270 * t288;
t334 = pkin(4) * t351 + t284 * t356;
t335 = t246 - t266;
t385 = -t287 * t247 + (-t270 * t277 + t272 * t335) * t284 + t334 + rSges(7,3) * t358 - t311;
t331 = t277 - t288;
t352 = t282 * t284;
t336 = -pkin(4) * t352 - t266 * t354;
t384 = -t331 * t357 + t336 + t396;
t383 = (t331 - rSges(7,3)) * t272 + (rSges(7,1) * t271 - rSges(7,2) * t269 + t335) * t270;
t305 = Icges(4,5) * t272 - Icges(4,6) * t270;
t209 = -Icges(4,3) * t287 + t284 * t305;
t210 = Icges(4,3) * t284 + t287 * t305;
t279 = t287 ^ 2;
t361 = Icges(4,4) * t272;
t307 = -Icges(4,2) * t270 + t361;
t212 = Icges(4,6) * t284 + t287 * t307;
t362 = Icges(4,4) * t270;
t309 = Icges(4,1) * t272 - t362;
t214 = Icges(4,5) * t284 + t287 * t309;
t303 = -t212 * t270 + t214 * t272;
t211 = -Icges(4,6) * t287 + t284 * t307;
t213 = -Icges(4,5) * t287 + t284 * t309;
t304 = t211 * t270 - t213 * t272;
t349 = t285 * t287;
t234 = -t272 * t352 - t349;
t350 = t284 * t285;
t235 = t272 * t350 - t351;
t168 = Icges(5,5) * t235 + Icges(5,6) * t234 + Icges(5,3) * t358;
t170 = Icges(5,4) * t235 + Icges(5,2) * t234 + Icges(5,6) * t358;
t172 = Icges(5,1) * t235 + Icges(5,4) * t234 + Icges(5,5) * t358;
t72 = t168 * t358 + t170 * t234 + t172 * t235;
t236 = -t272 * t351 + t350;
t237 = t272 * t349 + t352;
t169 = Icges(5,5) * t237 + Icges(5,6) * t236 + Icges(5,3) * t357;
t171 = Icges(5,4) * t237 + Icges(5,2) * t236 + Icges(5,6) * t357;
t173 = Icges(5,1) * t237 + Icges(5,4) * t236 + Icges(5,5) * t357;
t73 = t169 * t358 + t171 * t234 + t173 * t235;
t37 = t284 * t73 - t287 * t72;
t382 = -t37 - t279 * t209 - (t303 * t284 + (-t210 + t304) * t287) * t284 - t392;
t278 = t284 ^ 2;
t381 = -t272 / 0.2e1;
t380 = t284 / 0.2e1;
t379 = -t287 / 0.2e1;
t283 = sin(qJ(2));
t378 = pkin(2) * t283;
t377 = pkin(3) * t272;
t376 = pkin(9) * t270;
t375 = -pkin(3) + t266;
t374 = t399 + (t284 * t390 + t287 * t389) * t270;
t286 = cos(qJ(2));
t373 = rSges(3,1) * t286;
t372 = rSges(3,2) * t283;
t371 = t287 * rSges(3,3);
t370 = t65 * t287;
t369 = t66 * t284;
t368 = t67 * t287;
t367 = t68 * t284;
t79 = -t272 * t168 + (-t170 * t282 + t172 * t285) * t270;
t366 = t79 * t287;
t80 = -t272 * t169 + (-t171 * t282 + t173 * t285) * t270;
t365 = t80 * t284;
t364 = Icges(3,4) * t283;
t363 = Icges(3,4) * t286;
t204 = -Icges(5,6) * t272 + (Icges(5,4) * t285 - Icges(5,2) * t282) * t270;
t353 = t282 * t204;
t289 = -pkin(8) - pkin(7);
t348 = t287 * t289;
t346 = t385 * t357;
t156 = t224 * rSges(6,1) + t223 * rSges(6,2) + rSges(6,3) * t357;
t297 = -t287 * t356 - t336;
t333 = pkin(3) * t354 + pkin(9) * t357;
t167 = t297 - t333;
t344 = -t156 - t167;
t166 = (t272 * t375 - t376) * t284 - t334;
t188 = (pkin(9) + t288) * t272 + t375 * t270;
t343 = t272 * t166 + t188 * t358;
t312 = -t222 * rSges(6,1) - t221 * rSges(6,2);
t154 = rSges(6,3) * t358 - t312;
t196 = -rSges(6,3) * t272 + (rSges(6,1) * t271 - rSges(6,2) * t269) * t270;
t112 = t272 * t154 + t196 * t358;
t341 = -t188 - t196;
t267 = pkin(2) * t286 + pkin(1);
t256 = t287 * t267;
t276 = t287 * pkin(7);
t340 = t284 * (t348 + t276 + (-pkin(1) + t267) * t284) + t287 * (-t287 * pkin(1) + t256 + (-pkin(7) - t289) * t284);
t298 = rSges(4,1) * t354 - rSges(4,2) * t357 + t284 * rSges(4,3);
t314 = rSges(4,1) * t272 - rSges(4,2) * t270;
t159 = t284 * (-t287 * rSges(4,3) + t284 * t314) + t287 * t298;
t206 = -t272 * rSges(5,3) + (rSges(5,1) * t285 - rSges(5,2) * t282) * t270;
t245 = t270 * pkin(3) - t272 * pkin(9);
t339 = -t206 - t245;
t338 = t278 * (t376 + t377) + t287 * t333;
t332 = t284 * rSges(3,3) + t287 * t373;
t330 = t278 + t279;
t329 = t393 * t357 + t358 * t394;
t328 = -t167 - t384;
t327 = -t188 - t383;
t326 = -t245 + t341;
t175 = t237 * rSges(5,1) + t236 * rSges(5,2) + rSges(5,3) * t357;
t325 = t358 / 0.2e1;
t324 = t357 / 0.2e1;
t244 = rSges(4,1) * t270 + rSges(4,2) * t272;
t323 = -t244 - t378;
t322 = -t245 - t378;
t74 = t168 * t357 + t170 * t236 + t172 * t237;
t75 = t169 * t357 + t171 * t236 + t173 * t237;
t38 = t284 * t75 - t287 * t74;
t321 = (t278 * t210 + t38 + (t304 * t287 + (-t209 + t303) * t284) * t287 + t391) * t284;
t320 = -t284 * t289 + t256;
t46 = t272 * t385 + t358 * t383;
t319 = t284 * t166 + t287 * t167 + t338;
t318 = -t245 + t327;
t313 = -t235 * rSges(5,1) - t234 * rSges(5,2);
t174 = rSges(5,3) * t358 - t313;
t85 = t284 * t174 + t287 * t175 + t338;
t317 = -t188 + t322;
t316 = -t206 + t322;
t315 = -t372 + t373;
t310 = Icges(3,1) * t286 - t364;
t308 = -Icges(3,2) * t283 + t363;
t306 = Icges(3,5) * t286 - Icges(3,6) * t283;
t241 = Icges(4,2) * t272 + t362;
t242 = Icges(4,1) * t270 + t361;
t300 = -t241 * t270 + t242 * t272;
t299 = -t196 + t317;
t44 = t284 * t154 + t287 * t156 + t319;
t296 = t317 - t383;
t295 = t272 * t374 + t329;
t294 = (t388 - t390) * t325 + (t387 - t389) * t324;
t30 = t284 * t385 + t287 * t384 + t319;
t293 = (t369 - t370 + t367 - t368) * t381 + t393 * t380 + t394 * t379 + t392 * t325 + t391 * t324;
t292 = t287 * t382 + t321;
t203 = -Icges(5,3) * t272 + (Icges(5,5) * t285 - Icges(5,6) * t282) * t270;
t205 = -Icges(5,5) * t272 + (Icges(5,1) * t285 - Icges(5,4) * t282) * t270;
t101 = t203 * t358 + t204 * t234 + t205 * t235;
t17 = -t101 * t272 + (t284 * t72 + t287 * t73) * t270;
t102 = t203 * t357 + t204 * t236 + t205 * t237;
t18 = -t102 * t272 + (t284 * t74 + t287 * t75) * t270;
t291 = t17 * t379 + t18 * t380 + t38 * t324 + t37 * t325 + t293 + (t365 - t366) * t381;
t240 = Icges(4,5) * t270 + Icges(4,6) * t272;
t290 = -t370 / 0.2e1 + t369 / 0.2e1 - t368 / 0.2e1 + t367 / 0.2e1 - t366 / 0.2e1 + t365 / 0.2e1 + (t212 * t272 + t214 * t270 + t284 * t240 + t287 * t300 + t102 + t387) * t380 + (t211 * t272 + t213 * t270 - t287 * t240 + t284 * t300 + t101 + t388) * t379;
t255 = rSges(2,1) * t287 - rSges(2,2) * t284;
t254 = -rSges(2,1) * t284 - rSges(2,2) * t287;
t253 = rSges(3,1) * t283 + rSges(3,2) * t286;
t227 = Icges(3,3) * t284 + t287 * t306;
t226 = -Icges(3,3) * t287 + t284 * t306;
t208 = t323 * t287;
t207 = t323 * t284;
t198 = t284 * pkin(7) + (pkin(1) - t372) * t287 + t332;
t197 = t371 + t276 + (-pkin(1) - t315) * t284;
t187 = t270 * t285 * t205;
t184 = t298 + t320;
t183 = (rSges(4,3) - t289) * t287 + (-t267 - t314) * t284;
t179 = t339 * t287;
t178 = t339 * t284;
t177 = t287 * (-t287 * t372 + t332) + (t284 * t315 - t371) * t284;
t163 = t316 * t287;
t162 = t316 * t284;
t140 = t166 * t357;
t131 = t154 * t357;
t125 = t320 + t175 + t333;
t124 = -t348 + (-t377 - t267 + (-rSges(5,3) - pkin(9)) * t270) * t284 + t313;
t120 = t326 * t287;
t119 = t326 * t284;
t117 = -t175 * t272 - t206 * t357;
t116 = t174 * t272 + t206 * t358;
t115 = t299 * t287;
t114 = t299 * t284;
t113 = -t272 * t156 - t196 * t357;
t111 = -t272 * t203 - t270 * t353 + t187;
t110 = t297 + t320 + t156;
t109 = -t348 + (-rSges(6,3) * t270 - t266 * t272 - t267) * t284 + t312 + t334;
t108 = t159 + t340;
t107 = (t174 * t287 - t175 * t284) * t270;
t104 = -t277 * t357 + t320 + t396;
t103 = (t247 - t289) * t287 + (-t246 * t272 - t267 + (-rSges(7,3) + t277) * t270) * t284 + t311;
t98 = -t156 * t358 + t131;
t97 = t318 * t287;
t96 = t318 * t284;
t91 = t296 * t287;
t90 = t296 * t284;
t70 = t272 * t344 + t341 * t357;
t69 = t112 + t343;
t60 = t85 + t340;
t47 = -t272 * t384 - t357 * t383;
t45 = t344 * t358 + t131 + t140;
t43 = -t358 * t384 + t346;
t42 = t44 + t340;
t41 = t272 * t328 + t327 * t357;
t40 = t46 + t343;
t31 = t328 * t358 + t140 + t346;
t19 = t30 + t340;
t1 = [t286 * (Icges(3,2) * t286 + t364) + t283 * (Icges(3,1) * t283 + t363) + Icges(2,3) + t187 + (-t203 + t241 - t406) * t272 + (-t353 + t242 - t397) * t270 + m(7) * (t103 ^ 2 + t104 ^ 2) + m(6) * (t109 ^ 2 + t110 ^ 2) + m(5) * (t124 ^ 2 + t125 ^ 2) + m(4) * (t183 ^ 2 + t184 ^ 2) + m(3) * (t197 ^ 2 + t198 ^ 2) + m(2) * (t254 ^ 2 + t255 ^ 2) + t395; t290 + (t286 * (-Icges(3,6) * t287 + t284 * t308) + t283 * (-Icges(3,5) * t287 + t284 * t310)) * t379 + (t286 * (Icges(3,6) * t284 + t287 * t308) + t283 * (Icges(3,5) * t284 + t287 * t310)) * t380 + m(7) * (t103 * t91 + t104 * t90) + m(6) * (t109 * t115 + t110 * t114) + m(5) * (t124 * t163 + t125 * t162) + m(4) * (t183 * t208 + t184 * t207) + (t279 / 0.2e1 + t278 / 0.2e1) * (Icges(3,5) * t283 + Icges(3,6) * t286) + m(3) * (-t197 * t287 - t198 * t284) * t253; m(7) * (t19 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(6) * (t114 ^ 2 + t115 ^ 2 + t42 ^ 2) + m(5) * (t162 ^ 2 + t163 ^ 2 + t60 ^ 2) + m(4) * (t108 ^ 2 + t207 ^ 2 + t208 ^ 2) + m(3) * (t253 ^ 2 * t330 + t177 ^ 2) + t284 * t278 * t227 + t321 + (-t279 * t226 + (-t284 * t226 + t287 * t227) * t284 + t382) * t287; t290 + m(7) * (t103 * t97 + t104 * t96) + m(6) * (t109 * t120 + t110 * t119) + m(5) * (t124 * t179 + t125 * t178) + m(4) * (-t183 * t287 - t184 * t284) * t244; m(7) * (t30 * t19 + t90 * t96 + t91 * t97) + m(6) * (t114 * t119 + t115 * t120 + t44 * t42) + m(5) * (t162 * t178 + t163 * t179 + t60 * t85) + m(4) * (t159 * t108 + (-t207 * t284 - t208 * t287) * t244) + t292; m(7) * (t30 ^ 2 + t96 ^ 2 + t97 ^ 2) + m(6) * (t119 ^ 2 + t120 ^ 2 + t44 ^ 2) + m(5) * (t178 ^ 2 + t179 ^ 2 + t85 ^ 2) + m(4) * (t244 ^ 2 * t330 + t159 ^ 2) + t292; (-t111 - t386) * t272 + m(7) * (t103 * t40 + t104 * t41) + m(6) * (t109 * t69 + t110 * t70) + m(5) * (t116 * t124 + t117 * t125) + ((t102 / 0.2e1 + t80 / 0.2e1) * t287 + (t101 / 0.2e1 + t79 / 0.2e1) * t284) * t270 + t294; m(7) * (t31 * t19 + t40 * t91 + t41 * t90) + m(6) * (t114 * t70 + t115 * t69 + t45 * t42) + m(5) * (t107 * t60 + t116 * t163 + t117 * t162) + t291; m(7) * (t31 * t30 + t40 * t97 + t41 * t96) + m(6) * (t119 * t70 + t120 * t69 + t45 * t44) + m(5) * (t107 * t85 + t116 * t179 + t117 * t178) + t291; (t111 * t272 + t374) * t272 + (t287 * t18 + t284 * t17 - t272 * (t284 * t79 + t287 * t80)) * t270 + m(7) * (t31 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t45 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t107 ^ 2 + t116 ^ 2 + t117 ^ 2) + t329; -t399 + m(7) * (t103 * t46 + t104 * t47) + m(6) * (t109 * t112 + t110 * t113) + t294; m(7) * (t43 * t19 + t46 * t91 + t47 * t90) + m(6) * (t112 * t115 + t113 * t114 + t42 * t98) + t293; m(7) * (t43 * t30 + t46 * t97 + t47 * t96) + m(6) * (t112 * t120 + t113 * t119 + t44 * t98) + t293; m(7) * (t31 * t43 + t40 * t46 + t41 * t47) + m(6) * (t112 * t69 + t113 * t70 + t45 * t98) + t295; m(7) * (t43 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(6) * (t112 ^ 2 + t113 ^ 2 + t98 ^ 2) + t295; m(7) * (t103 * t287 + t104 * t284) * t270; m(7) * (-t272 * t19 + (t284 * t90 + t287 * t91) * t270); m(7) * (-t272 * t30 + (t284 * t96 + t287 * t97) * t270); m(7) * (-t272 * t31 + (t284 * t41 + t287 * t40) * t270); m(7) * (-t272 * t43 + (t284 * t47 + t287 * t46) * t270); m(7) * (t270 ^ 2 * t330 + t272 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
