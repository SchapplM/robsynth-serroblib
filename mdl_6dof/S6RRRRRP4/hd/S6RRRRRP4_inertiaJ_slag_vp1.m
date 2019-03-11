% Calculate joint inertia matrix for
% S6RRRRRP4
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
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:12:07
% EndTime: 2019-03-10 01:12:23
% DurationCPUTime: 7.07s
% Computational Cost: add. (17499->505), mult. (17033->696), div. (0->0), fcn. (18641->10), ass. (0->262)
t410 = Icges(7,4) + Icges(6,5);
t276 = qJ(4) + qJ(5);
t269 = cos(t276);
t277 = qJ(2) + qJ(3);
t270 = cos(t277);
t283 = cos(qJ(1));
t267 = sin(t276);
t280 = sin(qJ(1));
t355 = t267 * t280;
t223 = t269 * t283 + t270 * t355;
t346 = t280 * t269;
t224 = -t267 * t283 + t270 * t346;
t268 = sin(t277);
t353 = t268 * t280;
t135 = Icges(7,5) * t224 + Icges(7,6) * t353 + Icges(7,3) * t223;
t141 = Icges(6,4) * t224 - Icges(6,2) * t223 + Icges(6,6) * t353;
t409 = t135 - t141;
t350 = t270 * t283;
t225 = t267 * t350 - t346;
t226 = t269 * t350 + t355;
t352 = t268 * t283;
t136 = Icges(7,5) * t226 + Icges(7,6) * t352 + Icges(7,3) * t225;
t142 = Icges(6,4) * t226 - Icges(6,2) * t225 + Icges(6,6) * t352;
t408 = t136 - t142;
t137 = Icges(6,5) * t224 - Icges(6,6) * t223 + Icges(6,3) * t353;
t139 = Icges(7,4) * t224 + Icges(7,2) * t353 + Icges(7,6) * t223;
t407 = t137 + t139;
t138 = Icges(6,5) * t226 - Icges(6,6) * t225 + Icges(6,3) * t352;
t140 = Icges(7,4) * t226 + Icges(7,2) * t352 + Icges(7,6) * t225;
t406 = t138 + t140;
t143 = Icges(7,1) * t224 + Icges(7,4) * t353 + Icges(7,5) * t223;
t145 = Icges(6,1) * t224 - Icges(6,4) * t223 + Icges(6,5) * t353;
t405 = t143 + t145;
t144 = Icges(7,1) * t226 + Icges(7,4) * t352 + Icges(7,5) * t225;
t146 = Icges(6,1) * t226 - Icges(6,4) * t225 + Icges(6,5) * t352;
t404 = t144 + t146;
t189 = -Icges(7,6) * t270 + (Icges(7,5) * t269 + Icges(7,3) * t267) * t268;
t192 = -Icges(6,6) * t270 + (Icges(6,4) * t269 - Icges(6,2) * t267) * t268;
t403 = t189 - t192;
t402 = (-Icges(6,3) - Icges(7,2)) * t270 + (t410 * t269 + (-Icges(6,6) + Icges(7,6)) * t267) * t268;
t401 = -t410 * t270 + ((Icges(6,1) + Icges(7,1)) * t269 + (-Icges(6,4) + Icges(7,5)) * t267) * t268;
t394 = rSges(7,3) + qJ(6);
t399 = rSges(7,1) + pkin(5);
t400 = -t394 * t223 - t399 * t224;
t398 = t223 * t409 + t405 * t224 + t407 * t353;
t397 = t223 * t408 + t224 * t404 + t353 * t406;
t396 = t225 * t409 + t405 * t226 + t407 * t352;
t395 = t225 * t408 + t226 * t404 + t352 * t406;
t384 = t223 * t403 + t224 * t401 + t353 * t402;
t383 = t225 * t403 + t226 * t401 + t352 * t402;
t356 = t267 * t268;
t391 = t268 * t269 * t401 + t189 * t356;
t382 = -t192 * t356 - t270 * t402 + t391;
t393 = t382 * t270;
t390 = -t384 * t270 + (t280 * t398 + t283 * t397) * t268;
t389 = -t383 * t270 + (t280 * t396 + t283 * t395) * t268;
t388 = t280 * t397 - t283 * t398;
t387 = t280 * t395 - t283 * t396;
t63 = -t139 * t270 + (t135 * t267 + t143 * t269) * t268;
t65 = -t137 * t270 + (-t141 * t267 + t145 * t269) * t268;
t386 = -t63 - t65;
t64 = -t140 * t270 + (t136 * t267 + t144 * t269) * t268;
t66 = -t138 * t270 + (-t142 * t267 + t146 * t269) * t268;
t385 = -t64 - t66;
t381 = rSges(7,2) * t353 - t400;
t380 = rSges(7,2) * t352 + t394 * t225 + t226 * t399;
t379 = -rSges(7,2) * t270 + (t394 * t267 + t269 * t399) * t268;
t302 = Icges(4,5) * t270 - Icges(4,6) * t268;
t209 = -Icges(4,3) * t283 + t280 * t302;
t210 = Icges(4,3) * t280 + t283 * t302;
t275 = t283 ^ 2;
t357 = Icges(4,4) * t270;
t304 = -Icges(4,2) * t268 + t357;
t212 = Icges(4,6) * t280 + t283 * t304;
t358 = Icges(4,4) * t268;
t306 = Icges(4,1) * t270 - t358;
t214 = Icges(4,5) * t280 + t283 * t306;
t300 = -t212 * t268 + t214 * t270;
t211 = -Icges(4,6) * t283 + t280 * t304;
t213 = -Icges(4,5) * t283 + t280 * t306;
t301 = t211 * t268 - t213 * t270;
t281 = cos(qJ(4));
t344 = t281 * t283;
t278 = sin(qJ(4));
t348 = t278 * t280;
t235 = -t270 * t348 - t344;
t345 = t280 * t281;
t347 = t278 * t283;
t236 = t270 * t345 - t347;
t165 = Icges(5,5) * t236 + Icges(5,6) * t235 + Icges(5,3) * t353;
t167 = Icges(5,4) * t236 + Icges(5,2) * t235 + Icges(5,6) * t353;
t169 = Icges(5,1) * t236 + Icges(5,4) * t235 + Icges(5,5) * t353;
t70 = t165 * t353 + t167 * t235 + t169 * t236;
t237 = -t270 * t347 + t345;
t238 = t270 * t344 + t348;
t166 = Icges(5,5) * t238 + Icges(5,6) * t237 + Icges(5,3) * t352;
t168 = Icges(5,4) * t238 + Icges(5,2) * t237 + Icges(5,6) * t352;
t170 = Icges(5,1) * t238 + Icges(5,4) * t237 + Icges(5,5) * t352;
t71 = t166 * t353 + t168 * t235 + t170 * t236;
t35 = t280 * t71 - t283 * t70;
t378 = -t35 - t275 * t209 - (t300 * t280 + (-t210 + t301) * t283) * t280 - t388;
t274 = t280 ^ 2;
t377 = -t270 / 0.2e1;
t376 = t280 / 0.2e1;
t375 = -t283 / 0.2e1;
t279 = sin(qJ(2));
t374 = pkin(2) * t279;
t373 = pkin(3) * t270;
t372 = pkin(9) * t268;
t265 = pkin(4) * t281 + pkin(3);
t371 = -pkin(3) + t265;
t370 = t393 + (t280 * t386 + t283 * t385) * t268;
t282 = cos(qJ(2));
t369 = rSges(3,1) * t282;
t368 = rSges(3,2) * t279;
t367 = t283 * rSges(3,3);
t366 = t63 * t283;
t365 = t64 * t280;
t364 = t65 * t283;
t363 = t66 * t280;
t79 = -t165 * t270 + (-t167 * t278 + t169 * t281) * t268;
t362 = t79 * t283;
t80 = -t166 * t270 + (-t168 * t278 + t170 * t281) * t268;
t361 = t80 * t280;
t360 = Icges(3,4) * t279;
t359 = Icges(3,4) * t282;
t284 = -pkin(10) - pkin(9);
t351 = t268 * t284;
t204 = -Icges(5,6) * t270 + (Icges(5,4) * t281 - Icges(5,2) * t278) * t268;
t349 = t278 * t204;
t285 = -pkin(8) - pkin(7);
t343 = t283 * t285;
t341 = t381 * t352;
t150 = t226 * rSges(6,1) - t225 * rSges(6,2) + rSges(6,3) * t352;
t294 = pkin(4) * t348 + t265 * t350 - t283 * t351;
t330 = pkin(3) * t350 + pkin(9) * t352;
t164 = t294 - t330;
t339 = -t150 - t164;
t331 = -pkin(4) * t347 - t280 * t351;
t163 = (t270 * t371 - t372) * t280 + t331;
t187 = (pkin(9) + t284) * t270 + t371 * t268;
t338 = t270 * t163 + t187 * t353;
t308 = -t224 * rSges(6,1) + t223 * rSges(6,2);
t148 = rSges(6,3) * t353 - t308;
t196 = -rSges(6,3) * t270 + (rSges(6,1) * t269 - rSges(6,2) * t267) * t268;
t112 = t270 * t148 + t196 * t353;
t336 = -t187 - t196;
t266 = pkin(2) * t282 + pkin(1);
t254 = t283 * t266;
t273 = t283 * pkin(7);
t334 = t280 * (t343 + t273 + (-pkin(1) + t266) * t280) + t283 * (-t283 * pkin(1) + t254 + (-pkin(7) - t285) * t280);
t295 = rSges(4,1) * t350 - rSges(4,2) * t352 + t280 * rSges(4,3);
t310 = rSges(4,1) * t270 - rSges(4,2) * t268;
t153 = t280 * (-t283 * rSges(4,3) + t280 * t310) + t283 * t295;
t206 = -rSges(5,3) * t270 + (rSges(5,1) * t281 - rSges(5,2) * t278) * t268;
t245 = pkin(3) * t268 - pkin(9) * t270;
t333 = -t206 - t245;
t332 = t274 * (t372 + t373) + t283 * t330;
t329 = t280 * rSges(3,3) + t283 * t369;
t328 = t274 + t275;
t327 = t389 * t352 + t353 * t390;
t326 = -t164 - t380;
t325 = -t187 - t379;
t324 = -t245 + t336;
t172 = t238 * rSges(5,1) + t237 * rSges(5,2) + rSges(5,3) * t352;
t323 = t353 / 0.2e1;
t322 = t352 / 0.2e1;
t244 = rSges(4,1) * t268 + rSges(4,2) * t270;
t321 = -t244 - t374;
t320 = -t245 - t374;
t72 = t165 * t352 + t167 * t237 + t169 * t238;
t73 = t166 * t352 + t168 * t237 + t170 * t238;
t36 = t280 * t73 - t283 * t72;
t319 = (t274 * t210 + t36 + (t301 * t283 + (-t209 + t300) * t280) * t283 + t387) * t280;
t318 = -t280 * t285 + t254;
t317 = -t265 * t270 - t266;
t316 = t280 * t163 + t283 * t164 + t332;
t315 = -t245 + t325;
t75 = t270 * t381 + t353 * t379;
t309 = -t236 * rSges(5,1) - t235 * rSges(5,2);
t171 = rSges(5,3) * t353 - t309;
t87 = t280 * t171 + t283 * t172 + t332;
t314 = -t187 + t320;
t313 = -t206 + t320;
t312 = -t331 - t343;
t311 = -t368 + t369;
t307 = Icges(3,1) * t282 - t360;
t305 = -Icges(3,2) * t279 + t359;
t303 = Icges(3,5) * t282 - Icges(3,6) * t279;
t241 = Icges(4,2) * t270 + t358;
t242 = Icges(4,1) * t268 + t357;
t297 = -t241 * t268 + t242 * t270;
t296 = -t196 + t314;
t43 = t280 * t148 + t283 * t150 + t316;
t293 = t314 - t379;
t292 = t270 * t370 + t327;
t291 = (t384 - t386) * t323 + (t383 - t385) * t322;
t37 = t280 * t381 + t283 * t380 + t316;
t290 = (t365 - t366 + t363 - t364) * t377 + t389 * t376 + t390 * t375 + t388 * t323 + t387 * t322;
t289 = t283 * t378 + t319;
t288 = t294 + t318;
t203 = -Icges(5,3) * t270 + (Icges(5,5) * t281 - Icges(5,6) * t278) * t268;
t205 = -Icges(5,5) * t270 + (Icges(5,1) * t281 - Icges(5,4) * t278) * t268;
t99 = t203 * t353 + t204 * t235 + t205 * t236;
t17 = -t270 * t99 + (t280 * t70 + t283 * t71) * t268;
t100 = t203 * t352 + t204 * t237 + t205 * t238;
t18 = -t100 * t270 + (t280 * t72 + t283 * t73) * t268;
t287 = t17 * t375 + t18 * t376 + t36 * t322 + t35 * t323 + t290 + (t361 - t362) * t377;
t240 = Icges(4,5) * t268 + Icges(4,6) * t270;
t286 = -t366 / 0.2e1 + t365 / 0.2e1 - t364 / 0.2e1 + t363 / 0.2e1 - t362 / 0.2e1 + t361 / 0.2e1 + (t212 * t270 + t214 * t268 + t240 * t280 + t283 * t297 + t100 + t383) * t376 + (t211 * t270 + t213 * t268 - t240 * t283 + t280 * t297 + t384 + t99) * t375;
t253 = rSges(2,1) * t283 - rSges(2,2) * t280;
t252 = -rSges(2,1) * t280 - rSges(2,2) * t283;
t251 = rSges(3,1) * t279 + rSges(3,2) * t282;
t228 = Icges(3,3) * t280 + t283 * t303;
t227 = -Icges(3,3) * t283 + t280 * t303;
t208 = t321 * t283;
t207 = t321 * t280;
t198 = t280 * pkin(7) + (pkin(1) - t368) * t283 + t329;
t197 = t367 + t273 + (-pkin(1) - t311) * t280;
t186 = t268 * t281 * t205;
t183 = t295 + t318;
t182 = (rSges(4,3) - t285) * t283 + (-t266 - t310) * t280;
t177 = t333 * t283;
t176 = t333 * t280;
t175 = t283 * (-t283 * t368 + t329) + (t280 * t311 - t367) * t280;
t157 = t313 * t283;
t156 = t313 * t280;
t134 = t163 * t352;
t125 = t148 * t352;
t121 = t318 + t172 + t330;
t120 = -t343 + (-t373 - t266 + (-rSges(5,3) - pkin(9)) * t268) * t280 + t309;
t119 = t324 * t283;
t118 = t324 * t280;
t117 = -t172 * t270 - t206 * t352;
t116 = t171 * t270 + t206 * t353;
t115 = t296 * t283;
t114 = t296 * t280;
t113 = -t150 * t270 - t196 * t352;
t111 = -t270 * t203 - t268 * t349 + t186;
t110 = t288 + t150;
t109 = (-rSges(6,3) * t268 + t317) * t280 + t308 + t312;
t108 = t153 + t334;
t107 = t315 * t283;
t106 = t315 * t280;
t105 = (t171 * t283 - t172 * t280) * t268;
t102 = t293 * t283;
t101 = t293 * t280;
t96 = -t150 * t353 + t125;
t86 = t288 + t380;
t85 = (-rSges(7,2) * t268 + t317) * t280 + t312 + t400;
t76 = -t270 * t380 - t352 * t379;
t68 = t270 * t339 + t336 * t352;
t67 = t112 + t338;
t58 = t87 + t334;
t45 = -t353 * t380 + t341;
t44 = t339 * t353 + t125 + t134;
t42 = t270 * t326 + t325 * t352;
t41 = t75 + t338;
t40 = t43 + t334;
t38 = t326 * t353 + t134 + t341;
t29 = t37 + t334;
t1 = [t282 * (Icges(3,2) * t282 + t360) + t279 * (Icges(3,1) * t279 + t359) + Icges(2,3) + t186 + (-t267 * t192 + t242 - t349) * t268 + (-t203 + t241 - t402) * t270 + m(7) * (t85 ^ 2 + t86 ^ 2) + m(6) * (t109 ^ 2 + t110 ^ 2) + m(5) * (t120 ^ 2 + t121 ^ 2) + m(4) * (t182 ^ 2 + t183 ^ 2) + m(3) * (t197 ^ 2 + t198 ^ 2) + m(2) * (t252 ^ 2 + t253 ^ 2) + t391; (t274 / 0.2e1 + t275 / 0.2e1) * (Icges(3,5) * t279 + Icges(3,6) * t282) + m(7) * (t101 * t86 + t102 * t85) + m(6) * (t109 * t115 + t110 * t114) + m(5) * (t120 * t157 + t121 * t156) + m(4) * (t182 * t208 + t183 * t207) + t286 + ((-Icges(3,6) * t283 + t280 * t305) * t282 + (-Icges(3,5) * t283 + t280 * t307) * t279) * t375 + ((Icges(3,6) * t280 + t283 * t305) * t282 + (Icges(3,5) * t280 + t283 * t307) * t279) * t376 + m(3) * (-t197 * t283 - t198 * t280) * t251; m(7) * (t101 ^ 2 + t102 ^ 2 + t29 ^ 2) + m(6) * (t114 ^ 2 + t115 ^ 2 + t40 ^ 2) + m(5) * (t156 ^ 2 + t157 ^ 2 + t58 ^ 2) + m(4) * (t108 ^ 2 + t207 ^ 2 + t208 ^ 2) + t280 * t274 * t228 + m(3) * (t251 ^ 2 * t328 + t175 ^ 2) + t319 + (-t275 * t227 + (-t280 * t227 + t283 * t228) * t280 + t378) * t283; m(7) * (t106 * t86 + t107 * t85) + m(6) * (t109 * t119 + t110 * t118) + m(5) * (t120 * t177 + t121 * t176) + t286 + m(4) * (-t182 * t283 - t183 * t280) * t244; m(7) * (t101 * t106 + t102 * t107 + t29 * t37) + m(6) * (t114 * t118 + t115 * t119 + t40 * t43) + m(5) * (t156 * t176 + t157 * t177 + t58 * t87) + m(4) * (t108 * t153 + (-t207 * t280 - t208 * t283) * t244) + t289; m(7) * (t106 ^ 2 + t107 ^ 2 + t37 ^ 2) + m(6) * (t118 ^ 2 + t119 ^ 2 + t43 ^ 2) + m(5) * (t176 ^ 2 + t177 ^ 2 + t87 ^ 2) + m(4) * (t244 ^ 2 * t328 + t153 ^ 2) + t289; (-t111 - t382) * t270 + m(7) * (t41 * t85 + t42 * t86) + m(6) * (t109 * t67 + t110 * t68) + m(5) * (t116 * t120 + t117 * t121) + ((t80 / 0.2e1 + t100 / 0.2e1) * t283 + (t79 / 0.2e1 + t99 / 0.2e1) * t280) * t268 + t291; m(7) * (t101 * t42 + t102 * t41 + t29 * t38) + m(6) * (t114 * t68 + t115 * t67 + t40 * t44) + m(5) * (t105 * t58 + t116 * t157 + t117 * t156) + t287; m(7) * (t106 * t42 + t107 * t41 + t37 * t38) + m(6) * (t118 * t68 + t119 * t67 + t43 * t44) + m(5) * (t105 * t87 + t116 * t177 + t117 * t176) + t287; (t111 * t270 + t370) * t270 + (t283 * t18 + t280 * t17 - t270 * (t79 * t280 + t80 * t283)) * t268 + m(7) * (t38 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(6) * (t44 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(5) * (t105 ^ 2 + t116 ^ 2 + t117 ^ 2) + t327; -t393 + m(7) * (t75 * t85 + t76 * t86) + m(6) * (t109 * t112 + t110 * t113) + t291; m(7) * (t101 * t76 + t102 * t75 + t29 * t45) + m(6) * (t112 * t115 + t113 * t114 + t40 * t96) + t290; m(7) * (t106 * t76 + t107 * t75 + t37 * t45) + m(6) * (t112 * t119 + t113 * t118 + t43 * t96) + t290; m(7) * (t38 * t45 + t41 * t75 + t42 * t76) + m(6) * (t112 * t67 + t113 * t68 + t44 * t96) + t292; m(7) * (t45 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(6) * (t112 ^ 2 + t113 ^ 2 + t96 ^ 2) + t292; m(7) * (t223 * t86 + t225 * t85); m(7) * (t101 * t223 + t102 * t225 + t29 * t356); m(7) * (t106 * t223 + t107 * t225 + t356 * t37); m(7) * (t223 * t42 + t225 * t41 + t356 * t38); m(7) * (t223 * t76 + t225 * t75 + t356 * t45); m(7) * (t267 ^ 2 * t268 ^ 2 + t223 ^ 2 + t225 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
