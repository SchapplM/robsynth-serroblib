% Calculate time derivative of joint inertia matrix for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:34
% EndTime: 2019-12-05 17:09:51
% DurationCPUTime: 7.54s
% Computational Cost: add. (29688->530), mult. (30076->811), div. (0->0), fcn. (28582->10), ass. (0->304)
t252 = qJD(2) + qJD(3);
t255 = sin(pkin(9));
t249 = t255 ^ 2;
t256 = cos(pkin(9));
t250 = t256 ^ 2;
t417 = t249 + t250;
t411 = t252 * t417;
t253 = qJ(4) + qJ(5);
t245 = sin(t253);
t247 = cos(t253);
t254 = qJ(2) + qJ(3);
t246 = sin(t254);
t251 = qJD(4) + qJD(5);
t365 = t252 * t256;
t280 = t246 * t365 - t251 * t255;
t248 = cos(t254);
t367 = t251 * t256;
t339 = t248 * t367;
t177 = t280 * t245 - t247 * t339;
t178 = -t245 * t339 - t280 * t247;
t337 = t248 * t365;
t102 = rSges(6,1) * t178 + rSges(6,2) * t177 + rSges(6,3) * t337;
t259 = cos(qJ(4));
t394 = t259 * pkin(4);
t188 = -pkin(8) * t248 + t394 * t246;
t257 = sin(qJ(4));
t348 = qJD(4) * t257;
t344 = pkin(4) * t348;
t263 = -t188 * t252 - t248 * t344;
t347 = qJD(4) * t259;
t343 = pkin(4) * t347;
t418 = t255 * t343 + t256 * t263 + t102;
t416 = -Icges(4,5) * t246 - Icges(4,6) * t248;
t364 = t255 * t245;
t336 = t248 * t364;
t219 = -t247 * t256 - t336;
t363 = t255 * t247;
t335 = t248 * t363;
t220 = -t245 * t256 + t335;
t373 = t246 * t255;
t152 = rSges(6,1) * t220 + rSges(6,2) * t219 + rSges(6,3) * t373;
t266 = pkin(8) * t246 + t394 * t248;
t360 = t256 * t257;
t163 = -pkin(4) * t360 + t255 * t266;
t415 = t152 + t163;
t368 = t248 * t256;
t221 = -t245 * t368 + t363;
t222 = t247 * t368 + t364;
t372 = t246 * t256;
t153 = t222 * rSges(6,1) + t221 * rSges(6,2) + rSges(6,3) * t372;
t362 = t255 * t257;
t164 = pkin(4) * t362 + t256 * t266;
t414 = t153 + t164;
t145 = Icges(6,4) * t220 + Icges(6,2) * t219 + Icges(6,6) * t373;
t147 = Icges(6,1) * t220 + Icges(6,4) * t219 + Icges(6,5) * t373;
t366 = t252 * t255;
t279 = t246 * t366 + t367;
t175 = t279 * t245 - t251 * t335;
t176 = -t279 * t247 - t251 * t336;
t143 = Icges(6,5) * t220 + Icges(6,6) * t219 + Icges(6,3) * t373;
t369 = t248 * t252;
t338 = t248 * t366;
t95 = Icges(6,5) * t176 + Icges(6,6) * t175 + Icges(6,3) * t338;
t285 = t143 * t369 + t246 * t95;
t97 = Icges(6,4) * t176 + Icges(6,2) * t175 + Icges(6,6) * t338;
t99 = Icges(6,1) * t176 + Icges(6,4) * t175 + Icges(6,5) * t338;
t27 = t145 * t175 + t147 * t176 + t219 * t97 + t220 * t99 + t285 * t255;
t100 = Icges(6,1) * t178 + Icges(6,4) * t177 + Icges(6,5) * t337;
t146 = Icges(6,4) * t222 + Icges(6,2) * t221 + Icges(6,6) * t372;
t148 = Icges(6,1) * t222 + Icges(6,4) * t221 + Icges(6,5) * t372;
t144 = Icges(6,5) * t222 + Icges(6,6) * t221 + Icges(6,3) * t372;
t96 = Icges(6,5) * t178 + Icges(6,6) * t177 + Icges(6,3) * t337;
t284 = t144 * t369 + t246 * t96;
t98 = Icges(6,4) * t178 + Icges(6,2) * t177 + Icges(6,6) * t337;
t28 = t100 * t220 + t146 * t175 + t148 * t176 + t219 * t98 + t284 * t255;
t16 = t255 * t28 - t256 * t27;
t272 = t252 * t416;
t201 = t255 * t272;
t202 = t256 * t272;
t361 = t255 * t259;
t236 = t248 * t361 - t360;
t374 = t246 * t252;
t341 = t257 * t374;
t182 = -t236 * qJD(4) + t255 * t341;
t359 = t256 * t259;
t235 = -t248 * t362 - t359;
t340 = t259 * t374;
t183 = t235 * qJD(4) - t255 * t340;
t116 = Icges(5,4) * t183 + Icges(5,2) * t182 + Icges(5,6) * t338;
t118 = Icges(5,1) * t183 + Icges(5,4) * t182 + Icges(5,5) * t338;
t167 = Icges(5,4) * t236 + Icges(5,2) * t235 + Icges(5,6) * t373;
t169 = Icges(5,1) * t236 + Icges(5,4) * t235 + Icges(5,5) * t373;
t114 = Icges(5,5) * t183 + Icges(5,6) * t182 + Icges(5,3) * t338;
t165 = Icges(5,5) * t236 + Icges(5,6) * t235 + Icges(5,3) * t373;
t282 = t114 * t246 + t165 * t369;
t37 = t116 * t235 + t118 * t236 + t167 * t182 + t169 * t183 + t282 * t255;
t238 = t248 * t359 + t362;
t184 = -t238 * qJD(4) + t256 * t341;
t237 = -t248 * t360 + t361;
t185 = t237 * qJD(4) - t256 * t340;
t117 = Icges(5,4) * t185 + Icges(5,2) * t184 + Icges(5,6) * t337;
t119 = Icges(5,1) * t185 + Icges(5,4) * t184 + Icges(5,5) * t337;
t168 = Icges(5,4) * t238 + Icges(5,2) * t237 + Icges(5,6) * t372;
t170 = Icges(5,1) * t238 + Icges(5,4) * t237 + Icges(5,5) * t372;
t115 = Icges(5,5) * t185 + Icges(5,6) * t184 + Icges(5,3) * t337;
t166 = Icges(5,5) * t238 + Icges(5,6) * t237 + Icges(5,3) * t372;
t281 = t115 * t246 + t166 * t369;
t38 = t117 * t235 + t119 * t236 + t168 * t182 + t170 * t183 + t281 * t255;
t22 = t255 * t38 - t256 * t37;
t403 = t416 * t411;
t413 = -t16 - t22 - t201 * t250 - (-t202 * t256 + t403) * t255;
t120 = rSges(5,1) * t183 + rSges(5,2) * t182 + rSges(5,3) * t338;
t121 = rSges(5,1) * t185 + rSges(5,2) * t184 + rSges(5,3) * t337;
t408 = t255 * t120 + t256 * t121;
t258 = sin(qJ(2));
t260 = cos(qJ(2));
t407 = qJD(2) * (t258 * rSges(3,1) + t260 * rSges(3,2));
t101 = rSges(6,1) * t176 + rSges(6,2) * t175 + rSges(6,3) * t338;
t128 = t255 * t263 - t256 * t343;
t405 = t418 * t256 + (t101 + t128) * t255;
t402 = 2 * m(4);
t401 = 2 * m(5);
t400 = 2 * m(6);
t399 = -t248 / 0.2e1;
t398 = t255 / 0.2e1;
t397 = -t256 / 0.2e1;
t396 = pkin(2) * t258;
t391 = pkin(2) * qJD(2);
t65 = t143 * t372 + t145 * t221 + t147 * t222;
t390 = t255 * t65;
t74 = t165 * t372 + t167 * t237 + t169 * t238;
t389 = t255 * t74;
t64 = t144 * t373 + t146 * t219 + t148 * t220;
t388 = t256 * t64;
t73 = t166 * t373 + t168 * t235 + t170 * t236;
t387 = t256 * t73;
t386 = t101 * t372 + t152 * t337;
t381 = Icges(5,4) * t257;
t380 = Icges(5,4) * t259;
t379 = Icges(6,4) * t245;
t378 = Icges(6,4) * t247;
t299 = Icges(6,5) * t247 - Icges(6,6) * t245;
t377 = (t299 * t369 + (Icges(6,3) * t252 + (-Icges(6,5) * t245 - Icges(6,6) * t247) * t251) * t246) * t248;
t193 = -Icges(6,3) * t248 + t299 * t246;
t376 = t193 * t248;
t300 = Icges(5,5) * t259 - Icges(5,6) * t257;
t371 = t248 * (t300 * t369 + (Icges(5,3) * t252 + (-Icges(5,5) * t257 - Icges(5,6) * t259) * qJD(4)) * t246);
t207 = -Icges(5,3) * t248 + t300 * t246;
t370 = t248 * t207;
t313 = rSges(6,1) * t247 - rSges(6,2) * t245;
t136 = t313 * t369 + (rSges(6,3) * t252 + (-rSges(6,1) * t245 - rSges(6,2) * t247) * t251) * t246;
t156 = -t246 * t344 + t266 * t252;
t357 = -t136 - t156;
t314 = rSges(5,1) * t259 - rSges(5,2) * t257;
t154 = t314 * t369 + (rSges(5,3) * t252 + (-rSges(5,1) * t257 - rSges(5,2) * t259) * qJD(4)) * t246;
t317 = pkin(3) * t248 + pkin(7) * t246;
t228 = t317 * t252;
t355 = -t154 - t228;
t198 = -rSges(6,3) * t248 + t313 * t246;
t103 = t248 * t152 + t198 * t373;
t354 = -t188 - t198;
t239 = rSges(4,1) * t246 + rSges(4,2) * t248;
t140 = t239 * t411;
t240 = t246 * pkin(3) - t248 * pkin(7);
t353 = t240 * t411;
t352 = t417 * t260 * pkin(2);
t315 = rSges(4,1) * t248 - rSges(4,2) * t246;
t155 = t417 * t315;
t210 = -t248 * rSges(5,3) + t314 * t246;
t351 = -t210 - t240;
t350 = t417 * t317;
t29 = t177 * t145 + t178 * t147 + t221 * t97 + t222 * t99 + t285 * t256;
t30 = t222 * t100 + t177 * t146 + t178 * t148 + t221 * t98 + t284 * t256;
t17 = t255 * t30 - t256 * t29;
t39 = t237 * t116 + t238 * t118 + t184 * t167 + t185 * t169 + t282 * t256;
t40 = t237 * t117 + t238 * t119 + t184 * t168 + t185 * t170 + t281 * t256;
t23 = t255 * t40 - t256 * t39;
t346 = (t249 * t202 + t17 + t23 + (-t255 * t201 + t403) * t256) * t255;
t345 = t260 * t391;
t342 = t248 * t101 + t136 * t373 + t198 * t338;
t334 = -t228 + t357;
t333 = -t240 + t354;
t332 = t374 / 0.2e1;
t331 = t373 / 0.2e1;
t330 = t372 / 0.2e1;
t329 = t369 / 0.2e1;
t328 = -t239 - t396;
t327 = -t240 - t396;
t326 = t256 * t354;
t325 = t414 * t248;
t171 = rSges(5,1) * t236 + rSges(5,2) * t235 + rSges(5,3) * t373;
t172 = t238 * rSges(5,1) + t237 * rSges(5,2) + rSges(5,3) * t372;
t82 = t255 * t171 + t256 * t172 + t350;
t322 = t255 * t329;
t321 = t256 * t329;
t320 = -t210 + t327;
t227 = t315 * t252;
t319 = -t227 - t345;
t318 = -t228 - t345;
t298 = -t145 * t245 + t147 * t247;
t68 = -t143 * t248 + t298 * t246;
t297 = -t146 * t245 + t148 * t247;
t69 = -t144 * t248 + t297 * t246;
t312 = t68 * t255 + t69 * t256;
t296 = -t167 * t257 + t169 * t259;
t76 = -t248 * t165 + t296 * t246;
t295 = -t168 * t257 + t170 * t259;
t77 = -t248 * t166 + t295 * t246;
t311 = t76 * t255 + t77 * t256;
t306 = Icges(5,1) * t259 - t381;
t305 = Icges(6,1) * t247 - t379;
t302 = -Icges(5,2) * t257 + t380;
t301 = -Icges(6,2) * t245 + t378;
t294 = t171 * t256 - t172 * t255;
t194 = -Icges(6,6) * t248 + t301 * t246;
t195 = -Icges(6,5) * t248 + t305 * t246;
t293 = t194 * t245 - t195 * t247;
t208 = -Icges(5,6) * t248 + t302 * t246;
t209 = -Icges(5,5) * t248 + t306 * t246;
t292 = t208 * t257 - t209 * t259;
t288 = t327 + t354;
t287 = t417 * t258 * t391;
t286 = -t154 + t318;
t57 = t415 * t255 + t414 * t256 + t350;
t133 = t301 * t369 + (Icges(6,6) * t252 + (-Icges(6,2) * t247 - t379) * t251) * t246;
t134 = t305 * t369 + (Icges(6,5) * t252 + (-Icges(6,1) * t245 - t378) * t251) * t246;
t24 = (t252 * t298 - t95) * t248 + (t143 * t252 + (-t145 * t251 + t99) * t247 + (-t147 * t251 - t97) * t245) * t246;
t25 = (t252 * t297 - t96) * t248 + (t144 * t252 + (-t146 * t251 + t100) * t247 + (-t148 * t251 - t98) * t245) * t246;
t63 = t143 * t373 + t145 * t219 + t147 * t220;
t83 = t193 * t373 + t194 * t219 + t195 * t220;
t5 = (-t133 * t219 - t134 * t220 - t175 * t194 - t176 * t195 + (t388 + (t63 - t376) * t255) * t252) * t248 + (t83 * t252 + t28 * t256 + (t27 - t377) * t255) * t246;
t66 = t144 * t372 + t146 * t221 + t148 * t222;
t84 = t193 * t372 + t194 * t221 + t195 * t222;
t6 = (-t221 * t133 - t222 * t134 - t177 * t194 - t178 * t195 + (t390 + (t66 - t376) * t256) * t252) * t248 + (t84 * t252 + t29 * t255 + (t30 - t377) * t256) * t246;
t91 = -t293 * t246 - t376;
t278 = -t248 * ((t377 + (t293 * t248 + t312) * t252) * t248 + (t25 * t256 + t24 * t255 - (t193 * t252 + (-t194 * t251 + t134) * t247 + (-t195 * t251 - t133) * t245) * t248 + t91 * t252) * t246) + t6 * t372 + t5 * t373 + (-t83 * t248 + (t255 * t63 + t388) * t246) * t338 + (-t84 * t248 + (t256 * t66 + t390) * t246) * t337 + (t312 * t246 - t91 * t248) * t374;
t276 = t318 + t357;
t275 = t5 * t397 + t6 * t398 + (-t24 * t256 + t25 * t255) * t399 + t16 * t331 + t17 * t330 + (t255 * t64 - t256 * t63) * t322 + (t255 * t66 - t256 * t65) * t321 + (t255 * t69 - t256 * t68) * t332;
t268 = qJD(2) * (-Icges(3,5) * t258 - Icges(3,6) * t260);
t267 = t413 * t256 + t346;
t265 = -t287 - t353;
t150 = t302 * t369 + (Icges(5,6) * t252 + (-Icges(5,2) * t259 - t381) * qJD(4)) * t246;
t151 = t306 * t369 + (Icges(5,5) * t252 + (-Icges(5,1) * t257 - t380) * qJD(4)) * t246;
t72 = t165 * t373 + t167 * t235 + t169 * t236;
t86 = t207 * t373 + t208 * t235 + t209 * t236;
t10 = (-t150 * t235 - t151 * t236 - t182 * t208 - t183 * t209 + (t387 + (t72 - t370) * t255) * t252) * t248 + (t86 * t252 + t38 * t256 + (t37 - t371) * t255) * t246;
t75 = t166 * t372 + t168 * t237 + t170 * t238;
t87 = t207 * t372 + t208 * t237 + t209 * t238;
t11 = (-t237 * t150 - t238 * t151 - t184 * t208 - t185 * t209 + (t389 + (t75 - t370) * t256) * t252) * t248 + (t87 * t252 + t39 * t255 + (t40 - t371) * t256) * t246;
t33 = (t296 * t252 - t114) * t248 + (-t116 * t257 + t118 * t259 + t165 * t252 + (-t167 * t259 - t169 * t257) * qJD(4)) * t246;
t34 = (t295 * t252 - t115) * t248 + (-t117 * t257 + t119 * t259 + t166 * t252 + (-t168 * t259 - t170 * t257) * qJD(4)) * t246;
t264 = t10 * t397 + t11 * t398 + (t255 * t34 - t256 * t33) * t399 + t22 * t331 + t23 * t330 + t275 + (t255 * t73 - t256 * t72) * t322 + (t255 * t75 - t256 * t74) * t321 + (t255 * t77 - t256 * t76) * t332;
t230 = t256 * t268;
t229 = t255 * t268;
t212 = t328 * t256;
t211 = t328 * t255;
t187 = t319 * t256;
t186 = t319 * t255;
t179 = t417 * t407;
t174 = t351 * t256;
t173 = t351 * t255;
t160 = t320 * t256;
t159 = t320 * t255;
t138 = t152 * t372;
t137 = t153 * t374;
t130 = -t287 - t140;
t125 = t355 * t256;
t124 = t355 * t255;
t123 = t333 * t256;
t122 = t333 * t255;
t113 = -t172 * t248 - t210 * t372;
t112 = t171 * t248 + t210 * t373;
t109 = t288 * t256;
t108 = t288 * t255;
t107 = t286 * t256;
t106 = t286 * t255;
t105 = -t292 * t246 - t370;
t104 = -t248 * t153 - t198 * t372;
t90 = t155 + t352;
t88 = t294 * t246;
t85 = -t153 * t373 + t138;
t81 = t334 * t256;
t80 = t334 * t255;
t79 = t276 * t256;
t78 = t276 * t255;
t71 = t246 * t326 - t325;
t70 = t163 * t248 + t188 * t373 + t103;
t67 = t82 + t352;
t62 = -t353 + t408;
t61 = t138 + (t163 * t256 - t255 * t414) * t246;
t60 = (-t210 * t365 - t121) * t248 + (-t154 * t256 + t172 * t252) * t246;
t59 = (t210 * t366 + t120) * t248 + (t154 * t255 - t171 * t252) * t246;
t58 = t265 + t408;
t56 = -t136 * t372 + t137 + (-t198 * t365 - t102) * t248;
t55 = -t152 * t374 + t342;
t54 = t294 * t369 + (t120 * t256 - t121 * t255) * t246;
t53 = t57 + t352;
t49 = (-t102 * t246 - t153 * t369) * t255 + t386;
t45 = -t353 + t405;
t44 = t265 + t405;
t36 = t137 + (t164 * t252 + t357 * t256) * t246 + (t252 * t326 - t418) * t248;
t35 = (t188 * t366 + t128) * t248 + (t255 * t156 - t415 * t252) * t246 + t342;
t26 = (t128 * t246 + t163 * t369) * t256 + (-t246 * t418 - t252 * t325) * t255 + t386;
t1 = [0; -m(3) * t179 + m(4) * t130 + m(5) * t58 + m(6) * t44; (t108 * t78 + t109 * t79 + t44 * t53) * t400 + (t106 * t159 + t107 * t160 + t58 * t67) * t401 + (t130 * t90 + t186 * t211 + t187 * t212) * t402 + t255 * t249 * t230 + t346 + 0.2e1 * m(3) * (-t179 + t407) * t417 * (rSges(3,1) * t260 - rSges(3,2) * t258) + (-t250 * t229 + (-t255 * t229 + t256 * t230) * t255 + t413) * t256; -m(4) * t140 + m(5) * t62 + m(6) * t45; m(6) * (t108 * t80 + t109 * t81 + t122 * t78 + t123 * t79 + t44 * t57 + t45 * t53) + m(5) * (t106 * t173 + t107 * t174 + t124 * t159 + t125 * t160 + t58 * t82 + t62 * t67) + m(4) * (t155 * t130 - t140 * t90 + (-t186 * t255 - t187 * t256) * t239 + (-t211 * t255 - t212 * t256) * t227) + t267; (t122 * t80 + t123 * t81 + t45 * t57) * t400 + (t124 * t173 + t125 * t174 + t62 * t82) * t401 + (t227 * t239 * t417 - t140 * t155) * t402 + t267; m(5) * t54 + m(6) * t26; m(6) * (t108 * t36 + t109 * t35 + t26 * t53 + t44 * t61 + t70 * t79 + t71 * t78) + m(5) * (t106 * t113 + t107 * t112 + t159 * t60 + t160 * t59 + t54 * t67 + t58 * t88) + t264; m(6) * (t122 * t36 + t123 * t35 + t26 * t57 + t45 * t61 + t70 * t81 + t71 * t80) + m(5) * (t112 * t125 + t113 * t124 + t173 * t60 + t174 * t59 + t54 * t82 + t62 * t88) + t264; (-t86 * t248 + (t255 * t72 + t387) * t246) * t338 + t10 * t373 + (-t87 * t248 + (t256 * t75 + t389) * t246) * t337 + t11 * t372 + (t26 * t61 + t35 * t70 + t36 * t71) * t400 + (-t105 * t248 + t246 * t311) * t374 - t248 * ((t371 + (t248 * t292 + t311) * t252) * t248 + (t34 * t256 + t33 * t255 - (-t150 * t257 + t151 * t259 + t207 * t252 - t208 * t347 - t209 * t348) * t248 + t105 * t252) * t246) + (t112 * t59 + t113 * t60 + t54 * t88) * t401 + t278; m(6) * t49; m(6) * (t103 * t79 + t104 * t78 + t108 * t56 + t109 * t55 + t44 * t85 + t49 * t53) + t275; m(6) * (t103 * t81 + t104 * t80 + t122 * t56 + t123 * t55 + t45 * t85 + t49 * t57) + t275; m(6) * (t103 * t35 + t104 * t36 + t26 * t85 + t49 * t61 + t55 * t70 + t56 * t71) + t278; (t103 * t55 + t104 * t56 + t49 * t85) * t400 + t278;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
