% Calculate joint inertia matrix for
% S6RRRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:03:58
% EndTime: 2019-03-09 17:04:11
% DurationCPUTime: 6.65s
% Computational Cost: add. (26508->652), mult. (45618->897), div. (0->0), fcn. (57636->12), ass. (0->300)
t281 = cos(pkin(6));
t286 = sin(qJ(1));
t288 = cos(qJ(2));
t346 = t286 * t288;
t285 = sin(qJ(2));
t289 = cos(qJ(1));
t347 = t285 * t289;
t263 = t281 * t347 + t346;
t328 = qJ(3) + pkin(11);
t278 = sin(t328);
t309 = cos(t328);
t280 = sin(pkin(6));
t350 = t280 * t289;
t232 = t263 * t309 - t278 * t350;
t345 = t288 * t289;
t348 = t285 * t286;
t262 = -t281 * t345 + t348;
t283 = sin(qJ(5));
t359 = cos(qJ(5));
t193 = t232 * t283 - t262 * t359;
t194 = t232 * t359 + t262 * t283;
t301 = t280 * t309;
t231 = t263 * t278 + t289 * t301;
t362 = rSges(7,3) + qJ(6);
t363 = rSges(7,1) + pkin(5);
t344 = rSges(7,2) * t231 + t362 * t193 + t194 * t363;
t249 = t281 * t278 + t285 * t301;
t351 = t280 * t288;
t229 = t249 * t283 + t351 * t359;
t230 = t249 * t359 - t283 * t351;
t354 = t280 * t285;
t248 = t278 * t354 - t281 * t309;
t141 = Icges(7,5) * t230 + Icges(7,6) * t248 + Icges(7,3) * t229;
t143 = Icges(7,4) * t230 + Icges(7,2) * t248 + Icges(7,6) * t229;
t145 = Icges(7,1) * t230 + Icges(7,4) * t248 + Icges(7,5) * t229;
t70 = t229 * t141 + t248 * t143 + t230 * t145;
t142 = Icges(6,5) * t230 - Icges(6,6) * t229 + Icges(6,3) * t248;
t144 = Icges(6,4) * t230 - Icges(6,2) * t229 + Icges(6,6) * t248;
t146 = Icges(6,1) * t230 - Icges(6,4) * t229 + Icges(6,5) * t248;
t71 = t248 * t142 - t229 * t144 + t230 * t146;
t361 = -t70 - t71;
t360 = t280 ^ 2;
t287 = cos(qJ(3));
t277 = pkin(3) * t287 + pkin(2);
t358 = -pkin(2) + t277;
t209 = Icges(3,5) * t263 - Icges(3,6) * t262 - Icges(3,3) * t350;
t356 = t209 * t289;
t282 = -qJ(4) - pkin(9);
t355 = t262 * t282;
t353 = t280 * t286;
t352 = t280 * t287;
t284 = sin(qJ(3));
t349 = t281 * t284;
t265 = -t281 * t348 + t345;
t234 = t265 * t309 + t278 * t353;
t264 = t281 * t346 + t347;
t195 = t234 * t283 - t264 * t359;
t196 = t234 * t359 + t264 * t283;
t233 = t265 * t278 - t286 * t301;
t343 = t233 * rSges(7,2) + t362 * t195 + t196 * t363;
t258 = t262 * pkin(9);
t317 = t284 * t350;
t269 = pkin(3) * t317;
t162 = t263 * t358 - t258 - t269 - t355;
t137 = t264 * t162;
t182 = pkin(4) * t232 + t231 * pkin(10);
t342 = t264 * t182 + t137;
t341 = rSges(7,2) * t248 + t362 * t229 + t230 * t363;
t207 = pkin(3) * t349 + ((pkin(9) + t282) * t288 + t358 * t285) * t280;
t340 = t162 * t351 + t262 * t207;
t222 = t265 * pkin(2) + pkin(9) * t264;
t318 = t284 * t353;
t311 = pkin(3) * t318 - t264 * t282 + t265 * t277;
t163 = -t222 + t311;
t220 = t281 * t222;
t339 = t281 * t163 + t220;
t160 = t234 * rSges(5,1) - t233 * rSges(5,2) + t264 * rSges(5,3);
t338 = -t160 - t163;
t221 = pkin(2) * t263 + t258;
t337 = -t162 - t221;
t183 = t234 * pkin(4) + pkin(10) * t233;
t336 = -t163 - t183;
t237 = -t263 * t284 - t287 * t350;
t238 = t263 * t287 - t317;
t171 = rSges(4,1) * t238 + rSges(4,2) * t237 + rSges(4,3) * t262;
t335 = -t171 - t221;
t200 = Icges(5,4) * t249 - Icges(5,2) * t248 - Icges(5,6) * t351;
t201 = Icges(5,1) * t249 - Icges(5,4) * t248 - Icges(5,5) * t351;
t334 = -t248 * t200 + t249 * t201;
t260 = t281 * t287 - t284 * t354;
t261 = t285 * t352 + t349;
t205 = Icges(4,4) * t261 + Icges(4,2) * t260 - Icges(4,6) * t351;
t206 = Icges(4,1) * t261 + Icges(4,4) * t260 - Icges(4,5) * t351;
t333 = t260 * t205 + t261 * t206;
t202 = rSges(5,1) * t249 - rSges(5,2) * t248 - rSges(5,3) * t351;
t332 = -t202 - t207;
t203 = pkin(4) * t249 + pkin(10) * t248;
t331 = -t203 - t207;
t330 = t221 * t353 + t222 * t350;
t329 = t289 * pkin(1) + pkin(8) * t353;
t105 = Icges(7,5) * t194 + Icges(7,6) * t231 + Icges(7,3) * t193;
t109 = Icges(7,4) * t194 + Icges(7,2) * t231 + Icges(7,6) * t193;
t113 = Icges(7,1) * t194 + Icges(7,4) * t231 + Icges(7,5) * t193;
t40 = t105 * t193 + t109 * t231 + t113 * t194;
t106 = Icges(7,5) * t196 + Icges(7,6) * t233 + Icges(7,3) * t195;
t110 = Icges(7,4) * t196 + Icges(7,2) * t233 + Icges(7,6) * t195;
t114 = Icges(7,1) * t196 + Icges(7,4) * t233 + Icges(7,5) * t195;
t41 = t106 * t193 + t110 * t231 + t114 * t194;
t59 = t141 * t193 + t143 * t231 + t145 * t194;
t1 = t231 * t40 + t233 * t41 + t248 * t59;
t107 = Icges(6,5) * t194 - Icges(6,6) * t193 + Icges(6,3) * t231;
t111 = Icges(6,4) * t194 - Icges(6,2) * t193 + Icges(6,6) * t231;
t115 = Icges(6,1) * t194 - Icges(6,4) * t193 + Icges(6,5) * t231;
t42 = t107 * t231 - t111 * t193 + t115 * t194;
t108 = Icges(6,5) * t196 - Icges(6,6) * t195 + Icges(6,3) * t233;
t112 = Icges(6,4) * t196 - Icges(6,2) * t195 + Icges(6,6) * t233;
t116 = Icges(6,1) * t196 - Icges(6,4) * t195 + Icges(6,5) * t233;
t43 = t108 * t231 - t112 * t193 + t116 * t194;
t60 = t142 * t231 - t144 * t193 + t146 * t194;
t2 = t231 * t42 + t233 * t43 + t248 * t60;
t327 = t2 / 0.2e1 + t1 / 0.2e1;
t44 = t105 * t195 + t109 * t233 + t113 * t196;
t45 = t106 * t195 + t110 * t233 + t114 * t196;
t61 = t141 * t195 + t143 * t233 + t145 * t196;
t3 = t231 * t44 + t233 * t45 + t248 * t61;
t46 = t107 * t233 - t111 * t195 + t115 * t196;
t47 = t108 * t233 - t112 * t195 + t116 * t196;
t62 = t142 * t233 - t144 * t195 + t146 * t196;
t4 = t231 * t46 + t233 * t47 + t248 * t62;
t326 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = t262 * t40 + t264 * t41 - t351 * t59;
t6 = t262 * t42 + t264 * t43 - t351 * t60;
t325 = t5 / 0.2e1 + t6 / 0.2e1;
t7 = t262 * t44 + t264 * t45 - t351 * t61;
t8 = t262 * t46 + t264 * t47 - t351 * t62;
t324 = t7 / 0.2e1 + t8 / 0.2e1;
t10 = t60 * t281 + (t286 * t43 - t289 * t42) * t280;
t9 = t59 * t281 + (t286 * t41 - t289 * t40) * t280;
t323 = t9 / 0.2e1 + t10 / 0.2e1;
t11 = t61 * t281 + (t286 * t45 - t289 * t44) * t280;
t12 = t62 * t281 + (t286 * t47 - t289 * t46) * t280;
t322 = t11 / 0.2e1 + t12 / 0.2e1;
t52 = t105 * t229 + t109 * t248 + t113 * t230;
t53 = t106 * t229 + t110 * t248 + t114 * t230;
t66 = t70 * t248;
t13 = t52 * t231 + t53 * t233 + t66;
t54 = t107 * t248 - t111 * t229 + t115 * t230;
t55 = t108 * t248 - t112 * t229 + t116 * t230;
t67 = t71 * t248;
t14 = t54 * t231 + t55 * t233 + t67;
t321 = t13 / 0.2e1 + t14 / 0.2e1;
t15 = t52 * t262 + t53 * t264 - t351 * t70;
t16 = t54 * t262 + t55 * t264 - t351 * t71;
t320 = t16 / 0.2e1 + t15 / 0.2e1;
t68 = t70 * t281;
t17 = t68 + (t53 * t286 - t52 * t289) * t280;
t69 = t71 * t281;
t18 = t69 + (t55 * t286 - t54 * t289) * t280;
t319 = t17 / 0.2e1 + t18 / 0.2e1;
t120 = t196 * rSges(6,1) - t195 * rSges(6,2) + t233 * rSges(6,3);
t316 = -t120 + t336;
t148 = rSges(6,1) * t230 - rSges(6,2) * t229 + rSges(6,3) * t248;
t315 = -t148 + t331;
t314 = t281 * t183 + t339;
t313 = -t182 + t337;
t239 = -t265 * t284 + t286 * t352;
t240 = t265 * t287 + t318;
t172 = t240 * rSges(4,1) + t239 * rSges(4,2) + t264 * rSges(4,3);
t245 = Icges(3,3) * t281 + (Icges(3,5) * t285 + Icges(3,6) * t288) * t280;
t246 = Icges(3,6) * t281 + (Icges(3,4) * t285 + Icges(3,2) * t288) * t280;
t247 = Icges(3,5) * t281 + (Icges(3,1) * t285 + Icges(3,4) * t288) * t280;
t312 = t281 * t245 + t246 * t351 + t247 * t354;
t216 = t265 * rSges(3,1) - t264 * rSges(3,2) + rSges(3,3) * t353;
t310 = -t286 * pkin(1) + pkin(8) * t350;
t208 = rSges(4,1) * t261 + rSges(4,2) * t260 - rSges(4,3) * t351;
t266 = (pkin(2) * t285 - pkin(9) * t288) * t280;
t308 = t280 * (-t208 - t266);
t307 = t336 - t343;
t306 = t331 - t341;
t305 = t162 * t353 + t163 * t350 + t330;
t304 = t182 * t351 + t262 * t203 + t340;
t303 = t280 * (-t266 + t332);
t302 = -rSges(5,1) * t232 + rSges(5,2) * t231;
t300 = t311 + t329;
t299 = t280 * (-t266 + t315);
t298 = t54 / 0.2e1 + t52 / 0.2e1 + t60 / 0.2e1 + t59 / 0.2e1;
t297 = t61 / 0.2e1 + t62 / 0.2e1 + t55 / 0.2e1 + t53 / 0.2e1;
t296 = t182 * t353 + t183 * t350 + t305;
t295 = t280 * (-t266 + t306);
t294 = -t263 * t277 + t269 + t310;
t118 = rSges(6,1) * t194 - rSges(6,2) * t193 + rSges(6,3) * t231;
t215 = t263 * rSges(3,1) - t262 * rSges(3,2) - rSges(3,3) * t350;
t293 = t183 + t300;
t153 = Icges(5,5) * t232 - Icges(5,6) * t231 + Icges(5,3) * t262;
t155 = Icges(5,4) * t232 - Icges(5,2) * t231 + Icges(5,6) * t262;
t157 = Icges(5,1) * t232 - Icges(5,4) * t231 + Icges(5,5) * t262;
t81 = -t153 * t351 - t155 * t248 + t157 * t249;
t165 = Icges(4,5) * t238 + Icges(4,6) * t237 + Icges(4,3) * t262;
t167 = Icges(4,4) * t238 + Icges(4,2) * t237 + Icges(4,6) * t262;
t169 = Icges(4,1) * t238 + Icges(4,4) * t237 + Icges(4,5) * t262;
t91 = -t165 * t351 + t167 * t260 + t169 * t261;
t199 = Icges(5,5) * t249 - Icges(5,6) * t248 - Icges(5,3) * t351;
t94 = t199 * t262 - t200 * t231 + t201 * t232;
t204 = Icges(4,5) * t261 + Icges(4,6) * t260 - Icges(4,3) * t351;
t98 = t204 * t262 + t205 * t237 + t206 * t238;
t292 = t98 / 0.2e1 + t94 / 0.2e1 + t91 / 0.2e1 + t81 / 0.2e1 + t298;
t154 = Icges(5,5) * t234 - Icges(5,6) * t233 + Icges(5,3) * t264;
t156 = Icges(5,4) * t234 - Icges(5,2) * t233 + Icges(5,6) * t264;
t158 = Icges(5,1) * t234 - Icges(5,4) * t233 + Icges(5,5) * t264;
t82 = -t154 * t351 - t156 * t248 + t158 * t249;
t166 = Icges(4,5) * t240 + Icges(4,6) * t239 + Icges(4,3) * t264;
t168 = Icges(4,4) * t240 + Icges(4,2) * t239 + Icges(4,6) * t264;
t170 = Icges(4,1) * t240 + Icges(4,4) * t239 + Icges(4,5) * t264;
t92 = -t166 * t351 + t168 * t260 + t170 * t261;
t95 = t199 * t264 - t200 * t233 + t201 * t234;
t99 = t204 * t264 + t205 * t239 + t206 * t240;
t291 = t95 / 0.2e1 + t92 / 0.2e1 + t82 / 0.2e1 + t99 / 0.2e1 + t297;
t290 = -t182 + t294 + t355;
t271 = rSges(2,1) * t289 - t286 * rSges(2,2);
t270 = -t286 * rSges(2,1) - rSges(2,2) * t289;
t250 = rSges(3,3) * t281 + (rSges(3,1) * t285 + rSges(3,2) * t288) * t280;
t214 = Icges(3,1) * t265 - Icges(3,4) * t264 + Icges(3,5) * t353;
t213 = Icges(3,1) * t263 - Icges(3,4) * t262 - Icges(3,5) * t350;
t212 = Icges(3,4) * t265 - Icges(3,2) * t264 + Icges(3,6) * t353;
t211 = Icges(3,4) * t263 - Icges(3,2) * t262 - Icges(3,6) * t350;
t210 = Icges(3,5) * t265 - Icges(3,6) * t264 + Icges(3,3) * t353;
t198 = t216 + t329;
t197 = -t215 + t310;
t180 = -t281 * t215 - t250 * t350;
t179 = t216 * t281 - t250 * t353;
t164 = t312 * t281;
t159 = rSges(5,3) * t262 - t302;
t140 = (t215 * t286 + t216 * t289) * t280;
t139 = t245 * t353 - t246 * t264 + t247 * t265;
t138 = -t245 * t350 - t262 * t246 + t263 * t247;
t132 = t222 + t172 + t329;
t131 = t310 + t335;
t126 = -t172 * t351 - t208 * t264;
t125 = t171 * t351 + t208 * t262;
t124 = t210 * t281 + (t212 * t288 + t214 * t285) * t280;
t123 = t209 * t281 + (t211 * t288 + t213 * t285) * t280;
t122 = t300 + t160;
t121 = (-rSges(5,3) + t282) * t262 + t294 + t302;
t104 = -t204 * t351 + t333;
t103 = t104 * t281;
t102 = t171 * t264 - t172 * t262;
t101 = t281 * t335 + t289 * t308;
t100 = t172 * t281 + t286 * t308 + t220;
t97 = -t199 * t351 + t334;
t96 = t97 * t281;
t93 = (t171 * t286 + t172 * t289) * t280 + t330;
t90 = t293 + t120;
t89 = -t118 + t290;
t88 = t120 * t248 - t148 * t233;
t87 = -t118 * t248 + t148 * t231;
t86 = t166 * t264 + t168 * t239 + t170 * t240;
t85 = t165 * t264 + t167 * t239 + t169 * t240;
t84 = t166 * t262 + t168 * t237 + t170 * t238;
t83 = t165 * t262 + t167 * t237 + t169 * t238;
t80 = t264 * t332 + t338 * t351;
t79 = t159 * t351 + t202 * t262 + t340;
t78 = t154 * t264 - t156 * t233 + t158 * t234;
t77 = t153 * t264 - t155 * t233 + t157 * t234;
t76 = t154 * t262 - t156 * t231 + t158 * t232;
t75 = t153 * t262 - t155 * t231 + t157 * t232;
t74 = (-t159 + t337) * t281 + t289 * t303;
t73 = t160 * t281 + t286 * t303 + t339;
t72 = t118 * t233 - t120 * t231;
t65 = t293 + t343;
t64 = t290 - t344;
t63 = t159 * t264 + t262 * t338 + t137;
t58 = (t159 * t286 + t160 * t289) * t280 + t305;
t57 = -t233 * t341 + t248 * t343;
t56 = t231 * t341 - t248 * t344;
t51 = t264 * t315 + t316 * t351;
t50 = t118 * t351 + t148 * t262 + t304;
t49 = (-t118 + t313) * t281 + t289 * t299;
t48 = t120 * t281 + t286 * t299 + t314;
t39 = -t231 * t343 + t233 * t344;
t38 = t118 * t264 + t262 * t316 + t342;
t37 = (t118 * t286 + t120 * t289) * t280 + t296;
t36 = t103 + (t92 * t286 - t91 * t289) * t280;
t35 = -t104 * t351 + t91 * t262 + t92 * t264;
t34 = t264 * t306 + t307 * t351;
t33 = t262 * t341 + t344 * t351 + t304;
t32 = (t313 - t344) * t281 + t289 * t295;
t31 = t281 * t343 + t286 * t295 + t314;
t30 = t99 * t281 + (t286 * t86 - t289 * t85) * t280;
t29 = t98 * t281 + (t286 * t84 - t289 * t83) * t280;
t28 = t96 + (t82 * t286 - t81 * t289) * t280;
t27 = t262 * t85 + t264 * t86 - t351 * t99;
t26 = t262 * t83 + t264 * t84 - t351 * t98;
t25 = t81 * t262 + t82 * t264 - t351 * t97;
t24 = t95 * t281 + (t286 * t78 - t289 * t77) * t280;
t23 = t94 * t281 + (t286 * t76 - t289 * t75) * t280;
t22 = t262 * t307 + t264 * t344 + t342;
t21 = (t286 * t344 + t289 * t343) * t280 + t296;
t20 = t262 * t77 + t264 * t78 - t351 * t95;
t19 = t262 * t75 + t264 * t76 - t351 * t94;
t117 = [t312 + m(7) * (t64 ^ 2 + t65 ^ 2) + m(6) * (t89 ^ 2 + t90 ^ 2) + m(5) * (t121 ^ 2 + t122 ^ 2) + m(4) * (t131 ^ 2 + t132 ^ 2) + m(3) * (t197 ^ 2 + t198 ^ 2) + m(2) * (t270 ^ 2 + t271 ^ 2) + (-t199 - t204) * t351 + Icges(2,3) + t333 + t334 - t361; t96 + t103 + t68 + t69 + t164 + m(7) * (t31 * t65 + t32 * t64) + m(6) * (t48 * t90 + t49 * t89) + m(5) * (t121 * t74 + t122 * t73) + m(4) * (t100 * t132 + t101 * t131) + m(3) * (t179 * t198 + t180 * t197) + ((-t123 / 0.2e1 - t138 / 0.2e1 - t292) * t289 + (t124 / 0.2e1 + t139 / 0.2e1 + t291) * t286) * t280; (t17 + t18 + t28 + t36 + t164) * t281 + m(7) * (t21 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t37 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t58 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(4) * (t100 ^ 2 + t101 ^ 2 + t93 ^ 2) + m(3) * (t140 ^ 2 + t179 ^ 2 + t180 ^ 2) + ((-t9 - t10 - t29 - t23 + ((-t262 * t211 + t263 * t213) * t280 - t360 * t356) * t289) * t289 + (t11 + t12 + t30 + t24 + ((-t212 * t264 + t214 * t265 + (t210 * t286 - t356) * t280) * t286 + (t210 * t350 + t211 * t264 + t262 * t212 - t213 * t265 - t263 * t214) * t289) * t280) * t286 + ((-t123 - t138) * t289 + (t124 + t139) * t286) * t281) * t280; (-t104 - t97 + t361) * t351 + m(7) * (t33 * t64 + t34 * t65) + m(6) * (t50 * t89 + t51 * t90) + m(5) * (t121 * t79 + t122 * t80) + m(4) * (t125 * t131 + t126 * t132) + t291 * t264 + t292 * t262; (t35 / 0.2e1 + t25 / 0.2e1 + t320) * t281 + (t30 / 0.2e1 + t24 / 0.2e1 + t322) * t264 + (t29 / 0.2e1 + t23 / 0.2e1 + t323) * t262 + m(7) * (t21 * t22 + t31 * t34 + t32 * t33) + m(6) * (t37 * t38 + t48 * t51 + t49 * t50) + m(5) * (t63 * t58 + t73 * t80 + t74 * t79) + m(4) * (t100 * t126 + t101 * t125 + t102 * t93) + ((-t19 / 0.2e1 - t26 / 0.2e1 - t325) * t289 + (-t28 / 0.2e1 - t36 / 0.2e1 - t319) * t288 + (t20 / 0.2e1 + t27 / 0.2e1 + t324) * t286) * t280; (-t15 - t16 - t25 - t35) * t351 + (t7 + t8 + t20 + t27) * t264 + (t6 + t5 + t19 + t26) * t262 + m(7) * (t22 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(6) * (t38 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t63 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(4) * (t102 ^ 2 + t125 ^ 2 + t126 ^ 2); m(7) * (t262 * t65 + t264 * t64) + m(6) * (t262 * t90 + t264 * t89) + m(5) * (t121 * t264 + t122 * t262); m(7) * (-t21 * t351 + t262 * t31 + t264 * t32) + m(6) * (t262 * t48 + t264 * t49 - t351 * t37) + m(5) * (t262 * t73 + t264 * t74 - t351 * t58); m(7) * (-t22 * t351 + t262 * t34 + t264 * t33) + m(6) * (t262 * t51 + t264 * t50 - t351 * t38) + m(5) * (t262 * t80 + t264 * t79 - t351 * t63); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t288 ^ 2 * t360 + t262 ^ 2 + t264 ^ 2); t67 + t66 + m(7) * (t56 * t64 + t57 * t65) + m(6) * (t87 * t89 + t88 * t90) + t297 * t233 + t298 * t231; t321 * t281 + t319 * t248 + t322 * t233 + t323 * t231 + m(7) * (t21 * t39 + t31 * t57 + t32 * t56) + m(6) * (t37 * t72 + t48 * t88 + t49 * t87) + (t286 * t326 - t289 * t327) * t280; -t321 * t351 + t326 * t264 + t327 * t262 + t320 * t248 + t324 * t233 + t325 * t231 + m(7) * (t22 * t39 + t33 * t56 + t34 * t57) + m(6) * (t38 * t72 + t50 * t87 + t51 * t88); m(6) * (t262 * t88 + t264 * t87 - t351 * t72) + m(7) * (t262 * t57 + t264 * t56 - t351 * t39); (t14 + t13) * t248 + (t3 + t4) * t233 + (t1 + t2) * t231 + m(7) * (t39 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(6) * (t72 ^ 2 + t87 ^ 2 + t88 ^ 2); m(7) * (t193 * t65 + t195 * t64); m(7) * (t193 * t31 + t195 * t32 + t21 * t229); m(7) * (t193 * t34 + t195 * t33 + t22 * t229); m(7) * (t193 * t262 + t195 * t264 - t229 * t351); m(7) * (t193 * t57 + t195 * t56 + t229 * t39); m(7) * (t193 ^ 2 + t195 ^ 2 + t229 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t117(1) t117(2) t117(4) t117(7) t117(11) t117(16); t117(2) t117(3) t117(5) t117(8) t117(12) t117(17); t117(4) t117(5) t117(6) t117(9) t117(13) t117(18); t117(7) t117(8) t117(9) t117(10) t117(14) t117(19); t117(11) t117(12) t117(13) t117(14) t117(15) t117(20); t117(16) t117(17) t117(18) t117(19) t117(20) t117(21);];
Mq  = res;
