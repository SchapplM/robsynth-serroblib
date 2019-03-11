% Calculate joint inertia matrix for
% S6RRRPRP6
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
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:53:56
% EndTime: 2019-03-09 16:54:12
% DurationCPUTime: 6.73s
% Computational Cost: add. (27343->663), mult. (46203->907), div. (0->0), fcn. (58174->12), ass. (0->305)
t368 = rSges(7,3) + qJ(6) + pkin(10);
t282 = cos(pkin(6));
t288 = sin(qJ(1));
t291 = cos(qJ(2));
t350 = t288 * t291;
t287 = sin(qJ(2));
t292 = cos(qJ(1));
t351 = t287 * t292;
t263 = t282 * t351 + t350;
t332 = qJ(3) + pkin(11);
t279 = sin(t332);
t311 = cos(t332);
t281 = sin(pkin(6));
t354 = t281 * t292;
t233 = t263 * t311 - t279 * t354;
t349 = t291 * t292;
t352 = t287 * t288;
t262 = -t282 * t349 + t352;
t285 = sin(qJ(5));
t289 = cos(qJ(5));
t192 = -t233 * t285 + t262 * t289;
t360 = t262 * t285;
t193 = t233 * t289 + t360;
t302 = t281 * t311;
t232 = t263 * t279 + t292 * t302;
t369 = -rSges(7,1) * t193 - rSges(7,2) * t192 - t368 * t232;
t250 = t282 * t279 + t287 * t302;
t355 = t281 * t291;
t230 = -t250 * t285 - t289 * t355;
t322 = t285 * t355;
t231 = t250 * t289 - t322;
t358 = t281 * t287;
t249 = t279 * t358 - t282 * t311;
t142 = Icges(7,5) * t231 + Icges(7,6) * t230 + Icges(7,3) * t249;
t144 = Icges(7,4) * t231 + Icges(7,2) * t230 + Icges(7,6) * t249;
t146 = Icges(7,1) * t231 + Icges(7,4) * t230 + Icges(7,5) * t249;
t68 = t249 * t142 + t230 * t144 + t231 * t146;
t143 = Icges(6,5) * t231 + Icges(6,6) * t230 + Icges(6,3) * t249;
t145 = Icges(6,4) * t231 + Icges(6,2) * t230 + Icges(6,6) * t249;
t147 = Icges(6,1) * t231 + Icges(6,4) * t230 + Icges(6,5) * t249;
t69 = t249 * t143 + t230 * t145 + t231 * t147;
t367 = -t68 - t69;
t265 = -t282 * t352 + t349;
t357 = t281 * t288;
t235 = t265 * t311 + t279 * t357;
t264 = t282 * t350 + t351;
t194 = -t235 * t285 + t264 * t289;
t359 = t264 * t285;
t195 = t235 * t289 + t359;
t234 = t265 * t279 - t288 * t302;
t277 = pkin(5) * t289 + pkin(4);
t366 = t195 * rSges(7,1) + t194 * rSges(7,2) + pkin(5) * t359 + t234 * t368 + t235 * t277;
t365 = t281 ^ 2;
t290 = cos(qJ(3));
t278 = pkin(3) * t290 + pkin(2);
t364 = -pkin(2) + t278;
t363 = -pkin(4) + t277;
t208 = Icges(3,5) * t263 - Icges(3,6) * t262 - Icges(3,3) * t354;
t362 = t208 * t292;
t284 = -qJ(4) - pkin(9);
t361 = t262 * t284;
t356 = t281 * t290;
t286 = sin(qJ(3));
t353 = t282 * t286;
t228 = t232 * pkin(10);
t348 = pkin(5) * t360 + t233 * t363 - t228 - t369;
t183 = t235 * pkin(4) + pkin(10) * t234;
t347 = -t183 + t366;
t346 = rSges(7,1) * t231 + rSges(7,2) * t230 - pkin(5) * t322 + t363 * t250 + (-pkin(10) + t368) * t249;
t258 = t262 * pkin(9);
t320 = t286 * t354;
t269 = pkin(3) * t320;
t163 = t263 * t364 - t258 - t269 - t361;
t138 = t264 * t163;
t182 = pkin(4) * t233 + t228;
t345 = t264 * t182 + t138;
t206 = pkin(3) * t353 + ((pkin(9) + t284) * t291 + t364 * t287) * t281;
t344 = t163 * t355 + t262 * t206;
t223 = t265 * pkin(2) + pkin(9) * t264;
t321 = t286 * t357;
t313 = pkin(3) * t321 - t264 * t284 + t265 * t278;
t164 = -t223 + t313;
t220 = t282 * t223;
t343 = t282 * t164 + t220;
t161 = t235 * rSges(5,1) - t234 * rSges(5,2) + t264 * rSges(5,3);
t342 = -t161 - t164;
t222 = pkin(2) * t263 + t258;
t341 = -t163 - t222;
t340 = -t164 - t183;
t238 = -t263 * t286 - t290 * t354;
t239 = t263 * t290 - t320;
t172 = rSges(4,1) * t239 + rSges(4,2) * t238 + rSges(4,3) * t262;
t339 = -t172 - t222;
t199 = Icges(5,4) * t250 - Icges(5,2) * t249 - Icges(5,6) * t355;
t200 = Icges(5,1) * t250 - Icges(5,4) * t249 - Icges(5,5) * t355;
t338 = -t249 * t199 + t250 * t200;
t260 = t282 * t290 - t286 * t358;
t261 = t287 * t356 + t353;
t204 = Icges(4,4) * t261 + Icges(4,2) * t260 - Icges(4,6) * t355;
t205 = Icges(4,1) * t261 + Icges(4,4) * t260 - Icges(4,5) * t355;
t337 = t260 * t204 + t261 * t205;
t201 = rSges(5,1) * t250 - rSges(5,2) * t249 - rSges(5,3) * t355;
t336 = -t201 - t206;
t202 = pkin(4) * t250 + pkin(10) * t249;
t335 = -t202 - t206;
t334 = t222 * t357 + t223 * t354;
t333 = t292 * pkin(1) + pkin(8) * t357;
t107 = Icges(7,5) * t193 + Icges(7,6) * t192 + Icges(7,3) * t232;
t111 = Icges(7,4) * t193 + Icges(7,2) * t192 + Icges(7,6) * t232;
t115 = Icges(7,1) * t193 + Icges(7,4) * t192 + Icges(7,5) * t232;
t42 = t107 * t232 + t111 * t192 + t115 * t193;
t108 = Icges(7,5) * t195 + Icges(7,6) * t194 + Icges(7,3) * t234;
t112 = Icges(7,4) * t195 + Icges(7,2) * t194 + Icges(7,6) * t234;
t116 = Icges(7,1) * t195 + Icges(7,4) * t194 + Icges(7,5) * t234;
t43 = t108 * t232 + t112 * t192 + t116 * t193;
t59 = t142 * t232 + t144 * t192 + t146 * t193;
t1 = t232 * t42 + t234 * t43 + t249 * t59;
t109 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t232;
t113 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t232;
t117 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t232;
t44 = t109 * t232 + t113 * t192 + t117 * t193;
t110 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t234;
t114 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t234;
t118 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t234;
t45 = t110 * t232 + t114 * t192 + t118 * t193;
t60 = t143 * t232 + t145 * t192 + t147 * t193;
t2 = t232 * t44 + t234 * t45 + t249 * t60;
t331 = t2 / 0.2e1 + t1 / 0.2e1;
t46 = t107 * t234 + t111 * t194 + t115 * t195;
t47 = t108 * t234 + t112 * t194 + t116 * t195;
t61 = t142 * t234 + t144 * t194 + t146 * t195;
t3 = t232 * t46 + t234 * t47 + t249 * t61;
t48 = t109 * t234 + t113 * t194 + t117 * t195;
t49 = t110 * t234 + t114 * t194 + t118 * t195;
t62 = t143 * t234 + t145 * t194 + t147 * t195;
t4 = t232 * t48 + t234 * t49 + t249 * t62;
t330 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = t262 * t42 + t264 * t43 - t355 * t59;
t6 = t262 * t44 + t264 * t45 - t355 * t60;
t329 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t262 * t46 + t264 * t47 - t355 * t61;
t8 = t262 * t48 + t264 * t49 - t355 * t62;
t328 = t7 / 0.2e1 + t8 / 0.2e1;
t10 = t60 * t282 + (t288 * t45 - t292 * t44) * t281;
t9 = t59 * t282 + (t288 * t43 - t292 * t42) * t281;
t327 = t10 / 0.2e1 + t9 / 0.2e1;
t11 = t61 * t282 + (t288 * t47 - t292 * t46) * t281;
t12 = t62 * t282 + (t288 * t49 - t292 * t48) * t281;
t326 = t12 / 0.2e1 + t11 / 0.2e1;
t54 = t107 * t249 + t111 * t230 + t115 * t231;
t55 = t108 * t249 + t112 * t230 + t116 * t231;
t64 = t68 * t249;
t13 = t54 * t232 + t55 * t234 + t64;
t56 = t109 * t249 + t113 * t230 + t117 * t231;
t57 = t110 * t249 + t114 * t230 + t118 * t231;
t65 = t69 * t249;
t14 = t56 * t232 + t57 * t234 + t65;
t325 = -t14 / 0.2e1 - t13 / 0.2e1;
t15 = t54 * t262 + t55 * t264 - t355 * t68;
t16 = t56 * t262 + t57 * t264 - t355 * t69;
t324 = t16 / 0.2e1 + t15 / 0.2e1;
t66 = t68 * t282;
t17 = t66 + (t55 * t288 - t54 * t292) * t281;
t67 = t69 * t282;
t18 = t67 + (t57 * t288 - t56 * t292) * t281;
t323 = t17 / 0.2e1 + t18 / 0.2e1;
t122 = t195 * rSges(6,1) + t194 * rSges(6,2) + t234 * rSges(6,3);
t319 = -t122 + t340;
t149 = rSges(6,1) * t231 + rSges(6,2) * t230 + rSges(6,3) * t249;
t318 = -t149 + t335;
t317 = t282 * t183 + t343;
t316 = -t182 + t341;
t240 = -t265 * t286 + t288 * t356;
t241 = t265 * t290 + t321;
t173 = t241 * rSges(4,1) + t240 * rSges(4,2) + t264 * rSges(4,3);
t246 = Icges(3,3) * t282 + (Icges(3,5) * t287 + Icges(3,6) * t291) * t281;
t247 = Icges(3,6) * t282 + (Icges(3,4) * t287 + Icges(3,2) * t291) * t281;
t248 = Icges(3,5) * t282 + (Icges(3,1) * t287 + Icges(3,4) * t291) * t281;
t314 = t282 * t246 + t247 * t355 + t248 * t358;
t215 = t265 * rSges(3,1) - t264 * rSges(3,2) + rSges(3,3) * t357;
t312 = -t288 * pkin(1) + pkin(8) * t354;
t207 = rSges(4,1) * t261 + rSges(4,2) * t260 - rSges(4,3) * t355;
t266 = (pkin(2) * t287 - pkin(9) * t291) * t281;
t310 = t281 * (-t207 - t266);
t309 = t340 - t347;
t308 = t335 - t346;
t307 = t163 * t357 + t164 * t354 + t334;
t306 = t182 * t355 + t262 * t202 + t344;
t305 = t281 * (-t266 + t336);
t304 = -rSges(5,1) * t233 + rSges(5,2) * t232;
t301 = t313 + t333;
t300 = t281 * (-t266 + t318);
t299 = t54 / 0.2e1 + t56 / 0.2e1 + t59 / 0.2e1 + t60 / 0.2e1;
t298 = t55 / 0.2e1 + t57 / 0.2e1 + t62 / 0.2e1 + t61 / 0.2e1;
t297 = t182 * t357 + t183 * t354 + t307;
t296 = t281 * (-t266 + t308);
t295 = -t263 * t278 + t269 + t312;
t120 = rSges(6,1) * t193 + rSges(6,2) * t192 + rSges(6,3) * t232;
t214 = t263 * rSges(3,1) - t262 * rSges(3,2) - rSges(3,3) * t354;
t154 = Icges(5,5) * t233 - Icges(5,6) * t232 + Icges(5,3) * t262;
t156 = Icges(5,4) * t233 - Icges(5,2) * t232 + Icges(5,6) * t262;
t158 = Icges(5,1) * t233 - Icges(5,4) * t232 + Icges(5,5) * t262;
t79 = -t154 * t355 - t156 * t249 + t158 * t250;
t166 = Icges(4,5) * t239 + Icges(4,6) * t238 + Icges(4,3) * t262;
t168 = Icges(4,4) * t239 + Icges(4,2) * t238 + Icges(4,6) * t262;
t170 = Icges(4,1) * t239 + Icges(4,4) * t238 + Icges(4,5) * t262;
t91 = -t166 * t355 + t168 * t260 + t170 * t261;
t198 = Icges(5,5) * t250 - Icges(5,6) * t249 - Icges(5,3) * t355;
t94 = t198 * t262 - t199 * t232 + t200 * t233;
t203 = Icges(4,5) * t261 + Icges(4,6) * t260 - Icges(4,3) * t355;
t98 = t203 * t262 + t204 * t238 + t205 * t239;
t294 = t98 / 0.2e1 + t94 / 0.2e1 + t91 / 0.2e1 + t79 / 0.2e1 + t299;
t155 = Icges(5,5) * t235 - Icges(5,6) * t234 + Icges(5,3) * t264;
t157 = Icges(5,4) * t235 - Icges(5,2) * t234 + Icges(5,6) * t264;
t159 = Icges(5,1) * t235 - Icges(5,4) * t234 + Icges(5,5) * t264;
t80 = -t155 * t355 - t157 * t249 + t159 * t250;
t167 = Icges(4,5) * t241 + Icges(4,6) * t240 + Icges(4,3) * t264;
t169 = Icges(4,4) * t241 + Icges(4,2) * t240 + Icges(4,6) * t264;
t171 = Icges(4,1) * t241 + Icges(4,4) * t240 + Icges(4,5) * t264;
t92 = -t167 * t355 + t169 * t260 + t171 * t261;
t95 = t198 * t264 - t199 * t234 + t200 * t235;
t99 = t203 * t264 + t204 * t240 + t205 * t241;
t293 = t95 / 0.2e1 + t92 / 0.2e1 + t80 / 0.2e1 + t99 / 0.2e1 + t298;
t271 = rSges(2,1) * t292 - t288 * rSges(2,2);
t270 = -t288 * rSges(2,1) - rSges(2,2) * t292;
t251 = rSges(3,3) * t282 + (rSges(3,1) * t287 + rSges(3,2) * t291) * t281;
t213 = Icges(3,1) * t265 - Icges(3,4) * t264 + Icges(3,5) * t357;
t212 = Icges(3,1) * t263 - Icges(3,4) * t262 - Icges(3,5) * t354;
t211 = Icges(3,4) * t265 - Icges(3,2) * t264 + Icges(3,6) * t357;
t210 = Icges(3,4) * t263 - Icges(3,2) * t262 - Icges(3,6) * t354;
t209 = Icges(3,5) * t265 - Icges(3,6) * t264 + Icges(3,3) * t357;
t197 = t215 + t333;
t196 = -t214 + t312;
t181 = -t282 * t214 - t251 * t354;
t180 = t215 * t282 - t251 * t357;
t165 = t314 * t282;
t160 = rSges(5,3) * t262 - t304;
t141 = (t214 * t288 + t215 * t292) * t281;
t140 = t246 * t357 - t247 * t264 + t248 * t265;
t139 = -t246 * t354 - t262 * t247 + t263 * t248;
t134 = t223 + t173 + t333;
t133 = t312 + t339;
t128 = -t173 * t355 - t207 * t264;
t127 = t172 * t355 + t207 * t262;
t126 = t209 * t282 + (t211 * t291 + t213 * t287) * t281;
t125 = t208 * t282 + (t210 * t291 + t212 * t287) * t281;
t124 = t301 + t161;
t123 = (-rSges(5,3) + t284) * t262 + t295 + t304;
t104 = -t203 * t355 + t337;
t103 = t104 * t282;
t102 = t172 * t264 - t173 * t262;
t101 = t282 * t339 + t292 * t310;
t100 = t173 * t282 + t288 * t310 + t220;
t97 = -t198 * t355 + t338;
t96 = t97 * t282;
t93 = (t172 * t288 + t173 * t292) * t281 + t334;
t90 = t183 + t301 + t122;
t89 = -t120 - t182 + t295 + t361;
t88 = t122 * t249 - t149 * t234;
t87 = -t120 * t249 + t149 * t232;
t86 = t301 + t366;
t85 = -t233 * t277 + (-pkin(5) * t285 + t284) * t262 + t295 + t369;
t84 = t167 * t264 + t169 * t240 + t171 * t241;
t83 = t166 * t264 + t168 * t240 + t170 * t241;
t82 = t167 * t262 + t169 * t238 + t171 * t239;
t81 = t166 * t262 + t168 * t238 + t170 * t239;
t78 = t264 * t336 + t342 * t355;
t77 = t160 * t355 + t201 * t262 + t344;
t76 = t155 * t264 - t157 * t234 + t159 * t235;
t75 = t154 * t264 - t156 * t234 + t158 * t235;
t74 = t155 * t262 - t157 * t232 + t159 * t233;
t73 = t154 * t262 - t156 * t232 + t158 * t233;
t72 = (-t160 + t341) * t282 + t292 * t305;
t71 = t161 * t282 + t288 * t305 + t343;
t70 = t120 * t234 - t122 * t232;
t63 = t160 * t264 + t262 * t342 + t138;
t58 = (t160 * t288 + t161 * t292) * t281 + t307;
t53 = t264 * t318 + t319 * t355;
t52 = t120 * t355 + t149 * t262 + t306;
t51 = (-t120 + t316) * t282 + t292 * t300;
t50 = t122 * t282 + t288 * t300 + t317;
t41 = -t234 * t346 + t249 * t347;
t40 = t232 * t346 - t249 * t348;
t39 = t120 * t264 + t262 * t319 + t345;
t38 = (t120 * t288 + t122 * t292) * t281 + t297;
t37 = t103 + (t92 * t288 - t91 * t292) * t281;
t36 = -t104 * t355 + t91 * t262 + t92 * t264;
t35 = -t232 * t347 + t234 * t348;
t34 = t99 * t282 + (t288 * t84 - t292 * t83) * t281;
t33 = t98 * t282 + (t288 * t82 - t292 * t81) * t281;
t32 = t264 * t308 + t309 * t355;
t31 = t262 * t346 + t348 * t355 + t306;
t30 = t96 + (t80 * t288 - t79 * t292) * t281;
t29 = (t316 - t348) * t282 + t292 * t296;
t28 = t282 * t347 + t288 * t296 + t317;
t27 = t262 * t83 + t264 * t84 - t355 * t99;
t26 = t262 * t81 + t264 * t82 - t355 * t98;
t25 = t79 * t262 + t80 * t264 - t355 * t97;
t24 = t95 * t282 + (t288 * t76 - t292 * t75) * t281;
t23 = t94 * t282 + (t288 * t74 - t292 * t73) * t281;
t22 = t262 * t75 + t264 * t76 - t355 * t95;
t21 = t262 * t73 + t264 * t74 - t355 * t94;
t20 = t262 * t309 + t264 * t348 + t345;
t19 = (t288 * t348 + t292 * t347) * t281 + t297;
t105 = [t314 + m(7) * (t85 ^ 2 + t86 ^ 2) + m(6) * (t89 ^ 2 + t90 ^ 2) + m(5) * (t123 ^ 2 + t124 ^ 2) + m(4) * (t133 ^ 2 + t134 ^ 2) + m(3) * (t196 ^ 2 + t197 ^ 2) + m(2) * (t270 ^ 2 + t271 ^ 2) + (-t198 - t203) * t355 + Icges(2,3) + t337 + t338 - t367; t67 + t96 + t103 + t66 + t165 + m(3) * (t180 * t197 + t181 * t196) + m(4) * (t100 * t134 + t101 * t133) + m(7) * (t28 * t86 + t29 * t85) + m(6) * (t50 * t90 + t51 * t89) + m(5) * (t123 * t72 + t124 * t71) + ((-t125 / 0.2e1 - t139 / 0.2e1 - t294) * t292 + (t126 / 0.2e1 + t140 / 0.2e1 + t293) * t288) * t281; (t17 + t18 + t37 + t30 + t165) * t282 + m(7) * (t19 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t38 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t58 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(4) * (t100 ^ 2 + t101 ^ 2 + t93 ^ 2) + m(3) * (t141 ^ 2 + t180 ^ 2 + t181 ^ 2) + ((-t9 - t10 - t33 - t23 + ((-t262 * t210 + t263 * t212) * t281 - t365 * t362) * t292) * t292 + (t12 + t11 + t34 + t24 + ((-t264 * t211 + t265 * t213 + (t209 * t288 - t362) * t281) * t288 + (t209 * t354 + t264 * t210 + t262 * t211 - t265 * t212 - t263 * t213) * t292) * t281) * t288 + ((-t125 - t139) * t292 + (t126 + t140) * t288) * t282) * t281; (-t104 - t97 + t367) * t355 + m(7) * (t31 * t85 + t32 * t86) + m(6) * (t52 * t89 + t53 * t90) + m(5) * (t123 * t77 + t124 * t78) + m(4) * (t127 * t133 + t128 * t134) + t293 * t264 + t294 * t262; (t25 / 0.2e1 + t36 / 0.2e1 + t324) * t282 + (t24 / 0.2e1 + t34 / 0.2e1 + t326) * t264 + (t23 / 0.2e1 + t33 / 0.2e1 + t327) * t262 + m(4) * (t100 * t128 + t101 * t127 + t102 * t93) + m(7) * (t19 * t20 + t28 * t32 + t29 * t31) + m(6) * (t38 * t39 + t50 * t53 + t51 * t52) + m(5) * (t58 * t63 + t71 * t78 + t72 * t77) + ((-t21 / 0.2e1 - t26 / 0.2e1 - t329) * t292 + (-t30 / 0.2e1 - t37 / 0.2e1 - t323) * t291 + (t22 / 0.2e1 + t27 / 0.2e1 + t328) * t288) * t281; (-t15 - t16 - t25 - t36) * t355 + (t8 + t7 + t22 + t27) * t264 + (t6 + t5 + t26 + t21) * t262 + m(7) * (t20 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t39 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(5) * (t63 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(4) * (t102 ^ 2 + t127 ^ 2 + t128 ^ 2); m(7) * (t262 * t86 + t264 * t85) + m(6) * (t262 * t90 + t264 * t89) + m(5) * (t123 * t264 + t124 * t262); m(7) * (-t19 * t355 + t262 * t28 + t264 * t29) + m(6) * (t262 * t50 + t264 * t51 - t355 * t38) + m(5) * (t262 * t71 + t264 * t72 - t355 * t58); m(7) * (-t20 * t355 + t262 * t32 + t264 * t31) + m(6) * (t262 * t53 + t264 * t52 - t355 * t39) + m(5) * (t262 * t78 + t264 * t77 - t355 * t63); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t291 ^ 2 * t365 + t262 ^ 2 + t264 ^ 2); t64 + t65 + m(7) * (t40 * t85 + t41 * t86) + m(6) * (t87 * t89 + t88 * t90) + t298 * t234 + t299 * t232; -t325 * t282 + t323 * t249 + t326 * t234 + t327 * t232 + m(7) * (t19 * t35 + t28 * t41 + t29 * t40) + m(6) * (t38 * t70 + t50 * t88 + t51 * t87) + (t288 * t330 - t292 * t331) * t281; t325 * t355 + t330 * t264 + t331 * t262 + t324 * t249 + t328 * t234 + t329 * t232 + m(7) * (t20 * t35 + t31 * t40 + t32 * t41) + m(6) * (t39 * t70 + t52 * t87 + t53 * t88); m(6) * (t262 * t88 + t264 * t87 - t355 * t70) + m(7) * (t262 * t41 + t264 * t40 - t35 * t355); (t13 + t14) * t249 + (t4 + t3) * t234 + (t1 + t2) * t232 + m(7) * (t35 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t70 ^ 2 + t87 ^ 2 + t88 ^ 2); m(7) * (t232 * t86 + t234 * t85); m(7) * (t19 * t249 + t232 * t28 + t234 * t29); m(7) * (t20 * t249 + t232 * t32 + t234 * t31); m(7) * (t232 * t262 + t234 * t264 - t249 * t355); m(7) * (t232 * t41 + t234 * t40 + t249 * t35); m(7) * (t232 ^ 2 + t234 ^ 2 + t249 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t105(1) t105(2) t105(4) t105(7) t105(11) t105(16); t105(2) t105(3) t105(5) t105(8) t105(12) t105(17); t105(4) t105(5) t105(6) t105(9) t105(13) t105(18); t105(7) t105(8) t105(9) t105(10) t105(14) t105(19); t105(11) t105(12) t105(13) t105(14) t105(15) t105(20); t105(16) t105(17) t105(18) t105(19) t105(20) t105(21);];
Mq  = res;
