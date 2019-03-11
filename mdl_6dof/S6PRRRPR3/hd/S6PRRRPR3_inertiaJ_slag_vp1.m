% Calculate joint inertia matrix for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:11:04
% EndTime: 2019-03-08 23:11:21
% DurationCPUTime: 7.64s
% Computational Cost: add. (24767->506), mult. (42315->724), div. (0->0), fcn. (53111->12), ass. (0->248)
t363 = Icges(5,1) + Icges(6,2);
t362 = Icges(6,1) + Icges(5,3);
t361 = -Icges(5,4) - Icges(6,6);
t360 = Icges(6,4) - Icges(5,5);
t359 = Icges(6,5) - Icges(5,6);
t358 = Icges(5,2) + Icges(6,3);
t260 = sin(pkin(11));
t262 = cos(pkin(11));
t269 = cos(qJ(2));
t263 = cos(pkin(6));
t266 = sin(qJ(2));
t315 = t263 * t266;
t249 = t260 * t269 + t262 * t315;
t261 = sin(pkin(6));
t313 = qJ(3) + qJ(4);
t286 = cos(t313);
t278 = t261 * t286;
t285 = sin(t313);
t231 = t249 * t285 + t262 * t278;
t277 = t261 * t285;
t232 = t249 * t286 - t262 * t277;
t314 = t263 * t269;
t248 = t260 * t266 - t262 * t314;
t357 = t358 * t231 + t361 * t232 + t359 * t248;
t251 = -t260 * t315 + t262 * t269;
t233 = t251 * t285 - t260 * t278;
t234 = t251 * t286 + t260 * t277;
t250 = t260 * t314 + t262 * t266;
t356 = t358 * t233 + t361 * t234 + t359 * t250;
t355 = t359 * t231 - t360 * t232 + t362 * t248;
t354 = t359 * t233 - t360 * t234 + t362 * t250;
t353 = t361 * t231 + t363 * t232 - t360 * t248;
t352 = t361 * t233 + t363 * t234 - t360 * t250;
t246 = -t263 * t286 + t266 * t277;
t247 = t263 * t285 + t266 * t278;
t317 = t261 * t269;
t351 = t358 * t246 + t361 * t247 - t359 * t317;
t350 = t361 * t246 + t363 * t247 + t360 * t317;
t349 = t359 * t246 - t360 * t247 - t362 * t317;
t348 = t357 * t231 + t353 * t232 + t355 * t248;
t347 = t356 * t231 + t352 * t232 + t354 * t248;
t346 = t357 * t233 + t353 * t234 + t355 * t250;
t345 = t356 * t233 + t352 * t234 + t354 * t250;
t344 = t357 * t246 + t353 * t247 - t355 * t317;
t343 = t356 * t246 + t352 * t247 - t354 * t317;
t342 = t351 * t231 + t350 * t232 + t349 * t248;
t341 = t351 * t233 + t350 * t234 + t349 * t250;
t340 = t351 * t246 + t350 * t247 - t349 * t317;
t332 = m(6) + m(7);
t264 = sin(qJ(6));
t267 = cos(qJ(6));
t197 = t231 * t267 - t248 * t264;
t198 = t231 * t264 + t248 * t267;
t136 = rSges(7,1) * t198 + rSges(7,2) * t197 + rSges(7,3) * t232;
t312 = pkin(5) * t248 + pkin(10) * t232 + t136;
t235 = t246 * t267 + t264 * t317;
t236 = t246 * t264 - t267 * t317;
t154 = rSges(7,1) * t236 + rSges(7,2) * t235 + rSges(7,3) * t247;
t339 = -pkin(5) * t317 + pkin(10) * t247 + t154;
t130 = Icges(7,5) * t198 + Icges(7,6) * t197 + Icges(7,3) * t232;
t132 = Icges(7,4) * t198 + Icges(7,2) * t197 + Icges(7,6) * t232;
t134 = Icges(7,1) * t198 + Icges(7,4) * t197 + Icges(7,5) * t232;
t66 = t130 * t232 + t132 * t197 + t134 * t198;
t199 = t233 * t267 - t250 * t264;
t200 = t233 * t264 + t250 * t267;
t131 = Icges(7,5) * t200 + Icges(7,6) * t199 + Icges(7,3) * t234;
t133 = Icges(7,4) * t200 + Icges(7,2) * t199 + Icges(7,6) * t234;
t135 = Icges(7,1) * t200 + Icges(7,4) * t199 + Icges(7,5) * t234;
t67 = t131 * t232 + t133 * t197 + t135 * t198;
t151 = Icges(7,5) * t236 + Icges(7,6) * t235 + Icges(7,3) * t247;
t152 = Icges(7,4) * t236 + Icges(7,2) * t235 + Icges(7,6) * t247;
t153 = Icges(7,1) * t236 + Icges(7,4) * t235 + Icges(7,5) * t247;
t80 = t151 * t232 + t152 * t197 + t153 * t198;
t11 = t248 * t66 + t250 * t67 - t317 * t80;
t338 = t348 * t248 + t347 * t250 - t342 * t317 + t11;
t68 = t130 * t234 + t132 * t199 + t134 * t200;
t69 = t131 * t234 + t133 * t199 + t135 * t200;
t81 = t151 * t234 + t152 * t199 + t153 * t200;
t12 = t248 * t68 + t250 * t69 - t317 * t81;
t337 = t346 * t248 + t345 * t250 - t341 * t317 + t12;
t15 = t80 * t263 + (t260 * t67 - t262 * t66) * t261;
t336 = t15 + t342 * t263 + (t347 * t260 - t348 * t262) * t261;
t16 = t81 * t263 + (t260 * t69 - t262 * t68) * t261;
t335 = t16 + t341 * t263 + (t345 * t260 - t346 * t262) * t261;
t71 = t130 * t247 + t132 * t235 + t134 * t236;
t72 = t131 * t247 + t133 * t235 + t135 * t236;
t85 = t151 * t247 + t152 * t235 + t153 * t236;
t21 = t248 * t71 + t250 * t72 - t317 * t85;
t334 = t344 * t248 + t343 * t250 - t340 * t317 + t21;
t23 = t85 * t263 + (t260 * t72 - t262 * t71) * t261;
t333 = t23 + t340 * t263 + (t343 * t260 - t344 * t262) * t261;
t331 = t232 / 0.2e1;
t330 = t234 / 0.2e1;
t329 = t247 / 0.2e1;
t328 = t248 / 0.2e1;
t327 = t250 / 0.2e1;
t326 = t260 / 0.2e1;
t325 = -t262 / 0.2e1;
t324 = t263 / 0.2e1;
t268 = cos(qJ(3));
t323 = pkin(3) * t268;
t321 = t260 * t261;
t320 = t261 * t262;
t265 = sin(qJ(3));
t319 = t261 * t265;
t318 = t261 * t268;
t316 = t263 * t265;
t137 = rSges(7,1) * t200 + rSges(7,2) * t199 + rSges(7,3) * t234;
t311 = pkin(5) * t250 + pkin(10) * t234 + t137;
t168 = rSges(6,1) * t248 - rSges(6,2) * t232 + rSges(6,3) * t231;
t189 = pkin(4) * t232 + qJ(5) * t231;
t174 = t250 * t189;
t310 = t250 * t168 + t174;
t297 = t262 * t319;
t172 = -pkin(3) * t297 + pkin(9) * t248 + t249 * t323;
t225 = pkin(3) * t316 + (-pkin(9) * t269 + t266 * t323) * t261;
t309 = t172 * t317 + t248 * t225;
t298 = t260 * t319;
t173 = pkin(3) * t298 + pkin(9) * t250 + t251 * t323;
t230 = t251 * pkin(2) + t250 * pkin(8);
t228 = t263 * t230;
t308 = t263 * t173 + t228;
t169 = rSges(6,1) * t250 - rSges(6,2) * t234 + rSges(6,3) * t233;
t190 = pkin(4) * t234 + qJ(5) * t233;
t307 = -t169 - t190;
t171 = rSges(5,1) * t234 - rSges(5,2) * t233 + rSges(5,3) * t250;
t306 = -t171 - t173;
t229 = t249 * pkin(2) + t248 * pkin(8);
t305 = -t172 - t229;
t223 = pkin(4) * t247 + qJ(5) * t246;
t304 = t189 * t317 + t248 * t223;
t170 = rSges(5,1) * t232 - rSges(5,2) * t231 + rSges(5,3) * t248;
t210 = rSges(5,1) * t247 - rSges(5,2) * t246 - rSges(5,3) * t317;
t127 = t170 * t317 + t248 * t210;
t209 = -rSges(6,1) * t317 - rSges(6,2) * t247 + rSges(6,3) * t246;
t303 = -t209 - t223;
t302 = -t210 - t225;
t301 = t229 * t321 + t230 * t320;
t296 = t312 * t250 + t174;
t295 = -t190 - t311;
t294 = t263 * t190 + t308;
t293 = -t223 - t339;
t292 = -t173 + t307;
t291 = -t189 + t305;
t290 = -t225 + t303;
t287 = -t317 / 0.2e1;
t252 = t263 * t268 - t266 * t319;
t253 = t266 * t318 + t316;
t224 = rSges(4,1) * t253 + rSges(4,2) * t252 - rSges(4,3) * t317;
t254 = (pkin(2) * t266 - pkin(8) * t269) * t261;
t284 = (-t224 - t254) * t261;
t283 = -t173 + t295;
t282 = t172 * t321 + t173 * t320 + t301;
t281 = -t225 + t293;
t106 = t168 * t317 + t248 * t209 + t304;
t280 = (-t254 + t302) * t261;
t18 = t232 * t71 + t234 * t72 + t247 * t85;
t3 = t232 * t66 + t234 * t67 + t247 * t80;
t4 = t232 * t68 + t234 * t69 + t247 * t81;
t279 = t11 * t331 + t12 * t330 + t18 * t287 + t21 * t329 + t3 * t328 + t4 * t327;
t276 = t338 * t248 + t337 * t250;
t275 = (-t254 + t290) * t261;
t274 = t189 * t321 + t190 * t320 + t282;
t73 = t339 * t248 + t312 * t317 + t304;
t273 = (-t254 + t281) * t261;
t272 = -t317 * t334 + t276;
t271 = t336 * t328 + t335 * t327 + t334 * t324 + t337 * t321 / 0.2e1 - t338 * t320 / 0.2e1 + t333 * t287;
t245 = t263 * rSges(3,3) + (rSges(3,1) * t266 + rSges(3,2) * t269) * t261;
t244 = Icges(3,5) * t263 + (Icges(3,1) * t266 + Icges(3,4) * t269) * t261;
t243 = Icges(3,6) * t263 + (Icges(3,4) * t266 + Icges(3,2) * t269) * t261;
t242 = Icges(3,3) * t263 + (Icges(3,5) * t266 + Icges(3,6) * t269) * t261;
t240 = t251 * t268 + t298;
t239 = -t251 * t265 + t260 * t318;
t238 = t249 * t268 - t297;
t237 = -t249 * t265 - t262 * t318;
t222 = Icges(4,1) * t253 + Icges(4,4) * t252 - Icges(4,5) * t317;
t221 = Icges(4,4) * t253 + Icges(4,2) * t252 - Icges(4,6) * t317;
t220 = Icges(4,5) * t253 + Icges(4,6) * t252 - Icges(4,3) * t317;
t219 = rSges(3,1) * t251 - rSges(3,2) * t250 + rSges(3,3) * t321;
t218 = rSges(3,1) * t249 - rSges(3,2) * t248 - rSges(3,3) * t320;
t217 = Icges(3,1) * t251 - Icges(3,4) * t250 + Icges(3,5) * t321;
t216 = Icges(3,1) * t249 - Icges(3,4) * t248 - Icges(3,5) * t320;
t215 = Icges(3,4) * t251 - Icges(3,2) * t250 + Icges(3,6) * t321;
t214 = Icges(3,4) * t249 - Icges(3,2) * t248 - Icges(3,6) * t320;
t213 = Icges(3,5) * t251 - Icges(3,6) * t250 + Icges(3,3) * t321;
t212 = Icges(3,5) * t249 - Icges(3,6) * t248 - Icges(3,3) * t320;
t187 = -t218 * t263 - t245 * t320;
t186 = t219 * t263 - t245 * t321;
t182 = rSges(4,1) * t240 + rSges(4,2) * t239 + rSges(4,3) * t250;
t181 = rSges(4,1) * t238 + rSges(4,2) * t237 + rSges(4,3) * t248;
t180 = Icges(4,1) * t240 + Icges(4,4) * t239 + Icges(4,5) * t250;
t179 = Icges(4,1) * t238 + Icges(4,4) * t237 + Icges(4,5) * t248;
t178 = Icges(4,4) * t240 + Icges(4,2) * t239 + Icges(4,6) * t250;
t177 = Icges(4,4) * t238 + Icges(4,2) * t237 + Icges(4,6) * t248;
t176 = Icges(4,5) * t240 + Icges(4,6) * t239 + Icges(4,3) * t250;
t175 = Icges(4,5) * t238 + Icges(4,6) * t237 + Icges(4,3) * t248;
t144 = (t218 * t260 + t219 * t262) * t261;
t143 = t250 * t172;
t142 = t250 * t170;
t139 = -t182 * t317 - t224 * t250;
t138 = t181 * t317 + t224 * t248;
t128 = -t171 * t317 - t250 * t210;
t125 = -t220 * t317 + t221 * t252 + t222 * t253;
t124 = t181 * t250 - t182 * t248;
t123 = (-t181 - t229) * t263 + t262 * t284;
t122 = t263 * t182 + t260 * t284 + t228;
t119 = t220 * t250 + t221 * t239 + t222 * t240;
t118 = t220 * t248 + t221 * t237 + t222 * t238;
t117 = -t171 * t248 + t142;
t112 = (t181 * t260 + t182 * t262) * t261 + t301;
t111 = -t176 * t317 + t178 * t252 + t180 * t253;
t110 = -t175 * t317 + t177 * t252 + t179 * t253;
t109 = t137 * t247 - t154 * t234;
t108 = -t136 * t247 + t154 * t232;
t107 = t250 * t303 + t307 * t317;
t101 = t250 * t302 + t306 * t317;
t100 = t127 + t309;
t99 = t176 * t250 + t178 * t239 + t180 * t240;
t98 = t175 * t250 + t177 * t239 + t179 * t240;
t97 = t176 * t248 + t178 * t237 + t180 * t238;
t96 = t175 * t248 + t177 * t237 + t179 * t238;
t87 = (-t170 + t305) * t263 + t262 * t280;
t86 = t263 * t171 + t260 * t280 + t308;
t84 = t136 * t234 - t137 * t232;
t83 = t248 * t307 + t310;
t82 = t248 * t306 + t142 + t143;
t79 = t250 * t290 + t292 * t317;
t78 = t106 + t309;
t77 = (t170 * t260 + t171 * t262) * t261 + t282;
t76 = (-t168 + t291) * t263 + t262 * t275;
t75 = t263 * t169 + t260 * t275 + t294;
t74 = t250 * t293 + t295 * t317;
t70 = t248 * t292 + t143 + t310;
t65 = (t168 * t260 + t169 * t262) * t261 + t274;
t64 = t248 * t295 + t296;
t63 = t250 * t281 + t283 * t317;
t62 = t73 + t309;
t61 = (t291 - t312) * t263 + t262 * t273;
t60 = t260 * t273 + t263 * t311 + t294;
t59 = t125 * t263 + (-t110 * t262 + t111 * t260) * t261;
t58 = t110 * t248 + t111 * t250 - t125 * t317;
t53 = t248 * t283 + t143 + t296;
t52 = (t260 * t312 + t262 * t311) * t261 + t274;
t51 = t119 * t263 + (t260 * t99 - t262 * t98) * t261;
t50 = t118 * t263 + (t260 * t97 - t262 * t96) * t261;
t45 = -t119 * t317 + t248 * t98 + t250 * t99;
t44 = -t118 * t317 + t248 * t96 + t250 * t97;
t1 = [m(2) + m(3) + m(4) + m(5) + t332; m(3) * t144 + m(4) * t112 + m(5) * t77 + m(6) * t65 + m(7) * t52; m(4) * (t112 ^ 2 + t122 ^ 2 + t123 ^ 2) + m(3) * (t144 ^ 2 + t186 ^ 2 + t187 ^ 2) + m(7) * (t52 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(6) * (t65 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(5) * (t77 ^ 2 + t86 ^ 2 + t87 ^ 2) + (t51 + (t213 * t321 - t215 * t250 + t217 * t251) * t321 + t335) * t321 + (-t50 + (-t212 * t320 - t214 * t248 + t216 * t249) * t320 + (-t212 * t321 + t213 * t320 + t214 * t250 + t215 * t248 - t216 * t251 - t217 * t249) * t321 - t336) * t320 + (t59 + ((t215 * t269 + t217 * t266) * t260 - (t214 * t269 + t216 * t266) * t262) * t261 ^ 2 - (-t242 * t320 - t248 * t243 + t249 * t244) * t320 + (t242 * t321 - t250 * t243 + t251 * t244) * t321 + ((-t212 * t262 + t213 * t260 + t243 * t269 + t244 * t266) * t261 + t263 * t242) * t263 + t333) * t263; m(4) * t124 + m(5) * t82 + m(6) * t70 + m(7) * t53; (-t269 * t59 / 0.2e1 + t44 * t325 + t45 * t326) * t261 + t271 + m(7) * (t52 * t53 + t60 * t63 + t61 * t62) + m(6) * (t65 * t70 + t75 * t79 + t76 * t78) + m(5) * (t100 * t87 + t101 * t86 + t77 * t82) + m(4) * (t112 * t124 + t122 * t139 + t123 * t138) + t50 * t328 + t51 * t327 + t58 * t324; t248 * t44 + t250 * t45 + (-t58 - t334) * t317 + m(5) * (t100 ^ 2 + t101 ^ 2 + t82 ^ 2) + m(4) * (t124 ^ 2 + t138 ^ 2 + t139 ^ 2) + m(7) * (t53 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(6) * (t70 ^ 2 + t78 ^ 2 + t79 ^ 2) + t276; m(5) * t117 + m(6) * t83 + m(7) * t64; m(7) * (t52 * t64 + t60 * t74 + t61 * t73) + m(6) * (t106 * t76 + t107 * t75 + t65 * t83) + m(5) * (t117 * t77 + t127 * t87 + t128 * t86) + t271; m(7) * (t53 * t64 + t62 * t73 + t63 * t74) + m(6) * (t106 * t78 + t107 * t79 + t70 * t83) + m(5) * (t100 * t127 + t101 * t128 + t117 * t82) + t272; m(7) * (t64 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(6) * (t106 ^ 2 + t107 ^ 2 + t83 ^ 2) + m(5) * (t117 ^ 2 + t127 ^ 2 + t128 ^ 2) + t272; t246 * t332; m(7) * (t231 * t60 + t233 * t61 + t246 * t52) + m(6) * (t231 * t75 + t233 * t76 + t246 * t65); m(7) * (t231 * t63 + t233 * t62 + t246 * t53) + m(6) * (t231 * t79 + t233 * t78 + t246 * t70); m(7) * (t231 * t74 + t233 * t73 + t246 * t64) + m(6) * (t106 * t233 + t107 * t231 + t246 * t83); (t231 ^ 2 + t233 ^ 2 + t246 ^ 2) * t332; m(7) * t84; t15 * t331 + t18 * t324 + t23 * t329 + t16 * t330 + m(7) * (t108 * t61 + t109 * t60 + t52 * t84) + (t3 * t325 + t326 * t4) * t261; m(7) * (t108 * t62 + t109 * t63 + t53 * t84) + t279; m(7) * (t108 * t73 + t109 * t74 + t64 * t84) + t279; m(7) * (t108 * t233 + t109 * t231 + t246 * t84); t234 * t4 + t232 * t3 + t247 * t18 + m(7) * (t108 ^ 2 + t109 ^ 2 + t84 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
