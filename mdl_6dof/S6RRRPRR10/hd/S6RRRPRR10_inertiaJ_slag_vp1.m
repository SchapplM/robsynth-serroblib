% Calculate joint inertia matrix for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR10_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR10_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:15:41
% EndTime: 2019-03-09 19:15:54
% DurationCPUTime: 5.41s
% Computational Cost: add. (13609->601), mult. (26940->833), div. (0->0), fcn. (32450->10), ass. (0->287)
t259 = sin(qJ(2));
t368 = Icges(3,5) * t259;
t367 = t368 / 0.2e1;
t258 = sin(qJ(3));
t262 = cos(qJ(3));
t263 = cos(qJ(2));
t201 = -Icges(4,3) * t263 + (Icges(4,5) * t262 - Icges(4,6) * t258) * t259;
t204 = -Icges(5,2) * t263 + (Icges(5,4) * t262 + Icges(5,6) * t258) * t259;
t366 = -t201 - t204;
t264 = cos(qJ(1));
t260 = sin(qJ(1));
t334 = t260 * t263;
t227 = t258 * t334 + t262 * t264;
t228 = -t258 * t264 + t262 * t334;
t256 = qJ(5) + qJ(6);
t250 = sin(t256);
t251 = cos(t256);
t165 = t227 * t251 - t228 * t250;
t166 = t227 * t250 + t228 * t251;
t289 = -t166 * rSges(7,1) - t165 * rSges(7,2);
t338 = t259 * t260;
t112 = -rSges(7,3) * t338 - t289;
t245 = pkin(9) * t338;
t257 = sin(qJ(5));
t342 = t227 * t257;
t316 = pkin(5) * t342;
t265 = -pkin(10) - pkin(9);
t335 = t259 * t265;
t261 = cos(qJ(5));
t249 = pkin(5) * t261 + pkin(4);
t352 = -pkin(4) + t249;
t332 = t228 * t352 + t260 * t335 + t112 + t245 + t316;
t365 = t264 * t332;
t205 = -Icges(4,6) * t263 + (Icges(4,4) * t262 - Icges(4,2) * t258) * t259;
t339 = t258 * t259;
t200 = -Icges(5,6) * t263 + (Icges(5,5) * t262 + Icges(5,3) * t258) * t259;
t208 = -Icges(5,4) * t263 + (Icges(5,1) * t262 + Icges(5,5) * t258) * t259;
t209 = -Icges(4,5) * t263 + (Icges(4,1) * t262 - Icges(4,4) * t258) * t259;
t337 = t259 * t262;
t363 = t200 * t339 + (t208 + t209) * t337;
t364 = (-t205 * t339 + t263 * t366 + t363) * t263;
t196 = (-t250 * t262 + t251 * t258) * t259;
t197 = (t250 * t258 + t251 * t262) * t259;
t140 = rSges(7,1) * t197 + rSges(7,2) * t196 + rSges(7,3) * t263;
t340 = t257 * t258;
t172 = (-pkin(9) - t265) * t263 + (pkin(5) * t340 + t262 * t352) * t259;
t330 = t140 + t172;
t54 = -t263 * t332 - t330 * t338;
t362 = t260 ^ 2;
t361 = t264 ^ 2;
t344 = Icges(7,3) * t263;
t345 = Icges(7,6) * t196;
t346 = Icges(7,5) * t197;
t276 = t344 + t345 + t346;
t271 = t259 * t276;
t278 = Icges(7,1) * t197 + Icges(7,4) * t196 + Icges(7,5) * t263;
t274 = t166 * t278;
t277 = Icges(7,4) * t197 + Icges(7,2) * t196 + Icges(7,6) * t263;
t275 = t165 * t277;
t269 = -t260 * t271 + t274 + t275;
t106 = t166 * Icges(7,5) + t165 * Icges(7,6) - Icges(7,3) * t338;
t108 = t166 * Icges(7,4) + t165 * Icges(7,2) - Icges(7,6) * t338;
t110 = t166 * Icges(7,1) + t165 * Icges(7,4) - Icges(7,5) * t338;
t33 = -t106 * t338 + t108 * t165 + t110 * t166;
t333 = t263 * t264;
t229 = t258 * t333 - t260 * t262;
t230 = t258 * t260 + t262 * t333;
t167 = t229 * t251 - t230 * t250;
t168 = t229 * t250 + t230 * t251;
t336 = t259 * t264;
t107 = t168 * Icges(7,5) + t167 * Icges(7,6) - Icges(7,3) * t336;
t109 = t168 * Icges(7,4) + t167 * Icges(7,2) - Icges(7,6) * t336;
t111 = t168 * Icges(7,1) + t167 * Icges(7,4) - Icges(7,5) * t336;
t34 = -t107 * t338 + t109 * t165 + t111 * t166;
t7 = -t269 * t263 + (t260 * t33 + t264 * t34) * t259;
t177 = t227 * t261 - t228 * t257;
t178 = t228 * t261 + t342;
t116 = Icges(6,5) * t178 + Icges(6,6) * t177 - Icges(6,3) * t338;
t118 = Icges(6,4) * t178 + Icges(6,2) * t177 - Icges(6,6) * t338;
t120 = Icges(6,1) * t178 + Icges(6,4) * t177 - Icges(6,5) * t338;
t39 = -t116 * t338 + t118 * t177 + t120 * t178;
t179 = t229 * t261 - t230 * t257;
t341 = t229 * t257;
t180 = t230 * t261 + t341;
t117 = Icges(6,5) * t180 + Icges(6,6) * t179 - Icges(6,3) * t336;
t119 = Icges(6,4) * t180 + Icges(6,2) * t179 - Icges(6,6) * t336;
t121 = Icges(6,1) * t180 + Icges(6,4) * t179 - Icges(6,5) * t336;
t40 = -t117 * t338 + t119 * t177 + t121 * t178;
t215 = (-t257 * t262 + t258 * t261) * t259;
t216 = (t261 * t262 + t340) * t259;
t143 = Icges(6,5) * t216 + Icges(6,6) * t215 + Icges(6,3) * t263;
t144 = Icges(6,4) * t216 + Icges(6,2) * t215 + Icges(6,6) * t263;
t145 = Icges(6,1) * t216 + Icges(6,4) * t215 + Icges(6,5) * t263;
t64 = -t143 * t338 + t144 * t177 + t145 * t178;
t9 = -t64 * t263 + (t260 * t39 + t264 * t40) * t259;
t359 = t7 + t9;
t358 = -t260 / 0.2e1;
t357 = t263 / 0.2e1;
t356 = t264 / 0.2e1;
t41 = -t116 * t336 + t118 * t179 + t120 * t180;
t42 = -t117 * t336 + t119 * t179 + t121 * t180;
t65 = -t143 * t336 + t144 * t179 + t145 * t180;
t10 = -t65 * t263 + (t260 * t41 + t264 * t42) * t259;
t270 = t264 * t271;
t272 = t168 * t278;
t273 = t167 * t277;
t268 = t272 + t273 - t270;
t35 = -t106 * t336 + t108 * t167 + t110 * t168;
t36 = -t107 * t336 + t109 * t167 + t111 * t168;
t8 = -t268 * t263 + (t260 * t35 + t264 * t36) * t259;
t355 = t10 + t8;
t237 = rSges(3,1) * t259 + rSges(3,2) * t263;
t354 = m(3) * t237;
t353 = pkin(2) * t263;
t351 = rSges(5,3) * t227;
t49 = t106 * t263 + t108 * t196 + t110 * t197;
t50 = t107 * t263 + t109 * t196 + t111 * t197;
t310 = t196 * t277 + t197 * t278 + t263 * t276;
t66 = t310 * t263;
t14 = -t66 + (t49 * t260 + t50 * t264) * t259;
t13 = t263 * t14;
t52 = t116 * t263 + t118 * t215 + t120 * t216;
t53 = t117 * t263 + t119 * t215 + t121 * t216;
t308 = t263 * t143 + t215 * t144 + t216 * t145;
t72 = t308 * t263;
t15 = -t72 + (t52 * t260 + t53 * t264) * t259;
t350 = t263 * t15;
t349 = t264 * rSges(3,3);
t113 = t168 * rSges(7,1) + t167 * rSges(7,2) - rSges(7,3) * t336;
t85 = t263 * t113 + t140 * t336;
t347 = Icges(3,4) * t263;
t290 = -rSges(6,1) * t178 - rSges(6,2) * t177;
t122 = -rSges(6,3) * t338 - t290;
t343 = t122 * t264;
t163 = t230 * rSges(5,1) + rSges(5,2) * t336 + t229 * rSges(5,3);
t184 = t230 * pkin(3) + qJ(4) * t229;
t329 = -t163 - t184;
t220 = t227 * qJ(4);
t183 = t228 * pkin(3) + t220;
t169 = t183 * t336;
t194 = t228 * pkin(4) - t245;
t328 = t194 * t336 + t169;
t327 = t180 * rSges(6,1) + t179 * rSges(6,2);
t231 = (pkin(3) * t262 + qJ(4) * t258) * t259;
t326 = t263 * t183 + t231 * t338;
t225 = t230 * pkin(4);
t195 = -pkin(9) * t336 + t225;
t325 = -t184 - t195;
t214 = -rSges(4,3) * t263 + (rSges(4,1) * t262 - rSges(4,2) * t258) * t259;
t240 = pkin(2) * t259 - pkin(8) * t263;
t323 = -t214 - t240;
t319 = pkin(2) * t333 + pkin(8) * t336;
t322 = t362 * (pkin(8) * t259 + t353) + t264 * t319;
t232 = pkin(4) * t337 + t263 * pkin(9);
t321 = -t231 - t232;
t320 = -t231 - t240;
t318 = t264 * pkin(1) + t260 * pkin(7);
t254 = t264 * pkin(7);
t317 = t254 - t220;
t315 = -t52 / 0.2e1 - t64 / 0.2e1;
t314 = -t53 / 0.2e1 - t65 / 0.2e1;
t149 = Icges(5,5) * t228 + Icges(5,6) * t338 + Icges(5,3) * t227;
t153 = Icges(5,4) * t228 + Icges(5,2) * t338 + Icges(5,6) * t227;
t157 = Icges(5,1) * t228 + Icges(5,4) * t338 + Icges(5,5) * t227;
t88 = -t153 * t263 + (t149 * t258 + t157 * t262) * t259;
t151 = Icges(4,5) * t228 - Icges(4,6) * t227 + Icges(4,3) * t338;
t155 = Icges(4,4) * t228 - Icges(4,2) * t227 + Icges(4,6) * t338;
t159 = Icges(4,1) * t228 - Icges(4,4) * t227 + Icges(4,5) * t338;
t90 = -t151 * t263 + (-t155 * t258 + t159 * t262) * t259;
t313 = t90 / 0.2e1 + t88 / 0.2e1;
t150 = Icges(5,5) * t230 + Icges(5,6) * t336 + Icges(5,3) * t229;
t154 = Icges(5,4) * t230 + Icges(5,2) * t336 + Icges(5,6) * t229;
t158 = Icges(5,1) * t230 + Icges(5,4) * t336 + Icges(5,5) * t229;
t89 = -t154 * t263 + (t150 * t258 + t158 * t262) * t259;
t152 = Icges(4,5) * t230 - Icges(4,6) * t229 + Icges(4,3) * t336;
t156 = Icges(4,4) * t230 - Icges(4,2) * t229 + Icges(4,6) * t336;
t160 = Icges(4,1) * t230 - Icges(4,4) * t229 + Icges(4,5) * t336;
t91 = -t152 * t263 + (-t156 * t258 + t160 * t262) * t259;
t312 = -t91 / 0.2e1 - t89 / 0.2e1;
t123 = -rSges(6,3) * t336 + t327;
t311 = -t123 + t325;
t307 = pkin(5) * t341 + t230 * t249 + t264 * t335;
t129 = -t195 + t307;
t309 = -t129 + t325;
t213 = -rSges(5,2) * t263 + (rSges(5,1) * t262 + rSges(5,3) * t258) * t259;
t306 = -t213 + t320;
t164 = t230 * rSges(4,1) - t229 * rSges(4,2) + rSges(4,3) * t336;
t305 = -t232 + t320;
t304 = -pkin(1) - t353;
t303 = -t338 / 0.2e1;
t302 = -t336 / 0.2e1;
t146 = rSges(6,1) * t216 + rSges(6,2) * t215 + rSges(6,3) * t263;
t299 = -t146 + t305;
t298 = t260 * t183 + t264 * t184 + t322;
t297 = t263 * t194 + t232 * t338 + t326;
t296 = t318 + t319;
t18 = t260 * t34 - t264 * t33;
t19 = t260 * t36 - t264 * t35;
t23 = t50 * t260 - t49 * t264;
t295 = t18 * t303 + t19 * t302 + t23 * t357 + t7 * t356 + t8 * t358;
t294 = t66 + (t269 + t49) * t303 + (t268 + t50) * t302;
t293 = 0.2e1 * t8 * t302 + 0.2e1 * t7 * t303 + t13;
t292 = rSges(3,1) * t263 - rSges(3,2) * t259;
t291 = -rSges(4,1) * t228 + rSges(4,2) * t227;
t288 = t305 - t330;
t286 = -Icges(3,2) * t259 + t347;
t285 = Icges(3,5) * t263 - Icges(3,6) * t259;
t282 = rSges(3,1) * t333 - rSges(3,2) * t336 + t260 * rSges(3,3);
t281 = t260 * t194 + t264 * t195 + t298;
t92 = -t122 * t263 - t146 * t338;
t280 = (t260 * t7 + t264 * t8) * t259 - t13;
t279 = t184 + t296;
t101 = t200 * t227 + t204 * t338 + t208 * t228;
t102 = t201 * t338 - t205 * t227 + t209 * t228;
t267 = t101 / 0.2e1 + t102 / 0.2e1 + t49 / 0.2e1 + t275 / 0.2e1 + t274 / 0.2e1 + t313 - t315;
t103 = t200 * t229 + t204 * t336 + t208 * t230;
t104 = t201 * t336 - t205 * t229 + t209 * t230;
t266 = t103 / 0.2e1 + t104 / 0.2e1 + t50 / 0.2e1 + t273 / 0.2e1 + t272 / 0.2e1 - t312 - t314;
t239 = rSges(2,1) * t264 - rSges(2,2) * t260;
t238 = -rSges(2,1) * t260 - rSges(2,2) * t264;
t234 = Icges(3,6) * t263 + t368;
t203 = Icges(3,3) * t260 + t264 * t285;
t202 = -Icges(3,3) * t264 + t260 * t285;
t190 = t282 + t318;
t189 = t349 + t254 + (-pkin(1) - t292) * t260;
t182 = t323 * t264;
t181 = t323 * t260;
t162 = rSges(4,3) * t338 - t291;
t161 = rSges(5,1) * t228 + rSges(5,2) * t338 + t351;
t142 = t264 * t282 + (t260 * t292 - t349) * t260;
t139 = t306 * t264;
t138 = t306 * t260;
t133 = t296 + t164;
t132 = t254 + ((-rSges(4,3) - pkin(8)) * t259 + t304) * t260 + t291;
t131 = -t164 * t263 - t214 * t336;
t130 = t162 * t263 + t214 * t338;
t115 = t299 * t264;
t114 = t299 * t260;
t105 = (t162 * t264 - t164 * t260) * t259;
t100 = t279 + t163;
t99 = -t351 + (-rSges(5,1) - pkin(3)) * t228 + ((-rSges(5,2) - pkin(8)) * t259 + t304) * t260 + t317;
t97 = t113 * t338;
t96 = t162 * t260 + t164 * t264 + t322;
t95 = t329 * t263 + (-t213 - t231) * t336;
t94 = t161 * t263 + t213 * t338 + t326;
t93 = t123 * t263 + t146 * t336;
t87 = t288 * t264;
t86 = t288 * t260;
t84 = -t112 * t263 - t140 * t338;
t83 = t152 * t336 - t156 * t229 + t160 * t230;
t82 = t151 * t336 - t155 * t229 + t159 * t230;
t81 = t150 * t229 + t154 * t336 + t158 * t230;
t80 = t149 * t229 + t153 * t336 + t157 * t230;
t79 = t152 * t338 - t156 * t227 + t160 * t228;
t78 = t151 * t338 - t155 * t227 + t159 * t228;
t77 = t150 * t227 + t154 * t338 + t158 * t228;
t76 = t149 * t227 + t153 * t338 + t157 * t228;
t75 = t225 + (-rSges(6,3) - pkin(9)) * t336 + t279 + t327;
t74 = t245 + (-pkin(3) - pkin(4)) * t228 + ((rSges(6,3) - pkin(8)) * t259 + t304) * t260 + t290 + t317;
t73 = t169 + (t161 * t264 + t260 * t329) * t259;
t71 = (t123 * t260 - t343) * t259;
t70 = t113 + t279 + t307;
t69 = -t316 + (-pkin(3) - t249) * t228 + ((rSges(7,3) - pkin(8) - t265) * t259 + t304) * t260 + t289 + t317;
t68 = -t112 * t336 + t97;
t67 = t161 * t260 + t163 * t264 + t298;
t61 = t311 * t263 + (-t146 + t321) * t336;
t60 = -t92 + t297;
t55 = t129 * t263 + t172 * t336 + t85;
t51 = (t260 * t311 + t343) * t259 + t328;
t48 = t260 * t83 - t264 * t82;
t47 = t260 * t81 - t264 * t80;
t46 = t260 * t79 - t264 * t78;
t45 = t260 * t77 - t264 * t76;
t38 = t122 * t260 + t123 * t264 + t281;
t37 = t97 + (t129 * t260 - t365) * t259;
t32 = t309 * t263 + (-t172 + t321) * t336 - t85;
t31 = t297 - t54;
t30 = -t104 * t263 + (t260 * t82 + t264 * t83) * t259;
t29 = -t103 * t263 + (t260 * t80 + t264 * t81) * t259;
t28 = -t102 * t263 + (t260 * t78 + t264 * t79) * t259;
t27 = -t101 * t263 + (t260 * t76 + t264 * t77) * t259;
t26 = -t97 + (t260 * t309 + t365) * t259 + t328;
t25 = (t113 + t129) * t264 + t332 * t260 + t281;
t24 = t53 * t260 - t52 * t264;
t21 = t260 * t42 - t264 * t41;
t20 = t260 * t40 - t264 * t39;
t1 = [Icges(2,3) + (Icges(3,1) * t259 - t258 * t205 + t347) * t259 + (Icges(3,4) * t259 + Icges(3,2) * t263 + t366) * t263 + m(7) * (t69 ^ 2 + t70 ^ 2) + m(6) * (t74 ^ 2 + t75 ^ 2) + m(5) * (t100 ^ 2 + t99 ^ 2) + m(4) * (t132 ^ 2 + t133 ^ 2) + m(3) * (t189 ^ 2 + t190 ^ 2) + m(2) * (t238 ^ 2 + t239 ^ 2) + t308 + t310 + t363; m(7) * (t69 * t87 + t70 * t86) + m(6) * (t114 * t75 + t115 * t74) + m(5) * (t100 * t138 + t139 * t99) + m(4) * (t132 * t182 + t133 * t181) + (-t189 * t354 - (-Icges(3,6) * t264 + t260 * t286) * t263 / 0.2e1 + t264 * t367 + t234 * t356 - t267) * t264 + (t264 * t286 * t357 - t190 * t354 + t266 + (Icges(3,6) * t357 + t367 + t234 / 0.2e1) * t260) * t260; m(7) * (t25 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(6) * (t114 ^ 2 + t115 ^ 2 + t38 ^ 2) + m(5) * (t138 ^ 2 + t139 ^ 2 + t67 ^ 2) + m(4) * (t181 ^ 2 + t182 ^ 2 + t96 ^ 2) + m(3) * (t142 ^ 2 + (t361 + t362) * t237 ^ 2) + (-t361 * t202 - t18 - t20 - t45 - t46) * t264 + (t362 * t203 + t19 + t21 + t47 + t48 + (-t260 * t202 + t264 * t203) * t264) * t260; -t66 - t72 - t364 + m(7) * (t31 * t69 + t32 * t70) + m(6) * (t60 * t74 + t61 * t75) + m(5) * (t100 * t95 + t94 * t99) + m(4) * (t130 * t132 + t131 * t133) + ((-t270 / 0.2e1 + t266) * t264 + ((-t346 / 0.2e1 - t345 / 0.2e1 - t344 / 0.2e1) * t338 + t267) * t260) * t259; (-t23 / 0.2e1 - t24 / 0.2e1) * t263 + (-t7 / 0.2e1 - t9 / 0.2e1 - t28 / 0.2e1 - t27 / 0.2e1 + t313 * t263) * t264 + (t8 / 0.2e1 + t10 / 0.2e1 + t30 / 0.2e1 + t29 / 0.2e1 + t312 * t263) * t260 + m(7) * (t26 * t25 + t31 * t87 + t32 * t86) + m(6) * (t114 * t61 + t115 * t60 + t51 * t38) + m(5) * (t138 * t95 + t139 * t94 + t67 * t73) + m(4) * (t105 * t96 + t130 * t182 + t131 * t181) + ((t19 / 0.2e1 + t21 / 0.2e1 + t48 / 0.2e1 + t47 / 0.2e1) * t264 + (t18 / 0.2e1 + t20 / 0.2e1 + t46 / 0.2e1 + t45 / 0.2e1) * t260) * t259; (-t14 - t15 + t364) * t263 + m(7) * (t26 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t51 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t73 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(4) * (t105 ^ 2 + t130 ^ 2 + t131 ^ 2) + ((t29 + t30 + (-t89 - t91) * t263 + t355) * t264 + (t27 + t28 + (-t88 - t90) * t263 + t359) * t260) * t259; m(7) * (t227 * t70 + t229 * t69) + m(6) * (t227 * t75 + t229 * t74) + m(5) * (t100 * t227 + t229 * t99); m(7) * (t227 * t86 + t229 * t87 + t25 * t339) + m(6) * (t114 * t227 + t115 * t229 + t339 * t38) + m(5) * (t138 * t227 + t139 * t229 + t339 * t67); m(7) * (t227 * t32 + t229 * t31 + t26 * t339) + m(6) * (t227 * t61 + t229 * t60 + t339 * t51) + m(5) * (t227 * t95 + t229 * t94 + t339 * t73); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t258 ^ 2 * t259 ^ 2 + t227 ^ 2 + t229 ^ 2); t72 + m(7) * (t54 * t69 + t55 * t70) + m(6) * (t74 * t92 + t75 * t93) + (t260 * t315 + t264 * t314) * t259 + t294; t9 * t356 + t10 * t358 + t24 * t357 + (t20 * t358 - t264 * t21 / 0.2e1) * t259 + m(7) * (t37 * t25 + t54 * t87 + t55 * t86) + m(6) * (t114 * t93 + t115 * t92 + t38 * t71) + t295; t350 + (-t264 * t10 - t260 * t9) * t259 + m(7) * (t26 * t37 + t31 * t54 + t32 * t55) + m(6) * (t51 * t71 + t60 * t92 + t61 * t93) + t293; m(6) * (t227 * t93 + t229 * t92 + t339 * t71) + m(7) * (t227 * t55 + t229 * t54 + t339 * t37); -t350 - t13 + m(7) * (t37 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(6) * (t71 ^ 2 + t92 ^ 2 + t93 ^ 2) + (t260 * t359 + t264 * t355) * t259; m(7) * (t69 * t84 + t70 * t85) + t294; m(7) * (t68 * t25 + t84 * t87 + t85 * t86) + t295; m(7) * (t68 * t26 + t31 * t84 + t32 * t85) + t293; m(7) * (t227 * t85 + t229 * t84 + t339 * t68); m(7) * (t68 * t37 + t54 * t84 + t55 * t85) + t280; m(7) * (t68 ^ 2 + t84 ^ 2 + t85 ^ 2) + t280;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
