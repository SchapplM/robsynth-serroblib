% Calculate joint inertia matrix for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:20:19
% EndTime: 2019-03-09 00:20:38
% DurationCPUTime: 8.06s
% Computational Cost: add. (60619->664), mult. (167408->938), div. (0->0), fcn. (221084->14), ass. (0->310)
t357 = rSges(7,3) + qJ(6);
t254 = sin(pkin(12));
t256 = cos(pkin(12));
t262 = sin(qJ(2));
t257 = cos(pkin(6));
t264 = cos(qJ(2));
t334 = t257 * t264;
t246 = -t254 * t262 + t256 * t334;
t255 = sin(pkin(6));
t343 = cos(pkin(7));
t298 = t255 * t343;
t342 = sin(pkin(7));
t235 = -t246 * t342 - t256 * t298;
t335 = t257 * t262;
t247 = t254 * t264 + t256 * t335;
t221 = t247 * pkin(2) + t235 * pkin(9);
t248 = -t254 * t334 - t256 * t262;
t236 = -t248 * t342 + t254 * t298;
t249 = -t254 * t335 + t256 * t264;
t222 = t249 * pkin(2) + t236 * pkin(9);
t337 = t255 * t256;
t338 = t254 * t255;
t318 = t221 * t338 + t222 * t337;
t261 = sin(qJ(3));
t350 = cos(qJ(3));
t275 = t350 * t342;
t273 = t255 * t275;
t276 = t343 * t350;
t216 = -t248 * t276 + t249 * t261 - t254 * t273;
t297 = t255 * t342;
t217 = t249 * t350 + (t343 * t248 + t254 * t297) * t261;
t180 = rSges(4,1) * t217 - rSges(4,2) * t216 + rSges(4,3) * t236;
t260 = sin(qJ(4));
t349 = cos(qJ(4));
t201 = t217 * t260 - t236 * t349;
t202 = t217 * t349 + t236 * t260;
t147 = rSges(5,1) * t202 - rSges(5,2) * t201 + rSges(5,3) * t216;
t259 = sin(qJ(5));
t263 = cos(qJ(5));
t168 = -t202 * t259 + t216 * t263;
t340 = t216 * t259;
t169 = t202 * t263 + t340;
t124 = rSges(6,1) * t169 + rSges(6,2) * t168 + rSges(6,3) * t201;
t345 = pkin(5) * t263;
t332 = rSges(7,1) * t169 + rSges(7,2) * t168 + pkin(5) * t340 + t201 * t357 + t345 * t202;
t268 = -m(6) * t124 - m(7) * t332;
t266 = -m(5) * t147 + t268;
t356 = -m(4) * t180 + t266;
t355 = m(6) / 0.2e1 + m(7) / 0.2e1;
t354 = -0.2e1 * t235;
t353 = m(5) / 0.2e1;
t214 = -t246 * t276 + t247 * t261 + t256 * t273;
t341 = t214 * t259;
t336 = t255 * t262;
t233 = -t255 * t264 * t276 - t257 * t275 + t261 * t336;
t339 = t233 * t259;
t215 = t247 * t350 + (t343 * t246 - t256 * t297) * t261;
t200 = t215 * t349 + t235 * t260;
t166 = -t200 * t259 + t214 * t263;
t167 = t200 * t263 + t341;
t199 = t215 * t260 - t235 * t349;
t333 = rSges(7,1) * t167 + rSges(7,2) * t166 + pkin(5) * t341 + t199 * t357 + t345 * t200;
t122 = rSges(6,1) * t167 + rSges(6,2) * t166 + rSges(6,3) * t199;
t163 = t200 * pkin(4) + t199 * pkin(11);
t331 = -t122 - t163;
t164 = t202 * pkin(4) + t201 * pkin(11);
t330 = -t124 - t164;
t234 = t257 * t342 * t261 + (t343 * t261 * t264 + t350 * t262) * t255;
t245 = t257 * t343 - t264 * t297;
t219 = t234 * t349 + t245 * t260;
t197 = -t219 * t259 + t233 * t263;
t198 = t219 * t263 + t339;
t218 = t234 * t260 - t245 * t349;
t329 = rSges(7,1) * t198 + rSges(7,2) * t197 + pkin(5) * t339 + t218 * t357 + t345 * t219;
t146 = rSges(5,1) * t200 - rSges(5,2) * t199 + rSges(5,3) * t214;
t192 = pkin(3) * t215 + pkin(10) * t214;
t328 = -t146 - t192;
t149 = rSges(6,1) * t198 + rSges(6,2) * t197 + rSges(6,3) * t218;
t194 = t219 * pkin(4) + t218 * pkin(11);
t327 = -t149 - t194;
t156 = t236 * t163;
t184 = t236 * t192;
t325 = t156 + t184;
t193 = pkin(3) * t217 + pkin(10) * t216;
t186 = t245 * t193;
t324 = t245 * t164 + t186;
t181 = rSges(5,1) * t219 - rSges(5,2) * t218 + rSges(5,3) * t233;
t209 = pkin(3) * t234 + pkin(10) * t233;
t323 = -t181 - t209;
t322 = t193 * t354 + 0.2e1 * t184;
t196 = t235 * t209;
t321 = t235 * t194 + t196;
t220 = t257 * t222;
t320 = t257 * t193 + t220;
t319 = 0.2e1 * t318;
t109 = Icges(7,5) * t167 + Icges(7,6) * t166 + Icges(7,3) * t199;
t113 = Icges(7,4) * t167 + Icges(7,2) * t166 + Icges(7,6) * t199;
t117 = Icges(7,1) * t167 + Icges(7,4) * t166 + Icges(7,5) * t199;
t48 = t109 * t199 + t113 * t166 + t117 * t167;
t110 = Icges(7,5) * t169 + Icges(7,6) * t168 + Icges(7,3) * t201;
t114 = Icges(7,4) * t169 + Icges(7,2) * t168 + Icges(7,6) * t201;
t118 = Icges(7,1) * t169 + Icges(7,4) * t168 + Icges(7,5) * t201;
t49 = t110 * t199 + t114 * t166 + t118 * t167;
t134 = Icges(7,5) * t198 + Icges(7,6) * t197 + Icges(7,3) * t218;
t138 = Icges(7,4) * t198 + Icges(7,2) * t197 + Icges(7,6) * t218;
t142 = Icges(7,1) * t198 + Icges(7,4) * t197 + Icges(7,5) * t218;
t71 = t134 * t199 + t138 * t166 + t142 * t167;
t1 = t199 * t48 + t201 * t49 + t218 * t71;
t111 = Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t199;
t115 = Icges(6,4) * t167 + Icges(6,2) * t166 + Icges(6,6) * t199;
t119 = Icges(6,1) * t167 + Icges(6,4) * t166 + Icges(6,5) * t199;
t50 = t111 * t199 + t115 * t166 + t119 * t167;
t112 = Icges(6,5) * t169 + Icges(6,6) * t168 + Icges(6,3) * t201;
t116 = Icges(6,4) * t169 + Icges(6,2) * t168 + Icges(6,6) * t201;
t120 = Icges(6,1) * t169 + Icges(6,4) * t168 + Icges(6,5) * t201;
t51 = t112 * t199 + t116 * t166 + t120 * t167;
t135 = Icges(6,5) * t198 + Icges(6,6) * t197 + Icges(6,3) * t218;
t139 = Icges(6,4) * t198 + Icges(6,2) * t197 + Icges(6,6) * t218;
t143 = Icges(6,1) * t198 + Icges(6,4) * t197 + Icges(6,5) * t218;
t72 = t135 * t199 + t139 * t166 + t143 * t167;
t2 = t199 * t50 + t201 * t51 + t218 * t72;
t316 = -t2 / 0.2e1 - t1 / 0.2e1;
t52 = t109 * t201 + t113 * t168 + t117 * t169;
t53 = t110 * t201 + t114 * t168 + t118 * t169;
t73 = t134 * t201 + t138 * t168 + t142 * t169;
t3 = t199 * t52 + t201 * t53 + t218 * t73;
t54 = t111 * t201 + t115 * t168 + t119 * t169;
t55 = t112 * t201 + t116 * t168 + t120 * t169;
t74 = t135 * t201 + t139 * t168 + t143 * t169;
t4 = t199 * t54 + t201 * t55 + t218 * t74;
t315 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = t214 * t48 + t216 * t49 + t233 * t71;
t6 = t214 * t50 + t216 * t51 + t233 * t72;
t314 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t214 * t52 + t216 * t53 + t233 * t73;
t8 = t214 * t54 + t216 * t55 + t233 * t74;
t313 = t8 / 0.2e1 + t7 / 0.2e1;
t10 = t235 * t50 + t236 * t51 + t245 * t72;
t9 = t235 * t48 + t236 * t49 + t245 * t71;
t312 = t9 / 0.2e1 + t10 / 0.2e1;
t11 = t235 * t52 + t236 * t53 + t245 * t73;
t12 = t235 * t54 + t236 * t55 + t245 * t74;
t311 = t12 / 0.2e1 + t11 / 0.2e1;
t13 = t71 * t257 + (t254 * t49 - t256 * t48) * t255;
t14 = t72 * t257 + (t254 * t51 - t256 * t50) * t255;
t310 = t13 / 0.2e1 + t14 / 0.2e1;
t15 = t73 * t257 + (t254 * t53 - t256 * t52) * t255;
t16 = t74 * t257 + (t254 * t55 - t256 * t54) * t255;
t309 = t15 / 0.2e1 + t16 / 0.2e1;
t60 = t109 * t218 + t113 * t197 + t117 * t198;
t61 = t110 * t218 + t114 * t197 + t118 * t198;
t82 = t134 * t218 + t138 * t197 + t142 * t198;
t17 = t199 * t60 + t201 * t61 + t218 * t82;
t62 = t111 * t218 + t115 * t197 + t119 * t198;
t63 = t112 * t218 + t116 * t197 + t120 * t198;
t83 = t135 * t218 + t139 * t197 + t143 * t198;
t18 = t199 * t62 + t201 * t63 + t218 * t83;
t308 = t18 / 0.2e1 + t17 / 0.2e1;
t19 = t214 * t60 + t216 * t61 + t233 * t82;
t20 = t214 * t62 + t216 * t63 + t233 * t83;
t307 = t19 / 0.2e1 + t20 / 0.2e1;
t21 = t235 * t60 + t236 * t61 + t245 * t82;
t22 = t235 * t62 + t236 * t63 + t245 * t83;
t306 = t22 / 0.2e1 + t21 / 0.2e1;
t23 = t82 * t257 + (t254 * t61 - t256 * t60) * t255;
t24 = t83 * t257 + (t254 * t63 - t256 * t62) * t255;
t305 = t23 / 0.2e1 + t24 / 0.2e1;
t304 = -t163 - t333;
t303 = -t164 - t332;
t302 = -t192 + t331;
t301 = -t194 - t329;
t300 = -t209 + t327;
t299 = t257 * t164 + t320;
t206 = rSges(4,1) * t234 - rSges(4,2) * t233 + rSges(4,3) * t245;
t237 = pkin(2) * t336 + t245 * pkin(9);
t292 = (-t206 - t237) * t255;
t136 = Icges(5,5) * t200 - Icges(5,6) * t199 + Icges(5,3) * t214;
t140 = Icges(5,4) * t200 - Icges(5,2) * t199 + Icges(5,6) * t214;
t144 = Icges(5,1) * t200 - Icges(5,4) * t199 + Icges(5,5) * t214;
t78 = t136 * t214 - t140 * t199 + t144 * t200;
t137 = Icges(5,5) * t202 - Icges(5,6) * t201 + Icges(5,3) * t216;
t141 = Icges(5,4) * t202 - Icges(5,2) * t201 + Icges(5,6) * t216;
t145 = Icges(5,1) * t202 - Icges(5,4) * t201 + Icges(5,5) * t216;
t79 = t137 * t214 - t141 * t199 + t145 * t200;
t176 = Icges(5,5) * t219 - Icges(5,6) * t218 + Icges(5,3) * t233;
t177 = Icges(5,4) * t219 - Icges(5,2) * t218 + Icges(5,6) * t233;
t178 = Icges(5,1) * t219 - Icges(5,4) * t218 + Icges(5,5) * t233;
t90 = t176 * t214 - t177 * t199 + t178 * t200;
t25 = t214 * t78 + t216 * t79 + t233 * t90;
t291 = -t25 / 0.2e1 - t314;
t80 = t136 * t216 - t140 * t201 + t144 * t202;
t81 = t137 * t216 - t141 * t201 + t145 * t202;
t91 = t176 * t216 - t177 * t201 + t178 * t202;
t26 = t214 * t80 + t216 * t81 + t233 * t91;
t290 = t26 / 0.2e1 + t313;
t289 = -t192 + t304;
t288 = -t209 + t301;
t189 = t192 * t338;
t190 = t193 * t337;
t286 = 0.2e1 * t189 + 0.2e1 * t190 + t319;
t285 = t189 + t190 + t318;
t27 = t235 * t78 + t236 * t79 + t245 * t90;
t284 = t27 / 0.2e1 + t312;
t28 = t235 * t80 + t236 * t81 + t245 * t91;
t283 = t28 / 0.2e1 + t311;
t29 = t90 * t257 + (t254 * t79 - t256 * t78) * t255;
t282 = t29 / 0.2e1 + t310;
t30 = t91 * t257 + (t254 * t81 - t256 * t80) * t255;
t281 = t30 / 0.2e1 + t309;
t101 = t176 * t233 - t177 * t218 + t178 * t219;
t84 = t136 * t233 - t140 * t218 + t144 * t219;
t85 = t137 * t233 - t141 * t218 + t145 * t219;
t32 = t101 * t233 + t214 * t84 + t216 * t85;
t280 = t32 / 0.2e1 + t307;
t34 = t101 * t245 + t235 * t84 + t236 * t85;
t279 = t34 / 0.2e1 + t306;
t36 = t101 * t257 + (t254 * t85 - t256 * t84) * t255;
t278 = t36 / 0.2e1 + t305;
t277 = (-t237 + t323) * t255;
t274 = (-t237 + t300) * t255;
t160 = t163 * t338;
t161 = t164 * t337;
t271 = t160 + t161 + t285;
t270 = (-t237 + t288) * t255;
t269 = m(6) * t122 + m(7) * t333;
t267 = m(5) * t146 + t269;
t179 = rSges(4,1) * t215 - rSges(4,2) * t214 + rSges(4,3) * t235;
t265 = m(4) * t179 + t267;
t243 = t257 * rSges(3,3) + (rSges(3,1) * t262 + rSges(3,2) * t264) * t255;
t242 = Icges(3,5) * t257 + (Icges(3,1) * t262 + Icges(3,4) * t264) * t255;
t241 = Icges(3,6) * t257 + (Icges(3,4) * t262 + Icges(3,2) * t264) * t255;
t240 = Icges(3,3) * t257 + (Icges(3,5) * t262 + Icges(3,6) * t264) * t255;
t230 = rSges(3,1) * t249 + rSges(3,2) * t248 + rSges(3,3) * t338;
t229 = rSges(3,1) * t247 + rSges(3,2) * t246 - rSges(3,3) * t337;
t228 = Icges(3,1) * t249 + Icges(3,4) * t248 + Icges(3,5) * t338;
t227 = Icges(3,1) * t247 + Icges(3,4) * t246 - Icges(3,5) * t337;
t226 = Icges(3,4) * t249 + Icges(3,2) * t248 + Icges(3,6) * t338;
t225 = Icges(3,4) * t247 + Icges(3,2) * t246 - Icges(3,6) * t337;
t224 = Icges(3,5) * t249 + Icges(3,6) * t248 + Icges(3,3) * t338;
t223 = Icges(3,5) * t247 + Icges(3,6) * t246 - Icges(3,3) * t337;
t208 = -t229 * t257 - t243 * t337;
t207 = t230 * t257 - t243 * t338;
t205 = Icges(4,1) * t234 - Icges(4,4) * t233 + Icges(4,5) * t245;
t204 = Icges(4,4) * t234 - Icges(4,2) * t233 + Icges(4,6) * t245;
t203 = Icges(4,5) * t234 - Icges(4,6) * t233 + Icges(4,3) * t245;
t195 = (t229 * t254 + t230 * t256) * t255;
t175 = Icges(4,1) * t217 - Icges(4,4) * t216 + Icges(4,5) * t236;
t174 = Icges(4,1) * t215 - Icges(4,4) * t214 + Icges(4,5) * t235;
t173 = Icges(4,4) * t217 - Icges(4,2) * t216 + Icges(4,6) * t236;
t172 = Icges(4,4) * t215 - Icges(4,2) * t214 + Icges(4,6) * t235;
t171 = Icges(4,5) * t217 - Icges(4,6) * t216 + Icges(4,3) * t236;
t170 = Icges(4,5) * t215 - Icges(4,6) * t214 + Icges(4,3) * t235;
t165 = t214 * t194;
t155 = t233 * t164;
t152 = t216 * t163;
t132 = t180 * t245 - t206 * t236;
t131 = -t179 * t245 + t206 * t235;
t130 = (-t179 - t221) * t257 + t256 * t292;
t129 = t257 * t180 + t254 * t292 + t220;
t128 = t203 * t245 - t204 * t233 + t205 * t234;
t127 = t179 * t236 - t180 * t235;
t126 = t203 * t236 - t204 * t216 + t205 * t217;
t125 = t203 * t235 - t204 * t214 + t205 * t215;
t108 = (t179 * t254 + t180 * t256) * t255 + t318;
t105 = t147 * t233 - t181 * t216;
t104 = -t146 * t233 + t181 * t214;
t103 = t171 * t245 - t173 * t233 + t175 * t234;
t102 = t170 * t245 - t172 * t233 + t174 * t234;
t100 = t171 * t236 - t173 * t216 + t175 * t217;
t99 = t170 * t236 - t172 * t216 + t174 * t217;
t98 = t171 * t235 - t173 * t214 + t175 * t215;
t97 = t170 * t235 - t172 * t214 + t174 * t215;
t96 = t146 * t216 - t147 * t214;
t95 = t245 * t147 + t323 * t236 + t186;
t94 = t235 * t181 + t328 * t245 + t196;
t93 = (-t221 + t328) * t257 + t256 * t277;
t92 = t257 * t147 + t254 * t277 + t320;
t89 = t124 * t218 - t149 * t201;
t88 = -t122 * t218 + t149 * t199;
t87 = t236 * t146 + t184 + (-t147 - t193) * t235;
t86 = (t146 * t254 + t147 * t256) * t255 + t285;
t77 = t122 * t201 - t124 * t199;
t76 = t124 * t233 + t327 * t216 + t155;
t75 = t149 * t214 + t331 * t233 + t165;
t70 = (-t221 + t302) * t257 + t256 * t274;
t69 = t257 * t124 + t254 * t274 + t299;
t68 = t245 * t124 + t300 * t236 + t324;
t67 = t235 * t149 + t302 * t245 + t321;
t66 = t128 * t257 + (-t102 * t256 + t103 * t254) * t255;
t65 = t122 * t216 + t330 * t214 + t152;
t64 = t102 * t235 + t103 * t236 + t128 * t245;
t59 = (t122 * t254 + t124 * t256) * t255 + t271;
t58 = t126 * t257 + (t100 * t254 - t256 * t99) * t255;
t57 = t125 * t257 + (t254 * t98 - t256 * t97) * t255;
t56 = t236 * t122 + (-t193 + t330) * t235 + t325;
t47 = t100 * t236 + t126 * t245 + t235 * t99;
t46 = t125 * t245 + t235 * t97 + t236 * t98;
t45 = -t329 * t201 + t332 * t218;
t44 = t329 * t199 - t333 * t218;
t43 = t301 * t216 + t332 * t233 + t155;
t42 = t329 * t214 + t304 * t233 + t165;
t41 = (-t221 + t289) * t257 + t256 * t270;
t40 = t254 * t270 + t332 * t257 + t299;
t39 = t288 * t236 + t332 * t245 + t324;
t38 = t329 * t235 + t289 * t245 + t321;
t37 = -t332 * t199 + t333 * t201;
t35 = t303 * t214 + t333 * t216 + t152;
t33 = (t333 * t254 + t332 * t256) * t255 + t271;
t31 = t333 * t236 + (-t193 + t303) * t235 + t325;
t106 = [m(2) + m(3) + m(4) + m(5) + m(6) + m(7); m(4) * t319 / 0.2e1 + t286 * t353 + (m(3) * t230 - t356) * t337 + (m(3) * t229 + t265) * t338 + t355 * (0.2e1 * t160 + 0.2e1 * t161 + t286); (t33 ^ 2 + t40 ^ 2 + t41 ^ 2) * m(7) + (t59 ^ 2 + t69 ^ 2 + t70 ^ 2) * m(6) + (t86 ^ 2 + t92 ^ 2 + t93 ^ 2) * m(5) + (t108 ^ 2 + t129 ^ 2 + t130 ^ 2) * m(4) + m(3) * (t195 ^ 2 + t207 ^ 2 + t208 ^ 2) + (t15 + t16 + t30 + t58 + (t224 * t338 + t226 * t248 + t228 * t249) * t338) * t338 + (-t14 - t13 - t29 - t57 + (-t223 * t337 + t225 * t246 + t227 * t247) * t337 + (-t223 * t338 + t224 * t337 - t225 * t248 - t226 * t246 - t227 * t249 - t228 * t247) * t338) * t337 + ((t240 * t338 + t248 * t241 + t249 * t242) * t338 - (-t240 * t337 + t246 * t241 + t247 * t242) * t337 + t23 + t24 + t36 + t66 + ((t226 * t264 + t228 * t262) * t254 - (t225 * t264 + t227 * t262) * t256) * t255 ^ 2 + ((-t223 * t256 + t224 * t254 + t241 * t264 + t242 * t262) * t255 + t257 * t240) * t257) * t257; t322 * t353 + t265 * t236 + t356 * t235 + t355 * (t164 * t354 + 0.2e1 * t156 + t322); (t108 * t127 + t129 * t132 + t130 * t131) * m(4) + (t86 * t87 + t92 * t95 + t93 * t94) * m(5) + (t56 * t59 + t67 * t70 + t68 * t69) * m(6) + (t31 * t33 + t38 * t41 + t39 * t40) * m(7) + (t64 / 0.2e1 + t279) * t257 + (t66 / 0.2e1 + t278) * t245 + (t58 / 0.2e1 + t281) * t236 + (t57 / 0.2e1 + t282) * t235 + ((-t46 / 0.2e1 - t284) * t256 + (t47 / 0.2e1 + t283) * t254) * t255; (t31 ^ 2 + t38 ^ 2 + t39 ^ 2) * m(7) + (t56 ^ 2 + t67 ^ 2 + t68 ^ 2) * m(6) + (t87 ^ 2 + t94 ^ 2 + t95 ^ 2) * m(5) + (t127 ^ 2 + t131 ^ 2 + t132 ^ 2) * m(4) + (t22 + t21 + t34 + t64) * t245 + (t11 + t12 + t28 + t47) * t236 + (t9 + t10 + t27 + t46) * t235; t266 * t214 + t267 * t216 + 0.2e1 * t355 * (-t214 * t164 + t152); (t33 * t35 + t40 * t43 + t41 * t42) * m(7) + (t59 * t65 + t69 * t76 + t70 * t75) * m(6) + (t104 * t93 + t105 * t92 + t86 * t96) * m(5) + t280 * t257 + t278 * t233 + t281 * t216 + t282 * t214 + (t290 * t254 + t291 * t256) * t255; (t31 * t35 + t38 * t42 + t39 * t43) * m(7) + (t56 * t65 + t67 * t75 + t68 * t76) * m(6) + (t104 * t94 + t105 * t95 + t87 * t96) * m(5) + t280 * t245 + t290 * t236 - t291 * t235 + t279 * t233 + t283 * t216 + t284 * t214; (t35 ^ 2 + t42 ^ 2 + t43 ^ 2) * m(7) + (t65 ^ 2 + t75 ^ 2 + t76 ^ 2) * m(6) + (t104 ^ 2 + t105 ^ 2 + t96 ^ 2) * m(5) + (t19 + t20 + t32) * t233 + (t7 + t8 + t26) * t216 + (t6 + t5 + t25) * t214; t268 * t199 + t269 * t201; (t33 * t37 + t40 * t45 + t41 * t44) * m(7) + (t59 * t77 + t69 * t89 + t70 * t88) * m(6) + t308 * t257 + t305 * t218 + t309 * t201 + t310 * t199 + (t315 * t254 + t316 * t256) * t255; (t31 * t37 + t38 * t44 + t39 * t45) * m(7) + (t56 * t77 + t67 * t88 + t68 * t89) * m(6) + t308 * t245 + t315 * t236 - t316 * t235 + t306 * t218 + t311 * t201 + t312 * t199; (t35 * t37 + t42 * t44 + t43 * t45) * m(7) + (t65 * t77 + t75 * t88 + t76 * t89) * m(6) + t308 * t233 + t307 * t218 + t315 * t216 - t316 * t214 + t313 * t201 + t314 * t199; (t37 ^ 2 + t44 ^ 2 + t45 ^ 2) * m(7) + (t77 ^ 2 + t88 ^ 2 + t89 ^ 2) * m(6) + (t17 + t18) * t218 + (t3 + t4) * t201 + (t2 + t1) * t199; t218 * m(7); (t199 * t40 + t201 * t41 + t218 * t33) * m(7); (t199 * t39 + t201 * t38 + t218 * t31) * m(7); (t199 * t43 + t201 * t42 + t218 * t35) * m(7); (t199 * t45 + t201 * t44 + t218 * t37) * m(7); (t199 ^ 2 + t201 ^ 2 + t218 ^ 2) * m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t106(1) t106(2) t106(4) t106(7) t106(11) t106(16); t106(2) t106(3) t106(5) t106(8) t106(12) t106(17); t106(4) t106(5) t106(6) t106(9) t106(13) t106(18); t106(7) t106(8) t106(9) t106(10) t106(14) t106(19); t106(11) t106(12) t106(13) t106(14) t106(15) t106(20); t106(16) t106(17) t106(18) t106(19) t106(20) t106(21);];
Mq  = res;
