% Calculate joint inertia matrix for
% S6PRRRRP6
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
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:28:42
% EndTime: 2019-03-09 00:28:58
% DurationCPUTime: 7.92s
% Computational Cost: add. (59183->659), mult. (164330->930), div. (0->0), fcn. (217210->14), ass. (0->308)
t355 = rSges(7,1) + pkin(5);
t354 = rSges(7,3) + qJ(6);
t256 = sin(pkin(12));
t259 = cos(pkin(6));
t263 = sin(qJ(2));
t258 = cos(pkin(12));
t264 = cos(qJ(2));
t334 = t264 * t258;
t249 = -t256 * t263 + t259 * t334;
t257 = sin(pkin(6));
t341 = cos(pkin(7));
t298 = t257 * t341;
t340 = sin(pkin(7));
t238 = -t249 * t340 - t258 * t298;
t335 = t264 * t256;
t336 = t259 * t263;
t250 = t258 * t336 + t335;
t223 = t250 * pkin(2) + t238 * pkin(9);
t251 = -t258 * t263 - t259 * t335;
t239 = -t251 * t340 + t256 * t298;
t252 = -t256 * t336 + t334;
t224 = t252 * pkin(2) + t239 * pkin(9);
t338 = t257 * t258;
t339 = t256 * t257;
t318 = t223 * t339 + t224 * t338;
t262 = sin(qJ(3));
t347 = cos(qJ(3));
t275 = t347 * t340;
t273 = t257 * t275;
t276 = t341 * t347;
t218 = -t251 * t276 + t252 * t262 - t256 * t273;
t297 = t257 * t340;
t219 = t252 * t347 + (t341 * t251 + t256 * t297) * t262;
t180 = rSges(4,1) * t219 - rSges(4,2) * t218 + rSges(4,3) * t239;
t261 = sin(qJ(4));
t346 = cos(qJ(4));
t201 = t219 * t261 - t239 * t346;
t202 = t219 * t346 + t239 * t261;
t146 = rSges(5,1) * t202 - rSges(5,2) * t201 + rSges(5,3) * t218;
t260 = sin(qJ(5));
t345 = cos(qJ(5));
t168 = t202 * t260 - t218 * t345;
t169 = t202 * t345 + t218 * t260;
t122 = rSges(6,1) * t169 - rSges(6,2) * t168 + rSges(6,3) * t201;
t331 = rSges(7,2) * t201 + t354 * t168 + t169 * t355;
t268 = -m(6) * t122 - m(7) * t331;
t266 = -m(5) * t146 + t268;
t353 = -m(4) * t180 + t266;
t352 = m(6) / 0.2e1 + m(7) / 0.2e1;
t351 = -0.2e1 * t238;
t350 = m(5) / 0.2e1;
t337 = t257 * t263;
t217 = t250 * t347 + (t341 * t249 - t258 * t297) * t262;
t200 = t217 * t346 + t238 * t261;
t216 = -t249 * t276 + t250 * t262 + t258 * t273;
t166 = t200 * t260 - t216 * t345;
t167 = t200 * t345 + t216 * t260;
t199 = t217 * t261 - t238 * t346;
t333 = rSges(7,2) * t199 + t354 * t166 + t167 * t355;
t120 = rSges(6,1) * t167 - rSges(6,2) * t166 + rSges(6,3) * t199;
t163 = pkin(4) * t200 + pkin(11) * t199;
t332 = -t120 - t163;
t164 = pkin(4) * t202 + pkin(11) * t201;
t330 = -t122 - t164;
t145 = rSges(5,1) * t200 - rSges(5,2) * t199 + rSges(5,3) * t216;
t192 = pkin(3) * t217 + pkin(10) * t216;
t329 = -t145 - t192;
t237 = t259 * t340 * t262 + (t341 * t262 * t264 + t347 * t263) * t257;
t248 = t259 * t341 - t264 * t297;
t221 = t237 * t346 + t248 * t261;
t236 = -t257 * t264 * t276 - t259 * t275 + t262 * t337;
t197 = t221 * t260 - t236 * t345;
t198 = t221 * t345 + t236 * t260;
t220 = t237 * t261 - t248 * t346;
t328 = rSges(7,2) * t220 + t354 * t197 + t198 * t355;
t148 = rSges(6,1) * t198 - rSges(6,2) * t197 + rSges(6,3) * t220;
t194 = pkin(4) * t221 + pkin(11) * t220;
t327 = -t148 - t194;
t155 = t239 * t163;
t184 = t239 * t192;
t325 = t155 + t184;
t193 = pkin(3) * t219 + pkin(10) * t218;
t186 = t248 * t193;
t324 = t248 * t164 + t186;
t181 = rSges(5,1) * t221 - rSges(5,2) * t220 + rSges(5,3) * t236;
t209 = pkin(3) * t237 + pkin(10) * t236;
t323 = -t181 - t209;
t322 = t193 * t351 + 0.2e1 * t184;
t196 = t238 * t209;
t321 = t238 * t194 + t196;
t222 = t259 * t224;
t320 = t259 * t193 + t222;
t319 = 0.2e1 * t318;
t107 = Icges(7,5) * t167 + Icges(7,6) * t199 + Icges(7,3) * t166;
t111 = Icges(7,4) * t167 + Icges(7,2) * t199 + Icges(7,6) * t166;
t115 = Icges(7,1) * t167 + Icges(7,4) * t199 + Icges(7,5) * t166;
t46 = t107 * t166 + t111 * t199 + t115 * t167;
t108 = Icges(7,5) * t169 + Icges(7,6) * t201 + Icges(7,3) * t168;
t112 = Icges(7,4) * t169 + Icges(7,2) * t201 + Icges(7,6) * t168;
t116 = Icges(7,1) * t169 + Icges(7,4) * t201 + Icges(7,5) * t168;
t47 = t108 * t166 + t112 * t199 + t116 * t167;
t133 = Icges(7,5) * t198 + Icges(7,6) * t220 + Icges(7,3) * t197;
t137 = Icges(7,4) * t198 + Icges(7,2) * t220 + Icges(7,6) * t197;
t141 = Icges(7,1) * t198 + Icges(7,4) * t220 + Icges(7,5) * t197;
t71 = t133 * t166 + t137 * t199 + t141 * t167;
t1 = t199 * t46 + t201 * t47 + t220 * t71;
t109 = Icges(6,5) * t167 - Icges(6,6) * t166 + Icges(6,3) * t199;
t113 = Icges(6,4) * t167 - Icges(6,2) * t166 + Icges(6,6) * t199;
t117 = Icges(6,1) * t167 - Icges(6,4) * t166 + Icges(6,5) * t199;
t48 = t109 * t199 - t113 * t166 + t117 * t167;
t110 = Icges(6,5) * t169 - Icges(6,6) * t168 + Icges(6,3) * t201;
t114 = Icges(6,4) * t169 - Icges(6,2) * t168 + Icges(6,6) * t201;
t118 = Icges(6,1) * t169 - Icges(6,4) * t168 + Icges(6,5) * t201;
t49 = t110 * t199 - t114 * t166 + t118 * t167;
t134 = Icges(6,5) * t198 - Icges(6,6) * t197 + Icges(6,3) * t220;
t138 = Icges(6,4) * t198 - Icges(6,2) * t197 + Icges(6,6) * t220;
t142 = Icges(6,1) * t198 - Icges(6,4) * t197 + Icges(6,5) * t220;
t72 = t134 * t199 - t138 * t166 + t142 * t167;
t2 = t199 * t48 + t201 * t49 + t220 * t72;
t316 = t2 / 0.2e1 + t1 / 0.2e1;
t50 = t107 * t168 + t111 * t201 + t115 * t169;
t51 = t108 * t168 + t112 * t201 + t116 * t169;
t73 = t133 * t168 + t137 * t201 + t141 * t169;
t3 = t199 * t50 + t201 * t51 + t220 * t73;
t52 = t109 * t201 - t113 * t168 + t117 * t169;
t53 = t110 * t201 - t114 * t168 + t118 * t169;
t74 = t134 * t201 - t138 * t168 + t142 * t169;
t4 = t199 * t52 + t201 * t53 + t220 * t74;
t315 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = t216 * t46 + t218 * t47 + t236 * t71;
t6 = t216 * t48 + t218 * t49 + t236 * t72;
t314 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t216 * t50 + t218 * t51 + t236 * t73;
t8 = t216 * t52 + t218 * t53 + t236 * t74;
t313 = t7 / 0.2e1 + t8 / 0.2e1;
t10 = t238 * t48 + t239 * t49 + t248 * t72;
t9 = t238 * t46 + t239 * t47 + t248 * t71;
t312 = t10 / 0.2e1 + t9 / 0.2e1;
t11 = t238 * t50 + t239 * t51 + t248 * t73;
t12 = t238 * t52 + t239 * t53 + t248 * t74;
t311 = t12 / 0.2e1 + t11 / 0.2e1;
t13 = t259 * t71 + (t256 * t47 - t258 * t46) * t257;
t14 = t259 * t72 + (t256 * t49 - t258 * t48) * t257;
t310 = t14 / 0.2e1 + t13 / 0.2e1;
t15 = t259 * t73 + (t256 * t51 - t258 * t50) * t257;
t16 = t259 * t74 + (t256 * t53 - t258 * t52) * t257;
t309 = t15 / 0.2e1 + t16 / 0.2e1;
t58 = t107 * t197 + t111 * t220 + t115 * t198;
t59 = t108 * t197 + t112 * t220 + t116 * t198;
t82 = t133 * t197 + t137 * t220 + t141 * t198;
t17 = t199 * t58 + t201 * t59 + t220 * t82;
t60 = t109 * t220 - t113 * t197 + t117 * t198;
t61 = t110 * t220 - t114 * t197 + t118 * t198;
t83 = t134 * t220 - t138 * t197 + t142 * t198;
t18 = t199 * t60 + t201 * t61 + t220 * t83;
t308 = t17 / 0.2e1 + t18 / 0.2e1;
t19 = t216 * t58 + t218 * t59 + t236 * t82;
t20 = t216 * t60 + t218 * t61 + t236 * t83;
t307 = t19 / 0.2e1 + t20 / 0.2e1;
t21 = t238 * t58 + t239 * t59 + t248 * t82;
t22 = t238 * t60 + t239 * t61 + t248 * t83;
t306 = t22 / 0.2e1 + t21 / 0.2e1;
t23 = t259 * t82 + (t256 * t59 - t258 * t58) * t257;
t24 = t259 * t83 + (t256 * t61 - t258 * t60) * t257;
t305 = t23 / 0.2e1 + t24 / 0.2e1;
t304 = -t163 - t333;
t303 = -t192 + t332;
t302 = -t164 - t331;
t301 = -t194 - t328;
t300 = -t209 + t327;
t299 = t259 * t164 + t320;
t206 = rSges(4,1) * t237 - rSges(4,2) * t236 + rSges(4,3) * t248;
t240 = pkin(2) * t337 + t248 * pkin(9);
t292 = (-t206 - t240) * t257;
t135 = Icges(5,5) * t200 - Icges(5,6) * t199 + Icges(5,3) * t216;
t139 = Icges(5,4) * t200 - Icges(5,2) * t199 + Icges(5,6) * t216;
t143 = Icges(5,1) * t200 - Icges(5,4) * t199 + Icges(5,5) * t216;
t78 = t135 * t216 - t139 * t199 + t143 * t200;
t136 = Icges(5,5) * t202 - Icges(5,6) * t201 + Icges(5,3) * t218;
t140 = Icges(5,4) * t202 - Icges(5,2) * t201 + Icges(5,6) * t218;
t144 = Icges(5,1) * t202 - Icges(5,4) * t201 + Icges(5,5) * t218;
t79 = t136 * t216 - t140 * t199 + t144 * t200;
t176 = Icges(5,5) * t221 - Icges(5,6) * t220 + Icges(5,3) * t236;
t177 = Icges(5,4) * t221 - Icges(5,2) * t220 + Icges(5,6) * t236;
t178 = Icges(5,1) * t221 - Icges(5,4) * t220 + Icges(5,5) * t236;
t90 = t176 * t216 - t177 * t199 + t178 * t200;
t25 = t216 * t78 + t218 * t79 + t236 * t90;
t291 = t25 / 0.2e1 + t314;
t80 = t135 * t218 - t139 * t201 + t143 * t202;
t81 = t136 * t218 - t140 * t201 + t144 * t202;
t91 = t176 * t218 - t177 * t201 + t178 * t202;
t26 = t216 * t80 + t218 * t81 + t236 * t91;
t290 = t26 / 0.2e1 + t313;
t289 = -t192 + t304;
t288 = -t209 + t301;
t189 = t192 * t339;
t190 = t193 * t338;
t286 = 0.2e1 * t189 + 0.2e1 * t190 + t319;
t285 = t189 + t190 + t318;
t27 = t238 * t78 + t239 * t79 + t248 * t90;
t284 = t27 / 0.2e1 + t312;
t28 = t238 * t80 + t239 * t81 + t248 * t91;
t283 = t28 / 0.2e1 + t311;
t29 = t259 * t90 + (t256 * t79 - t258 * t78) * t257;
t282 = t29 / 0.2e1 + t310;
t30 = t259 * t91 + (t256 * t81 - t258 * t80) * t257;
t281 = t30 / 0.2e1 + t309;
t101 = t176 * t236 - t177 * t220 + t178 * t221;
t84 = t135 * t236 - t139 * t220 + t143 * t221;
t85 = t136 * t236 - t140 * t220 + t144 * t221;
t31 = t101 * t236 + t216 * t84 + t218 * t85;
t280 = t31 / 0.2e1 + t307;
t32 = t101 * t248 + t238 * t84 + t239 * t85;
t279 = t32 / 0.2e1 + t306;
t33 = t101 * t259 + (t256 * t85 - t258 * t84) * t257;
t278 = t33 / 0.2e1 + t305;
t277 = (-t240 + t323) * t257;
t274 = (-t240 + t300) * t257;
t159 = t163 * t339;
t160 = t164 * t338;
t271 = t159 + t160 + t285;
t270 = (-t240 + t288) * t257;
t269 = m(6) * t120 + m(7) * t333;
t267 = m(5) * t145 + t269;
t179 = rSges(4,1) * t217 - rSges(4,2) * t216 + rSges(4,3) * t238;
t265 = m(4) * t179 + t267;
t246 = t259 * rSges(3,3) + (rSges(3,1) * t263 + rSges(3,2) * t264) * t257;
t245 = Icges(3,5) * t259 + (Icges(3,1) * t263 + Icges(3,4) * t264) * t257;
t244 = Icges(3,6) * t259 + (Icges(3,4) * t263 + Icges(3,2) * t264) * t257;
t243 = Icges(3,3) * t259 + (Icges(3,5) * t263 + Icges(3,6) * t264) * t257;
t232 = rSges(3,1) * t252 + rSges(3,2) * t251 + rSges(3,3) * t339;
t231 = rSges(3,1) * t250 + rSges(3,2) * t249 - rSges(3,3) * t338;
t230 = Icges(3,1) * t252 + Icges(3,4) * t251 + Icges(3,5) * t339;
t229 = Icges(3,1) * t250 + Icges(3,4) * t249 - Icges(3,5) * t338;
t228 = Icges(3,4) * t252 + Icges(3,2) * t251 + Icges(3,6) * t339;
t227 = Icges(3,4) * t250 + Icges(3,2) * t249 - Icges(3,6) * t338;
t226 = Icges(3,5) * t252 + Icges(3,6) * t251 + Icges(3,3) * t339;
t225 = Icges(3,5) * t250 + Icges(3,6) * t249 - Icges(3,3) * t338;
t208 = -t231 * t259 - t246 * t338;
t207 = t232 * t259 - t246 * t339;
t205 = Icges(4,1) * t237 - Icges(4,4) * t236 + Icges(4,5) * t248;
t204 = Icges(4,4) * t237 - Icges(4,2) * t236 + Icges(4,6) * t248;
t203 = Icges(4,5) * t237 - Icges(4,6) * t236 + Icges(4,3) * t248;
t195 = (t231 * t256 + t232 * t258) * t257;
t175 = Icges(4,1) * t219 - Icges(4,4) * t218 + Icges(4,5) * t239;
t174 = Icges(4,1) * t217 - Icges(4,4) * t216 + Icges(4,5) * t238;
t173 = Icges(4,4) * t219 - Icges(4,2) * t218 + Icges(4,6) * t239;
t172 = Icges(4,4) * t217 - Icges(4,2) * t216 + Icges(4,6) * t238;
t171 = Icges(4,5) * t219 - Icges(4,6) * t218 + Icges(4,3) * t239;
t170 = Icges(4,5) * t217 - Icges(4,6) * t216 + Icges(4,3) * t238;
t165 = t216 * t194;
t154 = t236 * t164;
t151 = t218 * t163;
t132 = t180 * t248 - t206 * t239;
t131 = -t179 * t248 + t206 * t238;
t128 = (-t179 - t223) * t259 + t258 * t292;
t127 = t180 * t259 + t256 * t292 + t222;
t126 = t203 * t248 - t204 * t236 + t205 * t237;
t125 = t179 * t239 - t180 * t238;
t124 = t203 * t239 - t204 * t218 + t205 * t219;
t123 = t203 * t238 - t204 * t216 + t205 * t217;
t106 = (t179 * t256 + t180 * t258) * t257 + t318;
t105 = t146 * t236 - t181 * t218;
t104 = -t145 * t236 + t181 * t216;
t103 = t171 * t248 - t173 * t236 + t175 * t237;
t102 = t170 * t248 - t172 * t236 + t174 * t237;
t100 = t171 * t239 - t173 * t218 + t175 * t219;
t99 = t170 * t239 - t172 * t218 + t174 * t219;
t98 = t171 * t238 - t173 * t216 + t175 * t217;
t97 = t170 * t238 - t172 * t216 + t174 * t217;
t96 = t145 * t218 - t146 * t216;
t95 = t146 * t248 + t323 * t239 + t186;
t94 = t181 * t238 + t329 * t248 + t196;
t93 = (-t223 + t329) * t259 + t258 * t277;
t92 = t146 * t259 + t256 * t277 + t320;
t89 = t122 * t220 - t148 * t201;
t88 = -t120 * t220 + t148 * t199;
t87 = t145 * t239 + t184 + (-t146 - t193) * t238;
t86 = (t145 * t256 + t146 * t258) * t257 + t285;
t77 = t120 * t201 - t122 * t199;
t76 = t122 * t236 + t327 * t218 + t154;
t75 = t148 * t216 + t332 * t236 + t165;
t70 = (-t223 + t303) * t259 + t258 * t274;
t69 = t122 * t259 + t256 * t274 + t299;
t68 = t122 * t248 + t300 * t239 + t324;
t67 = t148 * t238 + t303 * t248 + t321;
t66 = t126 * t259 + (-t102 * t258 + t103 * t256) * t257;
t65 = t120 * t218 + t330 * t216 + t151;
t64 = -t328 * t201 + t331 * t220;
t63 = t328 * t199 - t333 * t220;
t62 = t102 * t238 + t103 * t239 + t126 * t248;
t57 = (t120 * t256 + t122 * t258) * t257 + t271;
t56 = t124 * t259 + (t100 * t256 - t258 * t99) * t257;
t55 = t123 * t259 + (t256 * t98 - t258 * t97) * t257;
t54 = t120 * t239 + (-t193 + t330) * t238 + t325;
t45 = t100 * t239 + t124 * t248 + t238 * t99;
t44 = t123 * t248 + t238 * t97 + t239 * t98;
t43 = t301 * t218 + t331 * t236 + t154;
t42 = t328 * t216 + t304 * t236 + t165;
t41 = (-t223 + t289) * t259 + t258 * t270;
t40 = t256 * t270 + t331 * t259 + t299;
t39 = -t331 * t199 + t333 * t201;
t38 = t288 * t239 + t331 * t248 + t324;
t37 = t328 * t238 + t289 * t248 + t321;
t36 = t302 * t216 + t333 * t218 + t151;
t35 = (t333 * t256 + t331 * t258) * t257 + t271;
t34 = t333 * t239 + (-t193 + t302) * t238 + t325;
t119 = [m(2) + m(3) + m(4) + m(5) + m(6) + m(7); m(4) * t319 / 0.2e1 + t286 * t350 + (m(3) * t232 - t353) * t338 + (m(3) * t231 + t265) * t339 + t352 * (0.2e1 * t159 + 0.2e1 * t160 + t286); (t35 ^ 2 + t40 ^ 2 + t41 ^ 2) * m(7) + (t57 ^ 2 + t69 ^ 2 + t70 ^ 2) * m(6) + (t86 ^ 2 + t92 ^ 2 + t93 ^ 2) * m(5) + (t106 ^ 2 + t127 ^ 2 + t128 ^ 2) * m(4) + m(3) * (t195 ^ 2 + t207 ^ 2 + t208 ^ 2) + (t15 + t16 + t30 + t56 + (t226 * t339 + t228 * t251 + t230 * t252) * t339) * t339 + (-t13 - t14 - t29 - t55 + (-t225 * t338 + t227 * t249 + t229 * t250) * t338 + (-t225 * t339 + t226 * t338 - t227 * t251 - t228 * t249 - t229 * t252 - t230 * t250) * t339) * t338 + ((t243 * t339 + t244 * t251 + t245 * t252) * t339 - (-t243 * t338 + t244 * t249 + t245 * t250) * t338 + t23 + t24 + t33 + t66 + ((t228 * t264 + t230 * t263) * t256 - (t227 * t264 + t229 * t263) * t258) * t257 ^ 2 + ((-t225 * t258 + t226 * t256 + t244 * t264 + t245 * t263) * t257 + t259 * t243) * t259) * t259; t322 * t350 + t265 * t239 + t353 * t238 + t352 * (t164 * t351 + 0.2e1 * t155 + t322); (t34 * t35 + t37 * t41 + t38 * t40) * m(7) + (t54 * t57 + t67 * t70 + t68 * t69) * m(6) + (t86 * t87 + t92 * t95 + t93 * t94) * m(5) + (t106 * t125 + t127 * t132 + t128 * t131) * m(4) + (t62 / 0.2e1 + t279) * t259 + (t66 / 0.2e1 + t278) * t248 + (t56 / 0.2e1 + t281) * t239 + (t55 / 0.2e1 + t282) * t238 + ((-t44 / 0.2e1 - t284) * t258 + (t45 / 0.2e1 + t283) * t256) * t257; (t34 ^ 2 + t37 ^ 2 + t38 ^ 2) * m(7) + (t54 ^ 2 + t67 ^ 2 + t68 ^ 2) * m(6) + (t87 ^ 2 + t94 ^ 2 + t95 ^ 2) * m(5) + (t125 ^ 2 + t131 ^ 2 + t132 ^ 2) * m(4) + (t21 + t22 + t32 + t62) * t248 + (t12 + t11 + t28 + t45) * t239 + (t10 + t9 + t27 + t44) * t238; t266 * t216 + t267 * t218 + 0.2e1 * t352 * (-t216 * t164 + t151); (t35 * t36 + t40 * t43 + t41 * t42) * m(7) + (t57 * t65 + t69 * t76 + t70 * t75) * m(6) + (t104 * t93 + t105 * t92 + t86 * t96) * m(5) + t280 * t259 + t278 * t236 + t281 * t218 + t282 * t216 + (t290 * t256 - t291 * t258) * t257; (t34 * t36 + t37 * t42 + t38 * t43) * m(7) + (t54 * t65 + t67 * t75 + t68 * t76) * m(6) + (t104 * t94 + t105 * t95 + t87 * t96) * m(5) + t280 * t248 + t290 * t239 + t291 * t238 + t279 * t236 + t283 * t218 + t284 * t216; (t36 ^ 2 + t42 ^ 2 + t43 ^ 2) * m(7) + (t65 ^ 2 + t75 ^ 2 + t76 ^ 2) * m(6) + (t104 ^ 2 + t105 ^ 2 + t96 ^ 2) * m(5) + (t20 + t19 + t31) * t236 + (t7 + t8 + t26) * t218 + (t6 + t5 + t25) * t216; t268 * t199 + t269 * t201; (t35 * t39 + t40 * t64 + t41 * t63) * m(7) + (t57 * t77 + t69 * t89 + t70 * t88) * m(6) + t308 * t259 + t305 * t220 + t309 * t201 + t310 * t199 + (t315 * t256 - t316 * t258) * t257; (t34 * t39 + t37 * t63 + t38 * t64) * m(7) + (t54 * t77 + t67 * t88 + t68 * t89) * m(6) + t308 * t248 + t315 * t239 + t316 * t238 + t306 * t220 + t311 * t201 + t312 * t199; (t36 * t39 + t42 * t63 + t43 * t64) * m(7) + (t65 * t77 + t75 * t88 + t76 * t89) * m(6) + t308 * t236 + t307 * t220 + t315 * t218 + t316 * t216 + t313 * t201 + t314 * t199; (t39 ^ 2 + t63 ^ 2 + t64 ^ 2) * m(7) + (t77 ^ 2 + t88 ^ 2 + t89 ^ 2) * m(6) + (t18 + t17) * t220 + (t4 + t3) * t201 + (t1 + t2) * t199; t197 * m(7); (t166 * t40 + t168 * t41 + t197 * t35) * m(7); (t166 * t38 + t168 * t37 + t197 * t34) * m(7); (t166 * t43 + t168 * t42 + t197 * t36) * m(7); (t166 * t64 + t168 * t63 + t197 * t39) * m(7); (t166 ^ 2 + t168 ^ 2 + t197 ^ 2) * m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t119(1) t119(2) t119(4) t119(7) t119(11) t119(16); t119(2) t119(3) t119(5) t119(8) t119(12) t119(17); t119(4) t119(5) t119(6) t119(9) t119(13) t119(18); t119(7) t119(8) t119(9) t119(10) t119(14) t119(19); t119(11) t119(12) t119(13) t119(14) t119(15) t119(20); t119(16) t119(17) t119(18) t119(19) t119(20) t119(21);];
Mq  = res;
