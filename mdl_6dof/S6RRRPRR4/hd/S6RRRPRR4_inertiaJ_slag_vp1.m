% Calculate joint inertia matrix for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:15:25
% EndTime: 2019-03-09 18:15:35
% DurationCPUTime: 4.42s
% Computational Cost: add. (15457->519), mult. (12837->735), div. (0->0), fcn. (13769->12), ass. (0->249)
t242 = pkin(11) + qJ(5);
t234 = qJ(6) + t242;
t227 = sin(t234);
t228 = cos(t234);
t250 = sin(qJ(1));
t245 = qJ(2) + qJ(3);
t236 = cos(t245);
t252 = cos(qJ(1));
t310 = t236 * t252;
t169 = -t227 * t310 + t228 * t250;
t170 = t227 * t250 + t228 * t310;
t235 = sin(t245);
t312 = t235 * t252;
t111 = t170 * rSges(7,1) + t169 * rSges(7,2) + rSges(7,3) * t312;
t247 = cos(pkin(11));
t229 = t247 * pkin(4) + pkin(3);
t233 = cos(t242);
t207 = pkin(5) * t233 + t229;
t232 = sin(t242);
t246 = sin(pkin(11));
t208 = pkin(4) * t246 + pkin(5) * t232;
t342 = t207 * t310 + t250 * t208 + t111;
t248 = -pkin(9) - qJ(4);
t241 = -pkin(10) + t248;
t293 = t241 - t248;
t309 = t246 * t250;
t298 = -pkin(4) * t309 - t229 * t310;
t341 = -t293 * t312 + t298 + t342;
t306 = t247 * t252;
t195 = -t236 * t309 - t306;
t307 = t247 * t250;
t308 = t246 * t252;
t196 = t236 * t307 - t308;
t313 = t235 * t250;
t129 = Icges(5,5) * t196 + Icges(5,6) * t195 + Icges(5,3) * t313;
t197 = -t236 * t308 + t307;
t198 = t236 * t306 + t309;
t130 = Icges(5,5) * t198 + Icges(5,6) * t197 + Icges(5,3) * t312;
t131 = Icges(5,4) * t196 + Icges(5,2) * t195 + Icges(5,6) * t313;
t132 = Icges(5,4) * t198 + Icges(5,2) * t197 + Icges(5,6) * t312;
t133 = Icges(5,1) * t196 + Icges(5,4) * t195 + Icges(5,5) * t313;
t134 = Icges(5,1) * t198 + Icges(5,4) * t197 + Icges(5,5) * t312;
t311 = t236 * t250;
t166 = -t227 * t311 - t228 * t252;
t167 = -t227 * t252 + t228 * t311;
t104 = Icges(7,5) * t167 + Icges(7,6) * t166 + Icges(7,3) * t313;
t106 = Icges(7,4) * t167 + Icges(7,2) * t166 + Icges(7,6) * t313;
t108 = Icges(7,1) * t167 + Icges(7,4) * t166 + Icges(7,5) * t313;
t33 = t104 * t313 + t106 * t166 + t108 * t167;
t105 = Icges(7,5) * t170 + Icges(7,6) * t169 + Icges(7,3) * t312;
t107 = Icges(7,4) * t170 + Icges(7,2) * t169 + Icges(7,6) * t312;
t109 = Icges(7,1) * t170 + Icges(7,4) * t169 + Icges(7,5) * t312;
t34 = t105 * t313 + t107 * t166 + t109 * t167;
t16 = t250 * t34 - t252 * t33;
t266 = Icges(4,5) * t236 - Icges(4,6) * t235;
t171 = -Icges(4,3) * t252 + t250 * t266;
t172 = Icges(4,3) * t250 + t252 * t266;
t182 = -t232 * t311 - t233 * t252;
t183 = -t232 * t252 + t233 * t311;
t114 = Icges(6,5) * t183 + Icges(6,6) * t182 + Icges(6,3) * t313;
t116 = Icges(6,4) * t183 + Icges(6,2) * t182 + Icges(6,6) * t313;
t118 = Icges(6,1) * t183 + Icges(6,4) * t182 + Icges(6,5) * t313;
t37 = t114 * t313 + t116 * t182 + t118 * t183;
t184 = -t232 * t310 + t233 * t250;
t185 = t232 * t250 + t233 * t310;
t115 = Icges(6,5) * t185 + Icges(6,6) * t184 + Icges(6,3) * t312;
t117 = Icges(6,4) * t185 + Icges(6,2) * t184 + Icges(6,6) * t312;
t119 = Icges(6,1) * t185 + Icges(6,4) * t184 + Icges(6,5) * t312;
t38 = t115 * t313 + t117 * t182 + t119 * t183;
t21 = t250 * t38 - t252 * t37;
t244 = t252 ^ 2;
t317 = Icges(4,4) * t236;
t268 = -Icges(4,2) * t235 + t317;
t174 = Icges(4,6) * t250 + t252 * t268;
t318 = Icges(4,4) * t235;
t270 = Icges(4,1) * t236 - t318;
t176 = Icges(4,5) * t250 + t252 * t270;
t264 = -t174 * t235 + t176 * t236;
t173 = -Icges(4,6) * t252 + t250 * t268;
t175 = -Icges(4,5) * t252 + t250 * t270;
t265 = t173 * t235 - t175 * t236;
t340 = -t16 - t21 + (t129 * t313 + t131 * t195 + t133 * t196) * t252 - t244 * t171 + (-t130 * t313 - t132 * t195 - t134 * t196 - (-t172 + t265) * t252 - t264 * t250) * t250;
t243 = t250 ^ 2;
t339 = m(5) / 0.2e1;
t338 = m(6) / 0.2e1;
t337 = m(7) / 0.2e1;
t145 = -Icges(7,3) * t236 + (Icges(7,5) * t228 - Icges(7,6) * t227) * t235;
t146 = -Icges(7,6) * t236 + (Icges(7,4) * t228 - Icges(7,2) * t227) * t235;
t147 = -Icges(7,5) * t236 + (Icges(7,1) * t228 - Icges(7,4) * t227) * t235;
t60 = t145 * t313 + t146 * t166 + t147 * t167;
t5 = -t60 * t236 + (t250 * t33 + t252 * t34) * t235;
t35 = t104 * t312 + t106 * t169 + t108 * t170;
t36 = t105 * t312 + t107 * t169 + t109 * t170;
t61 = t145 * t312 + t146 * t169 + t147 * t170;
t6 = -t61 * t236 + (t250 * t35 + t252 * t36) * t235;
t336 = t6 * t312 + t5 * t313;
t335 = -t236 / 0.2e1;
t334 = t250 / 0.2e1;
t333 = -t252 / 0.2e1;
t249 = sin(qJ(2));
t332 = pkin(2) * t249;
t331 = pkin(3) * t236;
t330 = -pkin(3) + t229;
t251 = cos(qJ(2));
t329 = rSges(3,1) * t251;
t328 = rSges(3,2) * t249;
t327 = t252 * rSges(3,3);
t45 = -t236 * t104 + (-t106 * t227 + t108 * t228) * t235;
t326 = t45 * t252;
t46 = -t236 * t105 + (-t107 * t227 + t109 * t228) * t235;
t325 = t46 * t250;
t50 = -t236 * t114 + (-t116 * t232 + t118 * t233) * t235;
t324 = t50 * t252;
t51 = -t236 * t115 + (-t117 * t232 + t119 * t233) * t235;
t323 = t51 * t250;
t139 = t235 * t228 * t147;
t315 = t146 * t227;
t75 = -t236 * t145 - t235 * t315 + t139;
t322 = t75 * t236;
t272 = -t167 * rSges(7,1) - t166 * rSges(7,2);
t110 = rSges(7,3) * t313 - t272;
t148 = -t236 * rSges(7,3) + (rSges(7,1) * t228 - rSges(7,2) * t227) * t235;
t80 = t236 * t110 + t148 * t313;
t320 = Icges(3,4) * t249;
t319 = Icges(3,4) * t251;
t316 = qJ(4) * t235;
t150 = -Icges(6,6) * t236 + (Icges(6,4) * t233 - Icges(6,2) * t232) * t235;
t314 = t150 * t232;
t253 = -pkin(8) - pkin(7);
t305 = t252 * t253;
t297 = t207 - t229;
t135 = t235 * t297 + t236 * t293;
t304 = -t135 - t148;
t144 = (qJ(4) + t248) * t236 + t330 * t235;
t205 = t235 * pkin(3) - t236 * qJ(4);
t303 = -t144 - t205;
t230 = pkin(2) * t251 + pkin(1);
t218 = t252 * t230;
t240 = t252 * pkin(7);
t302 = t250 * (t305 + t240 + (-pkin(1) + t230) * t250) + t252 * (-t252 * pkin(1) + t218 + (-pkin(7) - t253) * t250);
t259 = rSges(4,1) * t310 - rSges(4,2) * t312 + t250 * rSges(4,3);
t275 = rSges(4,1) * t236 - rSges(4,2) * t235;
t124 = t250 * (-t252 * rSges(4,3) + t250 * t275) + t252 * t259;
t164 = -t236 * rSges(5,3) + (rSges(5,1) * t247 - rSges(5,2) * t246) * t235;
t301 = -t164 - t205;
t295 = pkin(3) * t310 + qJ(4) * t312;
t300 = t243 * (t316 + t331) + t252 * t295;
t296 = pkin(4) * t308 + t248 * t313;
t294 = t250 * rSges(3,3) + t252 * t329;
t292 = t243 + t244;
t152 = -t236 * rSges(6,3) + (rSges(6,1) * t233 - rSges(6,2) * t232) * t235;
t291 = -t152 + t303;
t121 = t185 * rSges(6,1) + t184 * rSges(6,2) + rSges(6,3) * t312;
t290 = t198 * rSges(5,1) + t197 * rSges(5,2) + rSges(5,3) * t312;
t289 = t313 / 0.2e1;
t288 = t312 / 0.2e1;
t287 = -t205 - t332;
t206 = rSges(4,1) * t235 + rSges(4,2) * t236;
t286 = -t206 - t332;
t17 = t250 * t36 - t252 * t35;
t39 = t114 * t312 + t116 * t184 + t118 * t185;
t40 = t115 * t312 + t117 * t184 + t119 * t185;
t22 = t250 * t40 - t252 * t39;
t285 = (t17 + t22 + t243 * t172 + (t130 * t312 + t132 * t197 + t134 * t198) * t250 + (-t129 * t312 - t131 * t197 - t133 * t198 + (-t171 + t264) * t250 + t265 * t252) * t252) * t250;
t284 = (t45 + t60) * t289 + (t46 + t61) * t288;
t283 = -t250 * t253 + t218;
t11 = -t322 + (t250 * t45 + t252 * t46) * t235;
t282 = -t236 * t11 + t336;
t258 = -t248 * t312 - t298;
t281 = t250 * ((t330 * t236 - t316) * t250 - t296) + t252 * (t258 - t295) + t300;
t280 = t303 + t304;
t274 = -t196 * rSges(5,1) - t195 * rSges(5,2);
t62 = t250 * (rSges(5,3) * t313 - t274) + t252 * t290 + t300;
t279 = t16 * t289 + t17 * t288 + t5 * t333 + t6 * t334 + (t325 - t326) * t335;
t278 = -t144 + t287;
t277 = -t164 + t287;
t276 = -t328 + t329;
t273 = -t183 * rSges(6,1) - t182 * rSges(6,2);
t271 = Icges(3,1) * t251 - t320;
t269 = -Icges(3,2) * t249 + t319;
t267 = Icges(3,5) * t251 - Icges(3,6) * t249;
t203 = Icges(4,2) * t236 + t318;
t204 = Icges(4,1) * t235 + t317;
t261 = -t203 * t235 + t204 * t236;
t260 = -t152 + t278;
t120 = rSges(6,3) * t313 - t273;
t30 = t250 * t120 + t252 * t121 + t281;
t257 = t278 + t304;
t94 = -t252 * t208 + (-t235 * t241 + t236 * t297) * t250 + t296;
t23 = t281 + t341 * t252 + (t110 + t94) * t250;
t149 = -Icges(6,3) * t236 + (Icges(6,5) * t233 - Icges(6,6) * t232) * t235;
t151 = -Icges(6,5) * t236 + (Icges(6,1) * t233 - Icges(6,4) * t232) * t235;
t68 = t149 * t312 + t150 * t184 + t151 * t185;
t10 = -t68 * t236 + (t250 * t39 + t252 * t40) * t235;
t67 = t149 * t313 + t150 * t182 + t151 * t183;
t9 = -t67 * t236 + (t250 * t37 + t252 * t38) * t235;
t256 = t10 * t334 + t21 * t289 + t22 * t288 + t9 * t333 + (t323 - t324) * t335 + t279;
t255 = t340 * t252 + t285;
t159 = -Icges(5,3) * t236 + (Icges(5,5) * t247 - Icges(5,6) * t246) * t235;
t160 = -Icges(5,6) * t236 + (Icges(5,4) * t247 - Icges(5,2) * t246) * t235;
t161 = -Icges(5,5) * t236 + (Icges(5,1) * t247 - Icges(5,4) * t246) * t235;
t202 = Icges(4,5) * t235 + Icges(4,6) * t236;
t254 = -t326 / 0.2e1 + t325 / 0.2e1 - t324 / 0.2e1 + t323 / 0.2e1 + (t159 * t312 + t160 * t197 + t161 * t198 + t250 * t202 + t252 * t261 + t61 + t68 + (t174 - t130) * t236 + (-t132 * t246 + t134 * t247 + t176) * t235) * t334 + (t159 * t313 + t160 * t195 + t161 * t196 - t252 * t202 + t250 * t261 + t60 + t67 + (t173 - t129) * t236 + (-t131 * t246 + t133 * t247 + t175) * t235) * t333;
t216 = rSges(2,1) * t252 - rSges(2,2) * t250;
t215 = -rSges(2,1) * t250 - rSges(2,2) * t252;
t214 = rSges(3,1) * t249 + rSges(3,2) * t251;
t188 = Icges(3,3) * t250 + t252 * t267;
t187 = -Icges(3,3) * t252 + t250 * t267;
t168 = t286 * t252;
t165 = t286 * t250;
t154 = t250 * pkin(7) + (pkin(1) - t328) * t252 + t294;
t153 = t327 + t240 + (-pkin(1) - t276) * t250;
t143 = t259 + t283;
t142 = (rSges(4,3) - t253) * t252 + (-t230 - t275) * t250;
t140 = t235 * t233 * t151;
t138 = t252 * (-t252 * t328 + t294) + (t250 * t276 - t327) * t250;
t137 = t301 * t252;
t136 = t301 * t250;
t126 = t277 * t252;
t125 = t277 * t250;
t96 = t110 * t312;
t93 = t283 + t290 + t295;
t92 = -t305 + (-t331 - t230 + (-rSges(5,3) - qJ(4)) * t235) * t250 + t274;
t89 = t291 * t252;
t88 = t291 * t250;
t87 = t260 * t252;
t86 = t260 * t250;
t85 = -t121 * t236 - t152 * t312;
t84 = t120 * t236 + t152 * t313;
t83 = t258 + t283 + t121;
t82 = -t305 + (-rSges(6,3) * t235 - t229 * t236 - t230) * t250 + t273 + t296;
t81 = -t236 * t111 - t148 * t312;
t79 = t124 + t302;
t78 = -t236 * t149 - t235 * t314 + t140;
t77 = -t241 * t312 + t283 + t342;
t76 = (t208 - t253) * t252 + (-t207 * t236 - t230 + (-rSges(7,3) + t241) * t235) * t250 + t272;
t72 = (t120 * t252 - t121 * t250) * t235;
t71 = t280 * t252;
t70 = t280 * t250;
t69 = -t111 * t313 + t96;
t66 = t257 * t252;
t65 = t257 * t250;
t47 = t62 + t302;
t32 = -t236 * t341 + t304 * t312;
t31 = t135 * t313 + t236 * t94 + t80;
t29 = t96 + (-t250 * t341 + t252 * t94) * t235;
t28 = t30 + t302;
t14 = t23 + t302;
t1 = [t251 * (Icges(3,2) * t251 + t320) + t249 * (Icges(3,1) * t249 + t319) + Icges(2,3) + t139 + t140 + (-t145 - t149 - t159 + t203) * t236 + (-t160 * t246 + t161 * t247 + t204 - t314 - t315) * t235 + m(7) * (t76 ^ 2 + t77 ^ 2) + m(6) * (t82 ^ 2 + t83 ^ 2) + m(5) * (t92 ^ 2 + t93 ^ 2) + m(4) * (t142 ^ 2 + t143 ^ 2) + m(3) * (t153 ^ 2 + t154 ^ 2) + m(2) * (t215 ^ 2 + t216 ^ 2); m(7) * (t65 * t77 + t66 * t76) + m(6) * (t82 * t87 + t83 * t86) + m(5) * (t125 * t93 + t126 * t92) + m(4) * (t142 * t168 + t143 * t165) + m(3) * (-t153 * t252 - t154 * t250) * t214 + (t243 / 0.2e1 + t244 / 0.2e1) * (Icges(3,5) * t249 + Icges(3,6) * t251) + t254 + (t251 * (Icges(3,6) * t250 + t252 * t269) + t249 * (Icges(3,5) * t250 + t252 * t271)) * t334 + (t251 * (-Icges(3,6) * t252 + t250 * t269) + t249 * (-Icges(3,5) * t252 + t250 * t271)) * t333; m(7) * (t14 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(6) * (t28 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(5) * (t125 ^ 2 + t126 ^ 2 + t47 ^ 2) + m(4) * (t165 ^ 2 + t168 ^ 2 + t79 ^ 2) + m(3) * (t214 ^ 2 * t292 + t138 ^ 2) + t250 * t243 * t188 + t285 + (-t244 * t187 + (-t250 * t187 + t252 * t188) * t250 + t340) * t252; m(7) * (t70 * t77 + t71 * t76) + m(6) * (t82 * t89 + t83 * t88) + m(5) * (t136 * t93 + t137 * t92) + m(4) * (-t142 * t252 - t143 * t250) * t206 + t254; m(7) * (t14 * t23 + t65 * t70 + t66 * t71) + m(6) * (t28 * t30 + t86 * t88 + t87 * t89) + m(5) * (t125 * t136 + t126 * t137 + t47 * t62) + m(4) * (t124 * t79 + (-t165 * t250 - t168 * t252) * t206) + t255; m(7) * (t23 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(6) * (t30 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(5) * (t136 ^ 2 + t137 ^ 2 + t62 ^ 2) + m(4) * (t206 ^ 2 * t292 + t124 ^ 2) + t255; 0.2e1 * ((t250 * t77 + t252 * t76) * t337 + (t250 * t83 + t252 * t82) * t338 + (t250 * t93 + t252 * t92) * t339) * t235; m(7) * (-t236 * t14 + (t250 * t65 + t252 * t66) * t235) + m(6) * (-t236 * t28 + (t250 * t86 + t252 * t87) * t235) + m(5) * (-t236 * t47 + (t125 * t250 + t126 * t252) * t235); m(7) * (-t236 * t23 + (t250 * t70 + t252 * t71) * t235) + m(6) * (-t236 * t30 + (t250 * t88 + t252 * t89) * t235) + m(5) * (-t236 * t62 + (t136 * t250 + t137 * t252) * t235); 0.2e1 * (t339 + t338 + t337) * (t235 ^ 2 * t292 + t236 ^ 2); (-t75 - t78) * t236 + m(7) * (t31 * t76 + t32 * t77) + m(6) * (t82 * t84 + t83 * t85) + ((t51 / 0.2e1 + t68 / 0.2e1) * t252 + (t50 / 0.2e1 + t67 / 0.2e1) * t250) * t235 + t284; m(7) * (t14 * t29 + t31 * t66 + t32 * t65) + m(6) * (t28 * t72 + t84 * t87 + t85 * t86) + t256; m(7) * (t23 * t29 + t31 * t71 + t32 * t70) + m(6) * (t30 * t72 + t84 * t89 + t85 * t88) + t256; m(6) * (-t72 * t236 + (t250 * t85 + t252 * t84) * t235) + m(7) * (-t29 * t236 + (t250 * t32 + t252 * t31) * t235); (t78 * t236 - t11) * t236 + m(7) * (t29 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t72 ^ 2 + t84 ^ 2 + t85 ^ 2) + (t252 * t10 + t250 * t9 - t236 * (t250 * t50 + t252 * t51)) * t235 + t336; -t322 + m(7) * (t76 * t80 + t77 * t81) + t284; m(7) * (t14 * t69 + t65 * t81 + t66 * t80) + t279; m(7) * (t23 * t69 + t70 * t81 + t71 * t80) + t279; m(7) * (-t69 * t236 + (t250 * t81 + t252 * t80) * t235); m(7) * (t29 * t69 + t31 * t80 + t32 * t81) + t282; m(7) * (t69 ^ 2 + t80 ^ 2 + t81 ^ 2) + t282;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
