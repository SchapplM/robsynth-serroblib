% Calculate joint inertia matrix for
% S6RRRPRR2
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
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:06:29
% EndTime: 2019-03-09 18:06:40
% DurationCPUTime: 4.92s
% Computational Cost: add. (14320->470), mult. (11159->665), div. (0->0), fcn. (11913->12), ass. (0->242)
t346 = Icges(4,3) + Icges(5,3);
t234 = qJ(2) + qJ(3);
t220 = pkin(11) + t234;
t216 = sin(t220);
t217 = cos(t220);
t222 = sin(t234);
t224 = cos(t234);
t345 = Icges(4,5) * t224 + Icges(5,5) * t217 - Icges(4,6) * t222 - Icges(5,6) * t216;
t233 = qJ(5) + qJ(6);
t221 = sin(t233);
t240 = cos(qJ(1));
t300 = t240 * t221;
t223 = cos(t233);
t237 = sin(qJ(1));
t303 = t237 * t223;
t172 = -t217 * t300 + t303;
t299 = t240 * t223;
t304 = t237 * t221;
t173 = t217 * t299 + t304;
t308 = t216 * t240;
t105 = rSges(7,1) * t173 + rSges(7,2) * t172 + rSges(7,3) * t308;
t238 = cos(qJ(5));
t218 = pkin(5) * t238 + pkin(4);
t241 = -pkin(10) - pkin(9);
t235 = sin(qJ(5));
t302 = t237 * t235;
t307 = t217 * t240;
t344 = pkin(5) * t302 + t218 * t307 - t241 * t308 + t105;
t343 = t237 * t346 + t240 * t345;
t342 = -t237 * t345 + t240 * t346;
t310 = Icges(5,4) * t217;
t265 = -Icges(5,2) * t216 + t310;
t153 = Icges(5,6) * t237 + t240 * t265;
t311 = Icges(5,4) * t216;
t268 = Icges(5,1) * t217 - t311;
t155 = Icges(5,5) * t237 + t240 * t268;
t312 = Icges(4,4) * t224;
t266 = -Icges(4,2) * t222 + t312;
t167 = Icges(4,6) * t237 + t240 * t266;
t313 = Icges(4,4) * t222;
t269 = Icges(4,1) * t224 - t313;
t169 = Icges(4,5) * t237 + t240 * t269;
t341 = t153 * t216 - t155 * t217 + t167 * t222 - t169 * t224;
t152 = -Icges(5,6) * t240 + t237 * t265;
t154 = -Icges(5,5) * t240 + t237 * t268;
t166 = -Icges(4,6) * t240 + t237 * t266;
t168 = -Icges(4,5) * t240 + t237 * t269;
t340 = t152 * t216 - t154 * t217 + t166 * t222 - t168 * t224;
t231 = t237 ^ 2;
t339 = t237 * pkin(7);
t292 = pkin(4) * t307 + pkin(9) * t308;
t338 = -t292 + t344;
t337 = Icges(4,5) * t222 + Icges(5,5) * t216 + Icges(4,6) * t224 + Icges(5,6) * t217;
t187 = Icges(5,2) * t217 + t311;
t188 = Icges(5,1) * t216 + t310;
t195 = Icges(4,2) * t224 + t313;
t196 = Icges(4,1) * t222 + t312;
t336 = -t187 * t216 + t188 * t217 - t195 * t222 + t196 * t224;
t274 = rSges(4,1) * t224 - rSges(4,2) * t222;
t170 = -t217 * t304 - t299;
t171 = t217 * t303 - t300;
t309 = t216 * t237;
t100 = Icges(7,4) * t171 + Icges(7,2) * t170 + Icges(7,6) * t309;
t102 = Icges(7,1) * t171 + Icges(7,4) * t170 + Icges(7,5) * t309;
t98 = Icges(7,5) * t171 + Icges(7,6) * t170 + Icges(7,3) * t309;
t29 = t100 * t170 + t102 * t171 + t309 * t98;
t101 = Icges(7,4) * t173 + Icges(7,2) * t172 + Icges(7,6) * t308;
t103 = Icges(7,1) * t173 + Icges(7,4) * t172 + Icges(7,5) * t308;
t99 = Icges(7,5) * t173 + Icges(7,6) * t172 + Icges(7,3) * t308;
t30 = t101 * t170 + t103 * t171 + t309 * t99;
t15 = t237 * t30 - t240 * t29;
t297 = t240 * t238;
t182 = -t217 * t302 - t297;
t298 = t240 * t235;
t301 = t237 * t238;
t183 = t217 * t301 - t298;
t113 = Icges(6,5) * t183 + Icges(6,6) * t182 + Icges(6,3) * t309;
t115 = Icges(6,4) * t183 + Icges(6,2) * t182 + Icges(6,6) * t309;
t117 = Icges(6,1) * t183 + Icges(6,4) * t182 + Icges(6,5) * t309;
t43 = t113 * t309 + t115 * t182 + t117 * t183;
t184 = -t217 * t298 + t301;
t185 = t217 * t297 + t302;
t114 = Icges(6,5) * t185 + Icges(6,6) * t184 + Icges(6,3) * t308;
t116 = Icges(6,4) * t185 + Icges(6,2) * t184 + Icges(6,6) * t308;
t118 = Icges(6,1) * t185 + Icges(6,4) * t184 + Icges(6,5) * t308;
t44 = t114 * t309 + t116 * t182 + t118 * t183;
t22 = t237 * t44 - t240 * t43;
t232 = t240 ^ 2;
t335 = -t15 - t22 + t342 * t232 + (t341 * t237 + (-t340 + t343) * t240) * t237;
t242 = -pkin(8) - pkin(7);
t134 = -Icges(7,3) * t217 + (Icges(7,5) * t223 - Icges(7,6) * t221) * t216;
t135 = -Icges(7,6) * t217 + (Icges(7,4) * t223 - Icges(7,2) * t221) * t216;
t136 = -Icges(7,5) * t217 + (Icges(7,1) * t223 - Icges(7,4) * t221) * t216;
t58 = t134 * t309 + t135 * t170 + t136 * t171;
t5 = -t58 * t217 + (t237 * t29 + t240 * t30) * t216;
t31 = t100 * t172 + t102 * t173 + t308 * t98;
t32 = t101 * t172 + t103 * t173 + t308 * t99;
t59 = t134 * t308 + t135 * t172 + t136 * t173;
t6 = -t59 * t217 + (t237 * t31 + t240 * t32) * t216;
t334 = t308 * t6 + t309 * t5;
t333 = -t217 / 0.2e1;
t332 = t237 / 0.2e1;
t331 = -t240 / 0.2e1;
t236 = sin(qJ(2));
t330 = pkin(2) * t236;
t329 = pkin(3) * t222;
t328 = pkin(4) * t217;
t327 = -pkin(4) + t218;
t326 = pkin(9) + t241;
t239 = cos(qJ(2));
t219 = pkin(2) * t239 + pkin(1);
t325 = rSges(3,1) * t239;
t323 = rSges(3,2) * t236;
t321 = t240 * rSges(3,3);
t41 = -t217 * t98 + (-t100 * t221 + t102 * t223) * t216;
t320 = t41 * t240;
t42 = -t217 * t99 + (-t101 * t221 + t103 * t223) * t216;
t319 = t42 * t237;
t52 = -t217 * t113 + (-t115 * t235 + t117 * t238) * t216;
t318 = t52 * t240;
t53 = -t217 * t114 + (-t116 * t235 + t118 * t238) * t216;
t317 = t53 * t237;
t124 = t216 * t223 * t136;
t306 = t221 * t135;
t66 = -t217 * t134 - t216 * t306 + t124;
t316 = t66 * t217;
t271 = -t171 * rSges(7,1) - t170 * rSges(7,2);
t104 = rSges(7,3) * t309 - t271;
t137 = -t217 * rSges(7,3) + (rSges(7,1) * t223 - rSges(7,2) * t221) * t216;
t73 = t104 * t217 + t137 * t309;
t315 = Icges(3,4) * t236;
t314 = Icges(3,4) * t239;
t139 = -Icges(6,6) * t217 + (Icges(6,4) * t238 - Icges(6,2) * t235) * t216;
t305 = t235 * t139;
t198 = pkin(3) * t224 + t219;
t192 = t240 * t198;
t211 = t240 * t219;
t295 = t240 * (t192 - t211) + (t198 - t219) * t231;
t129 = t216 * t327 + t217 * t326;
t294 = -t129 - t137;
t229 = t240 * pkin(7);
t293 = t237 * (t229 + (-pkin(1) + t219) * t237) + t240 * (-t240 * pkin(1) + t211 - t339);
t252 = t237 * rSges(4,3) + t240 * t274;
t112 = t237 * (-t240 * rSges(4,3) + t237 * t274) + t240 * t252;
t291 = rSges(3,3) * t237 + t240 * t325;
t289 = t231 + t232;
t120 = rSges(6,1) * t185 + rSges(6,2) * t184 + rSges(6,3) * t308;
t288 = t309 / 0.2e1;
t287 = t308 / 0.2e1;
t197 = rSges(4,1) * t222 + rSges(4,2) * t224;
t286 = -t197 - t330;
t189 = rSges(5,1) * t216 + rSges(5,2) * t217;
t285 = -t189 - t329;
t190 = pkin(4) * t216 - t217 * pkin(9);
t284 = -t190 - t329;
t16 = t237 * t32 - t240 * t31;
t45 = t113 * t308 + t115 * t184 + t117 * t185;
t46 = t114 * t308 + t116 * t184 + t118 * t185;
t23 = t237 * t46 - t240 * t45;
t283 = (t16 + t23 + t343 * t231 + ((-t341 + t342) * t237 + t340 * t240) * t240) * t237;
t282 = (t41 + t58) * t288 + (t42 + t59) * t287;
t230 = -qJ(4) + t242;
t281 = -t230 * t237 + t192;
t9 = -t316 + (t237 * t41 + t240 * t42) * t216;
t280 = -t217 * t9 + t334;
t251 = rSges(5,1) * t307 - rSges(5,2) * t308 + rSges(5,3) * t237;
t273 = rSges(5,1) * t217 - rSges(5,2) * t216;
t62 = t237 * (-t240 * rSges(5,3) + t237 * t273) + t240 * t251 + t295;
t279 = t231 * (pkin(9) * t216 + t328) + t240 * t292 + t295;
t278 = t15 * t288 + t16 * t287 + t5 * t331 + t6 * t332 + (t319 - t320) * t333;
t141 = -t217 * rSges(6,3) + (rSges(6,1) * t238 - rSges(6,2) * t235) * t216;
t277 = -t141 + t284;
t276 = -t329 - t330;
t275 = -t323 + t325;
t272 = -t183 * rSges(6,1) - t182 * rSges(6,2);
t270 = Icges(3,1) * t239 - t315;
t267 = -Icges(3,2) * t236 + t314;
t264 = Icges(3,5) * t239 - Icges(3,6) * t236;
t253 = t284 + t294;
t250 = -t189 + t276;
t249 = -t190 + t276;
t119 = rSges(6,3) * t309 - t272;
t28 = t119 * t237 + t120 * t240 + t279;
t247 = -t141 + t249;
t246 = t249 + t294;
t110 = -pkin(5) * t298 + (-t216 * t326 + t217 * t327) * t237;
t25 = t279 + t338 * t240 + (t104 + t110) * t237;
t138 = -Icges(6,3) * t217 + (Icges(6,5) * t238 - Icges(6,6) * t235) * t216;
t140 = -Icges(6,5) * t217 + (Icges(6,1) * t238 - Icges(6,4) * t235) * t216;
t63 = t138 * t309 + t139 * t182 + t140 * t183;
t10 = -t63 * t217 + (t237 * t43 + t240 * t44) * t216;
t64 = t138 * t308 + t139 * t184 + t140 * t185;
t11 = -t64 * t217 + (t237 * t45 + t240 * t46) * t216;
t245 = t10 * t331 + t11 * t332 + t22 * t288 + t23 * t287 + (t317 - t318) * t333 + t278;
t244 = t240 * t335 + t283;
t243 = -t320 / 0.2e1 + t319 / 0.2e1 - t318 / 0.2e1 + t317 / 0.2e1 + (t153 * t217 + t155 * t216 + t224 * t167 + t222 * t169 + t237 * t337 + t240 * t336 + t59 + t64) * t332 + (t152 * t217 + t154 * t216 + t224 * t166 + t222 * t168 + t237 * t336 - t240 * t337 + t58 + t63) * t331;
t210 = rSges(2,1) * t240 - rSges(2,2) * t237;
t209 = -rSges(2,1) * t237 - rSges(2,2) * t240;
t208 = rSges(3,1) * t236 + rSges(3,2) * t239;
t177 = Icges(3,3) * t237 + t240 * t264;
t176 = -Icges(3,3) * t240 + t237 * t264;
t163 = t286 * t240;
t162 = t286 * t237;
t147 = t339 + (pkin(1) - t323) * t240 + t291;
t146 = t321 + t229 + (-pkin(1) - t275) * t237;
t145 = t285 * t240;
t144 = t285 * t237;
t133 = t250 * t240;
t132 = t250 * t237;
t131 = -t237 * t242 + t211 + t252;
t130 = (rSges(4,3) - t242) * t240 + (-t219 - t274) * t237;
t128 = t216 * t238 * t140;
t123 = t240 * (-t240 * t323 + t291) + (t237 * t275 - t321) * t237;
t122 = t251 + t281;
t121 = (rSges(5,3) - t230) * t240 + (-t198 - t273) * t237;
t97 = t277 * t240;
t96 = t277 * t237;
t88 = t247 * t240;
t87 = t247 * t237;
t86 = t104 * t308;
t81 = t281 + t120 + t292;
t80 = -t240 * t230 + (-t328 - t198 + (-rSges(6,3) - pkin(9)) * t216) * t237 + t272;
t79 = -t120 * t217 - t141 * t308;
t78 = t119 * t217 + t141 * t309;
t77 = t253 * t240;
t76 = t253 * t237;
t75 = t112 + t293;
t74 = -t217 * t105 - t137 * t308;
t72 = t246 * t240;
t71 = t246 * t237;
t70 = t281 + t344;
t69 = (pkin(5) * t235 - t230) * t240 + (-t217 * t218 - t198 + (-rSges(7,3) + t241) * t216) * t237 + t271;
t68 = -t217 * t138 - t216 * t305 + t128;
t67 = (t119 * t240 - t120 * t237) * t216;
t65 = -t105 * t309 + t86;
t47 = t62 + t293;
t36 = -t217 * t338 + t294 * t308;
t35 = t110 * t217 + t129 * t309 + t73;
t27 = t86 + (t110 * t240 - t237 * t338) * t216;
t26 = t28 + t293;
t20 = t25 + t293;
t1 = [t224 * t195 + t222 * t196 + t239 * (Icges(3,2) * t239 + t315) + t236 * (Icges(3,1) * t236 + t314) + Icges(2,3) + t124 + t128 + (-t134 - t138 + t187) * t217 + (t188 - t305 - t306) * t216 + m(7) * (t69 ^ 2 + t70 ^ 2) + m(6) * (t80 ^ 2 + t81 ^ 2) + m(5) * (t121 ^ 2 + t122 ^ 2) + m(4) * (t130 ^ 2 + t131 ^ 2) + m(3) * (t146 ^ 2 + t147 ^ 2) + m(2) * (t209 ^ 2 + t210 ^ 2); (t239 * (-Icges(3,6) * t240 + t237 * t267) + t236 * (-Icges(3,5) * t240 + t237 * t270)) * t331 + (t239 * (Icges(3,6) * t237 + t240 * t267) + t236 * (Icges(3,5) * t237 + t240 * t270)) * t332 + m(3) * (-t146 * t240 - t147 * t237) * t208 + m(7) * (t69 * t72 + t70 * t71) + m(6) * (t80 * t88 + t81 * t87) + m(5) * (t121 * t133 + t122 * t132) + m(4) * (t130 * t163 + t131 * t162) + t243 + (t231 / 0.2e1 + t232 / 0.2e1) * (Icges(3,5) * t236 + Icges(3,6) * t239); m(7) * (t20 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(6) * (t26 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(5) * (t132 ^ 2 + t133 ^ 2 + t47 ^ 2) + m(4) * (t162 ^ 2 + t163 ^ 2 + t75 ^ 2) + m(3) * (t208 ^ 2 * t289 + t123 ^ 2) + t237 * t231 * t177 + t283 + (-t232 * t176 + (-t237 * t176 + t240 * t177) * t237 + t335) * t240; m(7) * (t69 * t77 + t70 * t76) + m(6) * (t80 * t97 + t81 * t96) + m(5) * (t121 * t145 + t122 * t144) + t243 + m(4) * (-t130 * t240 - t131 * t237) * t197; m(7) * (t20 * t25 + t71 * t76 + t72 * t77) + m(6) * (t26 * t28 + t87 * t96 + t88 * t97) + m(5) * (t132 * t144 + t133 * t145 + t47 * t62) + m(4) * (t112 * t75 + (-t162 * t237 - t163 * t240) * t197) + t244; m(7) * (t25 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(6) * (t28 ^ 2 + t96 ^ 2 + t97 ^ 2) + m(5) * (t144 ^ 2 + t145 ^ 2 + t62 ^ 2) + m(4) * (t197 ^ 2 * t289 + t112 ^ 2) + t244; m(7) * (t237 * t69 - t240 * t70) + m(6) * (t237 * t80 - t240 * t81) + m(5) * (t121 * t237 - t122 * t240); m(7) * (t237 * t72 - t240 * t71) + m(6) * (t237 * t88 - t240 * t87) + m(5) * (-t132 * t240 + t133 * t237); m(7) * (t237 * t77 - t240 * t76) + m(6) * (t237 * t97 - t240 * t96) + m(5) * (-t144 * t240 + t145 * t237); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t289; (-t66 - t68) * t217 + m(7) * (t35 * t69 + t36 * t70) + m(6) * (t78 * t80 + t79 * t81) + ((t53 / 0.2e1 + t64 / 0.2e1) * t240 + (t52 / 0.2e1 + t63 / 0.2e1) * t237) * t216 + t282; m(7) * (t20 * t27 + t35 * t72 + t36 * t71) + m(6) * (t26 * t67 + t78 * t88 + t79 * t87) + t245; m(7) * (t25 * t27 + t35 * t77 + t36 * t76) + m(6) * (t28 * t67 + t78 * t97 + t79 * t96) + t245; m(6) * (t237 * t78 - t240 * t79) + m(7) * (t237 * t35 - t240 * t36); (t68 * t217 - t9) * t217 + m(7) * (t27 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(6) * (t67 ^ 2 + t78 ^ 2 + t79 ^ 2) + (t237 * t10 + t240 * t11 - t217 * (t237 * t52 + t240 * t53)) * t216 + t334; m(7) * (t69 * t73 + t70 * t74) - t316 + t282; m(7) * (t20 * t65 + t71 * t74 + t72 * t73) + t278; m(7) * (t25 * t65 + t73 * t77 + t74 * t76) + t278; m(7) * (t237 * t73 - t240 * t74); m(7) * (t27 * t65 + t35 * t73 + t36 * t74) + t280; m(7) * (t65 ^ 2 + t73 ^ 2 + t74 ^ 2) + t280;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
