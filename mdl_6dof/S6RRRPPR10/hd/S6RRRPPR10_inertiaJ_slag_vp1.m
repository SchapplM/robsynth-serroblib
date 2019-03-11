% Calculate joint inertia matrix for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR10_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR10_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:21:31
% EndTime: 2019-03-09 16:21:42
% DurationCPUTime: 4.44s
% Computational Cost: add. (16878->641), mult. (37638->881), div. (0->0), fcn. (47451->12), ass. (0->290)
t315 = m(6) / 0.2e1 + m(7) / 0.2e1;
t351 = 0.2e1 * t315;
t279 = cos(pkin(6));
t281 = sin(qJ(2));
t277 = sin(pkin(6));
t343 = sin(qJ(3));
t304 = t277 * t343;
t344 = cos(qJ(3));
t255 = -t279 * t344 + t281 * t304;
t275 = pkin(11) + qJ(6);
t272 = sin(t275);
t273 = cos(t275);
t283 = cos(qJ(2));
t336 = t277 * t283;
t220 = t255 * t273 + t272 * t336;
t221 = t255 * t272 - t273 * t336;
t305 = t277 * t344;
t256 = t279 * t343 + t281 * t305;
t134 = Icges(7,5) * t221 + Icges(7,6) * t220 + Icges(7,3) * t256;
t135 = Icges(7,4) * t221 + Icges(7,2) * t220 + Icges(7,6) * t256;
t136 = Icges(7,1) * t221 + Icges(7,4) * t220 + Icges(7,5) * t256;
t59 = t134 * t256 + t135 * t220 + t136 * t221;
t276 = sin(pkin(11));
t278 = cos(pkin(11));
t233 = t255 * t278 + t276 * t336;
t338 = t255 * t276;
t234 = -t278 * t336 + t338;
t139 = Icges(6,5) * t234 + Icges(6,6) * t233 + Icges(6,3) * t256;
t140 = Icges(6,4) * t234 + Icges(6,2) * t233 + Icges(6,6) * t256;
t141 = Icges(6,1) * t234 + Icges(6,4) * t233 + Icges(6,5) * t256;
t65 = t139 * t256 + t140 * t233 + t141 * t234;
t350 = -t59 - t65;
t284 = cos(qJ(1));
t335 = t277 * t284;
t331 = t283 * t284;
t282 = sin(qJ(1));
t334 = t281 * t282;
t260 = -t279 * t334 + t331;
t237 = t260 * t343 - t282 * t305;
t332 = t282 * t283;
t333 = t281 * t284;
t259 = t279 * t332 + t333;
t180 = t237 * t273 - t259 * t272;
t181 = t237 * t272 + t259 * t273;
t238 = t260 * t344 + t282 * t304;
t110 = rSges(7,1) * t181 + rSges(7,2) * t180 + rSges(7,3) * t238;
t271 = pkin(5) * t278 + pkin(4);
t280 = -pkin(10) - qJ(5);
t339 = t237 * t276;
t349 = pkin(5) * t339 - t238 * t280 + t259 * t271 + t110;
t258 = t279 * t333 + t332;
t236 = t258 * t344 - t284 * t304;
t235 = t258 * t343 + t284 * t305;
t257 = -t279 * t331 + t334;
t178 = t235 * t273 - t257 * t272;
t179 = t235 * t272 + t257 * t273;
t101 = Icges(7,5) * t179 + Icges(7,6) * t178 + Icges(7,3) * t236;
t103 = Icges(7,4) * t179 + Icges(7,2) * t178 + Icges(7,6) * t236;
t105 = Icges(7,1) * t179 + Icges(7,4) * t178 + Icges(7,5) * t236;
t38 = t101 * t238 + t103 * t180 + t105 * t181;
t102 = Icges(7,5) * t181 + Icges(7,6) * t180 + Icges(7,3) * t238;
t104 = Icges(7,4) * t181 + Icges(7,2) * t180 + Icges(7,6) * t238;
t106 = Icges(7,1) * t181 + Icges(7,4) * t180 + Icges(7,5) * t238;
t39 = t102 * t238 + t104 * t180 + t106 * t181;
t53 = t134 * t238 + t135 * t180 + t136 * t181;
t2 = t236 * t38 + t238 * t39 + t256 * t53;
t348 = t2 / 0.2e1;
t347 = t236 / 0.2e1;
t346 = t238 / 0.2e1;
t345 = t256 / 0.2e1;
t342 = pkin(4) - t271;
t205 = Icges(3,5) * t258 - Icges(3,6) * t257 - Icges(3,3) * t335;
t341 = t205 * t284;
t340 = t235 * t276;
t337 = t277 * t282;
t294 = -t179 * rSges(7,1) - t178 * rSges(7,2);
t109 = rSges(7,3) * t236 - t294;
t225 = t236 * qJ(5);
t314 = pkin(5) * t340;
t330 = -t236 * t280 - t257 * t342 + t109 - t225 + t314;
t196 = pkin(4) * t259 + qJ(5) * t238;
t329 = -t196 + t349;
t137 = rSges(7,1) * t221 + rSges(7,2) * t220 + rSges(7,3) * t256;
t328 = t137 + pkin(5) * t338 + t342 * t336 + (-qJ(5) - t280) * t256;
t156 = rSges(5,1) * t259 - rSges(5,2) * t238 + rSges(5,3) * t237;
t171 = pkin(3) * t238 + qJ(4) * t237;
t327 = -t156 - t171;
t224 = t235 * qJ(4);
t170 = pkin(3) * t236 + t224;
t160 = t259 * t170;
t195 = t257 * pkin(4) + t225;
t326 = t195 * t259 + t160;
t217 = pkin(3) * t256 + qJ(4) * t255;
t325 = t170 * t336 + t217 * t257;
t219 = pkin(2) * t260 + pkin(9) * t259;
t216 = t279 * t219;
t324 = t171 * t279 + t216;
t218 = t258 * pkin(2) + pkin(9) * t257;
t323 = -t170 - t218;
t322 = -t171 - t196;
t197 = -Icges(5,5) * t336 - Icges(5,6) * t256 + Icges(5,3) * t255;
t198 = -Icges(5,4) * t336 - Icges(5,2) * t256 + Icges(5,6) * t255;
t321 = t197 * t255 - t198 * t256;
t201 = Icges(4,4) * t256 - Icges(4,2) * t255 - Icges(4,6) * t336;
t202 = Icges(4,1) * t256 - Icges(4,4) * t255 - Icges(4,5) * t336;
t320 = -t201 * t255 + t202 * t256;
t203 = -rSges(5,1) * t336 - rSges(5,2) * t256 + rSges(5,3) * t255;
t319 = -t203 - t217;
t318 = t218 * t337 + t219 * t335;
t241 = -pkin(4) * t336 + t256 * qJ(5);
t317 = -t217 - t241;
t316 = pkin(1) * t284 + pkin(8) * t337;
t44 = t101 * t256 + t103 * t220 + t105 * t221;
t52 = t134 * t236 + t135 * t178 + t136 * t179;
t313 = t44 / 0.2e1 + t52 / 0.2e1;
t45 = t102 * t256 + t104 * t220 + t106 * t221;
t312 = t45 / 0.2e1 + t53 / 0.2e1;
t192 = t237 * t278 - t259 * t276;
t193 = t259 * t278 + t339;
t120 = rSges(6,1) * t193 + rSges(6,2) * t192 + rSges(6,3) * t238;
t311 = -t120 + t322;
t142 = rSges(6,1) * t234 + rSges(6,2) * t233 + rSges(6,3) * t256;
t310 = -t142 + t317;
t309 = t196 * t279 + t324;
t308 = -t195 + t323;
t158 = rSges(4,1) * t238 - rSges(4,2) * t237 + rSges(4,3) * t259;
t244 = Icges(3,3) * t279 + (Icges(3,5) * t281 + Icges(3,6) * t283) * t277;
t245 = Icges(3,6) * t279 + (Icges(3,4) * t281 + Icges(3,2) * t283) * t277;
t246 = Icges(3,5) * t279 + (Icges(3,1) * t281 + Icges(3,4) * t283) * t277;
t306 = t246 * t277 * t281 + t244 * t279 + t245 * t336;
t212 = rSges(3,1) * t260 - rSges(3,2) * t259 + rSges(3,3) * t337;
t303 = -t282 * pkin(1) + pkin(8) * t335;
t204 = rSges(4,1) * t256 - rSges(4,2) * t255 - rSges(4,3) * t336;
t261 = (pkin(2) * t281 - pkin(9) * t283) * t277;
t302 = t277 * (-t204 - t261);
t301 = t322 - t329;
t300 = t317 - t328;
t299 = t170 * t337 + t171 * t335 + t318;
t298 = t195 * t336 + t241 * t257 + t325;
t297 = t277 * (-t261 + t319);
t296 = -t257 * rSges(5,1) - t235 * rSges(5,3);
t190 = t235 * t278 - t257 * t276;
t191 = t257 * t278 + t340;
t295 = -t191 * rSges(6,1) - t190 * rSges(6,2);
t293 = t219 + t316;
t292 = t277 * (-t261 + t310);
t291 = t195 * t337 + t196 * t335 + t299;
t290 = t277 * (-t261 + t300);
t289 = -t218 + t303;
t157 = rSges(4,1) * t236 - rSges(4,2) * t235 + rSges(4,3) * t257;
t288 = -t224 + t289;
t211 = rSges(3,1) * t258 - rSges(3,2) * t257 - rSges(3,3) * t335;
t287 = t171 + t293;
t113 = Icges(6,5) * t191 + Icges(6,6) * t190 + Icges(6,3) * t236;
t115 = Icges(6,4) * t191 + Icges(6,2) * t190 + Icges(6,6) * t236;
t117 = Icges(6,1) * t191 + Icges(6,4) * t190 + Icges(6,5) * t236;
t46 = t113 * t256 + t115 * t233 + t117 * t234;
t55 = t139 * t236 + t140 * t190 + t141 * t191;
t143 = Icges(5,5) * t257 - Icges(5,6) * t236 + Icges(5,3) * t235;
t147 = Icges(5,4) * t257 - Icges(5,2) * t236 + Icges(5,6) * t235;
t151 = Icges(5,1) * t257 - Icges(5,4) * t236 + Icges(5,5) * t235;
t80 = t143 * t255 - t147 * t256 - t151 * t336;
t145 = Icges(4,5) * t236 - Icges(4,6) * t235 + Icges(4,3) * t257;
t149 = Icges(4,4) * t236 - Icges(4,2) * t235 + Icges(4,6) * t257;
t153 = Icges(4,1) * t236 - Icges(4,4) * t235 + Icges(4,5) * t257;
t82 = -t145 * t336 - t149 * t255 + t153 * t256;
t199 = -Icges(5,1) * t336 - Icges(5,4) * t256 + Icges(5,5) * t255;
t87 = t197 * t235 - t198 * t236 + t199 * t257;
t200 = Icges(4,5) * t256 - Icges(4,6) * t255 - Icges(4,3) * t336;
t89 = t200 * t257 - t201 * t235 + t202 * t236;
t286 = t46 / 0.2e1 + t89 / 0.2e1 + t87 / 0.2e1 + t55 / 0.2e1 + t82 / 0.2e1 + t80 / 0.2e1 + t313;
t114 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t238;
t116 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t238;
t118 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t238;
t47 = t114 * t256 + t116 * t233 + t118 * t234;
t56 = t139 * t238 + t140 * t192 + t141 * t193;
t144 = Icges(5,5) * t259 - Icges(5,6) * t238 + Icges(5,3) * t237;
t148 = Icges(5,4) * t259 - Icges(5,2) * t238 + Icges(5,6) * t237;
t152 = Icges(5,1) * t259 - Icges(5,4) * t238 + Icges(5,5) * t237;
t81 = t144 * t255 - t148 * t256 - t152 * t336;
t146 = Icges(4,5) * t238 - Icges(4,6) * t237 + Icges(4,3) * t259;
t150 = Icges(4,4) * t238 - Icges(4,2) * t237 + Icges(4,6) * t259;
t154 = Icges(4,1) * t238 - Icges(4,4) * t237 + Icges(4,5) * t259;
t83 = -t146 * t336 - t150 * t255 + t154 * t256;
t88 = t197 * t237 - t198 * t238 + t199 * t259;
t90 = t200 * t259 - t201 * t237 + t202 * t238;
t285 = t56 / 0.2e1 + t83 / 0.2e1 + t81 / 0.2e1 + t47 / 0.2e1 + t90 / 0.2e1 + t88 / 0.2e1 + t312;
t263 = rSges(2,1) * t284 - rSges(2,2) * t282;
t262 = -rSges(2,1) * t282 - rSges(2,2) * t284;
t247 = t279 * rSges(3,3) + (rSges(3,1) * t281 + rSges(3,2) * t283) * t277;
t210 = Icges(3,1) * t260 - Icges(3,4) * t259 + Icges(3,5) * t337;
t209 = Icges(3,1) * t258 - Icges(3,4) * t257 - Icges(3,5) * t335;
t208 = Icges(3,4) * t260 - Icges(3,2) * t259 + Icges(3,6) * t337;
t207 = Icges(3,4) * t258 - Icges(3,2) * t257 - Icges(3,6) * t335;
t206 = Icges(3,5) * t260 - Icges(3,6) * t259 + Icges(3,3) * t337;
t186 = t212 + t316;
t185 = -t211 + t303;
t164 = -t211 * t279 - t247 * t335;
t163 = t212 * t279 - t247 * t337;
t155 = -rSges(5,2) * t236 - t296;
t138 = t306 * t279;
t132 = (t211 * t282 + t212 * t284) * t277;
t131 = t244 * t337 - t245 * t259 + t246 * t260;
t130 = -t244 * t335 - t245 * t257 + t246 * t258;
t124 = t293 + t158;
t123 = -t157 + t289;
t119 = rSges(6,3) * t236 - t295;
t112 = -t158 * t336 - t204 * t259;
t111 = t157 * t336 + t204 * t257;
t108 = t279 * t206 + (t208 * t283 + t210 * t281) * t277;
t107 = t279 * t205 + (t207 * t283 + t209 * t281) * t277;
t99 = -t200 * t336 + t320;
t98 = -t199 * t336 + t321;
t97 = t99 * t279;
t96 = t98 * t279;
t95 = t287 + t156;
t94 = (rSges(5,2) - pkin(3)) * t236 + t288 + t296;
t93 = t157 * t259 - t158 * t257;
t92 = (-t157 - t218) * t279 + t284 * t302;
t91 = t279 * t158 + t282 * t302 + t216;
t86 = (t157 * t282 + t158 * t284) * t277 + t318;
t85 = t259 * t319 + t327 * t336;
t84 = t155 * t336 + t203 * t257 + t325;
t79 = t110 * t256 - t137 * t238;
t78 = -t109 * t256 + t137 * t236;
t77 = t196 + t287 + t120;
t76 = (-rSges(6,3) - pkin(3)) * t236 + t288 + t295 - t195;
t75 = (-t155 + t323) * t279 + t284 * t297;
t74 = t279 * t156 + t282 * t297 + t324;
t73 = t146 * t259 - t150 * t237 + t154 * t238;
t72 = t145 * t259 - t149 * t237 + t153 * t238;
t71 = t146 * t257 - t150 * t235 + t154 * t236;
t70 = t145 * t257 - t149 * t235 + t153 * t236;
t69 = t144 * t237 - t148 * t238 + t152 * t259;
t68 = t143 * t237 - t147 * t238 + t151 * t259;
t67 = t144 * t235 - t148 * t236 + t152 * t257;
t66 = t143 * t235 - t147 * t236 + t151 * t257;
t64 = t65 * t279;
t63 = t287 + t349;
t62 = -t314 - t257 * t271 + (-rSges(7,3) - pkin(3) + t280) * t236 + t288 + t294;
t61 = t259 * t155 + t257 * t327 + t160;
t60 = t109 * t238 - t110 * t236;
t58 = t59 * t279;
t57 = t59 * t256;
t54 = (t155 * t282 + t156 * t284) * t277 + t299;
t51 = t259 * t310 + t311 * t336;
t50 = t119 * t336 + t142 * t257 + t298;
t49 = (-t119 + t308) * t279 + t284 * t292;
t48 = t279 * t120 + t282 * t292 + t309;
t43 = t114 * t238 + t116 * t192 + t118 * t193;
t42 = t113 * t238 + t115 * t192 + t117 * t193;
t41 = t114 * t236 + t116 * t190 + t118 * t191;
t40 = t113 * t236 + t115 * t190 + t117 * t191;
t37 = t102 * t236 + t104 * t178 + t106 * t179;
t36 = t101 * t236 + t103 * t178 + t105 * t179;
t35 = t259 * t119 + t257 * t311 + t326;
t34 = (t119 * t282 + t120 * t284) * t277 + t291;
t33 = t259 * t300 + t301 * t336;
t32 = t257 * t328 + t330 * t336 + t298;
t31 = (t308 - t330) * t279 + t284 * t290;
t30 = t279 * t329 + t282 * t290 + t309;
t29 = t97 + (t83 * t282 - t82 * t284) * t277;
t28 = t96 + (t81 * t282 - t80 * t284) * t277;
t27 = t82 * t257 + t83 * t259 - t336 * t99;
t26 = t80 * t257 + t81 * t259 - t336 * t98;
t25 = t90 * t279 + (t282 * t73 - t284 * t72) * t277;
t24 = t89 * t279 + (t282 * t71 - t284 * t70) * t277;
t23 = t88 * t279 + (t282 * t69 - t284 * t68) * t277;
t22 = t87 * t279 + (t282 * t67 - t284 * t66) * t277;
t21 = t257 * t72 + t259 * t73 - t336 * t90;
t20 = t257 * t70 + t259 * t71 - t336 * t89;
t19 = t257 * t68 + t259 * t69 - t336 * t88;
t18 = t257 * t66 + t259 * t67 - t336 * t87;
t17 = t257 * t301 + t259 * t330 + t326;
t16 = (t282 * t330 + t284 * t329) * t277 + t291;
t15 = t64 + (t47 * t282 - t46 * t284) * t277;
t14 = t46 * t257 + t47 * t259 - t336 * t65;
t13 = t58 + (t45 * t282 - t44 * t284) * t277;
t12 = t44 * t257 + t45 * t259 - t336 * t59;
t11 = t44 * t236 + t45 * t238 + t57;
t10 = t56 * t279 + (t282 * t43 - t284 * t42) * t277;
t9 = t55 * t279 + (t282 * t41 - t284 * t40) * t277;
t8 = t257 * t42 + t259 * t43 - t336 * t56;
t7 = t257 * t40 + t259 * t41 - t336 * t55;
t6 = t53 * t279 + (t282 * t39 - t284 * t38) * t277;
t5 = t52 * t279 + (t282 * t37 - t284 * t36) * t277;
t4 = t257 * t38 + t259 * t39 - t336 * t53;
t3 = t257 * t36 + t259 * t37 - t336 * t52;
t1 = t236 * t36 + t238 * t37 + t256 * t52;
t100 = [t306 + (-t199 - t200) * t336 + m(7) * (t62 ^ 2 + t63 ^ 2) + m(6) * (t76 ^ 2 + t77 ^ 2) + m(5) * (t94 ^ 2 + t95 ^ 2) + m(4) * (t123 ^ 2 + t124 ^ 2) + m(3) * (t185 ^ 2 + t186 ^ 2) + m(2) * (t262 ^ 2 + t263 ^ 2) + Icges(2,3) + t320 + t321 - t350; t58 + t97 + t64 + t96 + t138 + m(7) * (t30 * t63 + t31 * t62) + m(6) * (t48 * t77 + t49 * t76) + m(5) * (t74 * t95 + t75 * t94) + m(4) * (t123 * t92 + t124 * t91) + m(3) * (t163 * t186 + t164 * t185) + ((-t107 / 0.2e1 - t130 / 0.2e1 - t286) * t284 + (t108 / 0.2e1 + t131 / 0.2e1 + t285) * t282) * t277; (t13 + t28 + t29 + t15 + t138) * t279 + m(7) * (t16 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t34 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t54 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(4) * (t86 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(3) * (t132 ^ 2 + t163 ^ 2 + t164 ^ 2) + ((-t5 - t9 - t24 - t22 + (-t207 * t257 + t209 * t258 - t277 * t341) * t335) * t284 + (t6 + t25 + t23 + t10 + ((-t208 * t259 + t210 * t260 + (t206 * t282 - t341) * t277) * t282 + (t206 * t335 + t207 * t259 + t208 * t257 - t209 * t260 - t210 * t258) * t284) * t277) * t282 + ((-t107 - t130) * t284 + (t108 + t131) * t282) * t279) * t277; (-t98 - t99 + t350) * t336 + m(7) * (t32 * t62 + t33 * t63) + m(6) * (t50 * t76 + t51 * t77) + m(5) * (t84 * t94 + t85 * t95) + m(4) * (t111 * t123 + t112 * t124) + t285 * t259 + t286 * t257; (t12 / 0.2e1 + t14 / 0.2e1 + t26 / 0.2e1 + t27 / 0.2e1) * t279 + (t6 / 0.2e1 + t10 / 0.2e1 + t23 / 0.2e1 + t25 / 0.2e1) * t259 + (t5 / 0.2e1 + t9 / 0.2e1 + t22 / 0.2e1 + t24 / 0.2e1) * t257 + m(7) * (t16 * t17 + t30 * t33 + t31 * t32) + m(6) * (t34 * t35 + t48 * t51 + t49 * t50) + m(5) * (t54 * t61 + t74 * t85 + t75 * t84) + m(4) * (t111 * t92 + t112 * t91 + t86 * t93) + ((-t3 / 0.2e1 - t7 / 0.2e1 - t18 / 0.2e1 - t20 / 0.2e1) * t284 + (-t13 / 0.2e1 - t15 / 0.2e1 - t28 / 0.2e1 - t29 / 0.2e1) * t283 + (t4 / 0.2e1 + t8 / 0.2e1 + t19 / 0.2e1 + t21 / 0.2e1) * t282) * t277; (-t12 - t14 - t26 - t27) * t336 + (t4 + t19 + t8 + t21) * t259 + (t3 + t7 + t20 + t18) * t257 + m(7) * (t17 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t35 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t61 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(4) * (t111 ^ 2 + t112 ^ 2 + t93 ^ 2); m(7) * (t235 * t63 + t237 * t62) + m(6) * (t235 * t77 + t237 * t76) + m(5) * (t235 * t95 + t237 * t94); m(7) * (t16 * t255 + t235 * t30 + t237 * t31) + m(6) * (t235 * t48 + t237 * t49 + t255 * t34) + m(5) * (t235 * t74 + t237 * t75 + t255 * t54); m(7) * (t17 * t255 + t235 * t33 + t237 * t32) + m(6) * (t235 * t51 + t237 * t50 + t255 * t35) + m(5) * (t235 * t85 + t237 * t84 + t255 * t61); 0.2e1 * (m(5) / 0.2e1 + t315) * (t235 ^ 2 + t237 ^ 2 + t255 ^ 2); m(7) * (t236 * t63 + t238 * t62) + m(6) * (t236 * t77 + t238 * t76); m(7) * (t16 * t256 + t236 * t30 + t238 * t31) + m(6) * (t236 * t48 + t238 * t49 + t256 * t34); m(7) * (t17 * t256 + t236 * t33 + t238 * t32) + m(6) * (t236 * t51 + t238 * t50 + t256 * t35); (t235 * t236 + t237 * t238 + t255 * t256) * t351; (t236 ^ 2 + t238 ^ 2 + t256 ^ 2) * t351; t57 + m(7) * (t62 * t78 + t63 * t79) + t312 * t238 + t313 * t236; m(7) * (t16 * t60 + t30 * t79 + t31 * t78) + t6 * t346 + t279 * t11 / 0.2e1 + t5 * t347 + t13 * t345 + (t282 * t348 - t284 * t1 / 0.2e1) * t277; -t11 * t336 / 0.2e1 + t12 * t345 + t4 * t346 + t259 * t348 + t257 * t1 / 0.2e1 + m(7) * (t17 * t60 + t32 * t78 + t33 * t79) + t3 * t347; m(7) * (t235 * t79 + t237 * t78 + t255 * t60); m(7) * (t236 * t79 + t238 * t78 + t256 * t60); t238 * t2 + t236 * t1 + t256 * t11 + m(7) * (t60 ^ 2 + t78 ^ 2 + t79 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t100(1) t100(2) t100(4) t100(7) t100(11) t100(16); t100(2) t100(3) t100(5) t100(8) t100(12) t100(17); t100(4) t100(5) t100(6) t100(9) t100(13) t100(18); t100(7) t100(8) t100(9) t100(10) t100(14) t100(19); t100(11) t100(12) t100(13) t100(14) t100(15) t100(20); t100(16) t100(17) t100(18) t100(19) t100(20) t100(21);];
Mq  = res;
