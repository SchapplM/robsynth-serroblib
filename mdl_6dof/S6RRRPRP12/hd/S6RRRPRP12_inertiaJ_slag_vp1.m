% Calculate joint inertia matrix for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP12_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP12_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:51:44
% EndTime: 2019-03-09 17:52:00
% DurationCPUTime: 6.07s
% Computational Cost: add. (17480->641), mult. (44807->878), div. (0->0), fcn. (56964->10), ass. (0->283)
t273 = cos(pkin(6));
t276 = sin(qJ(1));
t277 = cos(qJ(2));
t333 = t276 * t277;
t275 = sin(qJ(2));
t278 = cos(qJ(1));
t334 = t275 * t278;
t258 = t273 * t334 + t333;
t272 = sin(pkin(6));
t342 = cos(qJ(3));
t302 = t272 * t342;
t340 = sin(qJ(3));
t236 = t258 * t340 + t278 * t302;
t332 = t277 * t278;
t335 = t275 * t276;
t257 = -t273 * t332 + t335;
t274 = sin(qJ(5));
t341 = cos(qJ(5));
t193 = -t236 * t341 + t257 * t274;
t194 = t236 * t274 + t257 * t341;
t344 = rSges(7,3) + qJ(6);
t345 = rSges(7,1) + pkin(5);
t346 = -t193 * t344 - t194 * t345;
t301 = t272 * t340;
t255 = -t273 * t342 + t275 * t301;
t337 = t272 * t277;
t233 = t255 * t341 + t274 * t337;
t235 = t255 * t274 - t337 * t341;
t256 = t273 * t340 + t275 * t302;
t141 = Icges(7,5) * t235 + Icges(7,6) * t256 - Icges(7,3) * t233;
t143 = Icges(7,4) * t235 + Icges(7,2) * t256 - Icges(7,6) * t233;
t145 = Icges(7,1) * t235 + Icges(7,4) * t256 - Icges(7,5) * t233;
t71 = -t141 * t233 + t143 * t256 + t145 * t235;
t142 = Icges(6,5) * t235 + Icges(6,6) * t233 + Icges(6,3) * t256;
t144 = Icges(6,4) * t235 + Icges(6,2) * t233 + Icges(6,6) * t256;
t146 = Icges(6,1) * t235 + Icges(6,4) * t233 + Icges(6,5) * t256;
t72 = t142 * t256 + t144 * t233 + t146 * t235;
t343 = -t71 - t72;
t336 = t272 * t278;
t207 = Icges(3,5) * t258 - Icges(3,6) * t257 - Icges(3,3) * t336;
t339 = t207 * t278;
t338 = t272 * t276;
t237 = t258 * t342 - t278 * t301;
t331 = rSges(7,2) * t237 - t346;
t260 = -t273 * t335 + t332;
t238 = t260 * t340 - t276 * t302;
t259 = t273 * t333 + t334;
t195 = -t238 * t341 + t259 * t274;
t196 = t238 * t274 + t259 * t341;
t239 = t260 * t342 + t276 * t301;
t330 = t239 * rSges(7,2) + t195 * t344 + t196 * t345;
t329 = rSges(7,2) * t256 - t233 * t344 + t235 * t345;
t162 = rSges(5,1) * t259 - rSges(5,2) * t239 + rSges(5,3) * t238;
t176 = pkin(3) * t239 + qJ(4) * t238;
t328 = -t162 - t176;
t224 = t236 * qJ(4);
t175 = pkin(3) * t237 + t224;
t165 = t259 * t175;
t197 = pkin(4) * t257 + pkin(10) * t237;
t327 = t197 * t259 + t165;
t219 = pkin(3) * t256 + qJ(4) * t255;
t326 = t175 * t337 + t219 * t257;
t221 = pkin(2) * t260 + pkin(9) * t259;
t218 = t273 * t221;
t325 = t176 * t273 + t218;
t220 = pkin(2) * t258 + pkin(9) * t257;
t324 = -t175 - t220;
t198 = pkin(4) * t259 + pkin(10) * t239;
t323 = -t176 - t198;
t199 = -Icges(5,5) * t337 - Icges(5,6) * t256 + Icges(5,3) * t255;
t200 = -Icges(5,4) * t337 - Icges(5,2) * t256 + Icges(5,6) * t255;
t322 = t199 * t255 - t200 * t256;
t203 = Icges(4,4) * t256 - Icges(4,2) * t255 - Icges(4,6) * t337;
t204 = Icges(4,1) * t256 - Icges(4,4) * t255 - Icges(4,5) * t337;
t321 = -t203 * t255 + t204 * t256;
t205 = -rSges(5,1) * t337 - rSges(5,2) * t256 + rSges(5,3) * t255;
t320 = -t205 - t219;
t319 = t220 * t338 + t221 * t336;
t242 = -pkin(4) * t337 + pkin(10) * t256;
t318 = -t219 - t242;
t317 = pkin(1) * t278 + pkin(8) * t338;
t111 = Icges(7,5) * t194 + Icges(7,6) * t237 + Icges(7,3) * t193;
t115 = Icges(7,4) * t194 + Icges(7,2) * t237 + Icges(7,6) * t193;
t119 = Icges(7,1) * t194 + Icges(7,4) * t237 + Icges(7,5) * t193;
t40 = t111 * t193 + t115 * t237 + t119 * t194;
t112 = Icges(7,5) * t196 + Icges(7,6) * t239 + Icges(7,3) * t195;
t116 = Icges(7,4) * t196 + Icges(7,2) * t239 + Icges(7,6) * t195;
t120 = Icges(7,1) * t196 + Icges(7,4) * t239 + Icges(7,5) * t195;
t41 = t112 * t193 + t116 * t237 + t120 * t194;
t61 = t141 * t193 + t143 * t237 + t145 * t194;
t1 = t237 * t40 + t239 * t41 + t256 * t61;
t113 = Icges(6,5) * t194 - Icges(6,6) * t193 + Icges(6,3) * t237;
t117 = Icges(6,4) * t194 - Icges(6,2) * t193 + Icges(6,6) * t237;
t121 = Icges(6,1) * t194 - Icges(6,4) * t193 + Icges(6,5) * t237;
t42 = t113 * t237 - t117 * t193 + t121 * t194;
t114 = Icges(6,5) * t196 - Icges(6,6) * t195 + Icges(6,3) * t239;
t118 = Icges(6,4) * t196 - Icges(6,2) * t195 + Icges(6,6) * t239;
t122 = Icges(6,1) * t196 - Icges(6,4) * t195 + Icges(6,5) * t239;
t43 = t114 * t237 - t118 * t193 + t122 * t194;
t62 = t142 * t237 - t144 * t193 + t146 * t194;
t2 = t237 * t42 + t239 * t43 + t256 * t62;
t316 = -t2 / 0.2e1 - t1 / 0.2e1;
t44 = t111 * t195 + t115 * t239 + t119 * t196;
t45 = t112 * t195 + t116 * t239 + t120 * t196;
t63 = t141 * t195 + t143 * t239 + t145 * t196;
t3 = t237 * t44 + t239 * t45 + t256 * t63;
t46 = t113 * t239 - t117 * t195 + t121 * t196;
t47 = t114 * t239 - t118 * t195 + t122 * t196;
t64 = t142 * t239 - t144 * t195 + t146 * t196;
t4 = t237 * t46 + t239 * t47 + t256 * t64;
t315 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = t257 * t40 + t259 * t41 - t337 * t61;
t6 = t257 * t42 + t259 * t43 - t337 * t62;
t314 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t257 * t44 + t259 * t45 - t337 * t63;
t8 = t257 * t46 + t259 * t47 - t337 * t64;
t313 = t7 / 0.2e1 + t8 / 0.2e1;
t10 = t62 * t273 + (t276 * t43 - t278 * t42) * t272;
t9 = t61 * t273 + (t276 * t41 - t278 * t40) * t272;
t312 = t9 / 0.2e1 + t10 / 0.2e1;
t11 = t63 * t273 + (t276 * t45 - t278 * t44) * t272;
t12 = t64 * t273 + (t276 * t47 - t278 * t46) * t272;
t311 = t11 / 0.2e1 + t12 / 0.2e1;
t48 = -t111 * t233 + t115 * t256 + t119 * t235;
t49 = -t112 * t233 + t116 * t256 + t120 * t235;
t67 = t71 * t256;
t13 = t48 * t237 + t49 * t239 + t67;
t50 = t113 * t256 + t117 * t233 + t121 * t235;
t51 = t114 * t256 + t118 * t233 + t122 * t235;
t68 = t72 * t256;
t14 = t50 * t237 + t51 * t239 + t68;
t310 = t14 / 0.2e1 + t13 / 0.2e1;
t15 = t48 * t257 + t49 * t259 - t337 * t71;
t16 = t50 * t257 + t51 * t259 - t337 * t72;
t309 = t15 / 0.2e1 + t16 / 0.2e1;
t69 = t71 * t273;
t17 = t69 + (t49 * t276 - t48 * t278) * t272;
t70 = t72 * t273;
t18 = t70 + (t51 * t276 - t50 * t278) * t272;
t308 = t18 / 0.2e1 + t17 / 0.2e1;
t126 = rSges(6,1) * t196 - rSges(6,2) * t195 + rSges(6,3) * t239;
t307 = -t126 + t323;
t148 = rSges(6,1) * t235 + rSges(6,2) * t233 + rSges(6,3) * t256;
t306 = -t148 + t318;
t305 = t198 * t273 + t325;
t304 = -t197 + t324;
t164 = rSges(4,1) * t239 - rSges(4,2) * t238 + rSges(4,3) * t259;
t244 = Icges(3,3) * t273 + (Icges(3,5) * t275 + Icges(3,6) * t277) * t272;
t245 = Icges(3,6) * t273 + (Icges(3,4) * t275 + Icges(3,2) * t277) * t272;
t246 = Icges(3,5) * t273 + (Icges(3,1) * t275 + Icges(3,4) * t277) * t272;
t303 = t246 * t272 * t275 + t244 * t273 + t245 * t337;
t214 = rSges(3,1) * t260 - rSges(3,2) * t259 + rSges(3,3) * t338;
t300 = -t276 * pkin(1) + pkin(8) * t336;
t206 = rSges(4,1) * t256 - rSges(4,2) * t255 - rSges(4,3) * t337;
t261 = (pkin(2) * t275 - pkin(9) * t277) * t272;
t299 = t272 * (-t206 - t261);
t298 = t323 - t330;
t297 = t318 - t329;
t296 = t175 * t338 + t176 * t336 + t319;
t295 = t197 * t337 + t242 * t257 + t326;
t294 = t272 * (-t261 + t320);
t293 = -rSges(5,1) * t257 - rSges(5,3) * t236;
t292 = -rSges(6,1) * t194 + rSges(6,2) * t193;
t291 = t221 + t317;
t290 = t272 * (-t261 + t306);
t289 = t61 / 0.2e1 + t50 / 0.2e1 + t48 / 0.2e1 + t62 / 0.2e1;
t288 = t63 / 0.2e1 + t51 / 0.2e1 + t49 / 0.2e1 + t64 / 0.2e1;
t287 = t197 * t338 + t198 * t336 + t296;
t286 = t272 * (-t261 + t297);
t285 = -t220 + t300;
t163 = rSges(4,1) * t237 - rSges(4,2) * t236 + rSges(4,3) * t257;
t284 = -t224 + t285;
t213 = rSges(3,1) * t258 - rSges(3,2) * t257 - rSges(3,3) * t336;
t283 = t176 + t291;
t282 = -t197 + t284;
t149 = Icges(5,5) * t257 - Icges(5,6) * t237 + Icges(5,3) * t236;
t153 = Icges(5,4) * t257 - Icges(5,2) * t237 + Icges(5,6) * t236;
t157 = Icges(5,1) * t257 - Icges(5,4) * t237 + Icges(5,5) * t236;
t85 = t149 * t255 - t153 * t256 - t157 * t337;
t151 = Icges(4,5) * t237 - Icges(4,6) * t236 + Icges(4,3) * t257;
t155 = Icges(4,4) * t237 - Icges(4,2) * t236 + Icges(4,6) * t257;
t159 = Icges(4,1) * t237 - Icges(4,4) * t236 + Icges(4,5) * t257;
t87 = -t151 * t337 - t155 * t255 + t159 * t256;
t201 = -Icges(5,1) * t337 - Icges(5,4) * t256 + Icges(5,5) * t255;
t94 = t199 * t236 - t200 * t237 + t201 * t257;
t202 = Icges(4,5) * t256 - Icges(4,6) * t255 - Icges(4,3) * t337;
t96 = t202 * t257 - t203 * t236 + t204 * t237;
t281 = t87 / 0.2e1 + t85 / 0.2e1 + t96 / 0.2e1 + t94 / 0.2e1 + t289;
t150 = Icges(5,5) * t259 - Icges(5,6) * t239 + Icges(5,3) * t238;
t154 = Icges(5,4) * t259 - Icges(5,2) * t239 + Icges(5,6) * t238;
t158 = Icges(5,1) * t259 - Icges(5,4) * t239 + Icges(5,5) * t238;
t86 = t150 * t255 - t154 * t256 - t158 * t337;
t152 = Icges(4,5) * t239 - Icges(4,6) * t238 + Icges(4,3) * t259;
t156 = Icges(4,4) * t239 - Icges(4,2) * t238 + Icges(4,6) * t259;
t160 = Icges(4,1) * t239 - Icges(4,4) * t238 + Icges(4,5) * t259;
t88 = -t152 * t337 - t156 * t255 + t160 * t256;
t95 = t199 * t238 - t200 * t239 + t201 * t259;
t97 = t202 * t259 - t203 * t238 + t204 * t239;
t280 = t86 / 0.2e1 + t88 / 0.2e1 + t95 / 0.2e1 + t97 / 0.2e1 + t288;
t279 = t198 + t283;
t263 = rSges(2,1) * t278 - rSges(2,2) * t276;
t262 = -rSges(2,1) * t276 - rSges(2,2) * t278;
t247 = rSges(3,3) * t273 + (rSges(3,1) * t275 + rSges(3,2) * t277) * t272;
t212 = Icges(3,1) * t260 - Icges(3,4) * t259 + Icges(3,5) * t338;
t211 = Icges(3,1) * t258 - Icges(3,4) * t257 - Icges(3,5) * t336;
t210 = Icges(3,4) * t260 - Icges(3,2) * t259 + Icges(3,6) * t338;
t209 = Icges(3,4) * t258 - Icges(3,2) * t257 - Icges(3,6) * t336;
t208 = Icges(3,5) * t260 - Icges(3,6) * t259 + Icges(3,3) * t338;
t185 = t214 + t317;
t184 = -t213 + t300;
t168 = -t213 * t273 - t247 * t336;
t167 = t214 * t273 - t247 * t338;
t161 = -rSges(5,2) * t237 - t293;
t140 = t303 * t273;
t137 = (t213 * t276 + t214 * t278) * t272;
t136 = t244 * t338 - t245 * t259 + t246 * t260;
t135 = -t244 * t336 - t245 * t257 + t246 * t258;
t128 = t291 + t164;
t127 = -t163 + t285;
t124 = rSges(6,3) * t237 - t292;
t110 = -t164 * t337 - t206 * t259;
t109 = t163 * t337 + t206 * t257;
t108 = t208 * t273 + (t210 * t277 + t212 * t275) * t272;
t107 = t207 * t273 + (t209 * t277 + t211 * t275) * t272;
t106 = -t202 * t337 + t321;
t105 = -t201 * t337 + t322;
t104 = t106 * t273;
t103 = t105 * t273;
t102 = t283 + t162;
t101 = (rSges(5,2) - pkin(3)) * t237 + t284 + t293;
t100 = t163 * t259 - t164 * t257;
t99 = (-t163 - t220) * t273 + t278 * t299;
t98 = t164 * t273 + t276 * t299 + t218;
t93 = (t163 * t276 + t164 * t278) * t272 + t319;
t92 = t126 * t256 - t148 * t239;
t91 = -t124 * t256 + t148 * t237;
t90 = t259 * t320 + t328 * t337;
t89 = t161 * t337 + t205 * t257 + t326;
t84 = t279 + t126;
t83 = (-rSges(6,3) - pkin(3)) * t237 + t282 + t292;
t82 = (-t161 + t324) * t273 + t278 * t294;
t81 = t162 * t273 + t276 * t294 + t325;
t80 = t152 * t259 - t156 * t238 + t160 * t239;
t79 = t151 * t259 - t155 * t238 + t159 * t239;
t78 = t152 * t257 - t156 * t236 + t160 * t237;
t77 = t151 * t257 - t155 * t236 + t159 * t237;
t76 = t150 * t238 - t154 * t239 + t158 * t259;
t75 = t149 * t238 - t153 * t239 + t157 * t259;
t74 = t150 * t236 - t154 * t237 + t158 * t257;
t73 = t149 * t236 - t153 * t237 + t157 * t257;
t66 = t124 * t239 - t126 * t237;
t65 = t161 * t259 + t257 * t328 + t165;
t60 = (t161 * t276 + t162 * t278) * t272 + t296;
t59 = t279 + t330;
t58 = (-rSges(7,2) - pkin(3)) * t237 + t282 + t346;
t57 = t259 * t306 + t307 * t337;
t56 = t124 * t337 + t148 * t257 + t295;
t55 = (-t124 + t304) * t273 + t278 * t290;
t54 = t126 * t273 + t276 * t290 + t305;
t53 = -t239 * t329 + t256 * t330;
t52 = t237 * t329 - t256 * t331;
t39 = t124 * t259 + t257 * t307 + t327;
t38 = -t237 * t330 + t239 * t331;
t37 = (t124 * t276 + t126 * t278) * t272 + t287;
t36 = t259 * t297 + t298 * t337;
t35 = t257 * t329 + t331 * t337 + t295;
t34 = (t304 - t331) * t273 + t278 * t286;
t33 = t273 * t330 + t276 * t286 + t305;
t32 = t104 + (t88 * t276 - t87 * t278) * t272;
t31 = t103 + (t86 * t276 - t85 * t278) * t272;
t30 = -t106 * t337 + t87 * t257 + t88 * t259;
t29 = -t105 * t337 + t85 * t257 + t86 * t259;
t28 = t257 * t298 + t259 * t331 + t327;
t27 = (t276 * t331 + t278 * t330) * t272 + t287;
t26 = t97 * t273 + (t276 * t80 - t278 * t79) * t272;
t25 = t96 * t273 + (t276 * t78 - t278 * t77) * t272;
t24 = t95 * t273 + (t276 * t76 - t278 * t75) * t272;
t23 = t94 * t273 + (t276 * t74 - t278 * t73) * t272;
t22 = t257 * t79 + t259 * t80 - t337 * t97;
t21 = t257 * t77 + t259 * t78 - t337 * t96;
t20 = t257 * t75 + t259 * t76 - t337 * t95;
t19 = t257 * t73 + t259 * t74 - t337 * t94;
t123 = [(-t201 - t202) * t337 + t303 + Icges(2,3) + m(7) * (t58 ^ 2 + t59 ^ 2) + m(6) * (t83 ^ 2 + t84 ^ 2) + m(5) * (t101 ^ 2 + t102 ^ 2) + m(4) * (t127 ^ 2 + t128 ^ 2) + m(3) * (t184 ^ 2 + t185 ^ 2) + m(2) * (t262 ^ 2 + t263 ^ 2) + t321 + t322 - t343; t70 + t69 + t104 + t103 + t140 + m(7) * (t33 * t59 + t34 * t58) + m(6) * (t54 * t84 + t55 * t83) + m(5) * (t101 * t82 + t102 * t81) + m(4) * (t127 * t99 + t128 * t98) + m(3) * (t167 * t185 + t168 * t184) + ((-t107 / 0.2e1 - t135 / 0.2e1 - t281) * t278 + (t108 / 0.2e1 + t136 / 0.2e1 + t280) * t276) * t272; (t18 + t17 + t31 + t32 + t140) * t273 + m(7) * (t27 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(6) * (t37 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t60 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(4) * (t93 ^ 2 + t98 ^ 2 + t99 ^ 2) + m(3) * (t137 ^ 2 + t167 ^ 2 + t168 ^ 2) + ((-t10 - t9 - t23 - t25 + (-t209 * t257 + t211 * t258 - t272 * t339) * t336) * t278 + (t12 + t11 + t24 + t26 + ((-t210 * t259 + t212 * t260 + (t208 * t276 - t339) * t272) * t276 + (t208 * t336 + t209 * t259 + t210 * t257 - t211 * t260 - t212 * t258) * t278) * t272) * t276 + ((-t107 - t135) * t278 + (t108 + t136) * t276) * t273) * t272; (-t105 - t106 + t343) * t337 + m(7) * (t35 * t58 + t36 * t59) + m(6) * (t56 * t83 + t57 * t84) + m(5) * (t101 * t89 + t102 * t90) + m(4) * (t109 * t127 + t110 * t128) + t280 * t259 + t281 * t257; (t29 / 0.2e1 + t30 / 0.2e1 + t309) * t273 + (t24 / 0.2e1 + t26 / 0.2e1 + t311) * t259 + (t23 / 0.2e1 + t25 / 0.2e1 + t312) * t257 + m(6) * (t37 * t39 + t54 * t57 + t55 * t56) + m(7) * (t27 * t28 + t33 * t36 + t34 * t35) + m(5) * (t60 * t65 + t81 * t90 + t82 * t89) + m(4) * (t100 * t93 + t109 * t99 + t110 * t98) + ((-t19 / 0.2e1 - t21 / 0.2e1 - t314) * t278 + (-t31 / 0.2e1 - t32 / 0.2e1 - t308) * t277 + (t20 / 0.2e1 + t22 / 0.2e1 + t313) * t276) * t272; (-t15 - t16 - t29 - t30) * t337 + (t7 + t8 + t20 + t22) * t259 + (t6 + t5 + t19 + t21) * t257 + m(7) * (t28 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(6) * (t39 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t65 ^ 2 + t89 ^ 2 + t90 ^ 2) + m(4) * (t100 ^ 2 + t109 ^ 2 + t110 ^ 2); m(7) * (t236 * t59 + t238 * t58) + m(6) * (t236 * t84 + t238 * t83) + m(5) * (t101 * t238 + t102 * t236); m(7) * (t236 * t33 + t238 * t34 + t255 * t27) + m(6) * (t236 * t54 + t238 * t55 + t255 * t37) + m(5) * (t236 * t81 + t238 * t82 + t255 * t60); m(7) * (t236 * t36 + t238 * t35 + t255 * t28) + m(6) * (t236 * t57 + t238 * t56 + t255 * t39) + m(5) * (t236 * t90 + t238 * t89 + t255 * t65); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t236 ^ 2 + t238 ^ 2 + t255 ^ 2); t68 + t67 + m(7) * (t52 * t58 + t53 * t59) + m(6) * (t83 * t91 + t84 * t92) + t288 * t239 + t289 * t237; t310 * t273 + t308 * t256 + t311 * t239 + t312 * t237 + m(7) * (t27 * t38 + t33 * t53 + t34 * t52) + m(6) * (t37 * t66 + t54 * t92 + t55 * t91) + (t276 * t315 + t278 * t316) * t272; -t310 * t337 + t315 * t259 - t316 * t257 + t309 * t256 + t313 * t239 + t314 * t237 + m(7) * (t28 * t38 + t35 * t52 + t36 * t53) + m(6) * (t39 * t66 + t56 * t91 + t57 * t92); m(6) * (t236 * t92 + t238 * t91 + t255 * t66) + m(7) * (t236 * t53 + t238 * t52 + t255 * t38); (t14 + t13) * t256 + (t3 + t4) * t239 + (t2 + t1) * t237 + m(7) * (t38 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(6) * (t66 ^ 2 + t91 ^ 2 + t92 ^ 2); m(7) * (t193 * t59 + t195 * t58); m(7) * (t193 * t33 + t195 * t34 - t233 * t27); m(7) * (t193 * t36 + t195 * t35 - t233 * t28); m(7) * (t193 * t236 + t195 * t238 - t233 * t255); m(7) * (t193 * t53 + t195 * t52 - t233 * t38); m(7) * (t193 ^ 2 + t195 ^ 2 + t233 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t123(1) t123(2) t123(4) t123(7) t123(11) t123(16); t123(2) t123(3) t123(5) t123(8) t123(12) t123(17); t123(4) t123(5) t123(6) t123(9) t123(13) t123(18); t123(7) t123(8) t123(9) t123(10) t123(14) t123(19); t123(11) t123(12) t123(13) t123(14) t123(15) t123(20); t123(16) t123(17) t123(18) t123(19) t123(20) t123(21);];
Mq  = res;
