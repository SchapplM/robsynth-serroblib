% Calculate joint inertia matrix for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:56:35
% EndTime: 2019-03-09 15:56:45
% DurationCPUTime: 4.56s
% Computational Cost: add. (9177->545), mult. (18372->767), div. (0->0), fcn. (21757->10), ass. (0->259)
t242 = sin(qJ(2));
t333 = Icges(3,5) * t242;
t332 = t333 / 0.2e1;
t241 = sin(qJ(3));
t244 = cos(qJ(3));
t245 = cos(qJ(2));
t179 = -Icges(4,3) * t245 + (Icges(4,5) * t244 - Icges(4,6) * t241) * t242;
t182 = -Icges(5,2) * t245 + (Icges(5,4) * t244 + Icges(5,6) * t241) * t242;
t331 = t179 + t182;
t243 = sin(qJ(1));
t246 = cos(qJ(1));
t296 = t245 * t246;
t208 = t243 * t241 + t244 * t296;
t239 = cos(pkin(10));
t226 = pkin(5) * t239 + pkin(4);
t240 = -pkin(9) - qJ(5);
t298 = t242 * t246;
t207 = t241 * t296 - t243 * t244;
t238 = sin(pkin(10));
t303 = t207 * t238;
t234 = pkin(10) + qJ(6);
t228 = sin(t234);
t229 = cos(t234);
t128 = t207 * t229 - t208 * t228;
t129 = t207 * t228 + t208 * t229;
t86 = t129 * rSges(7,1) + t128 * rSges(7,2) - rSges(7,3) * t298;
t330 = pkin(5) * t303 + t208 * t226 + t240 * t298 + t86;
t183 = -Icges(4,6) * t245 + (Icges(4,4) * t244 - Icges(4,2) * t241) * t242;
t301 = t241 * t242;
t192 = (-t238 * t244 + t239 * t241) * t242;
t302 = t238 * t241;
t193 = (t239 * t244 + t302) * t242;
t119 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t245;
t120 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t245;
t121 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t245;
t178 = -Icges(5,6) * t245 + (Icges(5,5) * t244 + Icges(5,3) * t241) * t242;
t186 = -Icges(5,4) * t245 + (Icges(5,1) * t244 + Icges(5,5) * t241) * t242;
t187 = -Icges(4,5) * t245 + (Icges(4,1) * t244 - Icges(4,4) * t241) * t242;
t299 = t242 * t244;
t324 = t245 * t119 + t192 * t120 + t193 * t121 + t178 * t301 + (t186 + t187) * t299;
t329 = (t183 * t301 + t331 * t245 - t324) * t245;
t297 = t243 * t245;
t205 = t241 * t297 + t244 * t246;
t206 = -t241 * t246 + t244 * t297;
t151 = t205 * t239 - t206 * t238;
t304 = t205 * t238;
t152 = t206 * t239 + t304;
t300 = t242 * t243;
t90 = Icges(6,5) * t152 + Icges(6,6) * t151 - Icges(6,3) * t300;
t92 = Icges(6,4) * t152 + Icges(6,2) * t151 - Icges(6,6) * t300;
t94 = Icges(6,1) * t152 + Icges(6,4) * t151 - Icges(6,5) * t300;
t23 = t151 * t92 + t152 * t94 - t300 * t90;
t153 = t207 * t239 - t208 * t238;
t154 = t208 * t239 + t303;
t91 = Icges(6,5) * t154 + Icges(6,6) * t153 - Icges(6,3) * t298;
t93 = Icges(6,4) * t154 + Icges(6,2) * t153 - Icges(6,6) * t298;
t95 = Icges(6,1) * t154 + Icges(6,4) * t153 - Icges(6,5) * t298;
t24 = t151 * t93 + t152 * t95 - t300 * t91;
t43 = -t119 * t300 + t120 * t151 + t121 * t152;
t130 = Icges(5,5) * t206 + Icges(5,6) * t300 + Icges(5,3) * t205;
t134 = Icges(5,4) * t206 + Icges(5,2) * t300 + Icges(5,6) * t205;
t138 = Icges(5,1) * t206 + Icges(5,4) * t300 + Icges(5,5) * t205;
t54 = t130 * t205 + t134 * t300 + t138 * t206;
t131 = Icges(5,5) * t208 + Icges(5,6) * t298 + Icges(5,3) * t207;
t135 = Icges(5,4) * t208 + Icges(5,2) * t298 + Icges(5,6) * t207;
t139 = Icges(5,1) * t208 + Icges(5,4) * t298 + Icges(5,5) * t207;
t55 = t131 * t205 + t135 * t300 + t139 * t206;
t132 = Icges(4,5) * t206 - Icges(4,6) * t205 + Icges(4,3) * t300;
t136 = Icges(4,4) * t206 - Icges(4,2) * t205 + Icges(4,6) * t300;
t140 = Icges(4,1) * t206 - Icges(4,4) * t205 + Icges(4,5) * t300;
t56 = t132 * t300 - t136 * t205 + t140 * t206;
t133 = Icges(4,5) * t208 - Icges(4,6) * t207 + Icges(4,3) * t298;
t137 = Icges(4,4) * t208 - Icges(4,2) * t207 + Icges(4,6) * t298;
t141 = Icges(4,1) * t208 - Icges(4,4) * t207 + Icges(4,5) * t298;
t57 = t133 * t300 - t137 * t205 + t141 * t206;
t75 = t178 * t205 + t182 * t300 + t186 * t206;
t76 = t179 * t300 - t183 * t205 + t187 * t206;
t328 = (-t75 - t43 - t76) * t245 + ((t24 + t55 + t57) * t246 + (t23 + t54 + t56) * t243) * t242;
t25 = t153 * t92 + t154 * t94 - t298 * t90;
t26 = t153 * t93 + t154 * t95 - t298 * t91;
t44 = -t119 * t298 + t153 * t120 + t154 * t121;
t58 = t207 * t130 + t134 * t298 + t208 * t138;
t59 = t207 * t131 + t135 * t298 + t208 * t139;
t60 = t132 * t298 - t207 * t136 + t208 * t140;
t61 = t133 * t298 - t207 * t137 + t208 * t141;
t77 = t207 * t178 + t182 * t298 + t208 * t186;
t78 = t179 * t298 - t207 * t183 + t208 * t187;
t327 = (-t77 - t78 - t44) * t245 + ((t26 + t59 + t61) * t246 + (t25 + t58 + t60) * t243) * t242;
t35 = t192 * t92 + t193 * t94 + t245 * t90;
t66 = -t134 * t245 + (t130 * t241 + t138 * t244) * t242;
t68 = -t132 * t245 + (-t136 * t241 + t140 * t244) * t242;
t326 = -t35 - t66 - t68;
t36 = t192 * t93 + t193 * t95 + t245 * t91;
t67 = -t135 * t245 + (t131 * t241 + t139 * t244) * t242;
t69 = -t133 * t245 + (-t137 * t241 + t141 * t244) * t242;
t325 = t36 + t67 + t69;
t126 = t205 * t229 - t206 * t228;
t127 = t205 * t228 + t206 * t229;
t79 = Icges(7,5) * t127 + Icges(7,6) * t126 - Icges(7,3) * t300;
t81 = Icges(7,4) * t127 + Icges(7,2) * t126 - Icges(7,6) * t300;
t83 = Icges(7,1) * t127 + Icges(7,4) * t126 - Icges(7,5) * t300;
t19 = t126 * t81 + t127 * t83 - t300 * t79;
t80 = Icges(7,5) * t129 + Icges(7,6) * t128 - Icges(7,3) * t298;
t82 = Icges(7,4) * t129 + Icges(7,2) * t128 - Icges(7,6) * t298;
t84 = Icges(7,1) * t129 + Icges(7,4) * t128 - Icges(7,5) * t298;
t20 = t126 * t82 + t127 * t84 - t300 * t80;
t173 = (-t228 * t244 + t229 * t241) * t242;
t174 = (t228 * t241 + t229 * t244) * t242;
t111 = Icges(7,5) * t174 + Icges(7,6) * t173 + Icges(7,3) * t245;
t112 = Icges(7,4) * t174 + Icges(7,2) * t173 + Icges(7,6) * t245;
t113 = Icges(7,1) * t174 + Icges(7,4) * t173 + Icges(7,5) * t245;
t39 = -t111 * t300 + t112 * t126 + t113 * t127;
t1 = -t39 * t245 + (t19 * t243 + t20 * t246) * t242;
t21 = t128 * t81 + t129 * t83 - t298 * t79;
t22 = t128 * t82 + t129 * t84 - t298 * t80;
t40 = -t111 * t298 + t128 * t112 + t129 * t113;
t2 = -t40 * t245 + (t21 * t243 + t22 * t246) * t242;
t28 = t173 * t81 + t174 * t83 + t245 * t79;
t29 = t173 * t82 + t174 * t84 + t245 * t80;
t278 = t245 * t111 + t173 * t112 + t174 * t113;
t45 = t278 * t245;
t5 = -t45 + (t28 * t243 + t29 * t246) * t242;
t323 = (t243 * t1 + t246 * t2) * t242 - t245 * t5;
t236 = t243 ^ 2;
t237 = t246 ^ 2;
t322 = 0.2e1 * t242;
t321 = m(6) / 0.2e1;
t320 = m(7) / 0.2e1;
t319 = -t243 / 0.2e1;
t318 = t243 / 0.2e1;
t317 = -t245 / 0.2e1;
t316 = t245 / 0.2e1;
t315 = -t246 / 0.2e1;
t314 = t246 / 0.2e1;
t313 = pkin(2) * t245;
t311 = -pkin(4) + t226;
t310 = rSges(5,3) * t205;
t309 = t246 * rSges(3,3);
t220 = qJ(5) * t300;
t279 = pkin(5) * t304;
t259 = -rSges(7,1) * t127 - rSges(7,2) * t126;
t85 = -rSges(7,3) * t300 - t259;
t308 = t206 * t311 + t240 * t300 + t220 + t279 + t85;
t202 = t208 * pkin(4);
t172 = -qJ(5) * t298 + t202;
t307 = -t172 + t330;
t305 = Icges(3,4) * t245;
t116 = rSges(7,1) * t174 + rSges(7,2) * t173 + rSges(7,3) * t245;
t295 = t116 + (-qJ(5) - t240) * t245 + (pkin(5) * t302 + t244 * t311) * t242;
t144 = t208 * rSges(5,1) + rSges(5,2) * t298 + t207 * rSges(5,3);
t161 = t208 * pkin(3) + t207 * qJ(4);
t294 = -t144 - t161;
t197 = t205 * qJ(4);
t160 = pkin(3) * t206 + t197;
t146 = t160 * t298;
t171 = pkin(4) * t206 - t220;
t293 = t171 * t298 + t146;
t292 = t154 * rSges(6,1) + t153 * rSges(6,2);
t209 = (pkin(3) * t244 + qJ(4) * t241) * t242;
t291 = t245 * t160 + t209 * t300;
t290 = -t161 - t172;
t190 = -rSges(5,2) * t245 + (rSges(5,1) * t244 + rSges(5,3) * t241) * t242;
t288 = -t190 - t209;
t191 = -rSges(4,3) * t245 + (rSges(4,1) * t244 - rSges(4,2) * t241) * t242;
t218 = pkin(2) * t242 - pkin(8) * t245;
t287 = -t191 - t218;
t284 = pkin(2) * t296 + pkin(8) * t298;
t286 = t236 * (pkin(8) * t242 + t313) + t246 * t284;
t210 = pkin(4) * t299 + qJ(5) * t245;
t285 = -t209 - t210;
t283 = t246 * pkin(1) + t243 * pkin(7);
t232 = t246 * pkin(7);
t282 = t232 - t197;
t281 = t236 + t237;
t280 = t321 + t320;
t277 = -t28 / 0.2e1 - t39 / 0.2e1;
t276 = -t29 / 0.2e1 - t40 / 0.2e1;
t97 = -rSges(6,3) * t298 + t292;
t275 = -t97 + t290;
t122 = rSges(6,1) * t193 + rSges(6,2) * t192 + rSges(6,3) * t245;
t273 = -t122 + t285;
t271 = -t218 + t288;
t145 = t208 * rSges(4,1) - t207 * rSges(4,2) + rSges(4,3) * t298;
t270 = -pkin(1) - t313;
t269 = t290 - t307;
t268 = t285 - t295;
t267 = -t218 + t273;
t266 = t243 * t160 + t246 * t161 + t286;
t265 = t245 * t171 + t210 * t300 + t291;
t264 = t283 + t284;
t262 = rSges(3,1) * t245 - rSges(3,2) * t242;
t261 = -rSges(4,1) * t206 + rSges(4,2) * t205;
t260 = -rSges(6,1) * t152 - rSges(6,2) * t151;
t258 = -t218 + t268;
t256 = -Icges(3,2) * t242 + t305;
t255 = Icges(3,5) * t245 - Icges(3,6) * t242;
t252 = rSges(3,1) * t296 - rSges(3,2) * t298 + t243 * rSges(3,3);
t251 = t243 * t171 + t246 * t172 + t266;
t250 = t161 + t264;
t249 = t1 * t314 + (t29 * t243 - t28 * t246) * t316 + t2 * t319;
t248 = t76 / 0.2e1 + t75 / 0.2e1 + t43 / 0.2e1 + t68 / 0.2e1 + t35 / 0.2e1 + t66 / 0.2e1 - t277;
t247 = t36 / 0.2e1 + t78 / 0.2e1 + t77 / 0.2e1 + t44 / 0.2e1 + t69 / 0.2e1 + t67 / 0.2e1 - t276;
t235 = t242 ^ 2;
t217 = rSges(2,1) * t246 - t243 * rSges(2,2);
t216 = -t243 * rSges(2,1) - rSges(2,2) * t246;
t215 = rSges(3,1) * t242 + rSges(3,2) * t245;
t212 = Icges(3,6) * t245 + t333;
t181 = Icges(3,3) * t243 + t246 * t255;
t180 = -Icges(3,3) * t246 + t243 * t255;
t167 = t252 + t283;
t166 = t309 + t232 + (-pkin(1) - t262) * t243;
t159 = t287 * t246;
t158 = t287 * t243;
t143 = rSges(4,3) * t300 - t261;
t142 = rSges(5,1) * t206 + rSges(5,2) * t300 + t310;
t123 = t246 * t252 + (t243 * t262 - t309) * t243;
t115 = t271 * t246;
t114 = t271 * t243;
t109 = t264 + t145;
t108 = t232 + ((-rSges(4,3) - pkin(8)) * t242 + t270) * t243 + t261;
t105 = -t245 * t145 - t191 * t298;
t104 = t143 * t245 + t191 * t300;
t96 = -rSges(6,3) * t300 - t260;
t89 = t267 * t246;
t88 = t267 * t243;
t87 = (t143 * t246 - t145 * t243) * t242;
t74 = t250 + t144;
t73 = -t310 + (-rSges(5,1) - pkin(3)) * t206 + ((-rSges(5,2) - pkin(8)) * t242 + t270) * t243 + t282;
t72 = t243 * t143 + t145 * t246 + t286;
t71 = t245 * t294 + t288 * t298;
t70 = t142 * t245 + t190 * t300 + t291;
t65 = t258 * t246;
t64 = t258 * t243;
t63 = t116 * t298 + t245 * t86;
t62 = -t116 * t300 - t245 * t85;
t53 = t202 + (-rSges(6,3) - qJ(5)) * t298 + t250 + t292;
t52 = t220 + (-pkin(3) - pkin(4)) * t206 + ((rSges(6,3) - pkin(8)) * t242 + t270) * t243 + t260 + t282;
t51 = t146 + (t142 * t246 + t243 * t294) * t242;
t49 = t243 * t142 + t144 * t246 + t266;
t48 = t250 + t330;
t47 = -t279 + (-pkin(3) - t226) * t206 + ((rSges(7,3) - pkin(8) - t240) * t242 + t270) * t243 + t259 + t282;
t46 = (t243 * t86 - t246 * t85) * t242;
t42 = t245 * t275 + t273 * t298;
t41 = t122 * t300 + t245 * t96 + t265;
t34 = (t243 * t275 + t246 * t96) * t242 + t293;
t33 = t61 * t243 - t246 * t60;
t32 = t59 * t243 - t246 * t58;
t31 = t57 * t243 - t246 * t56;
t30 = t55 * t243 - t246 * t54;
t27 = t243 * t96 + t246 * t97 + t251;
t18 = t245 * t269 + t268 * t298;
t17 = t245 * t308 + t295 * t300 + t265;
t12 = (t269 * t243 + t246 * t308) * t242 + t293;
t11 = t243 * t308 + t246 * t307 + t251;
t9 = t26 * t243 - t246 * t25;
t8 = -t23 * t246 + t24 * t243;
t7 = -t21 * t246 + t22 * t243;
t6 = -t19 * t246 + t20 * t243;
t3 = [Icges(2,3) + (Icges(3,1) * t242 - t183 * t241 + t305) * t242 + (Icges(3,4) * t242 + Icges(3,2) * t245 - t331) * t245 + m(7) * (t47 ^ 2 + t48 ^ 2) + m(6) * (t52 ^ 2 + t53 ^ 2) + m(5) * (t73 ^ 2 + t74 ^ 2) + m(4) * (t108 ^ 2 + t109 ^ 2) + m(3) * (t166 ^ 2 + t167 ^ 2) + m(2) * (t216 ^ 2 + t217 ^ 2) + t278 + t324; ((-Icges(3,6) * t246 + t243 * t256) * t317 + t246 * t332 + t212 * t314 - t248) * t246 + ((Icges(3,6) * t243 + t246 * t256) * t316 + t243 * t332 + t212 * t318 + t247) * t243 + m(7) * (t47 * t65 + t48 * t64) + m(6) * (t52 * t89 + t53 * t88) + m(5) * (t114 * t74 + t115 * t73) + m(4) * (t108 * t159 + t109 * t158) + m(3) * (-t166 * t246 - t167 * t243) * t215; m(6) * (t27 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(5) * (t114 ^ 2 + t115 ^ 2 + t49 ^ 2) + m(4) * (t158 ^ 2 + t159 ^ 2 + t72 ^ 2) + m(3) * (t215 ^ 2 * t281 + t123 ^ 2) + m(7) * (t11 ^ 2 + t64 ^ 2 + t65 ^ 2) + (-t237 * t180 - t30 - t31 - t6 - t8) * t246 + (t236 * t181 + t32 + t33 + t7 + t9 + (-t243 * t180 + t246 * t181) * t246) * t243; -t45 + t329 + m(7) * (t17 * t47 + t18 * t48) + m(6) * (t41 * t52 + t42 * t53) + m(5) * (t70 * t73 + t71 * t74) + m(4) * (t104 * t108 + t105 * t109) + (t243 * t248 + t246 * t247) * t242; m(7) * (t11 * t12 + t17 * t65 + t18 * t64) + m(6) * (t27 * t34 + t41 * t89 + t42 * t88) + m(5) * (t114 * t71 + t115 * t70 + t49 * t51) + m(4) * (t104 * t159 + t105 * t158 + t72 * t87) + ((t7 / 0.2e1 + t33 / 0.2e1 + t32 / 0.2e1 + t9 / 0.2e1) * t246 + (t6 / 0.2e1 + t31 / 0.2e1 + t30 / 0.2e1 + t8 / 0.2e1) * t243) * t242 - t249 + t327 * t318 + (t325 * t243 + t326 * t246) * t317 + t328 * t315; m(6) * (t34 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(5) * (t51 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(4) * (t104 ^ 2 + t105 ^ 2 + t87 ^ 2) + m(7) * (t12 ^ 2 + t17 ^ 2 + t18 ^ 2) + (-t5 - t329) * t245 + ((-t325 * t245 + t2 + t327) * t246 + (t326 * t245 + t1 + t328) * t243) * t242; m(7) * (t205 * t48 + t207 * t47) + m(6) * (t205 * t53 + t207 * t52) + m(5) * (t205 * t74 + t207 * t73); m(7) * (t11 * t301 + t205 * t64 + t207 * t65) + m(6) * (t205 * t88 + t207 * t89 + t27 * t301) + m(5) * (t114 * t205 + t115 * t207 + t301 * t49); m(7) * (t12 * t301 + t17 * t207 + t18 * t205) + m(6) * (t205 * t42 + t207 * t41 + t301 * t34) + m(5) * (t205 * t71 + t207 * t70 + t301 * t51); 0.2e1 * (m(5) / 0.2e1 + t280) * (t235 * t241 ^ 2 + t205 ^ 2 + t207 ^ 2); ((-t243 * t48 - t246 * t47) * t320 + (-t243 * t53 - t246 * t52) * t321) * t322; m(7) * (t245 * t11 + (-t243 * t64 - t246 * t65) * t242) + m(6) * (t245 * t27 + (-t243 * t88 - t246 * t89) * t242); m(7) * (t245 * t12 + (-t17 * t246 - t18 * t243) * t242) + m(6) * (t245 * t34 + (-t243 * t42 - t246 * t41) * t242); t280 * (-t205 * t243 - t207 * t246 + t241 * t245) * t322; 0.2e1 * t280 * (t235 * t281 + t245 ^ 2); m(7) * (t47 * t62 + t48 * t63) + t45 + (t243 * t277 + t246 * t276) * t242; m(7) * (t11 * t46 + t62 * t65 + t63 * t64) + (t315 * t7 + t319 * t6) * t242 + t249; m(7) * (t12 * t46 + t17 * t62 + t18 * t63) - t323; m(7) * (t205 * t63 + t207 * t62 + t301 * t46); m(7) * (t46 * t245 + (-t243 * t63 - t246 * t62) * t242); m(7) * (t46 ^ 2 + t62 ^ 2 + t63 ^ 2) + t323;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
