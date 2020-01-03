% Calculate joint inertia matrix for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR11_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR11_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:38:50
% EndTime: 2019-12-31 22:39:00
% DurationCPUTime: 3.56s
% Computational Cost: add. (20128->526), mult. (43534->745), div. (0->0), fcn. (55725->12), ass. (0->257)
t253 = cos(pkin(5));
t257 = sin(qJ(1));
t259 = cos(qJ(2));
t300 = t257 * t259;
t256 = sin(qJ(2));
t260 = cos(qJ(1));
t302 = t256 * t260;
t235 = t253 * t302 + t300;
t255 = sin(qJ(3));
t252 = sin(pkin(5));
t311 = cos(qJ(3));
t276 = t252 * t311;
t215 = t235 * t255 + t260 * t276;
t317 = t215 / 0.2e1;
t298 = t259 * t260;
t303 = t256 * t257;
t237 = -t253 * t303 + t298;
t217 = t237 * t255 - t257 * t276;
t316 = t217 / 0.2e1;
t305 = t252 * t256;
t232 = -t253 * t311 + t255 * t305;
t315 = t232 / 0.2e1;
t234 = -t253 * t298 + t303;
t314 = t234 / 0.2e1;
t236 = t253 * t300 + t302;
t313 = t236 / 0.2e1;
t312 = t253 / 0.2e1;
t258 = cos(qJ(4));
t247 = pkin(4) * t258 + pkin(3);
t310 = -pkin(3) + t247;
t297 = t260 * t252;
t216 = t235 * t311 - t255 * t297;
t251 = qJ(4) + qJ(5);
t248 = sin(t251);
t249 = cos(t251);
t172 = -t216 * t248 + t234 * t249;
t173 = t216 * t249 + t234 * t248;
t111 = Icges(6,5) * t173 + Icges(6,6) * t172 + Icges(6,3) * t215;
t113 = Icges(6,4) * t173 + Icges(6,2) * t172 + Icges(6,6) * t215;
t115 = Icges(6,1) * t173 + Icges(6,4) * t172 + Icges(6,5) * t215;
t233 = t253 * t255 + t256 * t276;
t304 = t252 * t259;
t204 = -t233 * t248 - t249 * t304;
t205 = t233 * t249 - t248 * t304;
t59 = t111 * t232 + t113 * t204 + t115 * t205;
t309 = t59 * t215;
t301 = t257 * t252;
t218 = t237 * t311 + t255 * t301;
t174 = -t218 * t248 + t236 * t249;
t175 = t218 * t249 + t236 * t248;
t112 = Icges(6,5) * t175 + Icges(6,6) * t174 + Icges(6,3) * t217;
t114 = Icges(6,4) * t175 + Icges(6,2) * t174 + Icges(6,6) * t217;
t116 = Icges(6,1) * t175 + Icges(6,4) * t174 + Icges(6,5) * t217;
t60 = t112 * t232 + t114 * t204 + t116 * t205;
t308 = t60 * t217;
t254 = sin(qJ(4));
t307 = t234 * t254;
t306 = t236 * t254;
t299 = t257 * t260;
t211 = t215 * pkin(9);
t261 = -pkin(10) - pkin(9);
t285 = pkin(4) * t307;
t107 = -t215 * t261 + t310 * t216 - t211 + t285;
t267 = -t173 * rSges(6,1) - t172 * rSges(6,2);
t119 = t215 * rSges(6,3) - t267;
t296 = t107 + t119;
t167 = t218 * pkin(3) + pkin(9) * t217;
t278 = pkin(4) * t306 - t217 * t261 + t218 * t247;
t108 = -t167 + t278;
t120 = t175 * rSges(6,1) + t174 * rSges(6,2) + t217 * rSges(6,3);
t295 = t108 + t120;
t183 = -t218 * t254 + t236 * t258;
t184 = t218 * t258 + t306;
t128 = t184 * rSges(5,1) + t183 * rSges(5,2) + t217 * rSges(5,3);
t294 = -t128 - t167;
t144 = rSges(6,1) * t205 + rSges(6,2) * t204 + rSges(6,3) * t232;
t281 = t254 * t304;
t145 = -pkin(4) * t281 + t310 * t233 + (-pkin(9) - t261) * t232;
t293 = t144 + t145;
t213 = -t233 * t254 - t258 * t304;
t214 = t233 * t258 - t281;
t150 = rSges(5,1) * t214 + rSges(5,2) * t213 + rSges(5,3) * t232;
t200 = t233 * pkin(3) + t232 * pkin(9);
t292 = -t150 - t200;
t166 = t216 * pkin(3) + t211;
t291 = t166 * t304 + t234 * t200;
t202 = t237 * pkin(2) + pkin(8) * t236;
t199 = t253 * t202;
t290 = t253 * t167 + t199;
t201 = t235 * pkin(2) + t234 * pkin(8);
t289 = -t166 - t201;
t288 = t201 * t301 + t202 * t297;
t287 = t260 * pkin(1) + pkin(7) * t301;
t141 = Icges(6,5) * t205 + Icges(6,6) * t204 + Icges(6,3) * t232;
t142 = Icges(6,4) * t205 + Icges(6,2) * t204 + Icges(6,6) * t232;
t143 = Icges(6,1) * t205 + Icges(6,4) * t204 + Icges(6,5) * t232;
t78 = t232 * t141 + t204 * t142 + t205 * t143;
t73 = t78 * t232;
t26 = t308 + t73 + t309;
t47 = t111 * t215 + t113 * t172 + t115 * t173;
t48 = t112 * t215 + t114 * t172 + t116 * t173;
t67 = t141 * t215 + t142 * t172 + t143 * t173;
t7 = t215 * t47 + t217 * t48 + t232 * t67;
t49 = t111 * t217 + t113 * t174 + t115 * t175;
t50 = t112 * t217 + t114 * t174 + t116 * t175;
t68 = t141 * t217 + t142 * t174 + t143 * t175;
t8 = t215 * t49 + t217 * t50 + t232 * t68;
t286 = t215 * t7 + t217 * t8 + t232 * t26;
t185 = Icges(4,5) * t233 - Icges(4,6) * t232 - Icges(4,3) * t304;
t186 = Icges(4,4) * t233 - Icges(4,2) * t232 - Icges(4,6) * t304;
t187 = Icges(4,1) * t233 - Icges(4,4) * t232 - Icges(4,5) * t304;
t105 = -t185 * t304 - t232 * t186 + t233 * t187;
t147 = Icges(5,5) * t214 + Icges(5,6) * t213 + Icges(5,3) * t232;
t148 = Icges(5,4) * t214 + Icges(5,2) * t213 + Icges(5,6) * t232;
t149 = Icges(5,1) * t214 + Icges(5,4) * t213 + Icges(5,5) * t232;
t82 = t232 * t147 + t213 * t148 + t214 * t149;
t284 = -t105 - t78 - t82;
t122 = Icges(5,5) * t184 + Icges(5,6) * t183 + Icges(5,3) * t217;
t124 = Icges(5,4) * t184 + Icges(5,2) * t183 + Icges(5,6) * t217;
t126 = Icges(5,1) * t184 + Icges(5,4) * t183 + Icges(5,5) * t217;
t63 = t122 * t232 + t124 * t213 + t126 * t214;
t72 = t147 * t217 + t148 * t183 + t149 * t184;
t283 = t63 / 0.2e1 + t72 / 0.2e1;
t181 = -t216 * t254 + t234 * t258;
t182 = t216 * t258 + t307;
t121 = Icges(5,5) * t182 + Icges(5,6) * t181 + Icges(5,3) * t215;
t123 = Icges(5,4) * t182 + Icges(5,2) * t181 + Icges(5,6) * t215;
t125 = Icges(5,1) * t182 + Icges(5,4) * t181 + Icges(5,5) * t215;
t62 = t121 * t232 + t123 * t213 + t125 * t214;
t71 = t147 * t215 + t148 * t181 + t149 * t182;
t282 = t71 / 0.2e1 + t62 / 0.2e1;
t280 = -t167 - t295;
t279 = -t200 - t293;
t158 = t218 * rSges(4,1) - t217 * rSges(4,2) + t236 * rSges(4,3);
t222 = Icges(3,3) * t253 + (Icges(3,5) * t256 + Icges(3,6) * t259) * t252;
t223 = Icges(3,6) * t253 + (Icges(3,4) * t256 + Icges(3,2) * t259) * t252;
t224 = Icges(3,5) * t253 + (Icges(3,1) * t256 + Icges(3,4) * t259) * t252;
t277 = t253 * t222 + t223 * t304 + t224 * t305;
t196 = t237 * rSges(3,1) - t236 * rSges(3,2) + rSges(3,3) * t301;
t275 = -t304 / 0.2e1;
t274 = -t257 * pkin(1) + pkin(7) * t297;
t188 = rSges(4,1) * t233 - rSges(4,2) * t232 - rSges(4,3) * t304;
t238 = (pkin(2) * t256 - pkin(8) * t259) * t252;
t273 = t252 * (-t188 - t238);
t272 = t166 * t301 + t167 * t297 + t288;
t271 = t309 / 0.2e1 + t308 / 0.2e1 + t67 * t317 + t68 * t316 + t73;
t270 = t252 * (-t238 + t292);
t13 = t234 * t47 + t236 * t48 - t67 * t304;
t14 = t234 * t49 + t236 * t50 - t68 * t304;
t28 = t59 * t234 + t60 * t236 - t78 * t304;
t269 = t13 * t317 + t14 * t316 + t26 * t275 + t28 * t315 + t8 * t313 + t7 * t314;
t15 = t253 * t67 + (t257 * t48 - t260 * t47) * t252;
t16 = t253 * t68 + (t257 * t50 - t260 * t49) * t252;
t76 = t78 * t253;
t30 = t76 + (t60 * t257 - t59 * t260) * t252;
t268 = t15 * t317 + t16 * t316 + t26 * t312 + t30 * t315 + t8 * t301 / 0.2e1 - t7 * t297 / 0.2e1;
t266 = t202 + t287;
t265 = t252 * (-t238 + t279);
t264 = -t201 + t274;
t157 = rSges(4,1) * t216 - rSges(4,2) * t215 + rSges(4,3) * t234;
t127 = rSges(5,1) * t182 + rSges(5,2) * t181 + rSges(5,3) * t215;
t195 = rSges(3,1) * t235 - rSges(3,2) * t234 - rSges(3,3) * t297;
t151 = Icges(4,5) * t216 - Icges(4,6) * t215 + Icges(4,3) * t234;
t153 = Icges(4,4) * t216 - Icges(4,2) * t215 + Icges(4,6) * t234;
t155 = Icges(4,1) * t216 - Icges(4,4) * t215 + Icges(4,5) * t234;
t89 = -t151 * t304 - t153 * t232 + t155 * t233;
t98 = t185 * t234 - t186 * t215 + t187 * t216;
t263 = t67 / 0.2e1 + t59 / 0.2e1 + t98 / 0.2e1 + t89 / 0.2e1 + t282;
t152 = Icges(4,5) * t218 - Icges(4,6) * t217 + Icges(4,3) * t236;
t154 = Icges(4,4) * t218 - Icges(4,2) * t217 + Icges(4,6) * t236;
t156 = Icges(4,1) * t218 - Icges(4,4) * t217 + Icges(4,5) * t236;
t90 = -t152 * t304 - t154 * t232 + t156 * t233;
t99 = t185 * t236 - t186 * t217 + t187 * t218;
t262 = t68 / 0.2e1 + t60 / 0.2e1 + t99 / 0.2e1 + t90 / 0.2e1 + t283;
t240 = rSges(2,1) * t260 - rSges(2,2) * t257;
t239 = -rSges(2,1) * t257 - rSges(2,2) * t260;
t225 = rSges(3,3) * t253 + (rSges(3,1) * t256 + rSges(3,2) * t259) * t252;
t194 = Icges(3,1) * t237 - Icges(3,4) * t236 + Icges(3,5) * t301;
t193 = Icges(3,1) * t235 - Icges(3,4) * t234 - Icges(3,5) * t297;
t192 = Icges(3,4) * t237 - Icges(3,2) * t236 + Icges(3,6) * t301;
t191 = Icges(3,4) * t235 - Icges(3,2) * t234 - Icges(3,6) * t297;
t190 = Icges(3,5) * t237 - Icges(3,6) * t236 + Icges(3,3) * t301;
t189 = Icges(3,5) * t235 - Icges(3,6) * t234 - Icges(3,3) * t297;
t177 = t196 + t287;
t176 = -t195 + t274;
t161 = -t195 * t253 - t225 * t297;
t160 = t196 * t253 - t225 * t301;
t159 = t236 * t166;
t146 = t277 * t253;
t139 = (t195 * t257 + t196 * t260) * t252;
t138 = t222 * t301 - t223 * t236 + t224 * t237;
t137 = -t222 * t297 - t223 * t234 + t224 * t235;
t133 = t215 * t144;
t130 = t266 + t158;
t129 = -t157 + t264;
t118 = -t158 * t304 - t188 * t236;
t117 = t157 * t304 + t188 * t234;
t110 = t190 * t253 + (t192 * t259 + t194 * t256) * t252;
t109 = t189 * t253 + (t191 * t259 + t193 * t256) * t252;
t106 = t232 * t120;
t104 = t217 * t119;
t103 = t105 * t253;
t102 = t157 * t236 - t158 * t234;
t101 = (-t157 - t201) * t253 + t260 * t273;
t100 = t158 * t253 + t257 * t273 + t199;
t97 = t266 - t294;
t96 = -t127 - t166 + t264;
t95 = (t157 * t257 + t158 * t260) * t252 + t288;
t94 = t128 * t232 - t150 * t217;
t93 = -t127 * t232 + t150 * t215;
t92 = t266 + t278 + t120;
t91 = -t285 - t216 * t247 + (-rSges(6,3) + t261) * t215 + t264 + t267;
t88 = -t144 * t217 + t106;
t87 = -t119 * t232 + t133;
t86 = t152 * t236 - t154 * t217 + t156 * t218;
t85 = t151 * t236 - t153 * t217 + t155 * t218;
t84 = t152 * t234 - t154 * t215 + t156 * t216;
t83 = t151 * t234 - t153 * t215 + t155 * t216;
t81 = t82 * t253;
t80 = t82 * t232;
t79 = t127 * t217 - t128 * t215;
t77 = -t120 * t215 + t104;
t75 = t292 * t236 + t294 * t304;
t74 = t127 * t304 + t150 * t234 + t291;
t70 = (-t127 + t289) * t253 + t260 * t270;
t69 = t128 * t253 + t257 * t270 + t290;
t64 = t127 * t236 + t294 * t234 + t159;
t61 = (t127 * t257 + t128 * t260) * t252 + t272;
t58 = t122 * t217 + t124 * t183 + t126 * t184;
t57 = t121 * t217 + t123 * t183 + t125 * t184;
t56 = t122 * t215 + t124 * t181 + t126 * t182;
t55 = t121 * t215 + t123 * t181 + t125 * t182;
t52 = t108 * t232 - t293 * t217 + t106;
t51 = t145 * t215 - t296 * t232 + t133;
t46 = t279 * t236 + t280 * t304;
t45 = t293 * t234 + t296 * t304 + t291;
t44 = (t289 - t296) * t253 + t260 * t265;
t43 = t295 * t253 + t257 * t265 + t290;
t42 = t107 * t217 - t295 * t215 + t104;
t41 = t103 + (t90 * t257 - t89 * t260) * t252;
t40 = -t105 * t304 + t89 * t234 + t90 * t236;
t39 = t280 * t234 + t296 * t236 + t159;
t38 = (t296 * t257 + t295 * t260) * t252 + t272;
t37 = t253 * t99 + (t257 * t86 - t260 * t85) * t252;
t36 = t253 * t98 + (t257 * t84 - t260 * t83) * t252;
t35 = t234 * t85 + t236 * t86 - t99 * t304;
t34 = t234 * t83 + t236 * t84 - t98 * t304;
t33 = t81 + (t63 * t257 - t62 * t260) * t252;
t32 = t62 * t234 + t63 * t236 - t82 * t304;
t31 = t62 * t215 + t63 * t217 + t80;
t22 = t253 * t72 + (t257 * t58 - t260 * t57) * t252;
t21 = t253 * t71 + (t257 * t56 - t260 * t55) * t252;
t20 = t234 * t57 + t236 * t58 - t72 * t304;
t19 = t234 * t55 + t236 * t56 - t71 * t304;
t18 = t215 * t57 + t217 * t58 + t232 * t72;
t17 = t215 * t55 + t217 * t56 + t232 * t71;
t1 = [Icges(2,3) + m(6) * (t91 ^ 2 + t92 ^ 2) + m(5) * (t96 ^ 2 + t97 ^ 2) + m(4) * (t129 ^ 2 + t130 ^ 2) + m(3) * (t176 ^ 2 + t177 ^ 2) + m(2) * (t239 ^ 2 + t240 ^ 2) + t277 - t284; t76 + t81 + t103 + t146 + m(6) * (t43 * t92 + t44 * t91) + m(5) * (t69 * t97 + t70 * t96) + m(4) * (t100 * t130 + t101 * t129) + m(3) * (t160 * t177 + t161 * t176) + ((-t109 / 0.2e1 - t137 / 0.2e1 - t263) * t260 + (t138 / 0.2e1 + t110 / 0.2e1 + t262) * t257) * t252; (t30 + t33 + t41 + t146) * t253 + m(6) * (t38 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t61 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(4) * (t100 ^ 2 + t101 ^ 2 + t95 ^ 2) + m(3) * (t139 ^ 2 + t160 ^ 2 + t161 ^ 2) + (-t260 * t15 + t257 * t16 + t257 * t22 - t260 * t21 + t257 * t37 - t260 * t36 + (-t260 * ((-t192 * t234 + t194 * t235) * t257 - (-t191 * t234 + t193 * t235) * t260) + t257 * ((-t192 * t236 + t194 * t237) * t257 - (-t191 * t236 + t193 * t237) * t260) + (-t260 * (t189 * t260 ^ 2 - t190 * t299) + t257 * (t190 * t257 ^ 2 - t189 * t299)) * t252) * t252 + ((-t109 - t137) * t260 + (t110 + t138) * t257) * t253) * t252; t284 * t304 + m(6) * (t45 * t91 + t46 * t92) + m(5) * (t74 * t96 + t75 * t97) + m(4) * (t117 * t129 + t118 * t130) + t262 * t236 + t263 * t234; (t28 / 0.2e1 + t32 / 0.2e1 + t40 / 0.2e1) * t253 + (t16 / 0.2e1 + t22 / 0.2e1 + t37 / 0.2e1) * t236 + (t15 / 0.2e1 + t21 / 0.2e1 + t36 / 0.2e1) * t234 + m(6) * (t38 * t39 + t43 * t46 + t44 * t45) + m(5) * (t61 * t64 + t69 * t75 + t70 * t74) + m(4) * (t100 * t118 + t101 * t117 + t102 * t95) + ((-t13 / 0.2e1 - t19 / 0.2e1 - t34 / 0.2e1) * t260 + (-t30 / 0.2e1 - t33 / 0.2e1 - t41 / 0.2e1) * t259 + (t14 / 0.2e1 + t20 / 0.2e1 + t35 / 0.2e1) * t257) * t252; (-t28 - t32 - t40) * t304 + (t14 + t20 + t35) * t236 + (t13 + t19 + t34) * t234 + m(6) * (t39 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(5) * (t64 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(4) * (t102 ^ 2 + t117 ^ 2 + t118 ^ 2); t80 + t283 * t217 + t282 * t215 + m(6) * (t51 * t91 + t52 * t92) + m(5) * (t93 * t96 + t94 * t97) + t271; t33 * t315 + t22 * t316 + t21 * t317 + t31 * t312 + (t257 * t18 / 0.2e1 - t260 * t17 / 0.2e1) * t252 + m(6) * (t38 * t42 + t43 * t52 + t44 * t51) + m(5) * (t61 * t79 + t69 * t94 + t70 * t93) + t268; t31 * t275 + t32 * t315 + t18 * t313 + t20 * t316 + t17 * t314 + t19 * t317 + m(6) * (t39 * t42 + t45 * t51 + t46 * t52) + m(5) * (t64 * t79 + t74 * t93 + t75 * t94) + t269; t215 * t17 + t217 * t18 + t232 * t31 + m(6) * (t42 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(5) * (t79 ^ 2 + t93 ^ 2 + t94 ^ 2) + t286; m(6) * (t87 * t91 + t88 * t92) + t271; m(6) * (t38 * t77 + t43 * t88 + t44 * t87) + t268; m(6) * (t39 * t77 + t45 * t87 + t46 * t88) + t269; m(6) * (t42 * t77 + t51 * t87 + t52 * t88) + t286; m(6) * (t77 ^ 2 + t87 ^ 2 + t88 ^ 2) + t286;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
