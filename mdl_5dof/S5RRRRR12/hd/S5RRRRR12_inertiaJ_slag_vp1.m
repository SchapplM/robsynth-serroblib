% Calculate joint inertia matrix for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR12_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR12_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:40
% EndTime: 2019-12-31 22:47:00
% DurationCPUTime: 6.46s
% Computational Cost: add. (37313->541), mult. (103395->775), div. (0->0), fcn. (136008->14), ass. (0->263)
t243 = cos(pkin(5));
t248 = sin(qJ(1));
t250 = cos(qJ(2));
t297 = t248 * t250;
t247 = sin(qJ(2));
t251 = cos(qJ(1));
t299 = t247 * t251;
t229 = -t243 * t297 - t299;
t242 = sin(pkin(5));
t305 = cos(pkin(6));
t271 = t242 * t305;
t304 = sin(pkin(6));
t210 = -t229 * t304 + t248 * t271;
t295 = t250 * t251;
t298 = t248 * t247;
t227 = t243 * t295 - t298;
t209 = -t227 * t304 - t251 * t271;
t228 = t243 * t299 + t297;
t246 = sin(qJ(3));
t270 = t242 * t304;
t310 = cos(qJ(3));
t190 = t228 * t310 + (t305 * t227 - t251 * t270) * t246;
t245 = sin(qJ(4));
t309 = cos(qJ(4));
t166 = t190 * t245 - t209 * t309;
t230 = -t243 * t298 + t295;
t192 = t230 * t310 + (t305 * t229 + t248 * t270) * t246;
t168 = t192 * t245 - t210 * t309;
t208 = t243 * t304 * t246 + (t305 * t246 * t250 + t310 * t247) * t242;
t226 = t243 * t305 - t250 * t270;
t183 = t208 * t245 - t226 * t309;
t167 = t190 * t309 + t209 * t245;
t263 = t310 * t304;
t260 = t242 * t263;
t264 = t305 * t310;
t189 = -t227 * t264 + t228 * t246 + t251 * t260;
t244 = sin(qJ(5));
t249 = cos(qJ(5));
t127 = -t167 * t244 + t189 * t249;
t128 = t167 * t249 + t189 * t244;
t85 = Icges(6,5) * t128 + Icges(6,6) * t127 + Icges(6,3) * t166;
t87 = Icges(6,4) * t128 + Icges(6,2) * t127 + Icges(6,6) * t166;
t89 = Icges(6,1) * t128 + Icges(6,4) * t127 + Icges(6,5) * t166;
t24 = t127 * t87 + t128 * t89 + t166 * t85;
t169 = t192 * t309 + t210 * t245;
t191 = -t229 * t264 + t230 * t246 - t248 * t260;
t129 = -t169 * t244 + t191 * t249;
t130 = t169 * t249 + t191 * t244;
t86 = Icges(6,5) * t130 + Icges(6,6) * t129 + Icges(6,3) * t168;
t88 = Icges(6,4) * t130 + Icges(6,2) * t129 + Icges(6,6) * t168;
t90 = Icges(6,1) * t130 + Icges(6,4) * t129 + Icges(6,5) * t168;
t25 = t127 * t88 + t128 * t90 + t166 * t86;
t184 = t208 * t309 + t226 * t245;
t301 = t242 * t250;
t303 = t242 * t247;
t207 = -t243 * t263 + t246 * t303 - t264 * t301;
t160 = -t184 * t244 + t207 * t249;
t161 = t184 * t249 + t207 * t244;
t106 = Icges(6,5) * t161 + Icges(6,6) * t160 + Icges(6,3) * t183;
t107 = Icges(6,4) * t161 + Icges(6,2) * t160 + Icges(6,6) * t183;
t108 = Icges(6,1) * t161 + Icges(6,4) * t160 + Icges(6,5) * t183;
t41 = t106 * t166 + t107 * t127 + t108 * t128;
t1 = t166 * t24 + t168 * t25 + t183 * t41;
t316 = t1 / 0.2e1;
t26 = t129 * t87 + t130 * t89 + t168 * t85;
t27 = t129 * t88 + t130 * t90 + t168 * t86;
t42 = t106 * t168 + t107 * t129 + t108 * t130;
t2 = t166 * t26 + t168 * t27 + t183 * t42;
t315 = t2 / 0.2e1;
t32 = t160 * t87 + t161 * t89 + t183 * t85;
t33 = t160 * t88 + t161 * t90 + t183 * t86;
t49 = t183 * t106 + t160 * t107 + t161 * t108;
t45 = t49 * t183;
t9 = t32 * t166 + t33 * t168 + t45;
t314 = t9 / 0.2e1;
t313 = t166 / 0.2e1;
t312 = t168 / 0.2e1;
t311 = t183 / 0.2e1;
t308 = pkin(4) * t167;
t262 = -rSges(6,1) * t128 - rSges(6,2) * t127;
t91 = rSges(6,3) * t166 - t262;
t307 = pkin(11) * t166 + t308 + t91;
t92 = t130 * rSges(6,1) + t129 * rSges(6,2) + t168 * rSges(6,3);
t306 = t169 * pkin(4) + t168 * pkin(11) + t92;
t302 = t242 * t248;
t300 = t242 * t251;
t296 = t248 * t251;
t109 = rSges(6,1) * t161 + rSges(6,2) * t160 + rSges(6,3) * t183;
t294 = pkin(4) * t184 + pkin(11) * t183 + t109;
t116 = rSges(5,1) * t167 - rSges(5,2) * t166 + rSges(5,3) * t189;
t153 = pkin(3) * t190 + t189 * pkin(10);
t293 = -t116 - t153;
t136 = rSges(5,1) * t184 - rSges(5,2) * t183 + rSges(5,3) * t207;
t177 = pkin(3) * t208 + pkin(10) * t207;
t292 = -t136 - t177;
t154 = t192 * pkin(3) + t191 * pkin(10);
t195 = t230 * pkin(2) + t210 * pkin(9);
t193 = t243 * t195;
t291 = t243 * t154 + t193;
t194 = pkin(2) * t228 + t209 * pkin(9);
t290 = t194 * t302 + t195 * t300;
t289 = t251 * pkin(1) + pkin(8) * t302;
t110 = Icges(5,5) * t167 - Icges(5,6) * t166 + Icges(5,3) * t189;
t112 = Icges(5,4) * t167 - Icges(5,2) * t166 + Icges(5,6) * t189;
t114 = Icges(5,1) * t167 - Icges(5,4) * t166 + Icges(5,5) * t189;
t51 = t110 * t189 - t112 * t166 + t114 * t167;
t111 = Icges(5,5) * t169 - Icges(5,6) * t168 + Icges(5,3) * t191;
t113 = Icges(5,4) * t169 - Icges(5,2) * t168 + Icges(5,6) * t191;
t115 = Icges(5,1) * t169 - Icges(5,4) * t168 + Icges(5,5) * t191;
t52 = t111 * t189 - t113 * t166 + t115 * t167;
t131 = Icges(5,5) * t184 - Icges(5,6) * t183 + Icges(5,3) * t207;
t132 = Icges(5,4) * t184 - Icges(5,2) * t183 + Icges(5,6) * t207;
t133 = Icges(5,1) * t184 - Icges(5,4) * t183 + Icges(5,5) * t207;
t63 = t131 * t189 - t132 * t166 + t133 * t167;
t13 = t189 * t51 + t191 * t52 + t207 * t63;
t3 = t189 * t24 + t191 * t25 + t207 * t41;
t288 = t3 / 0.2e1 + t13 / 0.2e1;
t53 = t110 * t191 - t112 * t168 + t114 * t169;
t54 = t111 * t191 - t113 * t168 + t115 * t169;
t64 = t131 * t191 - t132 * t168 + t133 * t169;
t14 = t189 * t53 + t191 * t54 + t207 * t64;
t4 = t189 * t26 + t191 * t27 + t207 * t42;
t287 = t4 / 0.2e1 + t14 / 0.2e1;
t15 = t209 * t51 + t210 * t52 + t226 * t63;
t5 = t209 * t24 + t210 * t25 + t226 * t41;
t286 = t5 / 0.2e1 + t15 / 0.2e1;
t16 = t209 * t53 + t210 * t54 + t226 * t64;
t6 = t209 * t26 + t210 * t27 + t226 * t42;
t285 = t6 / 0.2e1 + t16 / 0.2e1;
t17 = t63 * t243 + (t248 * t52 - t251 * t51) * t242;
t7 = t41 * t243 + (-t24 * t251 + t248 * t25) * t242;
t284 = t7 / 0.2e1 + t17 / 0.2e1;
t18 = t64 * t243 + (t248 * t54 - t251 * t53) * t242;
t8 = t42 * t243 + (t248 * t27 - t251 * t26) * t242;
t283 = t8 / 0.2e1 + t18 / 0.2e1;
t46 = t49 * t207;
t10 = t32 * t189 + t33 * t191 + t46;
t55 = t110 * t207 - t112 * t183 + t114 * t184;
t56 = t111 * t207 - t113 * t183 + t115 * t184;
t73 = t207 * t131 - t183 * t132 + t184 * t133;
t69 = t73 * t207;
t19 = t55 * t189 + t56 * t191 + t69;
t282 = t10 / 0.2e1 + t19 / 0.2e1;
t47 = t49 * t226;
t11 = t32 * t209 + t33 * t210 + t47;
t70 = t73 * t226;
t20 = t55 * t209 + t56 * t210 + t70;
t281 = t11 / 0.2e1 + t20 / 0.2e1;
t48 = t49 * t243;
t12 = t48 + (t33 * t248 - t32 * t251) * t242;
t72 = t73 * t243;
t21 = t72 + (t56 * t248 - t55 * t251) * t242;
t280 = t12 / 0.2e1 + t21 / 0.2e1;
t279 = t33 / 0.2e1 + t42 / 0.2e1;
t278 = t41 / 0.2e1 + t32 / 0.2e1;
t277 = -t153 - t307;
t276 = -t177 - t294;
t170 = Icges(4,5) * t208 - Icges(4,6) * t207 + Icges(4,3) * t226;
t171 = Icges(4,4) * t208 - Icges(4,2) * t207 + Icges(4,6) * t226;
t172 = Icges(4,1) * t208 - Icges(4,4) * t207 + Icges(4,5) * t226;
t100 = t226 * t170 - t207 * t171 + t208 * t172;
t117 = t169 * rSges(5,1) - t168 * rSges(5,2) + t191 * rSges(5,3);
t144 = t192 * rSges(4,1) - t191 * rSges(4,2) + t210 * rSges(4,3);
t217 = Icges(3,3) * t243 + (Icges(3,5) * t247 + Icges(3,6) * t250) * t242;
t218 = Icges(3,6) * t243 + (Icges(3,4) * t247 + Icges(3,2) * t250) * t242;
t219 = Icges(3,5) * t243 + (Icges(3,1) * t247 + Icges(3,4) * t250) * t242;
t275 = t243 * t217 + t218 * t301 + t219 * t303;
t203 = t230 * rSges(3,1) + t229 * rSges(3,2) + rSges(3,3) * t302;
t274 = -t248 * pkin(1) + pkin(8) * t300;
t173 = rSges(4,1) * t208 - rSges(4,2) * t207 + rSges(4,3) * t226;
t213 = pkin(2) * t303 + t226 * pkin(9);
t269 = t242 * (-t173 - t213);
t268 = t153 * t302 + t154 * t300 + t290;
t265 = t242 * (-t213 + t292);
t261 = t242 * (-t213 + t276);
t259 = t56 / 0.2e1 + t64 / 0.2e1 + t279;
t258 = t55 / 0.2e1 + t63 / 0.2e1 + t278;
t143 = rSges(4,1) * t190 - rSges(4,2) * t189 + rSges(4,3) * t209;
t257 = -t194 + t274;
t202 = t228 * rSges(3,1) + t227 * rSges(3,2) - rSges(3,3) * t300;
t138 = Icges(4,5) * t192 - Icges(4,6) * t191 + Icges(4,3) * t210;
t140 = Icges(4,4) * t192 - Icges(4,2) * t191 + Icges(4,6) * t210;
t142 = Icges(4,1) * t192 - Icges(4,4) * t191 + Icges(4,5) * t210;
t79 = t138 * t226 - t140 * t207 + t142 * t208;
t94 = t170 * t210 - t171 * t191 + t172 * t192;
t256 = t94 / 0.2e1 + t79 / 0.2e1 + t259;
t137 = Icges(4,5) * t190 - Icges(4,6) * t189 + Icges(4,3) * t209;
t139 = Icges(4,4) * t190 - Icges(4,2) * t189 + Icges(4,6) * t209;
t141 = Icges(4,1) * t190 - Icges(4,4) * t189 + Icges(4,5) * t209;
t78 = t137 * t226 - t139 * t207 + t141 * t208;
t93 = t170 * t209 - t171 * t189 + t172 * t190;
t255 = t93 / 0.2e1 + t78 / 0.2e1 + t258;
t254 = t195 + t289;
t253 = -t153 + t257;
t252 = t154 + t254;
t236 = rSges(2,1) * t251 - t248 * rSges(2,2);
t235 = -t248 * rSges(2,1) - rSges(2,2) * t251;
t220 = rSges(3,3) * t243 + (rSges(3,1) * t247 + rSges(3,2) * t250) * t242;
t201 = Icges(3,1) * t230 + Icges(3,4) * t229 + Icges(3,5) * t302;
t200 = Icges(3,1) * t228 + Icges(3,4) * t227 - Icges(3,5) * t300;
t199 = Icges(3,4) * t230 + Icges(3,2) * t229 + Icges(3,6) * t302;
t198 = Icges(3,4) * t228 + Icges(3,2) * t227 - Icges(3,6) * t300;
t197 = Icges(3,5) * t230 + Icges(3,6) * t229 + Icges(3,3) * t302;
t196 = Icges(3,5) * t228 + Icges(3,6) * t227 - Icges(3,3) * t300;
t188 = t203 + t289;
t187 = -t202 + t274;
t176 = -t243 * t202 - t220 * t300;
t175 = t203 * t243 - t220 * t302;
t174 = t275 * t243;
t159 = t209 * t177;
t158 = (t202 * t248 + t203 * t251) * t242;
t157 = t217 * t302 + t218 * t229 + t219 * t230;
t156 = -t217 * t300 + t227 * t218 + t228 * t219;
t146 = t226 * t154;
t145 = t210 * t153;
t135 = t197 * t243 + (t199 * t250 + t201 * t247) * t242;
t134 = t196 * t243 + (t198 * t250 + t200 * t247) * t242;
t121 = t254 + t144;
t120 = -t143 + t257;
t105 = t144 * t226 - t173 * t210;
t104 = -t143 * t226 + t173 * t209;
t102 = (-t143 - t194) * t243 + t251 * t269;
t101 = t144 * t243 + t248 * t269 + t193;
t99 = t100 * t243;
t96 = t143 * t210 - t144 * t209;
t95 = t100 * t226;
t84 = (t143 * t248 + t144 * t251) * t242 + t290;
t83 = t252 + t117;
t82 = -t116 + t253;
t81 = t117 * t207 - t136 * t191;
t80 = -t116 * t207 + t136 * t189;
t77 = t138 * t210 - t140 * t191 + t142 * t192;
t76 = t137 * t210 - t139 * t191 + t141 * t192;
t75 = t138 * t209 - t140 * t189 + t142 * t190;
t74 = t137 * t209 - t139 * t189 + t141 * t190;
t71 = t116 * t191 - t117 * t189;
t68 = t117 * t226 + t292 * t210 + t146;
t67 = t136 * t209 + t293 * t226 + t159;
t66 = (-t194 + t293) * t243 + t251 * t265;
t65 = t117 * t243 + t248 * t265 + t291;
t62 = t252 + t306;
t61 = -t308 + (-rSges(6,3) - pkin(11)) * t166 + t253 + t262;
t60 = -t109 * t168 + t183 * t92;
t59 = t109 * t166 - t183 * t91;
t58 = t116 * t210 + t145 + (-t117 - t154) * t209;
t57 = (t116 * t248 + t117 * t251) * t242 + t268;
t50 = -t166 * t92 + t168 * t91;
t44 = -t294 * t191 + t306 * t207;
t43 = t294 * t189 - t307 * t207;
t40 = (-t194 + t277) * t243 + t251 * t261;
t39 = t306 * t243 + t248 * t261 + t291;
t38 = t276 * t210 + t306 * t226 + t146;
t37 = t294 * t209 + t277 * t226 + t159;
t36 = t99 + (t79 * t248 - t78 * t251) * t242;
t35 = -t306 * t189 + t307 * t191;
t34 = t78 * t209 + t79 * t210 + t95;
t31 = (t307 * t248 + t306 * t251) * t242 + t268;
t30 = t94 * t243 + (t248 * t77 - t251 * t76) * t242;
t29 = t93 * t243 + (t248 * t75 - t251 * t74) * t242;
t28 = t145 + t307 * t210 + (-t154 - t306) * t209;
t23 = t209 * t76 + t210 * t77 + t226 * t94;
t22 = t209 * t74 + t210 * t75 + t226 * t93;
t97 = [Icges(2,3) + m(6) * (t61 ^ 2 + t62 ^ 2) + m(5) * (t82 ^ 2 + t83 ^ 2) + m(4) * (t120 ^ 2 + t121 ^ 2) + m(3) * (t187 ^ 2 + t188 ^ 2) + m(2) * (t235 ^ 2 + t236 ^ 2) + t275 + t100 + t73 + t49; t48 + t72 + t99 + t174 + m(6) * (t39 * t62 + t40 * t61) + m(5) * (t65 * t83 + t66 * t82) + m(4) * (t101 * t121 + t102 * t120) + m(3) * (t175 * t188 + t176 * t187) + ((-t156 / 0.2e1 - t134 / 0.2e1 - t255) * t251 + (t157 / 0.2e1 + t135 / 0.2e1 + t256) * t248) * t242; (t12 + t21 + t36 + t174) * t243 + m(6) * (t31 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(5) * (t57 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(4) * (t101 ^ 2 + t102 ^ 2 + t84 ^ 2) + m(3) * (t158 ^ 2 + t175 ^ 2 + t176 ^ 2) + (-t251 * t7 + t248 * t8 + t248 * t18 - t251 * t17 + t248 * t30 - t251 * t29 + (-t251 * ((t227 * t199 + t228 * t201) * t248 - (t227 * t198 + t228 * t200) * t251) + t248 * ((t199 * t229 + t201 * t230) * t248 - (t198 * t229 + t200 * t230) * t251) + (-t251 * (t196 * t251 ^ 2 - t197 * t296) + t248 * (t197 * t248 ^ 2 - t196 * t296)) * t242) * t242 + ((-t134 - t156) * t251 + (t135 + t157) * t248) * t243) * t242; t47 + t70 + t95 + m(6) * (t37 * t61 + t38 * t62) + m(5) * (t67 * t82 + t68 * t83) + m(4) * (t104 * t120 + t105 * t121) + t256 * t210 + t255 * t209; (t34 / 0.2e1 + t281) * t243 + (t36 / 0.2e1 + t280) * t226 + (t30 / 0.2e1 + t283) * t210 + (t29 / 0.2e1 + t284) * t209 + m(6) * (t28 * t31 + t37 * t40 + t38 * t39) + m(5) * (t58 * t57 + t65 * t68 + t67 * t66) + m(4) * (t101 * t105 + t102 * t104 + t84 * t96) + ((-t22 / 0.2e1 - t286) * t251 + (t23 / 0.2e1 + t285) * t248) * t242; (t11 + t20 + t34) * t226 + (t6 + t16 + t23) * t210 + (t5 + t15 + t22) * t209 + m(6) * (t28 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(5) * (t58 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(4) * (t104 ^ 2 + t105 ^ 2 + t96 ^ 2); t46 + t69 + m(6) * (t43 * t61 + t44 * t62) + m(5) * (t80 * t82 + t81 * t83) + t259 * t191 + t258 * t189; t282 * t243 + t280 * t207 + t283 * t191 + t284 * t189 + m(6) * (t31 * t35 + t39 * t44 + t40 * t43) + m(5) * (t57 * t71 + t65 * t81 + t66 * t80) + (t287 * t248 - t288 * t251) * t242; t282 * t226 + t287 * t210 + t288 * t209 + t281 * t207 + t285 * t191 + t286 * t189 + m(6) * (t28 * t35 + t37 * t43 + t38 * t44) + m(5) * (t58 * t71 + t67 * t80 + t68 * t81); (t10 + t19) * t207 + (t4 + t14) * t191 + (t3 + t13) * t189 + m(6) * (t35 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t71 ^ 2 + t80 ^ 2 + t81 ^ 2); m(6) * (t59 * t61 + t60 * t62) + t45 + t279 * t168 + t278 * t166; m(6) * (t31 * t50 + t39 * t60 + t40 * t59) + t243 * t314 + t8 * t312 + t7 * t313 + t12 * t311 + (-t251 * t1 / 0.2e1 + t248 * t315) * t242; m(6) * (t28 * t50 + t37 * t59 + t38 * t60) + t226 * t314 + t210 * t315 + t209 * t316 + t11 * t311 + t6 * t312 + t5 * t313; m(6) * (t35 * t50 + t43 * t59 + t44 * t60) + t10 * t311 + t207 * t314 + t189 * t316 + t3 * t313 + t4 * t312 + t191 * t315; m(6) * (t50 ^ 2 + t59 ^ 2 + t60 ^ 2) + t168 * t2 + t166 * t1 + t183 * t9;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t97(1), t97(2), t97(4), t97(7), t97(11); t97(2), t97(3), t97(5), t97(8), t97(12); t97(4), t97(5), t97(6), t97(9), t97(13); t97(7), t97(8), t97(9), t97(10), t97(14); t97(11), t97(12), t97(13), t97(14), t97(15);];
Mq = res;
