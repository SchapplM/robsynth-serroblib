% Calculate joint inertia matrix for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:53:02
% EndTime: 2019-03-09 21:53:11
% DurationCPUTime: 3.97s
% Computational Cost: add. (11908->418), mult. (8751->606), div. (0->0), fcn. (8968->12), ass. (0->215)
t308 = Icges(5,3) + Icges(6,3);
t199 = qJ(2) + qJ(3);
t189 = qJ(4) + t199;
t182 = pkin(11) + t189;
t177 = sin(t182);
t178 = cos(t182);
t183 = sin(t189);
t184 = cos(t189);
t307 = Icges(5,5) * t184 + Icges(6,5) * t178 - Icges(5,6) * t183 - Icges(6,6) * t177;
t202 = sin(qJ(1));
t205 = cos(qJ(1));
t306 = t308 * t202 + t307 * t205;
t305 = -t307 * t202 + t308 * t205;
t271 = Icges(6,4) * t178;
t235 = -Icges(6,2) * t177 + t271;
t103 = Icges(6,6) * t202 + t235 * t205;
t272 = Icges(6,4) * t177;
t239 = Icges(6,1) * t178 - t272;
t105 = Icges(6,5) * t202 + t239 * t205;
t273 = Icges(5,4) * t184;
t236 = -Icges(5,2) * t183 + t273;
t121 = Icges(5,6) * t202 + t236 * t205;
t274 = Icges(5,4) * t183;
t240 = Icges(5,1) * t184 - t274;
t123 = Icges(5,5) * t202 + t240 * t205;
t304 = t103 * t177 - t105 * t178 + t121 * t183 - t123 * t184;
t102 = -Icges(6,6) * t205 + t235 * t202;
t104 = -Icges(6,5) * t205 + t239 * t202;
t120 = -Icges(5,6) * t205 + t236 * t202;
t122 = -Icges(5,5) * t205 + t240 * t202;
t303 = t102 * t177 - t104 * t178 + t120 * t183 - t122 * t184;
t197 = t202 ^ 2;
t302 = t202 * pkin(7);
t268 = t178 * t205;
t269 = t177 * t205;
t203 = cos(qJ(6));
t265 = t202 * t203;
t200 = sin(qJ(6));
t266 = t200 * t205;
t137 = -t178 * t266 + t265;
t264 = t203 * t205;
t267 = t200 * t202;
t138 = t178 * t264 + t267;
t76 = t138 * rSges(7,1) + t137 * rSges(7,2) + rSges(7,3) * t269;
t301 = pkin(5) * t268 + pkin(10) * t269 + t76;
t300 = Icges(5,5) * t183 + Icges(6,5) * t177 + Icges(5,6) * t184 + Icges(6,6) * t178;
t147 = Icges(6,2) * t178 + t272;
t148 = Icges(6,1) * t177 + t271;
t153 = Icges(5,2) * t184 + t274;
t154 = Icges(5,1) * t183 + t273;
t299 = -t147 * t177 + t148 * t178 - t153 * t183 + t154 * t184;
t245 = rSges(5,1) * t184 - rSges(5,2) * t183;
t186 = sin(t199);
t187 = cos(t199);
t246 = rSges(4,1) * t187 - rSges(4,2) * t186;
t95 = -t178 * rSges(7,3) + (rSges(7,1) * t203 - rSges(7,2) * t200) * t177;
t298 = -pkin(5) * t177 + pkin(10) * t178 - t95;
t233 = Icges(4,5) * t187 - Icges(4,6) * t186;
t126 = -Icges(4,3) * t205 + t233 * t202;
t127 = Icges(4,3) * t202 + t233 * t205;
t198 = t205 ^ 2;
t275 = Icges(4,4) * t187;
t237 = -Icges(4,2) * t186 + t275;
t129 = Icges(4,6) * t202 + t237 * t205;
t276 = Icges(4,4) * t186;
t241 = Icges(4,1) * t187 - t276;
t131 = Icges(4,5) * t202 + t241 * t205;
t225 = -t129 * t186 + t131 * t187;
t128 = -Icges(4,6) * t205 + t237 * t202;
t130 = -Icges(4,5) * t205 + t241 * t202;
t226 = t128 * t186 - t130 * t187;
t135 = -t178 * t267 - t264;
t136 = t178 * t265 - t266;
t270 = t177 * t202;
t67 = Icges(7,5) * t136 + Icges(7,6) * t135 + Icges(7,3) * t270;
t69 = Icges(7,4) * t136 + Icges(7,2) * t135 + Icges(7,6) * t270;
t71 = Icges(7,1) * t136 + Icges(7,4) * t135 + Icges(7,5) * t270;
t18 = t135 * t69 + t136 * t71 + t67 * t270;
t68 = Icges(7,5) * t138 + Icges(7,6) * t137 + Icges(7,3) * t269;
t70 = Icges(7,4) * t138 + Icges(7,2) * t137 + Icges(7,6) * t269;
t72 = Icges(7,1) * t138 + Icges(7,4) * t137 + Icges(7,5) * t269;
t19 = t135 * t70 + t136 * t72 + t68 * t270;
t8 = -t18 * t205 + t19 * t202;
t256 = -t8 + t305 * t198 + (t304 * t202 + (-t303 + t306) * t205) * t202;
t297 = -t198 * t126 - (t225 * t202 + (-t127 + t226) * t205) * t202 + t256;
t206 = -pkin(8) - pkin(7);
t296 = t202 / 0.2e1;
t295 = -t205 / 0.2e1;
t201 = sin(qJ(2));
t294 = pkin(2) * t201;
t293 = pkin(3) * t186;
t292 = pkin(4) * t183;
t291 = pkin(5) * t178;
t204 = cos(qJ(2));
t185 = t204 * pkin(2) + pkin(1);
t162 = pkin(3) * t187 + t185;
t145 = pkin(4) * t184 + t162;
t132 = t205 * t145;
t156 = t205 * t162;
t290 = t205 * (t132 - t156) + (t145 - t162) * t197;
t175 = t205 * t185;
t289 = t205 * (t156 - t175) + (t162 - t185) * t197;
t288 = rSges(3,1) * t204;
t285 = rSges(3,2) * t201;
t93 = -Icges(7,6) * t178 + (Icges(7,4) * t203 - Icges(7,2) * t200) * t177;
t282 = t200 * t93;
t281 = t205 * rSges(3,3);
t26 = -t178 * t67 + (-t200 * t69 + t203 * t71) * t177;
t280 = t26 * t205;
t27 = -t178 * t68 + (-t200 * t70 + t203 * t72) * t177;
t279 = t27 * t202;
t278 = Icges(3,4) * t201;
t277 = Icges(3,4) * t204;
t218 = t202 * rSges(5,3) + t245 * t205;
t66 = t202 * (-t205 * rSges(5,3) + t245 * t202) + t205 * t218;
t195 = t205 * pkin(7);
t263 = t202 * (t195 + (-pkin(1) + t185) * t202) + t205 * (-t205 * pkin(1) + t175 - t302);
t219 = t202 * rSges(4,3) + t246 * t205;
t77 = t202 * (-t205 * rSges(4,3) + t246 * t202) + t205 * t219;
t261 = t202 * rSges(3,3) + t205 * t288;
t258 = t197 + t198;
t196 = -pkin(9) + t206;
t20 = t137 * t69 + t138 * t71 + t67 * t269;
t21 = t137 * t70 + t138 * t72 + t68 * t269;
t9 = -t20 * t205 + t202 * t21;
t257 = (t9 + t306 * t197 + ((-t304 + t305) * t202 + t303 * t205) * t205) * t202;
t255 = t202 * (t197 * t127 + (t226 * t205 + (-t126 + t225) * t202) * t205) + t257;
t161 = rSges(4,1) * t186 + rSges(4,2) * t187;
t254 = -t161 - t294;
t155 = rSges(5,1) * t183 + rSges(5,2) * t184;
t253 = -t155 - t293;
t149 = rSges(6,1) * t177 + rSges(6,2) * t178;
t252 = -t149 - t292;
t217 = rSges(6,1) * t268 - rSges(6,2) * t269 + t202 * rSges(6,3);
t244 = rSges(6,1) * t178 - rSges(6,2) * t177;
t31 = t202 * (-t205 * rSges(6,3) + t244 * t202) + t205 * t217 + t290;
t188 = -qJ(5) + t196;
t251 = -t188 * t202 + t132;
t36 = t66 + t289;
t92 = -Icges(7,3) * t178 + (Icges(7,5) * t203 - Icges(7,6) * t200) * t177;
t94 = -Icges(7,5) * t178 + (Icges(7,1) * t203 - Icges(7,4) * t200) * t177;
t34 = t135 * t93 + t136 * t94 + t92 * t270;
t3 = -t34 * t178 + (t18 * t202 + t19 * t205) * t177;
t35 = t137 * t93 + t138 * t94 + t92 * t269;
t4 = -t35 * t178 + (t20 * t202 + t205 * t21) * t177;
t250 = t3 * t295 + t4 * t296 - t178 * (t279 - t280) / 0.2e1 + t8 * t270 / 0.2e1 + t9 * t269 / 0.2e1;
t249 = -t292 + t298;
t248 = -t292 - t293;
t247 = -t285 + t288;
t243 = -rSges(7,1) * t136 - rSges(7,2) * t135;
t15 = t31 + t289;
t242 = Icges(3,1) * t204 - t278;
t238 = -Icges(3,2) * t201 + t277;
t234 = Icges(3,5) * t204 - Icges(3,6) * t201;
t159 = Icges(4,2) * t187 + t276;
t160 = Icges(4,1) * t186 + t275;
t220 = -t159 * t186 + t160 * t187;
t75 = rSges(7,3) * t270 - t243;
t14 = t202 * t75 + t197 * (pkin(10) * t177 + t291) + t290 + t301 * t205;
t216 = t253 - t294;
t215 = -t149 + t248;
t214 = t248 + t298;
t213 = t248 - t294;
t212 = t256 * t205 + t257;
t12 = t14 + t289;
t211 = -t149 + t213;
t210 = t205 * t297 + t255;
t209 = t213 + t298;
t208 = -t280 / 0.2e1 + t279 / 0.2e1 + (t103 * t178 + t105 * t177 + t121 * t184 + t123 * t183 + t300 * t202 + t299 * t205 + t35) * t296 + (t102 * t178 + t104 * t177 + t120 * t184 + t122 * t183 + t299 * t202 - t300 * t205 + t34) * t295;
t158 = Icges(4,5) * t186 + Icges(4,6) * t187;
t207 = t208 + (t129 * t187 + t131 * t186 + t202 * t158 + t220 * t205) * t296 + (t128 * t187 + t130 * t186 - t205 * t158 + t220 * t202) * t295;
t174 = rSges(2,1) * t205 - rSges(2,2) * t202;
t173 = -rSges(2,1) * t202 - rSges(2,2) * t205;
t172 = rSges(3,1) * t201 + rSges(3,2) * t204;
t140 = Icges(3,3) * t202 + t234 * t205;
t139 = -Icges(3,3) * t205 + t234 * t202;
t125 = t254 * t205;
t124 = t254 * t202;
t111 = t302 + (pkin(1) - t285) * t205 + t261;
t110 = t281 + t195 + (-pkin(1) - t247) * t202;
t109 = t253 * t205;
t108 = t253 * t202;
t97 = t252 * t205;
t96 = t252 * t202;
t91 = t216 * t205;
t90 = t216 * t202;
t89 = -t202 * t206 + t175 + t219;
t88 = (rSges(4,3) - t206) * t205 + (-t185 - t246) * t202;
t85 = t215 * t205;
t84 = t215 * t202;
t83 = t177 * t203 * t94;
t82 = t205 * (-t205 * t285 + t261) + (t247 * t202 - t281) * t202;
t81 = -t202 * t196 + t156 + t218;
t80 = (rSges(5,3) - t196) * t205 + (-t162 - t245) * t202;
t79 = t211 * t205;
t78 = t211 * t202;
t65 = t217 + t251;
t64 = (rSges(6,3) - t188) * t205 + (-t145 - t244) * t202;
t57 = t249 * t205;
t56 = t249 * t202;
t51 = t214 * t205;
t50 = t214 * t202;
t47 = t209 * t205;
t46 = t209 * t202;
t43 = t77 + t263;
t42 = -t178 * t76 - t95 * t269;
t41 = t178 * t75 + t95 * t270;
t40 = t251 + t301;
t39 = -t188 * t205 + (-t291 - t145 + (-rSges(7,3) - pkin(10)) * t177) * t202 + t243;
t38 = -t177 * t282 - t178 * t92 + t83;
t37 = (-t202 * t76 + t205 * t75) * t177;
t28 = t36 + t263;
t13 = t15 + t263;
t10 = t12 + t263;
t1 = [t184 * t153 + t183 * t154 + t187 * t159 + t186 * t160 + t204 * (Icges(3,2) * t204 + t278) + t201 * (Icges(3,1) * t201 + t277) + Icges(2,3) + t83 + (-t92 + t147) * t178 + (t148 - t282) * t177 + m(7) * (t39 ^ 2 + t40 ^ 2) + m(6) * (t64 ^ 2 + t65 ^ 2) + m(5) * (t80 ^ 2 + t81 ^ 2) + m(4) * (t88 ^ 2 + t89 ^ 2) + m(3) * (t110 ^ 2 + t111 ^ 2) + m(2) * (t173 ^ 2 + t174 ^ 2); m(7) * (t39 * t47 + t40 * t46) + m(6) * (t64 * t79 + t65 * t78) + m(5) * (t80 * t91 + t81 * t90) + m(4) * (t124 * t89 + t125 * t88) + (t204 * (-Icges(3,6) * t205 + t238 * t202) + t201 * (-Icges(3,5) * t205 + t242 * t202)) * t295 + (t204 * (Icges(3,6) * t202 + t238 * t205) + t201 * (Icges(3,5) * t202 + t242 * t205)) * t296 + m(3) * (-t110 * t205 - t111 * t202) * t172 + t207 + (t198 / 0.2e1 + t197 / 0.2e1) * (Icges(3,5) * t201 + Icges(3,6) * t204); m(7) * (t10 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(6) * (t13 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(5) * (t28 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(4) * (t124 ^ 2 + t125 ^ 2 + t43 ^ 2) + m(3) * (t258 * t172 ^ 2 + t82 ^ 2) + t202 * t197 * t140 + t255 + (-t198 * t139 + (-t202 * t139 + t205 * t140) * t202 + t297) * t205; m(7) * (t39 * t51 + t40 * t50) + m(6) * (t64 * t85 + t65 * t84) + m(5) * (t108 * t81 + t109 * t80) + m(4) * (-t202 * t89 - t205 * t88) * t161 + t207; m(7) * (t10 * t12 + t46 * t50 + t47 * t51) + m(6) * (t15 * t13 + t78 * t84 + t79 * t85) + m(5) * (t108 * t90 + t109 * t91 + t36 * t28) + m(4) * (t77 * t43 + (-t124 * t202 - t125 * t205) * t161) + t210; m(7) * (t12 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t15 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(5) * (t108 ^ 2 + t109 ^ 2 + t36 ^ 2) + m(4) * (t258 * t161 ^ 2 + t77 ^ 2) + t210; m(5) * (-t202 * t81 - t205 * t80) * t155 + m(7) * (t39 * t57 + t40 * t56) + m(6) * (t64 * t97 + t65 * t96) + t208; m(7) * (t10 * t14 + t46 * t56 + t47 * t57) + m(6) * (t31 * t13 + t78 * t96 + t79 * t97) + m(5) * (t66 * t28 + (-t202 * t90 - t205 * t91) * t155) + t212; m(7) * (t12 * t14 + t50 * t56 + t51 * t57) + m(6) * (t31 * t15 + t84 * t96 + t85 * t97) + m(5) * (t66 * t36 + (-t108 * t202 - t109 * t205) * t155) + t212; m(7) * (t14 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(6) * (t31 ^ 2 + t96 ^ 2 + t97 ^ 2) + m(5) * (t258 * t155 ^ 2 + t66 ^ 2) + t212; m(7) * (t202 * t39 - t205 * t40) + m(6) * (t202 * t64 - t205 * t65); m(7) * (t202 * t47 - t205 * t46) + m(6) * (t202 * t79 - t205 * t78); m(7) * (t202 * t51 - t205 * t50) + m(6) * (t202 * t85 - t205 * t84); m(7) * (t202 * t57 - t205 * t56) + m(6) * (t202 * t97 - t205 * t96); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t258; m(7) * (t39 * t41 + t40 * t42) - t38 * t178 + ((t27 / 0.2e1 + t35 / 0.2e1) * t205 + (t26 / 0.2e1 + t34 / 0.2e1) * t202) * t177; m(7) * (t10 * t37 + t41 * t47 + t42 * t46) + t250; m(7) * (t12 * t37 + t41 * t51 + t42 * t50) + t250; m(7) * (t14 * t37 + t41 * t57 + t42 * t56) + t250; m(7) * (t202 * t41 - t205 * t42); t178 ^ 2 * t38 + m(7) * (t37 ^ 2 + t41 ^ 2 + t42 ^ 2) + (t205 * t4 + t202 * t3 - t178 * (t202 * t26 + t205 * t27)) * t177;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
