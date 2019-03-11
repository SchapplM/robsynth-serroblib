% Calculate joint inertia matrix for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:22:54
% EndTime: 2019-03-09 09:23:03
% DurationCPUTime: 4.20s
% Computational Cost: add. (7445->474), mult. (14400->677), div. (0->0), fcn. (17172->10), ass. (0->219)
t289 = Icges(4,1) + Icges(5,1);
t288 = -Icges(4,4) + Icges(5,5);
t287 = Icges(5,4) + Icges(4,5);
t286 = Icges(4,2) + Icges(5,3);
t285 = Icges(4,6) - Icges(5,6);
t202 = sin(qJ(1));
t284 = -t202 / 0.2e1;
t271 = t202 / 0.2e1;
t205 = cos(qJ(1));
t283 = t205 / 0.2e1;
t198 = sin(pkin(10));
t199 = cos(pkin(10));
t204 = cos(qJ(2));
t253 = t202 * t204;
t165 = t198 * t253 + t199 * t205;
t166 = -t198 * t205 + t199 * t253;
t201 = sin(qJ(2));
t256 = t201 * t202;
t280 = t286 * t165 + t288 * t166 - t285 * t256;
t252 = t204 * t205;
t167 = t198 * t252 - t202 * t199;
t168 = t198 * t202 + t199 * t252;
t255 = t201 * t205;
t279 = t286 * t167 + t288 * t168 - t285 * t255;
t105 = Icges(4,5) * t166 - Icges(4,6) * t165 + Icges(4,3) * t256;
t107 = Icges(5,4) * t166 + Icges(5,2) * t256 + Icges(5,6) * t165;
t282 = t105 + t107;
t106 = Icges(4,5) * t168 - Icges(4,6) * t167 + Icges(4,3) * t255;
t108 = Icges(5,4) * t168 + Icges(5,2) * t255 + Icges(5,6) * t167;
t281 = t108 + t106;
t278 = t288 * t165 + t289 * t166 + t287 * t256;
t277 = t288 * t167 + t289 * t168 + t287 * t255;
t195 = t202 ^ 2;
t196 = t205 ^ 2;
t276 = 0.2e1 * t201;
t275 = m(4) / 0.2e1;
t274 = m(5) / 0.2e1;
t273 = m(6) / 0.2e1;
t272 = m(7) / 0.2e1;
t270 = t204 / 0.2e1;
t269 = -t205 / 0.2e1;
t268 = pkin(2) * t204;
t203 = cos(qJ(5));
t187 = pkin(5) * t203 + pkin(4);
t267 = -pkin(4) + t187;
t184 = pkin(8) * t256;
t200 = sin(qJ(5));
t260 = t165 * t200;
t241 = pkin(5) * t260;
t206 = -pkin(9) - pkin(8);
t254 = t201 * t206;
t197 = qJ(5) + qJ(6);
t188 = sin(t197);
t189 = cos(t197);
t115 = t165 * t189 - t166 * t188;
t116 = t165 * t188 + t166 * t189;
t216 = -t116 * rSges(7,1) - t115 * rSges(7,2);
t66 = -rSges(7,3) * t256 - t216;
t266 = t166 * t267 + t202 * t254 + t184 + t241 + t66;
t117 = t167 * t189 - t168 * t188;
t118 = t167 * t188 + t168 * t189;
t67 = t118 * rSges(7,1) + t117 * rSges(7,2) - rSges(7,3) * t255;
t135 = (-t188 * t199 + t189 * t198) * t201;
t136 = (t188 * t198 + t189 * t199) * t201;
t94 = rSges(7,1) * t136 + rSges(7,2) * t135 + rSges(7,3) * t204;
t50 = t204 * t67 + t94 * t255;
t265 = t165 * rSges(5,3);
t264 = t205 * rSges(3,3);
t258 = t198 * t200;
t122 = (-pkin(8) - t206) * t204 + (pkin(5) * t258 + t199 * t267) * t201;
t263 = -t122 - t94;
t262 = Icges(3,4) * t201;
t261 = Icges(3,4) * t204;
t259 = t167 * t200;
t257 = t198 * t201;
t127 = t167 * t203 - t168 * t200;
t128 = t168 * t203 + t259;
t251 = t128 * rSges(6,1) + t127 * rSges(6,2);
t175 = pkin(2) * t201 - qJ(3) * t204;
t250 = t204 * rSges(4,3) - (rSges(4,1) * t199 - rSges(4,2) * t198) * t201 - t175;
t247 = pkin(2) * t252 + qJ(3) * t255;
t249 = t195 * (qJ(3) * t201 + t268) + t205 * t247;
t248 = -(pkin(3) * t199 + qJ(4) * t198) * t201 - t175;
t246 = t205 * pkin(1) + t202 * pkin(7);
t157 = t165 * qJ(4);
t192 = t205 * pkin(7);
t245 = t192 - t157;
t244 = t195 + t196;
t90 = Icges(7,5) * t136 + Icges(7,6) * t135 + Icges(7,3) * t204;
t91 = Icges(7,4) * t136 + Icges(7,2) * t135 + Icges(7,6) * t204;
t92 = Icges(7,1) * t136 + Icges(7,4) * t135 + Icges(7,5) * t204;
t243 = t135 * t91 + t136 * t92 + t204 * t90;
t153 = (t198 * t203 - t199 * t200) * t201;
t154 = (t199 * t203 + t258) * t201;
t96 = Icges(6,5) * t154 + Icges(6,6) * t153 + Icges(6,3) * t204;
t97 = Icges(6,4) * t154 + Icges(6,2) * t153 + Icges(6,6) * t204;
t98 = Icges(6,1) * t154 + Icges(6,4) * t153 + Icges(6,5) * t204;
t242 = t153 * t97 + t154 * t98 + t204 * t96;
t125 = t165 * t203 - t166 * t200;
t126 = t166 * t203 + t260;
t70 = Icges(6,5) * t126 + Icges(6,6) * t125 - Icges(6,3) * t256;
t72 = Icges(6,4) * t126 + Icges(6,2) * t125 - Icges(6,6) * t256;
t74 = Icges(6,1) * t126 + Icges(6,4) * t125 - Icges(6,5) * t256;
t30 = t153 * t72 + t154 * t74 + t204 * t70;
t38 = t125 * t97 + t126 * t98 - t256 * t96;
t240 = -t30 / 0.2e1 - t38 / 0.2e1;
t71 = Icges(6,5) * t128 + Icges(6,6) * t127 - Icges(6,3) * t255;
t73 = Icges(6,4) * t128 + Icges(6,2) * t127 - Icges(6,6) * t255;
t75 = Icges(6,1) * t128 + Icges(6,4) * t127 - Icges(6,5) * t255;
t31 = t153 * t73 + t154 * t75 + t204 * t71;
t39 = t127 * t97 + t128 * t98 - t255 * t96;
t239 = -t31 / 0.2e1 - t39 / 0.2e1;
t238 = pkin(5) * t259 + t168 * t187 + t205 * t254;
t237 = t204 * rSges(5,2) - (rSges(5,1) * t199 + rSges(5,3) * t198) * t201 + t248;
t236 = t168 * rSges(5,1) + rSges(5,2) * t255 + t167 * rSges(5,3);
t235 = t168 * rSges(4,1) - t167 * rSges(4,2) + rSges(4,3) * t255;
t234 = -t201 * t199 * pkin(4) - t204 * pkin(8) + t248;
t233 = -pkin(1) - t268;
t232 = -t256 / 0.2e1;
t231 = -t255 / 0.2e1;
t138 = -Icges(5,6) * t204 + (Icges(5,5) * t199 + Icges(5,3) * t198) * t201;
t141 = -Icges(4,6) * t204 + (Icges(4,4) * t199 - Icges(4,2) * t198) * t201;
t230 = t138 / 0.2e1 - t141 / 0.2e1;
t142 = -Icges(5,4) * t204 + (Icges(5,1) * t199 + Icges(5,5) * t198) * t201;
t143 = -Icges(4,5) * t204 + (Icges(4,1) * t199 - Icges(4,4) * t198) * t201;
t229 = t142 / 0.2e1 + t143 / 0.2e1;
t228 = t168 * pkin(3) + t167 * qJ(4);
t60 = Icges(7,5) * t116 + Icges(7,6) * t115 - Icges(7,3) * t256;
t62 = Icges(7,4) * t116 + Icges(7,2) * t115 - Icges(7,6) * t256;
t64 = Icges(7,1) * t116 + Icges(7,4) * t115 - Icges(7,5) * t256;
t16 = t115 * t62 + t116 * t64 - t256 * t60;
t61 = Icges(7,5) * t118 + Icges(7,6) * t117 - Icges(7,3) * t255;
t63 = Icges(7,4) * t118 + Icges(7,2) * t117 - Icges(7,6) * t255;
t65 = Icges(7,1) * t118 + Icges(7,4) * t117 - Icges(7,5) * t255;
t17 = t115 * t63 + t116 * t65 - t256 * t61;
t10 = -t16 * t205 + t17 * t202;
t18 = t117 * t62 + t118 * t64 - t255 * t60;
t19 = t117 * t63 + t118 * t65 - t255 * t61;
t11 = -t18 * t205 + t19 * t202;
t28 = t135 * t62 + t136 * t64 + t204 * t60;
t29 = t135 * t63 + t136 * t65 + t204 * t61;
t36 = t115 * t91 + t116 * t92 - t256 * t90;
t3 = t36 * t204 + (-t16 * t202 - t17 * t205) * t201;
t37 = t117 * t91 + t118 * t92 - t255 * t90;
t4 = t37 * t204 + (-t18 * t202 - t19 * t205) * t201;
t227 = t10 * t232 + t11 * t231 + t3 * t269 + t4 * t271 + (t29 * t202 - t28 * t205) * t270;
t226 = t274 + t273 + t272;
t99 = rSges(6,1) * t154 + rSges(6,2) * t153 + rSges(6,3) * t204;
t225 = -t99 + t234;
t224 = t202 * (t166 * pkin(3) + t157) + t205 * t228 + t249;
t223 = t246 + t247;
t40 = t243 * t204;
t222 = t40 + (t28 + t36) * t232 + (t29 + t37) * t231;
t162 = t168 * pkin(4);
t221 = -pkin(8) * t255 + t162;
t220 = t234 + t263;
t219 = rSges(3,1) * t204 - rSges(3,2) * t201;
t218 = -t166 * rSges(4,1) + t165 * rSges(4,2);
t217 = -t126 * rSges(6,1) - t125 * rSges(6,2);
t215 = Icges(3,1) * t204 - t262;
t214 = -Icges(3,2) * t201 + t261;
t213 = Icges(3,5) * t204 - Icges(3,6) * t201;
t210 = rSges(3,1) * t252 - rSges(3,2) * t255 + t202 * rSges(3,3);
t209 = t202 * (t166 * pkin(4) - t184) + t205 * t221 + t224;
t7 = t204 * (t40 + (-t202 * t28 - t205 * t29) * t201);
t208 = t7 + (-t202 * t3 - t205 * t4) * t201;
t207 = t223 + t228;
t194 = t201 ^ 2;
t178 = rSges(2,1) * t205 - rSges(2,2) * t202;
t177 = -rSges(2,1) * t202 - rSges(2,2) * t205;
t176 = rSges(3,1) * t201 + rSges(3,2) * t204;
t172 = Icges(3,5) * t201 + Icges(3,6) * t204;
t148 = Icges(3,3) * t202 + t205 * t213;
t147 = -Icges(3,3) * t205 + t202 * t213;
t134 = t210 + t246;
t133 = t264 + t192 + (-pkin(1) - t219) * t202;
t130 = t250 * t205;
t129 = t250 * t202;
t100 = t205 * t210 + (t202 * t219 - t264) * t202;
t89 = t237 * t205;
t88 = t237 * t202;
t83 = t223 + t235;
t82 = t192 + ((-rSges(4,3) - qJ(3)) * t201 + t233) * t202 + t218;
t81 = -t221 + t238;
t77 = -rSges(6,3) * t255 + t251;
t76 = -rSges(6,3) * t256 - t217;
t69 = t225 * t205;
t68 = t225 * t202;
t58 = t207 + t236;
t57 = -t265 + (-rSges(5,1) - pkin(3)) * t166 + ((-rSges(5,2) - qJ(3)) * t201 + t233) * t202 + t245;
t56 = t67 * t256;
t55 = t202 * (rSges(4,3) * t256 - t218) + t205 * t235 + t249;
t54 = t204 * t77 + t255 * t99;
t53 = -t204 * t76 - t256 * t99;
t52 = t220 * t205;
t51 = t220 * t202;
t49 = -t204 * t66 - t256 * t94;
t48 = t162 + (-rSges(6,3) - pkin(8)) * t255 + t207 + t251;
t47 = t184 + (-pkin(3) - pkin(4)) * t166 + ((rSges(6,3) - qJ(3)) * t201 + t233) * t202 + t217 + t245;
t46 = t242 * t204;
t45 = (t202 * t77 - t205 * t76) * t201;
t44 = t207 + t67 + t238;
t43 = -t241 + (-pkin(3) - t187) * t166 + ((rSges(7,3) - qJ(3) - t206) * t201 + t233) * t202 + t216 + t245;
t42 = -t255 * t66 + t56;
t41 = t202 * (t166 * rSges(5,1) + rSges(5,2) * t256 + t265) + t205 * t236 + t224;
t33 = t122 * t255 + t204 * t81 + t50;
t32 = -t204 * t266 + t256 * t263;
t25 = t127 * t73 + t128 * t75 - t255 * t71;
t24 = t127 * t72 + t128 * t74 - t255 * t70;
t23 = t125 * t73 + t126 * t75 - t256 * t71;
t22 = t125 * t72 + t126 * t74 - t256 * t70;
t21 = t202 * t76 + t205 * t77 + t209;
t20 = t56 + (t202 * t81 - t205 * t266) * t201;
t15 = (t67 + t81) * t205 + t266 * t202 + t209;
t13 = t202 * t25 - t205 * t24;
t12 = t202 * t23 - t205 * t22;
t6 = t39 * t204 + (-t202 * t24 - t205 * t25) * t201;
t5 = t38 * t204 + (-t202 * t22 - t205 * t23) * t201;
t1 = [Icges(2,3) + (t262 + (t285 * t198 - t287 * t199) * t201 + (Icges(5,2) + Icges(3,2) + Icges(4,3)) * t204) * t204 + (Icges(3,1) * t201 + t261 + (t142 + t143) * t199 + (t138 - t141) * t198) * t201 + m(7) * (t43 ^ 2 + t44 ^ 2) + m(6) * (t47 ^ 2 + t48 ^ 2) + m(5) * (t57 ^ 2 + t58 ^ 2) + m(4) * (t82 ^ 2 + t83 ^ 2) + m(3) * (t133 ^ 2 + t134 ^ 2) + m(2) * (t177 ^ 2 + t178 ^ 2) + t242 + t243; (-t36 / 0.2e1 - t28 / 0.2e1 + t172 * t283 - t229 * t166 - t230 * t165 + t240) * t205 + (t37 / 0.2e1 + t29 / 0.2e1 + t172 * t271 + t229 * t168 + t230 * t167 - t239) * t202 + m(7) * (t43 * t52 + t44 * t51) + m(6) * (t47 * t69 + t48 * t68) + m(5) * (t57 * t89 + t58 * t88) + m(4) * (t129 * t83 + t130 * t82) + m(3) * (-t133 * t205 - t134 * t202) * t176 + ((t105 / 0.2e1 + Icges(3,6) * t283 + t214 * t284 + t107 / 0.2e1) * t205 + (-t108 / 0.2e1 - t106 / 0.2e1 + Icges(3,6) * t271 + t214 * t283) * t202) * t204 + ((Icges(3,5) * t202 + t279 * t198 + t277 * t199 + t205 * t215) * t271 + (-Icges(3,5) * t205 + t280 * t198 + t278 * t199 + t202 * t215) * t269) * t201; m(7) * (t15 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(6) * (t21 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(5) * (t41 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(4) * (t129 ^ 2 + t130 ^ 2 + t55 ^ 2) + m(3) * (t176 ^ 2 * t244 + t100 ^ 2) + (-t196 * t147 - t10 - t12 + (t280 * t165 + t278 * t166 + t282 * t256) * t205) * t205 + (t11 + t13 + t195 * t148 + (t279 * t167 + t277 * t168 + t281 * t255) * t202 + (-t202 * t147 + t205 * t148 - t279 * t165 - t277 * t166 - t280 * t167 - t278 * t168 - t282 * t255 - t281 * t256) * t205) * t202; ((t202 * t44 + t205 * t43) * t272 + (t202 * t48 + t205 * t47) * t273 + (t202 * t58 + t205 * t57) * t274 + (t202 * t83 + t205 * t82) * t275) * t276; m(7) * (-t204 * t15 + (t202 * t51 + t205 * t52) * t201) + m(6) * (-t204 * t21 + (t202 * t68 + t205 * t69) * t201) + m(5) * (-t204 * t41 + (t202 * t88 + t205 * t89) * t201) + m(4) * (-t204 * t55 + (t129 * t202 + t130 * t205) * t201); 0.2e1 * (t275 + t226) * (t194 * t244 + t204 ^ 2); m(7) * (t165 * t44 + t167 * t43) + m(6) * (t165 * t48 + t167 * t47) + m(5) * (t165 * t58 + t167 * t57); m(7) * (t15 * t257 + t165 * t51 + t167 * t52) + m(6) * (t165 * t68 + t167 * t69 + t21 * t257) + m(5) * (t165 * t88 + t167 * t89 + t257 * t41); t226 * (t165 * t202 + t167 * t205 - t198 * t204) * t276; 0.2e1 * t226 * (t194 * t198 ^ 2 + t165 ^ 2 + t167 ^ 2); t46 + m(7) * (t32 * t43 + t33 * t44) + m(6) * (t47 * t53 + t48 * t54) + (t202 * t240 + t205 * t239) * t201 + t222; t5 * t269 + (t31 * t202 - t30 * t205) * t270 + t6 * t271 + (t12 * t284 + t13 * t269) * t201 + m(7) * (t15 * t20 + t32 * t52 + t33 * t51) + m(6) * (t21 * t45 + t53 * t69 + t54 * t68) + t227; m(6) * (-t45 * t204 + (t202 * t54 + t205 * t53) * t201) + m(7) * (-t20 * t204 + (t202 * t33 + t205 * t32) * t201); m(6) * (t165 * t54 + t167 * t53 + t257 * t45) + m(7) * (t165 * t33 + t167 * t32 + t20 * t257); t204 * t46 + t7 + m(7) * (t20 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t45 ^ 2 + t53 ^ 2 + t54 ^ 2) + ((-t204 * t31 - t4 - t6) * t205 + (-t204 * t30 - t3 - t5) * t202) * t201; m(7) * (t43 * t49 + t44 * t50) + t222; m(7) * (t15 * t42 + t49 * t52 + t50 * t51) + t227; m(7) * (-t42 * t204 + (t202 * t50 + t205 * t49) * t201); m(7) * (t165 * t50 + t167 * t49 + t257 * t42); m(7) * (t20 * t42 + t32 * t49 + t33 * t50) + t208; m(7) * (t42 ^ 2 + t49 ^ 2 + t50 ^ 2) + t208;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
