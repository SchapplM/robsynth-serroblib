% Calculate joint inertia matrix for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR9_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR9_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR9_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:22:16
% EndTime: 2019-12-31 21:22:25
% DurationCPUTime: 3.20s
% Computational Cost: add. (7263->415), mult. (8808->592), div. (0->0), fcn. (9497->10), ass. (0->218)
t202 = sin(qJ(2));
t283 = Icges(3,5) * t202;
t282 = t283 / 0.2e1;
t196 = qJ(3) + pkin(9);
t187 = sin(t196);
t188 = cos(t196);
t205 = cos(qJ(2));
t127 = -Icges(5,3) * t205 + (Icges(5,5) * t188 - Icges(5,6) * t187) * t202;
t201 = sin(qJ(3));
t204 = cos(qJ(3));
t143 = -Icges(4,3) * t205 + (Icges(4,5) * t204 - Icges(4,6) * t201) * t202;
t281 = -t127 - t143;
t128 = -Icges(5,6) * t205 + (Icges(5,4) * t188 - Icges(5,2) * t187) * t202;
t146 = -Icges(4,6) * t205 + (Icges(4,4) * t204 - Icges(4,2) * t201) * t202;
t280 = -t187 * t128 - t201 * t146;
t203 = sin(qJ(1));
t206 = cos(qJ(1));
t244 = t206 * t188;
t250 = t203 * t205;
t139 = -t187 * t250 - t244;
t245 = t206 * t187;
t140 = t188 * t250 - t245;
t253 = t202 * t203;
t89 = Icges(5,5) * t140 + Icges(5,6) * t139 + Icges(5,3) * t253;
t91 = Icges(5,4) * t140 + Icges(5,2) * t139 + Icges(5,6) * t253;
t93 = Icges(5,1) * t140 + Icges(5,4) * t139 + Icges(5,5) * t253;
t31 = t139 * t91 + t140 * t93 + t253 * t89;
t141 = t188 * t203 - t205 * t245;
t142 = t187 * t203 + t205 * t244;
t252 = t202 * t206;
t90 = Icges(5,5) * t142 + Icges(5,6) * t141 + Icges(5,3) * t252;
t92 = Icges(5,4) * t142 + Icges(5,2) * t141 + Icges(5,6) * t252;
t94 = Icges(5,1) * t142 + Icges(5,4) * t141 + Icges(5,5) * t252;
t32 = t139 * t92 + t140 * t94 + t253 * t90;
t242 = t206 * t204;
t159 = -t201 * t250 - t242;
t243 = t206 * t201;
t160 = t204 * t250 - t243;
t102 = Icges(4,5) * t160 + Icges(4,6) * t159 + Icges(4,3) * t253;
t104 = Icges(4,4) * t160 + Icges(4,2) * t159 + Icges(4,6) * t253;
t106 = Icges(4,1) * t160 + Icges(4,4) * t159 + Icges(4,5) * t253;
t41 = t102 * t253 + t104 * t159 + t106 * t160;
t161 = t203 * t204 - t205 * t243;
t251 = t203 * t201;
t162 = t205 * t242 + t251;
t103 = Icges(4,5) * t162 + Icges(4,6) * t161 + Icges(4,3) * t252;
t105 = Icges(4,4) * t162 + Icges(4,2) * t161 + Icges(4,6) * t252;
t107 = Icges(4,1) * t162 + Icges(4,4) * t161 + Icges(4,5) * t252;
t42 = t103 * t253 + t105 * t159 + t107 * t160;
t129 = -Icges(5,5) * t205 + (Icges(5,1) * t188 - Icges(5,4) * t187) * t202;
t55 = t127 * t253 + t128 * t139 + t129 * t140;
t149 = -Icges(4,5) * t205 + (Icges(4,1) * t204 - Icges(4,4) * t201) * t202;
t62 = t143 * t253 + t146 * t159 + t149 * t160;
t279 = (-t62 - t55) * t205 + ((t32 + t42) * t206 + (t31 + t41) * t203) * t202;
t33 = t141 * t91 + t142 * t93 + t252 * t89;
t34 = t141 * t92 + t142 * t94 + t252 * t90;
t43 = t102 * t252 + t104 * t161 + t106 * t162;
t44 = t103 * t252 + t105 * t161 + t107 * t162;
t56 = t127 * t252 + t128 * t141 + t129 * t142;
t63 = t143 * t252 + t146 * t161 + t149 * t162;
t278 = (-t63 - t56) * t205 + ((t34 + t44) * t206 + (t33 + t43) * t203) * t202;
t39 = -t205 * t89 + (-t187 * t91 + t188 * t93) * t202;
t47 = -t205 * t102 + (-t104 * t201 + t106 * t204) * t202;
t277 = -t39 - t47;
t40 = -t205 * t90 + (-t187 * t92 + t188 * t94) * t202;
t48 = -t205 * t103 + (-t105 * t201 + t107 * t204) * t202;
t276 = t40 + t48;
t186 = pkin(3) * t204 + pkin(2);
t164 = pkin(4) * t188 + t186;
t166 = pkin(3) * t201 + pkin(4) * t187;
t249 = t205 * t206;
t189 = qJ(5) + t196;
t185 = cos(t189);
t184 = sin(t189);
t247 = t206 * t184;
t135 = t185 * t203 - t205 * t247;
t246 = t206 * t185;
t136 = t184 * t203 + t205 * t246;
t88 = rSges(6,1) * t136 + rSges(6,2) * t135 + rSges(6,3) * t252;
t275 = t164 * t249 + t166 * t203 + t88;
t274 = (t129 * t188 + t149 * t204) * t202;
t198 = t203 ^ 2;
t199 = t206 ^ 2;
t273 = m(5) / 0.2e1;
t272 = m(6) / 0.2e1;
t133 = -t184 * t250 - t246;
t134 = t185 * t250 - t247;
t81 = Icges(6,5) * t134 + Icges(6,6) * t133 + Icges(6,3) * t253;
t83 = Icges(6,4) * t134 + Icges(6,2) * t133 + Icges(6,6) * t253;
t85 = Icges(6,1) * t134 + Icges(6,4) * t133 + Icges(6,5) * t253;
t26 = t133 * t83 + t134 * t85 + t253 * t81;
t82 = Icges(6,5) * t136 + Icges(6,6) * t135 + Icges(6,3) * t252;
t84 = Icges(6,4) * t136 + Icges(6,2) * t135 + Icges(6,6) * t252;
t86 = Icges(6,1) * t136 + Icges(6,4) * t135 + Icges(6,5) * t252;
t27 = t133 * t84 + t134 * t86 + t253 * t82;
t121 = -Icges(6,3) * t205 + (Icges(6,5) * t185 - Icges(6,6) * t184) * t202;
t122 = -Icges(6,6) * t205 + (Icges(6,4) * t185 - Icges(6,2) * t184) * t202;
t123 = -Icges(6,5) * t205 + (Icges(6,1) * t185 - Icges(6,4) * t184) * t202;
t51 = t121 * t253 + t122 * t133 + t123 * t134;
t5 = -t51 * t205 + (t203 * t26 + t206 * t27) * t202;
t28 = t135 * t83 + t136 * t85 + t252 * t81;
t29 = t135 * t84 + t136 * t86 + t252 * t82;
t52 = t121 * t252 + t122 * t135 + t123 * t136;
t6 = -t52 * t205 + (t203 * t28 + t206 * t29) * t202;
t271 = t252 * t6 + t253 * t5;
t270 = t203 / 0.2e1;
t269 = -t205 / 0.2e1;
t268 = -t206 / 0.2e1;
t171 = rSges(3,1) * t202 + rSges(3,2) * t205;
t267 = m(3) * t171;
t266 = pkin(2) * t205;
t265 = pkin(7) * t202;
t264 = -pkin(2) + t186;
t200 = -qJ(4) - pkin(7);
t263 = t202 * t280 + t205 * t281 + t274;
t195 = -pkin(8) + t200;
t231 = t195 - t200;
t235 = -pkin(3) * t251 - t186 * t249;
t262 = -t231 * t252 + t235 + t275;
t261 = t206 * rSges(3,3);
t113 = t202 * t185 * t123;
t256 = t184 * t122;
t61 = -t205 * t121 - t202 * t256 + t113;
t260 = t61 * t205;
t209 = -t200 * t252 - t235;
t233 = pkin(2) * t249 + pkin(7) * t252;
t109 = t209 - t233;
t96 = rSges(5,1) * t142 + rSges(5,2) * t141 + rSges(5,3) * t252;
t259 = -t109 - t96;
t124 = -t205 * rSges(6,3) + (rSges(6,1) * t185 - rSges(6,2) * t184) * t202;
t216 = -t134 * rSges(6,1) - t133 * rSges(6,2);
t87 = rSges(6,3) * t253 - t216;
t66 = t124 * t253 + t205 * t87;
t257 = Icges(3,4) * t205;
t248 = t206 * t166;
t234 = pkin(3) * t243 + t200 * t253;
t108 = (t205 * t264 - t265) * t203 - t234;
t126 = (pkin(7) + t200) * t205 + t264 * t202;
t241 = t108 * t205 + t126 * t253;
t132 = -t205 * rSges(5,3) + (rSges(5,1) * t188 - rSges(5,2) * t187) * t202;
t240 = -t126 - t132;
t152 = -t205 * rSges(4,3) + (rSges(4,1) * t204 - rSges(4,2) * t201) * t202;
t174 = t202 * pkin(2) - t205 * pkin(7);
t239 = -t152 - t174;
t237 = t198 * (t265 + t266) + t206 * t233;
t236 = t164 - t186;
t232 = t206 * pkin(1) + pkin(6) * t203;
t230 = t198 + t199;
t229 = -t109 - t262;
t112 = t202 * t236 + t205 * t231;
t228 = -t112 - t124 - t126;
t227 = -t174 + t240;
t111 = rSges(4,1) * t162 + rSges(4,2) * t161 + rSges(4,3) * t252;
t226 = t253 / 0.2e1;
t225 = t252 / 0.2e1;
t37 = -t205 * t81 + (-t184 * t83 + t185 * t85) * t202;
t38 = -t205 * t82 + (-t184 * t84 + t185 * t86) * t202;
t224 = (t37 + t51) * t226 + (t38 + t52) * t225;
t9 = -t260 + (t203 * t37 + t206 * t38) * t202;
t223 = -t205 * t9 + t271;
t222 = t108 * t203 + t109 * t206 + t237;
t221 = -t174 + t228;
t12 = t203 * t27 - t206 * t26;
t13 = t203 * t29 - t206 * t28;
t220 = t12 * t226 + t13 * t225 + t5 * t268 + t6 * t270 + (t38 * t203 - t37 * t206) * t269;
t219 = rSges(3,1) * t205 - rSges(3,2) * t202;
t218 = -t160 * rSges(4,1) - t159 * rSges(4,2);
t217 = -t140 * rSges(5,1) - t139 * rSges(5,2);
t214 = -Icges(3,2) * t202 + t257;
t213 = Icges(3,5) * t205 - Icges(3,6) * t202;
t210 = rSges(3,1) * t249 - rSges(3,2) * t252 + rSges(3,3) * t203;
t208 = t48 / 0.2e1 + t40 / 0.2e1 + t63 / 0.2e1 + t56 / 0.2e1;
t207 = t62 / 0.2e1 + t55 / 0.2e1 + t47 / 0.2e1 + t39 / 0.2e1;
t193 = t206 * pkin(6);
t173 = rSges(2,1) * t206 - rSges(2,2) * t203;
t172 = -rSges(2,1) * t203 - rSges(2,2) * t206;
t168 = Icges(3,6) * t205 + t283;
t145 = Icges(3,3) * t203 + t206 * t213;
t144 = -Icges(3,3) * t206 + t203 * t213;
t120 = t210 + t232;
t119 = t261 + t193 + (-pkin(1) - t219) * t203;
t115 = t239 * t206;
t114 = t239 * t203;
t110 = rSges(4,3) * t253 - t218;
t98 = t206 * t210 + (t203 * t219 - t261) * t203;
t97 = t108 * t252;
t95 = rSges(5,3) * t253 - t217;
t79 = t87 * t252;
t78 = t111 + t232 + t233;
t77 = t193 + (-t266 - pkin(1) + (-rSges(4,3) - pkin(7)) * t202) * t203 + t218;
t76 = t227 * t206;
t75 = t227 * t203;
t73 = -t248 + (-t195 * t202 + t205 * t236) * t203 + t234;
t72 = -t111 * t205 - t152 * t252;
t71 = t110 * t205 + t152 * t253;
t69 = t209 + t96 + t232;
t68 = t193 + (-rSges(5,3) * t202 - t186 * t205 - pkin(1)) * t203 + t217 + t234;
t67 = -t124 * t252 - t205 * t88;
t64 = (t110 * t206 - t111 * t203) * t202;
t60 = -t195 * t252 + t232 + t275;
t59 = t248 + t193 + (-t164 * t205 - pkin(1) + (-rSges(6,3) + t195) * t202) * t203 + t216;
t58 = t221 * t206;
t57 = t221 * t203;
t54 = -t253 * t88 + t79;
t53 = t110 * t203 + t111 * t206 + t237;
t46 = t205 * t259 + t240 * t252;
t45 = t132 * t253 + t205 * t95 + t241;
t30 = t97 + (t203 * t259 + t206 * t95) * t202;
t25 = t203 * t95 + t206 * t96 + t222;
t24 = t205 * t229 + t228 * t252;
t23 = t112 * t253 + t205 * t73 + t241 + t66;
t22 = t203 * t44 - t206 * t43;
t21 = t203 * t42 - t206 * t41;
t20 = t79 + t97 + (t203 * t229 + t206 * t73) * t202;
t19 = t262 * t206 + (t73 + t87) * t203 + t222;
t17 = t203 * t34 - t206 * t33;
t16 = t203 * t32 - t206 * t31;
t1 = [Icges(2,3) + t113 + (Icges(3,4) * t202 + Icges(3,2) * t205 - t121 + t281) * t205 + (Icges(3,1) * t202 - t256 + t257 + t280) * t202 + m(6) * (t59 ^ 2 + t60 ^ 2) + m(5) * (t68 ^ 2 + t69 ^ 2) + m(4) * (t77 ^ 2 + t78 ^ 2) + m(3) * (t119 ^ 2 + t120 ^ 2) + m(2) * (t172 ^ 2 + t173 ^ 2) + t274; m(6) * (t57 * t60 + t58 * t59) + m(5) * (t68 * t76 + t69 * t75) + m(4) * (t114 * t78 + t115 * t77) + (-t37 / 0.2e1 - t51 / 0.2e1 + t203 * t214 * t269 - t119 * t267 - t207 + (-Icges(3,6) * t269 + t282 + t168 / 0.2e1) * t206) * t206 + (t52 / 0.2e1 + t38 / 0.2e1 + t205 * (Icges(3,6) * t203 + t206 * t214) / 0.2e1 + t203 * t282 - t120 * t267 + t168 * t270 + t208) * t203; m(6) * (t19 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(5) * (t25 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(4) * (t114 ^ 2 + t115 ^ 2 + t53 ^ 2) + m(3) * (t171 ^ 2 * t230 + t98 ^ 2) + (-t199 * t144 - t12 - t16 - t21) * t206 + (t198 * t145 + t13 + t17 + t22 + (-t203 * t144 + t206 * t145) * t206) * t203; (-t61 - t263) * t205 + m(6) * (t23 * t59 + t24 * t60) + m(5) * (t45 * t68 + t46 * t69) + m(4) * (t71 * t77 + t72 * t78) + (t203 * t207 + t206 * t208) * t202 + t224; m(6) * (t19 * t20 + t23 * t58 + t24 * t57) + m(5) * (t25 * t30 + t45 * t76 + t46 * t75) + m(4) * (t114 * t72 + t115 * t71 + t53 * t64) + ((t17 / 0.2e1 + t22 / 0.2e1) * t206 + (t16 / 0.2e1 + t21 / 0.2e1) * t203) * t202 + t220 + t278 * t270 + (t203 * t276 + t206 * t277) * t269 + t279 * t268; m(6) * (t20 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(5) * (t30 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(4) * (t64 ^ 2 + t71 ^ 2 + t72 ^ 2) + (t263 * t205 - t9) * t205 + ((-t205 * t276 + t278) * t206 + (t205 * t277 + t279) * t203) * t202 + t271; 0.2e1 * ((t203 * t60 + t206 * t59) * t272 + (t203 * t69 + t206 * t68) * t273) * t202; m(6) * (-t205 * t19 + (t203 * t57 + t206 * t58) * t202) + m(5) * (-t205 * t25 + (t203 * t75 + t206 * t76) * t202); m(6) * (-t205 * t20 + (t203 * t24 + t206 * t23) * t202) + m(5) * (-t205 * t30 + (t203 * t46 + t206 * t45) * t202); 0.2e1 * (t273 + t272) * (t202 ^ 2 * t230 + t205 ^ 2); -t260 + m(6) * (t59 * t66 + t60 * t67) + t224; m(6) * (t19 * t54 + t57 * t67 + t58 * t66) + t220; m(6) * (t20 * t54 + t23 * t66 + t24 * t67) + t223; m(6) * (-t54 * t205 + (t203 * t67 + t206 * t66) * t202); m(6) * (t54 ^ 2 + t66 ^ 2 + t67 ^ 2) + t223;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
