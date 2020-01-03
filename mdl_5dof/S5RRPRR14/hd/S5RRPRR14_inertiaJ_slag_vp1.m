% Calculate joint inertia matrix for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR14_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR14_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR14_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:35:48
% EndTime: 2019-12-31 20:35:58
% DurationCPUTime: 3.74s
% Computational Cost: add. (13078->462), mult. (22250->660), div. (0->0), fcn. (27945->12), ass. (0->213)
t201 = sin(pkin(10));
t203 = cos(pkin(10));
t204 = cos(pkin(5));
t202 = sin(pkin(5));
t207 = sin(qJ(2));
t246 = t202 * t207;
t181 = -t201 * t246 + t204 * t203;
t242 = t204 * t201;
t182 = t203 * t246 + t242;
t210 = cos(qJ(2));
t244 = t202 * t210;
t129 = Icges(4,4) * t182 + Icges(4,2) * t181 - Icges(4,6) * t244;
t130 = Icges(4,1) * t182 + Icges(4,4) * t181 - Icges(4,5) * t244;
t168 = Icges(3,3) * t204 + (Icges(3,5) * t207 + Icges(3,6) * t210) * t202;
t169 = Icges(3,6) * t204 + (Icges(3,4) * t207 + Icges(3,2) * t210) * t202;
t170 = Icges(3,5) * t204 + (Icges(3,1) * t207 + Icges(3,4) * t210) * t202;
t260 = t181 * t129 + t182 * t130 + t204 * t168 + t169 * t244 + t170 * t246;
t128 = Icges(4,5) * t182 + Icges(4,6) * t181 - Icges(4,3) * t244;
t259 = (-t128 * t244 + t260) * t204;
t211 = cos(qJ(1));
t239 = t211 * t207;
t208 = sin(qJ(1));
t240 = t208 * t210;
t184 = t204 * t239 + t240;
t231 = pkin(10) + qJ(4);
t199 = sin(t231);
t223 = cos(t231);
t217 = t202 * t223;
t154 = t184 * t199 + t211 * t217;
t238 = t211 * t210;
t241 = t208 * t207;
t186 = -t204 * t241 + t238;
t156 = t186 * t199 - t208 * t217;
t171 = t199 * t246 - t204 * t223;
t245 = t202 * t208;
t157 = t186 * t223 + t199 * t245;
t185 = t204 * t240 + t239;
t206 = sin(qJ(5));
t209 = cos(qJ(5));
t120 = -t157 * t206 + t185 * t209;
t121 = t157 * t209 + t185 * t206;
t243 = t202 * t211;
t155 = t184 * t223 - t199 * t243;
t183 = -t204 * t238 + t241;
t118 = -t155 * t206 + t183 * t209;
t119 = t155 * t209 + t183 * t206;
t62 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t154;
t64 = Icges(6,4) * t119 + Icges(6,2) * t118 + Icges(6,6) * t154;
t66 = Icges(6,1) * t119 + Icges(6,4) * t118 + Icges(6,5) * t154;
t19 = t120 * t64 + t121 * t66 + t156 * t62;
t63 = Icges(6,5) * t121 + Icges(6,6) * t120 + Icges(6,3) * t156;
t65 = Icges(6,4) * t121 + Icges(6,2) * t120 + Icges(6,6) * t156;
t67 = Icges(6,1) * t121 + Icges(6,4) * t120 + Icges(6,5) * t156;
t20 = t120 * t65 + t121 * t67 + t156 * t63;
t172 = t204 * t199 + t207 * t217;
t152 = -t172 * t206 - t209 * t244;
t153 = t172 * t209 - t206 * t244;
t82 = Icges(6,5) * t153 + Icges(6,6) * t152 + Icges(6,3) * t171;
t83 = Icges(6,4) * t153 + Icges(6,2) * t152 + Icges(6,6) * t171;
t84 = Icges(6,1) * t153 + Icges(6,4) * t152 + Icges(6,5) * t171;
t28 = t120 * t83 + t121 * t84 + t156 * t82;
t2 = t19 * t154 + t20 * t156 + t28 * t171;
t258 = t2 / 0.2e1;
t257 = t154 / 0.2e1;
t256 = t156 / 0.2e1;
t255 = t171 / 0.2e1;
t254 = t155 * pkin(4);
t198 = t203 * pkin(3) + pkin(2);
t253 = -pkin(2) + t198;
t218 = -t119 * rSges(6,1) - t118 * rSges(6,2);
t68 = t154 * rSges(6,3) - t218;
t252 = t154 * pkin(9) + t254 + t68;
t69 = t121 * rSges(6,1) + t120 * rSges(6,2) + t156 * rSges(6,3);
t251 = t157 * pkin(4) + t156 * pkin(9) + t69;
t85 = t153 * rSges(6,1) + t152 * rSges(6,2) + t171 * rSges(6,3);
t250 = t172 * pkin(4) + t171 * pkin(9) + t85;
t147 = t186 * pkin(2) + t185 * qJ(3);
t145 = t204 * t147;
t205 = -pkin(8) - qJ(3);
t228 = t201 * t245;
t225 = pkin(3) * t228 - t185 * t205 + t186 * t198;
t98 = -t147 + t225;
t249 = t204 * t98 + t145;
t175 = t183 * qJ(3);
t146 = t184 * pkin(2) + t175;
t227 = t201 * t243;
t190 = pkin(3) * t227;
t247 = t183 * t205;
t97 = t184 * t253 - t175 - t190 - t247;
t248 = -t146 - t97;
t160 = -t184 * t201 - t203 * t243;
t161 = t184 * t203 - t227;
t106 = t161 * rSges(4,1) + t160 * rSges(4,2) + t183 * rSges(4,3);
t237 = -t106 - t146;
t125 = Icges(5,4) * t172 - Icges(5,2) * t171 - Icges(5,6) * t244;
t126 = Icges(5,1) * t172 - Icges(5,4) * t171 - Icges(5,5) * t244;
t236 = -t171 * t125 + t172 * t126;
t187 = (pkin(2) * t207 - qJ(3) * t210) * t202;
t234 = -pkin(3) * t242 - ((qJ(3) + t205) * t210 + t253 * t207) * t202 - t187;
t233 = t146 * t245 + t147 * t243;
t232 = t211 * pkin(1) + pkin(7) * t245;
t33 = t152 * t83 + t153 * t84 + t171 * t82;
t23 = t152 * t64 + t153 * t66 + t171 * t62;
t27 = t118 * t83 + t119 * t84 + t154 * t82;
t230 = t23 / 0.2e1 + t27 / 0.2e1;
t24 = t152 * t65 + t153 * t67 + t171 * t63;
t229 = t24 / 0.2e1 + t28 / 0.2e1;
t96 = t157 * rSges(5,1) - t156 * rSges(5,2) + t185 * rSges(5,3);
t162 = -t186 * t201 + t203 * t245;
t163 = t186 * t203 + t228;
t107 = t163 * rSges(4,1) + t162 * rSges(4,2) + t185 * rSges(4,3);
t141 = t186 * rSges(3,1) - t185 * rSges(3,2) + rSges(3,3) * t245;
t224 = -t208 * pkin(1) + pkin(7) * t243;
t222 = t202 * (-t182 * rSges(4,1) - t181 * rSges(4,2) + rSges(4,3) * t244 - t187);
t221 = t98 * t243 + t97 * t245 + t233;
t127 = t172 * rSges(5,1) - t171 * rSges(5,2) - rSges(5,3) * t244;
t220 = t202 * (-t127 + t234);
t219 = -t155 * rSges(5,1) + t154 * rSges(5,2);
t216 = t225 + t232;
t215 = t202 * (t234 - t250);
t89 = Icges(5,5) * t155 - Icges(5,6) * t154 + Icges(5,3) * t183;
t91 = Icges(5,4) * t155 - Icges(5,2) * t154 + Icges(5,6) * t183;
t93 = Icges(5,1) * t155 - Icges(5,4) * t154 + Icges(5,5) * t183;
t41 = -t171 * t91 + t172 * t93 - t244 * t89;
t124 = Icges(5,5) * t172 - Icges(5,6) * t171 - Icges(5,3) * t244;
t50 = t183 * t124 - t154 * t125 + t155 * t126;
t214 = t41 / 0.2e1 + t50 / 0.2e1 + t230;
t90 = Icges(5,5) * t157 - Icges(5,6) * t156 + Icges(5,3) * t185;
t92 = Icges(5,4) * t157 - Icges(5,2) * t156 + Icges(5,6) * t185;
t94 = Icges(5,1) * t157 - Icges(5,4) * t156 + Icges(5,5) * t185;
t42 = -t171 * t92 + t172 * t94 - t244 * t90;
t51 = t185 * t124 - t156 * t125 + t157 * t126;
t213 = t42 / 0.2e1 + t51 / 0.2e1 + t229;
t212 = -t184 * t198 + t190 + t224;
t140 = t184 * rSges(3,1) - t183 * rSges(3,2) - rSges(3,3) * t243;
t192 = t211 * rSges(2,1) - t208 * rSges(2,2);
t191 = -t208 * rSges(2,1) - t211 * rSges(2,2);
t173 = t204 * rSges(3,3) + (rSges(3,1) * t207 + rSges(3,2) * t210) * t202;
t139 = Icges(3,1) * t186 - Icges(3,4) * t185 + Icges(3,5) * t245;
t138 = Icges(3,1) * t184 - Icges(3,4) * t183 - Icges(3,5) * t243;
t137 = Icges(3,4) * t186 - Icges(3,2) * t185 + Icges(3,6) * t245;
t136 = Icges(3,4) * t184 - Icges(3,2) * t183 - Icges(3,6) * t243;
t135 = Icges(3,5) * t186 - Icges(3,6) * t185 + Icges(3,3) * t245;
t134 = Icges(3,5) * t184 - Icges(3,6) * t183 - Icges(3,3) * t243;
t123 = t141 + t232;
t122 = -t140 + t224;
t111 = -t204 * t140 - t173 * t243;
t110 = t204 * t141 - t173 * t245;
t105 = Icges(4,1) * t163 + Icges(4,4) * t162 + Icges(4,5) * t185;
t104 = Icges(4,1) * t161 + Icges(4,4) * t160 + Icges(4,5) * t183;
t103 = Icges(4,4) * t163 + Icges(4,2) * t162 + Icges(4,6) * t185;
t102 = Icges(4,4) * t161 + Icges(4,2) * t160 + Icges(4,6) * t183;
t101 = Icges(4,5) * t163 + Icges(4,6) * t162 + Icges(4,3) * t185;
t100 = Icges(4,5) * t161 + Icges(4,6) * t160 + Icges(4,3) * t183;
t95 = t183 * rSges(5,3) - t219;
t81 = (t140 * t208 + t141 * t211) * t202;
t80 = t168 * t245 - t185 * t169 + t186 * t170;
t79 = -t168 * t243 - t183 * t169 + t184 * t170;
t77 = t147 + t107 + t232;
t76 = t224 + t237;
t73 = t204 * t135 + (t137 * t210 + t139 * t207) * t202;
t72 = t204 * t134 + (t136 * t210 + t138 * t207) * t202;
t71 = t216 + t96;
t70 = (-rSges(5,3) + t205) * t183 + t212 + t219;
t61 = -t185 * t127 - t244 * t96;
t60 = t183 * t127 + t244 * t95;
t58 = t204 * t237 + t211 * t222;
t57 = t204 * t107 + t208 * t222 + t145;
t56 = -t124 * t244 + t236;
t55 = -t183 * t96 + t185 * t95;
t54 = t56 * t204;
t53 = t185 * t128 + t162 * t129 + t163 * t130;
t52 = t183 * t128 + t160 * t129 + t161 * t130;
t49 = (t106 * t208 + t107 * t211) * t202 + t233;
t48 = -t101 * t244 + t181 * t103 + t182 * t105;
t47 = -t100 * t244 + t181 * t102 + t182 * t104;
t46 = t216 + t251;
t45 = -t254 + t247 + (-rSges(6,3) - pkin(9)) * t154 + t212 + t218;
t44 = -t156 * t85 + t171 * t69;
t43 = t154 * t85 - t171 * t68;
t40 = -t156 * t92 + t157 * t94 + t185 * t90;
t39 = -t156 * t91 + t157 * t93 + t185 * t89;
t38 = -t154 * t92 + t155 * t94 + t183 * t90;
t37 = -t154 * t91 + t155 * t93 + t183 * t89;
t36 = (-t95 + t248) * t204 + t211 * t220;
t35 = t204 * t96 + t208 * t220 + t249;
t34 = -t154 * t69 + t156 * t68;
t32 = t33 * t204;
t31 = t33 * t171;
t30 = -t185 * t250 - t244 * t251;
t29 = t183 * t250 + t244 * t252;
t26 = (t208 * t95 + t211 * t96) * t202 + t221;
t25 = -t183 * t251 + t185 * t252;
t22 = (t248 - t252) * t204 + t211 * t215;
t21 = t204 * t251 + t208 * t215 + t249;
t18 = t118 * t65 + t119 * t67 + t154 * t63;
t17 = t118 * t64 + t119 * t66 + t154 * t62;
t16 = (t208 * t252 + t211 * t251) * t202 + t221;
t15 = t54 + (t42 * t208 - t41 * t211) * t202;
t14 = t41 * t183 + t42 * t185 - t244 * t56;
t13 = t51 * t204 + (t208 * t40 - t211 * t39) * t202;
t12 = t50 * t204 + (t208 * t38 - t211 * t37) * t202;
t11 = t39 * t183 + t40 * t185 - t244 * t51;
t10 = t37 * t183 + t38 * t185 - t244 * t50;
t9 = t32 + (t24 * t208 - t23 * t211) * t202;
t8 = t23 * t183 + t24 * t185 - t244 * t33;
t7 = t23 * t154 + t24 * t156 + t31;
t6 = t28 * t204 + (-t19 * t211 + t20 * t208) * t202;
t5 = t27 * t204 + (-t17 * t211 + t18 * t208) * t202;
t4 = t19 * t183 + t20 * t185 - t244 * t28;
t3 = t17 * t183 + t18 * t185 - t244 * t27;
t1 = t17 * t154 + t18 * t156 + t27 * t171;
t59 = [Icges(2,3) + (-t124 - t128) * t244 + m(6) * (t45 ^ 2 + t46 ^ 2) + m(5) * (t70 ^ 2 + t71 ^ 2) + m(4) * (t76 ^ 2 + t77 ^ 2) + m(3) * (t122 ^ 2 + t123 ^ 2) + m(2) * (t191 ^ 2 + t192 ^ 2) + t33 + t236 + t260; t32 + t54 + m(6) * (t21 * t46 + t22 * t45) + m(5) * (t35 * t71 + t36 * t70) + m(4) * (t57 * t77 + t58 * t76) + m(3) * (t110 * t123 + t111 * t122) + ((-t79 / 0.2e1 - t52 / 0.2e1 - t72 / 0.2e1 - t47 / 0.2e1 - t214) * t211 + (t80 / 0.2e1 + t53 / 0.2e1 + t73 / 0.2e1 + t48 / 0.2e1 + t213) * t208) * t202 + t259; m(6) * (t16 ^ 2 + t21 ^ 2 + t22 ^ 2) + m(5) * (t26 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(4) * (t49 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(3) * (t110 ^ 2 + t111 ^ 2 + t81 ^ 2) + (t6 + t13 + ((t185 * t101 + t162 * t103 + t163 * t105) * t208 - (t185 * t100 + t162 * t102 + t163 * t104) * t211) * t202 + (t135 * t245 - t185 * t137 + t186 * t139) * t245) * t245 + (-t5 - t12 - ((t183 * t101 + t160 * t103 + t161 * t105) * t208 - (t183 * t100 + t160 * t102 + t161 * t104) * t211) * t202 + (-t134 * t243 - t183 * t136 + t184 * t138) * t243 + (-t134 * t245 + t135 * t243 + t185 * t136 + t183 * t137 - t186 * t138 - t184 * t139) * t245) * t243 + (t9 + t15 + (t80 + t53) * t245 + (-t79 - t52) * t243 + ((-t47 - t72) * t211 + (t48 + t73) * t208) * t202 + t259) * t204; m(6) * (t183 * t46 + t185 * t45) + m(5) * (t183 * t71 + t185 * t70) + m(4) * (t183 * t77 + t185 * t76); m(6) * (-t16 * t244 + t183 * t21 + t185 * t22) + m(5) * (t183 * t35 + t185 * t36 - t244 * t26) + m(4) * (t183 * t57 + t185 * t58 - t244 * t49); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t202 ^ 2 * t210 ^ 2 + t183 ^ 2 + t185 ^ 2); (-t33 - t56) * t244 + m(6) * (t29 * t45 + t30 * t46) + m(5) * (t60 * t70 + t61 * t71) + t213 * t185 + t214 * t183; (t8 / 0.2e1 + t14 / 0.2e1) * t204 + (t6 / 0.2e1 + t13 / 0.2e1) * t185 + (t5 / 0.2e1 + t12 / 0.2e1) * t183 + m(6) * (t25 * t16 + t30 * t21 + t29 * t22) + m(5) * (t55 * t26 + t61 * t35 + t60 * t36) + ((-t3 / 0.2e1 - t10 / 0.2e1) * t211 + (-t9 / 0.2e1 - t15 / 0.2e1) * t210 + (t4 / 0.2e1 + t11 / 0.2e1) * t208) * t202; m(5) * (t61 * t183 + t60 * t185 - t244 * t55) + m(6) * (t30 * t183 + t29 * t185 - t244 * t25); (-t14 - t8) * t244 + (t4 + t11) * t185 + (t3 + t10) * t183 + m(6) * (t25 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(5) * (t55 ^ 2 + t60 ^ 2 + t61 ^ 2); m(6) * (t43 * t45 + t44 * t46) + t31 + t229 * t156 + t230 * t154; t9 * t255 + m(6) * (t34 * t16 + t44 * t21 + t43 * t22) + t6 * t256 + t204 * t7 / 0.2e1 + t5 * t257 + (t208 * t258 - t211 * t1 / 0.2e1) * t202; m(6) * (t44 * t183 + t43 * t185 - t244 * t34); m(6) * (t34 * t25 + t43 * t29 + t44 * t30) + t8 * t255 + t3 * t257 + t4 * t256 + t185 * t258 + t183 * t1 / 0.2e1 - t7 * t244 / 0.2e1; m(6) * (t34 ^ 2 + t43 ^ 2 + t44 ^ 2) + t156 * t2 + t154 * t1 + t171 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t59(1), t59(2), t59(4), t59(7), t59(11); t59(2), t59(3), t59(5), t59(8), t59(12); t59(4), t59(5), t59(6), t59(9), t59(13); t59(7), t59(8), t59(9), t59(10), t59(14); t59(11), t59(12), t59(13), t59(14), t59(15);];
Mq = res;
