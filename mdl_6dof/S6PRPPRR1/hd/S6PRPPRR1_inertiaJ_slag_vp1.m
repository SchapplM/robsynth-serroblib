% Calculate joint inertia matrix for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPPRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPPRR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:14:01
% EndTime: 2019-03-08 19:14:11
% DurationCPUTime: 5.35s
% Computational Cost: add. (18387->480), mult. (37322->723), div. (0->0), fcn. (48646->14), ass. (0->206)
t186 = sin(pkin(6));
t224 = sin(pkin(11));
t225 = cos(pkin(11));
t231 = sin(qJ(2));
t232 = cos(qJ(2));
t194 = t232 * t224 + t231 * t225;
t166 = t194 * t186;
t184 = sin(pkin(12));
t187 = cos(pkin(12));
t189 = cos(pkin(6));
t158 = -t166 * t184 + t187 * t189;
t223 = t184 * t189;
t159 = t166 * t187 + t223;
t175 = -t231 * t224 + t232 * t225;
t165 = t175 * t186;
t241 = Icges(4,4) * t166 - Icges(5,5) * t159 + Icges(4,6) * t189 - Icges(5,6) * t158 + (Icges(4,2) + Icges(5,3)) * t165;
t203 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t242 = 0.2e1 * t203;
t185 = sin(pkin(10));
t188 = cos(pkin(10));
t193 = t189 * t175;
t150 = -t185 * t194 + t188 * t193;
t167 = t194 * t189;
t152 = t167 * t188 + t175 * t185;
t221 = t186 * t188;
t101 = Icges(4,5) * t152 + Icges(4,6) * t150 - Icges(4,3) * t221;
t208 = t189 * t232;
t170 = -t185 * t231 + t188 * t208;
t207 = t189 * t231;
t171 = t185 * t232 + t188 * t207;
t139 = Icges(3,5) * t171 + Icges(3,6) * t170 - Icges(3,3) * t221;
t240 = -t139 - t101;
t153 = -t185 * t193 - t188 * t194;
t155 = -t167 * t185 + t175 * t188;
t222 = t185 * t186;
t102 = Icges(4,5) * t155 + Icges(4,6) * t153 + Icges(4,3) * t222;
t172 = -t185 * t208 - t188 * t231;
t173 = -t185 * t207 + t188 * t232;
t140 = Icges(3,5) * t173 + Icges(3,6) * t172 + Icges(3,3) * t222;
t239 = t140 + t102;
t238 = Icges(4,5) * t166 + Icges(4,6) * t165 + (t231 * Icges(3,5) + t232 * Icges(3,6)) * t186 + (Icges(4,3) + Icges(3,3)) * t189;
t215 = pkin(12) + qJ(5);
t182 = sin(t215);
t205 = cos(t215);
t198 = t186 * t205;
t124 = t152 * t182 + t188 * t198;
t126 = t155 * t182 - t185 * t198;
t156 = t166 * t182 - t189 * t205;
t125 = t152 * t205 - t182 * t221;
t191 = sin(qJ(6));
t192 = cos(qJ(6));
t91 = -t125 * t191 - t150 * t192;
t92 = t125 * t192 - t150 * t191;
t53 = Icges(7,5) * t92 + Icges(7,6) * t91 + Icges(7,3) * t124;
t55 = Icges(7,4) * t92 + Icges(7,2) * t91 + Icges(7,6) * t124;
t57 = Icges(7,1) * t92 + Icges(7,4) * t91 + Icges(7,5) * t124;
t19 = t124 * t53 + t55 * t91 + t57 * t92;
t127 = t155 * t205 + t182 * t222;
t93 = -t127 * t191 - t153 * t192;
t94 = t127 * t192 - t153 * t191;
t54 = Icges(7,5) * t94 + Icges(7,6) * t93 + Icges(7,3) * t126;
t56 = Icges(7,4) * t94 + Icges(7,2) * t93 + Icges(7,6) * t126;
t58 = Icges(7,1) * t94 + Icges(7,4) * t93 + Icges(7,5) * t126;
t20 = t124 * t54 + t56 * t91 + t58 * t92;
t157 = t166 * t205 + t189 * t182;
t122 = -t157 * t191 - t165 * t192;
t123 = t157 * t192 - t165 * t191;
t76 = Icges(7,5) * t123 + Icges(7,6) * t122 + Icges(7,3) * t156;
t77 = Icges(7,4) * t123 + Icges(7,2) * t122 + Icges(7,6) * t156;
t78 = Icges(7,1) * t123 + Icges(7,4) * t122 + Icges(7,5) * t156;
t29 = t124 * t76 + t77 * t91 + t78 * t92;
t1 = t124 * t19 + t126 * t20 + t156 * t29;
t237 = -t1 / 0.2e1;
t24 = t122 * t55 + t123 * t57 + t156 * t53;
t25 = t122 * t56 + t123 * t58 + t156 * t54;
t38 = t122 * t77 + t123 * t78 + t156 * t76;
t7 = t124 * t24 + t126 * t25 + t156 * t38;
t236 = t7 / 0.2e1;
t235 = t124 / 0.2e1;
t234 = t126 / 0.2e1;
t233 = t156 / 0.2e1;
t230 = t232 * pkin(2);
t229 = pkin(4) * t187;
t59 = rSges(7,1) * t92 + rSges(7,2) * t91 + rSges(7,3) * t124;
t228 = pkin(5) * t125 + pkin(9) * t124 + t59;
t60 = rSges(7,1) * t94 + rSges(7,2) * t93 + rSges(7,3) * t126;
t227 = pkin(5) * t127 + pkin(9) * t126 + t60;
t79 = rSges(7,1) * t123 + rSges(7,2) * t122 + rSges(7,3) * t156;
t226 = pkin(5) * t157 + pkin(9) * t156 + t79;
t118 = pkin(3) * t155 - qJ(4) * t153;
t206 = pkin(2) * t207 - qJ(3) * t186;
t148 = -t206 * t185 + t230 * t188;
t138 = t189 * t148;
t219 = t189 * t118 + t138;
t117 = pkin(3) * t152 - qJ(4) * t150;
t147 = t230 * t185 + t206 * t188;
t218 = -t117 - t147;
t217 = t147 * t222 + t148 * t221;
t176 = t186 * t231 * pkin(2) + t189 * qJ(3);
t216 = -pkin(3) * t166 + qJ(4) * t165 - t176;
t210 = t184 * t222;
t67 = pkin(4) * t210 - pkin(8) * t153 + t229 * t155;
t214 = t189 * t67 + t219;
t209 = t184 * t221;
t66 = -pkin(4) * t209 - pkin(8) * t150 + t229 * t152;
t213 = -t66 + t218;
t212 = -pkin(4) * t223 + pkin(8) * t165 - t229 * t166 + t216;
t211 = m(4) + m(5) + m(6) + m(7);
t204 = (-rSges(4,1) * t166 - rSges(4,2) * t165 - rSges(4,3) * t189 - t176) * t186;
t202 = t117 * t222 + t118 * t221 + t217;
t201 = (-rSges(5,1) * t159 - rSges(5,2) * t158 + rSges(5,3) * t165 + t216) * t186;
t107 = rSges(6,1) * t157 - rSges(6,2) * t156 - rSges(6,3) * t165;
t197 = (-t107 + t212) * t186;
t196 = t67 * t221 + t66 * t222 + t202;
t195 = (t212 - t226) * t186;
t164 = t189 * rSges(3,3) + (t231 * rSges(3,1) + t232 * rSges(3,2)) * t186;
t163 = Icges(3,5) * t189 + (t231 * Icges(3,1) + t232 * Icges(3,4)) * t186;
t162 = Icges(3,6) * t189 + (t231 * Icges(3,4) + t232 * Icges(3,2)) * t186;
t146 = rSges(3,1) * t173 + rSges(3,2) * t172 + rSges(3,3) * t222;
t145 = rSges(3,1) * t171 + rSges(3,2) * t170 - rSges(3,3) * t221;
t144 = Icges(3,1) * t173 + Icges(3,4) * t172 + Icges(3,5) * t222;
t143 = Icges(3,1) * t171 + Icges(3,4) * t170 - Icges(3,5) * t221;
t142 = Icges(3,4) * t173 + Icges(3,2) * t172 + Icges(3,6) * t222;
t141 = Icges(3,4) * t171 + Icges(3,2) * t170 - Icges(3,6) * t221;
t136 = Icges(4,1) * t166 + Icges(4,4) * t165 + Icges(4,5) * t189;
t131 = t155 * t187 + t210;
t130 = -t155 * t184 + t187 * t222;
t129 = t152 * t187 - t209;
t128 = -t152 * t184 - t187 * t221;
t120 = -t145 * t189 - t164 * t221;
t119 = t146 * t189 - t164 * t222;
t115 = Icges(5,1) * t159 + Icges(5,4) * t158 - Icges(5,5) * t165;
t114 = Icges(5,4) * t159 + Icges(5,2) * t158 - Icges(5,6) * t165;
t109 = rSges(4,1) * t155 + rSges(4,2) * t153 + rSges(4,3) * t222;
t108 = rSges(4,1) * t152 + rSges(4,2) * t150 - rSges(4,3) * t221;
t106 = Icges(4,1) * t155 + Icges(4,4) * t153 + Icges(4,5) * t222;
t105 = Icges(4,1) * t152 + Icges(4,4) * t150 - Icges(4,5) * t221;
t104 = Icges(4,4) * t155 + Icges(4,2) * t153 + Icges(4,6) * t222;
t103 = Icges(4,4) * t152 + Icges(4,2) * t150 - Icges(4,6) * t221;
t100 = Icges(6,1) * t157 - Icges(6,4) * t156 - Icges(6,5) * t165;
t99 = Icges(6,4) * t157 - Icges(6,2) * t156 - Icges(6,6) * t165;
t98 = Icges(6,5) * t157 - Icges(6,6) * t156 - Icges(6,3) * t165;
t96 = (t145 * t185 + t146 * t188) * t186;
t87 = rSges(5,1) * t131 + rSges(5,2) * t130 - rSges(5,3) * t153;
t86 = rSges(5,1) * t129 + rSges(5,2) * t128 - rSges(5,3) * t150;
t85 = Icges(5,1) * t131 + Icges(5,4) * t130 - Icges(5,5) * t153;
t84 = Icges(5,1) * t129 + Icges(5,4) * t128 - Icges(5,5) * t150;
t83 = Icges(5,4) * t131 + Icges(5,2) * t130 - Icges(5,6) * t153;
t82 = Icges(5,4) * t129 + Icges(5,2) * t128 - Icges(5,6) * t150;
t81 = Icges(5,5) * t131 + Icges(5,6) * t130 - Icges(5,3) * t153;
t80 = Icges(5,5) * t129 + Icges(5,6) * t128 - Icges(5,3) * t150;
t75 = rSges(6,1) * t127 - rSges(6,2) * t126 - rSges(6,3) * t153;
t74 = rSges(6,1) * t125 - rSges(6,2) * t124 - rSges(6,3) * t150;
t73 = Icges(6,1) * t127 - Icges(6,4) * t126 - Icges(6,5) * t153;
t72 = Icges(6,1) * t125 - Icges(6,4) * t124 - Icges(6,5) * t150;
t71 = Icges(6,4) * t127 - Icges(6,2) * t126 - Icges(6,6) * t153;
t70 = Icges(6,4) * t125 - Icges(6,2) * t124 - Icges(6,6) * t150;
t69 = Icges(6,5) * t127 - Icges(6,6) * t126 - Icges(6,3) * t153;
t68 = Icges(6,5) * t125 - Icges(6,6) * t124 - Icges(6,3) * t150;
t62 = (-t108 - t147) * t189 + t188 * t204;
t61 = t109 * t189 + t185 * t204 + t138;
t52 = (t108 * t185 + t109 * t188) * t186 + t217;
t51 = t107 * t153 - t165 * t75;
t50 = -t107 * t150 + t165 * t74;
t49 = t100 * t157 - t156 * t99 - t165 * t98;
t48 = t150 * t75 - t153 * t74;
t47 = t100 * t127 - t126 * t99 - t153 * t98;
t46 = t100 * t125 - t124 * t99 - t150 * t98;
t45 = (-t86 + t218) * t189 + t188 * t201;
t44 = t185 * t201 + t189 * t87 + t219;
t43 = -t126 * t79 + t156 * t60;
t42 = t124 * t79 - t156 * t59;
t41 = -t156 * t71 + t157 * t73 - t165 * t69;
t40 = -t156 * t70 + t157 * t72 - t165 * t68;
t39 = (t185 * t86 + t188 * t87) * t186 + t202;
t37 = -t126 * t71 + t127 * t73 - t153 * t69;
t36 = -t126 * t70 + t127 * t72 - t153 * t68;
t35 = -t124 * t71 + t125 * t73 - t150 * t69;
t34 = -t124 * t70 + t125 * t72 - t150 * t68;
t33 = -t124 * t60 + t126 * t59;
t32 = t226 * t153 - t227 * t165;
t31 = -t226 * t150 + t228 * t165;
t30 = t126 * t76 + t77 * t93 + t78 * t94;
t28 = (-t74 + t213) * t189 + t188 * t197;
t27 = t185 * t197 + t189 * t75 + t214;
t26 = t227 * t150 - t228 * t153;
t23 = (t185 * t74 + t188 * t75) * t186 + t196;
t22 = t126 * t54 + t56 * t93 + t58 * t94;
t21 = t126 * t53 + t55 * t93 + t57 * t94;
t18 = (t213 - t228) * t189 + t188 * t195;
t17 = t185 * t195 + t227 * t189 + t214;
t16 = (t228 * t185 + t227 * t188) * t186 + t196;
t15 = t189 * t49 + (t185 * t41 - t188 * t40) * t186;
t14 = -t150 * t40 - t153 * t41 - t165 * t49;
t13 = t189 * t47 + (t185 * t37 - t188 * t36) * t186;
t12 = t189 * t46 + (t185 * t35 - t188 * t34) * t186;
t11 = -t150 * t36 - t153 * t37 - t165 * t47;
t10 = -t150 * t34 - t153 * t35 - t165 * t46;
t9 = t189 * t38 + (t185 * t25 - t188 * t24) * t186;
t8 = -t150 * t24 - t153 * t25 - t165 * t38;
t6 = t189 * t30 + (t185 * t22 - t188 * t21) * t186;
t5 = t189 * t29 + (t185 * t20 - t188 * t19) * t186;
t4 = -t150 * t21 - t153 * t22 - t165 * t30;
t3 = -t150 * t19 - t153 * t20 - t165 * t29;
t2 = t124 * t21 + t126 * t22 + t156 * t30;
t63 = [m(2) + m(3) + t211; m(3) * t96 + m(4) * t52 + m(5) * t39 + m(6) * t23 + m(7) * t16; m(7) * (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) + m(6) * (t23 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(5) * (t39 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(4) * (t52 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(3) * (t119 ^ 2 + t120 ^ 2 + t96 ^ 2) + (t6 + t13 + ((t130 * t83 + t131 * t85 - t153 * t81) * t185 - (t130 * t82 + t131 * t84 - t153 * t80) * t188) * t186 + (t104 * t153 + t106 * t155 + t142 * t172 + t144 * t173 + t222 * t239) * t222) * t222 + (-t5 - t12 - ((t128 * t83 + t129 * t85 - t150 * t81) * t185 - (t128 * t82 + t129 * t84 - t150 * t80) * t188) * t186 + (t103 * t150 + t105 * t152 + t141 * t170 + t143 * t171 + t221 * t240) * t221 + (-t103 * t153 - t104 * t150 - t105 * t155 - t106 * t152 - t141 * t172 - t142 * t170 - t143 * t173 - t144 * t171 + t221 * t239 + t222 * t240) * t222) * t221 + (t9 + t15 + (t158 * t114 + t159 * t115 + t166 * t136 + t241 * t165 + t238 * t189) * t189 + (t114 * t130 + t115 * t131 + t136 * t155 + t189 * t140 + t241 * t153 + t162 * t172 + t163 * t173 + t238 * t222) * t222 + (-t114 * t128 - t115 * t129 - t136 * t152 - t189 * t139 - t241 * t150 - t162 * t170 - t163 * t171 + t238 * t221) * t221 + ((t232 * t142 + t231 * t144) * t222 - (t232 * t141 + t231 * t143) * t221 + (t232 * t162 + t231 * t163) * t189 + (-t101 * t189 - t105 * t166 - t158 * t82 - t159 * t84 - (t103 - t80) * t165) * t188 + (t102 * t189 + t106 * t166 + t158 * t83 + t159 * t85 - (-t104 + t81) * t165) * t185) * t186) * t189; t211 * t189; m(7) * (t16 * t189 + (-t17 * t188 + t18 * t185) * t186) + m(6) * (t189 * t23 + (t185 * t28 - t188 * t27) * t186) + m(5) * (t189 * t39 + (t185 * t45 - t188 * t44) * t186) + m(4) * (t189 * t52 + (t185 * t62 - t188 * t61) * t186); 0.2e1 * (m(4) / 0.2e1 + t203) * (t189 ^ 2 + (t185 ^ 2 + t188 ^ 2) * t186 ^ 2); -t165 * t242; m(7) * (-t150 * t17 - t153 * t18 - t16 * t165) + m(6) * (-t150 * t27 - t153 * t28 - t165 * t23) + m(5) * (-t150 * t44 - t153 * t45 - t165 * t39); (-t189 * t165 + (t150 * t188 - t153 * t185) * t186) * t242; (t150 ^ 2 + t153 ^ 2 + t165 ^ 2) * t242; m(6) * t48 + m(7) * t26; (t8 / 0.2e1 + t14 / 0.2e1) * t189 - (t9 / 0.2e1 + t15 / 0.2e1) * t165 + (-t6 / 0.2e1 - t13 / 0.2e1) * t153 + (-t5 / 0.2e1 - t12 / 0.2e1) * t150 + m(7) * (t16 * t26 + t17 * t32 + t18 * t31) + m(6) * (t23 * t48 + t27 * t51 + t28 * t50) + ((-t3 / 0.2e1 - t10 / 0.2e1) * t188 + (t4 / 0.2e1 + t11 / 0.2e1) * t185) * t186; m(6) * (t189 * t48 + (t185 * t50 - t188 * t51) * t186) + m(7) * (t189 * t26 + (t185 * t31 - t188 * t32) * t186); m(6) * (-t150 * t51 - t153 * t50 - t165 * t48) + m(7) * (-t150 * t32 - t153 * t31 - t165 * t26); -(t8 + t14) * t165 + (-t4 - t11) * t153 + (-t3 - t10) * t150 + m(7) * (t26 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t48 ^ 2 + t50 ^ 2 + t51 ^ 2); m(7) * t33; t189 * t236 + t9 * t233 + t5 * t235 + t6 * t234 + m(7) * (t16 * t33 + t17 * t43 + t18 * t42) + (t185 * t2 / 0.2e1 + t188 * t237) * t186; m(7) * (t189 * t33 + (t185 * t42 - t188 * t43) * t186); m(7) * (-t150 * t43 - t153 * t42 - t165 * t33); t150 * t237 - t165 * t236 + t4 * t234 + m(7) * (t26 * t33 + t31 * t42 + t32 * t43) + t8 * t233 - t153 * t2 / 0.2e1 + t3 * t235; t126 * t2 + t124 * t1 + t156 * t7 + m(7) * (t33 ^ 2 + t42 ^ 2 + t43 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t63(1) t63(2) t63(4) t63(7) t63(11) t63(16); t63(2) t63(3) t63(5) t63(8) t63(12) t63(17); t63(4) t63(5) t63(6) t63(9) t63(13) t63(18); t63(7) t63(8) t63(9) t63(10) t63(14) t63(19); t63(11) t63(12) t63(13) t63(14) t63(15) t63(20); t63(16) t63(17) t63(18) t63(19) t63(20) t63(21);];
Mq  = res;
