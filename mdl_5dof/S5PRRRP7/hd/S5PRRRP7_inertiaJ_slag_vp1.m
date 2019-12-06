% Calculate joint inertia matrix for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP7_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP7_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:09
% EndTime: 2019-12-05 16:54:19
% DurationCPUTime: 3.16s
% Computational Cost: add. (11927->424), mult. (30357->620), div. (0->0), fcn. (38827->10), ass. (0->194)
t212 = rSges(6,3) + qJ(5);
t211 = cos(qJ(3));
t173 = cos(qJ(4));
t210 = pkin(4) * t173;
t165 = sin(pkin(9));
t167 = cos(pkin(9));
t174 = cos(qJ(2));
t168 = cos(pkin(5));
t172 = sin(qJ(2));
t198 = t172 * t168;
t155 = t165 * t174 + t167 * t198;
t171 = sin(qJ(3));
t166 = sin(pkin(5));
t199 = t167 * t166;
t145 = t155 * t211 - t171 * t199;
t197 = t174 * t168;
t154 = t165 * t172 - t167 * t197;
t170 = sin(qJ(4));
t121 = -t145 * t170 + t154 * t173;
t204 = t154 * t170;
t122 = t145 * t173 + t204;
t179 = t166 * t211;
t144 = t155 * t171 + t167 * t179;
t208 = rSges(6,1) * t122 + rSges(6,2) * t121 + pkin(4) * t204 + t144 * t212 + t210 * t145;
t157 = -t165 * t198 + t167 * t174;
t201 = t166 * t171;
t147 = t157 * t211 + t165 * t201;
t156 = t165 * t197 + t167 * t172;
t123 = -t147 * t170 + t156 * t173;
t203 = t156 * t170;
t124 = t147 * t173 + t203;
t146 = t157 * t171 - t165 * t179;
t207 = rSges(6,1) * t124 + rSges(6,2) * t123 + pkin(4) * t203 + t146 * t212 + t210 * t147;
t159 = t168 * t171 + t172 * t179;
t200 = t166 * t174;
t148 = -t159 * t170 - t173 * t200;
t180 = t170 * t200;
t149 = t159 * t173 - t180;
t158 = -t168 * t211 + t172 * t201;
t206 = rSges(6,1) * t149 + rSges(6,2) * t148 - pkin(4) * t180 + t158 * t212 + t210 * t159;
t120 = pkin(3) * t147 + pkin(8) * t146;
t93 = rSges(5,1) * t124 + rSges(5,2) * t123 + rSges(5,3) * t146;
t205 = -t120 - t93;
t202 = t165 * t166;
t112 = rSges(5,1) * t149 + rSges(5,2) * t148 + rSges(5,3) * t158;
t143 = t159 * pkin(3) + t158 * pkin(8);
t196 = -t112 - t143;
t119 = pkin(3) * t145 + pkin(8) * t144;
t195 = t119 * t200 + t154 * t143;
t142 = pkin(2) * t157 + pkin(7) * t156;
t140 = t168 * t142;
t194 = t168 * t120 + t140;
t141 = pkin(2) * t155 + pkin(7) * t154;
t193 = -t119 - t141;
t192 = t141 * t202 + t142 * t199;
t78 = Icges(6,5) * t122 + Icges(6,6) * t121 + Icges(6,3) * t144;
t82 = Icges(6,4) * t122 + Icges(6,2) * t121 + Icges(6,6) * t144;
t86 = Icges(6,1) * t122 + Icges(6,4) * t121 + Icges(6,5) * t144;
t32 = t121 * t82 + t122 * t86 + t144 * t78;
t79 = Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t146;
t83 = Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t146;
t87 = Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t146;
t33 = t121 * t83 + t122 * t87 + t144 * t79;
t104 = Icges(6,5) * t149 + Icges(6,6) * t148 + Icges(6,3) * t158;
t106 = Icges(6,4) * t149 + Icges(6,2) * t148 + Icges(6,6) * t158;
t108 = Icges(6,1) * t149 + Icges(6,4) * t148 + Icges(6,5) * t158;
t50 = t104 * t144 + t106 * t121 + t108 * t122;
t1 = t144 * t32 + t146 * t33 + t158 * t50;
t80 = Icges(5,5) * t122 + Icges(5,6) * t121 + Icges(5,3) * t144;
t84 = Icges(5,4) * t122 + Icges(5,2) * t121 + Icges(5,6) * t144;
t88 = Icges(5,1) * t122 + Icges(5,4) * t121 + Icges(5,5) * t144;
t34 = t121 * t84 + t122 * t88 + t144 * t80;
t81 = Icges(5,5) * t124 + Icges(5,6) * t123 + Icges(5,3) * t146;
t85 = Icges(5,4) * t124 + Icges(5,2) * t123 + Icges(5,6) * t146;
t89 = Icges(5,1) * t124 + Icges(5,4) * t123 + Icges(5,5) * t146;
t35 = t121 * t85 + t122 * t89 + t144 * t81;
t105 = Icges(5,5) * t149 + Icges(5,6) * t148 + Icges(5,3) * t158;
t107 = Icges(5,4) * t149 + Icges(5,2) * t148 + Icges(5,6) * t158;
t109 = Icges(5,1) * t149 + Icges(5,4) * t148 + Icges(5,5) * t158;
t51 = t105 * t144 + t107 * t121 + t109 * t122;
t2 = t144 * t34 + t146 * t35 + t158 * t51;
t191 = t2 / 0.2e1 + t1 / 0.2e1;
t36 = t123 * t82 + t124 * t86 + t146 * t78;
t37 = t123 * t83 + t124 * t87 + t146 * t79;
t52 = t104 * t146 + t106 * t123 + t108 * t124;
t3 = t144 * t36 + t146 * t37 + t158 * t52;
t38 = t123 * t84 + t124 * t88 + t146 * t80;
t39 = t123 * t85 + t124 * t89 + t146 * t81;
t53 = t105 * t146 + t107 * t123 + t109 * t124;
t4 = t144 * t38 + t146 * t39 + t158 * t53;
t190 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = t154 * t32 + t156 * t33 - t200 * t50;
t6 = t154 * t34 + t156 * t35 - t200 * t51;
t189 = t5 / 0.2e1 + t6 / 0.2e1;
t7 = t154 * t36 + t156 * t37 - t200 * t52;
t8 = t154 * t38 + t156 * t39 - t200 * t53;
t188 = t7 / 0.2e1 + t8 / 0.2e1;
t10 = t168 * t51 + (t165 * t35 - t167 * t34) * t166;
t9 = t168 * t50 + (t165 * t33 - t167 * t32) * t166;
t187 = t9 / 0.2e1 + t10 / 0.2e1;
t186 = -t120 - t207;
t11 = t168 * t52 + (t165 * t37 - t167 * t36) * t166;
t12 = t168 * t53 + (t165 * t39 - t167 * t38) * t166;
t185 = t12 / 0.2e1 + t11 / 0.2e1;
t43 = t148 * t82 + t149 * t86 + t158 * t78;
t44 = t148 * t83 + t149 * t87 + t158 * t79;
t61 = t104 * t158 + t106 * t148 + t108 * t149;
t13 = t144 * t43 + t146 * t44 + t158 * t61;
t45 = t148 * t84 + t149 * t88 + t158 * t80;
t46 = t148 * t85 + t149 * t89 + t158 * t81;
t62 = t105 * t158 + t107 * t148 + t109 * t149;
t14 = t144 * t45 + t146 * t46 + t158 * t62;
t184 = t14 / 0.2e1 + t13 / 0.2e1;
t15 = t154 * t43 + t156 * t44 - t200 * t61;
t16 = t154 * t45 + t156 * t46 - t200 * t62;
t183 = t15 / 0.2e1 + t16 / 0.2e1;
t17 = t168 * t61 + (t165 * t44 - t167 * t43) * t166;
t18 = t168 * t62 + (t165 * t46 - t167 * t45) * t166;
t182 = t17 / 0.2e1 + t18 / 0.2e1;
t181 = -t143 - t206;
t137 = rSges(4,1) * t159 - rSges(4,2) * t158 - rSges(4,3) * t200;
t160 = (pkin(2) * t172 - pkin(7) * t174) * t166;
t178 = (-t137 - t160) * t166;
t177 = t119 * t202 + t120 * t199 + t192;
t176 = (-t160 + t196) * t166;
t175 = (-t160 + t181) * t166;
t153 = t168 * rSges(3,3) + (rSges(3,1) * t172 + rSges(3,2) * t174) * t166;
t152 = Icges(3,5) * t168 + (Icges(3,1) * t172 + Icges(3,4) * t174) * t166;
t151 = Icges(3,6) * t168 + (Icges(3,4) * t172 + Icges(3,2) * t174) * t166;
t150 = Icges(3,3) * t168 + (Icges(3,5) * t172 + Icges(3,6) * t174) * t166;
t136 = Icges(4,1) * t159 - Icges(4,4) * t158 - Icges(4,5) * t200;
t135 = Icges(4,4) * t159 - Icges(4,2) * t158 - Icges(4,6) * t200;
t134 = Icges(4,5) * t159 - Icges(4,6) * t158 - Icges(4,3) * t200;
t133 = rSges(3,1) * t157 - rSges(3,2) * t156 + rSges(3,3) * t202;
t132 = rSges(3,1) * t155 - rSges(3,2) * t154 - rSges(3,3) * t199;
t131 = Icges(3,1) * t157 - Icges(3,4) * t156 + Icges(3,5) * t202;
t130 = Icges(3,1) * t155 - Icges(3,4) * t154 - Icges(3,5) * t199;
t129 = Icges(3,4) * t157 - Icges(3,2) * t156 + Icges(3,6) * t202;
t128 = Icges(3,4) * t155 - Icges(3,2) * t154 - Icges(3,6) * t199;
t127 = Icges(3,5) * t157 - Icges(3,6) * t156 + Icges(3,3) * t202;
t126 = Icges(3,5) * t155 - Icges(3,6) * t154 - Icges(3,3) * t199;
t114 = -t132 * t168 - t153 * t199;
t113 = t133 * t168 - t153 * t202;
t110 = t156 * t119;
t103 = rSges(4,1) * t147 - rSges(4,2) * t146 + rSges(4,3) * t156;
t102 = rSges(4,1) * t145 - rSges(4,2) * t144 + rSges(4,3) * t154;
t101 = Icges(4,1) * t147 - Icges(4,4) * t146 + Icges(4,5) * t156;
t100 = Icges(4,1) * t145 - Icges(4,4) * t144 + Icges(4,5) * t154;
t99 = Icges(4,4) * t147 - Icges(4,2) * t146 + Icges(4,6) * t156;
t98 = Icges(4,4) * t145 - Icges(4,2) * t144 + Icges(4,6) * t154;
t97 = Icges(4,5) * t147 - Icges(4,6) * t146 + Icges(4,3) * t156;
t96 = Icges(4,5) * t145 - Icges(4,6) * t144 + Icges(4,3) * t154;
t94 = (t132 * t165 + t133 * t167) * t166;
t91 = rSges(5,1) * t122 + rSges(5,2) * t121 + rSges(5,3) * t144;
t77 = -t103 * t200 - t137 * t156;
t76 = t102 * t200 + t137 * t154;
t73 = -t134 * t200 - t135 * t158 + t136 * t159;
t72 = t102 * t156 - t103 * t154;
t71 = (-t102 - t141) * t168 + t167 * t178;
t70 = t103 * t168 + t165 * t178 + t140;
t69 = t134 * t156 - t135 * t146 + t136 * t147;
t68 = t134 * t154 - t135 * t144 + t136 * t145;
t67 = (t102 * t165 + t103 * t167) * t166 + t192;
t66 = -t112 * t146 + t158 * t93;
t65 = t112 * t144 - t158 * t91;
t64 = t101 * t159 - t158 * t99 - t200 * t97;
t63 = t100 * t159 - t158 * t98 - t200 * t96;
t60 = t101 * t147 - t146 * t99 + t156 * t97;
t59 = t100 * t147 - t146 * t98 + t156 * t96;
t58 = t101 * t145 - t144 * t99 + t154 * t97;
t57 = t100 * t145 - t144 * t98 + t154 * t96;
t56 = -t144 * t93 + t146 * t91;
t55 = t156 * t196 + t200 * t205;
t54 = t112 * t154 + t200 * t91 + t195;
t49 = (-t91 + t193) * t168 + t167 * t176;
t48 = t165 * t176 + t168 * t93 + t194;
t47 = t154 * t205 + t156 * t91 + t110;
t42 = (t165 * t91 + t167 * t93) * t166 + t177;
t41 = -t146 * t206 + t158 * t207;
t40 = t144 * t206 - t158 * t208;
t31 = t156 * t181 + t186 * t200;
t30 = t154 * t206 + t200 * t208 + t195;
t29 = (t193 - t208) * t168 + t167 * t175;
t28 = t165 * t175 + t168 * t207 + t194;
t27 = -t144 * t207 + t146 * t208;
t26 = t168 * t73 + (t165 * t64 - t167 * t63) * t166;
t25 = t154 * t63 + t156 * t64 - t200 * t73;
t24 = t154 * t186 + t156 * t208 + t110;
t23 = (t165 * t208 + t167 * t207) * t166 + t177;
t22 = t168 * t69 + (t165 * t60 - t167 * t59) * t166;
t21 = t168 * t68 + (t165 * t58 - t167 * t57) * t166;
t20 = t154 * t59 + t156 * t60 - t200 * t69;
t19 = t154 * t57 + t156 * t58 - t200 * t68;
t74 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t94 + m(4) * t67 + m(5) * t42 + m(6) * t23; m(6) * (t23 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(5) * (t42 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(4) * (t67 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(3) * (t113 ^ 2 + t114 ^ 2 + t94 ^ 2) + (t12 + t11 + t22 + (t127 * t202 - t129 * t156 + t131 * t157) * t202) * t202 + (-t9 - t10 - t21 + (-t126 * t199 - t128 * t154 + t130 * t155) * t199 + (-t126 * t202 + t127 * t199 + t128 * t156 + t129 * t154 - t130 * t157 - t131 * t155) * t202) * t199 + (-(-t150 * t199 - t151 * t154 + t152 * t155) * t199 + (t150 * t202 - t151 * t156 + t152 * t157) * t202 + t18 + t17 + t26 + ((t129 * t174 + t131 * t172) * t165 - (t128 * t174 + t130 * t172) * t167) * t166 ^ 2 + ((-t126 * t167 + t127 * t165 + t151 * t174 + t152 * t172) * t166 + t168 * t150) * t168) * t168; m(4) * t72 + m(5) * t47 + m(6) * t24; (t25 / 0.2e1 + t183) * t168 + (t22 / 0.2e1 + t185) * t156 + (t21 / 0.2e1 + t187) * t154 + m(6) * (t23 * t24 + t28 * t31 + t29 * t30) + m(5) * (t42 * t47 + t48 * t55 + t49 * t54) + m(4) * (t67 * t72 + t70 * t77 + t71 * t76) + ((-t26 / 0.2e1 - t182) * t174 + (-t19 / 0.2e1 - t189) * t167 + (t20 / 0.2e1 + t188) * t165) * t166; (-t15 - t16 - t25) * t200 + (t7 + t8 + t20) * t156 + (t6 + t5 + t19) * t154 + m(6) * (t24 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(5) * (t47 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(4) * (t72 ^ 2 + t76 ^ 2 + t77 ^ 2); m(5) * t56 + m(6) * t27; t184 * t168 + t182 * t158 + t185 * t146 + t187 * t144 + m(6) * (t23 * t27 + t28 * t41 + t29 * t40) + m(5) * (t42 * t56 + t48 * t66 + t49 * t65) + (t165 * t190 - t167 * t191) * t166; -t184 * t200 + t183 * t158 + t190 * t156 + t191 * t154 + t188 * t146 + t189 * t144 + m(6) * (t24 * t27 + t30 * t40 + t31 * t41) + m(5) * (t47 * t56 + t54 * t65 + t55 * t66); (t14 + t13) * t158 + (t3 + t4) * t146 + (t2 + t1) * t144 + m(6) * (t27 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(5) * (t56 ^ 2 + t65 ^ 2 + t66 ^ 2); m(6) * t158; m(6) * (t144 * t28 + t146 * t29 + t158 * t23); m(6) * (t144 * t31 + t146 * t30 + t158 * t24); m(6) * (t144 * t41 + t146 * t40 + t158 * t27); m(6) * (t144 ^ 2 + t146 ^ 2 + t158 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t74(1), t74(2), t74(4), t74(7), t74(11); t74(2), t74(3), t74(5), t74(8), t74(12); t74(4), t74(5), t74(6), t74(9), t74(13); t74(7), t74(8), t74(9), t74(10), t74(14); t74(11), t74(12), t74(13), t74(14), t74(15);];
Mq = res;
