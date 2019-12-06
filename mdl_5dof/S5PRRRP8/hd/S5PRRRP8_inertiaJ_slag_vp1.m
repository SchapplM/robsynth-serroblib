% Calculate joint inertia matrix for
% S5PRRRP8
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP8_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP8_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:16
% EndTime: 2019-12-05 16:58:26
% DurationCPUTime: 3.20s
% Computational Cost: add. (11504->419), mult. (29728->613), div. (0->0), fcn. (38165->10), ass. (0->191)
t208 = rSges(6,1) + pkin(4);
t207 = rSges(6,3) + qJ(5);
t206 = cos(qJ(3));
t205 = cos(qJ(4));
t166 = sin(pkin(9));
t168 = cos(pkin(9));
t173 = cos(qJ(2));
t169 = cos(pkin(5));
t172 = sin(qJ(2));
t197 = t172 * t169;
t157 = t173 * t166 + t168 * t197;
t171 = sin(qJ(3));
t167 = sin(pkin(5));
t198 = t168 * t167;
t145 = t157 * t206 - t171 * t198;
t196 = t173 * t169;
t156 = t172 * t166 - t168 * t196;
t170 = sin(qJ(4));
t121 = t145 * t170 - t156 * t205;
t122 = t145 * t205 + t156 * t170;
t178 = t167 * t206;
t144 = t157 * t171 + t168 * t178;
t204 = rSges(6,2) * t144 + t207 * t121 + t208 * t122;
t159 = -t166 * t197 + t173 * t168;
t200 = t167 * t171;
t147 = t159 * t206 + t166 * t200;
t158 = t166 * t196 + t168 * t172;
t123 = t147 * t170 - t158 * t205;
t124 = t147 * t205 + t158 * t170;
t146 = t159 * t171 - t166 * t178;
t203 = rSges(6,2) * t146 + t207 * t123 + t208 * t124;
t119 = pkin(3) * t147 + pkin(8) * t146;
t91 = rSges(5,1) * t124 - rSges(5,2) * t123 + rSges(5,3) * t146;
t202 = -t119 - t91;
t201 = t166 * t167;
t199 = t167 * t173;
t161 = t169 * t171 + t172 * t178;
t148 = t161 * t170 + t205 * t199;
t149 = t161 * t205 - t170 * t199;
t160 = -t169 * t206 + t172 * t200;
t195 = rSges(6,2) * t160 + t207 * t148 + t208 * t149;
t111 = rSges(5,1) * t149 - rSges(5,2) * t148 + rSges(5,3) * t160;
t143 = pkin(3) * t161 + pkin(8) * t160;
t194 = -t111 - t143;
t118 = pkin(3) * t145 + pkin(8) * t144;
t193 = t118 * t199 + t156 * t143;
t142 = pkin(2) * t159 + pkin(7) * t158;
t140 = t169 * t142;
t192 = t169 * t119 + t140;
t141 = pkin(2) * t157 + pkin(7) * t156;
t191 = -t118 - t141;
t190 = t141 * t201 + t142 * t198;
t76 = Icges(6,5) * t122 + Icges(6,6) * t144 + Icges(6,3) * t121;
t80 = Icges(6,4) * t122 + Icges(6,2) * t144 + Icges(6,6) * t121;
t84 = Icges(6,1) * t122 + Icges(6,4) * t144 + Icges(6,5) * t121;
t30 = t121 * t76 + t122 * t84 + t144 * t80;
t77 = Icges(6,5) * t124 + Icges(6,6) * t146 + Icges(6,3) * t123;
t81 = Icges(6,4) * t124 + Icges(6,2) * t146 + Icges(6,6) * t123;
t85 = Icges(6,1) * t124 + Icges(6,4) * t146 + Icges(6,5) * t123;
t31 = t121 * t77 + t122 * t85 + t144 * t81;
t103 = Icges(6,5) * t149 + Icges(6,6) * t160 + Icges(6,3) * t148;
t105 = Icges(6,4) * t149 + Icges(6,2) * t160 + Icges(6,6) * t148;
t107 = Icges(6,1) * t149 + Icges(6,4) * t160 + Icges(6,5) * t148;
t50 = t103 * t121 + t105 * t144 + t107 * t122;
t1 = t144 * t30 + t146 * t31 + t160 * t50;
t78 = Icges(5,5) * t122 - Icges(5,6) * t121 + Icges(5,3) * t144;
t82 = Icges(5,4) * t122 - Icges(5,2) * t121 + Icges(5,6) * t144;
t86 = Icges(5,1) * t122 - Icges(5,4) * t121 + Icges(5,5) * t144;
t32 = -t121 * t82 + t122 * t86 + t144 * t78;
t79 = Icges(5,5) * t124 - Icges(5,6) * t123 + Icges(5,3) * t146;
t83 = Icges(5,4) * t124 - Icges(5,2) * t123 + Icges(5,6) * t146;
t87 = Icges(5,1) * t124 - Icges(5,4) * t123 + Icges(5,5) * t146;
t33 = -t121 * t83 + t122 * t87 + t144 * t79;
t104 = Icges(5,5) * t149 - Icges(5,6) * t148 + Icges(5,3) * t160;
t106 = Icges(5,4) * t149 - Icges(5,2) * t148 + Icges(5,6) * t160;
t108 = Icges(5,1) * t149 - Icges(5,4) * t148 + Icges(5,5) * t160;
t51 = t104 * t144 - t106 * t121 + t108 * t122;
t2 = t144 * t32 + t146 * t33 + t160 * t51;
t189 = -t2 / 0.2e1 - t1 / 0.2e1;
t34 = t123 * t76 + t124 * t84 + t146 * t80;
t35 = t123 * t77 + t124 * t85 + t146 * t81;
t52 = t103 * t123 + t105 * t146 + t107 * t124;
t3 = t144 * t34 + t146 * t35 + t160 * t52;
t36 = -t123 * t82 + t124 * t86 + t146 * t78;
t37 = -t123 * t83 + t124 * t87 + t146 * t79;
t53 = t104 * t146 - t106 * t123 + t108 * t124;
t4 = t144 * t36 + t146 * t37 + t160 * t53;
t188 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = t30 * t156 + t31 * t158 - t50 * t199;
t6 = t32 * t156 + t33 * t158 - t51 * t199;
t187 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t34 * t156 + t35 * t158 - t52 * t199;
t8 = t36 * t156 + t37 * t158 - t53 * t199;
t186 = t7 / 0.2e1 + t8 / 0.2e1;
t10 = t169 * t51 + (t166 * t33 - t168 * t32) * t167;
t9 = t169 * t50 + (t166 * t31 - t168 * t30) * t167;
t185 = t10 / 0.2e1 + t9 / 0.2e1;
t184 = -t119 - t203;
t11 = t169 * t52 + (t166 * t35 - t168 * t34) * t167;
t12 = t169 * t53 + (t166 * t37 - t168 * t36) * t167;
t183 = t12 / 0.2e1 + t11 / 0.2e1;
t41 = t148 * t76 + t149 * t84 + t160 * t80;
t42 = t148 * t77 + t149 * t85 + t160 * t81;
t61 = t103 * t148 + t105 * t160 + t107 * t149;
t13 = t144 * t41 + t146 * t42 + t160 * t61;
t43 = -t148 * t82 + t149 * t86 + t160 * t78;
t44 = -t148 * t83 + t149 * t87 + t160 * t79;
t62 = t104 * t160 - t106 * t148 + t108 * t149;
t14 = t144 * t43 + t146 * t44 + t160 * t62;
t182 = t14 / 0.2e1 + t13 / 0.2e1;
t15 = t41 * t156 + t42 * t158 - t61 * t199;
t16 = t43 * t156 + t44 * t158 - t62 * t199;
t181 = t16 / 0.2e1 + t15 / 0.2e1;
t17 = t169 * t61 + (t166 * t42 - t168 * t41) * t167;
t18 = t169 * t62 + (t166 * t44 - t168 * t43) * t167;
t180 = t18 / 0.2e1 + t17 / 0.2e1;
t179 = -t143 - t195;
t137 = t161 * rSges(4,1) - t160 * rSges(4,2) - rSges(4,3) * t199;
t162 = (pkin(2) * t172 - pkin(7) * t173) * t167;
t177 = (-t137 - t162) * t167;
t176 = t118 * t201 + t119 * t198 + t190;
t175 = (-t162 + t194) * t167;
t174 = (-t162 + t179) * t167;
t153 = t169 * rSges(3,3) + (rSges(3,1) * t172 + rSges(3,2) * t173) * t167;
t152 = Icges(3,5) * t169 + (Icges(3,1) * t172 + Icges(3,4) * t173) * t167;
t151 = Icges(3,6) * t169 + (Icges(3,4) * t172 + Icges(3,2) * t173) * t167;
t150 = Icges(3,3) * t169 + (Icges(3,5) * t172 + Icges(3,6) * t173) * t167;
t136 = Icges(4,1) * t161 - Icges(4,4) * t160 - Icges(4,5) * t199;
t135 = Icges(4,4) * t161 - Icges(4,2) * t160 - Icges(4,6) * t199;
t134 = Icges(4,5) * t161 - Icges(4,6) * t160 - Icges(4,3) * t199;
t133 = rSges(3,1) * t159 - rSges(3,2) * t158 + rSges(3,3) * t201;
t132 = rSges(3,1) * t157 - rSges(3,2) * t156 - rSges(3,3) * t198;
t131 = Icges(3,1) * t159 - Icges(3,4) * t158 + Icges(3,5) * t201;
t130 = Icges(3,1) * t157 - Icges(3,4) * t156 - Icges(3,5) * t198;
t129 = Icges(3,4) * t159 - Icges(3,2) * t158 + Icges(3,6) * t201;
t128 = Icges(3,4) * t157 - Icges(3,2) * t156 - Icges(3,6) * t198;
t127 = Icges(3,5) * t159 - Icges(3,6) * t158 + Icges(3,3) * t201;
t126 = Icges(3,5) * t157 - Icges(3,6) * t156 - Icges(3,3) * t198;
t113 = -t132 * t169 - t153 * t198;
t112 = t133 * t169 - t153 * t201;
t109 = t158 * t118;
t102 = rSges(4,1) * t147 - rSges(4,2) * t146 + rSges(4,3) * t158;
t101 = rSges(4,1) * t145 - rSges(4,2) * t144 + rSges(4,3) * t156;
t100 = Icges(4,1) * t147 - Icges(4,4) * t146 + Icges(4,5) * t158;
t99 = Icges(4,1) * t145 - Icges(4,4) * t144 + Icges(4,5) * t156;
t98 = Icges(4,4) * t147 - Icges(4,2) * t146 + Icges(4,6) * t158;
t97 = Icges(4,4) * t145 - Icges(4,2) * t144 + Icges(4,6) * t156;
t96 = Icges(4,5) * t147 - Icges(4,6) * t146 + Icges(4,3) * t158;
t95 = Icges(4,5) * t145 - Icges(4,6) * t144 + Icges(4,3) * t156;
t94 = (t132 * t166 + t133 * t168) * t167;
t89 = rSges(5,1) * t122 - rSges(5,2) * t121 + rSges(5,3) * t144;
t75 = -t102 * t199 - t158 * t137;
t74 = t101 * t199 + t156 * t137;
t73 = -t134 * t199 - t160 * t135 + t161 * t136;
t72 = t101 * t158 - t102 * t156;
t71 = (-t101 - t141) * t169 + t168 * t177;
t70 = t102 * t169 + t166 * t177 + t140;
t69 = t134 * t158 - t135 * t146 + t136 * t147;
t68 = t134 * t156 - t135 * t144 + t136 * t145;
t67 = (t101 * t166 + t102 * t168) * t167 + t190;
t66 = -t111 * t146 + t160 * t91;
t65 = t111 * t144 - t160 * t89;
t64 = t161 * t100 - t160 * t98 - t96 * t199;
t63 = -t160 * t97 + t161 * t99 - t95 * t199;
t60 = t100 * t147 - t146 * t98 + t158 * t96;
t59 = -t146 * t97 + t147 * t99 + t158 * t95;
t58 = t100 * t145 - t144 * t98 + t156 * t96;
t57 = -t144 * t97 + t145 * t99 + t156 * t95;
t56 = -t144 * t91 + t146 * t89;
t55 = t194 * t158 + t202 * t199;
t54 = t156 * t111 + t89 * t199 + t193;
t49 = (-t89 + t191) * t169 + t168 * t175;
t48 = t166 * t175 + t169 * t91 + t192;
t47 = t202 * t156 + t158 * t89 + t109;
t46 = -t195 * t146 + t203 * t160;
t45 = t195 * t144 - t204 * t160;
t40 = (t166 * t89 + t168 * t91) * t167 + t176;
t39 = t179 * t158 + t184 * t199;
t38 = t195 * t156 + t204 * t199 + t193;
t29 = (t191 - t204) * t169 + t168 * t174;
t28 = t166 * t174 + t203 * t169 + t192;
t27 = -t203 * t144 + t204 * t146;
t26 = t184 * t156 + t204 * t158 + t109;
t25 = (t204 * t166 + t203 * t168) * t167 + t176;
t24 = t169 * t73 + (t166 * t64 - t168 * t63) * t167;
t23 = t63 * t156 + t64 * t158 - t73 * t199;
t22 = t169 * t69 + (t166 * t60 - t168 * t59) * t167;
t21 = t169 * t68 + (t166 * t58 - t168 * t57) * t167;
t20 = t59 * t156 + t60 * t158 - t69 * t199;
t19 = t57 * t156 + t58 * t158 - t68 * t199;
t88 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t94 + m(4) * t67 + m(5) * t40 + m(6) * t25; m(6) * (t25 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(5) * (t40 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(4) * (t67 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(3) * (t112 ^ 2 + t113 ^ 2 + t94 ^ 2) + (t12 + t11 + t22 + (t127 * t201 - t129 * t158 + t131 * t159) * t201) * t201 + (-t10 - t9 - t21 + (-t126 * t198 - t128 * t156 + t130 * t157) * t198 + (-t126 * t201 + t127 * t198 + t128 * t158 + t129 * t156 - t130 * t159 - t131 * t157) * t201) * t198 + ((t150 * t201 - t151 * t158 + t152 * t159) * t201 - (-t150 * t198 - t151 * t156 + t152 * t157) * t198 + t18 + t17 + t24 + ((t129 * t173 + t131 * t172) * t166 - (t128 * t173 + t130 * t172) * t168) * t167 ^ 2 + ((-t126 * t168 + t127 * t166 + t151 * t173 + t152 * t172) * t167 + t169 * t150) * t169) * t169; m(4) * t72 + m(5) * t47 + m(6) * t26; (t23 / 0.2e1 + t181) * t169 + (t22 / 0.2e1 + t183) * t158 + (t21 / 0.2e1 + t185) * t156 + m(6) * (t25 * t26 + t28 * t39 + t29 * t38) + m(5) * (t40 * t47 + t48 * t55 + t49 * t54) + m(4) * (t67 * t72 + t70 * t75 + t71 * t74) + ((-t24 / 0.2e1 - t180) * t173 + (-t19 / 0.2e1 - t187) * t168 + (t20 / 0.2e1 + t186) * t166) * t167; (-t15 - t16 - t23) * t199 + (t7 + t8 + t20) * t158 + (t5 + t6 + t19) * t156 + m(6) * (t26 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(5) * (t47 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(4) * (t72 ^ 2 + t74 ^ 2 + t75 ^ 2); m(5) * t56 + m(6) * t27; t182 * t169 + t180 * t160 + t183 * t146 + t185 * t144 + m(6) * (t25 * t27 + t28 * t46 + t29 * t45) + m(5) * (t40 * t56 + t48 * t66 + t49 * t65) + (t188 * t166 + t189 * t168) * t167; -t182 * t199 + t181 * t160 + t188 * t158 - t189 * t156 + t186 * t146 + t187 * t144 + m(6) * (t26 * t27 + t38 * t45 + t39 * t46) + m(5) * (t47 * t56 + t54 * t65 + t55 * t66); (t14 + t13) * t160 + (t4 + t3) * t146 + (t2 + t1) * t144 + m(6) * (t27 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(5) * (t56 ^ 2 + t65 ^ 2 + t66 ^ 2); m(6) * t148; m(6) * (t121 * t28 + t123 * t29 + t148 * t25); m(6) * (t121 * t39 + t123 * t38 + t148 * t26); m(6) * (t121 * t46 + t123 * t45 + t148 * t27); m(6) * (t121 ^ 2 + t123 ^ 2 + t148 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t88(1), t88(2), t88(4), t88(7), t88(11); t88(2), t88(3), t88(5), t88(8), t88(12); t88(4), t88(5), t88(6), t88(9), t88(13); t88(7), t88(8), t88(9), t88(10), t88(14); t88(11), t88(12), t88(13), t88(14), t88(15);];
Mq = res;
