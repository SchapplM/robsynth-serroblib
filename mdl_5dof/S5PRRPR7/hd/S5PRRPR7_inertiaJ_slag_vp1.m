% Calculate joint inertia matrix for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR7_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR7_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:35:03
% EndTime: 2019-12-05 16:35:16
% DurationCPUTime: 3.22s
% Computational Cost: add. (11924->410), mult. (30921->613), div. (0->0), fcn. (40309->12), ass. (0->190)
t206 = m(5) + m(6);
t168 = sin(pkin(9));
t170 = cos(pkin(9));
t176 = cos(qJ(2));
t171 = cos(pkin(5));
t174 = sin(qJ(2));
t191 = t174 * t171;
t158 = t168 * t176 + t170 * t191;
t173 = sin(qJ(3));
t169 = sin(pkin(5));
t192 = t170 * t169;
t201 = cos(qJ(3));
t148 = t158 * t201 - t173 * t192;
t190 = t176 * t171;
t157 = t168 * t174 - t170 * t190;
t167 = sin(pkin(10));
t196 = cos(pkin(10));
t120 = t148 * t167 - t157 * t196;
t160 = -t168 * t191 + t170 * t176;
t194 = t169 * t173;
t150 = t160 * t201 + t168 * t194;
t159 = t168 * t190 + t170 * t174;
t122 = t150 * t167 - t159 * t196;
t181 = t169 * t201;
t162 = t171 * t173 + t174 * t181;
t193 = t169 * t176;
t145 = t162 * t167 + t193 * t196;
t121 = t148 * t196 + t157 * t167;
t147 = t158 * t173 + t170 * t181;
t172 = sin(qJ(5));
t175 = cos(qJ(5));
t93 = -t121 * t172 + t147 * t175;
t94 = t121 * t175 + t147 * t172;
t62 = Icges(6,5) * t94 + Icges(6,6) * t93 + Icges(6,3) * t120;
t64 = Icges(6,4) * t94 + Icges(6,2) * t93 + Icges(6,6) * t120;
t66 = Icges(6,1) * t94 + Icges(6,4) * t93 + Icges(6,5) * t120;
t123 = t150 * t196 + t159 * t167;
t149 = t160 * t173 - t168 * t181;
t95 = -t123 * t172 + t149 * t175;
t96 = t123 * t175 + t149 * t172;
t18 = t122 * t62 + t64 * t95 + t66 * t96;
t63 = Icges(6,5) * t96 + Icges(6,6) * t95 + Icges(6,3) * t122;
t65 = Icges(6,4) * t96 + Icges(6,2) * t95 + Icges(6,6) * t122;
t67 = Icges(6,1) * t96 + Icges(6,4) * t95 + Icges(6,5) * t122;
t19 = t122 * t63 + t65 * t95 + t67 * t96;
t146 = t162 * t196 - t167 * t193;
t161 = -t171 * t201 + t174 * t194;
t125 = -t146 * t172 + t161 * t175;
t126 = t146 * t175 + t161 * t172;
t86 = Icges(6,5) * t126 + Icges(6,6) * t125 + Icges(6,3) * t145;
t87 = Icges(6,4) * t126 + Icges(6,2) * t125 + Icges(6,6) * t145;
t88 = Icges(6,1) * t126 + Icges(6,4) * t125 + Icges(6,5) * t145;
t35 = t122 * t86 + t87 * t95 + t88 * t96;
t2 = t120 * t18 + t122 * t19 + t145 * t35;
t205 = t2 / 0.2e1;
t204 = t120 / 0.2e1;
t203 = t122 / 0.2e1;
t202 = t145 / 0.2e1;
t68 = rSges(6,1) * t94 + rSges(6,2) * t93 + rSges(6,3) * t120;
t200 = pkin(4) * t121 + pkin(8) * t120 + t68;
t69 = rSges(6,1) * t96 + rSges(6,2) * t95 + rSges(6,3) * t122;
t199 = pkin(4) * t123 + pkin(8) * t122 + t69;
t118 = pkin(3) * t150 + qJ(4) * t149;
t85 = rSges(5,1) * t123 - rSges(5,2) * t122 + rSges(5,3) * t149;
t198 = -t118 - t85;
t89 = rSges(6,1) * t126 + rSges(6,2) * t125 + rSges(6,3) * t145;
t197 = pkin(4) * t146 + pkin(8) * t145 + t89;
t195 = t168 * t169;
t108 = rSges(5,1) * t146 - rSges(5,2) * t145 + rSges(5,3) * t161;
t144 = pkin(3) * t162 + qJ(4) * t161;
t189 = -t108 - t144;
t117 = pkin(3) * t148 + qJ(4) * t147;
t188 = t117 * t193 + t157 * t144;
t143 = pkin(2) * t160 + pkin(7) * t159;
t141 = t171 * t143;
t187 = t171 * t118 + t141;
t142 = pkin(2) * t158 + pkin(7) * t157;
t186 = -t117 - t142;
t185 = t142 * t195 + t143 * t192;
t183 = -t118 - t199;
t182 = -t144 - t197;
t138 = rSges(4,1) * t162 - rSges(4,2) * t161 - rSges(4,3) * t193;
t163 = (pkin(2) * t174 - pkin(7) * t176) * t169;
t180 = (-t138 - t163) * t169;
t179 = t117 * t195 + t118 * t192 + t185;
t178 = (-t163 + t189) * t169;
t177 = (-t163 + t182) * t169;
t154 = t171 * rSges(3,3) + (rSges(3,1) * t174 + rSges(3,2) * t176) * t169;
t153 = Icges(3,5) * t171 + (Icges(3,1) * t174 + Icges(3,4) * t176) * t169;
t152 = Icges(3,6) * t171 + (Icges(3,4) * t174 + Icges(3,2) * t176) * t169;
t151 = Icges(3,3) * t171 + (Icges(3,5) * t174 + Icges(3,6) * t176) * t169;
t137 = Icges(4,1) * t162 - Icges(4,4) * t161 - Icges(4,5) * t193;
t136 = Icges(4,4) * t162 - Icges(4,2) * t161 - Icges(4,6) * t193;
t135 = Icges(4,5) * t162 - Icges(4,6) * t161 - Icges(4,3) * t193;
t134 = rSges(3,1) * t160 - rSges(3,2) * t159 + rSges(3,3) * t195;
t133 = rSges(3,1) * t158 - rSges(3,2) * t157 - rSges(3,3) * t192;
t132 = Icges(3,1) * t160 - Icges(3,4) * t159 + Icges(3,5) * t195;
t131 = Icges(3,1) * t158 - Icges(3,4) * t157 - Icges(3,5) * t192;
t130 = Icges(3,4) * t160 - Icges(3,2) * t159 + Icges(3,6) * t195;
t129 = Icges(3,4) * t158 - Icges(3,2) * t157 - Icges(3,6) * t192;
t128 = Icges(3,5) * t160 - Icges(3,6) * t159 + Icges(3,3) * t195;
t127 = Icges(3,5) * t158 - Icges(3,6) * t157 - Icges(3,3) * t192;
t112 = -t133 * t171 - t154 * t192;
t111 = t134 * t171 - t154 * t195;
t109 = t159 * t117;
t107 = rSges(4,1) * t150 - rSges(4,2) * t149 + rSges(4,3) * t159;
t106 = rSges(4,1) * t148 - rSges(4,2) * t147 + rSges(4,3) * t157;
t105 = Icges(5,1) * t146 - Icges(5,4) * t145 + Icges(5,5) * t161;
t104 = Icges(5,4) * t146 - Icges(5,2) * t145 + Icges(5,6) * t161;
t103 = Icges(5,5) * t146 - Icges(5,6) * t145 + Icges(5,3) * t161;
t102 = Icges(4,1) * t150 - Icges(4,4) * t149 + Icges(4,5) * t159;
t101 = Icges(4,1) * t148 - Icges(4,4) * t147 + Icges(4,5) * t157;
t100 = Icges(4,4) * t150 - Icges(4,2) * t149 + Icges(4,6) * t159;
t99 = Icges(4,4) * t148 - Icges(4,2) * t147 + Icges(4,6) * t157;
t98 = Icges(4,5) * t150 - Icges(4,6) * t149 + Icges(4,3) * t159;
t97 = Icges(4,5) * t148 - Icges(4,6) * t147 + Icges(4,3) * t157;
t92 = (t133 * t168 + t134 * t170) * t169;
t84 = rSges(5,1) * t121 - rSges(5,2) * t120 + rSges(5,3) * t147;
t83 = Icges(5,1) * t123 - Icges(5,4) * t122 + Icges(5,5) * t149;
t82 = Icges(5,1) * t121 - Icges(5,4) * t120 + Icges(5,5) * t147;
t81 = Icges(5,4) * t123 - Icges(5,2) * t122 + Icges(5,6) * t149;
t80 = Icges(5,4) * t121 - Icges(5,2) * t120 + Icges(5,6) * t147;
t79 = Icges(5,5) * t123 - Icges(5,6) * t122 + Icges(5,3) * t149;
t78 = Icges(5,5) * t121 - Icges(5,6) * t120 + Icges(5,3) * t147;
t77 = -t107 * t193 - t138 * t159;
t76 = t106 * t193 + t138 * t157;
t75 = -t135 * t193 - t136 * t161 + t137 * t162;
t74 = t106 * t159 - t107 * t157;
t73 = (-t106 - t142) * t171 + t170 * t180;
t72 = t171 * t107 + t168 * t180 + t141;
t71 = t135 * t159 - t136 * t149 + t137 * t150;
t70 = t135 * t157 - t136 * t147 + t137 * t148;
t61 = (t106 * t168 + t107 * t170) * t169 + t185;
t60 = -t100 * t161 + t102 * t162 - t193 * t98;
t59 = t101 * t162 - t161 * t99 - t193 * t97;
t58 = t103 * t161 - t104 * t145 + t105 * t146;
t57 = -t100 * t149 + t102 * t150 + t159 * t98;
t56 = t101 * t150 - t149 * t99 + t159 * t97;
t55 = -t100 * t147 + t102 * t148 + t157 * t98;
t54 = t101 * t148 - t147 * t99 + t157 * t97;
t53 = t159 * t189 + t193 * t198;
t52 = t108 * t157 + t193 * t84 + t188;
t51 = t103 * t149 - t104 * t122 + t105 * t123;
t50 = t103 * t147 - t104 * t120 + t105 * t121;
t49 = (-t84 + t186) * t171 + t170 * t178;
t48 = t168 * t178 + t171 * t85 + t187;
t47 = -t122 * t89 + t145 * t69;
t46 = t120 * t89 - t145 * t68;
t45 = t157 * t198 + t159 * t84 + t109;
t44 = -t145 * t81 + t146 * t83 + t161 * t79;
t43 = -t145 * t80 + t146 * t82 + t161 * t78;
t42 = (t168 * t84 + t170 * t85) * t169 + t179;
t41 = t125 * t87 + t126 * t88 + t145 * t86;
t40 = -t122 * t81 + t123 * t83 + t149 * t79;
t39 = -t122 * t80 + t123 * t82 + t149 * t78;
t38 = -t120 * t81 + t121 * t83 + t147 * t79;
t37 = -t120 * t80 + t121 * t82 + t147 * t78;
t36 = -t120 * t69 + t122 * t68;
t34 = t120 * t86 + t87 * t93 + t88 * t94;
t33 = t159 * t182 + t183 * t193;
t32 = t157 * t197 + t193 * t200 + t188;
t31 = (t186 - t200) * t171 + t170 * t177;
t30 = t168 * t177 + t171 * t199 + t187;
t29 = t75 * t171 + (t168 * t60 - t170 * t59) * t169;
t28 = t157 * t59 + t159 * t60 - t193 * t75;
t27 = t125 * t65 + t126 * t67 + t145 * t63;
t26 = t125 * t64 + t126 * t66 + t145 * t62;
t25 = t71 * t171 + (t168 * t57 - t170 * t56) * t169;
t24 = t70 * t171 + (t168 * t55 - t170 * t54) * t169;
t23 = t157 * t56 + t159 * t57 - t193 * t71;
t22 = t157 * t54 + t159 * t55 - t193 * t70;
t21 = t157 * t183 + t159 * t200 + t109;
t20 = (t168 * t200 + t170 * t199) * t169 + t179;
t17 = t120 * t63 + t65 * t93 + t67 * t94;
t16 = t120 * t62 + t64 * t93 + t66 * t94;
t15 = t58 * t171 + (t168 * t44 - t170 * t43) * t169;
t14 = t157 * t43 + t159 * t44 - t193 * t58;
t13 = t51 * t171 + (t168 * t40 - t170 * t39) * t169;
t12 = t50 * t171 + (t168 * t38 - t170 * t37) * t169;
t11 = t157 * t39 + t159 * t40 - t193 * t51;
t10 = t157 * t37 + t159 * t38 - t193 * t50;
t9 = t41 * t171 + (t168 * t27 - t170 * t26) * t169;
t8 = t157 * t26 + t159 * t27 - t193 * t41;
t7 = t120 * t26 + t122 * t27 + t145 * t41;
t6 = t35 * t171 + (t168 * t19 - t170 * t18) * t169;
t5 = t34 * t171 + (-t16 * t170 + t168 * t17) * t169;
t4 = t157 * t18 + t159 * t19 - t193 * t35;
t3 = t157 * t16 + t159 * t17 - t193 * t34;
t1 = t120 * t16 + t122 * t17 + t145 * t34;
t90 = [m(2) + m(3) + m(4) + t206; m(3) * t92 + m(4) * t61 + m(5) * t42 + m(6) * t20; m(6) * (t20 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(5) * (t42 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(4) * (t61 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(3) * (t111 ^ 2 + t112 ^ 2 + t92 ^ 2) + (t6 + t13 + t25 + (t128 * t195 - t130 * t159 + t132 * t160) * t195) * t195 + (-t5 - t12 - t24 + (-t127 * t192 - t129 * t157 + t131 * t158) * t192 + (-t127 * t195 + t128 * t192 + t129 * t159 + t130 * t157 - t131 * t160 - t132 * t158) * t195) * t192 + (-(-t151 * t192 - t157 * t152 + t158 * t153) * t192 + (t151 * t195 - t159 * t152 + t160 * t153) * t195 + t9 + t29 + t15 + ((t130 * t176 + t132 * t174) * t168 - (t129 * t176 + t131 * t174) * t170) * t169 ^ 2 + ((-t127 * t170 + t128 * t168 + t152 * t176 + t153 * t174) * t169 + t171 * t151) * t171) * t171; m(4) * t74 + m(5) * t45 + m(6) * t21; (t8 / 0.2e1 + t28 / 0.2e1 + t14 / 0.2e1) * t171 + (t6 / 0.2e1 + t25 / 0.2e1 + t13 / 0.2e1) * t159 + (t5 / 0.2e1 + t12 / 0.2e1 + t24 / 0.2e1) * t157 + m(6) * (t20 * t21 + t30 * t33 + t31 * t32) + m(5) * (t42 * t45 + t48 * t53 + t49 * t52) + m(4) * (t61 * t74 + t72 * t77 + t73 * t76) + ((-t9 / 0.2e1 - t29 / 0.2e1 - t15 / 0.2e1) * t176 + (-t3 / 0.2e1 - t10 / 0.2e1 - t22 / 0.2e1) * t170 + (t4 / 0.2e1 + t11 / 0.2e1 + t23 / 0.2e1) * t168) * t169; (-t14 - t28 - t8) * t193 + (t4 + t11 + t23) * t159 + (t3 + t22 + t10) * t157 + m(6) * (t21 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(5) * (t45 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(4) * (t74 ^ 2 + t76 ^ 2 + t77 ^ 2); t161 * t206; m(6) * (t147 * t30 + t149 * t31 + t161 * t20) + m(5) * (t147 * t48 + t149 * t49 + t161 * t42); m(6) * (t147 * t33 + t149 * t32 + t161 * t21) + m(5) * (t147 * t53 + t149 * t52 + t161 * t45); (t147 ^ 2 + t149 ^ 2 + t161 ^ 2) * t206; m(6) * t36; m(6) * (t20 * t36 + t30 * t47 + t31 * t46) + t6 * t203 + t5 * t204 + t171 * t7 / 0.2e1 + t9 * t202 + (t168 * t205 - t170 * t1 / 0.2e1) * t169; m(6) * (t21 * t36 + t32 * t46 + t33 * t47) + t8 * t202 + t157 * t1 / 0.2e1 + t3 * t204 + t4 * t203 + t159 * t205 - t7 * t193 / 0.2e1; m(6) * (t147 * t47 + t149 * t46 + t161 * t36); m(6) * (t36 ^ 2 + t46 ^ 2 + t47 ^ 2) + t122 * t2 + t120 * t1 + t145 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t90(1), t90(2), t90(4), t90(7), t90(11); t90(2), t90(3), t90(5), t90(8), t90(12); t90(4), t90(5), t90(6), t90(9), t90(13); t90(7), t90(8), t90(9), t90(10), t90(14); t90(11), t90(12), t90(13), t90(14), t90(15);];
Mq = res;
