% Calculate joint inertia matrix for
% S5PRRPR6
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:08
% EndTime: 2019-12-05 16:30:20
% DurationCPUTime: 3.04s
% Computational Cost: add. (10904->414), mult. (23721->617), div. (0->0), fcn. (30221->12), ass. (0->195)
t212 = m(5) + m(6);
t170 = sin(pkin(9));
t173 = cos(pkin(9));
t178 = cos(qJ(2));
t174 = cos(pkin(5));
t177 = sin(qJ(2));
t195 = t174 * t177;
t156 = t170 * t178 + t173 * t195;
t176 = sin(qJ(3));
t171 = sin(pkin(5));
t207 = cos(qJ(3));
t183 = t171 * t207;
t147 = t156 * t176 + t173 * t183;
t158 = -t170 * t195 + t173 * t178;
t149 = t158 * t176 - t170 * t183;
t197 = t171 * t176;
t159 = -t174 * t207 + t177 * t197;
t150 = t158 * t207 + t170 * t197;
t194 = t174 * t178;
t157 = t170 * t194 + t173 * t177;
t168 = pkin(10) + qJ(5);
t166 = sin(t168);
t167 = cos(t168);
t118 = -t150 * t166 + t157 * t167;
t119 = t150 * t167 + t157 * t166;
t148 = t156 * t207 - t173 * t197;
t155 = t170 * t177 - t173 * t194;
t116 = -t148 * t166 + t155 * t167;
t117 = t148 * t167 + t155 * t166;
t70 = Icges(6,5) * t117 + Icges(6,6) * t116 + Icges(6,3) * t147;
t72 = Icges(6,4) * t117 + Icges(6,2) * t116 + Icges(6,6) * t147;
t74 = Icges(6,1) * t117 + Icges(6,4) * t116 + Icges(6,5) * t147;
t30 = t118 * t72 + t119 * t74 + t149 * t70;
t71 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t149;
t73 = Icges(6,4) * t119 + Icges(6,2) * t118 + Icges(6,6) * t149;
t75 = Icges(6,1) * t119 + Icges(6,4) * t118 + Icges(6,5) * t149;
t31 = t118 * t73 + t119 * t75 + t149 * t71;
t160 = t174 * t176 + t177 * t183;
t196 = t171 * t178;
t143 = -t160 * t166 - t167 * t196;
t144 = t160 * t167 - t166 * t196;
t89 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t159;
t90 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t159;
t91 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t159;
t43 = t118 * t90 + t119 * t91 + t149 * t89;
t2 = t30 * t147 + t31 * t149 + t43 * t159;
t211 = t2 / 0.2e1;
t210 = t147 / 0.2e1;
t209 = t149 / 0.2e1;
t208 = t159 / 0.2e1;
t172 = cos(pkin(10));
t206 = t172 * pkin(4);
t169 = sin(pkin(10));
t201 = t155 * t169;
t76 = t117 * rSges(6,1) + t116 * rSges(6,2) + t147 * rSges(6,3);
t205 = pkin(4) * t201 + pkin(8) * t147 + t206 * t148 + t76;
t200 = t157 * t169;
t77 = t119 * rSges(6,1) + t118 * rSges(6,2) + t149 * rSges(6,3);
t204 = pkin(4) * t200 + pkin(8) * t149 + t206 * t150 + t77;
t184 = t169 * t196;
t92 = t144 * rSges(6,1) + t143 * rSges(6,2) + t159 * rSges(6,3);
t203 = -pkin(4) * t184 + pkin(8) * t159 + t206 * t160 + t92;
t115 = t150 * pkin(3) + t149 * qJ(4);
t122 = -t150 * t169 + t157 * t172;
t123 = t150 * t172 + t200;
t87 = t123 * rSges(5,1) + t122 * rSges(5,2) + t149 * rSges(5,3);
t202 = -t115 - t87;
t199 = t170 * t171;
t198 = t171 * t173;
t145 = -t160 * t169 - t172 * t196;
t146 = t160 * t172 - t184;
t105 = t146 * rSges(5,1) + t145 * rSges(5,2) + t159 * rSges(5,3);
t142 = t160 * pkin(3) + t159 * qJ(4);
t192 = -t105 - t142;
t114 = t148 * pkin(3) + t147 * qJ(4);
t191 = t114 * t196 + t155 * t142;
t141 = t158 * pkin(2) + t157 * pkin(7);
t139 = t174 * t141;
t190 = t174 * t115 + t139;
t140 = t156 * pkin(2) + t155 * pkin(7);
t189 = -t114 - t140;
t188 = t140 * t199 + t141 * t198;
t186 = -t115 - t204;
t185 = -t142 - t203;
t136 = t160 * rSges(4,1) - t159 * rSges(4,2) - rSges(4,3) * t196;
t161 = (pkin(2) * t177 - pkin(7) * t178) * t171;
t182 = (-t136 - t161) * t171;
t181 = t114 * t199 + t115 * t198 + t188;
t180 = (-t161 + t192) * t171;
t179 = (-t161 + t185) * t171;
t154 = t174 * rSges(3,3) + (rSges(3,1) * t177 + rSges(3,2) * t178) * t171;
t153 = Icges(3,5) * t174 + (Icges(3,1) * t177 + Icges(3,4) * t178) * t171;
t152 = Icges(3,6) * t174 + (Icges(3,4) * t177 + Icges(3,2) * t178) * t171;
t151 = Icges(3,3) * t174 + (Icges(3,5) * t177 + Icges(3,6) * t178) * t171;
t135 = Icges(4,1) * t160 - Icges(4,4) * t159 - Icges(4,5) * t196;
t134 = Icges(4,4) * t160 - Icges(4,2) * t159 - Icges(4,6) * t196;
t133 = Icges(4,5) * t160 - Icges(4,6) * t159 - Icges(4,3) * t196;
t132 = t158 * rSges(3,1) - t157 * rSges(3,2) + rSges(3,3) * t199;
t131 = t156 * rSges(3,1) - t155 * rSges(3,2) - rSges(3,3) * t198;
t130 = Icges(3,1) * t158 - Icges(3,4) * t157 + Icges(3,5) * t199;
t129 = Icges(3,1) * t156 - Icges(3,4) * t155 - Icges(3,5) * t198;
t128 = Icges(3,4) * t158 - Icges(3,2) * t157 + Icges(3,6) * t199;
t127 = Icges(3,4) * t156 - Icges(3,2) * t155 - Icges(3,6) * t198;
t126 = Icges(3,5) * t158 - Icges(3,6) * t157 + Icges(3,3) * t199;
t125 = Icges(3,5) * t156 - Icges(3,6) * t155 - Icges(3,3) * t198;
t121 = t148 * t172 + t201;
t120 = -t148 * t169 + t155 * t172;
t109 = -t174 * t131 - t154 * t198;
t108 = t174 * t132 - t154 * t199;
t106 = t157 * t114;
t104 = t150 * rSges(4,1) - t149 * rSges(4,2) + t157 * rSges(4,3);
t103 = t148 * rSges(4,1) - t147 * rSges(4,2) + t155 * rSges(4,3);
t102 = Icges(5,1) * t146 + Icges(5,4) * t145 + Icges(5,5) * t159;
t101 = Icges(5,4) * t146 + Icges(5,2) * t145 + Icges(5,6) * t159;
t100 = Icges(5,5) * t146 + Icges(5,6) * t145 + Icges(5,3) * t159;
t99 = Icges(4,1) * t150 - Icges(4,4) * t149 + Icges(4,5) * t157;
t98 = Icges(4,1) * t148 - Icges(4,4) * t147 + Icges(4,5) * t155;
t97 = Icges(4,4) * t150 - Icges(4,2) * t149 + Icges(4,6) * t157;
t96 = Icges(4,4) * t148 - Icges(4,2) * t147 + Icges(4,6) * t155;
t95 = Icges(4,5) * t150 - Icges(4,6) * t149 + Icges(4,3) * t157;
t94 = Icges(4,5) * t148 - Icges(4,6) * t147 + Icges(4,3) * t155;
t88 = (t131 * t170 + t132 * t173) * t171;
t86 = t121 * rSges(5,1) + t120 * rSges(5,2) + t147 * rSges(5,3);
t85 = Icges(5,1) * t123 + Icges(5,4) * t122 + Icges(5,5) * t149;
t84 = Icges(5,1) * t121 + Icges(5,4) * t120 + Icges(5,5) * t147;
t83 = Icges(5,4) * t123 + Icges(5,2) * t122 + Icges(5,6) * t149;
t82 = Icges(5,4) * t121 + Icges(5,2) * t120 + Icges(5,6) * t147;
t81 = Icges(5,5) * t123 + Icges(5,6) * t122 + Icges(5,3) * t149;
t80 = Icges(5,5) * t121 + Icges(5,6) * t120 + Icges(5,3) * t147;
t79 = -t104 * t196 - t157 * t136;
t78 = t103 * t196 + t155 * t136;
t67 = -t133 * t196 - t159 * t134 + t160 * t135;
t66 = t157 * t103 - t155 * t104;
t65 = (-t103 - t140) * t174 + t173 * t182;
t64 = t174 * t104 + t170 * t182 + t139;
t63 = t157 * t133 - t149 * t134 + t150 * t135;
t62 = t155 * t133 - t147 * t134 + t148 * t135;
t61 = (t103 * t170 + t104 * t173) * t171 + t188;
t60 = -t159 * t97 + t160 * t99 - t95 * t196;
t59 = -t159 * t96 + t160 * t98 - t94 * t196;
t58 = -t149 * t92 + t159 * t77;
t57 = t147 * t92 - t159 * t76;
t56 = t159 * t100 + t145 * t101 + t146 * t102;
t55 = -t149 * t97 + t150 * t99 + t157 * t95;
t54 = -t149 * t96 + t150 * t98 + t157 * t94;
t53 = -t147 * t97 + t148 * t99 + t155 * t95;
t52 = -t147 * t96 + t148 * t98 + t155 * t94;
t51 = t143 * t90 + t144 * t91 + t159 * t89;
t50 = -t147 * t77 + t149 * t76;
t49 = t192 * t157 + t202 * t196;
t48 = t155 * t105 + t86 * t196 + t191;
t47 = t149 * t100 + t122 * t101 + t123 * t102;
t46 = t147 * t100 + t120 * t101 + t121 * t102;
t45 = (-t86 + t189) * t174 + t173 * t180;
t44 = t170 * t180 + t174 * t87 + t190;
t42 = t116 * t90 + t117 * t91 + t147 * t89;
t41 = t202 * t155 + t157 * t86 + t106;
t40 = t145 * t83 + t146 * t85 + t159 * t81;
t39 = t145 * t82 + t146 * t84 + t159 * t80;
t38 = (t170 * t86 + t173 * t87) * t171 + t181;
t37 = t143 * t73 + t144 * t75 + t159 * t71;
t36 = t143 * t72 + t144 * t74 + t159 * t70;
t35 = t122 * t83 + t123 * t85 + t149 * t81;
t34 = t122 * t82 + t123 * t84 + t149 * t80;
t33 = t120 * t83 + t121 * t85 + t147 * t81;
t32 = t120 * t82 + t121 * t84 + t147 * t80;
t29 = t116 * t73 + t117 * t75 + t147 * t71;
t28 = t116 * t72 + t117 * t74 + t147 * t70;
t27 = t185 * t157 + t186 * t196;
t26 = t203 * t155 + t205 * t196 + t191;
t25 = (t189 - t205) * t174 + t173 * t179;
t24 = t170 * t179 + t204 * t174 + t190;
t23 = t67 * t174 + (t170 * t60 - t173 * t59) * t171;
t22 = t59 * t155 + t60 * t157 - t67 * t196;
t21 = t186 * t155 + t205 * t157 + t106;
t20 = (t205 * t170 + t204 * t173) * t171 + t181;
t19 = t63 * t174 + (t170 * t55 - t173 * t54) * t171;
t18 = t62 * t174 + (t170 * t53 - t173 * t52) * t171;
t17 = t54 * t155 + t55 * t157 - t63 * t196;
t16 = t52 * t155 + t53 * t157 - t62 * t196;
t15 = t56 * t174 + (t170 * t40 - t173 * t39) * t171;
t14 = t39 * t155 + t40 * t157 - t56 * t196;
t13 = t51 * t174 + (t170 * t37 - t173 * t36) * t171;
t12 = t36 * t155 + t37 * t157 - t51 * t196;
t11 = t36 * t147 + t37 * t149 + t51 * t159;
t10 = t47 * t174 + (t170 * t35 - t173 * t34) * t171;
t9 = t46 * t174 + (t170 * t33 - t173 * t32) * t171;
t8 = t34 * t155 + t35 * t157 - t47 * t196;
t7 = t32 * t155 + t33 * t157 - t46 * t196;
t6 = t43 * t174 + (t170 * t31 - t173 * t30) * t171;
t5 = t42 * t174 + (t170 * t29 - t173 * t28) * t171;
t4 = t30 * t155 + t31 * t157 - t43 * t196;
t3 = t28 * t155 + t29 * t157 - t42 * t196;
t1 = t28 * t147 + t29 * t149 + t42 * t159;
t68 = [m(2) + m(3) + m(4) + t212; m(3) * t88 + m(4) * t61 + m(5) * t38 + m(6) * t20; m(6) * (t20 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t38 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(4) * (t61 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(3) * (t108 ^ 2 + t109 ^ 2 + t88 ^ 2) + (t6 + t10 + t19 + (t126 * t199 - t157 * t128 + t158 * t130) * t199) * t199 + (-t5 - t18 - t9 + (-t125 * t198 - t155 * t127 + t156 * t129) * t198 + (-t125 * t199 + t126 * t198 + t157 * t127 + t155 * t128 - t158 * t129 - t156 * t130) * t199) * t198 + (-(-t151 * t198 - t155 * t152 + t156 * t153) * t198 + (t151 * t199 - t157 * t152 + t158 * t153) * t199 + t13 + t23 + t15 + ((t128 * t178 + t130 * t177) * t170 - (t127 * t178 + t129 * t177) * t173) * t171 ^ 2 + ((-t125 * t173 + t126 * t170 + t152 * t178 + t153 * t177) * t171 + t174 * t151) * t174) * t174; m(4) * t66 + m(5) * t41 + m(6) * t21; (t12 / 0.2e1 + t22 / 0.2e1 + t14 / 0.2e1) * t174 + (t6 / 0.2e1 + t10 / 0.2e1 + t19 / 0.2e1) * t157 + (t5 / 0.2e1 + t18 / 0.2e1 + t9 / 0.2e1) * t155 + m(6) * (t21 * t20 + t27 * t24 + t26 * t25) + m(5) * (t41 * t38 + t49 * t44 + t48 * t45) + m(4) * (t66 * t61 + t79 * t64 + t78 * t65) + ((-t13 / 0.2e1 - t23 / 0.2e1 - t15 / 0.2e1) * t178 + (-t3 / 0.2e1 - t16 / 0.2e1 - t7 / 0.2e1) * t173 + (t4 / 0.2e1 + t17 / 0.2e1 + t8 / 0.2e1) * t170) * t171; (-t12 - t14 - t22) * t196 + (t4 + t8 + t17) * t157 + (t3 + t7 + t16) * t155 + m(6) * (t21 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(5) * (t41 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(4) * (t66 ^ 2 + t78 ^ 2 + t79 ^ 2); t159 * t212; m(6) * (t147 * t24 + t149 * t25 + t159 * t20) + m(5) * (t147 * t44 + t149 * t45 + t159 * t38); m(6) * (t147 * t27 + t149 * t26 + t159 * t21) + m(5) * (t147 * t49 + t149 * t48 + t159 * t41); (t147 ^ 2 + t149 ^ 2 + t159 ^ 2) * t212; m(6) * t50; m(6) * (t50 * t20 + t58 * t24 + t57 * t25) + t13 * t208 + t174 * t11 / 0.2e1 + t6 * t209 + t5 * t210 + (-t173 * t1 / 0.2e1 + t170 * t211) * t171; t155 * t1 / 0.2e1 + t3 * t210 + t12 * t208 + t4 * t209 + t157 * t211 + m(6) * (t50 * t21 + t57 * t26 + t58 * t27) - t11 * t196 / 0.2e1; m(6) * (t58 * t147 + t57 * t149 + t50 * t159); m(6) * (t50 ^ 2 + t57 ^ 2 + t58 ^ 2) + t149 * t2 + t147 * t1 + t159 * t11;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t68(1), t68(2), t68(4), t68(7), t68(11); t68(2), t68(3), t68(5), t68(8), t68(12); t68(4), t68(5), t68(6), t68(9), t68(13); t68(7), t68(8), t68(9), t68(10), t68(14); t68(11), t68(12), t68(13), t68(14), t68(15);];
Mq = res;
