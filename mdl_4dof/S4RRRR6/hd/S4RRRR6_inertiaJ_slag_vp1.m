% Calculate joint inertia matrix for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR6_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR6_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:23
% EndTime: 2019-12-31 17:29:28
% DurationCPUTime: 1.65s
% Computational Cost: add. (7177->336), mult. (18401->508), div. (0->0), fcn. (23366->10), ass. (0->166)
t156 = cos(pkin(4));
t159 = sin(qJ(2));
t163 = cos(qJ(1));
t179 = t163 * t159;
t160 = sin(qJ(1));
t162 = cos(qJ(2));
t182 = t160 * t162;
t142 = t156 * t179 + t182;
t158 = sin(qJ(3));
t155 = sin(pkin(4));
t192 = cos(qJ(3));
t172 = t155 * t192;
t123 = t142 * t158 + t163 * t172;
t178 = t163 * t162;
t183 = t160 * t159;
t144 = -t156 * t183 + t178;
t125 = t144 * t158 - t160 * t172;
t186 = t155 * t159;
t139 = -t156 * t192 + t158 * t186;
t180 = t163 * t155;
t124 = t142 * t192 - t158 * t180;
t141 = -t156 * t178 + t183;
t157 = sin(qJ(4));
t161 = cos(qJ(4));
t95 = -t124 * t157 + t141 * t161;
t96 = t124 * t161 + t141 * t157;
t56 = Icges(5,5) * t96 + Icges(5,6) * t95 + Icges(5,3) * t123;
t58 = Icges(5,4) * t96 + Icges(5,2) * t95 + Icges(5,6) * t123;
t60 = Icges(5,1) * t96 + Icges(5,4) * t95 + Icges(5,5) * t123;
t184 = t160 * t155;
t126 = t144 * t192 + t158 * t184;
t143 = t156 * t182 + t179;
t97 = -t126 * t157 + t143 * t161;
t98 = t126 * t161 + t143 * t157;
t18 = t125 * t56 + t97 * t58 + t98 * t60;
t57 = Icges(5,5) * t98 + Icges(5,6) * t97 + Icges(5,3) * t125;
t59 = Icges(5,4) * t98 + Icges(5,2) * t97 + Icges(5,6) * t125;
t61 = Icges(5,1) * t98 + Icges(5,4) * t97 + Icges(5,5) * t125;
t19 = t125 * t57 + t97 * t59 + t98 * t61;
t140 = t156 * t158 + t159 * t172;
t185 = t155 * t162;
t121 = -t140 * t157 - t161 * t185;
t122 = t140 * t161 - t157 * t185;
t73 = Icges(5,5) * t122 + Icges(5,6) * t121 + Icges(5,3) * t139;
t74 = Icges(5,4) * t122 + Icges(5,2) * t121 + Icges(5,6) * t139;
t75 = Icges(5,1) * t122 + Icges(5,4) * t121 + Icges(5,5) * t139;
t27 = t125 * t73 + t97 * t74 + t98 * t75;
t2 = t18 * t123 + t19 * t125 + t27 * t139;
t196 = t2 / 0.2e1;
t195 = t123 / 0.2e1;
t194 = t125 / 0.2e1;
t193 = t139 / 0.2e1;
t191 = t124 * pkin(3);
t33 = t121 * t74 + t122 * t75 + t139 * t73;
t100 = Icges(4,4) * t140 - Icges(4,2) * t139 - Icges(4,6) * t185;
t101 = Icges(4,1) * t140 - Icges(4,4) * t139 - Icges(4,5) * t185;
t99 = Icges(4,5) * t140 - Icges(4,6) * t139 - Icges(4,3) * t185;
t51 = -t139 * t100 + t140 * t101 - t99 * t185;
t190 = -t33 - t51;
t168 = -t96 * rSges(5,1) - t95 * rSges(5,2);
t62 = t123 * rSges(5,3) - t168;
t189 = t123 * pkin(8) + t191 + t62;
t63 = t98 * rSges(5,1) + t97 * rSges(5,2) + t125 * rSges(5,3);
t188 = t126 * pkin(3) + t125 * pkin(8) + t63;
t76 = t122 * rSges(5,1) + t121 * rSges(5,2) + t139 * rSges(5,3);
t187 = t140 * pkin(3) + t139 * pkin(8) + t76;
t181 = t160 * t163;
t115 = t142 * pkin(2) + t141 * pkin(7);
t116 = t144 * pkin(2) + t143 * pkin(7);
t177 = t115 * t184 + t116 * t180;
t176 = t163 * pkin(1) + pkin(6) * t184;
t21 = t121 * t58 + t122 * t60 + t139 * t56;
t26 = t123 * t73 + t95 * t74 + t96 * t75;
t175 = t21 / 0.2e1 + t26 / 0.2e1;
t22 = t121 * t59 + t122 * t61 + t139 * t57;
t174 = t27 / 0.2e1 + t22 / 0.2e1;
t84 = t126 * rSges(4,1) - t125 * rSges(4,2) + t143 * rSges(4,3);
t130 = Icges(3,3) * t156 + (Icges(3,5) * t159 + Icges(3,6) * t162) * t155;
t131 = Icges(3,6) * t156 + (Icges(3,4) * t159 + Icges(3,2) * t162) * t155;
t132 = Icges(3,5) * t156 + (Icges(3,1) * t159 + Icges(3,4) * t162) * t155;
t173 = t156 * t130 + t131 * t185 + t132 * t186;
t110 = t144 * rSges(3,1) - t143 * rSges(3,2) + rSges(3,3) * t184;
t171 = -t160 * pkin(1) + pkin(6) * t180;
t102 = t140 * rSges(4,1) - t139 * rSges(4,2) - rSges(4,3) * t185;
t145 = (pkin(2) * t159 - pkin(7) * t162) * t155;
t170 = t155 * (-t102 - t145);
t169 = t155 * (-t145 - t187);
t167 = t116 + t176;
t77 = Icges(4,5) * t124 - Icges(4,6) * t123 + Icges(4,3) * t141;
t79 = Icges(4,4) * t124 - Icges(4,2) * t123 + Icges(4,6) * t141;
t81 = Icges(4,1) * t124 - Icges(4,4) * t123 + Icges(4,5) * t141;
t38 = -t139 * t79 + t140 * t81 - t77 * t185;
t45 = -t123 * t100 + t124 * t101 + t141 * t99;
t166 = t45 / 0.2e1 + t38 / 0.2e1 + t175;
t78 = Icges(4,5) * t126 - Icges(4,6) * t125 + Icges(4,3) * t143;
t80 = Icges(4,4) * t126 - Icges(4,2) * t125 + Icges(4,6) * t143;
t82 = Icges(4,1) * t126 - Icges(4,4) * t125 + Icges(4,5) * t143;
t39 = -t139 * t80 + t140 * t82 - t78 * t185;
t46 = -t125 * t100 + t126 * t101 + t143 * t99;
t165 = t39 / 0.2e1 + t46 / 0.2e1 + t174;
t164 = -t115 + t171;
t83 = t124 * rSges(4,1) - t123 * rSges(4,2) + t141 * rSges(4,3);
t109 = t142 * rSges(3,1) - t141 * rSges(3,2) - rSges(3,3) * t180;
t147 = t163 * rSges(2,1) - t160 * rSges(2,2);
t146 = -t160 * rSges(2,1) - t163 * rSges(2,2);
t133 = t156 * rSges(3,3) + (rSges(3,1) * t159 + rSges(3,2) * t162) * t155;
t113 = t156 * t116;
t108 = Icges(3,1) * t144 - Icges(3,4) * t143 + Icges(3,5) * t184;
t107 = Icges(3,1) * t142 - Icges(3,4) * t141 - Icges(3,5) * t180;
t106 = Icges(3,4) * t144 - Icges(3,2) * t143 + Icges(3,6) * t184;
t105 = Icges(3,4) * t142 - Icges(3,2) * t141 - Icges(3,6) * t180;
t104 = Icges(3,5) * t144 - Icges(3,6) * t143 + Icges(3,3) * t184;
t103 = Icges(3,5) * t142 - Icges(3,6) * t141 - Icges(3,3) * t180;
t92 = t110 + t176;
t91 = -t109 + t171;
t86 = -t156 * t109 - t133 * t180;
t85 = t156 * t110 - t133 * t184;
t72 = t173 * t156;
t70 = (t109 * t160 + t110 * t163) * t155;
t69 = t130 * t184 - t143 * t131 + t144 * t132;
t68 = -t130 * t180 - t141 * t131 + t142 * t132;
t65 = t167 + t84;
t64 = t164 - t83;
t55 = -t143 * t102 - t84 * t185;
t54 = t141 * t102 + t83 * t185;
t53 = t156 * t104 + (t106 * t162 + t108 * t159) * t155;
t52 = t156 * t103 + (t105 * t162 + t107 * t159) * t155;
t50 = t51 * t156;
t49 = -t141 * t84 + t143 * t83;
t48 = (-t115 - t83) * t156 + t163 * t170;
t47 = t156 * t84 + t160 * t170 + t113;
t44 = t167 + t188;
t43 = -t191 + (-rSges(5,3) - pkin(8)) * t123 + t164 + t168;
t42 = (t160 * t83 + t163 * t84) * t155 + t177;
t41 = -t125 * t76 + t139 * t63;
t40 = t123 * t76 - t139 * t62;
t37 = -t125 * t80 + t126 * t82 + t143 * t78;
t36 = -t125 * t79 + t126 * t81 + t143 * t77;
t35 = -t123 * t80 + t124 * t82 + t141 * t78;
t34 = -t123 * t79 + t124 * t81 + t141 * t77;
t32 = t33 * t156;
t31 = t33 * t139;
t30 = -t123 * t63 + t125 * t62;
t29 = -t187 * t143 - t188 * t185;
t28 = t187 * t141 + t189 * t185;
t25 = (-t115 - t189) * t156 + t163 * t169;
t24 = t188 * t156 + t160 * t169 + t113;
t23 = -t188 * t141 + t189 * t143;
t20 = (t189 * t160 + t188 * t163) * t155 + t177;
t17 = t123 * t57 + t95 * t59 + t96 * t61;
t16 = t123 * t56 + t95 * t58 + t96 * t60;
t15 = t50 + (t39 * t160 - t38 * t163) * t155;
t14 = t38 * t141 + t39 * t143 - t51 * t185;
t13 = t46 * t156 + (t160 * t37 - t163 * t36) * t155;
t12 = t45 * t156 + (t160 * t35 - t163 * t34) * t155;
t11 = t36 * t141 + t37 * t143 - t46 * t185;
t10 = t34 * t141 + t35 * t143 - t45 * t185;
t9 = t32 + (t22 * t160 - t21 * t163) * t155;
t8 = t21 * t141 + t22 * t143 - t33 * t185;
t7 = t21 * t123 + t22 * t125 + t31;
t6 = t27 * t156 + (t160 * t19 - t163 * t18) * t155;
t5 = t26 * t156 + (-t16 * t163 + t160 * t17) * t155;
t4 = t18 * t141 + t19 * t143 - t27 * t185;
t3 = t16 * t141 + t17 * t143 - t26 * t185;
t1 = t16 * t123 + t17 * t125 + t26 * t139;
t66 = [Icges(2,3) + m(5) * (t43 ^ 2 + t44 ^ 2) + m(4) * (t64 ^ 2 + t65 ^ 2) + m(3) * (t91 ^ 2 + t92 ^ 2) + m(2) * (t146 ^ 2 + t147 ^ 2) + t173 - t190; t32 + t50 + t72 + m(5) * (t24 * t44 + t25 * t43) + m(4) * (t47 * t65 + t48 * t64) + m(3) * (t85 * t92 + t86 * t91) + ((-t68 / 0.2e1 - t52 / 0.2e1 - t166) * t163 + (t53 / 0.2e1 + t69 / 0.2e1 + t165) * t160) * t155; (t9 + t15 + t72) * t156 + m(5) * (t20 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(4) * (t42 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(3) * (t70 ^ 2 + t85 ^ 2 + t86 ^ 2) + (t160 * t6 - t163 * t5 + t160 * t13 - t163 * t12 + (t160 * ((-t143 * t106 + t144 * t108) * t160 - (-t143 * t105 + t144 * t107) * t163) - t163 * ((-t141 * t106 + t142 * t108) * t160 - (-t141 * t105 + t142 * t107) * t163) + (t160 * (t104 * t160 ^ 2 - t103 * t181) - t163 * (t103 * t163 ^ 2 - t104 * t181)) * t155) * t155 + ((-t52 - t68) * t163 + (t53 + t69) * t160) * t156) * t155; t190 * t185 + m(5) * (t28 * t43 + t29 * t44) + m(4) * (t54 * t64 + t55 * t65) + t165 * t143 + t166 * t141; (t8 / 0.2e1 + t14 / 0.2e1) * t156 + (t6 / 0.2e1 + t13 / 0.2e1) * t143 + (t5 / 0.2e1 + t12 / 0.2e1) * t141 + m(5) * (t23 * t20 + t29 * t24 + t28 * t25) + m(4) * (t49 * t42 + t55 * t47 + t54 * t48) + ((-t3 / 0.2e1 - t10 / 0.2e1) * t163 + (-t9 / 0.2e1 - t15 / 0.2e1) * t162 + (t4 / 0.2e1 + t11 / 0.2e1) * t160) * t155; (-t14 - t8) * t185 + (t4 + t11) * t143 + (t3 + t10) * t141 + m(5) * (t23 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(4) * (t49 ^ 2 + t54 ^ 2 + t55 ^ 2); t31 + m(5) * (t40 * t43 + t41 * t44) + t174 * t125 + t175 * t123; t6 * t194 + t5 * t195 + t156 * t7 / 0.2e1 + t9 * t193 + m(5) * (t30 * t20 + t41 * t24 + t40 * t25) + (t160 * t196 - t163 * t1 / 0.2e1) * t155; -t7 * t185 / 0.2e1 + t141 * t1 / 0.2e1 + t4 * t194 + t3 * t195 + t8 * t193 + t143 * t196 + m(5) * (t30 * t23 + t40 * t28 + t41 * t29); m(5) * (t30 ^ 2 + t40 ^ 2 + t41 ^ 2) + t125 * t2 + t123 * t1 + t139 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t66(1), t66(2), t66(4), t66(7); t66(2), t66(3), t66(5), t66(8); t66(4), t66(5), t66(6), t66(9); t66(7), t66(8), t66(9), t66(10);];
Mq = res;
