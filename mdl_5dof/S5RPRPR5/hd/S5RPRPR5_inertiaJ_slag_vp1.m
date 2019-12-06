% Calculate joint inertia matrix for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:55:53
% EndTime: 2019-12-05 17:56:01
% DurationCPUTime: 2.01s
% Computational Cost: add. (3926->332), mult. (4672->468), div. (0->0), fcn. (4933->10), ass. (0->159)
t133 = qJ(3) + pkin(9);
t126 = sin(t133);
t127 = cos(t133);
t136 = sin(pkin(8));
t137 = cos(pkin(8));
t86 = -Icges(5,3) * t137 + (Icges(5,5) * t127 - Icges(5,6) * t126) * t136;
t139 = sin(qJ(3));
t141 = cos(qJ(3));
t98 = -Icges(4,3) * t137 + (Icges(4,5) * t141 - Icges(4,6) * t139) * t136;
t193 = -t86 - t98;
t87 = -Icges(5,6) * t137 + (Icges(5,4) * t127 - Icges(5,2) * t126) * t136;
t99 = -Icges(4,6) * t137 + (Icges(4,4) * t141 - Icges(4,2) * t139) * t136;
t192 = -t126 * t87 - t139 * t99;
t100 = -Icges(4,5) * t137 + (Icges(4,1) * t141 - Icges(4,4) * t139) * t136;
t88 = -Icges(5,5) * t137 + (Icges(5,1) * t127 - Icges(5,4) * t126) * t136;
t191 = (t100 * t141 + t127 * t88) * t136;
t125 = t141 * pkin(3) + pkin(2);
t114 = pkin(4) * t127 + t125;
t158 = t114 - t125;
t190 = t158 * t137;
t183 = pkin(2) - t125;
t189 = pkin(6) * t136 + t137 * t183;
t188 = t136 ^ 2;
t187 = m(5) / 0.2e1;
t186 = m(6) / 0.2e1;
t185 = pkin(3) * t139;
t138 = -qJ(4) - pkin(6);
t182 = t136 * t192 + t137 * t193 + t191;
t142 = cos(qJ(1));
t170 = t136 * t142;
t128 = qJ(5) + t133;
t123 = sin(t128);
t124 = cos(t128);
t140 = sin(qJ(1));
t166 = t140 * t124;
t169 = t137 * t142;
t94 = -t123 * t169 + t166;
t167 = t140 * t123;
t95 = t124 * t169 + t167;
t149 = -t95 * rSges(6,1) - t94 * rSges(6,2);
t51 = rSges(6,3) * t170 - t149;
t83 = -rSges(6,3) * t137 + (rSges(6,1) * t124 - rSges(6,2) * t123) * t136;
t33 = t137 * t51 + t170 * t83;
t165 = t140 * t126;
t102 = t127 * t142 + t137 * t165;
t164 = t140 * t127;
t103 = t126 * t142 - t137 * t164;
t171 = t136 * t140;
t52 = Icges(5,5) * t103 + Icges(5,6) * t102 - Icges(5,3) * t171;
t161 = t141 * t142;
t163 = t140 * t139;
t109 = t137 * t163 + t161;
t162 = t140 * t141;
t168 = t139 * t142;
t110 = -t137 * t162 + t168;
t61 = Icges(4,5) * t110 + Icges(4,6) * t109 - Icges(4,3) * t171;
t181 = -t52 - t61;
t176 = t103 * rSges(5,1) + t102 * rSges(5,2);
t157 = pkin(3) * t168 + t138 * t171;
t67 = t140 * t189 + t157;
t180 = rSges(5,3) * t171 - t176 - t67;
t120 = t138 * t170;
t68 = pkin(3) * t163 - t142 * t189 - t120;
t85 = (pkin(6) + t138) * t137 - t183 * t136;
t179 = t137 * t68 + t170 * t85;
t104 = -t126 * t169 + t164;
t105 = t127 * t169 + t165;
t53 = Icges(5,5) * t105 + Icges(5,6) * t104 + Icges(5,3) * t170;
t111 = -t137 * t168 + t162;
t112 = t137 * t161 + t163;
t62 = Icges(4,5) * t112 + Icges(4,6) * t111 + Icges(4,3) * t170;
t178 = t62 + t53;
t92 = t124 * t142 + t137 * t167;
t93 = t123 * t142 - t137 * t166;
t177 = rSges(6,1) * t93 + rSges(6,2) * t92;
t81 = -Icges(6,6) * t137 + (Icges(6,4) * t124 - Icges(6,2) * t123) * t136;
t175 = t123 * t81;
t82 = -Icges(6,5) * t137 + (Icges(6,1) * t124 - Icges(6,4) * t123) * t136;
t72 = t136 * t124 * t82;
t80 = -Icges(6,3) * t137 + (Icges(6,5) * t124 - Icges(6,6) * t123) * t136;
t27 = -t136 * t175 - t137 * t80 + t72;
t172 = t27 * t137;
t160 = t110 * rSges(4,1) + t109 * rSges(4,2);
t115 = pkin(4) * t126 + t185;
t131 = -pkin(7) + t138;
t159 = t142 * t115 + t131 * t171;
t156 = t140 ^ 2 + t142 ^ 2;
t155 = t187 + t186;
t50 = -rSges(6,3) * t171 + t177;
t154 = t140 * t190 + t157 - t159 - t50 - t67;
t44 = Icges(6,5) * t93 + Icges(6,6) * t92 - Icges(6,3) * t171;
t46 = Icges(6,4) * t93 + Icges(6,2) * t92 - Icges(6,6) * t171;
t48 = Icges(6,1) * t93 + Icges(6,4) * t92 - Icges(6,5) * t171;
t10 = -t137 * t44 + (-t123 * t46 + t124 * t48) * t136;
t45 = Icges(6,5) * t95 + Icges(6,6) * t94 + Icges(6,3) * t170;
t47 = Icges(6,4) * t95 + Icges(6,2) * t94 + Icges(6,6) * t170;
t49 = Icges(6,1) * t95 + Icges(6,4) * t94 + Icges(6,5) * t170;
t11 = -t137 * t45 + (-t123 * t47 + t124 * t49) * t136;
t20 = -t171 * t80 + t81 * t92 + t82 * t93;
t21 = t170 * t80 + t81 * t94 + t82 * t95;
t151 = -(t10 + t20) * t171 / 0.2e1 + (t11 + t21) * t170 / 0.2e1;
t150 = -t114 * t137 - pkin(1);
t148 = -t112 * rSges(4,1) - t111 * rSges(4,2);
t147 = -t105 * rSges(5,1) - t104 * rSges(5,2);
t146 = -rSges(3,1) * t137 + rSges(3,2) * t136 - pkin(1);
t145 = -rSges(5,3) * t136 - t125 * t137 - pkin(1);
t1 = (-t21 * t137 - (t170 * t44 + t46 * t94 + t48 * t95) * t171 + (t170 * t45 + t47 * t94 + t49 * t95) * t170) * t170;
t2 = -t20 * t137 - (-t171 * t44 + t46 * t92 + t48 * t93) * t171 + (-t171 * t45 + t47 * t92 + t49 * t93) * t170;
t3 = -t172 + (-t10 * t140 + t11 * t142) * t136;
t144 = -t137 * t3 - t171 * t2 + t1;
t143 = -pkin(2) * t137 - pkin(1) + (-rSges(4,3) - pkin(6)) * t136;
t129 = t142 * qJ(2);
t118 = -rSges(2,1) * t142 + rSges(2,2) * t140;
t117 = -rSges(2,1) * t140 - rSges(2,2) * t142;
t101 = -rSges(4,3) * t137 + (rSges(4,1) * t141 - rSges(4,2) * t139) * t136;
t89 = -rSges(5,3) * t137 + (rSges(5,1) * t127 - rSges(5,2) * t126) * t136;
t79 = (-rSges(3,3) - qJ(2)) * t140 + t146 * t142;
t78 = rSges(3,3) * t142 + t140 * t146 + t129;
t76 = t85 * t171;
t74 = t83 * t171;
t71 = (t131 - t138) * t137 + t158 * t136;
t70 = rSges(4,3) * t170 - t148;
t69 = -rSges(4,3) * t171 + t160;
t66 = Icges(4,1) * t112 + Icges(4,4) * t111 + Icges(4,5) * t170;
t65 = Icges(4,1) * t110 + Icges(4,4) * t109 - Icges(4,5) * t171;
t64 = Icges(4,4) * t112 + Icges(4,2) * t111 + Icges(4,6) * t170;
t63 = Icges(4,4) * t110 + Icges(4,2) * t109 - Icges(4,6) * t171;
t59 = rSges(5,3) * t170 - t147;
t57 = Icges(5,1) * t105 + Icges(5,4) * t104 + Icges(5,5) * t170;
t56 = Icges(5,1) * t103 + Icges(5,4) * t102 - Icges(5,5) * t171;
t55 = Icges(5,4) * t105 + Icges(5,2) * t104 + Icges(5,6) * t170;
t54 = Icges(5,4) * t103 + Icges(5,2) * t102 - Icges(5,6) * t171;
t42 = -t140 * qJ(2) + t142 * t143 + t148;
t41 = t140 * t143 + t129 + t160;
t40 = t120 + (t115 - t185) * t140 + (-t131 * t136 + t190) * t142;
t38 = t101 * t170 + t137 * t70;
t37 = t101 * t171 - t137 * t69;
t35 = t120 + (-qJ(2) - t185) * t140 + t145 * t142 + t147;
t34 = t140 * t145 + t129 + t157 + t176;
t32 = -t137 * t50 + t74;
t30 = (-t140 * t70 - t142 * t69) * t136;
t29 = t100 * t112 + t111 * t99 + t170 * t98;
t28 = t100 * t110 + t109 * t99 - t171 * t98;
t26 = (-qJ(2) - t115) * t140 + ((-rSges(6,3) + t131) * t136 + t150) * t142 + t149;
t25 = t129 + (-rSges(6,3) * t136 + t150) * t140 + t159 + t177;
t24 = t104 * t87 + t105 * t88 + t170 * t86;
t23 = t102 * t87 + t103 * t88 - t171 * t86;
t22 = (-t140 * t51 - t142 * t50) * t136;
t17 = -t137 * t62 + (-t139 * t64 + t141 * t66) * t136;
t16 = -t137 * t61 + (-t139 * t63 + t141 * t65) * t136;
t15 = t137 * t59 + t170 * t89 + t179;
t14 = t137 * t180 + t171 * t89 + t76;
t13 = -t137 * t53 + (-t126 * t55 + t127 * t57) * t136;
t12 = -t137 * t52 + (-t126 * t54 + t127 * t56) * t136;
t7 = (t180 * t142 + (-t59 - t68) * t140) * t136;
t6 = t137 * t40 + t170 * t71 + t179 + t33;
t5 = t137 * t154 + t171 * t71 + t74 + t76;
t4 = (t154 * t142 + (-t40 - t51 - t68) * t140) * t136;
t8 = [Icges(2,3) + t72 + (Icges(3,2) * t137 + t193 - t80) * t137 + (Icges(3,1) * t136 + 0.2e1 * Icges(3,4) * t137 - t175 + t192) * t136 + m(2) * (t117 ^ 2 + t118 ^ 2) + m(3) * (t78 ^ 2 + t79 ^ 2) + m(4) * (t41 ^ 2 + t42 ^ 2) + m(5) * (t34 ^ 2 + t35 ^ 2) + m(6) * (t25 ^ 2 + t26 ^ 2) + t191; m(3) * (t140 * t78 + t142 * t79) + m(4) * (t140 * t41 + t142 * t42) + m(5) * (t140 * t34 + t142 * t35) + m(6) * (t140 * t25 + t142 * t26); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t155) * t156; (-t27 - t182) * t137 + m(4) * (t37 * t41 + t38 * t42) + m(5) * (t14 * t34 + t15 * t35) + m(6) * (t25 * t5 + t26 * t6) + ((t17 / 0.2e1 + t13 / 0.2e1 + t29 / 0.2e1 + t24 / 0.2e1) * t142 + (-t16 / 0.2e1 - t12 / 0.2e1 - t28 / 0.2e1 - t23 / 0.2e1) * t140) * t136 + t151; m(4) * (t140 * t37 + t142 * t38) + m(5) * (t14 * t140 + t142 * t15) + m(6) * (t140 * t5 + t142 * t6); t1 + m(6) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t14 ^ 2 + t15 ^ 2 + t7 ^ 2) + m(4) * (t30 ^ 2 + t37 ^ 2 + t38 ^ 2) + (t182 * t137 - t3) * t137 + (((t104 * t55 + t105 * t57 + t111 * t64 + t112 * t66 + t170 * t178) * t170 + (-t13 - t17 - t24 - t29) * t137) * t142 + (-t2 + ((t102 * t54 + t103 * t56 + t109 * t63 + t110 * t65) * t136 + t181 * t188 * t140) * t140 + (t23 + t28 + t16 + t12) * t137 + (-t102 * t55 - t103 * t57 - t104 * t54 - t105 * t56 - t109 * t64 - t110 * t66 - t111 * t63 - t112 * t65 + (t140 * t178 + t142 * t181) * t136) * t170) * t140) * t136; 0.2e1 * ((-t140 * t35 + t142 * t34) * t187 + (-t140 * t26 + t142 * t25) * t186) * t136; 0; m(6) * (-t137 * t4 + (-t140 * t6 + t142 * t5) * t136) + m(5) * (-t137 * t7 + (t14 * t142 - t140 * t15) * t136); 0.2e1 * t155 * (t137 ^ 2 + t156 * t188); -t172 + m(6) * (t25 * t32 + t26 * t33) + t151; m(6) * (t140 * t32 + t142 * t33); m(6) * (t22 * t4 + t32 * t5 + t33 * t6) + t144; m(6) * (-t22 * t137 + (-t140 * t33 + t142 * t32) * t136); m(6) * (t22 ^ 2 + t32 ^ 2 + t33 ^ 2) + t144;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;
