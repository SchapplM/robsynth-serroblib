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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:41:27
% EndTime: 2020-01-03 11:41:39
% DurationCPUTime: 2.46s
% Computational Cost: add. (3926->333), mult. (4672->468), div. (0->0), fcn. (4933->10), ass. (0->156)
t139 = qJ(3) + pkin(9);
t130 = sin(t139);
t131 = cos(t139);
t142 = sin(pkin(8));
t143 = cos(pkin(8));
t86 = -Icges(5,3) * t143 + (Icges(5,5) * t131 - Icges(5,6) * t130) * t142;
t145 = sin(qJ(3));
t147 = cos(qJ(3));
t98 = -Icges(4,3) * t143 + (Icges(4,5) * t147 - Icges(4,6) * t145) * t142;
t194 = -t86 - t98;
t87 = -Icges(5,6) * t143 + (Icges(5,4) * t131 - Icges(5,2) * t130) * t142;
t99 = -Icges(4,6) * t143 + (Icges(4,4) * t147 - Icges(4,2) * t145) * t142;
t193 = -t130 * t87 - t145 * t99;
t100 = -Icges(4,5) * t143 + (Icges(4,1) * t147 - Icges(4,4) * t145) * t142;
t88 = -Icges(5,5) * t143 + (Icges(5,1) * t131 - Icges(5,4) * t130) * t142;
t192 = (t100 * t147 + t131 * t88) * t142;
t191 = t142 ^ 2;
t190 = m(5) / 0.2e1;
t189 = m(6) / 0.2e1;
t188 = pkin(2) * t143;
t187 = pkin(3) * t145;
t144 = -qJ(4) - pkin(6);
t186 = pkin(6) + t144;
t129 = t147 * pkin(3) + pkin(2);
t185 = t193 * t142 + t194 * t143 + t192;
t148 = cos(qJ(1));
t175 = t142 * t148;
t146 = sin(qJ(1));
t176 = t142 * t146;
t132 = qJ(5) + t139;
t128 = cos(t132);
t127 = sin(t132);
t170 = t146 * t127;
t92 = -t128 * t148 - t143 * t170;
t169 = t146 * t128;
t93 = -t127 * t148 + t143 * t169;
t52 = t93 * rSges(6,1) + t92 * rSges(6,2) + rSges(6,3) * t176;
t173 = t143 * t148;
t94 = t127 * t173 - t169;
t95 = -t128 * t173 - t170;
t152 = -t95 * rSges(6,1) - t94 * rSges(6,2);
t53 = -rSges(6,3) * t175 - t152;
t22 = t52 * t175 + t53 * t176;
t168 = t146 * t130;
t102 = -t131 * t148 - t143 * t168;
t167 = t146 * t131;
t103 = -t130 * t148 + t143 * t167;
t54 = Icges(5,5) * t103 + Icges(5,6) * t102 + Icges(5,3) * t176;
t164 = t147 * t148;
t166 = t146 * t145;
t110 = -t143 * t166 - t164;
t165 = t146 * t147;
t172 = t145 * t148;
t111 = t143 * t165 - t172;
t65 = Icges(4,5) * t111 + Icges(4,6) * t110 + Icges(4,3) * t176;
t184 = t54 + t65;
t174 = t143 * t146;
t115 = t129 * t174;
t71 = -pkin(3) * t172 + t115 + (-t186 * t142 - t188) * t146;
t162 = pkin(2) * t173 + pkin(6) * t175;
t163 = pkin(3) * t166 + t129 * t173;
t72 = t144 * t175 + t162 - t163;
t183 = t71 * t175 + t72 * t176;
t104 = t130 * t173 - t167;
t105 = -t131 * t173 - t168;
t55 = Icges(5,5) * t105 + Icges(5,6) * t104 - Icges(5,3) * t175;
t112 = t143 * t172 - t165;
t113 = -t143 * t164 - t166;
t66 = Icges(4,5) * t113 + Icges(4,6) * t112 - Icges(4,3) * t175;
t182 = -t66 - t55;
t81 = -Icges(6,6) * t143 + (Icges(6,4) * t128 - Icges(6,2) * t127) * t142;
t181 = t127 * t81;
t82 = -Icges(6,5) * t143 + (Icges(6,1) * t128 - Icges(6,4) * t127) * t142;
t76 = t142 * t128 * t82;
t80 = -Icges(6,3) * t143 + (Icges(6,5) * t128 - Icges(6,6) * t127) * t142;
t27 = -t142 * t181 - t143 * t80 + t76;
t178 = t27 * t143;
t114 = pkin(4) * t131 + t129;
t177 = t114 * t143;
t117 = pkin(4) * t130 + t187;
t171 = t146 * t117;
t161 = t148 * pkin(1) + t146 * qJ(2);
t137 = -pkin(7) + t144;
t160 = t137 - t144;
t159 = t146 ^ 2 + t148 ^ 2;
t158 = t190 + t189;
t60 = t103 * rSges(5,1) + t102 * rSges(5,2) + rSges(5,3) * t176;
t73 = t111 * rSges(4,1) + t110 * rSges(4,2) + rSges(4,3) * t176;
t46 = Icges(6,5) * t93 + Icges(6,6) * t92 + Icges(6,3) * t176;
t48 = Icges(6,4) * t93 + Icges(6,2) * t92 + Icges(6,6) * t176;
t50 = Icges(6,1) * t93 + Icges(6,4) * t92 + Icges(6,5) * t176;
t10 = -t143 * t46 + (-t127 * t48 + t128 * t50) * t142;
t47 = Icges(6,5) * t95 + Icges(6,6) * t94 - Icges(6,3) * t175;
t49 = Icges(6,4) * t95 + Icges(6,2) * t94 - Icges(6,6) * t175;
t51 = Icges(6,1) * t95 + Icges(6,4) * t94 - Icges(6,5) * t175;
t11 = -t143 * t47 + (-t127 * t49 + t128 * t51) * t142;
t20 = t80 * t176 + t81 * t92 + t82 * t93;
t21 = -t80 * t175 + t94 * t81 + t95 * t82;
t155 = (t10 + t20) * t176 / 0.2e1 - (t11 + t21) * t175 / 0.2e1;
t85 = t186 * t143 + (-pkin(2) + t129) * t142;
t154 = t142 * (-t85 + rSges(5,3) * t143 - (rSges(5,1) * t131 - rSges(5,2) * t130) * t142);
t83 = -rSges(6,3) * t143 + (rSges(6,1) * t128 - rSges(6,2) * t127) * t142;
t153 = t142 * (-t160 * t143 - (t114 - t129) * t142 - t83 - t85);
t151 = rSges(3,1) * t143 - rSges(3,2) * t142;
t150 = -t105 * rSges(5,1) - t104 * rSges(5,2);
t1 = (-t20 * t143 + (t46 * t176 + t48 * t92 + t50 * t93) * t176 - (t47 * t176 + t49 * t92 + t51 * t93) * t175) * t176;
t2 = -t21 * t143 + (-t46 * t175 + t94 * t48 + t95 * t50) * t176 - (-t47 * t175 + t94 * t49 + t95 * t51) * t175;
t3 = -t178 + (t10 * t146 - t11 * t148) * t142;
t149 = -t143 * t3 - t2 * t175 + t1;
t74 = t113 * rSges(4,1) + t112 * rSges(4,2) - rSges(4,3) * t175;
t134 = t146 * pkin(1);
t119 = rSges(2,1) * t148 - t146 * rSges(2,2);
t118 = t146 * rSges(2,1) + rSges(2,2) * t148;
t106 = t114 * t174;
t101 = -rSges(4,3) * t143 + (rSges(4,1) * t147 - rSges(4,2) * t145) * t142;
t79 = t146 * rSges(3,3) + t151 * t148 + t161;
t78 = t134 + (-rSges(3,3) - qJ(2)) * t148 + t151 * t146;
t70 = Icges(4,1) * t113 + Icges(4,4) * t112 - Icges(4,5) * t175;
t69 = Icges(4,1) * t111 + Icges(4,4) * t110 + Icges(4,5) * t176;
t68 = Icges(4,4) * t113 + Icges(4,2) * t112 - Icges(4,6) * t175;
t67 = Icges(4,4) * t111 + Icges(4,2) * t110 + Icges(4,6) * t176;
t64 = t143 * t72;
t61 = -rSges(5,3) * t175 - t150;
t59 = Icges(5,1) * t105 + Icges(5,4) * t104 - Icges(5,5) * t175;
t58 = Icges(5,1) * t103 + Icges(5,4) * t102 + Icges(5,5) * t176;
t57 = Icges(5,4) * t105 + Icges(5,2) * t104 - Icges(5,6) * t175;
t56 = Icges(5,4) * t103 + Icges(5,2) * t102 + Icges(5,6) * t176;
t45 = t143 * t53;
t42 = -t74 + t161 + t162;
t41 = -qJ(2) * t148 + t134 + (pkin(6) * t142 + t188) * t146 + t73;
t40 = -t171 + (t160 * t142 - t177) * t148 + t163;
t39 = t106 - t115 + (-t117 + t187) * t148 - t160 * t176;
t38 = -t101 * t175 + t143 * t74;
t37 = -t101 * t176 - t143 * t73;
t35 = (rSges(5,3) - t144) * t175 + t150 + t161 + t163;
t34 = -t144 * t176 + t115 + t134 + (-qJ(2) - t187) * t148 + t60;
t33 = -t83 * t175 + t45;
t32 = -t143 * t52 - t83 * t176;
t30 = (t146 * t74 + t148 * t73) * t142;
t29 = t113 * t100 + t112 * t99 - t98 * t175;
t28 = t100 * t111 + t110 * t99 + t98 * t176;
t26 = t171 + (t177 + (rSges(6,3) - t137) * t142) * t148 + t152 + t161;
t25 = -t137 * t176 + t106 + t134 + (-qJ(2) - t117) * t148 + t52;
t24 = t104 * t87 + t105 * t88 - t86 * t175;
t23 = t102 * t87 + t103 * t88 + t86 * t176;
t17 = -t143 * t66 + (-t145 * t68 + t147 * t70) * t142;
t16 = -t143 * t65 + (-t145 * t67 + t147 * t69) * t142;
t15 = t143 * t61 + t148 * t154 + t64;
t14 = (-t60 - t71) * t143 + t146 * t154;
t13 = -t143 * t55 + (-t130 * t57 + t131 * t59) * t142;
t12 = -t143 * t54 + (-t130 * t56 + t131 * t58) * t142;
t7 = (t146 * t61 + t148 * t60) * t142 + t183;
t6 = t143 * t40 + t148 * t153 + t45 + t64;
t5 = (-t39 - t52 - t71) * t143 + t146 * t153;
t4 = (t146 * t40 + t148 * t39) * t142 + t183 + t22;
t8 = [Icges(2,3) + t76 + (Icges(3,2) * t143 + t194 - t80) * t143 + (Icges(3,1) * t142 + 0.2e1 * Icges(3,4) * t143 - t181 + t193) * t142 + m(2) * (t118 ^ 2 + t119 ^ 2) + m(3) * (t78 ^ 2 + t79 ^ 2) + m(4) * (t41 ^ 2 + t42 ^ 2) + m(5) * (t34 ^ 2 + t35 ^ 2) + m(6) * (t25 ^ 2 + t26 ^ 2) + t192; m(3) * (-t146 * t78 - t148 * t79) + m(4) * (-t146 * t41 - t148 * t42) + m(5) * (-t146 * t34 - t148 * t35) + m(6) * (-t146 * t25 - t148 * t26); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t158) * t159; (-t27 - t185) * t143 + m(4) * (t37 * t41 + t38 * t42) + m(5) * (t14 * t34 + t15 * t35) + m(6) * (t25 * t5 + t26 * t6) + ((-t17 / 0.2e1 - t13 / 0.2e1 - t29 / 0.2e1 - t24 / 0.2e1) * t148 + (t16 / 0.2e1 + t12 / 0.2e1 + t28 / 0.2e1 + t23 / 0.2e1) * t146) * t142 + t155; m(4) * (-t37 * t146 - t148 * t38) + m(5) * (-t14 * t146 - t148 * t15) + m(6) * (-t5 * t146 - t148 * t6); t1 + m(6) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t14 ^ 2 + t15 ^ 2 + t7 ^ 2) + m(4) * (t30 ^ 2 + t37 ^ 2 + t38 ^ 2) + (t185 * t143 - t3) * t143 + (((t102 * t56 + t103 * t58 + t110 * t67 + t111 * t69 + t184 * t176) * t176 + (-t12 - t16 - t23 - t28) * t143) * t146 + (-t2 + ((t104 * t57 + t105 * t59 + t112 * t68 + t113 * t70) * t142 + t182 * t191 * t148) * t148 + (t29 + t24 + t17 + t13) * t143 + (-t102 * t57 - t103 * t59 - t104 * t56 - t105 * t58 - t110 * t68 - t111 * t70 - t112 * t67 - t113 * t69 + (t182 * t146 + t184 * t148) * t142) * t176) * t148) * t142; 0.2e1 * ((t146 * t35 - t148 * t34) * t190 + (t146 * t26 - t148 * t25) * t189) * t142; 0; m(6) * (-t143 * t4 + (t146 * t6 - t148 * t5) * t142) + m(5) * (-t143 * t7 + (-t14 * t148 + t146 * t15) * t142); 0.2e1 * t158 * (t143 ^ 2 + t159 * t191); -t178 + m(6) * (t25 * t32 + t26 * t33) + t155; m(6) * (-t32 * t146 - t148 * t33); m(6) * (t22 * t4 + t32 * t5 + t33 * t6) + t149; m(6) * (-t22 * t143 + (t146 * t33 - t148 * t32) * t142); m(6) * (t22 ^ 2 + t32 ^ 2 + t33 ^ 2) + t149;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;
