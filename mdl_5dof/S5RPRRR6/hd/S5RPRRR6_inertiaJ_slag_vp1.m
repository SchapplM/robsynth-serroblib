% Calculate joint inertia matrix for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:00:57
% EndTime: 2019-12-31 19:01:00
% DurationCPUTime: 1.21s
% Computational Cost: add. (4774->253), mult. (4019->371), div. (0->0), fcn. (4235->10), ass. (0->134)
t121 = qJ(1) + pkin(9);
t117 = cos(t121);
t122 = qJ(3) + qJ(4);
t119 = cos(t122);
t157 = t117 * t119;
t118 = sin(t122);
t158 = t117 * t118;
t116 = sin(t121);
t126 = cos(qJ(5));
t123 = sin(qJ(5));
t155 = t119 * t123;
t87 = t116 * t126 - t117 * t155;
t154 = t119 * t126;
t88 = t116 * t123 + t117 * t154;
t51 = t88 * rSges(6,1) + t87 * rSges(6,2) + rSges(6,3) * t158;
t179 = pkin(4) * t157 + pkin(8) * t158 + t51;
t115 = t117 ^ 2;
t160 = Icges(5,4) * t119;
t135 = -Icges(5,2) * t118 + t160;
t66 = Icges(5,6) * t116 + t135 * t117;
t161 = Icges(5,4) * t118;
t137 = Icges(5,1) * t119 - t161;
t68 = Icges(5,5) * t116 + t137 * t117;
t142 = -t118 * t66 + t119 * t68;
t65 = -Icges(5,6) * t117 + t135 * t116;
t67 = -Icges(5,5) * t117 + t137 * t116;
t143 = t118 * t65 - t119 * t67;
t133 = Icges(5,5) * t119 - Icges(5,6) * t118;
t63 = -Icges(5,3) * t117 + t133 * t116;
t64 = Icges(5,3) * t116 + t133 * t117;
t159 = t116 * t118;
t85 = -t116 * t155 - t117 * t126;
t86 = t116 * t154 - t117 * t123;
t42 = Icges(6,5) * t86 + Icges(6,6) * t85 + Icges(6,3) * t159;
t44 = Icges(6,4) * t86 + Icges(6,2) * t85 + Icges(6,6) * t159;
t46 = Icges(6,1) * t86 + Icges(6,4) * t85 + Icges(6,5) * t159;
t14 = t42 * t159 + t44 * t85 + t46 * t86;
t43 = Icges(6,5) * t88 + Icges(6,6) * t87 + Icges(6,3) * t158;
t45 = Icges(6,4) * t88 + Icges(6,2) * t87 + Icges(6,6) * t158;
t47 = Icges(6,1) * t88 + Icges(6,4) * t87 + Icges(6,5) * t158;
t15 = t43 * t159 + t45 * t85 + t47 * t86;
t8 = t116 * t15 - t117 * t14;
t178 = -t115 * t63 - (t142 * t116 + (t143 - t64) * t117) * t116 - t8;
t114 = t116 ^ 2;
t177 = t116 / 0.2e1;
t176 = -t117 / 0.2e1;
t16 = t42 * t158 + t44 * t87 + t46 * t88;
t17 = t43 * t158 + t45 * t87 + t47 * t88;
t9 = t116 * t17 - t117 * t16;
t175 = (t114 * t64 + t9 + (t143 * t117 + (t142 - t63) * t116) * t117) * t116;
t124 = sin(qJ(3));
t174 = pkin(3) * t124;
t173 = pkin(4) * t119;
t125 = sin(qJ(1));
t172 = t125 * pkin(1);
t112 = t117 * pkin(6);
t127 = cos(qJ(3));
t113 = pkin(3) * t127 + pkin(2);
t129 = -pkin(7) - pkin(6);
t156 = t117 * t129;
t97 = t117 * t113;
t171 = t116 * (t156 + t112 + (-pkin(2) + t113) * t116) + t117 * (-t117 * pkin(2) + t97 + (-pkin(6) - t129) * t116);
t131 = rSges(5,1) * t157 - rSges(5,2) * t158 + t116 * rSges(5,3);
t144 = rSges(5,1) * t119 - rSges(5,2) * t118;
t38 = t116 * (-t117 * rSges(5,3) + t144 * t116) + t117 * t131;
t82 = -t119 * rSges(6,3) + (rSges(6,1) * t126 - rSges(6,2) * t123) * t118;
t170 = -pkin(4) * t118 + pkin(8) * t119 - t82;
t169 = rSges(4,1) * t127;
t168 = rSges(4,2) * t124;
t167 = t117 * rSges(4,3);
t80 = -Icges(6,6) * t119 + (Icges(6,4) * t126 - Icges(6,2) * t123) * t118;
t166 = t123 * t80;
t20 = -t119 * t42 + (-t123 * t44 + t126 * t46) * t118;
t165 = t20 * t117;
t21 = -t119 * t43 + (-t123 * t45 + t126 * t47) * t118;
t164 = t21 * t116;
t163 = Icges(4,4) * t124;
t162 = Icges(4,4) * t127;
t152 = t116 * rSges(4,3) + t117 * t169;
t151 = t114 + t115;
t95 = rSges(5,1) * t118 + rSges(5,2) * t119;
t150 = -t95 - t174;
t146 = -t86 * rSges(6,1) - t85 * rSges(6,2);
t50 = rSges(6,3) * t159 - t146;
t22 = t116 * t50 + t114 * (pkin(8) * t118 + t173) + t179 * t117;
t79 = -Icges(6,3) * t119 + (Icges(6,5) * t126 - Icges(6,6) * t123) * t118;
t81 = -Icges(6,5) * t119 + (Icges(6,1) * t126 - Icges(6,4) * t123) * t118;
t27 = t79 * t159 + t80 * t85 + t81 * t86;
t3 = -t27 * t119 + (t116 * t14 + t117 * t15) * t118;
t28 = t79 * t158 + t80 * t87 + t81 * t88;
t4 = -t28 * t119 + (t116 * t16 + t117 * t17) * t118;
t149 = t8 * t159 / 0.2e1 + t3 * t176 + t4 * t177 - t119 * (t164 - t165) / 0.2e1 + t9 * t158 / 0.2e1;
t148 = t170 - t174;
t128 = cos(qJ(1));
t120 = t128 * pkin(1);
t147 = -t116 * t129 + t120 + t97;
t145 = -t168 + t169;
t93 = Icges(5,2) * t119 + t161;
t94 = Icges(5,1) * t118 + t160;
t141 = -t118 * t93 + t119 * t94;
t138 = Icges(4,1) * t127 - t163;
t136 = -Icges(4,2) * t124 + t162;
t134 = Icges(4,5) * t127 - Icges(4,6) * t124;
t132 = t178 * t117 + t175;
t92 = Icges(5,5) * t118 + Icges(5,6) * t119;
t130 = -t165 / 0.2e1 + t164 / 0.2e1 + (t116 * t92 + t141 * t117 + t118 * t68 + t119 * t66 + t28) * t177 + (t141 * t116 - t117 * t92 + t118 * t67 + t119 * t65 + t27) * t176;
t108 = rSges(2,1) * t128 - rSges(2,2) * t125;
t107 = -rSges(2,1) * t125 - rSges(2,2) * t128;
t106 = rSges(4,1) * t124 + rSges(4,2) * t127;
t90 = rSges(3,1) * t117 - rSges(3,2) * t116 + t120;
t89 = -rSges(3,1) * t116 - rSges(3,2) * t117 - t172;
t74 = Icges(4,3) * t116 + t134 * t117;
t73 = -Icges(4,3) * t117 + t134 * t116;
t70 = t150 * t117;
t69 = t150 * t116;
t62 = t118 * t126 * t81;
t57 = t116 * pkin(6) + t120 + (pkin(2) - t168) * t117 + t152;
t56 = t167 - t172 + t112 + (-pkin(2) - t145) * t116;
t55 = t170 * t117;
t54 = t170 * t116;
t53 = t131 + t147;
t52 = -t172 + (rSges(5,3) - t129) * t117 + (-t113 - t144) * t116;
t49 = t148 * t117;
t48 = t148 * t116;
t41 = t117 * (-t117 * t168 + t152) + (t145 * t116 - t167) * t116;
t33 = -t118 * t166 - t119 * t79 + t62;
t32 = -t119 * t51 - t82 * t158;
t31 = t119 * t50 + t82 * t159;
t30 = t147 + t179;
t29 = -t172 - t156 + (-t173 - t113 + (-rSges(6,3) - pkin(8)) * t118) * t116 + t146;
t26 = (-t116 * t51 + t117 * t50) * t118;
t23 = t38 + t171;
t11 = t22 + t171;
t1 = [t127 * (Icges(4,2) * t127 + t163) + t124 * (Icges(4,1) * t124 + t162) + Icges(2,3) + Icges(3,3) + t62 + (-t79 + t93) * t119 + (t94 - t166) * t118 + m(6) * (t29 ^ 2 + t30 ^ 2) + m(5) * (t52 ^ 2 + t53 ^ 2) + m(4) * (t56 ^ 2 + t57 ^ 2) + m(3) * (t89 ^ 2 + t90 ^ 2) + m(2) * (t107 ^ 2 + t108 ^ 2); 0; m(3) + m(4) + m(5) + m(6); (t124 * (-Icges(4,5) * t117 + t138 * t116) + t127 * (-Icges(4,6) * t117 + t136 * t116)) * t176 + (t124 * (Icges(4,5) * t116 + t138 * t117) + t127 * (Icges(4,6) * t116 + t136 * t117)) * t177 + m(6) * (t29 * t49 + t30 * t48) + m(5) * (t52 * t70 + t53 * t69) + m(4) * (-t116 * t57 - t117 * t56) * t106 + (t115 / 0.2e1 + t114 / 0.2e1) * (Icges(4,5) * t124 + Icges(4,6) * t127) + t130; m(4) * t41 + m(5) * t23 + m(6) * t11; m(6) * (t11 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t23 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(4) * (t151 * t106 ^ 2 + t41 ^ 2) + t116 * t114 * t74 + t175 + (-t115 * t73 + (-t116 * t73 + t117 * t74) * t116 + t178) * t117; m(6) * (t29 * t55 + t30 * t54) + m(5) * (-t116 * t53 - t117 * t52) * t95 + t130; m(5) * t38 + m(6) * t22; m(6) * (t11 * t22 + t48 * t54 + t49 * t55) + m(5) * (t23 * t38 + (-t116 * t69 - t117 * t70) * t95) + t132; m(5) * (t151 * t95 ^ 2 + t38 ^ 2) + m(6) * (t22 ^ 2 + t54 ^ 2 + t55 ^ 2) + t132; m(6) * (t29 * t31 + t30 * t32) - t33 * t119 + ((t21 / 0.2e1 + t28 / 0.2e1) * t117 + (t20 / 0.2e1 + t27 / 0.2e1) * t116) * t118; m(6) * t26; m(6) * (t11 * t26 + t31 * t49 + t32 * t48) + t149; m(6) * (t22 * t26 + t31 * t55 + t32 * t54) + t149; m(6) * (t26 ^ 2 + t31 ^ 2 + t32 ^ 2) + t119 ^ 2 * t33 + (t117 * t4 + t116 * t3 - t119 * (t116 * t20 + t117 * t21)) * t118;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
