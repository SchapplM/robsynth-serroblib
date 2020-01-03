% Calculate joint inertia matrix for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP10_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP10_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:55
% EndTime: 2019-12-31 18:51:00
% DurationCPUTime: 1.65s
% Computational Cost: add. (3601->280), mult. (4498->419), div. (0->0), fcn. (4777->8), ass. (0->135)
t121 = -qJ(5) - pkin(7);
t180 = rSges(6,3) - t121;
t116 = pkin(8) + qJ(3);
t113 = sin(t116);
t179 = Icges(4,5) * t113;
t178 = t179 / 0.2e1;
t114 = cos(t116);
t123 = sin(qJ(4));
t125 = cos(qJ(4));
t69 = -Icges(6,3) * t114 + (Icges(6,5) * t125 - Icges(6,6) * t123) * t113;
t70 = -Icges(5,3) * t114 + (Icges(5,5) * t125 - Icges(5,6) * t123) * t113;
t177 = -t69 - t70;
t71 = -Icges(6,6) * t114 + (Icges(6,4) * t125 - Icges(6,2) * t123) * t113;
t72 = -Icges(5,6) * t114 + (Icges(5,4) * t125 - Icges(5,2) * t123) * t113;
t176 = (-t71 - t72) * t123;
t124 = sin(qJ(1));
t151 = t113 * t124;
t126 = cos(qJ(1));
t144 = t125 * t126;
t146 = t124 * t123;
t89 = -t114 * t146 - t144;
t145 = t124 * t125;
t147 = t123 * t126;
t90 = t114 * t145 - t147;
t44 = Icges(6,5) * t90 + Icges(6,6) * t89 + Icges(6,3) * t151;
t48 = Icges(6,4) * t90 + Icges(6,2) * t89 + Icges(6,6) * t151;
t52 = Icges(6,1) * t90 + Icges(6,4) * t89 + Icges(6,5) * t151;
t11 = t44 * t151 + t48 * t89 + t52 * t90;
t149 = t113 * t126;
t91 = -t114 * t147 + t145;
t92 = t114 * t144 + t146;
t45 = Icges(6,5) * t92 + Icges(6,6) * t91 + Icges(6,3) * t149;
t49 = Icges(6,4) * t92 + Icges(6,2) * t91 + Icges(6,6) * t149;
t53 = Icges(6,1) * t92 + Icges(6,4) * t91 + Icges(6,5) * t149;
t12 = t45 * t151 + t49 * t89 + t53 * t90;
t46 = Icges(5,5) * t90 + Icges(5,6) * t89 + Icges(5,3) * t151;
t50 = Icges(5,4) * t90 + Icges(5,2) * t89 + Icges(5,6) * t151;
t54 = Icges(5,1) * t90 + Icges(5,4) * t89 + Icges(5,5) * t151;
t13 = t46 * t151 + t50 * t89 + t54 * t90;
t47 = Icges(5,5) * t92 + Icges(5,6) * t91 + Icges(5,3) * t149;
t51 = Icges(5,4) * t92 + Icges(5,2) * t91 + Icges(5,6) * t149;
t55 = Icges(5,1) * t92 + Icges(5,4) * t91 + Icges(5,5) * t149;
t14 = t47 * t151 + t51 * t89 + t55 * t90;
t73 = -Icges(6,5) * t114 + (Icges(6,1) * t125 - Icges(6,4) * t123) * t113;
t26 = t69 * t151 + t71 * t89 + t73 * t90;
t74 = -Icges(5,5) * t114 + (Icges(5,1) * t125 - Icges(5,4) * t123) * t113;
t27 = t70 * t151 + t72 * t89 + t74 * t90;
t175 = (-t26 - t27) * t114 + ((t12 + t14) * t126 + (t11 + t13) * t124) * t113;
t15 = t44 * t149 + t91 * t48 + t92 * t52;
t16 = t45 * t149 + t91 * t49 + t92 * t53;
t17 = t46 * t149 + t91 * t50 + t92 * t54;
t18 = t47 * t149 + t91 * t51 + t92 * t55;
t28 = t69 * t149 + t91 * t71 + t92 * t73;
t29 = t70 * t149 + t91 * t72 + t92 * t74;
t174 = (-t28 - t29) * t114 + ((t16 + t18) * t126 + (t15 + t17) * t124) * t113;
t21 = -t114 * t44 + (-t123 * t48 + t125 * t52) * t113;
t23 = -t114 * t46 + (-t123 * t50 + t125 * t54) * t113;
t173 = -t21 - t23;
t22 = -t114 * t45 + (-t123 * t49 + t125 * t53) * t113;
t24 = -t114 * t47 + (-t123 * t51 + t125 * t55) * t113;
t172 = t22 + t24;
t171 = (t73 + t74) * t113 * t125;
t111 = pkin(4) * t125 + pkin(3);
t148 = t114 * t126;
t170 = t92 * rSges(6,1) + t91 * rSges(6,2) + pkin(4) * t146 + t111 * t148 + t180 * t149;
t169 = t114 ^ 2;
t117 = t124 ^ 2;
t118 = t126 ^ 2;
t97 = rSges(4,1) * t113 + rSges(4,2) * t114;
t168 = m(4) * t97;
t167 = -t114 / 0.2e1;
t166 = t124 / 0.2e1;
t164 = pkin(3) * t114;
t163 = -pkin(3) + t111;
t162 = pkin(7) + t121;
t161 = t113 * t176 + t177 * t114 + t171;
t138 = -t90 * rSges(6,1) - t89 * rSges(6,2);
t160 = -pkin(4) * t147 + (-t162 * t113 + t163 * t114) * t124 + rSges(6,3) * t151 - t138;
t143 = pkin(3) * t148 + pkin(7) * t149;
t159 = -t143 + t170;
t158 = (t162 - rSges(6,3)) * t114 + (rSges(6,1) * t125 - rSges(6,2) * t123 + t163) * t113;
t76 = -rSges(5,3) * t114 + (rSges(5,1) * t125 - rSges(5,2) * t123) * t113;
t99 = pkin(3) * t113 - pkin(7) * t114;
t157 = -t76 - t99;
t156 = t117 * (pkin(7) * t113 + t164) + t126 * t143;
t155 = rSges(3,3) + qJ(2);
t153 = Icges(4,4) * t114;
t142 = t117 + t118;
t141 = -t99 - t158;
t59 = t92 * rSges(5,1) + t91 * rSges(5,2) + rSges(5,3) * t149;
t120 = cos(pkin(8));
t110 = pkin(2) * t120 + pkin(1);
t122 = -pkin(6) - qJ(2);
t140 = t126 * t110 - t124 * t122;
t139 = -t90 * rSges(5,1) - t89 * rSges(5,2);
t137 = rSges(4,1) * t114 - rSges(4,2) * t113;
t133 = -Icges(4,2) * t113 + t153;
t132 = Icges(4,5) * t114 - Icges(4,6) * t113;
t131 = rSges(4,1) * t148 - rSges(4,2) * t149 + t124 * rSges(4,3);
t119 = sin(pkin(8));
t129 = rSges(3,1) * t120 - rSges(3,2) * t119 + pkin(1);
t128 = t21 / 0.2e1 + t27 / 0.2e1 + t26 / 0.2e1 + t23 / 0.2e1;
t127 = t24 / 0.2e1 + t22 / 0.2e1 + t29 / 0.2e1 + t28 / 0.2e1;
t101 = rSges(2,1) * t126 - t124 * rSges(2,2);
t100 = -t124 * rSges(2,1) - rSges(2,2) * t126;
t94 = Icges(4,6) * t114 + t179;
t78 = Icges(4,3) * t124 + t132 * t126;
t77 = -Icges(4,3) * t126 + t132 * t124;
t68 = t155 * t124 + t129 * t126;
t67 = -t129 * t124 + t155 * t126;
t63 = t131 + t140;
t62 = (rSges(4,3) - t122) * t126 + (-t110 - t137) * t124;
t61 = t157 * t126;
t60 = t157 * t124;
t57 = rSges(5,3) * t151 - t139;
t41 = t126 * t131 + (-t126 * rSges(4,3) + t137 * t124) * t124;
t40 = t140 + t59 + t143;
t39 = -t122 * t126 + (-t164 - t110 + (-rSges(5,3) - pkin(7)) * t113) * t124 + t139;
t38 = t141 * t126;
t37 = t141 * t124;
t36 = -t114 * t59 - t76 * t149;
t35 = t114 * t57 + t76 * t151;
t34 = t140 + t170;
t33 = (pkin(4) * t123 - t122) * t126 + (-t111 * t114 - t180 * t113 - t110) * t124 + t138;
t30 = (-t124 * t59 + t126 * t57) * t113;
t25 = t124 * t57 + t126 * t59 + t156;
t20 = -t159 * t114 - t158 * t149;
t19 = t160 * t114 + t158 * t151;
t10 = (-t159 * t124 + t160 * t126) * t113;
t9 = t160 * t124 + t159 * t126 + t156;
t8 = t18 * t124 - t126 * t17;
t7 = t16 * t124 - t126 * t15;
t6 = t14 * t124 - t126 * t13;
t5 = -t11 * t126 + t12 * t124;
t1 = [Icges(3,2) * t120 ^ 2 + Icges(2,3) + (Icges(3,1) * t119 + 0.2e1 * Icges(3,4) * t120) * t119 + (Icges(4,4) * t113 + Icges(4,2) * t114 + t177) * t114 + (Icges(4,1) * t113 + t153 + t176) * t113 + m(6) * (t33 ^ 2 + t34 ^ 2) + m(5) * (t39 ^ 2 + t40 ^ 2) + m(4) * (t62 ^ 2 + t63 ^ 2) + m(3) * (t67 ^ 2 + t68 ^ 2) + m(2) * (t100 ^ 2 + t101 ^ 2) + t171; m(6) * (t124 * t33 - t126 * t34) + m(5) * (t124 * t39 - t126 * t40) + m(4) * (t124 * t62 - t126 * t63) + m(3) * (t124 * t67 - t126 * t68); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t142; m(6) * (t33 * t38 + t34 * t37) + m(5) * (t39 * t61 + t40 * t60) + (t133 * t124 * t167 - t62 * t168 - t128 + (t178 - Icges(4,6) * t167 + t94 / 0.2e1) * t126) * t126 + (-t63 * t168 + t124 * t178 + t114 * (Icges(4,6) * t124 + t133 * t126) / 0.2e1 + t94 * t166 + t127) * t124; m(5) * (t61 * t124 - t126 * t60) + m(6) * (t38 * t124 - t126 * t37); m(6) * (t37 ^ 2 + t38 ^ 2 + t9 ^ 2) + m(5) * (t25 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(4) * (t142 * t97 ^ 2 + t41 ^ 2) + (t117 * t78 + t7 + t8) * t124 + (-t118 * t77 - t5 - t6 + (-t124 * t77 + t126 * t78) * t124) * t126; -t161 * t114 + m(6) * (t19 * t33 + t20 * t34) + m(5) * (t35 * t39 + t36 * t40) + (t128 * t124 + t127 * t126) * t113; m(5) * (t35 * t124 - t126 * t36) + m(6) * (t19 * t124 - t126 * t20); m(6) * (t10 * t9 + t19 * t38 + t20 * t37) + m(5) * (t25 * t30 + t35 * t61 + t36 * t60) + ((t8 / 0.2e1 + t7 / 0.2e1) * t126 + (t5 / 0.2e1 + t6 / 0.2e1) * t124) * t113 + (t172 * t124 + t173 * t126) * t167 + t174 * t166 - t175 * t126 / 0.2e1; m(6) * (t10 ^ 2 + t19 ^ 2 + t20 ^ 2) + m(5) * (t30 ^ 2 + t35 ^ 2 + t36 ^ 2) + t161 * t169 + (t174 * t126 + t175 * t124 + (t173 * t124 - t172 * t126) * t114) * t113; m(6) * (t124 * t34 + t126 * t33) * t113; 0; m(6) * (-t114 * t9 + (t124 * t37 + t126 * t38) * t113); m(6) * (-t114 * t10 + (t124 * t20 + t126 * t19) * t113); m(6) * (t142 * t113 ^ 2 + t169);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
