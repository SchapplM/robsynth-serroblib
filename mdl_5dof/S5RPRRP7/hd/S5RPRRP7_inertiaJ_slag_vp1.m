% Calculate joint inertia matrix for
% S5RPRRP7
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP7_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP7_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:34
% EndTime: 2019-12-31 18:44:38
% DurationCPUTime: 1.43s
% Computational Cost: add. (3546->257), mult. (4324->380), div. (0->0), fcn. (4696->8), ass. (0->129)
t118 = sin(qJ(3));
t174 = Icges(4,5) * t118;
t173 = t174 / 0.2e1;
t169 = rSges(6,3) + qJ(5);
t171 = rSges(6,1) + pkin(4);
t116 = qJ(1) + pkin(8);
t113 = sin(t116);
t114 = cos(t116);
t120 = cos(qJ(4));
t117 = sin(qJ(4));
t121 = cos(qJ(3));
t141 = t117 * t121;
t80 = t113 * t141 + t114 * t120;
t139 = t120 * t121;
t81 = t113 * t139 - t114 * t117;
t172 = -t169 * t80 - t171 * t81;
t85 = -Icges(5,3) * t121 + (Icges(5,5) * t120 - Icges(5,6) * t117) * t118;
t86 = -Icges(6,2) * t121 + (Icges(6,4) * t120 + Icges(6,6) * t117) * t118;
t170 = -t85 - t86;
t145 = t113 * t118;
t42 = Icges(6,5) * t81 + Icges(6,6) * t145 + Icges(6,3) * t80;
t46 = Icges(6,4) * t81 + Icges(6,2) * t145 + Icges(6,6) * t80;
t50 = Icges(6,1) * t81 + Icges(6,4) * t145 + Icges(6,5) * t80;
t11 = t145 * t46 + t42 * t80 + t50 * t81;
t144 = t114 * t118;
t82 = -t113 * t120 + t114 * t141;
t83 = t113 * t117 + t114 * t139;
t43 = Icges(6,5) * t83 + Icges(6,6) * t144 + Icges(6,3) * t82;
t47 = Icges(6,4) * t83 + Icges(6,2) * t144 + Icges(6,6) * t82;
t51 = Icges(6,1) * t83 + Icges(6,4) * t144 + Icges(6,5) * t82;
t12 = t145 * t47 + t43 * t80 + t51 * t81;
t44 = Icges(5,5) * t81 - Icges(5,6) * t80 + Icges(5,3) * t145;
t48 = Icges(5,4) * t81 - Icges(5,2) * t80 + Icges(5,6) * t145;
t52 = Icges(5,1) * t81 - Icges(5,4) * t80 + Icges(5,5) * t145;
t13 = t145 * t44 - t48 * t80 + t52 * t81;
t45 = Icges(5,5) * t83 - Icges(5,6) * t82 + Icges(5,3) * t144;
t49 = Icges(5,4) * t83 - Icges(5,2) * t82 + Icges(5,6) * t144;
t53 = Icges(5,1) * t83 - Icges(5,4) * t82 + Icges(5,5) * t144;
t14 = t145 * t45 - t49 * t80 + t53 * t81;
t84 = -Icges(6,6) * t121 + (Icges(6,5) * t120 + Icges(6,3) * t117) * t118;
t88 = -Icges(6,4) * t121 + (Icges(6,1) * t120 + Icges(6,5) * t117) * t118;
t29 = t145 * t86 + t80 * t84 + t81 * t88;
t87 = -Icges(5,6) * t121 + (Icges(5,4) * t120 - Icges(5,2) * t117) * t118;
t89 = -Icges(5,5) * t121 + (Icges(5,1) * t120 - Icges(5,4) * t117) * t118;
t30 = t145 * t85 - t80 * t87 + t81 * t89;
t168 = (-t29 - t30) * t121 + ((t12 + t14) * t114 + (t11 + t13) * t113) * t118;
t15 = t144 * t46 + t42 * t82 + t50 * t83;
t16 = t144 * t47 + t43 * t82 + t51 * t83;
t17 = t144 * t44 - t48 * t82 + t52 * t83;
t18 = t144 * t45 - t49 * t82 + t53 * t83;
t31 = t144 * t86 + t82 * t84 + t83 * t88;
t32 = t144 * t85 - t82 * t87 + t83 * t89;
t167 = (-t31 - t32) * t121 + ((t16 + t18) * t114 + (t15 + t17) * t113) * t118;
t19 = -t121 * t46 + (t117 * t42 + t120 * t50) * t118;
t21 = -t121 * t44 + (-t117 * t48 + t120 * t52) * t118;
t166 = -t19 - t21;
t20 = -t121 * t47 + (t117 * t43 + t120 * t51) * t118;
t22 = -t121 * t45 + (-t117 * t49 + t120 * t53) * t118;
t165 = t20 + t22;
t142 = t117 * t118;
t164 = t84 * t142 + (t88 + t89) * t118 * t120;
t163 = t113 ^ 2;
t162 = t114 ^ 2;
t99 = rSges(4,1) * t118 + rSges(4,2) * t121;
t161 = m(4) * t99;
t160 = t113 / 0.2e1;
t158 = -t121 / 0.2e1;
t157 = pkin(3) * t121;
t119 = sin(qJ(1));
t156 = t119 * pkin(1);
t155 = t121 * t170 - t87 * t142 + t164;
t154 = rSges(6,2) * t145 - t172;
t153 = rSges(6,2) * t144 + t169 * t82 + t171 * t83;
t143 = t114 * t121;
t138 = pkin(3) * t143 + pkin(7) * t144;
t151 = t163 * (pkin(7) * t118 + t157) + t114 * t138;
t150 = -t121 * rSges(6,2) + (t117 * t169 + t120 * t171) * t118;
t149 = t114 * rSges(4,3);
t105 = pkin(3) * t118 - pkin(7) * t121;
t91 = -t121 * rSges(5,3) + (rSges(5,1) * t120 - rSges(5,2) * t117) * t118;
t148 = -t105 - t91;
t146 = Icges(4,4) * t121;
t137 = -t105 - t150;
t57 = rSges(5,1) * t83 - rSges(5,2) * t82 + rSges(5,3) * t144;
t122 = cos(qJ(1));
t115 = t122 * pkin(1);
t136 = t114 * pkin(2) + t113 * pkin(6) + t115;
t135 = -pkin(2) - t157;
t134 = t114 * pkin(6) - t156;
t133 = -t81 * rSges(5,1) + t80 * rSges(5,2);
t132 = rSges(4,1) * t121 - rSges(4,2) * t118;
t129 = t136 + t138;
t127 = -Icges(4,2) * t118 + t146;
t126 = Icges(4,5) * t121 - Icges(4,6) * t118;
t125 = rSges(4,1) * t143 - rSges(4,2) * t144 + t113 * rSges(4,3);
t124 = t19 / 0.2e1 + t30 / 0.2e1 + t29 / 0.2e1 + t21 / 0.2e1;
t123 = t20 / 0.2e1 + t32 / 0.2e1 + t31 / 0.2e1 + t22 / 0.2e1;
t101 = rSges(2,1) * t122 - rSges(2,2) * t119;
t100 = -rSges(2,1) * t119 - rSges(2,2) * t122;
t96 = Icges(4,6) * t121 + t174;
t93 = rSges(3,1) * t114 - rSges(3,2) * t113 + t115;
t92 = -rSges(3,1) * t113 - rSges(3,2) * t114 - t156;
t65 = Icges(4,3) * t113 + t114 * t126;
t64 = -Icges(4,3) * t114 + t113 * t126;
t63 = t148 * t114;
t62 = t148 * t113;
t61 = t125 + t136;
t60 = t149 + (-pkin(2) - t132) * t113 + t134;
t55 = rSges(5,3) * t145 - t133;
t41 = t137 * t114;
t40 = t137 * t113;
t39 = t114 * t125 + (t113 * t132 - t149) * t113;
t36 = -t121 * t57 - t144 * t91;
t35 = t121 * t55 + t145 * t91;
t34 = t129 + t57;
t33 = ((-rSges(5,3) - pkin(7)) * t118 + t135) * t113 + t133 + t134;
t28 = (-t113 * t57 + t114 * t55) * t118;
t27 = t129 + t153;
t26 = ((-rSges(6,2) - pkin(7)) * t118 + t135) * t113 + t134 + t172;
t25 = -t121 * t153 - t144 * t150;
t24 = t121 * t154 + t145 * t150;
t23 = t113 * t55 + t114 * t57 + t151;
t10 = (-t113 * t153 + t114 * t154) * t118;
t9 = t113 * t154 + t114 * t153 + t151;
t8 = t113 * t18 - t114 * t17;
t7 = t113 * t16 - t114 * t15;
t6 = t113 * t14 - t114 * t13;
t5 = -t11 * t114 + t113 * t12;
t1 = [Icges(2,3) + Icges(3,3) + (Icges(4,1) * t118 - t117 * t87 + t146) * t118 + (Icges(4,4) * t118 + Icges(4,2) * t121 + t170) * t121 + m(6) * (t26 ^ 2 + t27 ^ 2) + m(5) * (t33 ^ 2 + t34 ^ 2) + m(4) * (t60 ^ 2 + t61 ^ 2) + m(3) * (t92 ^ 2 + t93 ^ 2) + m(2) * (t100 ^ 2 + t101 ^ 2) + t164; 0; m(3) + m(4) + m(5) + m(6); m(6) * (t26 * t41 + t27 * t40) + m(5) * (t33 * t63 + t34 * t62) + (t127 * t113 * t158 - t60 * t161 - t124 + (t173 - Icges(4,6) * t158 + t96 / 0.2e1) * t114) * t114 + (-t61 * t161 + t113 * t173 + t121 * (Icges(4,6) * t113 + t114 * t127) / 0.2e1 + t96 * t160 + t123) * t113; m(4) * t39 + m(5) * t23 + m(6) * t9; m(6) * (t40 ^ 2 + t41 ^ 2 + t9 ^ 2) + m(5) * (t23 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(4) * (t39 ^ 2 + (t162 + t163) * t99 ^ 2) + (t163 * t65 + t7 + t8) * t113 + (-t162 * t64 - t5 - t6 + (-t113 * t64 + t114 * t65) * t113) * t114; -t155 * t121 + m(6) * (t24 * t26 + t25 * t27) + m(5) * (t33 * t35 + t34 * t36) + (t113 * t124 + t114 * t123) * t118; m(5) * t28 + m(6) * t10; m(6) * (t10 * t9 + t24 * t41 + t25 * t40) + m(5) * (t23 * t28 + t35 * t63 + t36 * t62) + ((t8 / 0.2e1 + t7 / 0.2e1) * t114 + (t5 / 0.2e1 + t6 / 0.2e1) * t113) * t118 + t167 * t160 - t168 * t114 / 0.2e1 + (t113 * t165 + t114 * t166) * t158; m(6) * (t10 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t28 ^ 2 + t35 ^ 2 + t36 ^ 2) + t155 * t121 ^ 2 + ((-t121 * t165 + t167) * t114 + (t121 * t166 + t168) * t113) * t118; m(6) * (t26 * t82 + t27 * t80); m(6) * t142; m(6) * (t142 * t9 + t40 * t80 + t41 * t82); m(6) * (t10 * t142 + t24 * t82 + t25 * t80); m(6) * (t117 ^ 2 * t118 ^ 2 + t80 ^ 2 + t82 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
