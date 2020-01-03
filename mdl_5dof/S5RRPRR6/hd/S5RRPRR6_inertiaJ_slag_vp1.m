% Calculate joint inertia matrix for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:30
% EndTime: 2020-01-03 12:05:34
% DurationCPUTime: 1.35s
% Computational Cost: add. (4216->253), mult. (4188->370), div. (0->0), fcn. (4406->10), ass. (0->121)
t125 = sin(pkin(9));
t126 = cos(pkin(9));
t127 = sin(qJ(4));
t129 = cos(qJ(4));
t124 = qJ(1) + qJ(2);
t118 = sin(t124);
t153 = t118 * t125;
t120 = cos(t124);
t146 = t126 * t127;
t83 = -t118 * t146 - t120 * t129;
t145 = t126 * t129;
t148 = t120 * t127;
t84 = t118 * t145 - t148;
t50 = Icges(5,5) * t84 + Icges(5,6) * t83 + Icges(5,3) * t153;
t52 = Icges(5,4) * t84 + Icges(5,2) * t83 + Icges(5,6) * t153;
t54 = Icges(5,1) * t84 + Icges(5,4) * t83 + Icges(5,5) * t153;
t79 = -Icges(5,3) * t126 + (Icges(5,5) * t129 - Icges(5,6) * t127) * t125;
t80 = -Icges(5,6) * t126 + (Icges(5,4) * t129 - Icges(5,2) * t127) * t125;
t81 = -Icges(5,5) * t126 + (Icges(5,1) * t129 - Icges(5,4) * t127) * t125;
t161 = -t126 * t50 + (-t127 * t52 + t129 * t54) * t125 + t79 * t153 + t80 * t83 + t81 * t84;
t150 = t120 * t125;
t85 = -t118 * t129 + t120 * t146;
t151 = t118 * t127;
t86 = -t120 * t145 - t151;
t51 = Icges(5,5) * t86 + Icges(5,6) * t85 - Icges(5,3) * t150;
t53 = Icges(5,4) * t86 + Icges(5,2) * t85 - Icges(5,6) * t150;
t55 = Icges(5,1) * t86 + Icges(5,4) * t85 - Icges(5,5) * t150;
t160 = -t126 * t51 + (-t127 * t53 + t129 * t55) * t125 - t79 * t150 + t80 * t85 + t81 * t86;
t159 = t118 ^ 2;
t158 = t120 ^ 2;
t123 = qJ(4) + qJ(5);
t117 = sin(t123);
t119 = cos(t123);
t152 = t118 * t126;
t73 = -t117 * t152 - t120 * t119;
t74 = -t120 * t117 + t119 * t152;
t48 = t74 * rSges(6,1) + t73 * rSges(6,2) + rSges(6,3) * t153;
t149 = t120 * t126;
t75 = t117 * t149 - t118 * t119;
t76 = -t118 * t117 - t119 * t149;
t49 = t76 * rSges(6,1) + t75 * rSges(6,2) - rSges(6,3) * t150;
t17 = t48 * t150 + t49 * t153;
t70 = -Icges(6,6) * t126 + (Icges(6,4) * t119 - Icges(6,2) * t117) * t125;
t157 = t117 * t70;
t156 = t127 * t80;
t71 = -Icges(6,5) * t126 + (Icges(6,1) * t119 - Icges(6,4) * t117) * t125;
t64 = t125 * t119 * t71;
t69 = -Icges(6,3) * t126 + (Icges(6,5) * t119 - Icges(6,6) * t117) * t125;
t31 = -t125 * t157 - t126 * t69 + t64;
t155 = t31 * t126;
t154 = t118 * t120;
t131 = -pkin(8) - pkin(7);
t147 = t125 * t131;
t144 = -pkin(3) * t152 - pkin(7) * t153;
t143 = pkin(3) * t149 + pkin(7) * t150;
t142 = t120 * pkin(2) + t118 * qJ(3);
t93 = t118 * rSges(3,1) + t120 * rSges(3,2);
t58 = t84 * rSges(5,1) + t83 * rSges(5,2) + rSges(5,3) * t153;
t141 = t153 / 0.2e1;
t140 = -t150 / 0.2e1;
t18 = t69 * t153 + t70 * t73 + t71 * t74;
t19 = -t69 * t150 + t70 * t75 + t71 * t76;
t42 = Icges(6,5) * t74 + Icges(6,6) * t73 + Icges(6,3) * t153;
t44 = Icges(6,4) * t74 + Icges(6,2) * t73 + Icges(6,6) * t153;
t46 = Icges(6,1) * t74 + Icges(6,4) * t73 + Icges(6,5) * t153;
t7 = -t126 * t42 + (-t117 * t44 + t119 * t46) * t125;
t43 = Icges(6,5) * t76 + Icges(6,6) * t75 - Icges(6,3) * t150;
t45 = Icges(6,4) * t76 + Icges(6,2) * t75 - Icges(6,6) * t150;
t47 = Icges(6,1) * t76 + Icges(6,4) * t75 - Icges(6,5) * t150;
t8 = -t126 * t43 + (-t117 * t45 + t119 * t47) * t125;
t139 = (t18 + t7) * t141 + (t19 + t8) * t140;
t116 = t129 * pkin(4) + pkin(3);
t72 = -t126 * rSges(6,3) + (rSges(6,1) * t119 - rSges(6,2) * t117) * t125;
t138 = (-(pkin(7) + t131) * t126 - (-pkin(3) + t116) * t125 - t72) * t125;
t94 = t120 * rSges(3,1) - t118 * rSges(3,2);
t137 = t116 * t152 - t118 * t147;
t136 = -pkin(4) * t151 - t116 * t149 + t120 * t147;
t59 = t86 * rSges(5,1) + t85 * rSges(5,2) - rSges(5,3) * t150;
t135 = t139 - t155;
t1 = (-t18 * t126 + (t42 * t153 + t44 * t73 + t46 * t74) * t153 - (t43 * t153 + t45 * t73 + t47 * t74) * t150) * t153;
t2 = -t19 * t126 + (-t42 * t150 + t44 * t75 + t46 * t76) * t153 - (-t43 * t150 + t45 * t75 + t47 * t76) * t150;
t3 = -t155 + (t118 * t7 - t120 * t8) * t125;
t134 = -t126 * t3 - t2 * t150 + t1;
t63 = rSges(4,1) * t149 - rSges(4,2) * t150 + t118 * rSges(4,3) + t142;
t114 = t118 * pkin(2);
t36 = -t120 * qJ(3) + t114 - t144 + t58;
t62 = -rSges(4,2) * t153 + rSges(4,1) * t152 + t114 + (-rSges(4,3) - qJ(3)) * t120;
t37 = -t59 + t142 + t143;
t65 = t125 * t129 * t81;
t38 = -t125 * t156 - t126 * t79 + t65;
t133 = (-t38 - t31) * t126 + t139 + t161 * t141 + t160 * t140;
t28 = -t136 - t49 + t142;
t27 = t114 + (-pkin(4) * t127 - qJ(3)) * t120 + t137 + t48;
t132 = Icges(3,3) + t64 + t65 + (Icges(4,2) * t126 - t69 - t79) * t126 + (Icges(4,1) * t125 + 0.2e1 * Icges(4,4) * t126 - t156 - t157) * t125;
t130 = cos(qJ(1));
t128 = sin(qJ(1));
t122 = t130 * pkin(1);
t121 = t128 * pkin(1);
t97 = t130 * rSges(2,1) - t128 * rSges(2,2);
t96 = t128 * rSges(2,1) + t130 * rSges(2,2);
t88 = t122 + t94;
t87 = t121 + t93;
t82 = -t126 * rSges(5,3) + (rSges(5,1) * t129 - rSges(5,2) * t127) * t125;
t61 = t122 + t63;
t60 = t121 + t62;
t57 = t136 + t143;
t56 = -pkin(4) * t148 + t137 + t144;
t41 = t126 * t49;
t35 = t126 * t59 - t82 * t150;
t34 = -t126 * t58 - t82 * t153;
t33 = t122 + t37;
t32 = t121 + t36;
t30 = -t72 * t150 + t41;
t29 = -t126 * t48 - t72 * t153;
t26 = t122 + t28;
t25 = t121 + t27;
t20 = (t118 * t59 + t120 * t58) * t125;
t10 = t120 * t138 + t126 * t57 + t41;
t9 = (-t48 - t56) * t126 + t118 * t138;
t4 = (t118 * t57 + t120 * t56) * t125 + t17;
t5 = [Icges(2,3) + m(2) * (t96 ^ 2 + t97 ^ 2) + m(3) * (t87 ^ 2 + t88 ^ 2) + m(4) * (t60 ^ 2 + t61 ^ 2) + m(5) * (t32 ^ 2 + t33 ^ 2) + m(6) * (t25 ^ 2 + t26 ^ 2) + t132; m(3) * (t87 * t93 + t88 * t94) + m(4) * (t60 * t62 + t61 * t63) + m(5) * (t32 * t36 + t33 * t37) + m(6) * (t25 * t27 + t26 * t28) + t132; m(6) * (t27 ^ 2 + t28 ^ 2) + m(5) * (t36 ^ 2 + t37 ^ 2) + m(4) * (t62 ^ 2 + t63 ^ 2) + m(3) * (t93 ^ 2 + t94 ^ 2) + t132; m(4) * (-t118 * t60 - t120 * t61) + m(5) * (-t118 * t32 - t120 * t33) + m(6) * (-t118 * t25 - t120 * t26); m(6) * (-t118 * t27 - t120 * t28) + m(5) * (-t118 * t36 - t120 * t37) + m(4) * (-t118 * t62 - t120 * t63); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t158 + t159); m(5) * (t32 * t34 + t33 * t35) + m(6) * (t10 * t26 + t25 * t9) + t133; m(6) * (t10 * t28 + t27 * t9) + m(5) * (t34 * t36 + t35 * t37) + t133; m(5) * (-t118 * t34 - t120 * t35) + m(6) * (-t10 * t120 - t118 * t9); t1 + (t38 * t126 - t3) * t126 + (-t120 * t2 + (t118 * ((t52 * t83 + t54 * t84) * t118 - (t53 * t83 + t55 * t84) * t120) - t120 * ((t52 * t85 + t54 * t86) * t118 - (t53 * t85 + t55 * t86) * t120) + (t118 * (-t51 * t154 + t159 * t50) - t120 * (-t50 * t154 + t158 * t51)) * t125) * t125 + (-t161 * t118 + t160 * t120) * t126) * t125 + m(5) * (t20 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(6) * (t10 ^ 2 + t4 ^ 2 + t9 ^ 2); m(6) * (t25 * t29 + t26 * t30) + t135; m(6) * (t27 * t29 + t28 * t30) + t135; m(6) * (-t118 * t29 - t120 * t30); m(6) * (t10 * t30 + t17 * t4 + t29 * t9) + t134; m(6) * (t17 ^ 2 + t29 ^ 2 + t30 ^ 2) + t134;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
