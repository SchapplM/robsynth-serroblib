% Calculate joint inertia matrix for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:40
% EndTime: 2019-12-05 15:32:45
% DurationCPUTime: 1.38s
% Computational Cost: add. (2807->206), mult. (3645->328), div. (0->0), fcn. (3886->8), ass. (0->105)
t101 = sin(qJ(4));
t103 = cos(qJ(4));
t96 = qJ(2) + pkin(8);
t92 = sin(t96);
t93 = cos(t96);
t58 = -Icges(6,6) * t93 + (Icges(6,4) * t103 - Icges(6,2) * t101) * t92;
t59 = -Icges(5,6) * t93 + (Icges(5,4) * t103 - Icges(5,2) * t101) * t92;
t160 = -t58 - t59;
t60 = -Icges(6,5) * t93 + (Icges(6,1) * t103 - Icges(6,4) * t101) * t92;
t61 = -Icges(5,5) * t93 + (Icges(5,1) * t103 - Icges(5,4) * t101) * t92;
t159 = -t60 - t61;
t97 = sin(pkin(7));
t94 = t97 ^ 2;
t98 = cos(pkin(7));
t95 = t98 ^ 2;
t134 = t94 + t95;
t158 = Icges(3,3) + Icges(4,3);
t102 = sin(qJ(2));
t104 = cos(qJ(2));
t157 = Icges(3,5) * t104 + Icges(4,5) * t93 - Icges(3,6) * t102 - Icges(4,6) * t92;
t56 = -Icges(6,3) * t93 + (Icges(6,5) * t103 - Icges(6,6) * t101) * t92;
t57 = -Icges(5,3) * t93 + (Icges(5,5) * t103 - Icges(5,6) * t101) * t92;
t156 = (-t57 - t56) * t93;
t143 = t92 * t97;
t128 = Icges(6,3) * t92;
t124 = t98 * t103;
t127 = t97 * t101;
t80 = -t93 * t127 - t124;
t125 = t98 * t101;
t126 = t97 * t103;
t81 = t93 * t126 - t125;
t36 = Icges(6,5) * t81 + Icges(6,6) * t80 + t97 * t128;
t130 = Icges(6,6) * t92;
t40 = Icges(6,4) * t81 + Icges(6,2) * t80 + t97 * t130;
t132 = Icges(6,5) * t92;
t44 = Icges(6,1) * t81 + Icges(6,4) * t80 + t97 * t132;
t12 = t36 * t143 + t80 * t40 + t81 * t44;
t82 = -t93 * t125 + t126;
t83 = t93 * t124 + t127;
t37 = Icges(6,5) * t83 + Icges(6,6) * t82 + t98 * t128;
t41 = Icges(6,4) * t83 + Icges(6,2) * t82 + t98 * t130;
t45 = Icges(6,1) * t83 + Icges(6,4) * t82 + t98 * t132;
t13 = t37 * t143 + t80 * t41 + t81 * t45;
t129 = Icges(5,3) * t92;
t38 = Icges(5,5) * t81 + Icges(5,6) * t80 + t97 * t129;
t131 = Icges(5,6) * t92;
t42 = Icges(5,4) * t81 + Icges(5,2) * t80 + t97 * t131;
t133 = Icges(5,5) * t92;
t46 = Icges(5,1) * t81 + Icges(5,4) * t80 + t97 * t133;
t14 = t38 * t143 + t80 * t42 + t81 * t46;
t39 = Icges(5,5) * t83 + Icges(5,6) * t82 + t98 * t129;
t43 = Icges(5,4) * t83 + Icges(5,2) * t82 + t98 * t131;
t47 = Icges(5,1) * t83 + Icges(5,4) * t82 + t98 * t133;
t15 = t39 * t143 + t80 * t43 + t81 * t47;
t155 = (t159 * t81 + t160 * t80) * t93 + ((t13 + t15) * t98 + (t12 + t14 + t156) * t97) * t92;
t142 = t92 * t98;
t16 = t36 * t142 + t82 * t40 + t83 * t44;
t17 = t37 * t142 + t82 * t41 + t83 * t45;
t18 = t38 * t142 + t82 * t42 + t83 * t46;
t19 = t39 * t142 + t82 * t43 + t83 * t47;
t154 = (t159 * t83 + t160 * t82) * t93 + ((t17 + t19 + t156) * t98 + (t16 + t18) * t97) * t92;
t153 = -t157 * t97 + t158 * t98;
t152 = t157 * t98 + t158 * t97;
t151 = t93 ^ 2;
t146 = t103 * pkin(4);
t144 = pkin(2) * t102;
t105 = qJ(5) * t92 + t146 * t93;
t138 = t81 * rSges(6,1) + t80 * rSges(6,2) + rSges(6,3) * t143 - pkin(4) * t125 + t105 * t97;
t137 = t83 * rSges(6,1) + t82 * rSges(6,2) + rSges(6,3) * t142 + pkin(4) * t127 + t105 * t98;
t136 = (-qJ(5) - rSges(6,3)) * t93 + (rSges(6,1) * t103 - rSges(6,2) * t101 + t146) * t92;
t135 = t134 * t104 * pkin(2);
t123 = -t92 * rSges(4,1) - t93 * rSges(4,2) - t144;
t122 = -t92 * pkin(3) + t93 * pkin(6) - t144;
t121 = t135 + t134 * (pkin(3) * t93 + pkin(6) * t92);
t63 = -t93 * rSges(5,3) + (rSges(5,1) * t103 - rSges(5,2) * t101) * t92;
t120 = t122 - t63;
t109 = t122 - t136;
t87 = t102 * rSges(3,1) + t104 * rSges(3,2);
t65 = t123 * t98;
t64 = t123 * t97;
t52 = t134 * (rSges(3,1) * t104 - rSges(3,2) * t102);
t51 = t83 * rSges(5,1) + t82 * rSges(5,2) + rSges(5,3) * t142;
t49 = t81 * rSges(5,1) + t80 * rSges(5,2) + rSges(5,3) * t143;
t33 = t120 * t98;
t32 = t120 * t97;
t31 = -t63 * t142 - t93 * t51;
t30 = t63 * t143 + t93 * t49;
t29 = t109 * t98;
t28 = t109 * t97;
t27 = t135 + t134 * (rSges(4,1) * t93 - rSges(4,2) * t92);
t26 = (t49 * t98 - t51 * t97) * t92;
t25 = -t93 * t39 + (-t101 * t43 + t103 * t47) * t92;
t24 = -t93 * t38 + (-t101 * t42 + t103 * t46) * t92;
t23 = -t93 * t37 + (-t101 * t41 + t103 * t45) * t92;
t22 = -t93 * t36 + (-t101 * t40 + t103 * t44) * t92;
t21 = -t136 * t142 - t137 * t93;
t20 = t136 * t143 + t138 * t93;
t11 = t97 * t49 + t98 * t51 + t121;
t10 = (-t137 * t97 + t138 * t98) * t92;
t9 = t137 * t98 + t138 * t97 + t121;
t8 = -t18 * t98 + t19 * t97;
t7 = -t16 * t98 + t17 * t97;
t6 = -t14 * t98 + t15 * t97;
t5 = -t12 * t98 + t13 * t97;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t52 + m(4) * t27 + m(5) * t11 + m(6) * t9; m(6) * (t28 ^ 2 + t29 ^ 2 + t9 ^ 2) + m(5) * (t11 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(4) * (t27 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(3) * (t134 * t87 ^ 2 + t52 ^ 2) + (t152 * t94 + t7 + t8) * t97 + (-t5 - t6 + t153 * t95 + (t152 * t98 + t153 * t97) * t97) * t98; 0; m(6) * (-t98 * t28 + t97 * t29) + m(5) * (-t98 * t32 + t97 * t33) + m(4) * (-t98 * t64 + t97 * t65); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t134; m(5) * t26 + m(6) * t10; m(6) * (t10 * t9 + t20 * t29 + t21 * t28) + m(5) * (t26 * t11 + t30 * t33 + t31 * t32) + ((t8 / 0.2e1 + t7 / 0.2e1) * t98 + (t5 / 0.2e1 + t6 / 0.2e1) * t97) * t92 - ((-t22 - t24) * t98 + (t23 + t25) * t97) * t93 / 0.2e1 + t154 * t97 / 0.2e1 - t155 * t98 / 0.2e1; m(5) * (t30 * t97 - t31 * t98) + m(6) * (t20 * t97 - t21 * t98); -t93 * (t151 * t57 + (t25 * t98 + t24 * t97 - (-t101 * t59 + t103 * t61) * t93) * t92) + m(6) * (t10 ^ 2 + t20 ^ 2 + t21 ^ 2) + m(5) * (t26 ^ 2 + t30 ^ 2 + t31 ^ 2) - t93 * (t151 * t56 + (t23 * t98 + t22 * t97 - (-t101 * t58 + t103 * t60) * t93) * t92) + t155 * t143 + t154 * t142; -m(6) * t93; m(6) * (-t93 * t9 + (t28 * t97 + t29 * t98) * t92); 0; m(6) * (-t93 * t10 + (t20 * t98 + t21 * t97) * t92); m(6) * (t134 * t92 ^ 2 + t151);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
