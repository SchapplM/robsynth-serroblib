% Calculate joint inertia matrix for
% S5PRPRP5
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:35
% EndTime: 2019-12-05 15:37:41
% DurationCPUTime: 1.64s
% Computational Cost: add. (2856->243), mult. (4112->390), div. (0->0), fcn. (4462->8), ass. (0->113)
t107 = sin(qJ(2));
t108 = cos(qJ(2));
t100 = pkin(8) + qJ(4);
t96 = sin(t100);
t97 = cos(t100);
t64 = -Icges(6,6) * t108 + (Icges(6,5) * t97 + Icges(6,3) * t96) * t107;
t67 = -Icges(5,6) * t108 + (Icges(5,4) * t97 - Icges(5,2) * t96) * t107;
t149 = -t64 + t67;
t68 = -Icges(6,4) * t108 + (Icges(6,1) * t97 + Icges(6,5) * t96) * t107;
t69 = -Icges(5,5) * t108 + (Icges(5,1) * t97 - Icges(5,4) * t96) * t107;
t148 = -t68 - t69;
t65 = -Icges(5,3) * t108 + (Icges(5,5) * t97 - Icges(5,6) * t96) * t107;
t66 = -Icges(6,2) * t108 + (Icges(6,4) * t97 + Icges(6,6) * t96) * t107;
t147 = (-t65 - t66) * t108;
t146 = rSges(6,1) + pkin(4);
t103 = sin(pkin(7));
t98 = t103 ^ 2;
t105 = cos(pkin(7));
t99 = t105 ^ 2;
t131 = t98 + t99;
t145 = rSges(6,3) + qJ(5);
t126 = t103 * t107;
t125 = t103 * t108;
t78 = t105 * t97 + t96 * t125;
t79 = -t105 * t96 + t97 * t125;
t34 = Icges(6,5) * t79 + Icges(6,6) * t126 + Icges(6,3) * t78;
t38 = Icges(6,4) * t79 + Icges(6,2) * t126 + Icges(6,6) * t78;
t42 = Icges(6,1) * t79 + Icges(6,4) * t126 + Icges(6,5) * t78;
t12 = t126 * t38 + t34 * t78 + t42 * t79;
t123 = t105 * t107;
t122 = t105 * t108;
t80 = -t103 * t97 + t96 * t122;
t81 = t103 * t96 + t97 * t122;
t35 = Icges(6,5) * t81 + Icges(6,6) * t123 + Icges(6,3) * t80;
t39 = Icges(6,4) * t81 + Icges(6,2) * t123 + Icges(6,6) * t80;
t43 = Icges(6,1) * t81 + Icges(6,4) * t123 + Icges(6,5) * t80;
t13 = t126 * t39 + t35 * t78 + t43 * t79;
t36 = Icges(5,5) * t79 - Icges(5,6) * t78 + Icges(5,3) * t126;
t40 = Icges(5,4) * t79 - Icges(5,2) * t78 + Icges(5,6) * t126;
t44 = Icges(5,1) * t79 - Icges(5,4) * t78 + Icges(5,5) * t126;
t14 = t126 * t36 - t40 * t78 + t44 * t79;
t37 = Icges(5,5) * t81 - Icges(5,6) * t80 + Icges(5,3) * t123;
t41 = Icges(5,4) * t81 - Icges(5,2) * t80 + Icges(5,6) * t123;
t45 = Icges(5,1) * t81 - Icges(5,4) * t80 + Icges(5,5) * t123;
t15 = t126 * t37 - t41 * t78 + t45 * t79;
t144 = (t148 * t79 + t149 * t78) * t108 + ((t13 + t15) * t105 + (t12 + t14 + t147) * t103) * t107;
t16 = t123 * t38 + t34 * t80 + t42 * t81;
t17 = t123 * t39 + t35 * t80 + t43 * t81;
t18 = t123 * t36 - t40 * t80 + t44 * t81;
t19 = t123 * t37 - t41 * t80 + t45 * t81;
t143 = (t148 * t81 + t149 * t80) * t108 + ((t17 + t19 + t147) * t105 + (t16 + t18) * t103) * t107;
t142 = t108 ^ 2;
t104 = cos(pkin(8));
t138 = t104 * pkin(3);
t137 = rSges(6,2) * t126 + t145 * t78 + t146 * t79;
t136 = rSges(6,2) * t123 + t145 * t80 + t146 * t81;
t92 = t107 * pkin(2) - t108 * qJ(3);
t135 = pkin(6) * t108 - t138 * t107 - t92;
t134 = -t108 * rSges(6,2) + (t145 * t96 + t146 * t97) * t107;
t102 = sin(pkin(8));
t133 = t108 * rSges(4,3) - (rSges(4,1) * t104 - rSges(4,2) * t102) * t107 - t92;
t132 = t131 * (pkin(2) * t108 + qJ(3) * t107);
t130 = t107 * t96;
t127 = t103 * t102;
t124 = t105 * t102;
t120 = -m(4) - m(5) - m(6);
t71 = -t108 * rSges(5,3) + (rSges(5,1) * t97 - rSges(5,2) * t96) * t107;
t119 = -t71 + t135;
t109 = pkin(6) * t107 + t138 * t108;
t118 = t103 * (-pkin(3) * t124 + t109 * t103) + t105 * (pkin(3) * t127 + t109 * t105) + t132;
t117 = -t134 + t135;
t110 = Icges(3,5) * t108 - Icges(3,6) * t107;
t101 = t107 ^ 2;
t93 = t107 * rSges(3,1) + t108 * rSges(3,2);
t90 = t104 * t122 + t127;
t89 = -t102 * t122 + t103 * t104;
t88 = t104 * t125 - t124;
t87 = -t102 * t125 - t105 * t104;
t73 = Icges(3,3) * t103 + t110 * t105;
t72 = -Icges(3,3) * t105 + t110 * t103;
t62 = t133 * t105;
t61 = t133 * t103;
t58 = Icges(4,1) * t90 + Icges(4,4) * t89 + Icges(4,5) * t123;
t57 = Icges(4,1) * t88 + Icges(4,4) * t87 + Icges(4,5) * t126;
t56 = Icges(4,4) * t90 + Icges(4,2) * t89 + Icges(4,6) * t123;
t55 = Icges(4,4) * t88 + Icges(4,2) * t87 + Icges(4,6) * t126;
t54 = Icges(4,5) * t90 + Icges(4,6) * t89 + Icges(4,3) * t123;
t53 = Icges(4,5) * t88 + Icges(4,6) * t87 + Icges(4,3) * t126;
t50 = t131 * (rSges(3,1) * t108 - rSges(3,2) * t107);
t49 = t81 * rSges(5,1) - t80 * rSges(5,2) + rSges(5,3) * t123;
t47 = t79 * rSges(5,1) - t78 * rSges(5,2) + rSges(5,3) * t126;
t33 = t119 * t105;
t32 = t119 * t103;
t31 = -t108 * t49 - t123 * t71;
t30 = t108 * t47 + t126 * t71;
t29 = t117 * t105;
t28 = t117 * t103;
t27 = (-t103 * t49 + t105 * t47) * t107;
t26 = t103 * (rSges(4,1) * t88 + rSges(4,2) * t87 + rSges(4,3) * t126) + t105 * (rSges(4,1) * t90 + rSges(4,2) * t89 + rSges(4,3) * t123) + t132;
t25 = -t108 * t136 - t123 * t134;
t24 = t108 * t137 + t126 * t134;
t23 = -t108 * t37 + (-t41 * t96 + t45 * t97) * t107;
t22 = -t108 * t36 + (-t40 * t96 + t44 * t97) * t107;
t21 = -t108 * t39 + (t35 * t96 + t43 * t97) * t107;
t20 = -t108 * t38 + (t34 * t96 + t42 * t97) * t107;
t11 = (-t103 * t136 + t105 * t137) * t107;
t10 = t103 * t47 + t105 * t49 + t118;
t9 = t103 * t137 + t105 * t136 + t118;
t8 = t103 * t19 - t105 * t18;
t7 = t103 * t17 - t105 * t16;
t6 = t103 * t15 - t105 * t14;
t5 = t103 * t13 - t105 * t12;
t1 = [m(2) + m(3) - t120; m(3) * t50 + m(4) * t26 + m(5) * t10 + m(6) * t9; m(6) * (t28 ^ 2 + t29 ^ 2 + t9 ^ 2) + m(5) * (t10 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(4) * (t26 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(3) * (t131 * t93 ^ 2 + t50 ^ 2) + (-t99 * t72 - t5 - t6 + (t53 * t126 + t87 * t55 + t88 * t57) * t105) * t105 + (t8 + t7 + t98 * t73 + (t54 * t123 + t89 * t56 + t90 * t58) * t103 + (-t103 * t72 + t105 * t73 - t123 * t53 - t126 * t54 - t55 * t89 - t56 * t87 - t57 * t90 - t58 * t88) * t105) * t103; t120 * t108; m(6) * (-t108 * t9 + (t103 * t28 + t105 * t29) * t107) + m(5) * (-t108 * t10 + (t103 * t32 + t105 * t33) * t107) + m(4) * (-t108 * t26 + (t103 * t61 + t105 * t62) * t107); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t101 * t131 + t142); m(5) * t27 + m(6) * t11; m(6) * (t11 * t9 + t24 * t29 + t25 * t28) + m(5) * (t10 * t27 + t30 * t33 + t31 * t32) + ((t8 / 0.2e1 + t7 / 0.2e1) * t105 + (t6 / 0.2e1 + t5 / 0.2e1) * t103) * t107 + t143 * t103 / 0.2e1 - t144 * t105 / 0.2e1 - ((-t20 - t22) * t105 + (t21 + t23) * t103) * t108 / 0.2e1; m(5) * (-t27 * t108 + (t103 * t31 + t105 * t30) * t107) + m(6) * (-t11 * t108 + (t103 * t25 + t105 * t24) * t107); m(6) * (t11 ^ 2 + t24 ^ 2 + t25 ^ 2) - t108 * (t142 * t65 + (t23 * t105 + t22 * t103 - (-t67 * t96 + t69 * t97) * t108) * t107) - t108 * (t142 * t66 + (t21 * t105 + t20 * t103 - (t64 * t96 + t68 * t97) * t108) * t107) + m(5) * (t27 ^ 2 + t30 ^ 2 + t31 ^ 2) + t144 * t126 + t143 * t123; m(6) * t130; m(6) * (t130 * t9 + t28 * t78 + t29 * t80); m(6) * (t103 * t78 + t105 * t80 - t108 * t96) * t107; m(6) * (t11 * t130 + t24 * t80 + t25 * t78); m(6) * (t101 * t96 ^ 2 + t78 ^ 2 + t80 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
