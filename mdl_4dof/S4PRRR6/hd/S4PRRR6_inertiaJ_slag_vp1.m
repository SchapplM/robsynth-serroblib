% Calculate joint inertia matrix for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR6_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR6_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:40
% EndTime: 2019-12-31 16:34:43
% DurationCPUTime: 0.98s
% Computational Cost: add. (2818->198), mult. (4462->331), div. (0->0), fcn. (4933->8), ass. (0->116)
t99 = cos(pkin(7));
t142 = t99 ^ 2;
t98 = sin(pkin(7));
t144 = t98 ^ 2 + t142;
t103 = cos(qJ(2));
t141 = t103 ^ 2;
t140 = t98 / 0.2e1;
t139 = -t99 / 0.2e1;
t138 = -t103 / 0.2e1;
t102 = cos(qJ(3));
t137 = t102 * pkin(3);
t101 = sin(qJ(2));
t105 = pkin(6) * t101 + t137 * t103;
t100 = sin(qJ(3));
t125 = t98 * t100;
t130 = t101 * t99;
t126 = t103 * t99;
t97 = qJ(3) + qJ(4);
t95 = sin(t97);
t96 = cos(t97);
t79 = -t95 * t126 + t98 * t96;
t80 = t96 * t126 + t98 * t95;
t51 = t80 * rSges(5,1) + t79 * rSges(5,2) + rSges(5,3) * t130;
t135 = pkin(3) * t125 + t105 * t99 + t51;
t131 = t101 * t98;
t127 = t103 * t98;
t77 = -t95 * t127 - t99 * t96;
t78 = t96 * t127 - t99 * t95;
t50 = t78 * rSges(5,1) + t77 * rSges(5,2) + rSges(5,3) * t131;
t70 = -t103 * rSges(5,3) + (rSges(5,1) * t96 - rSges(5,2) * t95) * t101;
t36 = t103 * t50 + t70 * t131;
t66 = -pkin(6) * t103 + t137 * t101;
t134 = -t66 - t70;
t84 = -t103 * rSges(4,3) + (rSges(4,1) * t102 - rSges(4,2) * t100) * t101;
t93 = t101 * pkin(2) - t103 * pkin(5);
t133 = -t84 - t93;
t132 = t144 * (pkin(2) * t103 + pkin(5) * t101);
t67 = -Icges(5,3) * t103 + (Icges(5,5) * t96 - Icges(5,6) * t95) * t101;
t129 = t103 * t67;
t81 = -Icges(4,3) * t103 + (Icges(4,5) * t102 - Icges(4,6) * t100) * t101;
t128 = t103 * t81;
t124 = t99 * t100;
t123 = Icges(4,5) * t101;
t122 = Icges(5,5) * t101;
t121 = Icges(4,6) * t101;
t120 = Icges(5,6) * t101;
t119 = Icges(4,3) * t101;
t118 = Icges(5,3) * t101;
t117 = t100 * t103;
t116 = t102 * t103;
t115 = -t93 + t134;
t44 = Icges(5,5) * t78 + Icges(5,6) * t77 + t98 * t118;
t46 = Icges(5,4) * t78 + Icges(5,2) * t77 + t98 * t120;
t48 = Icges(5,1) * t78 + Icges(5,4) * t77 + t98 * t122;
t23 = -t103 * t44 + (-t46 * t95 + t48 * t96) * t101;
t45 = Icges(5,5) * t80 + Icges(5,6) * t79 + t99 * t118;
t47 = Icges(5,4) * t80 + Icges(5,2) * t79 + t99 * t120;
t49 = Icges(5,1) * t80 + Icges(5,4) * t79 + t99 * t122;
t24 = -t103 * t45 + (-t47 * t95 + t49 * t96) * t101;
t19 = t44 * t131 + t77 * t46 + t78 * t48;
t20 = t45 * t131 + t77 * t47 + t78 * t49;
t68 = -Icges(5,6) * t103 + (Icges(5,4) * t96 - Icges(5,2) * t95) * t101;
t69 = -Icges(5,5) * t103 + (Icges(5,1) * t96 - Icges(5,4) * t95) * t101;
t5 = -(t77 * t68 + t78 * t69) * t103 + (t20 * t99 + (t19 - t129) * t98) * t101;
t21 = t44 * t130 + t79 * t46 + t80 * t48;
t22 = t45 * t130 + t79 * t47 + t80 * t49;
t6 = -(t79 * t68 + t80 * t69) * t103 + (t21 * t98 + (t22 - t129) * t99) * t101;
t114 = -t103 * (t141 * t67 + (t24 * t99 + t23 * t98 - (-t68 * t95 + t69 * t96) * t103) * t101) + t6 * t130 + t5 * t131;
t12 = -t19 * t99 + t20 * t98;
t13 = -t21 * t99 + t22 * t98;
t113 = t12 * t131 / 0.2e1 + t5 * t139 + t6 * t140 + t13 * t130 / 0.2e1 + (-t23 * t99 + t24 * t98) * t138;
t106 = Icges(3,5) * t103 - Icges(3,6) * t101;
t92 = t101 * rSges(3,1) + t103 * rSges(3,2);
t90 = t99 * t116 + t125;
t89 = t98 * t102 - t99 * t117;
t88 = t98 * t116 - t124;
t87 = -t99 * t102 - t98 * t117;
t83 = -Icges(4,5) * t103 + (Icges(4,1) * t102 - Icges(4,4) * t100) * t101;
t82 = -Icges(4,6) * t103 + (Icges(4,4) * t102 - Icges(4,2) * t100) * t101;
t71 = -Icges(3,3) * t99 + t106 * t98;
t64 = t133 * t99;
t63 = t133 * t98;
t61 = -pkin(3) * t124 + t105 * t98;
t60 = t90 * rSges(4,1) + t89 * rSges(4,2) + rSges(4,3) * t130;
t59 = t88 * rSges(4,1) + t87 * rSges(4,2) + rSges(4,3) * t131;
t58 = Icges(4,1) * t90 + Icges(4,4) * t89 + t99 * t123;
t57 = Icges(4,1) * t88 + Icges(4,4) * t87 + t98 * t123;
t56 = Icges(4,4) * t90 + Icges(4,2) * t89 + t99 * t121;
t55 = Icges(4,4) * t88 + Icges(4,2) * t87 + t98 * t121;
t54 = Icges(4,5) * t90 + Icges(4,6) * t89 + t99 * t119;
t53 = Icges(4,5) * t88 + Icges(4,6) * t87 + t98 * t119;
t52 = t144 * (rSges(3,1) * t103 - rSges(3,2) * t101);
t42 = t50 * t130;
t41 = t115 * t99;
t40 = t115 * t98;
t39 = -t103 * t60 - t84 * t130;
t38 = t103 * t59 + t84 * t131;
t37 = -t103 * t51 - t70 * t130;
t35 = (t59 * t99 - t60 * t98) * t101;
t34 = -t51 * t131 + t42;
t33 = t98 * t59 + t99 * t60 + t132;
t32 = -t103 * t54 + (-t100 * t56 + t102 * t58) * t101;
t31 = -t103 * t53 + (-t100 * t55 + t102 * t57) * t101;
t30 = -t135 * t103 + t134 * t130;
t29 = t103 * t61 + t66 * t131 + t36;
t28 = t54 * t130 + t89 * t56 + t90 * t58;
t27 = t53 * t130 + t89 * t55 + t90 * t57;
t26 = t54 * t131 + t87 * t56 + t88 * t58;
t25 = t53 * t131 + t87 * t55 + t88 * t57;
t18 = t42 + (-t135 * t98 + t61 * t99) * t101;
t17 = t135 * t99 + (t50 + t61) * t98 + t132;
t16 = -t27 * t99 + t28 * t98;
t15 = -t25 * t99 + t26 * t98;
t8 = -(t89 * t82 + t90 * t83) * t103 + (t27 * t98 + (t28 - t128) * t99) * t101;
t7 = -(t87 * t82 + t88 * t83) * t103 + (t26 * t99 + (t25 - t128) * t98) * t101;
t1 = [m(2) + m(3) + m(4) + m(5); m(3) * t52 + m(4) * t33 + m(5) * t17; m(5) * (t17 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(4) * (t33 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(3) * (t144 * t92 ^ 2 + t52 ^ 2) + (-t142 * t71 - t12 - t15) * t99 + (-t98 * t71 * t99 + t13 + t16 + t144 * (Icges(3,3) * t98 + t106 * t99)) * t98; m(4) * t35 + m(5) * t18; (-t31 * t99 + t32 * t98) * t138 + t8 * t140 + t7 * t139 + (t99 * t16 / 0.2e1 + t15 * t140) * t101 + m(5) * (t18 * t17 + t29 * t41 + t30 * t40) + m(4) * (t35 * t33 + t38 * t64 + t39 * t63) + t113; t8 * t130 + t7 * t131 + m(5) * (t18 ^ 2 + t29 ^ 2 + t30 ^ 2) - t103 * (t141 * t81 + (t32 * t99 + t31 * t98 - (-t100 * t82 + t102 * t83) * t103) * t101) + m(4) * (t35 ^ 2 + t38 ^ 2 + t39 ^ 2) + t114; m(5) * t34; m(5) * (t34 * t17 + t36 * t41 + t37 * t40) + t113; m(5) * (t34 * t18 + t36 * t29 + t37 * t30) + t114; m(5) * (t34 ^ 2 + t36 ^ 2 + t37 ^ 2) + t114;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
