% Calculate joint inertia matrix for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:09
% EndTime: 2019-12-05 18:58:11
% DurationCPUTime: 0.68s
% Computational Cost: add. (3060->184), mult. (1852->258), div. (0->0), fcn. (1594->10), ass. (0->102)
t89 = qJ(1) + qJ(2);
t87 = qJ(3) + t89;
t80 = sin(t87);
t81 = cos(t87);
t136 = t80 * t81;
t77 = t80 ^ 2;
t78 = t81 ^ 2;
t135 = t80 / 0.2e1;
t134 = t81 / 0.2e1;
t88 = qJ(4) + qJ(5);
t83 = sin(t88);
t85 = cos(t88);
t97 = Icges(6,5) * t85 - Icges(6,6) * t83;
t27 = Icges(6,3) * t81 - t80 * t97;
t28 = Icges(6,3) * t80 + t97 * t81;
t133 = t80 * (t27 * t136 + t77 * t28) + t81 * (t28 * t136 + t78 * t27);
t90 = sin(qJ(4));
t92 = cos(qJ(4));
t69 = t90 * rSges(5,1) + t92 * rSges(5,2);
t132 = m(5) * t69;
t56 = t83 * rSges(6,1) + t85 * rSges(6,2);
t131 = m(6) * t56;
t84 = sin(t89);
t130 = pkin(2) * t84;
t86 = cos(t89);
t129 = pkin(2) * t86;
t91 = sin(qJ(1));
t128 = t91 * pkin(1);
t93 = cos(qJ(1));
t127 = t93 * pkin(1);
t82 = t92 * pkin(4) + pkin(3);
t126 = -pkin(3) + t82;
t125 = rSges(5,1) * t92;
t124 = rSges(6,1) * t85;
t123 = rSges(5,2) * t90;
t122 = rSges(6,2) * t83;
t121 = t81 * rSges(6,3) + t80 * t122;
t120 = t81 * rSges(5,3) + t80 * t123;
t119 = t77 + t78;
t118 = Icges(5,4) * t90;
t117 = Icges(5,4) * t92;
t116 = Icges(6,4) * t83;
t115 = Icges(6,4) * t85;
t101 = Icges(6,1) * t85 - t116;
t54 = Icges(6,2) * t85 + t116;
t55 = Icges(6,1) * t83 + t115;
t104 = t54 * t83 - t55 * t85;
t53 = Icges(6,5) * t83 + Icges(6,6) * t85;
t99 = -Icges(6,2) * t83 + t115;
t114 = (-t104 * t81 + t85 * (Icges(6,6) * t80 + t99 * t81) + t83 * (Icges(6,5) * t80 + t101 * t81) + t80 * t53) * t135 + (t104 * t80 + t85 * (Icges(6,6) * t81 - t80 * t99) + t83 * (Icges(6,5) * t81 - t101 * t80) + t81 * t53) * t134;
t113 = -pkin(3) - t125;
t112 = pkin(4) * t90 + t56;
t58 = -t86 * rSges(3,1) + t84 * rSges(3,2);
t51 = -t81 * rSges(4,1) + t80 * rSges(4,2);
t111 = -t82 - t124;
t110 = -t80 * rSges(6,3) + t81 * t122;
t67 = Icges(5,2) * t92 + t118;
t68 = Icges(5,1) * t90 + t117;
t109 = t85 * t54 + t83 * t55 + t92 * t67 + t90 * t68 + Icges(4,3);
t57 = -t84 * rSges(3,1) - t86 * rSges(3,2);
t50 = -t80 * rSges(4,1) - t81 * rSges(4,2);
t103 = t67 * t90 - t68 * t92;
t102 = Icges(5,1) * t92 - t118;
t100 = -Icges(5,2) * t90 + t117;
t98 = Icges(5,5) * t92 - Icges(5,6) * t90;
t45 = t51 - t129;
t96 = Icges(3,3) + t109;
t66 = Icges(5,5) * t90 + Icges(5,6) * t92;
t95 = t114 + (-t103 * t81 + t92 * (Icges(5,6) * t80 + t100 * t81) + t90 * (Icges(5,5) * t80 + t102 * t81) + t80 * t66) * t135 + (t103 * t80 + t92 * (Icges(5,6) * t81 - t100 * t80) + t90 * (Icges(5,5) * t81 - t102 * t80) + t81 * t66) * t134;
t44 = t50 - t130;
t76 = t81 * pkin(8);
t24 = t113 * t80 + t120 + t76;
t94 = -pkin(9) - pkin(8);
t72 = t80 * t94;
t21 = t111 * t81 + t110 + t72;
t20 = t111 * t80 - t81 * t94 + t121;
t64 = t81 * t123;
t25 = t64 + t113 * t81 + (-rSges(5,3) - pkin(8)) * t80;
t22 = t24 - t130;
t16 = t20 - t130;
t17 = t21 - t129;
t23 = t25 - t129;
t71 = -t93 * rSges(2,1) + t91 * rSges(2,2);
t70 = -t91 * rSges(2,1) - t93 * rSges(2,2);
t49 = t58 - t127;
t48 = t57 - t128;
t43 = t45 - t127;
t42 = t44 - t128;
t37 = Icges(5,3) * t80 + t98 * t81;
t36 = Icges(5,3) * t81 - t80 * t98;
t35 = t112 * t81;
t34 = t112 * t80;
t33 = -t80 * t124 + t121;
t26 = t81 * (t81 * t124 - t110);
t19 = t23 - t127;
t18 = t22 - t128;
t13 = t17 - t127;
t12 = t16 - t128;
t11 = t81 * (t80 * rSges(5,3) + t81 * t125 - t64) - t80 * (-t80 * t125 + t120);
t6 = -t80 * t33 + t26;
t3 = t26 + (t126 * t81 - t72) * t81 + (-t33 + t76 + t126 * t80 + (-pkin(8) + t94) * t81) * t80;
t1 = [Icges(2,3) + m(2) * (t70 ^ 2 + t71 ^ 2) + m(3) * (t48 ^ 2 + t49 ^ 2) + m(4) * (t42 ^ 2 + t43 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + m(6) * (t12 ^ 2 + t13 ^ 2) + t96; m(3) * (t57 * t48 + t58 * t49) + m(4) * (t44 * t42 + t45 * t43) + m(5) * (t22 * t18 + t23 * t19) + m(6) * (t16 * t12 + t17 * t13) + t96; m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t44 ^ 2 + t45 ^ 2) + m(3) * (t57 ^ 2 + t58 ^ 2) + t96; m(4) * (t50 * t42 + t51 * t43) + m(5) * (t24 * t18 + t25 * t19) + m(6) * (t20 * t12 + t21 * t13) + t109; m(6) * (t20 * t16 + t21 * t17) + m(5) * (t24 * t22 + t25 * t23) + m(4) * (t50 * t44 + t51 * t45) + t109; m(6) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t24 ^ 2 + t25 ^ 2) + m(4) * (t50 ^ 2 + t51 ^ 2) + t109; m(6) * (-t35 * t12 + t34 * t13) + (-t18 * t81 + t19 * t80) * t132 + t95; m(6) * (-t35 * t16 + t34 * t17) + (-t22 * t81 + t23 * t80) * t132 + t95; m(6) * (-t35 * t20 + t34 * t21) + (-t24 * t81 + t25 * t80) * t132 + t95; m(5) * (t119 * t69 ^ 2 + t11 ^ 2) + t81 * (t37 * t136 + t78 * t36) + t80 * (t36 * t136 + t77 * t37) + m(6) * (t3 ^ 2 + t34 ^ 2 + t35 ^ 2) + t133; (-t12 * t81 + t13 * t80) * t131 + t114; (-t16 * t81 + t17 * t80) * t131 + t114; (-t20 * t81 + t21 * t80) * t131 + t114; m(6) * (t6 * t3 + (t34 * t80 + t35 * t81) * t56) + t133; m(6) * (t119 * t56 ^ 2 + t6 ^ 2) + t133;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
