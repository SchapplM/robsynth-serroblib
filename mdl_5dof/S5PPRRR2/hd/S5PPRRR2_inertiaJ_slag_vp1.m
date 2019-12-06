% Calculate joint inertia matrix for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:11
% EndTime: 2019-12-05 15:14:14
% DurationCPUTime: 1.08s
% Computational Cost: add. (4624->211), mult. (4710->343), div. (0->0), fcn. (5157->8), ass. (0->115)
t104 = sin(pkin(8));
t100 = t104 ^ 2;
t105 = cos(pkin(8));
t101 = t105 ^ 2;
t120 = t100 + t101;
t102 = pkin(9) + qJ(3);
t97 = cos(t102);
t142 = t97 ^ 2;
t141 = -t97 / 0.2e1;
t140 = t104 / 0.2e1;
t139 = -t105 / 0.2e1;
t107 = cos(qJ(4));
t138 = pkin(4) * t107;
t96 = sin(t102);
t103 = qJ(4) + qJ(5);
t98 = sin(t103);
t99 = cos(t103);
t67 = -Icges(6,3) * t97 + (Icges(6,5) * t99 - Icges(6,6) * t98) * t96;
t137 = t67 * t97;
t106 = sin(qJ(4));
t71 = -Icges(5,3) * t97 + (Icges(5,5) * t107 - Icges(5,6) * t106) * t96;
t136 = t71 * t97;
t109 = pkin(7) * t96 + t138 * t97;
t124 = t104 * t106;
t127 = t105 * t96;
t126 = t105 * t98;
t128 = t104 * t99;
t85 = -t97 * t126 + t128;
t125 = t105 * t99;
t129 = t104 * t98;
t86 = t97 * t125 + t129;
t51 = rSges(6,1) * t86 + rSges(6,2) * t85 + rSges(6,3) * t127;
t134 = pkin(4) * t124 + t109 * t105 + t51;
t130 = t104 * t96;
t83 = -t97 * t129 - t125;
t84 = t97 * t128 - t126;
t50 = rSges(6,1) * t84 + rSges(6,2) * t83 + rSges(6,3) * t130;
t70 = -t97 * rSges(6,3) + (rSges(6,1) * t99 - rSges(6,2) * t98) * t96;
t36 = t70 * t130 + t97 * t50;
t66 = -pkin(7) * t97 + t138 * t96;
t133 = -t66 - t70;
t74 = -t97 * rSges(5,3) + (rSges(5,1) * t107 - rSges(5,2) * t106) * t96;
t93 = t96 * pkin(3) - t97 * pkin(6);
t132 = -t74 - t93;
t131 = t120 * (pkin(3) * t97 + pkin(6) * t96);
t123 = t104 * t107;
t122 = t105 * t106;
t121 = t105 * t107;
t119 = -t93 + t133;
t44 = Icges(6,5) * t84 + Icges(6,6) * t83 + Icges(6,3) * t130;
t46 = Icges(6,4) * t84 + Icges(6,2) * t83 + Icges(6,6) * t130;
t48 = Icges(6,1) * t84 + Icges(6,4) * t83 + Icges(6,5) * t130;
t25 = -t44 * t97 + (-t46 * t98 + t48 * t99) * t96;
t45 = Icges(6,5) * t86 + Icges(6,6) * t85 + Icges(6,3) * t127;
t47 = Icges(6,4) * t86 + Icges(6,2) * t85 + Icges(6,6) * t127;
t49 = Icges(6,1) * t86 + Icges(6,4) * t85 + Icges(6,5) * t127;
t26 = -t45 * t97 + (-t47 * t98 + t49 * t99) * t96;
t19 = t44 * t130 + t46 * t83 + t48 * t84;
t20 = t45 * t130 + t47 * t83 + t49 * t84;
t68 = -Icges(6,6) * t97 + (Icges(6,4) * t99 - Icges(6,2) * t98) * t96;
t69 = -Icges(6,5) * t97 + (Icges(6,1) * t99 - Icges(6,4) * t98) * t96;
t5 = -(t68 * t83 + t69 * t84) * t97 + (t20 * t105 + (t19 - t137) * t104) * t96;
t21 = t44 * t127 + t46 * t85 + t48 * t86;
t22 = t45 * t127 + t47 * t85 + t49 * t86;
t6 = -(t68 * t85 + t69 * t86) * t97 + (t21 * t104 + (t22 - t137) * t105) * t96;
t118 = t6 * t127 + t5 * t130 - t97 * (t142 * t67 + (t26 * t105 + t25 * t104 - (-t68 * t98 + t69 * t99) * t97) * t96);
t12 = t104 * t20 - t105 * t19;
t13 = t104 * t22 - t105 * t21;
t117 = t12 * t130 / 0.2e1 + t5 * t139 + t6 * t140 + t13 * t127 / 0.2e1 + (t104 * t26 - t105 * t25) * t141;
t110 = Icges(4,5) * t97 - Icges(4,6) * t96;
t92 = rSges(4,1) * t96 + rSges(4,2) * t97;
t90 = t97 * t121 + t124;
t89 = -t97 * t122 + t123;
t88 = t97 * t123 - t122;
t87 = -t97 * t124 - t121;
t76 = Icges(4,3) * t104 + t110 * t105;
t75 = -Icges(4,3) * t105 + t110 * t104;
t73 = -Icges(5,5) * t97 + (Icges(5,1) * t107 - Icges(5,4) * t106) * t96;
t72 = -Icges(5,6) * t97 + (Icges(5,4) * t107 - Icges(5,2) * t106) * t96;
t64 = t132 * t105;
t63 = t132 * t104;
t62 = rSges(5,1) * t90 + rSges(5,2) * t89 + rSges(5,3) * t127;
t61 = rSges(5,1) * t88 + rSges(5,2) * t87 + rSges(5,3) * t130;
t60 = Icges(5,1) * t90 + Icges(5,4) * t89 + Icges(5,5) * t127;
t59 = Icges(5,1) * t88 + Icges(5,4) * t87 + Icges(5,5) * t130;
t58 = Icges(5,4) * t90 + Icges(5,2) * t89 + Icges(5,6) * t127;
t57 = Icges(5,4) * t88 + Icges(5,2) * t87 + Icges(5,6) * t130;
t56 = Icges(5,5) * t90 + Icges(5,6) * t89 + Icges(5,3) * t127;
t55 = Icges(5,5) * t88 + Icges(5,6) * t87 + Icges(5,3) * t130;
t53 = -pkin(4) * t122 + t109 * t104;
t52 = t120 * (rSges(4,1) * t97 - rSges(4,2) * t96);
t42 = t50 * t127;
t41 = t119 * t105;
t40 = t119 * t104;
t39 = -t74 * t127 - t62 * t97;
t38 = t74 * t130 + t61 * t97;
t37 = -t70 * t127 - t51 * t97;
t35 = (-t104 * t62 + t105 * t61) * t96;
t34 = -t51 * t130 + t42;
t33 = t104 * t61 + t105 * t62 + t131;
t32 = -t56 * t97 + (-t106 * t58 + t107 * t60) * t96;
t31 = -t55 * t97 + (-t106 * t57 + t107 * t59) * t96;
t30 = t56 * t127 + t58 * t89 + t60 * t90;
t29 = t55 * t127 + t57 * t89 + t59 * t90;
t28 = t56 * t130 + t58 * t87 + t60 * t88;
t27 = t55 * t130 + t57 * t87 + t59 * t88;
t24 = t133 * t127 - t134 * t97;
t23 = t66 * t130 + t53 * t97 + t36;
t18 = t42 + (-t134 * t104 + t105 * t53) * t96;
t17 = t134 * t105 + (t50 + t53) * t104 + t131;
t16 = t104 * t30 - t105 * t29;
t15 = t104 * t28 - t105 * t27;
t8 = -(t72 * t89 + t73 * t90) * t97 + (t29 * t104 + (t30 - t136) * t105) * t96;
t7 = -(t72 * t87 + t73 * t88) * t97 + (t28 * t105 + (t27 - t136) * t104) * t96;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t120; m(4) * t52 + m(5) * t33 + m(6) * t17; m(5) * (t104 * t64 - t105 * t63) + m(6) * (t104 * t41 - t105 * t40); m(6) * (t17 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(5) * (t33 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(4) * (t120 * t92 ^ 2 + t52 ^ 2) + (t100 * t76 + t13 + t16) * t104 + (-t101 * t75 - t12 - t15 + (-t104 * t75 + t105 * t76) * t104) * t105; m(5) * t35 + m(6) * t18; m(5) * (t104 * t38 - t105 * t39) + m(6) * (t104 * t23 - t105 * t24); (t104 * t32 - t105 * t31) * t141 + t8 * t140 + t7 * t139 + (t15 * t140 + t105 * t16 / 0.2e1) * t96 + m(6) * (t17 * t18 + t23 * t41 + t24 * t40) + m(5) * (t33 * t35 + t38 * t64 + t39 * t63) + t117; m(6) * (t18 ^ 2 + t23 ^ 2 + t24 ^ 2) - t97 * (t142 * t71 + (t32 * t105 + t31 * t104 - (-t106 * t72 + t107 * t73) * t97) * t96) + t7 * t130 + m(5) * (t35 ^ 2 + t38 ^ 2 + t39 ^ 2) + t8 * t127 + t118; m(6) * t34; m(6) * (t104 * t36 - t105 * t37); m(6) * (t17 * t34 + t36 * t41 + t37 * t40) + t117; m(6) * (t18 * t34 + t23 * t36 + t24 * t37) + t118; m(6) * (t34 ^ 2 + t36 ^ 2 + t37 ^ 2) + t118;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
