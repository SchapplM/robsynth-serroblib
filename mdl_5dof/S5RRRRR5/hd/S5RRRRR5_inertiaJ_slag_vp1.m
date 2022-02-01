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
% m [6x1]
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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 12:01:39
% EndTime: 2022-01-20 12:01:40
% DurationCPUTime: 0.69s
% Computational Cost: add. (3060->180), mult. (1852->254), div. (0->0), fcn. (1594->10), ass. (0->98)
t95 = qJ(1) + qJ(2);
t92 = qJ(3) + t95;
t85 = sin(t92);
t86 = cos(t92);
t139 = t85 * t86;
t81 = t85 ^ 2;
t82 = t86 ^ 2;
t138 = t85 / 0.2e1;
t137 = -t86 / 0.2e1;
t96 = sin(qJ(4));
t98 = cos(qJ(4));
t70 = t96 * rSges(5,1) + t98 * rSges(5,2);
t136 = m(5) * t70;
t94 = qJ(4) + qJ(5);
t88 = sin(t94);
t90 = cos(t94);
t56 = t88 * rSges(6,1) + t90 * rSges(6,2);
t135 = m(6) * t56;
t89 = sin(t95);
t134 = pkin(2) * t89;
t97 = sin(qJ(1));
t133 = t97 * pkin(1);
t132 = rSges(5,1) * t98;
t131 = rSges(6,1) * t90;
t130 = rSges(5,2) * t96;
t129 = rSges(6,2) * t88;
t103 = t85 * rSges(6,3) + (-t129 + t131) * t86;
t128 = t86 * rSges(6,3) + t85 * t129;
t6 = t85 * (t85 * t131 - t128) + t86 * t103;
t127 = t86 * rSges(5,3) + t85 * t130;
t126 = -t86 * pkin(3) - t85 * pkin(8);
t125 = t81 + t82;
t124 = Icges(5,4) * t96;
t123 = Icges(5,4) * t98;
t122 = Icges(6,4) * t88;
t121 = Icges(6,4) * t90;
t107 = -Icges(6,2) * t88 + t121;
t109 = Icges(6,1) * t90 - t122;
t54 = Icges(6,2) * t90 + t122;
t55 = Icges(6,1) * t88 + t121;
t112 = -t54 * t88 + t55 * t90;
t53 = Icges(6,5) * t88 + Icges(6,6) * t90;
t120 = (t112 * t86 + t90 * (Icges(6,6) * t85 + t107 * t86) + t88 * (Icges(6,5) * t85 + t109 * t86) + t85 * t53) * t138 + (t112 * t85 + t90 * (-Icges(6,6) * t86 + t107 * t85) + t88 * (-Icges(6,5) * t86 + t109 * t85) - t86 * t53) * t137;
t105 = Icges(6,5) * t90 - Icges(6,6) * t88;
t28 = -Icges(6,3) * t86 + t105 * t85;
t29 = Icges(6,3) * t85 + t105 * t86;
t119 = -t86 * (-t29 * t139 + t82 * t28) + t85 * (-t28 * t139 + t81 * t29);
t118 = -pkin(4) * t96 - t56;
t91 = cos(t95);
t58 = t91 * rSges(3,1) - t89 * rSges(3,2);
t51 = t86 * rSges(4,1) - t85 * rSges(4,2);
t84 = pkin(2) * t91;
t45 = t51 + t84;
t68 = Icges(5,2) * t98 + t124;
t69 = Icges(5,1) * t96 + t123;
t117 = t90 * t54 + t88 * t55 + t98 * t68 + t96 * t69 + Icges(4,3);
t57 = -t89 * rSges(3,1) - t91 * rSges(3,2);
t50 = -t85 * rSges(4,1) - t86 * rSges(4,2);
t111 = -t68 * t96 + t69 * t98;
t110 = Icges(5,1) * t98 - t124;
t108 = -Icges(5,2) * t96 + t123;
t106 = Icges(5,5) * t98 - Icges(5,6) * t96;
t104 = t85 * rSges(5,3) + (-t130 + t132) * t86;
t102 = Icges(3,3) + t117;
t67 = Icges(5,5) * t96 + Icges(5,6) * t98;
t101 = t120 + (t111 * t86 + t98 * (Icges(5,6) * t85 + t108 * t86) + t96 * (Icges(5,5) * t85 + t110 * t86) + t85 * t67) * t138 + (t111 * t85 + t98 * (-Icges(5,6) * t86 + t108 * t85) + t96 * (-Icges(5,5) * t86 + t110 * t85) - t86 * t67) * t137;
t25 = t104 - t126;
t44 = t50 - t134;
t23 = t25 + t84;
t79 = t86 * pkin(8);
t24 = t79 + (-pkin(3) - t132) * t85 + t127;
t100 = -pkin(9) - pkin(8);
t87 = t98 * pkin(4) + pkin(3);
t61 = t86 * t87;
t21 = -t85 * t100 + t103 + t61;
t17 = t21 + t84;
t20 = -t86 * t100 + (-t87 - t131) * t85 + t128;
t22 = t24 - t134;
t16 = t20 - t134;
t99 = cos(qJ(1));
t93 = t99 * pkin(1);
t72 = t99 * rSges(2,1) - t97 * rSges(2,2);
t71 = -t97 * rSges(2,1) - t99 * rSges(2,2);
t49 = t58 + t93;
t48 = t57 - t133;
t43 = t45 + t93;
t42 = t44 - t133;
t37 = Icges(5,3) * t85 + t106 * t86;
t36 = -Icges(5,3) * t86 + t106 * t85;
t35 = t118 * t86;
t34 = t118 * t85;
t19 = t23 + t93;
t18 = t22 - t133;
t13 = t17 + t93;
t12 = t16 - t133;
t11 = t85 * (t85 * t132 - t127) + t86 * t104;
t3 = t86 * (t61 + t126) + (t79 + (-pkin(3) + t87) * t85) * t85 + t6;
t1 = [Icges(2,3) + m(6) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + m(4) * (t42 ^ 2 + t43 ^ 2) + m(3) * (t48 ^ 2 + t49 ^ 2) + m(2) * (t71 ^ 2 + t72 ^ 2) + t102; m(6) * (t16 * t12 + t17 * t13) + m(5) * (t22 * t18 + t23 * t19) + m(4) * (t44 * t42 + t45 * t43) + m(3) * (t57 * t48 + t58 * t49) + t102; m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t44 ^ 2 + t45 ^ 2) + m(3) * (t57 ^ 2 + t58 ^ 2) + t102; m(6) * (t20 * t12 + t21 * t13) + m(5) * (t24 * t18 + t25 * t19) + m(4) * (t50 * t42 + t51 * t43) + t117; m(6) * (t20 * t16 + t21 * t17) + m(5) * (t24 * t22 + t25 * t23) + m(4) * (t50 * t44 + t51 * t45) + t117; m(6) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t24 ^ 2 + t25 ^ 2) + m(4) * (t50 ^ 2 + t51 ^ 2) + t117; m(6) * (t35 * t12 + t34 * t13) + (-t18 * t86 - t19 * t85) * t136 + t101; m(6) * (t35 * t16 + t34 * t17) + (-t22 * t86 - t23 * t85) * t136 + t101; m(6) * (t35 * t20 + t34 * t21) + (-t24 * t86 - t25 * t85) * t136 + t101; m(5) * (t125 * t70 ^ 2 + t11 ^ 2) + t85 * (-t36 * t139 + t81 * t37) - t86 * (-t37 * t139 + t82 * t36) + m(6) * (t3 ^ 2 + t34 ^ 2 + t35 ^ 2) + t119; (-t12 * t86 - t13 * t85) * t135 + t120; (-t16 * t86 - t17 * t85) * t135 + t120; (-t20 * t86 - t21 * t85) * t135 + t120; m(6) * (t6 * t3 + (-t34 * t85 - t35 * t86) * t56) + t119; m(6) * (t125 * t56 ^ 2 + t6 ^ 2) + t119;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
