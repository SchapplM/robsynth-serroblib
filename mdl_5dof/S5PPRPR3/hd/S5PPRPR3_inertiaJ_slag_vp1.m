% Calculate joint inertia matrix for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:47
% EndTime: 2019-12-05 15:04:52
% DurationCPUTime: 1.51s
% Computational Cost: add. (3321->242), mult. (4944->381), div. (0->0), fcn. (5783->10), ass. (0->112)
t132 = Icges(5,3) + Icges(4,3);
t97 = sin(pkin(8));
t98 = sin(pkin(7));
t124 = t98 * t97;
t100 = cos(pkin(7));
t99 = cos(pkin(8));
t123 = t99 * t98;
t96 = qJ(3) + pkin(9);
t91 = sin(t96);
t92 = cos(t96);
t72 = t92 * t100 + t91 * t123;
t73 = -t91 * t100 + t92 * t123;
t105 = cos(qJ(3));
t109 = t105 * t100;
t103 = sin(qJ(3));
t112 = t98 * t103;
t84 = -t99 * t112 - t109;
t110 = t100 * t103;
t113 = t105 * t98;
t85 = t99 * t113 - t110;
t131 = Icges(4,5) * t85 + Icges(5,5) * t73 + Icges(4,6) * t84 - Icges(5,6) * t72 + t132 * t124;
t116 = t100 * t97;
t111 = t99 * t100;
t74 = t91 * t111 - t98 * t92;
t75 = t92 * t111 + t91 * t98;
t86 = -t99 * t110 + t113;
t87 = t99 * t109 + t112;
t130 = Icges(4,5) * t87 + Icges(5,5) * t75 + Icges(4,6) * t86 - Icges(5,6) * t74 + t132 * t116;
t129 = t132 * t99 + (-Icges(4,5) * t105 - Icges(5,5) * t92 + Icges(4,6) * t103 + Icges(5,6) * t91) * t97;
t128 = t99 ^ 2;
t127 = -m(5) - m(6);
t126 = pkin(3) * t105;
t125 = t91 * t97;
t102 = sin(qJ(5));
t104 = cos(qJ(5));
t114 = t104 * t97;
t62 = -t102 * t73 + t98 * t114;
t115 = t102 * t97;
t63 = t104 * t73 + t98 * t115;
t32 = rSges(6,1) * t63 + rSges(6,2) * t62 + rSges(6,3) * t72;
t121 = pkin(4) * t73 + pkin(6) * t72 + t32;
t106 = qJ(4) * t97 + t126 * t99;
t59 = pkin(3) * t112 + t106 * t100;
t120 = -rSges(5,1) * t75 + rSges(5,2) * t74 - rSges(5,3) * t116 - t59;
t81 = -t104 * t99 - t92 * t115;
t82 = -t102 * t99 + t92 * t114;
t48 = rSges(6,1) * t82 + rSges(6,2) * t81 + rSges(6,3) * t125;
t119 = t48 + (pkin(4) * t92 + pkin(6) * t91) * t97;
t58 = -pkin(3) * t110 + t106 * t98;
t67 = -qJ(4) * t99 + t126 * t97;
t118 = t67 * t124 + t99 * t58;
t117 = t100 ^ 2 + t98 ^ 2;
t108 = m(5) / 0.2e1 + m(6) / 0.2e1;
t64 = t100 * t114 - t102 * t75;
t65 = t100 * t115 + t104 * t75;
t33 = rSges(6,1) * t65 + rSges(6,2) * t64 + rSges(6,3) * t74;
t107 = -pkin(4) * t75 - pkin(6) * t74 - t33 - t59;
t79 = -t99 * rSges(4,3) + (rSges(4,1) * t105 - rSges(4,2) * t103) * t97;
t78 = -Icges(4,5) * t99 + (Icges(4,1) * t105 - Icges(4,4) * t103) * t97;
t77 = -Icges(4,6) * t99 + (Icges(4,4) * t105 - Icges(4,2) * t103) * t97;
t71 = -rSges(5,3) * t99 + (rSges(5,1) * t92 - rSges(5,2) * t91) * t97;
t70 = -Icges(5,5) * t99 + (Icges(5,1) * t92 - Icges(5,4) * t91) * t97;
t69 = -Icges(5,6) * t99 + (Icges(5,4) * t92 - Icges(5,2) * t91) * t97;
t61 = rSges(4,1) * t87 + rSges(4,2) * t86 + rSges(4,3) * t116;
t60 = rSges(4,1) * t85 + rSges(4,2) * t84 + rSges(4,3) * t124;
t57 = Icges(4,1) * t87 + Icges(4,4) * t86 + Icges(4,5) * t116;
t56 = Icges(4,1) * t85 + Icges(4,4) * t84 + Icges(4,5) * t124;
t55 = Icges(4,4) * t87 + Icges(4,2) * t86 + Icges(4,6) * t116;
t54 = Icges(4,4) * t85 + Icges(4,2) * t84 + Icges(4,6) * t124;
t47 = Icges(6,1) * t82 + Icges(6,4) * t81 + Icges(6,5) * t125;
t46 = Icges(6,4) * t82 + Icges(6,2) * t81 + Icges(6,6) * t125;
t45 = Icges(6,5) * t82 + Icges(6,6) * t81 + Icges(6,3) * t125;
t44 = t58 * t116;
t42 = rSges(5,1) * t73 - rSges(5,2) * t72 + rSges(5,3) * t124;
t41 = Icges(5,1) * t75 - Icges(5,4) * t74 + Icges(5,5) * t116;
t40 = Icges(5,1) * t73 - Icges(5,4) * t72 + Icges(5,5) * t124;
t39 = Icges(5,4) * t75 - Icges(5,2) * t74 + Icges(5,6) * t116;
t38 = Icges(5,4) * t73 - Icges(5,2) * t72 + Icges(5,6) * t124;
t35 = -t79 * t116 - t61 * t99;
t34 = t79 * t124 + t60 * t99;
t31 = Icges(6,1) * t65 + Icges(6,4) * t64 + Icges(6,5) * t74;
t30 = Icges(6,1) * t63 + Icges(6,4) * t62 + Icges(6,5) * t72;
t29 = Icges(6,4) * t65 + Icges(6,2) * t64 + Icges(6,6) * t74;
t28 = Icges(6,4) * t63 + Icges(6,2) * t62 + Icges(6,6) * t72;
t27 = Icges(6,5) * t65 + Icges(6,6) * t64 + Icges(6,3) * t74;
t26 = Icges(6,5) * t63 + Icges(6,6) * t62 + Icges(6,3) * t72;
t25 = (t100 * t60 - t61 * t98) * t97;
t24 = t120 * t99 + (-t67 - t71) * t116;
t23 = t71 * t124 + t42 * t99 + t118;
t22 = t33 * t125 - t48 * t74;
t21 = -t32 * t125 + t48 * t72;
t20 = t45 * t125 + t46 * t81 + t47 * t82;
t19 = t44 + (t100 * t42 + t120 * t98) * t97;
t18 = t32 * t74 - t33 * t72;
t17 = t45 * t74 + t46 * t64 + t47 * t65;
t16 = t45 * t72 + t46 * t62 + t47 * t63;
t15 = t107 * t99 + (-t67 - t119) * t116;
t14 = t119 * t124 + t121 * t99 + t118;
t13 = t27 * t125 + t29 * t81 + t31 * t82;
t12 = t26 * t125 + t28 * t81 + t30 * t82;
t11 = t27 * t74 + t29 * t64 + t31 * t65;
t10 = t26 * t74 + t28 * t64 + t30 * t65;
t9 = t27 * t72 + t29 * t62 + t31 * t63;
t8 = t26 * t72 + t28 * t62 + t30 * t63;
t7 = t44 + (t121 * t100 + t107 * t98) * t97;
t6 = -t20 * t99 + (t100 * t13 + t12 * t98) * t97;
t5 = t12 * t72 + t20 * t125 + t13 * t74;
t4 = -t17 * t99 + (t10 * t98 + t100 * t11) * t97;
t3 = -t16 * t99 + (t100 * t9 + t8 * t98) * t97;
t2 = t10 * t72 + t11 * t74 + t17 * t125;
t1 = t16 * t125 + t72 * t8 + t74 * t9;
t36 = [m(2) + m(3) + m(4) - t127; 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t108) * t117; m(4) * t25 + m(5) * t19 + m(6) * t7; m(4) * (-t100 * t35 + t34 * t98) + m(5) * (-t100 * t24 + t23 * t98) + m(6) * (-t100 * t15 + t14 * t98); m(6) * (t14 ^ 2 + t15 ^ 2 + t7 ^ 2) + m(5) * (t19 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(4) * (t25 ^ 2 + t34 ^ 2 + t35 ^ 2) + (-t6 + t129 * t128 + ((-t103 * t77 + t105 * t78 - t69 * t91 + t70 * t92) * t99 + (t131 * t99 + (t103 * t54 - t105 * t56 + t38 * t91 - t40 * t92) * t97) * t98 + (t130 * t99 + (t103 * t55 - t105 * t57 + t39 * t91 - t41 * t92) * t97) * t100) * t97) * t99 + (t3 + (t131 * t124 - t38 * t72 + t40 * t73 + t54 * t84 + t56 * t85) * t124 + (t129 * t124 + t69 * t72 - t70 * t73 - t77 * t84 - t78 * t85) * t99) * t124 + (t4 + (t130 * t116 - t39 * t74 + t41 * t75 + t55 * t86 + t57 * t87) * t116 + (t129 * t116 + t69 * t74 - t70 * t75 - t77 * t86 - t78 * t87) * t99 + (t131 * t116 + t130 * t124 - t38 * t74 - t39 * t72 + t40 * t75 + t41 * t73 + t54 * t86 + t55 * t84 + t56 * t87 + t57 * t85) * t124) * t116; t127 * t99; 0; m(6) * (-t7 * t99 + (t100 * t14 + t15 * t98) * t97) + m(5) * (-t19 * t99 + (t100 * t23 + t24 * t98) * t97); 0.2e1 * t108 * (t117 * t97 ^ 2 + t128); m(6) * t18; m(6) * (-t100 * t22 + t21 * t98); -t99 * t5 / 0.2e1 + t72 * t3 / 0.2e1 + t74 * t4 / 0.2e1 + m(6) * (t14 * t21 + t15 * t22 + t18 * t7) + (t100 * t2 / 0.2e1 + t98 * t1 / 0.2e1 + t91 * t6 / 0.2e1) * t97; m(6) * (-t18 * t99 + (t100 * t21 + t22 * t98) * t97); m(6) * (t18 ^ 2 + t21 ^ 2 + t22 ^ 2) + t74 * t2 + t72 * t1 + t5 * t125;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t36(1), t36(2), t36(4), t36(7), t36(11); t36(2), t36(3), t36(5), t36(8), t36(12); t36(4), t36(5), t36(6), t36(9), t36(13); t36(7), t36(8), t36(9), t36(10), t36(14); t36(11), t36(12), t36(13), t36(14), t36(15);];
Mq = res;
