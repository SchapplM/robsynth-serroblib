% Calculate joint inertia matrix for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:34
% EndTime: 2019-12-05 15:10:38
% DurationCPUTime: 1.30s
% Computational Cost: add. (3060->266), mult. (7884->413), div. (0->0), fcn. (9596->8), ass. (0->115)
t126 = rSges(6,1) + pkin(4);
t125 = rSges(6,3) + qJ(5);
t124 = cos(qJ(4));
t106 = sin(pkin(7));
t109 = sin(qJ(4));
t105 = sin(pkin(8));
t112 = t105 * t124;
t108 = cos(pkin(7));
t110 = sin(qJ(3));
t107 = cos(pkin(8));
t111 = cos(qJ(3));
t114 = t111 * t107;
t95 = t106 * t114 - t108 * t110;
t85 = -t106 * t112 + t95 * t109;
t117 = t106 * t105;
t86 = t109 * t117 + t95 * t124;
t115 = t110 * t107;
t94 = t106 * t115 + t108 * t111;
t123 = t94 * rSges(6,2) + t125 * t85 + t126 * t86;
t97 = t106 * t110 + t108 * t114;
t87 = -t108 * t112 + t97 * t109;
t116 = t108 * t105;
t88 = t109 * t116 + t97 * t124;
t96 = -t106 * t111 + t108 * t115;
t122 = t96 * rSges(6,2) + t125 * t87 + t126 * t88;
t59 = t88 * rSges(5,1) - t87 * rSges(5,2) + t96 * rSges(5,3);
t83 = t97 * pkin(3) + t96 * pkin(6);
t121 = -t59 - t83;
t118 = t105 * t110;
t98 = t105 * t111 * t109 + t107 * t124;
t99 = -t107 * t109 + t111 * t112;
t120 = rSges(6,2) * t118 + t125 * t98 + t126 * t99;
t100 = (pkin(3) * t111 + pkin(6) * t110) * t105;
t82 = t95 * pkin(3) + t94 * pkin(6);
t119 = t100 * t117 + t107 * t82;
t113 = -t83 - t122;
t93 = -t107 * rSges(4,3) + (rSges(4,1) * t111 - rSges(4,2) * t110) * t105;
t92 = -Icges(4,5) * t107 + (Icges(4,1) * t111 - Icges(4,4) * t110) * t105;
t91 = -Icges(4,6) * t107 + (Icges(4,4) * t111 - Icges(4,2) * t110) * t105;
t90 = -Icges(4,3) * t107 + (Icges(4,5) * t111 - Icges(4,6) * t110) * t105;
t80 = t82 * t116;
t79 = t99 * rSges(5,1) - t98 * rSges(5,2) + rSges(5,3) * t118;
t77 = Icges(5,1) * t99 - Icges(5,4) * t98 + Icges(5,5) * t118;
t76 = Icges(6,1) * t99 + Icges(6,4) * t118 + Icges(6,5) * t98;
t75 = Icges(5,4) * t99 - Icges(5,2) * t98 + Icges(5,6) * t118;
t74 = Icges(6,4) * t99 + Icges(6,2) * t118 + Icges(6,6) * t98;
t73 = Icges(5,5) * t99 - Icges(5,6) * t98 + Icges(5,3) * t118;
t72 = Icges(6,5) * t99 + Icges(6,6) * t118 + Icges(6,3) * t98;
t71 = t97 * rSges(4,1) - t96 * rSges(4,2) + rSges(4,3) * t116;
t70 = t95 * rSges(4,1) - t94 * rSges(4,2) + rSges(4,3) * t117;
t69 = Icges(4,1) * t97 - Icges(4,4) * t96 + Icges(4,5) * t116;
t68 = Icges(4,1) * t95 - Icges(4,4) * t94 + Icges(4,5) * t117;
t67 = Icges(4,4) * t97 - Icges(4,2) * t96 + Icges(4,6) * t116;
t66 = Icges(4,4) * t95 - Icges(4,2) * t94 + Icges(4,6) * t117;
t65 = Icges(4,5) * t97 - Icges(4,6) * t96 + Icges(4,3) * t116;
t64 = Icges(4,5) * t95 - Icges(4,6) * t94 + Icges(4,3) * t117;
t61 = -t107 * t71 - t93 * t116;
t60 = t107 * t70 + t93 * t117;
t57 = t86 * rSges(5,1) - t85 * rSges(5,2) + t94 * rSges(5,3);
t55 = Icges(5,1) * t88 - Icges(5,4) * t87 + Icges(5,5) * t96;
t54 = Icges(5,1) * t86 - Icges(5,4) * t85 + Icges(5,5) * t94;
t53 = Icges(6,1) * t88 + Icges(6,4) * t96 + Icges(6,5) * t87;
t52 = Icges(6,1) * t86 + Icges(6,4) * t94 + Icges(6,5) * t85;
t51 = Icges(5,4) * t88 - Icges(5,2) * t87 + Icges(5,6) * t96;
t50 = Icges(5,4) * t86 - Icges(5,2) * t85 + Icges(5,6) * t94;
t49 = Icges(6,4) * t88 + Icges(6,2) * t96 + Icges(6,6) * t87;
t48 = Icges(6,4) * t86 + Icges(6,2) * t94 + Icges(6,6) * t85;
t47 = Icges(5,5) * t88 - Icges(5,6) * t87 + Icges(5,3) * t96;
t46 = Icges(5,5) * t86 - Icges(5,6) * t85 + Icges(5,3) * t94;
t45 = Icges(6,5) * t88 + Icges(6,6) * t96 + Icges(6,3) * t87;
t44 = Icges(6,5) * t86 + Icges(6,6) * t94 + Icges(6,3) * t85;
t43 = (-t106 * t71 + t108 * t70) * t105;
t42 = t59 * t118 - t96 * t79;
t41 = -t57 * t118 + t94 * t79;
t40 = t73 * t118 - t98 * t75 + t99 * t77;
t39 = t74 * t118 + t98 * t72 + t99 * t76;
t38 = t96 * t57 - t94 * t59;
t37 = t121 * t107 + (-t100 - t79) * t116;
t36 = t107 * t57 + t79 * t117 + t119;
t35 = t96 * t73 - t87 * t75 + t88 * t77;
t34 = t87 * t72 + t96 * t74 + t88 * t76;
t33 = t94 * t73 - t85 * t75 + t86 * t77;
t32 = t85 * t72 + t94 * t74 + t86 * t76;
t31 = t80 + (t121 * t106 + t108 * t57) * t105;
t30 = t122 * t118 - t120 * t96;
t29 = -t123 * t118 + t120 * t94;
t28 = t47 * t118 - t98 * t51 + t99 * t55;
t27 = t46 * t118 - t98 * t50 + t99 * t54;
t26 = t49 * t118 + t98 * t45 + t99 * t53;
t25 = t48 * t118 + t98 * t44 + t99 * t52;
t24 = t113 * t107 + (-t100 - t120) * t116;
t23 = t123 * t107 + t120 * t117 + t119;
t22 = t96 * t47 - t87 * t51 + t88 * t55;
t21 = t96 * t46 - t87 * t50 + t88 * t54;
t20 = t87 * t45 + t96 * t49 + t88 * t53;
t19 = t87 * t44 + t96 * t48 + t88 * t52;
t18 = t94 * t47 - t85 * t51 + t86 * t55;
t17 = t94 * t46 - t85 * t50 + t86 * t54;
t16 = t85 * t45 + t94 * t49 + t86 * t53;
t15 = t85 * t44 + t94 * t48 + t86 * t52;
t14 = -t122 * t94 + t123 * t96;
t13 = t80 + (t113 * t106 + t123 * t108) * t105;
t12 = -t40 * t107 + (t106 * t27 + t108 * t28) * t105;
t11 = -t39 * t107 + (t106 * t25 + t108 * t26) * t105;
t10 = t40 * t118 + t27 * t94 + t28 * t96;
t9 = t39 * t118 + t25 * t94 + t26 * t96;
t8 = -t35 * t107 + (t106 * t21 + t108 * t22) * t105;
t7 = -t34 * t107 + (t106 * t19 + t108 * t20) * t105;
t6 = -t33 * t107 + (t106 * t17 + t108 * t18) * t105;
t5 = -t32 * t107 + (t106 * t15 + t108 * t16) * t105;
t4 = t35 * t118 + t21 * t94 + t22 * t96;
t3 = t34 * t118 + t19 * t94 + t20 * t96;
t2 = t33 * t118 + t17 * t94 + t18 * t96;
t1 = t32 * t118 + t15 * t94 + t16 * t96;
t56 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t106 ^ 2 + t108 ^ 2); m(4) * t43 + m(5) * t31 + m(6) * t13; m(4) * (t60 * t106 - t61 * t108) + m(5) * (t36 * t106 - t37 * t108) + m(6) * (t23 * t106 - t24 * t108); m(6) * (t13 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(5) * (t31 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(4) * (t43 ^ 2 + t60 ^ 2 + t61 ^ 2) + (t5 + t6 + (t64 * t117 - t94 * t66 + t95 * t68) * t117) * t117 + (t8 + t7 + (t65 * t116 - t96 * t67 + t97 * t69) * t116 + (t64 * t116 + t65 * t117 - t96 * t66 - t94 * t67 + t97 * t68 + t95 * t69) * t117) * t116 + ((-t90 * t116 + t96 * t91 - t97 * t92) * t116 + (-t90 * t117 + t94 * t91 - t95 * t92) * t117 - t11 - t12 - ((-t110 * t67 + t111 * t69) * t108 + (-t110 * t66 + t111 * t68) * t106) * t105 ^ 2 + (-(-t106 * t64 - t65 * t108 + t110 * t91 - t111 * t92) * t105 - t107 * t90) * t107) * t107; m(5) * t38 + m(6) * t14; m(5) * (t41 * t106 - t42 * t108) + m(6) * (t29 * t106 - t30 * t108); (t8 / 0.2e1 + t7 / 0.2e1) * t96 + (t6 / 0.2e1 + t5 / 0.2e1) * t94 + (-t9 / 0.2e1 - t10 / 0.2e1) * t107 + m(6) * (t14 * t13 + t29 * t23 + t30 * t24) + m(5) * (t38 * t31 + t41 * t36 + t42 * t37) + ((t12 / 0.2e1 + t11 / 0.2e1) * t110 + (t4 / 0.2e1 + t3 / 0.2e1) * t108 + (t1 / 0.2e1 + t2 / 0.2e1) * t106) * t105; (t3 + t4) * t96 + (t1 + t2) * t94 + (t10 + t9) * t118 + m(6) * (t14 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(5) * (t38 ^ 2 + t41 ^ 2 + t42 ^ 2); m(6) * t98; m(6) * (t87 * t106 - t85 * t108); m(6) * (t98 * t13 + t87 * t23 + t85 * t24); m(6) * (t98 * t14 + t87 * t29 + t85 * t30); m(6) * (t85 ^ 2 + t87 ^ 2 + t98 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t56(1), t56(2), t56(4), t56(7), t56(11); t56(2), t56(3), t56(5), t56(8), t56(12); t56(4), t56(5), t56(6), t56(9), t56(13); t56(7), t56(8), t56(9), t56(10), t56(14); t56(11), t56(12), t56(13), t56(14), t56(15);];
Mq = res;
