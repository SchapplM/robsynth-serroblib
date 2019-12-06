% Calculate time derivative of joint inertia matrix for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPPRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPPRR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:51
% EndTime: 2019-12-05 14:57:57
% DurationCPUTime: 1.91s
% Computational Cost: add. (8349->313), mult. (13055->509), div. (0->0), fcn. (13891->8), ass. (0->134)
t118 = cos(pkin(8));
t115 = pkin(9) + qJ(4);
t113 = sin(t115);
t114 = cos(t115);
t116 = sin(pkin(8));
t126 = t116 * qJD(4);
t99 = (-Icges(5,5) * t113 - Icges(5,6) * t114) * t126;
t144 = t118 * t99;
t119 = cos(pkin(7));
t117 = sin(pkin(7));
t130 = t118 * t117;
t105 = -t113 * t119 + t114 * t130;
t143 = 2 * m(6);
t120 = sin(qJ(5));
t121 = cos(qJ(5));
t131 = t117 * t116;
t87 = t105 * t121 + t120 * t131;
t104 = t113 * t130 + t114 * t119;
t95 = t104 * qJD(4);
t60 = -t87 * qJD(5) + t95 * t120;
t86 = -t105 * t120 + t121 * t131;
t61 = t86 * qJD(5) - t121 * t95;
t96 = t105 * qJD(4);
t43 = rSges(6,1) * t61 + rSges(6,2) * t60 + rSges(6,3) * t96;
t142 = -pkin(4) * t95 + pkin(6) * t96 + t43;
t129 = t118 * t119;
t107 = t113 * t117 + t114 * t129;
t128 = t119 * t116;
t89 = t107 * t121 + t120 * t128;
t106 = t113 * t129 - t117 * t114;
t97 = t106 * qJD(4);
t62 = -t89 * qJD(5) + t97 * t120;
t88 = -t107 * t120 + t121 * t128;
t63 = t88 * qJD(5) - t121 * t97;
t98 = t107 * qJD(4);
t44 = rSges(6,1) * t63 + rSges(6,2) * t62 + rSges(6,3) * t98;
t141 = pkin(4) * t97 - pkin(6) * t98 - t44;
t52 = rSges(6,1) * t87 + rSges(6,2) * t86 + rSges(6,3) * t104;
t140 = pkin(4) * t105 + pkin(6) * t104 + t52;
t53 = rSges(6,1) * t89 + rSges(6,2) * t88 + rSges(6,3) * t106;
t139 = -pkin(4) * t107 - pkin(6) * t106 - t53;
t123 = t114 * t126;
t132 = t114 * t116;
t110 = -t118 * t120 + t121 * t132;
t124 = t113 * t126;
t90 = -t110 * qJD(5) + t120 * t124;
t109 = -t118 * t121 - t120 * t132;
t91 = t109 * qJD(5) - t121 * t124;
t59 = rSges(6,1) * t91 + rSges(6,2) * t90 + rSges(6,3) * t123;
t138 = (-pkin(4) * t113 + pkin(6) * t114) * t126 + t59;
t134 = t113 * t116;
t83 = rSges(6,1) * t110 + rSges(6,2) * t109 + rSges(6,3) * t134;
t137 = (pkin(4) * t114 + pkin(6) * t113) * t116 + t83;
t136 = Icges(5,4) * t113;
t135 = Icges(5,4) * t114;
t127 = qJD(4) * t114;
t102 = (-rSges(5,1) * t113 - rSges(5,2) * t114) * t126;
t76 = -rSges(5,1) * t95 - rSges(5,2) * t96;
t54 = t102 * t131 + t118 * t76;
t77 = -rSges(5,1) * t97 - rSges(5,2) * t98;
t55 = -t102 * t128 - t118 * t77;
t122 = t117 * t54 - t119 * t55;
t101 = (-Icges(5,1) * t113 - t135) * t126;
t100 = (-Icges(5,2) * t114 - t136) * t126;
t93 = -Icges(5,5) * t118 + (Icges(5,1) * t114 - t136) * t116;
t92 = -Icges(5,6) * t118 + (-Icges(5,2) * t113 + t135) * t116;
t82 = Icges(6,1) * t110 + Icges(6,4) * t109 + Icges(6,5) * t134;
t81 = Icges(6,4) * t110 + Icges(6,2) * t109 + Icges(6,6) * t134;
t80 = Icges(6,5) * t110 + Icges(6,6) * t109 + Icges(6,3) * t134;
t75 = -Icges(5,1) * t97 - Icges(5,4) * t98;
t74 = -Icges(5,1) * t95 - Icges(5,4) * t96;
t73 = -Icges(5,4) * t97 - Icges(5,2) * t98;
t72 = -Icges(5,4) * t95 - Icges(5,2) * t96;
t71 = -Icges(5,5) * t97 - Icges(5,6) * t98;
t70 = -Icges(5,5) * t95 - Icges(5,6) * t96;
t69 = rSges(5,1) * t107 - rSges(5,2) * t106 + rSges(5,3) * t128;
t68 = rSges(5,1) * t105 - rSges(5,2) * t104 + rSges(5,3) * t131;
t67 = Icges(5,1) * t107 - Icges(5,4) * t106 + Icges(5,5) * t128;
t66 = Icges(5,1) * t105 - Icges(5,4) * t104 + Icges(5,5) * t131;
t65 = Icges(5,4) * t107 - Icges(5,2) * t106 + Icges(5,6) * t128;
t64 = Icges(5,4) * t105 - Icges(5,2) * t104 + Icges(5,6) * t131;
t58 = Icges(6,1) * t91 + Icges(6,4) * t90 + Icges(6,5) * t123;
t57 = Icges(6,4) * t91 + Icges(6,2) * t90 + Icges(6,6) * t123;
t56 = Icges(6,5) * t91 + Icges(6,6) * t90 + Icges(6,3) * t123;
t51 = Icges(6,1) * t89 + Icges(6,4) * t88 + Icges(6,5) * t106;
t50 = Icges(6,1) * t87 + Icges(6,4) * t86 + Icges(6,5) * t104;
t49 = Icges(6,4) * t89 + Icges(6,2) * t88 + Icges(6,6) * t106;
t48 = Icges(6,4) * t87 + Icges(6,2) * t86 + Icges(6,6) * t104;
t47 = Icges(6,5) * t89 + Icges(6,6) * t88 + Icges(6,3) * t106;
t46 = Icges(6,5) * t87 + Icges(6,6) * t86 + Icges(6,3) * t104;
t45 = (-t117 * t77 + t119 * t76) * t116;
t42 = Icges(6,1) * t63 + Icges(6,4) * t62 + Icges(6,5) * t98;
t41 = Icges(6,1) * t61 + Icges(6,4) * t60 + Icges(6,5) * t96;
t40 = Icges(6,4) * t63 + Icges(6,2) * t62 + Icges(6,6) * t98;
t39 = Icges(6,4) * t61 + Icges(6,2) * t60 + Icges(6,6) * t96;
t38 = Icges(6,5) * t63 + Icges(6,6) * t62 + Icges(6,3) * t98;
t37 = Icges(6,5) * t61 + Icges(6,6) * t60 + Icges(6,3) * t96;
t36 = -t106 * t83 + t53 * t134;
t35 = t104 * t83 - t52 * t134;
t34 = t109 * t81 + t110 * t82 + t80 * t134;
t33 = t139 * t118 - t137 * t128;
t32 = t140 * t118 + t137 * t131;
t31 = -t104 * t53 + t106 * t52;
t30 = t106 * t80 + t81 * t88 + t82 * t89;
t29 = t104 * t80 + t81 * t86 + t82 * t87;
t28 = (t139 * t117 + t140 * t119) * t116;
t27 = t109 * t49 + t110 * t51 + t47 * t134;
t26 = t109 * t48 + t110 * t50 + t46 * t134;
t25 = t141 * t118 - t138 * t128;
t24 = t142 * t118 + t138 * t131;
t23 = t106 * t47 + t49 * t88 + t51 * t89;
t22 = t106 * t46 + t48 * t88 + t50 * t89;
t21 = t104 * t47 + t49 * t86 + t51 * t87;
t20 = t104 * t46 + t48 * t86 + t50 * t87;
t19 = (t141 * t117 + t142 * t119) * t116;
t18 = -t106 * t59 - t83 * t98 + (t113 * t44 + t53 * t127) * t116;
t17 = t104 * t59 + t83 * t96 + (-t113 * t43 - t52 * t127) * t116;
t16 = t109 * t57 + t110 * t58 + t81 * t90 + t82 * t91 + (t113 * t56 + t80 * t127) * t116;
t15 = -t104 * t44 + t106 * t43 + t52 * t98 - t53 * t96;
t14 = t106 * t56 + t57 * t88 + t58 * t89 + t62 * t81 + t63 * t82 + t80 * t98;
t13 = t104 * t56 + t57 * t86 + t58 * t87 + t60 * t81 + t61 * t82 + t80 * t96;
t12 = t109 * t40 + t110 * t42 + t49 * t90 + t51 * t91 + (t113 * t38 + t47 * t127) * t116;
t11 = t109 * t39 + t110 * t41 + t48 * t90 + t50 * t91 + (t113 * t37 + t46 * t127) * t116;
t10 = t106 * t38 + t40 * t88 + t42 * t89 + t47 * t98 + t49 * t62 + t51 * t63;
t9 = t106 * t37 + t39 * t88 + t41 * t89 + t46 * t98 + t48 * t62 + t50 * t63;
t8 = t104 * t38 + t40 * t86 + t42 * t87 + t47 * t96 + t49 * t60 + t51 * t61;
t7 = t104 * t37 + t39 * t86 + t41 * t87 + t46 * t96 + t48 * t60 + t50 * t61;
t6 = -t118 * t16 + (t11 * t117 + t119 * t12) * t116;
t5 = -t118 * t14 + (t10 * t119 + t117 * t9) * t116;
t4 = -t118 * t13 + (t117 * t7 + t119 * t8) * t116;
t3 = t104 * t11 + t106 * t12 + t26 * t96 + t27 * t98 + (t113 * t16 + t34 * t127) * t116;
t2 = t10 * t106 + t104 * t9 + t22 * t96 + t98 * t23 + (t113 * t14 + t30 * t127) * t116;
t1 = t104 * t7 + t106 * t8 + t96 * t20 + t21 * t98 + (t113 * t13 + t29 * t127) * t116;
t78 = [0; 0; 0; 0; 0; 0; m(5) * t45 + m(6) * t19; m(5) * t122 + m(6) * (t117 * t24 - t119 * t25); m(5) * (-t118 * t45 + (t117 * t55 + t119 * t54) * t116) + m(6) * (-t118 * t19 + (t117 * t25 + t119 * t24) * t116); 0.2e1 * m(5) * ((t68 * t54 - t69 * t55) * t118 + ((-t117 * t69 + t119 * t68) * t45 + t122 * (-rSges(5,3) * t118 + (rSges(5,1) * t114 - rSges(5,2) * t113) * t116)) * t116) - t118 * (t118 ^ 2 * t99 + (((-t113 * t73 + t114 * t75) * t119 + (-t113 * t72 + t114 * t74) * t117 + ((-t113 * t67 - t114 * t65) * t119 + (-t113 * t66 - t114 * t64) * t117) * qJD(4)) * t116 + (t100 * t113 - t101 * t114 - t70 * t117 - t71 * t119 + (t113 * t93 + t114 * t92) * qJD(4)) * t118) * t116) + (t19 * t28 + t24 * t32 + t25 * t33) * t143 - t118 * t6 + (-(-t100 * t104 + t101 * t105 - t92 * t96 - t93 * t95) * t118 + t4 + (-t104 * t72 + t105 * t74 + t70 * t131 - t64 * t96 - t66 * t95 - t144) * t131) * t131 + (-(-t100 * t106 + t101 * t107 - t92 * t98 - t93 * t97) * t118 + t5 + (-t106 * t73 + t107 * t75 + t71 * t128 - t65 * t98 - t67 * t97 - t144) * t128 + (-t104 * t73 + t105 * t75 - t106 * t72 + t107 * t74 + t70 * t128 + t71 * t131 - t64 * t98 - t65 * t96 - t66 * t97 - t67 * t95) * t131) * t128; m(6) * t15; m(6) * (t117 * t17 - t119 * t18); m(6) * (-t118 * t15 + (t117 * t18 + t119 * t17) * t116); m(6) * (t15 * t28 + t17 * t32 + t18 * t33 + t19 * t31 + t24 * t35 + t25 * t36) + t106 * t5 / 0.2e1 + t104 * t4 / 0.2e1 + (-t98 * t30 / 0.2e1 - t96 * t29 / 0.2e1 - t3 / 0.2e1) * t118 + (t119 * t2 / 0.2e1 + t98 * (t117 * t22 + t119 * t23) / 0.2e1 + t117 * t1 / 0.2e1 + t96 * (t117 * t20 + t119 * t21) / 0.2e1 + t113 * t6 / 0.2e1 + (-t118 * t34 / 0.2e1 + (t117 * t26 + t119 * t27) * t116 / 0.2e1) * t127) * t116; (t15 * t31 + t17 * t35 + t18 * t36) * t143 + t98 * (t104 * t22 + t106 * t23 + t30 * t134) + t106 * t2 + t96 * (t104 * t20 + t106 * t21 + t29 * t134) + t104 * t1 + (t104 * t26 + t106 * t27 + t34 * t134) * t123 + t3 * t134;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t78(1), t78(2), t78(4), t78(7), t78(11); t78(2), t78(3), t78(5), t78(8), t78(12); t78(4), t78(5), t78(6), t78(9), t78(13); t78(7), t78(8), t78(9), t78(10), t78(14); t78(11), t78(12), t78(13), t78(14), t78(15);];
Mq = res;
