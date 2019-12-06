% Calculate time derivative of joint inertia matrix for
% S5PPPRR2
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPPRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPPRR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:26
% EndTime: 2019-12-05 14:59:30
% DurationCPUTime: 1.87s
% Computational Cost: add. (7615->319), mult. (21581->490), div. (0->0), fcn. (25308->10), ass. (0->135)
t125 = cos(pkin(9));
t126 = cos(pkin(8));
t129 = sin(qJ(4));
t123 = sin(pkin(8));
t146 = cos(qJ(4));
t134 = t123 * t146;
t118 = t125 * t134 - t126 * t129;
t147 = 2 * m(6);
t122 = sin(pkin(9));
t127 = cos(pkin(7));
t124 = sin(pkin(7));
t137 = t124 * t126;
t114 = -t122 * t127 + t125 * t137;
t132 = -t114 * t129 + t124 * t134;
t101 = t132 * qJD(4);
t138 = t123 * t129;
t106 = t114 * t146 + t124 * t138;
t102 = t106 * qJD(4);
t128 = sin(qJ(5));
t113 = t122 * t137 + t125 * t127;
t130 = cos(qJ(5));
t89 = t106 * t130 + t113 * t128;
t60 = -t89 * qJD(5) - t101 * t128;
t88 = -t106 * t128 + t113 * t130;
t61 = t88 * qJD(5) + t101 * t130;
t43 = rSges(6,1) * t61 + rSges(6,2) * t60 + rSges(6,3) * t102;
t145 = pkin(4) * t101 + pkin(6) * t102 + t43;
t136 = t126 * t127;
t116 = t122 * t124 + t125 * t136;
t131 = -t116 * t129 + t127 * t134;
t103 = t131 * qJD(4);
t108 = t116 * t146 + t127 * t138;
t104 = t108 * qJD(4);
t115 = t122 * t136 - t124 * t125;
t91 = t108 * t130 + t115 * t128;
t62 = -t91 * qJD(5) - t103 * t128;
t90 = -t108 * t128 + t115 * t130;
t63 = t90 * qJD(5) + t103 * t130;
t44 = rSges(6,1) * t63 + rSges(6,2) * t62 + rSges(6,3) * t104;
t144 = pkin(4) * t103 + pkin(6) * t104 + t44;
t52 = rSges(6,1) * t89 + rSges(6,2) * t88 - rSges(6,3) * t132;
t143 = pkin(4) * t106 - pkin(6) * t132 + t52;
t53 = rSges(6,1) * t91 + rSges(6,2) * t90 - rSges(6,3) * t131;
t142 = pkin(4) * t108 - pkin(6) * t131 + t53;
t117 = t125 * t138 + t126 * t146;
t111 = t117 * qJD(4);
t112 = t118 * qJD(4);
t139 = t123 * t122;
t110 = t118 * t130 + t128 * t139;
t86 = -t110 * qJD(5) + t111 * t128;
t109 = -t118 * t128 + t130 * t139;
t87 = t109 * qJD(5) - t111 * t130;
t59 = rSges(6,1) * t87 + rSges(6,2) * t86 + rSges(6,3) * t112;
t141 = -pkin(4) * t111 + pkin(6) * t112 + t59;
t73 = rSges(6,1) * t110 + rSges(6,2) * t109 + rSges(6,3) * t117;
t140 = pkin(4) * t118 + pkin(6) * t117 + t73;
t98 = -rSges(5,1) * t111 - rSges(5,2) * t112;
t97 = -Icges(5,1) * t111 - Icges(5,4) * t112;
t96 = -Icges(5,4) * t111 - Icges(5,2) * t112;
t95 = -Icges(5,5) * t111 - Icges(5,6) * t112;
t94 = rSges(5,1) * t118 - rSges(5,2) * t117 + rSges(5,3) * t139;
t93 = Icges(5,1) * t118 - Icges(5,4) * t117 + Icges(5,5) * t139;
t92 = Icges(5,4) * t118 - Icges(5,2) * t117 + Icges(5,6) * t139;
t81 = rSges(5,1) * t103 - rSges(5,2) * t104;
t80 = rSges(5,1) * t101 - rSges(5,2) * t102;
t79 = Icges(5,1) * t103 - Icges(5,4) * t104;
t78 = Icges(5,1) * t101 - Icges(5,4) * t102;
t77 = Icges(5,4) * t103 - Icges(5,2) * t104;
t76 = Icges(5,4) * t101 - Icges(5,2) * t102;
t75 = Icges(5,5) * t103 - Icges(5,6) * t104;
t74 = Icges(5,5) * t101 - Icges(5,6) * t102;
t72 = Icges(6,1) * t110 + Icges(6,4) * t109 + Icges(6,5) * t117;
t71 = Icges(6,4) * t110 + Icges(6,2) * t109 + Icges(6,6) * t117;
t70 = Icges(6,5) * t110 + Icges(6,6) * t109 + Icges(6,3) * t117;
t69 = rSges(5,1) * t108 + rSges(5,2) * t131 + rSges(5,3) * t115;
t68 = rSges(5,1) * t106 + rSges(5,2) * t132 + rSges(5,3) * t113;
t67 = Icges(5,1) * t108 + Icges(5,4) * t131 + Icges(5,5) * t115;
t66 = Icges(5,1) * t106 + Icges(5,4) * t132 + Icges(5,5) * t113;
t65 = Icges(5,4) * t108 + Icges(5,2) * t131 + Icges(5,6) * t115;
t64 = Icges(5,4) * t106 + Icges(5,2) * t132 + Icges(5,6) * t113;
t58 = Icges(6,1) * t87 + Icges(6,4) * t86 + Icges(6,5) * t112;
t57 = Icges(6,4) * t87 + Icges(6,2) * t86 + Icges(6,6) * t112;
t56 = Icges(6,5) * t87 + Icges(6,6) * t86 + Icges(6,3) * t112;
t55 = -t115 * t98 + t81 * t139;
t54 = t113 * t98 - t80 * t139;
t51 = Icges(6,1) * t91 + Icges(6,4) * t90 - Icges(6,5) * t131;
t50 = Icges(6,1) * t89 + Icges(6,4) * t88 - Icges(6,5) * t132;
t49 = Icges(6,4) * t91 + Icges(6,2) * t90 - Icges(6,6) * t131;
t48 = Icges(6,4) * t89 + Icges(6,2) * t88 - Icges(6,6) * t132;
t47 = Icges(6,5) * t91 + Icges(6,6) * t90 - Icges(6,3) * t131;
t46 = Icges(6,5) * t89 + Icges(6,6) * t88 - Icges(6,3) * t132;
t45 = -t113 * t81 + t115 * t80;
t42 = Icges(6,1) * t63 + Icges(6,4) * t62 + Icges(6,5) * t104;
t41 = Icges(6,1) * t61 + Icges(6,4) * t60 + Icges(6,5) * t102;
t40 = Icges(6,4) * t63 + Icges(6,2) * t62 + Icges(6,6) * t104;
t39 = Icges(6,4) * t61 + Icges(6,2) * t60 + Icges(6,6) * t102;
t38 = Icges(6,5) * t63 + Icges(6,6) * t62 + Icges(6,3) * t104;
t37 = Icges(6,5) * t61 + Icges(6,6) * t60 + Icges(6,3) * t102;
t36 = t117 * t53 + t131 * t73;
t35 = -t117 * t52 - t132 * t73;
t34 = t109 * t71 + t110 * t72 + t117 * t70;
t33 = -t131 * t52 + t132 * t53;
t32 = -t140 * t115 + t142 * t139;
t31 = t140 * t113 - t143 * t139;
t30 = -t131 * t70 + t71 * t90 + t72 * t91;
t29 = -t132 * t70 + t71 * t88 + t72 * t89;
t28 = -t142 * t113 + t143 * t115;
t27 = t109 * t49 + t110 * t51 + t117 * t47;
t26 = t109 * t48 + t110 * t50 + t117 * t46;
t25 = -t141 * t115 + t144 * t139;
t24 = t141 * t113 - t145 * t139;
t23 = -t131 * t47 + t49 * t90 + t51 * t91;
t22 = -t131 * t46 + t48 * t90 + t50 * t91;
t21 = -t132 * t47 + t49 * t88 + t51 * t89;
t20 = -t132 * t46 + t48 * t88 + t50 * t89;
t19 = -t144 * t113 + t145 * t115;
t18 = -t104 * t73 + t112 * t53 + t117 * t44 + t131 * t59;
t17 = t102 * t73 - t112 * t52 - t117 * t43 - t132 * t59;
t16 = t109 * t57 + t110 * t58 + t112 * t70 + t117 * t56 + t71 * t86 + t72 * t87;
t15 = -t102 * t53 + t104 * t52 - t131 * t43 + t132 * t44;
t14 = t104 * t70 - t131 * t56 + t57 * t90 + t58 * t91 + t62 * t71 + t63 * t72;
t13 = t102 * t70 - t132 * t56 + t57 * t88 + t58 * t89 + t60 * t71 + t61 * t72;
t12 = t109 * t40 + t110 * t42 + t112 * t47 + t117 * t38 + t49 * t86 + t51 * t87;
t11 = t109 * t39 + t110 * t41 + t112 * t46 + t117 * t37 + t48 * t86 + t50 * t87;
t10 = t104 * t47 - t131 * t38 + t40 * t90 + t42 * t91 + t49 * t62 + t51 * t63;
t9 = t104 * t46 - t131 * t37 + t39 * t90 + t41 * t91 + t48 * t62 + t50 * t63;
t8 = t102 * t47 - t132 * t38 + t40 * t88 + t42 * t89 + t49 * t60 + t51 * t61;
t7 = t102 * t46 - t132 * t37 + t39 * t88 + t41 * t89 + t48 * t60 + t50 * t61;
t6 = t11 * t113 + t115 * t12 + t16 * t139;
t5 = t10 * t115 + t113 * t9 + t14 * t139;
t4 = t113 * t7 + t115 * t8 + t13 * t139;
t3 = t102 * t26 + t104 * t27 - t11 * t132 + t112 * t34 + t117 * t16 - t12 * t131;
t2 = -t10 * t131 + t102 * t22 + t104 * t23 + t112 * t30 + t117 * t14 - t132 * t9;
t1 = t102 * t20 + t104 * t21 + t112 * t29 + t117 * t13 - t131 * t8 - t132 * t7;
t82 = [0; 0; 0; 0; 0; 0; m(5) * t45 + m(6) * t19; m(5) * (t124 * t54 - t127 * t55) + m(6) * (t124 * t24 - t127 * t25); m(5) * (-t126 * t45 + (t124 * t55 + t127 * t54) * t123) + m(6) * (-t126 * t19 + (t124 * t25 + t127 * t24) * t123); t115 * ((t103 * t67 - t104 * t65 + t108 * t79 + t115 * t75 + t131 * t77) * t115 + (t103 * t66 - t104 * t64 + t108 * t78 + t115 * t74 + t131 * t76) * t113 + (t103 * t93 - t104 * t92 + t108 * t97 + t115 * t95 + t131 * t96) * t139) + t113 * ((t101 * t67 - t102 * t65 + t106 * t79 + t113 * t75 + t132 * t77) * t115 + (t101 * t66 - t102 * t64 + t106 * t78 + t113 * t74 + t132 * t76) * t113 + (t101 * t93 - t102 * t92 + t106 * t97 + t113 * t95 + t132 * t96) * t139) + ((-t111 * t67 - t112 * t65 - t117 * t77 + t118 * t79 + t75 * t139) * t115 + (-t111 * t66 - t112 * t64 - t117 * t76 + t118 * t78 + t74 * t139) * t113 + (-t111 * t93 - t112 * t92 - t117 * t96 + t118 * t97 + t95 * t139) * t139) * t139 + t115 * t5 + t113 * t4 + t6 * t139 + 0.2e1 * m(5) * ((t113 * t94 - t68 * t139) * t54 + (-t115 * t94 + t69 * t139) * t55 + (-t113 * t69 + t115 * t68) * t45) + (t19 * t28 + t24 * t31 + t25 * t32) * t147; m(6) * t15; m(6) * (t124 * t17 - t127 * t18); m(6) * (-t126 * t15 + (t124 * t18 + t127 * t17) * t123); m(6) * (t15 * t28 + t17 * t31 + t18 * t32 + t19 * t33 + t24 * t35 + t25 * t36) + t115 * t2 / 0.2e1 + t104 * (t113 * t22 + t115 * t23 + t30 * t139) / 0.2e1 - t131 * t5 / 0.2e1 + t113 * t1 / 0.2e1 + t102 * (t113 * t20 + t115 * t21 + t29 * t139) / 0.2e1 - t132 * t4 / 0.2e1 + t3 * t139 / 0.2e1 + t112 * (t113 * t26 + t115 * t27 + t34 * t139) / 0.2e1 + t117 * t6 / 0.2e1; (t15 * t33 + t17 * t35 + t18 * t36) * t147 + t104 * (t117 * t30 - t131 * t23 - t132 * t22) - t131 * t2 + t102 * (t117 * t29 - t131 * t21 - t132 * t20) - t132 * t1 + t112 * (t117 * t34 - t131 * t27 - t132 * t26) + t117 * t3;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t82(1), t82(2), t82(4), t82(7), t82(11); t82(2), t82(3), t82(5), t82(8), t82(12); t82(4), t82(5), t82(6), t82(9), t82(13); t82(7), t82(8), t82(9), t82(10), t82(14); t82(11), t82(12), t82(13), t82(14), t82(15);];
Mq = res;
