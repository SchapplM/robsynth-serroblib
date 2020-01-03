% Calculate time derivative of joint inertia matrix for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR4_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR4_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:27
% EndTime: 2019-12-31 16:18:31
% DurationCPUTime: 1.46s
% Computational Cost: add. (3557->208), mult. (5486->357), div. (0->0), fcn. (5253->6), ass. (0->121)
t83 = cos(pkin(6));
t141 = t83 ^ 2;
t82 = sin(pkin(6));
t147 = t82 ^ 2 + t141;
t145 = t147 * qJD(3);
t143 = 2 * m(5);
t140 = t82 / 0.2e1;
t139 = t83 / 0.2e1;
t81 = pkin(7) + qJ(3);
t79 = sin(t81);
t134 = t79 * t82;
t118 = Icges(5,3) * t79;
t84 = sin(qJ(4));
t128 = t83 * t84;
t85 = cos(qJ(4));
t129 = t82 * t85;
t80 = cos(t81);
t73 = -t128 * t80 + t129;
t127 = t83 * t85;
t130 = t82 * t84;
t74 = t127 * t80 + t130;
t44 = Icges(5,5) * t74 + Icges(5,6) * t73 + t118 * t83;
t119 = Icges(5,6) * t79;
t46 = Icges(5,4) * t74 + Icges(5,2) * t73 + t119 * t83;
t120 = Icges(5,5) * t79;
t48 = Icges(5,1) * t74 + Icges(5,4) * t73 + t120 * t83;
t71 = -t130 * t80 - t127;
t72 = t129 * t80 - t128;
t16 = t134 * t44 + t46 * t71 + t48 * t72;
t138 = t16 * t83;
t133 = t79 * t83;
t43 = Icges(5,5) * t72 + Icges(5,6) * t71 + t118 * t82;
t45 = Icges(5,4) * t72 + Icges(5,2) * t71 + t119 * t82;
t47 = Icges(5,1) * t72 + Icges(5,4) * t71 + t120 * t82;
t17 = t133 * t43 + t45 * t73 + t47 * t74;
t137 = t17 * t82;
t105 = rSges(5,1) * t85 - rSges(5,2) * t84;
t60 = -t80 * rSges(5,3) + t105 * t79;
t136 = t60 * t80;
t115 = qJD(4) * t79;
t92 = Icges(5,5) * t85 - Icges(5,6) * t84;
t132 = t80 * ((-Icges(5,5) * t84 - Icges(5,6) * t85) * t115 + (t80 * t92 + t118) * qJD(3));
t57 = -Icges(5,3) * t80 + t79 * t92;
t131 = t80 * t57;
t107 = pkin(3) * t80 + pkin(5) * t79;
t42 = (-rSges(5,1) * t84 - rSges(5,2) * t85) * t115 + (rSges(5,3) * t79 + t105 * t80) * qJD(3);
t126 = -qJD(3) * t107 - t42;
t78 = pkin(3) * t79 - pkin(5) * t80;
t125 = -t60 - t78;
t122 = Icges(5,4) * t84;
t121 = Icges(5,4) * t85;
t117 = qJD(3) * t79;
t116 = qJD(3) * t80;
t114 = t84 * t117;
t113 = t85 * t117;
t112 = t82 * t116;
t111 = t83 * t116;
t110 = Icges(5,5) * t116;
t109 = Icges(5,6) * t116;
t108 = Icges(5,3) * t116;
t77 = rSges(4,1) * t79 + rSges(4,2) * t80;
t103 = -t45 * t84 + t47 * t85;
t19 = t103 * t79 - t80 * t43;
t102 = -t46 * t84 + t48 * t85;
t20 = t102 * t79 - t80 * t44;
t104 = t19 * t82 + t20 * t83;
t49 = rSges(5,1) * t72 + rSges(5,2) * t71 + rSges(5,3) * t134;
t50 = rSges(5,1) * t74 + rSges(5,2) * t73 + rSges(5,3) * t133;
t101 = t49 * t83 - t50 * t82;
t93 = -Icges(5,2) * t84 + t121;
t58 = -Icges(5,6) * t80 + t79 * t93;
t95 = Icges(5,1) * t85 - t122;
t59 = -Icges(5,5) * t80 + t79 * t95;
t100 = t58 * t84 - t59 * t85;
t53 = -qJD(4) * t72 + t114 * t82;
t54 = qJD(4) * t71 - t113 * t82;
t28 = Icges(5,5) * t54 + Icges(5,6) * t53 + t108 * t82;
t91 = t116 * t43 + t28 * t79;
t55 = -qJD(4) * t74 + t114 * t83;
t56 = qJD(4) * t73 - t113 * t83;
t29 = Icges(5,5) * t56 + Icges(5,6) * t55 + t108 * t83;
t90 = t116 * t44 + t29 * t79;
t86 = qJD(3) * (-Icges(4,5) * t79 - Icges(4,6) * t80);
t65 = t82 * t86;
t52 = t125 * t83;
t51 = t125 * t82;
t41 = (-Icges(5,1) * t84 - t121) * t115 + (t80 * t95 + t120) * qJD(3);
t40 = (-Icges(5,2) * t85 - t122) * t115 + (t80 * t93 + t119) * qJD(3);
t38 = t77 * t145;
t37 = t126 * t83;
t36 = t126 * t82;
t35 = rSges(5,1) * t56 + rSges(5,2) * t55 + rSges(5,3) * t111;
t34 = rSges(5,1) * t54 + rSges(5,2) * t53 + rSges(5,3) * t112;
t33 = Icges(5,1) * t56 + Icges(5,4) * t55 + t110 * t83;
t32 = Icges(5,1) * t54 + Icges(5,4) * t53 + t110 * t82;
t31 = Icges(5,4) * t56 + Icges(5,2) * t55 + t109 * t83;
t30 = Icges(5,4) * t54 + Icges(5,2) * t53 + t109 * t82;
t27 = -t133 * t60 - t50 * t80;
t26 = t134 * t60 + t49 * t80;
t25 = -t100 * t79 - t131;
t24 = t101 * t79;
t23 = t133 * t57 + t58 * t73 + t59 * t74;
t22 = t134 * t57 + t58 * t71 + t59 * t72;
t21 = (t107 * t83 + t50) * t83 + (t107 * t82 + t49) * t82;
t18 = t133 * t44 + t46 * t73 + t48 * t74;
t15 = t134 * t43 + t45 * t71 + t47 * t72;
t14 = -t145 * t78 + t82 * t34 + t83 * t35;
t13 = -t42 * t133 - t35 * t80 + (-t136 * t83 + t50 * t79) * qJD(3);
t12 = t42 * t134 + t34 * t80 + (t136 * t82 - t49 * t79) * qJD(3);
t11 = (t34 * t83 - t35 * t82) * t79 + t101 * t116;
t10 = t31 * t73 + t33 * t74 + t46 * t55 + t48 * t56 + t83 * t90;
t9 = t30 * t73 + t32 * t74 + t45 * t55 + t47 * t56 + t83 * t91;
t8 = t31 * t71 + t33 * t72 + t46 * t53 + t48 * t54 + t82 * t90;
t7 = t30 * t71 + t32 * t72 + t45 * t53 + t47 * t54 + t82 * t91;
t6 = (qJD(3) * t102 - t29) * t80 + (qJD(3) * t44 - t31 * t84 + t33 * t85 + (-t46 * t85 - t48 * t84) * qJD(4)) * t79;
t5 = (qJD(3) * t103 - t28) * t80 + (qJD(3) * t43 - t30 * t84 + t32 * t85 + (-t45 * t85 - t47 * t84) * qJD(4)) * t79;
t4 = t10 * t82 - t83 * t9;
t3 = -t7 * t83 + t8 * t82;
t2 = -(t40 * t73 + t41 * t74 + t55 * t58 + t56 * t59) * t80 + (t9 * t82 + (t10 - t132) * t83) * t79 + (t23 * t79 + (t137 + (t18 - t131) * t83) * t80) * qJD(3);
t1 = -(t40 * t71 + t41 * t72 + t53 * t58 + t54 * t59) * t80 + (t8 * t83 + (t7 - t132) * t82) * t79 + (t22 * t79 + (t138 + (t15 - t131) * t82) * t80) * qJD(3);
t39 = [0; 0; 0; -m(4) * t38 + m(5) * t14; m(5) * (-t36 * t83 + t37 * t82); (t14 * t21 + t36 * t51 + t37 * t52) * t143 + (-t141 * t65 - t3) * t83 + (t4 + (t147 * t86 - t82 * t65) * t83) * t82 + 0.2e1 * m(4) * (qJD(3) * t77 - t38) * t147 * (rSges(4,1) * t80 - rSges(4,2) * t79); m(5) * t11; m(5) * (t12 * t82 - t13 * t83); m(5) * (t11 * t21 + t12 * t52 + t13 * t51 + t14 * t24 + t26 * t37 + t27 * t36) + t2 * t140 - t83 * t1 / 0.2e1 - t80 * (-t5 * t83 + t6 * t82) / 0.2e1 + (t4 * t139 + t3 * t140) * t79 + (t79 * (-t19 * t83 + t20 * t82) / 0.2e1 + ((-t17 * t83 + t18 * t82) * t139 + (-t15 * t83 + t82 * t16) * t140) * t80) * qJD(3); (t11 * t24 + t12 * t26 + t13 * t27) * t143 + (-t23 * t80 + (t18 * t83 + t137) * t79) * t111 + t2 * t133 + (-t22 * t80 + (t15 * t82 + t138) * t79) * t112 + t1 * t134 + (t104 * t79 - t25 * t80) * t117 - t80 * ((t132 + (t100 * t80 + t104) * qJD(3)) * t80 + (t6 * t83 + t5 * t82 - (qJD(3) * t57 - t40 * t84 + t41 * t85 + (-t58 * t85 - t59 * t84) * qJD(4)) * t80 + t25 * qJD(3)) * t79);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t39(1), t39(2), t39(4), t39(7); t39(2), t39(3), t39(5), t39(8); t39(4), t39(5), t39(6), t39(9); t39(7), t39(8), t39(9), t39(10);];
Mq = res;
