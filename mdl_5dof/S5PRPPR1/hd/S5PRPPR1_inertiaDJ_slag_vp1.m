% Calculate time derivative of joint inertia matrix for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:47
% EndTime: 2019-12-05 15:21:52
% DurationCPUTime: 1.29s
% Computational Cost: add. (3741->235), mult. (3750->381), div. (0->0), fcn. (3483->8), ass. (0->117)
t81 = pkin(9) + qJ(5);
t79 = cos(t81);
t117 = Icges(6,4) * t79;
t77 = sin(t81);
t84 = sin(pkin(8));
t86 = cos(pkin(8));
t57 = -Icges(6,6) * t86 + (-Icges(6,2) * t77 + t117) * t84;
t118 = Icges(6,4) * t77;
t58 = -Icges(6,5) * t86 + (Icges(6,1) * t79 - t118) * t84;
t134 = -t57 * t79 - t58 * t77;
t83 = sin(pkin(9));
t85 = cos(pkin(9));
t133 = (rSges(5,3) + qJ(4)) * t84 + (rSges(5,1) * t85 - rSges(5,2) * t83 + pkin(3)) * t86;
t132 = 2 * m(6);
t131 = pkin(4) * t83;
t108 = qJD(5) * t84;
t64 = (-Icges(6,5) * t77 - Icges(6,6) * t79) * t108;
t66 = (-Icges(6,1) * t77 - t117) * t108;
t101 = t84 * t79 * t66 - t86 * t64;
t65 = (-Icges(6,2) * t79 - t118) * t108;
t126 = t77 * t65;
t130 = ((qJD(5) * t134 - t126) * t84 + t101) * t86;
t127 = (pkin(4) * t85 + pkin(3)) * t86;
t82 = pkin(7) + qJ(2);
t78 = sin(t82);
t125 = t78 * t84;
t124 = t78 * t86;
t80 = cos(t82);
t123 = t80 * t84;
t122 = t80 * t86;
t63 = t79 * t122 + t77 * t78;
t93 = t77 * t124 + t79 * t80;
t45 = t93 * qJD(2) - t63 * qJD(5);
t61 = t79 * t124 - t77 * t80;
t62 = -t77 * t122 + t78 * t79;
t46 = -t61 * qJD(2) + t62 * qJD(5);
t121 = t46 * rSges(6,1) + t45 * rSges(6,2);
t111 = qJD(2) * t80;
t120 = qJ(3) * t111 + qJD(3) * t78;
t119 = t80 * pkin(2) + t78 * qJ(3);
t116 = Icges(6,5) * t84;
t115 = Icges(6,6) * t84;
t114 = Icges(6,3) * t84;
t112 = qJD(2) * t78;
t110 = qJD(2) * t84;
t109 = qJD(4) * t84;
t42 = t63 * rSges(6,1) + t62 * rSges(6,2) + rSges(6,3) * t123;
t107 = t80 * t109 + t120;
t106 = t78 * t110;
t105 = t80 * t110;
t104 = Icges(6,5) * t110;
t103 = Icges(6,6) * t110;
t102 = Icges(6,3) * t110;
t72 = qJD(3) * t80;
t100 = -t78 * t109 + t72;
t99 = rSges(4,1) * t86 - rSges(4,2) * t84;
t98 = t83 * rSges(5,1) + t85 * rSges(5,2);
t47 = t62 * qJD(2) - t61 * qJD(5);
t48 = t63 * qJD(2) - t93 * qJD(5);
t97 = -t48 * rSges(6,1) - t47 * rSges(6,2);
t96 = -t61 * rSges(6,1) + rSges(6,2) * t93;
t95 = -pkin(2) - t99;
t87 = -pkin(6) - qJ(4);
t92 = -t127 - pkin(2) + (-rSges(6,3) + t87) * t84;
t91 = rSges(4,3) * t80 + t95 * t78;
t90 = t80 * t131 + t92 * t78;
t89 = -pkin(2) - t133;
t88 = t89 * t78 + t98 * t80;
t74 = t80 * qJ(3);
t67 = (-rSges(6,1) * t77 - rSges(6,2) * t79) * t108;
t59 = -rSges(6,3) * t86 + (rSges(6,1) * t79 - rSges(6,2) * t77) * t84;
t56 = -Icges(6,3) * t86 + (Icges(6,5) * t79 - Icges(6,6) * t77) * t84;
t52 = rSges(4,3) * t78 + t99 * t80 + t119;
t51 = t74 + t91;
t50 = t72 + ((-rSges(4,3) - qJ(3)) * t78 + t95 * t80) * qJD(2);
t49 = t91 * qJD(2) + t120;
t41 = rSges(6,3) * t125 - t96;
t40 = Icges(6,1) * t63 + Icges(6,4) * t62 + t80 * t116;
t39 = Icges(6,1) * t61 - Icges(6,4) * t93 + t78 * t116;
t38 = Icges(6,4) * t63 + Icges(6,2) * t62 + t80 * t115;
t37 = Icges(6,4) * t61 - Icges(6,2) * t93 + t78 * t115;
t36 = Icges(6,5) * t63 + Icges(6,6) * t62 + t80 * t114;
t35 = Icges(6,5) * t61 - Icges(6,6) * t93 + t78 * t114;
t34 = t133 * t80 + t98 * t78 + t119;
t33 = t74 + t88;
t32 = -t59 * t123 - t42 * t86;
t31 = t59 * t125 + t41 * t86;
t30 = t78 * t131 + (-t84 * t87 + t127) * t80 + t42 + t119;
t29 = t74 + t90 + t96;
t28 = ((-qJ(3) - t98) * t78 + t89 * t80) * qJD(2) + t100;
t27 = t88 * qJD(2) + t107;
t26 = rSges(6,3) * t105 - t97;
t25 = -rSges(6,3) * t106 + t121;
t24 = Icges(6,1) * t48 + Icges(6,4) * t47 + t80 * t104;
t23 = Icges(6,1) * t46 + Icges(6,4) * t45 - t78 * t104;
t22 = Icges(6,4) * t48 + Icges(6,2) * t47 + t80 * t103;
t21 = Icges(6,4) * t46 + Icges(6,2) * t45 - t78 * t103;
t20 = Icges(6,5) * t48 + Icges(6,6) * t47 + t80 * t102;
t19 = Icges(6,5) * t46 + Icges(6,6) * t45 - t78 * t102;
t18 = t56 * t123 + t57 * t62 + t58 * t63;
t17 = t56 * t125 - t57 * t93 + t58 * t61;
t15 = t26 * t86 + (t59 * t111 + t67 * t78) * t84;
t14 = -t25 * t86 + (t59 * t112 - t67 * t80) * t84;
t13 = ((-qJ(3) - t131) * t78 + t92 * t80) * qJD(2) + t97 + t100;
t12 = t90 * qJD(2) + t107 + t121;
t11 = -t36 * t86 + (-t38 * t77 + t40 * t79) * t84;
t10 = -t35 * t86 + (-t37 * t77 + t39 * t79) * t84;
t9 = t36 * t123 + t38 * t62 + t40 * t63;
t8 = t35 * t123 + t37 * t62 + t39 * t63;
t7 = t36 * t125 - t38 * t93 + t40 * t61;
t6 = t35 * t125 - t37 * t93 + t39 * t61;
t5 = t47 * t57 + t48 * t58 - t93 * t65 + t61 * t66 + (t56 * t111 + t64 * t78) * t84;
t4 = t45 * t57 + t46 * t58 + t62 * t65 + t63 * t66 + (-t56 * t112 + t64 * t80) * t84;
t3 = (-t25 * t78 + t26 * t80 + (-t41 * t78 - t42 * t80) * qJD(2)) * t84;
t2 = -t19 * t86 + (-t21 * t77 + t23 * t79 + (-t38 * t79 - t40 * t77) * qJD(5)) * t84;
t1 = -t20 * t86 + (-t22 * t77 + t24 * t79 + (-t37 * t79 - t39 * t77) * qJD(5)) * t84;
t16 = [0; 0; -t84 * t126 + (t12 * t30 + t13 * t29) * t132 + 0.2e1 * m(5) * (t27 * t34 + t28 * t33) + 0.2e1 * m(4) * (t49 * t52 + t50 * t51) + t101 + t134 * t108; 0; m(6) * (-t12 * t80 + t13 * t78 + (t29 * t80 + t30 * t78) * qJD(2)) + m(5) * (-t27 * t80 + t28 * t78 + (t33 * t80 + t34 * t78) * qJD(2)) + m(4) * (-t49 * t80 + t50 * t78 + (t51 * t80 + t52 * t78) * qJD(2)); 0; 0; 0.2e1 * (m(6) * (t30 * t111 - t29 * t112 + t12 * t78 + t13 * t80) / 0.2e1 + m(5) * (t34 * t111 - t33 * t112 + t27 * t78 + t28 * t80) / 0.2e1) * t84; 0; 0; m(6) * t3; -t130 + m(6) * (t12 * t32 + t13 * t31 + t14 * t30 + t15 * t29) + ((t2 / 0.2e1 + t4 / 0.2e1) * t80 + (t1 / 0.2e1 + t5 / 0.2e1) * t78 + ((t10 / 0.2e1 + t17 / 0.2e1) * t80 + (-t11 / 0.2e1 - t18 / 0.2e1) * t78) * qJD(2)) * t84; m(6) * (-t14 * t80 + t15 * t78 + (t31 * t80 + t32 * t78) * qJD(2)); m(6) * (-t3 * t86 + (t14 * t78 + t15 * t80 + (-t31 * t78 + t32 * t80) * qJD(2)) * t84); (t32 * t14 + t31 * t15 + (t41 * t80 - t42 * t78) * t3 * t84) * t132 - (-t18 * t86 + (t78 * t8 + t80 * t9) * t84) * t106 + (-t4 * t86 + ((t62 * t21 + t63 * t23 + t45 * t38 + t46 * t40) * t80 - t9 * t112 + (t62 * t22 + t63 * t24 + t45 * t37 + t46 * t39) * t78 + t8 * t111 + ((-t36 * t112 + t19 * t80) * t80 + (-t35 * t112 + t20 * t80) * t78) * t84) * t84) * t123 + (-t17 * t86 + (t6 * t78 + t7 * t80) * t84) * t105 + (-t5 * t86 + ((-t21 * t93 + t61 * t23 + t47 * t38 + t48 * t40) * t80 - t7 * t112 + (-t22 * t93 + t61 * t24 + t47 * t37 + t48 * t39) * t78 + t6 * t111 + ((t36 * t111 + t19 * t78) * t80 + (t35 * t111 + t20 * t78) * t78) * t84) * t84) * t125 - t86 * (-t130 + (t1 * t78 + t2 * t80 + (t10 * t80 - t11 * t78) * qJD(2)) * t84);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t16(1), t16(2), t16(4), t16(7), t16(11); t16(2), t16(3), t16(5), t16(8), t16(12); t16(4), t16(5), t16(6), t16(9), t16(13); t16(7), t16(8), t16(9), t16(10), t16(14); t16(11), t16(12), t16(13), t16(14), t16(15);];
Mq = res;
