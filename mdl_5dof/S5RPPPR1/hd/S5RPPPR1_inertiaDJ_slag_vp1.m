% Calculate time derivative of joint inertia matrix for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:13
% EndTime: 2022-01-20 09:12:19
% DurationCPUTime: 1.60s
% Computational Cost: add. (3787->242), mult. (3828->383), div. (0->0), fcn. (3531->10), ass. (0->119)
t82 = pkin(9) + qJ(5);
t79 = cos(t82);
t123 = Icges(6,4) * t79;
t77 = sin(t82);
t85 = sin(pkin(8));
t87 = cos(pkin(8));
t57 = -Icges(6,6) * t87 + (-Icges(6,2) * t77 + t123) * t85;
t124 = Icges(6,4) * t77;
t58 = -Icges(6,5) * t87 + (Icges(6,1) * t79 - t124) * t85;
t140 = -t57 * t79 - t58 * t77;
t84 = sin(pkin(9));
t86 = cos(pkin(9));
t139 = (rSges(5,3) + qJ(4)) * t85 + (rSges(5,1) * t86 - rSges(5,2) * t84 + pkin(3)) * t87;
t138 = 2 * m(6);
t137 = pkin(4) * t84;
t136 = sin(qJ(1)) * pkin(1);
t81 = cos(qJ(1)) * pkin(1);
t114 = qJD(5) * t85;
t64 = (-Icges(6,5) * t77 - Icges(6,6) * t79) * t114;
t66 = (-Icges(6,1) * t77 - t123) * t114;
t106 = t85 * t79 * t66 - t87 * t64;
t65 = (-Icges(6,2) * t79 - t124) * t114;
t131 = t77 * t65;
t135 = ((t140 * qJD(5) - t131) * t85 + t106) * t87;
t132 = (pkin(4) * t86 + pkin(3)) * t87;
t83 = qJ(1) + pkin(7);
t78 = sin(t83);
t130 = t78 * t85;
t129 = t78 * t87;
t80 = cos(t83);
t128 = t80 * t85;
t127 = t80 * t87;
t63 = t79 * t127 + t77 * t78;
t97 = t77 * t129 + t79 * t80;
t47 = t97 * qJD(1) - t63 * qJD(5);
t61 = t79 * t129 - t77 * t80;
t62 = -t77 * t127 + t78 * t79;
t48 = -t61 * qJD(1) + t62 * qJD(5);
t126 = t48 * rSges(6,1) + t47 * rSges(6,2);
t117 = qJD(1) * t80;
t125 = qJ(3) * t117 + qJD(3) * t78;
t122 = Icges(6,5) * t85;
t121 = Icges(6,6) * t85;
t120 = Icges(6,3) * t85;
t118 = qJD(1) * t78;
t116 = qJD(1) * t85;
t115 = qJD(4) * t85;
t42 = t63 * rSges(6,1) + t62 * rSges(6,2) + rSges(6,3) * t128;
t113 = t80 * t115 + t125;
t112 = t80 * pkin(2) + t78 * qJ(3) + t81;
t111 = t78 * t116;
t110 = t80 * t116;
t109 = Icges(6,5) * t116;
t108 = Icges(6,6) * t116;
t107 = Icges(6,3) * t116;
t72 = qJD(3) * t80;
t105 = -t78 * t115 + t72;
t103 = rSges(4,1) * t87 - rSges(4,2) * t85;
t102 = t84 * rSges(5,1) + t86 * rSges(5,2);
t49 = t62 * qJD(1) - t61 * qJD(5);
t50 = t63 * qJD(1) - t97 * qJD(5);
t101 = -t50 * rSges(6,1) - t49 * rSges(6,2);
t100 = -rSges(6,1) * t61 + rSges(6,2) * t97;
t99 = -pkin(2) - t103;
t88 = -pkin(6) - qJ(4);
t96 = -t132 - pkin(2) + (-rSges(6,3) + t88) * t85;
t94 = rSges(4,3) * t80 + t99 * t78 - t136;
t93 = -pkin(2) - t139;
t92 = t80 * t137 + t96 * t78 - t136;
t91 = t102 * t80 + t93 * t78 - t136;
t74 = t80 * qJ(3);
t67 = (-rSges(6,1) * t77 - rSges(6,2) * t79) * t114;
t59 = -rSges(6,3) * t87 + (rSges(6,1) * t79 - rSges(6,2) * t77) * t85;
t56 = -Icges(6,3) * t87 + (Icges(6,5) * t79 - Icges(6,6) * t77) * t85;
t52 = rSges(4,3) * t78 + t103 * t80 + t112;
t51 = t74 + t94;
t44 = t72 + (-t81 + (-rSges(4,3) - qJ(3)) * t78 + t99 * t80) * qJD(1);
t43 = t94 * qJD(1) + t125;
t41 = rSges(6,3) * t130 - t100;
t40 = Icges(6,1) * t63 + Icges(6,4) * t62 + t80 * t122;
t39 = Icges(6,1) * t61 - Icges(6,4) * t97 + t78 * t122;
t38 = Icges(6,4) * t63 + Icges(6,2) * t62 + t80 * t121;
t37 = Icges(6,4) * t61 - Icges(6,2) * t97 + t78 * t121;
t36 = Icges(6,5) * t63 + Icges(6,6) * t62 + t80 * t120;
t35 = Icges(6,5) * t61 - Icges(6,6) * t97 + t78 * t120;
t34 = t102 * t78 + t139 * t80 + t112;
t33 = t74 + t91;
t32 = -t59 * t128 - t42 * t87;
t31 = t59 * t130 + t41 * t87;
t30 = t78 * t137 + (-t85 * t88 + t132) * t80 + t112 + t42;
t29 = t100 + t74 + t92;
t28 = (-t81 + (-qJ(3) - t102) * t78 + t93 * t80) * qJD(1) + t105;
t27 = t91 * qJD(1) + t113;
t26 = rSges(6,3) * t110 - t101;
t25 = -rSges(6,3) * t111 + t126;
t24 = Icges(6,1) * t50 + Icges(6,4) * t49 + t80 * t109;
t23 = Icges(6,1) * t48 + Icges(6,4) * t47 - t78 * t109;
t22 = Icges(6,4) * t50 + Icges(6,2) * t49 + t80 * t108;
t21 = Icges(6,4) * t48 + Icges(6,2) * t47 - t78 * t108;
t20 = Icges(6,5) * t50 + Icges(6,6) * t49 + t80 * t107;
t19 = Icges(6,5) * t48 + Icges(6,6) * t47 - t78 * t107;
t18 = t56 * t128 + t57 * t62 + t58 * t63;
t17 = t56 * t130 - t57 * t97 + t58 * t61;
t15 = t26 * t87 + (t59 * t117 + t67 * t78) * t85;
t14 = -t25 * t87 + (t59 * t118 - t67 * t80) * t85;
t13 = (-t81 + (-qJ(3) - t137) * t78 + t96 * t80) * qJD(1) + t101 + t105;
t12 = t92 * qJD(1) + t113 + t126;
t11 = -t36 * t87 + (-t38 * t77 + t40 * t79) * t85;
t10 = -t35 * t87 + (-t37 * t77 + t39 * t79) * t85;
t9 = t36 * t128 + t38 * t62 + t40 * t63;
t8 = t35 * t128 + t37 * t62 + t39 * t63;
t7 = t36 * t130 - t38 * t97 + t40 * t61;
t6 = t35 * t130 - t37 * t97 + t39 * t61;
t5 = t49 * t57 + t50 * t58 - t97 * t65 + t61 * t66 + (t56 * t117 + t64 * t78) * t85;
t4 = t47 * t57 + t48 * t58 + t62 * t65 + t63 * t66 + (-t56 * t118 + t64 * t80) * t85;
t3 = (-t25 * t78 + t26 * t80 + (-t41 * t78 - t42 * t80) * qJD(1)) * t85;
t2 = -t19 * t87 + (-t21 * t77 + t23 * t79 + (-t38 * t79 - t40 * t77) * qJD(5)) * t85;
t1 = -t20 * t87 + (-t22 * t77 + t24 * t79 + (-t37 * t79 - t39 * t77) * qJD(5)) * t85;
t16 = [(t12 * t30 + t13 * t29) * t138 + 0.2e1 * m(5) * (t27 * t34 + t28 * t33) - t85 * t131 + 0.2e1 * m(4) * (t43 * t52 + t44 * t51) + t106 + t140 * t114; 0; 0; m(6) * (-t12 * t80 + t13 * t78 + (t29 * t80 + t30 * t78) * qJD(1)) + m(5) * (-t27 * t80 + t28 * t78 + (t33 * t80 + t34 * t78) * qJD(1)) + m(4) * (-t43 * t80 + t44 * t78 + (t51 * t80 + t52 * t78) * qJD(1)); 0; 0; 0.2e1 * (m(6) * (t30 * t117 - t29 * t118 + t12 * t78 + t13 * t80) / 0.2e1 + m(5) * (t34 * t117 - t33 * t118 + t27 * t78 + t28 * t80) / 0.2e1) * t85; 0; 0; 0; m(6) * (t12 * t32 + t13 * t31 + t14 * t30 + t15 * t29) - t135 + ((t2 / 0.2e1 + t4 / 0.2e1) * t80 + (t1 / 0.2e1 + t5 / 0.2e1) * t78 + ((t10 / 0.2e1 + t17 / 0.2e1) * t80 + (-t11 / 0.2e1 - t18 / 0.2e1) * t78) * qJD(1)) * t85; m(6) * t3; m(6) * (-t14 * t80 + t15 * t78 + (t31 * t80 + t32 * t78) * qJD(1)); m(6) * (-t3 * t87 + (t14 * t78 + t15 * t80 + (-t31 * t78 + t32 * t80) * qJD(1)) * t85); (t32 * t14 + t31 * t15 + (t41 * t80 - t42 * t78) * t3 * t85) * t138 - (-t18 * t87 + (t78 * t8 + t80 * t9) * t85) * t111 + (-t4 * t87 + ((t62 * t21 + t63 * t23 + t47 * t38 + t48 * t40) * t80 - t9 * t118 + (t62 * t22 + t63 * t24 + t47 * t37 + t48 * t39) * t78 + t8 * t117 + ((-t36 * t118 + t19 * t80) * t80 + (-t35 * t118 + t20 * t80) * t78) * t85) * t85) * t128 + (-t17 * t87 + (t6 * t78 + t7 * t80) * t85) * t110 + (-t5 * t87 + ((-t21 * t97 + t61 * t23 + t49 * t38 + t50 * t40) * t80 - t7 * t118 + (-t22 * t97 + t61 * t24 + t49 * t37 + t50 * t39) * t78 + t6 * t117 + ((t36 * t117 + t19 * t78) * t80 + (t35 * t117 + t20 * t78) * t78) * t85) * t85) * t130 - t87 * (-t135 + (t1 * t78 + t2 * t80 + (t10 * t80 - t11 * t78) * qJD(1)) * t85);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t16(1), t16(2), t16(4), t16(7), t16(11); t16(2), t16(3), t16(5), t16(8), t16(12); t16(4), t16(5), t16(6), t16(9), t16(13); t16(7), t16(8), t16(9), t16(10), t16(14); t16(11), t16(12), t16(13), t16(14), t16(15);];
Mq = res;
