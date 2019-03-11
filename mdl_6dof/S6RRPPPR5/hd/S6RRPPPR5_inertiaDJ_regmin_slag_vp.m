% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:23:44
% EndTime: 2019-03-09 08:23:49
% DurationCPUTime: 1.42s
% Computational Cost: add. (840->197), mult. (2129->357), div. (0->0), fcn. (1682->6), ass. (0->109)
t79 = sin(pkin(9));
t77 = t79 ^ 2;
t80 = cos(pkin(9));
t134 = (t80 ^ 2 + t77) * qJD(3);
t81 = sin(qJ(6));
t83 = cos(qJ(6));
t44 = t81 * t79 - t83 * t80;
t82 = sin(qJ(2));
t31 = t44 * t82;
t133 = pkin(4) + qJ(3);
t124 = t80 * t82;
t84 = cos(qJ(2));
t131 = pkin(4) * t124 + t84 * qJ(5);
t120 = pkin(3) + qJ(5);
t130 = t120 * t80;
t70 = t82 * qJD(2);
t129 = -qJ(4) * t70 + t84 * qJD(4);
t42 = -0.2e1 * t134;
t108 = t79 * qJD(4);
t47 = t80 * qJD(5) + t108;
t128 = 0.2e1 * t47;
t127 = t79 * t82;
t126 = t79 * t84;
t33 = -t82 * qJD(3) + (pkin(2) * t82 - qJ(3) * t84) * qJD(2);
t125 = t80 * t33;
t123 = t80 * t84;
t119 = -pkin(5) - qJ(4);
t106 = t84 * qJD(2);
t58 = t80 * t106;
t118 = -qJ(4) * t58 - qJD(4) * t124;
t117 = pkin(3) * t127 + t82 * pkin(7);
t113 = t82 * qJ(3);
t93 = -t84 * pkin(2) - t113;
t48 = -pkin(1) + t93;
t63 = pkin(7) * t123;
t27 = t79 * t48 + t63;
t116 = qJ(3) * t134;
t50 = t133 * t79;
t51 = t133 * t80;
t115 = qJ(4) * t80;
t114 = qJ(5) * t79;
t112 = qJD(6) * t81;
t111 = qJD(6) * t83;
t110 = qJD(6) * t84;
t109 = t77 * qJD(4);
t67 = t79 * qJD(3);
t107 = t80 * qJD(3);
t61 = pkin(7) * t126;
t103 = -0.2e1 * pkin(1) * qJD(2);
t102 = pkin(7) * t70;
t101 = pkin(7) * t106;
t100 = t79 * (-pkin(4) - pkin(8));
t56 = t79 * t106;
t99 = t82 * t106;
t98 = -pkin(7) * t79 - pkin(3);
t97 = t79 * qJ(4) + pkin(2);
t26 = t80 * t48 - t61;
t96 = pkin(4) * t58 + t84 * qJD(5) - t125;
t76 = t84 * pkin(3);
t25 = -t26 + t76;
t24 = t84 * qJ(4) - t27;
t28 = t79 * t33;
t95 = t28 - t129;
t94 = t120 * t79 + pkin(7);
t11 = t61 + t76 + (pkin(8) * t82 - t48) * t80 + t131;
t12 = t82 * t100 + t119 * t84 + t27;
t92 = t83 * t11 + t81 * t12;
t91 = t81 * t11 - t83 * t12;
t21 = t79 * t102 + t125;
t22 = -t80 * t102 + t28;
t90 = -t21 * t79 + t22 * t80;
t39 = t79 * pkin(8) + t50;
t40 = t80 * pkin(8) + t51;
t89 = t83 * t39 + t81 * t40;
t88 = t81 * t39 - t83 * t40;
t43 = t83 * t79 + t81 * t80;
t87 = qJD(5) * t127 + t118;
t86 = (-qJ(5) + t98) * t82;
t34 = t43 * qJD(6);
t45 = -t80 * pkin(3) - t97;
t85 = (-t45 * t84 + t113) * qJD(2);
t55 = t84 * t67;
t38 = t83 * t110 - t81 * t70;
t37 = t81 * t110 + t83 * t70;
t36 = t97 + t130;
t35 = t80 * t111 - t79 * t112;
t32 = t43 * t82;
t30 = -t82 * t115 + t117;
t29 = t119 * t79 - pkin(2) - t130;
t23 = (-t114 + t115) * t82 - t117;
t20 = (pkin(3) * t79 + pkin(7)) * t106 + t118;
t19 = -pkin(4) * t127 - t24;
t18 = (t119 * t80 + t114) * t82 + t117;
t17 = t98 * t70 - t125;
t16 = t25 + t131;
t15 = -t22 + t129;
t14 = -qJD(6) * t31 + t43 * t106;
t13 = t82 * t34 + t81 * t56 - t83 * t58;
t10 = t94 * t106 + t87;
t9 = (-pkin(4) * t126 - pkin(7) * t124) * qJD(2) + t95;
t8 = qJD(2) * t86 + t96;
t7 = -t44 * qJD(3) - t89 * qJD(6);
t6 = -t43 * qJD(3) + t88 * qJD(6);
t5 = (-pkin(5) * t80 + t94) * t106 + t87;
t4 = ((-pkin(7) * t80 + pkin(5)) * t82 + t84 * t100) * qJD(2) + t95;
t3 = (pkin(8) * t123 + t86) * qJD(2) + t96;
t2 = -t92 * qJD(6) - t81 * t3 + t83 * t4;
t1 = t91 * qJD(6) - t83 * t3 - t81 * t4;
t41 = [0, 0, 0, 0.2e1 * t99, 0.2e1 * (-t82 ^ 2 + t84 ^ 2) * qJD(2), 0, 0, 0, t82 * t103, t84 * t103, -0.2e1 * t21 * t84 + 0.2e1 * (t26 + 0.2e1 * t61) * t70, 0.2e1 * t22 * t84 + 0.2e1 * (-t27 + 0.2e1 * t63) * t70, 0.2e1 * (-t21 * t80 - t22 * t79) * t82 + 0.2e1 * (-t26 * t80 - t27 * t79) * t106, 0.2e1 * pkin(7) ^ 2 * t99 + 0.2e1 * t26 * t21 + 0.2e1 * t27 * t22, 0.2e1 * (t15 * t79 + t17 * t80) * t82 + 0.2e1 * (t24 * t79 + t25 * t80) * t106, -0.2e1 * t20 * t127 - 0.2e1 * t17 * t84 + 0.2e1 * (-t30 * t126 + t25 * t82) * qJD(2), -0.2e1 * t20 * t124 + 0.2e1 * t15 * t84 + 0.2e1 * (-t30 * t123 - t24 * t82) * qJD(2), 0.2e1 * t24 * t15 + 0.2e1 * t25 * t17 + 0.2e1 * t30 * t20, -0.2e1 * t10 * t124 - 0.2e1 * t9 * t84 + 0.2e1 * (t23 * t123 + t19 * t82) * qJD(2), 0.2e1 * (t79 * t9 - t8 * t80) * t82 + 0.2e1 * (-t16 * t80 + t19 * t79) * t106, 0.2e1 * t10 * t127 + 0.2e1 * t8 * t84 + 0.2e1 * (-t23 * t126 - t16 * t82) * qJD(2), -0.2e1 * t23 * t10 + 0.2e1 * t16 * t8 + 0.2e1 * t19 * t9, 0.2e1 * t32 * t14, -0.2e1 * t32 * t13 - 0.2e1 * t14 * t31, -0.2e1 * t14 * t84 + 0.2e1 * t32 * t70, 0.2e1 * t13 * t84 - 0.2e1 * t31 * t70, -0.2e1 * t99, 0.2e1 * t18 * t13 - 0.2e1 * t2 * t84 + 0.2e1 * t5 * t31 - 0.2e1 * t70 * t91, -0.2e1 * t1 * t84 + 0.2e1 * t18 * t14 + 0.2e1 * t5 * t32 - 0.2e1 * t70 * t92; 0, 0, 0, 0, 0, t106, -t70, 0, -t101, t102, t55 + (t93 * t79 - t63) * qJD(2), t84 * t107 + (t93 * t80 + t61) * qJD(2), t90, -pkin(2) * t101 + (-t26 * t79 + t27 * t80) * qJD(3) + t90 * qJ(3), -t15 * t80 + t17 * t79, t82 * t109 + t20 * t80 + t79 * t85 - t55, -t20 * t79 + (-qJD(3) * t84 + t108 * t82 + t85) * t80, t20 * t45 + (-qJ(3) * t15 - qJD(3) * t24) * t80 + (qJ(3) * t17 + qJD(3) * t25 - qJD(4) * t30) * t79, t51 * t70 - t10 * t79 + (t47 * t82 + (qJD(2) * t36 - qJD(3)) * t84) * t80, -t8 * t79 - t9 * t80 + (-t50 * t80 + t51 * t79) * t106, -t47 * t127 - t10 * t80 + t55 + (-t36 * t126 - t50 * t82) * qJD(2), -t10 * t36 + t23 * t47 + t8 * t50 + t9 * t51 + (t16 * t79 + t19 * t80) * qJD(3), t14 * t44 + t32 * t34, -t44 * t13 + t14 * t43 - t34 * t31 + t32 * t35, -t34 * t84 + t44 * t70, -t35 * t84 + t43 * t70, 0, t29 * t13 - t18 * t35 - t47 * t31 - t5 * t43 - t7 * t84 - t70 * t88, t29 * t14 + t18 * t34 - t47 * t32 + t5 * t44 - t6 * t84 - t70 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, 0.2e1 * t116, -t42, -0.2e1 * t80 * t108, 0.2e1 * t109, -0.2e1 * t108 * t45 + 0.2e1 * t116, t79 * t128, t42, t80 * t128, 0.2e1 * t36 * t47 + 0.2e1 * (t50 * t79 + t51 * t80) * qJD(3), 0.2e1 * t44 * t34, 0.2e1 * t34 * t43 + 0.2e1 * t44 * t35, 0, 0, 0, -0.2e1 * t29 * t35 + 0.2e1 * t47 * t43, 0.2e1 * t29 * t34 - 0.2e1 * t47 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t58, 0, t101, 0, -t56, -t58, t20, -t58, 0, t56, t10, 0, 0, 0, 0, 0, t13, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, 0, 0, 0, -t47, 0, 0, 0, 0, 0, -t35, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t70, 0, t17, 0, -t58, -t70, t8, 0, 0, 0, 0, 0, t38, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t56, 0, t9, 0, 0, 0, 0, 0, t37, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, t70, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t35, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t41;
