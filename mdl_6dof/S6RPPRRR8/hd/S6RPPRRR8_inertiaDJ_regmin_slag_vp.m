% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x31]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:36:31
% EndTime: 2019-03-09 02:36:33
% DurationCPUTime: 1.01s
% Computational Cost: add. (1557->160), mult. (3472->288), div. (0->0), fcn. (3478->8), ass. (0->107)
t128 = sin(qJ(4));
t129 = cos(qJ(4));
t73 = sin(pkin(10));
t74 = cos(pkin(10));
t48 = t128 * t74 + t129 * t73;
t64 = t73 * pkin(3) + qJ(2);
t80 = -t128 * t73 + t129 * t74;
t32 = t48 * pkin(4) - pkin(8) * t80 + t64;
t75 = -pkin(1) - qJ(3);
t130 = -pkin(7) + t75;
t54 = t130 * t73;
t55 = t130 * t74;
t35 = t128 * t55 + t129 * t54;
t79 = cos(qJ(5));
t33 = t79 * t35;
t77 = sin(qJ(5));
t116 = t77 * t32 + t33;
t113 = qJD(5) * t77;
t93 = qJD(4) * t128;
t94 = qJD(4) * t129;
t43 = -t73 * t94 - t74 * t93;
t117 = t79 * t43;
t83 = -t113 * t80 + t117;
t72 = t79 ^ 2;
t114 = t77 ^ 2 - t72;
t92 = qJD(5) * t114;
t56 = (t73 ^ 2 + t74 ^ 2) * qJD(3);
t108 = qJD(5) + qJD(6);
t46 = t80 ^ 2;
t135 = 2 * qJD(2);
t134 = pkin(8) + pkin(9);
t44 = -t73 * t93 + t74 * t94;
t133 = t44 * pkin(5);
t132 = t48 * pkin(5);
t112 = qJD(5) * t79;
t19 = t48 * qJD(3) + t54 * t93 - t55 * t94;
t31 = t44 * pkin(4) - t43 * pkin(8) + qJD(2);
t6 = -t32 * t112 + t35 * t113 + t79 * t19 - t77 * t31;
t119 = t77 * t43;
t84 = -t112 * t80 - t119;
t5 = t84 * pkin(9) - t6;
t78 = cos(qJ(6));
t131 = t78 * t5;
t127 = t80 * t43;
t126 = t80 * t77;
t125 = t80 * t79;
t124 = t48 * t44;
t76 = sin(qJ(6));
t120 = t76 * t77;
t50 = -t78 * t79 + t120;
t123 = t50 * t44;
t51 = t76 * t79 + t78 * t77;
t122 = t51 * t44;
t13 = -pkin(9) * t126 + t116;
t121 = t76 * t13;
t118 = t78 * t13;
t111 = qJD(6) * t76;
t110 = qJD(6) * t78;
t109 = qJ(2) * qJD(2);
t107 = -0.2e1 * pkin(4) * qJD(5);
t106 = pkin(5) * t113;
t105 = pkin(5) * t111;
t104 = pkin(5) * t110;
t102 = t77 * t112;
t101 = t48 ^ 2 + t46;
t96 = t77 * t19 + t79 * t31;
t4 = -pkin(9) * t117 + t133 + (-t33 + (pkin(9) * t80 - t32) * t77) * qJD(5) + t96;
t100 = t78 * t4 - t76 * t5;
t95 = t79 * t32 - t77 * t35;
t12 = -pkin(9) * t125 + t132 + t95;
t99 = -t12 - t132;
t98 = qJD(5) * t134;
t97 = -0.4e1 * t77 * t125;
t91 = -pkin(4) * t43 - pkin(8) * t44;
t90 = -pkin(4) * t80 - pkin(8) * t48;
t34 = t128 * t54 - t129 * t55;
t89 = t78 * t12 - t121;
t88 = t76 * t12 + t118;
t37 = t108 * t120 - t79 * t110 - t78 * t112;
t87 = t37 * t48 - t122;
t38 = t108 * t51;
t9 = -t38 * t80 - t43 * t50;
t15 = -t38 * t48 - t123;
t57 = t134 * t77;
t58 = t134 * t79;
t86 = -t78 * t57 - t76 * t58;
t85 = -t76 * t57 + t78 * t58;
t30 = t48 * t112 + t77 * t44;
t29 = t48 * t113 - t79 * t44;
t81 = -0.2e1 * t124 - 0.2e1 * t127;
t20 = t80 * qJD(3) + t35 * qJD(4);
t67 = -t79 * pkin(5) - pkin(4);
t53 = t79 * t98;
t52 = t77 * t98;
t36 = 0.2e1 * t124;
t28 = t50 * t80;
t27 = t51 * t80;
t21 = pkin(5) * t126 + t34;
t17 = -t85 * qJD(6) + t76 * t52 - t78 * t53;
t16 = -t86 * qJD(6) + t78 * t52 + t76 * t53;
t14 = -t84 * pkin(5) + t20;
t11 = -t111 * t126 + (t108 * t125 + t119) * t78 + t83 * t76;
t10 = t108 * t48 * t50 - t122;
t7 = -t116 * qJD(5) + t96;
t2 = -t88 * qJD(6) + t100;
t1 = -t89 * qJD(6) - t76 * t4 - t131;
t3 = [0, 0, 0, 0, t135, 0.2e1 * t109, t73 * t135, t74 * t135, 0.2e1 * t56, -0.2e1 * t75 * t56 + 0.2e1 * t109, 0.2e1 * t127, -0.2e1 * t43 * t48 - 0.2e1 * t44 * t80, 0, 0, 0, 0.2e1 * qJD(2) * t48 + 0.2e1 * t64 * t44, 0.2e1 * qJD(2) * t80 + 0.2e1 * t64 * t43, -0.2e1 * t46 * t102 + 0.2e1 * t72 * t127, t43 * t97 + 0.2e1 * t46 * t92, 0.2e1 * t48 * t117 - 0.2e1 * t29 * t80, -0.2e1 * t48 * t119 - 0.2e1 * t30 * t80, t36, 0.2e1 * t20 * t126 - 0.2e1 * t84 * t34 + 0.2e1 * t95 * t44 + 0.2e1 * t7 * t48, -0.2e1 * t116 * t44 + 0.2e1 * t20 * t125 + 0.2e1 * t83 * t34 + 0.2e1 * t6 * t48, -0.2e1 * t28 * t9, 0.2e1 * t28 * t11 - 0.2e1 * t9 * t27, -0.2e1 * t28 * t44 + 0.2e1 * t9 * t48, -0.2e1 * t11 * t48 - 0.2e1 * t27 * t44, t36, 0.2e1 * t21 * t11 + 0.2e1 * t14 * t27 + 0.2e1 * t2 * t48 + 0.2e1 * t89 * t44, 0.2e1 * t1 * t48 - 0.2e1 * t14 * t28 + 0.2e1 * t21 * t9 - 0.2e1 * t88 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101 * t112 + t77 * t81, t101 * t113 + t79 * t81, 0, 0, 0, 0, 0, -t80 * t11 - t43 * t27 + (t10 - t122) * t48, t43 * t28 - t80 * t9 + (-t15 + t123) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, t44, t43, 0, 0, 0, 0, 0, -t29, -t30, 0, 0, 0, 0, 0, t15, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t44, 0, -t20, t19, t77 * t117 - t80 * t92, qJD(5) * t97 - t114 * t43, t30, -t29, 0, -t20 * t79 + t91 * t77 + (t34 * t77 + t90 * t79) * qJD(5), t20 * t77 + t91 * t79 + (t34 * t79 - t90 * t77) * qJD(5), t28 * t37 + t9 * t51, -t51 * t11 + t37 * t27 + t28 * t38 - t9 * t50, -t87, t15, 0, t27 * t106 + t67 * t11 + t14 * t50 + t17 * t48 + t21 * t38 + t86 * t44, -t28 * t106 + t14 * t51 + t16 * t48 - t21 * t37 - t85 * t44 + t67 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t44, 0, 0, 0, 0, 0, t83, t84, 0, 0, 0, 0, 0, t9, t37 * t80 - t43 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t102, -0.2e1 * t92, 0, 0, 0, t77 * t107, t79 * t107, -0.2e1 * t51 * t37, 0.2e1 * t37 * t50 - 0.2e1 * t51 * t38, 0, 0, 0, 0.2e1 * t50 * t106 + 0.2e1 * t67 * t38, 0.2e1 * t51 * t106 - 0.2e1 * t67 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t84, t44, t7, t6, 0, 0, t9, -t11, t44, t78 * t133 + (t99 * t76 - t118) * qJD(6) + t100, -t131 + (-t4 - t133) * t76 + (t99 * t78 + t121) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t29, 0, 0, 0, 0, 0, t10, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t112, 0, 0, 0, 0, 0, -t38, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, -t113, 0, -pkin(8) * t112, pkin(8) * t113, 0, 0, -t37, -t38, 0, t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t105, -0.2e1 * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t11, t44, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t38, 0, t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, -t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
