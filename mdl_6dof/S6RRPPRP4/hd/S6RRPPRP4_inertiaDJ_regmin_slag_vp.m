% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:39:41
% EndTime: 2019-03-09 08:39:46
% DurationCPUTime: 1.76s
% Computational Cost: add. (1674->232), mult. (3971->406), div. (0->0), fcn. (3334->6), ass. (0->102)
t86 = cos(qJ(2));
t111 = t86 * qJD(2);
t82 = cos(pkin(9));
t84 = sin(qJ(2));
t128 = t82 * t84;
t64 = t82 * t111;
t122 = qJ(4) * t64 + qJD(4) * t128;
t81 = sin(pkin(9));
t134 = pkin(7) + (pkin(3) + pkin(4)) * t81;
t23 = -t134 * t111 + t122;
t79 = t81 ^ 2;
t138 = (t82 ^ 2 + t79) * qJD(3);
t120 = t84 * qJ(3);
t93 = -t86 * pkin(2) - t120;
t54 = -pkin(1) + t93;
t130 = t81 * t86;
t69 = pkin(7) * t130;
t78 = t86 * pkin(3);
t24 = t86 * pkin(4) + t69 + t78 + (-pkin(8) * t84 - t54) * t82;
t131 = t81 * t84;
t127 = t82 * t86;
t70 = pkin(7) * t127;
t37 = t81 * t54 + t70;
t34 = -t86 * qJ(4) + t37;
t27 = pkin(8) * t131 + t34;
t83 = sin(qJ(5));
t85 = cos(qJ(5));
t137 = t83 * t24 + t85 * t27;
t112 = t84 * qJD(2);
t136 = qJ(4) * t112 - t86 * qJD(4);
t53 = -t82 * pkin(3) - t81 * qJ(4) - pkin(2);
t135 = (-t53 * t86 + t120) * qJD(2);
t118 = qJD(5) * t85;
t119 = qJD(5) * t83;
t44 = t81 * t118 - t82 * t119;
t52 = t85 * t81 - t83 * t82;
t100 = -pkin(7) * t81 - pkin(3);
t42 = -t84 * qJD(3) + (pkin(2) * t84 - qJ(3) * t86) * qJD(2);
t129 = t82 * t42;
t13 = -t129 + (-pkin(8) * t127 + (-pkin(4) + t100) * t84) * qJD(2);
t38 = t81 * t42;
t14 = t38 + (-pkin(7) * t128 + pkin(8) * t130) * qJD(2) + t136;
t4 = -qJD(5) * t137 + t85 * t13 - t83 * t14;
t133 = 2 * qJD(6);
t124 = -pkin(8) + qJ(3);
t121 = qJ(3) * t138;
t117 = qJD(5) * t86;
t116 = t79 * qJD(4);
t115 = t81 * qJD(3);
t114 = t81 * qJD(4);
t113 = t82 * qJD(3);
t109 = t86 * qJD(6);
t107 = -0.2e1 * pkin(1) * qJD(2);
t106 = pkin(5) * t112;
t105 = pkin(7) * t112;
t104 = pkin(7) * t111;
t61 = t81 * t111;
t101 = t84 * t111;
t99 = pkin(3) * t81 + pkin(7);
t98 = qJ(6) * t112;
t97 = t124 * t81;
t36 = t82 * t54 - t69;
t45 = t82 * pkin(4) - t53;
t95 = t85 * t97;
t92 = t85 * t24 - t83 * t27;
t29 = t81 * t105 + t129;
t30 = -t82 * t105 + t38;
t90 = -t29 * t81 + t30 * t82;
t51 = t83 * t81 + t85 * t82;
t57 = t124 * t82;
t15 = -qJD(5) * t95 - t85 * t113 - t83 * t115 + t57 * t119;
t32 = t85 * t57 + t83 * t97;
t88 = -t32 * t112 - t15 * t86;
t16 = t32 * qJD(5) + t83 * t113 - t85 * t115;
t31 = t83 * t57 - t95;
t87 = t31 * t112 - t16 * t86;
t3 = -t24 * t118 + t27 * t119 - t83 * t13 - t85 * t14;
t43 = t51 * qJD(5);
t68 = qJ(4) * t128;
t33 = -t134 * t84 + t68;
t60 = t86 * t115;
t50 = 0.2e1 * t138;
t47 = -t83 * t112 + t85 * t117;
t46 = -t85 * t112 - t83 * t117;
t41 = t51 * t84;
t40 = t52 * t84;
t39 = t99 * t84 - t68;
t35 = -t36 + t78;
t28 = t99 * t111 - t122;
t26 = t100 * t112 - t129;
t22 = t30 + t136;
t20 = t51 * t111 + t44 * t84;
t19 = t84 * t43 - t85 * t61 + t83 * t64;
t17 = t51 * pkin(5) - t52 * qJ(6) + t45;
t9 = t44 * pkin(5) + t43 * qJ(6) - t52 * qJD(6) + t114;
t8 = -t40 * pkin(5) - t41 * qJ(6) + t33;
t7 = -t86 * pkin(5) - t92;
t6 = t86 * qJ(6) + t137;
t5 = -t19 * pkin(5) + t20 * qJ(6) + t41 * qJD(6) - t23;
t2 = t106 - t4;
t1 = -t3 - t98 + t109;
t10 = [0, 0, 0, 0.2e1 * t101, 0.2e1 * (-t84 ^ 2 + t86 ^ 2) * qJD(2), 0, 0, 0, t84 * t107, t86 * t107, -0.2e1 * t29 * t86 + 0.2e1 * (t36 + 0.2e1 * t69) * t112, 0.2e1 * t30 * t86 + 0.2e1 * (-t37 + 0.2e1 * t70) * t112, 0.2e1 * (-t29 * t82 - t30 * t81) * t84 + 0.2e1 * (-t36 * t82 - t37 * t81) * t111, 0.2e1 * pkin(7) ^ 2 * t101 + 0.2e1 * t36 * t29 + 0.2e1 * t37 * t30, 0.2e1 * t28 * t131 + 0.2e1 * t26 * t86 + 0.2e1 * (t39 * t130 - t35 * t84) * qJD(2), 0.2e1 * (-t22 * t81 + t26 * t82) * t84 + 0.2e1 * (-t34 * t81 + t35 * t82) * t111, -0.2e1 * t28 * t128 - 0.2e1 * t22 * t86 + 0.2e1 * (-t39 * t127 + t34 * t84) * qJD(2), 0.2e1 * t34 * t22 + 0.2e1 * t35 * t26 + 0.2e1 * t39 * t28, 0.2e1 * t41 * t20, -0.2e1 * t41 * t19 + 0.2e1 * t20 * t40, -0.2e1 * t41 * t112 + 0.2e1 * t20 * t86, -0.2e1 * t40 * t112 - 0.2e1 * t19 * t86, -0.2e1 * t101, -0.2e1 * t92 * t112 + 0.2e1 * t33 * t19 - 0.2e1 * t23 * t40 + 0.2e1 * t4 * t86, 0.2e1 * t112 * t137 + 0.2e1 * t33 * t20 + 0.2e1 * t23 * t41 + 0.2e1 * t3 * t86, 0.2e1 * t7 * t112 + 0.2e1 * t8 * t19 - 0.2e1 * t2 * t86 + 0.2e1 * t5 * t40, 0.2e1 * t1 * t40 - 0.2e1 * t6 * t19 + 0.2e1 * t2 * t41 + 0.2e1 * t7 * t20, 0.2e1 * t1 * t86 - 0.2e1 * t112 * t6 - 0.2e1 * t8 * t20 + 0.2e1 * t5 * t41, 0.2e1 * t6 * t1 + 0.2e1 * t7 * t2 - 0.2e1 * t8 * t5; 0, 0, 0, 0, 0, t111, -t112, 0, -t104, t105, t60 + (t93 * t81 - t70) * qJD(2), t86 * t113 + (t93 * t82 + t69) * qJD(2), t90, -pkin(2) * t104 + (-t36 * t81 + t37 * t82) * qJD(3) + t90 * qJ(3), -t84 * t116 - t81 * t135 - t28 * t82 + t60, t22 * t82 + t26 * t81, -t28 * t81 + (-qJD(3) * t86 + t84 * t114 + t135) * t82, t28 * t53 + (qJ(3) * t22 + qJD(3) * t34) * t82 + (qJ(3) * t26 + qJD(3) * t35 - qJD(4) * t39) * t81, t20 * t52 - t41 * t43, -t52 * t19 - t20 * t51 - t43 * t40 - t41 * t44, -t52 * t112 - t43 * t86, t51 * t112 - t44 * t86, 0, -t40 * t114 + t45 * t19 + t23 * t51 + t33 * t44 + t87, t41 * t114 + t45 * t20 + t23 * t52 - t33 * t43 - t88, t17 * t19 - t9 * t40 + t8 * t44 - t5 * t51 + t87, -t1 * t51 - t15 * t40 + t16 * t41 - t32 * t19 + t2 * t52 + t31 * t20 - t7 * t43 - t6 * t44, -t17 * t20 - t9 * t41 + t8 * t43 + t5 * t52 + t88, t1 * t32 - t6 * t15 + t7 * t16 - t5 * t17 + t2 * t31 + t8 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0.2e1 * t121, 0.2e1 * t82 * t114, t50, 0.2e1 * t116, -0.2e1 * t53 * t114 + 0.2e1 * t121, -0.2e1 * t52 * t43, 0.2e1 * t43 * t51 - 0.2e1 * t52 * t44, 0, 0, 0, 0.2e1 * t51 * t114 + 0.2e1 * t45 * t44, 0.2e1 * t52 * t114 - 0.2e1 * t45 * t43, 0.2e1 * t17 * t44 + 0.2e1 * t9 * t51, 0.2e1 * t15 * t51 + 0.2e1 * t16 * t52 - 0.2e1 * t31 * t43 - 0.2e1 * t32 * t44, 0.2e1 * t17 * t43 - 0.2e1 * t9 * t52, -0.2e1 * t32 * t15 + 0.2e1 * t31 * t16 + 0.2e1 * t17 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t64, 0, t104, t61, 0, -t64, t28, 0, 0, 0, 0, 0, -t19, -t20, -t19, 0, t20, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, 0, 0, 0, 0, 0, -t44, t43, -t44, 0, -t43, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, t64, 0, t26, 0, 0, 0, 0, 0, t46, -t47, t46, -t83 * t19 - t85 * t20 + (t40 * t85 + t41 * t83) * qJD(5), t47, t1 * t83 - t2 * t85 + (t6 * t85 + t7 * t83) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, 0, 0, 0, 0, 0, 0, 0, 0, t85 * t43 - t83 * t44 + (-t51 * t85 + t52 * t83) * qJD(5), 0, -t15 * t83 - t16 * t85 + (t31 * t83 + t32 * t85) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, -t112, t4, t3, t4 - 0.2e1 * t106, -t20 * pkin(5) - t19 * qJ(6) + t40 * qJD(6), -t3 - 0.2e1 * t98 + 0.2e1 * t109, -t2 * pkin(5) + t1 * qJ(6) + t6 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t44, 0, -t16, t15, -t16, pkin(5) * t43 - t44 * qJ(6) - t51 * qJD(6), -t15, -t16 * pkin(5) - t15 * qJ(6) + t32 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, -t118, -t119, 0, t118, t83 * qJD(6) + (-pkin(5) * t83 + qJ(6) * t85) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, qJ(6) * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t20, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
