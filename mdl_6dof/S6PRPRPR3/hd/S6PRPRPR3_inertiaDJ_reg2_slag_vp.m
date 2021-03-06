% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRPR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:37:39
% EndTime: 2019-03-08 19:37:44
% DurationCPUTime: 1.63s
% Computational Cost: add. (1138->156), mult. (3008->275), div. (0->0), fcn. (2908->10), ass. (0->113)
t60 = cos(qJ(4));
t113 = qJ(5) * t60;
t122 = pkin(4) + pkin(9);
t57 = sin(qJ(4));
t99 = t122 * t57;
t129 = t99 - t113;
t128 = t122 * t60;
t111 = cos(pkin(6));
t53 = sin(pkin(11));
t55 = cos(pkin(11));
t54 = sin(pkin(6));
t110 = qJD(2) * t54;
t117 = cos(qJ(2));
t80 = t117 * t110;
t58 = sin(qJ(2));
t95 = t58 * t110;
t24 = t53 * t95 - t55 * t80;
t127 = t111 * qJD(4) - t24;
t52 = t60 ^ 2;
t86 = qJD(4) * (t57 ^ 2 - t52);
t56 = sin(qJ(6));
t49 = t56 ^ 2;
t59 = cos(qJ(6));
t51 = t59 ^ 2;
t116 = t49 - t51;
t85 = t116 * qJD(6);
t107 = qJD(5) * t60;
t112 = t57 * qJ(5);
t44 = t53 * pkin(2) + pkin(8);
t119 = pkin(5) + t44;
t38 = t119 * t60;
t126 = (t112 + t128) * qJD(4) - qJD(6) * t38 - t107;
t102 = t57 * qJD(5);
t106 = qJD(6) * t56;
t45 = -t55 * pkin(2) - pkin(3);
t70 = t45 - t112;
t28 = t70 - t128;
t89 = t119 * t57;
t81 = t59 * t89;
t109 = qJD(4) * t56;
t98 = t38 * t109;
t5 = t28 * t106 - t59 * (t129 * qJD(4) - t102) - t98 - qJD(6) * t81;
t105 = qJD(6) * t59;
t6 = -t28 * t105 + (-qJD(6) * t119 + qJD(5)) * t56 * t57 + (-t56 * t99 + (t56 * qJ(5) + t59 * t119) * t60) * qJD(4);
t16 = -t56 * t28 + t81;
t17 = t59 * t28 + t56 * t89;
t74 = t16 * t56 - t17 * t59;
t2 = -t74 * qJD(6) - t5 * t56 + t6 * t59;
t67 = t54 * (t117 * t53 + t55 * t58);
t19 = -t111 * t60 + t57 * t67;
t27 = (-t117 * t55 + t53 * t58) * t54;
t12 = t19 * t56 + t27 * t59;
t25 = qJD(2) * t67;
t65 = qJD(4) * t67;
t9 = t127 * t57 + t60 * t65;
t3 = -t12 * qJD(6) - t25 * t56 + t9 * t59;
t11 = t19 * t59 - t27 * t56;
t4 = t11 * qJD(6) + t25 * t59 + t9 * t56;
t76 = t11 * t56 - t12 * t59;
t125 = t76 * qJD(6) - t3 * t59 - t4 * t56;
t15 = t27 * t25;
t10 = t127 * t60 - t57 * t65;
t20 = t111 * t57 + t60 * t67;
t7 = t20 * t10;
t124 = 0.2e1 * t19 * t9 + 0.2e1 * t15 + 0.2e1 * t7;
t123 = 0.2e1 * qJD(5);
t121 = t60 * pkin(4);
t120 = t9 * t60;
t47 = qJD(4) * t60;
t8 = t10 * t57;
t118 = t20 * t47 + t8;
t115 = t49 + t51;
t46 = qJD(4) * t57;
t108 = qJD(4) * t59;
t104 = qJD(6) * t60;
t103 = qJD(6) * t122;
t101 = qJ(5) * qJD(6);
t100 = 0.2e1 * qJD(4) * t45;
t97 = t56 * t104;
t96 = t59 * t104;
t94 = t57 * t47;
t93 = t44 * t46;
t92 = t57 * t108;
t91 = t56 * t105;
t90 = t44 * t47;
t87 = t57 * t115;
t41 = 0.2e1 * t94;
t83 = t56 * t92;
t82 = t52 * t91;
t77 = t11 * t59 + t12 * t56;
t75 = t16 * t59 + t17 * t56;
t72 = t10 * qJ(5) + t20 * qJD(5);
t69 = t10 * t56 + t20 * t105;
t68 = t10 * t59 - t20 * t106;
t30 = t119 * t46;
t66 = t129 * qJD(6) - t30;
t64 = t107 + (-t112 - t121) * qJD(4);
t63 = t10 * t60 + t9 * t57 + (t19 * t60 - t20 * t57) * qJD(4);
t62 = t63 * t44;
t48 = qJ(5) * t123;
t42 = -0.2e1 * t94;
t39 = -0.2e1 * t86;
t37 = t70 - t121;
t36 = -t56 * t46 + t96;
t35 = t57 * t105 + t56 * t47;
t34 = t92 + t97;
t33 = -t57 * t106 + t59 * t47;
t32 = qJD(4) * t87;
t31 = t102 + (-pkin(4) * t57 + t113) * qJD(4);
t22 = -t60 * t85 - t83;
t14 = t25 * t60 - t27 * t46;
t13 = t25 * t57 + t27 * t47;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t24 * t67 + 0.2e1 * t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t11 * t3 + 0.2e1 * t12 * t4 + 0.2e1 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t80, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t24, 0 (-t24 * t53 - t25 * t55) * pkin(2), 0, 0, 0, 0, 0, 0, -t14, t13, t63, t25 * t45 + t62, 0, 0, 0, 0, 0, 0, t63, t14, -t13, t25 * t37 - t27 * t31 + t62, 0, 0, 0, 0, 0, 0 (-t20 * t108 + t3) * t57 + (qJD(4) * t11 + t68) * t60 (t20 * t109 - t4) * t57 + (-qJD(4) * t12 - t69) * t60, -t76 * t46 + (t77 * qJD(6) + t3 * t56 - t4 * t59) * t60, t10 * t38 + t11 * t6 - t12 * t5 + t3 * t16 + t4 * t17 - t20 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t39, 0, t42, 0, 0, t57 * t100, t60 * t100, 0, 0, 0, 0, 0, t41, t39, t42, 0, -0.2e1 * t31 * t60 - 0.2e1 * t37 * t46, 0.2e1 * t31 * t57 - 0.2e1 * t37 * t47, -0.2e1 * t37 * t31, -0.2e1 * t49 * t94 + 0.2e1 * t82, -0.2e1 * t52 * t85 - 0.4e1 * t60 * t83, 0.2e1 * t56 * t86 - 0.2e1 * t57 * t96, -0.2e1 * t51 * t94 - 0.2e1 * t82, 0.2e1 * t57 * t97 + 0.2e1 * t59 * t86, t41, 0.2e1 * (-t38 * t108 + t6) * t57 + 0.2e1 * (qJD(4) * t16 - t38 * t106 - t30 * t59) * t60, 0.2e1 * (t5 + t98) * t57 + 0.2e1 * (-qJD(4) * t17 - t38 * t105 + t30 * t56) * t60, -0.2e1 * t74 * t46 + 0.2e1 * (t75 * qJD(6) + t5 * t59 + t56 * t6) * t60, 0.2e1 * t16 * t6 - 0.2e1 * t17 * t5 - 0.2e1 * t38 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t46 + t118 - t120, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 - t120 + (t19 * t57 + t20 * t60) * qJD(4), 0, 0, 0, 0, 0, 0, 0, 0, 0, t125 * t60 + t77 * t46 + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t75 * qJD(4) - t30) * t57 + (qJD(4) * t38 - t2) * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (0.1e1 - t115) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, -t9 * pkin(4) + t72, 0, 0, 0, 0, 0, 0, t69, t68, t125, t122 * t125 + t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, -t46, 0, -t90, t93, 0, 0, 0, -t47, t46, 0, 0, 0, t64, t90, -t93, t64 * t44, -t22, -t116 * t46 + 0.4e1 * t60 * t91, t33, t22, -t35, 0, -t126 * t59 + t66 * t56, t126 * t56 + t66 * t59, -t2, -t30 * qJ(5) + t38 * qJD(5) - t122 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t47, t31, 0, 0, 0, 0, 0, 0, t35, t33, -t32, t102 + (-t122 * t87 + t113) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, t48, -0.2e1 * t91, 0.2e1 * t85, 0, 0.2e1 * t91, 0, 0, 0.2e1 * qJD(5) * t56 + 0.2e1 * t59 * t101, 0.2e1 * qJD(5) * t59 - 0.2e1 * t56 * t101, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, t90, 0, 0, 0, 0, 0, 0, t33, -t35, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, 0, t34, t47, t6, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t36, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, 0, -t105, 0, t56 * t103, t59 * t103, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t105, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t1;
