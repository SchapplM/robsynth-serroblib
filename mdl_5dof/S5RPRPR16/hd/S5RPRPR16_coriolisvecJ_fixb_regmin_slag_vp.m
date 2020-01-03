% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR16_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:32
% EndTime: 2019-12-31 18:39:35
% DurationCPUTime: 0.94s
% Computational Cost: add. (616->175), mult. (1293->255), div. (0->0), fcn. (672->4), ass. (0->108)
t59 = cos(qJ(3));
t97 = t59 * qJD(1);
t45 = qJD(5) + t97;
t132 = qJD(5) - t45;
t57 = sin(qJ(3));
t104 = qJD(1) * t57;
t58 = cos(qJ(5));
t56 = sin(qJ(5));
t99 = t56 * qJD(3);
t131 = t58 * t104 - t99;
t61 = -pkin(1) - pkin(6);
t41 = t61 * qJD(1) + qJD(2);
t18 = (pkin(4) * qJD(1) - t41) * t59;
t95 = qJD(4) + t18;
t113 = t59 * t41;
t31 = t57 * t41;
t92 = qJD(3) * qJ(4);
t22 = -t31 - t92;
t120 = t22 * t59;
t19 = (-qJD(4) - t113) * qJD(3);
t107 = qJD(3) * pkin(3);
t74 = -qJD(4) + t107;
t20 = -t74 - t113;
t130 = ((-t20 + t113) * t57 + t120) * qJD(3) + t19 * t57;
t70 = (t45 + t97) * t57;
t129 = t57 * pkin(3) + qJ(2);
t98 = t58 * qJD(3);
t28 = t56 * t104 + t98;
t91 = qJD(1) * qJD(3);
t81 = t59 * t91;
t7 = qJD(5) * t28 - t58 * t81;
t52 = qJD(1) * qJD(2);
t127 = 0.2e1 * t52;
t60 = -pkin(3) - pkin(7);
t6 = qJD(5) * t131 + t56 * t81;
t126 = t57 * t6;
t125 = t6 * t58;
t9 = (qJD(4) - t18) * qJD(3);
t124 = t9 * t56;
t123 = t9 * t58;
t122 = pkin(4) - t61;
t119 = t131 * t45;
t118 = t131 * t57;
t117 = t28 * t45;
t116 = t45 * t59;
t115 = t45 * t60;
t114 = t58 * t59;
t62 = qJD(3) ^ 2;
t112 = t62 * t57;
t111 = t62 * t59;
t93 = qJ(4) * qJD(1);
t30 = pkin(3) * t97 + t57 * t93;
t110 = pkin(3) * t104 + qJD(1) * qJ(2);
t54 = t57 ^ 2;
t55 = t59 ^ 2;
t109 = t54 - t55;
t63 = qJD(1) ^ 2;
t108 = t62 + t63;
t106 = t59 * qJ(4);
t105 = t63 * qJ(2);
t103 = qJD(3) * t57;
t102 = qJD(3) * t59;
t101 = qJD(5) * t58;
t17 = -pkin(4) * t104 + t31;
t12 = t17 + t92;
t100 = t12 * qJD(5);
t96 = t59 * qJD(4);
t94 = qJ(2) * qJD(3);
t90 = t56 * t116;
t89 = t59 * t63 * t57;
t80 = t57 * t91;
t88 = pkin(3) * t81 + qJ(4) * t80 + t52;
t87 = 0.2e1 * qJD(1);
t86 = qJD(5) * t56 * t45;
t85 = t57 * t101;
t83 = pkin(3) * t102 + t57 * t92 + qJD(2);
t29 = t41 * t103;
t13 = -pkin(4) * t80 + t29;
t68 = (qJD(3) * pkin(7) - qJD(4)) * t59;
t3 = qJD(1) * t68 + t88;
t82 = t58 * t13 - t56 * t3;
t78 = qJD(3) * t122;
t21 = -t59 * t93 + t110;
t32 = -t106 + t129;
t76 = qJD(1) * t32 + t21;
t72 = t57 * pkin(7) - t106;
t11 = t72 * qJD(1) + t110;
t4 = qJD(3) * t60 + t95;
t2 = t58 * t11 + t56 * t4;
t73 = t56 * t11 - t58 * t4;
t23 = t72 + t129;
t34 = t122 * t59;
t71 = t58 * t23 + t56 * t34;
t69 = -qJD(1) * t54 + t116;
t67 = t45 * (qJD(5) * t59 + qJD(1));
t15 = t83 - t96;
t5 = -qJD(1) * t96 + t88;
t65 = -qJD(1) * t15 + t61 * t62 - t5;
t37 = t56 * t80;
t36 = t108 * t59;
t35 = t108 * t57;
t33 = t122 * t57;
t25 = t59 * t78;
t24 = t57 * t78;
t16 = pkin(7) * t97 + t30;
t14 = t21 * t97;
t8 = t68 + t83;
t1 = [0, 0, 0, 0, t127, qJ(2) * t127, -0.2e1 * t59 * t80, 0.2e1 * t109 * t91, -t112, -t111, 0, -t61 * t112 + (qJD(2) * t57 + t59 * t94) * t87, -t61 * t111 + (qJD(2) * t59 - t57 * t94) * t87, t130, -t76 * t102 + t65 * t57, t76 * t103 + t65 * t59, -t130 * t61 + t21 * t15 + t5 * t32, t56 * t126 + (t59 * t99 + t85) * t28, (t131 * t56 + t28 * t58) * t102 + (-t56 * t7 + t125 + (t131 * t58 - t28 * t56) * qJD(5)) * t57, t45 * t85 + t6 * t59 + (-t28 * t57 + t69 * t56) * qJD(3), -t57 * t86 - t7 * t59 + (t69 * t58 - t118) * qJD(3), -qJD(3) * t70, (-t58 * t24 - t56 * t8) * t45 + t25 * t131 - t33 * t7 + (-t12 * t98 + t82) * t59 + (-t2 * t59 - t71 * t45) * qJD(5) + (t56 * t100 - t123 + (-(-t56 * t23 + t58 * t34) * qJD(1) + t73) * qJD(3)) * t57, -t25 * t28 - t33 * t6 + (-(qJD(5) * t34 + t8) * t45 - (qJD(5) * t4 + t3) * t59) * t58 + (-(-qJD(5) * t23 - t24) * t45 + (t12 * qJD(3) + qJD(5) * t11 - t13) * t59) * t56 + (t58 * t100 + t124 + (t71 * qJD(1) + t2) * qJD(3)) * t57; 0, 0, 0, 0, -t63, -t105, 0, 0, 0, 0, 0, -t35, -t36, 0, t35, t36, -t21 * qJD(1) - t130, 0, 0, 0, 0, 0, t57 * t7 + t56 * t67 + (-t131 * t59 + t58 * t70) * qJD(3), t126 + t58 * t67 + (t28 * t59 - t56 * t70) * qJD(3); 0, 0, 0, 0, 0, 0, t89, -t109 * t63, 0, 0, 0, -t59 * t105, t57 * t105, ((-t22 - t92) * t59 + (t20 + t74) * t57) * qJD(1), t30 * t104 + t14, 0.2e1 * qJD(3) * qJD(4) + (-t21 * t57 + t30 * t59) * qJD(1), -t19 * qJ(4) - t22 * qJD(4) - t21 * t30 + (t120 + (-t20 - t107) * t57) * t41, -t117 * t56 + t125, (-t7 - t117) * t58 + (-t6 - t119) * t56, -t86 + (-t90 + (t28 - t98) * t57) * qJD(1), -t45 * t101 + t37 + (-t45 * t114 + t118) * qJD(1), t45 * t104, qJ(4) * t7 + t124 - (-t56 * t16 + t58 * t17) * t45 - t95 * t131 + (-t56 * t115 + t12 * t58) * qJD(5) + (t12 * t114 + (-t60 * t98 - t73) * t57) * qJD(1), qJ(4) * t6 + t123 + (t58 * t16 + t56 * t17) * t45 + t95 * t28 + (-t58 * t115 - t12 * t56) * qJD(5) + (-t2 * t57 + (t60 * t103 - t12 * t59) * t56) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t55 * t63 - t62, t22 * qJD(3) + t14 + t29, 0, 0, 0, 0, 0, -t86 + qJD(3) * t131 + (-t57 * t98 - t90) * qJD(1), -t45 ^ 2 * t58 - qJD(3) * t28 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * t131, -t131 ^ 2 + t28 ^ 2, t6 - t119, t117 - t7, -t80, -t12 * t28 - t132 * t2 + t82, -t12 * t131 - t56 * t13 + t132 * t73 - t58 * t3;];
tauc_reg = t1;
