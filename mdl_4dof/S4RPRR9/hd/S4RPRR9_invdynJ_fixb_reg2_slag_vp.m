% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:31
% EndTime: 2019-12-31 16:56:34
% DurationCPUTime: 1.10s
% Computational Cost: add. (1144->224), mult. (2252->308), div. (0->0), fcn. (1273->6), ass. (0->121)
t53 = sin(qJ(1));
t56 = cos(qJ(1));
t143 = g(1) * t53 - g(2) * t56;
t59 = qJD(1) ^ 2;
t79 = -t59 * qJ(2) - t143;
t52 = sin(qJ(3));
t35 = t52 * qJD(1) + qJD(4);
t105 = qJD(3) * t52;
t57 = -pkin(1) - pkin(5);
t32 = t57 * qJDD(1) + qJDD(2);
t33 = t57 * qJD(1) + qJD(2);
t55 = cos(qJ(3));
t13 = -qJDD(3) * pkin(3) + t33 * t105 - t55 * t32;
t140 = g(3) * t52;
t63 = -t143 * t55 - t13 + t140;
t142 = -pkin(6) * qJD(4) * t35 + t63;
t49 = t52 ^ 2;
t50 = t55 ^ 2;
t112 = t49 + t50;
t86 = t112 * t32;
t108 = qJD(1) * t55;
t54 = cos(qJ(4));
t51 = sin(qJ(4));
t99 = t51 * qJD(3);
t28 = t54 * t108 + t99;
t103 = qJD(4) * t28;
t90 = t52 * t99;
t95 = t55 * qJDD(1);
t8 = -qJD(1) * t90 - t54 * qJDD(3) + t51 * t95 + t103;
t77 = t52 * pkin(3) - t55 * pkin(6);
t100 = qJD(4) * t55;
t98 = t54 * qJD(3);
t67 = t51 * t100 + t52 * t98;
t7 = t67 * qJD(1) - qJD(4) * t98 - t51 * qJDD(3) - t54 * t95;
t141 = 0.2e1 * qJ(2);
t139 = g(3) * t55;
t30 = qJ(2) + t77;
t17 = t30 * qJD(1);
t124 = t52 * t33;
t22 = qJD(3) * pkin(6) + t124;
t5 = t54 * t17 - t51 * t22;
t138 = t5 * t35;
t135 = t55 * t7;
t134 = t55 * t8;
t6 = t51 * t17 + t54 * t22;
t133 = t6 * t35;
t132 = t7 * t51;
t131 = t8 * t54;
t26 = t51 * t108 - t98;
t130 = t26 * t35;
t129 = t28 * t26;
t128 = t28 * t35;
t127 = t35 * t51;
t126 = t35 * t54;
t94 = qJD(1) * qJD(3);
t85 = t55 * t94;
t96 = t52 * qJDD(1);
t24 = qJDD(4) + t85 + t96;
t125 = t51 * t24;
t123 = t52 * t57;
t122 = t53 * t51;
t121 = t53 * t54;
t120 = t54 * t24;
t119 = t55 * t33;
t118 = t56 * t51;
t117 = t56 * t54;
t88 = 0.2e1 * qJD(1) * qJD(2);
t116 = (qJDD(1) * qJ(2) + t88) * qJ(2);
t115 = t56 * pkin(1) + t53 * qJ(2);
t113 = t49 - t50;
t58 = qJD(3) ^ 2;
t111 = -t58 - t59;
t109 = pkin(1) * qJDD(1);
t107 = qJD(3) * t26;
t106 = qJD(3) * t28;
t104 = qJD(3) * t55;
t102 = qJD(4) * t51;
t101 = qJD(4) * t54;
t97 = qJDD(3) * t52;
t93 = t55 * t59 * t52;
t89 = t57 * t104;
t84 = t112 * qJDD(1);
t83 = qJD(4) * t52 + qJD(1);
t82 = g(2) * (t56 * pkin(5) + t115);
t81 = qJDD(2) - t109;
t80 = t52 * t85;
t78 = pkin(3) * t55 + pkin(6) * t52;
t76 = g(1) * t56 + g(2) * t53;
t74 = t5 * t54 + t51 * t6;
t73 = t5 * t51 - t54 * t6;
t16 = t54 * t123 + t51 * t30;
t15 = -t51 * t123 + t54 * t30;
t14 = qJDD(3) * pkin(6) + t33 * t104 + t52 * t32;
t71 = -qJD(4) * t17 + t139 - t14;
t70 = t35 * t101 + t125;
t69 = -t35 * t102 + t120;
t68 = qJDD(3) * t57 + t94 * t141;
t66 = -t32 - t79;
t25 = t78 * qJD(3) + qJD(2);
t23 = -qJD(3) * pkin(3) - t119;
t65 = -pkin(6) * t24 + t35 * t23;
t64 = -t143 * t52 - t139;
t62 = qJDD(1) * t141 - t76 + t88;
t10 = t25 * qJD(1) + t30 * qJDD(1);
t1 = qJD(4) * t5 + t51 * t10 + t54 * t14;
t9 = t54 * t10;
t2 = -qJD(4) * t6 - t51 * t14 + t9;
t61 = -t74 * qJD(4) + t1 * t54 - t2 * t51;
t60 = -t57 * t58 + t62;
t43 = t56 * qJ(2);
t40 = qJDD(3) * t55;
t29 = t78 * qJD(1);
t21 = t52 * t117 - t122;
t20 = t52 * t118 + t121;
t19 = t52 * t121 + t118;
t18 = -t52 * t122 + t117;
t12 = t54 * t119 + t51 * t29;
t11 = -t51 * t119 + t54 * t29;
t4 = -t16 * qJD(4) + t54 * t25 - t51 * t89;
t3 = t15 * qJD(4) + t51 * t25 + t54 * t89;
t27 = [0, 0, 0, 0, 0, qJDD(1), t143, t76, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - 0.2e1 * t109 - t143, t62, -t81 * pkin(1) - g(1) * (-t53 * pkin(1) + t43) - g(2) * t115 + t116, t50 * qJDD(1) - 0.2e1 * t80, 0.2e1 * t113 * t94 - 0.2e1 * t52 * t95, -t58 * t52 + t40, t49 * qJDD(1) + 0.2e1 * t80, -t58 * t55 - t97, 0, t60 * t52 + t68 * t55, -t68 * t52 + t60 * t55, -t57 * t84 + t143 - t86, -g(1) * (t57 * t53 + t43) - t82 + t57 * t86 + t116, -t54 * t135 - t28 * t67, (t26 * t54 + t28 * t51) * t105 + (t132 - t131 + (t26 * t51 - t28 * t54) * qJD(4)) * t55, (-t35 * t98 - t7) * t52 + (t69 + t106) * t55, t51 * t134 + (t100 * t54 - t90) * t26, (t35 * t99 - t8) * t52 + (-t70 - t107) * t55, t104 * t35 + t24 * t52, -g(1) * t21 - g(2) * t19 + t15 * t24 + t4 * t35 + (t2 + (-t23 * t51 + t26 * t57) * qJD(3)) * t52 + (qJD(3) * t5 + t101 * t23 + t13 * t51 - t57 * t8) * t55, g(1) * t20 - g(2) * t18 - t16 * t24 - t3 * t35 + (-t1 + (-t23 * t54 + t28 * t57) * qJD(3)) * t52 + (-qJD(3) * t6 - t102 * t23 + t13 * t54 + t57 * t7) * t55, t15 * t7 - t16 * t8 - t3 * t26 - t4 * t28 + t74 * t105 + (qJD(4) * t73 - t1 * t51 - t2 * t54 + t76) * t55, t1 * t16 + t6 * t3 + t2 * t15 + t5 * t4 - g(1) * (t77 * t56 + t43) - t82 + (t105 * t23 - t13 * t55) * t57 + (-g(1) * t57 - g(2) * t77) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t59, t79 + t81, 0, 0, 0, 0, 0, 0, t111 * t52 + t40, t111 * t55 - t97, -t84, t86 + t79, 0, 0, 0, 0, 0, 0, -t134 + (t107 - t125) * t52 + (-t54 * t83 - t55 * t99) * t35, t135 + (t106 - t120) * t52 + (t51 * t83 - t55 * t98) * t35, (-t26 * t104 + t28 * t83 - t8 * t52) * t54 + (t104 * t28 + t83 * t26 - t7 * t52) * t51, -t74 * qJD(1) + (-qJD(3) * t73 - t13) * t55 + (qJD(3) * t23 + t61) * t52 - t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, -t113 * t59, t95, -t93, -t96, qJDD(3), -t66 * t55 + t140, t66 * t52 + t139, 0, 0, t126 * t28 - t132, (-t7 - t130) * t54 + (-t8 - t128) * t51, (t126 * t52 - t28 * t55) * qJD(1) + t70, t26 * t127 - t131, (-t127 * t52 + t26 * t55) * qJD(1) + t69, -t35 * t108, -pkin(3) * t8 - t5 * t108 - t11 * t35 - t26 * t124 + t142 * t54 + t65 * t51, pkin(3) * t7 + t6 * t108 + t12 * t35 - t28 * t124 - t142 * t51 + t65 * t54, t11 * t28 + t12 * t26 + (t1 - t138 + (-t8 + t103) * pkin(6)) * t54 + (-t2 - t133 + (qJD(4) * t26 - t7) * pkin(6)) * t51 + t64, -t23 * t124 - t5 * t11 - t6 * t12 + t63 * pkin(3) + (t61 + t64) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, -t26 ^ 2 + t28 ^ 2, t130 - t7, -t129, t128 - t8, t24, -g(1) * t18 - g(2) * t20 - t101 * t22 - t23 * t28 + t51 * t71 + t133 + t9, g(1) * t19 - g(2) * t21 + t23 * t26 + t138 + (qJD(4) * t22 - t10) * t51 + t71 * t54, 0, 0;];
tau_reg = t27;
