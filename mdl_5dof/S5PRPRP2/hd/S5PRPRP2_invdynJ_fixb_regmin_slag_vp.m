% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:07
% EndTime: 2019-12-05 15:31:11
% DurationCPUTime: 0.80s
% Computational Cost: add. (692->182), mult. (1423->247), div. (0->0), fcn. (925->6), ass. (0->105)
t54 = sin(pkin(8));
t55 = cos(pkin(8));
t31 = -pkin(3) * t55 - pkin(6) * t54 - pkin(2);
t24 = t31 * qJD(2) + qJD(3);
t92 = qJ(3) * qJD(2);
t29 = qJD(1) * t54 + t55 * t92;
t57 = sin(qJ(4));
t58 = cos(qJ(4));
t101 = qJD(2) * t54;
t80 = qJ(5) * t101;
t7 = (t24 - t80) * t57 + t29 * t58;
t129 = qJD(4) * t7;
t51 = pkin(7) + qJ(2);
t47 = sin(t51);
t48 = cos(t51);
t126 = g(1) * t47 - g(2) * t48;
t94 = qJDD(2) * pkin(2);
t127 = t94 + t126;
t86 = qJD(2) * qJD(3);
t87 = qJ(3) * qJDD(2);
t63 = t86 + t87;
t111 = t55 * t57;
t17 = t47 * t111 + t48 * t58;
t19 = -t48 * t111 + t47 * t58;
t125 = -g(1) * t19 + g(2) * t17;
t122 = pkin(4) * t57;
t85 = qJD(2) * qJD(4);
t76 = t58 * t85;
t89 = t54 * qJDD(2);
t124 = t54 * pkin(4) * t76 + t89 * t122 + qJDD(5);
t100 = qJD(2) * t55;
t38 = -qJD(4) + t100;
t74 = t58 * t24 - t29 * t57;
t6 = -t58 * t80 + t74;
t3 = -pkin(4) * t38 + t6;
t123 = -t6 + t3;
t45 = t55 * qJDD(1);
t21 = t63 * t54 - t45;
t118 = t21 * t54;
t88 = t55 * qJDD(2);
t36 = -qJDD(4) + t88;
t116 = t36 * t55;
t115 = (pkin(4) * t58 + pkin(3)) * t55;
t114 = t47 * t57;
t113 = t48 * t57;
t49 = t54 ^ 2;
t59 = qJD(2) ^ 2;
t112 = t49 * t59;
t110 = t55 * t58;
t96 = qJD(4) * t58;
t98 = qJD(3) * t55;
t109 = t31 * t96 + t58 * t98;
t37 = qJ(3) * t110;
t108 = t57 * t31 + t37;
t107 = t48 * pkin(2) + t47 * qJ(3);
t106 = t55 ^ 2 + t49;
t52 = t57 ^ 2;
t53 = t58 ^ 2;
t105 = -t52 - t53;
t104 = t52 - t53;
t103 = qJ(3) * t57;
t102 = qJ(5) * t54;
t99 = qJD(2) * t57;
t97 = qJD(4) * t57;
t95 = qJD(5) * t54;
t93 = -qJD(4) - t38;
t91 = qJDD(2) * t57;
t90 = qJDD(2) * t58;
t84 = qJD(2) * qJD(5);
t83 = t55 * t103;
t82 = t58 * t102;
t81 = t38 * t97;
t79 = -pkin(2) * t47 + t48 * qJ(3);
t77 = t57 * t85;
t75 = qJ(3) + t122;
t73 = t36 - t88;
t72 = t36 + t88;
t71 = qJD(2) * t93;
t70 = g(1) * t48 + g(2) * t47;
t68 = -t3 * t58 - t57 * t7;
t67 = t3 * t57 - t58 * t7;
t22 = qJDD(1) * t54 + t63 * t55;
t66 = t22 * t55 + t118;
t46 = t55 * qJD(1);
t28 = t54 * t92 - t46;
t65 = t28 * t54 + t29 * t55;
t23 = t31 * qJDD(2) + qJDD(3);
t62 = -t58 * t22 - t57 * t23 - t24 * t96 + t29 * t97;
t44 = qJDD(3) - t94;
t61 = -t44 + t127;
t60 = -t38 ^ 2 - t112;
t30 = t54 * t55 * t77;
t27 = t58 * t31;
t20 = t48 * t110 + t114;
t18 = -t47 * t110 + t113;
t15 = t58 * t23;
t12 = t75 * t101 + qJD(5) - t46;
t10 = -t57 * t102 + t108;
t9 = -t82 + t27 + (-pkin(4) - t103) * t55;
t8 = t21 + t124;
t5 = -t57 * t98 - t58 * t95 + (-t37 + (-t31 + t102) * t57) * qJD(4);
t4 = -t57 * t95 + (-t82 - t83) * qJD(4) + t109;
t2 = (-t57 * t84 + (-t76 - t91) * qJ(5)) * t54 - t62;
t1 = -pkin(4) * t36 - t57 * t22 + t15 + (-qJ(5) * qJDD(2) - t84) * t58 * t54 - t129;
t11 = [qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, -t21 * t55 + t22 * t54 - g(3), 0, 0, 0, 0, 0, (t73 * t57 + (t38 - t100) * t96) * t54, t30 + (t73 * t58 - t81) * t54, 0, -t55 * t8 - g(3) + (t68 * qJD(4) - t1 * t57 + t2 * t58) * t54; 0, qJDD(2), t126, t70, t61 * t55, -t61 * t54, t63 * t106 + t66 - t70, -t44 * pkin(2) - g(1) * t79 - g(2) * t107 + t66 * qJ(3) + t65 * qJD(3), (qJDD(2) * t53 - 0.2e1 * t57 * t76) * t49, 0.2e1 * (t104 * t85 - t57 * t90) * t49, t30 + (-t72 * t58 + t81) * t54, (t72 * t57 + (t38 + t100) * t96) * t54, t116, -g(1) * t18 - g(2) * t20 - t15 * t55 - t27 * t36 + ((qJD(2) * t49 + t38 * t55) * qJ(3) + t65) * t96 + (-(-qJD(4) * t31 - t98) * t38 - (-qJD(4) * t24 - t22) * t55 + t49 * t86 + t118 + (t49 * qJDD(2) + t116) * qJ(3)) * t57, (-qJD(4) * t83 + t109) * t38 + t108 * t36 - t62 * t55 - g(1) * t17 - g(2) * t19 + (t21 * t58 - t28 * t97) * t54 + (t58 * t86 + (-t77 + t90) * qJ(3)) * t49, ((-t129 - qJDD(2) * t9 - t1 + (-qJD(4) * t10 - t5) * qJD(2)) * t58 + (qJD(4) * t3 - qJDD(2) * t10 - t2 + (qJD(4) * t9 - t4) * qJD(2)) * t57 + t126) * t54, t2 * t10 + t7 * t4 + t1 * t9 + t3 * t5 - g(1) * (pkin(4) * t113 - t47 * t115 + t79) - g(2) * (pkin(4) * t114 + t48 * t115 + t107) + (t8 * t75 + t12 * (pkin(4) * t96 + qJD(3)) - t126 * (-qJ(5) - pkin(6))) * t54; 0, 0, 0, 0, -t88, t89, -t106 * t59, -t65 * qJD(2) + qJDD(3) - t127, 0, 0, 0, 0, 0, -t36 * t58 + t60 * t57, t36 * t57 + t60 * t58, t105 * t89, t1 * t58 + t2 * t57 - t67 * qJD(4) + (-t12 * t54 + t67 * t55) * qJD(2) - t126; 0, 0, 0, 0, 0, 0, 0, 0, t58 * t57 * t112, -t104 * t112, (t57 * t71 + t90) * t54, (t58 * t71 - t91) * t54, -t36, t15 + (-t28 * t101 + t93 * t29) * t58 + (g(3) * t54 + t93 * t24 - t22) * t57 + t125, -t74 * t38 + g(1) * t20 - g(2) * t18 + (g(3) * t58 + t28 * t99) * t54 + t62, (-pkin(4) * t90 + (pkin(4) * qJD(4) - t123) * t99) * t54, t123 * t7 + (t1 + (-qJD(2) * t12 * t58 + g(3) * t57) * t54 + t125) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105 * t112, g(3) * t55 - t45 + (t87 + (qJD(3) - t68) * qJD(2) - t70) * t54 + t124;];
tau_reg = t11;
