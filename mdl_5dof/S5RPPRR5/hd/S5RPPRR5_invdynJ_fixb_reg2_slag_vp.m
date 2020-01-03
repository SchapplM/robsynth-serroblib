% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:43
% EndTime: 2019-12-31 17:56:44
% DurationCPUTime: 0.68s
% Computational Cost: add. (1514->162), mult. (2214->192), div. (0->0), fcn. (1159->10), ass. (0->103)
t106 = qJD(1) - qJD(4);
t133 = t106 ^ 2;
t66 = cos(pkin(8));
t51 = -t66 * pkin(1) - pkin(2);
t43 = -pkin(3) + t51;
t32 = t43 * qJD(1) + qJD(3);
t65 = sin(pkin(8));
t46 = t65 * pkin(1) + qJ(3);
t38 = t46 * qJD(1);
t68 = sin(qJ(4));
t71 = cos(qJ(4));
t11 = t71 * t32 - t68 * t38;
t126 = t106 * pkin(4);
t9 = -t11 + t126;
t132 = t106 * t9;
t125 = t106 * t11;
t12 = t68 * t32 + t71 * t38;
t124 = t12 * t106;
t110 = qJD(5) * t106;
t59 = qJDD(1) - qJDD(4);
t119 = t71 * t59;
t73 = qJD(5) ^ 2;
t131 = -t119 + (-t73 - t133) * t68;
t130 = t71 * t133;
t19 = t68 * t43 + t71 * t46;
t10 = -pkin(7) * t106 + t12;
t67 = sin(qJ(5));
t70 = cos(qJ(5));
t7 = -t70 * qJD(2) - t67 * t10;
t114 = t7 * qJD(5);
t111 = qJD(4) * t71;
t112 = qJD(4) * t68;
t30 = t43 * qJDD(1) + qJDD(3);
t108 = qJD(3) * qJD(1);
t113 = pkin(1) * qJDD(1);
t52 = t65 * t113;
t117 = qJDD(1) * qJ(3) + t52;
t33 = t108 + t117;
t5 = t32 * t111 - t38 * t112 + t68 * t30 + t71 * t33;
t3 = -t59 * pkin(7) + t5;
t1 = -t67 * qJDD(2) + t70 * t3 + t114;
t109 = qJD(2) * qJD(5);
t53 = t67 * t109;
t98 = -qJD(5) * t10 - qJDD(2);
t2 = -t67 * t3 + t98 * t70 + t53;
t8 = -t67 * qJD(2) + t70 * t10;
t79 = -(t67 * t8 + t7 * t70) * qJD(5) + t1 * t70 - t2 * t67;
t107 = qJ(1) + pkin(8);
t100 = sin(t107);
t54 = cos(t107);
t26 = -t100 * t68 - t54 * t71;
t27 = -t100 * t71 + t54 * t68;
t95 = g(1) * t26 + g(2) * t27;
t75 = t79 + t95;
t128 = t51 * qJDD(1);
t127 = t59 * pkin(4);
t18 = t71 * t43 - t68 * t46;
t13 = t71 * qJD(3) + t18 * qJD(4);
t123 = t13 * t106;
t14 = t68 * qJD(3) + t19 * qJD(4);
t122 = t14 * t106;
t121 = t67 * t70;
t120 = t70 * t59;
t118 = g(1) * t100 - g(2) * t54;
t62 = t67 ^ 2;
t63 = t70 ^ 2;
t116 = t62 - t63;
t115 = t62 + t63;
t105 = t133 * t121;
t72 = cos(qJ(1));
t104 = t72 * pkin(1) + t54 * pkin(2) + t100 * qJ(3);
t103 = t115 * t59;
t102 = t38 * t111 + t32 * t112 - t71 * t30 + t68 * t33;
t101 = t54 * pkin(3) + t104;
t97 = -0.2e1 * t110 * t121;
t96 = g(1) * t27 - g(2) * t26;
t69 = sin(qJ(1));
t94 = g(1) * t69 - g(2) * t72;
t91 = t7 * t67 - t8 * t70;
t90 = -g(3) - t98;
t4 = t102 + t127;
t89 = -t4 + t96;
t88 = qJDD(3) + t128;
t87 = -t3 - t95 + t132;
t86 = g(1) * t54 + g(2) * t100;
t85 = t96 - t102;
t84 = -t69 * pkin(1) - t100 * pkin(2) + t54 * qJ(3);
t83 = -pkin(7) * qJDD(5) + (t11 + t9 + t126) * qJD(5);
t16 = pkin(4) - t18;
t17 = -pkin(7) + t19;
t82 = -qJDD(5) * t17 + (-t106 * t16 - t13 - t9) * qJD(5);
t81 = -qJDD(5) * t68 + 0.2e1 * t110 * t71;
t80 = -t5 - t95;
t78 = pkin(7) * t73 + t124 + t127 - t89;
t77 = -t16 * t59 + t17 * t73 - t122 + t89;
t76 = -t100 * pkin(3) + t84;
t64 = qJDD(2) - g(3);
t36 = qJDD(5) * t70 - t73 * t67;
t35 = qJDD(5) * t67 + t73 * t70;
t21 = t63 * t59 + t97;
t20 = -t62 * t59 + t97;
t15 = t116 * t110 - t67 * t120;
t6 = [0, 0, 0, 0, 0, qJDD(1), t94, g(1) * t72 + g(2) * t69, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t66 * t113 + t118, -0.2e1 * t52 + t86, 0, (t94 + (t65 ^ 2 + t66 ^ 2) * t113) * pkin(1), 0, 0, 0, qJDD(1), 0, 0, -qJDD(3) + t118 - 0.2e1 * t128, 0, t46 * qJDD(1) + 0.2e1 * t108 + t117 - t86, -g(1) * t84 - g(2) * t104 + t38 * qJD(3) + t33 * t46 + t88 * t51, 0, 0, 0, 0, 0, t59, -t18 * t59 + t122 - t85, t19 * t59 + t123 - t80, 0, -g(1) * t76 - g(2) * t101 - t102 * t18 - t11 * t14 + t12 * t13 + t5 * t19, -t20, -0.2e1 * t15, -t35, t21, -t36, 0, t82 * t67 - t77 * t70, t77 * t67 + t82 * t70, -t17 * t103 - t115 * t123 - t75, t4 * t16 + t9 * t14 - g(1) * (t27 * pkin(4) + t26 * pkin(7) + t76) - g(2) * (-t26 * pkin(4) + t27 * pkin(7) + t101) - t91 * t13 + t79 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, 0, 0, 0, 0, 0, -t36, t35, 0, t91 * qJD(5) - t1 * t67 - t2 * t70 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -qJD(1) ^ 2, -t38 * qJD(1) - t118 + t88, 0, 0, 0, 0, 0, 0, -t133 * t68 - t119, t68 * t59 - t130, 0, (-t102 - t124) * t71 + (t5 + t125) * t68 - t118, 0, 0, 0, 0, 0, 0, t131 * t70 + t81 * t67, -t131 * t67 + t81 * t70, -t68 * t103 + t115 * t130, (t106 * t91 - t4) * t71 + (t79 - t132) * t68 - t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t85 - t124, t80 - t125, 0, 0, t20, 0.2e1 * t15, t35, -t21, t36, 0, t83 * t67 - t78 * t70, t78 * t67 + t83 * t70, -pkin(7) * t103 + t115 * t125 + t75, t89 * pkin(4) + t75 * pkin(7) + t91 * t11 - t9 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t116 * t133, -t67 * t59, t105, -t120, qJDD(5), t8 * qJD(5) + t87 * t67 - t90 * t70 + t53, t114 + t90 * t67 + (t87 + t109) * t70, 0, 0;];
tau_reg = t6;
