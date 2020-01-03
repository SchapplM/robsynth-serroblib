% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRR7
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
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:56
% EndTime: 2019-12-31 17:59:58
% DurationCPUTime: 0.87s
% Computational Cost: add. (755->182), mult. (1380->254), div. (0->0), fcn. (862->10), ass. (0->105)
t53 = sin(pkin(8));
t38 = t53 * pkin(1) + qJ(3);
t109 = t38 * qJDD(1);
t48 = qJ(1) + pkin(8);
t43 = sin(t48);
t44 = cos(t48);
t78 = g(1) * t43 - g(2) * t44;
t33 = qJD(1) * t38;
t131 = -t33 * qJD(1) - t78;
t54 = cos(pkin(8));
t40 = -t54 * pkin(1) - pkin(2);
t130 = t40 * qJDD(1);
t36 = -pkin(6) + t40;
t23 = t36 * qJD(1) + qJD(3);
t56 = sin(qJ(4));
t59 = cos(qJ(4));
t12 = -t56 * qJD(2) + t59 * t23;
t10 = -qJD(4) * pkin(4) - t12;
t104 = qJD(4) * t59;
t79 = t56 * pkin(4) - t59 * pkin(7);
t20 = t38 + t79;
t91 = qJD(1) * qJD(4);
t86 = t59 * t91;
t93 = t56 * qJDD(1);
t25 = qJDD(5) + t86 + t93;
t37 = t56 * qJD(1) + qJD(5);
t21 = t36 * qJDD(1) + qJDD(3);
t13 = t59 * qJD(2) + t56 * t23;
t99 = t13 * qJD(4);
t4 = -qJDD(4) * pkin(4) + t56 * qJDD(2) - t59 * t21 + t99;
t100 = t12 * qJD(4);
t14 = t20 * qJD(1);
t83 = -qJDD(4) * pkin(7) - qJD(5) * t14 - t59 * qJDD(2) - t56 * t21 - t100;
t129 = -(qJD(5) * t20 + t36 * t104) * t37 + t4 * t59 + (-t10 * qJD(4) - t36 * t25 + t83) * t56;
t58 = cos(qJ(5));
t110 = t58 * t59;
t101 = qJD(5) * t59;
t55 = sin(qJ(5));
t97 = t58 * qJD(4);
t69 = t55 * t101 + t56 * t97;
t128 = t25 * t110 - t69 * t37;
t127 = 0.2e1 * qJD(4) * t33 + qJDD(4) * t36;
t106 = qJD(1) * t59;
t98 = t55 * qJD(4);
t29 = t58 * t106 + t98;
t88 = t56 * t98;
t92 = t59 * qJDD(1);
t9 = -qJD(1) * t88 + t29 * qJD(5) - t58 * qJDD(4) + t55 * t92;
t80 = pkin(4) * t59 + pkin(7) * t56;
t126 = (pkin(7) * qJD(5) + t80 * qJD(1)) * t37 + t78 * t59 - g(3) * t56 + t4;
t123 = t56 * t9;
t8 = -t69 * qJD(1) + qJD(5) * t97 + t55 * qJDD(4) + t58 * t92;
t122 = t8 * t55;
t121 = t29 * t104 + t8 * t56;
t120 = t20 * t25;
t27 = t55 * t106 - t97;
t119 = t27 * t37;
t118 = t29 * t37;
t117 = t29 * t58;
t116 = t36 * t56;
t115 = t36 * t59;
t114 = t55 * t25;
t113 = t55 * t56;
t112 = t56 * t58;
t111 = t58 * t25;
t51 = t59 ^ 2;
t108 = t56 ^ 2 - t51;
t61 = qJD(4) ^ 2;
t62 = qJD(1) ^ 2;
t107 = -t61 - t62;
t105 = qJD(4) * t27;
t103 = qJD(5) * t56;
t102 = qJD(5) * t58;
t52 = qJDD(2) - g(3);
t95 = qJDD(4) * t56;
t94 = qJDD(4) * t59;
t90 = qJD(3) * qJD(1);
t11 = qJD(4) * pkin(7) + t13;
t26 = t80 * qJD(4) + qJD(3);
t7 = t26 * qJD(1) + t79 * qJDD(1) + t109;
t84 = qJD(5) * t11 - t7;
t81 = qJD(1) + t103;
t57 = sin(qJ(1));
t60 = cos(qJ(1));
t77 = g(1) * t57 - g(2) * t60;
t76 = -qJD(4) * t23 - t52;
t73 = g(3) * t59 + t83;
t72 = t37 * t102 + t114;
t70 = qJDD(3) + t130;
t68 = -g(1) * t44 - g(2) * t43 + t109;
t67 = -pkin(7) * t25 + (t10 + t12) * t37;
t66 = qJD(2) * qJD(4) - t131 - t21;
t24 = t90 + t109;
t64 = -t36 * t61 + t24 + t68 + t90;
t32 = -t61 * t56 + t94;
t31 = -t61 * t59 - t95;
t22 = t37 * t88;
t18 = t44 * t112 - t43 * t55;
t17 = t44 * t113 + t43 * t58;
t16 = t43 * t112 + t44 * t55;
t15 = -t43 * t113 + t44 * t58;
t5 = t58 * t7;
t2 = t58 * t11 + t55 * t14;
t1 = -t55 * t11 + t58 * t14;
t3 = [qJDD(1), t77, g(1) * t60 + g(2) * t57, (t77 + (t53 ^ 2 + t54 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(3) + 0.2e1 * t130 - t78, t68 + 0.2e1 * t90 + t109, t24 * t38 + t33 * qJD(3) + t70 * t40 - g(1) * (-t57 * pkin(1) - t43 * pkin(2) + t44 * qJ(3)) - g(2) * (t60 * pkin(1) + t44 * pkin(2) + t43 * qJ(3)), t51 * qJDD(1) - 0.2e1 * t56 * t86, 0.2e1 * t108 * t91 - 0.2e1 * t56 * t92, t32, t31, 0, t127 * t59 + t64 * t56, -t127 * t56 + t64 * t59, t8 * t110 - t69 * t29, (t27 * t58 + t29 * t55) * t56 * qJD(4) + (-t122 - t58 * t9 + (t27 * t55 - t117) * qJD(5)) * t59, t121 + t128, -t123 + t22 + (-t72 - t105) * t59, t37 * t104 + t25 * t56, -t9 * t115 - g(1) * t18 - g(2) * t16 + t5 * t56 + (t1 * t59 + t27 * t116) * qJD(4) + (t120 + t26 * t37 + (t10 * t59 + (-t36 * t37 - t11) * t56) * qJD(5)) * t58 + t129 * t55, -t8 * t115 + g(1) * t17 - g(2) * t15 + (t29 * t116 - t2 * t59) * qJD(4) + (-(-t36 * t103 + t26) * t37 - t120 + t84 * t56 - t10 * t101) * t55 + t129 * t58; 0, 0, 0, t52, 0, 0, t52, 0, 0, 0, 0, 0, t31, -t32, 0, 0, 0, 0, 0, t123 + t22 + (-t72 + t105) * t59, t121 - t128; 0, 0, 0, 0, qJDD(1), -t62, t70 + t131, 0, 0, 0, 0, 0, t107 * t56 + t94, t107 * t59 - t95, 0, 0, 0, 0, 0, -t59 * t9 + (t105 - t114) * t56 + (-t81 * t58 - t59 * t98) * t37, -t59 * t8 + (qJD(4) * t29 - t111) * t56 + (t81 * t55 - t59 * t97) * t37; 0, 0, 0, 0, 0, 0, 0, t59 * t62 * t56, -t108 * t62, t92, -t93, qJDD(4), t76 * t56 - t66 * t59 + t99, t66 * t56 + t76 * t59 + t100, t37 * t117 + t122, (t8 - t119) * t58 + (-t9 - t118) * t55, (t37 * t112 - t29 * t59) * qJD(1) + t72, -qJD(5) * t55 * t37 + t111 + (-t37 * t113 + t27 * t59) * qJD(1), -t37 * t106, -pkin(4) * t9 - t1 * t106 - t126 * t58 - t13 * t27 + t67 * t55, -pkin(4) * t8 + t2 * t106 + t126 * t55 - t13 * t29 + t67 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t27, -t27 ^ 2 + t29 ^ 2, t8 + t119, t118 - t9, t25, -g(1) * t15 - g(2) * t17 - t10 * t29 - t11 * t102 + t2 * t37 + t73 * t55 + t5, g(1) * t16 - g(2) * t18 + t1 * t37 + t10 * t27 + t84 * t55 + t73 * t58;];
tau_reg = t3;
