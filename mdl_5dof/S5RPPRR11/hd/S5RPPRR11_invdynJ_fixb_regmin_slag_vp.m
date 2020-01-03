% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRR11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR11_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:56
% EndTime: 2019-12-31 18:05:58
% DurationCPUTime: 0.80s
% Computational Cost: add. (652->194), mult. (1193->247), div. (0->0), fcn. (713->6), ass. (0->113)
t54 = sin(qJ(1));
t57 = cos(qJ(1));
t111 = g(1) * t54 - g(2) * t57;
t51 = pkin(1) + qJ(3);
t135 = qJD(1) * t51;
t56 = cos(qJ(4));
t104 = qJD(4) * t56;
t31 = qJ(2) * qJD(1) + qJD(3);
t26 = -pkin(6) * qJD(1) + t31;
t115 = t56 * t26;
t11 = -qJD(4) * pkin(4) - t115;
t92 = qJD(1) * qJD(4);
t83 = t56 * t92;
t53 = sin(qJ(4));
t95 = t53 * qJDD(1);
t16 = qJDD(5) + t83 + t95;
t75 = t53 * pkin(4) - t56 * pkin(7);
t24 = t75 + t51;
t29 = t53 * qJD(1) + qJD(5);
t50 = -pkin(6) + qJ(2);
t105 = qJD(4) * t53;
t44 = qJDD(1) * qJ(2);
t45 = qJD(1) * qJD(2);
t85 = qJDD(3) + t44 + t45;
t21 = -pkin(6) * qJDD(1) + t85;
t7 = -qJDD(4) * pkin(4) + t26 * t105 - t56 * t21;
t9 = t24 * qJD(1) - qJD(2);
t82 = -qJDD(4) * pkin(7) - qJD(5) * t9 - t26 * t104 - t53 * t21;
t134 = -(qJD(5) * t24 + t50 * t104) * t29 + t7 * t56 + (-qJD(2) * t29 - t11 * qJD(4) - t50 * t16 + t82) * t53;
t27 = -qJD(2) + t135;
t133 = (qJD(2) + t27 + t135) * qJD(4) + qJDD(4) * t50;
t52 = sin(qJ(5));
t100 = t52 * qJD(4);
t108 = qJD(1) * t56;
t55 = cos(qJ(5));
t20 = t55 * t108 + t100;
t84 = t53 * t92;
t94 = t56 * qJDD(1);
t6 = t20 * qJD(5) - t55 * qJDD(4) + (-t84 + t94) * t52;
t131 = g(3) * t53;
t74 = g(1) * t57 + g(2) * t54;
t76 = pkin(4) * t56 + pkin(7) * t53;
t132 = (pkin(7) * qJD(5) + t76 * qJD(1)) * t29 + t74 * t56 + t7 - t131;
t39 = 0.2e1 * t45;
t130 = g(3) * t56;
t101 = qJD(5) * t56;
t99 = t55 * qJD(4);
t65 = -t52 * t101 - t53 * t99;
t5 = t65 * qJD(1) + qJD(5) * t99 + t52 * qJDD(4) + t55 * t94;
t129 = t5 * t52;
t128 = t56 * t5;
t18 = t52 * t108 - t99;
t127 = t18 * t29;
t126 = t20 * t29;
t125 = t24 * t16;
t124 = t29 * t52;
t123 = t29 * t55;
t122 = t52 * t16;
t121 = t53 * t26;
t120 = t54 * t52;
t119 = t54 * t55;
t118 = t55 * t16;
t117 = t56 * t18;
t116 = t56 * t20;
t114 = t57 * t52;
t113 = t57 * t55;
t112 = t57 * pkin(1) + t54 * qJ(2);
t48 = t56 ^ 2;
t110 = t53 ^ 2 - t48;
t58 = qJD(4) ^ 2;
t59 = qJD(1) ^ 2;
t109 = -t58 - t59;
t107 = qJD(4) * t18;
t106 = qJD(4) * t20;
t103 = qJD(5) * t53;
t102 = qJD(5) * t55;
t97 = qJDD(4) * t53;
t96 = t51 * qJDD(1);
t49 = qJDD(1) * pkin(1);
t93 = t49 - qJDD(2);
t91 = qJD(3) * qJD(1);
t90 = t53 * t124;
t89 = t53 * t123;
t87 = qJDD(2) - t111;
t43 = qJDD(1) * qJ(3);
t86 = -t43 - t93;
t10 = qJD(4) * pkin(7) + t121;
t17 = t76 * qJD(4) + qJD(3);
t4 = t17 * qJD(1) + t75 * qJDD(1) - t86;
t80 = qJD(5) * t10 - t4;
t79 = -0.2e1 * t83;
t78 = qJD(1) + t103;
t77 = -t49 + t87;
t72 = -t43 + t77;
t71 = t82 + t130;
t69 = t29 * t102 + t122;
t68 = qJD(5) * t124 - t118;
t66 = 0.2e1 * t44 + t39 - t74;
t64 = t27 * qJD(1) + t74;
t63 = -t21 + t64;
t62 = -pkin(7) * t16 + (t11 + t115) * t29;
t22 = -t86 + t91;
t61 = -t50 * t58 + t111 + t22 + t91 + t96;
t37 = t57 * qJ(2);
t34 = qJDD(4) * t56;
t15 = t53 * t113 - t120;
t14 = -t53 * t114 - t119;
t13 = -t53 * t119 - t114;
t12 = t53 * t120 - t113;
t3 = t55 * t4;
t2 = t55 * t10 + t52 * t9;
t1 = -t52 * t10 + t55 * t9;
t8 = [qJDD(1), t111, t74, -0.2e1 * t49 + t87, t66, t93 * pkin(1) - g(1) * (-t54 * pkin(1) + t37) - g(2) * t112 + (t44 + t39) * qJ(2), qJDD(3) + t66, -t72 + 0.2e1 * t91 + t96, t22 * t51 + t27 * qJD(3) + t85 * qJ(2) + t31 * qJD(2) - g(1) * (-t51 * t54 + t37) - g(2) * (t57 * qJ(3) + t112), t48 * qJDD(1) + t53 * t79, 0.2e1 * t110 * t92 - 0.2e1 * t53 * t94, -t58 * t53 + t34, -t58 * t56 - t97, 0, t133 * t56 + t61 * t53, -t133 * t53 + t61 * t56, t55 * t128 + t65 * t20, (t18 * t55 + t20 * t52) * t105 + (-t129 - t55 * t6 + (t18 * t52 - t20 * t55) * qJD(5)) * t56, (-t29 * t99 + t5) * t53 + (-t68 + t106) * t56, (t29 * t100 - t6) * t53 + (-t69 - t107) * t56, t29 * t104 + t16 * t53, -g(1) * t13 - g(2) * t15 + (t50 * t107 + t3) * t53 + (-qJD(2) * t18 + t1 * qJD(4) - t50 * t6) * t56 + (t125 + t17 * t29 + (t11 * t56 + (-t29 * t50 - t10) * t53) * qJD(5)) * t55 + t134 * t52, t50 * t20 * t105 - g(1) * t12 - g(2) * t14 + (-qJD(2) * t20 - t2 * qJD(4) - t50 * t5) * t56 + (-(-t50 * t103 + t17) * t29 - t125 + t80 * t53 - t11 * t101) * t52 + t134 * t55; 0, 0, 0, qJDD(1), -t59, -t59 * qJ(2) + t77, -t59, -qJDD(1), (-qJD(3) - t31) * qJD(1) + t72, 0, 0, 0, 0, 0, t79 - t95, 0.2e1 * t84 - t94, 0, 0, 0, 0, 0, (t90 + t117) * qJD(1) + t68, (t89 + t116) * qJD(1) + t69; 0, 0, 0, 0, 0, 0, qJDD(1), -t59, -t64 + t85, 0, 0, 0, 0, 0, t109 * t53 + t34, t109 * t56 - t97, 0, 0, 0, 0, 0, -t56 * t6 + (t107 - t122) * t53 + (-t56 * t100 - t78 * t55) * t29, -t128 + (t106 - t118) * t53 + (t78 * t52 - t56 * t99) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t59 * t53, -t110 * t59, t94, -t95, qJDD(4), -t63 * t56 + t131, t63 * t53 + t130, t20 * t123 + t129, (t5 - t127) * t55 + (-t6 - t126) * t52, (t89 - t116) * qJD(1) + t69, (-t90 + t117) * qJD(1) - t68, -t29 * t108, -pkin(4) * t6 - t1 * t108 - t18 * t121 - t132 * t55 + t62 * t52, -pkin(4) * t5 + t2 * t108 - t20 * t121 + t132 * t52 + t62 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t18, -t18 ^ 2 + t20 ^ 2, t5 + t127, t126 - t6, t16, -g(1) * t14 + g(2) * t12 - t10 * t102 - t11 * t20 + t2 * t29 + t71 * t52 + t3, g(1) * t15 - g(2) * t13 + t1 * t29 + t11 * t18 + t80 * t52 + t71 * t55;];
tau_reg = t8;
