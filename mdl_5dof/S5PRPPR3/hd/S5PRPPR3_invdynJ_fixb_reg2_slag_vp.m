% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPPR3
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:27:01
% EndTime: 2019-12-05 15:27:02
% DurationCPUTime: 0.77s
% Computational Cost: add. (776->152), mult. (1396->182), div. (0->0), fcn. (993->10), ass. (0->101)
t56 = sin(qJ(2));
t110 = qJD(1) * t56;
t58 = cos(qJ(2));
t105 = t58 * qJD(1);
t32 = qJD(2) * pkin(2) + t105;
t51 = sin(pkin(8));
t53 = cos(pkin(8));
t15 = t53 * t110 + t51 * t32;
t13 = qJD(2) * qJ(4) + t15;
t48 = qJ(2) + pkin(8);
t43 = sin(t48);
t52 = sin(pkin(7));
t54 = cos(pkin(7));
t81 = g(1) * t54 + g(2) * t52;
t130 = t81 * t43;
t44 = cos(t48);
t91 = g(3) * t44 - t130;
t133 = -t13 * qJD(2) + t91;
t96 = qJD(1) * qJD(2);
t98 = t56 * qJDD(1);
t132 = t58 * t96 + t98;
t131 = g(1) * t52 - g(2) * t54;
t69 = -qJDD(3) + t131;
t129 = -t51 * t56 + t53 * t58;
t64 = -g(3) * t43 - t81 * t44;
t125 = pkin(3) + pkin(6);
t34 = t51 * t110;
t14 = t53 * t32 - t34;
t75 = qJD(4) - t14;
t11 = -t125 * qJD(2) + t75;
t55 = sin(qJ(5));
t57 = cos(qJ(5));
t7 = -t55 * qJD(3) + t57 * t11;
t111 = t7 * qJD(5);
t109 = qJDD(2) * pkin(2);
t45 = t58 * qJDD(1);
t21 = -t56 * t96 + t109 + t45;
t93 = t132 * t51 - t53 * t21;
t84 = qJDD(4) + t93;
t4 = -t125 * qJDD(2) + t84;
t1 = t57 * qJDD(3) + t55 * t4 + t111;
t3 = t57 * t4;
t8 = t57 * qJD(3) + t55 * t11;
t2 = -t8 * qJD(5) - t55 * qJDD(3) + t3;
t62 = -(t55 * t7 - t57 * t8) * qJD(5) + t1 * t55 + t2 * t57;
t23 = t51 * t58 + t53 * t56;
t18 = t23 * qJD(1);
t42 = -t53 * pkin(2) - pkin(3);
t36 = -pkin(6) + t42;
t38 = t51 * pkin(2) + qJ(4);
t127 = qJDD(5) * t36 + (qJD(2) * t38 + t13 - t18) * qJD(5);
t17 = qJD(2) * t23;
t126 = 0.2e1 * qJD(5) * t17 - qJDD(5) * t129;
t124 = pkin(2) * t56;
t119 = t51 * t21;
t116 = t55 * t57;
t49 = t55 ^ 2;
t50 = t57 ^ 2;
t115 = t49 - t50;
t114 = t49 + t50;
t59 = qJD(5) ^ 2;
t60 = qJD(2) ^ 2;
t113 = -t59 - t60;
t112 = qJ(4) * t44;
t107 = t17 * qJD(2);
t106 = t18 * qJD(2);
t90 = t53 * t105;
t20 = -t34 + t90;
t104 = qJD(4) - t20;
t103 = qJDD(1) - g(3);
t100 = qJDD(5) * t55;
t99 = qJDD(5) * t57;
t97 = t57 * qJDD(2);
t95 = qJD(2) * qJD(5);
t94 = t60 * t116;
t92 = t58 * pkin(2) + t44 * pkin(3) + t43 * qJ(4);
t89 = t53 * t98;
t85 = t114 * qJDD(2);
t83 = t95 * t116;
t82 = -pkin(3) * t43 - t124;
t80 = t20 - t90;
t78 = t8 * t55 + t7 * t57;
t19 = t129 * qJD(2);
t68 = t89 + t119;
t5 = qJDD(2) * qJ(4) + (qJD(4) + t90) * qJD(2) + t68;
t73 = t13 * t19 + t5 * t23 - g(3);
t71 = -qJDD(2) * t129 + t107;
t70 = t19 * qJD(2) + t23 * qJDD(2);
t67 = t104 * t13 + t5 * t38;
t6 = -qJDD(2) * pkin(3) + t84;
t66 = t129 * t59 + t70;
t65 = t91 + t93 - t106;
t63 = -g(3) * t58 + t81 * t56;
t61 = t104 * qJD(2) + t38 * qJDD(2) - t36 * t59 + t5 + t64;
t28 = t54 * t112;
t27 = t52 * t112;
t26 = -t59 * t55 + t99;
t25 = -t59 * t57 - t100;
t12 = -qJD(2) * pkin(3) + t75;
t10 = t132 * t53 + t119;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, 0, 0, 0, 0, 0, t58 * qJDD(2) - t60 * t56, -qJDD(2) * t56 - t60 * t58, 0, -g(3) + (t56 ^ 2 + t58 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t71, -t70, 0, t10 * t23 - t129 * t93 - t14 * t17 + t15 * t19 - g(3), 0, 0, 0, 0, 0, 0, 0, t71, t70, t12 * t17 - t129 * t6 + t73, 0, 0, 0, 0, 0, 0, t126 * t57 + t66 * t55, -t126 * t55 + t66 * t57, -t114 * t107 + t129 * t85, -t129 * t62 + t78 * t17 + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t45 + t63, -t103 * t56 + t81 * t58, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t53 * t109 - t65, -t89 + (-t21 - t109) * t51 + t80 * qJD(2) - t64, 0, t14 * t18 - t15 * t20 + (t10 * t51 - t53 * t93 + t63) * pkin(2), qJDD(2), 0, 0, 0, 0, 0, 0, qJDD(4) + (-pkin(3) + t42) * qJDD(2) + t65, (qJ(4) + t38) * qJDD(2) + (0.2e1 * qJD(4) - t80) * qJD(2) + t64 + t68, t6 * t42 - t12 * t18 - g(1) * (t82 * t54 + t28) - g(2) * (t82 * t52 + t27) - g(3) * t92 + t67, t50 * qJDD(2) - 0.2e1 * t83, 0.2e1 * t115 * t95 - 0.2e1 * t55 * t97, t26, t49 * qJDD(2) + 0.2e1 * t83, t25, 0, t127 * t57 + t61 * t55, -t127 * t55 + t61 * t57, t114 * t106 - t36 * t85 - t62 - t91, -g(1) * (-t54 * t124 + t28) - g(2) * (-t52 * t124 + t27) - g(3) * (t44 * pkin(6) + t92) - t78 * t18 + t62 * t36 + t67 + t125 * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, 0, 0, 0, 0, 0, t25, -t26, 0, -t78 * qJD(5) + t1 * t57 - t2 * t55 - t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t60, t6 + t133, 0, 0, 0, 0, 0, 0, t113 * t55 + t99, t113 * t57 - t100, -t85, t62 + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t115 * t60, t97, -t94, -t55 * qJDD(2), qJDD(5), t133 * t57 + t69 * t55 + t3, t111 + (-qJD(5) * t11 + t69) * t57 + (qJD(5) * qJD(3) - t133 - t4) * t55, 0, 0;];
tau_reg = t9;
