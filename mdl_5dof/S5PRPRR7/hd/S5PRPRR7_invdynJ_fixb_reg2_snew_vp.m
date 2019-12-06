% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPRR7
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPRR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:45
% EndTime: 2019-12-05 16:00:49
% DurationCPUTime: 0.91s
% Computational Cost: add. (1891->168), mult. (3628->221), div. (0->0), fcn. (2351->8), ass. (0->111)
t90 = sin(qJ(5));
t91 = sin(qJ(4));
t93 = cos(qJ(5));
t94 = cos(qJ(4));
t53 = (-t94 * t90 - t91 * t93) * qJD(2);
t111 = qJD(2) * t94;
t55 = -t90 * t91 * qJD(2) + t93 * t111;
t126 = t55 * t53;
t79 = qJDD(4) + qJDD(5);
t132 = t79 + t126;
t135 = t132 * t90;
t134 = t132 * t93;
t80 = qJD(4) + qJD(5);
t125 = t80 * t53;
t109 = qJD(2) * qJD(4);
t104 = t94 * t109;
t110 = t91 * qJDD(2);
t61 = -t104 - t110;
t105 = t91 * t109;
t76 = t94 * qJDD(2);
t62 = t76 - t105;
t29 = t53 * qJD(5) + t90 * t61 + t93 * t62;
t133 = t29 + t125;
t97 = qJD(2) ^ 2;
t113 = t94 * t97;
t128 = t62 * pkin(7);
t86 = sin(pkin(8));
t87 = cos(pkin(8));
t69 = -t87 * g(1) - t86 * g(2);
t84 = -g(3) + qJDD(1);
t92 = sin(qJ(2));
t95 = cos(qJ(2));
t47 = -t92 * t69 + t95 * t84;
t102 = qJDD(3) - t47;
t85 = qJDD(2) * pkin(2);
t41 = -t97 * qJ(3) + t102 - t85;
t38 = -qJDD(2) * pkin(6) + t41;
t35 = t94 * t38;
t68 = -t86 * g(1) + t87 * g(2);
t131 = qJDD(4) * pkin(4) - t128 + t35 + (-pkin(4) * t113 - pkin(7) * t109 - t68) * t91;
t51 = t53 ^ 2;
t52 = t55 ^ 2;
t78 = t80 ^ 2;
t130 = 2 * qJD(3);
t129 = -pkin(6) - pkin(2);
t82 = t91 ^ 2;
t122 = t82 * t97;
t27 = t91 * t38 + t94 * t68;
t72 = qJD(4) * pkin(4) - pkin(7) * t111;
t16 = -pkin(4) * t122 + t61 * pkin(7) - qJD(4) * t72 + t27;
t6 = -t93 * t131 + t90 * t16;
t117 = t93 * t16;
t7 = t131 * t90 + t117;
t2 = -t93 * t6 + t90 * t7;
t127 = t94 * t2;
t124 = t80 * t90;
t123 = t80 * t93;
t83 = t94 ^ 2;
t121 = t83 * t97;
t88 = t97 * pkin(6);
t108 = qJDD(2) * qJ(3);
t48 = t95 * t69 + t92 * t84;
t99 = -t97 * pkin(2) + t108 + t48;
t25 = -t61 * pkin(4) - pkin(7) * t122 - t88 + (t72 * t94 + t130) * qJD(2) + t99;
t120 = t90 * t25;
t33 = -t126 + t79;
t119 = t90 * t33;
t107 = t91 * t113;
t118 = t91 * (qJDD(4) + t107);
t116 = t93 * t25;
t115 = t93 * t33;
t71 = qJDD(4) - t107;
t114 = t94 * t71;
t112 = t82 + t83;
t106 = qJD(2) * t130;
t3 = t90 * t6 + t93 * t7;
t26 = t91 * t68 - t35;
t103 = -t93 * t61 + t90 * t62;
t12 = -t94 * t26 + t91 * t27;
t100 = (-qJD(5) + t80) * t55 - t103;
t40 = t99 + t106;
t96 = qJD(4) ^ 2;
t67 = t112 * t97;
t66 = -t95 * qJDD(2) + t92 * t97;
t65 = t92 * qJDD(2) + t95 * t97;
t64 = t112 * qJDD(2);
t63 = t76 - 0.2e1 * t105;
t60 = 0.2e1 * t104 + t110;
t46 = -t52 + t78;
t45 = t51 - t78;
t44 = -t52 - t78;
t43 = -t118 + t94 * (-t96 - t121);
t42 = t91 * (-t96 - t122) + t114;
t37 = t52 - t51;
t36 = t40 - t88;
t31 = -t78 - t51;
t30 = -t51 - t52;
t28 = -t55 * qJD(5) - t103;
t24 = -t90 * t44 - t115;
t23 = t93 * t44 - t119;
t22 = -t125 + t29;
t17 = (qJD(5) + t80) * t55 + t103;
t15 = t93 * t31 - t135;
t14 = t90 * t31 + t134;
t11 = t94 * t23 + t91 * t24;
t10 = t100 * t93 + t90 * t22;
t9 = t100 * t90 - t93 * t22;
t8 = t94 * t14 + t91 * t15;
t4 = t91 * t10 + t94 * t9;
t1 = t91 * t3 + t127;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, -t66, -t65, 0, t95 * t47 + t92 * t48, 0, 0, 0, 0, 0, 0, 0, t66, t65, t92 * t40 - t95 * t41, 0, 0, 0, 0, 0, 0, -t95 * t42 + t92 * t60, -t95 * t43 + t92 * t63, t95 * t64 - t92 * t67, -t95 * t12 + t92 * t36, 0, 0, 0, 0, 0, 0, t92 * t17 - t95 * t8, -t95 * t11 + t133 * t92, t92 * t30 - t95 * t4, -t95 * t1 + t92 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t47, -t48, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, t102 - 0.2e1 * t85, t106 + t48 + 0.2e1 * t108, -pkin(2) * t41 + qJ(3) * t40, (t62 - t105) * t94, -t94 * t60 - t91 * t63, t114 - t91 * (t96 - t121), (-t61 + t104) * t91, t94 * (-t96 + t122) - t118, 0, qJ(3) * t60 + t129 * t42 + t91 * t36, qJ(3) * t63 + t129 * t43 + t94 * t36, -qJ(3) * t67 - t129 * t64 - t12, qJ(3) * t36 + t129 * t12, t94 * (-t55 * t124 + t93 * t29) - t91 * (t55 * t123 + t90 * t29), t94 * (-t133 * t90 - t93 * t17) - t91 * (t133 * t93 - t90 * t17), t94 * (-t90 * t46 + t134) - t91 * (t93 * t46 + t135), t94 * (-t53 * t123 - t90 * t28) - t91 * (-t53 * t124 + t93 * t28), t94 * (t93 * t45 - t119) - t91 * (t90 * t45 + t115), (t94 * (t53 * t93 + t55 * t90) - t91 * (t53 * t90 - t55 * t93)) * t80, t94 * (-pkin(7) * t14 + t120) - t91 * (-pkin(4) * t17 + pkin(7) * t15 - t116) + qJ(3) * t17 + t129 * t8, t94 * (-pkin(7) * t23 + t116) - t91 * (-pkin(4) * t133 + pkin(7) * t24 + t120) + qJ(3) * t133 + t129 * t11, t94 * (-pkin(7) * t9 - t2) - t91 * (-pkin(4) * t30 + pkin(7) * t10 + t3) + qJ(3) * t30 + t129 * t4, -pkin(7) * t127 - t91 * (-pkin(4) * t25 + pkin(7) * t3) + qJ(3) * t25 + t129 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t97, t41, 0, 0, 0, 0, 0, 0, t42, t43, -t64, t12, 0, 0, 0, 0, 0, 0, t8, t11, t4, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, (-t82 + t83) * t97, t76, -t107, -t110, qJDD(4), -t26, -t27, 0, 0, -t126, t37, t22, t126, t100, t79, pkin(4) * t14 - t6, -t117 - t90 * (-pkin(7) * t105 - t128 - t26) + (-t90 * t71 + t23) * pkin(4), pkin(4) * t9, pkin(4) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, t37, t22, t126, t100, t79, -t6, -t7, 0, 0;];
tauJ_reg = t5;
