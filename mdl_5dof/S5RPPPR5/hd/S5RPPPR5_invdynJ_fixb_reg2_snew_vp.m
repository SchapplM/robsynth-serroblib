% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPPR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:30
% EndTime: 2019-12-31 17:46:34
% DurationCPUTime: 1.02s
% Computational Cost: add. (2519->169), mult. (4938->240), div. (0->0), fcn. (2569->8), ass. (0->106)
t85 = sin(pkin(8));
t87 = cos(pkin(8));
t89 = sin(qJ(5));
t91 = cos(qJ(5));
t60 = (t85 * t89 - t87 * t91) * qJD(1);
t103 = t85 * t91 + t87 * t89;
t61 = t103 * qJD(1);
t125 = t60 * t61;
t134 = qJDD(5) - t125;
t136 = t134 * t89;
t135 = t134 * t91;
t111 = qJD(1) * qJD(4);
t129 = qJD(1) ^ 2;
t90 = sin(qJ(1));
t92 = cos(qJ(1));
t105 = t92 * g(1) + t90 * g(2);
t102 = 0.2e1 * qJD(2) * qJD(1) - t105;
t82 = qJDD(1) * qJ(2);
t100 = t102 + t82;
t128 = pkin(1) + pkin(2);
t50 = -t128 * t129 + t100;
t86 = sin(pkin(7));
t88 = cos(pkin(7));
t109 = t90 * g(1) - t92 * g(2);
t104 = qJDD(2) - t109;
t99 = -t129 * qJ(2) + t104;
t97 = -t128 * qJDD(1) + t99;
t33 = t88 * t50 + t86 * t97;
t30 = -t129 * pkin(3) - qJDD(1) * qJ(4) + t33;
t133 = t30 - 0.2e1 * t111;
t81 = t87 ^ 2;
t94 = t85 ^ 2;
t119 = t81 + t94;
t68 = t119 * t129;
t126 = t87 * pkin(4);
t130 = 0.2e1 * t85;
t84 = g(3) + qJDD(3);
t75 = t87 * t84;
t132 = (pkin(6) * qJDD(1) + t129 * t126 - t30) * t85 + t111 * t130 + t75;
t59 = t103 * qJDD(1);
t107 = t86 * t50 - t88 * t97;
t101 = -t129 * qJ(4) + qJDD(4) + t107;
t29 = qJDD(1) * pkin(3) + t101;
t131 = t29 + (qJ(2) * t86 + pkin(3)) * qJDD(1);
t57 = t60 ^ 2;
t58 = t61 ^ 2;
t113 = t87 * qJDD(1);
t115 = t81 * t129;
t25 = t133 * t87 + t85 * t84;
t18 = -pkin(4) * t115 - pkin(6) * t113 + t25;
t10 = t132 * t89 + t91 * t18;
t9 = -t91 * t132 + t89 * t18;
t3 = t89 * t10 - t91 * t9;
t127 = t85 * t3;
t23 = (pkin(3) + t126) * qJDD(1) + t101 + (-t94 * t129 - t115) * pkin(6);
t124 = t89 * t23;
t37 = qJDD(5) + t125;
t123 = t89 * t37;
t122 = t91 * t23;
t121 = t91 * t37;
t114 = t85 * qJDD(1);
t31 = -t91 * t113 + t89 * t114;
t118 = qJD(5) * t61;
t117 = qJDD(1) * pkin(1);
t116 = t60 * qJD(5);
t112 = t88 * qJDD(1);
t4 = t91 * t10 + t89 * t9;
t106 = qJ(2) * t88 - qJ(4);
t24 = t133 * t85 - t75;
t14 = t85 * t24 + t87 * t25;
t93 = qJD(5) ^ 2;
t67 = t86 * t129 + t112;
t66 = -t86 * qJDD(1) + t88 * t129;
t65 = t119 * qJDD(1);
t64 = t87 * t68;
t63 = t85 * t68;
t56 = -t99 + t117;
t53 = -t58 - t93;
t52 = -t58 + t93;
t51 = t57 - t93;
t46 = -t87 * t112 - t86 * t64;
t45 = t85 * t112 + t86 * t63;
t44 = -t86 * t65 + t88 * t68;
t42 = t116 - t59;
t41 = t59 - 0.2e1 * t116;
t40 = 0.2e1 * t118 + t31;
t39 = t118 + t31;
t35 = -t93 - t57;
t34 = -t57 - t58;
t27 = -t89 * t53 - t121;
t26 = t91 * t53 - t123;
t22 = t91 * t31 - t89 * t59;
t21 = t89 * t31 + t91 * t59;
t20 = t91 * t35 - t136;
t19 = t89 * t35 + t135;
t16 = -t107 * t88 + t86 * t33;
t15 = -t85 * t26 + t87 * t27;
t13 = -t85 * t21 + t87 * t22;
t12 = -t85 * t19 + t87 * t20;
t11 = t86 * t15 + t88 * t41;
t7 = t86 * t12 + t88 * t40;
t6 = t86 * t13 - t88 * t34;
t5 = t86 * t14 - t88 * t29;
t2 = t87 * t4 - t127;
t1 = t86 * t2 - t88 * t23;
t8 = [0, 0, 0, 0, 0, qJDD(1), t109, t105, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -t104 + 0.2e1 * t117, 0, t102 + 0.2e1 * t82, pkin(1) * t56 + qJ(2) * (-t129 * pkin(1) + t100), 0, 0, 0, 0, 0, qJDD(1), -qJ(2) * t66 + t128 * t67 + t107, qJ(2) * t67 + t128 * t66 + t33, 0, qJ(2) * (t107 * t86 + t88 * t33) - t128 * t16, t94 * qJDD(1), t113 * t130, 0, t81 * qJDD(1), 0, 0, -t106 * t64 - t128 * t46 + t131 * t87, t106 * t63 - t128 * t45 - t131 * t85, qJ(2) * (-t88 * t65 - t86 * t68) - pkin(3) * t68 + qJ(4) * t65 - t128 * t44 - t14, qJ(2) * (t88 * t14 + t86 * t29) + pkin(3) * t29 - qJ(4) * t14 - t128 * t5, -t85 * (t89 * t118 + t91 * t42) - t87 * (-t91 * t118 + t89 * t42), -t85 * (t91 * t40 + t89 * t41) - t87 * (t89 * t40 - t91 * t41), -t85 * (-t89 * t52 + t135) - t87 * (t91 * t52 + t136), -t85 * (-t91 * t116 - t89 * t39) - t87 * (-t89 * t116 + t91 * t39), -t85 * (t91 * t51 - t123) - t87 * (t89 * t51 + t121), (-t85 * (t60 * t91 - t61 * t89) - t87 * (t60 * t89 + t61 * t91)) * qJD(5), qJ(2) * (t88 * t12 - t86 * t40) - t85 * (-pkin(6) * t19 + t124) - t87 * (pkin(4) * t40 + pkin(6) * t20 - t122) - pkin(3) * t40 - qJ(4) * t12 - t128 * t7, qJ(2) * (t88 * t15 - t86 * t41) - t85 * (-pkin(6) * t26 + t122) - t87 * (pkin(4) * t41 + pkin(6) * t27 + t124) - pkin(3) * t41 - qJ(4) * t15 - t128 * t11, qJ(2) * (t88 * t13 + t86 * t34) - t85 * (-pkin(6) * t21 - t3) - t87 * (-pkin(4) * t34 + pkin(6) * t22 + t4) + pkin(3) * t34 - qJ(4) * t13 - t128 * t6, qJ(2) * (t88 * t2 + t86 * t23) + pkin(6) * t127 - t87 * (-pkin(4) * t23 + pkin(6) * t4) + pkin(3) * t23 - qJ(4) * t2 - t128 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t129, -t56, 0, 0, 0, 0, 0, 0, -t67, -t66, 0, t16, 0, 0, 0, 0, 0, 0, t46, t45, t44, t5, 0, 0, 0, 0, 0, 0, t7, t11, t6, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87 * t24 + t85 * t25, 0, 0, 0, 0, 0, 0, t87 * t19 + t85 * t20, t87 * t26 + t85 * t27, t87 * t21 + t85 * t22, t87 * t3 + t85 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, -t114, -t68, t29, 0, 0, 0, 0, 0, 0, -t40, -t41, t34, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t58 - t57, -t59, -t125, t31, qJDD(5), -t9, -t10, 0, 0;];
tauJ_reg = t8;
