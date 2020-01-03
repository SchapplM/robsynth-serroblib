% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPPR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:15
% EndTime: 2019-12-31 17:45:18
% DurationCPUTime: 0.74s
% Computational Cost: add. (1727->132), mult. (3502->188), div. (0->0), fcn. (2124->8), ass. (0->95)
t74 = sin(pkin(8));
t76 = cos(pkin(8));
t79 = sin(qJ(5));
t81 = cos(qJ(5));
t94 = t74 * t81 + t76 * t79;
t45 = t94 * qJD(1);
t47 = (-t74 * t79 + t76 * t81) * qJD(1);
t115 = t47 * t45;
t123 = qJDD(5) - t115;
t125 = t123 * t79;
t124 = t123 * t81;
t68 = t74 ^ 2;
t69 = t76 ^ 2;
t108 = t68 + t69;
t84 = qJD(1) ^ 2;
t122 = t108 * t84;
t103 = qJD(1) * qJD(4);
t73 = qJDD(1) * pkin(2);
t80 = sin(qJ(1));
t82 = cos(qJ(1));
t101 = t80 * g(1) - t82 * g(2);
t53 = qJDD(1) * pkin(1) + t101;
t95 = t82 * g(1) + t80 * g(2);
t54 = -t84 * pkin(1) - t95;
t75 = sin(pkin(7));
t77 = cos(pkin(7));
t97 = t77 * t53 - t75 * t54;
t28 = -t84 * qJ(3) + qJDD(3) - t73 - t97;
t88 = -qJDD(1) * qJ(4) + t28;
t121 = (-0.2e1 * t103 + t88) * t76;
t22 = t94 * qJDD(1);
t120 = -t84 * qJ(4) + qJDD(4);
t42 = t45 ^ 2;
t43 = t47 ^ 2;
t118 = -0.2e1 * t74;
t117 = pkin(1) * (t75 * qJDD(1) + t77 * t84);
t116 = pkin(1) * t75;
t114 = t68 * t84;
t109 = t75 * t53 + t77 * t54;
t93 = 0.2e1 * qJD(3) * qJD(1) + t109;
t90 = -t84 * pkin(2) + t93;
t98 = t74 * pkin(4) + qJ(3);
t17 = t98 * qJDD(1) + t90 + (-t69 * t84 - t114) * pkin(6) + t120;
t113 = t79 * t17;
t31 = qJDD(5) + t115;
t112 = t79 * t31;
t111 = t81 * t17;
t110 = t81 * t31;
t107 = t45 * qJD(5);
t106 = t47 * qJD(5);
t105 = t74 * qJDD(1);
t104 = t76 * qJDD(1);
t102 = qJDD(1) * qJ(3);
t72 = -g(3) + qJDD(2);
t21 = t103 * t118 + t76 * t72 + t74 * t88;
t16 = -pkin(4) * t114 - pkin(6) * t105 + t21;
t87 = -pkin(6) * t104 + (-pkin(4) * t76 * t84 - t72) * t74 + t121;
t7 = t79 * t16 - t81 * t87;
t8 = t81 * t16 + t79 * t87;
t3 = t79 * t7 + t81 * t8;
t99 = qJ(3) + t116;
t96 = pkin(1) * t77 + pkin(2) + qJ(4);
t2 = -t81 * t7 + t79 * t8;
t1 = t76 * t2 + t74 * t3;
t44 = t81 * t104 - t79 * t105;
t20 = t74 * t72 - t121;
t10 = -t76 * t20 + t74 * t21;
t92 = t98 + t116;
t91 = -pkin(1) * (-t77 * qJDD(1) + t75 * t84) + t97;
t27 = t90 + t102;
t26 = t27 + t120;
t89 = t99 * qJDD(1) + t26;
t83 = qJD(5) ^ 2;
t55 = t108 * qJDD(1);
t52 = t74 * t122;
t51 = t76 * t122;
t40 = -t43 - t83;
t39 = -t43 + t83;
t38 = t42 - t83;
t36 = t44 - t107;
t35 = t44 - 0.2e1 * t107;
t34 = -t22 - t106;
t33 = 0.2e1 * t106 + t22;
t29 = -t83 - t42;
t23 = -t42 - t43;
t19 = -t79 * t40 - t110;
t18 = t81 * t40 - t112;
t15 = -t81 * t22 + t79 * t44;
t14 = -t79 * t22 - t81 * t44;
t13 = t81 * t29 - t125;
t12 = t79 * t29 + t124;
t9 = t76 * t18 + t74 * t19;
t5 = t76 * t14 + t74 * t15;
t4 = t76 * t12 + t74 * t13;
t6 = [0, 0, 0, 0, 0, qJDD(1), t101, t95, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t91, -t109 - t117, 0, pkin(1) * (t75 * t109 + t77 * t97), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - 0.2e1 * t73 - t91, 0.2e1 * t102 + t93 + t117, pkin(1) * (t75 * t27 - t77 * t28) - pkin(2) * t28 + qJ(3) * t27, t69 * qJDD(1), t104 * t118, 0, t68 * qJDD(1), 0, 0, t96 * t52 + t89 * t74, t96 * t51 + t89 * t76, -t122 * t99 + t96 * t55 - t10, -t96 * t10 + t99 * t26, -t74 * (t81 * t106 + t79 * t36) + t76 * (-t79 * t106 + t81 * t36), -t74 * (-t79 * t33 + t81 * t35) + t76 * (-t81 * t33 - t79 * t35), -t74 * (t81 * t39 + t125) + t76 * (-t79 * t39 + t124), -t74 * (t79 * t107 + t81 * t34) + t76 * (t81 * t107 - t79 * t34), -t74 * (t79 * t38 + t110) + t76 * (t81 * t38 - t112), (-t74 * (-t45 * t79 - t47 * t81) + t76 * (-t45 * t81 + t47 * t79)) * qJD(5), -t74 * (pkin(6) * t13 - t111) + t76 * (-pkin(6) * t12 + t113) + t92 * t33 - t96 * t4, -t74 * (pkin(6) * t19 + t113) + t76 * (-pkin(6) * t18 + t111) + t92 * t35 - t96 * t9, -t74 * (pkin(6) * t15 + t3) + t76 * (-pkin(6) * t14 - t2) + t92 * t23 - t96 * t5, t92 * t17 + (-pkin(6) - t96) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t20 + t76 * t21, 0, 0, 0, 0, 0, 0, -t74 * t12 + t76 * t13, -t74 * t18 + t76 * t19, -t74 * t14 + t76 * t15, -t74 * t2 + t76 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t84, t28, 0, 0, 0, 0, 0, 0, -t52, -t51, -t55, t10, 0, 0, 0, 0, 0, 0, t4, t9, t5, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t104, -t122, t26, 0, 0, 0, 0, 0, 0, t33, t35, t23, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, t43 - t42, t44, -t115, -t22, qJDD(5), -t7, -t8, 0, 0;];
tauJ_reg = t6;
