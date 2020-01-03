% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:01
% EndTime: 2019-12-31 17:44:04
% DurationCPUTime: 0.90s
% Computational Cost: add. (1648->158), mult. (3624->215), div. (0->0), fcn. (2181->8), ass. (0->101)
t85 = sin(pkin(8));
t87 = cos(pkin(8));
t89 = sin(qJ(5));
t90 = cos(qJ(5));
t103 = t85 * t89 + t87 * t90;
t51 = t103 * qJD(1);
t119 = qJD(1) * t87;
t120 = qJD(1) * t85;
t53 = -t89 * t119 + t90 * t120;
t128 = t53 * t51;
t139 = qJDD(5) - t128;
t141 = t139 * t89;
t140 = t139 * t90;
t133 = 2 * qJD(3);
t129 = sin(qJ(1));
t130 = cos(qJ(1));
t100 = t130 * g(1) + t129 * g(2);
t92 = qJD(1) ^ 2;
t63 = -t92 * pkin(1) - t100;
t86 = sin(pkin(7));
t88 = cos(pkin(7));
t99 = t129 * g(1) - t130 * g(2);
t96 = qJDD(1) * pkin(1) + t99;
t122 = t88 * t63 + t86 * t96;
t30 = -t92 * pkin(2) + qJDD(1) * qJ(3) + t122;
t105 = qJD(1) * t133 + t30;
t127 = t87 * t92;
t80 = t85 ^ 2;
t81 = t87 ^ 2;
t138 = t80 + t81;
t107 = pkin(1) * t86 + qJ(3);
t65 = t138 * t92;
t137 = t107 * t85 * t65;
t108 = qJD(4) * t120;
t78 = t85 * qJDD(1);
t73 = qJ(4) * t78;
t136 = 0.2e1 * t73 + 0.2e1 * t108;
t121 = t85 * qJ(4);
t131 = t87 * pkin(3);
t61 = (-t121 - t131) * qJD(1);
t102 = t30 + (t133 + t61) * qJD(1);
t83 = -g(3) + qJDD(2);
t72 = t87 * t83;
t114 = qJDD(4) - t72;
t135 = (-pkin(4) * t127 - pkin(6) * qJDD(1) + t102) * t85 + t114;
t26 = t103 * qJDD(1);
t110 = pkin(1) * t88 + pkin(2);
t101 = t110 + t121;
t97 = t87 * (pkin(3) + pkin(4)) + t101;
t48 = t51 ^ 2;
t49 = t53 ^ 2;
t116 = t87 * qJDD(1);
t106 = -t86 * t63 + t88 * t96;
t29 = -qJDD(1) * pkin(2) - t92 * qJ(3) + qJDD(3) - t106;
t75 = pkin(3) * t116;
t98 = -t29 + t75;
t23 = -t73 - t98 - 0.2e1 * t108;
t15 = -pkin(4) * t116 + pkin(6) * t65 + t23;
t126 = t89 * t15;
t33 = qJDD(5) + t128;
t125 = t89 * t33;
t124 = t90 * t15;
t123 = t90 * t33;
t118 = t51 * qJD(5);
t117 = t53 * qJD(5);
t115 = t88 * qJDD(1);
t25 = t105 * t87 + t85 * t83;
t77 = t80 * qJDD(1);
t79 = t81 * qJDD(1);
t64 = t79 + t77;
t113 = pkin(1) * (t86 * t64 + t88 * t65) + qJ(3) * t64 + pkin(2) * t65;
t59 = t138 * t127;
t112 = pkin(2) * t116 - qJ(3) * t59 + pkin(1) * (t87 * t115 - t86 * t59);
t109 = t85 * t116;
t24 = t105 * t85 - t72;
t10 = t85 * t24 + t87 * t25;
t19 = t61 * t119 + t25;
t12 = -t81 * t92 * pkin(4) - pkin(6) * t116 + t19;
t4 = t89 * t12 - t90 * t135;
t5 = t90 * t12 + t135 * t89;
t2 = -t90 * t4 + t89 * t5;
t3 = t89 * t4 + t90 * t5;
t50 = -t89 * t116 + t90 * t78;
t91 = qJD(5) ^ 2;
t44 = -t49 - t91;
t43 = -t49 + t91;
t42 = t48 - t91;
t38 = t50 - t118;
t37 = t50 - 0.2e1 * t118;
t36 = -t26 - t117;
t35 = 0.2e1 * t117 + t26;
t31 = -t91 - t48;
t27 = -t48 - t49;
t21 = -t89 * t44 - t123;
t20 = t90 * t44 - t125;
t18 = t102 * t85 + t114;
t17 = -t90 * t26 + t89 * t50;
t16 = -t89 * t26 - t90 * t50;
t14 = t90 * t31 - t141;
t13 = t89 * t31 + t140;
t1 = [0, 0, 0, 0, 0, qJDD(1), t99, t100, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (-t86 * t92 + t115) + t106, pkin(1) * (-t86 * qJDD(1) - t88 * t92) - t122, 0, pkin(1) * (t88 * t106 + t86 * t122), t77, 0.2e1 * t109, 0, t79, 0, 0, -t87 * t29 + t112, t137 + (-t110 * qJDD(1) + t29) * t85, t10 + t113, -pkin(2) * t29 + qJ(3) * t10 + pkin(1) * (t86 * t10 - t88 * t29), t77, 0, -0.2e1 * t109, 0, 0, t79, (-t29 + 0.2e1 * t75 + t136) * t87 + t112, t87 * (pkin(3) * t65 + t19) + (qJ(4) * t65 + (qJD(1) * t61 + t105) * t85 + t114) * t85 + t113, -t137 + ((t110 + t131) * qJDD(1) + t98 + t136) * t85, t107 * (t85 * t18 + t87 * t19) + (-t101 - t131) * t23, t85 * (-t89 * t117 + t90 * t38) + t87 * (-t90 * t117 - t89 * t38), t85 * (-t90 * t35 - t89 * t37) + t87 * (t89 * t35 - t90 * t37), t85 * (-t89 * t43 + t140) + t87 * (-t90 * t43 - t141), t85 * (t90 * t118 - t89 * t36) + t87 * (-t89 * t118 - t90 * t36), t85 * (t90 * t42 - t125) + t87 * (-t89 * t42 - t123), (t85 * (-t51 * t90 + t53 * t89) + t87 * (t51 * t89 + t53 * t90)) * qJD(5), t85 * (-pkin(6) * t13 - t126) + t87 * (-pkin(6) * t14 - t124) + t107 * (t85 * t13 + t87 * t14) + t97 * t35, t85 * (-pkin(6) * t20 - t124) + t87 * (-pkin(6) * t21 + t126) + t107 * (t85 * t20 + t87 * t21) + t97 * t37, t85 * (-pkin(6) * t16 - t2) + t87 * (-pkin(6) * t17 - t3) + t107 * (t85 * t16 + t87 * t17) + t97 * t27, -t97 * t15 + (t107 - pkin(6)) * (t85 * t2 + t87 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87 * t24 + t85 * t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87 * t18 + t85 * t19, 0, 0, 0, 0, 0, 0, -t87 * t13 + t85 * t14, -t87 * t20 + t85 * t21, -t87 * t16 + t85 * t17, -t87 * t2 + t85 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, t78, -t65, t29, 0, 0, 0, 0, 0, 0, -t116, -t65, -t78, t23, 0, 0, 0, 0, 0, 0, -t35, -t37, -t27, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85 * t127, t78, -t80 * t92, t18, 0, 0, 0, 0, 0, 0, t13, t20, t16, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t49 - t48, t50, -t128, -t26, qJDD(5), -t4, -t5, 0, 0;];
tauJ_reg = t1;
