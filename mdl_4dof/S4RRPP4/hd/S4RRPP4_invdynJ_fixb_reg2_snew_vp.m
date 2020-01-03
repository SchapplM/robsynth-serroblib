% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRPP4
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRPP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:18
% EndTime: 2019-12-31 16:59:20
% DurationCPUTime: 0.52s
% Computational Cost: add. (643->148), mult. (1480->146), div. (0->0), fcn. (645->4), ass. (0->94)
t69 = sin(qJ(2));
t71 = cos(qJ(2));
t74 = qJD(1) ^ 2;
t100 = t71 * t74 * t69;
t42 = qJDD(2) + t100;
t136 = pkin(3) * t42;
t105 = t71 * qJDD(1);
t103 = qJD(1) * qJD(2);
t57 = t69 * t103;
t32 = -t57 + t105;
t135 = t32 - t57;
t43 = qJDD(2) - t100;
t116 = t71 * t43;
t106 = t69 * qJDD(1);
t95 = t71 * t103;
t30 = 0.2e1 * t95 + t106;
t65 = t69 ^ 2;
t123 = t65 * t74;
t73 = qJD(2) ^ 2;
t47 = t73 + t123;
t134 = pkin(1) * t30 + pkin(5) * (-t69 * t47 + t116);
t108 = qJD(1) * t69;
t41 = -qJD(2) * pkin(3) - qJ(4) * t108;
t133 = -t32 * qJ(4) + qJD(2) * t41;
t66 = t71 ^ 2;
t122 = t66 * t74;
t50 = -t73 - t122;
t132 = pkin(2) * t42 + qJ(3) * t50;
t107 = qJD(2) * t71;
t31 = t95 + t106;
t70 = sin(qJ(1));
t72 = cos(qJ(1));
t114 = t70 * g(1) - t72 * g(2);
t19 = qJDD(1) * pkin(1) + t74 * pkin(5) + t114;
t81 = t135 * pkin(2) + t19;
t77 = t32 * pkin(3) - qJ(4) * t122 + qJDD(4) + t81;
t128 = 2 * qJD(3);
t94 = (t128 + t41) * t69;
t1 = t31 * qJ(3) + t77 + (qJ(3) * t107 + t94) * qJD(1);
t131 = qJDD(1) * pkin(5);
t13 = t116 + t69 * (-t73 + t122);
t127 = pkin(2) + pkin(3);
t129 = -t71 * t127 - pkin(1);
t125 = t71 * g(3);
t110 = t69 * qJ(3);
t89 = -t71 * pkin(2) - t110;
t28 = t89 * qJD(1);
t109 = qJD(1) * t28;
t91 = t72 * g(1) + t70 * g(2);
t20 = -t74 * pkin(1) + t131 - t91;
t92 = t20 + t109;
t97 = qJDD(2) * pkin(2) + t73 * qJ(3) - qJDD(3);
t6 = t92 * t69 + t125 - t97;
t33 = -0.2e1 * t57 + t105;
t126 = pkin(2) * t33;
t113 = t65 + t66;
t29 = t113 * t131;
t120 = t69 * t42;
t124 = pkin(1) * t33 + pkin(5) * (t71 * t50 - t120);
t121 = t69 * t33;
t117 = t71 * t30;
t16 = -t69 * g(3) + t71 * t20;
t39 = t113 * t74;
t115 = pkin(1) * t39 + t29;
t112 = qJ(3) * t39;
t111 = qJ(3) * t71;
t104 = qJ(4) * qJDD(1);
t102 = qJD(1) * qJD(4);
t101 = qJD(3) * qJD(2);
t99 = 0.2e1 * t102;
t98 = qJ(4) * t107;
t15 = t69 * t20 + t125;
t96 = t69 * t15 + t71 * t16;
t93 = qJDD(2) * qJ(3) + t71 * t109 + t16;
t88 = t31 + t95;
t87 = t117 + t121;
t14 = -t71 * (-t73 + t123) + t120;
t85 = -t31 * qJ(4) - t136 - t97;
t84 = t73 * pkin(2) - t93;
t59 = 0.2e1 * t101;
t5 = t59 - t84;
t83 = pkin(2) * t47 + qJ(3) * t43 + t5;
t80 = t15 + t85;
t79 = (t30 + t88) * t110 + t134;
t78 = pkin(3) * t122 - t133 + t84;
t76 = t108 * t128 + t81;
t55 = -0.2e1 * t71 * t102;
t54 = t69 * t99;
t40 = (t65 - t66) * t74;
t12 = t88 * t69;
t11 = t135 * t71;
t3 = (t98 + (-0.2e1 * qJD(4) + t28) * t69) * qJD(1) + t80;
t2 = t55 + t59 - t78;
t4 = [0, 0, 0, 0, 0, qJDD(1), t114, t91, 0, 0, t12, t87, t14, t11, t13, 0, t71 * t19 + t124, -t69 * t19 - t134, t96 + t115, pkin(1) * t19 + pkin(5) * t96, t12, t14, -t87, 0, -t13, t11, t71 * (t76 + t126) + (t71 * t88 + t121) * qJ(3) + t124, t71 * (t59 + (t39 - t73) * pkin(2) + t93) + (t6 + t112) * t69 + t115, pkin(2) * t117 + t69 * t76 + t79, pkin(5) * (t71 * t5 + t69 * t6) + (pkin(1) - t89) * (t88 * qJ(3) + t76), t12, -t87, -t14, t11, t13, 0, t69 * (qJ(3) * t33 + qJ(4) * t42) + t124 + (pkin(3) * t33 - qJ(4) * t50 + t1 + t126) * t71, t71 * (-qJ(4) * t43 + t127 * t30) + (qJ(4) * t47 + qJD(1) * t94 + t77) * t69 + t79, -t29 + (-0.2e1 * t101 + (t99 + t104) * t71 + t78) * t71 + t129 * t39 + (-t112 + t54 + (-t92 + t104) * t69 + (-qJ(4) * t103 - g(3)) * t71 - t85) * t69, (t110 - t129) * t1 + (pkin(5) - qJ(4)) * (t71 * t2 + t69 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, t40, t106, t100, t105, qJDD(2), -t15, -t16, 0, 0, -t100, t106, -t40, qJDD(2), -t105, t100, -t6 + t132, (-pkin(2) * t69 + t111) * qJDD(1), t83, -pkin(2) * t6 + qJ(3) * t5, -t100, -t40, -t106, t100, t105, qJDD(2), t136 + t54 + (-t28 * t69 - t98) * qJD(1) - t80 + t132, t55 + (t47 - t122) * pkin(3) + t83 + t133, (t127 * t69 - t111) * qJDD(1), qJ(3) * t2 - t127 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t106, -t47, t6, 0, 0, 0, 0, 0, 0, -t42, -t47, -t106, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t30, -t39, t1;];
tauJ_reg = t4;
