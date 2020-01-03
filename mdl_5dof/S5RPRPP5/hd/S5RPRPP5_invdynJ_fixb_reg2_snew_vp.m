% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:42
% EndTime: 2019-12-31 18:16:44
% DurationCPUTime: 0.78s
% Computational Cost: add. (944->172), mult. (1989->153), div. (0->0), fcn. (833->4), ass. (0->106)
t130 = -pkin(6) - pkin(1);
t70 = sin(qJ(3));
t72 = cos(qJ(3));
t75 = qJD(1) ^ 2;
t55 = t70 * t75 * t72;
t45 = qJDD(3) + t55;
t121 = t70 * t45;
t69 = t72 ^ 2;
t123 = t69 * t75;
t74 = qJD(3) ^ 2;
t52 = t74 + t123;
t17 = t72 * t52 + t121;
t106 = t72 * qJDD(1);
t104 = qJD(1) * qJD(3);
t99 = t70 * t104;
t37 = -0.2e1 * t99 + t106;
t95 = -qJ(2) * t37 + t130 * t17;
t68 = t70 ^ 2;
t117 = t68 + t69;
t41 = t117 * qJDD(1);
t42 = t117 * t75;
t94 = -qJ(2) * t42 - t130 * t41;
t124 = t68 * t75;
t132 = 2 * qJD(4);
t107 = t70 * qJDD(1);
t98 = t72 * t104;
t35 = t98 + t107;
t111 = qJD(1) * t72;
t44 = -qJD(3) * pkin(4) - qJ(5) * t111;
t36 = -t99 + t106;
t116 = qJ(4) * t36;
t128 = pkin(3) * t35;
t135 = -pkin(3) * t98 - qJ(4) * t99;
t67 = qJDD(1) * qJ(2);
t71 = sin(qJ(1));
t73 = cos(qJ(1));
t93 = g(1) * t73 + g(2) * t71;
t88 = -t67 + t93;
t81 = -t116 - t88 + t128 - t135;
t140 = -t35 * pkin(4) - qJ(5) * t124 + qJDD(5) - t81 + (-(2 * qJD(2)) + (t132 + t44) * t72) * qJD(1);
t139 = 0.2e1 * qJD(1);
t131 = pkin(3) + pkin(4);
t138 = t35 * qJ(5) + qJD(3) * t44;
t137 = pkin(3) * t52 + qJ(4) * t45;
t46 = qJDD(3) - t55;
t51 = -t74 - t124;
t136 = pkin(3) * t46 + qJ(4) * t51;
t129 = pkin(1) * t75;
t134 = -pkin(6) * t75 - t129 - t140;
t102 = t130 * t75;
t100 = g(1) * t71 - t73 * g(2);
t92 = qJDD(2) - t100;
t87 = -qJ(2) * t75 + t92;
t24 = t130 * qJDD(1) + t87;
t22 = t70 * t24;
t133 = qJDD(3) * qJ(4) + qJD(3) * t132 + t22;
t127 = g(3) * t70;
t126 = g(3) * t72;
t125 = t24 * t72;
t122 = t70 * t37;
t120 = t72 * t46;
t118 = t42 - t74;
t115 = qJ(4) * t70;
t114 = qJ(4) * t74;
t113 = t72 * qJ(4);
t91 = t70 * pkin(3) - t113;
t31 = t91 * qJD(1);
t112 = qJD(1) * t31;
t110 = qJD(4) * t72;
t109 = qJDD(1) * pkin(1);
t108 = qJDD(3) * pkin(3);
t105 = -t31 * t111 - qJDD(4);
t103 = qJD(2) * qJD(1);
t101 = qJD(5) * t139;
t97 = t131 * qJDD(3);
t16 = t70 * t51 + t120;
t34 = 0.2e1 * t98 + t107;
t96 = qJ(2) * t34 + t130 * t16;
t85 = -pkin(3) * t74 - t70 * t112 + t133;
t8 = t85 - t126;
t80 = t70 * t101 + t138 + t8;
t4 = -pkin(4) * t124 + t80;
t83 = -pkin(4) * t55 + (t24 + t101) * t72 + t105;
t90 = t36 + t99;
t76 = t90 * qJ(5) + t114 + t127 + t83;
t5 = t97 + t76;
t1 = t4 * t70 + t5 * t72;
t11 = t125 + t127;
t12 = -t22 + t126;
t7 = t11 * t72 - t12 * t70;
t10 = t72 * t34 + t122;
t89 = -(-t74 + t124) * t72 + t121;
t21 = (-t74 + t123) * t70 + t120;
t84 = -t102 + t88;
t82 = -t105 - t108 - t114 - t125;
t9 = -t82 + t127;
t78 = (-qJD(2) + t110) * t139 + t84 + t135;
t64 = 0.2e1 * t103;
t43 = (-t68 + t69) * t75;
t25 = -t87 + t109;
t23 = t84 - 0.2e1 * t103;
t20 = (t36 - t99) * t72;
t15 = (t35 + t98) * t70;
t3 = -t102 + t140;
t2 = t70 * t8 + t72 * t9;
t6 = [0, 0, 0, 0, 0, qJDD(1), t100, t93, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t92 - 0.2e1 * t109, t64 + 0.2e1 * t67 - t93, pkin(1) * t25 + qJ(2) * (t64 - t88 - t129), t20, -t10, t21, t15, -t89, 0, -t23 * t70 + t96, -t23 * t72 - t95, -t7 + t94, -qJ(2) * t23 + t130 * t7, t20, t21, t10, 0, t89, t15, -t34 * t113 - t70 * (t116 + (-t34 - t35) * pkin(3) + t78) + t96, t72 * (qJ(4) * t42 + t82) - t70 * (pkin(3) * t42 + t85) + t94, t72 * (-t128 + (t36 + t37) * qJ(4) + t78) - pkin(3) * t122 + t95, t130 * t2 + (qJ(2) + t91) * (-0.2e1 * qJD(1) * t110 + t102 + t64 + t81), t20, t10, -t21, t15, -t89, 0, t72 * (-qJ(4) * t34 + qJ(5) * t46) + t96 + (qJ(5) * t51 + t131 * t34 + t134) * t70, -t70 * (-qJ(5) * t45 + t131 * t37) + t95 + (qJ(4) * t37 + qJ(5) * t52 - t134) * t72, ((qJ(5) * qJDD(1) + t101 - t112) * t70 + (t42 - t124) * pkin(4) + t118 * pkin(3) + t133 + t138) * t70 + (t97 + (t90 + t106) * qJ(5) - t118 * qJ(4) + t83) * t72 - t94, (-t70 * t131 - qJ(2) + t113) * t3 + (qJ(5) + t130) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t75, -t25, 0, 0, 0, 0, 0, 0, t16, -t17, -t41, t7, 0, 0, 0, 0, 0, 0, t16, -t41, t17, t2, 0, 0, 0, 0, 0, 0, t16, t17, t41, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t43, t106, -t55, -t107, qJDD(3), t11, t12, 0, 0, t55, t106, -t43, qJDD(3), t107, -t55, t9 + t136, (-pkin(3) * t72 - t115) * qJDD(1), t8 + t137, pkin(3) * t9 + qJ(4) * t8, t55, -t43, -t106, -t55, -t107, qJDD(3), t108 + (qJDD(3) + t46) * pkin(4) + t76 + t136, (t52 - t124) * pkin(4) + t80 + t137, (t72 * t131 + t115) * qJDD(1), qJ(4) * t4 + t131 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t106, -t52, -t9, 0, 0, 0, 0, 0, 0, -t46, -t52, -t106, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t37, -t42, t3;];
tauJ_reg = t6;
