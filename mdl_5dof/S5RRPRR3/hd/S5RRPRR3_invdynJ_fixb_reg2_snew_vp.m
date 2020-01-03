% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:36
% EndTime: 2020-01-03 12:00:38
% DurationCPUTime: 0.70s
% Computational Cost: add. (4357->125), mult. (5715->187), div. (0->0), fcn. (3326->10), ass. (0->99)
t108 = sin(pkin(9));
t109 = cos(pkin(9));
t103 = qJDD(1) + qJDD(2);
t100 = qJDD(4) + t103;
t111 = sin(qJ(4));
t115 = cos(qJ(4));
t104 = qJD(1) + qJD(2);
t101 = qJD(4) + t104;
t99 = t101 ^ 2;
t123 = -t111 * t100 - t115 * t99;
t76 = t115 * t100 - t111 * t99;
t47 = -t108 * t76 + t109 * t123;
t48 = t108 * t123 + t109 * t76;
t110 = sin(qJ(5));
t114 = cos(qJ(5));
t107 = -g(1) + qJDD(3);
t112 = sin(qJ(2));
t116 = cos(qJ(2));
t113 = sin(qJ(1));
t117 = cos(qJ(1));
t125 = -t117 * g(2) - t113 * g(3);
t87 = qJDD(1) * pkin(1) + t125;
t124 = t113 * g(2) - t117 * g(3);
t88 = -qJD(1) ^ 2 * pkin(1) - t124;
t60 = -t112 * t88 + t116 * t87;
t54 = t103 * pkin(2) + t60;
t102 = t104 ^ 2;
t61 = t112 * t87 + t116 * t88;
t55 = -t102 * pkin(2) + t61;
t35 = -t108 * t55 + t109 * t54;
t33 = t103 * pkin(3) + t35;
t36 = t108 * t54 + t109 * t55;
t34 = -t102 * pkin(3) + t36;
t23 = t111 * t33 + t115 * t34;
t21 = -t99 * pkin(4) + t100 * pkin(8) + t23;
t15 = -t114 * t107 + t110 * t21;
t16 = t110 * t107 + t114 * t21;
t9 = t110 * t15 + t114 * t16;
t22 = -t111 * t34 + t115 * t33;
t20 = -t100 * pkin(4) - t99 * pkin(8) - t22;
t139 = -pkin(4) * t20 + pkin(8) * t9;
t11 = t111 * t23 + t115 * t22;
t12 = -t111 * t22 + t115 * t23;
t4 = t108 * t12 + t109 * t11;
t138 = pkin(2) * t4 + pkin(3) * t11;
t91 = t114 * t99 * t110;
t85 = qJDD(5) + t91;
t137 = t110 * t85;
t86 = qJDD(5) - t91;
t136 = t114 * t86;
t135 = t110 * t100;
t134 = qJD(5) * t101;
t118 = qJD(5) ^ 2;
t105 = t110 ^ 2;
t94 = t105 * t99;
t89 = -t94 - t118;
t65 = -t110 * t89 - t136;
t71 = 0.2e1 * t114 * t134 + t135;
t133 = -pkin(4) * t71 + pkin(8) * t65 + t110 * t20;
t106 = t114 ^ 2;
t95 = t106 * t99;
t90 = -t95 - t118;
t64 = t114 * t90 - t137;
t93 = t114 * t100;
t72 = -0.2e1 * t110 * t134 + t93;
t132 = pkin(4) * t72 + pkin(8) * t64 - t114 * t20;
t83 = -t109 * t102 - t108 * t103;
t131 = pkin(2) * t83 - t36;
t6 = t111 * t9 - t115 * t20;
t7 = t111 * t20 + t115 * t9;
t2 = t108 * t7 + t109 * t6;
t130 = pkin(2) * t2 + pkin(3) * t6 + t139;
t74 = (t105 + t106) * t100;
t79 = t94 + t95;
t129 = pkin(4) * t79 + pkin(8) * t74 + t9;
t120 = t108 * t102 - t109 * t103;
t128 = -pkin(2) * t120 + t35;
t39 = t111 * t64 + t115 * t72;
t41 = -t111 * t72 + t115 * t64;
t28 = t108 * t41 + t109 * t39;
t127 = pkin(2) * t28 + pkin(3) * t39 + t132;
t40 = t111 * t65 - t115 * t71;
t42 = t111 * t71 + t115 * t65;
t29 = t108 * t42 + t109 * t40;
t126 = pkin(2) * t29 + pkin(3) * t40 + t133;
t122 = pkin(2) * t48 + pkin(3) * t76 + t22;
t49 = t111 * t74 + t115 * t79;
t50 = -t111 * t79 + t115 * t74;
t31 = t108 * t50 + t109 * t49;
t121 = pkin(2) * t31 + pkin(3) * t49 + t129;
t119 = pkin(2) * t47 + pkin(3) * t123 - t23;
t63 = t137 + t114 * (-t94 + t118);
t62 = t110 * (t95 - t118) + t136;
t57 = t71 * t110;
t56 = t72 * t114;
t43 = t110 * t72 + t114 * t71;
t25 = t108 * t36 + t109 * t35;
t24 = pkin(2) * t25;
t1 = [0, 0, 0, 0, 0, qJDD(1), t125, t124, 0, 0, 0, 0, 0, 0, 0, t103, pkin(1) * (-t112 * t102 + t116 * t103) + t60, pkin(1) * (-t116 * t102 - t112 * t103) - t61, 0, pkin(1) * (t112 * t61 + t116 * t60), 0, 0, 0, 0, 0, t103, pkin(1) * (t112 * t83 - t116 * t120) + t128, pkin(1) * (t112 * t120 + t116 * t83) + t131, 0, pkin(1) * (t112 * (-t108 * t35 + t109 * t36) + t116 * t25) + t24, 0, 0, 0, 0, 0, t100, pkin(1) * (t112 * t47 + t116 * t48) + t122, pkin(1) * (-t112 * t48 + t116 * t47) + t119, 0, pkin(1) * (t112 * (-t108 * t11 + t109 * t12) + t116 * t4) + t138, t57, t43, t63, t56, t62, 0, pkin(1) * (t112 * (-t108 * t39 + t109 * t41) + t116 * t28) + t127, pkin(1) * (t112 * (-t108 * t40 + t109 * t42) + t116 * t29) + t126, pkin(1) * (t112 * (-t108 * t49 + t109 * t50) + t116 * t31) + t121, pkin(1) * (t112 * (-t108 * t6 + t109 * t7) + t116 * t2) + t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, t60, -t61, 0, 0, 0, 0, 0, 0, 0, t103, t128, t131, 0, t24, 0, 0, 0, 0, 0, t100, t122, t119, 0, t138, t57, t43, t63, t56, t62, 0, t127, t126, t121, t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, 0, 0, 0, 0, 0, 0, t110 * t90 + t114 * t85, -t110 * t86 + t114 * t89, 0, t110 * t16 - t114 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t22, -t23, 0, 0, t57, t43, t63, t56, t62, 0, t132, t133, t129, t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t94 - t95, t135, t91, t93, qJDD(5), -t15, -t16, 0, 0;];
tauJ_reg = t1;
