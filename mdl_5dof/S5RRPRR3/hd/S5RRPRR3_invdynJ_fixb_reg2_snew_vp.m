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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:30:32
% EndTime: 2019-12-05 18:30:35
% DurationCPUTime: 0.70s
% Computational Cost: add. (4357->125), mult. (5715->187), div. (0->0), fcn. (3326->10), ass. (0->99)
t110 = sin(pkin(9));
t111 = cos(pkin(9));
t105 = qJDD(1) + qJDD(2);
t100 = qJDD(4) + t105;
t113 = sin(qJ(4));
t117 = cos(qJ(4));
t106 = qJD(1) + qJD(2);
t101 = qJD(4) + t106;
t99 = t101 ^ 2;
t125 = -t113 * t100 - t117 * t99;
t76 = t117 * t100 - t113 * t99;
t47 = -t110 * t76 + t111 * t125;
t48 = t110 * t125 + t111 * t76;
t112 = sin(qJ(5));
t116 = cos(qJ(5));
t109 = -g(1) + qJDD(3);
t114 = sin(qJ(2));
t118 = cos(qJ(2));
t115 = sin(qJ(1));
t119 = cos(qJ(1));
t136 = t119 * g(2) + t115 * g(3);
t87 = qJDD(1) * pkin(1) + t136;
t126 = -t115 * g(2) + t119 * g(3);
t88 = -qJD(1) ^ 2 * pkin(1) - t126;
t60 = -t114 * t88 + t118 * t87;
t54 = t105 * pkin(2) + t60;
t104 = t106 ^ 2;
t61 = t114 * t87 + t118 * t88;
t55 = -t104 * pkin(2) + t61;
t35 = -t110 * t55 + t111 * t54;
t33 = t105 * pkin(3) + t35;
t36 = t110 * t54 + t111 * t55;
t34 = -t104 * pkin(3) + t36;
t23 = t113 * t33 + t117 * t34;
t21 = -t99 * pkin(4) + t100 * pkin(8) + t23;
t15 = -t116 * t109 + t112 * t21;
t16 = t112 * t109 + t116 * t21;
t9 = t112 * t15 + t116 * t16;
t22 = -t113 * t34 + t117 * t33;
t20 = -t100 * pkin(4) - t99 * pkin(8) - t22;
t141 = -pkin(4) * t20 + pkin(8) * t9;
t11 = t113 * t23 + t117 * t22;
t12 = -t113 * t22 + t117 * t23;
t4 = t111 * t11 + t110 * t12;
t140 = pkin(2) * t4 + pkin(3) * t11;
t91 = t116 * t99 * t112;
t85 = qJDD(5) + t91;
t139 = t112 * t85;
t86 = qJDD(5) - t91;
t138 = t116 * t86;
t137 = t112 * t100;
t135 = qJD(5) * t101;
t120 = qJD(5) ^ 2;
t107 = t112 ^ 2;
t94 = t107 * t99;
t89 = -t94 - t120;
t65 = -t112 * t89 - t138;
t71 = 0.2e1 * t116 * t135 + t137;
t134 = -pkin(4) * t71 + pkin(8) * t65 + t112 * t20;
t108 = t116 ^ 2;
t95 = t108 * t99;
t90 = -t95 - t120;
t64 = t116 * t90 - t139;
t93 = t116 * t100;
t72 = -0.2e1 * t112 * t135 + t93;
t133 = pkin(4) * t72 + pkin(8) * t64 - t116 * t20;
t83 = -t111 * t104 - t110 * t105;
t132 = pkin(2) * t83 - t36;
t6 = t113 * t9 - t117 * t20;
t7 = t113 * t20 + t117 * t9;
t2 = t110 * t7 + t111 * t6;
t131 = pkin(2) * t2 + pkin(3) * t6 + t141;
t74 = (t107 + t108) * t100;
t79 = t94 + t95;
t130 = pkin(4) * t79 + pkin(8) * t74 + t9;
t122 = t110 * t104 - t111 * t105;
t129 = -pkin(2) * t122 + t35;
t39 = t113 * t64 + t117 * t72;
t41 = -t113 * t72 + t117 * t64;
t28 = t110 * t41 + t111 * t39;
t128 = pkin(2) * t28 + pkin(3) * t39 + t133;
t40 = t113 * t65 - t117 * t71;
t42 = t113 * t71 + t117 * t65;
t29 = t110 * t42 + t111 * t40;
t127 = pkin(2) * t29 + pkin(3) * t40 + t134;
t124 = pkin(2) * t48 + pkin(3) * t76 + t22;
t49 = t113 * t74 + t117 * t79;
t50 = -t113 * t79 + t117 * t74;
t31 = t110 * t50 + t111 * t49;
t123 = pkin(2) * t31 + pkin(3) * t49 + t130;
t121 = pkin(2) * t47 + pkin(3) * t125 - t23;
t63 = t139 + t116 * (-t94 + t120);
t62 = t112 * (t95 - t120) + t138;
t57 = t71 * t112;
t56 = t72 * t116;
t43 = t112 * t72 + t116 * t71;
t25 = t110 * t36 + t111 * t35;
t24 = pkin(2) * t25;
t1 = [0, 0, 0, 0, 0, qJDD(1), t136, t126, 0, 0, 0, 0, 0, 0, 0, t105, pkin(1) * (-t114 * t104 + t118 * t105) + t60, pkin(1) * (-t118 * t104 - t114 * t105) - t61, 0, pkin(1) * (t114 * t61 + t118 * t60), 0, 0, 0, 0, 0, t105, pkin(1) * (t114 * t83 - t118 * t122) + t129, pkin(1) * (t114 * t122 + t118 * t83) + t132, 0, pkin(1) * (t114 * (-t110 * t35 + t111 * t36) + t118 * t25) + t24, 0, 0, 0, 0, 0, t100, pkin(1) * (t114 * t47 + t118 * t48) + t124, pkin(1) * (-t114 * t48 + t118 * t47) + t121, 0, pkin(1) * (t114 * (-t110 * t11 + t111 * t12) + t118 * t4) + t140, t57, t43, t63, t56, t62, 0, pkin(1) * (t114 * (-t110 * t39 + t111 * t41) + t118 * t28) + t128, pkin(1) * (t114 * (-t110 * t40 + t111 * t42) + t118 * t29) + t127, pkin(1) * (t114 * (-t110 * t49 + t111 * t50) + t118 * t31) + t123, pkin(1) * (t114 * (-t110 * t6 + t111 * t7) + t118 * t2) + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t60, -t61, 0, 0, 0, 0, 0, 0, 0, t105, t129, t132, 0, t24, 0, 0, 0, 0, 0, t100, t124, t121, 0, t140, t57, t43, t63, t56, t62, 0, t128, t127, t123, t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, 0, 0, 0, 0, 0, 0, t112 * t90 + t116 * t85, -t112 * t86 + t116 * t89, 0, t112 * t16 - t116 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t22, -t23, 0, 0, t57, t43, t63, t56, t62, 0, t133, t134, t130, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t94 - t95, t137, t91, t93, qJDD(5), -t15, -t16, 0, 0;];
tauJ_reg = t1;
