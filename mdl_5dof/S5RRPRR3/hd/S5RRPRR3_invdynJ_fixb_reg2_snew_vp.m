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
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:34:21
% EndTime: 2022-01-20 10:34:24
% DurationCPUTime: 0.75s
% Computational Cost: add. (4357->125), mult. (5715->187), div. (0->0), fcn. (3326->10), ass. (0->99)
t109 = sin(pkin(9));
t110 = cos(pkin(9));
t104 = qJDD(1) + qJDD(2);
t100 = qJDD(4) + t104;
t112 = sin(qJ(4));
t116 = cos(qJ(4));
t105 = qJD(1) + qJD(2);
t101 = qJD(4) + t105;
t99 = t101 ^ 2;
t124 = -t112 * t100 - t116 * t99;
t76 = t116 * t100 - t112 * t99;
t47 = -t109 * t76 + t110 * t124;
t48 = t109 * t124 + t110 * t76;
t111 = sin(qJ(5));
t115 = cos(qJ(5));
t108 = -g(3) + qJDD(3);
t113 = sin(qJ(2));
t117 = cos(qJ(2));
t114 = sin(qJ(1));
t118 = cos(qJ(1));
t130 = t114 * g(1) - t118 * g(2);
t87 = qJDD(1) * pkin(1) + t130;
t125 = t118 * g(1) + t114 * g(2);
t88 = -qJD(1) ^ 2 * pkin(1) - t125;
t60 = -t113 * t88 + t117 * t87;
t54 = t104 * pkin(2) + t60;
t103 = t105 ^ 2;
t61 = t113 * t87 + t117 * t88;
t55 = -t103 * pkin(2) + t61;
t35 = -t109 * t55 + t110 * t54;
t33 = t104 * pkin(3) + t35;
t36 = t109 * t54 + t110 * t55;
t34 = -t103 * pkin(3) + t36;
t23 = t112 * t33 + t116 * t34;
t21 = -t99 * pkin(4) + t100 * pkin(8) + t23;
t15 = -t115 * t108 + t111 * t21;
t16 = t111 * t108 + t115 * t21;
t9 = t111 * t15 + t115 * t16;
t22 = -t112 * t34 + t116 * t33;
t20 = -t100 * pkin(4) - t99 * pkin(8) - t22;
t140 = -pkin(4) * t20 + pkin(8) * t9;
t11 = t112 * t23 + t116 * t22;
t12 = -t112 * t22 + t116 * t23;
t4 = t109 * t12 + t110 * t11;
t139 = pkin(2) * t4 + pkin(3) * t11;
t91 = t115 * t99 * t111;
t85 = qJDD(5) + t91;
t138 = t111 * t85;
t86 = qJDD(5) - t91;
t137 = t115 * t86;
t136 = t111 * t100;
t135 = qJD(5) * t101;
t119 = qJD(5) ^ 2;
t106 = t111 ^ 2;
t94 = t106 * t99;
t89 = -t94 - t119;
t65 = -t111 * t89 - t137;
t71 = 0.2e1 * t115 * t135 + t136;
t134 = -pkin(4) * t71 + pkin(8) * t65 + t111 * t20;
t107 = t115 ^ 2;
t95 = t107 * t99;
t90 = -t95 - t119;
t64 = t115 * t90 - t138;
t93 = t115 * t100;
t72 = -0.2e1 * t111 * t135 + t93;
t133 = pkin(4) * t72 + pkin(8) * t64 - t115 * t20;
t83 = -t110 * t103 - t109 * t104;
t132 = pkin(2) * t83 - t36;
t6 = t112 * t9 - t116 * t20;
t7 = t112 * t20 + t116 * t9;
t2 = t109 * t7 + t110 * t6;
t131 = pkin(2) * t2 + pkin(3) * t6 + t140;
t74 = (t106 + t107) * t100;
t79 = t94 + t95;
t129 = pkin(4) * t79 + pkin(8) * t74 + t9;
t121 = t109 * t103 - t110 * t104;
t128 = -pkin(2) * t121 + t35;
t39 = t112 * t64 + t116 * t72;
t41 = -t112 * t72 + t116 * t64;
t28 = t109 * t41 + t110 * t39;
t127 = pkin(2) * t28 + pkin(3) * t39 + t133;
t40 = t112 * t65 - t116 * t71;
t42 = t112 * t71 + t116 * t65;
t29 = t109 * t42 + t110 * t40;
t126 = pkin(2) * t29 + pkin(3) * t40 + t134;
t123 = pkin(2) * t48 + pkin(3) * t76 + t22;
t49 = t112 * t74 + t116 * t79;
t50 = -t112 * t79 + t116 * t74;
t31 = t109 * t50 + t110 * t49;
t122 = pkin(2) * t31 + pkin(3) * t49 + t129;
t120 = pkin(2) * t47 + pkin(3) * t124 - t23;
t63 = t138 + t115 * (-t94 + t119);
t62 = t111 * (t95 - t119) + t137;
t57 = t71 * t111;
t56 = t72 * t115;
t43 = t111 * t72 + t115 * t71;
t25 = t109 * t36 + t110 * t35;
t24 = pkin(2) * t25;
t1 = [0, 0, 0, 0, 0, qJDD(1), t130, t125, 0, 0, 0, 0, 0, 0, 0, t104, pkin(1) * (-t113 * t103 + t117 * t104) + t60, pkin(1) * (-t117 * t103 - t113 * t104) - t61, 0, pkin(1) * (t113 * t61 + t117 * t60), 0, 0, 0, 0, 0, t104, pkin(1) * (t113 * t83 - t117 * t121) + t128, pkin(1) * (t113 * t121 + t117 * t83) + t132, 0, pkin(1) * (t113 * (-t109 * t35 + t110 * t36) + t117 * t25) + t24, 0, 0, 0, 0, 0, t100, pkin(1) * (t113 * t47 + t117 * t48) + t123, pkin(1) * (-t113 * t48 + t117 * t47) + t120, 0, pkin(1) * (t113 * (-t109 * t11 + t110 * t12) + t117 * t4) + t139, t57, t43, t63, t56, t62, 0, pkin(1) * (t113 * (-t109 * t39 + t110 * t41) + t117 * t28) + t127, pkin(1) * (t113 * (-t109 * t40 + t110 * t42) + t117 * t29) + t126, pkin(1) * (t113 * (-t109 * t49 + t110 * t50) + t117 * t31) + t122, pkin(1) * (t113 * (-t109 * t6 + t110 * t7) + t117 * t2) + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t60, -t61, 0, 0, 0, 0, 0, 0, 0, t104, t128, t132, 0, t24, 0, 0, 0, 0, 0, t100, t123, t120, 0, t139, t57, t43, t63, t56, t62, 0, t127, t126, t122, t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, 0, 0, 0, 0, 0, 0, t111 * t90 + t115 * t85, -t111 * t86 + t115 * t89, 0, t111 * t16 - t115 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t22, -t23, 0, 0, t57, t43, t63, t56, t62, 0, t133, t134, t129, t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t94 - t95, t136, t91, t93, qJDD(5), -t15, -t16, 0, 0;];
tauJ_reg = t1;
