% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:37
% EndTime: 2020-01-03 12:07:40
% DurationCPUTime: 0.72s
% Computational Cost: add. (4831->125), mult. (6090->187), div. (0->0), fcn. (3542->10), ass. (0->102)
t110 = sin(qJ(3));
t114 = cos(qJ(3));
t103 = qJDD(1) + qJDD(2);
t100 = qJDD(3) + t103;
t107 = sin(pkin(9));
t108 = cos(pkin(9));
t104 = qJD(1) + qJD(2);
t101 = qJD(3) + t104;
t99 = t101 ^ 2;
t121 = -t100 * t107 - t108 * t99;
t77 = t100 * t108 - t107 * t99;
t49 = -t110 * t77 + t114 * t121;
t50 = t110 * t121 + t114 * t77;
t109 = sin(qJ(5));
t113 = cos(qJ(5));
t136 = -g(1) + qJDD(4);
t111 = sin(qJ(2));
t115 = cos(qJ(2));
t112 = sin(qJ(1));
t116 = cos(qJ(1));
t123 = -g(2) * t116 - g(3) * t112;
t87 = qJDD(1) * pkin(1) + t123;
t122 = g(2) * t112 - g(3) * t116;
t88 = -qJD(1) ^ 2 * pkin(1) - t122;
t60 = -t111 * t88 + t115 * t87;
t54 = pkin(2) * t103 + t60;
t102 = t104 ^ 2;
t61 = t111 * t87 + t115 * t88;
t55 = -pkin(2) * t102 + t61;
t37 = -t110 * t55 + t114 * t54;
t35 = pkin(3) * t100 + t37;
t38 = t110 * t54 + t114 * t55;
t36 = -pkin(3) * t99 + t38;
t23 = t107 * t35 + t108 * t36;
t21 = -pkin(4) * t99 + pkin(8) * t100 + t23;
t15 = t109 * t21 - t113 * t136;
t16 = t109 * t136 + t113 * t21;
t9 = t109 * t15 + t113 * t16;
t22 = -t107 * t36 + t108 * t35;
t11 = t107 * t23 + t108 * t22;
t10 = pkin(3) * t11;
t12 = -t107 * t22 + t108 * t23;
t4 = t11 * t114 + t110 * t12;
t140 = pkin(2) * t4 + t10;
t91 = t113 * t99 * t109;
t85 = qJDD(5) + t91;
t139 = t109 * t85;
t86 = qJDD(5) - t91;
t138 = t113 * t86;
t137 = t109 * t100;
t135 = qJD(5) * t101;
t20 = -pkin(4) * t100 - pkin(8) * t99 - t22;
t6 = t107 * t9 - t108 * t20;
t134 = pkin(3) * t6 - pkin(4) * t20 + pkin(8) * t9;
t133 = pkin(3) * t121 - t23;
t7 = t107 * t20 + t108 * t9;
t2 = t110 * t7 + t114 * t6;
t132 = pkin(2) * t2 + t134;
t117 = qJD(5) ^ 2;
t105 = t109 ^ 2;
t94 = t105 * t99;
t89 = -t94 - t117;
t65 = -t109 * t89 - t138;
t73 = 0.2e1 * t113 * t135 + t137;
t42 = t107 * t65 - t108 * t73;
t131 = pkin(3) * t42 - pkin(4) * t73 + pkin(8) * t65 + t109 * t20;
t106 = t113 ^ 2;
t95 = t106 * t99;
t90 = -t95 - t117;
t64 = t113 * t90 - t139;
t93 = t113 * t100;
t74 = -0.2e1 * t109 * t135 + t93;
t41 = t107 * t64 + t108 * t74;
t130 = pkin(3) * t41 + pkin(4) * t74 + pkin(8) * t64 - t113 * t20;
t129 = pkin(2) * t49 + t133;
t128 = pkin(3) * t77 + t22;
t82 = t100 * t114 - t110 * t99;
t127 = pkin(2) * t82 + t37;
t80 = (t105 + t106) * t100;
t83 = t94 + t95;
t51 = t107 * t80 + t108 * t83;
t126 = pkin(3) * t51 + pkin(4) * t83 + pkin(8) * t80 + t9;
t43 = -t107 * t74 + t108 * t64;
t28 = t110 * t43 + t114 * t41;
t125 = pkin(2) * t28 + t130;
t44 = t107 * t73 + t108 * t65;
t29 = t110 * t44 + t114 * t42;
t124 = pkin(2) * t29 + t131;
t81 = -t100 * t110 - t114 * t99;
t120 = pkin(2) * t50 + t128;
t52 = -t107 * t83 + t108 * t80;
t31 = t110 * t52 + t114 * t51;
t119 = pkin(2) * t31 + t126;
t118 = pkin(2) * t81 - t38;
t63 = t139 + t113 * (-t94 + t117);
t62 = t109 * (t95 - t117) + t138;
t57 = t73 * t109;
t56 = t74 * t113;
t47 = t109 * t74 + t113 * t73;
t25 = t110 * t38 + t114 * t37;
t24 = pkin(2) * t25;
t1 = [0, 0, 0, 0, 0, qJDD(1), t123, t122, 0, 0, 0, 0, 0, 0, 0, t103, pkin(1) * (-t102 * t111 + t103 * t115) + t60, pkin(1) * (-t102 * t115 - t103 * t111) - t61, 0, pkin(1) * (t111 * t61 + t115 * t60), 0, 0, 0, 0, 0, t100, pkin(1) * (t111 * t81 + t115 * t82) + t127, pkin(1) * (-t111 * t82 + t115 * t81) + t118, 0, pkin(1) * (t111 * (-t110 * t37 + t114 * t38) + t115 * t25) + t24, 0, 0, 0, 0, 0, t100, pkin(1) * (t111 * t49 + t115 * t50) + t120, pkin(1) * (-t111 * t50 + t115 * t49) + t129, 0, pkin(1) * (t111 * (-t11 * t110 + t114 * t12) + t115 * t4) + t140, t57, t47, t63, t56, t62, 0, pkin(1) * (t111 * (-t110 * t41 + t114 * t43) + t115 * t28) + t125, pkin(1) * (t111 * (-t110 * t42 + t114 * t44) + t115 * t29) + t124, pkin(1) * (t111 * (-t110 * t51 + t114 * t52) + t115 * t31) + t119, pkin(1) * (t111 * (-t110 * t6 + t114 * t7) + t115 * t2) + t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, t60, -t61, 0, 0, 0, 0, 0, 0, 0, t100, t127, t118, 0, t24, 0, 0, 0, 0, 0, t100, t120, t129, 0, t140, t57, t47, t63, t56, t62, 0, t125, t124, t119, t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t37, -t38, 0, 0, 0, 0, 0, 0, 0, t100, t128, t133, 0, t10, t57, t47, t63, t56, t62, 0, t130, t131, t126, t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, 0, 0, 0, 0, 0, 0, t109 * t90 + t113 * t85, -t109 * t86 + t113 * t89, 0, t109 * t16 - t113 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t94 - t95, t137, t91, t93, qJDD(5), -t15, -t16, 0, 0;];
tauJ_reg = t1;
