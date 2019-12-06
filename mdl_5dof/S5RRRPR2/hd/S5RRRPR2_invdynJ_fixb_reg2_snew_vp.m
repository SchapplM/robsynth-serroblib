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
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:40:59
% EndTime: 2019-12-05 18:41:02
% DurationCPUTime: 0.72s
% Computational Cost: add. (4831->125), mult. (6090->187), div. (0->0), fcn. (3542->10), ass. (0->102)
t112 = sin(qJ(3));
t116 = cos(qJ(3));
t105 = qJDD(1) + qJDD(2);
t100 = qJDD(3) + t105;
t109 = sin(pkin(9));
t110 = cos(pkin(9));
t106 = qJD(1) + qJD(2);
t101 = qJD(3) + t106;
t99 = t101 ^ 2;
t123 = -t109 * t100 - t110 * t99;
t77 = t110 * t100 - t109 * t99;
t49 = -t112 * t77 + t116 * t123;
t50 = t112 * t123 + t116 * t77;
t111 = sin(qJ(5));
t115 = cos(qJ(5));
t138 = -g(1) + qJDD(4);
t113 = sin(qJ(2));
t117 = cos(qJ(2));
t114 = sin(qJ(1));
t118 = cos(qJ(1));
t137 = t118 * g(2) + t114 * g(3);
t87 = qJDD(1) * pkin(1) + t137;
t124 = -t114 * g(2) + t118 * g(3);
t88 = -qJD(1) ^ 2 * pkin(1) - t124;
t60 = -t113 * t88 + t117 * t87;
t54 = t105 * pkin(2) + t60;
t104 = t106 ^ 2;
t61 = t113 * t87 + t117 * t88;
t55 = -t104 * pkin(2) + t61;
t37 = -t112 * t55 + t116 * t54;
t35 = t100 * pkin(3) + t37;
t38 = t112 * t54 + t116 * t55;
t36 = -t99 * pkin(3) + t38;
t23 = t109 * t35 + t110 * t36;
t21 = -t99 * pkin(4) + t100 * pkin(8) + t23;
t15 = t111 * t21 - t115 * t138;
t16 = t111 * t138 + t115 * t21;
t9 = t111 * t15 + t115 * t16;
t22 = -t109 * t36 + t110 * t35;
t11 = t109 * t23 + t110 * t22;
t10 = pkin(3) * t11;
t12 = -t109 * t22 + t110 * t23;
t4 = t116 * t11 + t112 * t12;
t142 = pkin(2) * t4 + t10;
t91 = t115 * t99 * t111;
t85 = qJDD(5) + t91;
t141 = t111 * t85;
t86 = qJDD(5) - t91;
t140 = t115 * t86;
t139 = t111 * t100;
t136 = qJD(5) * t101;
t20 = -t100 * pkin(4) - t99 * pkin(8) - t22;
t6 = t109 * t9 - t110 * t20;
t135 = pkin(3) * t6 - pkin(4) * t20 + pkin(8) * t9;
t134 = pkin(3) * t123 - t23;
t7 = t109 * t20 + t110 * t9;
t2 = t112 * t7 + t116 * t6;
t133 = pkin(2) * t2 + t135;
t119 = qJD(5) ^ 2;
t107 = t111 ^ 2;
t94 = t107 * t99;
t89 = -t94 - t119;
t65 = -t111 * t89 - t140;
t73 = 0.2e1 * t115 * t136 + t139;
t42 = t109 * t65 - t110 * t73;
t132 = pkin(3) * t42 - pkin(4) * t73 + pkin(8) * t65 + t111 * t20;
t108 = t115 ^ 2;
t95 = t108 * t99;
t90 = -t95 - t119;
t64 = t115 * t90 - t141;
t93 = t115 * t100;
t74 = -0.2e1 * t111 * t136 + t93;
t41 = t109 * t64 + t110 * t74;
t131 = pkin(3) * t41 + pkin(4) * t74 + pkin(8) * t64 - t115 * t20;
t130 = pkin(2) * t49 + t134;
t129 = pkin(3) * t77 + t22;
t82 = t116 * t100 - t112 * t99;
t128 = pkin(2) * t82 + t37;
t80 = (t107 + t108) * t100;
t83 = t94 + t95;
t51 = t109 * t80 + t110 * t83;
t127 = pkin(3) * t51 + pkin(4) * t83 + pkin(8) * t80 + t9;
t43 = -t109 * t74 + t110 * t64;
t28 = t112 * t43 + t116 * t41;
t126 = pkin(2) * t28 + t131;
t44 = t109 * t73 + t110 * t65;
t29 = t112 * t44 + t116 * t42;
t125 = pkin(2) * t29 + t132;
t81 = -t112 * t100 - t116 * t99;
t122 = pkin(2) * t50 + t129;
t52 = -t109 * t83 + t110 * t80;
t31 = t112 * t52 + t116 * t51;
t121 = pkin(2) * t31 + t127;
t120 = pkin(2) * t81 - t38;
t63 = t141 + t115 * (-t94 + t119);
t62 = t111 * (t95 - t119) + t140;
t57 = t73 * t111;
t56 = t74 * t115;
t47 = t111 * t74 + t115 * t73;
t25 = t112 * t38 + t116 * t37;
t24 = pkin(2) * t25;
t1 = [0, 0, 0, 0, 0, qJDD(1), t137, t124, 0, 0, 0, 0, 0, 0, 0, t105, pkin(1) * (-t113 * t104 + t117 * t105) + t60, pkin(1) * (-t117 * t104 - t113 * t105) - t61, 0, pkin(1) * (t113 * t61 + t117 * t60), 0, 0, 0, 0, 0, t100, pkin(1) * (t113 * t81 + t117 * t82) + t128, pkin(1) * (-t113 * t82 + t117 * t81) + t120, 0, pkin(1) * (t113 * (-t112 * t37 + t116 * t38) + t117 * t25) + t24, 0, 0, 0, 0, 0, t100, pkin(1) * (t113 * t49 + t117 * t50) + t122, pkin(1) * (-t113 * t50 + t117 * t49) + t130, 0, pkin(1) * (t113 * (-t112 * t11 + t116 * t12) + t117 * t4) + t142, t57, t47, t63, t56, t62, 0, pkin(1) * (t113 * (-t112 * t41 + t116 * t43) + t117 * t28) + t126, pkin(1) * (t113 * (-t112 * t42 + t116 * t44) + t117 * t29) + t125, pkin(1) * (t113 * (-t112 * t51 + t116 * t52) + t117 * t31) + t121, pkin(1) * (t113 * (-t112 * t6 + t116 * t7) + t117 * t2) + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t60, -t61, 0, 0, 0, 0, 0, 0, 0, t100, t128, t120, 0, t24, 0, 0, 0, 0, 0, t100, t122, t130, 0, t142, t57, t47, t63, t56, t62, 0, t126, t125, t121, t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t37, -t38, 0, 0, 0, 0, 0, 0, 0, t100, t129, t134, 0, t10, t57, t47, t63, t56, t62, 0, t131, t132, t127, t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, 0, 0, 0, 0, 0, 0, t111 * t90 + t115 * t85, -t111 * t86 + t115 * t89, 0, t111 * t16 - t115 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t94 - t95, t139, t91, t93, qJDD(5), -t15, -t16, 0, 0;];
tauJ_reg = t1;
