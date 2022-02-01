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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:30:53
% EndTime: 2022-01-20 11:30:56
% DurationCPUTime: 0.76s
% Computational Cost: add. (4831->125), mult. (6090->187), div. (0->0), fcn. (3542->10), ass. (0->102)
t111 = sin(qJ(3));
t115 = cos(qJ(3));
t104 = qJDD(1) + qJDD(2);
t100 = qJDD(3) + t104;
t108 = sin(pkin(9));
t109 = cos(pkin(9));
t105 = qJD(1) + qJD(2);
t101 = qJD(3) + t105;
t99 = t101 ^ 2;
t122 = -t108 * t100 - t109 * t99;
t77 = t109 * t100 - t108 * t99;
t49 = -t111 * t77 + t115 * t122;
t50 = t111 * t122 + t115 * t77;
t110 = sin(qJ(5));
t114 = cos(qJ(5));
t137 = -g(3) + qJDD(4);
t112 = sin(qJ(2));
t116 = cos(qJ(2));
t113 = sin(qJ(1));
t117 = cos(qJ(1));
t132 = t113 * g(1) - t117 * g(2);
t87 = qJDD(1) * pkin(1) + t132;
t123 = t117 * g(1) + t113 * g(2);
t88 = -qJD(1) ^ 2 * pkin(1) - t123;
t60 = -t112 * t88 + t116 * t87;
t54 = t104 * pkin(2) + t60;
t103 = t105 ^ 2;
t61 = t112 * t87 + t116 * t88;
t55 = -t103 * pkin(2) + t61;
t37 = -t111 * t55 + t115 * t54;
t35 = t100 * pkin(3) + t37;
t38 = t111 * t54 + t115 * t55;
t36 = -t99 * pkin(3) + t38;
t23 = t108 * t35 + t109 * t36;
t21 = -t99 * pkin(4) + t100 * pkin(8) + t23;
t15 = t110 * t21 - t114 * t137;
t16 = t110 * t137 + t114 * t21;
t9 = t110 * t15 + t114 * t16;
t22 = -t108 * t36 + t109 * t35;
t11 = t108 * t23 + t109 * t22;
t10 = pkin(3) * t11;
t12 = -t108 * t22 + t109 * t23;
t4 = t115 * t11 + t111 * t12;
t141 = pkin(2) * t4 + t10;
t91 = t114 * t99 * t110;
t85 = qJDD(5) + t91;
t140 = t110 * t85;
t86 = qJDD(5) - t91;
t139 = t114 * t86;
t138 = t110 * t100;
t136 = qJD(5) * t101;
t20 = -t100 * pkin(4) - t99 * pkin(8) - t22;
t6 = t108 * t9 - t109 * t20;
t135 = pkin(3) * t6 - pkin(4) * t20 + pkin(8) * t9;
t134 = pkin(3) * t122 - t23;
t7 = t108 * t20 + t109 * t9;
t2 = t111 * t7 + t115 * t6;
t133 = pkin(2) * t2 + t135;
t118 = qJD(5) ^ 2;
t106 = t110 ^ 2;
t94 = t106 * t99;
t89 = -t94 - t118;
t65 = -t110 * t89 - t139;
t73 = 0.2e1 * t114 * t136 + t138;
t42 = t108 * t65 - t109 * t73;
t131 = pkin(3) * t42 - pkin(4) * t73 + pkin(8) * t65 + t110 * t20;
t107 = t114 ^ 2;
t95 = t107 * t99;
t90 = -t95 - t118;
t64 = t114 * t90 - t140;
t93 = t114 * t100;
t74 = -0.2e1 * t110 * t136 + t93;
t41 = t108 * t64 + t109 * t74;
t130 = pkin(3) * t41 + pkin(4) * t74 + pkin(8) * t64 - t114 * t20;
t129 = pkin(2) * t49 + t134;
t128 = pkin(3) * t77 + t22;
t82 = t115 * t100 - t111 * t99;
t127 = pkin(2) * t82 + t37;
t80 = (t106 + t107) * t100;
t83 = t94 + t95;
t51 = t108 * t80 + t109 * t83;
t126 = pkin(3) * t51 + pkin(4) * t83 + pkin(8) * t80 + t9;
t43 = -t108 * t74 + t109 * t64;
t28 = t111 * t43 + t115 * t41;
t125 = pkin(2) * t28 + t130;
t44 = t108 * t73 + t109 * t65;
t29 = t111 * t44 + t115 * t42;
t124 = pkin(2) * t29 + t131;
t81 = -t111 * t100 - t115 * t99;
t121 = pkin(2) * t50 + t128;
t52 = -t108 * t83 + t109 * t80;
t31 = t111 * t52 + t115 * t51;
t120 = pkin(2) * t31 + t126;
t119 = pkin(2) * t81 - t38;
t63 = t140 + t114 * (-t94 + t118);
t62 = t110 * (t95 - t118) + t139;
t57 = t73 * t110;
t56 = t74 * t114;
t47 = t110 * t74 + t114 * t73;
t25 = t111 * t38 + t115 * t37;
t24 = pkin(2) * t25;
t1 = [0, 0, 0, 0, 0, qJDD(1), t132, t123, 0, 0, 0, 0, 0, 0, 0, t104, pkin(1) * (-t112 * t103 + t116 * t104) + t60, pkin(1) * (-t116 * t103 - t112 * t104) - t61, 0, pkin(1) * (t112 * t61 + t116 * t60), 0, 0, 0, 0, 0, t100, pkin(1) * (t112 * t81 + t116 * t82) + t127, pkin(1) * (-t112 * t82 + t116 * t81) + t119, 0, pkin(1) * (t112 * (-t111 * t37 + t115 * t38) + t116 * t25) + t24, 0, 0, 0, 0, 0, t100, pkin(1) * (t112 * t49 + t116 * t50) + t121, pkin(1) * (-t112 * t50 + t116 * t49) + t129, 0, pkin(1) * (t112 * (-t111 * t11 + t115 * t12) + t116 * t4) + t141, t57, t47, t63, t56, t62, 0, pkin(1) * (t112 * (-t111 * t41 + t115 * t43) + t116 * t28) + t125, pkin(1) * (t112 * (-t111 * t42 + t115 * t44) + t116 * t29) + t124, pkin(1) * (t112 * (-t111 * t51 + t115 * t52) + t116 * t31) + t120, pkin(1) * (t112 * (-t111 * t6 + t115 * t7) + t116 * t2) + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t60, -t61, 0, 0, 0, 0, 0, 0, 0, t100, t127, t119, 0, t24, 0, 0, 0, 0, 0, t100, t121, t129, 0, t141, t57, t47, t63, t56, t62, 0, t125, t124, t120, t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t37, -t38, 0, 0, 0, 0, 0, 0, 0, t100, t128, t134, 0, t10, t57, t47, t63, t56, t62, 0, t130, t131, t126, t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, 0, 0, 0, 0, 0, 0, t110 * t90 + t114 * t85, -t110 * t86 + t114 * t89, 0, t110 * t16 - t114 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t94 - t95, t138, t91, t93, qJDD(5), -t15, -t16, 0, 0;];
tauJ_reg = t1;
