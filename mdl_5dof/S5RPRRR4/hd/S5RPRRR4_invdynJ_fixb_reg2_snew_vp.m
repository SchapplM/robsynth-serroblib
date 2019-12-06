% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:36
% EndTime: 2019-12-05 18:14:40
% DurationCPUTime: 0.71s
% Computational Cost: add. (3491->125), mult. (4975->187), div. (0->0), fcn. (2890->10), ass. (0->97)
t102 = sin(qJ(3));
t106 = cos(qJ(3));
t101 = sin(qJ(4));
t105 = cos(qJ(4));
t94 = qJD(1) + qJD(3);
t89 = qJD(4) + t94;
t87 = t89 ^ 2;
t93 = qJDD(1) + qJDD(3);
t88 = qJDD(4) + t93;
t112 = -t101 * t88 - t105 * t87;
t65 = t101 * t87 - t105 * t88;
t40 = t102 * t65 + t106 * t112;
t132 = -t102 * t112 + t106 * t65;
t100 = sin(qJ(5));
t104 = cos(qJ(5));
t103 = sin(qJ(1));
t107 = cos(qJ(1));
t127 = t107 * g(2) + t103 * g(3);
t75 = qJDD(1) * pkin(1) + t127;
t109 = qJD(1) ^ 2;
t113 = -t103 * g(2) + t107 * g(3);
t76 = -t109 * pkin(1) - t113;
t98 = sin(pkin(9));
t99 = cos(pkin(9));
t119 = t99 * t75 - t98 * t76;
t44 = qJDD(1) * pkin(2) + t119;
t128 = t98 * t75 + t99 * t76;
t45 = -t109 * pkin(2) + t128;
t29 = -t102 * t45 + t106 * t44;
t27 = t93 * pkin(3) + t29;
t30 = t102 * t44 + t106 * t45;
t92 = t94 ^ 2;
t28 = -t92 * pkin(3) + t30;
t21 = t101 * t27 + t105 * t28;
t19 = -t87 * pkin(4) + t88 * pkin(8) + t21;
t97 = -g(1) + qJDD(2);
t13 = t100 * t19 - t104 * t97;
t14 = t100 * t97 + t104 * t19;
t7 = t100 * t13 + t104 * t14;
t20 = -t101 * t28 + t105 * t27;
t18 = -t88 * pkin(4) - t87 * pkin(8) - t20;
t129 = -pkin(4) * t18 + pkin(8) * t7;
t79 = t104 * t87 * t100;
t73 = qJDD(5) + t79;
t126 = t100 * t73;
t125 = t100 * t88;
t74 = qJDD(5) - t79;
t124 = t104 * t74;
t123 = qJD(5) * t89;
t4 = t101 * t7 - t105 * t18;
t122 = pkin(3) * t4 + t129;
t108 = qJD(5) ^ 2;
t95 = t100 ^ 2;
t82 = t95 * t87;
t77 = -t82 - t108;
t53 = -t100 * t77 - t124;
t59 = 0.2e1 * t104 * t123 + t125;
t121 = -pkin(4) * t59 + pkin(8) * t53 + t100 * t18;
t96 = t104 ^ 2;
t83 = t96 * t87;
t78 = -t83 - t108;
t52 = t104 * t78 - t126;
t81 = t104 * t88;
t60 = -0.2e1 * t100 * t123 + t81;
t120 = pkin(4) * t60 + pkin(8) * t52 - t104 * t18;
t62 = (t95 + t96) * t88;
t67 = t82 + t83;
t118 = pkin(4) * t67 + pkin(8) * t62 + t7;
t34 = t101 * t53 - t105 * t59;
t117 = pkin(3) * t34 + t121;
t33 = t101 * t52 + t105 * t60;
t116 = pkin(3) * t33 + t120;
t115 = -pkin(3) * t65 + t20;
t39 = t101 * t62 + t105 * t67;
t114 = pkin(3) * t39 + t118;
t71 = -t102 * t93 - t106 * t92;
t111 = t102 * t92 - t106 * t93;
t110 = pkin(3) * t112 - t21;
t51 = t126 + t104 * (-t82 + t108);
t50 = t100 * (t83 - t108) + t124;
t47 = t59 * t100;
t46 = t60 * t104;
t42 = -t101 * t67 + t105 * t62;
t37 = t100 * t60 + t104 * t59;
t36 = t101 * t59 + t105 * t53;
t35 = -t101 * t60 + t105 * t52;
t25 = t102 * t42 + t106 * t39;
t24 = t102 * t36 + t106 * t34;
t23 = t102 * t35 + t106 * t33;
t22 = t102 * t30 + t106 * t29;
t10 = -t101 * t20 + t105 * t21;
t9 = t101 * t21 + t105 * t20;
t8 = pkin(3) * t9;
t5 = t101 * t18 + t105 * t7;
t2 = t102 * t10 + t106 * t9;
t1 = t102 * t5 + t106 * t4;
t3 = [0, 0, 0, 0, 0, qJDD(1), t127, t113, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t99 * qJDD(1) - t98 * t109) + t119, pkin(1) * (-t98 * qJDD(1) - t99 * t109) - t128, 0, pkin(1) * (t99 * t119 + t98 * t128), 0, 0, 0, 0, 0, t93, pkin(1) * (-t111 * t99 + t98 * t71) - pkin(2) * t111 + t29, pkin(1) * (t98 * t111 + t99 * t71) + pkin(2) * t71 - t30, 0, pkin(1) * (t98 * (-t102 * t29 + t106 * t30) + t99 * t22) + pkin(2) * t22, 0, 0, 0, 0, 0, t88, pkin(1) * (-t132 * t99 + t98 * t40) - pkin(2) * t132 + t115, pkin(1) * (t98 * t132 + t99 * t40) + pkin(2) * t40 + t110, 0, pkin(1) * (t98 * (t106 * t10 - t102 * t9) + t99 * t2) + pkin(2) * t2 + t8, t47, t37, t51, t46, t50, 0, pkin(1) * (t98 * (-t102 * t33 + t106 * t35) + t99 * t23) + pkin(2) * t23 + t116, pkin(1) * (t98 * (-t102 * t34 + t106 * t36) + t99 * t24) + pkin(2) * t24 + t117, pkin(1) * (t98 * (-t102 * t39 + t106 * t42) + t99 * t25) + pkin(2) * t25 + t114, pkin(1) * (t98 * (-t102 * t4 + t106 * t5) + t99 * t1) + pkin(2) * t1 + t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, 0, 0, 0, t100 * t78 + t104 * t73, -t100 * t74 + t104 * t77, 0, t100 * t14 - t104 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t29, -t30, 0, 0, 0, 0, 0, 0, 0, t88, t115, t110, 0, t8, t47, t37, t51, t46, t50, 0, t116, t117, t114, t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t20, -t21, 0, 0, t47, t37, t51, t46, t50, 0, t120, t121, t118, t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t82 - t83, t125, t79, t81, qJDD(5), -t13, -t14, 0, 0;];
tauJ_reg = t3;
