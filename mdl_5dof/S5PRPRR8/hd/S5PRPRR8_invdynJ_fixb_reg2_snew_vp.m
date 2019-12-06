% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPRR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:27
% EndTime: 2019-12-05 16:04:32
% DurationCPUTime: 0.90s
% Computational Cost: add. (2166->176), mult. (4086->247), div. (0->0), fcn. (2790->10), ass. (0->121)
t100 = cos(qJ(4));
t121 = qJD(2) * t100;
t96 = sin(qJ(5));
t99 = cos(qJ(5));
t60 = -t99 * qJD(4) + t96 * t121;
t62 = t96 * qJD(4) + t99 * t121;
t44 = t62 * t60;
t119 = qJD(2) * qJD(4);
t114 = t100 * t119;
t97 = sin(qJ(4));
t120 = t97 * qJDD(2);
t65 = -t114 - t120;
t58 = qJDD(5) - t65;
t139 = -t44 + t58;
t142 = t139 * t96;
t141 = t139 * t99;
t91 = sin(pkin(9));
t93 = cos(pkin(9));
t70 = t91 * g(1) - t93 * g(2);
t88 = -g(3) + qJDD(1);
t92 = sin(pkin(5));
t94 = cos(pkin(5));
t140 = t70 * t94 + t88 * t92;
t116 = t97 * t119;
t81 = t100 * qJDD(2);
t66 = t81 - t116;
t115 = -t99 * qJDD(4) + t96 * t66;
t79 = t97 * qJD(2) + qJD(5);
t23 = (qJD(5) - t79) * t62 + t115;
t56 = t60 ^ 2;
t57 = t62 ^ 2;
t78 = t79 ^ 2;
t138 = -pkin(7) - pkin(2);
t136 = t79 * t96;
t135 = t79 * t99;
t102 = qJD(4) ^ 2;
t103 = qJD(2) ^ 2;
t101 = cos(qJ(2));
t71 = -t93 * g(1) - t91 * g(2);
t98 = sin(qJ(2));
t37 = t140 * t101 - t98 * t71;
t106 = qJDD(3) - t37;
t90 = qJDD(2) * pkin(2);
t32 = -t103 * qJ(3) + t106 - t90;
t104 = -qJDD(2) * pkin(7) + t32;
t50 = -t92 * t70 + t94 * t88;
t19 = -t100 * t104 + t97 * t50;
t110 = t97 * pkin(4) - t100 * pkin(8);
t63 = t110 * qJD(2);
t13 = -qJDD(4) * pkin(4) - t102 * pkin(8) + t63 * t121 + t19;
t133 = t96 * t13;
t35 = t44 + t58;
t132 = t96 * t35;
t117 = t100 * t103 * t97;
t72 = qJDD(4) + t117;
t131 = t97 * t72;
t130 = t99 * t13;
t129 = t99 * t35;
t86 = t97 ^ 2;
t87 = t100 ^ 2;
t128 = t86 + t87;
t127 = t100 * t50;
t73 = qJDD(4) - t117;
t126 = t100 * t73;
t125 = t86 * t103;
t124 = t87 * t103;
t122 = qJD(5) + t79;
t118 = t97 * t44;
t38 = t101 * t71 + t140 * t98;
t14 = -t102 * pkin(4) + qJDD(4) * pkin(8) + t127 + (-qJD(2) * t63 + t104) * t97;
t108 = -t65 + t114;
t109 = -t66 + t116;
t113 = 0.2e1 * qJD(3) * qJD(2) + t38;
t84 = qJDD(2) * qJ(3);
t112 = t84 + t113;
t30 = t138 * t103 + t112;
t16 = t108 * pkin(4) + t109 * pkin(8) + t30;
t5 = t96 * t14 - t99 * t16;
t6 = t99 * t14 + t96 * t16;
t3 = t96 * t5 + t99 * t6;
t111 = t99 * t5 - t96 * t6;
t20 = t104 * t97 + t127;
t8 = -t100 * t19 + t97 * t20;
t107 = -t96 * qJDD(4) - t99 * t66;
t105 = qJ(3) + t110;
t40 = -t60 * qJD(5) - t107;
t77 = -t102 - t124;
t76 = -t102 - t125;
t69 = t128 * t103;
t68 = t128 * qJDD(2);
t67 = t81 - 0.2e1 * t116;
t64 = 0.2e1 * t114 + t120;
t53 = (-qJDD(2) * t101 + t103 * t98) * t92;
t52 = (qJDD(2) * t98 + t101 * t103) * t92;
t51 = t79 * t60;
t49 = -t57 + t78;
t48 = t56 - t78;
t47 = t100 * t77 - t131;
t46 = t97 * t76 + t126;
t45 = t94 * t50;
t43 = t57 - t56;
t42 = -t57 - t78;
t41 = -t78 - t56;
t39 = -t62 * qJD(5) - t115;
t33 = t56 + t57;
t31 = -t103 * pkin(2) + t112;
t28 = t122 * t60 + t107;
t27 = t40 + t51;
t26 = t40 - t51;
t24 = -t122 * t62 - t115;
t22 = -t96 * t42 - t129;
t21 = t99 * t42 - t132;
t18 = t99 * t41 - t142;
t17 = t96 * t41 + t141;
t12 = -t23 * t99 + t96 * t27;
t11 = -t23 * t96 - t99 * t27;
t10 = t100 * t28 + t97 * t22;
t9 = t100 * t24 + t97 * t18;
t7 = t100 * t33 + t97 * t12;
t1 = -t100 * t13 + t97 * t3;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, 0, 0, 0, 0, 0, -t53, -t52, 0, t45 + (t101 * t37 + t38 * t98) * t92, 0, 0, 0, 0, 0, 0, 0, t53, t52, t45 + (-t101 * t32 + t31 * t98) * t92, 0, 0, 0, 0, 0, 0, t94 * (t100 * t76 - t97 * t73) + (-t101 * t46 + t98 * t64) * t92, t94 * (-t100 * t72 - t97 * t77) + (-t101 * t47 + t98 * t67) * t92, (t101 * t68 - t69 * t98) * t92, t94 * (t100 * t20 + t97 * t19) + (-t101 * t8 + t98 * t30) * t92, 0, 0, 0, 0, 0, 0, t94 * (t100 * t18 - t97 * t24) + (-t101 * t9 + t98 * t17) * t92, t94 * (t100 * t22 - t97 * t28) + (-t101 * t10 + t98 * t21) * t92, t94 * (t100 * t12 - t97 * t33) + (-t101 * t7 + t98 * t11) * t92, t94 * (t100 * t3 + t97 * t13) + (-t101 * t1 - t111 * t98) * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t37, -t38, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, t106 - 0.2e1 * t90, 0.2e1 * t84 + t113, -pkin(2) * t32 + qJ(3) * t31, -t109 * t100, -t100 * t64 - t97 * t67, t126 - t97 * (t102 - t124), t108 * t97, t100 * (-t102 + t125) - t131, 0, qJ(3) * t64 + t138 * t46 + t97 * t30, qJ(3) * t67 + t100 * t30 + t138 * t47, -qJ(3) * t69 - t138 * t68 - t8, qJ(3) * t30 + t138 * t8, t100 * (-t62 * t136 + t99 * t40) + t118, t100 * (t99 * t24 - t96 * t26) + t97 * t43, t100 * (-t96 * t49 + t141) + t97 * t27, t100 * (t60 * t135 - t96 * t39) - t118, t100 * (t99 * t48 - t132) - t97 * t23, t97 * t58 + t100 * (-t60 * t99 + t62 * t96) * t79, t100 * (-pkin(8) * t17 + t133) - t97 * (-pkin(4) * t17 + t5) + qJ(3) * t17 + t138 * t9, t100 * (-pkin(8) * t21 + t130) - t97 * (-pkin(4) * t21 + t6) + qJ(3) * t21 + t138 * t10, t100 * t111 + t105 * t11 + t138 * t7, t138 * t1 - t105 * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t103, t32, 0, 0, 0, 0, 0, 0, t46, t47, -t68, t8, 0, 0, 0, 0, 0, 0, t9, t10, t7, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, (-t86 + t87) * t103, t81, -t117, -t120, qJDD(4), -t19, -t20, 0, 0, t62 * t135 + t96 * t40, t96 * t24 + t99 * t26, t99 * t49 + t142, t60 * t136 + t99 * t39, t96 * t48 + t129, (-t60 * t96 - t62 * t99) * t79, pkin(4) * t24 + pkin(8) * t18 - t130, pkin(4) * t28 + pkin(8) * t22 + t133, pkin(4) * t33 + pkin(8) * t12 + t3, -pkin(4) * t13 + pkin(8) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t43, t27, -t44, -t23, t58, -t5, -t6, 0, 0;];
tauJ_reg = t2;
