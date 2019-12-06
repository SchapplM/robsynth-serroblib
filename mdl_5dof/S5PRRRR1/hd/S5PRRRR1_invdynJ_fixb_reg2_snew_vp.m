% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRR1
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:03:23
% EndTime: 2019-12-05 17:03:29
% DurationCPUTime: 1.44s
% Computational Cost: add. (2362->211), mult. (4760->312), div. (0->0), fcn. (3723->8), ass. (0->151)
t125 = sin(qJ(3));
t124 = sin(qJ(4));
t128 = cos(qJ(4));
t122 = -g(3) + qJDD(1);
t126 = sin(qJ(2));
t130 = cos(qJ(2));
t105 = t130 * g(1) - t126 * t122;
t129 = cos(qJ(3));
t138 = t129 * g(2) + t125 * t105;
t132 = qJD(2) ^ 2;
t108 = t125 * t132 * t129;
t147 = qJDD(3) + t108;
t133 = pkin(2) * t147 + t138;
t121 = t129 ^ 2;
t116 = t121 * t132;
t131 = qJD(3) ^ 2;
t169 = -t116 - t131;
t89 = t125 * g(2) - t129 * t105;
t78 = t169 * pkin(2) + t89;
t57 = t124 * t78 - t128 * t133;
t48 = t128 * t57;
t159 = t128 * t78;
t58 = t124 * t133 + t159;
t172 = -t124 * t58 + t48;
t178 = t125 * t172;
t123 = sin(qJ(5));
t119 = qJD(3) + qJD(4);
t127 = cos(qJ(5));
t154 = t125 * t128;
t97 = (t129 * t124 + t154) * qJD(2);
t81 = -t127 * t119 + t123 * t97;
t83 = t123 * t119 + t127 * t97;
t61 = t83 * t81;
t113 = t125 * qJDD(2);
t149 = qJD(2) * qJD(3);
t143 = t129 * t149;
t101 = t113 + t143;
t114 = t129 * qJDD(2);
t144 = t125 * t149;
t102 = t114 - t144;
t140 = t124 * t101 - t128 * t102;
t63 = -t97 * qJD(4) - t140;
t62 = qJDD(5) - t63;
t173 = -t61 + t62;
t177 = t123 * t173;
t118 = qJDD(3) + qJDD(4);
t155 = t125 * t124;
t95 = (-t129 * t128 + t155) * qJD(2);
t77 = t97 * t95;
t171 = -t77 + t118;
t176 = t124 * t171;
t175 = t127 * t173;
t174 = t128 * t171;
t170 = t102 - t144;
t135 = t128 * t101 + t124 * t102;
t64 = -t95 * qJD(4) + t135;
t139 = -t127 * t118 + t123 * t64;
t92 = qJD(5) + t95;
t26 = (qJD(5) - t92) * t83 + t139;
t79 = t81 ^ 2;
t80 = t83 ^ 2;
t91 = t92 ^ 2;
t93 = t95 ^ 2;
t94 = t97 ^ 2;
t117 = t119 ^ 2;
t168 = pkin(2) * t129;
t104 = t126 * g(1) + t130 * t122;
t76 = t170 * pkin(2) + t104;
t34 = t123 * t58 + t127 * t76;
t35 = -t123 * t76 + t127 * t58;
t11 = t123 * t34 + t127 * t35;
t37 = t61 + t62;
t167 = t123 * t37;
t166 = t123 * t92;
t165 = t124 * t57;
t71 = t77 + t118;
t163 = t124 * t71;
t162 = t127 * t37;
t47 = t127 * t57;
t161 = t127 * t92;
t160 = t128 * t71;
t158 = t119 * t124;
t157 = t119 * t128;
t156 = t125 * t147;
t153 = t129 * (qJDD(3) - t108);
t151 = qJD(5) + t92;
t150 = qJD(4) + t119;
t148 = t126 * qJDD(2);
t146 = t124 * t61;
t145 = t128 * t61;
t142 = t128 * t58 + t165;
t141 = -t125 * t138 + t129 * t89;
t9 = t123 * t35 - t127 * t34;
t136 = -t123 * t118 - t127 * t64;
t134 = (-qJD(4) + t119) * t97 - t140;
t42 = -t81 * qJD(5) - t136;
t120 = t125 ^ 2;
t115 = t120 * t132;
t103 = t114 - 0.2e1 * t144;
t100 = t113 + 0.2e1 * t143;
t99 = t130 * t104;
t90 = t119 * t95;
t87 = -t94 + t117;
t86 = t93 - t117;
t84 = -t94 - t117;
t75 = t94 - t93;
t69 = -t117 - t93;
t68 = t92 * t81;
t67 = -t80 + t91;
t66 = t79 - t91;
t65 = -t93 - t94;
t60 = t80 - t79;
t59 = -t80 - t91;
t56 = t128 * t84 - t163;
t55 = t64 + t90;
t54 = t64 - t90;
t53 = -t150 * t95 + t135;
t50 = t150 * t97 + t140;
t46 = t123 * t57;
t45 = -t91 - t79;
t44 = t124 * t69 + t174;
t43 = t79 + t80;
t41 = -t83 * qJD(5) - t139;
t40 = (t123 * t83 - t127 * t81) * t92;
t39 = (-t123 * t81 - t127 * t83) * t92;
t31 = t151 * t81 + t136;
t30 = t42 + t68;
t29 = t42 - t68;
t28 = -t151 * t83 - t139;
t25 = t127 * t42 - t83 * t166;
t24 = t123 * t42 + t83 * t161;
t23 = -t123 * t41 + t81 * t161;
t22 = t127 * t41 + t81 * t166;
t20 = t124 * t134 - t128 * t55;
t19 = t127 * t66 - t167;
t18 = -t123 * t67 + t175;
t17 = t123 * t66 + t162;
t16 = t127 * t67 + t177;
t15 = -t123 * t59 - t162;
t14 = t127 * t59 - t167;
t13 = t127 * t45 - t177;
t12 = t123 * t45 + t175;
t8 = -t123 * t29 + t127 * t28;
t7 = t123 * t30 - t127 * t26;
t6 = t123 * t28 + t127 * t29;
t5 = -t123 * t26 - t127 * t30;
t4 = t124 * t11 - t48;
t3 = t124 * t15 + t128 * t31;
t2 = t124 * t13 + t128 * t28;
t1 = t124 * t7 + t128 * t43;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t122, 0, 0, 0, 0, 0, 0, t130 * qJDD(2) - t126 * t132, -t130 * t132 - t148, 0, -t126 * t105 + t99, 0, 0, 0, 0, 0, 0, t126 * (t129 * t169 - t156) + t130 * t103, t126 * (-t153 - t125 * (-t115 - t131)) - t130 * t100, t130 * (t115 + t116) + (t120 + t121) * t148, t126 * t141 + t99, 0, 0, 0, 0, 0, 0, t126 * (t129 * (t128 * t69 - t176) - t125 * t44) - t130 * t50, t126 * (t129 * (-t124 * t84 - t160) - t125 * t56) - t130 * t53, t126 * (t129 * (t124 * t55 + t128 * t134) - t125 * t20) - t130 * t65, t126 * (t129 * t142 + t178) + t130 * t76, 0, 0, 0, 0, 0, 0, t126 * (t129 * (-t124 * t28 + t128 * t13) - t125 * t2) - t130 * t12, t126 * (t129 * (-t124 * t31 + t128 * t15) - t125 * t3) - t130 * t14, t126 * (t129 * (-t124 * t43 + t128 * t7) - t125 * t1) - t130 * t5, t126 * (t129 * (t128 * t11 + t165) - t125 * t4) - t130 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t104, t105, 0, 0, (t101 + t143) * t125, t129 * t100 + t125 * t103, t156 + t129 * (-t115 + t131), t170 * t129, t125 * (t116 - t131) + t153, 0, t129 * t104, -t125 * t104, t141, 0, t125 * (t128 * t64 - t97 * t158) + t129 * (t124 * t64 + t97 * t157), t125 * (-t124 * t54 - t128 * t50) + t129 * (-t124 * t50 + t128 * t54), t125 * (-t124 * t87 + t174) + t129 * (t128 * t87 + t176), t125 * (-t124 * t63 + t95 * t157) + t129 * (t128 * t63 + t95 * t158), t125 * (t128 * t86 - t163) + t129 * (t124 * t86 + t160), (t125 * (t124 * t97 - t128 * t95) + t129 * (-t124 * t95 - t128 * t97)) * t119, -t76 * t155 + t129 * (-pkin(2) * t50 + t128 * t76), -t76 * t154 + t129 * (-pkin(2) * t53 - t124 * t76), t178 + t129 * (-pkin(2) * t65 + t142), t76 * t168, t125 * (t128 * t25 + t146) + t129 * (t124 * t25 - t145), t125 * (t124 * t60 + t128 * t8) + t129 * (t124 * t8 - t128 * t60), t125 * (t124 * t30 + t128 * t18) + t129 * (t124 * t18 - t128 * t30), t125 * (t128 * t23 - t146) + t129 * (t124 * t23 + t145), t125 * (-t124 * t26 + t128 * t19) + t129 * (t124 * t19 + t128 * t26), t125 * (t124 * t62 + t128 * t40) + t129 * (t124 * t40 - t128 * t62), t125 * (t123 * t48 - t124 * t34) + t129 * (-pkin(2) * t12 + t123 * t165 + t128 * t34), t125 * (-t124 * t35 + t127 * t48) + t129 * (-pkin(2) * t14 + t124 * t47 + t128 * t35), -t9 * t154 + t129 * (-pkin(2) * t5 - t124 * t9), -t9 * t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, t115 - t116, t113, t108, t114, qJDD(3), t138, -t89, 0, 0, t77, t75, t55, -t77, t134, t118, pkin(2) * t44 - t57, -t159 - t124 * t138 + (-t124 * t147 + t56) * pkin(2), pkin(2) * t20, -pkin(2) * t172, t24, t6, t16, t22, t17, t39, pkin(2) * t2 - t47, pkin(2) * t3 + t46, pkin(2) * t1 + t11, pkin(2) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t75, t55, -t77, t134, t118, -t57, -t58, 0, 0, t24, t6, t16, t22, t17, t39, -t47, t46, t11, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t60, t30, -t61, -t26, t62, -t34, -t35, 0, 0;];
tauJ_reg = t10;
