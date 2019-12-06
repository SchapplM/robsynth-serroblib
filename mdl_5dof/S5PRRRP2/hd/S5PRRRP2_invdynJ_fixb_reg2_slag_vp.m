% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRP2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:42:06
% EndTime: 2019-12-05 16:42:09
% DurationCPUTime: 1.17s
% Computational Cost: add. (1216->211), mult. (1785->231), div. (0->0), fcn. (911->8), ass. (0->129)
t81 = qJD(2) + qJD(3);
t85 = sin(qJ(4));
t82 = t85 ^ 2;
t87 = cos(qJ(4));
t83 = t87 ^ 2;
t179 = (t82 + t83) * t81;
t79 = qJDD(2) + qJDD(3);
t153 = t83 * t79;
t154 = t82 * t79;
t177 = t153 + t154;
t86 = sin(qJ(3));
t132 = qJDD(2) * t86;
t88 = cos(qJ(3));
t140 = qJD(3) * t88;
t176 = qJD(1) * qJD(4) + t79 * pkin(7) + (qJD(2) * t140 + t132) * pkin(2);
t145 = pkin(2) * qJD(2);
t125 = t86 * t145;
t39 = t81 * pkin(7) + t125;
t152 = t85 * t39;
t21 = t87 * qJD(1) - t152;
t175 = qJD(5) - t21;
t109 = t87 * pkin(4) + t85 * qJ(5);
t128 = t85 * qJDD(1) + t176 * t87;
t129 = qJDD(4) * qJ(5);
t4 = t129 + (qJD(5) - t152) * qJD(4) + t128;
t137 = qJD(4) * t87;
t119 = -t87 * qJDD(1) + t39 * t137 + t176 * t85;
t136 = qJDD(4) * pkin(4);
t173 = qJDD(5) - t136;
t5 = t119 + t173;
t174 = t4 * t87 + t5 * t85;
t22 = t85 * qJD(1) + t87 * t39;
t138 = qJD(4) * t85;
t7 = -t39 * t138 + t128;
t92 = t7 * t87 + (-t21 * t87 - t22 * t85) * qJD(4) + t119 * t85;
t80 = pkin(8) + qJ(2);
t76 = qJ(3) + t80;
t64 = sin(t76);
t65 = cos(t76);
t148 = g(1) * t65 + g(2) * t64;
t15 = -qJD(4) * pkin(4) + t175;
t133 = qJD(4) * qJ(5);
t16 = t22 + t133;
t142 = qJD(2) * t88;
t124 = pkin(2) * t142;
t172 = t177 * pkin(7) - t124 * t179;
t126 = pkin(2) * t140;
t68 = t86 * pkin(2) + pkin(7);
t171 = t126 * t179 + t177 * t68;
t106 = t21 * t85 - t22 * t87;
t164 = t81 * pkin(3);
t40 = -t124 - t164;
t170 = t106 * t88 - t40 * t86;
t139 = qJD(4) * t81;
t146 = -t82 + t83;
t151 = t87 * t79;
t169 = 0.2e1 * t146 * t139 + 0.2e1 * t85 * t151;
t89 = qJD(4) ^ 2;
t168 = pkin(7) * t89;
t59 = g(1) * t64;
t72 = sin(t80);
t167 = g(1) * t72;
t166 = g(2) * t65;
t165 = t79 * pkin(3);
t162 = t88 * pkin(2);
t42 = -pkin(3) - t109;
t160 = t42 * t79;
t159 = t42 * t81;
t158 = t64 * t85;
t157 = t65 * t85;
t156 = t68 * t89;
t155 = t81 * t85;
t150 = g(1) * t158 - g(2) * t157;
t149 = t65 * pkin(3) + t64 * pkin(7);
t147 = -qJD(3) * t125 + qJDD(2) * t162;
t143 = pkin(7) * qJDD(4);
t141 = qJD(3) * t86;
t135 = t16 * qJD(4);
t134 = t85 * qJD(5);
t131 = qJDD(4) * t68;
t116 = t81 * t125;
t49 = t87 * t59;
t127 = t87 * t116 + t124 * t138 + t49;
t123 = t81 * t141;
t23 = -t147 - t165;
t122 = -t23 - t166;
t121 = t21 + t152;
t120 = t40 * t137 + t23 * t85 - t150;
t118 = t109 * t65 + t149;
t114 = t137 * t155;
t56 = t65 * pkin(7);
t113 = g(1) * (-t64 * pkin(3) + t56);
t112 = -t165 + t168;
t73 = cos(t80);
t110 = -g(2) * t73 + t167;
t108 = pkin(4) * t85 - qJ(5) * t87;
t107 = t15 * t85 + t16 * t87;
t105 = t147 + t59 - t166;
t1 = (t108 * qJD(4) - t134) * t81 + t160 - t147;
t104 = -t1 - t160 - t168;
t32 = t42 - t162;
t103 = t32 * t81 - t126;
t27 = pkin(4) * t138 - t87 * t133 - t134;
t101 = g(1) * t157 + g(2) * t158 - g(3) * t87 - t119;
t100 = -t85 * t135 + t15 * t137 - t148 + t174;
t17 = pkin(2) * t141 + t27;
t99 = -t17 * t81 - t32 * t79 - t1 - t156;
t69 = -pkin(3) - t162;
t98 = pkin(2) * t123 + t69 * t79 + t156;
t96 = t22 * qJD(4) + t101;
t95 = -g(1) * t56 - t42 * t59;
t94 = -t131 + (t69 * t81 - t126) * qJD(4);
t93 = (t15 * t87 - t16 * t85) * qJD(4) + t174;
t91 = -t148 + t92;
t84 = qJDD(1) - g(3);
t78 = t81 ^ 2;
t63 = pkin(2) * t73;
t61 = t85 * t79;
t52 = t85 * t78 * t87;
t44 = qJDD(4) * t87 - t89 * t85;
t43 = qJDD(4) * t85 + t89 * t87;
t33 = t146 * t78;
t29 = t40 * t138;
t28 = t108 * t81;
t26 = -0.2e1 * t114 + t153;
t25 = 0.2e1 * t114 + t154;
t13 = -t124 + t159;
t9 = t13 * t138;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, t44, -t43, 0, -t106 * qJD(4) - t119 * t87 + t7 * t85 - g(3), 0, 0, 0, 0, 0, 0, t44, 0, t43, t107 * qJD(4) + t4 * t85 - t5 * t87 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t110, g(1) * t73 + g(2) * t72, 0, 0, 0, 0, 0, 0, 0, t79, (t79 * t88 - t123) * pkin(2) + t105, ((-qJDD(2) - t79) * t86 + (-qJD(2) - t81) * t140) * pkin(2) + t148, 0, (t110 + (t86 ^ 2 + t88 ^ 2) * qJDD(2) * pkin(2)) * pkin(2), t25, t169, t43, t26, t44, 0, t29 + t49 + t94 * t85 + (t122 - t98) * t87, t98 * t85 + t94 * t87 + t120, t91 + t171, t23 * t69 - t113 - g(2) * (t63 + t149) + (-t170 * qJD(3) + t167) * pkin(2) + t92 * t68, t25, t43, -t169, 0, -t44, t26, t49 + t9 + (t103 * qJD(4) - t131) * t85 + (t99 - t166) * t87, t100 + t171, (t131 + (-t103 - t13) * qJD(4)) * t87 + t99 * t85 + t150, t1 * t32 + t13 * t17 - g(2) * (t63 + t118) + (t107 * t140 + t167) * pkin(2) + t93 * t68 + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t105 + t116, (-t132 + (-qJD(3) + t81) * t142) * pkin(2) + t148, 0, 0, t25, t169, t43, t26, t44, 0, t29 + (-pkin(3) * t139 - t143) * t85 + (-t112 + t122) * t87 + t127, (-t143 + (t124 - t164) * qJD(4)) * t87 + (t112 - t116) * t85 + t120, t91 + t172, -t23 * pkin(3) + t92 * pkin(7) - g(2) * t149 + t170 * t145 - t113, t25, t43, -t169, 0, -t44, t26, t9 + (t42 * t139 - t143) * t85 + (-t27 * t81 + t104 - t166) * t87 + t127, t100 + t172, (t143 + (-t124 - t13 - t159) * qJD(4)) * t87 + ((-t27 + t125) * t81 + t104) * t85 + t150, t1 * t42 + t13 * t27 - g(2) * t118 + (-t107 * t88 - t13 * t86) * t145 + t93 * pkin(7) + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t33, t61, t52, t151, qJDD(4), -t40 * t155 + t96, g(3) * t85 + t121 * qJD(4) + (-t40 * t81 + t148) * t87 - t128, 0, 0, -t52, t61, t33, qJDD(4), -t151, t52, 0.2e1 * t136 - qJDD(5) + (-t13 * t85 + t28 * t87) * t81 + t96, -t108 * t79, 0.2e1 * t129 + (t28 * t81 - g(3)) * t85 + (t13 * t81 - t148) * t87 + (0.2e1 * qJD(5) - t121) * qJD(4) + t128, -t5 * pkin(4) - g(3) * t109 + t4 * qJ(5) + t148 * t108 - t13 * t28 - t15 * t22 + t16 * t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t52, t61, -t82 * t78 - t89, t13 * t155 - t101 - t135 + t173;];
tau_reg = t2;
