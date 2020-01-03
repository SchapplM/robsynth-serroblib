% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:43
% EndTime: 2019-12-31 20:15:46
% DurationCPUTime: 1.00s
% Computational Cost: add. (1233->207), mult. (1639->266), div. (0->0), fcn. (1030->12), ass. (0->138)
t101 = qJ(1) + qJ(2);
t89 = sin(t101);
t91 = cos(t101);
t176 = g(1) * t89 - g(2) * t91;
t102 = sin(qJ(5));
t103 = sin(qJ(4));
t147 = qJD(5) * t102;
t149 = qJD(4) * t103;
t177 = -t102 * t149 - t103 * t147;
t106 = cos(qJ(5));
t107 = cos(qJ(4));
t50 = -t102 * t103 + t106 * t107;
t130 = -g(1) * t91 - g(2) * t89;
t92 = t103 * pkin(4);
t74 = qJ(3) + t92;
t96 = qJD(4) + qJD(5);
t104 = sin(qJ(2));
t165 = pkin(1) * qJD(1);
t142 = t104 * t165;
t97 = qJD(1) + qJD(2);
t158 = t97 * qJ(3);
t51 = t142 + t158;
t119 = -t51 * t97 - t176;
t93 = t97 ^ 2;
t110 = -pkin(2) - pkin(7);
t95 = qJDD(1) + qJDD(2);
t175 = t95 * pkin(2);
t108 = cos(qJ(2));
t82 = -t108 * pkin(1) - pkin(2);
t70 = -pkin(7) + t82;
t174 = -pkin(8) + t70;
t35 = t50 * t97;
t154 = t102 * t107;
t49 = t106 * t103 + t154;
t36 = t49 * t97;
t173 = t35 * t36;
t171 = -pkin(8) + t110;
t21 = t96 * t49;
t94 = qJDD(4) + qJDD(5);
t13 = -t21 * t96 + t50 * t94;
t148 = qJD(4) * t107;
t159 = t95 * qJ(3);
t140 = qJD(2) * t165;
t156 = pkin(1) * qJDD(1);
t169 = t104 * t156 + t108 * t140;
t32 = t97 * qJD(3) + t159 + t169;
t170 = t32 * t103 + t51 * t148;
t168 = t91 * pkin(2) + t89 * qJ(3);
t167 = -t104 * t140 + t108 * t156;
t99 = t107 ^ 2;
t166 = t103 ^ 2 - t99;
t164 = t103 * t95;
t163 = t104 * t96;
t141 = t108 * t165;
t129 = qJD(3) - t141;
t40 = t110 * t97 + t129;
t23 = (-pkin(8) * t97 + t40) * t103;
t162 = t106 * t23;
t161 = t107 * t95;
t160 = t107 * t97;
t111 = qJD(4) ^ 2;
t157 = -t111 - t93;
t153 = t103 * t107;
t151 = qJD(2) * t104;
t150 = qJD(2) * t108;
t146 = qJDD(4) * t70;
t145 = qJDD(4) * t103;
t144 = qJDD(4) * t110;
t143 = pkin(1) * t151;
t139 = t97 * t151;
t138 = t97 * t148;
t137 = qJDD(3) - t167;
t136 = -t89 * pkin(2) + t91 * qJ(3);
t44 = t174 * t107;
t55 = t171 * t107;
t73 = t104 * pkin(1) + qJ(3);
t133 = -t169 - t130;
t132 = t167 + t176;
t131 = t96 * t107;
t24 = -pkin(8) * t160 + t107 * t40;
t63 = pkin(1) * t150 + qJD(3);
t22 = t106 * t131 + t177;
t14 = -t22 * t96 - t49 * t94;
t20 = qJD(4) * pkin(4) + t24;
t128 = -t102 * t20 - t162;
t43 = t174 * t103;
t127 = t102 * t44 + t106 * t43;
t126 = -t102 * t43 + t106 * t44;
t54 = t171 * t103;
t125 = t102 * t55 + t106 * t54;
t124 = -t102 * t54 + t106 * t55;
t123 = t95 * t154 + t177 * t97;
t34 = t137 - t175;
t122 = t73 * t97 + t143;
t121 = -t142 + t158;
t120 = t138 + t164;
t118 = t97 * t142 + t132;
t17 = t120 * pkin(4) + t32;
t38 = t74 * t97 + t142;
t100 = qJ(4) + qJ(5);
t88 = sin(t100);
t117 = t130 * t88 + t17 * t49 + t38 * t22;
t90 = cos(t100);
t116 = t130 * t90 + t17 * t50 - t38 * t21;
t115 = -t111 * t70 + t63 * t97 + t73 * t95 + t130;
t29 = t110 * t95 + t137;
t26 = t107 * t29;
t6 = -t40 * t149 + qJDD(4) * pkin(4) + t26 + (t97 * t149 - t161) * pkin(8);
t114 = t23 * t147 + g(3) * t90 + (-t23 * t96 - t6) * t102 + t38 * t36 + t176 * t88;
t113 = -t110 * t111 + t129 * t97 + t130 + t159;
t7 = -t120 * pkin(8) + t103 * t29 + t40 * t148;
t112 = g(3) * t88 + t128 * qJD(5) - t102 * t7 + t106 * t6 - t176 * t90 - t38 * t35;
t8 = -t21 * t97 + t50 * t95;
t109 = cos(qJ(1));
t105 = sin(qJ(1));
t87 = qJDD(4) * t107;
t86 = pkin(4) * t148;
t85 = pkin(8) * t149;
t64 = qJD(3) + t86;
t60 = -t111 * t103 + t87;
t59 = -t111 * t107 - t145;
t58 = t73 + t92;
t48 = -t97 * pkin(2) + t129;
t47 = t63 + t86;
t46 = qJD(4) * t55;
t45 = -t110 * t149 + t85;
t39 = -0.2e1 * t103 * t138 + t99 * t95;
t31 = qJD(4) * t44 + t103 * t143;
t30 = t107 * t143 - t70 * t149 + t85;
t28 = t32 * t107;
t25 = 0.2e1 * t166 * t97 * qJD(4) - 0.2e1 * t95 * t153;
t10 = t35 ^ 2 - t36 ^ 2;
t9 = (t97 * t131 + t164) * t106 + t123;
t4 = t35 * t96 + (-t96 * t160 - t164) * t106 - t123;
t3 = t36 * t96 + t8;
t2 = -t35 * t21 + t8 * t50;
t1 = t21 * t36 - t35 * t22 - t8 * t49 - t50 * t9;
t5 = [qJDD(1), g(1) * t105 - g(2) * t109, g(1) * t109 + g(2) * t105, t95, (t108 * t95 - t139) * pkin(1) + t132, (-t104 * t95 - t97 * t150) * pkin(1) + t133, pkin(1) * t139 + qJDD(3) + (-pkin(2) + t82) * t95 - t132, (qJD(3) + t63) * t97 + (qJ(3) + t73) * t95 - t133, t32 * t73 + t51 * t63 + t34 * t82 + t48 * t143 - g(1) * (-t105 * pkin(1) + t136) - g(2) * (t109 * pkin(1) + t168), t39, t25, t60, t59, 0, (t122 * qJD(4) + t146) * t107 + t115 * t103 + t170, t28 + (-t146 + (-t122 - t51) * qJD(4)) * t103 + t115 * t107, t2, t1, t13, t14, 0, t47 * t36 + t58 * t9 + (-t127 * qJD(5) - t102 * t31 + t106 * t30) * t96 + t126 * t94 + t117, t47 * t35 + t58 * t8 - (t126 * qJD(5) + t102 * t30 + t106 * t31) * t96 - t127 * t94 + t116; 0, 0, 0, t95, t118, t97 * t141 + t133, qJDD(3) - t118 - 0.2e1 * t175, 0.2e1 * t159 + (0.2e1 * qJD(3) - t141) * t97 - t133, t32 * qJ(3) + t51 * qJD(3) - t34 * pkin(2) - g(1) * t136 - g(2) * t168 + (-t104 * t48 - t108 * t51) * t165, t39, t25, t60, t59, 0, (t121 * qJD(4) + t144) * t107 + t113 * t103 + t170, t28 + (-t144 + (-t121 - t51) * qJD(4)) * t103 + t113 * t107, t2, t1, t13, t14, 0, t64 * t36 + t74 * t9 + (-t125 * qJD(5) - t102 * t46 + t106 * t45) * t96 + t124 * t94 + (-t108 * t36 - t50 * t163) * t165 + t117, t64 * t35 + t74 * t8 - (t124 * qJD(5) + t102 * t45 + t106 * t46) * t96 - t125 * t94 + (-t108 * t35 + t49 * t163) * t165 + t116; 0, 0, 0, 0, 0, 0, t95, -t93, t34 + t119, 0, 0, 0, 0, 0, t157 * t103 + t87, t157 * t107 - t145, 0, 0, 0, 0, 0, -t36 * t97 + t13, -t35 * t97 + t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, t93 * t153, -t166 * t93, t161, -t164, qJDD(4), g(3) * t103 + t119 * t107 + t26, g(3) * t107 + (-t119 - t29) * t103, t173, t10, t3, t4, t94, -(-t102 * t24 - t162) * t96 + (t106 * t94 - t96 * t147 - t36 * t160) * pkin(4) + t112, (-qJD(5) * t20 + t24 * t96 - t7) * t106 + (-qJD(5) * t106 * t96 - t102 * t94 - t35 * t160) * pkin(4) + t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t10, t3, t4, t94, -t128 * t96 + t112, (-t7 + (-qJD(5) + t96) * t20) * t106 + t114;];
tau_reg = t5;
