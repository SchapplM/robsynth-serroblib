% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tau_reg [4x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:39
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:39:14
% EndTime: 2021-01-15 14:39:21
% DurationCPUTime: 1.53s
% Computational Cost: add. (1319->275), mult. (3035->369), div. (0->0), fcn. (1895->6), ass. (0->139)
t86 = sin(qJ(2));
t134 = qJD(3) * t86;
t177 = qJD(1) * t134 - qJDD(2);
t127 = t86 * qJDD(1);
t89 = cos(qJ(2));
t128 = t89 * qJD(1);
t85 = sin(qJ(3));
t88 = cos(qJ(3));
t11 = ((qJD(3) + t128) * qJD(2) + t127) * t85 + t177 * t88;
t126 = qJD(1) * qJD(2);
t119 = t89 * t126;
t130 = t88 * qJD(2);
t10 = -qJD(3) * t130 + (-t119 - t127) * t88 + t177 * t85;
t139 = qJD(1) * t86;
t121 = t85 * t139;
t47 = t121 - t130;
t67 = -qJD(3) + t128;
t160 = t47 * t67;
t176 = -t10 + t160;
t131 = t85 * qJD(2);
t49 = t88 * t139 + t131;
t159 = t49 * t67;
t175 = t11 - t159;
t117 = t86 * t126;
t80 = t89 * qJDD(1);
t174 = -t117 + t80;
t167 = g(3) * t86;
t90 = cos(qJ(1));
t151 = t90 * t88;
t87 = sin(qJ(1));
t155 = t87 * t85;
t31 = t89 * t155 + t151;
t152 = t90 * t85;
t154 = t87 * t88;
t33 = -t89 * t152 + t154;
t173 = -g(1) * t33 + g(2) * t31 + t85 * t167;
t171 = t49 ^ 2;
t54 = -t89 * pkin(2) - t86 * pkin(6) - pkin(1);
t39 = t54 * qJD(1);
t77 = pkin(5) * t128;
t60 = qJD(2) * pkin(6) + t77;
t16 = t88 * t39 - t85 * t60;
t7 = -t49 * qJ(4) + t16;
t6 = -t67 * pkin(3) + t7;
t170 = -t7 + t6;
t166 = g(3) * t89;
t44 = qJDD(3) - t174;
t165 = t44 * pkin(3);
t164 = t47 * pkin(3);
t163 = t85 * pkin(3);
t142 = qJ(4) * t89;
t103 = pkin(3) * t86 - t88 * t142;
t84 = qJ(4) + pkin(6);
t114 = qJD(3) * t84;
t110 = pkin(2) * t86 - pkin(6) * t89;
t51 = t110 * qJD(1);
t146 = pkin(5) * t121 + t88 * t51;
t162 = -t103 * qJD(1) - t85 * qJD(4) - t88 * t114 - t146;
t161 = t10 * t85;
t158 = t49 * t88;
t157 = t85 * t89;
t156 = t86 * t88;
t153 = t88 * t89;
t129 = t88 * qJD(4);
t35 = t85 * t51;
t150 = -t85 * t114 + t129 - t35 - (-pkin(5) * t156 - t85 * t142) * qJD(1);
t133 = qJD(3) * t88;
t52 = t110 * qJD(2);
t149 = t54 * t133 + t85 * t52;
t148 = t86 * pkin(5) * t131 + t88 * t52;
t147 = (g(1) * t151 + g(2) * t154) * t86;
t68 = pkin(5) * t153;
t145 = t85 * t54 + t68;
t82 = t86 ^ 2;
t144 = -t89 ^ 2 + t82;
t143 = qJ(4) * t86;
t141 = t10 * qJ(4);
t140 = t11 * qJ(4);
t138 = qJD(2) * t47;
t137 = qJD(2) * t49;
t136 = qJD(2) * t89;
t135 = qJD(3) * t85;
t59 = -qJD(2) * pkin(2) + pkin(5) * t139;
t132 = t59 * qJD(3);
t125 = t67 * t130;
t124 = t67 * t135;
t123 = t67 * t133;
t115 = -qJD(4) - t164;
t20 = -t115 + t59;
t122 = t20 * t133;
t73 = pkin(5) + t163;
t75 = pkin(5) * t127;
t29 = -qJDD(2) * pkin(2) + pkin(5) * t119 + t75;
t5 = t11 * pkin(3) + qJDD(4) + t29;
t120 = -t5 - t166;
t74 = t88 * pkin(3) + pkin(2);
t116 = t74 * t89 + t84 * t86;
t18 = qJD(1) * t52 + t54 * qJDD(1);
t28 = t174 * pkin(5) + qJDD(2) * pkin(6);
t113 = t39 * t133 - t60 * t135 + t85 * t18 + t88 * t28;
t111 = -qJD(3) * pkin(6) * t67 + t29;
t109 = -g(1) * t31 - g(2) * t33;
t32 = -t87 * t153 + t152;
t34 = t89 * t151 + t155;
t108 = -g(1) * t32 - g(2) * t34;
t107 = g(1) * t90 + g(2) * t87;
t106 = g(1) * t87 - g(2) * t90;
t17 = t85 * t39 + t88 * t60;
t8 = -t47 * qJ(4) + t17;
t105 = t6 * t88 + t8 * t85;
t104 = -pkin(6) * t44 + t132;
t101 = t107 * t86;
t100 = t85 * t44 - t123;
t99 = t88 * t44 + t124;
t98 = -0.2e1 * pkin(1) * t126 - pkin(5) * qJDD(2);
t92 = qJD(1) ^ 2;
t97 = pkin(1) * t92 + t107;
t96 = g(1) * t34 - g(2) * t32 + g(3) * t156 - t113;
t91 = qJD(2) ^ 2;
t95 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t91 + t106;
t14 = t88 * t18;
t94 = -t17 * qJD(3) - t85 * t28 + t14;
t93 = t94 + t173;
t71 = g(3) * t157;
t56 = t84 * t88;
t55 = t84 * t85;
t53 = t73 * t86;
t46 = t88 * t54;
t43 = t47 ^ 2;
t41 = t128 * t163 + t77;
t27 = pkin(1) + t116;
t21 = pkin(5) * t136 + (t89 * t131 + t86 * t133) * pkin(3);
t19 = -t85 * t143 + t145;
t15 = -t88 * t143 + t46 + (-pkin(5) * t85 - pkin(3)) * t89;
t4 = (-pkin(5) * qJD(2) - qJ(4) * qJD(3)) * t156 + (-qJD(4) * t86 + (-pkin(5) * qJD(3) - qJ(4) * qJD(2)) * t89) * t85 + t149;
t3 = -t86 * t129 + t103 * qJD(2) + (-t68 + (-t54 + t143) * t85) * qJD(3) + t148;
t2 = -t47 * qJD(4) + t113 - t140;
t1 = -t49 * qJD(4) + t141 + t165 + t94;
t9 = [qJDD(1), t106, t107, t82 * qJDD(1) + 0.2e1 * t89 * t117, -0.2e1 * t144 * t126 + 0.2e1 * t86 * t80, qJDD(2) * t86 + t91 * t89, qJDD(2) * t89 - t91 * t86, 0, t98 * t86 + t95 * t89, -t95 * t86 + t98 * t89, -t10 * t156 + (t89 * t130 - t85 * t134) * t49, (-t47 * t88 - t49 * t85) * t136 + (t161 - t11 * t88 + (t47 * t85 - t158) * qJD(3)) * t86, (t10 - t125) * t89 + (t99 + t137) * t86, (t67 * t131 + t11) * t89 + (-t100 - t138) * t86, -t67 * qJD(2) * t86 - t44 * t89, -(-t54 * t135 + t148) * t67 + t46 * t44 + (t60 * t133 - t14 + (t123 + t138) * pkin(5) + (-pkin(5) * t44 + qJD(2) * t59 + qJD(3) * t39 + t28) * t85) * t89 + (pkin(5) * t11 + t16 * qJD(2) + t88 * t132 + t29 * t85) * t86 + t108, t149 * t67 - t145 * t44 + (t59 * t130 + (-t124 + t137) * pkin(5) + t113) * t89 + (-t85 * t132 - t17 * qJD(2) + t29 * t88 + (-t10 - t125) * pkin(5)) * t86 + t109, t53 * t11 + t15 * t44 + t21 * t47 - t3 * t67 + (t20 * t131 - t1) * t89 + (qJD(2) * t6 + t5 * t85 + t122) * t86 + t108, -t53 * t10 - t19 * t44 + t21 * t49 + t4 * t67 + (t20 * t130 + t2) * t89 + (-qJD(2) * t8 - t20 * t135 + t5 * t88) * t86 + t109, t15 * t10 - t19 * t11 - t3 * t49 - t4 * t47 - t105 * t136 + (-t1 * t88 - t2 * t85 + (t6 * t85 - t8 * t88) * qJD(3) + t106) * t86, t2 * t19 + t8 * t4 + t1 * t15 + t6 * t3 + t5 * t53 + t20 * t21 - g(1) * (-t27 * t87 + t73 * t90) - g(2) * (t27 * t90 + t73 * t87); 0, 0, 0, -t86 * t92 * t89, t144 * t92, t127, t80, qJDD(2), t97 * t86 - t166 - t75, t167 + (-pkin(5) * qJDD(1) + t97) * t89, -t67 * t158 - t161, -t175 * t85 + t176 * t88, (t67 * t153 - t49 * t86) * qJD(1) + t100, (-t67 * t157 + t47 * t86) * qJD(1) + t99, t67 * t139, -pkin(2) * t11 + t146 * t67 + t104 * t85 + (-t111 - t166) * t88 + (-t16 * t86 + (-pkin(5) * t47 - t59 * t85) * t89) * qJD(1) + t147, pkin(2) * t10 - t35 * t67 + t71 + t104 * t88 + (-t59 * t153 + t17 * t86 + (t67 * t156 - t49 * t89) * pkin(5)) * qJD(1) + (-t101 + t111) * t85, -t6 * t139 - t74 * t11 - t41 * t47 - t55 * t44 + t120 * t88 - t162 * t67 + (-t20 * t128 + (t20 + t164) * qJD(3)) * t85 + t147, t122 + t74 * t10 - t41 * t49 - t56 * t44 + t71 + t150 * t67 + (-t20 * t153 + t8 * t86) * qJD(1) + (pkin(3) * qJD(3) * t49 - t101 + t5) * t85, -t167 - t1 * t85 - t55 * t10 - t56 * t11 + t2 * t88 - t162 * t49 - t150 * t47 - t105 * qJD(3) + (qJD(1) * t105 - t107) * t89, t2 * t56 - t1 * t55 - t5 * t74 - g(3) * t116 + t150 * t8 + t162 * t6 - t107 * (-t86 * t74 + t84 * t89) + (pkin(3) * t135 - t41) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t47, -t43 + t171, -t10 - t160, -t11 - t159, t44, -t17 * t67 - t59 * t49 + t93, -t16 * t67 + t59 * t47 + t96, 0.2e1 * t165 + t141 - t8 * t67 + (t115 - t20) * t49 + t93, -t171 * pkin(3) + t140 - t7 * t67 + (qJD(4) + t20) * t47 + t96, t10 * pkin(3) - t170 * t47, t170 * t8 + (-t20 * t49 + t1 + t173) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, t176, -t43 - t171, t8 * t47 + t6 * t49 - t101 - t120;];
tau_reg = t9;
