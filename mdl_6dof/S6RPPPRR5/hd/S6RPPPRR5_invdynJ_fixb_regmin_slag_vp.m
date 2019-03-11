% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPPRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:45
% EndTime: 2019-03-09 01:37:48
% DurationCPUTime: 1.36s
% Computational Cost: add. (1302->263), mult. (2107->345), div. (0->0), fcn. (1305->8), ass. (0->143)
t86 = sin(qJ(5));
t151 = qJD(6) * t86;
t181 = -qJD(1) * t151 + qJDD(5);
t140 = t86 * qJDD(1);
t89 = cos(qJ(5));
t144 = t89 * qJD(1);
t85 = sin(qJ(6));
t88 = cos(qJ(6));
t11 = ((qJD(6) + t144) * qJD(5) + t140) * t85 - t181 * t88;
t76 = qJD(1) * qJD(2);
t77 = qJ(2) * qJDD(1);
t133 = qJDD(3) + t76 + t77;
t41 = pkin(3) * qJDD(1) + t133;
t84 = -pkin(1) - qJ(3);
t179 = t84 * qJDD(1);
t113 = qJDD(2) + t179;
t138 = qJD(3) * qJD(1);
t42 = t113 - t138;
t81 = sin(pkin(9));
t82 = cos(pkin(9));
t161 = -t41 * t82 + t42 * t81;
t87 = sin(qJ(1));
t90 = cos(qJ(1));
t37 = t81 * t90 + t82 * t87;
t173 = g(2) * t37;
t36 = -t81 * t87 + t82 * t90;
t175 = g(1) * t36;
t49 = t82 * qJD(2) + t81 * qJD(3);
t180 = qJD(1) * t49 - t161 - t173 - t175;
t145 = t88 * qJD(5);
t103 = t145 * t89 - t151 * t85;
t139 = qJD(1) * qJD(5);
t128 = t86 * t139;
t66 = t89 * qJDD(1);
t35 = qJDD(6) - t66 + t128;
t163 = t88 * t35;
t57 = -qJD(6) + t144;
t178 = -t103 * t57 + t163 * t86;
t174 = g(1) * t37;
t120 = -g(2) * t36 + t174;
t121 = pkin(5) * t86 - pkin(8) * t89;
t63 = qJ(2) * qJD(1) + qJD(3);
t54 = pkin(3) * qJD(1) + t63;
t55 = qJD(1) * t84 + qJD(2);
t23 = t54 * t81 + t55 * t82;
t21 = qJD(1) * pkin(7) + t23;
t15 = qJD(4) * t86 + t21 * t89;
t147 = t15 * qJD(5);
t19 = t41 * t81 + t42 * t82;
t17 = qJDD(1) * pkin(7) + t19;
t2 = -qJDD(5) * pkin(5) - t89 * qJDD(4) + t86 * t17 + t147;
t176 = t120 * t86 + (pkin(8) * qJD(6) + qJD(1) * t121) * t57 - g(3) * t89 - t2;
t71 = 0.2e1 * t76;
t130 = t89 * t139;
t10 = qJD(6) * t145 + (t130 + t140) * t88 + t181 * t85;
t172 = t10 * t85;
t171 = t36 * t89;
t156 = qJD(1) * t86;
t38 = t156 * t85 - t145;
t170 = t38 * t57;
t146 = t85 * qJD(5);
t40 = t156 * t88 + t146;
t169 = t40 * t57;
t168 = t57 * t88;
t167 = t85 * t35;
t166 = t85 * t89;
t165 = t86 * t38;
t164 = t86 * t40;
t162 = t88 * t89;
t83 = pkin(3) + qJ(2);
t34 = t81 * t83 + t82 * t84;
t160 = pkin(1) * t90 + qJ(2) * t87;
t159 = g(1) * t87 - g(2) * t90;
t79 = t86 ^ 2;
t158 = -t89 ^ 2 + t79;
t91 = qJD(5) ^ 2;
t92 = qJD(1) ^ 2;
t157 = t91 + t92;
t155 = qJD(5) * t38;
t154 = qJD(5) * t40;
t153 = qJD(5) * t86;
t152 = qJD(6) * t85;
t150 = qJD(6) * t88;
t149 = qJDD(1) * pkin(1);
t14 = qJD(4) * t89 - t21 * t86;
t148 = t14 * qJD(5);
t143 = qJDD(4) - g(3);
t142 = qJDD(5) * t86;
t141 = qJDD(5) * t89;
t136 = qJ(3) * t90 + t160;
t135 = t57 * t146;
t132 = -qJDD(2) + t159;
t32 = pkin(7) + t34;
t8 = qJD(5) * pkin(8) + t15;
t131 = t32 * t57 + t8;
t114 = -pkin(5) * t89 - pkin(8) * t86 - pkin(4);
t22 = t82 * t54 - t55 * t81;
t9 = qJD(1) * t114 - t22;
t127 = -qJDD(5) * pkin(8) - qJD(6) * t9 - t86 * qJDD(4) - t89 * t17 - t148;
t126 = -t10 * t89 + t153 * t40;
t33 = -t81 * t84 + t82 * t83;
t124 = 0.2e1 * t128;
t122 = qJDD(2) - t149;
t119 = g(1) * t90 + g(2) * t87;
t118 = -t11 + t135;
t70 = t90 * qJ(2);
t117 = t84 * t87 + t70;
t116 = -qJD(6) * t8 + t173;
t115 = qJD(5) * t21 - t143;
t111 = t150 * t57 - t167;
t110 = t152 * t57 + t163;
t109 = t121 * qJD(5);
t108 = -t119 + t71 + 0.2e1 * t77;
t107 = -t157 * t89 - t142;
t106 = t157 * t86 - t141;
t105 = 0.2e1 * t130 + t140;
t104 = t124 - t66;
t101 = t111 + t155;
t100 = -g(2) * t171 + g(3) * t86 + t127;
t20 = -qJD(1) * pkin(4) - t22;
t31 = -pkin(4) - t33;
t50 = qJD(2) * t81 - qJD(3) * t82;
t99 = -qJDD(5) * t32 + (qJD(1) * t31 + t20 - t50) * qJD(5);
t7 = -qJD(5) * pkin(5) - t14;
t98 = -pkin(8) * t35 + (-t14 - t7) * t57;
t97 = qJD(5) * t7 - t32 * t35 + t50 * t57 - t127;
t96 = -qJD(1) * t20 - qJD(4) * qJD(5) + t120 - t17;
t95 = t32 * t91 - t180 + (-pkin(4) + t31) * qJDD(1);
t94 = (-t145 * t57 + t10) * t86 + (-t110 + t154) * t89;
t93 = t101 * t89 - t118 * t86;
t52 = -t86 * t91 + t141;
t51 = t89 * t91 + t142;
t48 = qJDD(1) * t82 - t81 * t92;
t47 = -qJDD(1) * t81 - t82 * t92;
t25 = t109 - t49;
t24 = t114 - t33;
t13 = t162 * t37 - t36 * t85;
t12 = -t166 * t37 - t36 * t88;
t6 = qJD(1) * t109 + qJDD(1) * t114 + t161;
t5 = t88 * t6;
t4 = t8 * t88 + t85 * t9;
t3 = -t8 * t85 + t88 * t9;
t1 = [qJDD(1), t159, t119, -t132 - 0.2e1 * t149, t108, -t122 * pkin(1) - g(1) * (-pkin(1) * t87 + t70) - g(2) * t160 + (t71 + t77) * qJ(2), qJDD(3) + t108, t132 + 0.2e1 * t138 - 0.2e1 * t179, -g(1) * t117 - g(2) * t136 + qJ(2) * t133 + t63 * qJD(2) - t55 * qJD(3) + t42 * t84, qJDD(1) * t33 + t180, -qJD(1) * t50 - qJDD(1) * t34 + t120 - t19, t19 * t34 + t23 * t50 - t161 * t33 + t22 * t49 - g(1) * (pkin(3) * t90 + t117) - g(2) * (pkin(3) * t87 + t136) qJDD(1) * t79 + t124 * t89, -0.2e1 * t139 * t158 + 0.2e1 * t66 * t86, t51, t52, 0, t86 * t99 - t89 * t95, t86 * t95 + t89 * t99, t10 * t88 * t86 + t103 * t40 (-t38 * t88 - t40 * t85) * t89 * qJD(5) + (-t172 - t11 * t88 + (t38 * t85 - t40 * t88) * qJD(6)) * t86, t126 + t178 (t11 + t135) * t89 + (t111 - t155) * t86, -t153 * t57 - t35 * t89, -g(2) * t13 + (t155 * t32 - t5) * t89 + (qJD(5) * t3 + t11 * t32 + t38 * t50) * t86 + (-g(1) * t171 + t24 * t35 - t25 * t57 + (t131 * t89 + t7 * t86) * qJD(6)) * t88 + (-(-qJD(6) * t24 + t153 * t32) * t57 + t2 * t86 - t174 + t97 * t89) * t85 (t150 * t24 + t85 * t25) * t57 - t24 * t167 - t88 * t174 - g(2) * t12 + (t32 * t154 + (-qJD(6) * t131 + t175 + t6) * t85 + t97 * t88) * t89 + (-t7 * t152 + t32 * t10 + t2 * t88 + t50 * t40 + (-t168 * t32 - t4) * qJD(5)) * t86; 0, 0, 0, qJDD(1), -t92, -qJ(2) * t92 + t122 - t159, -t92, -qJDD(1) (-qJD(3) - t63) * qJD(1) + t113 - t159, t47, -t48, t161 * t81 + t19 * t82 + (-t22 * t82 - t23 * t81) * qJD(1) - t159, 0, 0, 0, 0, 0, t104 * t81 + t107 * t82, t105 * t81 + t106 * t82, 0, 0, 0, 0, 0, t110 * t81 + t93 * t82 + ((-t166 * t81 - t82 * t88) * t57 - t81 * t165) * qJD(1), t111 * t81 + t94 * t82 + (-(t162 * t81 - t82 * t85) * t57 - t81 * t164) * qJD(1); 0, 0, 0, 0, 0, 0, qJDD(1), -t92, qJD(1) * t55 - t119 + t133, t48, t47, -t161 * t82 + t19 * t81 + (-t22 * t81 + t23 * t82) * qJD(1) - t119, 0, 0, 0, 0, 0, -t104 * t82 + t107 * t81, -t105 * t82 + t106 * t81, 0, 0, 0, 0, 0, -t110 * t82 + t93 * t81 + ((t166 * t82 - t81 * t88) * t57 + t82 * t165) * qJD(1), -t111 * t82 + t94 * t81 + (-(-t162 * t82 - t81 * t85) * t57 + t82 * t164) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, 0, 0, 0, 0, 0, t52, -t51, 0, 0, 0, 0, 0, t101 * t86 + t118 * t89, t126 - t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86 * t92 * t89, t158 * t92, t140, t66, qJDD(5), -t115 * t89 + t86 * t96 + t147, t115 * t86 + t89 * t96 + t148, -t168 * t40 + t172 (t10 + t170) * t88 + (-t11 + t169) * t85 (t162 * t57 - t164) * qJD(1) - t111 (-t166 * t57 + t165) * qJD(1) + t110, t57 * t156, -pkin(5) * t11 - t15 * t38 - t156 * t3 + t176 * t88 + t85 * t98, -pkin(5) * t10 - t15 * t40 + t156 * t4 - t176 * t85 + t88 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 * t38, -t38 ^ 2 + t40 ^ 2, t10 - t170, -t11 - t169, t35, -g(1) * t12 + t100 * t85 + t116 * t88 - t4 * t57 - t7 * t40 + t5, g(1) * t13 - t3 * t57 + t7 * t38 + (-t116 - t6) * t85 + t100 * t88;];
tau_reg  = t1;
