% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPP2
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:52:01
% EndTime: 2019-12-31 20:52:05
% DurationCPUTime: 1.77s
% Computational Cost: add. (1667->316), mult. (2337->342), div. (0->0), fcn. (1173->8), ass. (0->183)
t116 = cos(qJ(3));
t190 = qJ(4) * t116;
t113 = sin(qJ(3));
t211 = pkin(3) + pkin(4);
t212 = t211 * t113;
t217 = t190 - t212;
t97 = t113 * qJ(4);
t154 = pkin(2) + t97;
t106 = qJD(1) + qJD(2);
t175 = qJD(3) * t116;
t157 = t106 * t175;
t105 = qJDD(1) + qJDD(2);
t78 = t113 * t105;
t216 = t157 + t78;
t160 = t211 * qJD(3);
t186 = t106 * t113;
t114 = sin(qJ(2));
t201 = pkin(1) * qJD(1);
t166 = t114 * t201;
t53 = t106 * pkin(7) + t166;
t44 = t113 * t53;
t24 = qJ(5) * t186 - t44;
t184 = qJD(4) - t24;
t15 = -t160 + t184;
t156 = t211 * qJDD(3);
t172 = qJDD(1) * t114;
t117 = cos(qJ(2));
t178 = qJD(2) * t117;
t33 = t105 * pkin(7) + (qJD(1) * t178 + t172) * pkin(1);
t28 = t113 * t33;
t41 = t53 * t175;
t162 = qJDD(4) + t28 + t41;
t187 = qJDD(3) * pkin(3);
t10 = t162 - t187;
t107 = qJDD(3) * qJ(4);
t108 = qJD(3) * qJD(4);
t176 = qJD(3) * t113;
t29 = t116 * t33;
t9 = -t53 * t176 + t107 + t108 + t29;
t215 = t10 * t113 + t9 * t116;
t112 = qJ(1) + qJ(2);
t95 = sin(t112);
t96 = cos(t112);
t203 = g(1) * t96 + g(2) * t95;
t109 = qJD(3) * qJ(4);
t185 = t106 * t116;
t45 = t116 * t53;
t25 = -qJ(5) * t185 + t45;
t21 = t109 + t25;
t210 = g(2) * t96;
t85 = g(1) * t95;
t214 = t85 - t210;
t180 = qJD(1) * t117;
t165 = pkin(1) * t180;
t12 = t165 + qJD(5) + (t211 * t116 + t154) * t106;
t183 = qJD(5) + t12;
t213 = t183 * t113;
t120 = qJD(3) ^ 2;
t209 = pkin(7) * t120;
t115 = sin(qJ(1));
t208 = g(1) * t115;
t93 = t105 * pkin(2);
t207 = t106 * pkin(2);
t101 = t116 * pkin(3);
t100 = t116 * pkin(4);
t206 = t117 * pkin(1);
t205 = pkin(7) - qJ(5);
t198 = t113 * t96;
t199 = t113 * t95;
t204 = g(1) * t199 - g(2) * t198;
t202 = -qJD(2) * t166 + qJDD(1) * t206;
t192 = t101 + t97;
t164 = t100 + t192;
t47 = pkin(2) + t164;
t200 = t106 * t47;
t197 = t116 * t21;
t196 = t116 * t96;
t89 = t114 * pkin(1) + pkin(7);
t195 = t120 * t89;
t148 = -qJD(3) * pkin(3) + qJD(4);
t31 = t148 + t44;
t194 = t31 * t116;
t193 = -qJ(5) + t89;
t191 = pkin(7) * qJDD(3);
t189 = qJ(5) * t105;
t49 = t193 * t116;
t188 = qJD(3) * t49;
t79 = t116 * t105;
t110 = t113 ^ 2;
t111 = t116 ^ 2;
t182 = t110 - t111;
t181 = t110 + t111;
t179 = qJD(2) * t114;
t177 = qJD(3) * t106;
t174 = qJDD(3) * t89;
t94 = t113 * qJD(4);
t173 = t113 * qJD(5);
t171 = t21 * t176 + t203;
t170 = t29 + 0.2e1 * t107 + 0.2e1 * t108;
t146 = t106 * t166;
t72 = t116 * t85;
t169 = t116 * t146 + t165 * t176 + t72;
t32 = -t93 - t202;
t39 = pkin(3) * t176 - qJ(4) * t175 - t94;
t168 = pkin(1) * t178;
t167 = pkin(1) * t179;
t104 = t106 ^ 2;
t163 = t113 * t104 * t116;
t90 = -pkin(2) - t206;
t161 = -t32 - t210;
t159 = t106 * t179;
t158 = t106 * t176;
t87 = qJ(5) * t176;
t61 = t205 * t116;
t131 = pkin(3) * t79 + t216 * qJ(4) + t106 * t94 - t32;
t127 = pkin(4) * t79 + qJDD(5) + t131;
t3 = -t211 * t158 + t127;
t155 = t3 * t113 + t12 * t175 + t204;
t82 = t96 * pkin(7);
t153 = -t96 * qJ(5) + t82;
t54 = -t165 - t207;
t152 = t32 * t113 + t54 * t175 - t204;
t151 = pkin(3) * t196 + t95 * pkin(7) + t154 * t96;
t150 = t181 * t105;
t149 = g(1) * t198 + g(2) * t199 - g(3) * t116 - t28;
t147 = 0.2e1 * t157;
t145 = t181 * t206;
t118 = cos(qJ(1));
t143 = t118 * pkin(1) + t151;
t142 = -qJD(5) + t168;
t141 = -t93 + t209;
t140 = t202 + t214;
t139 = -qJDD(4) + t149;
t138 = pkin(3) * t113 - t190;
t26 = -pkin(4) * t176 - t39;
t19 = t26 - t167;
t46 = t90 - t192;
t34 = t100 - t46;
t137 = t105 * t34 + t106 * t19;
t136 = t105 * t47 + t106 * t26;
t36 = t45 + t109;
t135 = t113 * t31 + t116 * t36;
t134 = -t154 - t101;
t133 = -t41 + t139;
t132 = t106 * t46 - t168;
t130 = t31 * t175 - t36 * t176 - t203 + t215;
t55 = -pkin(2) - t192;
t6 = pkin(3) * t158 - t131;
t129 = -t105 * t55 - t106 * t39 - t209 - t6;
t27 = t167 + t39;
t128 = -t105 * t46 - t106 * t27 - t195 - t6;
t126 = pkin(1) * t159 + t105 * t90 + t195;
t125 = -qJ(5) * t78 - t133;
t124 = -g(1) * t82 - t134 * t85;
t123 = -t174 + (t106 * t90 - t168) * qJD(3);
t122 = (-t36 * t113 + t194) * qJD(3) + t215;
t121 = (-g(1) * (t134 - t100) + g(2) * qJ(5)) * t95;
t67 = pkin(4) * t196;
t66 = t96 * t190;
t64 = t95 * t190;
t62 = t106 * t87;
t60 = t205 * t113;
t58 = -t110 * t104 - t120;
t57 = qJDD(3) * t116 - t120 * t113;
t56 = qJDD(3) * t113 + t120 * t116;
t51 = t113 * t146;
t50 = -qJDD(3) - t163;
t48 = t193 * t113;
t42 = t54 * t176;
t40 = t138 * t106;
t38 = qJD(3) * t61 - t173;
t37 = -pkin(7) * t176 - t116 * qJD(5) + t87;
t35 = t110 * t105 + t113 * t147;
t23 = t217 * t106;
t20 = t134 * t106 - t165;
t18 = 0.2e1 * t113 * t79 - 0.2e1 * t182 * t177;
t17 = t142 * t113 + t188;
t16 = t142 * t116 - t89 * t176 + t87;
t13 = t20 * t176;
t5 = t62 + (-qJD(5) * t106 - t189) * t116 + t9;
t4 = -t216 * qJ(5) - t106 * t173 - t156 + t162;
t2 = t3 * t116;
t1 = [qJDD(1), -g(2) * t118 + t208, g(1) * t118 + g(2) * t115, t105, (t105 * t117 - t159) * pkin(1) + t140, ((-qJDD(1) - t105) * t114 + (-qJD(1) - t106) * t178) * pkin(1) + t203, t35, t18, t56, t57, 0, t42 + t72 + t123 * t113 + (-t126 + t161) * t116, t113 * t126 + t116 * t123 + t152, t13 + t72 + (qJD(3) * t132 - t174) * t113 + (t128 - t210) * t116, qJD(2) * t106 * t145 + t150 * t89 + t130, (t174 + (-t132 - t20) * qJD(3)) * t116 + t128 * t113 + t204, t6 * t46 + t20 * t27 - g(2) * t143 + (t135 * t178 + t208) * pkin(1) + t122 * t89 + t124, -t48 * qJDD(3) + t2 + t72 + (t137 - t210) * t116 + (-t17 + (-t106 * t34 - t12) * t113) * qJD(3), t49 * qJDD(3) + t137 * t113 + (t34 * t185 + t16) * qJD(3) + t155, (-t105 * t48 - t4 + (-t17 + t188) * t106) * t113 + (-t105 * t49 - t106 * t16 - t5 + (-t106 * t48 - t15) * qJD(3)) * t116 + t171, t5 * t49 + t21 * t16 + t4 * t48 + t15 * t17 + t3 * t34 + t12 * t19 - g(1) * (-t115 * pkin(1) + t153) - g(2) * (t67 + t143) + t121; 0, 0, 0, t105, t140 + t146, (-t172 + (-qJD(2) + t106) * t180) * pkin(1) + t203, t35, t18, t56, t57, 0, t42 + (-pkin(2) * t177 - t191) * t113 + (-t141 + t161) * t116 + t169, -t51 + t141 * t113 + (-t191 + (t165 - t207) * qJD(3)) * t116 + t152, t13 + (t55 * t177 - t191) * t113 + (t129 - t210) * t116 + t169, -t181 * t106 * t165 + pkin(7) * t150 + t130, t51 + (t191 + (-t106 * t55 - t165 - t20) * qJD(3)) * t116 + t129 * t113 + t204, t6 * t55 + t20 * t39 - g(2) * t151 + (-t114 * t20 - t117 * t135) * t201 + t122 * pkin(7) + t124, -t60 * qJDD(3) + t2 + (t136 - t210) * t116 + (-t38 + (-t12 - t200) * t113) * qJD(3) + t169, t61 * qJDD(3) + t51 + t136 * t113 + (t37 + (-t165 + t200) * t116) * qJD(3) + t155, (-t105 * t60 - t4) * t113 + (-qJD(3) * t15 - t105 * t61 - t5) * t116 + (-t113 * t38 - t116 * t37 + (t113 * t61 - t116 * t60) * qJD(3) + qJD(1) * t145) * t106 + t171, t5 * t61 + t21 * t37 + t4 * t60 + t15 * t38 + t3 * t47 + t12 * t26 - g(1) * t153 - g(2) * (t67 + t151) + t121 + (t114 * t12 + (-t113 * t15 - t197) * t117) * t201; 0, 0, 0, 0, 0, 0, -t163, t182 * t104, t78, t79, qJDD(3), -t54 * t186 + t149, g(3) * t113 - t29 + (-t106 * t54 + t203) * t116, 0.2e1 * t187 + (-t113 * t20 + t116 * t40) * t106 + t139, -t138 * t105 + ((t36 - t109) * t113 + (t148 - t31) * t116) * t106, (t106 * t40 - g(3)) * t113 + (t106 * t20 - t203) * t116 + t170, t9 * qJ(4) - t10 * pkin(3) - t20 * t40 - t53 * t194 - g(1) * (-pkin(3) * t198 + t66) - g(2) * (-pkin(3) * t199 + t64) - g(3) * t192 + (qJD(4) + t44) * t36, t25 * qJD(3) + 0.2e1 * t156 + ((qJ(5) * qJD(3) - t23) * t116 + t213) * t106 - t125, -t24 * qJD(3) + t62 + (-qJD(3) * t53 - t106 * t23 - g(3)) * t113 + (-t183 * t106 - t189 - t203) * t116 + t170, -t217 * t105, -g(1) * t66 - g(2) * t64 - g(3) * t164 + t5 * qJ(4) - t12 * t23 - t15 * t25 + t184 * t21 + t203 * t212 - t211 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t78, t58, -qJD(3) * t36 + t20 * t186 - t133 - t187, t50, t58, -t78, -t21 * qJD(3) - t156 + (-qJ(5) * t175 - t213) * t106 + t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79 - 0.2e1 * t158, t78 + t147, -t181 * t104, (t197 + (t15 - t160) * t113) * t106 + t127 + t214;];
tau_reg = t1;
