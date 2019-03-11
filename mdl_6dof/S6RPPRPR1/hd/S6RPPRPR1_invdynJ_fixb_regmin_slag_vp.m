% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:40:06
% EndTime: 2019-03-09 01:40:12
% DurationCPUTime: 2.98s
% Computational Cost: add. (4066->359), mult. (8888->469), div. (0->0), fcn. (6835->18), ass. (0->189)
t164 = pkin(10) + qJ(4);
t155 = sin(t164);
t158 = cos(t164);
t165 = qJ(1) + pkin(9);
t156 = sin(t165);
t159 = cos(t165);
t212 = g(1) * t159 + g(2) * t156;
t185 = -g(3) * t158 + t155 * t212;
t174 = sin(qJ(4));
t254 = cos(qJ(4));
t168 = sin(pkin(9));
t141 = pkin(1) * t168 + qJ(3);
t135 = t141 * qJD(1);
t170 = cos(pkin(10));
t153 = t170 * qJD(2);
t167 = sin(pkin(10));
t244 = pkin(7) * qJD(1);
t104 = t153 + (-t135 - t244) * t167;
t110 = t167 * qJD(2) + t170 * t135;
t105 = t170 * t244 + t110;
t51 = t174 * t104 + t254 * t105;
t262 = t51 * qJD(4);
t125 = qJD(1) * qJD(3) + t141 * qJDD(1);
t150 = t170 * qJDD(2);
t90 = t150 + (-pkin(7) * qJDD(1) - t125) * t167;
t107 = t167 * qJDD(2) + t170 * t125;
t222 = t170 * qJDD(1);
t91 = pkin(7) * t222 + t107;
t192 = t174 * t91 - t254 * t90 + t262;
t19 = -qJDD(4) * pkin(4) + qJDD(5) + t192;
t180 = t19 - t185;
t171 = cos(pkin(9));
t253 = pkin(1) * t171;
t146 = -pkin(2) - t253;
t224 = qJDD(1) * t146;
t131 = qJDD(3) + t224;
t216 = -g(1) * t156 + g(2) * t159;
t261 = -t131 - t216;
t218 = t254 * t170;
t140 = qJD(1) * t218;
t230 = t174 * t167;
t217 = qJD(1) * t230;
t116 = -t140 + t217;
t112 = qJD(6) + t116;
t130 = t254 * t167 + t174 * t170;
t118 = t130 * qJD(1);
t166 = sin(pkin(11));
t169 = cos(pkin(11));
t100 = -t169 * qJD(4) + t118 * t166;
t102 = qJD(4) * t166 + t118 * t169;
t173 = sin(qJ(6));
t176 = cos(qJ(6));
t47 = t176 * t100 + t102 * t173;
t260 = t112 * t47;
t129 = t166 * t176 + t169 * t173;
t120 = t129 * qJD(6);
t237 = t129 * t116 + t120;
t198 = t100 * t173 - t102 * t176;
t259 = t112 * t198;
t50 = t254 * t104 - t174 * t105;
t258 = t216 * t155;
t193 = t216 * t158;
t257 = -qJD(6) + t112;
t127 = t166 * t173 - t176 * t169;
t238 = t112 * t127;
t122 = t130 * qJD(4);
t214 = qJDD(1) * t254;
t223 = t167 * qJDD(1);
t205 = -t170 * t214 + t174 * t223;
t82 = qJD(1) * t122 + t205;
t80 = qJDD(6) + t82;
t256 = t112 * t238 - t129 * t80;
t249 = g(3) * t155;
t183 = -t212 * t158 - t249;
t221 = qJD(4) * t140 + t167 * t214 + t174 * t222;
t81 = -qJD(4) * t217 + t221;
t69 = -t169 * qJDD(4) + t166 * t81;
t70 = qJDD(4) * t166 + t169 * t81;
t12 = -qJD(6) * t198 + t173 * t70 + t176 * t69;
t115 = t116 ^ 2;
t226 = qJD(6) * t176;
t227 = qJD(6) * t173;
t11 = -t100 * t226 - t102 * t227 - t173 * t69 + t176 * t70;
t189 = t218 - t230;
t255 = -t11 * t189 - t122 * t198;
t252 = pkin(8) * t169;
t175 = sin(qJ(1));
t250 = g(1) * t175;
t247 = pkin(7) + t141;
t246 = pkin(8) + qJ(5);
t194 = t174 * t90 + t254 * t91;
t18 = qJDD(4) * qJ(5) + (qJD(5) + t50) * qJD(4) + t194;
t145 = pkin(3) * t170 + pkin(2);
t132 = -t145 - t253;
t111 = qJDD(1) * t132 + qJDD(3);
t28 = pkin(4) * t82 - qJ(5) * t81 - qJD(5) * t118 + t111;
t7 = t166 * t28 + t169 * t18;
t121 = t189 * qJD(4);
t233 = t130 * t169;
t234 = t130 * t166;
t35 = t121 * t129 + t226 * t233 - t227 * t234;
t75 = t129 * t130;
t245 = -t35 * t112 - t75 * t80;
t45 = qJD(4) * qJ(5) + t51;
t114 = qJD(1) * t132 + qJD(3);
t59 = pkin(4) * t116 - qJ(5) * t118 + t114;
t22 = t166 * t59 + t169 * t45;
t79 = pkin(4) * t118 + qJ(5) * t116;
t32 = t166 * t79 + t169 * t50;
t123 = t247 * t167;
t124 = t247 * t170;
t190 = -t254 * t123 - t174 * t124;
t52 = qJD(3) * t189 + qJD(4) * t190;
t60 = pkin(4) * t122 - qJ(5) * t121 - qJD(5) * t130;
t25 = t166 * t60 + t169 * t52;
t72 = -pkin(4) * t189 - qJ(5) * t130 + t132;
t78 = -t174 * t123 + t254 * t124;
t37 = t166 * t72 + t169 * t78;
t243 = t118 * t47;
t242 = t118 * t198;
t240 = t166 * t82;
t239 = t169 * t82;
t236 = t116 * t166;
t235 = t121 * t166;
t232 = t156 * t158;
t231 = t158 * t159;
t43 = -qJD(4) * pkin(4) + qJD(5) - t50;
t229 = -qJD(5) + t43;
t228 = t167 ^ 2 + t170 ^ 2;
t6 = -t166 * t18 + t169 * t28;
t2 = pkin(5) * t82 - pkin(8) * t70 + t6;
t5 = -pkin(8) * t69 + t7;
t220 = -t173 * t5 + t176 * t2;
t21 = -t166 * t45 + t169 * t59;
t31 = -t166 * t50 + t169 * t79;
t24 = -t166 * t52 + t169 * t60;
t36 = -t166 * t78 + t169 * t72;
t213 = -t237 * t112 - t127 * t80;
t177 = cos(qJ(1));
t210 = -g(2) * t177 + t250;
t209 = -t166 * t7 - t169 * t6;
t208 = -t6 * t166 + t7 * t169;
t207 = t173 * t2 + t176 * t5;
t34 = -t130 * t120 - t121 * t127;
t76 = t127 * t130;
t206 = -t112 * t34 + t76 * t80;
t10 = pkin(5) * t116 - pkin(8) * t102 + t21;
t14 = -pkin(8) * t100 + t22;
t3 = t10 * t176 - t14 * t173;
t4 = t10 * t173 + t14 * t176;
t204 = t12 * t189 - t122 * t47;
t203 = -t21 * t166 + t22 * t169;
t23 = -pkin(5) * t189 - pkin(8) * t233 + t36;
t30 = -pkin(8) * t234 + t37;
t202 = -t173 * t30 + t176 * t23;
t201 = t173 * t23 + t176 * t30;
t200 = -t116 * t121 - t130 * t82;
t199 = -t100 * t169 + t102 * t166;
t106 = -t125 * t167 + t150;
t197 = -t106 * t167 + t107 * t170;
t196 = (-t135 * t167 + t153) * t167 - t110 * t170;
t195 = pkin(4) * t158 + qJ(5) * t155 + t145;
t133 = t246 * t166;
t188 = pkin(8) * t236 - t169 * qJD(5) + qJD(6) * t133 + t32;
t134 = t246 * t169;
t187 = pkin(5) * t118 + t166 * qJD(5) + qJD(6) * t134 + t116 * t252 + t31;
t186 = -t224 + t261;
t182 = t43 * t121 + t19 * t130 - t212;
t53 = qJD(3) * t130 + qJD(4) * t78;
t172 = -pkin(7) - qJ(3);
t163 = pkin(11) + qJ(6);
t160 = t177 * pkin(1);
t157 = cos(t163);
t154 = sin(t163);
t147 = -pkin(5) * t169 - pkin(4);
t99 = t154 * t156 + t157 * t231;
t98 = -t154 * t231 + t156 * t157;
t97 = t154 * t159 - t157 * t232;
t96 = t154 * t232 + t157 * t159;
t84 = -qJD(4) * t122 + qJDD(4) * t189;
t83 = qJD(4) * t121 + qJDD(4) * t130;
t56 = pkin(5) * t234 - t190;
t39 = pkin(5) * t235 + t53;
t38 = -pkin(5) * t236 + t51;
t33 = t100 * pkin(5) + t43;
t17 = -pkin(8) * t235 + t25;
t13 = pkin(5) * t122 - t121 * t252 + t24;
t8 = t69 * pkin(5) + t19;
t1 = [qJDD(1), t210, g(1) * t177 + g(2) * t175 (t210 + (t168 ^ 2 + t171 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t186 * t170, -t186 * t167, t125 * t228 + t197 - t212, t131 * t146 - g(1) * (-pkin(1) * t175 - pkin(2) * t156 + qJ(3) * t159) - g(2) * (pkin(2) * t159 + qJ(3) * t156 + t160) + t197 * t141 - t196 * qJD(3), t118 * t121 + t130 * t81, -t118 * t122 + t189 * t81 + t200, t83, t84, 0, -qJD(4) * t53 + qJDD(4) * t190 - t111 * t189 + t114 * t122 + t132 * t82 - t193, -qJD(4) * t52 - qJDD(4) * t78 + t111 * t130 + t114 * t121 + t132 * t81 + t258, t53 * t100 + t24 * t116 + t21 * t122 + t166 * t182 - t169 * t193 - t189 * t6 - t190 * t69 + t36 * t82, t53 * t102 - t25 * t116 - t22 * t122 + t166 * t193 + t169 * t182 + t189 * t7 - t190 * t70 - t37 * t82, -t100 * t25 - t102 * t24 - t36 * t70 - t37 * t69 - t258 + t209 * t130 + (-t166 * t22 - t169 * t21) * t121, pkin(1) * t250 - g(2) * t160 - t19 * t190 + t21 * t24 + t22 * t25 + t6 * t36 + t7 * t37 + t43 * t53 + (g(1) * t172 - g(2) * t195) * t159 + (g(1) * t195 + g(2) * t172) * t156, -t11 * t76 - t198 * t34, -t11 * t75 + t12 * t76 + t198 * t35 - t34 * t47, -t206 + t255, t204 + t245, t112 * t122 - t189 * t80 (t13 * t176 - t17 * t173) * t112 + t202 * t80 - t220 * t189 + t3 * t122 + t39 * t47 + t56 * t12 + t8 * t75 + t33 * t35 - g(1) * t97 - g(2) * t99 + (-t112 * t201 + t189 * t4) * qJD(6) -(t13 * t173 + t17 * t176) * t112 - t201 * t80 + t207 * t189 - t4 * t122 - t39 * t198 + t56 * t11 - t8 * t76 + t33 * t34 - g(1) * t96 - g(2) * t98 + (-t112 * t202 + t189 * t3) * qJD(6); 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, t106 * t170 + t107 * t167 - g(3), 0, 0, 0, 0, 0, t84, -t83, t100 * t122 + t166 * t200 - t189 * t69, t102 * t122 + t169 * t200 - t189 * t70 (t166 * t70 - t169 * t69) * t130 + t199 * t121, t121 * t203 + t122 * t43 + t130 * t208 - t189 * t19 - g(3), 0, 0, 0, 0, 0, -t204 + t245, t206 + t255; 0, 0, 0, 0, -t222, t223, -t228 * qJD(1) ^ 2, qJD(1) * t196 - t261, 0, 0, 0, 0, 0, 0.2e1 * t118 * qJD(4) + t205 (-t116 - t217) * qJD(4) + t221, -t100 * t118 - t115 * t166 + t239, -t102 * t118 - t115 * t169 - t240, t116 * t199 - t166 * t69 - t169 * t70, t116 * t203 - t118 * t43 - t209 + t216, 0, 0, 0, 0, 0, t213 - t243, t242 + t256; 0, 0, 0, 0, 0, 0, 0, 0, t118 * t116, t118 ^ 2 - t115 (t116 - t217) * qJD(4) + t221, -t205, qJDD(4), -t114 * t118 + t185 - t192 + t262, t114 * t116 - t183 - t194, -qJ(5) * t240 - pkin(4) * t69 - t100 * t51 - t118 * t21 + (t229 * t166 - t31) * t116 - t180 * t169, -qJ(5) * t239 - pkin(4) * t70 - t102 * t51 + t118 * t22 + (t229 * t169 + t32) * t116 + t180 * t166, t100 * t32 + t102 * t31 + (-qJ(5) * t69 - qJD(5) * t100 - t116 * t21 + t7) * t169 + (qJ(5) * t70 + qJD(5) * t102 - t116 * t22 - t6) * t166 + t183, -t21 * t31 - t22 * t32 - t43 * t51 + t203 * qJD(5) - t180 * pkin(4) + (t183 + t208) * qJ(5), t11 * t129 + t198 * t238, -t11 * t127 - t12 * t129 + t198 * t237 + t238 * t47, t242 - t256, t213 + t243, -t112 * t118 (-t133 * t176 - t134 * t173) * t80 + t147 * t12 + t8 * t127 - t3 * t118 - t38 * t47 + t237 * t33 + (t173 * t188 - t176 * t187) * t112 + t185 * t157 -(-t133 * t173 + t134 * t176) * t80 + t147 * t11 + t8 * t129 + t4 * t118 + t38 * t198 - t238 * t33 + (t173 * t187 + t176 * t188) * t112 - t185 * t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102 * t116 + t69, -t100 * t116 + t70, -t100 ^ 2 - t102 ^ 2, t100 * t22 + t102 * t21 + t180, 0, 0, 0, 0, 0, t12 - t259, t11 - t260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198 * t47, t198 ^ 2 - t47 ^ 2, t11 + t260, -t12 - t259, t80, -g(1) * t98 + g(2) * t96 + t154 * t249 + t198 * t33 + t257 * t4 + t220, g(1) * t99 - g(2) * t97 + t157 * t249 + t257 * t3 + t33 * t47 - t207;];
tau_reg  = t1;
