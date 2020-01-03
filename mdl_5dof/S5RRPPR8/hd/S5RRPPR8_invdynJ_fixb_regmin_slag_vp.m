% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:17
% EndTime: 2019-12-31 19:39:24
% DurationCPUTime: 2.14s
% Computational Cost: add. (1672->330), mult. (3686->430), div. (0->0), fcn. (2526->10), ass. (0->185)
t153 = qJD(2) - qJD(5);
t163 = sin(qJ(5));
t166 = cos(qJ(5));
t160 = sin(pkin(8));
t161 = cos(pkin(8));
t164 = sin(qJ(2));
t167 = cos(qJ(2));
t93 = t164 * t160 + t167 * t161;
t81 = t93 * qJD(1);
t224 = qJD(1) * t167;
t211 = t160 * t224;
t225 = qJD(1) * t164;
t83 = t161 * t225 - t211;
t248 = t163 * t83 + t166 * t81;
t238 = t248 * t153;
t220 = qJD(5) * t166;
t221 = qJD(5) * t163;
t218 = qJD(1) * qJD(2);
t210 = t164 * t218;
t176 = qJDD(1) * t93 - t161 * t210;
t209 = t167 * t218;
t41 = t160 * t209 + t176;
t139 = t164 * qJDD(1);
t217 = t167 * qJDD(1);
t196 = t161 * t139 - t160 * t217;
t86 = t93 * qJD(2);
t42 = qJD(1) * t86 + t196;
t3 = -t163 * t41 + t166 * t42 - t81 * t220 - t83 * t221;
t257 = t3 - t238;
t34 = -t163 * t81 + t166 * t83;
t237 = t34 * t153;
t4 = t163 * t42 + t166 * t41 + t83 * t220 - t81 * t221;
t256 = -t4 - t237;
t255 = t34 * t248;
t254 = t209 + t139;
t253 = t248 ^ 2 - t34 ^ 2;
t129 = pkin(6) * t139;
t208 = pkin(6) * t209 + qJDD(3) + t129;
t219 = t164 * qJD(4);
t245 = pkin(2) + pkin(3);
t26 = -t254 * qJ(4) - qJD(1) * t219 - t245 * qJDD(2) + t208;
t130 = pkin(6) * t217;
t154 = qJDD(2) * qJ(3);
t155 = qJD(2) * qJD(3);
t215 = t130 + t154 + t155;
t223 = qJD(2) * t164;
t241 = pkin(6) - qJ(4);
t77 = -t167 * qJD(4) - t241 * t223;
t28 = -qJ(4) * t217 + qJD(1) * t77 + t215;
t206 = t160 * t28 - t161 * t26;
t1 = -qJDD(2) * pkin(4) - t42 * pkin(7) - t206;
t152 = pkin(8) + qJ(5);
t137 = sin(t152);
t138 = cos(t152);
t185 = t164 * t137 + t167 * t138;
t6 = t160 * t26 + t161 * t28;
t2 = -t41 * pkin(7) + t6;
t91 = -qJD(1) * pkin(1) - pkin(2) * t224 - qJ(3) * t225;
t57 = pkin(3) * t224 + qJD(4) - t91;
t29 = t81 * pkin(4) + t57;
t165 = sin(qJ(1));
t186 = t167 * t137 - t164 * t138;
t53 = t186 * t165;
t168 = cos(qJ(1));
t230 = t167 * t168;
t233 = t164 * t168;
t55 = t137 * t230 - t138 * t233;
t252 = -g(1) * t55 - g(2) * t53 - g(3) * t185 - t166 * t1 + t163 * t2 + t29 * t34;
t54 = t185 * t165;
t56 = t185 * t168;
t251 = -g(1) * t56 - g(2) * t54 + g(3) * t186 + t163 * t1 + t166 * t2 - t29 * t248;
t228 = t167 * pkin(2) + t164 * qJ(3);
t250 = -pkin(1) - t228;
t249 = g(1) * t168 + g(2) * t165;
t247 = qJD(5) + t153;
t236 = pkin(6) * qJDD(2);
t246 = (qJD(1) * t250 + t91) * qJD(2) - t236;
t244 = t81 * pkin(7);
t243 = t83 * pkin(7);
t149 = g(1) * t165;
t242 = g(2) * t168;
t143 = t167 * pkin(3);
t109 = t241 * t167;
t78 = qJD(2) * t109 - t219;
t23 = t160 * t78 + t161 * t77;
t212 = t245 * qJD(2);
t133 = pkin(6) * t225;
t99 = qJ(4) * t225 - t133;
t65 = qJD(3) - t212 - t99;
t134 = pkin(6) * t224;
t101 = -qJ(4) * t224 + t134;
t156 = qJD(2) * qJ(3);
t90 = t101 + t156;
t21 = t160 * t65 + t161 * t90;
t46 = t160 * t101 + t161 * t99;
t108 = t241 * t164;
t49 = t160 * t108 + t161 * t109;
t159 = qJDD(1) * pkin(1);
t235 = qJDD(2) * pkin(2);
t234 = t164 * t165;
t171 = qJD(1) ^ 2;
t232 = t164 * t171;
t231 = t165 * t167;
t140 = t164 * qJD(3);
t222 = qJD(2) * t167;
t229 = qJ(3) * t222 + t140;
t157 = t164 ^ 2;
t158 = t167 ^ 2;
t226 = t157 - t158;
t216 = t167 * t232;
t214 = t143 + t228;
t213 = -g(1) * t233 - g(2) * t234 + g(3) * t167;
t207 = t149 - t242;
t20 = -t160 * t90 + t161 * t65;
t22 = -t160 * t77 + t161 * t78;
t45 = t161 * t101 - t160 * t99;
t48 = t161 * t108 - t160 * t109;
t205 = -qJD(2) * pkin(2) + qJD(3);
t102 = -t160 * qJ(3) - t161 * t245;
t203 = t160 * qJD(3) + t45;
t202 = t161 * qJD(3) - t46;
t92 = pkin(1) + t214;
t201 = t168 * pkin(1) + pkin(2) * t230 + t165 * pkin(6) + qJ(3) * t233;
t200 = -t129 - t213;
t199 = t164 * t212;
t170 = qJD(2) ^ 2;
t198 = pkin(6) * t170 + t242;
t10 = -qJD(2) * pkin(4) + t20 - t243;
t11 = t21 - t244;
t195 = t166 * t10 - t163 * t11;
t194 = -t163 * t10 - t166 * t11;
t193 = t20 * t160 - t21 * t161;
t182 = t167 * t160 - t164 * t161;
t24 = pkin(7) * t182 + t48;
t25 = -t93 * pkin(7) + t49;
t192 = -t163 * t25 + t166 * t24;
t191 = t163 * t24 + t166 * t25;
t43 = -t163 * t182 + t166 * t93;
t44 = -t163 * t93 - t166 * t182;
t190 = pkin(2) * t217 + t254 * qJ(3) + qJD(1) * t140 + t159;
t103 = t161 * qJ(3) - t160 * t245;
t98 = -pkin(4) + t102;
t189 = t166 * t103 + t163 * t98;
t188 = -t163 * t103 + t166 * t98;
t105 = t133 + t205;
t107 = t134 + t156;
t187 = t105 * t167 - t107 * t164;
t184 = t166 * t160 + t163 * t161;
t183 = t163 * t160 - t166 * t161;
t127 = qJ(3) * t224;
t74 = -t245 * t225 + t127;
t73 = t208 - t235;
t180 = -0.2e1 * pkin(1) * t218 - t236;
t179 = t184 * t153;
t178 = t183 * t153;
t52 = -t199 + t229;
t177 = -t198 + 0.2e1 * t159;
t39 = pkin(2) * t210 - t190;
t79 = pkin(2) * t223 - t229;
t175 = -qJD(1) * t79 - qJDD(1) * t250 - t198 - t39;
t64 = -pkin(6) * t210 + t215;
t172 = qJD(2) * t187 + t73 * t164 + t64 * t167;
t14 = pkin(3) * t217 - qJD(1) * t199 + qJDD(4) + t190;
t151 = qJDD(2) - qJDD(5);
t145 = t168 * pkin(6);
t121 = g(1) * t231;
t117 = qJ(3) * t230;
t115 = qJ(3) * t231;
t100 = pkin(2) * t225 - t127;
t85 = t160 * t222 - t161 * t223;
t72 = t93 * t168;
t71 = t182 * t168;
t70 = t93 * t165;
t69 = t182 * t165;
t47 = t93 * pkin(4) + t92;
t40 = -t83 * pkin(4) + t74;
t27 = t85 * pkin(4) + t52;
t16 = t46 + t243;
t15 = t45 - t244;
t13 = -t85 * pkin(7) + t23;
t12 = -t86 * pkin(7) + t22;
t9 = t44 * qJD(5) + t163 * t86 + t166 * t85;
t8 = -t43 * qJD(5) - t163 * t85 + t166 * t86;
t7 = t41 * pkin(4) + t14;
t5 = [qJDD(1), t207, t249, t157 * qJDD(1) + 0.2e1 * t164 * t209, 0.2e1 * t164 * t217 - 0.2e1 * t226 * t218, qJDD(2) * t164 + t170 * t167, qJDD(2) * t167 - t170 * t164, 0, t164 * t180 + t167 * t177 + t121, t180 * t167 + (-t177 - t149) * t164, t246 * t164 + t175 * t167 + t121, (t157 + t158) * qJDD(1) * pkin(6) + t172 - t249, -t246 * t167 + (t175 + t149) * t164, pkin(6) * t172 - g(1) * t145 - g(2) * t201 + t91 * t79 + (-t149 + t39) * t250, g(1) * t70 - g(2) * t72 - t22 * qJD(2) - t48 * qJDD(2) + t14 * t93 + t92 * t41 + t52 * t81 + t57 * t85, -g(1) * t69 + g(2) * t71 + t23 * qJD(2) + t49 * qJDD(2) - t14 * t182 + t92 * t42 + t52 * t83 + t57 * t86, -t182 * t206 - t20 * t86 - t21 * t85 - t22 * t83 - t23 * t81 - t49 * t41 - t48 * t42 - t6 * t93 + t249, t6 * t49 + t21 * t23 - t206 * t48 + t20 * t22 + t14 * t92 + t57 * t52 - g(1) * (-t168 * qJ(4) + t145) - g(2) * (pkin(3) * t230 + t201) + (-g(1) * (t250 - t143) + g(2) * qJ(4)) * t165, t3 * t44 + t34 * t8, -t248 * t8 - t3 * t43 - t34 * t9 - t44 * t4, -t44 * t151 - t8 * t153, t43 * t151 + t9 * t153, 0, t27 * t248 + t47 * t4 + t7 * t43 + t29 * t9 - (-qJD(5) * t191 + t166 * t12 - t163 * t13) * t153 - t192 * t151 + g(1) * t54 - g(2) * t56, t27 * t34 + t47 * t3 + t7 * t44 + t29 * t8 + (qJD(5) * t192 + t163 * t12 + t166 * t13) * t153 + t191 * t151 - g(1) * t53 + g(2) * t55; 0, 0, 0, -t216, t226 * t171, t139, t217, qJDD(2), pkin(1) * t232 + t200, g(3) * t164 - t130 + (pkin(1) * t171 + t249) * t167, 0.2e1 * t235 - qJDD(3) + (t100 * t167 - t164 * t91) * qJD(1) + t200, (-pkin(2) * t164 + qJ(3) * t167) * qJDD(1) + ((t107 - t156) * t164 + (-t105 + t205) * t167) * qJD(1), t130 + 0.2e1 * t154 + 0.2e1 * t155 + (qJD(1) * t100 - g(3)) * t164 + (qJD(1) * t91 - t249) * t167, t64 * qJ(3) + t107 * qJD(3) - t73 * pkin(2) - t91 * t100 - g(1) * (-pkin(2) * t233 + t117) - g(2) * (-pkin(2) * t234 + t115) - g(3) * t228 - t187 * qJD(1) * pkin(6), -g(1) * t71 - g(2) * t69 - g(3) * t93 + qJD(2) * t203 - t102 * qJDD(2) + t57 * t83 - t74 * t81 + t206, -g(1) * t72 - g(2) * t70 + g(3) * t182 + qJD(2) * t202 + t103 * qJDD(2) - t57 * t81 - t74 * t83 + t6, -t102 * t42 - t103 * t41 + (t203 - t21) * t83 + (t20 - t202) * t81, t249 * t164 * t245 - g(1) * t117 - g(2) * t115 - g(3) * t214 - t193 * qJD(3) - t102 * t206 + t6 * t103 - t20 * t45 - t21 * t46 - t57 * t74, -t255, t253, -t257, -t256, t151, -t188 * t151 - t40 * t248 + (t166 * t15 - t163 * t16) * t153 + qJD(3) * t179 + (t153 * t189 - t194) * qJD(5) + t252, t189 * t151 - t40 * t34 - (t163 * t15 + t166 * t16) * t153 - qJD(3) * t178 + (t188 * t153 + t195) * qJD(5) + t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) - t216, t139, -t157 * t171 - t170, -t107 * qJD(2) + t225 * t91 + t213 + t73, -t161 * qJDD(2) - t160 * t170 - t225 * t81, t160 * qJDD(2) - t161 * t170 - t225 * t83, -t160 * t41 - t161 * t42 + (-t160 * t83 + t161 * t81) * qJD(2), qJD(2) * t193 + t6 * t160 - t161 * t206 - t225 * t57 + t213, 0, 0, 0, 0, 0, t151 * t183 - t153 * t179 - t225 * t248, t151 * t184 + t153 * t178 - t225 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t83 + t211) * qJD(2) + t176, 0.2e1 * t81 * qJD(2) + t196, -t81 ^ 2 - t83 ^ 2, t20 * t83 + t21 * t81 + t14 + t207, 0, 0, 0, 0, 0, t4 - t237, t3 + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t255, -t253, t257, t256, -t151, t247 * t194 - t252, -t247 * t195 - t251;];
tau_reg = t5;
