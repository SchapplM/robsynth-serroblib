% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:56:35
% EndTime: 2021-01-15 12:56:47
% DurationCPUTime: 2.97s
% Computational Cost: add. (2804->327), mult. (6757->425), div. (0->0), fcn. (4829->10), ass. (0->200)
t139 = sin(pkin(8));
t141 = sin(qJ(4));
t145 = cos(qJ(3));
t199 = qJDD(1) * t145;
t142 = sin(qJ(3));
t200 = qJDD(1) * t142;
t144 = cos(qJ(4));
t261 = t139 * t144;
t197 = qJD(3) + qJD(4);
t87 = t141 * t145 + t142 * t144;
t265 = t197 * t87;
t22 = t139 * (qJD(1) * t265 + t141 * t200) - t199 * t261;
t134 = t139 ^ 2;
t256 = 0.2e1 * t134;
t231 = qJDD(1) * pkin(1);
t122 = qJDD(2) - t231;
t143 = sin(qJ(1));
t146 = cos(qJ(1));
t259 = -g(2) * t146 - g(3) * t143;
t264 = t259 - t122;
t140 = cos(pkin(8));
t198 = t140 * qJDD(1);
t113 = -qJDD(3) + t198;
t100 = -qJDD(4) + t113;
t206 = qJD(4) * t144;
t211 = qJD(1) * t140;
t115 = -qJD(3) + t211;
t106 = -qJD(4) + t115;
t252 = pkin(3) * t106;
t263 = t141 * pkin(3) * t100 + t206 * t252;
t227 = t139 * t145;
t196 = pkin(7) * t227;
t234 = qJ(2) * t142;
t160 = -t140 * t234 - t196;
t93 = -pkin(2) * t140 - pkin(6) * t139 - pkin(1);
t73 = qJD(1) * t93 + qJD(2);
t65 = t145 * t73;
t42 = qJD(1) * t160 + t65;
t28 = -pkin(3) * t115 + t42;
t233 = qJ(2) * t145;
t114 = t140 * t233;
t212 = qJD(1) * t139;
t193 = t142 * t212;
t43 = -pkin(7) * t193 + qJD(1) * t114 + t142 * t73;
t35 = t141 * t43;
t184 = t144 * t28 - t35;
t172 = t141 * t193;
t210 = qJD(1) * t145;
t192 = t139 * t210;
t60 = t144 * t192 - t172;
t54 = t60 * qJ(5);
t8 = t184 - t54;
t241 = qJ(5) * t22;
t94 = t100 * pkin(4);
t262 = t241 - t94;
t203 = qJD(1) * qJD(2);
t204 = qJ(2) * qJDD(1);
t163 = t203 + t204;
t260 = t140 * t163;
t138 = qJ(3) + qJ(4);
t125 = sin(t138);
t250 = g(1) * t139;
t126 = cos(t138);
t226 = t140 * t143;
t66 = -t125 * t226 - t126 * t146;
t218 = t146 * t125;
t221 = t143 * t126;
t68 = t140 * t218 - t221;
t258 = -g(2) * t66 - g(3) * t68 + t125 * t250;
t257 = t60 ^ 2;
t135 = t140 ^ 2;
t5 = -pkin(4) * t106 + t8;
t255 = -t8 + t5;
t251 = pkin(7) * t139;
t249 = g(2) * t143;
t248 = g(3) * t146;
t161 = qJD(1) * t87;
t57 = t139 * t161;
t247 = t60 * t57;
t246 = t144 * t42 - t35;
t83 = t145 * t93;
t46 = -t196 + t83 + (-pkin(3) - t234) * t140;
t228 = t139 * t142;
t235 = t142 * t93 + t114;
t53 = -pkin(7) * t228 + t235;
t245 = t141 * t46 + t144 * t53;
t224 = t141 * t142;
t86 = t144 * t145 - t224;
t244 = (t197 - t211) * t86;
t243 = t140 * t161 - t265;
t242 = t197 * t172;
t176 = t197 * t145;
t187 = t141 * t199;
t23 = (t187 + (qJD(1) * t176 + t200) * t144) * t139 - t242;
t240 = qJ(5) * t23;
t239 = qJ(5) * t57;
t238 = t106 * t57;
t237 = t106 * t60;
t37 = t144 * t43;
t208 = qJD(3) * t145;
t209 = qJD(2) * t140;
t236 = t145 * t209 + t93 * t208;
t232 = qJD(3) * t73;
t230 = (-qJ(5) - pkin(7) - pkin(6)) * t139;
t147 = qJD(1) ^ 2;
t229 = t134 * t147;
t225 = t140 * t146;
t223 = t142 * t143;
t222 = t142 * t146;
t220 = t143 * t145;
t219 = t145 * t146;
t186 = -pkin(4) * t57 - qJD(5);
t74 = pkin(3) * t193 + qJ(2) * t212;
t40 = -t186 + t74;
t217 = qJD(5) + t40;
t92 = t145 * pkin(3) + pkin(4) * t126;
t109 = t139 * pkin(3) * t208;
t80 = t139 * qJD(2) + t109;
t116 = pkin(3) * t228;
t85 = t139 * qJ(2) + t116;
t216 = t146 * pkin(1) + t143 * qJ(2);
t214 = t134 + t135;
t137 = t145 ^ 2;
t213 = t142 ^ 2 - t137;
t207 = qJD(4) * t141;
t205 = qJD(3) + t115;
t202 = qJD(1) * qJD(3);
t201 = qJDD(1) * t139;
t194 = qJ(2) * qJD(3) * t140;
t191 = qJ(2) * t198;
t190 = t145 * t203;
t189 = t145 * t202;
t171 = qJD(1) * t194;
t72 = qJDD(1) * t93 + qJDD(2);
t64 = t145 * t72;
t15 = -pkin(3) * t113 + t64 + (-pkin(7) * t201 - t171) * t145 + (-t191 - t232 + (qJD(3) * t251 - t209) * qJD(1)) * t142;
t178 = t140 * t190 + t142 * t72 + t145 * t191 + t73 * t208;
t156 = -t142 * t171 + t178;
t21 = (-t189 - t200) * t251 + t156;
t185 = -t141 * t21 + t144 * t15;
t38 = qJD(3) * t160 + t236;
t39 = -t142 * t209 + (-t114 + (-t93 + t251) * t142) * qJD(3);
t183 = -t141 * t38 + t144 * t39;
t182 = -t141 * t42 - t37;
t181 = -t141 * t53 + t144 * t46;
t180 = -t141 * t15 - t144 * t21 - t28 * t206 + t43 * t207;
t179 = t214 * t147;
t177 = qJD(1) * t205;
t175 = pkin(3) * t192;
t49 = qJ(2) * t201 + qJD(1) * t109 + qJDD(1) * t116 + t139 * t203;
t174 = t113 + t198;
t173 = 0.2e1 * t214;
t170 = g(2) * t68 - g(3) * t66;
t67 = t140 * t221 - t218;
t69 = t125 * t143 + t126 * t225;
t169 = -g(2) * t69 - g(3) * t67;
t167 = -t248 + t249;
t166 = -t141 * t28 - t37;
t164 = qJD(3) * (t115 + t211);
t162 = t141 * t39 + t144 * t38 + t46 * t206 - t207 * t53;
t12 = pkin(4) * t23 + qJDD(5) + t49;
t159 = -t115 ^ 2 - t229;
t158 = g(2) * t67 - g(3) * t69 + t126 * t250 + t180;
t157 = t173 * t203 + t248;
t155 = qJD(4) * t166 + t185;
t154 = t74 * t57 + t158;
t152 = t217 * t57 + t158 + t240;
t151 = t155 + t258;
t150 = -t74 * t60 + t151;
t128 = t143 * pkin(1);
t121 = pkin(3) * t144 + pkin(4);
t91 = pkin(3) * t142 + pkin(4) * t125;
t88 = pkin(2) + t92;
t78 = t140 * t219 + t223;
t77 = t140 * t222 - t220;
t76 = t140 * t220 - t222;
t75 = -t140 * t223 - t219;
t71 = t86 * t139;
t70 = t87 * t139;
t55 = t57 ^ 2;
t50 = pkin(4) * t60 + t175;
t47 = pkin(4) * t70 + t85;
t34 = -t197 * t139 * t224 + t176 * t261;
t33 = t265 * t139;
t25 = pkin(4) * t34 + t80;
t24 = -t55 + t257;
t20 = -qJ(5) * t70 + t245;
t18 = -pkin(4) * t140 - qJ(5) * t71 + t181;
t17 = -t237 + (-t187 + (-t197 * t210 - t200) * t144) * t139 + t242;
t16 = -t22 - t238;
t11 = -t54 + t246;
t10 = t182 + t239;
t9 = -t166 - t239;
t7 = -t100 * t86 - t243 * t106 - t57 * t212;
t6 = t100 * t87 + t244 * t106 - t60 * t212;
t4 = qJ(5) * t33 - qJD(4) * t245 - qJD(5) * t71 + t183;
t3 = -qJ(5) * t34 - qJD(5) * t70 + t162;
t2 = -qJD(5) * t57 - t180 - t240;
t1 = -qJD(5) * t60 + t155 + t262;
t13 = [qJDD(1), t259, t167, (t231 + t264) * t140, t173 * t204 + t157 - t249, -t122 * pkin(1) - g(2) * t216 - g(3) * t128 + (t214 * t204 + t157) * qJ(2), (qJDD(1) * t137 - 0.2e1 * t142 * t189) * t134, (-t142 * t199 + t202 * t213) * t256, (t142 * t164 - t145 * t174) * t139, (t142 * t174 + t145 * t164) * t139, t113 * t140, -g(2) * t78 - g(3) * t76 - t83 * t113 - t64 * t140 + (t115 * t140 + (t256 + t135) * qJD(1)) * qJ(2) * t208 + (qJD(3) * t93 * t115 + t163 * t256 + (qJ(2) * t113 + qJD(2) * t115 + t232 + t260) * t140) * t142, (-t142 * t194 + t236) * t115 + t235 * t113 + t156 * t140 + g(2) * t77 - g(3) * t75 + (t190 + (-t142 * t202 + t199) * qJ(2)) * t256, -t22 * t71 - t33 * t60, t22 * t70 - t23 * t71 + t33 * t57 - t34 * t60, -t100 * t71 + t106 * t33 + t140 * t22, t100 * t70 + t106 * t34 + t140 * t23, t100 * t140, -t183 * t106 - t181 * t100 - t185 * t140 + t80 * t57 + t85 * t23 + t49 * t70 + t74 * t34 + (t106 * t245 - t140 * t166) * qJD(4) + t169, t100 * t245 + t106 * t162 - t140 * t180 - t85 * t22 - t74 * t33 + t49 * t71 + t80 * t60 + t170, -t1 * t140 - t100 * t18 - t106 * t4 + t12 * t70 + t23 * t47 + t25 * t57 + t34 * t40 + t169, t100 * t20 + t106 * t3 + t12 * t71 + t140 * t2 - t22 * t47 + t25 * t60 - t33 * t40 + t170, -t1 * t71 + t139 * t259 + t18 * t22 - t2 * t70 - t20 * t23 - t3 * t57 + t33 * t5 - t34 * t9 - t4 * t60, t2 * t20 + t9 * t3 + t1 * t18 + t5 * t4 + t12 * t47 + t40 * t25 - g(2) * (t143 * t91 + t216) - g(3) * (-t143 * t230 + t226 * t88 + t128) + (-g(2) * (t140 * t88 - t230) - g(3) * (-qJ(2) - t91)) * t146; 0, 0, 0, -t198, -t179, -qJ(2) * t179 - t264, 0, 0, 0, 0, 0, -t113 * t145 + t142 * t159, t113 * t142 + t145 * t159, 0, 0, 0, 0, 0, t7, t6, t7, t6, t22 * t86 - t23 * t87 - t243 * t60 - t244 * t57, t1 * t86 + t2 * t87 - t212 * t40 + t243 * t5 + t244 * t9 - t259; 0, 0, 0, 0, 0, 0, t145 * t142 * t229, -t213 * t229, (-t142 * t177 + t199) * t139, (-t145 * t177 - t200) * t139, -t113, -g(2) * t75 - g(3) * t77 + t64 + (-t140 * t177 - t229) * t233 + (-t205 * t73 + t250 - t260) * t142, g(1) * t227 + g(2) * t76 - g(3) * t78 - t65 * t115 + (t205 * t211 + t229) * t234 - t178, t247, t24, t16, t17, -t100, t182 * t106 + (-t100 * t144 + t106 * t207 - t192 * t57) * pkin(3) + t150, -t106 * t246 - t175 * t60 + t154 + t263, t10 * t106 - t100 * t121 - t50 * t57 - t217 * t60 + (-t37 + (-t28 + t252) * t141) * qJD(4) + t185 + t258 + t262, -t106 * t11 - t50 * t60 + t152 + t263, t121 * t22 + (t10 + t9) * t60 + (t11 - t5) * t57 + (-t141 * t23 + (t141 * t60 - t144 * t57) * qJD(4)) * pkin(3), t1 * t121 - t9 * t11 - t5 * t10 - t40 * t50 + t91 * t250 - g(2) * (-t146 * t92 - t226 * t91) - g(3) * (-t143 * t92 + t225 * t91) + (t2 * t141 + (-t141 * t5 + t144 * t9) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t247, t24, t16, t17, -t100, t106 * t166 + t150, -t106 * t184 + t154, t241 - t106 * t9 - 0.2e1 * t94 + (t186 - t40) * t60 + t151, -pkin(4) * t257 - t106 * t8 + t152, pkin(4) * t22 - t255 * t57, t255 * t9 + (-t40 * t60 + t1 + t258) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 - t237, -t22 + t238, -t55 - t257, g(1) * t140 - t139 * t167 + t5 * t60 + t57 * t9 + t12;];
tau_reg = t13;
