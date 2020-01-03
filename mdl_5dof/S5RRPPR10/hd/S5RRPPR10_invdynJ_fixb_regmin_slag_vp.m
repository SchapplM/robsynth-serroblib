% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPR10
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:43
% EndTime: 2019-12-31 19:44:52
% DurationCPUTime: 3.21s
% Computational Cost: add. (2056->409), mult. (4779->532), div. (0->0), fcn. (3272->8), ass. (0->209)
t150 = cos(qJ(2));
t141 = g(3) * t150;
t147 = sin(qJ(2));
t217 = t147 * qJDD(1);
t127 = pkin(6) * t217;
t218 = qJD(1) * qJD(2);
t203 = t150 * t218;
t174 = qJDD(2) * pkin(2) - pkin(6) * t203 - qJDD(3) - t127;
t206 = t174 - t141;
t148 = sin(qJ(1));
t151 = cos(qJ(1));
t191 = g(1) * t151 + g(2) * t148;
t271 = t191 * t147;
t154 = -t271 - t206;
t144 = sin(pkin(8));
t149 = cos(qJ(5));
t145 = cos(pkin(8));
t146 = sin(qJ(5));
t243 = t146 * t145;
t95 = t149 * t144 - t243;
t69 = t95 * t147;
t232 = t150 * pkin(2) + t147 * qJ(3);
t273 = -pkin(1) - t232;
t229 = qJD(1) * t147;
t209 = t145 * t229;
t224 = t144 * qJD(2);
t91 = t209 + t224;
t87 = t91 ^ 2;
t210 = t144 * t229;
t222 = t145 * qJD(2);
t89 = t210 - t222;
t272 = -t89 ^ 2 - t87;
t228 = qJD(2) * t147;
t270 = qJ(4) * t228 - t150 * qJD(4);
t134 = t150 * qJDD(1);
t202 = t147 * t218;
t164 = -t202 + t134;
t220 = t150 * qJD(1);
t120 = qJD(5) + t220;
t269 = -qJD(5) + t120;
t37 = t146 * t89 + t149 * t91;
t132 = t145 * qJDD(2);
t165 = t203 + t217;
t59 = t144 * t165 - t132;
t60 = t144 * qJDD(2) + t145 * t165;
t8 = qJD(5) * t37 + t146 * t60 - t149 * t59;
t233 = t91 * qJD(4);
t249 = t60 * qJ(4);
t156 = t174 + t233 + t249;
t265 = pkin(3) + pkin(4);
t6 = -t265 * t59 + t156;
t267 = t6 + t271;
t266 = -0.2e1 * pkin(1);
t264 = pkin(6) * t91;
t263 = t59 * pkin(3);
t262 = g(1) * t148;
t259 = g(2) * t151;
t258 = g(3) * t147;
t137 = t150 * pkin(3);
t257 = -pkin(7) + qJ(3);
t187 = pkin(2) * t147 - qJ(3) * t150;
t75 = qJD(2) * t187 - t147 * qJD(3);
t33 = qJD(1) * t75 + qJDD(1) * t273;
t66 = pkin(6) * t164 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t13 = t144 * t33 + t145 * t66;
t94 = t146 * t144 + t149 * t145;
t166 = t94 * t150;
t256 = -qJD(1) * t166 - t94 * qJD(5);
t211 = t144 * t220;
t215 = t150 * t243;
t225 = qJD(5) * t149;
t226 = qJD(5) * t146;
t255 = -qJD(1) * t215 + t144 * t225 - t145 * t226 + t149 * t211;
t130 = pkin(6) * t220;
t108 = qJD(2) * qJ(3) + t130;
t85 = t273 * qJD(1);
t40 = t145 * t108 + t144 * t85;
t254 = t13 * t145;
t253 = t145 * t75;
t98 = t187 * qJD(1);
t252 = t145 * t98;
t35 = t146 * t91 - t149 * t89;
t251 = t35 * t120;
t250 = t37 * t120;
t244 = t145 * t150;
t62 = pkin(6) * t244 + t144 * t273;
t248 = qJ(4) * t145;
t247 = t144 * t147;
t246 = t144 * t150;
t245 = t145 * t147;
t242 = t147 * t148;
t241 = t147 * t151;
t240 = t148 * t145;
t239 = t148 * t150;
t237 = t150 * qJ(4);
t236 = t150 * t151;
t235 = t151 * t144;
t234 = t151 * t145;
t142 = t147 ^ 2;
t143 = t150 ^ 2;
t231 = t142 - t143;
t230 = qJD(1) * t145;
t227 = qJD(2) * t150;
t223 = t144 * qJD(4);
t221 = t145 * qJD(3);
t216 = pkin(6) * t228;
t214 = -pkin(6) * t144 - pkin(3);
t12 = -t144 * t66 + t145 * t33;
t175 = pkin(3) * t134 + qJDD(4) - t12;
t11 = -pkin(3) * t202 + t175;
t4 = pkin(4) * t164 - t60 * pkin(7) + t11;
t9 = qJ(4) * t202 + (-qJ(4) * qJDD(1) - qJD(1) * qJD(4)) * t150 + t13;
t5 = t59 * pkin(7) + t9;
t213 = -t146 * t5 + t149 * t4;
t212 = qJ(3) * t228;
t208 = qJD(4) * t245;
t10 = -t156 + t263;
t207 = -t10 - t141;
t205 = g(3) * t232;
t204 = qJ(3) * t134;
t201 = t144 * qJ(4) + pkin(2);
t39 = -t144 * t108 + t145 * t85;
t172 = -t265 * t144 + t248;
t199 = -t172 * t220 + t130 + t223;
t114 = pkin(6) * t246;
t61 = t145 * t273 - t114;
t198 = qJD(3) * t211 + t144 * t204 + (g(1) * t234 + g(2) * t240) * t147;
t197 = t151 * pkin(1) + pkin(2) * t236 + t148 * pkin(6) + qJ(3) * t241;
t196 = -t145 * qJ(3) * t59 - t89 * t221 - t258;
t195 = t145 * t204;
t79 = t144 * t239 + t234;
t81 = t150 * t235 - t240;
t194 = -g(1) * t79 + g(2) * t81;
t80 = t145 * t239 - t235;
t82 = t148 * t144 + t150 * t234;
t193 = g(1) * t80 - g(2) * t82;
t49 = -t237 + t62;
t192 = t214 * t147;
t190 = t146 * t4 + t149 * t5;
t189 = qJD(2) * pkin(2) - pkin(6) * t229 - qJD(3);
t188 = pkin(6) * t89 - t144 * t189;
t186 = pkin(3) * t144 - t248;
t26 = pkin(3) * t220 + qJD(4) - t39;
t14 = pkin(4) * t220 - t91 * pkin(7) + t26;
t31 = -qJ(4) * t220 + t40;
t16 = t89 * pkin(7) + t31;
t1 = t149 * t14 - t146 * t16;
t2 = t146 * t14 + t149 * t16;
t32 = t150 * pkin(4) + t114 + t137 + (-pkin(7) * t147 - t273) * t145;
t38 = pkin(7) * t247 + t49;
t184 = t146 * t38 - t149 * t32;
t183 = t146 * t32 + t149 * t38;
t182 = t80 * t146 - t79 * t149;
t181 = t79 * t146 + t80 * t149;
t180 = qJ(3) * t60 + qJD(3) * t91;
t83 = t144 * t98;
t55 = -pkin(6) * t209 + t83;
t65 = t144 * t75;
t46 = -t145 * t216 + t65;
t179 = t120 ^ 2;
t177 = t145 * pkin(3) + t201;
t176 = pkin(6) + t186;
t171 = -pkin(6) * qJDD(2) + t218 * t266;
t170 = -t146 * t59 - t149 * t60 - t89 * t225 + t226 * t91;
t104 = t257 * t144;
t125 = qJ(4) * t229;
t167 = -pkin(6) * t245 + pkin(7) * t246;
t169 = qJD(1) * t167 - qJD(5) * t104 + t125 - t221 + t83;
t105 = t257 * t145;
t158 = -pkin(7) * t244 + (-pkin(4) + t214) * t147;
t168 = -qJD(1) * t158 + t144 * qJD(3) - qJD(5) * t105 + t252;
t70 = t94 * t147;
t163 = -pkin(6) + t172;
t153 = qJD(1) ^ 2;
t162 = pkin(1) * t153 + t191;
t161 = t91 * qJ(4) + t189;
t152 = qJD(2) ^ 2;
t160 = pkin(6) * t152 + qJDD(1) * t266 + t259;
t159 = t273 * t262;
t155 = t144 * t217 - t132 + (-t91 + t224) * t220;
t139 = t151 * pkin(6);
t122 = g(1) * t242;
t118 = qJ(3) * t236;
t115 = qJ(3) * t239;
t93 = -qJDD(5) - t164;
t84 = t265 * t145 + t201;
t73 = t89 * t220;
t67 = t176 * t147;
t58 = t186 * t220 + t130;
t54 = pkin(6) * t210 + t252;
t53 = t137 - t61;
t48 = t163 * t147;
t45 = t144 * t216 + t253;
t44 = qJD(1) * t192 - t252;
t43 = t125 + t55;
t42 = t176 * t227 - t208;
t34 = qJD(2) * t192 - t253;
t30 = t81 * t146 + t82 * t149;
t29 = -t82 * t146 + t81 * t149;
t25 = t163 * t227 + t208;
t24 = t46 + t270;
t23 = t73 + t60;
t22 = t89 * pkin(3) - t161;
t21 = qJD(2) * t166 + qJD(5) * t69;
t20 = -t149 * t150 * t224 + qJD(2) * t215 + qJD(5) * t70;
t18 = qJD(2) * t167 + t270 + t65;
t17 = qJD(2) * t158 - t253;
t15 = -t265 * t89 + t161;
t3 = [qJDD(1), -t259 + t262, t191, t142 * qJDD(1) + 0.2e1 * t150 * t202, 0.2e1 * t134 * t147 - 0.2e1 * t218 * t231, qJDD(2) * t147 + t152 * t150, qJDD(2) * t150 - t152 * t147, 0, t171 * t147 + (-t160 + t262) * t150, t147 * t160 + t150 * t171 - t122, (pkin(6) * t59 - t174 * t144 + (qJD(1) * t61 + t39) * qJD(2)) * t147 + (-qJD(1) * t45 + qJD(2) * t188 - qJDD(1) * t61 - t12) * t150 + t193, (pkin(6) * t60 - t174 * t145 + (-qJD(1) * t62 - t40) * qJD(2)) * t147 + (qJD(1) * t46 + qJDD(1) * t62 + t13 + (-t145 * t189 + t264) * qJD(2)) * t150 + t194, -t45 * t91 - t46 * t89 - t62 * t59 - t61 * t60 + t122 + (-t144 * t40 - t145 * t39) * t227 + (-t12 * t145 - t13 * t144 - t259) * t147, t13 * t62 + t40 * t46 + t12 * t61 + t39 * t45 - g(1) * t139 - g(2) * t197 - t159 + (-t147 * t174 - t189 * t227) * pkin(6), t42 * t89 + t67 * t59 + (t10 * t144 + (-qJD(1) * t53 - t26) * qJD(2)) * t147 + (qJD(1) * t34 + qJDD(1) * t53 + t22 * t224 + t11) * t150 + t193, -t24 * t89 + t34 * t91 - t49 * t59 + t53 * t60 + t122 + (-t144 * t31 + t145 * t26) * t227 + (t11 * t145 - t144 * t9 - t259) * t147, -t42 * t91 - t67 * t60 + (-t10 * t145 + (qJD(1) * t49 + t31) * qJD(2)) * t147 + (-qJD(1) * t24 - qJDD(1) * t49 - t22 * t222 - t9) * t150 - t194, t9 * t49 + t31 * t24 + t10 * t67 + t22 * t42 + t11 * t53 + t26 * t34 - g(1) * (-t80 * pkin(3) - t79 * qJ(4) + t139) - g(2) * (t82 * pkin(3) + t81 * qJ(4) + t197) - t159, -t170 * t70 + t37 * t21, -t170 * t69 - t37 * t20 - t21 * t35 - t70 * t8, t21 * t120 - t150 * t170 - t228 * t37 - t70 * t93, -t20 * t120 - t8 * t150 + t228 * t35 - t69 * t93, -t120 * t228 - t93 * t150, (-t146 * t18 + t149 * t17) * t120 + t184 * t93 + t213 * t150 - t1 * t228 + t25 * t35 + t48 * t8 - t6 * t69 + t15 * t20 + g(1) * t181 - g(2) * t30 + (-t120 * t183 - t150 * t2) * qJD(5), -(t146 * t17 + t149 * t18) * t120 + t183 * t93 - t190 * t150 + t2 * t228 + t25 * t37 - t48 * t170 + t6 * t70 + t15 * t21 - g(1) * t182 - g(2) * t29 + (-t1 * t150 + t120 * t184) * qJD(5); 0, 0, 0, -t147 * t153 * t150, t231 * t153, t217, t134, qJDD(2), t147 * t162 - t127 - t141, t258 + (-pkin(6) * qJDD(1) + t162) * t150, -pkin(2) * t59 + t206 * t145 + ((-qJ(3) * t224 - t39) * t147 + (-t188 + t54) * t150) * qJD(1) + t198, t195 - pkin(2) * t60 + t154 * t144 + ((-qJ(3) * t222 + t40) * t147 + (-t264 - t55 + (qJD(3) + t189) * t145) * t150) * qJD(1), t254 + t54 * t91 + t55 * t89 + (t230 * t39 - t191) * t150 + (t220 * t40 - t12 + t180) * t144 + t196, t174 * pkin(2) - t40 * t55 - t39 * t54 + t189 * t130 - g(1) * (-pkin(2) * t241 + t118) - g(2) * (-pkin(2) * t242 + t115) - t205 + (-t39 * t144 + t40 * t145) * qJD(3) + (-t12 * t144 + t254) * qJ(3), -t177 * t59 + (-t58 - t223) * t89 + t207 * t145 + (t147 * t26 - t150 * t44 + (-t150 * t22 - t212) * t144) * qJD(1) + t198, t9 * t145 + t43 * t89 - t44 * t91 + (-t230 * t26 - t191) * t150 + (t220 * t31 + t11 + t180) * t144 + t196, -t195 + t58 * t91 + t177 * t60 + (t207 + t233 + t271) * t144 + (-t147 * t31 + t150 * t43 + (t212 + (-qJD(3) + t22) * t150) * t145) * qJD(1), -t31 * t43 - t22 * t58 - t26 * t44 - g(1) * t118 - g(2) * t115 - t205 + (-g(3) * t137 + t9 * qJ(3) + t31 * qJD(3)) * t145 + (-g(3) * t237 + t11 * qJ(3) + t26 * qJD(3) - t22 * qJD(4)) * t144 + (-t10 + t271) * t177, -t170 * t95 + t256 * t37, t170 * t94 - t255 * t37 - t256 * t35 - t95 * t8, t120 * t256 + t229 * t37 - t95 * t93, -t120 * t255 - t229 * t35 + t94 * t93, t120 * t229, -(t149 * t104 - t146 * t105) * t93 + t84 * t8 + t199 * t35 - g(3) * t166 + t255 * t15 + (t146 * t169 + t149 * t168) * t120 + t1 * t229 + t267 * t94, (t146 * t104 + t149 * t105) * t93 - t84 * t170 + t199 * t37 + t256 * t15 + (-t146 * t168 + t149 * t169) * t120 - t2 * t229 + (-t141 + t267) * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t23, t272, t39 * t91 + t40 * t89 + t154, t155, t272, -t23, t263 - t249 + t31 * t89 + (-qJD(4) - t26) * t91 + t154, 0, 0, 0, 0, 0, -t8 - t250, t170 + t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91 * t89 + t164, -t73 + t60, -t143 * t153 - t87, -g(3) * t247 - g(1) * t81 - g(2) * t79 + t22 * t91 + (-pkin(3) * t228 + t150 * t31) * qJD(1) + t175, 0, 0, 0, 0, 0, -t146 * t179 - t149 * t93 - t91 * t35, t146 * t93 - t149 * t179 - t91 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t35, -t35 ^ 2 + t37 ^ 2, -t170 + t251, -t8 + t250, -t93, -g(1) * t29 + g(2) * t182 - g(3) * t69 - t15 * t37 + t2 * t269 + t213, g(1) * t30 + g(2) * t181 + g(3) * t70 + t1 * t269 + t15 * t35 - t190;];
tau_reg = t3;
