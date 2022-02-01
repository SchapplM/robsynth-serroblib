% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:49:28
% EndTime: 2022-01-20 11:49:32
% DurationCPUTime: 1.86s
% Computational Cost: add. (3288->303), mult. (4938->357), div. (0->0), fcn. (3237->12), ass. (0->204)
t275 = cos(qJ(4));
t161 = sin(qJ(4));
t165 = cos(qJ(3));
t156 = qJD(1) + qJD(2);
t163 = sin(qJ(2));
t274 = pkin(1) * t163;
t226 = qJD(1) * t274;
t102 = pkin(7) * t156 + t226;
t212 = pkin(8) * t156 + t102;
t63 = t212 * t165;
t55 = t161 * t63;
t162 = sin(qJ(3));
t62 = t212 * t162;
t58 = qJD(3) * pkin(3) - t62;
t211 = t275 * t58 - t55;
t220 = t275 * t162;
t94 = t161 * t165 + t220;
t72 = t94 * t156;
t251 = t72 * qJ(5);
t19 = t211 - t251;
t234 = qJD(2) * t163;
t143 = pkin(1) * t234;
t166 = cos(qJ(2));
t273 = pkin(1) * t166;
t237 = -qJD(1) * t143 + qJDD(1) * t273;
t160 = qJ(1) + qJ(2);
t148 = cos(t160);
t269 = g(2) * t148;
t153 = qJDD(1) + qJDD(2);
t272 = pkin(2) * t153;
t215 = t237 + t272 - t269;
t214 = qJD(4) * t275;
t219 = t275 * t165;
t288 = -qJD(3) * t219 - t165 * t214;
t150 = t165 * pkin(3);
t266 = pkin(2) + t150;
t277 = -pkin(7) - pkin(8);
t221 = qJD(3) * t277;
t100 = t165 * t221;
t235 = qJD(1) * t166;
t224 = pkin(1) * t235;
t120 = t277 * t162;
t149 = t165 * pkin(8);
t121 = pkin(7) * t165 + t149;
t258 = t161 * t120 + t275 * t121;
t99 = t162 * t221;
t287 = -t258 * qJD(4) + t275 * t100 - t161 * t99 + t94 * t224;
t243 = t161 * t162;
t185 = t219 - t243;
t230 = qJD(4) * t161;
t286 = -t161 * t100 - t120 * t214 + t121 * t230 + t185 * t224 - t275 * t99;
t231 = qJD(3) * t165;
t216 = t156 * t231;
t242 = t162 * t153;
t229 = qJDD(1) * t163;
t233 = qJD(2) * t166;
t76 = t153 * pkin(7) + (qJD(1) * t233 + t229) * pkin(1);
t29 = -t102 * t231 + qJDD(3) * pkin(3) - t162 * t76 + (-t216 - t242) * pkin(8);
t232 = qJD(3) * t162;
t217 = t156 * t232;
t241 = t165 * t153;
t30 = -t102 * t232 + t165 * t76 + (-t217 + t241) * pkin(8);
t285 = -t161 * t30 + t275 * t29;
t152 = qJDD(3) + qJDD(4);
t155 = qJD(3) + qJD(4);
t270 = pkin(3) * t155;
t284 = -t161 * pkin(3) * t152 - t214 * t270;
t146 = sin(t160);
t283 = g(1) * t148 + g(2) * t146;
t142 = pkin(3) * t232;
t282 = t142 - t226;
t144 = t152 * pkin(4);
t195 = t155 * t243;
t203 = -t153 * t220 + t288 * t156 - t161 * t241;
t25 = t156 * t195 + t203;
t255 = t25 * qJ(5);
t281 = t144 + t255;
t133 = g(1) * t146;
t280 = t269 - t133;
t159 = qJ(3) + qJ(4);
t145 = sin(t159);
t247 = t145 * t148;
t248 = t145 * t146;
t147 = cos(t159);
t268 = g(3) * t147;
t279 = g(1) * t247 + g(2) * t248 - t268;
t278 = t72 ^ 2;
t276 = pkin(4) * t185;
t271 = pkin(2) * t156;
t223 = t156 * t243;
t70 = -t156 * t219 + t223;
t267 = t72 * t70;
t137 = pkin(7) + t274;
t265 = -pkin(8) - t137;
t47 = t155 * t94;
t261 = -t47 * qJ(5) + qJD(5) * t185;
t264 = t261 - t286;
t46 = t195 + t288;
t194 = t46 * qJ(5) - t94 * qJD(5);
t263 = t194 + t287;
t18 = pkin(4) * t155 + t19;
t262 = t18 - t19;
t260 = -t275 * t62 - t55;
t90 = t265 * t162;
t91 = t137 * t165 + t149;
t259 = t161 * t90 + t275 * t91;
t257 = qJ(5) * t94;
t196 = -t153 * t219 + t161 * t242;
t26 = t47 * t156 + t196;
t254 = t26 * qJ(5);
t253 = t70 * qJ(5);
t252 = t70 * t155;
t103 = -t224 - t271;
t249 = t103 * t232 + t165 * t133;
t246 = t146 * t147;
t245 = t147 * t148;
t244 = t156 * t162;
t213 = t70 * pkin(4) + qJD(5);
t74 = -t156 * t266 - t224;
t40 = t213 + t74;
t240 = qJD(5) + t40;
t238 = pkin(4) * t147 + t150;
t157 = t162 ^ 2;
t236 = -t165 ^ 2 + t157;
t228 = pkin(3) * t244;
t227 = t103 * t231 - t215 * t162;
t225 = pkin(1) * t233;
t57 = t275 * t63;
t139 = -pkin(2) - t273;
t218 = t156 * t234;
t43 = pkin(4) * t47 + t142;
t210 = t161 * t62 - t57;
t209 = -t161 * t91 + t275 * t90;
t208 = qJD(3) * t265;
t207 = t275 * t120 - t121 * t161;
t154 = qJ(5) - t277;
t98 = pkin(2) + t238;
t206 = -t146 * t98 + t154 * t148;
t205 = t146 * t154 + t148 * t98;
t204 = t156 * t226;
t115 = t139 - t150;
t200 = t43 - t226;
t199 = -g(1) * t248 + g(2) * t247;
t198 = g(1) * t246 - g(2) * t245;
t190 = -t161 * t58 - t57;
t20 = -t190 - t253;
t177 = t190 * qJD(4) + t285;
t3 = -t72 * qJD(5) + t177 + t281;
t175 = t161 * t29 + t58 * t214 - t63 * t230 + t275 * t30;
t4 = -t70 * qJD(5) + t175 - t254;
t193 = t18 * t46 + t185 * t4 - t20 * t47 - t3 * t94 - t283;
t192 = -t237 + t280;
t44 = pkin(3) * t217 - t153 * t266 - t237;
t15 = t26 * pkin(4) + qJDD(5) + t44;
t191 = -t15 * t185 + t40 * t47 + t198;
t189 = t15 * t94 - t40 * t46 + t199;
t188 = t44 * t94 - t74 * t46 + t199;
t187 = -t185 * t44 + t74 * t47 + t198;
t59 = t162 * t208 + t165 * t225;
t60 = -t162 * t225 + t165 * t208;
t186 = t161 * t60 + t90 * t214 - t91 * t230 + t275 * t59;
t183 = -t103 * t156 + t283 - t76;
t168 = qJD(3) ^ 2;
t182 = pkin(7) * t168 - t204 - t272;
t181 = pkin(1) * t218 + t137 * t168 + t139 * t153;
t180 = -t155 * t223 - t203;
t179 = -pkin(7) * qJDD(3) + (t224 - t271) * qJD(3);
t178 = -qJDD(3) * t137 + (t139 * t156 - t225) * qJD(3);
t176 = -t259 * qJD(4) - t161 * t59 + t275 * t60;
t173 = g(1) * t245 + g(2) * t246 + g(3) * t145 - t175;
t172 = t177 + t279;
t171 = t74 * t70 + t173;
t170 = -t74 * t72 + t172;
t169 = t240 * t70 + t173 + t254;
t167 = cos(qJ(1));
t164 = sin(qJ(1));
t151 = t156 ^ 2;
t138 = t275 * pkin(3) + pkin(4);
t114 = qJDD(3) * t165 - t162 * t168;
t113 = qJDD(3) * t162 + t165 * t168;
t101 = t143 + t142;
t89 = t185 * qJ(5);
t77 = t153 * t157 + 0.2e1 * t162 * t216;
t69 = -t266 - t276;
t68 = t70 ^ 2;
t61 = t115 - t276;
t49 = -0.2e1 * t236 * t156 * qJD(3) + 0.2e1 * t162 * t241;
t48 = pkin(4) * t72 + t228;
t42 = t89 + t258;
t41 = t207 - t257;
t39 = t143 + t43;
t34 = t89 + t259;
t33 = t209 - t257;
t32 = t152 * t185 - t155 * t47;
t31 = t152 * t94 - t155 * t46;
t28 = -t68 + t278;
t22 = -t251 + t260;
t21 = t210 + t253;
t16 = t180 + t252;
t8 = -t25 * t94 - t46 * t72;
t7 = t176 + t194;
t6 = t186 + t261;
t5 = -t185 * t25 - t26 * t94 + t46 * t70 - t47 * t72;
t1 = [qJDD(1), g(1) * t164 - g(2) * t167, g(1) * t167 + g(2) * t164, t153, (t153 * t166 - t218) * pkin(1) - t192, ((-qJDD(1) - t153) * t163 + (-qJD(1) - t156) * t233) * pkin(1) + t283, t77, t49, t113, t114, 0, t178 * t162 + (-t181 + t215) * t165 + t249, t178 * t165 + (t181 - t133) * t162 + t227, t8, t5, t31, t32, 0, t101 * t70 + t115 * t26 + t209 * t152 + t176 * t155 + t187, t101 * t72 - t115 * t25 - t259 * t152 - t186 * t155 + t188, t152 * t33 + t155 * t7 + t26 * t61 + t39 * t70 + t191, -t152 * t34 - t155 * t6 - t25 * t61 + t39 * t72 + t189, t25 * t33 - t26 * t34 - t6 * t70 - t7 * t72 + t193, t4 * t34 + t20 * t6 + t3 * t33 + t18 * t7 + t15 * t61 + t40 * t39 - g(1) * (-pkin(1) * t164 + t206) - g(2) * (pkin(1) * t167 + t205); 0, 0, 0, t153, -t192 + t204, (-t229 + (-qJD(2) + t156) * t235) * pkin(1) + t283, t77, t49, t113, t114, 0, t179 * t162 + (-t182 + t215) * t165 + t249, t179 * t165 + (t182 - t133) * t162 + t227, t8, t5, t31, t32, 0, t207 * t152 + t287 * t155 - t26 * t266 + t282 * t70 + t187, -t258 * t152 + t286 * t155 + t25 * t266 + t282 * t72 + t188, t152 * t41 + t263 * t155 + t200 * t70 + t26 * t69 + t191, -t152 * t42 - t264 * t155 + t200 * t72 - t25 * t69 + t189, t25 * t41 - t26 * t42 - t263 * t72 - t264 * t70 + t193, -g(1) * t206 - g(2) * t205 + t15 * t69 + t263 * t18 + t264 * t20 + t200 * t40 + t3 * t41 + t4 * t42; 0, 0, 0, 0, 0, 0, -t162 * t151 * t165, t236 * t151, t242, t241, qJDD(3), -g(3) * t165 + t162 * t183, g(3) * t162 + t165 * t183, t267, t28, t16, -t196, t152, -t210 * t155 + (t275 * t152 - t155 * t230 - t70 * t244) * pkin(3) + t170, t260 * t155 - t72 * t228 + t171 + t284, t138 * t152 - t21 * t155 - t48 * t70 - t240 * t72 + (-t57 + (-t58 - t270) * t161) * qJD(4) + t279 + t281 + t285, t22 * t155 - t48 * t72 + t169 + t284, t138 * t25 + (t20 + t21) * t72 + (-t18 + t22) * t70 + (-t161 * t26 + (t161 * t72 - t275 * t70) * qJD(4)) * pkin(3), t3 * t138 - t20 * t22 - t18 * t21 - t40 * t48 - g(3) * t238 - t283 * (-pkin(3) * t162 - pkin(4) * t145) + (t4 * t161 + (-t161 * t18 + t275 * t20) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t267, t28, t16, -t196, t152, -t155 * t190 + t170, t211 * t155 + t171, t255 + t20 * t155 + 0.2e1 * t144 + (-t213 - t40) * t72 + t172, -t278 * pkin(4) + t19 * t155 + t169, t25 * pkin(4) - t262 * t70, t262 * t20 + (t145 * t283 - t40 * t72 - t268 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t155 + t26, t180 - t252, -t68 - t278, t18 * t72 + t20 * t70 + t15 + t280;];
tau_reg = t1;
