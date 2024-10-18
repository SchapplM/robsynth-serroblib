% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRR11
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR11_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:57
% EndTime: 2024-09-27 21:46:00
% DurationCPUTime: 2.09s
% Computational Cost: add. (2773->318), mult. (7355->446), div. (0->0), fcn. (5857->22), ass. (0->207)
t164 = qJ(3) + qJ(4);
t216 = pkin(5) - t164;
t199 = -qJ(5) + t216;
t186 = sin(t199);
t167 = sin(qJ(5));
t170 = cos(qJ(5));
t171 = cos(qJ(4));
t165 = sin(pkin(5));
t172 = cos(qJ(3));
t243 = qJD(2) * t172;
t223 = t165 * t243;
t201 = t171 * t223;
t168 = sin(qJ(4));
t169 = sin(qJ(3));
t245 = qJD(2) * t165;
t224 = t169 * t245;
t202 = t168 * t224;
t90 = -t201 + t202;
t92 = -t168 * t223 - t171 * t224;
t194 = t167 * t90 + t170 * t92;
t166 = cos(pkin(5));
t139 = qJDD(2) * t166 + qJDD(3);
t125 = qJDD(4) + t139;
t244 = qJD(2) * t166;
t143 = qJD(3) + t244;
t278 = pkin(7) + pkin(8);
t227 = t278 * t169;
t204 = t165 * t227;
t259 = t165 * t172;
t138 = qJD(1) * t259;
t257 = t166 * t172;
t146 = pkin(2) * t257;
t249 = qJD(2) * t146 + t138;
t63 = -qJD(2) * t204 + t249;
t52 = pkin(3) * t143 + t63;
t258 = t166 * t169;
t145 = pkin(2) * t258;
t246 = qJD(1) * t169;
t64 = t165 * t246 + (t278 * t259 + t145) * qJD(2);
t58 = t171 * t64;
t193 = -t168 * t52 - t58;
t276 = pkin(2) * t166;
t228 = qJD(3) * t276;
t203 = qJD(2) * t228;
t226 = t278 * t172;
t234 = qJDD(2) * t172;
t219 = t166 * t234;
t235 = qJDD(1) * t165;
t250 = pkin(2) * t219 + t172 * t235;
t26 = -t169 * t203 + pkin(3) * t139 + (-qJDD(2) * t227 + (-qJD(2) * t226 - t246) * qJD(3)) * t165 + t250;
t220 = t165 * t234;
t191 = pkin(7) * t220 + qJD(3) * t138 + qJDD(2) * t145 + t169 * t235 + t172 * t203;
t236 = qJD(2) * qJD(3);
t222 = t169 * t236;
t30 = (-pkin(7) * t222 + (-t222 + t234) * pkin(8)) * t165 + t191;
t213 = -t168 * t30 + t171 * t26;
t180 = t193 * qJD(4) + t213;
t279 = qJD(3) + qJD(4);
t233 = t169 * qJDD(2);
t217 = t165 * t233;
t221 = t172 * t236;
t282 = t165 * t221 + t217;
t35 = qJD(4) * t201 + t168 * t220 + t282 * t171 - t279 * t202;
t2 = pkin(4) * t125 - pkin(9) * t35 + t180;
t240 = qJD(4) * t171;
t241 = qJD(4) * t168;
t208 = -t168 * t26 - t171 * t30 - t52 * t240 + t64 * t241;
t189 = t168 * t172 + t169 * t171;
t176 = t279 * t189;
t36 = (t176 * qJD(2) + t168 * t233) * t165 - t171 * t220;
t3 = -pkin(9) * t36 - t208;
t225 = -t167 * t3 + t170 * t2;
t156 = pkin(5) + t164;
t148 = qJ(5) + t156;
t230 = sin(t148) / 0.2e1;
t108 = (-pkin(3) * t172 - pkin(2)) * t165;
t153 = t166 * qJD(1);
t93 = qJD(2) * t108 + t153;
t51 = pkin(4) * t90 + t93;
t159 = qJ(5) + t164;
t149 = sin(t159);
t161 = pkin(10) + qJ(2);
t154 = sin(t161);
t155 = cos(t161);
t187 = cos(t199);
t231 = cos(t148) / 0.2e1;
t182 = t187 / 0.2e1 + t231;
t59 = t149 * t154 - t155 * t182;
t61 = -t155 * t149 - t154 * t182;
t287 = -g(1) * t61 + g(2) * t59 - g(3) * (t230 + t186 / 0.2e1) + t225 + t51 * t194;
t277 = pkin(9) * t90;
t18 = -t193 - t277;
t239 = qJD(5) * t167;
t16 = t18 * t239;
t44 = t167 * t92 - t170 * t90;
t102 = t230 - t186 / 0.2e1;
t150 = cos(t159);
t60 = -t102 * t155 - t150 * t154;
t62 = t102 * t154 - t150 * t155;
t286 = -g(1) * t62 - g(2) * t60 - g(3) * (t231 - t187 / 0.2e1) + t16 - t51 * t44;
t133 = qJD(4) + t143;
t56 = t168 * t64;
t212 = t171 * t52 - t56;
t86 = t92 * pkin(9);
t17 = t212 + t86;
t15 = pkin(4) * t133 + t17;
t263 = t170 * t18;
t195 = -t15 * t167 - t263;
t285 = t195 * qJD(5) + t287;
t122 = qJD(5) + t133;
t284 = (-t18 * t122 - t2) * t167 + t286;
t272 = t194 * t44;
t99 = (-t168 * t169 + t171 * t172) * t165;
t6 = t194 ^ 2 - t44 ^ 2;
t238 = qJD(5) * t170;
t8 = -t167 * t36 + t170 * t35 - t90 * t238 + t92 * t239;
t4 = -t122 * t44 + t8;
t179 = t194 * qJD(5) - t167 * t35 - t170 * t36;
t5 = -t122 * t194 + t179;
t281 = (qJD(3) * t143 * t172 + t139 * t169) * t165;
t275 = pkin(3) * t168;
t274 = pkin(7) * t169;
t273 = t166 * t179;
t269 = t92 * t90;
t115 = qJDD(5) + t125;
t100 = t189 * t165;
t50 = t100 * t170 + t167 * t99;
t53 = t279 * t99;
t54 = t176 * t165;
t14 = t50 * qJD(5) + t167 * t53 + t170 * t54;
t49 = t100 * t167 - t170 * t99;
t268 = -t49 * t115 - t14 * t122;
t267 = t99 * t125 - t54 * t133;
t266 = t171 * t63 - t56;
t79 = pkin(3) * t166 + t146 - t204;
t144 = pkin(7) * t259;
t248 = t144 + t145;
t88 = pkin(8) * t259 + t248;
t265 = t168 * t79 + t171 * t88;
t264 = t166 * t36;
t262 = pkin(2) * qJDD(2);
t111 = -pkin(2) * t245 + t153;
t261 = t111 * t172;
t160 = t165 ^ 2;
t260 = t160 * qJD(2) ^ 2;
t256 = t167 * t115;
t254 = t168 * t170;
t253 = t170 * t115;
t251 = t282 * t166;
t162 = t169 ^ 2;
t247 = -t172 ^ 2 + t162;
t242 = qJD(3) * t169;
t237 = qJD(3) - t143;
t142 = cos(t156) / 0.2e1;
t147 = cos(t216);
t232 = t147 / 0.2e1 + t142;
t229 = sin(t156) / 0.2e1;
t215 = qJD(5) * t15 + t3;
t211 = -t168 * t63 - t58;
t210 = -t168 * t88 + t171 * t79;
t135 = t172 * t228;
t80 = -qJD(3) * t204 + t135;
t81 = (-t165 * t226 - t145) * qJD(3);
t209 = -t168 * t80 + t171 * t81;
t207 = t143 + t244;
t206 = t237 * qJD(2);
t205 = pkin(3) * t165 * t242;
t198 = sin(t216);
t197 = t169 * t206;
t13 = -t49 * qJD(5) - t167 * t54 + t170 * t53;
t196 = -t115 * t50 - t122 * t13;
t190 = -t100 * t125 - t133 * t53;
t105 = t229 - t198 / 0.2e1;
t158 = cos(t164);
t70 = -t105 * t155 - t154 * t158;
t72 = t105 * t154 - t155 * t158;
t185 = t168 * t81 + t171 * t80 + t79 * t240 - t88 * t241;
t152 = t166 * qJDD(1);
t68 = qJD(2) * t205 + qJDD(2) * t108 + t152;
t178 = -g(1) * t72 - g(2) * t70 - g(3) * (t142 - t147 / 0.2e1) + t93 * t90 + t208;
t157 = sin(t164);
t69 = t154 * t157 - t155 * t232;
t71 = -t154 * t232 - t155 * t157;
t174 = -g(1) * t71 + g(2) * t69 - g(3) * (t229 + t198 / 0.2e1) + t93 * t92 + t180;
t151 = pkin(3) * t171 + pkin(4);
t109 = -t165 * t262 + t152;
t107 = t139 * t259;
t97 = -t154 * t258 + t155 * t172;
t96 = -t154 * t257 - t155 * t169;
t95 = -t154 * t172 - t155 * t258;
t94 = t154 * t169 - t155 * t257;
t66 = pkin(3) * t224 - pkin(4) * t92;
t65 = -pkin(4) * t99 + t108;
t39 = pkin(4) * t54 + t205;
t37 = -t90 ^ 2 + t92 ^ 2;
t31 = t35 * t166;
t29 = pkin(9) * t99 + t265;
t27 = pkin(4) * t166 - pkin(9) * t100 + t210;
t23 = t86 + t266;
t22 = t211 + t277;
t21 = -t133 * t92 - t36;
t20 = t133 * t90 + t35;
t19 = pkin(4) * t36 + t68;
t11 = -pkin(9) * t53 - qJD(4) * t265 + t209;
t10 = -pkin(9) * t54 + t185;
t7 = t8 * t166;
t1 = [qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, t107 + (-t219 + (-t143 + t244) * t242) * t165, t251 - t281, 0, 0, 0, 0, 0, t264 + t267, t190 + t31, 0, 0, 0, 0, 0, t268 - t273, t196 + t7; 0, qJDD(2), g(1) * t154 - g(2) * t155, g(1) * t155 + g(2) * t154, (qJDD(2) * t162 + 0.2e1 * t169 * t221) * t160, 0.2e1 * (t172 * t233 - t247 * t236) * t160, t251 + t281, t107 + (-t207 * t242 + t219) * t165, t139 * t166, (-t165 * t274 + t146) * t139 + (-pkin(7) * t217 + t250) * t166 - g(1) * t95 - g(2) * t97 + (-t109 * t165 + t160 * t262) * t172 + (-t207 * t144 + ((t111 - t153) * t165 + (-t166 * t143 + (-t166 ^ 2 - t160) * qJD(2)) * pkin(2)) * t169) * qJD(3), -t135 * t143 - t248 * t139 - t191 * t166 - g(1) * t94 - g(2) * t96 + (t109 * t169 + (t207 * t274 + t261) * qJD(3)) * t165 + (-t221 - t233) * t160 * pkin(2), t100 * t35 - t53 * t92, -t100 * t36 + t35 * t99 - t53 * t90 + t54 * t92, -t190 + t31, -t264 + t267, t125 * t166, t209 * t133 + t210 * t125 + t213 * t166 + t90 * t205 + t108 * t36 - t68 * t99 + t93 * t54 - g(1) * t70 + g(2) * t72 + (-t133 * t265 + t193 * t166) * qJD(4), -g(1) * t69 - g(2) * t71 + t68 * t100 + t108 * t35 - t265 * t125 - t185 * t133 + t208 * t166 - t92 * t205 + t93 * t53, -t13 * t194 + t50 * t8, t13 * t44 + t14 * t194 + t179 * t50 - t49 * t8, -t196 + t7, t268 + t273, t115 * t166, (-t10 * t167 + t11 * t170) * t122 + (-t167 * t29 + t170 * t27) * t115 + t225 * t166 - t39 * t44 - t65 * t179 + t19 * t49 + t51 * t14 - g(1) * t60 + g(2) * t62 + ((-t167 * t27 - t170 * t29) * t122 + t195 * t166) * qJD(5), -g(1) * t59 - g(2) * t61 + t51 * t13 + t16 * t166 + t19 * t50 - t39 * t194 + t65 * t8 + (-(-qJD(5) * t29 + t11) * t122 - t27 * t115 - t2 * t166) * t167 + (-(qJD(5) * t27 + t10) * t122 - t29 * t115 - t215 * t166) * t170; 0, 0, 0, 0, -t169 * t172 * t260, t247 * t260, (t237 * t243 + t233) * t165, (-t197 + t234) * t165, t139, -g(1) * t96 + g(2) * t94 - t197 * t276 + ((-pkin(7) * t206 - g(3)) * t172 + (-qJDD(2) * pkin(7) - t237 * qJD(1) - t111 * qJD(2)) * t169) * t165 + t250, t249 * t143 + g(1) * t97 - g(2) * t95 + (g(3) * t169 + (t237 * t274 - t261) * qJD(2)) * t165 - t191, -t269, t37, t20, t21, t125, -t211 * t133 + (t171 * t125 - t133 * t241 - t90 * t224) * pkin(3) + t174, t266 * t133 + (-t168 * t125 - t133 * t240 + t92 * t224) * pkin(3) + t178, t272, t6, t4, t5, t115, t151 * t253 - (-t167 * t23 + t170 * t22) * t122 + t66 * t44 + (-t168 * t256 + (-t167 * t171 - t254) * t122 * qJD(4)) * pkin(3) + ((-pkin(3) * t254 - t151 * t167) * t122 + t195) * qJD(5) + t287, t66 * t194 + (-t151 * t115 - t2 + (t22 - (-qJD(4) - qJD(5)) * t275) * t122) * t167 + (-t115 * t275 + (-pkin(3) * t240 - qJD(5) * t151 + t23) * t122 - t215) * t170 + t286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t269, t37, t20, t21, t125, -t193 * t133 + t174, t212 * t133 + t178, t272, t6, t4, t5, t115, -(-t167 * t17 - t263) * t122 + (-t122 * t239 - t44 * t92 + t253) * pkin(4) + t285, (t17 * t122 - t215) * t170 + (-t122 * t238 - t194 * t92 - t256) * pkin(4) + t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272, t6, t4, t5, t115, -t195 * t122 + t285, (-t3 + (-qJD(5) + t122) * t15) * t170 + t284;];
tau_reg = t1;
