% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:07:15
% EndTime: 2019-03-09 05:07:22
% DurationCPUTime: 3.21s
% Computational Cost: add. (2954->397), mult. (6783->517), div. (0->0), fcn. (4334->8), ass. (0->215)
t159 = cos(qJ(3));
t238 = qJD(1) * t159;
t298 = qJD(4) - t238;
t297 = qJD(6) - qJD(4);
t155 = sin(qJ(4));
t158 = cos(qJ(4));
t227 = t158 * qJD(3);
t156 = sin(qJ(3));
t239 = qJD(1) * t156;
t107 = t155 * t239 - t227;
t237 = qJD(3) * t155;
t109 = t158 * t239 + t237;
t154 = sin(qJ(6));
t157 = cos(qJ(6));
t181 = -t157 * t107 + t109 * t154;
t51 = t107 * t154 + t109 * t157;
t296 = -t181 ^ 2 + t51 ^ 2;
t232 = qJD(4) * t158;
t210 = t156 * t232;
t235 = qJD(3) * t159;
t166 = t155 * t235 + t210;
t224 = qJD(3) * qJD(4);
t68 = t166 * qJD(1) + t155 * t224;
t128 = t298 * qJD(5);
t225 = qJD(1) * qJD(3);
t141 = t156 * t225;
t134 = qJ(5) * t141;
t139 = sin(pkin(10)) * pkin(1) + pkin(7);
t124 = t139 * qJD(1);
t289 = t159 * qJD(2) - t156 * t124;
t77 = t289 * qJD(3);
t140 = -cos(pkin(10)) * pkin(1) - pkin(2);
t100 = -pkin(3) * t159 - pkin(8) * t156 + t140;
t79 = t100 * qJD(1);
t192 = pkin(3) * t156 - pkin(8) * t159;
t118 = t192 * qJD(3);
t99 = qJD(1) * t118;
t222 = t155 * t99 + t158 * t77 + t79 * t232;
t233 = qJD(4) * t155;
t146 = t156 * qJD(2);
t86 = t159 * t124 + t146;
t76 = qJD(3) * pkin(8) + t86;
t173 = t76 * t233 - t222;
t9 = t128 + t134 - t173;
t6 = pkin(9) * t68 + t9;
t203 = t155 * t77 - t158 * t99 + t76 * t232 + t79 * t233;
t283 = pkin(4) + pkin(5);
t207 = t158 * t225;
t211 = t156 * t233;
t67 = qJD(1) * t211 - t158 * t224 - t159 * t207;
t7 = pkin(9) * t67 - t283 * t141 + t203;
t215 = -t154 * t6 + t157 * t7;
t194 = qJD(3) * pkin(3) + t289;
t170 = qJ(5) * t109 + t194;
t24 = -t283 * t107 + t170;
t295 = t24 * t51 - t215;
t293 = t24 * t181;
t292 = t51 * t181;
t197 = pkin(4) * t141;
t10 = -t197 + t203;
t131 = t298 * qJ(5);
t36 = t155 * t79 + t158 * t76;
t30 = t131 + t36;
t291 = -t298 * t30 + t10;
t290 = qJD(5) * t155 + t86;
t35 = -t155 * t76 + t158 * t79;
t241 = qJD(5) - t35;
t208 = t159 * t227;
t167 = t208 - t211;
t226 = -qJD(6) + t298;
t288 = -qJD(6) - t226;
t228 = qJD(6) * t157;
t229 = qJD(6) * t154;
t12 = t107 * t228 - t109 * t229 + t154 * t68 - t157 * t67;
t286 = t181 * t226 - t12;
t258 = qJ(5) * t155;
t285 = -t283 * t158 - t258;
t284 = t109 ^ 2;
t282 = pkin(8) - pkin(9);
t281 = pkin(9) * t156;
t257 = qJ(5) * t158;
t188 = pkin(4) * t155 - t257;
t280 = -t298 * t188 + t290;
t246 = t156 * t158;
t279 = -t107 * t208 - t68 * t246;
t111 = t154 * t155 + t157 * t158;
t172 = t159 * t111;
t278 = -qJD(1) * t172 - t111 * t297;
t214 = t155 * t238;
t249 = t154 * t158;
t277 = t154 * t232 + t155 * t228 - t158 * t229 - t238 * t249 + (t214 - t233) * t157;
t174 = -t283 * t155 + t257;
t276 = t298 * t174 + t290;
t115 = t192 * qJD(1);
t275 = t155 * t115 + t158 * t289;
t274 = t100 * t232 + t155 * t118;
t34 = pkin(4) * t107 - t170;
t273 = t109 * t34;
t204 = -t154 * t67 - t157 * t68;
t13 = t51 * qJD(6) + t204;
t272 = t13 * t159;
t248 = t155 * t157;
t180 = -t248 + t249;
t28 = -t156 * t180 * t297 + qJD(3) * t172;
t271 = t226 * t28;
t26 = pkin(9) * t107 + t36;
t20 = t131 + t26;
t269 = t154 * t20;
t78 = qJD(3) * t146 + t124 * t235;
t168 = -qJ(5) * t67 + qJD(5) * t109 - t78;
t16 = pkin(4) * t68 - t168;
t268 = t155 * t16;
t267 = t155 * t194;
t266 = t158 * t16;
t265 = t158 * t194;
t264 = t159 * t67;
t263 = t159 * t68;
t262 = t78 * t155;
t261 = t78 * t158;
t245 = t158 * t159;
t114 = t139 * t245;
t260 = t155 * t100 + t114;
t259 = qJ(5) * t107;
t89 = t154 * t246 - t156 * t248;
t256 = qJD(1) * t89;
t90 = t111 * t156;
t255 = qJD(1) * t90;
t254 = t107 * t298;
t253 = t109 * t107;
t252 = t109 * t298;
t251 = t298 * t158;
t250 = t139 * t155;
t247 = t155 * t159;
t161 = qJD(3) ^ 2;
t244 = t161 * t156;
t243 = t161 * t159;
t242 = pkin(9) * t109 - t241;
t150 = t156 ^ 2;
t240 = -t159 ^ 2 + t150;
t125 = qJD(1) * t140;
t236 = qJD(3) * t156;
t234 = qJD(4) * t107;
t230 = qJD(5) * t158;
t18 = -t283 * t298 - t242;
t223 = t154 * t7 + t157 * t6 + t18 * t228;
t221 = pkin(9) * t245;
t220 = t155 * pkin(8) * t298;
t219 = pkin(8) * t251;
t218 = qJ(5) * t236 + t274;
t40 = qJ(5) * t239 + t275;
t217 = pkin(8) * t236;
t216 = pkin(8) * t227;
t127 = t282 * t158;
t213 = t109 * t235;
t212 = t298 * t233;
t205 = -pkin(4) - t250;
t202 = t115 * t158 - t155 * t289;
t129 = t150 * t207;
t201 = t129 - t264;
t200 = -t67 + t234;
t113 = t139 * t247;
t199 = t100 * t158 - t113;
t198 = t226 ^ 2;
t196 = t109 * t210;
t195 = t298 * t210;
t45 = -qJ(5) * t159 + t260;
t193 = -qJD(4) * t114 - t100 * t233 + t118 * t158;
t126 = t282 * t155;
t191 = pkin(9) * t214 - qJD(6) * t126 + t282 * t233 + t40;
t190 = (-t283 * t156 - t221) * qJD(1) - t202 + t297 * t127;
t189 = pkin(4) * t158 + t258;
t2 = t154 * t18 + t157 * t20;
t149 = t159 * pkin(4);
t37 = pkin(5) * t159 + t113 + t149 + (-t100 - t281) * t158;
t43 = t155 * t281 + t45;
t187 = -t154 * t43 + t157 * t37;
t186 = t154 * t37 + t157 * t43;
t29 = -pkin(4) * t298 + t241;
t185 = -t155 * t30 + t158 * t29;
t184 = t155 * t29 + t158 * t30;
t183 = qJ(5) * t157 - t154 * t283;
t182 = qJ(5) * t154 + t157 * t283;
t178 = 0.2e1 * qJD(3) * t125;
t177 = t139 + t188;
t176 = -t20 * t229 + t223;
t175 = t298 * t36 - t203;
t171 = (-qJD(1) * t150 - t159 * t298) * t155;
t169 = -t139 + t174;
t165 = t298 * t35 + t173;
t164 = t185 * qJD(4) + t10 * t155 + t158 * t9;
t163 = qJD(3) * t171 + t107 * t236 - t195 - t263;
t162 = qJD(1) ^ 2;
t119 = -pkin(3) - t189;
t103 = t298 * t208;
t101 = pkin(3) - t285;
t92 = t109 * t236;
t69 = t177 * t156;
t56 = pkin(4) * t109 + t259;
t52 = t169 * t156;
t46 = t149 - t199;
t42 = -pkin(4) * t239 - t202;
t39 = -t283 * t109 - t259;
t38 = -t67 + t254;
t32 = (t189 * qJD(4) - t230) * t156 + t177 * t235;
t27 = qJD(6) * t90 + t167 * t154 - t166 * t157;
t23 = t205 * t236 - t193;
t22 = (t285 * qJD(4) + t230) * t156 + t169 * t235;
t21 = t27 * t226;
t19 = -qJD(5) * t159 + (-t156 * t227 - t159 * t233) * t139 + t218;
t15 = (pkin(9) * qJD(4) - qJD(3) * t139) * t246 + (-qJD(5) + (pkin(9) * qJD(3) - qJD(4) * t139) * t155) * t159 + t218;
t14 = pkin(9) * t211 + (-t221 + (-pkin(5) + t205) * t156) * qJD(3) - t193;
t11 = t12 * t159;
t8 = -t283 * t68 + t168;
t1 = t157 * t18 - t269;
t3 = [0, 0, 0, 0, 0.2e1 * t159 * t141, -0.2e1 * t240 * t225, t243, -t244, 0, -t139 * t243 + t156 * t178, t139 * t244 + t159 * t178, t167 * t109 - t67 * t246, -t196 + (-t213 + (t67 + t234) * t156) * t155 + t279, -t211 * t298 + t103 + t129 + t264 + t92, -t195 + t263 + (-t107 * t156 + t171) * qJD(3) (t298 - t238) * t236, t193 * t298 + ((t107 * t139 - t267) * qJD(3) + t203) * t159 + (-t194 * t232 + t139 * t68 + t262 + (t199 * qJD(1) + t250 * t298 + t35) * qJD(3)) * t156, -t274 * t298 + ((t139 * t298 - t76) * t233 + (t109 * t139 - t265) * qJD(3) + t222) * t159 + (t194 * t233 - t139 * t67 + t261 + (-t260 * qJD(1) + t139 * t251 - t36) * qJD(3)) * t156, t107 * t32 - t298 * t23 + t68 * t69 + (t34 * t237 + t10) * t159 + (t34 * t232 + t268 + (-qJD(1) * t46 - t29) * qJD(3)) * t156, -t107 * t19 + t109 * t23 - t45 * t68 - t46 * t67 + t185 * t235 + (-qJD(4) * t184 + t10 * t158 - t155 * t9) * t156, -t109 * t32 + t298 * t19 + t67 * t69 + (-t34 * t227 - t9) * t159 + (t34 * t233 - t266 + (qJD(1) * t45 + t30) * qJD(3)) * t156, t10 * t46 + t16 * t69 + t19 * t30 + t23 * t29 + t32 * t34 + t45 * t9, t12 * t90 + t28 * t51, -t12 * t89 - t13 * t90 - t181 * t28 - t27 * t51, -t271 + t11 + (-t51 - t255) * t236, -t272 + t21 + (t181 + t256) * t236 (t226 - t238) * t236 -(t14 * t157 - t15 * t154) * t226 + t215 * t159 + t22 * t181 + t52 * t13 + t8 * t89 + t24 * t27 + (-t159 * t2 + t186 * t226) * qJD(6) + (-qJD(1) * t187 - t1) * t236 (qJD(6) * t187 + t14 * t154 + t15 * t157) * t226 - t176 * t159 + t22 * t51 + t52 * t12 + t8 * t90 + t24 * t28 + (qJD(1) * t186 + t2) * t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244, -t243, 0, 0, 0, 0, 0, t163, -t167 * t298 - t201 + t92, t163, t196 + (t156 * t200 + t213) * t155 + t279, t103 + (-qJD(3) * t109 - t212) * t156 + t201 (qJD(3) * t184 - t16) * t159 + (qJD(3) * t34 + t164) * t156, 0, 0, 0, 0, 0, t272 + t21 + (-t181 + t256) * t236, t271 + t11 + (-t51 + t255) * t236; 0, 0, 0, 0, -t156 * t162 * t159, t240 * t162, 0, 0, 0, qJD(3) * t86 - t125 * t239 - t78, -t125 * t238, t109 * t251 - t155 * t67 (-t67 - t254) * t158 + (-t252 - t68) * t155, t298 * t232 + (-t298 * t245 + (-t109 + t237) * t156) * qJD(1), -t212 + (t298 * t247 + (t107 + t227) * t156) * qJD(1), -t298 * t239, -pkin(3) * t68 - t261 - t202 * t298 - t86 * t107 + (-t219 - t267) * qJD(4) + (-t35 * t156 + (t159 * t194 - t217) * t155) * qJD(1), pkin(3) * t67 + t262 + t275 * t298 - t86 * t109 + (t220 - t265) * qJD(4) + (t194 * t245 + (t36 - t216) * t156) * qJD(1), t119 * t68 + t298 * t42 - t266 - t280 * t107 + (t155 * t34 - t219) * qJD(4) + (t156 * t29 + (-t159 * t34 - t217) * t155) * qJD(1), t107 * t40 - t109 * t42 + (t9 + t298 * t29 + (qJD(4) * t109 - t68) * pkin(8)) * t158 + (pkin(8) * t200 + t291) * t155, t119 * t67 - t298 * t40 - t268 + t280 * t109 + (-t158 * t34 - t220) * qJD(4) + (t34 * t245 + (-t30 + t216) * t156) * qJD(1), t164 * pkin(8) + t119 * t16 - t280 * t34 - t29 * t42 - t30 * t40, -t12 * t180 + t278 * t51, -t111 * t12 + t13 * t180 - t181 * t278 - t277 * t51, -t278 * t226 + (qJD(3) * t180 + t51) * t239, t277 * t226 + (qJD(3) * t111 - t181) * t239, -t226 * t239, t101 * t13 + t8 * t111 + t276 * t181 + t277 * t24 - (t154 * t191 - t157 * t190) * t226 + (-(t126 * t157 - t127 * t154) * qJD(3) + t1) * t239, t101 * t12 - t8 * t180 + t276 * t51 + t278 * t24 - (t154 * t190 + t157 * t191) * t226 + ((t126 * t154 + t127 * t157) * qJD(3) - t2) * t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, -t107 ^ 2 + t284, t38, -t68 + t252, t141, t109 * t194 + t175, -t107 * t194 + t165, -t107 * t56 + t175 + 0.2e1 * t197 - t273, pkin(4) * t67 - qJ(5) * t68 + (t30 - t36) * t109 + (t29 - t241) * t107, -t107 * t34 + t109 * t56 + 0.2e1 * t128 + 0.2e1 * t134 - t165, -pkin(4) * t10 + qJ(5) * t9 + t241 * t30 - t29 * t36 - t34 * t56, -t292, -t296, t286, t226 * t51 + t13, t141, t182 * t141 - t39 * t181 - (t154 * t242 - t157 * t26) * t226 + (t183 * t226 + t2) * qJD(6) + t295, t183 * t141 - t39 * t51 - t293 - (t154 * t26 + t157 * t242) * t226 + (-t182 * t226 - t269) * qJD(6) + t223; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141 + t253, t38, -t298 ^ 2 - t284, t273 + t291, 0, 0, 0, 0, 0, -t109 * t181 - t141 * t157 - t154 * t198, -t109 * t51 + t141 * t154 - t157 * t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t292, t296, -t286, t288 * t51 - t204, -t141, t288 * t2 - t295, -t1 * t226 - t176 + t293;];
tauc_reg  = t3;
