% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PPRRRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:01:33
% EndTime: 2019-03-08 19:01:43
% DurationCPUTime: 4.38s
% Computational Cost: add. (4111->374), mult. (10168->571), div. (0->0), fcn. (9436->18), ass. (0->218)
t147 = qJD(4) + qJD(5);
t157 = sin(qJ(5));
t161 = cos(qJ(4));
t283 = cos(qJ(5));
t229 = qJD(3) * t283;
t158 = sin(qJ(4));
t252 = qJD(3) * t158;
t296 = -t157 * t252 + t161 * t229;
t297 = t296 * t147;
t103 = qJD(6) - t296;
t295 = t103 - qJD(6);
t287 = pkin(9) + pkin(10);
t152 = sin(pkin(7));
t153 = sin(pkin(6));
t154 = cos(pkin(13));
t155 = cos(pkin(7));
t265 = cos(pkin(6));
t100 = -t153 * t154 * t152 + t265 * t155;
t150 = qJ(4) + qJ(5);
t143 = sin(t150);
t144 = cos(t150);
t151 = sin(pkin(13));
t264 = cos(pkin(12));
t204 = t265 * t264;
t263 = sin(pkin(12));
t101 = t151 * t204 + t263 * t154;
t159 = sin(qJ(3));
t162 = cos(qJ(3));
t173 = t263 * t151 - t154 * t204;
t220 = t153 * t264;
t291 = t152 * t220 + t173 * t155;
t58 = t101 * t162 - t159 * t291;
t203 = t265 * t263;
t102 = -t151 * t203 + t264 * t154;
t174 = t264 * t151 + t154 * t203;
t219 = t153 * t263;
t290 = -t152 * t219 + t174 * t155;
t60 = t102 * t162 - t159 * t290;
t256 = t154 * t155;
t191 = t151 * t162 + t159 * t256;
t180 = t191 * t153;
t221 = t152 * t265;
t82 = t159 * t221 + t180;
t83 = t173 * t152 - t155 * t220;
t84 = t174 * t152 + t155 * t219;
t187 = -g(3) * (t100 * t144 - t82 * t143) - g(2) * (-t58 * t143 + t83 * t144) - g(1) * (-t60 * t143 + t84 * t144);
t145 = qJDD(4) + qJDD(5);
t126 = t265 * qJDD(1) + qJDD(2);
t133 = t265 * qJD(1) + qJD(2);
t238 = t162 * t256;
t259 = t151 * t159;
t192 = t238 - t259;
t250 = qJD(3) * t162;
t42 = qJDD(3) * pkin(9) + (t126 * t159 + t133 * t250) * t152 + (t192 * qJD(3) * qJD(1) + t191 * qJDD(1)) * t153;
t218 = pkin(10) * qJDD(3) + t42;
t258 = t152 * t159;
t75 = qJD(1) * t180 + t133 * t258;
t224 = t287 * qJD(3) + t75;
t253 = qJD(1) * t153;
t235 = t154 * t253;
t99 = t155 * t133 - t152 * t235;
t48 = t158 * t99 + t224 * t161;
t244 = qJDD(1) * t153;
t225 = t154 * t244;
t98 = t155 * t126 - t152 * t225;
t89 = t161 * t98;
t11 = qJDD(4) * pkin(4) - qJD(4) * t48 - t218 * t158 + t89;
t47 = -t224 * t158 + t161 * t99;
t13 = qJD(4) * t47 + t158 * t98 + t218 * t161;
t240 = t283 * t48;
t46 = qJD(4) * pkin(4) + t47;
t20 = t157 * t46 + t240;
t289 = t20 * qJD(5) - t283 * t11 + t157 * t13;
t3 = -t145 * pkin(5) + t289;
t183 = t187 - t3;
t213 = t155 * t235;
t251 = qJD(3) * t159;
t234 = t152 * t251;
t236 = t151 * t253;
t239 = t153 * t259;
t169 = -(t126 * t152 + t155 * t225) * t162 + qJDD(1) * t239 + t133 * t234 + t213 * t251 + t236 * t250;
t57 = t101 * t159 + t162 * t291;
t59 = t102 * t159 + t162 * t290;
t208 = t162 * t221;
t81 = -t153 * t238 - t208 + t239;
t186 = g(1) * t59 + g(2) * t57 + g(3) * t81;
t294 = t75 * qJD(3) - t169 + t186;
t178 = t186 * t144;
t222 = qJDD(3) * t283;
t243 = t158 * qJDD(3);
t202 = t157 * t243 - t161 * t222;
t255 = t157 * t161;
t112 = t283 * t158 + t255;
t86 = t147 * t112;
t66 = t86 * qJD(3) + t202;
t64 = qJDD(6) + t66;
t141 = -t161 * pkin(4) - pkin(3);
t188 = -t157 * t158 + t283 * t161;
t80 = -pkin(5) * t188 - t112 * pkin(11) + t141;
t293 = t80 * t64 + t178;
t249 = qJD(4) * t158;
t212 = pkin(4) * t249 - t75;
t156 = sin(qJ(6));
t61 = t100 * t161 - t82 * t158;
t62 = t100 * t158 + t82 * t161;
t26 = t157 * t61 + t283 * t62;
t160 = cos(qJ(6));
t78 = t81 * t160;
t292 = -t156 * t26 + t78;
t275 = t157 * t48;
t19 = t283 * t46 - t275;
t17 = -t147 * pkin(5) - t19;
t210 = g(1) * t60 + g(2) * t58;
t185 = g(3) * t82 + t210;
t228 = qJD(5) * t283;
t248 = qJD(5) * t157;
t172 = t157 * t11 + t283 * t13 + t46 * t228 - t48 * t248;
t2 = t145 * pkin(11) + t172;
t237 = qJD(4) * t287;
t116 = t158 * t237;
t117 = t161 * t237;
t121 = t287 * t158;
t122 = t287 * t161;
t189 = -t283 * t121 - t157 * t122;
t74 = -t159 * t236 + (t133 * t152 + t213) * t162;
t281 = -t189 * qJD(5) + t283 * t116 + t157 * t117 + t188 * t74;
t108 = -qJD(3) * t255 - t158 * t229;
t68 = t141 * qJD(3) - t74;
t45 = -pkin(5) * t296 + t108 * pkin(11) + t68;
t85 = t147 * t188;
t95 = -t157 * t121 + t283 * t122;
t288 = (qJD(6) * t45 + t2) * t188 + t17 * t85 + t3 * t112 + (-qJD(6) * t80 + t281) * t103 - t95 * t64 - t185;
t198 = t160 * t108 - t156 * t147;
t242 = t161 * qJDD(3);
t65 = t157 * t242 + t158 * t222 + t297;
t39 = -t198 * qJD(6) - t160 * t145 + t156 * t65;
t280 = t86 * pkin(5) - t85 * pkin(11) + t212;
t279 = t95 * qJD(5) - t112 * t74 - t157 * t116 + t283 * t117;
t278 = qJD(3) * pkin(3);
t276 = t156 * t64;
t273 = t160 * t64;
t272 = t160 * t198;
t271 = t17 * t296;
t270 = t17 * t112;
t246 = qJD(6) * t160;
t247 = qJD(6) * t156;
t38 = t108 * t247 + t156 * t145 + t147 * t246 + t160 * t65;
t269 = t38 * t156;
t268 = t81 * t156;
t90 = -t156 * t108 - t160 * t147;
t267 = t90 * t103;
t266 = t198 * t103;
t261 = t103 * t108;
t260 = t108 * t296;
t257 = t152 * t162;
t148 = t158 ^ 2;
t254 = -t161 ^ 2 + t148;
t245 = qJD(3) * qJD(4);
t233 = t152 * t250;
t18 = t147 * pkin(11) + t20;
t200 = t156 * t18 - t160 * t45;
t231 = -t108 * t200 + t17 * t247;
t227 = t158 * t245;
t226 = t162 * t245;
t139 = t157 * pkin(4) + pkin(11);
t79 = -t108 * pkin(5) - pkin(11) * t296;
t216 = pkin(4) * t252 + qJD(6) * t139 + t79;
t215 = t103 * t160;
t21 = t157 * t47 + t240;
t211 = pkin(4) * t248 - t21;
t209 = -g(1) * t84 - g(2) * t83;
t201 = -t139 * t64 - t271;
t8 = t156 * t45 + t160 * t18;
t199 = t160 * t26 + t268;
t22 = t283 * t47 - t275;
t197 = -pkin(4) * t228 + t22;
t164 = qJD(3) ^ 2;
t196 = qJDD(3) * t162 - t159 * t164;
t195 = -t157 * t62 + t283 * t61;
t104 = t161 * t155 - t158 * t258;
t105 = t158 * t155 + t161 * t258;
t73 = t157 * t104 + t283 * t105;
t194 = -t156 * t73 - t160 * t257;
t193 = t156 * t257 - t160 * t73;
t190 = t283 * t104 - t157 * t105;
t184 = -t8 * t108 - t183 * t156 + t17 * t246;
t69 = -t74 - t278;
t179 = -t69 * qJD(3) + t210 - t42;
t177 = -pkin(9) * qJDD(4) + (t69 + t74 - t278) * qJD(4);
t163 = qJD(4) ^ 2;
t168 = 0.2e1 * qJDD(3) * pkin(3) - pkin(9) * t163 + t294;
t37 = pkin(4) * t227 + t141 * qJDD(3) + t169;
t166 = t68 * t108 + t187 - t289;
t32 = t83 * t143 + t58 * t144;
t34 = t84 * t143 + t60 * t144;
t54 = t100 * t143 + t82 * t144;
t165 = g(1) * t34 + g(2) * t32 + g(3) * t54 - t296 * t68 - t172;
t140 = -t283 * pkin(4) - pkin(5);
t88 = -t105 * qJD(4) - t158 * t233;
t87 = t104 * qJD(4) + t161 * t233;
t77 = t82 * qJD(3);
t76 = (t192 * t153 + t208) * qJD(3);
t67 = t108 ^ 2 - t296 ^ 2;
t50 = -t202 + (-qJD(3) * t112 - t108) * t147;
t49 = t65 - t297;
t36 = t61 * qJD(4) + t76 * t161;
t35 = -t62 * qJD(4) - t76 * t158;
t28 = t73 * qJD(5) + t157 * t87 - t283 * t88;
t27 = t190 * qJD(5) + t157 * t88 + t283 * t87;
t24 = t103 * t215 - t108 * t198 + t276;
t23 = -t103 ^ 2 * t156 - t90 * t108 + t273;
t16 = -t198 * t215 + t269;
t12 = t66 * pkin(5) - t65 * pkin(11) + t37;
t10 = t160 * t12;
t6 = t26 * qJD(5) + t157 * t36 - t283 * t35;
t5 = t195 * qJD(5) + t157 * t35 + t283 * t36;
t4 = (t38 - t267) * t160 + (-t39 + t266) * t156;
t1 = [qJDD(1) - g(3), t126 * t265 - g(3) + (t151 ^ 2 + t154 ^ 2) * t153 ^ 2 * qJDD(1), 0, -t77 * qJD(3) - t81 * qJDD(3), -t76 * qJD(3) - t82 * qJDD(3), 0, 0, 0, 0, 0, -t81 * t242 + t35 * qJD(4) + t61 * qJDD(4) + (-t161 * t77 + t249 * t81) * qJD(3), t81 * t243 - t36 * qJD(4) - t62 * qJDD(4) + (qJD(4) * t161 * t81 + t158 * t77) * qJD(3), 0, 0, 0, 0, 0, t145 * t195 - t6 * t147 - t296 * t77 + t81 * t66, -t77 * t108 - t26 * t145 - t5 * t147 + t81 * t65, 0, 0, 0, 0, 0 (-qJD(6) * t199 - t156 * t5 + t77 * t160) * t103 + t292 * t64 + t6 * t90 - t195 * t39 -(qJD(6) * t292 + t77 * t156 + t160 * t5) * t103 - t199 * t64 - t6 * t198 - t195 * t38; 0, -g(3) * t265 + (-g(1) * t263 + t264 * g(2)) * t153 + t126, 0, t196 * t152 (-qJDD(3) * t159 - t162 * t164) * t152, 0, 0, 0, 0, 0, t88 * qJD(4) + t104 * qJDD(4) + (-t158 * t226 + t161 * t196) * t152, -t87 * qJD(4) - t105 * qJDD(4) + (-t158 * t196 - t161 * t226) * t152, 0, 0, 0, 0, 0, t190 * t145 - t28 * t147 + (-t162 * t66 - t251 * t296) * t152, -t73 * t145 - t27 * t147 + (-t108 * t251 - t162 * t65) * t152, 0, 0, 0, 0, 0 (qJD(6) * t193 - t156 * t27 + t160 * t234) * t103 + t194 * t64 + t28 * t90 - t190 * t39 -(qJD(6) * t194 + t156 * t234 + t160 * t27) * t103 + t193 * t64 - t28 * t198 - t190 * t38; 0, 0, qJDD(3), t294, -t126 * t258 - t191 * t244 + (-t133 * t257 - t192 * t253 + t74) * qJD(3) + t185, t148 * qJDD(3) + 0.2e1 * t161 * t227, 0.2e1 * t158 * t242 - 0.2e1 * t245 * t254, qJDD(4) * t158 + t163 * t161, qJDD(4) * t161 - t163 * t158, 0, t158 * t177 + t161 * t168, -t158 * t168 + t161 * t177, -t108 * t85 + t65 * t112, t108 * t86 - t112 * t66 + t188 * t65 + t296 * t85, t112 * t145 + t85 * t147, t145 * t188 - t86 * t147, 0, t141 * t66 + t145 * t189 - t147 * t279 - t188 * t37 - t212 * t296 + t68 * t86 + t178, -t108 * t212 + t37 * t112 + t141 * t65 - t143 * t186 - t95 * t145 + t147 * t281 + t68 * t85, -t85 * t272 + (t38 * t160 + t198 * t247) * t112 (t156 * t198 - t160 * t90) * t85 + (-t269 - t160 * t39 + (t156 * t90 + t272) * qJD(6)) * t112, t112 * t273 - t38 * t188 - t198 * t86 + (-t112 * t247 + t160 * t85) * t103, -t112 * t276 + t39 * t188 - t90 * t86 + (-t112 * t246 - t156 * t85) * t103, t103 * t86 - t188 * t64, -t10 * t188 - t189 * t39 - t200 * t86 + t279 * t90 + (t280 * t103 + (-t95 * t103 + t18 * t188 + t270) * qJD(6) + t293) * t160 + t288 * t156, -t189 * t38 - t8 * t86 - t279 * t198 + ((-qJD(6) * t18 + t12) * t188 - qJD(6) * t270 + (qJD(6) * t95 - t280) * t103 - t293) * t156 + t288 * t160; 0, 0, 0, 0, 0, -t158 * t164 * t161, t254 * t164, t243, t242, qJDD(4), -g(3) * t61 + t158 * t179 + t161 * t209 + t89, g(3) * t62 + (-t209 - t98) * t158 + t179 * t161, t260, t67, t49, t50, t145, t21 * t147 + (t145 * t283 - t147 * t248 + t252 * t296) * pkin(4) + t166, t22 * t147 + (t108 * t252 - t145 * t157 - t147 * t228) * pkin(4) + t165, t16, t4, t24, t23, t261, t140 * t39 + t211 * t90 + (t103 * t197 + t201) * t156 + (-t103 * t216 + t183) * t160 + t231, t140 * t38 - t211 * t198 + t201 * t160 + (t156 * t216 + t160 * t197) * t103 + t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t260, t67, t49, t50, t145, t20 * t147 + t166, t19 * t147 + t165, t16, t4, t24, t23, t261, -pkin(5) * t39 - t20 * t90 + (-pkin(11) * t64 + t19 * t103 - t271) * t156 + ((-pkin(11) * qJD(6) - t79) * t103 + t183) * t160 + t231, -pkin(5) * t38 + (t156 * t79 + t160 * t19) * t103 + t20 * t198 - t160 * t271 + (t103 * t247 - t273) * pkin(11) + t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198 * t90, t198 ^ 2 - t90 ^ 2, t38 + t267, -t39 - t266, t64, -t156 * t2 + t10 + t17 * t198 - g(1) * (-t34 * t156 + t59 * t160) - g(2) * (-t32 * t156 + t57 * t160) - g(3) * (-t54 * t156 + t78) + t295 * t8, -t160 * t2 - t156 * t12 + t17 * t90 - g(1) * (-t59 * t156 - t34 * t160) - g(2) * (-t57 * t156 - t32 * t160) - g(3) * (-t54 * t160 - t268) - t295 * t200;];
tau_reg  = t1;
