% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRP4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% tau_reg [6x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:44:01
% EndTime: 2019-03-08 21:44:08
% DurationCPUTime: 3.57s
% Computational Cost: add. (2686->407), mult. (5981->540), div. (0->0), fcn. (4346->10), ass. (0->229)
t145 = sin(pkin(6));
t150 = sin(qJ(2));
t153 = cos(qJ(2));
t231 = qJD(1) * qJD(2);
t212 = t153 * t231;
t146 = cos(pkin(6));
t242 = qJD(3) * t146;
t267 = qJDD(2) * pkin(8);
t319 = qJD(1) * t242 + t267 + (qJDD(1) * t150 + t212) * t145;
t149 = sin(qJ(3));
t230 = qJD(2) * qJD(3);
t210 = t149 * t230;
t152 = cos(qJ(3));
t226 = t152 * qJDD(2);
t316 = -t210 + t226;
t294 = pkin(4) + pkin(8);
t289 = pkin(3) * t152;
t318 = pkin(2) + t289;
t248 = qJD(1) * t145;
t220 = t150 * t248;
t102 = qJD(2) * pkin(8) + t220;
t247 = qJD(1) * t152;
t56 = t149 * t102 - t146 * t247;
t317 = qJD(4) + t56;
t257 = t149 * qJ(4);
t315 = -t257 - t318;
t156 = qJD(2) ^ 2;
t186 = qJDD(2) * t153 - t150 * t156;
t208 = t153 * t230;
t260 = t145 * t153;
t218 = qJD(2) * t260;
t259 = t146 * t149;
t261 = t145 * t152;
t84 = t150 * t261 + t259;
t42 = qJD(3) * t84 + t149 * t218;
t262 = t145 * t150;
t222 = t149 * t262;
t83 = -t146 * t152 + t222;
t314 = -t42 * qJD(3) - t83 * qJDD(3) + t145 * (-t149 * t208 + t152 * t186);
t41 = -qJD(3) * t222 + (t218 + t242) * t152;
t313 = t41 * qJD(3) + t84 * qJDD(3) + t145 * (t149 * t186 + t152 * t208);
t154 = -pkin(3) - pkin(9);
t229 = qJDD(1) * t146;
t238 = qJD(3) * t152;
t303 = t102 * t238 + t319 * t149;
t173 = -t152 * t229 + t303;
t167 = qJDD(4) + t173;
t209 = t152 * t230;
t227 = t149 * qJDD(2);
t172 = t209 + t227;
t11 = t172 * pkin(4) + t154 * qJDD(3) + t167;
t151 = cos(qJ(5));
t10 = t151 * t11;
t148 = sin(qJ(5));
t245 = qJD(2) * t149;
t305 = pkin(4) * t245 + t317;
t31 = t154 * qJD(3) + t305;
t215 = t153 * t248;
t207 = -pkin(2) - t257;
t89 = t154 * t152 + t207;
t47 = qJD(2) * t89 - t215;
t16 = t148 * t31 + t151 * t47;
t272 = qJ(4) * t152;
t194 = pkin(9) * t149 - t272;
t233 = t149 * qJD(4);
t166 = qJD(3) * t194 - t233;
t213 = t150 * t231;
t193 = -qJDD(1) * t260 + t145 * t213;
t183 = pkin(3) * t210 + t193;
t19 = qJD(2) * t166 + qJDD(2) * t89 + t183;
t164 = -qJD(5) * t16 - t148 * t19 + t10;
t241 = qJD(3) * t148;
t243 = qJD(2) * t152;
t93 = t151 * t243 + t241;
t29 = qJD(5) * t93 - t151 * qJDD(3) + t316 * t148;
t92 = qJDD(5) + t172;
t214 = t148 * t243;
t239 = qJD(3) * t151;
t95 = -t214 + t239;
t1 = pkin(5) * t92 + qJ(6) * t29 - qJD(6) * t95 + t164;
t129 = qJD(5) + t245;
t274 = cos(pkin(10));
t205 = t145 * t274;
t204 = t274 * t150;
t144 = sin(pkin(10));
t263 = t144 * t153;
t78 = t146 * t204 + t263;
t37 = t78 * t149 + t152 * t205;
t203 = t274 * t153;
t264 = t144 * t150;
t80 = -t146 * t264 + t203;
t39 = -t144 * t261 + t149 * t80;
t176 = g(1) * t39 + g(2) * t37 + g(3) * t83;
t236 = qJD(5) * t151;
t225 = -t148 * t11 - t151 * t19 - t31 * t236;
t237 = qJD(5) * t148;
t181 = t237 * t47 + t225;
t30 = -qJD(5) * t214 + t148 * qJDD(3) + (qJD(3) * qJD(5) + t316) * t151;
t2 = -qJ(6) * t30 - qJD(6) * t93 - t181;
t15 = -t148 * t47 + t151 * t31;
t7 = -qJ(6) * t95 + t15;
t6 = pkin(5) * t129 + t7;
t8 = -qJ(6) * t93 + t16;
t312 = -(t129 * t6 - t2) * t148 + (t129 * t8 + t1) * t151 - t176;
t101 = t294 * t238;
t254 = t151 * t153;
t169 = t145 * (-t148 * t150 + t149 * t254);
t311 = -qJD(1) * t169 + t151 * t101;
t141 = qJD(3) * qJ(4);
t57 = qJD(1) * t259 + t152 * t102;
t50 = -t141 - t57;
t114 = t294 * t149;
t256 = t149 * t153;
t168 = (t148 * t256 + t150 * t151) * t145;
t240 = qJD(3) * t149;
t134 = pkin(3) * t240;
t59 = t134 + t166;
t310 = qJD(1) * t168 - t148 * t101 - t114 * t236 - t151 * t59 + t237 * t89;
t49 = -qJD(3) * pkin(3) + t317;
t280 = t129 * t93;
t309 = t29 - t280;
t279 = t129 * t95;
t308 = -t30 + t279;
t76 = t151 * t92;
t307 = -t129 * t237 + t76;
t43 = t148 * t260 + t151 * t83;
t77 = -t146 * t203 + t264;
t79 = t146 * t263 + t204;
t306 = -g(1) * (-t148 * t79 + t151 * t39) - g(2) * (-t148 * t77 + t151 * t37) - g(3) * t43;
t106 = t207 - t289;
t246 = qJD(2) * t106;
t58 = -t215 + t246;
t304 = t58 * t245 + qJDD(4);
t46 = pkin(4) * t243 + t57;
t32 = t141 + t46;
t300 = t129 * t32 + t154 * t92;
t197 = g(1) * t79 + g(2) * t77;
t155 = qJD(3) ^ 2;
t288 = pkin(8) * t155;
t299 = 0.2e1 * qJDD(2) * pkin(2) + t145 * (-g(3) * t153 + t213) - t193 + t197 - t288;
t170 = -g(3) * t260 + t197;
t228 = qJDD(2) * t106;
t177 = -qJ(4) * t238 - t233;
t23 = qJD(2) * t177 + t183 + t228;
t75 = t134 + t177;
t296 = qJD(2) * (-t75 + t220) + t170 - t228 - t23 - t288;
t295 = t95 ^ 2;
t293 = -t7 + t6;
t202 = qJ(6) * t152 - t89;
t291 = pkin(5) * t238 + t202 * t236 + (-qJ(6) * t240 - qJD(5) * t114 + qJD(6) * t152 - t59) * t148 + t311;
t235 = qJD(5) * t152;
t216 = t148 * t235;
t232 = t151 * qJD(6);
t290 = -t152 * t232 + (t149 * t239 + t216) * qJ(6) - t310;
t287 = g(3) * t150;
t131 = pkin(5) * t151 + pkin(4);
t286 = pkin(8) + t131;
t135 = pkin(3) * t245;
t70 = qJD(2) * t194 + t135;
t285 = t148 * t46 + t151 * t70;
t253 = qJ(6) - t154;
t258 = t148 * t149;
t34 = t151 * t46;
t284 = t237 * t253 - t232 + t148 * t70 - t34 - (pkin(5) * t152 - qJ(6) * t258) * qJD(2);
t105 = t253 * t151;
t283 = -qJ(6) * t151 * t245 - qJD(5) * t105 - t148 * qJD(6) - t285;
t282 = t148 * t114 + t151 * t89;
t281 = qJD(2) * pkin(2);
t278 = t148 * t92;
t277 = t151 * t29;
t276 = t151 * t95;
t273 = pkin(8) * qJDD(3);
t271 = qJD(3) * t57;
t270 = qJD(3) * t93;
t269 = qJD(3) * t95;
t266 = qJDD(3) * pkin(3);
t265 = t129 * t151;
t255 = t151 * t152;
t251 = qJDD(1) - g(3);
t115 = t294 * t152;
t142 = t149 ^ 2;
t143 = t152 ^ 2;
t250 = t142 - t143;
t249 = t142 + t143;
t244 = qJD(2) * t150;
t234 = qJD(5) * t154;
t224 = t315 * t77;
t223 = t315 * t79;
t221 = t149 * t156 * t152;
t219 = t145 * t244;
t206 = pkin(5) * t148 + qJ(4);
t200 = t102 * t240 - t149 * t229 - t319 * t152;
t196 = g(1) * t80 + g(2) * t78;
t191 = g(3) * (t145 * qJ(4) * t256 + pkin(8) * t262 + t318 * t260);
t187 = t129 * t148;
t139 = qJDD(3) * qJ(4);
t140 = qJD(3) * qJD(4);
t20 = -t139 - t140 + t200;
t147 = -qJ(6) - pkin(9);
t184 = pkin(5) * t258 - t147 * t152;
t44 = t145 * t254 - t83 * t148;
t179 = -t129 * t236 - t278;
t38 = -t149 * t205 + t152 * t78;
t40 = t144 * t145 * t149 + t152 * t80;
t175 = -g(1) * t40 - g(2) * t38 - g(3) * t84;
t12 = t316 * pkin(4) - t20;
t171 = t12 + t175;
t103 = -t215 - t281;
t163 = -t273 + (t103 + t215 - t281) * qJD(3);
t162 = t273 + (-t215 - t58 - t246) * qJD(3);
t161 = qJD(3) * t56 + t175 - t200;
t160 = -t173 + t176;
t3 = pkin(5) * t30 + qJDD(6) + t12;
t21 = t167 - t266;
t157 = t21 * t149 - t20 * t152 + (t149 * t50 + t152 * t49) * qJD(3) - t196;
t104 = t253 * t148;
t100 = t294 * t240;
t99 = -qJ(4) * t243 + t135;
t98 = t151 * t114;
t91 = t93 ^ 2;
t74 = t83 * pkin(3);
t36 = t39 * pkin(3);
t35 = t37 * pkin(3);
t28 = -qJ(6) * t255 + t282;
t25 = t149 * pkin(5) + t148 * t202 + t98;
t24 = pkin(5) * t93 + qJD(6) + t32;
t14 = qJD(5) * t43 + t42 * t148 + t151 * t219;
t13 = qJD(5) * t44 - t148 * t219 + t42 * t151;
t4 = [t251, 0, t186 * t145 (-qJDD(2) * t150 - t153 * t156) * t145, 0, 0, 0, 0, 0, t314, -t313 (t149 * t83 + t152 * t84) * qJDD(2) + (t149 * t42 + t152 * t41 + (-t149 * t84 + t152 * t83) * qJD(3)) * qJD(2), -t314, t313, -t20 * t84 + t21 * t83 - t41 * t50 + t42 * t49 - g(3) + (-t153 * t23 + t244 * t58) * t145, 0, 0, 0, 0, 0, t129 * t13 + t30 * t84 + t41 * t93 + t43 * t92, -t129 * t14 - t29 * t84 + t41 * t95 + t44 * t92, -t13 * t95 - t14 * t93 + t29 * t43 + t30 * t44, t1 * t43 + t13 * t6 + t14 * t8 - t2 * t44 + t24 * t41 + t3 * t84 - g(3); 0, qJDD(2), t251 * t260 + t197, -t251 * t262 + t196, qJDD(2) * t142 + 0.2e1 * t149 * t209, 0.2e1 * t149 * t226 - 0.2e1 * t230 * t250, qJDD(3) * t149 + t152 * t155, qJDD(3) * t152 - t149 * t155, 0, t163 * t149 + t152 * t299, -t149 * t299 + t163 * t152, t249 * t267 + (-t212 * t249 - t287) * t145 + t157, t162 * t149 - t152 * t296, t149 * t296 + t162 * t152, t23 * t106 + t58 * t75 - g(1) * t223 - g(2) * t224 - t191 + (-t150 * t58 + (-t149 * t49 + t152 * t50) * t153) * t248 + t157 * pkin(8), -t235 * t276 + (t152 * t29 + t240 * t95) * t148 (-t148 * t93 + t276) * t240 + (t148 * t30 + t277 + (t148 * t95 + t151 * t93) * qJD(5)) * t152 (t129 * t241 - t29) * t149 + (t179 + t269) * t152 (t129 * t239 - t30) * t149 + (-t270 - t307) * t152, t129 * t238 + t149 * t92 (-t148 * t89 + t98) * t92 - t100 * t93 + t115 * t30 - t196 * t151 + (-t32 * t239 + t10 + (-t19 + t197) * t148) * t149 - g(3) * t168 + (-t148 * t59 + t311) * t129 + (-t129 * t282 - t149 * t16) * qJD(5) + (qJD(3) * t15 + t12 * t151 - t215 * t93 - t237 * t32) * t152, -t282 * t92 - t100 * t95 - t115 * t29 + t196 * t148 + (t197 * t151 + (qJD(3) * t32 + qJD(5) * t47) * t148 + t225) * t149 - g(3) * t169 + t310 * t129 + (-qJD(3) * t16 - t12 * t148 - t215 * t95 - t236 * t32) * t152, t25 * t29 - t28 * t30 - t291 * t95 - t290 * t93 + (-t148 * t6 + t151 * t8) * t240 + (t1 * t148 - t151 * t2 + (t148 * t8 + t151 * t6) * qJD(5) + t170) * t152, t2 * t28 + t1 * t25 + t3 * (pkin(5) * t255 + t115) - t24 * pkin(5) * t216 - g(1) * (-t184 * t79 + t286 * t80 + t223) - g(2) * (-t184 * t77 + t286 * t78 + t224) - t191 + t290 * t8 + t291 * t6 - t24 * t286 * t240 + (-t131 * t287 + (-g(3) * t184 - t24 * t247) * t153) * t145; 0, 0, 0, 0, -t221, t250 * t156, t227, t226, qJDD(3), -t103 * t245 + t160 + t271, -t103 * t243 - t161 (-pkin(3) * t149 + t272) * qJDD(2), -0.2e1 * t266 - t271 + (-qJD(2) * t99 - t229) * t152 - t176 + t303 + t304, 0.2e1 * t139 + 0.2e1 * t140 + (t149 * t99 + t152 * t58) * qJD(2) + t161, -t20 * qJ(4) - t21 * pkin(3) - t58 * t99 - t49 * t57 - g(1) * (qJ(4) * t40 - t36) - g(2) * (qJ(4) * t38 - t35) - g(3) * (qJ(4) * t84 - t74) - t317 * t50, -t187 * t95 - t277 (-t30 - t279) * t151 + (t29 + t280) * t148 (-t129 * t258 - t152 * t95) * qJD(2) + t307 (-t149 * t265 + t152 * t93) * qJD(2) + t179, -t129 * t243, -t15 * t243 + qJ(4) * t30 - t34 * t129 + t305 * t93 + t300 * t151 + ((t70 - t234) * t129 + t171) * t148, -qJ(4) * t29 + t285 * t129 + t16 * t243 + t305 * t95 - t300 * t148 + (-t129 * t234 + t171) * t151, t104 * t30 - t105 * t29 - t283 * t93 - t284 * t95 - t312, -t2 * t104 - t1 * t105 + t3 * t206 - g(1) * (t39 * t147 + t206 * t40 - t36) - g(2) * (t37 * t147 + t206 * t38 - t35) - g(3) * (t83 * t147 + t206 * t84 - t74) + t283 * t8 + t284 * t6 + (pkin(5) * t265 + t305) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, qJDD(3) + t221, -t142 * t156 - t155, qJD(3) * t50 - t160 - t266 + t304, 0, 0, 0, 0, 0, -t129 * t187 - t270 + t76, -t129 ^ 2 * t151 - t269 - t278, t148 * t308 + t151 * t309, -qJD(3) * t24 + t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95 * t93, -t91 + t295, -t309, t308, t92, t16 * t129 - t32 * t95 + t164 + t306, t15 * t129 + t32 * t93 - g(1) * (-t148 * t39 - t151 * t79) - g(2) * (-t148 * t37 - t151 * t77) - g(3) * t44 + t181, pkin(5) * t29 - t293 * t93, t293 * t8 + (-t24 * t95 + t1 + t306) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91 - t295, t6 * t95 + t8 * t93 + t175 + t3;];
tau_reg  = t4;
