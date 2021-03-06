% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRP12_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_inertiaDJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:50:39
% EndTime: 2019-03-09 06:51:00
% DurationCPUTime: 8.68s
% Computational Cost: add. (15529->486), mult. (45352->842), div. (0->0), fcn. (48089->12), ass. (0->237)
t150 = sin(pkin(7));
t152 = cos(pkin(12));
t291 = cos(pkin(6));
t261 = pkin(1) * t291;
t149 = sin(pkin(12));
t151 = sin(pkin(6));
t288 = t149 * t151;
t290 = cos(pkin(7));
t168 = t152 * t261 + t291 * pkin(2) + (-t290 * pkin(9) - qJ(2)) * t288;
t187 = t151 * (-pkin(9) * t149 * t150 - pkin(2) * t152 - pkin(1));
t326 = t150 * t187 + t168 * t290;
t241 = t291 * t150;
t243 = t152 * t290;
t325 = t151 * t243 + t241;
t153 = sin(qJ(5));
t155 = cos(qJ(5));
t154 = sin(qJ(4));
t156 = cos(qJ(4));
t308 = sin(qJ(3));
t113 = t156 * t150 * t308 + t154 * t290;
t309 = cos(qJ(3));
t260 = t150 * t309;
t93 = t153 * t113 + t155 * t260;
t237 = t153 * t260;
t94 = t155 * t113 - t237;
t215 = -t153 * t94 + t155 * t93;
t242 = t156 * t290;
t259 = t154 * t308;
t112 = t150 * t259 - t242;
t246 = qJD(3) * t309;
t231 = t156 * t246;
t167 = -qJD(4) * t112 + t150 * t231;
t245 = qJD(3) * t308;
t232 = t150 * t245;
t64 = t93 * qJD(5) - t153 * t232 - t155 * t167;
t142 = qJD(5) * t155;
t65 = -qJD(5) * t237 + t113 * t142 + t153 * t167 - t155 * t232;
t178 = qJD(5) * t215 + t65 * t153 - t64 * t155;
t143 = qJD(4) * t154;
t276 = qJD(4) * t156;
t233 = t150 * t246;
t279 = qJD(4) * t113;
t92 = t154 * t233 + t279;
t196 = t112 * t276 + t92 * t154;
t324 = -t113 * t143 + t167 * t156 + t196;
t287 = t151 * t152;
t191 = t150 * t287 - t291 * t290;
t89 = t308 * t241 + (t309 * t149 + t308 * t243) * t151;
t70 = -t154 * t191 + t89 * t156;
t88 = t308 * t288 - t325 * t309;
t59 = t153 * t70 - t88 * t155;
t60 = t153 * t88 + t155 * t70;
t217 = t153 * t60 + t155 * t59;
t275 = qJD(5) * t153;
t83 = t88 * qJD(3);
t58 = -t89 * t143 + (-qJD(4) * t191 - t83) * t156;
t84 = t89 * qJD(3);
t33 = t84 * t153 - t70 * t275 + (qJD(5) * t88 + t58) * t155;
t296 = t33 * t153;
t32 = qJD(5) * t60 + t58 * t153 - t84 * t155;
t297 = t32 * t155;
t298 = t155 * t60;
t299 = t153 * t59;
t323 = t154 * (qJD(5) * (-t298 + t299) - t296 - t297) - t217 * t276;
t322 = 0.2e1 * t191;
t57 = qJD(4) * t70 - t83 * t154;
t69 = t154 * t89 + t156 * t191;
t36 = -t155 * t57 + t69 * t275;
t321 = pkin(11) * t36;
t145 = t153 ^ 2;
t147 = t155 ^ 2;
t283 = t145 - t147;
t130 = t283 * qJD(5);
t165 = t325 * pkin(9) + qJ(2) * t287 + t149 * t261;
t55 = t309 * t165 + t326 * t308;
t318 = t308 * t165 - t326 * t309;
t273 = qJD(5) * t156;
t254 = t153 * t273;
t277 = qJD(4) * t155;
t115 = t154 * t277 + t254;
t302 = t154 * pkin(11);
t228 = -pkin(4) * t156 - t302;
t209 = -pkin(3) + t228;
t189 = t155 * t209;
t304 = pkin(11) * t156;
t227 = pkin(4) * t154 - t304;
t195 = qJD(4) * t227;
t75 = pkin(10) * t115 - qJD(5) * t189 - t153 * t195;
t221 = pkin(5) * t153 - qJ(6) * t155;
t206 = pkin(10) + t221;
t110 = t206 * t154;
t111 = qJD(5) * t221 - t153 * qJD(6);
t222 = pkin(5) * t155 + qJ(6) * t153;
t126 = -pkin(4) - t222;
t317 = qJD(4) * (-t126 * t156 + t302) - qJD(5) * t110 - t111 * t154;
t295 = t33 * t155;
t315 = qJD(5) * t217 + t153 * t32 - t295;
t314 = t222 * qJD(5) - t155 * qJD(6);
t249 = t153 * t143;
t305 = pkin(10) * t156;
t140 = t155 * t305;
t284 = t153 * t209 + t140;
t76 = pkin(10) * t249 - qJD(5) * t284 + t155 * t195;
t278 = qJD(4) * t153;
t37 = t69 * t142 + t153 * t57;
t313 = t154 * (qJD(4) * t59 + t37) + t156 * (t69 * t278 - t32);
t312 = 0.2e1 * qJD(6);
t311 = pkin(11) * t33;
t310 = t57 * pkin(5);
t307 = pkin(4) * t155;
t306 = pkin(10) * t153;
t289 = qJD(6) * t69;
t281 = qJD(2) * t151;
t258 = t149 * t281;
t207 = t290 * t258;
t257 = t152 * t281;
t43 = qJD(3) * t55 + t309 * t207 + t308 * t257;
t158 = t57 * pkin(4) - t58 * pkin(11) + t43;
t52 = pkin(3) * t191 + t318;
t160 = t69 * pkin(4) - t70 * pkin(11) + t52;
t159 = qJD(5) * t160;
t236 = t150 * t258;
t182 = pkin(3) * t84 + pkin(10) * t83 + t236;
t42 = qJD(3) * t318 + t308 * t207 - t309 * t257;
t66 = -t150 * t168 + t290 * t187;
t47 = t88 * pkin(3) - t89 * pkin(10) + t66;
t53 = -pkin(10) * t191 + t55;
t14 = t53 * t143 - t154 * t182 + t156 * t42 - t47 * t276;
t186 = t84 * pkin(11) - t14;
t27 = t154 * t47 + t156 * t53;
t24 = pkin(11) * t88 + t27;
t3 = -t153 * t158 + t24 * t275 + (-t159 - t186) * t155;
t300 = qJ(6) * t57;
t1 = t289 - t3 + t300;
t303 = t1 * t155;
t301 = t3 * t155;
t10 = t153 * t160 + t155 * t24;
t77 = t112 * t92;
t294 = t57 * t156;
t293 = t58 * t154;
t146 = t154 ^ 2;
t282 = -t156 ^ 2 + t146;
t274 = qJD(5) * t154;
t271 = t156 * qJD(6);
t270 = 0.2e1 * t59 * t32;
t39 = 0.2e1 * t69 * t57;
t269 = 0.2e1 * t88 * t84;
t268 = -0.2e1 * pkin(3) * qJD(4);
t267 = -0.2e1 * pkin(4) * qJD(5);
t266 = t153 * t305;
t265 = pkin(5) * t143;
t264 = pkin(10) * t276;
t263 = pkin(11) * t275;
t262 = pkin(11) * t142;
t256 = qJ(6) * t143;
t255 = t155 * t276;
t253 = t154 * t142;
t252 = t155 * t273;
t248 = t153 * t142;
t247 = t154 * t276;
t244 = qJD(4) * t309;
t26 = -t154 * t53 + t156 * t47;
t240 = t282 * qJD(4);
t239 = 0.2e1 * t247;
t235 = t153 * t255;
t234 = t146 * t248;
t6 = qJ(6) * t69 + t10;
t9 = -t153 * t24 + t155 * t160;
t7 = -t69 * pkin(5) - t9;
t226 = -t153 * t6 + t155 * t7;
t225 = t32 * t60 + t33 * t59;
t224 = t32 * t69 + t57 * t59;
t223 = -t10 * t153 - t155 * t9;
t214 = -0.2e1 * t291 * t281;
t100 = -qJ(6) * t156 + t284;
t102 = -t155 * (-pkin(3) - t302) + (pkin(5) + t306 + t307) * t156;
t213 = -t100 * t153 + t102 * t155;
t107 = t189 - t266;
t212 = -t107 * t155 - t153 * t284;
t23 = -t88 * pkin(4) - t26;
t16 = t59 * pkin(5) - t60 * qJ(6) + t23;
t15 = -t47 * t143 + t154 * t42 + t156 * t182 - t53 * t276;
t12 = -t84 * pkin(4) - t15;
t5 = t32 * pkin(5) - t33 * qJ(6) - t60 * qJD(6) + t12;
t205 = -t16 * t142 - t5 * t153;
t204 = -t5 * t155 + t16 * t275;
t38 = t69 * t143 - t294;
t203 = t154 * t84 + t88 * t276;
t202 = t88 * t143 - t156 * t84;
t201 = t12 * t153 + t23 * t142;
t200 = -t12 * t155 + t23 * t275;
t199 = t59 * t275 - t297;
t198 = -0.2e1 * t94 * t64 + 0.2e1 * t65 * t93 + 0.2e1 * t77;
t197 = (t298 + t299) * pkin(11);
t72 = t112 * t275 - t155 * t92;
t194 = -t94 * t32 + t33 * t93 + t64 * t59 + t60 * t65;
t193 = t112 * t32 - t57 * t93 + t92 * t59 - t65 * t69;
t192 = 0.2e1 * (t149 ^ 2 + t152 ^ 2) * t151 ^ 2 * qJD(2);
t188 = t37 * pkin(11);
t183 = t112 * t253 - t93 * t143 + t196 * t153 + t156 * t65;
t181 = t112 * t33 - t57 * t94 + t60 * t92 + t64 * t69;
t175 = -t153 * t186 + t155 * t158;
t173 = -t14 * t156 - t15 * t154 + (-t154 * t27 - t156 * t26) * qJD(4);
t172 = t178 * pkin(11);
t68 = -t75 + t256 - t271;
t73 = -t265 - t76;
t171 = qJD(5) * t213 + t73 * t153 + t68 * t155;
t170 = qJD(5) * t212 - t76 * t153 - t75 * t155;
t169 = t59 * t253 + (t154 * t32 + t59 * t276) * t153;
t161 = t215 * t276 + (t153 * t64 + t155 * t65 + (-t153 * t93 - t155 * t94) * qJD(5)) * t154;
t157 = -t24 * t142 - t153 * t159 + t175;
t138 = -0.2e1 * t247;
t137 = -0.2e1 * t248;
t136 = 0.2e1 * t248;
t135 = pkin(11) * t252;
t117 = t249 - t252;
t116 = t153 * t276 + t253;
t114 = -t153 * t274 + t255;
t106 = 0.2e1 * t147 * t247 - 0.2e1 * t234;
t105 = 0.2e1 * t145 * t247 + 0.2e1 * t234;
t101 = t274 * t283 - t235;
t99 = -t153 * t240 + t154 * t252;
t98 = 0.4e1 * t154 * t248 + t283 * t276;
t96 = 0.2e1 * t154 * t254 + 0.2e1 * t282 * t277;
t95 = t130 * t146 - 0.2e1 * t154 * t235;
t74 = t154 * t314 + t206 * t276;
t71 = t112 * t142 + t153 * t92;
t34 = (t112 * t277 - t64) * t156 + (-qJD(4) * t94 - t72) * t154;
t31 = pkin(11) * t297;
t21 = 0.2e1 * t60 * t33;
t18 = t60 * t142 + t296;
t17 = t60 * t255 + (-t60 * t275 + t295) * t154;
t13 = 0.2e1 * t33 * t69 + 0.2e1 * t57 * t60;
t8 = (t69 * t277 - t33) * t156 + (qJD(4) * t60 - t36) * t154;
t4 = -t10 * qJD(5) + t175;
t2 = -t157 - t310;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149 * t214, t152 * t214, t192, qJ(2) * t192, -0.2e1 * t89 * t83, 0.2e1 * t83 * t88 - 0.2e1 * t84 * t89, t83 * t322, t269, t84 * t322, 0, 0.2e1 * t191 * t43 + 0.2e1 * t236 * t88 + 0.2e1 * t66 * t84, -0.2e1 * t191 * t42 + 0.2e1 * t236 * t89 - 0.2e1 * t66 * t83, -0.2e1 * t318 * t83 + 0.2e1 * t42 * t88 + 0.2e1 * t43 * t89 - 0.2e1 * t55 * t84, 0.2e1 * t236 * t66 + 0.2e1 * t318 * t43 - 0.2e1 * t42 * t55, 0.2e1 * t70 * t58, -0.2e1 * t57 * t70 - 0.2e1 * t58 * t69, 0.2e1 * t58 * t88 + 0.2e1 * t70 * t84, t39, -0.2e1 * t57 * t88 - 0.2e1 * t69 * t84, t269, 0.2e1 * t15 * t88 + 0.2e1 * t26 * t84 + 0.2e1 * t43 * t69 + 0.2e1 * t52 * t57, 0.2e1 * t14 * t88 - 0.2e1 * t27 * t84 + 0.2e1 * t43 * t70 + 0.2e1 * t52 * t58, 0.2e1 * t14 * t69 - 0.2e1 * t15 * t70 - 0.2e1 * t26 * t58 - 0.2e1 * t27 * t57, -0.2e1 * t14 * t27 + 0.2e1 * t15 * t26 + 0.2e1 * t43 * t52, t21, -0.2e1 * t225, t13, t270, -0.2e1 * t224, t39, 0.2e1 * t12 * t59 + 0.2e1 * t23 * t32 + 0.2e1 * t4 * t69 + 0.2e1 * t57 * t9, -0.2e1 * t10 * t57 + 0.2e1 * t12 * t60 + 0.2e1 * t23 * t33 + 0.2e1 * t3 * t69, -0.2e1 * t10 * t32 + 0.2e1 * t3 * t59 - 0.2e1 * t33 * t9 - 0.2e1 * t4 * t60, -0.2e1 * t10 * t3 + 0.2e1 * t12 * t23 + 0.2e1 * t4 * t9, t21, t13, 0.2e1 * t225, t39, 0.2e1 * t224, t270, 0.2e1 * t16 * t32 - 0.2e1 * t2 * t69 + 0.2e1 * t5 * t59 - 0.2e1 * t57 * t7, -0.2e1 * t1 * t59 + 0.2e1 * t2 * t60 - 0.2e1 * t32 * t6 + 0.2e1 * t33 * t7, 0.2e1 * t1 * t69 - 0.2e1 * t16 * t33 - 0.2e1 * t5 * t60 + 0.2e1 * t57 * t6, 0.2e1 * t1 * t6 + 0.2e1 * t16 * t5 + 0.2e1 * t2 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191 * t232 + t290 * t84, t191 * t233 - t290 * t83 (-t308 * t84 + t309 * t83 + (t308 * t89 - t309 * t88) * qJD(3)) * t150 (t207 - t308 * t42 - t309 * t43 + (t308 * t318 + t309 * t55) * qJD(3)) * t150, 0, 0, 0, 0, 0, 0, -t112 * t84 - t92 * t88 + (t69 * t245 - t309 * t57) * t150, -t113 * t84 - t167 * t88 + t70 * t232 - t58 * t260, t112 * t58 - t113 * t57 - t167 * t69 + t92 * t70, -t15 * t112 - t14 * t113 + t27 * t167 + t52 * t232 - t26 * t92 - t43 * t260, 0, 0, 0, 0, 0, 0, t193, t181, t194, -t10 * t64 + t112 * t12 + t23 * t92 - t3 * t94 - t4 * t93 - t65 * t9, 0, 0, 0, 0, 0, 0, t193, t194, -t181, t1 * t94 + t112 * t5 + t16 * t92 + t2 * t93 - t6 * t64 + t65 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t242 * t279 + 0.2e1 * t77 + 0.2e1 * (t113 * (-qJD(4) * t259 + t231) - t309 * t232) * t150, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, 0, -t84, 0, -t43, t42, 0, 0, t70 * t276 + t293, -t154 * t57 + t58 * t156 + (-t154 * t70 - t156 * t69) * qJD(4), t203, t38, -t202, 0, -pkin(3) * t57 - pkin(10) * t203 + t52 * t143 - t43 * t156, -pkin(3) * t58 + pkin(10) * t202 + t43 * t154 + t52 * t276 (t293 - t294 + (t154 * t69 + t156 * t70) * qJD(4)) * pkin(10) + t173, -t43 * pkin(3) + pkin(10) * t173, t17, t323, t8, t169, -t313, t38, t107 * t57 + t76 * t69 + (-t4 + (pkin(10) * t59 + t153 * t23) * qJD(4)) * t156 + (pkin(10) * t32 + qJD(4) * t9 + t201) * t154, -t284 * t57 + t75 * t69 + (-t3 + (pkin(10) * t60 + t155 * t23) * qJD(4)) * t156 + (pkin(10) * t33 - qJD(4) * t10 - t200) * t154, -t107 * t33 - t284 * t32 + t75 * t59 - t76 * t60 + t223 * t276 + (t153 * t3 - t155 * t4 + (-t10 * t155 + t153 * t9) * qJD(5)) * t154, -t10 * t75 + t4 * t107 - t3 * t284 + t9 * t76 + (t12 * t154 + t23 * t276) * pkin(10), t17, t8, -t323, t38, t313, t169, -t102 * t57 + t110 * t32 + t74 * t59 - t73 * t69 + (t16 * t278 + t2) * t156 + (-qJD(4) * t7 - t205) * t154, -t100 * t32 + t102 * t33 - t68 * t59 + t73 * t60 + t226 * t276 + (-t1 * t153 + t155 * t2 + (-t153 * t7 - t155 * t6) * qJD(5)) * t154, t100 * t57 - t110 * t33 - t74 * t60 + t68 * t69 + (-t16 * t277 - t1) * t156 + (qJD(4) * t6 + t204) * t154, t1 * t100 + t102 * t2 + t110 * t5 + t16 * t74 + t6 * t68 + t7 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t232, -t233, 0, 0, 0, 0, 0, 0, 0, 0 (-t154 * t244 - t156 * t245) * t150 (t154 * t245 - t156 * t244) * t150, t324, -pkin(3) * t232 + t324 * pkin(10), 0, 0, 0, 0, 0, 0, t183, t34, t161, pkin(10) * t196 - t65 * t107 - t284 * t64 - t94 * t75 - t93 * t76, 0, 0, 0, 0, 0, 0, t183, t161, -t34, -t100 * t64 + t102 * t65 + t110 * t92 + t112 * t74 + t68 * t94 + t73 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t239, -0.2e1 * t240, 0, t138, 0, 0, t154 * t268, t156 * t268, 0, 0, t106, 0.2e1 * t95, t96, t105, 0.2e1 * t99, t138, 0.2e1 * t107 * t143 - 0.2e1 * t76 * t156 + 0.2e1 * (t142 * t146 + t153 * t239) * pkin(10), -0.2e1 * t284 * t143 - 0.2e1 * t75 * t156 + 0.2e1 * (-t146 * t275 + t155 * t239) * pkin(10), 0.2e1 * t212 * t276 + 0.2e1 * (t153 * t75 - t155 * t76 + (t107 * t153 - t155 * t284) * qJD(5)) * t154, 0.2e1 * pkin(10) ^ 2 * t247 + 0.2e1 * t107 * t76 - 0.2e1 * t284 * t75, t106, t96, -0.2e1 * t95, t138, -0.2e1 * t99, t105, 0.2e1 * (t110 * t278 + t73) * t156 + 0.2e1 * (-qJD(4) * t102 + t110 * t142 + t74 * t153) * t154, 0.2e1 * t213 * t276 + 0.2e1 * (-t153 * t68 + t155 * t73 + (-t100 * t155 - t102 * t153) * qJD(5)) * t154, 0.2e1 * (-t110 * t277 - t68) * t156 + 0.2e1 * (qJD(4) * t100 + t110 * t275 - t74 * t155) * t154, 0.2e1 * t100 * t68 + 0.2e1 * t102 * t73 + 0.2e1 * t110 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, -t57, t84, t15, t14, 0, 0, t18, -t315, t37, t199, -t36, 0, -pkin(4) * t32 - t188 + t200, -pkin(4) * t33 + t201 + t321, -t301 - t31 + (-t4 + t311) * t153 + (t197 + t223) * qJD(5), -t12 * pkin(4) + (qJD(5) * t223 - t4 * t153 - t301) * pkin(11), t18, t37, t315, 0, t36, t199, t111 * t59 + t126 * t32 - t188 + t204, t303 - t31 + (t2 + t311) * t153 + (t197 + t226) * qJD(5), -t111 * t60 - t126 * t33 + t205 - t321, t16 * t111 + t5 * t126 + (qJD(5) * t226 + t2 * t153 + t303) * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t167, 0, 0, 0, 0, 0, 0, 0, 0, t72, t71, t178, -t92 * pkin(4) + t172, 0, 0, 0, 0, 0, 0, t72, t178, -t71, t112 * t111 + t92 * t126 + t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t276, 0, -t143, 0, -t264, pkin(10) * t143, 0, 0, -t101, -t98, t117, t101, t115, 0, t135 + (t306 - t307) * t274 + (t153 * t228 - t140) * qJD(4) (pkin(10) * t154 * t155 + t153 * t227) * qJD(5) + (t155 * t228 + t266) * qJD(4), t170, -pkin(4) * t264 + pkin(11) * t170, -t101, t117, t98, 0, -t115, t101, t135 + (t126 * t274 - t74) * t155 - t317 * t153, t171 (-t74 + (t126 * t154 + t304) * qJD(5)) * t153 + t317 * t155, pkin(11) * t171 + t110 * t111 + t74 * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, -0.2e1 * t130, 0, t137, 0, 0, t153 * t267, t155 * t267, 0, 0, t136, 0, 0.2e1 * t130, 0, 0, t137, -0.2e1 * t111 * t155 + 0.2e1 * t126 * t275, 0, -0.2e1 * t111 * t153 - 0.2e1 * t126 * t142, 0.2e1 * t126 * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, -t32, t57, t4, t3, 0, 0, 0, t33, 0, t57, t32, 0, t157 + 0.2e1 * t310, -pkin(5) * t33 - qJ(6) * t32 - qJD(6) * t59, 0.2e1 * t289 - t3 + 0.2e1 * t300, -pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t64, 0, 0, 0, 0, 0, 0, 0, 0, -t65, 0, -t64, -pkin(5) * t65 - qJ(6) * t64 + qJD(6) * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, 0, -t116, t143, t76, t75, 0, 0, 0, t114, 0, t143, t116, 0, t76 + 0.2e1 * t265 (-pkin(5) * t276 - qJ(6) * t274) * t155 + (-qJ(6) * t276 + (pkin(5) * qJD(5) - qJD(6)) * t154) * t153, -t75 + 0.2e1 * t256 - 0.2e1 * t271, -pkin(5) * t73 + qJ(6) * t68 + qJD(6) * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, 0, -t275, 0, -t262, t263, 0, 0, 0, t142, 0, 0, t275, 0, -t262, -t314, -t263, -t314 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t312, qJ(6) * t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t33, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t114, 0, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, 0, t262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t11;
