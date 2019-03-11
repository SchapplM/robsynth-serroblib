% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PPRRPR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PPRRPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:51:00
% EndTime: 2019-03-08 18:51:08
% DurationCPUTime: 4.76s
% Computational Cost: add. (5299->482), mult. (13436->653), div. (0->0), fcn. (12219->14), ass. (0->249)
t167 = sin(qJ(4));
t284 = qJD(4) * t167;
t312 = cos(pkin(6));
t148 = qJD(1) * t312 + qJD(2);
t162 = sin(pkin(12));
t164 = sin(pkin(6));
t168 = sin(qJ(3));
t165 = cos(pkin(12));
t311 = cos(pkin(7));
t248 = t165 * t311;
t238 = t168 * t248;
t325 = cos(qJ(3));
t189 = (t325 * t162 + t238) * t164;
t163 = sin(pkin(7));
t298 = t163 * t168;
t69 = qJD(1) * t189 + t148 * t298;
t347 = pkin(4) * t284 - t69;
t277 = qJD(3) * qJD(4);
t259 = t167 * t277;
t170 = cos(qJ(4));
t271 = t170 * qJDD(3);
t344 = -t259 + t271;
t307 = qJ(5) * t170;
t231 = pkin(10) * t167 - t307;
t278 = t167 * qJD(5);
t192 = qJD(4) * t231 - t278;
t346 = -t192 - t347;
t67 = qJD(3) * pkin(9) + t69;
t266 = t164 * t165 * t163;
t98 = -qJD(1) * t266 + t148 * t311;
t92 = t170 * t98;
t46 = t167 * t67 - t92;
t345 = -qJD(5) - t46;
t106 = t167 * t311 + t170 * t298;
t173 = qJD(3) ^ 2;
t205 = t325 * qJDD(3) - t168 * t173;
t260 = qJD(3) * t325;
t240 = qJD(4) * t260;
t105 = t167 * t298 - t170 * t311;
t243 = t163 * t260;
t84 = qJD(4) * t105 - t170 * t243;
t343 = -t84 * qJD(4) + t106 * qJDD(4) + t163 * (t167 * t205 + t170 * t240);
t85 = qJD(4) * t106 + t167 * t243;
t342 = -t85 * qJD(4) - t105 * qJDD(4) + t163 * (-t167 * t240 + t170 * t205);
t171 = -pkin(4) - pkin(10);
t296 = t167 * qJ(5);
t256 = -pkin(3) - t296;
t111 = t171 * t170 + t256;
t144 = t312 * qJDD(1) + qJDD(2);
t218 = t325 * t248;
t206 = t164 * t218;
t224 = t164 * t238;
t300 = t162 * t164;
t264 = qJD(1) * t300;
t242 = qJD(3) * t264;
t257 = qJDD(1) * t300;
t290 = qJD(3) * t168;
t263 = t163 * t290;
t265 = t163 * t325;
t215 = qJD(3) * qJD(1) * t224 - qJDD(1) * t206 - t144 * t265 + t148 * t263 + t168 * t257 + t325 * t242;
t203 = pkin(4) * t259 + t215;
t14 = qJD(3) * t192 + qJDD(3) * t111 + t203;
t166 = sin(qJ(6));
t169 = cos(qJ(6));
t294 = qJD(5) - t92 + (pkin(5) * qJD(3) + t67) * t167;
t34 = t171 * qJD(4) + t294;
t131 = t168 * t264;
t201 = qJD(1) * t206;
t181 = -t148 * t265 + t131 - t201;
t52 = qJD(3) * t111 + t181;
t229 = t166 * t52 - t169 * t34;
t258 = t170 * t277;
t272 = t167 * qJDD(3);
t204 = t258 + t272;
t282 = qJD(4) * t170;
t227 = -qJD(3) * t201 - qJDD(1) * t224 - t144 * t298 - t148 * t243 - t325 * t257;
t35 = -t168 * t242 - t227;
t32 = qJDD(3) * pkin(9) + t35;
t95 = -qJDD(1) * t266 + t144 * t311;
t245 = t167 * t32 - t170 * t95 + t67 * t282 + t98 * t284;
t235 = qJDD(5) + t245;
t6 = t204 * pkin(5) + t171 * qJDD(4) + t235;
t1 = -t229 * qJD(6) + t169 * t14 + t166 * t6;
t291 = qJD(3) * t167;
t150 = qJD(6) + t291;
t341 = t150 * t229 + t1;
t13 = t166 * t34 + t169 * t52;
t2 = -qJD(6) * t13 - t166 * t14 + t169 * t6;
t340 = t13 * t150 + t2;
t40 = -qJD(4) * pkin(4) - t345;
t47 = t167 * t98 + t170 * t67;
t43 = -qJD(4) * qJ(5) - t47;
t285 = qJD(4) * t166;
t289 = qJD(3) * t170;
t119 = t169 * t289 + t285;
t225 = t119 * t150;
t70 = qJD(6) * t119 - t169 * qJDD(4) + t166 * t344;
t339 = t70 - t225;
t261 = t166 * t289;
t283 = qJD(4) * t169;
t121 = -t261 + t283;
t301 = t121 * t150;
t71 = -qJD(6) * t261 + qJDD(4) * t166 + (qJD(4) * qJD(6) + t344) * t169;
t338 = -t71 + t301;
t118 = qJDD(6) + t204;
t104 = t169 * t118;
t281 = qJD(6) * t166;
t337 = -t150 * t281 + t104;
t310 = cos(pkin(11));
t233 = t312 * t310;
t309 = sin(pkin(11));
t187 = t162 * t309 - t165 * t233;
t250 = t164 * t310;
t336 = t163 * t250 + t187 * t311;
t232 = t312 * t309;
t188 = t162 * t310 + t165 * t232;
t249 = t164 * t309;
t335 = -t163 * t249 + t188 * t311;
t132 = -pkin(4) * t170 + t256;
t292 = qJD(3) * t132;
t63 = t181 + t292;
t334 = t63 * t291 + qJDD(5);
t160 = t167 ^ 2;
t161 = t170 ^ 2;
t274 = qJDD(3) * t161;
t275 = qJDD(3) * t160;
t333 = -(-t160 - t161) * t181 * qJD(3) + (t274 + t275) * pkin(9);
t154 = pkin(5) * t289;
t37 = t154 - t43;
t332 = t118 * t171 + t150 * t37;
t211 = -qJ(5) * t282 - t278;
t276 = qJDD(3) * t132;
t15 = qJD(3) * t211 + t203 + t276;
t102 = t162 * t233 + t165 * t309;
t57 = t102 * t168 + t325 * t336;
t103 = -t162 * t232 + t165 * t310;
t59 = t103 * t168 + t325 * t335;
t251 = t163 * t312;
t219 = t325 * t251;
t299 = t162 * t168;
t80 = t164 * t299 - t206 - t219;
t208 = g(1) * t59 + g(2) * t57 + g(3) * t80;
t313 = t211 + t347;
t172 = qJD(4) ^ 2;
t324 = pkin(9) * t172;
t329 = qJD(3) * t313 + t15 - t208 + t276 + t324;
t199 = t311 * t312 - t266;
t81 = t168 * t251 + t189;
t62 = t167 * t199 + t81 * t170;
t74 = (t219 + (t218 - t299) * t164) * qJD(3);
t18 = qJD(4) * t62 + t74 * t167;
t61 = t167 * t81 - t170 * t199;
t75 = t81 * qJD(3);
t328 = qJD(3) * (-t170 * t75 + t284 * t80) - qJD(4) * t18 - qJDD(4) * t61 - t271 * t80;
t19 = -t81 * t284 + (qJD(4) * t199 + t74) * t170;
t327 = qJD(3) * (t167 * t75 + t282 * t80) - qJD(4) * t19 - qJDD(4) * t62 + t272 * t80;
t326 = pkin(5) + pkin(9);
t139 = t326 * t170;
t128 = qJD(4) * t139;
t297 = t166 * t167;
t138 = t326 * t167;
t82 = -t111 * t166 + t138 * t169;
t323 = qJD(6) * t82 + t128 * t166 - t169 * t346 + t297 * t181;
t295 = t167 * t169;
t83 = t111 * t169 + t138 * t166;
t322 = -qJD(6) * t83 + t128 * t169 + t166 * t346 + t295 * t181;
t321 = qJD(3) * pkin(3);
t318 = t169 * t70;
t317 = t170 * t57;
t316 = t170 * t59;
t315 = t170 * t80;
t314 = t71 * t166;
t308 = pkin(9) * qJDD(4);
t306 = qJDD(3) * pkin(3);
t305 = qJDD(4) * pkin(4);
t304 = t118 * t166;
t302 = t121 * t119;
t293 = t160 - t161;
t287 = qJD(4) * t119;
t286 = qJD(4) * t121;
t280 = qJD(6) * t169;
t279 = qJD(6) * t170;
t273 = qJDD(4) * qJ(5);
t53 = t57 * pkin(3);
t270 = -pkin(4) * t317 - t57 * t296 - t53;
t54 = t59 * pkin(3);
t269 = -pkin(4) * t316 - t59 * t296 - t54;
t79 = t80 * pkin(3);
t268 = -pkin(4) * t315 - t80 * t296 - t79;
t176 = t163 * t187 - t250 * t311;
t58 = t102 * t325 - t168 * t336;
t22 = t167 * t58 - t170 * t176;
t23 = t167 * t176 + t58 * t170;
t254 = -t22 * pkin(4) + qJ(5) * t23;
t177 = t163 * t188 + t249 * t311;
t60 = t103 * t325 - t168 * t335;
t24 = t167 * t60 - t170 * t177;
t25 = t167 * t177 + t60 * t170;
t253 = -t24 * pkin(4) + qJ(5) * t25;
t252 = -t61 * pkin(4) + qJ(5) * t62;
t246 = -t167 * t95 - t170 * t32 - t98 * t282 + t67 * t284;
t241 = t167 * t258;
t230 = t13 * t169 + t166 * t229;
t26 = -t166 * t80 + t169 * t61;
t27 = t166 * t61 + t169 * t80;
t226 = t150 * t166;
t212 = -t150 * t280 - t304;
t210 = g(1) * t24 + g(2) * t22 + g(3) * t61;
t209 = -g(1) * t25 - g(2) * t23 - g(3) * t62;
t207 = g(1) * t60 + g(2) * t58 + g(3) * t81;
t86 = t169 * t105 + t166 * t265;
t202 = -t166 * t105 + t169 * t265;
t66 = t181 - t321;
t200 = -t308 + (t66 - t181 - t321) * qJD(4);
t196 = qJD(3) * t69 + t208;
t195 = qJD(4) * qJD(5) - t246 + t273;
t194 = t308 + (-t63 + t181 - t292) * qJD(4);
t193 = -g(1) * t249 + g(2) * t250 - g(3) * t312;
t191 = t210 - t245;
t190 = t209 - t246;
t7 = pkin(5) * t344 + t195;
t186 = -qJD(6) * t150 * t171 + t209 + t7;
t184 = qJD(4) * t47 + t191;
t33 = t215 - t306;
t180 = t196 + t306 - t33 - t324;
t9 = t235 - t305;
t179 = t9 * t167 + t195 * t170 + (t167 * t43 + t170 * t40) * qJD(4) - t207;
t178 = -t246 * t170 + t245 * t167 + (-t167 * t47 + t170 * t46) * qJD(4) - t207;
t175 = (t167 * t61 + t170 * t62) * qJDD(3) + (t167 * t18 + t170 * t19 + (-t167 * t62 + t170 * t61) * qJD(4)) * qJD(3);
t174 = (t105 * t167 + t106 * t170) * qJDD(3) + (t167 * t85 - t170 * t84 + (t105 * t170 - t106 * t167) * qJD(4)) * qJD(3);
t156 = pkin(4) * t291;
t147 = t167 * t173 * t170;
t136 = t293 * t173;
t134 = qJDD(4) * t170 - t167 * t172;
t133 = qJDD(4) * t167 + t170 * t172;
t127 = t326 * t284;
t122 = -qJ(5) * t289 + t156;
t113 = -0.2e1 * t241 + t274;
t112 = 0.2e1 * t241 + t275;
t100 = qJD(3) * t231 + t156;
t96 = 0.2e1 * t167 * t271 - 0.2e1 * t277 * t293;
t45 = qJD(6) * t202 - t166 * t263 + t169 * t85;
t44 = qJD(6) * t86 + t166 * t85 + t169 * t263;
t39 = t154 + t47;
t17 = t100 * t169 + t166 * t39;
t16 = -t100 * t166 + t169 * t39;
t4 = qJD(6) * t26 + t166 * t18 + t169 * t75;
t3 = -qJD(6) * t27 - t166 * t75 + t169 * t18;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t144 * t312 - g(3) + (t162 ^ 2 + t165 ^ 2) * t164 ^ 2 * qJDD(1), 0, 0, 0, 0, 0, 0, -qJD(3) * t75 - qJDD(3) * t80, -qJD(3) * t74 - qJDD(3) * t81, 0, t181 * t75 + t199 * t95 + t215 * t80 + t35 * t81 + t69 * t74 - g(3), 0, 0, 0, 0, 0, 0, t328, t327, t175, t18 * t46 + t19 * t47 + t245 * t61 - t246 * t62 + t33 * t80 + t66 * t75 - g(3), 0, 0, 0, 0, 0, 0, t175, -t328, -t327, t15 * t80 + t18 * t40 - t19 * t43 + t195 * t62 + t61 * t9 + t63 * t75 - g(3), 0, 0, 0, 0, 0, 0, t118 * t26 + t119 * t19 + t150 * t3 + t62 * t71, -t118 * t27 + t121 * t19 - t150 * t4 - t62 * t70, -t119 * t4 - t121 * t3 + t26 * t70 - t27 * t71, t1 * t27 + t13 * t4 + t19 * t37 + t2 * t26 - t229 * t3 + t62 * t7 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193 + t144, 0, 0, 0, 0, 0, 0, t205 * t163 (-qJDD(3) * t168 - t325 * t173) * t163, 0, t95 * t311 + (-t325 * t215 + t168 * t35 + (t168 * t181 + t325 * t69) * qJD(3)) * t163 + t193, 0, 0, 0, 0, 0, 0, t342, -t343, t174, -t246 * t106 + t245 * t105 + t46 * t85 - t47 * t84 + (t66 * t290 - t325 * t33) * t163 + t193, 0, 0, 0, 0, 0, 0, t174, -t342, t343, t9 * t105 + t195 * t106 + t40 * t85 + t43 * t84 + (-t325 * t15 + t63 * t290) * t163 + t193, 0, 0, 0, 0, 0, 0, t106 * t71 + t118 * t86 - t119 * t84 + t150 * t45, -t106 * t70 + t118 * t202 - t121 * t84 - t150 * t44, -t119 * t44 - t121 * t45 + t202 * t71 + t70 * t86, -t1 * t202 + t7 * t106 + t13 * t44 + t2 * t86 - t229 * t45 - t37 * t84 + t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t196 - t215 (-t181 + t131) * qJD(3) + t207 + t227, 0, 0, t112, t96, t133, t113, t134, 0, t167 * t200 + t170 * t180, -t167 * t180 + t170 * t200, t178 + t333, -t33 * pkin(3) + g(1) * t54 + g(2) * t53 + g(3) * t79 - t66 * t69 - (-t167 * t46 - t170 * t47) * t181 + t178 * pkin(9), 0, -t133, -t134, t112, t96, t113, t179 + t333, t194 * t167 + t170 * t329, -t167 * t329 + t194 * t170, t15 * t132 - g(1) * t269 - g(2) * t270 - g(3) * t268 - (-t167 * t40 + t170 * t43) * t181 + t313 * t63 + t179 * pkin(9), t70 * t166 * t170 + (t166 * t284 - t169 * t279) * t121 (-t119 * t166 + t121 * t169) * t284 + (t314 + t318 + (t119 * t169 + t121 * t166) * qJD(6)) * t170 (t150 * t285 - t70) * t167 + (t212 + t286) * t170, t71 * t169 * t170 + (-t166 * t279 - t167 * t283) * t119 (t150 * t283 - t71) * t167 + (-t287 - t337) * t170, t118 * t167 + t150 * t282, t82 * t118 - t127 * t119 + t139 * t71 - t207 * t169 + (t166 * t208 - t283 * t37 + t2) * t167 + t322 * t150 + (-qJD(4) * t229 + t119 * t181 + t169 * t7 - t281 * t37) * t170, -t83 * t118 - t127 * t121 - t139 * t70 + t207 * t166 + (t169 * t208 + t285 * t37 - t1) * t167 - t323 * t150 + (-qJD(4) * t13 + t121 * t181 - t166 * t7 - t280 * t37) * t170, t70 * t82 - t71 * t83 - t322 * t121 - t323 * t119 + t230 * t284 + (-t1 * t169 + t166 * t2 + (t13 * t166 - t169 * t229) * qJD(6) + t208) * t170, t1 * t83 + t2 * t82 + t7 * t139 - g(1) * (-pkin(10) * t316 + t326 * t60 + t269) - g(2) * (-pkin(10) * t317 + t326 * t58 + t270) - g(3) * (-pkin(10) * t315 + t326 * t81 + t268) + (t170 * t181 - t127) * t37 + t323 * t13 - t322 * t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t136, t272, t147, t271, qJDD(4), -t291 * t66 + t184, -qJD(4) * t46 - t289 * t66 - t190, 0, 0, qJDD(4), -t272, -t271, -t147, t136, t147 (-pkin(4) * t167 + t307) * qJDD(3), -t122 * t289 - t184 - 0.2e1 * t305 + t334, 0.2e1 * t273 + (0.2e1 * qJD(5) + t46) * qJD(4) + (t122 * t167 + t170 * t63) * qJD(3) + t190, -t9 * pkin(4) - g(1) * t253 - g(2) * t254 - g(3) * t252 + t195 * qJ(5) - t63 * t122 + t345 * t43 - t40 * t47, -t121 * t226 - t318 (-t71 - t301) * t169 + (t70 + t225) * t166 (-t121 * t170 - t150 * t297) * qJD(3) + t337, t169 * t225 + t314 (t119 * t170 - t150 * t295) * qJD(3) + t212, -t150 * t289, qJ(5) * t71 + t294 * t119 - t150 * t16 + t186 * t166 + t169 * t332 + t229 * t289, -qJ(5) * t70 + t294 * t121 + t13 * t289 + t150 * t17 - t166 * t332 + t186 * t169, t119 * t17 + t121 * t16 + (-t13 * t291 + t171 * t70 - t2 + (-t119 * t171 - t13) * qJD(6)) * t169 + (-t229 * t291 - t171 * t71 - t1 + (t121 * t171 - t229) * qJD(6)) * t166 + t210, t7 * qJ(5) - t13 * t17 + t229 * t16 - g(1) * (-pkin(10) * t24 + t253) - g(2) * (-pkin(10) * t22 + t254) - g(3) * (-pkin(10) * t61 + t252) + t294 * t37 + (qJD(6) * t230 + t1 * t166 + t2 * t169) * t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272, qJDD(4) + t147, -t160 * t173 - t172, qJD(4) * t43 - t191 - t305 + t334, 0, 0, 0, 0, 0, 0, -t150 * t226 + t104 - t287, -t150 ^ 2 * t169 - t286 - t304, t166 * t338 + t169 * t339, -qJD(4) * t37 + t166 * t341 + t340 * t169 - t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, -t119 ^ 2 + t121 ^ 2, -t339, -t302, t338, t118, -t37 * t121 - g(1) * (-t166 * t59 + t169 * t24) - g(2) * (-t166 * t57 + t169 * t22) - g(3) * t26 + t340, t37 * t119 - g(1) * (-t166 * t24 - t169 * t59) - g(2) * (-t166 * t22 - t169 * t57) + g(3) * t27 - t341, 0, 0;];
tau_reg  = t5;
