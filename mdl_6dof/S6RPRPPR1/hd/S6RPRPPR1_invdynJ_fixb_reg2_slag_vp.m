% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRPPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:39:42
% EndTime: 2019-03-09 02:39:53
% DurationCPUTime: 7.26s
% Computational Cost: add. (8933->524), mult. (19705->652), div. (0->0), fcn. (14348->18), ass. (0->265)
t212 = sin(pkin(10));
t219 = sin(qJ(3));
t221 = cos(qJ(3));
t324 = cos(pkin(10));
t166 = t212 * t221 + t219 * t324;
t156 = t166 * qJD(1);
t211 = sin(pkin(11));
t214 = cos(pkin(11));
t138 = qJD(3) * t211 + t156 * t214;
t218 = sin(qJ(6));
t349 = cos(qJ(6));
t305 = t211 * t156;
t367 = qJD(3) * t214 - t305;
t241 = t349 * t367;
t75 = -t218 * t138 + t241;
t370 = t75 ^ 2;
t279 = t324 * t221;
t181 = qJD(1) * t279;
t298 = qJD(1) * t219;
t153 = t212 * t298 - t181;
t150 = qJD(6) + t153;
t369 = t150 * t75;
t74 = t349 * t138 + t218 * t367;
t368 = t74 ^ 2;
t366 = 0.2e1 * qJD(3);
t285 = t349 * t214;
t242 = -t218 * t211 + t285;
t295 = qJD(6) * t218;
t355 = qJD(6) * t285 - t211 * t295;
t326 = -t242 * t153 - t355;
t167 = t349 * t211 + t218 * t214;
t160 = t167 * qJD(6);
t325 = t167 * t153 + t160;
t213 = sin(pkin(9));
t189 = pkin(1) * t213 + pkin(7);
t175 = t189 * qJDD(1);
t273 = -qJD(2) * qJD(3) - t175;
t207 = qJ(3) + pkin(10);
t196 = sin(t207);
t199 = cos(t207);
t208 = qJ(1) + pkin(9);
t197 = sin(t208);
t200 = cos(t208);
t267 = g(1) * t200 + g(2) * t197;
t363 = -g(3) * t199 + t267 * t196;
t202 = t221 * qJDD(2);
t292 = qJD(1) * qJD(4);
t177 = t189 * qJD(1);
t275 = qJ(4) * qJD(1) + t177;
t359 = t221 * t275;
t63 = qJDD(3) * pkin(3) + t202 - qJD(3) * t359 + (-qJ(4) * qJDD(1) + t273 - t292) * t219;
t286 = -t219 * qJDD(2) + t221 * t273;
t296 = qJD(3) * t219;
t110 = -t177 * t296 - t286;
t293 = qJD(1) * qJD(3);
t283 = t219 * t293;
t289 = t221 * qJDD(1);
t70 = t221 * t292 + (-t283 + t289) * qJ(4) + t110;
t31 = -t212 * t70 + t324 * t63;
t30 = -qJDD(3) * pkin(4) + qJDD(5) - t31;
t226 = t30 - t363;
t323 = pkin(1) * qJDD(1);
t362 = t156 * t367;
t229 = qJDD(1) * t166 - t212 * t283;
t120 = qJD(3) * t181 + t229;
t101 = qJDD(3) * t211 + t120 * t214;
t240 = -t212 * t219 + t279;
t158 = t240 * qJD(3);
t255 = t101 * t166 + t138 * t158;
t361 = t211 * t255;
t155 = t166 * qJD(3);
t290 = t219 * qJDD(1);
t261 = -qJDD(1) * t279 + t212 * t290;
t119 = qJD(1) * t155 + t261;
t327 = -t119 * t166 - t153 * t158;
t360 = t214 * t327;
t347 = g(1) * t197;
t281 = g(2) * t200 - t347;
t243 = t281 * t199;
t215 = cos(pkin(9));
t192 = -pkin(1) * t215 - pkin(2);
t204 = t221 * pkin(3);
t357 = t192 - t204;
t203 = t221 * qJD(2);
t141 = -t219 * t275 + t203;
t294 = t219 * qJD(2);
t142 = t294 + t359;
t280 = t324 * t142;
t78 = t141 * t212 + t280;
t356 = t78 * qJD(3) + t363;
t100 = -qJDD(3) * t214 + t120 * t211;
t21 = -qJD(6) * t241 + t218 * t100 - t349 * t101 + t138 * t295;
t354 = -t21 * t242 - t325 * t74;
t117 = qJDD(6) + t119;
t353 = t167 * t117 - t150 * t326;
t151 = t153 ^ 2;
t321 = t119 * t211;
t352 = -t151 * t214 - t321;
t351 = -t100 * t240 - t155 * t367;
t344 = g(3) * t196;
t232 = -t267 * t199 - t344;
t131 = qJD(3) * pkin(3) + t141;
t69 = t212 * t131 + t280;
t60 = qJD(3) * qJ(5) + t69;
t152 = qJD(1) * t357 + qJD(4);
t88 = t153 * pkin(4) - t156 * qJ(5) + t152;
t34 = -t211 * t60 + t214 * t88;
t20 = pkin(5) * t153 - pkin(8) * t138 + t34;
t35 = t211 * t88 + t214 * t60;
t24 = pkin(8) * t367 + t35;
t246 = -t349 * t20 + t218 * t24;
t32 = t212 * t63 + t324 * t70;
t29 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t32;
t140 = pkin(3) * t283 + qJDD(1) * t357 + qJDD(4);
t41 = t119 * pkin(4) - t120 * qJ(5) - t156 * qJD(5) + t140;
t12 = -t211 * t29 + t214 * t41;
t6 = pkin(5) * t119 - pkin(8) * t101 + t12;
t13 = t211 * t41 + t214 * t29;
t9 = -pkin(8) * t100 + t13;
t1 = -qJD(6) * t246 + t218 * t6 + t349 * t9;
t350 = t156 ^ 2;
t348 = pkin(3) * t219;
t343 = g(3) * t221;
t342 = t100 * pkin(5);
t341 = t214 * pkin(5);
t220 = sin(qJ(1));
t340 = t220 * pkin(1);
t339 = t74 * t75;
t184 = pkin(3) * t212 + qJ(5);
t338 = pkin(8) + t184;
t109 = t242 * t166;
t278 = t349 * t100 + t218 * t101;
t22 = qJD(6) * t74 + t278;
t47 = -t158 * t242 + t160 * t166;
t337 = -t109 * t22 - t47 * t75;
t108 = t167 * t166;
t48 = t158 * t167 + t166 * t355;
t336 = -t108 * t117 - t48 * t150;
t335 = t74 * t155 + t21 * t240;
t105 = pkin(3) * t298 + pkin(4) * t156 + qJ(5) * t153;
t128 = t212 * t142;
t79 = t141 * t324 - t128;
t43 = t211 * t105 + t214 * t79;
t161 = t338 * t211;
t162 = t338 * t214;
t113 = -t349 * t161 - t218 * t162;
t316 = t153 * t214;
t42 = t214 * t105 - t211 * t79;
t26 = pkin(5) * t156 + pkin(8) * t316 + t42;
t317 = t153 * t211;
t33 = pkin(8) * t317 + t43;
t334 = qJD(5) * t242 + qJD(6) * t113 - t218 * t26 - t349 * t33;
t114 = -t218 * t161 + t349 * t162;
t333 = -qJD(5) * t167 - qJD(6) * t114 + t218 * t33 - t349 * t26;
t288 = pkin(3) * t296;
t77 = pkin(4) * t155 - qJ(5) * t158 - qJD(5) * t166 + t288;
t303 = qJ(4) + t189;
t276 = qJD(3) * t303;
t144 = t221 * qJD(4) - t219 * t276;
t235 = -t219 * qJD(4) - t221 * t276;
t87 = t144 * t324 + t212 * t235;
t40 = t211 * t77 + t214 * t87;
t249 = t367 * t158;
t312 = t166 * t214;
t332 = -t100 * t312 + t214 * t249;
t331 = t156 * t75;
t329 = t74 * t156;
t328 = -t101 * t240 + t138 * t155;
t322 = t101 * t214;
t320 = t138 * t156;
t319 = t138 * t211;
t318 = t153 * t156;
t315 = t158 * t211;
t313 = t166 * t211;
t310 = t177 * t219;
t309 = t177 * t221;
t308 = t196 * t200;
t307 = t197 * t199;
t306 = t199 * t200;
t112 = t214 * t119;
t68 = t131 * t324 - t128;
t59 = -qJD(3) * pkin(4) + qJD(5) - t68;
t302 = -qJD(5) + t59;
t107 = -pkin(4) * t240 - qJ(5) * t166 + t357;
t163 = t303 * t221;
t277 = t303 * t219;
t116 = t163 * t324 - t212 * t277;
t50 = t211 * t107 + t214 * t116;
t301 = -t211 * t151 + t112;
t194 = t204 + pkin(2);
t222 = cos(qJ(1));
t205 = t222 * pkin(1);
t300 = t200 * t194 + t205;
t209 = t219 ^ 2;
t210 = t221 ^ 2;
t299 = t209 - t210;
t178 = qJD(1) * t192;
t176 = qJDD(1) * t192;
t224 = qJD(1) ^ 2;
t287 = t219 * t224 * t221;
t216 = -qJ(4) - pkin(7);
t282 = pkin(5) * t211 - t216;
t39 = -t211 * t87 + t214 * t77;
t49 = t214 * t107 - t116 * t211;
t86 = t144 * t212 - t324 * t235;
t115 = t163 * t212 + t324 * t277;
t272 = -t167 * t22 - t326 * t75;
t270 = t221 * t283;
t269 = g(2) * t308 - t196 * t347;
t268 = t242 * t117 - t150 * t325;
t191 = -pkin(3) * t324 - pkin(4);
t265 = g(1) * t220 - g(2) * t222;
t264 = -t108 * t21 + t48 * t74;
t263 = -t200 * t216 - t340;
t262 = -pkin(4) * t199 - qJ(5) * t196;
t260 = -t12 * t214 - t13 * t211;
t259 = -t12 * t211 + t13 * t214;
t258 = t155 * t75 + t22 * t240;
t257 = -t211 * t34 + t214 * t35;
t256 = -t109 * t117 + t150 * t47;
t254 = -t119 * t240 + t153 * t155;
t252 = t120 * t240 - t155 * t156;
t190 = pkin(4) + t341;
t217 = -pkin(8) - qJ(5);
t251 = t190 * t199 - t196 * t217;
t250 = t214 * t367;
t148 = t294 + t309;
t38 = -pkin(5) * t240 - pkin(8) * t312 + t49;
t45 = -pkin(8) * t313 + t50;
t14 = -t218 * t45 + t349 * t38;
t8 = t218 * t20 + t349 * t24;
t15 = t218 * t38 + t349 * t45;
t239 = t327 * t211;
t237 = -qJD(1) * t178 + t267;
t234 = -qJDD(3) * t189 + t178 * t366;
t230 = t158 * t59 + t166 * t30 - t267;
t223 = qJD(3) ^ 2;
t228 = -t189 * t223 - 0.2e1 * t176 - t281;
t2 = -qJD(6) * t8 - t218 * t9 + t349 * t6;
t111 = -t148 * qJD(3) - t219 * t175 + t202;
t147 = t203 - t310;
t227 = t110 * t221 - t111 * t219 + (-t147 * t221 - t148 * t219) * qJD(3);
t206 = pkin(11) + qJ(6);
t198 = cos(t206);
t195 = sin(t206);
t174 = qJDD(3) * t221 - t219 * t223;
t173 = qJDD(3) * t219 + t221 * t223;
t171 = t191 - t341;
t135 = t195 * t197 + t198 * t306;
t134 = -t195 * t306 + t197 * t198;
t133 = t195 * t200 - t198 * t307;
t132 = t195 * t307 + t198 * t200;
t122 = qJD(3) * t158 + qJDD(3) * t166;
t121 = -qJD(3) * t155 + qJDD(3) * t240;
t91 = t211 * t100;
t85 = pkin(5) * t313 + t115;
t55 = pkin(5) * t315 + t86;
t51 = -pkin(5) * t317 + t78;
t46 = -pkin(5) * t367 + t59;
t28 = -pkin(8) * t315 + t40;
t23 = -pkin(8) * t158 * t214 + pkin(5) * t155 + t39;
t17 = t30 + t342;
t4 = -t15 * qJD(6) - t218 * t28 + t349 * t23;
t3 = t14 * qJD(6) + t218 * t23 + t349 * t28;
t5 = [0, 0, 0, 0, 0, qJDD(1), t265, g(1) * t222 + g(2) * t220, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t215 * t323 - t281, -0.2e1 * t213 * t323 + t267, 0 (t265 + (t213 ^ 2 + t215 ^ 2) * t323) * pkin(1), qJDD(1) * t209 + 0.2e1 * t270, 0.2e1 * t219 * t289 - 0.2e1 * t293 * t299, t173, qJDD(1) * t210 - 0.2e1 * t270, t174, 0, t219 * t234 + t221 * t228, -t219 * t228 + t221 * t234 (t209 + t210) * t175 + t227 - t267, t176 * t192 - g(1) * (-pkin(2) * t197 + pkin(7) * t200 - t340) - g(2) * (pkin(2) * t200 + pkin(7) * t197 + t205) + t227 * t189, t120 * t166 + t156 * t158, t252 + t327, t122, t254, t121, 0, -t115 * qJDD(3) + t357 * t119 - t140 * t240 + t152 * t155 - t243 + (t153 * t348 - t86) * qJD(3), -t116 * qJDD(3) + t357 * t120 + t140 * t166 + t152 * t158 + (t156 * t348 - t87) * qJD(3) + t269, t115 * t120 - t116 * t119 - t153 * t87 - t155 * t69 + t156 * t86 - t158 * t68 - t166 * t31 + t240 * t32 - t267, t32 * t116 + t69 * t87 - t31 * t115 - t68 * t86 + t140 * t357 + t152 * t288 - g(1) * (-t197 * t194 + t263) - g(2) * (-t197 * t216 + t300) t255 * t214, t332 - t361, t328 - t360 (t100 * t166 - t249) * t211, t239 - t351, t254, t115 * t100 + t49 * t119 - t12 * t240 + t39 * t153 + t34 * t155 + t211 * t230 - t214 * t243 - t367 * t86, t115 * t101 - t50 * t119 + t13 * t240 + t86 * t138 - t40 * t153 - t35 * t155 + t211 * t243 + t214 * t230, t40 * t367 - t50 * t100 - t39 * t138 - t49 * t101 + t260 * t166 + (-t211 * t35 - t214 * t34) * t158 - t269, t13 * t50 + t35 * t40 + t12 * t49 + t34 * t39 + t30 * t115 + t59 * t86 - g(1) * t263 - g(2) * (pkin(4) * t306 + qJ(5) * t308 + t300) + (-g(1) * (-t194 + t262) + g(2) * t216) * t197, -t109 * t21 - t47 * t74, -t264 + t337, -t256 + t335, t108 * t22 - t48 * t75, t258 + t336, -t117 * t240 + t150 * t155, -g(1) * t133 - g(2) * t135 + t108 * t17 + t117 * t14 + t150 * t4 - t155 * t246 - t2 * t240 + t22 * t85 + t46 * t48 - t55 * t75, -g(1) * t132 - g(2) * t134 + t1 * t240 + t109 * t17 - t117 * t15 - t150 * t3 - t155 * t8 - t21 * t85 - t46 * t47 + t55 * t74, -t1 * t108 - t109 * t2 + t14 * t21 - t15 * t22 - t246 * t47 + t3 * t75 - t4 * t74 - t48 * t8 - t269, t1 * t15 + t8 * t3 + t2 * t14 - t246 * t4 + t17 * t85 + t46 * t55 + g(1) * t340 - g(2) * t300 + (-g(1) * t282 - g(2) * t251) * t200 + (-g(1) * (-t194 - t251) - g(2) * t282) * t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t174, -t173, 0, t110 * t219 + t111 * t221 - g(3) + (-t147 * t219 + t148 * t221) * qJD(3), 0, 0, 0, 0, 0, 0, t121, -t122, -t252 + t327, -t155 * t68 + t158 * t69 + t166 * t32 + t240 * t31 - g(3), 0, 0, 0, 0, 0, 0, t239 + t351, t328 + t360, t332 + t361, t59 * t155 + t158 * t257 + t166 * t259 - t240 * t30 - g(3), 0, 0, 0, 0, 0, 0, -t258 + t336, t256 + t335, t264 + t337, t1 * t109 - t108 * t2 + t155 * t46 - t17 * t240 + t246 * t48 - t47 * t8 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t287, t299 * t224, t290, t287, t289, qJDD(3), -t343 + t202 + (t148 - t309) * qJD(3) + (t237 + t273) * t219, g(3) * t219 + (t147 + t310) * qJD(3) + t237 * t221 + t286, 0, 0, t318, -t151 + t350 (t181 + t153) * qJD(3) + t229, -t318, -t261, qJDD(3), -t152 * t156 + (qJDD(3) * t324 - t153 * t298) * pkin(3) + t31 + t356, t79 * qJD(3) + t152 * t153 + (-qJDD(3) * t212 - t156 * t298) * pkin(3) - t32 - t232 (t69 - t78) * t156 + (-t68 + t79) * t153 + (-t119 * t212 - t120 * t324) * pkin(3), t68 * t78 - t69 * t79 + (t324 * t31 - t343 + t212 * t32 + (-qJD(1) * t152 + t267) * t219) * pkin(3), t101 * t211 + t138 * t316, t322 - t91 + (t250 - t319) * t153, -t320 - t352, -t100 * t214 - t317 * t367, t301 - t362, -t318, -t184 * t321 + t191 * t100 - t78 * t305 - t34 * t156 + (t211 * t302 - t42) * t153 + (-t30 + t356) * t214, -t184 * t112 + t191 * t101 - t78 * t138 + t35 * t156 + (t214 * t302 + t43) * t153 + t226 * t211, t42 * t138 + t43 * t305 + (-qJD(5) * t305 - t184 * t100 - t34 * t153 + t13 + (qJD(5) * t214 - t43) * qJD(3)) * t214 + (qJD(5) * t138 + t101 * t184 - t153 * t35 - t12) * t211 + t232, t30 * t191 - t35 * t43 - t34 * t42 - t59 * t78 - g(3) * (t204 - t262) + t259 * t184 + t257 * qJD(5) + t267 * (pkin(4) * t196 - qJ(5) * t199 + t348) -t21 * t167 - t326 * t74, t272 + t354, -t329 + t353, -t22 * t242 - t325 * t75, t268 - t331, -t150 * t156, t113 * t117 + t333 * t150 + t246 * t156 - t17 * t242 + t171 * t22 + t198 * t363 + t325 * t46 + t51 * t75, -t114 * t117 - t150 * t334 + t8 * t156 + t17 * t167 - t171 * t21 - t195 * t363 - t326 * t46 - t51 * t74, t1 * t242 + t113 * t21 - t114 * t22 - t2 * t167 - t246 * t326 - t325 * t8 - t333 * t74 + t334 * t75 + t232, t1 * t114 + t2 * t113 + t17 * t171 - t46 * t51 - g(3) * (t204 + t251) + t334 * t8 - t333 * t246 + t267 * (t190 * t196 + t199 * t217 + t348); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156 * t366 + t261 (t181 - t153) * qJD(3) + t229, -t151 - t350, t69 * t153 + t68 * t156 + t140 + t281, 0, 0, 0, 0, 0, 0, t301 + t362, -t320 + t352, -t322 - t91 + (t250 + t319) * t153, t153 * t257 - t59 * t156 - t260 + t281, 0, 0, 0, 0, 0, 0, t268 + t331, -t329 - t353, t272 - t354, t1 * t167 - t46 * t156 + t2 * t242 + t246 * t325 - t326 * t8 + t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138 * t153 + t100, t153 * t367 + t101, -t138 ^ 2 - t367 ^ 2, t138 * t34 - t35 * t367 + t226, 0, 0, 0, 0, 0, 0, t74 * t150 + t22, -t21 + t369, -t368 - t370, -t246 * t74 - t8 * t75 + t226 + t342; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t339, t368 - t370, -t21 - t369, t339, -t278 + (-qJD(6) + t150) * t74, t117, -g(1) * t134 + g(2) * t132 + t8 * t150 + t195 * t344 - t46 * t74 + t2, g(1) * t135 - g(2) * t133 - t150 * t246 + t198 * t344 - t46 * t75 - t1, 0, 0;];
tau_reg  = t5;
