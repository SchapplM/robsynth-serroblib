% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRP9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRP9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:29
% EndTime: 2019-12-31 20:07:42
% DurationCPUTime: 5.60s
% Computational Cost: add. (12107->359), mult. (26446->464), div. (0->0), fcn. (18219->8), ass. (0->231)
t219 = sin(qJ(2));
t222 = cos(qJ(2));
t216 = sin(pkin(8));
t217 = cos(pkin(8));
t260 = qJD(1) * t219;
t188 = -t217 * qJD(2) + t216 * t260;
t190 = t216 * qJD(2) + t217 * t260;
t218 = sin(qJ(4));
t221 = cos(qJ(4));
t160 = t221 * t188 + t218 * t190;
t208 = t219 * qJDD(1);
t257 = qJD(1) * qJD(2);
t252 = t222 * t257;
t194 = t208 + t252;
t169 = t216 * qJDD(2) + t217 * t194;
t246 = -t217 * qJDD(2) + t216 * t194;
t112 = -t160 * qJD(4) + t221 * t169 - t218 * t246;
t258 = t222 * qJD(1);
t203 = -qJD(4) + t258;
t277 = t160 * t203;
t308 = t112 + t277;
t205 = t219 * t257;
t256 = t222 * qJDD(1);
t195 = -t205 + t256;
t191 = -qJDD(4) + t195;
t162 = -t218 * t188 + t221 * t190;
t276 = t162 * t160;
t227 = t191 - t276;
t267 = t218 * t227;
t159 = t162 ^ 2;
t293 = t203 ^ 2;
t306 = -t159 - t293;
t63 = -t221 * t306 - t267;
t263 = t221 * t227;
t69 = -t218 * t306 + t263;
t41 = t216 * t69 - t217 * t63;
t43 = t216 * t63 + t217 * t69;
t355 = -pkin(6) * (t219 * t308 + t222 * t43) + pkin(1) * t41;
t353 = qJ(3) * t41;
t352 = -pkin(2) * t41 + pkin(3) * t63;
t351 = pkin(2) * t308 - qJ(3) * t43;
t295 = t160 ^ 2;
t139 = t295 - t293;
t79 = t218 * t139 - t263;
t83 = t221 * t139 + t267;
t148 = t162 * t203;
t247 = -t218 * t169 - t221 * t246;
t230 = t162 * qJD(4) - t247;
t89 = t148 + t230;
t349 = t219 * (t216 * t79 - t217 * t83) - t222 * t89;
t305 = t159 - t295;
t307 = -t148 + t230;
t50 = -t218 * t307 + t221 * t308;
t285 = t218 * t308;
t52 = t221 * t307 + t285;
t348 = t219 * (t216 * t50 + t217 * t52) + t222 * t305;
t346 = pkin(7) * t63;
t345 = pkin(7) * t69;
t344 = t216 * t52 - t217 * t50;
t342 = t216 * t83 + t217 * t79;
t115 = t191 + t276;
t266 = t218 * t115;
t303 = -t293 - t295;
t311 = t221 * t303 + t266;
t113 = t221 * t115;
t312 = t218 * t303 - t113;
t321 = t216 * t311 + t217 * t312;
t338 = qJ(3) * t321;
t337 = -pkin(2) * t321 - pkin(3) * t312;
t322 = -t216 * t312 + t217 * t311;
t336 = -pkin(2) * t307 + qJ(3) * t322;
t140 = -t159 + t293;
t323 = t221 * t140 - t266;
t324 = -t218 * t140 - t113;
t335 = t216 * t324 + t217 * t323;
t334 = pkin(6) * (t219 * t307 + t222 * t322) - pkin(1) * t321;
t304 = -t277 + t112;
t333 = t219 * (-t216 * t323 + t217 * t324) - t222 * t304;
t100 = -t295 - t159;
t332 = pkin(2) * t100;
t331 = pkin(3) * t100;
t329 = pkin(7) * t311;
t328 = pkin(7) * t312;
t326 = qJ(5) * t308;
t224 = qJD(1) ^ 2;
t220 = sin(qJ(1));
t223 = cos(qJ(1));
t240 = t223 * g(1) + t220 * g(2);
t278 = qJDD(1) * pkin(6);
t181 = -t224 * pkin(1) - t240 + t278;
t237 = -t222 * pkin(2) - t219 * qJ(3);
t245 = t224 * t237 + t181;
t289 = t222 * g(3);
t292 = qJD(2) ^ 2;
t135 = -qJDD(2) * pkin(2) - t292 * qJ(3) + t219 * t245 + qJDD(3) + t289;
t170 = -pkin(3) * t258 - t190 * pkin(7);
t294 = t188 ^ 2;
t98 = t246 * pkin(3) - t294 * pkin(7) + t190 * t170 + t135;
t327 = pkin(4) * t230 - t326 + t98;
t325 = t219 * t100;
t274 = t190 * t188;
t229 = -t195 - t274;
t310 = t216 * t229;
t309 = t217 * t229;
t174 = t188 * t258;
t152 = -t169 + t174;
t175 = t190 * t258;
t150 = -t246 - t175;
t123 = t160 * pkin(4) - t162 * qJ(5);
t250 = t220 * g(1) - t223 * g(2);
t180 = qJDD(1) * pkin(1) + t224 * pkin(6) + t250;
t236 = t194 + t252;
t131 = -t236 * qJ(3) + (-t195 + t205) * pkin(2) - t180;
t290 = t219 * g(3);
t136 = -t292 * pkin(2) + qJDD(2) * qJ(3) + t222 * t245 - t290;
t96 = 0.2e1 * qJD(3) * t190 - t217 * t131 + t216 * t136;
t60 = t229 * pkin(3) + t152 * pkin(7) - t96;
t97 = -0.2e1 * qJD(3) * t188 + t216 * t131 + t217 * t136;
t62 = -t294 * pkin(3) - pkin(7) * t246 + t170 * t258 + t97;
t34 = t218 * t60 + t221 * t62;
t249 = -t191 * qJ(5) - t160 * t123 + t34;
t300 = -pkin(4) * (t306 + t293) - qJ(5) * t227 + t249;
t228 = (t160 * t218 + t162 * t221) * t203;
t273 = t203 * t218;
t138 = t162 * t273;
t272 = t203 * t221;
t255 = t160 * t272;
t238 = -t138 + t255;
t299 = t216 * t238 + t217 * t228;
t232 = t218 * t230 - t255;
t239 = -t160 * t273 - t221 * t230;
t298 = t216 * t232 + t217 * t239;
t297 = t219 * (-t216 * t228 + t217 * t238) + t222 * t191;
t254 = t222 * t276;
t296 = t219 * (-t216 * t239 + t217 * t232) + t254;
t187 = t190 ^ 2;
t291 = pkin(4) * t221;
t33 = t218 * t62 - t221 * t60;
t18 = t218 * t34 - t221 * t33;
t288 = t216 * t18;
t287 = t217 * t18;
t284 = t218 * t304;
t283 = t218 * t98;
t281 = t221 * t98;
t279 = qJ(5) * t221;
t271 = t216 * t135;
t153 = t195 - t274;
t270 = t216 * t153;
t269 = t217 * t135;
t268 = t217 * t153;
t202 = t222 * t224 * t219;
t264 = t219 * (qJDD(2) + t202);
t262 = t222 * (-t202 + qJDD(2));
t259 = qJD(5) * t203;
t253 = t222 * t274;
t251 = -qJ(5) * t218 - pkin(3);
t57 = t216 * t96 + t217 * t97;
t19 = t218 * t33 + t221 * t34;
t166 = t219 * t181 + t289;
t167 = t222 * t181 - t290;
t248 = t219 * t166 + t222 * t167;
t75 = t218 * t112 - t162 * t272;
t76 = t221 * t112 + t138;
t244 = t219 * (-t216 * t75 + t217 * t76) - t254;
t199 = -0.2e1 * t259;
t243 = t199 + t249;
t24 = -pkin(4) * t293 + t243;
t25 = t191 * pkin(4) - qJ(5) * t293 + t162 * t123 + qJDD(5) + t33;
t242 = -pkin(4) * t25 + qJ(5) * t24;
t241 = -pkin(4) * t304 - qJ(5) * t89;
t235 = t216 * t97 - t217 * t96;
t234 = -pkin(1) + t237;
t226 = -pkin(4) * t115 + qJ(5) * t303 - t25;
t225 = 0.2e1 * qJD(5) * t162 - t327;
t214 = t222 ^ 2;
t213 = t219 ^ 2;
t210 = t214 * t224;
t209 = t213 * t224;
t196 = -0.2e1 * t205 + t256;
t193 = t208 + 0.2e1 * t252;
t183 = t222 * t195;
t173 = -t187 - t210;
t172 = -t187 + t210;
t171 = -t210 + t294;
t163 = -t210 - t294;
t151 = t169 + t174;
t149 = -t175 + t246;
t143 = -t187 - t294;
t128 = -t216 * t173 + t268;
t127 = t217 * t173 + t270;
t120 = t217 * t163 - t310;
t119 = t216 * t163 + t309;
t110 = t217 * t150 - t216 * t152;
t90 = (-qJD(4) - t203) * t162 + t247;
t85 = t221 * t304;
t55 = t221 * t90 + t284;
t53 = -t221 * t89 + t284;
t51 = t218 * t90 - t85;
t49 = -t218 * t89 - t85;
t47 = t281 + t346;
t46 = t216 * t76 + t217 * t75;
t40 = t283 - t328;
t35 = -pkin(3) * t308 + t283 + t345;
t31 = (-pkin(4) * t203 - 0.2e1 * qJD(5)) * t162 + t327;
t30 = -pkin(3) * t307 - t281 + t329;
t29 = -t216 * t51 + t217 * t55;
t28 = -t216 * t49 + t217 * t53;
t27 = t216 * t55 + t217 * t51;
t26 = t216 * t53 + t217 * t49;
t23 = t225 + (-t307 + t148) * pkin(4);
t22 = pkin(4) * t148 + t225 + t326;
t21 = -qJ(5) * t100 + t25;
t20 = (-t100 - t293) * pkin(4) + t243;
t17 = -t218 * t23 - t279 * t307 - t328;
t16 = -pkin(4) * t285 + t221 * t22 - t346;
t15 = -pkin(3) * t98 + pkin(7) * t19;
t14 = t221 * t23 + t251 * t307 + t329;
t13 = -t345 + t218 * t22 + (pkin(3) + t291) * t308;
t12 = -pkin(7) * t51 - t18;
t11 = t218 * t25 + t221 * t24;
t10 = t218 * t24 - t221 * t25;
t9 = pkin(7) * t55 + t19 - t331;
t8 = -pkin(7) * t49 - t218 * t20 + t221 * t21;
t7 = t217 * t19 - t288;
t6 = t216 * t19 + t287;
t5 = pkin(7) * t53 + t221 * t20 + t218 * t21 - t331;
t4 = -pkin(7) * t10 + (pkin(4) * t218 - t279) * t31;
t3 = -t216 * t10 + t217 * t11;
t2 = t217 * t10 + t216 * t11;
t1 = pkin(7) * t11 + (t251 - t291) * t31;
t32 = [0, 0, 0, 0, 0, qJDD(1), t250, t240, 0, 0, t236 * t219, t222 * t193 + t219 * t196, t264 + t222 * (-t209 + t292), -t219 * t252 + t183, t219 * (t210 - t292) + t262, 0, t222 * t180 + pkin(1) * t196 + pkin(6) * (t222 * (-t210 - t292) - t264), -t219 * t180 - pkin(1) * t193 + pkin(6) * (-t262 - t219 * (-t209 - t292)), pkin(1) * (t209 + t210) + (t213 + t214) * t278 + t248, pkin(1) * t180 + pkin(6) * t248, t219 * (t217 * t169 + t175 * t216) - t253, t219 * (-t217 * t149 - t216 * t151) + t222 * (-t187 + t294), t219 * (-t216 * t172 + t309) + t222 * t152, t219 * (-t174 * t217 + t216 * t246) + t253, t219 * (t217 * t171 + t270) - t222 * t150, t183 + t219 * (t188 * t217 - t190 * t216) * t258, t219 * (-qJ(3) * t119 + t271) + t222 * (-pkin(2) * t119 + t96) - pkin(1) * t119 + pkin(6) * (t222 * t120 + t219 * t149), t219 * (-qJ(3) * t127 + t269) + t222 * (-pkin(2) * t127 + t97) - pkin(1) * t127 + pkin(6) * (t222 * t128 + t219 * t151), -t219 * t235 + pkin(6) * (t222 * t110 + t219 * t143) + t234 * (t216 * t150 + t217 * t152), pkin(6) * (t219 * t135 + t222 * t57) + t234 * t235, t244, -t348, t333, t296, -t349, t297, t219 * (-t216 * t30 + t217 * t40 - t338) + t222 * (t33 + t337) + t334, t219 * (-t216 * t35 + t217 * t47 - t353) + t222 * (t34 + t352) - t355, t219 * (-qJ(3) * t27 + t217 * t12 - t216 * t9) + t222 * (-pkin(2) * t27 - pkin(3) * t51) - pkin(1) * t27 + pkin(6) * (t222 * t29 + t325), t219 * (-pkin(7) * t287 - qJ(3) * t6 - t216 * t15) + t222 * (-pkin(2) * t6 - pkin(3) * t18) - pkin(1) * t6 + pkin(6) * (t219 * t98 + t222 * t7), t244, t333, t348, t297, t349, t296, t219 * (-t216 * t14 + t217 * t17 - t338) + t222 * (-t226 + t337) + t334, t219 * (-qJ(3) * t26 - t216 * t5 + t217 * t8) + t222 * (-pkin(2) * t26 - pkin(3) * t49 - t241) - pkin(1) * t26 + pkin(6) * (t222 * t28 + t325), t219 * (-t216 * t13 + t217 * t16 + t353) + t222 * (0.2e1 * t259 - t300 - t352) + t355, t219 * (-qJ(3) * t2 - t216 * t1 + t217 * t4) + t222 * (-pkin(2) * t2 - pkin(3) * t10 - t242) - pkin(1) * t2 + pkin(6) * (t219 * t31 + t222 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, t209 - t210, t208, t202, t256, qJDD(2), -t166, -t167, 0, 0, t216 * t169 - t175 * t217, -t216 * t149 + t217 * t151, t217 * t172 + t310, -t174 * t216 - t217 * t246, t216 * t171 - t268, (t188 * t216 + t190 * t217) * t258, -pkin(2) * t149 + qJ(3) * t120 - t269, -pkin(2) * t151 + qJ(3) * t128 + t271, -pkin(2) * t143 + qJ(3) * t110 + t57, -pkin(2) * t135 + qJ(3) * t57, t46, -t344, t335, t298, t342, t299, t216 * t40 + t217 * t30 + t336, t216 * t47 + t217 * t35 - t351, qJ(3) * t29 + t216 * t12 + t217 * t9 - t332, -pkin(2) * t98 - pkin(7) * t288 + qJ(3) * t7 + t217 * t15, t46, t335, t344, t299, -t342, t298, t217 * t14 + t216 * t17 + t336, qJ(3) * t28 + t216 * t8 + t217 * t5 - t332, t217 * t13 + t216 * t16 + t351, -pkin(2) * t31 + qJ(3) * t3 + t217 * t1 + t216 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, t151, t143, t135, 0, 0, 0, 0, 0, 0, t307, t308, t100, t98, 0, 0, 0, 0, 0, 0, t307, t100, -t308, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t276, t305, t304, -t276, -t89, -t191, -t33, -t34, 0, 0, t276, t304, -t305, -t191, t89, -t276, t226, t241, t199 + t300, t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, t304, t306, t25;];
tauJ_reg = t32;
