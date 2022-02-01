% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x24]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:43:18
% EndTime: 2022-01-20 11:43:26
% DurationCPUTime: 4.07s
% Computational Cost: add. (3912->279), mult. (7856->343), div. (0->0), fcn. (8424->8), ass. (0->231)
t398 = qJD(3) + qJD(5);
t237 = sin(qJ(3));
t352 = -qJ(4) - pkin(7);
t216 = t352 * t237;
t240 = cos(qJ(3));
t217 = t352 * t240;
t234 = sin(pkin(9));
t235 = cos(pkin(9));
t149 = t234 * t216 - t235 * t217;
t319 = t235 * t240;
t324 = t234 * t237;
t200 = -t319 + t324;
t196 = t200 * pkin(8);
t115 = t196 - t149;
t236 = sin(qJ(5));
t239 = cos(qJ(5));
t320 = t235 * t237;
t323 = t234 * t240;
t202 = t320 + t323;
t197 = t202 * pkin(8);
t271 = t235 * t216 + t217 * t234;
t370 = -t197 + t271;
t412 = t398 * (t239 * t115 - t236 * t370);
t411 = t398 * (-t236 * t115 - t239 * t370);
t238 = sin(qJ(2));
t354 = t238 * pkin(1);
t226 = pkin(7) + t354;
t305 = qJ(4) + t226;
t198 = t305 * t237;
t199 = t305 * t240;
t272 = -t235 * t198 - t199 * t234;
t369 = -t197 + t272;
t134 = -t234 * t198 + t235 * t199;
t92 = -t196 + t134;
t404 = t398 * (-t236 * t369 - t239 * t92);
t403 = t398 * (t236 * t92 - t239 * t369);
t367 = t134 / 0.2e1;
t364 = t149 / 0.2e1;
t292 = qJD(1) + qJD(2);
t140 = t239 * t200 + t236 * t202;
t189 = t239 * t202;
t314 = t236 * t200;
t371 = t189 - t314;
t385 = t140 ^ 2 - t371 ^ 2;
t395 = t292 * t385;
t241 = cos(qJ(2));
t353 = t241 * pkin(1);
t167 = t200 * t353;
t358 = -t236 / 0.2e1;
t166 = t202 * t353;
t360 = -t166 / 0.2e1;
t249 = -t167 * t358 + t239 * t360;
t228 = -t240 * pkin(3) - pkin(2);
t165 = t200 * pkin(4) + t228;
t158 = t165 - t353;
t275 = t165 / 0.2e1 + t158 / 0.2e1;
t379 = t275 * t371;
t25 = t249 - t379;
t390 = t292 * t140;
t376 = t292 * t371;
t387 = t398 * t140;
t386 = t140 * qJD(4);
t250 = -t166 * t358 + t239 * t167 / 0.2e1;
t26 = t275 * t140 + t250;
t384 = -t234 / 0.2e1;
t280 = t189 / 0.2e1;
t88 = 0.2e1 * t280 - t314;
t380 = t292 * t88;
t251 = t200 * t384 - t235 * t202 / 0.2e1;
t357 = -t237 / 0.2e1;
t106 = (t357 + t251) * pkin(3);
t378 = t292 * t106;
t138 = t280 - t189 / 0.2e1;
t377 = t292 * t138;
t145 = t200 ^ 2 + t202 ^ 2;
t375 = t292 * t145;
t374 = t292 * t200;
t373 = t292 * t202;
t222 = -t237 ^ 2 + t240 ^ 2;
t372 = t292 * t222;
t227 = -pkin(2) - t353;
t287 = -t353 / 0.2e1;
t368 = t287 + t227 / 0.2e1;
t366 = -t272 / 0.2e1;
t363 = -t271 / 0.2e1;
t356 = pkin(3) * t234;
t355 = t237 * pkin(3);
t351 = pkin(1) * qJD(1);
t350 = pkin(1) * qJD(2);
t349 = pkin(2) * qJD(2);
t348 = qJD(3) * pkin(3);
t144 = t145 * qJD(4);
t54 = t166 * t202 + t167 * t200;
t339 = t54 * qJD(2) + t144;
t335 = t134 * t200;
t336 = t272 * t202;
t42 = -t335 - t336;
t338 = qJD(1) * t42;
t337 = qJD(1) * t54;
t334 = t271 * t202;
t333 = t149 * t200;
t332 = t158 * t371;
t331 = t158 * t140;
t330 = t165 * t371;
t329 = t165 * t140;
t215 = t228 - t353;
t328 = t215 * t200;
t327 = t215 * t202;
t326 = t228 * t200;
t325 = t228 * t202;
t30 = t215 * t355;
t309 = t30 * qJD(1);
t35 = -t134 * t167 - t166 * t272 + t215 * t354;
t308 = t35 * qJD(1);
t168 = pkin(4) * t202 + t355;
t77 = t168 * t140;
t37 = t77 + t332;
t307 = t37 * qJD(1);
t78 = t168 * t371;
t38 = t78 - t331;
t306 = t38 * qJD(1);
t184 = t200 * t355;
t113 = t184 + t327;
t303 = qJD(1) * t113;
t185 = t202 * t355;
t114 = t185 - t328;
t302 = qJD(1) * t114;
t301 = qJD(1) * t158;
t300 = qJD(1) * t227;
t299 = qJD(2) * t165;
t296 = qJD(5) * t158;
t295 = qJD(5) * t165;
t294 = t138 * qJD(5);
t293 = t237 * qJD(3);
t233 = t240 * qJD(3);
t291 = t238 * t350;
t290 = t238 * t351;
t289 = t355 / 0.2e1;
t288 = t354 / 0.2e1;
t284 = t237 * t300;
t283 = t240 * t300;
t282 = t140 * t301;
t281 = t371 * t301;
t277 = t271 / 0.2e1 + t272 / 0.2e1;
t274 = t228 / 0.2e1 + t215 / 0.2e1;
t273 = pkin(1) * t292;
t269 = t140 * t290;
t268 = t371 * t290;
t267 = t200 * t290;
t266 = t202 * t290;
t265 = t237 * t290;
t264 = t238 * t273;
t2 = (t363 + t366 + t277) * t200;
t263 = qJD(2) * t2;
t262 = qJD(1) * t2;
t242 = t134 * t363 + t149 * t366 + t271 * t367 + t272 * t364;
t252 = t167 * t384 + t235 * t360;
t3 = (-t274 * t237 + t252) * pkin(3) + t242;
t36 = t228 * t355;
t261 = -t3 * qJD(1) + t36 * qJD(2);
t12 = -t77 + t25;
t39 = t77 + t330;
t260 = -t12 * qJD(1) + t39 * qJD(2);
t13 = t26 - t78;
t40 = t78 - t329;
t259 = -t13 * qJD(1) + t40 * qJD(2);
t19 = (t364 + t367) * t200 + t277 * t202 + t288;
t52 = -t333 - t334;
t258 = -qJD(1) * t19 + qJD(2) * t52;
t122 = t184 + t325;
t244 = (-t323 / 0.2e1 - t320 / 0.2e1) * t353;
t48 = -t274 * t202 - t184 + t244;
t257 = qJD(1) * t48 - qJD(2) * t122;
t123 = t185 - t326;
t243 = (-t319 / 0.2e1 + t324 / 0.2e1) * t353;
t49 = t274 * t200 - t185 + t243;
t256 = qJD(1) * t49 - qJD(2) * t123;
t45 = qJD(3) * t371 + qJD(5) * t88;
t255 = t287 + pkin(2) / 0.2e1 - t227 / 0.2e1;
t159 = t255 * t237;
t254 = qJD(1) * t159 + t237 * t349;
t160 = t255 * t240;
t253 = qJD(1) * t160 + t240 * t349;
t248 = -t26 * qJD(1) - t140 * t299;
t247 = -t25 * qJD(1) + t299 * t371;
t225 = pkin(3) * t235 + pkin(4);
t186 = -t225 * t239 + t236 * t356;
t246 = qJD(3) * t186;
t187 = t225 * t236 + t239 * t356;
t245 = qJD(3) * t187;
t223 = t237 * t233;
t221 = t237 * t291;
t218 = t222 * qJD(3);
t195 = t202 * qJD(4);
t194 = t202 * qJD(3);
t193 = t200 * qJD(3);
t192 = t200 * qJD(4);
t191 = t292 * t240 * t237;
t176 = t202 * t291;
t175 = t200 * t291;
t172 = t187 * qJD(5);
t171 = t186 * qJD(5);
t162 = (-pkin(2) / 0.2e1 + t368) * t240;
t161 = pkin(2) * t357 + t368 * t237;
t135 = t371 * qJD(4);
t129 = t138 * qJD(3);
t128 = t138 * qJD(4);
t121 = (t200 * t235 - t202 * t234) * t348;
t120 = t371 * t291;
t119 = t140 * t291;
t105 = t251 * pkin(3) + t289;
t104 = t106 * qJD(3);
t103 = t106 * qJD(4);
t102 = t105 * qJD(3);
t101 = t105 * qJD(4);
t80 = t88 * qJD(4);
t51 = t185 - t326 / 0.2e1 - t328 / 0.2e1 + t243;
t50 = t184 + t325 / 0.2e1 + t327 / 0.2e1 + t244;
t44 = -qJD(3) * t88 - qJD(5) * t371;
t33 = t371 * t390;
t32 = t387 * t371;
t31 = t140 * t376;
t28 = t249 + t379;
t27 = t250 - (t158 + t165) * t140 / 0.2e1;
t20 = -t333 / 0.2e1 - t335 / 0.2e1 - t334 / 0.2e1 - t336 / 0.2e1 + t288;
t15 = t78 - t329 / 0.2e1 - t331 / 0.2e1 + t250;
t14 = t77 + t330 / 0.2e1 + t332 / 0.2e1 + t249;
t7 = t398 * t385;
t4 = t252 * pkin(3) - t242 + (t215 + t228) * t289;
t1 = t2 * qJD(3);
t5 = [0, 0, 0, 0, -t291, -t241 * t350, t223, t218, 0, 0, 0, t227 * t293 - t240 * t291, t227 * t233 + t221, qJD(3) * t113 + t175, qJD(3) * t114 + t176, t339, qJD(2) * t35 + qJD(3) * t30 + qJD(4) * t42, -t32, t7, 0, 0, 0, qJD(3) * t37 + t296 * t371 + t119, qJD(3) * t38 - t140 * t296 + t120; 0, 0, 0, 0, -t264, -t241 * t273, t223, t218, 0, 0, 0, t161 * qJD(3) - t240 * t264, qJD(3) * t162 + t221 + t265, qJD(3) * t50 + t175 + t267, qJD(3) * t51 + t176 + t266, t1 + t337 + t339, t308 + (-t149 * t167 - t166 * t271 + t228 * t354) * qJD(2) + t4 * qJD(3) + t20 * qJD(4), -t32, t7, 0, 0, 0, qJD(3) * t14 + qJD(5) * t28 + t119 + t269, qJD(3) * t15 + qJD(5) * t27 + t120 + t268; 0, 0, 0, 0, 0, 0, t191, t372, t233, -t293, 0, qJD(2) * t161 - t226 * t233 + t284, qJD(2) * t162 + t226 * t293 + t283, qJD(2) * t50 - qJD(3) * t134 + t303, qJD(2) * t51 - qJD(3) * t272 + t302, t121 + t263, t309 + t4 * qJD(2) + t101 + (-t134 * t235 + t234 * t272) * t348, -t33, t395, -t387, -t45, 0, t14 * qJD(2) + t307 + t404, t15 * qJD(2) + t306 + t403; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t375, qJD(2) * t20 + t102 + t338, 0, 0, 0, 0, 0, t294, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t395, -t387, t44, 0, t28 * qJD(2) + t128 + t281 + t404, t27 * qJD(2) - t282 + t403; 0, 0, 0, 0, t290, t241 * t351, t223, t218, 0, 0, 0, -qJD(3) * t159 + t240 * t290, -qJD(3) * t160 - t265, -qJD(3) * t48 - t267, -qJD(3) * t49 - t266, t1 + t144 - t337, -qJD(3) * t3 - qJD(4) * t19 - t308, -t32, t7, 0, 0, 0, -qJD(3) * t12 - qJD(5) * t25 - t269, -qJD(3) * t13 - qJD(5) * t26 - t268; 0, 0, 0, 0, 0, 0, t223, t218, 0, 0, 0, -pkin(2) * t293, -pkin(2) * t233, t122 * qJD(3), t123 * qJD(3), t144, qJD(3) * t36 + qJD(4) * t52, -t32, t7, 0, 0, 0, qJD(3) * t39 + t295 * t371, qJD(3) * t40 - t140 * t295; 0, 0, 0, 0, 0, 0, t191, t372, t233, -t293, 0, -pkin(7) * t233 - t254, pkin(7) * t293 - t253, -qJD(3) * t149 - t257, -qJD(3) * t271 - t256, t121 + t262, t101 + (-t149 * t235 + t234 * t271) * t348 + t261, -t33, t395, -t387, -t45, 0, t260 + t412, t259 + t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t375, t102 + t258, 0, 0, 0, 0, 0, t294, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t395, -t387, t44, 0, t128 + t247 + t412, t248 + t411; 0, 0, 0, 0, 0, 0, -t191, -t372, 0, 0, 0, qJD(2) * t159 - t284, qJD(2) * t160 - t283, qJD(2) * t48 - t195 - t303, qJD(2) * t49 + t192 - t302, -t263, qJD(2) * t3 + t103 - t309, t33, -t395, 0, -t294, 0, qJD(2) * t12 - t135 - t307, qJD(2) * t13 - t306 + t386; 0, 0, 0, 0, 0, 0, -t191, -t372, 0, 0, 0, t254, t253, -t195 + t257, t192 + t256, -t262, t103 - t261, t33, -t395, 0, -t294, 0, -t135 - t260, -t259 + t386; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t373, t374, 0, t378, 0, 0, 0, 0, 0, -t376, t390; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t377, 0, -t172 - t245, t171 + t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, -t193, -t375, qJD(2) * t19 - t104 - t338, 0, 0, 0, 0, 0, t45, -t387; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, -t193, -t375, -t104 - t258, 0, 0, 0, 0, 0, t45, -t387; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t373, -t374, 0, -t378, 0, 0, 0, 0, 0, t376, -t390; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t380, -t390; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t395, 0, t129, 0, qJD(2) * t25 - t281 - t80, qJD(2) * t26 + t282 + t386; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t395, 0, t129, 0, -t247 - t80, -t248 + t386; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t377, 0, t245, -t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t380, t390; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
