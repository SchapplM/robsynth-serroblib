% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 06:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRRPP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:37:42
% EndTime: 2019-05-05 06:38:04
% DurationCPUTime: 8.93s
% Computational Cost: add. (20015->417), mult. (40294->550), div. (0->0), fcn. (29425->12), ass. (0->276)
t240 = sin(pkin(6));
t242 = cos(pkin(6));
t246 = sin(qJ(3));
t247 = sin(qJ(2));
t249 = cos(qJ(3));
t250 = cos(qJ(2));
t245 = sin(qJ(4));
t248 = cos(qJ(4));
t302 = qJD(2) * t246;
t209 = -t248 * qJD(3) + t245 * t302;
t296 = qJD(2) * qJD(3);
t289 = t249 * t296;
t295 = t246 * qJDD(2);
t213 = t289 + t295;
t270 = -t245 * qJDD(3) - t248 * t213;
t183 = -t209 * qJD(4) - t270;
t239 = sin(pkin(11));
t241 = cos(pkin(11));
t210 = t245 * qJD(3) + t248 * t302;
t271 = t248 * qJDD(3) - t245 * t213;
t260 = -t210 * qJD(4) + t271;
t256 = t241 * t183 + t239 * t260;
t188 = t241 * t209 + t239 * t210;
t228 = t249 * qJD(2) - qJD(4);
t326 = t188 * t228;
t351 = t256 + t326;
t231 = t246 * t296;
t294 = t249 * qJDD(2);
t214 = -t231 + t294;
t208 = -qJDD(4) + t214;
t190 = -t239 * t209 + t241 * t210;
t325 = t190 * t188;
t262 = t208 - t325;
t316 = t239 * t262;
t187 = t190 ^ 2;
t340 = t228 ^ 2;
t353 = -t187 - t340;
t77 = -t241 * t353 - t316;
t312 = t241 * t262;
t85 = -t239 * t353 + t312;
t58 = t245 * t77 + t248 * t85;
t43 = t246 * t351 + t249 * t58;
t56 = t245 * t85 - t248 * t77;
t429 = (t247 * t43 - t250 * t56) * t240 + t242 * (t246 * t58 - t249 * t351);
t428 = pkin(2) * t56 - pkin(8) * t43;
t424 = pkin(3) * t56;
t423 = pkin(9) * t56;
t421 = pkin(3) * t351 - pkin(9) * t58;
t169 = t190 * t228;
t285 = t239 * t183 - t241 * t260;
t106 = -t285 - t169;
t350 = -t326 + t256;
t373 = t241 * t106 + t239 * t350;
t374 = t239 * t106 - t241 * t350;
t387 = t245 * t373 + t248 * t374;
t342 = t188 ^ 2;
t115 = -t342 - t187;
t388 = -t245 * t374 + t248 * t373;
t404 = t246 * t115 + t249 * t388;
t419 = -pkin(2) * t387 + pkin(8) * t404;
t162 = t342 - t340;
t90 = t239 * t162 - t312;
t94 = t241 * t162 + t316;
t418 = t249 * t106 + t246 * (t245 * t90 - t248 * t94);
t352 = t187 - t342;
t355 = -t169 + t285;
t66 = -t239 * t355 + t241 * t351;
t318 = t239 * t351;
t68 = t241 * t355 + t318;
t417 = t246 * (t245 * t66 + t248 * t68) + t249 * t352;
t416 = t242 * (-t249 * t115 + t246 * t388) + (t247 * t404 - t250 * t387) * t240;
t415 = pkin(4) * t77;
t413 = pkin(3) * t387;
t412 = pkin(9) * t387;
t411 = qJ(5) * t77;
t410 = qJ(5) * t85;
t132 = t208 + t325;
t315 = t239 * t132;
t349 = -t340 - t342;
t358 = t241 * t349 + t315;
t127 = t241 * t132;
t359 = t239 * t349 - t127;
t371 = t245 * t358 + t248 * t359;
t372 = -t245 * t359 + t248 * t358;
t390 = t246 * t355 + t249 * t372;
t408 = -pkin(2) * t371 + pkin(8) * t390;
t407 = -pkin(3) * t115 + pkin(9) * t388;
t406 = t245 * t68 - t248 * t66;
t403 = t245 * t94 + t248 * t90;
t400 = t242 * (t246 * t372 - t249 * t355) + (t247 * t390 - t250 * t371) * t240;
t398 = pkin(3) * t371;
t337 = pkin(4) * t374;
t397 = pkin(9) * t371;
t394 = qJ(5) * t374;
t392 = -pkin(3) * t355 + pkin(9) * t372;
t391 = -pkin(4) * t115 + qJ(5) * t373;
t163 = -t187 + t340;
t375 = t241 * t163 - t315;
t376 = -t239 * t163 - t127;
t386 = t245 * t376 + t248 * t375;
t385 = t246 * (-t245 * t375 + t248 * t376) - t249 * t350;
t336 = pkin(4) * t359;
t379 = qJ(6) * t351;
t328 = sin(pkin(10));
t329 = cos(pkin(10));
t263 = g(1) * t328 - g(2) * t329;
t304 = -g(3) + qJDD(1);
t254 = -t240 * t263 + t242 * t304;
t195 = t249 * t254;
t218 = -g(1) * t329 - g(2) * t328;
t257 = t242 * t263;
t360 = t240 * t304 + t257;
t176 = t250 * t218 + t247 * t360;
t251 = qJD(2) ^ 2;
t160 = -t251 * pkin(2) + qJDD(2) * pkin(8) + t176;
t279 = -t249 * pkin(3) - t246 * pkin(9);
t282 = t251 * t279 + t160;
t339 = qJD(3) ^ 2;
t118 = -qJDD(3) * pkin(3) - t339 * pkin(9) + t282 * t246 - t195;
t196 = -t228 * pkin(4) - t210 * qJ(5);
t341 = t209 ^ 2;
t75 = -t260 * pkin(4) - t341 * qJ(5) + t210 * t196 + qJDD(5) + t118;
t382 = pkin(5) * t285 - t379 + t75;
t381 = qJ(5) * t358;
t380 = qJ(5) * t359;
t324 = t210 * t209;
t261 = -t208 - t324;
t357 = t245 * t261;
t356 = t248 * t261;
t301 = qJD(5) * t188;
t180 = -0.2e1 * t301;
t299 = qJD(6) * t228;
t354 = t180 - 0.2e1 * t299;
t199 = t209 * t228;
t154 = t183 - t199;
t146 = pkin(5) * t188 - qJ(6) * t190;
t253 = t246 * t254;
t119 = -pkin(3) * t339 + qJDD(3) * pkin(9) + t249 * t282 + t253;
t276 = t247 * t218 - t250 * t360;
t159 = -qJDD(2) * pkin(2) - t251 * pkin(8) + t276;
t274 = -t214 + t231;
t275 = t213 + t289;
t126 = pkin(3) * t274 - pkin(9) * t275 + t159;
t73 = t245 * t119 - t248 * t126;
t55 = t261 * pkin(4) - qJ(5) * t154 - t73;
t74 = t248 * t119 + t245 * t126;
t61 = -pkin(4) * t341 + qJ(5) * t260 + t228 * t196 + t74;
t334 = t239 * t55 + t241 * t61;
t283 = -t208 * qJ(6) - t188 * t146 + t334;
t347 = t415 - pkin(5) * (t353 + t340) - qJ(6) * t262 + t283;
t150 = (qJD(4) + t228) * t210 - t271;
t266 = (t188 * t239 + t190 * t241) * t228;
t323 = t228 * t239;
t161 = t190 * t323;
t322 = t228 * t241;
t293 = t188 * t322;
t277 = -t161 + t293;
t346 = t245 * t277 + t248 * t266;
t268 = t239 * t285 - t293;
t278 = -t188 * t323 - t241 * t285;
t345 = t245 * t268 + t248 * t278;
t200 = t249 * t208;
t344 = t246 * (-t245 * t266 + t248 * t277) + t200;
t292 = t249 * t325;
t343 = t246 * (-t245 * t278 + t248 * t268) + t292;
t207 = t210 ^ 2;
t287 = t239 * t61 - t241 * t55;
t300 = qJD(5) * t190;
t29 = t287 + 0.2e1 * t300;
t30 = t180 + t334;
t16 = t239 * t30 - t241 * t29;
t338 = pkin(4) * t16;
t335 = pkin(5) * t241;
t333 = t239 * t75;
t332 = t241 * t75;
t331 = t245 * t16;
t330 = t248 * t16;
t327 = qJ(6) * t241;
t321 = t228 * t245;
t320 = t228 * t248;
t311 = t245 * t118;
t172 = t208 - t324;
t310 = t245 * t172;
t227 = t246 * t251 * t249;
t220 = qJDD(3) + t227;
t309 = t246 * t220;
t308 = t248 * t118;
t307 = t248 * t172;
t219 = -t227 + qJDD(3);
t305 = t249 * t219;
t298 = qJD(4) - t228;
t291 = t249 * t324;
t288 = -qJ(6) * t239 - pkin(4);
t17 = t239 * t29 + t241 * t30;
t46 = t245 * t73 + t248 * t74;
t144 = t246 * t160 - t195;
t145 = t249 * t160 + t253;
t81 = t246 * t144 + t249 * t145;
t284 = (0.2e1 * qJD(5) + t146) * t190;
t281 = -t334 - t415;
t100 = -t190 * t322 + t239 * t256;
t101 = t241 * t256 + t161;
t280 = t246 * (-t245 * t100 + t248 * t101) - t292;
t45 = t245 * t74 - t248 * t73;
t272 = -pkin(2) + t279;
t269 = t283 + t354;
t267 = t208 * pkin(5) - qJ(6) * t340 + qJDD(6) + t287;
t24 = -pkin(5) * t340 + t269;
t25 = t284 + t267;
t12 = t239 * t24 - t241 * t25;
t265 = pkin(4) * t12 - pkin(5) * t25 + qJ(6) * t24;
t264 = -pkin(5) * t350 + qJ(6) * t106 + t337;
t255 = pkin(5) * t132 - qJ(6) * t349 + t267 - t336;
t252 = 0.2e1 * qJD(6) * t190 - t382;
t236 = t249 ^ 2;
t235 = t246 ^ 2;
t234 = t236 * t251;
t232 = t235 * t251;
t225 = -t234 - t339;
t224 = -t232 - t339;
t217 = t232 + t234;
t216 = (t235 + t236) * qJDD(2);
t215 = -0.2e1 * t231 + t294;
t212 = 0.2e1 * t289 + t295;
t198 = -t207 + t340;
t197 = -t340 + t341;
t194 = -t246 * t224 - t305;
t193 = t249 * t225 - t309;
t192 = t207 - t341;
t191 = -t207 - t340;
t184 = -t340 - t341;
t182 = -0.2e1 * t300;
t181 = 0.2e1 * t301;
t171 = t207 + t341;
t155 = t209 * t298 + t270;
t153 = t183 + t199;
t151 = -t210 * t298 + t271;
t143 = -t245 * t191 + t307;
t142 = t248 * t191 + t310;
t130 = t248 * t184 - t357;
t129 = t245 * t184 + t356;
t113 = -t150 * t248 + t245 * t154;
t112 = -t150 * t245 - t248 * t154;
t87 = t249 * t143 - t246 * t155;
t82 = t249 * t130 - t246 * t151;
t76 = t249 * t113 - t246 * t171;
t63 = t248 * t100 + t245 * t101;
t52 = t332 + t411;
t47 = t333 - t380;
t42 = -pkin(4) * t351 + t333 + t410;
t39 = (-pkin(5) * t228 - 0.2e1 * qJD(6)) * t190 + t382;
t38 = -pkin(4) * t355 - t332 + t381;
t37 = t246 * t118 + t249 * t46;
t32 = t252 + (-t355 + t169) * pkin(5);
t31 = pkin(5) * t169 + t252 + t379;
t23 = -qJ(6) * t115 + t25;
t22 = (-t115 - t340) * pkin(5) + t269;
t21 = -t239 * t32 - t327 * t355 - t380;
t20 = -pkin(5) * t318 + t241 * t31 - t411;
t19 = t241 * t32 + t288 * t355 + t381;
t18 = -t410 + t239 * t31 + (pkin(4) + t335) * t351;
t15 = -pkin(4) * t75 + qJ(5) * t17;
t14 = -t16 - t394;
t13 = t239 * t25 + t241 * t24;
t11 = t17 + t391;
t10 = -t239 * t22 + t241 * t23 - t394;
t9 = t241 * t22 + t239 * t23 + t391;
t8 = -qJ(5) * t12 + (pkin(5) * t239 - t327) * t39;
t7 = t248 * t17 - t331;
t6 = t245 * t17 + t330;
t5 = t246 * t75 + t249 * t7;
t4 = qJ(5) * t13 + (t288 - t335) * t39;
t3 = -t245 * t12 + t248 * t13;
t2 = t248 * t12 + t245 * t13;
t1 = t246 * t39 + t249 * t3;
t26 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t304, 0, 0, 0, 0, 0, 0, (qJDD(2) * t250 - t247 * t251) * t240, (-qJDD(2) * t247 - t250 * t251) * t240, 0, t242 ^ 2 * t304 + (t247 * t176 - t250 * t276 - t257) * t240, 0, 0, 0, 0, 0, 0, t242 * (t249 * t220 + t246 * t225) + (t247 * t193 + t250 * t215) * t240, t242 * (-t246 * t219 + t249 * t224) + (t247 * t194 - t250 * t212) * t240, (t216 * t247 + t217 * t250) * t240, t242 * (-t249 * t144 + t246 * t145) + (-t250 * t159 + t247 * t81) * t240, 0, 0, 0, 0, 0, 0, t242 * (t246 * t130 + t249 * t151) + (-t250 * t129 + t247 * t82) * t240, t242 * (t246 * t143 + t249 * t155) + (-t250 * t142 + t247 * t87) * t240, t242 * (t246 * t113 + t249 * t171) + (-t250 * t112 + t247 * t76) * t240, t242 * (-t249 * t118 + t246 * t46) + (t247 * t37 - t250 * t45) * t240, 0, 0, 0, 0, 0, 0, t400, t429, t416, t242 * (t246 * t7 - t249 * t75) + (t247 * t5 - t250 * t6) * t240, 0, 0, 0, 0, 0, 0, t400, t416, -t429, t242 * (t246 * t3 - t249 * t39) + (t247 * t1 - t250 * t2) * t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t276, -t176, 0, 0, t275 * t246, t249 * t212 + t246 * t215, t309 + t249 * (-t232 + t339), -t274 * t249, t246 * (t234 - t339) + t305, 0, pkin(2) * t215 + pkin(8) * t193 - t249 * t159, -pkin(2) * t212 + pkin(8) * t194 + t246 * t159, pkin(2) * t217 + pkin(8) * t216 + t81, -pkin(2) * t159 + pkin(8) * t81, t246 * (t248 * t183 + t210 * t321) - t291, t246 * (t248 * t151 - t245 * t153) - t249 * t192, t246 * (-t245 * t198 + t356) - t249 * t154, t246 * (-t209 * t320 - t245 * t260) + t291, t246 * (t248 * t197 + t310) + t249 * t150, t200 + t246 * (t209 * t248 - t210 * t245) * t228, t246 * (-pkin(9) * t129 + t311) + t249 * (-pkin(3) * t129 + t73) - pkin(2) * t129 + pkin(8) * t82, t246 * (-pkin(9) * t142 + t308) + t249 * (-pkin(3) * t142 + t74) - pkin(2) * t142 + pkin(8) * t87, pkin(8) * t76 + t112 * t272 - t246 * t45, pkin(8) * t37 + t272 * t45, t280, -t417, t385, t343, -t418, t344, t246 * (-t245 * t38 + t248 * t47 - t397) + t249 * (t29 - t336 - t398) + t408, t246 * (-t245 * t42 + t248 * t52 - t423) + t249 * (t180 - t281 - t424) - t428, t246 * (-t245 * t11 + t248 * t14 - t412) + t249 * (-t337 - t413) + t419, t246 * (-pkin(9) * t6 - qJ(5) * t330 - t245 * t15) + t249 * (-pkin(3) * t6 - t338) - pkin(2) * t6 + pkin(8) * t5, t280, t385, t417, t344, t418, t343, t246 * (-t245 * t19 + t248 * t21 - t397) + t249 * (t255 + t284 - t398) + t408, t246 * (t248 * t10 - t245 * t9 - t412) + t249 * (-t264 - t413) + t419, t246 * (-t245 * t18 + t248 * t20 + t423) + t249 * (t181 + 0.2e1 * t299 - t347 + t424) + t428, t246 * (-pkin(9) * t2 - t245 * t4 + t248 * t8) + t249 * (-pkin(3) * t2 - t265) - pkin(2) * t2 + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t227, t232 - t234, t295, t227, t294, qJDD(3), -t144, -t145, 0, 0, t245 * t183 - t210 * t320, t245 * t151 + t248 * t153, t248 * t198 + t357, -t209 * t321 + t248 * t260, t245 * t197 - t307, (t209 * t245 + t210 * t248) * t228, pkin(3) * t151 + pkin(9) * t130 - t308, pkin(3) * t155 + pkin(9) * t143 + t311, pkin(3) * t171 + pkin(9) * t113 + t46, -pkin(3) * t118 + pkin(9) * t46, t63, -t406, t386, t345, t403, t346, t245 * t47 + t248 * t38 + t392, t245 * t52 + t248 * t42 - t421, t248 * t11 + t245 * t14 + t407, -pkin(3) * t75 + pkin(9) * t7 - qJ(5) * t331 + t248 * t15, t63, t386, t406, t346, -t403, t345, t248 * t19 + t245 * t21 + t392, t245 * t10 + t248 * t9 + t407, t248 * t18 + t245 * t20 + t421, -pkin(3) * t39 + pkin(9) * t3 + t245 * t8 + t248 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, t192, t154, -t324, -t150, -t208, -t73, -t74, 0, 0, t325, t352, t350, -t325, t106, -t208, t182 - t287 + t336, t181 + t281, t337, t338, t325, t350, -t352, -t208, -t106, -t325, -t190 * t146 + t182 - t255, t264, t347 + t354, t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t355, t351, t115, t75, 0, 0, 0, 0, 0, 0, t355, t115, -t351, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t350, t353, t25;];
tauJ_reg  = t26;
