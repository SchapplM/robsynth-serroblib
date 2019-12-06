% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PRPRP6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PRPRP6_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP6_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:21
% EndTime: 2019-12-05 15:41:30
% DurationCPUTime: 3.64s
% Computational Cost: add. (4244->320), mult. (8153->404), div. (0->0), fcn. (4491->6), ass. (0->230)
t318 = cos(qJ(4));
t320 = qJD(4) ^ 2;
t308 = t318 ^ 2;
t321 = qJD(2) ^ 2;
t368 = t308 * t321;
t334 = t320 + t368;
t316 = sin(qJ(4));
t343 = t316 * t318 * t321;
t290 = qJDD(4) + t343;
t360 = t316 * t290;
t239 = t318 * t334 + t360;
t352 = qJD(2) * qJD(4);
t340 = t316 * t352;
t347 = qJDD(2) * t318;
t279 = -0.2e1 * t340 + t347;
t317 = sin(qJ(2));
t319 = cos(qJ(2));
t193 = t239 * t317 - t279 * t319;
t371 = t290 * t318;
t248 = t316 * t334 - t371;
t312 = sin(pkin(7));
t313 = cos(pkin(7));
t163 = t193 * t312 + t248 * t313;
t419 = qJ(1) * t163;
t166 = t193 * t313 - t248 * t312;
t418 = qJ(1) * t166;
t189 = t239 * t319 + t279 * t317;
t417 = pkin(5) * t189;
t393 = pkin(6) + pkin(2);
t416 = pkin(1) * t189 + t393 * t239;
t415 = pkin(1) * t248 + pkin(5) * t193;
t339 = t318 * t352;
t349 = qJDD(2) * t316;
t276 = 0.2e1 * t339 + t349;
t215 = t276 * t316 - t279 * t318;
t307 = t316 ^ 2;
t287 = (-t307 + t308) * t321;
t183 = t215 * t317 - t287 * t319;
t217 = t276 * t318 + t279 * t316;
t414 = t183 * t312 - t217 * t313;
t413 = t183 * t313 + t217 * t312;
t369 = t307 * t321;
t292 = -t320 + t369;
t238 = t292 * t316 + t371;
t345 = t319 * qJDD(2);
t198 = t238 * t317 - t316 * t345;
t242 = -t292 * t318 + t360;
t412 = t198 * t312 + t242 * t313;
t411 = t198 * t313 - t242 * t312;
t348 = qJDD(2) * t317;
t284 = t319 * t321 + t348;
t288 = g(1) * t312 - t313 * g(2);
t230 = -pkin(5) * t284 + t319 * t288;
t407 = t312 * t230;
t406 = t313 * t230;
t404 = t393 * t248;
t403 = -pkin(3) * t239 - qJ(3) * t248;
t285 = t317 * t321 - t345;
t229 = pkin(5) * t285 - t288 * t317;
t398 = t312 * t229;
t397 = t313 * t229;
t396 = t215 * t319 + t287 * t317;
t395 = t238 * t319 + t316 * t348;
t289 = g(1) * t313 + g(2) * t312;
t309 = g(3) - qJDD(1);
t329 = t319 * t289 + t317 * t309;
t336 = pkin(1) * t284 + qJ(1) * t285 - t329;
t254 = -t289 * t317 + t319 * t309;
t335 = -pkin(1) * t285 + qJ(1) * t284 - t254;
t293 = -t320 - t369;
t291 = qJDD(4) - t343;
t370 = t291 * t318;
t237 = t293 * t316 + t370;
t273 = t316 * t291;
t243 = t293 * t318 - t273;
t394 = pkin(3) * t237 - qJ(3) * t243;
t265 = t313 * t288;
t223 = -t289 * t312 + t265;
t355 = t307 + t308;
t286 = t355 * t321;
t392 = pkin(3) * t286;
t391 = pkin(4) * t316;
t390 = pkin(4) * t318;
t188 = -t237 * t319 + t276 * t317;
t389 = pkin(5) * t188;
t283 = t355 * qJDD(2);
t376 = t283 * t319;
t221 = -t286 * t317 + t376;
t388 = pkin(5) * t221;
t191 = t237 * t317 + t276 * t319;
t161 = t191 * t312 - t243 * t313;
t386 = qJ(1) * t161;
t377 = t283 * t317;
t222 = -t286 * t319 - t377;
t385 = qJ(1) * t222;
t381 = qJ(5) * t316;
t380 = qJ(5) * t318;
t373 = t288 * t312;
t367 = t312 * t284;
t366 = t312 * t285;
t365 = t312 * t309;
t212 = t313 * t222;
t364 = t313 * t284;
t363 = t313 * t285;
t362 = t313 * t309;
t350 = qJDD(2) * qJ(3);
t325 = -t321 * pkin(2) - t329 + t350;
t351 = (qJD(3) * qJD(2));
t218 = t325 + (2 * t351);
t314 = t321 * pkin(6);
t207 = t218 - t314;
t361 = t316 * t207;
t359 = t318 * t207;
t357 = -pkin(1) * t243 + pkin(5) * t191;
t311 = qJDD(2) * pkin(2);
t220 = -t321 * qJ(3) + qJDD(3) + t254 - t311;
t213 = -qJDD(2) * pkin(6) + t220;
t180 = -t318 * t213 - t316 * t288;
t356 = -t316 * t213 + t318 * t288;
t330 = -t380 + t391;
t354 = t321 * t330;
t353 = qJD(5) * t318;
t346 = t313 * qJDD(2);
t344 = t393 * t243;
t303 = -2 * t351;
t338 = -t303 + 0.2e1 * t350 + t336;
t337 = -qJDD(3) + 0.2e1 * t311 + t335;
t176 = t319 * t218 + t220 * t317;
t186 = t254 * t317 - t319 * t329;
t224 = -t313 * t289 - t373;
t333 = t317 * t343;
t332 = t319 * t343;
t331 = t381 + t390;
t145 = -t180 * t318 - t316 * t356;
t146 = t316 * t180 - t318 * t356;
t174 = t218 * t317 - t220 * t319;
t185 = t254 * t319 + t317 * t329;
t328 = -pkin(1) * t188 + t393 * t237;
t327 = qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t316 * t354 - t356;
t326 = -pkin(1) * t221 + qJ(3) * t286 - t393 * t283;
t324 = -qJDD(4) * pkin(4) - t320 * qJ(5) + t318 * t354 + qJDD(5) + t180;
t277 = -t339 - t349;
t278 = -t340 + t347;
t323 = -t277 * pkin(4) - t278 * qJ(5) - t314 + t325;
t322 = 0.2e1 * qJD(2) * t353 + t303 - t323;
t299 = t312 * qJDD(2);
t294 = t320 - t368;
t280 = pkin(1) * t288;
t272 = t355 * t352;
t253 = qJDD(4) * t319 - t272 * t317;
t252 = qJDD(4) * t317 + t272 * t319;
t251 = -t278 * t316 - t308 * t352;
t250 = -t277 * t318 - t307 * t352;
t247 = t294 * t316 - t370;
t244 = (-t278 + t340) * t318;
t240 = -t294 * t318 - t273;
t235 = (t277 - t339) * t316;
t234 = -pkin(3) * t283 - t331 * qJDD(2);
t233 = t313 * t253;
t232 = t312 * t253;
t219 = pkin(5) * t222;
t211 = t312 * t222;
t208 = qJ(1) * t212;
t203 = -t250 * t317 - t332;
t202 = -t251 * t317 + t332;
t201 = t250 * t319 - t333;
t200 = t251 * t319 + t333;
t199 = -t240 * t317 + t318 * t345;
t196 = t240 * t319 + t317 * t347;
t178 = t186 * t313 - t373;
t177 = t186 * t312 + t265;
t173 = t203 * t313 - t235 * t312;
t172 = t202 * t313 - t244 * t312;
t171 = t203 * t312 + t235 * t313;
t170 = t202 * t312 + t244 * t313;
t169 = t199 * t313 - t247 * t312;
t168 = t199 * t312 + t247 * t313;
t167 = -t320 * pkin(4) + t327;
t164 = t191 * t313 + t243 * t312;
t159 = qJ(1) * t164;
t158 = (t331 * qJD(4) + (2 * qJD(3)) - 0.2e1 * t353) * qJD(2) + t323;
t157 = t176 * t313 - t373;
t156 = t176 * t312 + t265;
t155 = qJ(5) * t286 + t324;
t154 = (t286 - t320) * pkin(4) + t327;
t153 = -qJ(5) * t340 + (-t276 - t339) * pkin(4) + t322;
t152 = -pkin(4) * t339 + (t279 - t340) * qJ(5) + t322;
t151 = -pkin(5) * t174 + (-pkin(2) * t317 + qJ(3) * t319) * t288;
t150 = t356 + t403;
t149 = -t180 + t394;
t148 = pkin(3) * t276 - t344 + t359;
t147 = pkin(3) * t279 - t361 - t404;
t144 = -t146 - t392;
t143 = -pkin(1) * t174 + pkin(2) * t220 - qJ(3) * t218;
t142 = t145 * t317 + t207 * t319;
t141 = -t145 * t319 + t207 * t317;
t140 = pkin(4) * t291 + qJ(5) * t293 - t324 + t394;
t139 = -qJ(3) * t279 - t359 - t416;
t138 = -qJ(3) * t276 + t328 - t361;
t137 = pkin(4) * t368 + qJ(5) * t290 + t327 - t403;
t136 = t167 * t318 + t316 * t324;
t135 = t316 * t167 - t318 * t324;
t134 = -pkin(3) * t376 - t144 * t317 - t388;
t133 = t145 + t326;
t132 = -t316 * t152 + (-pkin(3) - t390) * t279 + t404;
t131 = -t318 * t153 + (pkin(3) + t381) * t276 - t344;
t130 = -t154 * t318 - t155 * t316 - t392;
t129 = pkin(3) * t145 - qJ(3) * t146;
t128 = -t318 * t152 + (qJ(3) + t391) * t279 + t416;
t127 = t316 * t153 + (-qJ(3) + t380) * t276 + t328;
t126 = t135 * t317 + t158 * t319;
t125 = -t135 * t319 + t158 * t317;
t124 = t154 * t316 - t155 * t318 + t326;
t123 = pkin(3) * t207 - t393 * t146;
t122 = -t130 * t317 + t234 * t319 - t388;
t121 = -t147 * t317 + t150 * t319 - t417;
t120 = -t148 * t317 + t149 * t319 - t389;
t119 = t142 * t313 + t146 * t312;
t118 = t142 * t312 - t146 * t313;
t117 = -t131 * t317 + t140 * t319 - t389;
t116 = -t132 * t317 + t137 * t319 + t417;
t115 = -pkin(1) * t141 - qJ(3) * t207 + t393 * t145;
t114 = t126 * t313 + t136 * t312;
t113 = t126 * t312 - t136 * t313;
t112 = pkin(3) * t135 - pkin(4) * t324 - qJ(3) * t136 + qJ(5) * t167;
t111 = -t393 * t136 + (pkin(3) + t331) * t158;
t110 = -pkin(5) * t141 - t123 * t317 + t129 * t319;
t109 = -pkin(1) * t125 + t393 * t135 + (-qJ(3) - t330) * t158;
t108 = -pkin(5) * t125 - t111 * t317 + t112 * t319;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, 0, 0, 0, 0, 0, 0, -t364, t363, 0, t178, 0, 0, 0, 0, 0, 0, 0, t364, -t363, t157, 0, 0, 0, 0, 0, 0, t164, -t166, t212, t119, 0, 0, 0, 0, 0, 0, t164, t212, t166, t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, 0, 0, 0, 0, 0, 0, -t367, t366, 0, t177, 0, 0, 0, 0, 0, 0, 0, t367, -t366, t156, 0, 0, 0, 0, 0, 0, t161, -t163, t211, t118, 0, 0, 0, 0, 0, 0, t161, t211, t163, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t309, 0, 0, 0, 0, 0, 0, -t285, -t284, 0, -t185, 0, 0, 0, 0, 0, 0, 0, t285, t284, t174, 0, 0, 0, 0, 0, 0, t188, t189, t221, t141, 0, 0, 0, 0, 0, 0, t188, t221, -t189, t125; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t365, -t362, -t223, -qJ(1) * t223, 0, 0, -t363, 0, -t364, t299, t312 * t335 + t397, -t312 * t336 - t406, t313 * t185, -qJ(1) * t177 - (pkin(1) * t312 - pkin(5) * t313) * t185, t299, t363, t364, 0, 0, 0, -t313 * t174, -t312 * t337 - t397, t312 * t338 + t406, -qJ(1) * t156 - t143 * t312 + t151 * t313, t172, -t413, t169, t173, t411, t233, t120 * t313 - t138 * t312 - t386, t121 * t313 - t139 * t312 + t419, t313 * t134 + (-t133 - t385) * t312, -qJ(1) * t118 + t110 * t313 - t115 * t312, t172, t169, t413, t233, -t411, t173, t117 * t313 - t127 * t312 - t386, t313 * t122 + (-t124 - t385) * t312, t116 * t313 - t128 * t312 - t419, -qJ(1) * t113 + t108 * t313 - t109 * t312; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t362, -t365, t224, qJ(1) * t224, 0, 0, -t366, 0, -t367, -t346, -t313 * t335 + t398, t313 * t336 - t407, t312 * t185, qJ(1) * t178 - (-pkin(1) * t313 - pkin(5) * t312) * t185, -t346, t366, t367, 0, 0, 0, -t312 * t174, t313 * t337 - t398, -t313 * t338 + t407, qJ(1) * t157 + t143 * t313 + t151 * t312, t170, -t414, t168, t171, t412, t232, t120 * t312 + t138 * t313 + t159, t121 * t312 + t139 * t313 - t418, t133 * t313 + t134 * t312 + t208, qJ(1) * t119 + t110 * t312 + t115 * t313, t170, t168, t414, t232, -t412, t171, t117 * t312 + t127 * t313 + t159, t122 * t312 + t124 * t313 + t208, t116 * t312 + t128 * t313 + t418, qJ(1) * t114 + t108 * t312 + t109 * t313; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t288, t289, 0, 0, 0, 0, t284, 0, -t285, 0, t230, t229, t186, pkin(5) * t186 + t280, 0, -t284, t285, 0, 0, 0, t176, -t230, -t229, pkin(5) * t176 + t280 + (pkin(2) * t319 + qJ(3) * t317) * t288, t200, t396, t196, t201, -t395, t252, t148 * t319 + t149 * t317 + t357, t147 * t319 + t150 * t317 - t415, -pkin(3) * t377 + t144 * t319 + t219, -pkin(1) * t146 + pkin(5) * t142 + t123 * t319 + t129 * t317, t200, t196, -t396, t252, t395, t201, t131 * t319 + t140 * t317 + t357, t130 * t319 + t234 * t317 + t219, t132 * t319 + t137 * t317 + t415, -pkin(1) * t136 + pkin(5) * t126 + t111 * t319 + t112 * t317;];
tauB_reg = t1;
