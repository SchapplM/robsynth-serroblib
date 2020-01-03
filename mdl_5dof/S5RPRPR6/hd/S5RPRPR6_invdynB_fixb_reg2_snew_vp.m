% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRPR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRPR6_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR6_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:51
% EndTime: 2019-12-31 18:17:57
% DurationCPUTime: 4.38s
% Computational Cost: add. (11729->320), mult. (17357->448), div. (0->0), fcn. (9986->8), ass. (0->212)
t336 = (qJD(1) + qJD(3));
t334 = t336 ^ 2;
t343 = sin(qJ(3));
t335 = qJDD(1) + qJDD(3);
t346 = cos(qJ(3));
t371 = t346 * t335;
t303 = t343 * t334 - t371;
t339 = g(3) - qJDD(2);
t279 = pkin(6) * t303 - t343 * t339;
t340 = sin(pkin(8));
t341 = cos(pkin(8));
t381 = t343 * t335;
t300 = t346 * t334 + t381;
t357 = t340 * t300 + t341 * t303;
t362 = pkin(6) * t300 - t346 * t339;
t201 = qJ(2) * t357 + t341 * t279 + t340 * t362;
t247 = t341 * t300 - t340 * t303;
t344 = sin(qJ(1));
t347 = cos(qJ(1));
t211 = t347 * t247 - t344 * t357;
t418 = qJ(2) * t247 - t340 * t279 + t341 * t362;
t431 = pkin(5) * t211 - t344 * t201 + t347 * t418;
t415 = t344 * t247 + t347 * t357;
t422 = pkin(5) * t415 + t347 * t201 + t344 * t418;
t321 = t344 * g(1) - t347 * g(2);
t309 = qJDD(1) * pkin(1) + t321;
t322 = t347 * g(1) + t344 * g(2);
t349 = qJD(1) ^ 2;
t310 = -t349 * pkin(1) - t322;
t260 = -t341 * t309 + t340 * t310;
t258 = qJDD(1) * pkin(2) - t260;
t261 = t340 * t309 + t341 * t310;
t259 = -t349 * pkin(2) + t261;
t223 = -t346 * t258 + t343 * t259;
t224 = t343 * t258 + t346 * t259;
t360 = t343 * t223 + t346 * t224;
t177 = t346 * t223 - t343 * t224;
t389 = t341 * t177;
t149 = -t340 * t360 + t389;
t392 = t340 * t177;
t427 = t341 * t360 + t392;
t134 = t344 * t149 + t347 * t427;
t437 = t347 * t149 - t344 * t427;
t434 = -pkin(1) * t247 - pkin(2) * t300;
t433 = -pkin(1) * t357 - pkin(2) * t303;
t350 = (2 * qJD(4) * t336) + t224;
t395 = t335 * qJ(4);
t197 = -t334 * pkin(3) + t350 + t395;
t327 = t335 * pkin(3);
t351 = qJDD(4) + t223 - t327;
t210 = -t334 * qJ(4) + t351;
t165 = t343 * t197 - t346 * t210;
t361 = t346 * t197 + t343 * t210;
t144 = t341 * t165 + t340 * t361;
t417 = -t340 * t165 + t341 * t361;
t128 = -t344 * t144 + t347 * t417;
t127 = t347 * t144 + t344 * t417;
t359 = t340 * t260 + t341 * t261;
t227 = t341 * t260 - t340 * t261;
t369 = t347 * t227;
t428 = -t344 * t359 + t369;
t379 = t344 * t227;
t180 = t347 * t359 + t379;
t313 = t340 * qJDD(1) + t341 * t349;
t314 = t341 * qJDD(1) - t340 * t349;
t262 = -t344 * t313 + t347 * t314;
t287 = qJ(2) * t313 - t341 * t339;
t352 = -qJ(2) * t314 - t340 * t339;
t426 = -pkin(5) * t262 + t344 * t287 + t347 * t352;
t401 = t347 * t313 + t344 * t314;
t413 = pkin(5) * t401 + t347 * t287 - t344 * t352;
t397 = pkin(3) + pkin(7);
t396 = pkin(1) * t339;
t342 = sin(qJ(5));
t337 = t342 ^ 2;
t394 = t337 * t334;
t345 = cos(qJ(5));
t338 = t345 ^ 2;
t393 = t338 * t334;
t189 = -t334 * pkin(7) + t197;
t386 = t342 * t189;
t366 = t342 * t334 * t345;
t311 = qJDD(5) + t366;
t385 = t342 * t311;
t312 = qJDD(5) - t366;
t384 = t342 * t312;
t383 = t342 * t335;
t368 = t337 + t338;
t299 = t368 * t335;
t382 = t343 * t299;
t376 = t345 * t189;
t375 = t345 * t311;
t374 = t345 * t312;
t373 = t345 * t335;
t372 = t346 * t299;
t367 = qJD(5) * t336;
t365 = t342 * t367;
t364 = t345 * t367;
t196 = -t335 * pkin(7) + t210;
t187 = -t345 * t196 - t342 * t339;
t277 = -t344 * t321 - t347 * t322;
t355 = t343 * t366;
t354 = t346 * t366;
t316 = t347 * qJDD(1) - t344 * t349;
t353 = -pkin(5) * t316 - t344 * g(3);
t188 = t342 * t196 - t345 * t339;
t160 = -t345 * t187 + t342 * t188;
t161 = t342 * t187 + t345 * t188;
t276 = t347 * t321 - t344 * t322;
t348 = qJD(5) ^ 2;
t320 = -t348 - t393;
t319 = t348 - t393;
t318 = -t348 - t394;
t317 = -t348 + t394;
t315 = t344 * qJDD(1) + t347 * t349;
t305 = (-t337 + t338) * t334;
t304 = t368 * t334;
t295 = -0.2e1 * t365 + t373;
t294 = -t365 + t373;
t293 = -t364 - t383;
t292 = 0.2e1 * t364 + t383;
t291 = -pkin(5) * t315 + t347 * g(3);
t290 = t368 * t367;
t275 = t346 * qJDD(5) - t343 * t290;
t274 = t343 * qJDD(5) + t346 * t290;
t273 = -t342 * t294 - t338 * t367;
t272 = -t345 * t293 - t337 * t367;
t271 = -t342 * t320 - t375;
t270 = t345 * t318 - t384;
t269 = t345 * t320 - t385;
t268 = -t345 * t319 - t384;
t267 = t342 * t318 + t374;
t266 = -t342 * t317 - t375;
t256 = -t346 * t304 - t382;
t255 = -t343 * t304 + t372;
t245 = t342 * t292 - t345 * t295;
t244 = -t343 * t268 + t345 * t371;
t243 = -t343 * t266 - t342 * t371;
t242 = t346 * t268 + t343 * t373;
t241 = t346 * t266 - t342 * t381;
t240 = -t343 * t272 - t354;
t239 = -t343 * t273 + t354;
t238 = t346 * t272 - t355;
t237 = t346 * t273 + t355;
t236 = t343 * t269 + t346 * t295;
t235 = t343 * t267 + t346 * t292;
t234 = -t346 * t269 + t343 * t295;
t233 = -t346 * t267 + t343 * t292;
t232 = -t340 * t274 + t341 * t275;
t231 = t341 * t274 + t340 * t275;
t230 = -t343 * t245 + t346 * t305;
t229 = t346 * t245 + t343 * t305;
t222 = qJ(2) * t359 + t396;
t220 = -t340 * t255 + t341 * t256;
t219 = t341 * t255 + t340 * t256;
t209 = -t340 * t242 + t341 * t244;
t208 = -t340 * t241 + t341 * t243;
t207 = t341 * t242 + t340 * t244;
t206 = t341 * t241 + t340 * t243;
t195 = -t340 * t238 + t341 * t240;
t194 = -t340 * t237 + t341 * t239;
t193 = t341 * t238 + t340 * t240;
t192 = t341 * t237 + t340 * t239;
t186 = -t340 * t234 + t341 * t236;
t185 = -t340 * t233 + t341 * t235;
t184 = t341 * t234 + t340 * t236;
t183 = t341 * t233 + t340 * t235;
t182 = -t340 * t229 + t341 * t230;
t181 = t341 * t229 + t340 * t230;
t174 = pkin(4) * t269 - qJ(4) * t271 - t188;
t173 = pkin(4) * t267 - qJ(4) * t270 - t187;
t172 = -t344 * t219 + t347 * t220;
t171 = t347 * t219 + t344 * t220;
t170 = pkin(2) * t339 + pkin(6) * t360;
t169 = pkin(4) * t292 - t397 * t270 + t376;
t168 = pkin(4) * t295 - t397 * t271 - t386;
t163 = -pkin(6) * t165 + (-pkin(3) * t343 + qJ(4) * t346) * t339;
t162 = pkin(6) * t361 + (pkin(3) * t346 + qJ(4) * t343 + pkin(2)) * t339;
t159 = -t344 * t184 + t347 * t186;
t158 = -t344 * t183 + t347 * t185;
t157 = t347 * t184 + t344 * t186;
t156 = t347 * t183 + t344 * t185;
t155 = -pkin(4) * t304 - t161;
t154 = -pkin(4) * t372 - pkin(6) * t255 - t343 * t155;
t153 = -pkin(4) * t382 + pkin(6) * t256 + t346 * t155;
t152 = t343 * t160 + t346 * t189;
t151 = -t346 * t160 + t343 * t189;
t142 = -pkin(6) * t234 - t343 * t168 + t346 * t174;
t141 = -pkin(6) * t233 - t343 * t169 + t346 * t173;
t140 = -pkin(2) * t271 + pkin(6) * t236 + t346 * t168 + t343 * t174;
t139 = -pkin(2) * t270 + pkin(6) * t235 + t346 * t169 + t343 * t173;
t138 = pkin(4) * t160 - qJ(4) * t161;
t137 = pkin(4) * t189 - t397 * t161;
t136 = -t340 * t151 + t341 * t152;
t135 = t341 * t151 + t340 * t152;
t132 = pkin(6) * t389 + qJ(2) * t149 - t340 * t170;
t131 = pkin(6) * t392 + qJ(2) * t427 + t341 * t170 + t396;
t130 = -qJ(2) * t219 - t340 * t153 + t341 * t154;
t129 = qJ(2) * t220 + t341 * t153 + t340 * t154;
t126 = -qJ(2) * t144 - t340 * t162 + t341 * t163;
t125 = qJ(2) * t417 + t341 * t162 + t340 * t163 + t396;
t124 = -qJ(2) * t184 - t340 * t140 + t341 * t142;
t123 = -qJ(2) * t183 - t340 * t139 + t341 * t141;
t122 = -pkin(1) * t271 + qJ(2) * t186 + t341 * t140 + t340 * t142;
t121 = -pkin(1) * t270 + qJ(2) * t185 + t341 * t139 + t340 * t141;
t120 = -pkin(6) * t151 - t343 * t137 + t346 * t138;
t119 = -t344 * t135 + t347 * t136;
t118 = t347 * t135 + t344 * t136;
t117 = -pkin(2) * t161 + pkin(6) * t152 + t346 * t137 + t343 * t138;
t116 = -qJ(2) * t135 - t340 * t117 + t341 * t120;
t115 = -pkin(1) * t161 + qJ(2) * t136 + t341 * t117 + t340 * t120;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t315, -t316, 0, t277, 0, 0, 0, 0, 0, 0, -t401, -t262, 0, t180, 0, 0, 0, 0, 0, 0, -t211, t415, 0, t134, 0, 0, 0, 0, 0, 0, 0, t211, -t415, t128, 0, 0, 0, 0, 0, 0, t158, t159, t172, t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t316, -t315, 0, t276, 0, 0, 0, 0, 0, 0, t262, -t401, 0, -t428, 0, 0, 0, 0, 0, 0, -t415, -t211, 0, -t437, 0, 0, 0, 0, 0, 0, 0, t415, t211, t127, 0, 0, 0, 0, 0, 0, t156, t157, t171, t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t339, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t339, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t339, 0, 0, 0, 0, 0, 0, t270, t271, 0, t161; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t316, 0, -t315, 0, t353, -t291, -t276, -pkin(5) * t276, 0, 0, t262, 0, -t401, 0, t426, t413, t428, pkin(5) * t428 + qJ(2) * t369 - t344 * t222, 0, 0, -t415, 0, -t211, 0, t422, t431, t437, pkin(5) * t437 - t344 * t131 + t347 * t132, 0, t415, t211, 0, 0, 0, -t127, -t422, -t431, -pkin(5) * t127 - t344 * t125 + t347 * t126, -t344 * t192 + t347 * t194, -t344 * t181 + t347 * t182, -t344 * t207 + t347 * t209, -t344 * t193 + t347 * t195, -t344 * t206 + t347 * t208, -t344 * t231 + t347 * t232, -pkin(5) * t156 - t344 * t121 + t347 * t123, -pkin(5) * t157 - t344 * t122 + t347 * t124, -pkin(5) * t171 - t344 * t129 + t347 * t130, -pkin(5) * t118 - t344 * t115 + t347 * t116; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t315, 0, t316, 0, t291, t353, t277, pkin(5) * t277, 0, 0, t401, 0, t262, 0, -t413, t426, t180, pkin(5) * t180 + qJ(2) * t379 + t347 * t222, 0, 0, t211, 0, -t415, 0, -t431, t422, t134, pkin(5) * t134 + t347 * t131 + t344 * t132, 0, -t211, t415, 0, 0, 0, t128, t431, -t422, pkin(5) * t128 + t347 * t125 + t344 * t126, t347 * t192 + t344 * t194, t347 * t181 + t344 * t182, t347 * t207 + t344 * t209, t347 * t193 + t344 * t195, t347 * t206 + t344 * t208, t347 * t231 + t344 * t232, pkin(5) * t158 + t347 * t121 + t344 * t123, pkin(5) * t159 + t347 * t122 + t344 * t124, pkin(5) * t172 + t347 * t129 + t344 * t130, pkin(5) * t119 + t347 * t115 + t344 * t116; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t321, t322, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t314 - t260, -pkin(1) * t313 - t261, 0, -pkin(1) * t227, 0, 0, 0, 0, 0, t335, -t223 + t433, -t224 + t434, 0, -pkin(1) * t149 - pkin(2) * t177, t335, 0, 0, 0, 0, 0, 0, -t327 + t351 - t433, t350 + 0.2e1 * t395 - t434, pkin(1) * t144 + pkin(2) * t165 - pkin(3) * t210 + qJ(4) * t197, (t294 - t365) * t345, -t345 * t292 - t342 * t295, -t342 * t319 + t374, (-t293 + t364) * t342, t345 * t317 - t385, 0, pkin(1) * t183 + pkin(2) * t233 + qJ(4) * t292 - t397 * t267 + t386, pkin(1) * t184 + pkin(2) * t234 + qJ(4) * t295 - t397 * t269 + t376, pkin(1) * t219 + pkin(2) * t255 - qJ(4) * t304 + t397 * t299 - t160, pkin(1) * t135 + pkin(2) * t151 + qJ(4) * t189 - t397 * t160;];
tauB_reg = t1;
