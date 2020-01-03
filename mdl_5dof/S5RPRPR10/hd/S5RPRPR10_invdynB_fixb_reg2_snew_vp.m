% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRPR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRPR10_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR10_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:07
% EndTime: 2019-12-31 18:26:13
% DurationCPUTime: 3.98s
% Computational Cost: add. (13254->311), mult. (18911->434), div. (0->0), fcn. (8076->8), ass. (0->195)
t311 = -qJD(1) + qJD(3);
t309 = t311 ^ 2;
t310 = qJDD(1) - qJDD(3);
t316 = sin(pkin(8));
t317 = cos(pkin(8));
t274 = t317 * t309 - t316 * t310;
t276 = t316 * t309 + t317 * t310;
t319 = sin(qJ(3));
t322 = cos(qJ(3));
t232 = t319 * t274 + t322 * t276;
t315 = g(3) + qJDD(4);
t255 = qJ(4) * t274 + t317 * t315;
t367 = qJ(4) * t276 + t316 * t315;
t191 = pkin(6) * t232 + t319 * t255 + t322 * t367;
t320 = sin(qJ(1));
t323 = cos(qJ(1));
t364 = -t322 * t274 + t319 * t276;
t375 = -pkin(6) * t364 + t322 * t255 - t319 * t367;
t376 = t323 * t232 + t320 * t364;
t381 = -pkin(5) * t376 + t191 * t323 - t320 * t375;
t197 = t320 * t232 - t323 * t364;
t380 = -pkin(5) * t197 + t191 * t320 + t323 * t375;
t325 = qJD(1) ^ 2;
t294 = t320 * g(1) - t323 * g(2);
t336 = -qJDD(2) + t294;
t330 = -t325 * qJ(2) - t336;
t359 = pkin(1) + pkin(2);
t328 = -t359 * qJDD(1) + t330;
t312 = qJDD(1) * qJ(2);
t295 = t323 * g(1) + t320 * g(2);
t332 = (2 * qJD(2) * qJD(1)) - t295;
t331 = t312 + t332;
t329 = -t359 * t325 + t331;
t216 = t319 * t328 + t322 * t329;
t210 = -t309 * pkin(3) + t216;
t327 = -t319 * t329 + t322 * t328;
t326 = -t310 * pkin(3) + t327;
t173 = t316 * t210 - t317 * t326;
t174 = t317 * t210 + t316 * t326;
t340 = t173 * t316 + t317 * t174;
t151 = t173 * t317 - t174 * t316;
t355 = t151 * t322;
t135 = -t319 * t340 + t355;
t356 = t151 * t319;
t136 = t322 * t340 + t356;
t125 = t135 * t323 + t136 * t320;
t377 = t135 * t320 - t136 * t323;
t333 = t319 * t309 + t322 * t310;
t259 = pkin(6) * t333 + t319 * g(3);
t339 = -t322 * t309 + t319 * t310;
t360 = t320 * t339 + t323 * t333;
t361 = -pkin(6) * t339 + t322 * g(3);
t371 = -pkin(5) * t360 + t259 * t323 - t320 * t361;
t236 = t320 * t333 - t323 * t339;
t369 = -pkin(5) * t236 + t259 * t320 + t323 * t361;
t183 = -t319 * t216 - t322 * t327;
t184 = t322 * t216 - t319 * t327;
t157 = t183 * t323 + t184 * t320;
t365 = t183 * t320 - t184 * t323;
t357 = qJDD(1) * pkin(1);
t318 = sin(qJ(5));
t313 = t318 ^ 2;
t354 = t313 * t309;
t170 = t310 * pkin(4) - t309 * pkin(7) + t173;
t351 = t318 * t170;
t321 = cos(qJ(5));
t293 = t321 * t309 * t318;
t285 = qJDD(5) + t293;
t350 = t318 * t285;
t286 = qJDD(5) - t293;
t349 = t318 * t286;
t348 = t318 * t310;
t347 = t321 * t170;
t346 = t321 * t285;
t345 = t321 * t286;
t297 = t321 * t310;
t171 = -t309 * pkin(4) - t310 * pkin(7) + t174;
t166 = t321 * t171 + t318 * t315;
t314 = t321 ^ 2;
t344 = t313 + t314;
t343 = qJD(5) * t311;
t342 = t318 * t343;
t341 = t321 * t343;
t165 = t171 * t318 - t321 * t315;
t263 = -t325 * pkin(1) + t331;
t265 = -t330 + t357;
t222 = t323 * t263 - t320 * t265;
t251 = -t320 * t294 - t323 * t295;
t338 = t316 * t293;
t337 = t317 * t293;
t287 = t320 * qJDD(1) + t323 * t325;
t268 = -pkin(5) * t287 + t323 * g(3);
t288 = t323 * qJDD(1) - t320 * t325;
t267 = pkin(5) * t288 + t320 * g(3);
t146 = t165 * t321 - t166 * t318;
t147 = t318 * t165 + t321 * t166;
t221 = t320 * t263 + t323 * t265;
t250 = t323 * t294 - t320 * t295;
t324 = qJD(5) ^ 2;
t299 = t314 * t309;
t292 = -t299 - t324;
t291 = t299 - t324;
t290 = -t324 - t354;
t289 = t324 - t354;
t283 = t299 - t354;
t282 = t299 + t354;
t277 = t344 * t310;
t272 = -t297 - 0.2e1 * t342;
t271 = -t297 - t342;
t270 = t341 - t348;
t269 = 0.2e1 * t341 - t348;
t266 = t344 * t343;
t249 = t316 * qJDD(5) + t317 * t266;
t248 = -t317 * qJDD(5) + t316 * t266;
t247 = t321 * t270 - t313 * t343;
t246 = -t318 * t271 - t314 * t343;
t245 = -t318 * t290 - t345;
t244 = -t318 * t289 + t346;
t243 = t321 * t292 - t350;
t242 = t321 * t291 - t349;
t241 = t321 * t290 - t349;
t240 = t318 * t292 + t346;
t235 = -t317 * t277 - t316 * t282;
t234 = -t316 * t277 + t317 * t282;
t227 = -t318 * t269 + t321 * t272;
t226 = t317 * t244 - t316 * t348;
t225 = t317 * t242 - t316 * t297;
t224 = t316 * t244 + t317 * t348;
t223 = t316 * t242 + t317 * t297;
t220 = t317 * t247 - t338;
t219 = t317 * t246 + t338;
t218 = t316 * t247 + t337;
t217 = t316 * t246 - t337;
t214 = t317 * t245 + t316 * t269;
t213 = t317 * t243 - t316 * t272;
t212 = t316 * t245 - t317 * t269;
t211 = t316 * t243 + t317 * t272;
t206 = -t319 * t248 + t322 * t249;
t205 = -t322 * t248 - t319 * t249;
t204 = t317 * t227 - t316 * t283;
t203 = t316 * t227 + t317 * t283;
t202 = -t319 * t234 + t322 * t235;
t201 = t322 * t234 + t319 * t235;
t196 = -t319 * t224 + t322 * t226;
t195 = -t319 * t223 + t322 * t225;
t194 = -t322 * t224 - t319 * t226;
t193 = -t322 * t223 - t319 * t225;
t188 = -t319 * t218 + t322 * t220;
t187 = -t319 * t217 + t322 * t219;
t186 = -t322 * t218 - t319 * t220;
t185 = -t322 * t217 - t319 * t219;
t180 = pkin(6) * t183 + qJ(2) * g(3);
t179 = -t319 * t212 + t322 * t214;
t178 = -t319 * t211 + t322 * t213;
t177 = t322 * t212 + t319 * t214;
t176 = t322 * t211 + t319 * t213;
t175 = -pkin(6) * t184 + t359 * g(3);
t169 = -t319 * t203 + t322 * t204;
t168 = -t322 * t203 - t319 * t204;
t164 = -pkin(7) * t241 + t347;
t163 = -pkin(7) * t240 + t351;
t162 = -pkin(4) * t241 + t166;
t161 = -pkin(4) * t240 + t165;
t160 = t201 * t320 + t202 * t323;
t159 = -t201 * t323 + t202 * t320;
t156 = t177 * t320 + t179 * t323;
t155 = t176 * t320 + t178 * t323;
t154 = -t177 * t323 + t179 * t320;
t153 = -t176 * t323 + t178 * t320;
t148 = -pkin(3) * t315 + qJ(4) * t340;
t144 = -qJ(4) * t234 + t146 * t317;
t143 = qJ(4) * t235 + t146 * t316;
t142 = -qJ(4) * t212 - t162 * t316 + t164 * t317;
t141 = -qJ(4) * t211 - t161 * t316 + t163 * t317;
t140 = t147 * t317 + t170 * t316;
t139 = t147 * t316 - t170 * t317;
t138 = -pkin(3) * t241 + qJ(4) * t214 + t162 * t317 + t164 * t316;
t137 = -pkin(3) * t240 + qJ(4) * t213 + t161 * t317 + t163 * t316;
t132 = -pkin(6) * t201 - t143 * t319 + t144 * t322;
t131 = -pkin(6) * t202 - t143 * t322 - t144 * t319;
t130 = -t139 * t319 + t140 * t322;
t129 = t139 * t322 + t140 * t319;
t128 = -pkin(6) * t177 + qJ(2) * t241 - t138 * t319 + t142 * t322;
t127 = -pkin(6) * t176 + qJ(2) * t240 - t137 * t319 + t141 * t322;
t124 = -pkin(6) * t179 - t322 * t138 - t319 * t142 + t359 * t241;
t123 = -pkin(6) * t178 - t322 * t137 - t319 * t141 + t359 * t240;
t122 = -qJ(4) * t139 - (pkin(4) * t316 - pkin(7) * t317) * t146;
t121 = pkin(6) * t135 + qJ(2) * t315 + qJ(4) * t355 - t148 * t319;
t120 = -pkin(6) * t136 - qJ(4) * t356 - t322 * t148 + t359 * t315;
t119 = qJ(4) * t140 - (-pkin(4) * t317 - pkin(7) * t316 - pkin(3)) * t146;
t118 = t129 * t320 + t130 * t323;
t117 = -t129 * t323 + t130 * t320;
t116 = -pkin(6) * t129 - qJ(2) * t146 - t119 * t319 + t122 * t322;
t115 = -pkin(6) * t130 - t322 * t119 - t319 * t122 - t146 * t359;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t287, -t288, 0, t251, 0, 0, 0, 0, 0, 0, -t287, 0, t288, t222, 0, 0, 0, 0, 0, 0, -t236, t360, 0, -t365, 0, 0, 0, 0, 0, 0, -t197, t376, 0, -t377, 0, 0, 0, 0, 0, 0, t155, t156, t160, t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t288, -t287, 0, t250, 0, 0, 0, 0, 0, 0, t288, 0, t287, t221, 0, 0, 0, 0, 0, 0, t360, t236, 0, t157, 0, 0, 0, 0, 0, 0, t376, t197, 0, t125, 0, 0, 0, 0, 0, 0, t153, t154, t159, t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t315, 0, 0, 0, 0, 0, 0, -t240, -t241, 0, t146; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t288, 0, -t287, 0, -t267, -t268, -t250, -pkin(5) * t250, 0, t288, 0, 0, t287, 0, -t267, -t221, t268, -pkin(5) * t221 + (-pkin(1) * t320 + qJ(2) * t323) * g(3), 0, 0, -t360, 0, -t236, 0, t371, t369, t157, -pkin(5) * t157 - t175 * t320 + t180 * t323, 0, 0, -t376, 0, -t197, 0, t381, t380, t125, -pkin(5) * t125 - t120 * t320 + t121 * t323, -t186 * t320 + t188 * t323, -t168 * t320 + t169 * t323, -t194 * t320 + t196 * t323, -t185 * t320 + t187 * t323, -t193 * t320 + t195 * t323, -t205 * t320 + t206 * t323, -pkin(5) * t153 - t123 * t320 + t127 * t323, -pkin(5) * t154 - t124 * t320 + t128 * t323, -pkin(5) * t159 - t131 * t320 + t132 * t323, -pkin(5) * t117 - t115 * t320 + t116 * t323; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t287, 0, t288, 0, t268, -t267, t251, pkin(5) * t251, 0, t287, 0, 0, -t288, 0, t268, t222, t267, pkin(5) * t222 + (pkin(1) * t323 + qJ(2) * t320) * g(3), 0, 0, -t236, 0, t360, 0, t369, -t371, t365, -pkin(5) * t365 + t175 * t323 + t180 * t320, 0, 0, -t197, 0, t376, 0, t380, -t381, t377, -pkin(5) * t377 + t120 * t323 + t121 * t320, t186 * t323 + t188 * t320, t168 * t323 + t169 * t320, t194 * t323 + t196 * t320, t185 * t323 + t187 * t320, t193 * t323 + t195 * t320, t205 * t323 + t206 * t320, pkin(5) * t155 + t123 * t323 + t127 * t320, pkin(5) * t156 + t124 * t323 + t128 * t320, pkin(5) * t160 + t131 * t323 + t132 * t320, pkin(5) * t118 + t115 * t323 + t116 * t320; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t294, t295, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t336 + 0.2e1 * t357, 0, 0.2e1 * t312 + t332, pkin(1) * t265 + qJ(2) * t263, 0, 0, 0, 0, 0, t310, qJ(2) * t339 + t333 * t359 - t327, qJ(2) * t333 - t339 * t359 + t216, 0, qJ(2) * t184 + t183 * t359, 0, 0, 0, 0, 0, t310, pkin(3) * t276 + qJ(2) * t364 + t232 * t359 + t173, pkin(3) * t274 + qJ(2) * t232 - t359 * t364 + t174, 0, pkin(3) * t151 + qJ(2) * t136 + t135 * t359, (-t270 - t341) * t318, -t269 * t321 - t272 * t318, -t289 * t321 - t350, (-t271 + t342) * t321, -t291 * t318 - t345, 0, -pkin(3) * t211 - pkin(4) * t272 - pkin(7) * t243 + qJ(2) * t178 - t359 * t176 + t347, -pkin(3) * t212 + pkin(4) * t269 - pkin(7) * t245 + qJ(2) * t179 - t359 * t177 - t351, -pkin(3) * t234 - pkin(4) * t282 + pkin(7) * t277 + qJ(2) * t202 - t359 * t201 - t147, -pkin(3) * t139 + pkin(4) * t170 - pkin(7) * t147 + qJ(2) * t130 - t359 * t129;];
tauB_reg = t1;
