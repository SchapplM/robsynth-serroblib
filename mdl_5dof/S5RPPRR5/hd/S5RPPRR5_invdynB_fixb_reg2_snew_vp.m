% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPPRR5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR5_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:40
% EndTime: 2019-12-31 17:56:45
% DurationCPUTime: 3.50s
% Computational Cost: add. (10629->317), mult. (16428->438), div. (0->0), fcn. (8290->8), ass. (0->203)
t322 = g(3) - qJDD(2);
t327 = sin(qJ(4));
t317 = -qJD(1) + qJD(4);
t315 = t317 ^ 2;
t316 = qJDD(1) - qJDD(4);
t330 = cos(qJ(4));
t356 = t330 * t316;
t335 = t327 * t315 + t356;
t256 = pkin(6) * t335 + t327 * t322;
t324 = sin(pkin(8));
t325 = cos(pkin(8));
t359 = t327 * t316;
t342 = -t330 * t315 + t359;
t378 = t324 * t342 + t325 * t335;
t380 = -pkin(6) * t342 + t330 * t322;
t179 = -qJ(2) * t378 + t256 * t325 - t324 * t380;
t328 = sin(qJ(1));
t331 = cos(qJ(1));
t227 = t324 * t335 - t325 * t342;
t395 = -qJ(2) * t227 + t256 * t324 + t325 * t380;
t396 = -t227 * t328 + t331 * t378;
t404 = -pkin(5) * t396 + t179 * t331 - t328 * t395;
t189 = t227 * t331 + t328 * t378;
t403 = -pkin(5) * t189 + t179 * t328 + t331 * t395;
t376 = qJD(1) ^ 2;
t290 = t324 * qJDD(1) + t325 * t376;
t291 = t325 * qJDD(1) - t324 * t376;
t238 = t331 * t290 + t291 * t328;
t259 = qJ(2) * t291 + t324 * t322;
t348 = -qJ(2) * t290 + t325 * t322;
t196 = -pkin(5) * t238 - t259 * t328 + t331 * t348;
t300 = t331 * g(1) + t328 * g(2);
t287 = -t376 * pkin(1) - t300;
t299 = t328 * g(1) - t331 * g(2);
t334 = qJDD(1) * pkin(1) + t299;
t236 = t325 * t287 + t324 * t334;
t313 = 2 * qJD(3) * qJD(1);
t318 = qJDD(1) * qJ(3);
t339 = t313 + t318 + t236;
t375 = pkin(2) + pkin(3);
t214 = -t375 * t376 + t339;
t323 = qJDD(1) * pkin(2);
t235 = t287 * t324 - t325 * t334;
t336 = -qJDD(3) - t235;
t224 = -t376 * qJ(3) - t323 - t336;
t333 = -qJDD(1) * pkin(3) + t224;
t175 = t327 * t214 - t330 * t333;
t176 = t330 * t214 + t327 * t333;
t146 = t175 * t330 - t176 * t327;
t147 = t175 * t327 + t176 * t330;
t125 = t146 * t324 - t147 * t325;
t388 = t146 * t325 + t147 * t324;
t117 = -t125 * t328 + t331 * t388;
t401 = t125 * t331 + t328 * t388;
t346 = t235 * t324 + t325 * t236;
t200 = t235 * t325 - t236 * t324;
t371 = t200 * t331;
t397 = -t328 * t346 + t371;
t372 = t200 * t328;
t153 = t331 * t346 + t372;
t343 = -t290 * t328 + t331 * t291;
t377 = pkin(5) * t343 + t259 * t331 + t328 * t348;
t223 = -t376 * pkin(2) + t339;
t186 = t223 * t324 - t224 * t325;
t347 = t325 * t223 + t224 * t324;
t149 = -t186 * t328 + t331 * t347;
t148 = t186 * t331 + t328 * t347;
t374 = pkin(1) * t322;
t326 = sin(qJ(5));
t329 = cos(qJ(5));
t298 = t329 * t315 * t326;
t288 = qJDD(5) + t298;
t368 = t288 * t326;
t367 = t288 * t329;
t289 = qJDD(5) - t298;
t366 = t289 * t326;
t365 = t289 * t329;
t320 = t326 ^ 2;
t362 = t320 * t315;
t169 = t316 * pkin(4) - t315 * pkin(7) + t175;
t361 = t326 * t169;
t360 = t326 * t316;
t357 = t329 * t169;
t302 = t329 * t316;
t170 = -pkin(4) * t315 - pkin(7) * t316 + t176;
t161 = t329 * t170 + t326 * t322;
t321 = t329 ^ 2;
t352 = t320 + t321;
t351 = qJD(5) * t317;
t350 = t326 * t351;
t349 = t329 * t351;
t160 = t170 * t326 - t329 * t322;
t254 = -t328 * t299 - t331 * t300;
t341 = t327 * t298;
t340 = t330 * t298;
t338 = -pkin(1) * t290 - t236;
t293 = t331 * qJDD(1) - t328 * t376;
t337 = -pkin(5) * t293 - g(3) * t328;
t136 = t160 * t329 - t161 * t326;
t137 = t326 * t160 + t329 * t161;
t253 = t299 * t331 - t300 * t328;
t332 = qJD(5) ^ 2;
t305 = t321 * t315;
t297 = -t305 - t332;
t296 = t305 - t332;
t295 = -t332 - t362;
t294 = t332 - t362;
t292 = t328 * qJDD(1) + t331 * t376;
t285 = pkin(1) * t291;
t280 = t305 - t362;
t279 = t305 + t362;
t274 = t352 * t316;
t270 = -t302 - 0.2e1 * t350;
t269 = -t302 - t350;
t268 = t349 - t360;
t267 = 0.2e1 * t349 - t360;
t266 = -pkin(5) * t292 + t331 * g(3);
t265 = t352 * t351;
t252 = qJDD(5) * t327 + t265 * t330;
t251 = qJDD(5) * t330 - t265 * t327;
t250 = t268 * t329 - t320 * t351;
t249 = -t269 * t326 - t321 * t351;
t248 = -t295 * t326 - t365;
t247 = -t294 * t326 + t367;
t246 = t297 * t329 - t368;
t245 = t296 * t329 - t366;
t244 = t295 * t329 - t366;
t243 = t297 * t326 + t367;
t233 = -t274 * t330 - t279 * t327;
t232 = -t274 * t327 + t279 * t330;
t225 = -t267 * t326 + t270 * t329;
t222 = t247 * t330 - t326 * t359;
t221 = t245 * t330 - t327 * t302;
t220 = -t247 * t327 - t326 * t356;
t219 = -t245 * t327 - t329 * t356;
t218 = t250 * t330 - t341;
t217 = t249 * t330 + t341;
t216 = -t250 * t327 - t340;
t215 = -t249 * t327 + t340;
t209 = t248 * t330 + t267 * t327;
t208 = t246 * t330 - t270 * t327;
t207 = t248 * t327 - t267 * t330;
t206 = t246 * t327 + t270 * t330;
t205 = -t251 * t324 + t252 * t325;
t204 = t251 * t325 + t252 * t324;
t203 = t225 * t330 - t280 * t327;
t202 = -t225 * t327 - t280 * t330;
t195 = qJ(2) * t346 + t374;
t194 = t232 * t324 + t233 * t325;
t193 = -t232 * t325 + t233 * t324;
t184 = -t220 * t324 + t222 * t325;
t183 = -t219 * t324 + t221 * t325;
t182 = t220 * t325 + t222 * t324;
t181 = t219 * t325 + t221 * t324;
t174 = -t216 * t324 + t218 * t325;
t173 = -t215 * t324 + t217 * t325;
t172 = t216 * t325 + t218 * t324;
t171 = t215 * t325 + t217 * t324;
t168 = -qJ(2) * t186 + (-pkin(2) * t324 + qJ(3) * t325) * t322;
t166 = t207 * t324 + t209 * t325;
t165 = t206 * t324 + t208 * t325;
t164 = -t207 * t325 + t209 * t324;
t163 = -t206 * t325 + t208 * t324;
t162 = qJ(2) * t347 + (pkin(2) * t325 + qJ(3) * t324 + pkin(1)) * t322;
t159 = -t202 * t324 + t203 * t325;
t158 = t202 * t325 + t203 * t324;
t157 = -pkin(7) * t244 + t357;
t156 = -pkin(7) * t243 + t361;
t155 = -pkin(4) * t244 + t161;
t154 = -pkin(4) * t243 + t160;
t151 = -t193 * t328 + t194 * t331;
t150 = t193 * t331 + t194 * t328;
t143 = pkin(6) * t146 + qJ(3) * t322;
t142 = -pkin(6) * t147 + t375 * t322;
t141 = -t164 * t328 + t166 * t331;
t140 = -t163 * t328 + t165 * t331;
t139 = t164 * t331 + t166 * t328;
t138 = t163 * t331 + t165 * t328;
t134 = -pkin(6) * t232 + t136 * t330;
t133 = -pkin(6) * t233 - t136 * t327;
t132 = t137 * t330 + t169 * t327;
t131 = t137 * t327 - t169 * t330;
t130 = -pkin(6) * t207 + qJ(3) * t244 - t155 * t327 + t157 * t330;
t129 = -pkin(6) * t206 + qJ(3) * t243 - t154 * t327 + t156 * t330;
t124 = -pkin(6) * t209 - t330 * t155 - t327 * t157 + t375 * t244;
t123 = -pkin(6) * t208 - t330 * t154 - t327 * t156 + t375 * t243;
t122 = -qJ(2) * t193 - t133 * t324 + t134 * t325;
t121 = qJ(2) * t194 + t133 * t325 + t134 * t324;
t120 = t131 * t324 + t132 * t325;
t119 = -t131 * t325 + t132 * t324;
t116 = -qJ(2) * t388 - t142 * t324 + t143 * t325;
t115 = -qJ(2) * t125 + t142 * t325 + t143 * t324 + t374;
t114 = -qJ(2) * t164 - t124 * t324 + t130 * t325;
t113 = -qJ(2) * t163 - t123 * t324 + t129 * t325;
t112 = pkin(1) * t244 + qJ(2) * t166 + t124 * t325 + t130 * t324;
t111 = pkin(1) * t243 + qJ(2) * t165 + t123 * t325 + t129 * t324;
t110 = -pkin(6) * t131 - (pkin(4) * t327 - pkin(7) * t330 + qJ(3)) * t136;
t109 = -pkin(6) * t132 - (pkin(4) * t330 + pkin(7) * t327 + t375) * t136;
t108 = -t119 * t328 + t120 * t331;
t107 = t119 * t331 + t120 * t328;
t106 = -qJ(2) * t119 - t109 * t324 + t110 * t325;
t105 = -pkin(1) * t136 + qJ(2) * t120 + t109 * t325 + t110 * t324;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t292, -t293, 0, t254, 0, 0, 0, 0, 0, 0, -t238, -t343, 0, t153, 0, 0, 0, 0, 0, 0, -t238, 0, t343, t149, 0, 0, 0, 0, 0, 0, -t189, t396, 0, -t401, 0, 0, 0, 0, 0, 0, t140, t141, t151, t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t293, -t292, 0, t253, 0, 0, 0, 0, 0, 0, t343, -t238, 0, -t397, 0, 0, 0, 0, 0, 0, t343, 0, t238, t148, 0, 0, 0, 0, 0, 0, t396, t189, 0, t117, 0, 0, 0, 0, 0, 0, t138, t139, t150, t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t322, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t322, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t322, 0, 0, 0, 0, 0, 0, -t243, -t244, 0, t136; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t293, 0, -t292, 0, t337, -t266, -t253, -pkin(5) * t253, 0, 0, t343, 0, -t238, 0, -t377, -t196, t397, pkin(5) * t397 + qJ(2) * t371 - t195 * t328, 0, t343, 0, 0, t238, 0, -t377, -t148, t196, -pkin(5) * t148 - t162 * t328 + t168 * t331, 0, 0, -t396, 0, -t189, 0, t404, t403, t117, -pkin(5) * t117 - t115 * t328 + t116 * t331, -t172 * t328 + t174 * t331, -t158 * t328 + t159 * t331, -t182 * t328 + t184 * t331, -t171 * t328 + t173 * t331, -t181 * t328 + t183 * t331, -t204 * t328 + t205 * t331, -pkin(5) * t138 - t111 * t328 + t113 * t331, -pkin(5) * t139 - t112 * t328 + t114 * t331, -pkin(5) * t150 - t121 * t328 + t122 * t331, -pkin(5) * t107 - t105 * t328 + t106 * t331; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t292, 0, t293, 0, t266, t337, t254, pkin(5) * t254, 0, 0, t238, 0, t343, 0, t196, -t377, t153, pkin(5) * t153 + qJ(2) * t372 + t195 * t331, 0, t238, 0, 0, -t343, 0, t196, t149, t377, pkin(5) * t149 + t162 * t331 + t168 * t328, 0, 0, -t189, 0, t396, 0, t403, -t404, t401, -pkin(5) * t401 + t115 * t331 + t116 * t328, t172 * t331 + t174 * t328, t158 * t331 + t159 * t328, t182 * t331 + t184 * t328, t171 * t331 + t173 * t328, t181 * t331 + t183 * t328, t204 * t331 + t205 * t328, pkin(5) * t140 + t111 * t331 + t113 * t328, pkin(5) * t141 + t112 * t331 + t114 * t328, pkin(5) * t151 + t121 * t331 + t122 * t328, pkin(5) * t108 + t105 * t331 + t106 * t328; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t299, t300, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t235 + t285, t338, 0, -pkin(1) * t200, 0, 0, 0, qJDD(1), 0, 0, t285 + 0.2e1 * t323 + t336, 0, t313 + 0.2e1 * t318 - t338, pkin(1) * t186 - pkin(2) * t224 + qJ(3) * t223, 0, 0, 0, 0, 0, t316, pkin(1) * t378 + qJ(3) * t342 + t335 * t375 + t175, pkin(1) * t227 + qJ(3) * t335 - t342 * t375 + t176, 0, pkin(1) * t388 + qJ(3) * t147 + t146 * t375, (-t268 - t349) * t326, -t267 * t329 - t270 * t326, -t294 * t329 - t368, (-t269 + t350) * t329, -t296 * t326 - t365, 0, pkin(1) * t163 - pkin(4) * t270 - pkin(7) * t246 + qJ(3) * t208 - t375 * t206 + t357, pkin(1) * t164 + pkin(4) * t267 - pkin(7) * t248 + qJ(3) * t209 - t375 * t207 - t361, pkin(1) * t193 - pkin(4) * t279 + pkin(7) * t274 + qJ(3) * t233 - t375 * t232 - t137, pkin(1) * t119 + pkin(4) * t169 - pkin(7) * t137 + qJ(3) * t132 - t375 * t131;];
tauB_reg = t1;
