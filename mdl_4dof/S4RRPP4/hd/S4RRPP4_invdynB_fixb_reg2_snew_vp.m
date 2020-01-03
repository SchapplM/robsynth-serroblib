% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRPP4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP4_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_invdynB_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:20
% EndTime: 2019-12-31 16:59:23
% DurationCPUTime: 2.06s
% Computational Cost: add. (2403->271), mult. (5533->315), div. (0->0), fcn. (2597->4), ass. (0->188)
t300 = qJD(2) ^ 2;
t296 = sin(qJ(2));
t291 = t296 ^ 2;
t301 = qJD(1) ^ 2;
t363 = t291 * t301;
t271 = t300 + t363;
t298 = cos(qJ(2));
t330 = t298 * t301 * t296;
t267 = qJDD(2) - t330;
t348 = t298 * t267;
t228 = -t296 * t271 + t348;
t337 = qJD(1) * qJD(2);
t325 = t298 * t337;
t335 = t296 * qJDD(1);
t255 = 0.2e1 * t325 + t335;
t297 = sin(qJ(1));
t299 = cos(qJ(1));
t189 = t297 * t228 + t299 * t255;
t373 = pkin(4) * t189;
t192 = t299 * t228 - t297 * t255;
t187 = pkin(4) * t192;
t344 = pkin(1) * t255 + pkin(5) * t228;
t266 = qJDD(2) + t330;
t248 = t298 * t266;
t272 = -t300 + t363;
t227 = t296 * t272 + t248;
t332 = t299 * qJDD(1);
t196 = t297 * t227 - t296 * t332;
t334 = t297 * qJDD(1);
t198 = t299 * t227 + t296 * t334;
t357 = t296 * t267;
t222 = t298 * t271 + t357;
t383 = pkin(1) * t222;
t292 = t298 ^ 2;
t341 = t291 + t292;
t260 = t341 * qJDD(1);
t263 = t341 * t301;
t207 = t297 * t260 + t299 * t263;
t372 = pkin(4) * t207;
t209 = t299 * t260 - t297 * t263;
t206 = pkin(4) * t209;
t370 = pkin(5) * t222;
t281 = t296 * t337;
t333 = t298 * qJDD(1);
t258 = -0.2e1 * t281 + t333;
t349 = t298 * t258;
t360 = t296 * t255;
t205 = -t349 + t360;
t264 = (t291 - t292) * t301;
t182 = t297 * t205 + t299 * t264;
t183 = t299 * t205 - t297 * t264;
t362 = t292 * t301;
t273 = -t300 + t362;
t226 = -t298 * t273 + t357;
t382 = t297 * t226 + t298 * t332;
t381 = t299 * t226 - t297 * t333;
t343 = pkin(1) * t263 + pkin(5) * t260;
t219 = t296 * t273 + t348;
t274 = -t300 - t362;
t220 = t296 * t274 + t248;
t376 = pkin(1) * t220;
t380 = -qJ(3) * t274 - t376;
t336 = (qJD(3) * qJD(2));
t286 = -2 * t336;
t379 = -qJ(3) * t267 + t286 - t383;
t257 = -t281 + t333;
t378 = t257 * pkin(3) - qJ(4) * t362 + qJDD(4);
t269 = t299 * g(1) + t297 * g(2);
t242 = -t301 * pkin(1) + qJDD(1) * pkin(5) - t269;
t315 = -pkin(2) * t298 - qJ(3) * t296;
t253 = t315 * qJD(1);
t320 = qJD(1) * t253 + t242;
t322 = qJDD(2) * pkin(2) + t300 * qJ(3) - qJDD(3);
t368 = t298 * g(3);
t305 = t320 * t296 - t322 + t368;
t377 = pkin(2) + pkin(3);
t358 = t296 * t266;
t225 = t298 * t274 - t358;
t188 = t297 * t225 + t299 * t258;
t374 = pkin(4) * t188;
t371 = pkin(5) * t220;
t369 = t257 * pkin(2);
t367 = qJ(3) * t263;
t364 = qJ(3) * t298;
t268 = t297 * g(1) - t299 * g(2);
t241 = qJDD(1) * pkin(1) + t301 * pkin(5) + t268;
t361 = t296 * t241;
t359 = t296 * t258;
t351 = t298 * t241;
t350 = t298 * t255;
t345 = pkin(1) * t258 + pkin(5) * t225;
t342 = t296 * g(3) - t298 * t242;
t340 = qJD(1) * t296;
t339 = qJD(1) * t298;
t338 = qJD(2) * t298;
t331 = 0.2e1 * t340;
t329 = (t257 + t258) * pkin(2);
t328 = qJD(4) * t339;
t231 = t296 * t242 + t368;
t177 = t296 * t231 - t298 * t342;
t216 = -t297 * t268 - t299 * t269;
t265 = -qJD(2) * pkin(3) - qJ(4) * t340;
t321 = ((2 * qJD(3)) + t265) * t296;
t319 = qJD(3) * t331;
t318 = t297 * t330;
t317 = t299 * t330;
t262 = -t297 * t301 + t332;
t316 = -pkin(4) * t262 - t297 * g(3);
t314 = pkin(2) * t296 - t364;
t256 = t325 + t335;
t313 = t256 + t325;
t312 = t300 * pkin(2) - qJDD(2) * qJ(3) - t253 * t339 + t342;
t176 = t298 * t231 + t296 * t342;
t311 = t350 + t359;
t221 = -t298 * t272 + t358;
t215 = t299 * t268 - t297 * t269;
t285 = 2 * t336;
t178 = t285 - t312;
t309 = pkin(3) * t266 + t256 * qJ(4) + t322;
t308 = pkin(3) * t362 - qJD(2) * t265 + t312;
t307 = -pkin(2) * t281 + t241;
t306 = t257 * qJ(4) + t308;
t304 = t369 + (t255 + t313) * qJ(3) + t307;
t303 = t313 * qJ(3) + t307 + t319;
t167 = (qJ(4) * t338 + (-0.2e1 * qJD(4) + t253) * t296) * qJD(1) - t309 + t231;
t302 = t256 * qJ(3) + (qJ(3) * t338 + t321) * qJD(1) + t307 + t378;
t279 = 0.2e1 * t328;
t261 = t299 * t301 + t334;
t251 = t314 * qJDD(1);
t247 = t341 * t337;
t240 = -pkin(4) * t261 + t299 * g(3);
t238 = (-t377 * t296 + t364) * qJDD(1);
t236 = t297 * qJDD(2) + t299 * t247;
t235 = t298 * t256 - t291 * t337;
t234 = -t299 * qJDD(2) + t297 * t247;
t233 = -t296 * t257 - t292 * t337;
t218 = t313 * t296;
t217 = (t257 - t281) * t298;
t214 = qJ(3) * t258 + qJ(4) * t266;
t202 = t299 * t235 - t318;
t201 = t299 * t233 + t318;
t200 = t297 * t235 + t317;
t199 = t297 * t233 - t317;
t194 = -qJ(4) * t267 + t377 * t255;
t191 = t299 * t225 - t297 * t258;
t186 = pkin(4) * t191;
t185 = -t351 + t370;
t184 = -t361 - t371;
t180 = -t342 + t383;
t179 = t231 - t376;
t174 = t305 + t367;
t173 = pkin(2) * t263 + t178;
t172 = t303 + t369;
t171 = t329 + t303;
t170 = t304 + t319;
t169 = t299 * t177 - t297 * t241;
t168 = t297 * t177 + t299 * t241;
t166 = t285 - t306 - 0.2e1 * t328;
t165 = -pkin(2) * t266 + t305 + t380;
t164 = -pkin(2) * t271 + t312 + t379;
t163 = t302 + t369;
t162 = -t367 + qJD(4) * t331 + (-qJ(4) * t337 - g(3)) * t298 + (qJ(4) * qJDD(1) - t320) * t296 + t309;
t161 = t298 * t178 + t296 * t305;
t160 = t296 * t178 - t298 * t305;
t159 = t279 + t286 - t377 * t263 + (t257 + t333) * qJ(4) + t308;
t158 = qJ(4) * t271 + qJD(1) * t321 + t304 + t378;
t157 = -pkin(2) * t360 + t298 * t170 - t370;
t156 = qJ(3) * t349 - t296 * t171 - t371;
t155 = -t296 * t173 + t298 * t174;
t154 = pkin(3) * t258 - qJ(4) * t274 + t302 + t329;
t153 = -t377 * t271 + t279 + t306 + t379;
t152 = -t377 * t266 + t167 + t380;
t151 = t299 * t161 - t297 * t172;
t150 = t297 * t161 + t299 * t172;
t149 = t298 * t166 + t296 * t167;
t148 = t296 * t166 - t298 * t167;
t147 = qJ(3) * t163 - qJ(4) * t167;
t146 = -t296 * t154 + t298 * t214 - t371;
t145 = t298 * t158 - t296 * t194 - t370;
t144 = -pkin(1) * t160 + pkin(2) * t305 - qJ(3) * t178;
t143 = -t296 * t159 + t298 * t162;
t142 = -pkin(5) * t160 - t314 * t172;
t141 = t299 * t149 - t297 * t163;
t140 = t297 * t149 + t299 * t163;
t139 = -qJ(4) * t166 + t377 * t163;
t138 = -pkin(1) * t148 - qJ(3) * t166 + t377 * t167;
t137 = -pkin(5) * t148 - t296 * t139 + t298 * t147;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t261, -t262, 0, t216, 0, 0, 0, 0, 0, 0, t191, -t192, t209, t169, 0, 0, 0, 0, 0, 0, t191, t209, t192, t151, 0, 0, 0, 0, 0, 0, t191, t192, -t209, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t262, -t261, 0, t215, 0, 0, 0, 0, 0, 0, t188, -t189, t207, t168, 0, 0, 0, 0, 0, 0, t188, t207, t189, t150, 0, 0, 0, 0, 0, 0, t188, t189, -t207, t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t220, -t222, 0, -t176, 0, 0, 0, 0, 0, 0, t220, 0, t222, t160, 0, 0, 0, 0, 0, 0, t220, t222, 0, t148; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t262, 0, -t261, 0, t316, -t240, -t215, -pkin(4) * t215, t202, -t183, t198, t201, -t381, t236, -t297 * t179 + t299 * t184 - t374, -t297 * t180 + t299 * t185 + t373, t299 * t176 - t372, -pkin(4) * t168 - (pkin(1) * t297 - pkin(5) * t299) * t176, t202, t198, t183, t236, t381, t201, t299 * t156 - t297 * t165 - t374, t299 * t155 - t297 * t251 - t372, t299 * t157 - t297 * t164 - t373, -pkin(4) * t150 + t299 * t142 - t297 * t144, t202, t183, -t198, t201, -t381, t236, t299 * t146 - t297 * t152 - t374, t299 * t145 - t297 * t153 - t373, t299 * t143 - t297 * t238 + t372, -pkin(4) * t140 + t299 * t137 - t297 * t138; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t261, 0, t262, 0, t240, t316, t216, pkin(4) * t216, t200, -t182, t196, t199, -t382, t234, t299 * t179 + t297 * t184 + t186, t299 * t180 + t297 * t185 - t187, t297 * t176 + t206, pkin(4) * t169 - (-pkin(1) * t299 - pkin(5) * t297) * t176, t200, t196, t182, t234, t382, t199, t297 * t156 + t299 * t165 + t186, t297 * t155 + t299 * t251 + t206, t297 * t157 + t299 * t164 + t187, pkin(4) * t151 + t297 * t142 + t299 * t144, t200, t182, -t196, t199, -t382, t234, t297 * t146 + t299 * t152 + t186, t297 * t145 + t299 * t153 + t187, t297 * t143 + t299 * t238 - t206, pkin(4) * t141 + t297 * t137 + t299 * t138; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t268, t269, 0, 0, t218, t311, t221, t217, t219, 0, t345 + t351, -t344 - t361, t177 + t343, pkin(1) * t241 + pkin(5) * t177, t218, t221, -t311, 0, -t219, t217, qJ(3) * t359 + t298 * t171 + t345, t298 * t173 + t296 * t174 + t343, pkin(2) * t350 + t296 * t170 + t344, pkin(5) * t161 + (pkin(1) - t315) * t172, t218, -t311, -t221, t217, t219, 0, t298 * t154 + t296 * t214 + t345, t296 * t158 + t298 * t194 + t344, t298 * t159 + t296 * t162 - t343, pkin(1) * t163 + pkin(5) * t149 + t298 * t139 + t296 * t147;];
tauB_reg = t1;
