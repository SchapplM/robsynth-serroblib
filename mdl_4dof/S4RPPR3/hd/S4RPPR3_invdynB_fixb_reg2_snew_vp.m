% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPPR3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPPR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:58
% EndTime: 2019-12-31 16:38:02
% DurationCPUTime: 3.21s
% Computational Cost: add. (8098->320), mult. (17700->485), div. (0->0), fcn. (11435->8), ass. (0->221)
t336 = sin(qJ(1));
t338 = cos(qJ(1));
t311 = t338 * g(1) + t336 * g(2);
t378 = qJD(1) ^ 2;
t300 = -pkin(1) * t378 - t311;
t332 = sin(pkin(6));
t334 = cos(pkin(6));
t310 = t336 * g(1) - t338 * g(2);
t344 = qJDD(1) * pkin(1) + t310;
t256 = t334 * t300 + t332 * t344;
t390 = -pkin(2) * t378 + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t256;
t255 = t332 * t300 - t334 * t344;
t351 = t332 * t255 + t334 * t256;
t212 = t334 * t255 - t332 * t256;
t365 = t338 * t212;
t394 = -t336 * t351 + t365;
t370 = t336 * t212;
t179 = t338 * t351 + t370;
t358 = t332 * qJDD(1);
t303 = t334 * t378 + t358;
t357 = t334 * qJDD(1);
t304 = -t332 * t378 + t357;
t261 = -t336 * t303 + t338 * t304;
t329 = g(3) - qJDD(2);
t278 = qJ(2) * t303 - t334 * t329;
t348 = -qJ(2) * t304 - t332 * t329;
t393 = -pkin(4) * t261 + t336 * t278 + t338 * t348;
t335 = sin(qJ(4));
t333 = cos(pkin(7));
t337 = cos(qJ(4));
t331 = sin(pkin(7));
t376 = t331 * t335;
t289 = (-t333 * t337 + t376) * qJD(1);
t345 = t331 * t337 + t333 * t335;
t291 = t345 * qJD(1);
t254 = t291 * t289;
t381 = qJDD(4) - t254;
t392 = t335 * t381;
t391 = t337 * t381;
t383 = t338 * t303 + t336 * t304;
t388 = pkin(4) * t383 + t338 * t278 - t336 * t348;
t327 = t331 ^ 2;
t341 = t333 ^ 2;
t380 = t378 * (t327 + t341);
t298 = t333 * t380;
t353 = t333 * t357;
t267 = -t332 * t298 + t353;
t269 = t334 * t298 + t333 * t358;
t226 = t338 * t267 - t336 * t269;
t387 = t336 * t267 + t338 * t269;
t324 = t341 * t378;
t362 = t327 * t378;
t306 = t324 + t362;
t382 = t345 * qJDD(1);
t285 = t289 ^ 2;
t286 = t291 ^ 2;
t330 = qJDD(1) * pkin(2);
t317 = t333 * t329;
t361 = t333 * t378;
t215 = -t317 + (pkin(3) * t361 - pkin(5) * qJDD(1) - t390) * t331;
t236 = -t331 * t329 + t390 * t333;
t359 = qJDD(1) * t333;
t216 = -pkin(3) * t324 + pkin(5) * t359 + t236;
t181 = -t337 * t215 + t335 * t216;
t182 = t335 * t215 + t337 * t216;
t148 = -t337 * t181 + t335 * t182;
t377 = t331 * t148;
t242 = -t378 * qJ(3) + qJDD(3) + t255 - t330;
t375 = t332 * t242;
t374 = t333 * t148;
t373 = t334 * t242;
t234 = -pkin(3) * t359 - t306 * pkin(5) + t242;
t372 = t335 * t234;
t247 = qJDD(4) + t254;
t371 = t335 * t247;
t367 = t337 * t234;
t366 = t337 * t247;
t364 = t289 * qJD(4);
t363 = t291 * qJD(4);
t355 = t332 * t254;
t354 = t334 * t254;
t352 = -t242 + t330;
t149 = t335 * t181 + t337 * t182;
t235 = t390 * t331 + t317;
t197 = t331 * t235 + t333 * t236;
t272 = -t336 * t310 - t338 * t311;
t309 = t338 * qJDD(1) - t336 * t378;
t349 = -pkin(4) * t309 - t336 * g(3);
t287 = qJDD(1) * t376 - t337 * t359;
t196 = t333 * t235 - t331 * t236;
t273 = t303 * t333 * t331;
t274 = (-t332 * t361 + t353) * t331;
t347 = t338 * t273 + t336 * t274;
t346 = t336 * t273 - t338 * t274;
t271 = t338 * t310 - t336 * t311;
t339 = qJD(4) ^ 2;
t322 = t341 * qJDD(1);
t321 = t327 * qJDD(1);
t308 = t336 * qJDD(1) + t338 * t378;
t307 = t324 - t362;
t302 = t322 - t321;
t301 = t322 + t321;
t297 = t331 * t380;
t284 = -pkin(4) * t308 + t338 * g(3);
t281 = -t286 - t339;
t280 = -t286 + t339;
t279 = t285 - t339;
t268 = t334 * t297 + t331 * t358;
t265 = t332 * t297 - t331 * t357;
t260 = t334 * t302 - t332 * t307;
t259 = t334 * t301 - t332 * t306;
t258 = t332 * t302 + t334 * t307;
t257 = t332 * t301 + t334 * t306;
t253 = -t286 + t285;
t252 = t382 - t364;
t251 = t382 - 0.2e1 * t364;
t250 = -t287 - t363;
t249 = t287 + 0.2e1 * t363;
t245 = -t339 - t285;
t241 = (-t289 * t337 + t291 * t335) * qJD(4);
t240 = (-t289 * t335 - t291 * t337) * qJD(4);
t238 = -t285 - t286;
t233 = t337 * t252 - t335 * t363;
t232 = t335 * t252 + t337 * t363;
t231 = -t335 * t250 + t337 * t364;
t230 = t337 * t250 + t335 * t364;
t227 = -t336 * t265 + t338 * t268;
t225 = t338 * t265 + t336 * t268;
t224 = -t335 * t281 - t366;
t223 = -t335 * t280 + t391;
t222 = t337 * t279 - t371;
t221 = t337 * t281 - t371;
t220 = t337 * t280 + t392;
t219 = t335 * t279 + t366;
t218 = -t336 * t257 + t338 * t259;
t217 = t338 * t257 + t336 * t259;
t207 = -t337 * t249 - t335 * t251;
t206 = -t287 * t337 + t335 * t382;
t205 = -t335 * t249 + t337 * t251;
t204 = -t287 * t335 - t337 * t382;
t203 = t337 * t245 - t392;
t202 = t335 * t245 + t391;
t201 = pkin(1) * t329 + qJ(2) * t351;
t200 = -t331 * t240 + t333 * t241;
t199 = t332 * qJDD(4) + t334 * t200;
t198 = -t334 * qJDD(4) + t332 * t200;
t194 = -pkin(5) * t221 + t367;
t193 = -t331 * t232 + t333 * t233;
t192 = -t331 * t230 + t333 * t231;
t191 = -t331 * t221 + t333 * t224;
t190 = -t331 * t220 + t333 * t223;
t189 = -t331 * t219 + t333 * t222;
t188 = t333 * t221 + t331 * t224;
t187 = -qJ(2) * t265 - t332 * t236 + t333 * t373;
t186 = -qJ(2) * t267 - t332 * t235 + t331 * t373;
t185 = qJ(2) * t268 + t334 * t236 + t333 * t375;
t184 = -qJ(2) * t269 + t334 * t235 + t331 * t375;
t183 = -pkin(5) * t202 + t372;
t177 = t334 * t190 + t332 * t382;
t176 = t334 * t189 - t332 * t287;
t175 = t332 * t190 - t334 * t382;
t174 = t332 * t189 + t334 * t287;
t173 = -qJ(2) * t257 + t334 * t196;
t172 = qJ(2) * t259 + t332 * t196;
t171 = -t331 * t205 + t333 * t207;
t170 = -t331 * t204 + t333 * t206;
t169 = t333 * t204 + t331 * t206;
t168 = -pkin(3) * t251 + pkin(5) * t224 + t372;
t167 = -t331 * t202 + t333 * t203;
t166 = t333 * t202 + t331 * t203;
t165 = t334 * t193 + t355;
t164 = t334 * t192 - t355;
t163 = t332 * t193 - t354;
t162 = t332 * t192 + t354;
t161 = t334 * t197 + t375;
t160 = t332 * t197 - t373;
t159 = t334 * t191 + t332 * t251;
t158 = t332 * t191 - t334 * t251;
t157 = -pkin(3) * t249 + pkin(5) * t203 - t367;
t156 = t334 * t171 - t332 * t253;
t155 = t332 * t171 + t334 * t253;
t154 = t334 * t167 + t332 * t249;
t153 = t332 * t167 - t334 * t249;
t152 = t334 * t170 + t332 * t238;
t151 = t332 * t170 - t334 * t238;
t150 = -pkin(2) * t169 - pkin(3) * t204;
t147 = -pkin(2) * t188 - pkin(3) * t221 + t182;
t146 = -t336 * t160 + t338 * t161;
t145 = t338 * t160 + t336 * t161;
t144 = -t336 * t158 + t338 * t159;
t143 = t338 * t158 + t336 * t159;
t142 = -pkin(3) * t234 + pkin(5) * t149;
t141 = -pkin(5) * t204 - t148;
t140 = -pkin(2) * t166 - pkin(3) * t202 + t181;
t139 = -t336 * t153 + t338 * t154;
t138 = t338 * t153 + t336 * t154;
t137 = -qJ(2) * t160 - (pkin(2) * t332 - qJ(3) * t334) * t196;
t136 = -qJ(3) * t188 - t331 * t168 + t333 * t194;
t135 = -pkin(3) * t238 + pkin(5) * t206 + t149;
t134 = -t336 * t151 + t338 * t152;
t133 = t338 * t151 + t336 * t152;
t132 = -qJ(3) * t166 - t331 * t157 + t333 * t183;
t131 = qJ(2) * t161 - (-pkin(2) * t334 - qJ(3) * t332 - pkin(1)) * t196;
t130 = t333 * t149 - t377;
t129 = t331 * t149 + t374;
t128 = t334 * t130 + t332 * t234;
t127 = t332 * t130 - t334 * t234;
t126 = -pkin(2) * t129 - pkin(3) * t148;
t125 = -qJ(2) * t158 + t334 * t136 - t332 * t147;
t124 = -qJ(3) * t169 - t331 * t135 + t333 * t141;
t123 = -pkin(1) * t188 + qJ(2) * t159 + t332 * t136 + t334 * t147;
t122 = -qJ(2) * t153 + t334 * t132 - t332 * t140;
t121 = -pkin(1) * t166 + qJ(2) * t154 + t332 * t132 + t334 * t140;
t120 = -pkin(5) * t374 - qJ(3) * t129 - t331 * t142;
t119 = -t336 * t127 + t338 * t128;
t118 = t338 * t127 + t336 * t128;
t117 = -qJ(2) * t151 + t334 * t124 - t332 * t150;
t116 = -pkin(1) * t169 + qJ(2) * t152 + t332 * t124 + t334 * t150;
t115 = -qJ(2) * t127 + t334 * t120 - t332 * t126;
t114 = -pkin(1) * t129 + qJ(2) * t128 + t332 * t120 + t334 * t126;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t308, -t309, 0, t272, 0, 0, 0, 0, 0, 0, -t383, -t261, 0, t179, 0, 0, 0, 0, 0, 0, -t387, t227, t218, t146, 0, 0, 0, 0, 0, 0, t139, t144, t134, t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t309, -t308, 0, t271, 0, 0, 0, 0, 0, 0, t261, -t383, 0, -t394, 0, 0, 0, 0, 0, 0, t226, t225, t217, t145, 0, 0, 0, 0, 0, 0, t138, t143, t133, t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t329, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t196, 0, 0, 0, 0, 0, 0, t166, t188, t169, t129; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t309, 0, -t308, 0, t349, -t284, -t271, -pkin(4) * t271, 0, 0, t261, 0, -t383, 0, t393, t388, t394, pkin(4) * t394 + qJ(2) * t365 - t336 * t201, -t346, -t336 * t258 + t338 * t260, t227, t346, t387, 0, -pkin(4) * t226 - t336 * t184 + t338 * t186, -pkin(4) * t225 - t336 * t185 + t338 * t187, -pkin(4) * t217 - t336 * t172 + t338 * t173, -pkin(4) * t145 - t336 * t131 + t338 * t137, -t336 * t163 + t338 * t165, -t336 * t155 + t338 * t156, -t336 * t175 + t338 * t177, -t336 * t162 + t338 * t164, -t336 * t174 + t338 * t176, -t336 * t198 + t338 * t199, -pkin(4) * t138 - t336 * t121 + t338 * t122, -pkin(4) * t143 - t336 * t123 + t338 * t125, -pkin(4) * t133 - t336 * t116 + t338 * t117, -pkin(4) * t118 - t336 * t114 + t338 * t115; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t308, 0, t309, 0, t284, t349, t272, pkin(4) * t272, 0, 0, t383, 0, t261, 0, -t388, t393, t179, pkin(4) * t179 + qJ(2) * t370 + t338 * t201, t347, t338 * t258 + t336 * t260, t225, -t347, -t226, 0, -pkin(4) * t387 + t338 * t184 + t336 * t186, pkin(4) * t227 + t338 * t185 + t336 * t187, pkin(4) * t218 + t338 * t172 + t336 * t173, pkin(4) * t146 + t338 * t131 + t336 * t137, t338 * t163 + t336 * t165, t338 * t155 + t336 * t156, t338 * t175 + t336 * t177, t338 * t162 + t336 * t164, t338 * t174 + t336 * t176, t338 * t198 + t336 * t199, pkin(4) * t139 + t338 * t121 + t336 * t122, pkin(4) * t144 + t338 * t123 + t336 * t125, pkin(4) * t134 + t338 * t116 + t336 * t117, pkin(4) * t119 + t338 * t114 + t336 * t115; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t310, t311, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t304 - t255, -pkin(1) * t303 - t256, 0, -pkin(1) * t212, t321, 0.2e1 * t331 * t359, 0, t322, 0, 0, pkin(1) * t267 - qJ(3) * t298 + t333 * t352, pkin(1) * t265 + qJ(3) * t297 - t331 * t352, pkin(1) * t257 + pkin(2) * t306 + qJ(3) * t301 + t197, pkin(1) * t160 - pkin(2) * t242 + qJ(3) * t197, t333 * t232 + t331 * t233, t333 * t205 + t331 * t207, t333 * t220 + t331 * t223, t333 * t230 + t331 * t231, t333 * t219 + t331 * t222, t333 * t240 + t331 * t241, pkin(1) * t153 - pkin(2) * t249 + qJ(3) * t167 + t333 * t157 + t331 * t183, pkin(1) * t158 - pkin(2) * t251 + qJ(3) * t191 + t333 * t168 + t331 * t194, pkin(1) * t151 - pkin(2) * t238 + qJ(3) * t170 + t333 * t135 + t331 * t141, pkin(1) * t127 - pkin(2) * t234 - pkin(5) * t377 + qJ(3) * t130 + t333 * t142;];
tauB_reg = t1;
