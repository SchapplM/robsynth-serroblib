% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PRPRR9
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PRPRR9_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR9_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:46
% EndTime: 2019-12-31 17:39:52
% DurationCPUTime: 3.41s
% Computational Cost: add. (9190->309), mult. (14246->430), div. (0->0), fcn. (8266->8), ass. (0->199)
t299 = g(3) - qJDD(1);
t304 = sin(qJ(4));
t294 = -qJD(2) + qJD(4);
t292 = t294 ^ 2;
t293 = qJDD(2) - qJDD(4);
t307 = cos(qJ(4));
t331 = t307 * t293;
t311 = t292 * t304 + t331;
t238 = pkin(6) * t311 + t299 * t304;
t305 = sin(qJ(2));
t308 = cos(qJ(2));
t336 = t304 * t293;
t317 = -t292 * t307 + t336;
t357 = t305 * t317 + t308 * t311;
t359 = -pkin(6) * t317 + t299 * t307;
t159 = -pkin(5) * t357 + t238 * t308 - t305 * t359;
t301 = sin(pkin(8));
t302 = cos(pkin(8));
t209 = t305 * t311 - t308 * t317;
t374 = -pkin(5) * t209 + t238 * t305 + t308 * t359;
t375 = -t209 * t301 + t302 * t357;
t383 = -qJ(1) * t375 + t159 * t302 - t301 * t374;
t167 = t209 * t302 + t301 * t357;
t382 = -qJ(1) * t167 + t159 * t301 + t302 * t374;
t355 = qJD(2) ^ 2;
t270 = qJDD(2) * t305 + t308 * t355;
t271 = qJDD(2) * t308 - t305 * t355;
t217 = t270 * t302 + t271 * t301;
t241 = pkin(5) * t271 + t299 * t305;
t324 = -pkin(5) * t270 + t299 * t308;
t177 = -qJ(1) * t217 - t241 * t301 + t302 * t324;
t291 = 2 * qJD(3) * qJD(2);
t295 = qJDD(2) * qJ(3);
t272 = g(1) * t301 - g(2) * t302;
t273 = g(1) * t302 + g(2) * t301;
t329 = -t272 * t305 + t273 * t308;
t314 = t291 + t295 - t329;
t354 = pkin(2) + pkin(3);
t204 = -t354 * t355 + t314;
t300 = qJDD(2) * pkin(2);
t318 = t272 * t308 + t273 * t305;
t312 = -qJDD(3) + t318;
t215 = -qJ(3) * t355 - t300 - t312;
t310 = -qJDD(2) * pkin(3) + t215;
t165 = t304 * t204 - t307 * t310;
t166 = t204 * t307 + t304 * t310;
t128 = t165 * t307 - t166 * t304;
t129 = t165 * t304 + t166 * t307;
t107 = t128 * t305 - t129 * t308;
t367 = t128 * t308 + t129 * t305;
t99 = -t107 * t301 + t302 * t367;
t380 = t107 * t302 + t301 * t367;
t322 = -t305 * t318 - t308 * t329;
t184 = t305 * t329 - t308 * t318;
t345 = t302 * t184;
t376 = -t301 * t322 + t345;
t350 = t301 * t184;
t137 = t302 * t322 + t350;
t320 = -t270 * t301 + t271 * t302;
t356 = qJ(1) * t320 + t241 * t302 + t301 * t324;
t206 = -pkin(2) * t355 + t314;
t174 = t206 * t305 - t215 * t308;
t323 = t206 * t308 + t215 * t305;
t133 = -t174 * t301 + t302 * t323;
t132 = t174 * t302 + t301 * t323;
t353 = pkin(1) * t299;
t303 = sin(qJ(5));
t297 = t303 ^ 2;
t351 = t297 * t292;
t346 = t301 * t299;
t341 = t302 * t299;
t155 = pkin(4) * t293 - pkin(7) * t292 + t165;
t340 = t303 * t155;
t306 = cos(qJ(5));
t278 = t303 * t292 * t306;
t268 = qJDD(5) + t278;
t339 = t303 * t268;
t269 = qJDD(5) - t278;
t338 = t303 * t269;
t337 = t303 * t293;
t334 = t306 * t155;
t333 = t306 * t268;
t332 = t306 * t269;
t280 = t306 * t293;
t156 = -pkin(4) * t292 - pkin(7) * t293 + t166;
t147 = t156 * t306 + t299 * t303;
t298 = t306 ^ 2;
t328 = t297 + t298;
t327 = qJD(5) * t294;
t326 = t303 * t327;
t325 = t306 * t327;
t146 = t156 * t303 - t299 * t306;
t232 = -t272 * t301 - t273 * t302;
t316 = t304 * t278;
t315 = t307 * t278;
t313 = -pkin(1) * t270 + t329;
t122 = t146 * t306 - t147 * t303;
t123 = t303 * t146 + t306 * t147;
t231 = t272 * t302 - t273 * t301;
t309 = qJD(5) ^ 2;
t282 = t298 * t292;
t277 = -t282 - t309;
t276 = t282 - t309;
t275 = -t309 - t351;
t274 = t309 - t351;
t267 = pkin(1) * t271;
t258 = t282 - t351;
t257 = t282 + t351;
t252 = t328 * t293;
t251 = -t280 - 0.2e1 * t326;
t250 = -t280 - t326;
t249 = t325 - t337;
t248 = 0.2e1 * t325 - t337;
t247 = t328 * t327;
t236 = qJDD(5) * t304 + t247 * t307;
t235 = qJDD(5) * t307 - t247 * t304;
t230 = t249 * t306 - t297 * t327;
t229 = -t250 * t303 - t298 * t327;
t228 = -t275 * t303 - t332;
t227 = -t274 * t303 + t333;
t226 = t277 * t306 - t339;
t225 = t276 * t306 - t338;
t221 = t275 * t306 - t338;
t220 = t277 * t303 + t333;
t211 = -t252 * t307 - t257 * t304;
t207 = -t252 * t304 + t257 * t307;
t205 = -t248 * t303 + t251 * t306;
t200 = t227 * t307 - t303 * t336;
t199 = t225 * t307 - t280 * t304;
t198 = -t227 * t304 - t303 * t331;
t197 = -t225 * t304 - t306 * t331;
t195 = t230 * t307 - t316;
t194 = t229 * t307 + t316;
t193 = -t230 * t304 - t315;
t192 = -t229 * t304 + t315;
t191 = t228 * t307 + t248 * t304;
t190 = t226 * t307 - t251 * t304;
t189 = t228 * t304 - t248 * t307;
t188 = t226 * t304 + t251 * t307;
t187 = -t235 * t305 + t236 * t308;
t186 = t235 * t308 + t236 * t305;
t183 = t205 * t307 - t258 * t304;
t180 = -t205 * t304 - t258 * t307;
t179 = pkin(5) * t322 + t353;
t172 = t207 * t305 + t211 * t308;
t171 = -t207 * t308 + t211 * t305;
t164 = -t198 * t305 + t200 * t308;
t163 = -t197 * t305 + t199 * t308;
t162 = t198 * t308 + t200 * t305;
t161 = t197 * t308 + t199 * t305;
t154 = -t193 * t305 + t195 * t308;
t153 = -t192 * t305 + t194 * t308;
t152 = t193 * t308 + t195 * t305;
t151 = t192 * t308 + t194 * t305;
t150 = -pkin(5) * t174 + (-pkin(2) * t305 + qJ(3) * t308) * t299;
t148 = pkin(5) * t323 + (pkin(2) * t308 + qJ(3) * t305 + pkin(1)) * t299;
t145 = t189 * t305 + t191 * t308;
t144 = t188 * t305 + t190 * t308;
t143 = -t189 * t308 + t191 * t305;
t142 = -t188 * t308 + t190 * t305;
t141 = -t180 * t305 + t183 * t308;
t140 = t180 * t308 + t183 * t305;
t139 = -pkin(7) * t221 + t334;
t138 = -pkin(7) * t220 + t340;
t135 = -pkin(4) * t221 + t147;
t134 = -pkin(4) * t220 + t146;
t131 = -t171 * t301 + t172 * t302;
t130 = t171 * t302 + t172 * t301;
t125 = pkin(6) * t128 + qJ(3) * t299;
t124 = -pkin(6) * t129 + t299 * t354;
t120 = -t143 * t301 + t145 * t302;
t119 = -t142 * t301 + t144 * t302;
t118 = t143 * t302 + t145 * t301;
t117 = t142 * t302 + t144 * t301;
t116 = -pkin(6) * t207 + t122 * t307;
t115 = -pkin(6) * t211 - t122 * t304;
t114 = t123 * t307 + t155 * t304;
t113 = t123 * t304 - t155 * t307;
t112 = -pkin(6) * t189 + qJ(3) * t221 - t135 * t304 + t139 * t307;
t111 = -pkin(6) * t188 + qJ(3) * t220 - t134 * t304 + t138 * t307;
t106 = -pkin(6) * t191 - t307 * t135 - t304 * t139 + t221 * t354;
t105 = -pkin(6) * t190 - t307 * t134 - t304 * t138 + t220 * t354;
t104 = -pkin(5) * t171 - t115 * t305 + t116 * t308;
t103 = pkin(5) * t172 + t115 * t308 + t116 * t305;
t102 = t113 * t305 + t114 * t308;
t101 = -t113 * t308 + t114 * t305;
t98 = -pkin(5) * t367 - t124 * t305 + t125 * t308;
t97 = -pkin(5) * t107 + t124 * t308 + t125 * t305 + t353;
t96 = -pkin(5) * t143 - t106 * t305 + t112 * t308;
t95 = -pkin(5) * t142 - t105 * t305 + t111 * t308;
t94 = pkin(1) * t221 + pkin(5) * t145 + t106 * t308 + t112 * t305;
t93 = pkin(1) * t220 + pkin(5) * t144 + t105 * t308 + t111 * t305;
t92 = -pkin(6) * t113 - (pkin(4) * t304 - pkin(7) * t307 + qJ(3)) * t122;
t91 = -pkin(6) * t114 - (pkin(4) * t307 + pkin(7) * t304 + t354) * t122;
t90 = -t101 * t301 + t102 * t302;
t89 = t101 * t302 + t102 * t301;
t88 = -pkin(5) * t101 - t305 * t91 + t308 * t92;
t87 = -pkin(1) * t122 + pkin(5) * t102 + t305 * t92 + t308 * t91;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, 0, 0, 0, 0, 0, 0, -t217, -t320, 0, t137, 0, 0, 0, 0, 0, 0, -t217, 0, t320, t133, 0, 0, 0, 0, 0, 0, -t167, t375, 0, -t380, 0, 0, 0, 0, 0, 0, t119, t120, t131, t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t231, 0, 0, 0, 0, 0, 0, t320, -t217, 0, -t376, 0, 0, 0, 0, 0, 0, t320, 0, t217, t132, 0, 0, 0, 0, 0, 0, t375, t167, 0, t99, 0, 0, 0, 0, 0, 0, t117, t118, t130, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, 0, 0, 0, 0, 0, 0, -t220, -t221, 0, t122; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t346, -t341, -t231, -qJ(1) * t231, 0, 0, t320, 0, -t217, 0, -t356, -t177, t376, pkin(5) * t345 + qJ(1) * t376 - t179 * t301, 0, t320, 0, 0, t217, 0, -t356, -t132, t177, -qJ(1) * t132 - t148 * t301 + t150 * t302, 0, 0, -t375, 0, -t167, 0, t383, t382, t99, -qJ(1) * t99 - t301 * t97 + t302 * t98, -t152 * t301 + t154 * t302, -t140 * t301 + t141 * t302, -t162 * t301 + t164 * t302, -t151 * t301 + t153 * t302, -t161 * t301 + t163 * t302, -t186 * t301 + t187 * t302, -qJ(1) * t117 - t301 * t93 + t302 * t95, -qJ(1) * t118 - t301 * t94 + t302 * t96, -qJ(1) * t130 - t103 * t301 + t104 * t302, -qJ(1) * t89 - t301 * t87 + t302 * t88; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t341, -t346, t232, qJ(1) * t232, 0, 0, t217, 0, t320, 0, t177, -t356, t137, pkin(5) * t350 + qJ(1) * t137 + t179 * t302, 0, t217, 0, 0, -t320, 0, t177, t133, t356, qJ(1) * t133 + t148 * t302 + t150 * t301, 0, 0, -t167, 0, t375, 0, t382, -t383, t380, -qJ(1) * t380 + t301 * t98 + t302 * t97, t152 * t302 + t154 * t301, t140 * t302 + t141 * t301, t162 * t302 + t164 * t301, t151 * t302 + t153 * t301, t161 * t302 + t163 * t301, t186 * t302 + t187 * t301, qJ(1) * t119 + t301 * t95 + t302 * t93, qJ(1) * t120 + t301 * t96 + t302 * t94, qJ(1) * t131 + t103 * t302 + t104 * t301, qJ(1) * t90 + t301 * t88 + t302 * t87; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t272, t273, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t267 + t318, t313, 0, -pkin(1) * t184, 0, 0, 0, qJDD(2), 0, 0, t267 + 0.2e1 * t300 + t312, 0, t291 + 0.2e1 * t295 - t313, pkin(1) * t174 - pkin(2) * t215 + qJ(3) * t206, 0, 0, 0, 0, 0, t293, pkin(1) * t357 + qJ(3) * t317 + t311 * t354 + t165, pkin(1) * t209 + qJ(3) * t311 - t317 * t354 + t166, 0, pkin(1) * t367 + qJ(3) * t129 + t128 * t354, (-t249 - t325) * t303, -t248 * t306 - t251 * t303, -t274 * t306 - t339, (-t250 + t326) * t306, -t276 * t303 - t332, 0, pkin(1) * t142 - pkin(4) * t251 - pkin(7) * t226 + qJ(3) * t190 - t188 * t354 + t334, pkin(1) * t143 + pkin(4) * t248 - pkin(7) * t228 + qJ(3) * t191 - t189 * t354 - t340, pkin(1) * t171 - pkin(4) * t257 + pkin(7) * t252 + qJ(3) * t211 - t207 * t354 - t123, pkin(1) * t101 + pkin(4) * t155 - pkin(7) * t123 + qJ(3) * t114 - t113 * t354;];
tauB_reg = t1;
