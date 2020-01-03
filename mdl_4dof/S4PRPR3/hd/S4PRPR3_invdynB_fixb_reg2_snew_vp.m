% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PRPR3
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PRPR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:57
% EndTime: 2019-12-31 16:21:01
% DurationCPUTime: 3.17s
% Computational Cost: add. (7129->312), mult. (16223->478), div. (0->0), fcn. (11411->8), ass. (0->218)
t314 = sin(pkin(6));
t316 = cos(pkin(6));
t294 = t316 * g(1) + t314 * g(2);
t318 = sin(qJ(2));
t320 = cos(qJ(2));
t334 = t314 * g(1) - t316 * g(2);
t347 = t320 * t294 - t318 * t334;
t365 = qJD(2) ^ 2;
t377 = -t365 * pkin(2) + qJDD(2) * qJ(3) + (2 * qJD(2) * qJD(3)) - t347;
t330 = t318 * t294 + t320 * t334;
t332 = -t318 * t330 - t320 * t347;
t211 = t318 * t347 - t320 * t330;
t355 = t316 * t211;
t381 = -t314 * t332 + t355;
t361 = t314 * t211;
t175 = t316 * t332 + t361;
t340 = t318 * qJDD(2);
t292 = t320 * t365 + t340;
t339 = t320 * qJDD(2);
t343 = t318 * t365;
t293 = t339 - t343;
t248 = -t314 * t292 + t316 * t293;
t311 = g(3) - qJDD(1);
t269 = pkin(4) * t292 - t320 * t311;
t329 = -pkin(4) * t293 - t318 * t311;
t380 = -qJ(1) * t248 + t314 * t269 + t316 * t329;
t317 = sin(qJ(4));
t315 = cos(pkin(7));
t319 = cos(qJ(4));
t313 = sin(pkin(7));
t362 = t313 * t317;
t276 = (-t315 * t319 + t362) * qJD(2);
t326 = t313 * t319 + t315 * t317;
t278 = t326 * qJD(2);
t241 = t278 * t276;
t368 = qJDD(4) - t241;
t379 = t317 * t368;
t378 = t319 * t368;
t370 = t316 * t292 + t314 * t293;
t375 = qJ(1) * t370 + t316 * t269 - t314 * t329;
t309 = t313 ^ 2;
t323 = t315 ^ 2;
t367 = t365 * (t309 + t323);
t282 = t315 * t367;
t335 = t315 * t339;
t257 = -t318 * t282 + t335;
t259 = t320 * t282 + t315 * t340;
t216 = t316 * t257 - t314 * t259;
t374 = t314 * t257 + t316 * t259;
t306 = t323 * t365;
t344 = t309 * t365;
t289 = t306 + t344;
t369 = t326 * qJDD(2);
t252 = -t316 * t294 - t314 * t334;
t251 = -t314 * t294 + t316 * t334;
t272 = t276 ^ 2;
t273 = t278 ^ 2;
t312 = qJDD(2) * pkin(2);
t300 = t315 * t311;
t201 = -t300 + (pkin(3) * t365 * t315 - pkin(5) * qJDD(2) - t377) * t313;
t225 = -t313 * t311 + t377 * t315;
t341 = qJDD(2) * t315;
t202 = -pkin(3) * t306 + pkin(5) * t341 + t225;
t168 = -t319 * t201 + t317 * t202;
t169 = t317 * t201 + t319 * t202;
t136 = -t319 * t168 + t317 * t169;
t364 = t313 * t136;
t363 = t313 * t315;
t357 = t314 * t311;
t356 = t315 * t136;
t354 = t316 * t311;
t234 = -t365 * qJ(3) + qJDD(3) - t312 - t330;
t222 = -pkin(3) * t341 - t289 * pkin(5) + t234;
t353 = t317 * t222;
t232 = qJDD(4) + t241;
t352 = t317 * t232;
t351 = t318 * t234;
t350 = t319 * t222;
t349 = t319 * t232;
t348 = t320 * t234;
t346 = t276 * qJD(4);
t345 = t278 * qJD(4);
t337 = t318 * t241;
t336 = t320 * t241;
t333 = -t234 + t312;
t137 = t317 * t168 + t319 * t169;
t224 = t377 * t313 + t300;
t185 = t313 * t224 + t315 * t225;
t274 = qJDD(2) * t362 - t319 * t341;
t184 = t315 * t224 - t313 * t225;
t261 = t292 * t363;
t262 = t313 * t335 - t343 * t363;
t328 = t316 * t261 + t314 * t262;
t327 = t314 * t261 - t316 * t262;
t321 = qJD(4) ^ 2;
t305 = t323 * qJDD(2);
t304 = t309 * qJDD(2);
t290 = t306 - t344;
t288 = t305 - t304;
t287 = t305 + t304;
t281 = t313 * t367;
t265 = -t273 - t321;
t264 = -t273 + t321;
t263 = t272 - t321;
t258 = t320 * t281 + t313 * t340;
t255 = t318 * t281 - t313 * t339;
t246 = t320 * t288 - t318 * t290;
t245 = t320 * t287 - t318 * t289;
t244 = t318 * t288 + t320 * t290;
t243 = t318 * t287 + t320 * t289;
t239 = -t273 + t272;
t238 = t369 - t346;
t237 = t369 - 0.2e1 * t346;
t236 = -t274 - t345;
t235 = t274 + 0.2e1 * t345;
t230 = -t321 - t272;
t228 = (-t276 * t319 + t278 * t317) * qJD(4);
t227 = (-t276 * t317 - t278 * t319) * qJD(4);
t226 = -t272 - t273;
t220 = t319 * t238 - t317 * t345;
t218 = -t314 * t255 + t316 * t258;
t217 = t317 * t238 + t319 * t345;
t215 = t316 * t255 + t314 * t258;
t214 = -t317 * t236 + t319 * t346;
t213 = t319 * t236 + t317 * t346;
t208 = -t317 * t265 - t349;
t207 = -t317 * t264 + t378;
t206 = t319 * t263 - t352;
t205 = t319 * t265 - t352;
t204 = t319 * t264 + t379;
t203 = t317 * t263 + t349;
t199 = -t314 * t243 + t316 * t245;
t198 = t316 * t243 + t314 * t245;
t197 = pkin(1) * t311 + pkin(4) * t332;
t194 = -t319 * t235 - t317 * t237;
t193 = -t274 * t319 + t317 * t369;
t192 = -t317 * t235 + t319 * t237;
t191 = -t274 * t317 - t319 * t369;
t190 = t319 * t230 - t379;
t189 = t317 * t230 + t378;
t188 = -t313 * t227 + t315 * t228;
t187 = t318 * qJDD(4) + t320 * t188;
t186 = -t320 * qJDD(4) + t318 * t188;
t182 = -pkin(5) * t205 + t350;
t181 = -pkin(4) * t255 - t318 * t225 + t315 * t348;
t180 = -pkin(4) * t257 - t318 * t224 + t313 * t348;
t179 = pkin(4) * t258 + t320 * t225 + t315 * t351;
t178 = -pkin(4) * t259 + t320 * t224 + t313 * t351;
t177 = -t313 * t217 + t315 * t220;
t176 = -t313 * t213 + t315 * t214;
t173 = -t313 * t205 + t315 * t208;
t172 = -t313 * t204 + t315 * t207;
t171 = -t313 * t203 + t315 * t206;
t170 = t315 * t205 + t313 * t208;
t167 = -pkin(5) * t189 + t353;
t165 = -pkin(4) * t243 + t320 * t184;
t164 = pkin(4) * t245 + t318 * t184;
t163 = t320 * t185 + t351;
t162 = t318 * t185 - t348;
t161 = t320 * t172 + t318 * t369;
t160 = t320 * t171 - t318 * t274;
t159 = t318 * t172 - t320 * t369;
t158 = t318 * t171 + t320 * t274;
t157 = -pkin(3) * t237 + pkin(5) * t208 + t353;
t156 = -t313 * t192 + t315 * t194;
t155 = -t313 * t191 + t315 * t193;
t154 = t315 * t191 + t313 * t193;
t153 = -t313 * t189 + t315 * t190;
t152 = t315 * t189 + t313 * t190;
t151 = t320 * t177 + t337;
t150 = t320 * t176 - t337;
t149 = t318 * t177 - t336;
t148 = t318 * t176 + t336;
t147 = -pkin(3) * t235 + pkin(5) * t190 - t350;
t146 = t320 * t173 + t318 * t237;
t145 = t318 * t173 - t320 * t237;
t144 = t320 * t156 - t318 * t239;
t143 = t318 * t156 + t320 * t239;
t142 = t320 * t153 + t318 * t235;
t141 = t318 * t153 - t320 * t235;
t140 = t320 * t155 + t318 * t226;
t139 = t318 * t155 - t320 * t226;
t138 = -pkin(2) * t154 - pkin(3) * t191;
t135 = -t314 * t162 + t316 * t163;
t134 = t316 * t162 + t314 * t163;
t133 = -pkin(2) * t170 - pkin(3) * t205 + t169;
t132 = -pkin(3) * t222 + pkin(5) * t137;
t131 = -t314 * t145 + t316 * t146;
t130 = t316 * t145 + t314 * t146;
t129 = -pkin(5) * t191 - t136;
t128 = -pkin(2) * t152 - pkin(3) * t189 + t168;
t127 = -pkin(4) * t162 - (pkin(2) * t318 - qJ(3) * t320) * t184;
t126 = -pkin(3) * t226 + pkin(5) * t193 + t137;
t125 = -qJ(3) * t170 - t313 * t157 + t315 * t182;
t124 = -t314 * t141 + t316 * t142;
t123 = t316 * t141 + t314 * t142;
t122 = -t314 * t139 + t316 * t140;
t121 = t316 * t139 + t314 * t140;
t120 = -qJ(3) * t152 - t313 * t147 + t315 * t167;
t119 = pkin(4) * t163 - (-pkin(2) * t320 - qJ(3) * t318 - pkin(1)) * t184;
t118 = t315 * t137 - t364;
t117 = t313 * t137 + t356;
t116 = t320 * t118 + t318 * t222;
t115 = t318 * t118 - t320 * t222;
t114 = -pkin(2) * t117 - pkin(3) * t136;
t113 = -pkin(4) * t145 + t320 * t125 - t318 * t133;
t112 = -qJ(3) * t154 - t313 * t126 + t315 * t129;
t111 = -pkin(1) * t170 + pkin(4) * t146 + t318 * t125 + t320 * t133;
t110 = -pkin(4) * t141 + t320 * t120 - t318 * t128;
t109 = -pkin(5) * t356 - qJ(3) * t117 - t313 * t132;
t108 = -pkin(1) * t152 + pkin(4) * t142 + t318 * t120 + t320 * t128;
t107 = -t314 * t115 + t316 * t116;
t106 = t316 * t115 + t314 * t116;
t105 = -pkin(4) * t139 + t320 * t112 - t318 * t138;
t104 = -pkin(1) * t154 + pkin(4) * t140 + t318 * t112 + t320 * t138;
t103 = -pkin(4) * t115 + t320 * t109 - t318 * t114;
t102 = -pkin(1) * t117 + pkin(4) * t116 + t318 * t109 + t320 * t114;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t252, 0, 0, 0, 0, 0, 0, -t370, -t248, 0, t175, 0, 0, 0, 0, 0, 0, -t374, t218, t199, t135, 0, 0, 0, 0, 0, 0, t124, t131, t122, t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, 0, 0, 0, 0, 0, 0, t248, -t370, 0, -t381, 0, 0, 0, 0, 0, 0, t216, t215, t198, t134, 0, 0, 0, 0, 0, 0, t123, t130, t121, t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t311, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t311, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t184, 0, 0, 0, 0, 0, 0, t152, t170, t154, t117; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t357, -t354, -t251, -qJ(1) * t251, 0, 0, t248, 0, -t370, 0, t380, t375, t381, pkin(4) * t355 + qJ(1) * t381 - t314 * t197, -t327, -t314 * t244 + t316 * t246, t218, t327, t374, 0, -qJ(1) * t216 - t314 * t178 + t316 * t180, -qJ(1) * t215 - t314 * t179 + t316 * t181, -qJ(1) * t198 - t314 * t164 + t316 * t165, -qJ(1) * t134 - t314 * t119 + t316 * t127, -t314 * t149 + t316 * t151, -t314 * t143 + t316 * t144, -t314 * t159 + t316 * t161, -t314 * t148 + t316 * t150, -t314 * t158 + t316 * t160, -t314 * t186 + t316 * t187, -qJ(1) * t123 - t314 * t108 + t316 * t110, -qJ(1) * t130 - t314 * t111 + t316 * t113, -qJ(1) * t121 - t314 * t104 + t316 * t105, -qJ(1) * t106 - t314 * t102 + t316 * t103; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t354, -t357, t252, qJ(1) * t252, 0, 0, t370, 0, t248, 0, -t375, t380, t175, pkin(4) * t361 + qJ(1) * t175 + t316 * t197, t328, t316 * t244 + t314 * t246, t215, -t328, -t216, 0, -qJ(1) * t374 + t316 * t178 + t314 * t180, qJ(1) * t218 + t316 * t179 + t314 * t181, qJ(1) * t199 + t316 * t164 + t314 * t165, qJ(1) * t135 + t316 * t119 + t314 * t127, t316 * t149 + t314 * t151, t316 * t143 + t314 * t144, t316 * t159 + t314 * t161, t316 * t148 + t314 * t150, t316 * t158 + t314 * t160, t316 * t186 + t314 * t187, qJ(1) * t124 + t316 * t108 + t314 * t110, qJ(1) * t131 + t316 * t111 + t314 * t113, qJ(1) * t122 + t316 * t104 + t314 * t105, qJ(1) * t107 + t316 * t102 + t314 * t103; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t334, t294, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(1) * t293 + t330, -pkin(1) * t292 + t347, 0, -pkin(1) * t211, t304, 0.2e1 * t313 * t341, 0, t305, 0, 0, pkin(1) * t257 - qJ(3) * t282 + t333 * t315, pkin(1) * t255 + qJ(3) * t281 - t333 * t313, pkin(1) * t243 + pkin(2) * t289 + qJ(3) * t287 + t185, pkin(1) * t162 - pkin(2) * t234 + qJ(3) * t185, t315 * t217 + t313 * t220, t315 * t192 + t313 * t194, t315 * t204 + t313 * t207, t315 * t213 + t313 * t214, t315 * t203 + t313 * t206, t315 * t227 + t313 * t228, pkin(1) * t141 - pkin(2) * t235 + qJ(3) * t153 + t315 * t147 + t313 * t167, pkin(1) * t145 - pkin(2) * t237 + qJ(3) * t173 + t315 * t157 + t313 * t182, pkin(1) * t139 - pkin(2) * t226 + qJ(3) * t155 + t315 * t126 + t313 * t129, pkin(1) * t115 - pkin(2) * t222 - pkin(5) * t364 + qJ(3) * t118 + t315 * t132;];
tauB_reg = t1;
