% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPRP3
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPRP3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP3_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:49
% EndTime: 2019-12-31 16:42:52
% DurationCPUTime: 1.97s
% Computational Cost: add. (3814->267), mult. (8010->360), div. (0->0), fcn. (4360->6), ass. (0->202)
t295 = sin(qJ(1));
t297 = cos(qJ(1));
t271 = t297 * g(1) + t295 * g(2);
t299 = qJD(1) ^ 2;
t253 = -t299 * pkin(1) - t271;
t291 = sin(pkin(6));
t292 = cos(pkin(6));
t270 = t295 * g(1) - t297 * g(2);
t301 = qJDD(1) * pkin(1) + t270;
t204 = t291 * t253 - t292 * t301;
t205 = t292 * t253 + t291 * t301;
t308 = t291 * t204 + t292 * t205;
t169 = t292 * t204 - t291 * t205;
t326 = t297 * t169;
t357 = -t295 * t308 + t326;
t332 = t295 * t169;
t129 = t297 * t308 + t332;
t260 = t291 * qJDD(1) + t292 * t299;
t261 = t292 * qJDD(1) - t291 * t299;
t210 = -t295 * t260 + t297 * t261;
t289 = g(3) - qJDD(2);
t238 = qJ(2) * t260 - t292 * t289;
t302 = -qJ(2) * t261 - t291 * t289;
t356 = -pkin(4) * t210 + t295 * t238 + t297 * t302;
t349 = t297 * t260 + t295 * t261;
t354 = pkin(4) * t349 + t297 * t238 - t295 * t302;
t294 = sin(qJ(3));
t296 = cos(qJ(3));
t276 = t296 * t299 * t294;
t320 = qJDD(3) + t276;
t353 = t320 * pkin(3);
t194 = -t299 * pkin(2) + qJDD(1) * pkin(5) + t205;
t180 = t294 * t194 + t296 * t289;
t323 = qJD(1) * qJD(3);
t311 = t296 * t323;
t322 = t294 * qJDD(1);
t256 = t311 + t322;
t245 = t256 * qJ(4);
t348 = -t245 - t180 + t353;
t298 = qJD(3) ^ 2;
t287 = t294 ^ 2;
t337 = t287 * t299;
t273 = -t298 - t337;
t269 = qJDD(3) - t276;
t333 = t294 * t269;
t226 = t296 * t273 - t333;
t346 = pkin(2) * t226;
t288 = t296 ^ 2;
t285 = t288 * t299;
t275 = -t285 - t298;
t334 = t294 * t320;
t228 = t296 * t275 - t334;
t312 = t294 * t323;
t321 = t296 * qJDD(1);
t258 = -0.2e1 * t312 + t321;
t187 = t291 * t228 + t292 * t258;
t189 = t292 * t228 - t291 * t258;
t147 = t297 * t187 + t295 * t189;
t345 = pkin(4) * t147;
t327 = t296 * t269;
t230 = -t294 * t273 - t327;
t255 = 0.2e1 * t311 + t322;
t188 = t291 * t230 - t292 * t255;
t190 = t292 * t230 + t291 * t255;
t148 = t297 * t188 + t295 * t190;
t344 = pkin(4) * t148;
t325 = t287 + t288;
t262 = t325 * qJDD(1);
t265 = t285 + t337;
t214 = t291 * t262 + t292 * t265;
t215 = t292 * t262 - t291 * t265;
t174 = t297 * t214 + t295 * t215;
t343 = pkin(4) * t174;
t251 = t296 * t320;
t224 = t294 * t275 + t251;
t342 = pkin(5) * t224;
t341 = pkin(5) * t226;
t340 = qJ(2) * t187;
t339 = qJ(2) * t188;
t338 = qJ(2) * t214;
t152 = (qJ(4) * qJD(3) * t296 - 0.2e1 * qJD(4) * t294) * qJD(1) + t348;
t336 = t294 * t152;
t193 = -qJDD(1) * pkin(2) - t299 * pkin(5) + t204;
t335 = t294 * t193;
t329 = t296 * t152;
t328 = t296 * t193;
t182 = t296 * t194 - t294 * t289;
t324 = qJD(1) * t294;
t319 = 0.2e1 * qJD(1) * qJD(4);
t318 = pkin(1) * t187 + pkin(2) * t258 + pkin(5) * t228;
t317 = pkin(1) * t188 - pkin(2) * t255 + pkin(5) * t230;
t315 = pkin(1) * t214 + pkin(2) * t265 + pkin(5) * t262;
t314 = t291 * t322;
t313 = t292 * t322;
t310 = -pkin(1) * t224 + qJ(2) * t189;
t309 = -pkin(1) * t226 + qJ(2) * t190;
t141 = t294 * t180 + t296 * t182;
t219 = -t295 * t270 - t297 * t271;
t306 = t291 * t276;
t305 = t292 * t276;
t164 = -pkin(2) * t224 + t180;
t264 = t297 * qJDD(1) - t295 * t299;
t304 = -pkin(4) * t264 - t295 * g(3);
t140 = t296 * t180 - t294 * t182;
t218 = t297 * t270 - t295 * t271;
t257 = -t312 + t321;
t267 = qJD(3) * pkin(3) - qJ(4) * t324;
t300 = t257 * qJ(4) - qJD(3) * t267 + t296 * t319 + t182;
t166 = -t257 * pkin(3) - qJ(4) * t285 + t267 * t324 + qJDD(4) + t193;
t278 = t294 * t319;
t274 = t285 - t298;
t272 = t298 - t337;
t266 = t285 - t337;
t263 = t295 * qJDD(1) + t297 * t299;
t250 = t325 * t323;
t239 = -pkin(4) * t263 + t297 * g(3);
t234 = t296 * t256 - t287 * t323;
t233 = -t294 * t257 - t288 * t323;
t232 = t291 * qJDD(3) + t292 * t250;
t231 = -t292 * qJDD(3) + t291 * t250;
t229 = -t294 * t272 + t251;
t227 = t296 * t274 - t333;
t225 = t296 * t272 + t334;
t223 = t294 * t274 + t327;
t222 = (t256 + t311) * t294;
t221 = (t257 - t312) * t296;
t220 = -pkin(3) * t255 - qJ(4) * t269;
t208 = qJ(2) * t215;
t207 = -t294 * t255 + t296 * t258;
t206 = t296 * t255 + t294 * t258;
t202 = t292 * t234 - t306;
t201 = t292 * t233 + t306;
t200 = t291 * t234 + t305;
t199 = t291 * t233 - t305;
t198 = t292 * t229 + t314;
t197 = t292 * t227 + t291 * t321;
t196 = t291 * t229 - t313;
t195 = t291 * t227 - t292 * t321;
t181 = t292 * t207 - t291 * t266;
t179 = t291 * t207 + t292 * t266;
t177 = -t295 * t231 + t297 * t232;
t176 = t297 * t231 + t295 * t232;
t175 = -t295 * t214 + t297 * t215;
t173 = pkin(4) * t175;
t172 = t328 - t341;
t171 = t335 - t342;
t165 = t182 - t346;
t163 = pkin(1) * t289 + qJ(2) * t308;
t162 = -t295 * t200 + t297 * t202;
t161 = -t295 * t199 + t297 * t201;
t160 = t297 * t200 + t295 * t202;
t159 = t297 * t199 + t295 * t201;
t158 = -t295 * t196 + t297 * t198;
t157 = -t295 * t195 + t297 * t197;
t156 = t297 * t196 + t295 * t198;
t155 = t297 * t195 + t295 * t197;
t154 = -qJ(4) * t273 + t166;
t153 = -pkin(3) * t285 + t300;
t151 = t278 + (-t311 + t322) * qJ(4) - t348;
t150 = -t295 * t188 + t297 * t190;
t149 = -t295 * t187 + t297 * t189;
t146 = pkin(3) * t258 + qJ(4) * t275 - t166;
t145 = pkin(4) * t150;
t144 = pkin(4) * t149;
t143 = qJ(4) * t321 + (t265 - t285) * pkin(3) + t300;
t142 = -t295 * t179 + t297 * t181;
t139 = t297 * t179 + t295 * t181;
t137 = -t346 + (-t273 - t285) * pkin(3) + t300;
t136 = -qJ(4) * t311 + t164 + t245 + t278 - 0.2e1 * t353;
t135 = t292 * t140 - t338;
t134 = t291 * t140 + t208;
t133 = -qJ(4) * t251 - t294 * t146 - t342;
t132 = t296 * t154 - t294 * t220 - t341;
t131 = t292 * t141 + t291 * t193;
t130 = t291 * t141 - t292 * t193;
t127 = -pkin(3) * t166 + qJ(4) * t153;
t126 = t296 * t153 - t336;
t125 = t294 * t153 + t329;
t124 = -t291 * t165 + t292 * t172 - t339;
t123 = -t291 * t164 + t292 * t171 - t340;
t122 = -t294 * t143 + t296 * t151;
t121 = t292 * t165 + t291 * t172 + t309;
t120 = t292 * t164 + t291 * t171 + t310;
t119 = -pkin(3) * t314 + t292 * t122 - t338;
t118 = pkin(3) * t313 + t291 * t122 + t208;
t117 = t292 * t126 + t291 * t166;
t116 = t291 * t126 - t292 * t166;
t115 = -pkin(2) * t125 - pkin(3) * t152;
t114 = -t295 * t130 + t297 * t131;
t113 = t297 * t130 + t295 * t131;
t112 = t292 * t132 - t291 * t137 - t339;
t111 = t292 * t133 - t291 * t136 - t340;
t110 = t291 * t132 + t292 * t137 + t309;
t109 = t291 * t133 + t292 * t136 + t310;
t108 = -qJ(2) * t130 - (pkin(2) * t291 - pkin(5) * t292) * t140;
t107 = -pkin(5) * t125 - qJ(4) * t329 - t294 * t127;
t106 = qJ(2) * t131 - (-pkin(2) * t292 - pkin(5) * t291 - pkin(1)) * t140;
t105 = -t295 * t116 + t297 * t117;
t104 = t297 * t116 + t295 * t117;
t103 = -qJ(2) * t116 + t292 * t107 - t291 * t115;
t102 = -pkin(1) * t125 + qJ(2) * t117 + t291 * t107 + t292 * t115;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t263, -t264, 0, t219, 0, 0, 0, 0, 0, 0, -t349, -t210, 0, t129, 0, 0, 0, 0, 0, 0, t149, t150, t175, t114, 0, 0, 0, 0, 0, 0, t149, t150, t175, t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t264, -t263, 0, t218, 0, 0, 0, 0, 0, 0, t210, -t349, 0, -t357, 0, 0, 0, 0, 0, 0, t147, t148, t174, t113, 0, 0, 0, 0, 0, 0, t147, t148, t174, t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t289, 0, 0, 0, 0, 0, 0, t224, t226, 0, -t140, 0, 0, 0, 0, 0, 0, t224, t226, 0, t125; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t264, 0, -t263, 0, t304, -t239, -t218, -pkin(4) * t218, 0, 0, t210, 0, -t349, 0, t356, t354, t357, pkin(4) * t357 + qJ(2) * t326 - t295 * t163, t162, t142, t158, t161, t157, t177, -t295 * t120 + t297 * t123 - t345, -t295 * t121 + t297 * t124 - t344, -t295 * t134 + t297 * t135 - t343, -pkin(4) * t113 - t295 * t106 + t297 * t108, t162, t142, t158, t161, t157, t177, -t295 * t109 + t297 * t111 - t345, -t295 * t110 + t297 * t112 - t344, -t295 * t118 + t297 * t119 - t343, -pkin(4) * t104 - t295 * t102 + t297 * t103; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t263, 0, t264, 0, t239, t304, t219, pkin(4) * t219, 0, 0, t349, 0, t210, 0, -t354, t356, t129, pkin(4) * t129 + qJ(2) * t332 + t297 * t163, t160, t139, t156, t159, t155, t176, t297 * t120 + t295 * t123 + t144, t297 * t121 + t295 * t124 + t145, t297 * t134 + t295 * t135 + t173, pkin(4) * t114 + t297 * t106 + t295 * t108, t160, t139, t156, t159, t155, t176, t297 * t109 + t295 * t111 + t144, t297 * t110 + t295 * t112 + t145, t297 * t118 + t295 * t119 + t173, pkin(4) * t105 + t297 * t102 + t295 * t103; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t270, t271, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t261 - t204, -pkin(1) * t260 - t205, 0, -pkin(1) * t169, t222, t206, t225, t221, t223, 0, t318 - t328, t317 + t335, t141 + t315, pkin(1) * t130 - pkin(2) * t193 + pkin(5) * t141, t222, t206, t225, t221, t223, 0, -qJ(4) * t334 + t296 * t146 + t318, t294 * t154 + t296 * t220 + t317, t296 * t143 + t294 * t151 + t315, pkin(1) * t116 - pkin(2) * t166 + pkin(5) * t126 - qJ(4) * t336 + t296 * t127;];
tauB_reg = t1;
