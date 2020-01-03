% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRPR3
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRPR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:35
% EndTime: 2019-12-31 17:01:39
% DurationCPUTime: 2.87s
% Computational Cost: add. (8256->277), mult. (12125->396), div. (0->0), fcn. (7328->8), ass. (0->185)
t287 = qJD(1) + qJD(2);
t285 = t287 ^ 2;
t286 = qJDD(1) + qJDD(2);
t291 = sin(pkin(7));
t292 = cos(pkin(7));
t249 = t285 * t292 + t286 * t291;
t252 = t285 * t291 - t286 * t292;
t294 = sin(qJ(2));
t297 = cos(qJ(2));
t204 = t249 * t297 - t252 * t294;
t290 = g(3) - qJDD(3);
t235 = qJ(3) * t249 - t290 * t292;
t336 = qJ(3) * t252 - t290 * t291;
t162 = pkin(5) * t204 + t235 * t297 - t294 * t336;
t208 = t249 * t294 + t252 * t297;
t295 = sin(qJ(1));
t298 = cos(qJ(1));
t173 = t204 * t295 + t208 * t298;
t344 = pkin(5) * t208 + t235 * t294 + t297 * t336;
t353 = pkin(4) * t173 + t162 * t295 + t298 * t344;
t335 = t204 * t298 - t208 * t295;
t352 = pkin(4) * t335 + t162 * t298 - t295 * t344;
t274 = g(1) * t295 - g(2) * t298;
t263 = qJDD(1) * pkin(1) + t274;
t275 = g(1) * t298 + g(2) * t295;
t300 = qJD(1) ^ 2;
t264 = -pkin(1) * t300 - t275;
t219 = t263 * t294 + t264 * t297;
t212 = -pkin(2) * t285 + t219;
t302 = t263 * t297 - t264 * t294;
t301 = pkin(2) * t286 + t302;
t177 = t212 * t291 - t292 * t301;
t178 = t212 * t292 + t291 * t301;
t309 = t177 * t291 + t178 * t292;
t137 = t177 * t292 - t178 * t291;
t315 = t297 * t137;
t118 = -t294 * t309 + t315;
t319 = t294 * t137;
t338 = t297 * t309 + t319;
t109 = t118 * t295 + t298 * t338;
t349 = t118 * t298 - t295 * t338;
t256 = t285 * t297 + t286 * t294;
t259 = t285 * t294 - t286 * t297;
t215 = t256 * t295 + t259 * t298;
t239 = pkin(5) * t256 - g(3) * t297;
t337 = pkin(5) * t259 - g(3) * t294;
t346 = pkin(4) * t215 + t239 * t295 + t298 * t337;
t303 = t256 * t298 - t259 * t295;
t345 = pkin(4) * t303 + t239 * t298 - t295 * t337;
t308 = t219 * t297 - t294 * t302;
t182 = -t219 * t294 - t297 * t302;
t314 = t298 * t182;
t339 = -t295 * t308 + t314;
t318 = t295 * t182;
t141 = t298 * t308 + t318;
t293 = sin(qJ(4));
t288 = t293 ^ 2;
t326 = t288 * t285;
t168 = -pkin(3) * t286 - pkin(6) * t285 + t177;
t323 = t293 * t168;
t296 = cos(qJ(4));
t273 = t296 * t285 * t293;
t265 = qJDD(4) + t273;
t322 = t293 * t265;
t266 = qJDD(4) - t273;
t321 = t293 * t266;
t320 = t293 * t286;
t317 = t296 * t168;
t316 = t296 * t266;
t277 = t296 * t286;
t169 = -pkin(3) * t285 + pkin(6) * t286 + t178;
t154 = t169 * t296 - t290 * t293;
t289 = t296 ^ 2;
t313 = t288 + t289;
t312 = qJD(4) * t287;
t311 = t293 * t312;
t310 = t296 * t312;
t153 = t169 * t293 + t290 * t296;
t132 = t153 * t293 + t154 * t296;
t231 = -t274 * t295 - t275 * t298;
t306 = t291 * t273;
t305 = t292 * t273;
t268 = qJDD(1) * t298 - t295 * t300;
t304 = -pkin(4) * t268 - g(3) * t295;
t131 = t153 * t296 - t154 * t293;
t230 = t274 * t298 - t275 * t295;
t299 = qJD(4) ^ 2;
t278 = t289 * t285;
t272 = -t278 - t299;
t271 = t278 - t299;
t270 = -t299 - t326;
t269 = t299 - t326;
t267 = qJDD(1) * t295 + t298 * t300;
t261 = t278 - t326;
t260 = t278 + t326;
t255 = t296 * t265;
t254 = t313 * t286;
t247 = t277 - 0.2e1 * t311;
t246 = t277 - t311;
t245 = t310 + t320;
t244 = 0.2e1 * t310 + t320;
t243 = -pkin(4) * t267 + g(3) * t298;
t242 = t313 * t312;
t229 = qJDD(4) * t291 + t242 * t292;
t228 = -qJDD(4) * t292 + t242 * t291;
t227 = t245 * t296 - t288 * t312;
t226 = -t246 * t293 - t289 * t312;
t225 = -t270 * t293 - t316;
t224 = -t269 * t293 + t255;
t223 = t272 * t296 - t322;
t222 = t271 * t296 - t321;
t221 = t270 * t296 - t321;
t220 = t272 * t293 + t255;
t211 = t254 * t292 - t260 * t291;
t210 = t254 * t291 + t260 * t292;
t203 = -t244 * t293 + t247 * t296;
t199 = t224 * t292 + t291 * t320;
t198 = t222 * t292 + t277 * t291;
t197 = t224 * t291 - t292 * t320;
t196 = t222 * t291 - t277 * t292;
t195 = t227 * t292 - t306;
t194 = t226 * t292 + t306;
t193 = t227 * t291 + t305;
t192 = t226 * t291 - t305;
t191 = t225 * t292 + t244 * t291;
t190 = t223 * t292 - t247 * t291;
t189 = t225 * t291 - t244 * t292;
t188 = t223 * t291 + t247 * t292;
t187 = -t228 * t294 + t229 * t297;
t186 = t228 * t297 + t229 * t294;
t185 = t203 * t292 - t261 * t291;
t184 = t203 * t291 + t261 * t292;
t179 = pkin(1) * g(3) + pkin(5) * t308;
t176 = -t210 * t294 + t211 * t297;
t175 = t210 * t297 + t211 * t294;
t167 = -t197 * t294 + t199 * t297;
t166 = -t196 * t294 + t198 * t297;
t165 = t197 * t297 + t199 * t294;
t164 = t196 * t297 + t198 * t294;
t158 = -t193 * t294 + t195 * t297;
t157 = -t192 * t294 + t194 * t297;
t156 = t193 * t297 + t195 * t294;
t155 = t192 * t297 + t194 * t294;
t151 = -t189 * t294 + t191 * t297;
t150 = -t188 * t294 + t190 * t297;
t149 = t189 * t297 + t191 * t294;
t148 = t188 * t297 + t190 * t294;
t147 = -pkin(6) * t221 + t317;
t146 = -pkin(6) * t220 + t323;
t145 = -t184 * t294 + t185 * t297;
t144 = t184 * t297 + t185 * t294;
t143 = -pkin(3) * t221 + t154;
t142 = -pkin(3) * t220 + t153;
t139 = -t175 * t295 + t176 * t298;
t136 = t175 * t298 + t176 * t295;
t133 = pkin(2) * t290 + qJ(3) * t309;
t129 = -t149 * t295 + t151 * t298;
t128 = -t148 * t295 + t150 * t298;
t127 = t149 * t298 + t151 * t295;
t126 = t148 * t298 + t150 * t295;
t125 = -qJ(3) * t210 + t131 * t292;
t124 = qJ(3) * t211 + t131 * t291;
t123 = t132 * t292 + t168 * t291;
t122 = t132 * t291 - t168 * t292;
t121 = -qJ(3) * t189 - t143 * t291 + t147 * t292;
t120 = -qJ(3) * t188 - t142 * t291 + t146 * t292;
t115 = -pkin(2) * t221 + qJ(3) * t191 + t143 * t292 + t147 * t291;
t114 = -pkin(2) * t220 + qJ(3) * t190 + t142 * t292 + t146 * t291;
t113 = -pkin(5) * t175 - t124 * t294 + t125 * t297;
t112 = pkin(5) * t176 + t124 * t297 + t125 * t294;
t111 = -t122 * t294 + t123 * t297;
t110 = t122 * t297 + t123 * t294;
t107 = pkin(5) * t118 + qJ(3) * t315 - t133 * t294;
t106 = pkin(1) * t290 + pkin(5) * t338 + qJ(3) * t319 + t133 * t297;
t105 = -qJ(3) * t122 - (pkin(3) * t291 - pkin(6) * t292) * t131;
t104 = -pkin(5) * t149 - t115 * t294 + t121 * t297;
t103 = -pkin(5) * t148 - t114 * t294 + t120 * t297;
t102 = -pkin(1) * t221 + pkin(5) * t151 + t115 * t297 + t121 * t294;
t101 = -pkin(1) * t220 + pkin(5) * t150 + t114 * t297 + t120 * t294;
t100 = qJ(3) * t123 - (-pkin(3) * t292 - pkin(6) * t291 - pkin(2)) * t131;
t99 = -t110 * t295 + t111 * t298;
t98 = t110 * t298 + t111 * t295;
t97 = -pkin(5) * t110 - t100 * t294 + t105 * t297;
t96 = pkin(1) * t131 + pkin(5) * t111 + t100 * t297 + t105 * t294;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t267, -t268, 0, t231, 0, 0, 0, 0, 0, 0, -t303, t215, 0, t141, 0, 0, 0, 0, 0, 0, -t335, t173, 0, t109, 0, 0, 0, 0, 0, 0, t128, t129, t139, t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t268, -t267, 0, t230, 0, 0, 0, 0, 0, 0, -t215, -t303, 0, -t339, 0, 0, 0, 0, 0, 0, -t173, -t335, 0, -t349, 0, 0, 0, 0, 0, 0, t126, t127, t136, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t290, 0, 0, 0, 0, 0, 0, t220, t221, 0, -t131; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t268, 0, -t267, 0, t304, -t243, -t230, -pkin(4) * t230, 0, 0, -t215, 0, -t303, 0, t346, t345, t339, pkin(4) * t339 + pkin(5) * t314 - t179 * t295, 0, 0, -t173, 0, -t335, 0, t353, t352, t349, pkin(4) * t349 - t106 * t295 + t107 * t298, -t156 * t295 + t158 * t298, -t144 * t295 + t145 * t298, -t165 * t295 + t167 * t298, -t155 * t295 + t157 * t298, -t164 * t295 + t166 * t298, -t186 * t295 + t187 * t298, -pkin(4) * t126 - t101 * t295 + t103 * t298, -pkin(4) * t127 - t102 * t295 + t104 * t298, -pkin(4) * t136 - t112 * t295 + t113 * t298, -pkin(4) * t98 - t295 * t96 + t298 * t97; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t267, 0, t268, 0, t243, t304, t231, pkin(4) * t231, 0, 0, t303, 0, -t215, 0, -t345, t346, t141, pkin(4) * t141 + pkin(5) * t318 + t179 * t298, 0, 0, t335, 0, -t173, 0, -t352, t353, t109, pkin(4) * t109 + t106 * t298 + t107 * t295, t156 * t298 + t158 * t295, t144 * t298 + t145 * t295, t165 * t298 + t167 * t295, t155 * t298 + t157 * t295, t164 * t298 + t166 * t295, t186 * t298 + t187 * t295, pkin(4) * t128 + t101 * t298 + t103 * t295, pkin(4) * t129 + t102 * t298 + t104 * t295, pkin(4) * t139 + t112 * t298 + t113 * t295, pkin(4) * t99 + t295 * t97 + t298 * t96; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t274, t275, 0, 0, 0, 0, 0, 0, 0, t286, -pkin(1) * t259 + t302, -pkin(1) * t256 - t219, 0, -pkin(1) * t182, 0, 0, 0, 0, 0, t286, -pkin(1) * t208 - pkin(2) * t252 - t177, -pkin(1) * t204 - pkin(2) * t249 - t178, 0, -pkin(1) * t118 - pkin(2) * t137, (t245 + t310) * t293, t244 * t296 + t247 * t293, t269 * t296 + t322, (t246 - t311) * t296, t271 * t293 + t316, 0, pkin(1) * t148 + pkin(2) * t188 + pkin(3) * t247 + pkin(6) * t223 - t317, pkin(1) * t149 + pkin(2) * t189 - pkin(3) * t244 + pkin(6) * t225 + t323, pkin(1) * t175 + pkin(2) * t210 + pkin(3) * t260 + pkin(6) * t254 + t132, pkin(1) * t110 + pkin(2) * t122 - pkin(3) * t168 + pkin(6) * t132;];
tauB_reg = t1;
