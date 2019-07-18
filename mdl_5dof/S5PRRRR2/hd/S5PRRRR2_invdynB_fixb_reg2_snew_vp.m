% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PRRRR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:20
% EndTime: 2019-07-18 13:30:26
% DurationCPUTime: 3.11s
% Computational Cost: add. (9149->293), mult. (12029->404), div. (0->0), fcn. (7518->8), ass. (0->191)
t284 = qJD(2) + qJD(3);
t279 = qJD(4) + t284;
t277 = t279 ^ 2;
t293 = cos(qJ(4));
t283 = qJDD(2) + qJDD(3);
t278 = qJDD(4) + t283;
t289 = sin(qJ(4));
t323 = t289 * t278;
t244 = t293 * t277 + t323;
t316 = t293 * t278;
t247 = t289 * t277 - t316;
t290 = sin(qJ(3));
t294 = cos(qJ(3));
t198 = t294 * t244 - t290 * t247;
t287 = g(3) - qJDD(1);
t229 = pkin(5) * t244 - t293 * t287;
t338 = pkin(5) * t247 - t289 * t287;
t156 = pkin(4) * t198 + t294 * t229 - t290 * t338;
t203 = t290 * t244 + t294 * t247;
t291 = sin(qJ(2));
t295 = cos(qJ(2));
t168 = t291 * t198 + t295 * t203;
t349 = pkin(4) * t203 + t290 * t229 + t294 * t338;
t356 = qJ(1) * t168 + t291 * t156 + t295 * t349;
t337 = t295 * t198 - t291 * t203;
t355 = qJ(1) * t337 + t295 * t156 - t291 * t349;
t271 = t291 * g(1) - t295 * g(2);
t261 = qJDD(2) * pkin(2) + t271;
t297 = qJD(2) ^ 2;
t302 = t295 * g(1) + t291 * g(2);
t262 = -t297 * pkin(2) - t302;
t217 = t290 * t261 + t294 * t262;
t282 = t284 ^ 2;
t208 = -t282 * pkin(3) + t217;
t216 = -t294 * t261 + t290 * t262;
t298 = t283 * pkin(3) - t216;
t173 = t289 * t208 - t293 * t298;
t174 = t293 * t208 + t289 * t298;
t307 = t289 * t173 + t293 * t174;
t137 = t293 * t173 - t289 * t174;
t314 = t294 * t137;
t115 = -t290 * t307 + t314;
t321 = t290 * t137;
t341 = t294 * t307 + t321;
t107 = t291 * t115 + t295 * t341;
t352 = t295 * t115 - t291 * t341;
t254 = t294 * t282 + t290 * t283;
t257 = t290 * t282 - t294 * t283;
t211 = t291 * t254 + t295 * t257;
t235 = pkin(4) * t254 - t294 * t287;
t339 = pkin(4) * t257 - t290 * t287;
t348 = qJ(1) * t211 + t291 * t235 + t295 * t339;
t299 = t295 * t254 - t291 * t257;
t347 = qJ(1) * t299 + t295 * t235 - t291 * t339;
t306 = t290 * t216 + t294 * t217;
t178 = t294 * t216 - t290 * t217;
t312 = t295 * t178;
t342 = -t291 * t306 + t312;
t319 = t291 * t178;
t140 = t295 * t306 + t319;
t170 = t278 * pkin(6) + t174;
t288 = sin(qJ(5));
t292 = cos(qJ(5));
t162 = t288 * t170 + t292 * t287;
t163 = t292 * t170 - t288 * t287;
t130 = t292 * t162 - t288 * t163;
t340 = t293 * t130;
t330 = pkin(1) * t287;
t329 = pkin(2) * t287;
t285 = t288 ^ 2;
t328 = t285 * t277;
t171 = -t277 * pkin(6) + t173;
t327 = t288 * t171;
t267 = t288 * t277 * t292;
t258 = qJDD(5) + t267;
t326 = t288 * t258;
t259 = qJDD(5) - t267;
t325 = t288 * t259;
t324 = t288 * t278;
t318 = t292 * t171;
t317 = t292 * t259;
t273 = t292 * t278;
t286 = t292 ^ 2;
t311 = t285 + t286;
t310 = qJD(5) * t279;
t309 = t288 * t310;
t308 = t292 * t310;
t131 = t288 * t162 + t292 * t163;
t231 = -t291 * t271 - t295 * t302;
t304 = t289 * t267;
t303 = t293 * t267;
t268 = t291 * qJDD(2) + t295 * t297;
t301 = qJ(1) * t268 - t295 * t287;
t269 = t295 * qJDD(2) - t291 * t297;
t300 = -qJ(1) * t269 - t291 * t287;
t230 = t295 * t271 - t291 * t302;
t296 = qJD(5) ^ 2;
t274 = t286 * t277;
t266 = -t274 - t296;
t265 = t274 - t296;
t264 = -t296 - t328;
t263 = t296 - t328;
t250 = t292 * t258;
t249 = t274 - t328;
t248 = t274 + t328;
t242 = t311 * t278;
t241 = t273 - 0.2e1 * t309;
t240 = t273 - t309;
t239 = t308 + t324;
t238 = 0.2e1 * t308 + t324;
t237 = t311 * t310;
t225 = t289 * qJDD(5) + t293 * t237;
t224 = -t293 * qJDD(5) + t289 * t237;
t223 = -t288 * t264 - t317;
t222 = -t288 * t263 + t250;
t221 = t292 * t266 - t326;
t220 = t292 * t265 - t325;
t219 = t292 * t264 - t325;
t218 = t288 * t266 + t250;
t215 = t292 * t239 - t285 * t310;
t214 = -t288 * t240 - t286 * t310;
t201 = t293 * t242 - t289 * t248;
t197 = t289 * t242 + t293 * t248;
t196 = -t288 * t238 + t292 * t241;
t195 = t293 * t222 + t288 * t323;
t194 = t293 * t220 + t289 * t273;
t193 = t289 * t222 - t288 * t316;
t192 = t289 * t220 - t292 * t316;
t191 = t293 * t215 - t304;
t190 = t293 * t214 + t304;
t189 = t289 * t215 + t303;
t188 = t289 * t214 - t303;
t187 = t293 * t223 + t289 * t238;
t186 = t293 * t221 - t289 * t241;
t185 = t289 * t223 - t293 * t238;
t184 = t289 * t221 + t293 * t241;
t183 = -t290 * t224 + t294 * t225;
t182 = t294 * t224 + t290 * t225;
t181 = t293 * t196 - t289 * t249;
t180 = t289 * t196 + t293 * t249;
t175 = pkin(4) * t306 + t329;
t167 = -t290 * t197 + t294 * t201;
t164 = t294 * t197 + t290 * t201;
t161 = -t290 * t193 + t294 * t195;
t160 = -t290 * t192 + t294 * t194;
t159 = t294 * t193 + t290 * t195;
t158 = t294 * t192 + t290 * t194;
t152 = -t290 * t189 + t294 * t191;
t151 = -t290 * t188 + t294 * t190;
t150 = t294 * t189 + t290 * t191;
t149 = t294 * t188 + t290 * t190;
t148 = -t290 * t185 + t294 * t187;
t147 = -t290 * t184 + t294 * t186;
t146 = t294 * t185 + t290 * t187;
t145 = t294 * t184 + t290 * t186;
t144 = -pkin(6) * t219 + t318;
t143 = -pkin(6) * t218 + t327;
t142 = -t290 * t180 + t294 * t181;
t141 = t294 * t180 + t290 * t181;
t134 = pkin(3) * t287 + pkin(5) * t307;
t133 = -t291 * t164 + t295 * t167;
t132 = t295 * t164 + t291 * t167;
t128 = -pkin(5) * t197 + t340;
t127 = pkin(5) * t201 + t289 * t130;
t126 = -t291 * t146 + t295 * t148;
t125 = -t291 * t145 + t295 * t147;
t124 = t295 * t146 + t291 * t148;
t123 = t295 * t145 + t291 * t147;
t122 = -pkin(5) * t185 + t293 * t144 - t289 * t163;
t121 = -pkin(5) * t184 + t293 * t143 - t289 * t162;
t120 = t293 * t131 + t289 * t171;
t119 = t289 * t131 - t293 * t171;
t118 = -pkin(3) * t219 + pkin(5) * t187 + t289 * t144 + t293 * t163;
t117 = -pkin(3) * t218 + pkin(5) * t186 + t289 * t143 + t293 * t162;
t112 = -pkin(5) * t119 + pkin(6) * t340;
t111 = -t290 * t119 + t294 * t120;
t110 = t294 * t119 + t290 * t120;
t109 = -pkin(4) * t164 - t290 * t127 + t294 * t128;
t108 = pkin(4) * t167 + t294 * t127 + t290 * t128;
t105 = pkin(4) * t115 + pkin(5) * t314 - t290 * t134;
t104 = pkin(4) * t341 + pkin(5) * t321 + t294 * t134 + t329;
t103 = pkin(5) * t120 - (-pkin(6) * t289 - pkin(3)) * t130;
t102 = -pkin(4) * t146 - t290 * t118 + t294 * t122;
t101 = -pkin(4) * t145 - t290 * t117 + t294 * t121;
t100 = -pkin(2) * t219 + pkin(4) * t148 + t294 * t118 + t290 * t122;
t99 = -pkin(2) * t218 + pkin(4) * t147 + t294 * t117 + t290 * t121;
t98 = -t291 * t110 + t295 * t111;
t97 = t295 * t110 + t291 * t111;
t96 = -pkin(4) * t110 - t290 * t103 + t294 * t112;
t95 = pkin(2) * t130 + pkin(4) * t111 + t294 * t103 + t290 * t112;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t268, -t269, 0, t231, 0, 0, 0, 0, 0, 0, -t299, t211, 0, t140, 0, 0, 0, 0, 0, 0, -t337, t168, 0, t107, 0, 0, 0, 0, 0, 0, t125, t126, t133, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t269, -t268, 0, t230, 0, 0, 0, 0, 0, 0, -t211, -t299, 0, -t342, 0, 0, 0, 0, 0, 0, -t168, -t337, 0, -t352, 0, 0, 0, 0, 0, 0, t123, t124, t132, t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t287, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t287, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t287, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t287, 0, 0, 0, 0, 0, 0, t218, t219, 0, -t130; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, 0, -t287, g(2), qJ(1) * g(2), 0, 0, t269, 0, -t268, 0, t300, t301, -t230, -qJ(1) * t230, 0, 0, -t211, 0, -t299, 0, t348, t347, t342, pkin(4) * t312 + qJ(1) * t342 - t291 * t175, 0, 0, -t168, 0, -t337, 0, t356, t355, t352, qJ(1) * t352 - t291 * t104 + t295 * t105, -t291 * t150 + t295 * t152, -t291 * t141 + t295 * t142, -t291 * t159 + t295 * t161, -t291 * t149 + t295 * t151, -t291 * t158 + t295 * t160, -t291 * t182 + t295 * t183, -qJ(1) * t123 + t295 * t101 - t291 * t99, -qJ(1) * t124 - t291 * t100 + t295 * t102, -qJ(1) * t132 - t291 * t108 + t295 * t109, -qJ(1) * t97 - t291 * t95 + t295 * t96; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t287, 0, -g(1), -qJ(1) * g(1), 0, 0, t268, 0, t269, 0, -t301, t300, t231, qJ(1) * t231 + t330, 0, 0, t299, 0, -t211, 0, -t347, t348, t140, pkin(4) * t319 + qJ(1) * t140 + t295 * t175 + t330, 0, 0, t337, 0, -t168, 0, -t355, t356, t107, qJ(1) * t107 + t295 * t104 + t291 * t105 + t330, t295 * t150 + t291 * t152, t295 * t141 + t291 * t142, t295 * t159 + t291 * t161, t295 * t149 + t291 * t151, t295 * t158 + t291 * t160, t295 * t182 + t291 * t183, -pkin(1) * t218 + qJ(1) * t125 + t291 * t101 + t295 * t99, -pkin(1) * t219 + qJ(1) * t126 + t295 * t100 + t291 * t102, qJ(1) * t133 + t295 * t108 + t291 * t109, pkin(1) * t130 + qJ(1) * t98 + t291 * t96 + t295 * t95; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(1) * t269 + t271, -pkin(1) * t268 + t302, 0, pkin(1) * t230, 0, 0, 0, 0, 0, t283, -pkin(1) * t211 - pkin(2) * t257 - t216, -pkin(1) * t299 - pkin(2) * t254 - t217, 0, -pkin(1) * t342 - pkin(2) * t178, 0, 0, 0, 0, 0, t278, -pkin(1) * t168 - pkin(2) * t203 - pkin(3) * t247 - t173, -pkin(1) * t337 - pkin(2) * t198 - pkin(3) * t244 - t174, 0, -pkin(1) * t352 - pkin(2) * t115 - pkin(3) * t137, (t239 + t308) * t288, t292 * t238 + t288 * t241, t292 * t263 + t326, (t240 - t309) * t292, t288 * t265 + t317, 0, pkin(1) * t123 + pkin(2) * t145 + pkin(3) * t184 + pkin(6) * t221 - t318, pkin(1) * t124 + pkin(2) * t146 + pkin(3) * t185 + pkin(6) * t223 + t327, pkin(1) * t132 + pkin(2) * t164 + pkin(3) * t197 + pkin(6) * t242 + t131, pkin(1) * t97 + pkin(2) * t110 + pkin(3) * t119 + pkin(6) * t131;];
tauB_reg  = t1;
