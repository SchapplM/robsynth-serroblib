% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPRR8
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPRR8_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR8_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:18
% EndTime: 2019-12-31 16:55:21
% DurationCPUTime: 2.41s
% Computational Cost: add. (5204->293), mult. (10483->415), div. (0->0), fcn. (6460->6), ass. (0->202)
t288 = sin(qJ(4));
t283 = qJDD(3) + qJDD(4);
t289 = sin(qJ(3));
t291 = cos(qJ(4));
t292 = cos(qJ(3));
t253 = (-t288 * t292 - t289 * t291) * qJD(1);
t320 = qJD(1) * t292;
t255 = -t288 * t289 * qJD(1) + t291 * t320;
t340 = t255 * t253;
t345 = t283 + t340;
t347 = t288 * t345;
t346 = t291 * t345;
t251 = t253 ^ 2;
t252 = t255 ^ 2;
t284 = qJD(3) + qJD(4);
t282 = t284 ^ 2;
t344 = pkin(5) + pkin(1);
t343 = qJDD(1) * pkin(1);
t318 = qJD(1) * qJD(3);
t309 = t289 * t318;
t314 = t292 * qJDD(1);
t262 = -t309 + t314;
t295 = qJD(1) ^ 2;
t323 = t292 * t295;
t290 = sin(qJ(1));
t293 = cos(qJ(1));
t271 = t290 * g(1) - t293 * g(2);
t302 = qJDD(2) - t271;
t298 = -t295 * qJ(2) + t302;
t243 = -t344 * qJDD(1) + t298;
t326 = t292 * t243;
t202 = qJDD(3) * pkin(3) - t262 * pkin(6) + t326 + (-pkin(3) * t323 - pkin(6) * t318 + g(3)) * t289;
t224 = -t292 * g(3) + t289 * t243;
t308 = t292 * t318;
t316 = t289 * qJDD(1);
t261 = -t308 - t316;
t299 = -qJD(3) * pkin(3) + pkin(6) * t320;
t286 = t289 ^ 2;
t337 = t286 * t295;
t203 = -pkin(3) * t337 + t261 * pkin(6) + qJD(3) * t299 + t224;
t171 = -t291 * t202 + t288 * t203;
t172 = t288 * t202 + t291 * t203;
t146 = -t291 * t171 + t288 * t172;
t342 = t146 * t289;
t341 = t146 * t292;
t339 = t284 * t288;
t338 = t284 * t291;
t287 = t292 ^ 2;
t336 = t287 * t295;
t272 = t293 * g(1) + t290 * g(2);
t285 = qJDD(1) * qJ(2);
t300 = t272 - t285;
t317 = qJD(2) * qJD(1);
t296 = t300 - 0.2e1 * t317;
t204 = t261 * pkin(3) + t299 * t320 + (pkin(6) * t286 + t344) * t295 + t296;
t335 = t288 * t204;
t218 = -t340 + t283;
t334 = t288 * t218;
t238 = t344 * t295 + t296;
t333 = t289 * t238;
t310 = t289 * t323;
t269 = qJDD(3) + t310;
t332 = t289 * t269;
t270 = qJDD(3) - t310;
t331 = t289 * t270;
t321 = t286 + t287;
t264 = t321 * qJDD(1);
t330 = t290 * t264;
t329 = t291 * t204;
t328 = t291 * t218;
t327 = t292 * t238;
t325 = t292 * t269;
t324 = t292 * t270;
t322 = t293 * t264;
t319 = qJD(4) + t284;
t315 = t290 * qJDD(1);
t313 = t293 * qJDD(1);
t312 = t290 * t340;
t311 = t293 * t340;
t147 = t288 * t171 + t291 * t172;
t280 = 0.2e1 * t317;
t245 = -t295 * pkin(1) + t280 - t300;
t246 = -t298 + t343;
t213 = t293 * t245 - t290 * t246;
t307 = -t291 * t261 + t288 * t262;
t228 = -t290 * t271 - t293 * t272;
t306 = t290 * t310;
t305 = t293 * t310;
t265 = -t290 * t295 + t313;
t304 = pkin(4) * t265 + t290 * g(3);
t266 = t293 * t295 + t315;
t303 = -pkin(4) * t266 + t293 * g(3);
t223 = t289 * g(3) + t326;
t195 = t292 * t223 + t289 * t224;
t196 = -t289 * t223 + t292 * t224;
t210 = t290 * t245 + t293 * t246;
t301 = t288 * t261 + t291 * t262;
t227 = t293 * t271 - t290 * t272;
t297 = (-qJD(4) + t284) * t255 - t307;
t206 = t253 * qJD(4) + t301;
t294 = qJD(3) ^ 2;
t276 = -t294 - t336;
t275 = t294 - t336;
t274 = -t294 - t337;
t273 = -t294 + t337;
t268 = (-t286 + t287) * t295;
t267 = t321 * t295;
t263 = -0.2e1 * t309 + t314;
t260 = 0.2e1 * t308 + t316;
t258 = t321 * t318;
t244 = t284 * t253;
t242 = -t252 + t282;
t241 = t251 - t282;
t240 = -t289 * t262 - t287 * t318;
t239 = -t292 * t261 - t286 * t318;
t235 = -t252 - t282;
t234 = -t289 * t276 - t325;
t233 = t292 * t274 - t331;
t232 = t292 * t276 - t332;
t231 = -t292 * t275 - t331;
t230 = t289 * t274 + t324;
t229 = -t289 * t273 - t325;
t226 = -t293 * t267 - t330;
t225 = -t290 * t267 + t322;
t222 = t289 * t260 - t292 * t263;
t220 = t252 - t251;
t216 = -t282 - t251;
t215 = t290 * t232 + t293 * t263;
t214 = t290 * t230 + t293 * t260;
t212 = -t293 * t232 + t290 * t263;
t211 = -t293 * t230 + t290 * t260;
t209 = (t253 * t291 + t255 * t288) * t284;
t208 = (t253 * t288 - t255 * t291) * t284;
t207 = -t251 - t252;
t205 = -t255 * qJD(4) - t307;
t200 = t291 * t241 - t334;
t199 = -t288 * t242 + t346;
t198 = t288 * t241 + t328;
t197 = t291 * t242 + t347;
t194 = -t288 * t235 - t328;
t193 = t291 * t235 - t334;
t192 = t206 - t244;
t191 = t206 + t244;
t190 = t319 * t253 + t301;
t187 = t319 * t255 + t307;
t186 = t291 * t206 - t255 * t339;
t185 = t288 * t206 + t255 * t338;
t184 = -t288 * t205 - t253 * t338;
t183 = t291 * t205 - t253 * t339;
t182 = -pkin(2) * t267 - t196;
t181 = t291 * t216 - t347;
t180 = t288 * t216 + t346;
t179 = pkin(2) * t232 - qJ(2) * t234 - t224;
t178 = pkin(2) * t230 - qJ(2) * t233 + t223;
t177 = pkin(2) * t260 - t344 * t233 - t327;
t176 = pkin(2) * t263 - t344 * t234 + t333;
t175 = t290 * t195 - t293 * t238;
t174 = -t293 * t195 - t290 * t238;
t173 = -t292 * t208 - t289 * t209;
t170 = -pkin(6) * t193 - t329;
t168 = pkin(2) * t195 - qJ(2) * t196;
t167 = -t292 * t198 - t289 * t200;
t166 = -t292 * t197 - t289 * t199;
t165 = -pkin(6) * t180 - t335;
t164 = -t289 * t193 + t292 * t194;
t163 = t292 * t193 + t289 * t194;
t162 = t288 * t192 + t291 * t297;
t161 = -t291 * t187 - t288 * t191;
t160 = -t291 * t192 + t288 * t297;
t159 = -t288 * t187 + t291 * t191;
t158 = -pkin(2) * t238 - t344 * t196;
t157 = -t292 * t185 - t289 * t186;
t156 = -t292 * t183 - t289 * t184;
t155 = -t289 * t180 + t292 * t181;
t154 = t292 * t180 + t289 * t181;
t153 = -pkin(3) * t190 + pkin(6) * t194 - t335;
t152 = -pkin(3) * t187 + pkin(6) * t181 + t329;
t151 = t290 * t163 + t293 * t190;
t150 = -t293 * t163 + t290 * t190;
t149 = t290 * t154 + t293 * t187;
t148 = -t293 * t154 + t290 * t187;
t145 = -t289 * t160 + t292 * t162;
t144 = t292 * t160 + t289 * t162;
t143 = -t292 * t159 - t289 * t161;
t142 = pkin(3) * t204 + pkin(6) * t147;
t141 = t290 * t144 + t293 * t207;
t140 = -t293 * t144 + t290 * t207;
t139 = -pkin(6) * t160 - t146;
t138 = -pkin(3) * t207 + pkin(6) * t162 + t147;
t137 = pkin(2) * t163 + pkin(3) * t193 - qJ(2) * t164 - t172;
t136 = t147 * t292 - t342;
t135 = t147 * t289 + t341;
t134 = pkin(2) * t154 + pkin(3) * t180 - qJ(2) * t155 - t171;
t133 = t135 * t290 - t204 * t293;
t132 = -t135 * t293 - t204 * t290;
t131 = pkin(2) * t190 - t292 * t153 - t344 * t164 - t289 * t170;
t130 = pkin(2) * t144 + pkin(3) * t160 - qJ(2) * t145;
t129 = pkin(2) * t187 - t292 * t152 - t344 * t155 - t289 * t165;
t128 = pkin(2) * t135 + pkin(3) * t146 - qJ(2) * t136;
t127 = pkin(2) * t207 - t292 * t138 - t289 * t139 - t344 * t145;
t126 = -pkin(2) * t204 + pkin(6) * t342 - t344 * t136 - t292 * t142;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t266, -t265, 0, t228, 0, 0, 0, 0, 0, 0, 0, t266, t265, t213, 0, 0, 0, 0, 0, 0, t214, t215, t226, t175, 0, 0, 0, 0, 0, 0, t149, t151, t141, t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t265, -t266, 0, t227, 0, 0, 0, 0, 0, 0, 0, -t265, t266, t210, 0, 0, 0, 0, 0, 0, t211, t212, t225, t174, 0, 0, 0, 0, 0, 0, t148, t150, t140, t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t233, t234, 0, t196, 0, 0, 0, 0, 0, 0, t155, t164, t145, t136; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t265, 0, -t266, 0, -t304, -t303, -t227, -pkin(4) * t227, 0, -t265, t266, 0, 0, 0, -t210, t304, t303, -pkin(4) * t210 + (-pkin(1) * t290 + qJ(2) * t293) * g(3), -t240 * t290 + t305, -t222 * t290 + t268 * t293, -t231 * t290 + t292 * t313, -t239 * t290 - t305, -t229 * t290 - t289 * t313, qJDD(3) * t293 - t258 * t290, -pkin(4) * t211 - t177 * t290 + t178 * t293, -pkin(4) * t212 - t176 * t290 + t179 * t293, -pkin(2) * t322 - pkin(4) * t225 - t182 * t290, -pkin(4) * t174 - t158 * t290 + t168 * t293, -t157 * t290 - t311, -t143 * t290 + t220 * t293, -t166 * t290 + t192 * t293, -t156 * t290 + t311, -t167 * t290 + t293 * t297, -t173 * t290 + t283 * t293, -pkin(4) * t148 - t129 * t290 + t134 * t293, -pkin(4) * t150 - t131 * t290 + t137 * t293, -pkin(4) * t140 - t127 * t290 + t130 * t293, -pkin(4) * t132 - t126 * t290 + t128 * t293; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t266, 0, t265, 0, t303, -t304, t228, pkin(4) * t228, 0, -t266, -t265, 0, 0, 0, t213, -t303, t304, pkin(4) * t213 + (pkin(1) * t293 + qJ(2) * t290) * g(3), t240 * t293 + t306, t222 * t293 + t268 * t290, t231 * t293 + t290 * t314, t239 * t293 - t306, t229 * t293 - t289 * t315, qJDD(3) * t290 + t258 * t293, pkin(4) * t214 + t177 * t293 + t178 * t290, pkin(4) * t215 + t176 * t293 + t179 * t290, -pkin(2) * t330 + pkin(4) * t226 + t182 * t293, pkin(4) * t175 + t158 * t293 + t168 * t290, t157 * t293 - t312, t143 * t293 + t220 * t290, t166 * t293 + t192 * t290, t156 * t293 + t312, t167 * t293 + t290 * t297, t173 * t293 + t283 * t290, pkin(4) * t149 + t129 * t293 + t134 * t290, pkin(4) * t151 + t131 * t293 + t137 * t290, pkin(4) * t141 + t127 * t293 + t130 * t290, pkin(4) * t133 + t126 * t293 + t128 * t290; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t271, t272, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t302 - 0.2e1 * t343, -t272 + t280 + 0.2e1 * t285, pkin(1) * t246 + qJ(2) * t245, (t262 - t309) * t292, -t260 * t292 - t263 * t289, -t275 * t289 + t324, (-t261 + t308) * t289, t273 * t292 - t332, 0, qJ(2) * t260 - t344 * t230 - t333, qJ(2) * t263 - t344 * t232 - t327, -qJ(2) * t267 + t344 * t264 - t195, -qJ(2) * t238 - t344 * t195, -t185 * t289 + t186 * t292, -t159 * t289 + t161 * t292, -t197 * t289 + t199 * t292, -t183 * t289 + t184 * t292, -t198 * t289 + t200 * t292, -t208 * t289 + t209 * t292, qJ(2) * t187 - t289 * t152 - t344 * t154 + t292 * t165, qJ(2) * t190 - t289 * t153 - t344 * t163 + t292 * t170, qJ(2) * t207 - t289 * t138 + t292 * t139 - t344 * t144, -pkin(6) * t341 - qJ(2) * t204 - t344 * t135 - t289 * t142;];
tauB_reg = t1;
