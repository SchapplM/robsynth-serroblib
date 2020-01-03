% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPRR9
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPRR9_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR9_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:30
% EndTime: 2019-12-31 16:56:34
% DurationCPUTime: 2.01s
% Computational Cost: add. (4579->288), mult. (8974->403), div. (0->0), fcn. (5332->6), ass. (0->199)
t287 = sin(qJ(4));
t290 = cos(qJ(4));
t291 = cos(qJ(3));
t323 = qJD(1) * t291;
t254 = -t290 * qJD(3) + t287 * t323;
t256 = t287 * qJD(3) + t290 * t323;
t223 = t256 * t254;
t321 = qJD(1) * qJD(3);
t309 = t291 * t321;
t288 = sin(qJ(3));
t319 = t288 * qJDD(1);
t259 = -t309 - t319;
t295 = qJDD(4) - t259;
t344 = -t223 + t295;
t347 = t287 * t344;
t346 = t290 * t344;
t294 = qJD(1) ^ 2;
t343 = pkin(5) + pkin(1);
t345 = t343 * t294;
t276 = t288 * qJD(1) + qJD(4);
t240 = t276 * t254;
t310 = t288 * t321;
t317 = t291 * qJDD(1);
t260 = -t310 + t317;
t312 = t254 * qJD(4) - t287 * qJDD(3) - t290 * t260;
t196 = t312 + t240;
t308 = -t290 * qJDD(3) + t287 * t260;
t192 = (qJD(4) - t276) * t256 + t308;
t250 = t254 ^ 2;
t251 = t256 ^ 2;
t275 = t276 ^ 2;
t342 = pkin(3) * t288;
t341 = qJDD(1) * pkin(1);
t340 = t276 * t287;
t339 = t276 * t290;
t285 = t288 ^ 2;
t338 = t285 * t294;
t286 = t291 ^ 2;
t337 = t286 * t294;
t289 = sin(qJ(1));
t292 = cos(qJ(1));
t269 = t289 * g(1) - t292 * g(2);
t302 = qJDD(2) - t269;
t296 = -t294 * qJ(2) + t302;
t239 = -t343 * qJDD(1) + t296;
t217 = t288 * g(3) + t291 * t239;
t293 = qJD(3) ^ 2;
t305 = -pkin(6) * t291 + t342;
t298 = t294 * t305;
t199 = qJDD(3) * pkin(3) + t293 * pkin(6) - t291 * t298 + t217;
t336 = t287 * t199;
t211 = t223 + t295;
t335 = t287 * t211;
t270 = t292 * g(1) + t289 * g(2);
t284 = qJDD(1) * qJ(2);
t299 = t270 - t284;
t320 = qJD(2) * qJD(1);
t234 = t299 - 0.2e1 * t320 + t345;
t334 = t288 * t234;
t314 = t288 * t294 * t291;
t267 = qJDD(3) + t314;
t333 = t288 * t267;
t268 = qJDD(3) - t314;
t332 = t288 * t268;
t324 = t285 + t286;
t262 = t324 * qJDD(1);
t331 = t289 * t262;
t330 = t290 * t199;
t329 = t290 * t211;
t328 = t291 * t234;
t327 = t291 * t267;
t326 = t291 * t268;
t325 = t292 * t262;
t282 = 0.2e1 * t320;
t297 = t282 - t299;
t300 = -t260 + t310;
t301 = -t259 + t309;
t187 = t301 * pkin(3) + t300 * pkin(6) + t297 - t345;
t218 = -t291 * g(3) + t288 * t239;
t200 = -t293 * pkin(3) + qJDD(3) * pkin(6) - t288 * t298 + t218;
t160 = t287 * t187 + t290 * t200;
t318 = t289 * qJDD(1);
t316 = t292 * qJDD(1);
t315 = t288 * t223;
t313 = t291 * t223;
t311 = pkin(3) * t291 + pkin(2);
t159 = -t290 * t187 + t287 * t200;
t138 = t287 * t159 + t290 * t160;
t241 = -t294 * pkin(1) + t297;
t242 = -t296 + t341;
t206 = t292 * t241 - t289 * t242;
t225 = -t289 * t269 - t292 * t270;
t307 = t289 * t314;
t306 = t292 * t314;
t263 = -t289 * t294 + t316;
t304 = pkin(4) * t263 + t289 * g(3);
t264 = t292 * t294 + t318;
t303 = -pkin(4) * t264 + t292 * g(3);
t137 = -t290 * t159 + t287 * t160;
t185 = t291 * t217 + t288 * t218;
t186 = -t288 * t217 + t291 * t218;
t203 = t289 * t241 + t292 * t242;
t224 = t292 * t269 - t289 * t270;
t274 = -t293 - t337;
t273 = t293 - t337;
t272 = -t293 - t338;
t271 = -t293 + t338;
t266 = (-t285 + t286) * t294;
t265 = t324 * t294;
t261 = -0.2e1 * t310 + t317;
t258 = 0.2e1 * t309 + t319;
t253 = t324 * t321;
t238 = -t251 + t275;
t237 = t250 - t275;
t236 = -t288 * t260 - t286 * t321;
t235 = -t291 * t259 - t285 * t321;
t231 = -t288 * t274 - t327;
t230 = t291 * t272 - t332;
t229 = t291 * t274 - t333;
t228 = -t291 * t273 - t332;
t227 = t288 * t272 + t326;
t226 = -t288 * t271 - t327;
t222 = -t251 + t250;
t221 = -t292 * t265 - t331;
t220 = -t289 * t265 + t325;
t219 = -t251 - t275;
t216 = t288 * t258 - t291 * t261;
t215 = -t275 - t250;
t213 = -t256 * qJD(4) - t308;
t209 = t250 + t251;
t208 = t289 * t229 + t292 * t261;
t207 = t289 * t227 + t292 * t258;
t205 = -t292 * t229 + t289 * t261;
t204 = -t292 * t227 + t289 * t258;
t202 = (-t254 * t290 + t256 * t287) * t276;
t201 = (-t254 * t287 - t256 * t290) * t276;
t197 = -t240 + t312;
t193 = (-qJD(4) - t276) * t256 - t308;
t191 = -t256 * t340 - t290 * t312;
t190 = t256 * t339 - t287 * t312;
t189 = -t287 * t213 + t254 * t339;
t188 = t290 * t213 + t254 * t340;
t182 = -t288 * t202 + t291 * t295;
t181 = t290 * t237 - t335;
t180 = -t287 * t238 + t346;
t179 = t287 * t237 + t329;
t178 = t290 * t238 + t347;
t177 = -t287 * t219 - t329;
t176 = t290 * t219 - t335;
t175 = -pkin(2) * t265 - t186;
t174 = pkin(2) * t229 - qJ(2) * t231 - t218;
t173 = pkin(2) * t227 - qJ(2) * t230 + t217;
t172 = t290 * t215 - t347;
t171 = t287 * t215 + t346;
t170 = pkin(2) * t258 - t343 * t230 - t328;
t169 = pkin(2) * t261 - t343 * t231 + t334;
t168 = t289 * t185 - t292 * t234;
t167 = -t292 * t185 - t289 * t234;
t166 = -t288 * t191 + t313;
t165 = -t288 * t189 - t313;
t164 = -t192 * t290 - t287 * t197;
t163 = t290 * t193 + t196 * t287;
t162 = -t192 * t287 + t290 * t197;
t161 = t287 * t193 - t196 * t290;
t157 = -t288 * t181 - t291 * t192;
t156 = -t288 * t180 - t291 * t197;
t155 = -pkin(6) * t176 - t330;
t154 = pkin(2) * t185 - qJ(2) * t186;
t153 = -pkin(6) * t171 - t336;
t152 = t291 * t177 - t288 * t196;
t151 = t288 * t177 + t291 * t196;
t150 = t291 * t172 - t288 * t193;
t149 = t288 * t172 + t291 * t193;
t148 = -t288 * t163 - t291 * t222;
t147 = -pkin(2) * t234 - t343 * t186;
t146 = t291 * t164 - t288 * t209;
t145 = t288 * t164 + t291 * t209;
t144 = -pkin(3) * t176 + t160;
t143 = -pkin(3) * t171 + t159;
t142 = t289 * t151 + t292 * t176;
t141 = -t292 * t151 + t289 * t176;
t140 = t289 * t149 + t292 * t171;
t139 = -t292 * t149 + t289 * t171;
t136 = t289 * t145 + t292 * t162;
t135 = -t292 * t145 + t289 * t162;
t134 = t291 * t138 - t288 * t199;
t133 = t288 * t138 + t291 * t199;
t132 = -pkin(6) * t162 - t137;
t131 = pkin(2) * t151 + pkin(3) * t196 + pkin(6) * t177 - qJ(2) * t152 - t336;
t130 = pkin(2) * t149 + pkin(3) * t193 + pkin(6) * t172 - qJ(2) * t150 + t330;
t129 = t289 * t133 + t292 * t137;
t128 = -t292 * t133 + t289 * t137;
t127 = pkin(2) * t176 - t291 * t144 - t343 * t152 - t288 * t155;
t126 = pkin(2) * t171 - t291 * t143 - t343 * t150 - t288 * t153;
t125 = pkin(2) * t145 + pkin(3) * t209 + pkin(6) * t164 - qJ(2) * t146 + t138;
t124 = pkin(2) * t133 + pkin(3) * t199 + pkin(6) * t138 - qJ(2) * t134;
t123 = -t288 * t132 - t343 * t146 + t311 * t162;
t122 = -t343 * t134 + (pkin(6) * t288 + t311) * t137;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t264, -t263, 0, t225, 0, 0, 0, 0, 0, 0, 0, t264, t263, t206, 0, 0, 0, 0, 0, 0, t207, t208, t221, t168, 0, 0, 0, 0, 0, 0, t140, t142, t136, t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t263, -t264, 0, t224, 0, 0, 0, 0, 0, 0, 0, -t263, t264, t203, 0, 0, 0, 0, 0, 0, t204, t205, t220, t167, 0, 0, 0, 0, 0, 0, t139, t141, t135, t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t230, t231, 0, t186, 0, 0, 0, 0, 0, 0, t150, t152, t146, t134; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t263, 0, -t264, 0, -t304, -t303, -t224, -pkin(4) * t224, 0, -t263, t264, 0, 0, 0, -t203, t304, t303, -pkin(4) * t203 + (-pkin(1) * t289 + qJ(2) * t292) * g(3), -t289 * t236 + t306, -t289 * t216 + t292 * t266, -t289 * t228 + t291 * t316, -t289 * t235 - t306, -t289 * t226 - t288 * t316, t292 * qJDD(3) - t289 * t253, -pkin(4) * t204 - t289 * t170 + t292 * t173, -pkin(4) * t205 - t289 * t169 + t292 * t174, -pkin(2) * t325 - pkin(4) * t220 - t289 * t175, -pkin(4) * t167 - t289 * t147 + t292 * t154, -t289 * t166 + t292 * t190, -t289 * t148 + t292 * t161, -t289 * t156 + t292 * t178, -t289 * t165 + t292 * t188, -t289 * t157 + t292 * t179, -t289 * t182 + t292 * t201, -pkin(4) * t139 - t289 * t126 + t292 * t130, -pkin(4) * t141 - t289 * t127 + t292 * t131, -pkin(4) * t135 - t289 * t123 + t292 * t125, -pkin(4) * t128 - t289 * t122 + t292 * t124; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t264, 0, t263, 0, t303, -t304, t225, pkin(4) * t225, 0, -t264, -t263, 0, 0, 0, t206, -t303, t304, pkin(4) * t206 + (pkin(1) * t292 + qJ(2) * t289) * g(3), t292 * t236 + t307, t292 * t216 + t289 * t266, t292 * t228 + t289 * t317, t292 * t235 - t307, t292 * t226 - t288 * t318, t289 * qJDD(3) + t292 * t253, pkin(4) * t207 + t292 * t170 + t289 * t173, pkin(4) * t208 + t292 * t169 + t289 * t174, -pkin(2) * t331 + pkin(4) * t221 + t292 * t175, pkin(4) * t168 + t292 * t147 + t289 * t154, t292 * t166 + t289 * t190, t292 * t148 + t289 * t161, t292 * t156 + t289 * t178, t292 * t165 + t289 * t188, t292 * t157 + t289 * t179, t292 * t182 + t289 * t201, pkin(4) * t140 + t292 * t126 + t289 * t130, pkin(4) * t142 + t292 * t127 + t289 * t131, pkin(4) * t136 + t292 * t123 + t289 * t125, pkin(4) * t129 + t292 * t122 + t289 * t124; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t269, t270, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t302 - 0.2e1 * t341, -t270 + t282 + 0.2e1 * t284, pkin(1) * t242 + qJ(2) * t241, -t300 * t291, -t291 * t258 - t288 * t261, -t288 * t273 + t326, t301 * t288, t291 * t271 - t333, 0, qJ(2) * t258 - t343 * t227 - t334, qJ(2) * t261 - t343 * t229 - t328, -qJ(2) * t265 + t343 * t262 - t185, -qJ(2) * t234 - t343 * t185, t291 * t191 + t315, t291 * t163 - t288 * t222, t291 * t180 - t288 * t197, t291 * t189 - t315, t291 * t181 - t288 * t192, t291 * t202 + t288 * t295, qJ(2) * t171 - t288 * t143 - t343 * t149 + t291 * t153, qJ(2) * t176 - t288 * t144 - t343 * t151 + t291 * t155, t291 * t132 + (qJ(2) + t342) * t162 - t343 * t145, -t343 * t133 + (qJ(2) + t305) * t137;];
tauB_reg = t1;
