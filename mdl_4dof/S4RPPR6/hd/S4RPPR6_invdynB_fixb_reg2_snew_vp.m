% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPPR6
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPPR6_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR6_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:46
% EndTime: 2019-12-31 16:40:49
% DurationCPUTime: 2.42s
% Computational Cost: add. (3928->287), mult. (9667->405), div. (0->0), fcn. (5953->6), ass. (0->201)
t291 = sin(qJ(4));
t289 = sin(pkin(6));
t290 = cos(pkin(6));
t293 = cos(qJ(4));
t305 = t289 * t291 + t290 * t293;
t252 = t305 * qJD(1);
t327 = qJD(1) * t290;
t328 = qJD(1) * t289;
t254 = -t291 * t327 + t293 * t328;
t343 = t254 * t252;
t354 = qJDD(4) - t343;
t356 = t291 * t354;
t355 = t293 * t354;
t287 = t289 ^ 2;
t288 = t290 ^ 2;
t296 = qJD(1) ^ 2;
t266 = (t287 + t288) * t296;
t257 = t289 * t266;
t292 = sin(qJ(1));
t294 = cos(qJ(1));
t320 = t294 * qJDD(1);
t228 = -t292 * t257 + t289 * t320;
t353 = pkin(4) * t228;
t321 = t292 * qJDD(1);
t230 = t294 * t257 + t289 * t321;
t352 = pkin(4) * t230;
t272 = t294 * g(1) + t292 * g(2);
t256 = -t296 * pkin(1) + qJDD(1) * qJ(2) - t272;
t324 = qJD(1) * qJD(2);
t351 = (-t256 - 0.2e1 * t324) * t290;
t301 = t305 * qJDD(1);
t248 = t252 ^ 2;
t249 = t254 ^ 2;
t350 = pkin(2) + pkin(3);
t349 = pkin(2) * t290;
t282 = t287 * qJDD(1);
t283 = t288 * qJDD(1);
t264 = t283 + t282;
t224 = t292 * t264 + t294 * t266;
t348 = pkin(4) * t224;
t258 = t290 * t266;
t313 = t290 * t320;
t229 = -t292 * t258 + t313;
t347 = pkin(4) * t229;
t346 = qJ(2) * t257;
t345 = qJ(3) * t289;
t344 = qJDD(1) * pkin(1);
t342 = t287 * t296;
t341 = t288 * t296;
t340 = t289 * t290;
t339 = t290 * t296;
t271 = t292 * g(1) - t294 * g(2);
t300 = t296 * qJ(2) - qJDD(2) + t271;
t299 = 0.2e1 * qJD(3) * t328 + t300;
t312 = pkin(1) + t345;
t198 = (t350 * t290 + t312) * qJDD(1) + t299 + (-t342 - t341) * pkin(5);
t338 = t291 * t198;
t212 = qJDD(4) + t343;
t337 = t291 * t212;
t246 = t300 + t344;
t336 = t292 * t246;
t335 = t292 * t296;
t334 = t293 * t198;
t333 = t293 * t212;
t332 = t294 * t246;
t331 = t290 * g(3) + t289 * t256;
t330 = pkin(1) * t266 + qJ(2) * t264;
t322 = qJDD(1) * t290;
t329 = pkin(1) * t322 - qJ(2) * t258;
t326 = t252 * qJD(4);
t325 = t254 * qJD(4);
t323 = qJDD(1) * t289;
t319 = t292 * t343;
t318 = t294 * t343;
t317 = pkin(1) + t349;
t316 = t289 * t324;
t314 = t289 * t322;
t276 = 0.2e1 * t316;
t309 = -t345 - t349;
t261 = t309 * qJD(1);
t311 = t261 * t328 + qJDD(3) + t331;
t304 = t276 + t311;
t183 = (-pkin(3) * t339 - pkin(5) * qJDD(1)) * t289 + t304;
t222 = -t289 * g(3) - t351;
t245 = t261 * t327;
t202 = t222 + t245;
t193 = -pkin(3) * t341 - pkin(5) * t322 + t202;
t158 = -t293 * t183 + t291 * t193;
t221 = t276 + t331;
t180 = t289 * t221 + t290 * t222;
t234 = -t292 * t271 - t294 * t272;
t270 = t320 - t335;
t310 = -pkin(4) * t270 - t292 * g(3);
t308 = pkin(2) * t289 - qJ(3) * t290;
t251 = -t291 * t322 + t293 * t323;
t159 = t291 * t183 + t293 * t193;
t137 = -t293 * t158 + t291 * t159;
t138 = t291 * t158 + t293 * t159;
t179 = t290 * t221 - t289 * t222;
t265 = t283 - t282;
t267 = (t287 - t288) * t296;
t307 = t294 * t265 + t292 * t267;
t306 = t292 * t265 - t294 * t267;
t233 = t294 * t271 - t292 * t272;
t269 = t294 * t296 + t321;
t303 = pkin(1) - t309;
t302 = t294 * t258 + t290 * t321;
t295 = qJD(4) ^ 2;
t259 = t308 * qJDD(1);
t247 = -pkin(4) * t269 + t294 * g(3);
t239 = -t249 - t295;
t238 = -t249 + t295;
t237 = t248 - t295;
t236 = t289 * t313 - t335 * t340;
t235 = t269 * t340;
t226 = pkin(4) * t302;
t225 = t294 * t264 - t292 * t266;
t223 = pkin(4) * t225;
t219 = t249 - t248;
t217 = t251 - t326;
t216 = t251 - 0.2e1 * t326;
t215 = -t301 - t325;
t214 = 0.2e1 * t325 + t301;
t210 = -t295 - t248;
t209 = qJDD(1) * t303 + t299;
t208 = (t312 + 0.2e1 * t349) * qJDD(1) + t299;
t207 = (t317 + 0.2e1 * t345) * qJDD(1) + t299;
t206 = (-t252 * t293 + t254 * t291) * qJD(4);
t205 = (t252 * t291 + t254 * t293) * qJD(4);
t204 = -pkin(2) * t282 + t290 * t207;
t203 = qJ(3) * t283 - t289 * t208;
t201 = -t248 - t249;
t200 = -t311 - 0.2e1 * t316;
t197 = t293 * t217 - t291 * t325;
t196 = -t291 * t217 - t293 * t325;
t195 = -t291 * t215 + t293 * t326;
t194 = -t293 * t215 - t291 * t326;
t192 = pkin(2) * t266 + t202;
t191 = -t291 * t239 - t333;
t190 = -t291 * t238 + t355;
t189 = t293 * t237 - t337;
t188 = t293 * t239 - t337;
t187 = -t293 * t238 - t356;
t186 = -t291 * t237 - t333;
t185 = qJ(3) * t266 + t304;
t184 = -pkin(2) * t342 - t245 + (qJ(3) * t339 + g(3)) * t289 + t351;
t182 = (-pkin(2) * t340 + qJ(3) * t288) * t296 + t304;
t177 = -t293 * t214 - t291 * t216;
t176 = t291 * t251 - t293 * t301;
t175 = t291 * t214 - t293 * t216;
t174 = -t293 * t251 - t291 * t301;
t173 = t293 * t210 - t356;
t172 = t291 * t210 + t355;
t171 = -t289 * t205 + t290 * t206;
t170 = t294 * t180 - t336;
t169 = t292 * t180 + t332;
t168 = -t289 * t200 + t290 * t202;
t167 = t290 * t200 + t289 * t202;
t166 = -t289 * t196 + t290 * t197;
t165 = -t289 * t194 + t290 * t195;
t164 = t289 * t188 + t290 * t191;
t163 = -t289 * t187 + t290 * t190;
t162 = -t289 * t186 + t290 * t189;
t161 = -t290 * t188 + t289 * t191;
t160 = t290 * t185 - t289 * t192;
t157 = t294 * t168 - t292 * t209;
t156 = t292 * t168 + t294 * t209;
t155 = -pkin(5) * t188 + qJ(3) * t216 + t334;
t154 = -t289 * t175 + t290 * t177;
t153 = t289 * t174 + t290 * t176;
t152 = -t290 * t174 + t289 * t176;
t151 = t289 * t172 + t290 * t173;
t150 = -t290 * t172 + t289 * t173;
t149 = t294 * t164 - t292 * t216;
t148 = t292 * t164 + t294 * t216;
t147 = -pkin(5) * t172 + qJ(3) * t214 + t338;
t146 = -qJ(2) * t167 - t308 * t209;
t145 = -pkin(5) * t191 + t350 * t216 - t338;
t144 = t294 * t151 - t292 * t214;
t143 = t292 * t151 + t294 * t214;
t142 = -pkin(5) * t173 + t350 * t214 + t334;
t141 = -pkin(1) * t167 - pkin(2) * t200 - qJ(3) * t202;
t140 = t294 * t153 - t292 * t201;
t139 = t292 * t153 + t294 * t201;
t136 = -pkin(5) * t137 + qJ(3) * t198;
t135 = -pkin(5) * t138 + t350 * t198;
t134 = -pkin(1) * t152 - qJ(3) * t176 + t350 * t174;
t133 = -pkin(5) * t174 + qJ(3) * t201 - t137;
t132 = -pkin(5) * t176 + t350 * t201 - t138;
t131 = -pkin(1) * t161 - qJ(3) * t191 + t350 * t188 - t159;
t130 = -qJ(2) * t161 - t289 * t145 + t290 * t155;
t129 = -pkin(1) * t150 - qJ(3) * t173 + t350 * t172 - t158;
t128 = t289 * t137 + t290 * t138;
t127 = -t290 * t137 + t289 * t138;
t126 = -qJ(2) * t150 - t289 * t142 + t290 * t147;
t125 = t294 * t128 - t292 * t198;
t124 = t292 * t128 + t294 * t198;
t123 = -qJ(2) * t152 - t289 * t132 + t290 * t133;
t122 = -qJ(2) * t127 - t289 * t135 + t290 * t136;
t121 = -pkin(1) * t127 - qJ(3) * t138 + t350 * t137;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t269, -t270, 0, t234, 0, 0, 0, 0, 0, 0, -t302, t230, t225, t170, 0, 0, 0, 0, 0, 0, -t302, t225, -t230, t157, 0, 0, 0, 0, 0, 0, t144, t149, t140, t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t270, -t269, 0, t233, 0, 0, 0, 0, 0, 0, t229, -t228, t224, t169, 0, 0, 0, 0, 0, 0, t229, t224, t228, t156, 0, 0, 0, 0, 0, 0, t143, t148, t139, t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, 0, 0, 0, 0, 0, 0, t150, t161, t152, t127; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t270, 0, -t269, 0, t310, -t247, -t233, -pkin(4) * t233, t236, t307, t230, -t236, t302, 0, -t292 * t221 - t289 * t332 - t347, -t292 * t222 - t290 * t332 + t353, t294 * t179 - t348, -pkin(4) * t169 - (pkin(1) * t292 - qJ(2) * t294) * t179, t236, t230, -t307, 0, -t302, -t236, -t292 * t182 + t294 * t203 - t347, t294 * t160 - t292 * t259 - t348, -t292 * t184 + t294 * t204 - t353, -pkin(4) * t156 - t292 * t141 + t294 * t146, t294 * t166 - t319, t294 * t154 - t292 * t219, t294 * t163 - t292 * t251, t294 * t165 + t319, t294 * t162 + t292 * t301, -t292 * qJDD(4) + t294 * t171, -pkin(4) * t143 + t294 * t126 - t292 * t129, -pkin(4) * t148 + t294 * t130 - t292 * t131, -pkin(4) * t139 + t294 * t123 - t292 * t134, -pkin(4) * t124 - t292 * t121 + t294 * t122; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t269, 0, t270, 0, t247, t310, t234, pkin(4) * t234, t235, t306, -t228, -t235, -t229, 0, t294 * t221 - t289 * t336 - t226, t294 * t222 - t290 * t336 + t352, t292 * t179 + t223, pkin(4) * t170 - (-pkin(1) * t294 - qJ(2) * t292) * t179, t235, -t228, -t306, 0, t229, -t235, t294 * t182 + t292 * t203 - t226, t292 * t160 + t294 * t259 + t223, t294 * t184 + t292 * t204 - t352, pkin(4) * t157 + t294 * t141 + t292 * t146, t292 * t166 + t318, t292 * t154 + t294 * t219, t292 * t163 + t294 * t251, t292 * t165 - t318, t292 * t162 - t294 * t301, t294 * qJDD(4) + t292 * t171, pkin(4) * t144 + t292 * t126 + t294 * t129, pkin(4) * t149 + t292 * t130 + t294 * t131, pkin(4) * t140 + t292 * t123 + t294 * t134, pkin(4) * t125 + t294 * t121 + t292 * t122; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t271, t272, 0, 0, t282, 0.2e1 * t314, 0, t283, 0, 0, t290 * t246 + t329, t346 + (-t246 - t344) * t289, t180 + t330, pkin(1) * t246 + qJ(2) * t180, t282, 0, -0.2e1 * t314, 0, 0, t283, (qJ(3) * t323 + t208) * t290 + t329, t289 * t185 + t290 * t192 + t330, -t346 + (t317 * qJDD(1) + t207) * t289, qJ(2) * t168 + t209 * t303, t290 * t196 + t289 * t197, t290 * t175 + t289 * t177, t290 * t187 + t289 * t190, t290 * t194 + t289 * t195, t290 * t186 + t289 * t189, t290 * t205 + t289 * t206, pkin(1) * t214 + qJ(2) * t151 + t290 * t142 + t289 * t147, pkin(1) * t216 + qJ(2) * t164 + t290 * t145 + t289 * t155, pkin(1) * t201 + qJ(2) * t153 + t290 * t132 + t289 * t133, pkin(1) * t198 + qJ(2) * t128 + t290 * t135 + t289 * t136;];
tauB_reg = t1;
