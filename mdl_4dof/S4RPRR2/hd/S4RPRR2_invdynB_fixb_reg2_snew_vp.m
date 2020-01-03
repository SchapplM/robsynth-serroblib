% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPRR2
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPRR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:14
% EndTime: 2019-12-31 16:48:18
% DurationCPUTime: 2.88s
% Computational Cost: add. (7609->275), mult. (12125->395), div. (0->0), fcn. (7328->8), ass. (0->188)
t279 = qJD(1) + qJD(3);
t277 = t279 ^ 2;
t289 = cos(qJ(3));
t278 = qJDD(1) + qJDD(3);
t286 = sin(qJ(3));
t310 = t278 * t286;
t243 = t277 * t289 + t310;
t309 = t278 * t289;
t246 = t277 * t286 - t309;
t283 = sin(pkin(7));
t284 = cos(pkin(7));
t196 = t243 * t284 - t246 * t283;
t282 = g(3) - qJDD(2);
t227 = pkin(5) * t243 - t282 * t289;
t337 = pkin(5) * t246 - t282 * t286;
t154 = qJ(2) * t196 + t227 * t284 - t283 * t337;
t200 = t243 * t283 + t246 * t284;
t287 = sin(qJ(1));
t290 = cos(qJ(1));
t164 = t196 * t287 + t200 * t290;
t343 = qJ(2) * t200 + t227 * t283 + t284 * t337;
t350 = pkin(4) * t164 + t154 * t287 + t290 * t343;
t336 = t196 * t290 - t200 * t287;
t349 = pkin(4) * t336 + t154 * t290 - t287 * t343;
t265 = g(1) * t287 - t290 * g(2);
t252 = qJDD(1) * pkin(1) + t265;
t266 = g(1) * t290 + g(2) * t287;
t292 = qJD(1) ^ 2;
t253 = -pkin(1) * t292 - t266;
t207 = t283 * t252 + t284 * t253;
t205 = -pkin(2) * t292 + t207;
t294 = t252 * t284 - t283 * t253;
t293 = qJDD(1) * pkin(2) + t294;
t170 = t205 * t286 - t289 * t293;
t171 = t289 * t205 + t286 * t293;
t302 = t170 * t286 + t289 * t171;
t130 = t170 * t289 - t171 * t286;
t322 = t130 * t284;
t110 = -t283 * t302 + t322;
t323 = t130 * t283;
t339 = t284 * t302 + t323;
t101 = t110 * t287 + t290 * t339;
t346 = t110 * t290 - t287 * t339;
t301 = t284 * t207 - t283 * t294;
t174 = -t207 * t283 - t284 * t294;
t318 = t174 * t290;
t340 = -t287 * t301 + t318;
t319 = t174 * t287;
t133 = t290 * t301 + t319;
t256 = qJDD(1) * t283 + t284 * t292;
t257 = qJDD(1) * t284 - t283 * t292;
t208 = -t256 * t287 + t290 * t257;
t231 = qJ(2) * t256 - t282 * t284;
t295 = -qJ(2) * t257 - t282 * t283;
t338 = -pkin(4) * t208 + t231 * t287 + t290 * t295;
t326 = t290 * t256 + t257 * t287;
t334 = pkin(4) * t326 + t231 * t290 - t287 * t295;
t324 = pkin(1) * t282;
t160 = -t278 * pkin(3) - t277 * pkin(6) + t170;
t285 = sin(qJ(4));
t321 = t160 * t285;
t288 = cos(qJ(4));
t320 = t160 * t288;
t264 = t288 * t277 * t285;
t254 = qJDD(4) + t264;
t317 = t254 * t285;
t255 = qJDD(4) - t264;
t316 = t255 * t285;
t315 = t255 * t288;
t280 = t285 ^ 2;
t312 = t277 * t280;
t311 = t278 * t285;
t268 = t288 * t278;
t161 = -pkin(3) * t277 + pkin(6) * t278 + t171;
t146 = t288 * t161 - t285 * t282;
t281 = t288 ^ 2;
t306 = t280 + t281;
t305 = qJD(4) * t279;
t304 = t285 * t305;
t303 = t288 * t305;
t145 = t161 * t285 + t288 * t282;
t124 = t145 * t285 + t288 * t146;
t223 = -t265 * t287 - t290 * t266;
t298 = t286 * t264;
t297 = t289 * t264;
t259 = qJDD(1) * t290 - t287 * t292;
t296 = -pkin(4) * t259 - g(3) * t287;
t123 = t145 * t288 - t146 * t285;
t222 = t265 * t290 - t266 * t287;
t291 = qJD(4) ^ 2;
t269 = t281 * t277;
t263 = -t269 - t291;
t262 = t269 - t291;
t261 = -t291 - t312;
t260 = t291 - t312;
t258 = qJDD(1) * t287 + t290 * t292;
t248 = t269 - t312;
t247 = t269 + t312;
t242 = t288 * t254;
t241 = t306 * t278;
t238 = t268 - 0.2e1 * t304;
t237 = t268 - t304;
t236 = t303 + t311;
t235 = 0.2e1 * t303 + t311;
t234 = -pkin(4) * t258 + g(3) * t290;
t233 = t306 * t305;
t221 = qJDD(4) * t286 + t233 * t289;
t220 = -qJDD(4) * t289 + t233 * t286;
t219 = t236 * t288 - t280 * t305;
t218 = -t237 * t285 - t281 * t305;
t217 = -t261 * t285 - t315;
t216 = -t260 * t285 + t242;
t215 = t263 * t288 - t317;
t214 = t262 * t288 - t316;
t213 = t261 * t288 - t316;
t212 = t263 * t285 + t242;
t203 = t241 * t289 - t247 * t286;
t202 = t241 * t286 + t247 * t289;
t192 = -t235 * t285 + t238 * t288;
t191 = t216 * t289 + t285 * t310;
t190 = t214 * t289 + t286 * t268;
t189 = t216 * t286 - t285 * t309;
t188 = t214 * t286 - t289 * t268;
t187 = t219 * t289 - t298;
t186 = t218 * t289 + t298;
t185 = t219 * t286 + t297;
t184 = t218 * t286 - t297;
t183 = t217 * t289 + t235 * t286;
t182 = t215 * t289 - t238 * t286;
t181 = t217 * t286 - t235 * t289;
t180 = t215 * t286 + t238 * t289;
t179 = -t220 * t283 + t221 * t284;
t178 = t220 * t284 + t221 * t283;
t177 = t192 * t289 - t248 * t286;
t176 = t192 * t286 + t248 * t289;
t169 = qJ(2) * t301 + t324;
t167 = -t202 * t283 + t203 * t284;
t166 = t202 * t284 + t203 * t283;
t158 = -t189 * t283 + t191 * t284;
t157 = -t188 * t283 + t190 * t284;
t156 = t189 * t284 + t191 * t283;
t155 = t188 * t284 + t190 * t283;
t150 = -t185 * t283 + t187 * t284;
t149 = -t184 * t283 + t186 * t284;
t148 = t185 * t284 + t187 * t283;
t147 = t184 * t284 + t186 * t283;
t143 = -t181 * t283 + t183 * t284;
t142 = -t180 * t283 + t182 * t284;
t141 = t181 * t284 + t183 * t283;
t140 = t180 * t284 + t182 * t283;
t139 = -pkin(6) * t213 + t320;
t138 = -pkin(6) * t212 + t321;
t137 = -pkin(3) * t213 + t146;
t136 = -pkin(3) * t212 + t145;
t135 = -t176 * t283 + t177 * t284;
t134 = t176 * t284 + t177 * t283;
t127 = -t166 * t287 + t167 * t290;
t126 = t166 * t290 + t167 * t287;
t125 = pkin(2) * t282 + pkin(5) * t302;
t121 = -t141 * t287 + t143 * t290;
t120 = -t140 * t287 + t142 * t290;
t119 = t141 * t290 + t143 * t287;
t118 = t140 * t290 + t142 * t287;
t117 = -pkin(5) * t202 + t123 * t289;
t116 = pkin(5) * t203 + t123 * t286;
t115 = t124 * t289 + t160 * t286;
t114 = t124 * t286 - t160 * t289;
t113 = -pkin(5) * t181 - t137 * t286 + t139 * t289;
t112 = -pkin(5) * t180 - t136 * t286 + t138 * t289;
t107 = -pkin(2) * t213 + pkin(5) * t183 + t137 * t289 + t139 * t286;
t106 = -pkin(2) * t212 + pkin(5) * t182 + t136 * t289 + t138 * t286;
t105 = -qJ(2) * t166 - t116 * t283 + t117 * t284;
t104 = qJ(2) * t167 + t116 * t284 + t117 * t283;
t103 = -t114 * t283 + t115 * t284;
t102 = t114 * t284 + t115 * t283;
t99 = pkin(5) * t322 + qJ(2) * t110 - t125 * t283;
t98 = pkin(5) * t323 + qJ(2) * t339 + t125 * t284 + t324;
t97 = -pkin(5) * t114 - (pkin(3) * t286 - pkin(6) * t289) * t123;
t96 = -qJ(2) * t141 - t107 * t283 + t113 * t284;
t95 = -qJ(2) * t140 - t106 * t283 + t112 * t284;
t94 = -pkin(1) * t213 + qJ(2) * t143 + t107 * t284 + t113 * t283;
t93 = -pkin(1) * t212 + qJ(2) * t142 + t106 * t284 + t112 * t283;
t92 = pkin(5) * t115 - (-pkin(3) * t289 - pkin(6) * t286 - pkin(2)) * t123;
t91 = -t102 * t287 + t103 * t290;
t90 = t102 * t290 + t103 * t287;
t89 = -qJ(2) * t102 - t283 * t92 + t284 * t97;
t88 = pkin(1) * t123 + qJ(2) * t103 + t283 * t97 + t284 * t92;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t258, -t259, 0, t223, 0, 0, 0, 0, 0, 0, -t326, -t208, 0, t133, 0, 0, 0, 0, 0, 0, -t336, t164, 0, t101, 0, 0, 0, 0, 0, 0, t120, t121, t127, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t259, -t258, 0, t222, 0, 0, 0, 0, 0, 0, t208, -t326, 0, -t340, 0, 0, 0, 0, 0, 0, -t164, -t336, 0, -t346, 0, 0, 0, 0, 0, 0, t118, t119, t126, t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t282, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t282, 0, 0, 0, 0, 0, 0, t212, t213, 0, -t123; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t259, 0, -t258, 0, t296, -t234, -t222, -pkin(4) * t222, 0, 0, t208, 0, -t326, 0, t338, t334, t340, pkin(4) * t340 + qJ(2) * t318 - t169 * t287, 0, 0, -t164, 0, -t336, 0, t350, t349, t346, pkin(4) * t346 - t287 * t98 + t290 * t99, -t148 * t287 + t150 * t290, -t134 * t287 + t135 * t290, -t156 * t287 + t158 * t290, -t147 * t287 + t149 * t290, -t155 * t287 + t157 * t290, -t178 * t287 + t179 * t290, -pkin(4) * t118 - t287 * t93 + t290 * t95, -pkin(4) * t119 - t287 * t94 + t290 * t96, -pkin(4) * t126 - t104 * t287 + t105 * t290, -pkin(4) * t90 - t287 * t88 + t290 * t89; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t258, 0, t259, 0, t234, t296, t223, pkin(4) * t223, 0, 0, t326, 0, t208, 0, -t334, t338, t133, pkin(4) * t133 + qJ(2) * t319 + t169 * t290, 0, 0, t336, 0, -t164, 0, -t349, t350, t101, pkin(4) * t101 + t287 * t99 + t290 * t98, t148 * t290 + t150 * t287, t134 * t290 + t135 * t287, t156 * t290 + t158 * t287, t147 * t290 + t149 * t287, t155 * t290 + t157 * t287, t178 * t290 + t179 * t287, pkin(4) * t120 + t287 * t95 + t290 * t93, pkin(4) * t121 + t287 * t96 + t290 * t94, pkin(4) * t127 + t104 * t290 + t105 * t287, pkin(4) * t91 + t287 * t89 + t290 * t88; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t265, t266, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t257 + t294, -pkin(1) * t256 - t207, 0, -pkin(1) * t174, 0, 0, 0, 0, 0, t278, -pkin(1) * t200 - pkin(2) * t246 - t170, -pkin(1) * t196 - pkin(2) * t243 - t171, 0, -pkin(1) * t110 - pkin(2) * t130, (t236 + t303) * t285, t235 * t288 + t238 * t285, t260 * t288 + t317, (t237 - t304) * t288, t262 * t285 + t315, 0, pkin(1) * t140 + pkin(2) * t180 + pkin(3) * t238 + pkin(6) * t215 - t320, pkin(1) * t141 + pkin(2) * t181 - pkin(3) * t235 + pkin(6) * t217 + t321, pkin(1) * t166 + pkin(2) * t202 + pkin(3) * t247 + pkin(6) * t241 + t124, pkin(1) * t102 + pkin(2) * t114 - pkin(3) * t160 + pkin(6) * t124;];
tauB_reg = t1;
