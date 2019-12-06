% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PPPRR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PPPRR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:39
% EndTime: 2019-12-05 14:59:45
% DurationCPUTime: 3.43s
% Computational Cost: add. (10083->351), mult. (16324->565), div. (0->0), fcn. (13340->10), ass. (0->236)
t290 = cos(qJ(4));
t292 = qJD(4) ^ 2;
t288 = sin(qJ(4));
t304 = qJDD(4) * t288;
t259 = t290 * t292 + t304;
t303 = t290 * qJDD(4);
t260 = t288 * t292 - t303;
t281 = sin(pkin(8));
t283 = cos(pkin(9));
t284 = cos(pkin(8));
t314 = t283 * t284;
t221 = -t259 * t281 + t260 * t314;
t282 = sin(pkin(7));
t280 = sin(pkin(9));
t285 = cos(pkin(7));
t319 = t280 * t285;
t193 = t221 * t282 - t260 * t319;
t320 = t280 * t282;
t195 = t221 * t285 + t260 * t320;
t308 = g(3) - qJDD(1);
t334 = t282 * t308;
t333 = t285 * t308;
t332 = pkin(1) * t280;
t264 = g(1) * t285 + g(2) * t282;
t249 = -t284 * t264 - t281 * t308;
t263 = g(1) * t282 - g(2) * t285;
t257 = -qJDD(2) + t263;
t210 = t249 * t283 - t257 * t280;
t248 = -t264 * t281 + t284 * t308;
t243 = qJDD(3) + t248;
t176 = t210 * t288 - t290 * t243;
t166 = -qJDD(4) * pkin(4) - t292 * pkin(6) + t176;
t287 = sin(qJ(5));
t331 = t166 * t287;
t289 = cos(qJ(5));
t330 = t166 * t289;
t209 = t249 * t280 + t283 * t257;
t329 = t209 * t280;
t328 = t209 * t283;
t327 = t243 * t284;
t271 = t287 * t292 * t289;
t265 = qJDD(5) + t271;
t326 = t265 * t287;
t325 = t265 * t289;
t266 = qJDD(5) - t271;
t324 = t266 * t287;
t323 = t266 * t289;
t277 = t287 ^ 2;
t322 = t277 * t292;
t321 = t280 * t281;
t318 = t281 * t257;
t317 = t281 * t283;
t316 = t282 * t283;
t315 = t282 * t284;
t313 = t283 * t285;
t312 = t284 * t257;
t311 = t284 * t285;
t310 = t288 * t209;
t309 = t290 * t209;
t177 = t290 * t210 + t288 * t243;
t167 = -pkin(4) * t292 + qJDD(4) * pkin(6) + t177;
t146 = t289 * t167 + t287 * t209;
t278 = t289 ^ 2;
t307 = t277 + t278;
t306 = qJD(4) * qJD(5);
t305 = qJDD(4) * t287;
t274 = t289 * qJDD(4);
t300 = pkin(2) * t280 + pkin(5);
t299 = t287 * t306;
t298 = t289 * t306;
t145 = t167 * t287 - t289 * t209;
t118 = t145 * t287 + t146 * t289;
t258 = t307 * qJDD(4);
t276 = t278 * t292;
t261 = t276 + t322;
t226 = t258 * t288 + t261 * t290;
t227 = t258 * t290 - t261 * t288;
t297 = pkin(3) * t226 + pkin(4) * t261 + pkin(6) * t258 - qJ(3) * t227 + t118;
t296 = -pkin(3) * t260 + qJ(3) * t259 - t176;
t295 = -pkin(3) * t259 - qJ(3) * t260 - t177;
t199 = t248 * t281 + t284 * t249;
t229 = -t263 * t282 - t285 * t264;
t294 = t288 * t271;
t293 = t290 * t271;
t117 = t145 * t289 - t146 * t287;
t138 = t176 * t290 - t177 * t288;
t139 = t176 * t288 + t177 * t290;
t159 = -t210 * t280 + t328;
t160 = t210 * t283 + t329;
t198 = t248 * t284 - t249 * t281;
t228 = t263 * t285 - t264 * t282;
t219 = t259 * t284 + t260 * t317;
t291 = qJD(5) ^ 2;
t270 = -t276 - t291;
t269 = t276 - t291;
t268 = -t291 - t322;
t267 = t291 - t322;
t262 = t276 - t322;
t256 = t274 - 0.2e1 * t299;
t255 = t274 - t299;
t254 = t298 + t305;
t253 = 0.2e1 * t298 + t305;
t252 = t307 * t306;
t247 = qJDD(5) * t288 + t252 * t290;
t246 = t254 * t289 - t277 * t306;
t245 = qJDD(5) * t290 - t252 * t288;
t244 = -t255 * t287 - t278 * t306;
t242 = -t268 * t287 - t323;
t241 = -t267 * t287 + t325;
t240 = t270 * t289 - t326;
t239 = t269 * t289 - t324;
t238 = t268 * t289 - t324;
t237 = -t267 * t289 - t326;
t236 = t270 * t287 + t325;
t235 = -t269 * t287 - t323;
t234 = (-t254 - t298) * t287;
t233 = (-t255 + t299) * t289;
t225 = -t253 * t287 + t256 * t289;
t224 = -t253 * t289 - t256 * t287;
t223 = -t259 * t314 - t260 * t281;
t220 = -t259 * t317 + t260 * t284;
t218 = t246 * t290 - t294;
t217 = t244 * t290 + t294;
t216 = -t246 * t288 - t293;
t215 = -t244 * t288 + t293;
t214 = t241 * t290 + t287 * t304;
t213 = t239 * t290 + t288 * t274;
t212 = -t241 * t288 + t287 * t303;
t211 = -t239 * t288 + t289 * t303;
t208 = t242 * t290 + t253 * t288;
t207 = t240 * t290 - t256 * t288;
t206 = t242 * t288 - t253 * t290;
t205 = t240 * t288 + t256 * t290;
t201 = t225 * t290 - t262 * t288;
t200 = -t225 * t288 - t262 * t290;
t196 = t223 * t285 - t259 * t320;
t194 = t223 * t282 + t259 * t319;
t192 = -t245 * t281 + t247 * t314;
t191 = pkin(5) * t259 + t309;
t190 = pkin(5) * t260 + t310;
t189 = t226 * t281 + t227 * t314;
t188 = -t226 * t284 + t227 * t317;
t187 = t218 * t283 - t234 * t280;
t186 = t217 * t283 - t233 * t280;
t185 = -t218 * t280 - t234 * t283;
t184 = -t217 * t280 - t233 * t283;
t183 = t214 * t283 - t237 * t280;
t182 = t213 * t283 - t235 * t280;
t181 = -t214 * t280 - t237 * t283;
t180 = -t213 * t280 - t235 * t283;
t179 = t199 * t285 - t257 * t282;
t178 = t199 * t282 + t257 * t285;
t175 = t300 * t259 + t309;
t174 = -t300 * t260 - t310;
t173 = -t210 * t281 + t243 * t314;
t172 = -t209 * t281 + t280 * t327;
t171 = t208 * t283 + t238 * t280;
t170 = t207 * t283 + t236 * t280;
t169 = t208 * t280 - t238 * t283;
t168 = t207 * t280 - t236 * t283;
t162 = t201 * t283 - t224 * t280;
t161 = -t201 * t280 - t224 * t283;
t156 = t189 * t285 + t227 * t320;
t155 = t189 * t282 - t227 * t319;
t154 = -pkin(6) * t238 + t330;
t153 = -pkin(6) * t236 + t331;
t152 = t187 * t284 - t216 * t281;
t151 = t186 * t284 - t215 * t281;
t150 = t183 * t284 - t212 * t281;
t149 = t182 * t284 - t211 * t281;
t148 = t160 * t284 + t243 * t281;
t147 = t160 * t281 - t327;
t144 = t171 * t284 + t206 * t281;
t143 = t170 * t284 + t205 * t281;
t142 = t171 * t281 - t206 * t284;
t141 = t170 * t281 - t205 * t284;
t140 = t162 * t284 - t200 * t281;
t135 = -pkin(4) * t238 + t146;
t134 = -pkin(4) * t236 + t145;
t133 = t283 * t191 + t295 * t280;
t132 = t283 * t190 + t296 * t280;
t131 = -pkin(3) * t205 - pkin(4) * t256 - pkin(6) * t240 + t330;
t130 = -pkin(3) * t206 + pkin(4) * t253 - pkin(6) * t242 - t331;
t129 = t139 * t283 + t329;
t128 = t139 * t280 - t328;
t127 = -pkin(1) * t219 - pkin(2) * t259 - t280 * t191 + t295 * t283;
t126 = -pkin(1) * t220 - pkin(2) * t260 - t280 * t190 + t296 * t283;
t125 = t144 * t285 + t169 * t282;
t124 = t143 * t285 + t168 * t282;
t123 = t144 * t282 - t169 * t285;
t122 = t143 * t282 - t168 * t285;
t121 = t148 * t285 - t159 * t282;
t120 = t148 * t282 + t159 * t285;
t119 = -pkin(1) * t147 + pkin(2) * t243 - qJ(3) * t160;
t115 = -qJ(2) * t219 + t133 * t284 - t174 * t281;
t114 = -qJ(2) * t220 + t132 * t284 - t175 * t281;
t113 = -pkin(5) * t226 + t117 * t290;
t112 = t138 * t314 + t139 * t281;
t110 = -pkin(5) * t206 - t135 * t288 + t154 * t290;
t109 = -pkin(5) * t205 - t134 * t288 + t153 * t290;
t108 = -qJ(2) * t147 - (pkin(2) * t281 - qJ(3) * t284) * t159;
t107 = -t288 * t117 - t300 * t227;
t106 = t118 * t290 + t166 * t288;
t105 = t118 * t288 - t166 * t290;
t104 = t129 * t284 - t138 * t281;
t103 = t129 * t281 + t138 * t284;
t102 = -pkin(2) * t128 + pkin(3) * t209 - pkin(5) * t139;
t101 = -pkin(2) * t169 + pkin(3) * t238 - pkin(5) * t208 - t135 * t290 - t154 * t288;
t100 = -pkin(2) * t168 + pkin(3) * t236 - pkin(5) * t207 - t134 * t290 - t153 * t288;
t99 = -qJ(3) * t128 - (pkin(3) * t280 - pkin(5) * t283) * t138;
t98 = t104 * t285 + t128 * t282;
t97 = t104 * t282 - t128 * t285;
t96 = -qJ(3) * t169 + t110 * t283 - t130 * t280;
t95 = -qJ(3) * t168 + t109 * t283 - t131 * t280;
t94 = t106 * t283 - t117 * t280;
t93 = t106 * t280 + t117 * t283;
t92 = t283 * t113 + t297 * t280;
t91 = -pkin(3) * t105 + pkin(4) * t166 - pkin(6) * t118;
t90 = -pkin(1) * t142 + pkin(2) * t206 - qJ(3) * t171 - t110 * t280 - t130 * t283;
t89 = -pkin(1) * t141 + pkin(2) * t205 - qJ(3) * t170 - t109 * t280 - t131 * t283;
t88 = -pkin(1) * t188 + pkin(2) * t226 - t280 * t113 + t297 * t283;
t87 = -pkin(5) * t105 - (pkin(4) * t288 - pkin(6) * t290) * t117;
t86 = t105 * t281 + t284 * t94;
t85 = -t105 * t284 + t281 * t94;
t84 = -pkin(1) * t103 - qJ(3) * t129 - (pkin(3) * t283 + pkin(5) * t280 + pkin(2)) * t138;
t83 = -qJ(2) * t188 - t107 * t281 + t284 * t92;
t82 = -qJ(2) * t142 - t101 * t281 + t284 * t96;
t81 = -qJ(2) * t141 - t100 * t281 + t284 * t95;
t80 = -qJ(2) * t103 - t102 * t281 + t284 * t99;
t79 = t282 * t93 + t285 * t86;
t78 = t282 * t86 - t285 * t93;
t77 = -pkin(2) * t93 - pkin(5) * t106 - (pkin(4) * t290 + pkin(6) * t288 + pkin(3)) * t117;
t76 = -qJ(3) * t93 - t280 * t91 + t283 * t87;
t75 = -pkin(1) * t85 + pkin(2) * t105 - qJ(3) * t94 - t280 * t87 - t283 * t91;
t74 = -qJ(2) * t85 - t281 * t77 + t284 * t76;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, 0, 0, 0, 0, 0, 0, t196, t195, 0, t98, 0, 0, 0, 0, 0, 0, t124, t125, t156, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, 0, 0, 0, 0, 0, 0, t194, t193, 0, t97, 0, 0, 0, 0, 0, 0, t122, t123, t155, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t308, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, 0, 0, 0, 0, 0, 0, t220, t219, 0, t103, 0, 0, 0, 0, 0, 0, t141, t142, t188, t85; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t334, -t333, -t228, -qJ(1) * t228, 0, 0, 0, 0, 0, 0, -t248 * t282 - t285 * t318, -t249 * t282 - t257 * t311, t285 * t198, -qJ(1) * t178 - (pkin(1) * t282 - qJ(2) * t285) * t198, 0, 0, 0, 0, 0, 0, t172 * t285 - t243 * t316, t173 * t285 + t243 * t320, t159 * t311 + t160 * t282, -qJ(1) * t120 + t108 * t285 - t119 * t282, 0, 0, -t195, 0, t196, (t280 * t311 - t316) * qJDD(4), -qJ(1) * t194 + t114 * t285 - t126 * t282, -qJ(1) * t193 + t115 * t285 - t127 * t282, t112 * t285 + t138 * t320, -qJ(1) * t97 - t282 * t84 + t285 * t80, t152 * t285 - t185 * t282, t140 * t285 - t161 * t282, t150 * t285 - t181 * t282, t151 * t285 - t184 * t282, t149 * t285 - t180 * t282, t192 * t285 + t247 * t320, -qJ(1) * t122 - t282 * t89 + t285 * t81, -qJ(1) * t123 - t282 * t90 + t285 * t82, -qJ(1) * t155 - t282 * t88 + t285 * t83, -qJ(1) * t78 - t282 * t75 + t285 * t74; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t333, -t334, t229, qJ(1) * t229, 0, 0, 0, 0, 0, 0, t248 * t285 - t282 * t318, t249 * t285 - t282 * t312, t282 * t198, qJ(1) * t179 - (-pkin(1) * t285 - qJ(2) * t282) * t198, 0, 0, 0, 0, 0, 0, t172 * t282 + t243 * t313, t173 * t282 - t243 * t319, t159 * t315 - t160 * t285, qJ(1) * t121 + t108 * t282 + t119 * t285, 0, 0, -t193, 0, t194, (t280 * t315 + t313) * qJDD(4), qJ(1) * t196 + t114 * t282 + t126 * t285, qJ(1) * t195 + t115 * t282 + t127 * t285, t112 * t282 - t138 * t319, qJ(1) * t98 + t282 * t80 + t285 * t84, t152 * t282 + t185 * t285, t140 * t282 + t161 * t285, t150 * t282 + t181 * t285, t151 * t282 + t184 * t285, t149 * t282 + t180 * t285, t192 * t282 - t247 * t319, qJ(1) * t124 + t282 * t81 + t285 * t89, qJ(1) * t125 + t282 * t82 + t285 * t90, qJ(1) * t156 + t282 * t83 + t285 * t88, qJ(1) * t79 + t282 * t74 + t285 * t75; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t263, t264, 0, 0, 0, 0, 0, 0, 0, 0, t312, -t318, t199, pkin(1) * t257 + qJ(2) * t199, 0, 0, 0, 0, 0, 0, t209 * t284 + t243 * t321, t210 * t284 + t243 * t317, t281 * t159, qJ(2) * t148 - (-pkin(2) * t284 - qJ(3) * t281 - pkin(1)) * t159, 0, 0, -t219, 0, t220, qJDD(4) * t321, qJ(2) * t223 + t132 * t281 + t175 * t284 + t259 * t332, qJ(2) * t221 + t133 * t281 + t174 * t284 - t260 * t332, t138 * t317 - t139 * t284, -pkin(1) * t128 + qJ(2) * t104 + t102 * t284 + t281 * t99, t187 * t281 + t216 * t284, t162 * t281 + t200 * t284, t183 * t281 + t212 * t284, t186 * t281 + t215 * t284, t182 * t281 + t211 * t284, t245 * t284 + t247 * t317, -pkin(1) * t168 + qJ(2) * t143 + t100 * t284 + t281 * t95, -pkin(1) * t169 + qJ(2) * t144 + t101 * t284 + t281 * t96, qJ(2) * t189 + t107 * t284 - t227 * t332 + t281 * t92, -pkin(1) * t93 + qJ(2) * t86 + t281 * t76 + t284 * t77;];
tauB_reg = t1;
