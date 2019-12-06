% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PPPRR1
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PPPRR1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR1_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:09
% EndTime: 2019-12-05 14:58:16
% DurationCPUTime: 4.37s
% Computational Cost: add. (11253->339), mult. (18403->545), div. (0->0), fcn. (15186->10), ass. (0->233)
t288 = sin(pkin(9));
t290 = sin(pkin(7));
t293 = cos(pkin(7));
t271 = t293 * g(1) + t290 * g(2);
t289 = sin(pkin(8));
t292 = cos(pkin(8));
t318 = g(3) - qJDD(1);
t257 = -t292 * t271 - t289 * t318;
t270 = t290 * g(1) - t293 * g(2);
t264 = -qJDD(2) + t270;
t291 = cos(pkin(9));
t210 = t288 * t257 + t291 * t264;
t211 = t291 * t257 - t288 * t264;
t296 = sin(qJ(4));
t298 = cos(qJ(4));
t176 = t298 * t210 + t296 * t211;
t177 = -t296 * t210 + t298 * t211;
t305 = t296 * t176 + t298 * t177;
t135 = t298 * t176 - t296 * t177;
t328 = t291 * t135;
t106 = -t288 * t305 + t328;
t332 = t288 * t135;
t107 = t291 * t305 + t332;
t300 = qJD(4) ^ 2;
t313 = t296 * qJDD(4);
t266 = t298 * t300 + t313;
t312 = t298 * qJDD(4);
t267 = -t296 * t300 + t312;
t231 = t288 * t266 - t291 * t267;
t347 = t289 * t231;
t256 = -t289 * t271 + t292 * t318;
t251 = qJDD(3) + t256;
t214 = pkin(5) * t266 + t298 * t251;
t302 = -pkin(5) * t267 + t296 * t251;
t346 = t288 * t214 + t291 * t302;
t345 = t291 * t214 - t288 * t302;
t301 = -t291 * t266 - t288 * t267;
t326 = t292 * t293;
t192 = t231 * t326 + t290 * t301;
t329 = t290 * t292;
t190 = t231 * t329 - t293 * t301;
t341 = t290 * t318;
t337 = t293 * t318;
t336 = pkin(2) * t251;
t335 = t251 * t290;
t334 = t251 * t293;
t295 = sin(qJ(5));
t285 = t295 ^ 2;
t333 = t285 * t300;
t331 = t289 * t251;
t330 = t290 * t264;
t237 = t292 * t251;
t327 = t292 * t264;
t325 = t293 * t264;
t169 = -qJDD(4) * pkin(4) - t300 * pkin(6) + t176;
t324 = t295 * t169;
t297 = cos(qJ(5));
t278 = t295 * t300 * t297;
t272 = qJDD(5) + t278;
t323 = t295 * t272;
t273 = qJDD(5) - t278;
t322 = t295 * t273;
t321 = t297 * t169;
t320 = t297 * t272;
t319 = t297 * t273;
t170 = -t300 * pkin(4) + qJDD(4) * pkin(6) + t177;
t153 = t297 * t170 + t295 * t251;
t286 = t297 ^ 2;
t317 = t285 + t286;
t316 = qJD(4) * qJD(5);
t315 = qJDD(4) * t289;
t314 = t295 * qJDD(4);
t281 = t297 * qJDD(4);
t311 = t295 * t316;
t310 = t297 * t316;
t309 = -pkin(1) * t289 - qJ(3);
t152 = t295 * t170 - t297 * t251;
t124 = t295 * t152 + t297 * t153;
t265 = t317 * qJDD(4);
t283 = t286 * t300;
t268 = t283 + t333;
t233 = t296 * t265 + t298 * t268;
t234 = t298 * t265 - t296 * t268;
t194 = t291 * t233 + t288 * t234;
t195 = -t288 * t233 + t291 * t234;
t308 = -pkin(2) * t194 - pkin(3) * t233 - pkin(4) * t268 - pkin(6) * t265 + qJ(2) * t195 - t124;
t307 = -pkin(2) * t301 + pkin(3) * t266 + qJ(2) * t231 + t177;
t306 = pkin(2) * t231 - pkin(3) * t267 + qJ(2) * t301 + t176;
t200 = t289 * t256 + t292 * t257;
t236 = -t290 * t270 - t293 * t271;
t304 = t296 * t278;
t303 = t298 * t278;
t123 = t297 * t152 - t295 * t153;
t174 = t291 * t210 - t288 * t211;
t175 = t288 * t210 + t291 * t211;
t199 = t292 * t256 - t289 * t257;
t235 = t293 * t270 - t290 * t271;
t299 = qJD(5) ^ 2;
t277 = -t283 - t299;
t276 = t283 - t299;
t275 = -t299 - t333;
t274 = t299 - t333;
t269 = t283 - t333;
t263 = t281 - 0.2e1 * t311;
t262 = t281 - t311;
t261 = t310 + t314;
t260 = 0.2e1 * t310 + t314;
t259 = t317 * t316;
t255 = t296 * qJDD(5) + t298 * t259;
t254 = t297 * t261 - t285 * t316;
t253 = -t298 * qJDD(5) + t296 * t259;
t252 = -t295 * t262 - t286 * t316;
t250 = -t295 * t275 - t319;
t249 = -t295 * t274 + t320;
t248 = t297 * t277 - t323;
t247 = t297 * t276 - t322;
t246 = t297 * t275 - t322;
t245 = -t297 * t274 - t323;
t244 = t295 * t277 + t320;
t243 = -t295 * t276 - t319;
t242 = (-t261 - t310) * t295;
t241 = (-t262 + t311) * t297;
t226 = -t295 * t260 + t297 * t263;
t225 = -t297 * t260 - t295 * t263;
t224 = t289 * t301;
t223 = t298 * t254 - t304;
t222 = t298 * t252 + t304;
t221 = t296 * t254 + t303;
t220 = t296 * t252 - t303;
t219 = t298 * t249 + t295 * t313;
t218 = t298 * t247 + t281 * t296;
t217 = t296 * t249 - t295 * t312;
t216 = t296 * t247 - t297 * t312;
t209 = t298 * t250 + t296 * t260;
t208 = t298 * t248 - t296 * t263;
t207 = t296 * t250 - t298 * t260;
t206 = t296 * t248 + t298 * t263;
t202 = t298 * t226 - t296 * t269;
t201 = t296 * t226 + t298 * t269;
t197 = -t288 * t253 + t291 * t255;
t196 = -t291 * t253 - t288 * t255;
t193 = -t231 * t290 + t301 * t326;
t191 = t231 * t293 + t301 * t329;
t189 = t293 * t200 - t330;
t188 = t290 * t200 + t325;
t187 = -t289 * t211 + t237 * t291;
t186 = -t289 * t210 + t237 * t288;
t185 = -t288 * t221 + t291 * t223;
t184 = -t288 * t220 + t291 * t222;
t183 = -t291 * t221 - t288 * t223;
t182 = -t291 * t220 - t288 * t222;
t181 = -t288 * t217 + t291 * t219;
t180 = -t288 * t216 + t291 * t218;
t179 = -t291 * t217 - t288 * t219;
t178 = -t291 * t216 - t288 * t218;
t168 = -t288 * t207 + t291 * t209;
t167 = -t288 * t206 + t291 * t208;
t166 = t291 * t207 + t288 * t209;
t165 = t291 * t206 + t288 * t208;
t163 = -t288 * t201 + t291 * t202;
t162 = -t291 * t201 - t288 * t202;
t161 = t292 * t185 - t289 * t242;
t160 = t292 * t184 - t289 * t241;
t159 = t292 * t181 - t289 * t245;
t158 = t292 * t180 - t289 * t243;
t157 = -qJ(3) * t301 + t345;
t156 = qJ(3) * t231 + t346;
t155 = t292 * t175 + t331;
t154 = t289 * t175 - t237;
t151 = -pkin(6) * t246 + t321;
t150 = -pkin(6) * t244 + t324;
t149 = t292 * t168 + t289 * t246;
t148 = t292 * t167 + t289 * t244;
t147 = t289 * t168 - t292 * t246;
t146 = t289 * t167 - t292 * t244;
t143 = t290 * t194 + t195 * t326;
t142 = -t293 * t194 + t195 * t329;
t141 = t292 * t163 - t289 * t225;
t140 = t231 * t309 - t346;
t139 = t301 * t309 + t345;
t138 = -pkin(4) * t246 + t153;
t137 = -pkin(4) * t244 + t152;
t132 = -pkin(3) * t251 + pkin(5) * t305;
t131 = t293 * t155 - t174 * t290;
t130 = t290 * t155 + t174 * t293;
t129 = t293 * t149 + t290 * t166;
t128 = t293 * t148 + t290 * t165;
t127 = t290 * t149 - t293 * t166;
t126 = t290 * t148 - t293 * t165;
t125 = -pkin(1) * t154 - qJ(3) * t175 + t336;
t121 = -pkin(2) * t166 - pkin(3) * t207 + pkin(4) * t260 - pkin(6) * t250 - t324;
t120 = -pkin(2) * t165 - pkin(3) * t206 - pkin(4) * t263 - pkin(6) * t248 + t321;
t119 = t292 * t157 - t289 * t307;
t118 = t292 * t156 - t289 * t306;
t117 = -pkin(5) * t233 + t298 * t123;
t116 = pkin(5) * t234 + t296 * t123;
t115 = -qJ(2) * t154 - (pkin(2) * t289 - qJ(3) * t292) * t174;
t114 = -pkin(5) * t207 - t296 * t138 + t298 * t151;
t113 = -pkin(5) * t206 - t296 * t137 + t298 * t150;
t112 = -pkin(3) * t246 + pkin(5) * t209 + t298 * t138 + t296 * t151;
t111 = -pkin(3) * t244 + pkin(5) * t208 + t298 * t137 + t296 * t150;
t110 = t298 * t124 + t296 * t169;
t109 = t296 * t124 - t298 * t169;
t103 = t292 * t107 + t331;
t102 = t289 * t107 - t237;
t101 = pkin(2) * t106 + pkin(3) * t135;
t100 = -qJ(3) * t194 - t288 * t116 + t291 * t117;
t99 = -t291 * t116 - t288 * t117 + t195 * t309;
t98 = -t288 * t109 + t291 * t110;
t97 = t291 * t109 + t288 * t110;
t96 = pkin(5) * t328 + qJ(3) * t106 - t288 * t132;
t95 = -qJ(3) * t166 - t288 * t112 + t291 * t114;
t94 = -qJ(3) * t165 - t288 * t111 + t291 * t113;
t93 = t293 * t103 - t106 * t290;
t92 = t290 * t103 + t106 * t293;
t91 = -pkin(5) * t109 - (pkin(4) * t296 - pkin(6) * t298) * t123;
t90 = -pkin(1) * t147 + pkin(2) * t246 - qJ(3) * t168 - t291 * t112 - t288 * t114;
t89 = -pkin(1) * t146 + pkin(2) * t244 - qJ(3) * t167 - t291 * t111 - t288 * t113;
t88 = -t123 * t289 + t292 * t98;
t87 = t123 * t292 + t289 * t98;
t86 = t292 * t100 - t289 * t308;
t85 = pkin(5) * t110 - (-pkin(4) * t298 - pkin(6) * t296 - pkin(3)) * t123;
t84 = -qJ(2) * t147 - t289 * t121 + t292 * t95;
t83 = -qJ(2) * t146 - t289 * t120 + t292 * t94;
t82 = -pkin(1) * t102 - pkin(5) * t332 - qJ(3) * t107 - t291 * t132 + t336;
t81 = -pkin(2) * t97 - pkin(3) * t109 + pkin(4) * t169 - pkin(6) * t124;
t80 = t290 * t97 + t293 * t88;
t79 = t290 * t88 - t293 * t97;
t78 = -qJ(2) * t102 - t289 * t101 + t292 * t96;
t77 = -qJ(3) * t97 - t288 * t85 + t291 * t91;
t76 = -pkin(1) * t87 - pkin(2) * t123 - qJ(3) * t98 - t288 * t91 - t291 * t85;
t75 = -qJ(2) * t87 - t289 * t81 + t292 * t77;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, 0, 0, 0, 0, 0, 0, t193, t192, 0, t93, 0, 0, 0, 0, 0, 0, t128, t129, t143, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t235, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, 0, 0, 0, 0, 0, 0, t191, t190, 0, t92, 0, 0, 0, 0, 0, 0, t126, t127, t142, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t318, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, 0, 0, 0, 0, 0, 0, t224, t347, 0, t102, 0, 0, 0, 0, 0, 0, t146, t147, t289 * t195, t87; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t341, -t337, -t235, -qJ(1) * t235, 0, 0, 0, 0, 0, 0, -t290 * t256 - t289 * t325, -t290 * t257 - t292 * t325, t293 * t199, -qJ(1) * t188 - (pkin(1) * t290 - qJ(2) * t293) * t199, 0, 0, 0, 0, 0, 0, t293 * t186 - t291 * t335, t293 * t187 + t288 * t335, t174 * t326 + t175 * t290, -qJ(1) * t130 + t293 * t115 - t290 * t125, 0, 0, -t192, 0, t193, t293 * t315, -qJ(1) * t191 + t293 * t118 - t290 * t139, -qJ(1) * t190 + t293 * t119 - t290 * t140, t106 * t326 + t107 * t290, -qJ(1) * t92 - t290 * t82 + t293 * t78, t293 * t161 - t290 * t183, t293 * t141 - t290 * t162, t293 * t159 - t290 * t179, t293 * t160 - t290 * t182, t293 * t158 - t290 * t178, -t290 * t196 + t197 * t326, -qJ(1) * t126 - t290 * t89 + t293 * t83, -qJ(1) * t127 - t290 * t90 + t293 * t84, -qJ(1) * t142 - t290 * t99 + t293 * t86, -qJ(1) * t79 - t290 * t76 + t293 * t75; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t337, -t341, t236, qJ(1) * t236, 0, 0, 0, 0, 0, 0, t293 * t256 - t289 * t330, t293 * t257 - t290 * t327, t290 * t199, qJ(1) * t189 - (-pkin(1) * t293 - qJ(2) * t290) * t199, 0, 0, 0, 0, 0, 0, t290 * t186 + t291 * t334, t290 * t187 - t288 * t334, t174 * t329 - t175 * t293, qJ(1) * t131 + t290 * t115 + t293 * t125, 0, 0, -t190, 0, t191, t290 * t315, qJ(1) * t193 + t290 * t118 + t293 * t139, qJ(1) * t192 + t290 * t119 + t293 * t140, t106 * t329 - t107 * t293, qJ(1) * t93 + t290 * t78 + t293 * t82, t290 * t161 + t293 * t183, t290 * t141 + t293 * t162, t290 * t159 + t293 * t179, t290 * t160 + t293 * t182, t290 * t158 + t293 * t178, t293 * t196 + t197 * t329, qJ(1) * t128 + t290 * t83 + t293 * t89, qJ(1) * t129 + t290 * t84 + t293 * t90, qJ(1) * t143 + t290 * t86 + t293 * t99, qJ(1) * t80 + t290 * t75 + t293 * t76; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t270, t271, 0, 0, 0, 0, 0, 0, 0, 0, t327, -t289 * t264, t200, pkin(1) * t264 + qJ(2) * t200, 0, 0, 0, 0, 0, 0, t292 * t210 + t288 * t331, t292 * t211 + t291 * t331, t289 * t174, qJ(2) * t155 - (-pkin(2) * t292 - qJ(3) * t289 - pkin(1)) * t174, 0, 0, -t347, 0, t224, -t292 * qJDD(4), pkin(1) * t231 + t289 * t156 + t292 * t306, -pkin(1) * t301 + t289 * t157 + t292 * t307, t289 * t106, pkin(1) * t106 + qJ(2) * t103 + t292 * t101 + t289 * t96, t289 * t185 + t292 * t242, t289 * t163 + t292 * t225, t289 * t181 + t292 * t245, t289 * t184 + t292 * t241, t289 * t180 + t292 * t243, t289 * t197, -pkin(1) * t165 + qJ(2) * t148 + t292 * t120 + t289 * t94, -pkin(1) * t166 + qJ(2) * t149 + t292 * t121 + t289 * t95, -pkin(1) * t194 + t289 * t100 + t292 * t308, -pkin(1) * t97 + qJ(2) * t88 + t289 * t77 + t292 * t81;];
tauB_reg = t1;
