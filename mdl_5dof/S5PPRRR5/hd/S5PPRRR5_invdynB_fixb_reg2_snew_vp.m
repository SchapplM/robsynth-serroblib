% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PPRRR5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PPRRR5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR5_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:47
% EndTime: 2019-12-31 17:35:52
% DurationCPUTime: 3.40s
% Computational Cost: add. (7704->295), mult. (11245->414), div. (0->0), fcn. (8028->8), ass. (0->190)
t270 = qJD(3) + qJD(4);
t268 = t270 ^ 2;
t280 = cos(qJ(4));
t269 = qJDD(3) + qJDD(4);
t277 = sin(qJ(4));
t301 = t269 * t277;
t234 = t268 * t280 + t301;
t300 = t269 * t280;
t237 = t268 * t277 - t300;
t278 = sin(qJ(3));
t281 = cos(qJ(3));
t194 = t234 * t278 + t237 * t281;
t273 = g(3) - qJDD(1);
t222 = pkin(6) * t234 + t273 * t280;
t324 = pkin(6) * t237 + t273 * t277;
t145 = pkin(5) * t194 + t222 * t278 + t281 * t324;
t274 = sin(pkin(8));
t275 = cos(pkin(8));
t320 = -t234 * t281 + t237 * t278;
t330 = -pkin(5) * t320 + t222 * t281 - t278 * t324;
t331 = t194 * t275 + t274 * t320;
t336 = -qJ(1) * t331 + t145 * t275 - t274 * t330;
t154 = t194 * t274 - t275 * t320;
t335 = -qJ(1) * t154 + t145 * t274 + t275 * t330;
t251 = g(1) * t274 - t275 * g(2);
t248 = -qJDD(2) + t251;
t252 = g(1) * t275 + g(2) * t274;
t213 = -t278 * t248 - t281 * t252;
t283 = qJD(3) ^ 2;
t199 = -pkin(3) * t283 + t213;
t212 = t248 * t281 - t252 * t278;
t284 = qJDD(3) * pkin(3) - t212;
t165 = t199 * t277 - t280 * t284;
t166 = t280 * t199 + t277 * t284;
t291 = t165 * t277 + t280 * t166;
t125 = t280 * t165 - t277 * t166;
t312 = t125 * t281;
t105 = -t278 * t291 + t312;
t313 = t125 * t278;
t106 = t281 * t291 + t313;
t95 = t105 * t275 + t106 * t274;
t332 = t105 * t274 - t106 * t275;
t249 = qJDD(3) * t278 + t281 * t283;
t250 = -t281 * qJDD(3) + t278 * t283;
t203 = t275 * t249 + t250 * t274;
t225 = pkin(5) * t250 + t273 * t278;
t287 = pkin(5) * t249 + t273 * t281;
t326 = -qJ(1) * t203 + t225 * t274 + t275 * t287;
t290 = -t249 * t274 + t275 * t250;
t322 = -qJ(1) * t290 + t225 * t275 - t274 * t287;
t169 = t212 * t281 - t213 * t278;
t170 = t212 * t278 + t213 * t281;
t127 = t169 * t275 + t170 * t274;
t321 = t169 * t274 - t170 * t275;
t315 = pkin(1) + pkin(2);
t314 = qJ(2) * t273;
t158 = -t269 * pkin(4) - t268 * pkin(7) + t165;
t276 = sin(qJ(5));
t311 = t158 * t276;
t279 = cos(qJ(5));
t310 = t158 * t279;
t257 = t276 * t268 * t279;
t246 = qJDD(5) + t257;
t309 = t246 * t276;
t308 = t246 * t279;
t247 = qJDD(5) - t257;
t307 = t247 * t276;
t306 = t247 * t279;
t271 = t276 ^ 2;
t303 = t268 * t271;
t302 = t269 * t276;
t297 = t274 * t273;
t260 = t275 * t273;
t259 = t279 * t269;
t159 = -pkin(4) * t268 + pkin(7) * t269 + t166;
t152 = t279 * t159 + t276 * t273;
t272 = t279 ^ 2;
t296 = t271 + t272;
t295 = qJD(5) * t270;
t294 = t315 * t273;
t293 = t276 * t295;
t292 = t279 * t295;
t151 = t159 * t276 - t279 * t273;
t243 = t274 * t252;
t200 = t248 * t275 - t243;
t216 = t251 * t275 - t243;
t244 = t275 * t252;
t201 = -t248 * t274 - t244;
t217 = -t251 * t274 - t244;
t289 = t277 * t257;
t288 = t280 * t257;
t118 = t151 * t279 - t152 * t276;
t119 = t151 * t276 + t152 * t279;
t282 = qJD(5) ^ 2;
t261 = t272 * t268;
t256 = -t261 - t282;
t255 = t261 - t282;
t254 = -t282 - t303;
t253 = t282 - t303;
t240 = t261 - t303;
t239 = t261 + t303;
t233 = t296 * t269;
t232 = t259 - 0.2e1 * t293;
t231 = t259 - t293;
t230 = t292 + t302;
t229 = 0.2e1 * t292 + t302;
t228 = t296 * t295;
t219 = qJDD(5) * t277 + t228 * t280;
t218 = -qJDD(5) * t280 + t228 * t277;
t215 = t230 * t279 - t271 * t295;
t214 = -t231 * t276 - t272 * t295;
t211 = -t254 * t276 - t306;
t210 = -t253 * t276 + t308;
t209 = t256 * t279 - t309;
t208 = t255 * t279 - t307;
t205 = t254 * t279 - t307;
t204 = t256 * t276 + t308;
t192 = t233 * t280 - t239 * t277;
t188 = t233 * t277 + t239 * t280;
t187 = -t229 * t276 + t232 * t279;
t186 = t210 * t280 + t276 * t301;
t185 = t208 * t280 + t277 * t259;
t184 = t210 * t277 - t276 * t300;
t183 = t208 * t277 - t280 * t259;
t182 = t215 * t280 - t289;
t181 = t214 * t280 + t289;
t180 = t215 * t277 + t288;
t179 = t214 * t277 - t288;
t178 = t211 * t280 + t229 * t277;
t177 = t209 * t280 - t232 * t277;
t176 = t211 * t277 - t229 * t280;
t175 = t209 * t277 + t232 * t280;
t174 = -t218 * t278 + t219 * t281;
t173 = -t218 * t281 - t219 * t278;
t172 = t187 * t280 - t240 * t277;
t171 = t187 * t277 + t240 * t280;
t164 = pkin(5) * t169 + t314;
t162 = -pkin(5) * t170 + t294;
t161 = -t188 * t278 + t192 * t281;
t160 = t188 * t281 + t192 * t278;
t150 = -t184 * t278 + t186 * t281;
t149 = -t183 * t278 + t185 * t281;
t148 = -t184 * t281 - t186 * t278;
t147 = -t183 * t281 - t185 * t278;
t142 = -t180 * t278 + t182 * t281;
t141 = -t179 * t278 + t181 * t281;
t140 = -t180 * t281 - t182 * t278;
t139 = -t179 * t281 - t181 * t278;
t138 = -t176 * t278 + t178 * t281;
t137 = -t175 * t278 + t177 * t281;
t136 = t176 * t281 + t178 * t278;
t135 = t175 * t281 + t177 * t278;
t134 = -pkin(7) * t205 + t310;
t133 = -pkin(7) * t204 + t311;
t132 = -pkin(4) * t205 + t152;
t131 = -pkin(4) * t204 + t151;
t130 = -t171 * t278 + t172 * t281;
t129 = -t171 * t281 - t172 * t278;
t122 = -pkin(3) * t273 + pkin(6) * t291;
t121 = t160 * t274 + t161 * t275;
t120 = -t160 * t275 + t161 * t274;
t116 = t136 * t274 + t138 * t275;
t115 = t135 * t274 + t137 * t275;
t114 = -t136 * t275 + t138 * t274;
t113 = -t135 * t275 + t137 * t274;
t112 = -pkin(6) * t188 + t118 * t280;
t111 = pkin(6) * t192 + t118 * t277;
t110 = t119 * t280 + t158 * t277;
t109 = t119 * t277 - t158 * t280;
t108 = -pkin(6) * t176 - t132 * t277 + t134 * t280;
t107 = -pkin(6) * t175 - t131 * t277 + t133 * t280;
t102 = -pkin(3) * t205 + pkin(6) * t178 + t132 * t280 + t134 * t277;
t101 = -pkin(3) * t204 + pkin(6) * t177 + t131 * t280 + t133 * t277;
t100 = -pkin(5) * t160 - t111 * t278 + t112 * t281;
t99 = -pkin(5) * t161 - t111 * t281 - t112 * t278;
t98 = -t109 * t278 + t110 * t281;
t97 = t109 * t281 + t110 * t278;
t94 = pkin(5) * t105 + pkin(6) * t312 - t122 * t278 + t314;
t93 = -pkin(5) * t106 - pkin(6) * t313 - t122 * t281 + t294;
t92 = -pkin(6) * t109 - (pkin(4) * t277 - pkin(7) * t280) * t118;
t91 = -pkin(5) * t136 + qJ(2) * t205 - t102 * t278 + t108 * t281;
t90 = -pkin(5) * t135 + qJ(2) * t204 - t101 * t278 + t107 * t281;
t89 = -pkin(5) * t138 - t102 * t281 - t108 * t278 + t315 * t205;
t88 = -pkin(5) * t137 - t101 * t281 - t107 * t278 + t315 * t204;
t87 = pkin(6) * t110 - (-pkin(4) * t280 - pkin(7) * t277 - pkin(3)) * t118;
t86 = t274 * t97 + t275 * t98;
t85 = t274 * t98 - t275 * t97;
t84 = -pkin(5) * t97 - qJ(2) * t118 - t278 * t87 + t281 * t92;
t83 = -pkin(5) * t98 - t118 * t315 - t278 * t92 - t281 * t87;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t217, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, 0, 0, 0, 0, 0, 0, -t203, t290, 0, -t321, 0, 0, 0, 0, 0, 0, -t154, t331, 0, -t332, 0, 0, 0, 0, 0, 0, t115, t116, t121, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t216, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200, 0, 0, 0, 0, 0, 0, t290, t203, 0, t127, 0, 0, 0, 0, 0, 0, t331, t154, 0, t95, 0, 0, 0, 0, 0, 0, t113, t114, t120, t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t273, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t273, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t273, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t273, 0, 0, 0, 0, 0, 0, -t204, -t205, 0, t118; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t297, -t260, -t216, -qJ(1) * t216, 0, 0, 0, 0, 0, 0, -t297, -t200, t260, -qJ(1) * t200 + (-pkin(1) * t274 + qJ(2) * t275) * t273, 0, 0, -t290, 0, -t203, 0, t322, t326, t127, -qJ(1) * t127 - t162 * t274 + t164 * t275, 0, 0, -t331, 0, -t154, 0, t336, t335, t95, -qJ(1) * t95 - t274 * t93 + t275 * t94, -t140 * t274 + t142 * t275, -t129 * t274 + t130 * t275, -t148 * t274 + t150 * t275, -t139 * t274 + t141 * t275, -t147 * t274 + t149 * t275, -t173 * t274 + t174 * t275, -qJ(1) * t113 - t274 * t88 + t275 * t90, -qJ(1) * t114 - t274 * t89 + t275 * t91, -qJ(1) * t120 + t100 * t275 - t274 * t99, -qJ(1) * t85 - t274 * t83 + t275 * t84; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t260, -t297, t217, qJ(1) * t217, 0, 0, 0, 0, 0, 0, t260, t201, t297, qJ(1) * t201 + (pkin(1) * t275 + qJ(2) * t274) * t273, 0, 0, -t203, 0, t290, 0, t326, -t322, t321, -qJ(1) * t321 + t162 * t275 + t164 * t274, 0, 0, -t154, 0, t331, 0, t335, -t336, t332, -qJ(1) * t332 + t274 * t94 + t275 * t93, t140 * t275 + t142 * t274, t129 * t275 + t130 * t274, t148 * t275 + t150 * t274, t139 * t275 + t141 * t274, t147 * t275 + t149 * t274, t173 * t275 + t174 * t274, qJ(1) * t115 + t274 * t90 + t275 * t88, qJ(1) * t116 + t274 * t91 + t275 * t89, qJ(1) * t121 + t100 * t274 + t275 * t99, qJ(1) * t86 + t274 * t84 + t275 * t83; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t251, t252, 0, 0, 0, 0, 0, 0, 0, 0, t248, 0, -t252, pkin(1) * t248 - qJ(2) * t252, 0, 0, 0, 0, 0, -qJDD(3), -qJ(2) * t249 + t315 * t250 + t212, qJ(2) * t250 + t315 * t249 + t213, 0, qJ(2) * t170 + t169 * t315, 0, 0, 0, 0, 0, -t269, pkin(3) * t237 + qJ(2) * t320 + t194 * t315 + t165, pkin(3) * t234 + qJ(2) * t194 - t315 * t320 + t166, 0, pkin(3) * t125 + qJ(2) * t106 + t105 * t315, (-t230 - t292) * t276, -t229 * t279 - t232 * t276, -t253 * t279 - t309, (-t231 + t293) * t279, -t255 * t276 - t306, 0, -pkin(3) * t175 - pkin(4) * t232 - pkin(7) * t209 + qJ(2) * t137 - t315 * t135 + t310, -pkin(3) * t176 + pkin(4) * t229 - pkin(7) * t211 + qJ(2) * t138 - t315 * t136 - t311, -pkin(3) * t188 - pkin(4) * t239 - pkin(7) * t233 + qJ(2) * t161 - t315 * t160 - t119, -pkin(3) * t109 + pkin(4) * t158 - pkin(7) * t119 + qJ(2) * t98 - t315 * t97;];
tauB_reg = t1;
