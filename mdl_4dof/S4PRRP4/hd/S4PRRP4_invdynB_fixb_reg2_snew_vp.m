% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PRRP4
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PRRP4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP4_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:59
% EndTime: 2019-12-31 16:28:03
% DurationCPUTime: 2.60s
% Computational Cost: add. (3259->267), mult. (6981->349), div. (0->0), fcn. (4192->6), ass. (0->200)
t290 = qJD(3) ^ 2;
t286 = sin(qJ(3));
t279 = t286 ^ 2;
t291 = qJD(2) ^ 2;
t348 = t279 * t291;
t264 = t290 + t348;
t288 = cos(qJ(3));
t269 = t286 * t291 * t288;
t263 = qJDD(3) - t269;
t330 = t288 * t263;
t225 = -t286 * t264 + t330;
t322 = qJD(2) * qJD(3);
t312 = t288 * t322;
t320 = t286 * qJDD(2);
t253 = 0.2e1 * t312 + t320;
t287 = sin(qJ(2));
t289 = cos(qJ(2));
t179 = t287 * t225 + t289 * t253;
t183 = t289 * t225 - t287 * t253;
t283 = sin(pkin(6));
t284 = cos(pkin(6));
t135 = t284 * t179 + t283 * t183;
t379 = qJ(1) * t135;
t139 = t283 * t179 - t284 * t183;
t378 = qJ(1) * t139;
t377 = pkin(4) * t179;
t376 = pkin(1) * t179 + pkin(5) * t225;
t336 = t286 * t263;
t219 = t288 * t264 + t336;
t375 = -pkin(1) * t219 + pkin(4) * t183;
t314 = t286 * t322;
t318 = t288 * qJDD(2);
t254 = -0.2e1 * t314 + t318;
t215 = t288 * t254;
t216 = t286 * t253;
t200 = -t215 + t216;
t280 = t288 ^ 2;
t260 = (t279 - t280) * t291;
t171 = t287 * t200 + t289 * t260;
t173 = t289 * t200 - t287 * t260;
t374 = t284 * t171 + t283 * t173;
t373 = t283 * t171 - t284 * t173;
t347 = t280 * t291;
t266 = -t290 + t347;
t223 = -t288 * t266 + t336;
t317 = t289 * qJDD(2);
t188 = t287 * t223 + t288 * t317;
t191 = t289 * t223 - t287 * t318;
t372 = t284 * t188 + t283 * t191;
t371 = t283 * t188 - t284 * t191;
t261 = t284 * g(1) + t283 * g(2);
t309 = t283 * g(1) - t284 * g(2);
t306 = t287 * t261 + t289 * t309;
t327 = t289 * t261 - t287 * t309;
t308 = -t287 * t306 - t289 * t327;
t166 = t287 * t327 - t289 * t306;
t341 = t284 * t166;
t369 = -t283 * t308 + t341;
t346 = t283 * t166;
t128 = t284 * t308 + t346;
t319 = t287 * qJDD(2);
t257 = t289 * t291 + t319;
t258 = -t287 * t291 + t317;
t205 = -t283 * t257 + t284 * t258;
t281 = g(3) - qJDD(1);
t234 = pkin(4) * t257 - t289 * t281;
t300 = -pkin(4) * t258 - t287 * t281;
t368 = -qJ(1) * t205 + t283 * t234 + t284 * t300;
t367 = 2 * qJD(4);
t365 = pkin(2) * t219;
t364 = pkin(5) * t219;
t358 = t284 * t257 + t283 * t258;
t362 = qJ(1) * t358 + t284 * t234 - t283 * t300;
t354 = pkin(3) * t288;
t302 = -qJ(4) * t286 - t354;
t251 = t302 * qJD(2);
t198 = -t291 * pkin(2) + qJDD(2) * pkin(5) - t327;
t325 = -t288 * t198 + t286 * t281;
t299 = t288 * qJD(2) * t251 + qJDD(3) * qJ(4) + (qJD(3) * t367) - t325;
t357 = t286 * t266 + t330;
t323 = qJD(2) * t286;
t356 = t251 * t323 + qJDD(4);
t211 = -t284 * t261 - t283 * t309;
t210 = -t283 * t261 + t284 * t309;
t267 = -t290 - t347;
t262 = qJDD(3) + t269;
t337 = t286 * t262;
t222 = t288 * t267 - t337;
t178 = t287 * t222 + t289 * t254;
t353 = pkin(4) * t178;
t324 = t279 + t280;
t256 = t324 * qJDD(2);
t259 = t324 * t291;
t208 = t287 * t256 + t289 * t259;
t352 = pkin(4) * t208;
t247 = t288 * t262;
t217 = t286 * t267 + t247;
t351 = pkin(5) * t217;
t182 = t289 * t222 - t287 * t254;
t134 = t284 * t178 + t283 * t182;
t350 = qJ(1) * t134;
t209 = t289 * t256 - t287 * t259;
t160 = t284 * t208 + t283 * t209;
t349 = qJ(1) * t160;
t342 = t283 * t281;
t340 = t284 * t281;
t297 = qJDD(2) * pkin(2) + t291 * pkin(5) + t306;
t339 = t286 * t297;
t338 = t286 * t254;
t331 = t288 * t297;
t177 = t286 * t198 + t288 * t281;
t326 = t259 - t290;
t316 = pkin(1) * t178 + pkin(2) * t254 + pkin(5) * t222;
t315 = pkin(1) * t208 + pkin(2) * t259 + pkin(5) * t256;
t310 = -pkin(1) * t217 + pkin(4) * t182;
t142 = t286 * t177 - t288 * t325;
t305 = t287 * t269;
t304 = t289 * t269;
t155 = -pkin(2) * t217 + t177;
t301 = pkin(3) * t286 - qJ(4) * t288;
t141 = t288 * t177 + t286 * t325;
t298 = t288 * t253 + t338;
t296 = t312 + t320;
t295 = -t314 + t318;
t294 = -qJDD(3) * pkin(3) + t177 + t356;
t293 = -t295 * pkin(3) - t297 + (-t296 - t312) * qJ(4);
t292 = t323 * t367 - t293;
t265 = t290 - t348;
t250 = t301 * qJDD(2);
t246 = t324 * t322;
t230 = t287 * qJDD(3) + t289 * t246;
t229 = -t279 * t322 + t288 * t296;
t228 = -t289 * qJDD(3) + t287 * t246;
t227 = -t280 * t322 - t286 * t295;
t224 = -t286 * t265 + t247;
t218 = t288 * t265 + t337;
t202 = pkin(4) * t209;
t196 = t289 * t229 - t305;
t195 = t289 * t227 + t305;
t194 = t287 * t229 + t304;
t193 = t287 * t227 - t304;
t192 = t289 * t224 + t286 * t319;
t189 = t287 * t224 - t286 * t317;
t169 = -t283 * t228 + t284 * t230;
t168 = t284 * t228 + t283 * t230;
t163 = -t331 + t364;
t162 = -t339 - t351;
t161 = -t283 * t208 + t284 * t209;
t159 = qJ(1) * t161;
t158 = pkin(1) * t281 + pkin(4) * t308;
t157 = t290 * qJ(4) - t294;
t156 = -t325 + t365;
t154 = -t290 * pkin(3) + t299;
t153 = -t283 * t194 + t284 * t196;
t152 = -t283 * t193 + t284 * t195;
t151 = t284 * t194 + t283 * t196;
t150 = t284 * t193 + t283 * t195;
t149 = -t283 * t189 + t284 * t192;
t148 = t284 * t189 + t283 * t192;
t147 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t323 + t293;
t146 = t326 * qJ(4) + t294;
t145 = t326 * pkin(3) + t299;
t144 = (t254 - t314) * pkin(3) + t292;
t143 = -pkin(3) * t314 + qJ(4) * t253 + t292;
t137 = -t283 * t178 + t284 * t182;
t133 = qJ(1) * t137;
t132 = (-t267 - t290) * qJ(4) + (-qJDD(3) - t262) * pkin(3) + t155 + t356;
t131 = -t365 - qJ(4) * t263 + (-t264 + t290) * pkin(3) - t299;
t130 = t289 * t141 - t352;
t129 = t287 * t141 + t202;
t126 = t289 * t142 - t287 * t297;
t125 = t287 * t142 + t289 * t297;
t124 = -pkin(3) * t216 + t288 * t143 - t364;
t123 = qJ(4) * t215 - t286 * t144 - t351;
t122 = t288 * t154 - t286 * t157;
t121 = t286 * t154 + t288 * t157;
t120 = -t286 * t145 + t288 * t146;
t119 = -t287 * t156 + t289 * t163 + t377;
t118 = -t287 * t155 + t289 * t162 - t353;
t117 = t289 * t156 + t287 * t163 - t375;
t116 = t289 * t155 + t287 * t162 + t310;
t115 = t289 * t120 - t287 * t250 - t352;
t114 = t287 * t120 + t289 * t250 + t202;
t113 = t289 * t122 + t287 * t147;
t112 = t287 * t122 - t289 * t147;
t111 = -t283 * t125 + t284 * t126;
t110 = t284 * t125 + t283 * t126;
t109 = -pkin(2) * t121 - pkin(3) * t157 - qJ(4) * t154;
t108 = t289 * t123 - t287 * t132 - t353;
t107 = t289 * t124 - t287 * t131 - t377;
t106 = -pkin(5) * t121 + t301 * t147;
t105 = t287 * t123 + t289 * t132 + t310;
t104 = -pkin(4) * t125 - (pkin(2) * t287 - pkin(5) * t289) * t141;
t103 = t287 * t124 + t289 * t131 + t375;
t102 = pkin(4) * t126 - (-pkin(2) * t289 - pkin(5) * t287 - pkin(1)) * t141;
t101 = -t283 * t112 + t284 * t113;
t100 = t284 * t112 + t283 * t113;
t99 = -pkin(4) * t112 + t289 * t106 - t287 * t109;
t98 = -pkin(1) * t121 + pkin(4) * t113 + t287 * t106 + t289 * t109;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, 0, 0, 0, 0, 0, 0, -t358, -t205, 0, t128, 0, 0, 0, 0, 0, 0, t137, t139, t161, t111, 0, 0, 0, 0, 0, 0, t137, t161, -t139, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, 0, 0, 0, 0, 0, 0, t205, -t358, 0, -t369, 0, 0, 0, 0, 0, 0, t134, -t135, t160, t110, 0, 0, 0, 0, 0, 0, t134, t160, t135, t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t281, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t281, 0, 0, 0, 0, 0, 0, t217, -t219, 0, -t141, 0, 0, 0, 0, 0, 0, t217, 0, t219, t121; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t342, -t340, -t210, -qJ(1) * t210, 0, 0, t205, 0, -t358, 0, t368, t362, t369, pkin(4) * t341 + qJ(1) * t369 - t283 * t158, t153, t373, t149, t152, t371, t169, -t283 * t116 + t284 * t118 - t350, -t283 * t117 + t284 * t119 + t379, -t283 * t129 + t284 * t130 - t349, -qJ(1) * t110 - t283 * t102 + t284 * t104, t153, t149, -t373, t169, -t371, t152, -t283 * t105 + t284 * t108 - t350, -t283 * t114 + t284 * t115 - t349, -t283 * t103 + t284 * t107 - t379, -qJ(1) * t100 - t283 * t98 + t284 * t99; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t340, -t342, t211, qJ(1) * t211, 0, 0, t358, 0, t205, 0, -t362, t368, t128, pkin(4) * t346 + qJ(1) * t128 + t284 * t158, t151, -t374, t148, t150, -t372, t168, t284 * t116 + t283 * t118 + t133, t284 * t117 + t283 * t119 + t378, t284 * t129 + t283 * t130 + t159, qJ(1) * t111 + t284 * t102 + t283 * t104, t151, t148, t374, t168, t372, t150, t284 * t105 + t283 * t108 + t133, t284 * t114 + t283 * t115 + t159, t284 * t103 + t283 * t107 - t378, qJ(1) * t101 + t283 * t99 + t284 * t98; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t309, t261, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(1) * t258 + t306, -pkin(1) * t257 + t327, 0, -pkin(1) * t166, t216, t298, t218, t215, t357, 0, t316 + t331, -pkin(2) * t253 - t339 - t376, t142 + t315, pkin(1) * t125 + pkin(2) * t297 + pkin(5) * t142, t216, t218, -t298, 0, -t357, t215, qJ(4) * t338 + t288 * t144 + t316, t288 * t145 + t286 * t146 + t315, t286 * t143 + (pkin(2) + t354) * t253 + t376, pkin(1) * t112 + pkin(5) * t122 + (-pkin(2) + t302) * t147;];
tauB_reg = t1;
