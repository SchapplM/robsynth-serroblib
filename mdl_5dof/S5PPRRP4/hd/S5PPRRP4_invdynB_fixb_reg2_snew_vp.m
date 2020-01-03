% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PPRRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PPRRP4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP4_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:43
% EndTime: 2019-12-31 17:34:47
% DurationCPUTime: 2.19s
% Computational Cost: add. (4102->289), mult. (8253->377), div. (0->0), fcn. (4924->6), ass. (0->205)
t283 = cos(qJ(3));
t285 = qJD(3) ^ 2;
t281 = sin(qJ(3));
t308 = t281 * qJDD(3);
t250 = t283 * t285 + t308;
t306 = t283 * qJDD(3);
t251 = t281 * t285 - t306;
t277 = sin(pkin(7));
t278 = cos(pkin(7));
t205 = t278 * t250 + t277 * t251;
t275 = g(3) - qJDD(1);
t230 = pkin(5) * t251 + t281 * t275;
t293 = pkin(5) * t250 + t283 * t275;
t342 = -qJ(1) * t205 + t277 * t230 + t278 * t293;
t334 = pkin(1) + pkin(2);
t298 = -t277 * t250 + t278 * t251;
t340 = -qJ(1) * t298 + t278 * t230 - t277 * t293;
t254 = t277 * g(1) - t278 * g(2);
t248 = -qJDD(2) + t254;
t255 = t278 * g(1) + t277 * g(2);
t208 = t283 * t248 - t281 * t255;
t209 = -t281 * t248 - t283 * t255;
t168 = t283 * t208 - t281 * t209;
t169 = t281 * t208 + t283 * t209;
t127 = t278 * t168 + t277 * t169;
t339 = t277 * t168 - t278 * t169;
t280 = sin(qJ(4));
t282 = cos(qJ(4));
t263 = t280 * t285 * t282;
t305 = qJDD(4) + t263;
t338 = t305 * pkin(4);
t310 = qJD(3) * qJD(4);
t300 = t282 * t310;
t309 = t280 * qJDD(3);
t245 = t300 + t309;
t237 = t245 * qJ(5);
t197 = -t285 * pkin(3) + qJDD(3) * pkin(6) + t209;
t313 = -t280 * t197 + t282 * t275;
t336 = -t237 + t313 + t338;
t284 = qJD(4) ^ 2;
t273 = t280 ^ 2;
t325 = t273 * t285;
t260 = -t284 - t325;
t258 = qJDD(4) - t263;
t318 = t280 * t258;
t220 = t282 * t260 - t318;
t333 = pkin(3) * t220;
t274 = t282 ^ 2;
t312 = t273 + t274;
t249 = t312 * qJDD(3);
t272 = t274 * t285;
t252 = t272 + t325;
t210 = t281 * t249 + t283 * t252;
t332 = pkin(5) * t210;
t211 = t283 * t249 - t281 * t252;
t331 = pkin(5) * t211;
t262 = -t272 - t284;
t315 = t282 * t305;
t218 = t280 * t262 + t315;
t330 = pkin(6) * t218;
t329 = pkin(6) * t220;
t319 = t280 * t305;
t222 = t282 * t262 - t319;
t302 = t280 * t310;
t307 = t282 * qJDD(3);
t247 = -0.2e1 * t302 + t307;
t182 = t281 * t222 + t283 * t247;
t184 = t283 * t222 - t281 * t247;
t141 = -t278 * t182 + t277 * t184;
t328 = qJ(1) * t141;
t314 = t282 * t258;
t224 = -t280 * t260 - t314;
t244 = 0.2e1 * t300 + t309;
t183 = t281 * t224 - t283 * t244;
t185 = t283 * t224 + t281 * t244;
t142 = -t278 * t183 + t277 * t185;
t327 = qJ(1) * t142;
t170 = -t278 * t210 + t277 * t211;
t326 = qJ(1) * t170;
t322 = t277 * t275;
t267 = t278 * t275;
t147 = (qJ(5) * qJD(4) * t282 - 0.2e1 * qJD(5) * t280) * qJD(3) + t336;
t321 = t280 * t147;
t196 = -qJDD(3) * pkin(3) - t285 * pkin(6) + t208;
t320 = t280 * t196;
t317 = t282 * t147;
t316 = t282 * t196;
t179 = t282 * t197 + t280 * t275;
t311 = qJD(3) * t280;
t304 = 0.2e1 * qJD(3) * qJD(5);
t301 = t280 * t308;
t299 = t280 * t306;
t239 = t277 * t255;
t202 = t278 * t248 - t239;
t212 = t278 * t254 - t239;
t240 = t278 * t255;
t203 = -t277 * t248 - t240;
t213 = -t277 * t254 - t240;
t297 = t281 * t263;
t296 = t283 * t263;
t160 = -pkin(3) * t218 - t313;
t295 = -pkin(5) * t182 + qJ(2) * t218;
t294 = -pkin(5) * t183 + qJ(2) * t220;
t137 = -t280 * t179 - t282 * t313;
t138 = t282 * t179 - t280 * t313;
t291 = -pkin(5) * t184 + t334 * t218;
t290 = -pkin(5) * t185 + t334 * t220;
t246 = -t302 + t307;
t256 = qJD(4) * pkin(4) - qJ(5) * t311;
t289 = t246 * qJ(5) - qJD(4) * t256 + t282 * t304 + t179;
t288 = -pkin(3) * t247 - pkin(6) * t222 + qJ(2) * t184 - t334 * t182;
t287 = pkin(3) * t244 - pkin(6) * t224 + qJ(2) * t185 - t334 * t183;
t286 = -pkin(3) * t252 - pkin(6) * t249 + qJ(2) * t211 - t334 * t210;
t162 = -t246 * pkin(4) - qJ(5) * t272 + t256 * t311 + qJDD(5) + t196;
t265 = t280 * t304;
t261 = t272 - t284;
t259 = t284 - t325;
t253 = t272 - t325;
t243 = t312 * t310;
t228 = t281 * qJDD(4) + t283 * t243;
t227 = t282 * t245 - t273 * t310;
t226 = t283 * qJDD(4) - t281 * t243;
t225 = -t280 * t246 - t274 * t310;
t223 = -t280 * t259 + t315;
t221 = t282 * t261 - t318;
t219 = -t282 * t259 - t319;
t217 = -t280 * t261 - t314;
t216 = (-t245 - t300) * t280;
t215 = (-t246 + t302) * t282;
t214 = -pkin(4) * t244 - qJ(5) * t258;
t199 = -t280 * t244 + t282 * t247;
t198 = -t282 * t244 - t280 * t247;
t195 = t283 * t227 - t297;
t194 = t283 * t225 + t297;
t193 = -t281 * t227 - t296;
t192 = -t281 * t225 + t296;
t191 = t283 * t223 + t301;
t190 = t283 * t221 + t281 * t307;
t189 = -t281 * t223 + t299;
t188 = -t281 * t221 + t282 * t306;
t175 = t283 * t199 - t281 * t253;
t174 = -t281 * t199 - t283 * t253;
t173 = -t277 * t226 + t278 * t228;
t172 = t278 * t226 + t277 * t228;
t171 = t277 * t210 + t278 * t211;
t165 = t316 - t329;
t164 = t320 - t330;
t163 = qJ(1) * t171;
t161 = t179 - t333;
t159 = pkin(5) * t168 + qJ(2) * t275;
t158 = -pkin(5) * t169 + t334 * t275;
t157 = -qJ(5) * t260 + t162;
t156 = -t277 * t193 + t278 * t195;
t155 = -t277 * t192 + t278 * t194;
t154 = t278 * t193 + t277 * t195;
t153 = t278 * t192 + t277 * t194;
t152 = -t277 * t189 + t278 * t191;
t151 = -t277 * t188 + t278 * t190;
t150 = t278 * t189 + t277 * t191;
t149 = t278 * t188 + t277 * t190;
t148 = -pkin(4) * t272 + t289;
t146 = t265 + (-t300 + t309) * qJ(5) - t336;
t145 = pkin(4) * t247 + qJ(5) * t262 - t162;
t144 = t277 * t183 + t278 * t185;
t143 = t277 * t182 + t278 * t184;
t140 = qJ(1) * t144;
t139 = qJ(1) * t143;
t135 = qJ(5) * t307 + (t252 - t272) * pkin(4) + t289;
t134 = -t277 * t174 + t278 * t175;
t133 = t278 * t174 + t277 * t175;
t132 = -t333 + (-t260 - t272) * pkin(4) + t289;
t131 = -qJ(5) * t300 + t160 + t237 + t265 - 0.2e1 * t338;
t130 = t283 * t137 - t332;
t129 = -t281 * t137 - t331;
t126 = -qJ(5) * t315 - t280 * t145 - t330;
t125 = t282 * t157 - t280 * t214 - t329;
t124 = t283 * t138 + t281 * t196;
t123 = t281 * t138 - t283 * t196;
t122 = -pkin(4) * t162 + qJ(5) * t148;
t121 = t282 * t148 - t321;
t120 = t280 * t148 + t317;
t119 = -t280 * t135 + t282 * t146;
t118 = -t281 * t161 + t283 * t165 + t294;
t117 = -t281 * t160 + t283 * t164 + t295;
t116 = -pkin(4) * t301 + t283 * t119 - t332;
t115 = -pkin(4) * t299 - t281 * t119 - t331;
t114 = -t283 * t161 - t281 * t165 + t290;
t113 = -t283 * t160 - t281 * t164 + t291;
t112 = t283 * t121 + t281 * t162;
t111 = t281 * t121 - t283 * t162;
t110 = -pkin(3) * t120 - pkin(4) * t147;
t109 = t277 * t123 + t278 * t124;
t108 = -t278 * t123 + t277 * t124;
t107 = t283 * t125 - t281 * t132 + t294;
t106 = t283 * t126 - t281 * t131 + t295;
t105 = -t281 * t125 - t283 * t132 + t290;
t104 = -t281 * t126 - t283 * t131 + t291;
t103 = -pkin(6) * t120 - qJ(5) * t317 - t280 * t122;
t102 = -pkin(5) * t123 - (pkin(3) * t281 - pkin(6) * t283 + qJ(2)) * t137;
t101 = t277 * t111 + t278 * t112;
t100 = -t278 * t111 + t277 * t112;
t99 = -pkin(5) * t124 - (pkin(3) * t283 + pkin(6) * t281 + t334) * t137;
t98 = -pkin(5) * t111 + qJ(2) * t120 + t283 * t103 - t281 * t110;
t97 = -pkin(5) * t112 - t281 * t103 - t283 * t110 + t334 * t120;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, 0, 0, 0, 0, 0, 0, -t205, t298, 0, -t339, 0, 0, 0, 0, 0, 0, t143, t144, t171, t109, 0, 0, 0, 0, 0, 0, t143, t144, t171, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, 0, 0, 0, 0, 0, 0, t298, t205, 0, t127, 0, 0, 0, 0, 0, 0, t141, t142, t170, t108, 0, 0, 0, 0, 0, 0, t141, t142, t170, t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t275, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t275, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t275, 0, 0, 0, 0, 0, 0, -t218, -t220, 0, t137, 0, 0, 0, 0, 0, 0, -t218, -t220, 0, -t120; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t322, -t267, -t212, -qJ(1) * t212, 0, 0, 0, 0, 0, 0, -t322, -t202, t267, -qJ(1) * t202 + (-pkin(1) * t277 + qJ(2) * t278) * t275, 0, 0, -t298, 0, -t205, 0, t340, t342, t127, -qJ(1) * t127 - t277 * t158 + t278 * t159, t156, t134, t152, t155, t151, t173, -t277 * t113 + t278 * t117 - t328, -t277 * t114 + t278 * t118 - t327, -t277 * t129 + t278 * t130 - t326, -qJ(1) * t108 + t278 * t102 - t277 * t99, t156, t134, t152, t155, t151, t173, -t277 * t104 + t278 * t106 - t328, -t277 * t105 + t278 * t107 - t327, -t277 * t115 + t278 * t116 - t326, -qJ(1) * t100 - t277 * t97 + t278 * t98; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t267, -t322, t213, qJ(1) * t213, 0, 0, 0, 0, 0, 0, t267, t203, t322, qJ(1) * t203 + (pkin(1) * t278 + qJ(2) * t277) * t275, 0, 0, -t205, 0, t298, 0, t342, -t340, t339, -qJ(1) * t339 + t278 * t158 + t277 * t159, t154, t133, t150, t153, t149, t172, t278 * t113 + t277 * t117 + t139, t278 * t114 + t277 * t118 + t140, t278 * t129 + t277 * t130 + t163, qJ(1) * t109 + t277 * t102 + t278 * t99, t154, t133, t150, t153, t149, t172, t278 * t104 + t277 * t106 + t139, t278 * t105 + t277 * t107 + t140, t278 * t115 + t277 * t116 + t163, qJ(1) * t101 + t277 * t98 + t278 * t97; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t254, t255, 0, 0, 0, 0, 0, 0, 0, 0, t248, 0, -t255, pkin(1) * t248 - qJ(2) * t255, 0, 0, 0, 0, 0, -qJDD(3), -qJ(2) * t250 + t334 * t251 + t208, qJ(2) * t251 + t334 * t250 + t209, 0, qJ(2) * t169 + t168 * t334, t216, t198, t219, t215, t217, 0, t288 + t316, t287 - t320, -t138 + t286, pkin(3) * t196 - pkin(6) * t138 + qJ(2) * t124 - t334 * t123, t216, t198, t219, t215, t217, 0, qJ(5) * t319 - t282 * t145 + t288, -t280 * t157 - t282 * t214 + t287, -t282 * t135 - t280 * t146 + t286, pkin(3) * t162 - pkin(6) * t121 + qJ(2) * t112 + qJ(5) * t321 - t334 * t111 - t282 * t122;];
tauB_reg = t1;
