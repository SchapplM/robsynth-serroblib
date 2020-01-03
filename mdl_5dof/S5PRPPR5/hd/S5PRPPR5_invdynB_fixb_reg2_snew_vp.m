% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PRPPR5
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PRPPR5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR5_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:39
% EndTime: 2019-12-31 17:38:44
% DurationCPUTime: 2.59s
% Computational Cost: add. (6895->311), mult. (11467->445), div. (0->0), fcn. (6254->8), ass. (0->201)
t279 = sin(pkin(7));
t278 = sin(pkin(8));
t280 = cos(pkin(8));
t288 = qJD(2) ^ 2;
t247 = -t278 * qJDD(2) + t280 * t288;
t248 = t280 * qJDD(2) + t278 * t288;
t284 = sin(qJ(2));
t286 = cos(qJ(2));
t294 = -t284 * t247 + t286 * t248;
t335 = t279 * t294;
t250 = t284 * qJDD(2) + t286 * t288;
t281 = cos(pkin(7));
t254 = t279 * g(1) - t281 * g(2);
t300 = -pkin(5) * t250 + t286 * t254;
t334 = t279 * t300;
t333 = t281 * t294;
t332 = t281 * t300;
t190 = t286 * t247 + t284 * t248;
t245 = qJDD(4) + t254;
t198 = qJ(4) * t248 + t278 * t245;
t290 = qJ(4) * t247 + t280 * t245;
t139 = -pkin(5) * t190 + t284 * t198 + t286 * t290;
t138 = -pkin(5) * t294 + t286 * t198 - t284 * t290;
t255 = t281 * g(1) + t279 * g(2);
t276 = g(3) - qJDD(1);
t226 = -t286 * t255 - t284 * t276;
t272 = qJDD(2) * qJ(3);
t306 = (qJD(3) * qJD(2));
t291 = (2 * t306) + t272 + t226;
t326 = pkin(2) + pkin(3);
t182 = -t288 * t326 + t291;
t225 = -t284 * t255 + t286 * t276;
t277 = qJDD(2) * pkin(2);
t203 = -t288 * qJ(3) + qJDD(3) + t225 - t277;
t289 = -qJDD(2) * pkin(3) + t203;
t149 = t278 * t182 - t280 * t289;
t150 = t280 * t182 + t278 * t289;
t116 = t280 * t149 - t278 * t150;
t117 = t278 * t149 + t280 * t150;
t330 = t286 * t116 + t284 * t117;
t329 = t284 * t116 - t286 * t117;
t251 = t286 * qJDD(2) - t284 * t288;
t295 = -pkin(1) * t250 + qJ(1) * t251 - t226;
t236 = t281 * t254;
t204 = -t279 * t255 + t236;
t324 = qJ(1) * t250;
t232 = t279 * t250;
t233 = t279 * t251;
t322 = t279 * t254;
t320 = t279 * t276;
t234 = t281 * t250;
t235 = t281 * t251;
t319 = t281 * t276;
t147 = qJDD(2) * pkin(4) - t288 * pkin(6) + t149;
t283 = sin(qJ(5));
t318 = t283 * t147;
t285 = cos(qJ(5));
t262 = t283 * t288 * t285;
t256 = qJDD(5) + t262;
t317 = t283 * t256;
t257 = qJDD(5) - t262;
t316 = t283 * t257;
t313 = t285 * t147;
t312 = t285 * t256;
t311 = t285 * t257;
t274 = t283 ^ 2;
t310 = t288 * t274;
t148 = -t288 * pkin(4) - qJDD(2) * pkin(6) + t150;
t132 = t285 * t148 + t283 * t245;
t206 = pkin(5) * t251 + t284 * t254;
t309 = -qJ(1) * t234 - t279 * t206;
t275 = t285 ^ 2;
t308 = -t274 - t275;
t307 = qJD(2) * qJD(5);
t305 = t281 * qJDD(2);
t304 = t283 * qJDD(2);
t303 = t285 * qJDD(2);
t302 = t283 * t307;
t301 = t285 * t307;
t131 = t283 * t148 - t285 * t245;
t111 = t283 * t131 + t285 * t132;
t249 = t308 * qJDD(2);
t269 = t275 * t288;
t252 = t269 + t310;
t199 = t278 * t249 + t280 * t252;
t202 = t280 * t249 - t278 * t252;
t156 = -t286 * t199 + t284 * t202;
t157 = t284 * t199 + t286 * t202;
t299 = -pkin(1) * t156 + pkin(4) * t252 + pkin(6) * t249 + qJ(1) * t157 - qJ(3) * t202 + t199 * t326 + t111;
t298 = -pkin(1) * t190 + qJ(1) * t294 - qJ(3) * t248 - t247 * t326 - t150;
t297 = -pkin(1) * t294 - qJ(1) * t190 + qJ(3) * t247 - t248 * t326 - t149;
t296 = -0.2e1 * t272 - (2 * t306) + t295;
t189 = -t288 * pkin(2) + t291;
t153 = t286 * t189 + t284 * t203;
t163 = t284 * t225 + t286 * t226;
t205 = -t281 * t255 - t322;
t293 = t278 * t262;
t292 = t280 * t262;
t110 = t285 * t131 - t283 * t132;
t151 = t284 * t189 - t286 * t203;
t162 = t286 * t225 - t284 * t226;
t186 = -pkin(1) * t251 + t225;
t287 = qJD(5) ^ 2;
t267 = t279 * qJDD(2);
t261 = -t269 - t287;
t260 = t269 - t287;
t259 = -t287 - t310;
t258 = t287 - t310;
t253 = t269 - t310;
t246 = pkin(1) * t254;
t244 = t302 - t303;
t243 = 0.2e1 * t302 - t303;
t242 = -t301 - t304;
t241 = 0.2e1 * t301 + t304;
t240 = t308 * t307;
t224 = t285 * t242 + t274 * t307;
t223 = -t283 * t244 + t275 * t307;
t222 = t278 * qJDD(5) + t280 * t240;
t221 = t280 * qJDD(5) - t278 * t240;
t219 = -t283 * t259 - t311;
t218 = -t283 * t258 + t312;
t217 = t285 * t261 - t317;
t216 = t285 * t260 - t316;
t215 = t285 * t259 - t316;
t214 = t285 * t258 + t317;
t213 = t283 * t261 + t312;
t212 = t283 * t260 + t311;
t211 = (t242 - t301) * t283;
t210 = (t244 + t302) * t285;
t196 = t281 * t206;
t188 = t283 * t241 + t285 * t243;
t187 = -t285 * t241 + t283 * t243;
t184 = t281 * t190;
t183 = t279 * t190;
t178 = t280 * t224 - t293;
t177 = t280 * t223 + t293;
t176 = -t278 * t224 - t292;
t175 = -t278 * t223 + t292;
t174 = t280 * t218 - t278 * t304;
t173 = t280 * t216 - t278 * t303;
t172 = -t278 * t218 - t280 * t304;
t171 = -t278 * t216 - t280 * t303;
t169 = qJDD(3) + t186 - 0.2e1 * t277;
t167 = t280 * t219 - t278 * t241;
t166 = t280 * t217 - t278 * t243;
t165 = t278 * t219 + t280 * t241;
t164 = t278 * t217 + t280 * t243;
t160 = t280 * t188 - t278 * t253;
t159 = -t278 * t188 - t280 * t253;
t158 = -t284 * t221 + t286 * t222;
t155 = t281 * t163 - t322;
t154 = t279 * t163 + t236;
t146 = -t284 * t176 + t286 * t178;
t145 = -t284 * t175 + t286 * t177;
t144 = -t284 * t172 + t286 * t174;
t143 = -t284 * t171 + t286 * t173;
t141 = t281 * t153 - t322;
t140 = t279 * t153 + t236;
t137 = t284 * t165 + t286 * t167;
t136 = t284 * t164 + t286 * t166;
t135 = -t286 * t165 + t284 * t167;
t134 = -t286 * t164 + t284 * t166;
t133 = -pkin(5) * t151 + (-pkin(2) * t284 + qJ(3) * t286) * t254;
t129 = -t284 * t159 + t286 * t160;
t128 = -pkin(6) * t215 + t313;
t127 = -pkin(6) * t213 + t318;
t126 = -pkin(4) * t215 + t132;
t125 = -pkin(4) * t213 + t131;
t124 = t281 * t137 - t279 * t215;
t123 = t281 * t136 - t279 * t213;
t122 = t279 * t137 + t281 * t215;
t121 = t279 * t136 + t281 * t213;
t120 = -pkin(1) * t151 + pkin(2) * t203 - qJ(3) * t189;
t113 = qJ(3) * t245 + qJ(4) * t116;
t112 = -qJ(4) * t117 + t245 * t326;
t108 = -qJ(4) * t199 + t280 * t110;
t107 = -qJ(4) * t202 - t278 * t110;
t106 = t280 * t111 + t278 * t147;
t105 = t278 * t111 - t280 * t147;
t104 = qJ(3) * t215 - qJ(4) * t165 - t278 * t126 + t280 * t128;
t103 = qJ(3) * t213 - qJ(4) * t164 - t278 * t125 + t280 * t127;
t99 = -t279 * t245 - t281 * t329;
t98 = t281 * t245 - t279 * t329;
t97 = -qJ(4) * t167 - t280 * t126 - t278 * t128 + t215 * t326;
t96 = -qJ(4) * t166 - t280 * t125 - t278 * t127 + t213 * t326;
t95 = -pkin(1) * t135 + pkin(4) * t241 + pkin(6) * t219 - qJ(3) * t167 + t165 * t326 + t318;
t94 = -pkin(1) * t134 + pkin(4) * t243 + pkin(6) * t217 - qJ(3) * t166 + t164 * t326 - t313;
t92 = -pkin(5) * t156 - t284 * t107 + t286 * t108;
t91 = t284 * t105 + t286 * t106;
t90 = -t286 * t105 + t284 * t106;
t89 = -pkin(5) * t330 - t284 * t112 + t286 * t113;
t88 = -pkin(5) * t135 + t286 * t104 - t284 * t97;
t87 = -pkin(5) * t134 + t286 * t103 - t284 * t96;
t86 = -pkin(1) * t330 - qJ(3) * t117 - t116 * t326;
t85 = t110 * t279 + t281 * t91;
t84 = -t110 * t281 + t279 * t91;
t83 = -qJ(4) * t105 - (pkin(4) * t278 - pkin(6) * t280 + qJ(3)) * t110;
t82 = -qJ(4) * t106 - (pkin(4) * t280 + pkin(6) * t278 + t326) * t110;
t81 = -pkin(1) * t90 - pkin(4) * t147 + pkin(6) * t111 - qJ(3) * t106 + t105 * t326;
t80 = -pkin(5) * t90 - t284 * t82 + t286 * t83;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, 0, 0, 0, 0, 0, 0, -t234, -t235, 0, t155, 0, 0, 0, 0, 0, 0, -t234, 0, t235, t141, 0, 0, 0, 0, 0, 0, -t184, t333, 0, t99, 0, 0, 0, 0, 0, 0, t123, t124, t281 * t157, t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, 0, 0, 0, 0, 0, 0, -t232, -t233, 0, t154, 0, 0, 0, 0, 0, 0, -t232, 0, t233, t140, 0, 0, 0, 0, 0, 0, -t183, t335, 0, t98, 0, 0, 0, 0, 0, 0, t121, t122, t279 * t157, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t276, 0, 0, 0, 0, 0, 0, t251, -t250, 0, -t162, 0, 0, 0, 0, 0, 0, t251, 0, t250, t151, 0, 0, 0, 0, 0, 0, t294, t190, 0, t330, 0, 0, 0, 0, 0, 0, t134, t135, t156, t90; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t320, -t319, -t204, -qJ(1) * t204, 0, 0, t235, 0, -t234, t267, -t196 + (-t186 + t324) * t279, t279 * t295 - t332, t281 * t162, -qJ(1) * t154 - (pkin(1) * t279 - pkin(5) * t281) * t162, 0, t235, 0, t267, t234, 0, -t196 + (-t169 + t324) * t279, -t281 * t151, -t279 * t296 + t332, -qJ(1) * t140 - t279 * t120 + t281 * t133, 0, 0, -t333, 0, -t184, t267, t281 * t138 - t279 * t297, t281 * t139 - t279 * t298, t281 * t330, -qJ(1) * t98 - t279 * t86 + t281 * t89, t281 * t146 - t279 * t211, t281 * t129 - t279 * t187, t281 * t144 - t279 * t214, t281 * t145 - t279 * t210, t281 * t143 - t279 * t212, t281 * t158, -qJ(1) * t121 - t279 * t94 + t281 * t87, -qJ(1) * t122 - t279 * t95 + t281 * t88, -t279 * t299 + t281 * t92, -qJ(1) * t84 - t279 * t81 + t281 * t80; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t319, -t320, t205, qJ(1) * t205, 0, 0, t233, 0, -t232, -t305, t281 * t186 + t309, -t281 * t295 - t334, t279 * t162, qJ(1) * t155 - (-pkin(1) * t281 - pkin(5) * t279) * t162, 0, t233, 0, -t305, t232, 0, t281 * t169 + t309, -t279 * t151, t281 * t296 + t334, qJ(1) * t141 + t281 * t120 + t279 * t133, 0, 0, -t335, 0, -t183, -t305, t279 * t138 + t281 * t297, t279 * t139 + t281 * t298, t279 * t330, qJ(1) * t99 + t279 * t89 + t281 * t86, t279 * t146 + t281 * t211, t279 * t129 + t281 * t187, t279 * t144 + t281 * t214, t279 * t145 + t281 * t210, t279 * t143 + t281 * t212, t279 * t158, qJ(1) * t123 + t279 * t87 + t281 * t94, qJ(1) * t124 + t279 * t88 + t281 * t95, t279 * t92 + t281 * t299, qJ(1) * t85 + t279 * t80 + t281 * t81; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t254, t255, 0, 0, 0, 0, t250, 0, t251, 0, t300, -t206, t163, pkin(5) * t163 + t246, 0, t250, 0, 0, -t251, 0, t300, t153, t206, pkin(5) * t153 + t246 + (pkin(2) * t286 + qJ(3) * t284) * t254, 0, 0, -t190, 0, t294, 0, t139, -t138, t329, pkin(1) * t245 - pkin(5) * t329 + t286 * t112 + t284 * t113, t286 * t176 + t284 * t178, t286 * t159 + t284 * t160, t286 * t172 + t284 * t174, t286 * t175 + t284 * t177, t286 * t171 + t284 * t173, t286 * t221 + t284 * t222, pkin(1) * t213 + pkin(5) * t136 + t284 * t103 + t286 * t96, pkin(1) * t215 + pkin(5) * t137 + t284 * t104 + t286 * t97, pkin(5) * t157 + t286 * t107 + t284 * t108, -pkin(1) * t110 + pkin(5) * t91 + t284 * t83 + t286 * t82;];
tauB_reg = t1;
