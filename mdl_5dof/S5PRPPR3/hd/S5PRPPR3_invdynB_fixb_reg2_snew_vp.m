% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PRPPR3
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PRPPR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:57
% EndTime: 2019-12-05 15:27:04
% DurationCPUTime: 4.22s
% Computational Cost: add. (6621->316), mult. (11130->453), div. (0->0), fcn. (7200->8), ass. (0->212)
t283 = sin(pkin(7));
t285 = cos(pkin(7));
t262 = g(1) * t283 - g(2) * t285;
t253 = -qJDD(3) + t262;
t282 = sin(pkin(8));
t284 = cos(pkin(8));
t293 = qJD(2) ^ 2;
t254 = qJDD(2) * t284 - t282 * t293;
t202 = qJ(3) * t254 + t253 * t282;
t289 = sin(qJ(2));
t291 = cos(qJ(2));
t255 = qJDD(2) * t282 + t284 * t293;
t300 = t254 * t289 + t255 * t291;
t306 = -qJ(3) * t255 + t253 * t284;
t143 = -pkin(5) * t300 - t202 * t289 + t291 * t306;
t374 = t283 * t143;
t373 = t285 * t143;
t302 = t254 * t291 - t255 * t289;
t353 = pkin(5) * t302 + t202 * t291 + t289 * t306;
t372 = t283 * t353;
t371 = t285 * t353;
t369 = t283 * t300;
t368 = t283 * t302;
t367 = t285 * t300;
t366 = t285 * t302;
t263 = g(1) * t285 + g(2) * t283;
t280 = g(3) - qJDD(1);
t232 = -t263 * t289 + t280 * t291;
t226 = qJDD(2) * pkin(2) - t232;
t233 = -t263 * t291 - t280 * t289;
t227 = -pkin(2) * t293 + t233;
t162 = -t226 * t284 + t227 * t282;
t163 = t226 * t282 + t227 * t284;
t303 = t162 * t282 + t163 * t284;
t130 = t162 * t284 - t163 * t282;
t327 = t291 * t130;
t106 = -t289 * t303 + t327;
t335 = t289 * t130;
t107 = t291 * t303 + t335;
t309 = -pkin(1) * t302 - pkin(2) * t254 - qJ(1) * t300 + t162;
t365 = pkin(1) * t300 + pkin(2) * t255 - qJ(1) * t302;
t294 = (2 * qJD(4) * qJD(2)) + t163;
t320 = qJDD(2) * qJ(4);
t157 = -pkin(3) * t293 + t294 + t320;
t281 = qJDD(2) * pkin(3);
t160 = -qJ(4) * t293 + qJDD(4) + t162 - t281;
t120 = t157 * t282 - t160 * t284;
t304 = t157 * t284 + t160 * t282;
t100 = -t120 * t289 + t291 * t304;
t98 = t120 * t291 + t289 * t304;
t243 = t285 * t262;
t210 = -t263 * t283 + t243;
t352 = pkin(3) + pkin(6);
t288 = sin(qJ(5));
t278 = t288 ^ 2;
t349 = t278 * t293;
t290 = cos(qJ(5));
t279 = t290 ^ 2;
t348 = t279 * t293;
t322 = t278 + t279;
t257 = t322 * qJDD(2);
t347 = t282 * t257;
t346 = t283 * t253;
t259 = qJDD(2) * t291 - t289 * t293;
t345 = t283 * t259;
t344 = t283 * t262;
t342 = t283 * t280;
t340 = t284 * t257;
t339 = t285 * t259;
t338 = t285 * t280;
t153 = -pkin(6) * t293 + t157;
t337 = t288 * t153;
t316 = t290 * t293 * t288;
t264 = qJDD(5) + t316;
t336 = t288 * t264;
t330 = t290 * t153;
t329 = t290 * t264;
t265 = qJDD(5) - t316;
t328 = t290 * t265;
t321 = qJD(2) * qJD(5);
t319 = t285 * qJDD(2);
t318 = t288 * qJDD(2);
t317 = t290 * qJDD(2);
t315 = t288 * t321;
t314 = t290 * t321;
t154 = -qJDD(2) * pkin(6) + t160;
t146 = -t154 * t290 - t253 * t288;
t147 = t154 * t288 - t253 * t290;
t112 = -t290 * t146 + t288 * t147;
t260 = t322 * t293;
t204 = -t260 * t282 + t340;
t209 = -t260 * t284 - t347;
t158 = t204 * t291 + t209 * t289;
t159 = -t204 * t289 + t209 * t291;
t313 = -pkin(1) * t158 - pkin(2) * t204 + qJ(1) * t159 + qJ(4) * t260 - t257 * t352 + t112;
t312 = -t294 - 0.2e1 * t320 - t365;
t311 = -qJDD(4) + 0.2e1 * t281 - t309;
t310 = t163 + t365;
t258 = qJDD(2) * t289 + t291 * t293;
t308 = -pkin(1) * t258 + qJ(1) * t259 - t233;
t307 = pkin(1) * t259 + qJ(1) * t258 - t232;
t169 = t232 * t289 + t233 * t291;
t211 = -t263 * t285 - t344;
t292 = qJD(5) ^ 2;
t298 = -t292 - t348;
t297 = t282 * t316;
t296 = t284 * t316;
t214 = pkin(5) * t258 - t262 * t291;
t213 = -pkin(5) * t259 - t262 * t289;
t113 = t146 * t288 + t147 * t290;
t168 = t232 * t291 - t233 * t289;
t272 = t283 * qJDD(2);
t268 = t292 - t348;
t267 = -t292 - t349;
t266 = -t292 + t349;
t261 = (-t278 + t279) * t293;
t252 = -0.2e1 * t315 + t317;
t251 = -t315 + t317;
t250 = -t314 - t318;
t249 = 0.2e1 * t314 + t318;
t248 = pkin(1) * t253;
t247 = t288 * t265;
t246 = t322 * t321;
t242 = t285 * t258;
t241 = t283 * t258;
t236 = t285 * t253;
t231 = -t251 * t288 - t279 * t321;
t230 = -t250 * t290 - t278 * t321;
t229 = qJDD(5) * t284 - t246 * t282;
t228 = qJDD(5) * t282 + t246 * t284;
t224 = -t288 * t298 - t329;
t223 = t268 * t288 - t328;
t222 = (-t251 + t315) * t290;
t221 = t267 * t290 - t247;
t220 = -t266 * t290 + t336;
t219 = t290 * t298 - t336;
t218 = -t268 * t290 - t247;
t217 = t267 * t288 + t328;
t216 = -t266 * t288 - t329;
t215 = (t250 - t314) * t288;
t193 = t249 * t290 + t252 * t288;
t192 = t249 * t288 - t252 * t290;
t181 = -t230 * t282 - t296;
t180 = -t231 * t282 + t296;
t179 = t230 * t284 - t297;
t178 = t231 * t284 + t297;
t177 = -t218 * t282 + t284 * t317;
t176 = -t216 * t282 - t284 * t318;
t175 = t218 * t284 + t282 * t317;
t174 = t216 * t284 - t282 * t318;
t173 = t219 * t282 + t252 * t284;
t172 = t217 * t282 + t249 * t284;
t171 = -t219 * t284 + t252 * t282;
t170 = -t217 * t284 + t249 * t282;
t166 = -t192 * t282 + t261 * t284;
t165 = t192 * t284 + t261 * t282;
t164 = -t228 * t289 + t229 * t291;
t156 = t169 * t285 - t344;
t155 = t169 * t283 + t243;
t151 = -t179 * t289 + t181 * t291;
t150 = -t178 * t289 + t180 * t291;
t149 = -t175 * t289 + t177 * t291;
t148 = -t174 * t289 + t176 * t291;
t140 = -t171 * t289 + t173 * t291;
t139 = -t170 * t289 + t172 * t291;
t138 = t171 * t291 + t173 * t289;
t137 = t170 * t291 + t172 * t289;
t133 = -t165 * t289 + t166 * t291;
t127 = t140 * t285 + t224 * t283;
t126 = t139 * t285 + t221 * t283;
t125 = t140 * t283 - t224 * t285;
t124 = t139 * t283 - t221 * t285;
t123 = pkin(2) * t253 + qJ(3) * t303;
t118 = pkin(4) * t219 - qJ(4) * t224 - t147;
t117 = pkin(4) * t249 - t221 * t352 + t330;
t116 = pkin(4) * t217 - qJ(4) * t221 - t146;
t115 = pkin(4) * t252 - t224 * t352 - t337;
t114 = -qJ(3) * t120 + (-pkin(3) * t282 + qJ(4) * t284) * t253;
t111 = qJ(3) * t304 + (pkin(3) * t284 + qJ(4) * t282 + pkin(2)) * t253;
t110 = -pkin(4) * t260 - t113;
t109 = -pkin(4) * t340 - qJ(3) * t204 - t110 * t282;
t108 = -pkin(4) * t347 + qJ(3) * t209 + t110 * t284;
t104 = t112 * t282 + t153 * t284;
t103 = -t112 * t284 + t153 * t282;
t102 = t107 * t285 - t346;
t101 = t107 * t283 + t236;
t97 = -pkin(1) * t138 - pkin(2) * t171 - qJ(4) * t252 + t219 * t352 - t330;
t96 = -pkin(1) * t137 - pkin(2) * t170 - qJ(4) * t249 + t217 * t352 - t337;
t95 = t100 * t285 - t346;
t94 = t100 * t283 + t236;
t92 = -qJ(3) * t171 - t115 * t282 + t118 * t284;
t91 = -qJ(3) * t170 + t116 * t284 - t117 * t282;
t90 = -pkin(2) * t224 + qJ(3) * t173 + t115 * t284 + t118 * t282;
t89 = -pkin(2) * t221 + qJ(3) * t172 + t116 * t282 + t117 * t284;
t88 = pkin(4) * t112 - qJ(4) * t113;
t87 = pkin(1) * t106 + pkin(2) * t130;
t86 = pkin(4) * t153 - t113 * t352;
t85 = pkin(5) * t106 + qJ(3) * t327 - t123 * t289;
t84 = -t103 * t289 + t104 * t291;
t83 = t103 * t291 + t104 * t289;
t82 = -pkin(1) * t98 - pkin(2) * t120 + pkin(3) * t160 - qJ(4) * t157;
t81 = -pkin(5) * t158 - t108 * t289 + t109 * t291;
t80 = -pkin(5) * t98 - t111 * t289 + t114 * t291;
t79 = t113 * t283 + t285 * t84;
t78 = -t113 * t285 + t283 * t84;
t77 = -pkin(5) * t138 - t289 * t90 + t291 * t92;
t76 = -pkin(5) * t137 - t289 * t89 + t291 * t91;
t75 = -qJ(3) * t103 - t282 * t86 + t284 * t88;
t74 = -pkin(1) * t83 - pkin(2) * t103 - qJ(4) * t153 + t112 * t352;
t73 = -pkin(2) * t113 + qJ(3) * t104 + t282 * t88 + t284 * t86;
t72 = -pkin(5) * t83 - t289 * t73 + t291 * t75;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, 0, 0, 0, 0, 0, 0, -t242, -t339, 0, t156, 0, 0, 0, 0, 0, 0, -t367, -t366, 0, t102, 0, 0, 0, 0, 0, 0, 0, t367, t366, t95, 0, 0, 0, 0, 0, 0, t126, t127, t285 * t159, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, 0, 0, 0, 0, 0, 0, -t241, -t345, 0, t155, 0, 0, 0, 0, 0, 0, -t369, -t368, 0, t101, 0, 0, 0, 0, 0, 0, 0, t369, t368, t94, 0, 0, 0, 0, 0, 0, t124, t125, t283 * t159, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t280, 0, 0, 0, 0, 0, 0, t259, -t258, 0, -t168, 0, 0, 0, 0, 0, 0, t302, -t300, 0, -t106, 0, 0, 0, 0, 0, 0, 0, -t302, t300, t98, 0, 0, 0, 0, 0, 0, t137, t138, t158, t83; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t342, -t338, -t210, -qJ(1) * t210, 0, 0, t339, 0, -t242, t272, t285 * t213 + t283 * t307, t285 * t214 + t283 * t308, t285 * t168, -qJ(1) * t155 - (pkin(1) * t283 - pkin(5) * t285) * t168, 0, 0, t366, 0, -t367, t272, -t283 * t309 - t371, -t283 * t310 - t373, t285 * t106, -qJ(1) * t101 - t283 * t87 + t285 * t85, t272, -t366, t367, 0, 0, 0, -t285 * t98, -t283 * t311 + t371, -t283 * t312 + t373, -qJ(1) * t94 - t283 * t82 + t285 * t80, t150 * t285 - t222 * t283, t133 * t285 - t193 * t283, t149 * t285 - t223 * t283, t151 * t285 - t215 * t283, t148 * t285 - t220 * t283, t285 * t164, -qJ(1) * t124 - t283 * t96 + t285 * t76, -qJ(1) * t125 - t283 * t97 + t285 * t77, -t283 * t313 + t285 * t81, -qJ(1) * t78 - t283 * t74 + t285 * t72; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t338, -t342, t211, qJ(1) * t211, 0, 0, t345, 0, -t241, -t319, t283 * t213 - t285 * t307, t283 * t214 - t285 * t308, t283 * t168, qJ(1) * t156 - (-pkin(1) * t285 - pkin(5) * t283) * t168, 0, 0, t368, 0, -t369, -t319, t285 * t309 - t372, t285 * t310 - t374, t283 * t106, qJ(1) * t102 + t283 * t85 + t285 * t87, -t319, -t368, t369, 0, 0, 0, -t283 * t98, t285 * t311 + t372, t285 * t312 + t374, qJ(1) * t95 + t283 * t80 + t285 * t82, t150 * t283 + t222 * t285, t133 * t283 + t193 * t285, t149 * t283 + t223 * t285, t151 * t283 + t215 * t285, t148 * t283 + t220 * t285, t283 * t164, qJ(1) * t126 + t283 * t76 + t285 * t96, qJ(1) * t127 + t283 * t77 + t285 * t97, t283 * t81 + t285 * t313, qJ(1) * t79 + t283 * t72 + t285 * t74; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t262, t263, 0, 0, 0, 0, t258, 0, t259, 0, -t214, t213, t169, pkin(1) * t262 + pkin(5) * t169, 0, 0, t300, 0, t302, 0, t143, -t353, t107, pkin(5) * t107 + qJ(3) * t335 + t123 * t291 + t248, 0, -t300, -t302, 0, 0, 0, t100, -t143, t353, pkin(5) * t100 + t111 * t291 + t114 * t289 + t248, t178 * t291 + t180 * t289, t165 * t291 + t166 * t289, t175 * t291 + t177 * t289, t179 * t291 + t181 * t289, t174 * t291 + t176 * t289, t228 * t291 + t229 * t289, -pkin(1) * t221 + pkin(5) * t139 + t289 * t91 + t291 * t89, -pkin(1) * t224 + pkin(5) * t140 + t289 * t92 + t291 * t90, pkin(5) * t159 + t108 * t291 + t109 * t289, -pkin(1) * t113 + pkin(5) * t84 + t289 * t75 + t291 * t73;];
tauB_reg = t1;
