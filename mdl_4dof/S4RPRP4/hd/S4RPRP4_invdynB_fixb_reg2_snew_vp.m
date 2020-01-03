% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPRP4
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPRP4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP4_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:55
% EndTime: 2019-12-31 16:43:58
% DurationCPUTime: 2.69s
% Computational Cost: add. (3746->275), mult. (7735->357), div. (0->0), fcn. (4216->6), ass. (0->202)
t302 = qJD(3) ^ 2;
t298 = sin(qJ(3));
t291 = t298 ^ 2;
t303 = qJD(1) ^ 2;
t353 = t291 * t303;
t275 = t302 + t353;
t300 = cos(qJ(3));
t280 = t300 * t303 * t298;
t272 = qJDD(3) - t280;
t337 = t300 * t272;
t231 = -t298 * t275 + t337;
t331 = qJD(1) * qJD(3);
t322 = t300 * t331;
t329 = t298 * qJDD(1);
t261 = 0.2e1 * t322 + t329;
t295 = sin(pkin(6));
t296 = cos(pkin(6));
t186 = t295 * t231 + t296 * t261;
t189 = t296 * t231 - t295 * t261;
t299 = sin(qJ(1));
t301 = cos(qJ(1));
t146 = t301 * t186 + t299 * t189;
t384 = pkin(4) * t146;
t150 = t299 * t186 - t301 * t189;
t383 = pkin(4) * t150;
t382 = qJ(2) * t186;
t381 = pkin(1) * t186 + pkin(5) * t231;
t344 = t298 * t272;
t225 = t300 * t275 + t344;
t380 = -pkin(1) * t225 + qJ(2) * t189;
t323 = t298 * t331;
t328 = t300 * qJDD(1);
t262 = -0.2e1 * t323 + t328;
t221 = t300 * t262;
t222 = t298 * t261;
t209 = -t221 + t222;
t292 = t300 ^ 2;
t270 = (t291 - t292) * t303;
t178 = t295 * t209 + t296 * t270;
t181 = t296 * t209 - t295 * t270;
t379 = t301 * t178 + t299 * t181;
t378 = t299 * t178 - t301 * t181;
t352 = t292 * t303;
t277 = -t302 + t352;
t229 = -t300 * t277 + t344;
t196 = t295 * t229 + t296 * t328;
t199 = t296 * t229 - t295 * t328;
t377 = t301 * t196 + t299 * t199;
t376 = t299 * t196 - t301 * t199;
t274 = t301 * g(1) + t299 * g(2);
t259 = -t303 * pkin(1) - t274;
t273 = t299 * g(1) - t301 * g(2);
t309 = qJDD(1) * pkin(1) + t273;
t206 = t295 * t259 - t296 * t309;
t207 = t296 * t259 + t295 * t309;
t320 = t295 * t206 + t296 * t207;
t167 = t296 * t206 - t295 * t207;
t336 = t301 * t167;
t374 = -t299 * t320 + t336;
t341 = t299 * t167;
t130 = t301 * t320 + t341;
t264 = t295 * qJDD(1) + t296 * t303;
t265 = t296 * qJDD(1) - t295 * t303;
t212 = -t299 * t264 + t301 * t265;
t293 = g(3) - qJDD(2);
t240 = qJ(2) * t264 - t296 * t293;
t311 = -qJ(2) * t265 - t295 * t293;
t373 = -pkin(4) * t212 + t299 * t240 + t301 * t311;
t372 = 2 * qJD(4);
t370 = pkin(2) * t225;
t369 = pkin(5) * t225;
t363 = t301 * t264 + t299 * t265;
t367 = pkin(4) * t363 + t301 * t240 - t299 * t311;
t359 = pkin(3) * t300;
t314 = -qJ(4) * t298 - t359;
t258 = t314 * qJD(1);
t194 = -t303 * pkin(2) + qJDD(1) * pkin(5) + t207;
t334 = -t300 * t194 + t298 * t293;
t312 = t300 * qJD(1) * t258 + qJDD(3) * qJ(4) + (qJD(3) * t372) - t334;
t362 = t298 * t277 + t337;
t332 = qJD(1) * t298;
t361 = t258 * t332 + qJDD(4);
t278 = -t302 - t352;
t271 = qJDD(3) + t280;
t345 = t298 * t271;
t228 = t300 * t278 - t345;
t185 = t295 * t228 + t296 * t262;
t188 = t296 * t228 - t295 * t262;
t145 = t301 * t185 + t299 * t188;
t358 = pkin(4) * t145;
t333 = t291 + t292;
t266 = t333 * qJDD(1);
t269 = t333 * t303;
t216 = t295 * t266 + t296 * t269;
t217 = t296 * t266 - t295 * t269;
t172 = t301 * t216 + t299 * t217;
t357 = pkin(4) * t172;
t253 = t300 * t271;
t223 = t298 * t278 + t253;
t356 = pkin(5) * t223;
t355 = qJ(2) * t185;
t354 = qJ(2) * t216;
t193 = -qJDD(1) * pkin(2) - t303 * pkin(5) + t206;
t347 = t298 * t193;
t346 = t298 * t262;
t338 = t300 * t193;
t179 = t298 * t194 + t300 * t293;
t335 = t269 - t302;
t327 = pkin(1) * t185 + pkin(2) * t262 + pkin(5) * t228;
t326 = pkin(1) * t216 + pkin(2) * t269 + pkin(5) * t266;
t321 = -pkin(1) * t223 + qJ(2) * t188;
t141 = t298 * t179 - t300 * t334;
t220 = -t299 * t273 - t301 * t274;
t318 = t295 * t280;
t317 = t296 * t280;
t162 = -pkin(2) * t223 + t179;
t268 = t301 * qJDD(1) - t299 * t303;
t315 = -pkin(4) * t268 - t299 * g(3);
t313 = pkin(3) * t298 - qJ(4) * t300;
t140 = t300 * t179 + t298 * t334;
t310 = t300 * t261 + t346;
t219 = t301 * t273 - t299 * t274;
t308 = t322 + t329;
t307 = -t323 + t328;
t306 = -qJDD(3) * pkin(3) + t179 + t361;
t305 = -t307 * pkin(3) + t193 + (-t308 - t322) * qJ(4);
t304 = t332 * t372 - t305;
t276 = t302 - t353;
t267 = t299 * qJDD(1) + t301 * t303;
t256 = t313 * qJDD(1);
t252 = t333 * t331;
t243 = -pkin(4) * t267 + t301 * g(3);
t236 = -t291 * t331 + t300 * t308;
t235 = -t292 * t331 - t298 * t307;
t234 = t295 * qJDD(3) + t296 * t252;
t233 = -t296 * qJDD(3) + t295 * t252;
t230 = -t298 * t276 + t253;
t224 = t300 * t276 + t345;
t210 = qJ(2) * t217;
t204 = t296 * t236 - t318;
t203 = t296 * t235 + t318;
t202 = t295 * t236 + t317;
t201 = t295 * t235 - t317;
t200 = t296 * t230 + t295 * t329;
t197 = t295 * t230 - t296 * t329;
t175 = -t299 * t233 + t301 * t234;
t174 = t301 * t233 + t299 * t234;
t173 = -t299 * t216 + t301 * t217;
t171 = pkin(4) * t173;
t170 = t338 + t369;
t169 = t347 - t356;
t164 = t302 * qJ(4) - t306;
t163 = -t334 + t370;
t161 = -t302 * pkin(3) + t312;
t160 = pkin(1) * t293 + qJ(2) * t320;
t159 = -t299 * t202 + t301 * t204;
t158 = -t299 * t201 + t301 * t203;
t157 = t301 * t202 + t299 * t204;
t156 = t301 * t201 + t299 * t203;
t155 = -t299 * t197 + t301 * t200;
t154 = t301 * t197 + t299 * t200;
t153 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t332 + t305;
t152 = t335 * qJ(4) + t306;
t151 = t335 * pkin(3) + t312;
t148 = -t299 * t185 + t301 * t188;
t144 = pkin(4) * t148;
t143 = (t262 - t323) * pkin(3) + t304;
t142 = -pkin(3) * t323 + qJ(4) * t261 + t304;
t138 = (-t278 - t302) * qJ(4) + (-qJDD(3) - t271) * pkin(3) + t162 + t361;
t137 = -t370 - qJ(4) * t272 + (-t275 + t302) * pkin(3) - t312;
t136 = t296 * t140 - t354;
t135 = t295 * t140 + t210;
t134 = -pkin(3) * t222 + t300 * t142 - t369;
t133 = qJ(4) * t221 - t298 * t143 - t356;
t132 = t296 * t141 + t295 * t193;
t131 = t295 * t141 - t296 * t193;
t128 = t300 * t161 - t298 * t164;
t127 = t298 * t161 + t300 * t164;
t126 = -t298 * t151 + t300 * t152;
t125 = -t295 * t163 + t296 * t170 + t382;
t124 = -t295 * t162 + t296 * t169 - t355;
t123 = t296 * t163 + t295 * t170 - t380;
t122 = t296 * t162 + t295 * t169 + t321;
t121 = t296 * t126 - t295 * t256 - t354;
t120 = t295 * t126 + t296 * t256 + t210;
t119 = t296 * t128 + t295 * t153;
t118 = t295 * t128 - t296 * t153;
t117 = -t299 * t131 + t301 * t132;
t116 = t301 * t131 + t299 * t132;
t115 = -pkin(2) * t127 - pkin(3) * t164 - qJ(4) * t161;
t114 = t296 * t133 - t295 * t138 - t355;
t113 = t296 * t134 - t295 * t137 - t382;
t112 = -pkin(5) * t127 + t313 * t153;
t111 = t295 * t133 + t296 * t138 + t321;
t110 = t295 * t134 + t296 * t137 + t380;
t109 = -qJ(2) * t131 - (pkin(2) * t295 - pkin(5) * t296) * t140;
t108 = qJ(2) * t132 - (-pkin(2) * t296 - pkin(5) * t295 - pkin(1)) * t140;
t107 = -t299 * t118 + t301 * t119;
t106 = t301 * t118 + t299 * t119;
t105 = -qJ(2) * t118 + t296 * t112 - t295 * t115;
t104 = -pkin(1) * t127 + qJ(2) * t119 + t295 * t112 + t296 * t115;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t267, -t268, 0, t220, 0, 0, 0, 0, 0, 0, -t363, -t212, 0, t130, 0, 0, 0, 0, 0, 0, t148, t150, t173, t117, 0, 0, 0, 0, 0, 0, t148, t173, -t150, t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t268, -t267, 0, t219, 0, 0, 0, 0, 0, 0, t212, -t363, 0, -t374, 0, 0, 0, 0, 0, 0, t145, -t146, t172, t116, 0, 0, 0, 0, 0, 0, t145, t172, t146, t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t293, 0, 0, 0, 0, 0, 0, t223, -t225, 0, -t140, 0, 0, 0, 0, 0, 0, t223, 0, t225, t127; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t268, 0, -t267, 0, t315, -t243, -t219, -pkin(4) * t219, 0, 0, t212, 0, -t363, 0, t373, t367, t374, pkin(4) * t374 + qJ(2) * t336 - t299 * t160, t159, t378, t155, t158, t376, t175, -t299 * t122 + t301 * t124 - t358, -t299 * t123 + t301 * t125 + t384, -t299 * t135 + t301 * t136 - t357, -pkin(4) * t116 - t299 * t108 + t301 * t109, t159, t155, -t378, t175, -t376, t158, -t299 * t111 + t301 * t114 - t358, -t299 * t120 + t301 * t121 - t357, -t299 * t110 + t301 * t113 - t384, -pkin(4) * t106 - t299 * t104 + t301 * t105; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t267, 0, t268, 0, t243, t315, t220, pkin(4) * t220, 0, 0, t363, 0, t212, 0, -t367, t373, t130, pkin(4) * t130 + qJ(2) * t341 + t301 * t160, t157, -t379, t154, t156, -t377, t174, t301 * t122 + t299 * t124 + t144, t301 * t123 + t299 * t125 + t383, t301 * t135 + t299 * t136 + t171, pkin(4) * t117 + t301 * t108 + t299 * t109, t157, t154, t379, t174, t377, t156, t301 * t111 + t299 * t114 + t144, t301 * t120 + t299 * t121 + t171, t301 * t110 + t299 * t113 - t383, pkin(4) * t107 + t301 * t104 + t299 * t105; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t273, t274, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t265 - t206, -pkin(1) * t264 - t207, 0, -pkin(1) * t167, t222, t310, t224, t221, t362, 0, t327 - t338, -pkin(2) * t261 + t347 - t381, t141 + t326, pkin(1) * t131 - pkin(2) * t193 + pkin(5) * t141, t222, t224, -t310, 0, -t362, t221, qJ(4) * t346 + t300 * t143 + t327, t300 * t151 + t298 * t152 + t326, t298 * t142 + (pkin(2) + t359) * t261 + t381, pkin(1) * t118 + pkin(5) * t128 + (-pkin(2) + t314) * t153;];
tauB_reg = t1;
