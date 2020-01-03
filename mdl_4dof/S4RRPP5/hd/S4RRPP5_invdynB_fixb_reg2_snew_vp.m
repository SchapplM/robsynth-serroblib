% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRPP5
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRPP5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP5_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_invdynB_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:39
% EndTime: 2019-12-31 17:00:42
% DurationCPUTime: 2.07s
% Computational Cost: add. (2403->274), mult. (5511->308), div. (0->0), fcn. (2575->4), ass. (0->185)
t294 = qJD(2) ^ 2;
t290 = sin(qJ(2));
t286 = t290 ^ 2;
t295 = qJD(1) ^ 2;
t353 = t286 * t295;
t265 = -t294 - t353;
t292 = cos(qJ(2));
t318 = t292 * t295 * t290;
t260 = qJDD(2) - t318;
t336 = t292 * t260;
t219 = t290 * t265 + t336;
t324 = qJD(1) * qJD(2);
t315 = t292 * t324;
t322 = t290 * qJDD(1);
t248 = 0.2e1 * t315 + t322;
t291 = sin(qJ(1));
t293 = cos(qJ(1));
t179 = t291 * t219 + t293 * t248;
t366 = pkin(4) * t179;
t183 = t293 * t219 - t291 * t248;
t176 = pkin(4) * t183;
t287 = t292 ^ 2;
t352 = t287 * t295;
t267 = -t294 - t352;
t259 = qJDD(2) + t318;
t349 = t290 * t259;
t216 = -t292 * t267 + t349;
t274 = t290 * t324;
t320 = t292 * qJDD(1);
t251 = -0.2e1 * t274 + t320;
t178 = t291 * t216 - t293 * t251;
t367 = pkin(4) * t178;
t182 = t293 * t216 + t291 * t251;
t175 = pkin(4) * t182;
t337 = t292 * t259;
t208 = t290 * t267 + t337;
t369 = pkin(1) * t208;
t379 = qJ(3) * t267 + t369;
t364 = pkin(5) * t208;
t331 = pkin(1) * t248 + pkin(5) * t219;
t264 = t294 - t353;
t218 = -t290 * t264 + t337;
t319 = t293 * qJDD(1);
t187 = t291 * t218 - t290 * t319;
t321 = t291 * qJDD(1);
t189 = t293 * t218 + t290 * t321;
t348 = t290 * t260;
t212 = -t292 * t265 + t348;
t378 = pkin(1) * t212;
t363 = pkin(5) * t212;
t377 = pkin(5) * t216;
t361 = t292 * g(3);
t304 = -qJDD(2) * pkin(2) - t294 * qJ(3) + qJDD(3) + t361;
t262 = t293 * g(1) + t291 * g(2);
t233 = -t295 * pkin(1) + qJDD(1) * pkin(5) - t262;
t355 = qJ(3) * t290;
t308 = -pkin(2) * t292 - t355;
t326 = t295 * t308;
t312 = t233 + t326;
t305 = t312 * t290;
t170 = -t304 - t305;
t338 = t292 * t251;
t350 = t290 * t248;
t196 = -t338 + t350;
t257 = (t286 - t287) * t295;
t376 = t291 * t196 + t293 * t257;
t375 = t293 * t196 - t291 * t257;
t266 = -t294 + t352;
t217 = -t292 * t266 + t348;
t186 = t291 * t217 + t292 * t319;
t188 = t293 * t217 - t291 * t320;
t250 = -t274 + t320;
t374 = (t250 + t251) * pkin(2);
t249 = t315 + t322;
t307 = t249 + t315;
t373 = t307 * qJ(3);
t325 = qJD(1) * t290;
t258 = pkin(3) * t325 - qJD(2) * qJ(4);
t372 = t250 * qJ(4) + ((qJ(3) * qJD(2) + (2 * qJD(4))) * t292 + ((2 * qJD(3)) + t258) * t290) * qJD(1);
t371 = t290 * t266 + t336;
t323 = qJD(3) * qJD(2);
t370 = -qJ(3) * t260 - 0.2e1 * t323 - t378;
t327 = t286 + t287;
t253 = t327 * qJDD(1);
t256 = t327 * t295;
t198 = t291 * t253 + t293 * t256;
t365 = pkin(4) * t198;
t362 = t250 * pkin(2);
t360 = t295 * pkin(5);
t359 = pkin(2) + qJ(4);
t358 = qJ(3) * t256;
t354 = qJ(3) * t292;
t261 = t291 * g(1) - t293 * g(2);
t303 = qJDD(1) * pkin(1) + t261;
t232 = t303 + t360;
t351 = t290 * t232;
t340 = t292 * t232;
t339 = t292 * t248;
t332 = pkin(1) * t251 - t377;
t330 = pkin(1) * t256 + pkin(5) * t253;
t328 = t290 * g(3) - t292 * t233;
t222 = t290 * t233 + t361;
t166 = t290 * t222 - t292 * t328;
t204 = -t291 * t261 - t293 * t262;
t311 = t291 * t318;
t310 = t293 * t318;
t255 = -t291 * t295 + t319;
t309 = -pkin(4) * t255 - t291 * g(3);
t306 = t294 * pkin(2) - qJDD(2) * qJ(3) - t292 * t326 + t328;
t165 = t292 * t222 + t290 * t328;
t194 = t290 * t251 + t339;
t211 = t292 * t264 + t349;
t203 = t293 * t261 - t291 * t262;
t279 = 0.2e1 * t323;
t167 = t279 - t306;
t268 = pkin(2) * t274;
t302 = t232 - t268;
t301 = t250 * pkin(3) - qJ(4) * t352 + qJD(2) * t258 + qJDD(4) - t306;
t156 = t279 + t301;
t299 = -0.2e1 * qJD(4) * qJD(2) + t304 - qJ(4) * t259 + (t249 - t315) * pkin(3);
t298 = -0.2e1 * qJD(3) * t325 - t302;
t297 = -t298 + t362;
t155 = t305 + t299;
t296 = t249 * qJ(3) - t268 + t303 + t372;
t276 = pkin(2) * t322;
t254 = t293 * t295 + t321;
t244 = -qJ(3) * t320 + t276;
t239 = t327 * t324;
t231 = -pkin(4) * t254 + t293 * g(3);
t228 = t276 + (qJ(4) * t290 - t354) * qJDD(1);
t227 = t291 * qJDD(2) + t293 * t239;
t226 = t292 * t249 - t286 * t324;
t225 = -t293 * qJDD(2) + t291 * t239;
t224 = -t290 * t250 - t287 * t324;
t207 = t307 * t290;
t206 = (t250 - t274) * t292;
t205 = -pkin(3) * t259 + qJ(3) * t251;
t199 = t293 * t253 - t291 * t256;
t197 = pkin(4) * t199;
t193 = t293 * t226 - t311;
t192 = t293 * t224 + t311;
t191 = t291 * t226 + t310;
t190 = t291 * t224 - t310;
t185 = pkin(3) * t260 + t359 * t248;
t174 = -t340 + t363;
t173 = -t351 - t364;
t169 = -t328 + t378;
t168 = t222 - t369;
t163 = t358 - t170;
t162 = pkin(2) * t256 + t167;
t161 = t297 + t373;
t160 = t298 - t373 - t374;
t159 = (t248 + t307) * qJ(3) + t297;
t158 = t293 * t166 - t291 * t232;
t157 = t291 * t166 + t293 * t232;
t154 = pkin(2) * t259 + t170 + t379;
t153 = pkin(2) * t265 + t306 + t370;
t152 = t358 + (pkin(3) * qJDD(1) + t312) * t290 + t299;
t151 = t362 + (pkin(3) * t287 + pkin(5)) * t295 + t296;
t150 = t292 * t167 - t290 * t170;
t149 = t290 * t167 + t292 * t170;
t148 = pkin(3) * t320 + t359 * t256 + t156;
t147 = -pkin(2) * t350 + t292 * t159 - t363;
t146 = -qJ(3) * t338 - t290 * t160 + t364;
t145 = t362 + (t248 + t249) * qJ(3) + (t265 + t352) * pkin(3) + t302 + t372;
t144 = -t290 * t162 + t292 * t163;
t143 = t360 + qJ(4) * t251 + (t267 + t352) * pkin(3) + t374 + t296;
t142 = t359 * t265 - t301 + t370;
t141 = -t359 * t259 + t155 - t379;
t140 = t290 * t155 + t292 * t156;
t139 = -t292 * t155 + t290 * t156;
t138 = t293 * t150 - t291 * t161;
t137 = t291 * t150 + t293 * t161;
t136 = pkin(3) * t155 + qJ(3) * t151;
t135 = -t290 * t143 + t292 * t205 - t364;
t134 = t292 * t145 - t290 * t185 - t363;
t133 = -pkin(1) * t149 - pkin(2) * t170 - qJ(3) * t167;
t132 = -t290 * t148 + t292 * t152;
t131 = -pkin(5) * t149 + (-pkin(2) * t290 + t354) * t161;
t130 = t293 * t140 - t291 * t151;
t129 = t291 * t140 + t293 * t151;
t128 = pkin(3) * t156 + t359 * t151;
t127 = -pkin(1) * t139 - qJ(3) * t156 + t359 * t155;
t126 = -pkin(5) * t139 - t290 * t128 + t292 * t136;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t254, -t255, 0, t204, 0, 0, 0, 0, 0, 0, -t182, -t183, t199, t158, 0, 0, 0, 0, 0, 0, t199, t182, t183, t138, 0, 0, 0, 0, 0, 0, t199, t183, -t182, t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t255, -t254, 0, t203, 0, 0, 0, 0, 0, 0, -t178, -t179, t198, t157, 0, 0, 0, 0, 0, 0, t198, t178, t179, t137, 0, 0, 0, 0, 0, 0, t198, t179, -t178, t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t208, -t212, 0, -t165, 0, 0, 0, 0, 0, 0, 0, -t208, t212, t149, 0, 0, 0, 0, 0, 0, 0, t212, t208, t139; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t255, 0, -t254, 0, t309, -t231, -t203, -pkin(4) * t203, t193, -t375, t189, t192, -t188, t227, -t291 * t168 + t293 * t173 + t367, -t291 * t169 + t293 * t174 + t366, t293 * t165 - t365, -pkin(4) * t157 - (pkin(1) * t291 - pkin(5) * t293) * t165, t227, -t189, t188, t193, -t375, t192, t293 * t144 - t291 * t244 - t365, t293 * t146 - t291 * t154 - t367, t293 * t147 - t291 * t153 - t366, -pkin(4) * t137 + t293 * t131 - t291 * t133, t227, t188, t189, t192, t375, t193, t293 * t132 - t291 * t228 - t365, t293 * t134 - t291 * t142 - t366, t293 * t135 - t291 * t141 + t367, -pkin(4) * t129 + t293 * t126 - t291 * t127; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t254, 0, t255, 0, t231, t309, t204, pkin(4) * t204, t191, -t376, t187, t190, -t186, t225, t293 * t168 + t291 * t173 - t175, t293 * t169 + t291 * t174 - t176, t291 * t165 + t197, pkin(4) * t158 - (-pkin(1) * t293 - pkin(5) * t291) * t165, t225, -t187, t186, t191, -t376, t190, t291 * t144 + t293 * t244 + t197, t291 * t146 + t293 * t154 + t175, t291 * t147 + t293 * t153 + t176, pkin(4) * t138 + t291 * t131 + t293 * t133, t225, t186, t187, t190, t376, t191, t291 * t132 + t293 * t228 + t197, t291 * t134 + t293 * t142 + t176, t291 * t135 + t293 * t141 - t175, pkin(4) * t130 + t291 * t126 + t293 * t127; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t261, t262, 0, 0, t207, t194, t211, t206, t371, 0, t332 + t340, -t331 - t351, t166 + t330, pkin(1) * t232 + pkin(5) * t166, 0, -t211, -t371, t207, t194, t206, t292 * t162 + t290 * t163 + t330, t377 + t292 * t160 + (-pkin(1) - t355) * t251, pkin(2) * t339 + t290 * t159 + t331, pkin(5) * t150 + (pkin(1) - t308) * t161, 0, -t371, t211, t206, -t194, t207, t292 * t148 + t290 * t152 + t330, t290 * t145 + t292 * t185 + t331, t292 * t143 + t290 * t205 + t332, pkin(1) * t151 + pkin(5) * t140 + t292 * t128 + t290 * t136;];
tauB_reg = t1;
