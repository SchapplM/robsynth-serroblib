% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPPRR1
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPPRR1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR1_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:19
% EndTime: 2019-12-05 17:38:26
% DurationCPUTime: 2.85s
% Computational Cost: add. (6898->321), mult. (13665->440), div. (0->0), fcn. (7647->6), ass. (0->218)
t316 = sin(qJ(5));
t310 = qJDD(4) + qJDD(5);
t317 = sin(qJ(4));
t319 = cos(qJ(5));
t320 = cos(qJ(4));
t274 = (-t316 * t320 - t317 * t319) * qJD(1);
t352 = qJD(1) * t320;
t276 = -t316 * t317 * qJD(1) + t319 * t352;
t374 = t276 * t274;
t383 = t310 + t374;
t386 = t316 * t383;
t385 = t319 * t383;
t351 = qJD(1) * qJD(4);
t341 = t317 * t351;
t347 = t320 * qJDD(1);
t282 = -t341 + t347;
t340 = t320 * t351;
t349 = t317 * qJDD(1);
t329 = -t340 - t349;
t218 = t274 * qJD(5) + t319 * t282 + t316 * t329;
t311 = qJD(4) + qJD(5);
t373 = t311 * t274;
t384 = t218 + t373;
t323 = qJD(1) ^ 2;
t318 = sin(qJ(1));
t321 = cos(qJ(1));
t293 = t321 * g(1) + t318 * g(2);
t312 = qJDD(1) * qJ(2);
t331 = t293 - t312;
t377 = pkin(1) + qJ(3);
t382 = t377 * t323 - qJDD(3) + t331;
t272 = t274 ^ 2;
t273 = t276 ^ 2;
t309 = t311 ^ 2;
t381 = 2 * qJD(3);
t380 = pkin(2) + pkin(3);
t346 = t321 * qJDD(1);
t285 = -t318 * t323 + t346;
t379 = pkin(5) * t285;
t348 = t318 * qJDD(1);
t286 = t321 * t323 + t348;
t378 = pkin(5) * t286;
t376 = qJ(2) - pkin(6);
t375 = qJDD(1) * pkin(1);
t372 = t311 * t316;
t371 = t311 * t319;
t313 = t317 ^ 2;
t370 = t313 * t323;
t314 = t320 ^ 2;
t369 = t314 * t323;
t291 = qJD(4) * pkin(4) - pkin(7) * t352;
t315 = t323 * pkin(6);
t292 = t318 * g(1) - t321 * g(2);
t332 = -qJDD(2) + t292;
t326 = t323 * qJ(2) + t332;
t338 = t377 * qJDD(1);
t324 = t338 + t326;
t214 = -t329 * pkin(4) - pkin(7) * t370 - t315 + (t291 * t320 + t381) * qJD(1) + t324;
t368 = t316 * t214;
t231 = -t374 + t310;
t367 = t316 * t231;
t354 = t320 * t323;
t350 = qJD(2) * qJD(1);
t307 = 0.2e1 * t350;
t245 = -qJDD(1) * pkin(6) + t307 - t382;
t357 = t320 * t245;
t208 = qJDD(4) * pkin(4) - t282 * pkin(7) + t357 + (-pkin(4) * t354 - pkin(7) * t351 + g(3)) * t317;
t339 = t320 * g(3) - t317 * t245;
t209 = -pkin(4) * t370 + t329 * pkin(7) - qJD(4) * t291 - t339;
t180 = -t319 * t208 + t316 * t209;
t181 = t316 * t208 + t319 * t209;
t157 = -t319 * t180 + t316 * t181;
t366 = t317 * t157;
t345 = qJD(1) * t381;
t257 = t324 + t345;
t244 = -t315 + t257;
t365 = t317 * t244;
t281 = 0.2e1 * t340 + t349;
t364 = t317 * t281;
t344 = t317 * t354;
t289 = qJDD(4) + t344;
t363 = t317 * t289;
t290 = qJDD(4) - t344;
t362 = t317 * t290;
t361 = t319 * t214;
t360 = t319 * t231;
t359 = t320 * t157;
t358 = t320 * t244;
t356 = t320 * t289;
t355 = t320 * t290;
t353 = t313 + t314;
t343 = t318 * t374;
t342 = t321 * t374;
t158 = t316 * t180 + t319 * t181;
t233 = t317 * g(3) + t357;
t194 = -t317 * t233 - t320 * t339;
t258 = -0.2e1 * t350 + t382;
t216 = -t318 * t257 - t321 * t258;
t266 = -t323 * pkin(1) + t307 - t331;
t267 = t326 + t375;
t225 = t321 * t266 - t318 * t267;
t337 = t316 * t282 - t319 * t329;
t248 = -t318 * t292 - t321 * t293;
t336 = t318 * t344;
t335 = t321 * t344;
t334 = t318 * g(3) + t379;
t333 = t321 * g(3) - t378;
t193 = t320 * t233 - t317 * t339;
t215 = t321 * t257 - t318 * t258;
t222 = t318 * t266 + t321 * t267;
t247 = t321 * t292 - t318 * t293;
t330 = -t293 + 0.2e1 * t312 + t307;
t304 = -pkin(2) * t323 + g(3);
t328 = -pkin(2) * t348 + t321 * t304 - t378;
t325 = (-qJD(5) + t311) * t276 - t337;
t322 = qJD(4) ^ 2;
t298 = -t322 - t369;
t297 = t322 - t369;
t296 = -t322 - t370;
t295 = -t322 + t370;
t288 = (-t313 + t314) * t323;
t287 = t353 * t323;
t284 = t353 * qJDD(1);
t283 = -0.2e1 * t341 + t347;
t279 = t353 * t351;
t263 = -t273 + t309;
t262 = t272 - t309;
t261 = t313 * t351 + t320 * t329;
t260 = t317 * t282 + t314 * t351;
t256 = -t273 - t309;
t255 = -t317 * t298 - t356;
t254 = t317 * t295 + t356;
t253 = t320 * t297 + t362;
t252 = t320 * t296 - t362;
t251 = t320 * t298 - t363;
t250 = t317 * t296 + t355;
t249 = t380 * t284;
t246 = -pkin(2) * t346 - t318 * t304 - t379;
t242 = -t321 * t284 + t318 * t287;
t241 = -t318 * t284 - t321 * t287;
t240 = -pkin(2) * t257 + qJ(2) * g(3);
t238 = t320 * t283 - t364;
t236 = -pkin(2) * t258 + t377 * g(3);
t235 = t273 - t272;
t229 = -t309 - t272;
t227 = t321 * t251 - t318 * t283;
t226 = t321 * t250 - t318 * t281;
t224 = t318 * t251 + t321 * t283;
t223 = t318 * t250 + t321 * t281;
t221 = (t274 * t319 + t276 * t316) * t311;
t220 = (t274 * t316 - t276 * t319) * t311;
t219 = -t272 - t273;
t217 = -t276 * qJD(5) - t337;
t213 = t319 * t262 - t367;
t212 = -t316 * t263 + t385;
t211 = t316 * t262 + t360;
t210 = t319 * t263 + t386;
t207 = -t316 * t256 - t360;
t206 = t319 * t256 - t367;
t205 = -t373 + t218;
t200 = (qJD(5) + t311) * t276 + t337;
t198 = t319 * t218 - t276 * t372;
t197 = t316 * t218 + t276 * t371;
t196 = -t316 * t217 - t274 * t371;
t195 = t319 * t217 - t274 * t372;
t192 = t319 * t229 - t386;
t191 = t316 * t229 + t385;
t190 = t320 * t220 + t317 * t221;
t189 = t380 * t287 + t194;
t188 = -t376 * t255 - t380 * t283 + t365;
t187 = -t376 * t252 - t380 * t281 - t358;
t186 = t321 * t193 - t318 * t244;
t185 = t318 * t193 + t321 * t244;
t184 = -pkin(7) * t206 + t361;
t183 = t320 * t211 + t317 * t213;
t182 = t320 * t210 + t317 * t212;
t179 = t380 * t251 - t377 * t255 + t339;
t178 = t380 * t250 - t377 * t252 + t233;
t177 = -t317 * t206 + t320 * t207;
t176 = t320 * t206 + t317 * t207;
t175 = t316 * t205 + t319 * t325;
t174 = -t319 * t200 - t316 * t384;
t173 = -t319 * t205 + t316 * t325;
t172 = -t316 * t200 + t319 * t384;
t170 = -pkin(7) * t191 + t368;
t169 = t320 * t197 + t317 * t198;
t168 = t320 * t195 + t317 * t196;
t167 = -t317 * t191 + t320 * t192;
t166 = t320 * t191 + t317 * t192;
t165 = -pkin(4) * t384 + pkin(7) * t207 + t368;
t164 = -pkin(4) * t200 + pkin(7) * t192 - t361;
t163 = t321 * t176 - t318 * t384;
t162 = t318 * t176 + t321 * t384;
t161 = -t376 * t194 - t380 * t244;
t160 = t321 * t166 - t318 * t200;
t159 = t318 * t166 + t321 * t200;
t156 = -t317 * t173 + t320 * t175;
t155 = t320 * t172 + t317 * t174;
t154 = t320 * t173 + t317 * t175;
t153 = t380 * t193 - t377 * t194;
t152 = t321 * t154 - t318 * t219;
t151 = t318 * t154 + t321 * t219;
t150 = -pkin(4) * t214 + pkin(7) * t158;
t149 = -pkin(7) * t173 - t157;
t148 = -pkin(4) * t219 + pkin(7) * t175 + t158;
t147 = t320 * t158 - t366;
t146 = t317 * t158 + t359;
t145 = t321 * t146 - t318 * t214;
t144 = t318 * t146 + t321 * t214;
t143 = t320 * t165 - t376 * t177 + t317 * t184 - t380 * t384;
t142 = pkin(4) * t206 + t380 * t176 - t377 * t177 - t181;
t141 = t320 * t164 - t376 * t167 + t317 * t170 - t380 * t200;
t140 = pkin(4) * t191 + t380 * t166 - t377 * t167 - t180;
t139 = pkin(4) * t173 + t380 * t154 - t377 * t156;
t138 = t320 * t148 + t317 * t149 - t376 * t156 - t380 * t219;
t137 = -pkin(7) * t366 - t376 * t147 + t320 * t150 - t380 * t214;
t136 = pkin(4) * t157 + t380 * t146 - t377 * t147;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t286, -t285, 0, t248, 0, 0, 0, 0, 0, 0, 0, t286, t285, t225, 0, 0, 0, 0, 0, 0, 0, t285, -t286, t216, 0, 0, 0, 0, 0, 0, t226, t227, t242, t186, 0, 0, 0, 0, 0, 0, t160, t163, t152, t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t285, -t286, 0, t247, 0, 0, 0, 0, 0, 0, 0, -t285, t286, t222, 0, 0, 0, 0, 0, 0, 0, t286, t285, t215, 0, 0, 0, 0, 0, 0, t223, t224, t241, t185, 0, 0, 0, 0, 0, 0, t159, t162, t151, t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t252, t255, 0, t194, 0, 0, 0, 0, 0, 0, t167, t177, t156, t147; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t285, 0, -t286, 0, -t334, -t333, -t247, -pkin(5) * t247, 0, -t285, t286, 0, 0, 0, -t222, t334, t333, -pkin(5) * t222 + (-pkin(1) * t318 + qJ(2) * t321) * g(3), 0, t286, t285, 0, 0, 0, -t215, t328, t246, -pkin(5) * t215 - t318 * t236 + t321 * t240, t321 * t260 - t336, t321 * t238 - t318 * t288, t321 * t253 - t318 * t347, t321 * t261 + t336, t321 * t254 + t317 * t348, -t318 * qJDD(4) - t321 * t279, -pkin(5) * t223 - t318 * t178 + t321 * t187, -pkin(5) * t224 - t318 * t179 + t321 * t188, -pkin(5) * t241 + t321 * t189 + t318 * t249, -pkin(5) * t185 - t318 * t153 + t321 * t161, t321 * t169 + t343, t321 * t155 - t318 * t235, t321 * t182 - t318 * t205, t321 * t168 - t343, t321 * t183 - t318 * t325, t321 * t190 - t318 * t310, -pkin(5) * t159 - t318 * t140 + t321 * t141, -pkin(5) * t162 - t318 * t142 + t321 * t143, -pkin(5) * t151 + t321 * t138 - t318 * t139, -pkin(5) * t144 - t318 * t136 + t321 * t137; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t286, 0, t285, 0, t333, -t334, t248, pkin(5) * t248, 0, -t286, -t285, 0, 0, 0, t225, -t333, t334, pkin(5) * t225 + (pkin(1) * t321 + qJ(2) * t318) * g(3), 0, -t285, t286, 0, 0, 0, t216, -t246, t328, pkin(5) * t216 + t321 * t236 + t318 * t240, t318 * t260 + t335, t318 * t238 + t321 * t288, t318 * t253 + t320 * t346, t318 * t261 - t335, t318 * t254 - t317 * t346, t321 * qJDD(4) - t318 * t279, pkin(5) * t226 + t321 * t178 + t318 * t187, pkin(5) * t227 + t321 * t179 + t318 * t188, pkin(5) * t242 + t318 * t189 - t321 * t249, pkin(5) * t186 + t321 * t153 + t318 * t161, t318 * t169 - t342, t318 * t155 + t321 * t235, t318 * t182 + t321 * t205, t318 * t168 + t342, t318 * t183 + t321 * t325, t318 * t190 + t321 * t310, pkin(5) * t160 + t321 * t140 + t318 * t141, pkin(5) * t163 + t321 * t142 + t318 * t143, pkin(5) * t152 + t318 * t138 + t321 * t139, pkin(5) * t145 + t321 * t136 + t318 * t137; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t292, t293, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t332 - 0.2e1 * t375, t330, pkin(1) * t267 + qJ(2) * t266, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) + t330, t332 + 0.2e1 * t338 + t345, -qJ(2) * t258 + t377 * t257, (t282 - t341) * t320, -t320 * t281 - t317 * t283, -t317 * t297 + t355, t364, t320 * t295 - t363, 0, t376 * t250 + t377 * t281 + t365, t376 * t251 + t377 * t283 + t358, -t376 * t284 - t377 * t287 - t193, t376 * t193 + t377 * t244, -t317 * t197 + t320 * t198, -t317 * t172 + t320 * t174, -t317 * t210 + t320 * t212, -t317 * t195 + t320 * t196, -t317 * t211 + t320 * t213, -t317 * t220 + t320 * t221, -t317 * t164 + t376 * t166 + t320 * t170 + t377 * t200, -t317 * t165 + t376 * t176 + t320 * t184 + t377 * t384, -t317 * t148 + t320 * t149 + t376 * t154 + t377 * t219, -pkin(7) * t359 + t376 * t146 - t317 * t150 + t377 * t214;];
tauB_reg = t1;
