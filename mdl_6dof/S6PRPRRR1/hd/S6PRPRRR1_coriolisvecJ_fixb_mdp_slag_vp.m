% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:25:12
% EndTime: 2019-03-08 20:25:19
% DurationCPUTime: 3.02s
% Computational Cost: add. (2514->294), mult. (6166->424), div. (0->0), fcn. (4856->12), ass. (0->145)
t314 = sin(qJ(5));
t318 = cos(qJ(4));
t394 = cos(qJ(5));
t355 = qJD(2) * t394;
t315 = sin(qJ(4));
t370 = qJD(2) * t315;
t397 = -t314 * t370 + t318 * t355;
t306 = qJD(4) + qJD(5);
t309 = sin(pkin(12));
t310 = sin(pkin(6));
t311 = cos(pkin(12));
t316 = sin(qJ(2));
t319 = cos(qJ(2));
t274 = (t309 * t316 - t311 * t319) * t310;
t396 = (t315 ^ 2 - t318 ^ 2) * MDP(7);
t371 = qJD(1) * t310;
t358 = t316 * t371;
t289 = t311 * t358;
t357 = t319 * t371;
t270 = t309 * t357 + t289;
t369 = qJD(4) * t315;
t345 = pkin(4) * t369 - t270;
t312 = cos(pkin(6));
t298 = qJD(1) * t312 + qJD(3);
t287 = qJD(2) * pkin(2) + t357;
t264 = t309 * t287 + t289;
t351 = t264 + (pkin(8) + pkin(9)) * qJD(2);
t245 = t318 * t298 - t351 * t315;
t272 = qJD(2) * t274;
t269 = qJD(1) * t272;
t224 = t245 * qJD(4) - t318 * t269;
t246 = t298 * t315 + t318 * t351;
t225 = -t246 * qJD(4) + t315 * t269;
t243 = qJD(4) * pkin(4) + t245;
t354 = qJD(5) * t394;
t368 = qJD(5) * t314;
t204 = t394 * t224 + t314 * t225 + t243 * t354 - t246 * t368;
t359 = t394 * t246;
t219 = t314 * t243 + t359;
t350 = t314 * t224 - t394 * t225;
t205 = qJD(5) * t219 + t350;
t381 = t314 * t246;
t218 = t394 * t243 - t381;
t215 = -t306 * pkin(5) - t218;
t288 = t309 * t358;
t263 = t287 * t311 - t288;
t360 = -pkin(4) * t318 - pkin(3);
t254 = qJD(2) * t360 - t263;
t380 = t314 * t318;
t281 = -qJD(2) * t380 - t315 * t355;
t230 = -pkin(5) * t397 + pkin(10) * t281 + t254;
t286 = t394 * t315 + t380;
t393 = pkin(2) * t311;
t291 = t360 - t393;
t329 = -t314 * t315 + t394 * t318;
t249 = -pkin(5) * t329 - pkin(10) * t286 + t291;
t301 = pkin(2) * t309 + pkin(8);
t392 = pkin(9) + t301;
t283 = t392 * t315;
t284 = t392 * t318;
t251 = -t314 * t283 + t394 * t284;
t258 = t306 * t286;
t253 = t258 * qJD(2);
t257 = t306 * t329;
t278 = qJD(6) - t397;
t273 = t311 * t357 - t288;
t352 = qJD(4) * t392;
t276 = t315 * t352;
t277 = t318 * t352;
t330 = -t394 * t283 - t314 * t284;
t377 = -t330 * qJD(5) + t329 * t273 + t394 * t276 + t314 * t277;
t395 = (qJD(6) * t230 + t204) * t329 + t205 * t286 + t215 * t257 + (-qJD(6) * t249 + t377) * t278 - t251 * t253;
t391 = t215 * t397;
t390 = t215 * t286;
t313 = sin(qJ(6));
t367 = qJD(6) * t313;
t252 = t397 * t306;
t317 = cos(qJ(6));
t366 = qJD(6) * t317;
t374 = t317 * t252 + t306 * t366;
t231 = t281 * t367 + t374;
t389 = t231 * t313;
t388 = t249 * t253;
t387 = t252 * t313;
t384 = t281 * t313;
t265 = -t317 * t306 - t384;
t386 = t265 * t278;
t336 = t281 * t317 - t306 * t313;
t385 = t336 * t278;
t382 = t313 * t253;
t379 = t317 * t253;
t378 = -t231 * t329 - t258 * t336;
t376 = t251 * qJD(5) - t286 * t273 - t314 * t276 + t394 * t277;
t375 = pkin(5) * t258 - pkin(10) * t257 + t345;
t372 = MDP(11) * t315;
t365 = qJD(2) * qJD(4);
t364 = pkin(4) * t370;
t362 = t286 * t382;
t361 = t286 * t379;
t353 = t315 * t365;
t349 = t278 * t317;
t261 = -qJD(2) * pkin(3) - t263;
t348 = -qJD(2) * t261 + t269;
t255 = -pkin(5) * t281 - pkin(10) * t397;
t304 = pkin(4) * t314 + pkin(10);
t346 = qJD(6) * t304 + t255 + t364;
t221 = t314 * t245 + t359;
t344 = pkin(4) * t368 - t221;
t216 = t306 * pkin(10) + t219;
t209 = t216 * t317 + t230 * t313;
t343 = t205 * t313 - t209 * t281 + t215 * t366;
t341 = -t253 * t304 - t391;
t340 = t216 * t313 - t230 * t317;
t232 = -qJD(6) * t336 + t387;
t339 = t232 * t329 - t258 * t265;
t275 = (t309 * t319 + t311 * t316) * t310;
t259 = -t275 * t315 + t312 * t318;
t260 = t275 * t318 + t312 * t315;
t234 = t314 * t259 + t394 * t260;
t338 = t234 * t317 + t274 * t313;
t337 = -t234 * t313 + t274 * t317;
t222 = t394 * t245 - t381;
t334 = -pkin(4) * t354 + t222;
t333 = -t205 * t317 + t215 * t367 - t281 * t340;
t332 = t254 * t281 - t350;
t331 = t394 * t259 - t314 * t260;
t328 = -t257 * t313 - t286 * t366;
t327 = -t257 * t317 + t286 * t367;
t271 = qJD(2) * t275;
t268 = qJD(1) * t271;
t320 = qJD(4) ^ 2;
t326 = qJD(2) * t270 - t301 * t320 - t268;
t325 = qJD(4) * (qJD(2) * (-pkin(3) - t393) + t261 + t273);
t256 = pkin(4) * t353 + t268;
t324 = ((t231 - t386) * t317 + (-t232 + t385) * t313) * MDP(21) + (-t336 * t349 + t389) * MDP(20) + (-t278 ^ 2 * t313 - t265 * t281 + t379) * MDP(23) + (t278 * t349 - t281 * t336 + t382) * MDP(22) + t252 * MDP(15) + (t281 ^ 2 - t397 ^ 2) * MDP(14) + (MDP(13) * t397 + t278 * MDP(24)) * t281 + (-t397 * MDP(15) + (-qJD(2) * t286 - t281) * MDP(16)) * t306;
t322 = -t254 * t397 - t204;
t321 = qJD(2) ^ 2;
t305 = -t394 * pkin(4) - pkin(5);
t236 = qJD(4) * t259 - t272 * t318;
t235 = -qJD(4) * t260 + t272 * t315;
t220 = pkin(5) * t253 - pkin(10) * t252 + t256;
t217 = t317 * t220;
t207 = t234 * qJD(5) - t394 * t235 + t314 * t236;
t206 = t331 * qJD(5) + t314 * t235 + t394 * t236;
t1 = [(-t263 * t271 - t264 * t272 + t268 * t274 - t269 * t275) * MDP(5) + (-t207 * t306 + t274 * t253 - t271 * t397) * MDP(18) + (-t206 * t306 + t252 * t274 - t271 * t281) * MDP(19) + ((-qJD(6) * t338 - t206 * t313 + t271 * t317) * t278 + t337 * t253 + t207 * t265 - t331 * t232) * MDP(25) + (-(qJD(6) * t337 + t206 * t317 + t271 * t313) * t278 - t338 * t253 - t207 * t336 - t331 * t231) * MDP(26) + (-MDP(3) * t316 - MDP(4) * t319) * t321 * t310 + (MDP(11) * t235 - MDP(12) * t236) * qJD(4) + ((-t271 * t318 + t274 * t369) * MDP(11) + (qJD(4) * t274 * t318 + t271 * t315) * MDP(12)) * qJD(2); (t263 * t270 - t264 * t273 + (-t268 * t311 - t269 * t309) * pkin(2)) * MDP(5) - 0.2e1 * t365 * t396 + (t252 * t286 - t257 * t281) * MDP(13) + (t252 * t329 - t253 * t286 + t257 * t397 + t258 * t281) * MDP(14) + (t253 * t291 + t254 * t258 - t256 * t329 - t345 * t397) * MDP(18) + (t252 * t291 + t254 * t257 + t256 * t286 - t281 * t345) * MDP(19) + (t231 * t286 * t317 + t327 * t336) * MDP(20) + ((-t265 * t317 + t313 * t336) * t257 + (-t389 - t232 * t317 + (t265 * t313 + t317 * t336) * qJD(6)) * t286) * MDP(21) + (-t278 * t327 + t361 + t378) * MDP(22) + (t278 * t328 + t339 - t362) * MDP(23) + (-t253 * t329 + t258 * t278) * MDP(24) + (-t340 * t258 - t217 * t329 - t330 * t232 + t376 * t265 + (t388 + t375 * t278 + (t216 * t329 - t251 * t278 + t390) * qJD(6)) * t317 + t395 * t313) * MDP(25) + (-t209 * t258 - t330 * t231 - t376 * t336 + (-t388 + (-qJD(6) * t216 + t220) * t329 - qJD(6) * t390 + (qJD(6) * t251 - t375) * t278) * t313 + t395 * t317) * MDP(26) + (t325 * MDP(11) - t326 * MDP(12) - t320 * MDP(9)) * t315 + (t326 * MDP(11) + t325 * MDP(12) + 0.2e1 * MDP(6) * t353 + t320 * MDP(8)) * t318 + (t257 * MDP(15) - t258 * MDP(16) - t376 * MDP(18) + t377 * MDP(19)) * t306; (-t339 - t362) * MDP(25) + (-t361 + t378) * MDP(26) + (-MDP(12) * t318 - t372) * t320 + (-MDP(18) * t258 - MDP(19) * t257) * t306 + (MDP(25) * t328 + MDP(26) * t327) * t278; t324 + t348 * t372 + (t305 * t232 + t341 * t313 + t344 * t265 + (t313 * t334 - t317 * t346) * t278 + t333) * MDP(25) + (t305 * t231 + t341 * t317 - t344 * t336 + (t313 * t346 + t317 * t334) * t278 + t343) * MDP(26) + (t397 * t364 + t221 * t306 + (-t359 + (-pkin(4) * t306 - t243) * t314) * qJD(5) + t332) * MDP(18) + (t222 * t306 + (t281 * t370 - t306 * t354) * pkin(4) + t322) * MDP(19) + t321 * t396 + (-t315 * t321 * MDP(6) + t348 * MDP(12)) * t318; (t332 + (-qJD(5) + t306) * t219) * MDP(18) + (t218 * t306 + t322) * MDP(19) + (-pkin(5) * t232 - (-t218 * t313 + t255 * t317) * t278 - t219 * t265 - t313 * t391 + (-t278 * t366 - t382) * pkin(10) + t333) * MDP(25) + (-pkin(5) * t231 + (t218 * t317 + t255 * t313) * t278 + t219 * t336 - t317 * t391 + (t278 * t367 - t379) * pkin(10) + t343) * MDP(26) + t324; -t336 * t265 * MDP(20) + (-t265 ^ 2 + t336 ^ 2) * MDP(21) + (t374 + t386) * MDP(22) + (-t385 - t387) * MDP(23) + t253 * MDP(24) + (-t204 * t313 + t209 * t278 + t215 * t336 + t217) * MDP(25) + (-t204 * t317 + t215 * t265 - t220 * t313 - t278 * t340) * MDP(26) + (MDP(22) * t384 + MDP(23) * t336 - MDP(25) * t209 + MDP(26) * t340) * qJD(6);];
tauc  = t1;
