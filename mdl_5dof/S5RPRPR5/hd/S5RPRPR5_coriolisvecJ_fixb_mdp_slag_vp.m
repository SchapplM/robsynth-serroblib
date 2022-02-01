% Calculate Coriolis joint torque vector for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:26:04
% EndTime: 2022-01-23 09:26:10
% DurationCPUTime: 2.98s
% Computational Cost: add. (1863->254), mult. (5179->367), div. (0->0), fcn. (3720->8), ass. (0->134)
t313 = cos(pkin(9));
t312 = sin(pkin(8));
t318 = cos(qJ(3));
t356 = qJD(3) * t318;
t345 = t312 * t356;
t335 = qJD(1) * t345;
t311 = sin(pkin(9));
t316 = sin(qJ(3));
t362 = qJD(1) * t312;
t348 = t316 * t362;
t336 = t311 * t348;
t266 = qJD(3) * t336 - t313 * t335;
t288 = t311 * t318 + t313 * t316;
t278 = t288 * t312;
t275 = qJD(3) * t278;
t267 = qJD(1) * t275;
t347 = t318 * t362;
t272 = t313 * t347 - t336;
t315 = sin(qJ(5));
t317 = cos(qJ(5));
t354 = qJD(5) * t315;
t325 = qJD(1) * t288;
t268 = t312 * t325;
t370 = t317 * t268;
t207 = -qJD(5) * t370 + t315 * t266 - t317 * t267 - t272 * t354;
t330 = -t268 * t315 + t317 * t272;
t208 = qJD(5) * t330 - t317 * t266 - t267 * t315;
t232 = t272 * t315 + t370;
t314 = cos(pkin(8));
t361 = qJD(1) * t314;
t396 = qJD(3) - t361;
t296 = -qJD(5) - t396;
t375 = t232 * t296;
t376 = t330 * t296;
t398 = t232 * MDP(18) * t330 + (-t208 - t376) * MDP(21) + (-t232 ^ 2 + t330 ^ 2) * MDP(19) + (t207 - t375) * MDP(20);
t397 = MDP(10) * t318 + MDP(9) * t316;
t307 = t312 ^ 2;
t308 = t314 ^ 2;
t395 = (t307 + t308) * (qJ(2) * MDP(6) + MDP(5));
t291 = -pkin(2) * t314 - pkin(6) * t312 - pkin(1);
t281 = qJD(1) * t291 + qJD(2);
t277 = t318 * t281;
t377 = qJ(4) * t312;
t379 = qJ(2) * t316;
t324 = -t314 * t379 - t318 * t377;
t249 = qJD(1) * t324 + t277;
t240 = pkin(3) * t396 + t249;
t378 = qJ(2) * t318;
t350 = t314 * t378;
t371 = t316 * t281;
t250 = -qJ(4) * t348 + qJD(1) * t350 + t371;
t372 = t313 * t250;
t217 = t311 * t240 + t372;
t382 = pkin(7) * t268;
t206 = t217 - t382;
t205 = t206 * t354;
t282 = pkin(3) * t348 + qJ(2) * t362 + qJD(4);
t248 = pkin(4) * t268 + t282;
t394 = t248 * t232 + t205;
t384 = 0.2e1 * t307;
t391 = MDP(8) * (t316 ^ 2 - t318 ^ 2);
t389 = t384 + t308;
t387 = qJD(5) + t296;
t355 = qJD(4) * t312;
t320 = qJD(3) * t324 - t316 * t355;
t352 = qJD(1) * qJD(2);
t344 = t314 * t352;
t367 = t281 * t356 + t318 * t344;
t226 = qJD(1) * t320 + t367;
t360 = qJD(2) * t316;
t346 = t314 * t360;
t323 = -t318 * t355 - t346;
t357 = qJD(3) * t316;
t373 = t312 * t316;
t227 = -t281 * t357 + ((qJ(4) * t373 - t350) * qJD(3) + t323) * qJD(1);
t200 = -t226 * t311 + t313 * t227;
t198 = pkin(7) * t267 + t200;
t201 = t313 * t226 + t311 * t227;
t199 = pkin(7) * t266 + t201;
t341 = t317 * t198 - t315 * t199;
t386 = -t248 * t330 + t341;
t342 = -t291 + t377;
t385 = t316 * t342 - t350;
t383 = pkin(3) * t311;
t381 = pkin(7) * t272;
t319 = qJD(1) ^ 2;
t374 = t307 * t319;
t244 = t311 * t250;
t359 = qJD(2) * t318;
t366 = t291 * t356 + t314 * t359;
t238 = t320 + t366;
t239 = qJD(3) * t385 + t323;
t212 = t313 * t238 + t311 * t239;
t220 = t313 * t249 - t244;
t254 = -t342 * t318 + (-pkin(3) - t379) * t314;
t222 = t311 * t254 - t313 * t385;
t369 = t288 * qJD(3) - t314 * t325;
t329 = t311 * t316 - t313 * t318;
t368 = t396 * t329;
t280 = pkin(3) * t335 + t312 * t352;
t286 = pkin(3) * t345 + t312 * qJD(2);
t289 = pkin(3) * t373 + t312 * qJ(2);
t358 = qJD(3) * t312;
t353 = qJD(3) - t396;
t349 = qJ(2) * t357;
t211 = -t238 * t311 + t313 * t239;
t216 = t313 * t240 - t244;
t219 = -t249 * t311 - t372;
t221 = t313 * t254 + t311 * t385;
t339 = qJD(1) * t353;
t338 = pkin(3) * t347;
t337 = t318 * t307 * t316 * MDP(7);
t334 = qJD(5) * t288 + t369;
t333 = -qJD(5) * t329 - t368;
t204 = pkin(4) * t396 + t216 - t381;
t331 = -t315 * t204 - t317 * t206;
t279 = t329 * t312;
t241 = t278 * t317 - t279 * t315;
t242 = -t278 * t315 - t279 * t317;
t304 = pkin(3) * t313 + pkin(4);
t271 = t329 * t358;
t257 = pkin(4) * t272 + t338;
t255 = pkin(4) * t278 + t289;
t251 = -pkin(4) * t271 + t286;
t243 = -pkin(4) * t266 + t280;
t218 = -pkin(7) * t278 + t222;
t215 = -pkin(4) * t314 + pkin(7) * t279 + t221;
t214 = qJD(5) * t242 - t317 * t271 - t275 * t315;
t213 = -qJD(5) * t241 + t271 * t315 - t275 * t317;
t210 = t220 - t381;
t209 = t219 + t382;
t203 = pkin(7) * t271 + t212;
t202 = pkin(7) * t275 + t211;
t1 = [(-t396 * t346 + ((-t291 * t316 - t350) * t396 + t314 * t371) * qJD(3) + t389 * qJD(1) * (qJ(2) * t356 + t360)) * MDP(12) + (-(-t314 * t349 + t366) * t396 + t367 * t314 + (-t349 * t389 + t359 * t384) * qJD(1)) * MDP(13) + (-t200 * t314 + t211 * t396 - t266 * t289 + t268 * t286 - t271 * t282 + t278 * t280) * MDP(14) + (t201 * t314 - t212 * t396 - t267 * t289 + t272 * t286 - t275 * t282 - t279 * t280) * MDP(15) + (t200 * t279 - t201 * t278 - t211 * t272 - t212 * t268 + t216 * t275 + t217 * t271 + t221 * t267 + t222 * t266) * MDP(16) + (t200 * t221 + t201 * t222 + t211 * t216 + t212 * t217 + t280 * t289 + t282 * t286) * MDP(17) + (t207 * t242 + t213 * t330) * MDP(18) + (-t207 * t241 - t208 * t242 - t213 * t232 - t214 * t330) * MDP(19) + (-t207 * t314 - t213 * t296) * MDP(20) + (t208 * t314 + t214 * t296) * MDP(21) + (-(t202 * t317 - t203 * t315) * t296 - t341 * t314 + t251 * t232 + t255 * t208 + t243 * t241 + t248 * t214 + (-(-t215 * t315 - t218 * t317) * t296 - t331 * t314) * qJD(5)) * MDP(23) + (-t205 * t314 + t255 * t207 + t248 * t213 + t251 * t330 + t243 * t242 + ((-qJD(5) * t218 + t202) * t296 + t198 * t314) * t315 + ((qJD(5) * t215 + t203) * t296 + (qJD(5) * t204 + t199) * t314) * t317) * MDP(24) + (t384 * t391 - 0.2e1 * t337) * qJD(1) * qJD(3) + 0.2e1 * t352 * t395 + t397 * (-t396 + t361) * t358; (-t268 * t362 - t369 * t396) * MDP(14) + (-t272 * t362 + t368 * t396) * MDP(15) + (t266 * t288 - t267 * t329 + t268 * t368 + t272 * t369) * MDP(16) + (-t200 * t329 + t201 * t288 - t216 * t369 - t217 * t368 - t282 * t362) * MDP(17) + (-t232 * t362 + (t315 * t333 + t317 * t334) * t296) * MDP(23) + (-t330 * t362 + (-t315 * t334 + t317 * t333) * t296) * MDP(24) - t319 * t395 + (MDP(12) * t316 + MDP(13) * t318) * (-t396 ^ 2 - t374); t319 * t337 - t374 * t391 + ((-t281 * t353 - t344) * t316 + (-t314 * t339 - t374) * t378) * MDP(12) + (t277 * t396 + (t353 * t361 + t374) * t379 - t367) * MDP(13) + (-t219 * t396 - t268 * t338 - t272 * t282 + t200) * MDP(14) + (t220 * t396 + t268 * t282 - t272 * t338 - t201) * MDP(15) + ((t217 + t219) * t272 + (-t216 + t220) * t268 + (t266 * t311 + t267 * t313) * pkin(3)) * MDP(16) + (-t216 * t219 - t217 * t220 + (t200 * t313 + t201 * t311 - t282 * t347) * pkin(3)) * MDP(17) + ((t209 * t317 - t210 * t315) * t296 - t257 * t232 + (-(-t304 * t315 - t317 * t383) * t296 + t331) * qJD(5) + t386) * MDP(23) + (-t317 * t199 - t315 * t198 - (t209 * t315 + t210 * t317) * t296 - t257 * t330 + ((t304 * t317 - t315 * t383) * t296 - t317 * t204) * qJD(5) + t394) * MDP(24) - t397 * t312 * t339 + t398; (t272 * t396 - t266) * MDP(14) - t272 ^ 2 * MDP(16) + (t216 * t272 + t280) * MDP(17) + (t208 - t376) * MDP(23) + (t207 + t375) * MDP(24) - ((qJD(3) + t396) * MDP(15) - t217 * MDP(17) + MDP(16) * t268) * t268; (t331 * t387 + t386) * MDP(23) + ((t206 * t296 - t198) * t315 + (-t204 * t387 - t199) * t317 + t394) * MDP(24) + t398;];
tauc = t1;
