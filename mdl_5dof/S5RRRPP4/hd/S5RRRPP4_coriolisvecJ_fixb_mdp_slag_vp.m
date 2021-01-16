% Calculate Coriolis joint torque vector for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:24:55
% EndTime: 2021-01-15 22:25:01
% DurationCPUTime: 2.23s
% Computational Cost: add. (3319->289), mult. (8747->363), div. (0->0), fcn. (6035->6), ass. (0->147)
t320 = cos(qJ(3));
t321 = cos(qJ(2));
t363 = qJD(1) * t321;
t355 = t320 * t363;
t318 = sin(qJ(3));
t319 = sin(qJ(2));
t364 = qJD(1) * t319;
t356 = t318 * t364;
t276 = -t355 + t356;
t278 = -t318 * t363 - t320 * t364;
t316 = sin(pkin(8));
t317 = cos(pkin(8));
t248 = -t276 * t317 + t278 * t316;
t313 = qJD(2) + qJD(3);
t394 = t248 * t313;
t340 = -t276 * t316 - t278 * t317;
t393 = t340 ^ 2;
t359 = qJD(1) * qJD(2);
t392 = -0.2e1 * t359;
t391 = MDP(4) * t319;
t390 = MDP(5) * (t319 ^ 2 - t321 ^ 2);
t309 = -pkin(2) * t321 - pkin(1);
t294 = t309 * qJD(1);
t259 = pkin(3) * t276 + qJD(4) + t294;
t218 = -pkin(4) * t248 - qJ(5) * t340 + t259;
t378 = t218 * t340;
t387 = pkin(6) + pkin(7);
t295 = t387 * t319;
t289 = qJD(1) * t295;
t385 = qJD(2) * pkin(2);
t283 = -t289 + t385;
t357 = qJD(2) * t387;
t346 = qJD(1) * t357;
t284 = t319 * t346;
t389 = (qJD(3) * t283 - t284) * t320;
t273 = t278 * qJ(4);
t296 = t387 * t321;
t291 = qJD(1) * t296;
t279 = t318 * t291;
t366 = -t289 * t320 - t279;
t241 = t273 + t366;
t372 = t320 * t291;
t338 = t289 * t318 - t372;
t384 = qJ(4) * t276;
t330 = t338 + t384;
t358 = pkin(2) * t316 * t318;
t361 = qJD(3) * t320;
t367 = qJD(3) * t358 + t316 * t330 + (-pkin(2) * t361 + t241) * t317;
t350 = t283 * t320 - t279;
t238 = t273 + t350;
t288 = t318 * t321 + t319 * t320;
t388 = qJD(1) * t288;
t386 = pkin(3) * t278;
t258 = t313 * t288;
t253 = t258 * qJD(1);
t285 = t321 * t346;
t362 = qJD(3) * t318;
t348 = -t318 * t285 - t291 * t362;
t209 = -qJ(4) * t253 - qJD(4) * t276 + t348 + t389;
t354 = t321 * t359;
t252 = qJD(3) * t355 - t313 * t356 + t320 * t354;
t339 = -t283 * t318 - t372;
t349 = t318 * t284 - t285 * t320;
t326 = qJD(3) * t339 + t349;
t210 = -qJ(4) * t252 + qJD(4) * t278 + t326;
t198 = t209 * t316 - t210 * t317;
t287 = t318 * t319 - t320 * t321;
t337 = t295 * t318 - t296 * t320;
t244 = -qJ(4) * t287 - t337;
t329 = -qJ(4) * t288 - t295 * t320 - t296 * t318;
t225 = t244 * t316 - t317 * t329;
t383 = t198 * t225;
t290 = t319 * t357;
t292 = t321 * t357;
t332 = -t290 * t320 - t292 * t318 - t295 * t361 - t296 * t362;
t219 = -qJ(4) * t258 - qJD(4) * t287 + t332;
t257 = t313 * t287;
t325 = qJD(3) * t337 + t290 * t318 - t292 * t320;
t324 = qJ(4) * t257 - qJD(4) * t288 + t325;
t202 = t219 * t316 - t317 * t324;
t382 = t202 * t313;
t203 = t317 * t219 + t316 * t324;
t381 = t203 * t313;
t239 = -t339 - t384;
t235 = t317 * t239;
t215 = t238 * t316 + t235;
t380 = t215 * t340;
t377 = t239 * t316;
t216 = t238 * t317 - t377;
t379 = t216 * t313;
t376 = t340 * t259;
t375 = t294 * t278;
t374 = t317 * t318;
t322 = qJD(2) ^ 2;
t373 = t319 * t322;
t371 = t321 * t322;
t323 = qJD(1) ^ 2;
t370 = t321 * t323;
t199 = t209 * t317 + t210 * t316;
t369 = t241 * t316 - t317 * t330 - (t316 * t320 + t374) * qJD(3) * pkin(2);
t368 = -qJD(5) + t367;
t234 = pkin(3) * t313 + t238;
t214 = t234 * t316 + t235;
t308 = pkin(2) * t320 + pkin(3);
t272 = pkin(2) * t374 + t308 * t316;
t360 = qJD(5) - t216;
t311 = t319 * t385;
t310 = pkin(2) * t364;
t353 = -pkin(2) * t313 - t283;
t242 = pkin(3) * t253 + qJD(2) * t310;
t251 = pkin(3) * t258 + t311;
t352 = t369 * t340;
t351 = pkin(1) * t392;
t227 = t252 * t316 + t253 * t317;
t345 = t218 * t248 + t199;
t344 = -t248 * t259 - t199;
t213 = t234 * t317 - t377;
t211 = -pkin(4) * t313 + qJD(5) - t213;
t212 = qJ(5) * t313 + t214;
t343 = -t211 * t248 + t212 * t340;
t342 = t213 * t248 + t214 * t340;
t228 = t252 * t317 - t253 * t316;
t261 = pkin(3) * t287 + t309;
t271 = t308 * t317 - t358;
t336 = t215 * t313 - t198;
t334 = t294 * t276 - t348;
t333 = -t278 * t276 * MDP(11) + t252 * MDP(13) + (-t276 ^ 2 + t278 ^ 2) * MDP(12) + (t276 * MDP(13) + (-t278 - t388) * MDP(14)) * t313;
t223 = pkin(4) * t340 - qJ(5) * t248 - t386;
t331 = t313 * t369 - t198;
t226 = t244 * t317 + t316 * t329;
t256 = -t287 * t316 + t288 * t317;
t327 = t198 * t256 + t202 * t340 + t203 * t248 + t225 * t228 - t226 * t227;
t201 = pkin(4) * t227 - qJ(5) * t228 - qJD(5) * t340 + t242;
t312 = t313 * qJD(5);
t306 = -pkin(3) * t317 - pkin(4);
t305 = pkin(3) * t316 + qJ(5);
t266 = -pkin(4) - t271;
t265 = qJ(5) + t272;
t260 = t310 - t386;
t255 = t287 * t317 + t288 * t316;
t230 = -t257 * t317 - t258 * t316;
t229 = -t257 * t316 + t258 * t317;
t224 = pkin(4) * t255 - qJ(5) * t256 + t261;
t222 = t223 + t310;
t204 = pkin(4) * t229 - qJ(5) * t230 - qJD(5) * t256 + t251;
t197 = t312 + t199;
t1 = [0.2e1 * t354 * t391 + t390 * t392 + MDP(6) * t371 - MDP(7) * t373 + (-pkin(6) * t371 + t319 * t351) * MDP(9) + (pkin(6) * t373 + t321 * t351) * MDP(10) + (t252 * t288 + t257 * t278) * MDP(11) + (-t252 * t287 - t253 * t288 + t257 * t276 + t258 * t278) * MDP(12) + (t309 * t253 + t294 * t258 + (qJD(1) * t287 + t276) * t311) * MDP(16) + (t309 * t252 - t294 * t257 + (-t278 + t388) * t311) * MDP(17) + (t227 * t261 + t229 * t259 + t242 * t255 - t248 * t251 - t382) * MDP(18) + (t228 * t261 + t230 * t259 + t242 * t256 + t251 * t340 - t381) * MDP(19) + (-t199 * t255 - t213 * t230 - t214 * t229 + t327) * MDP(20) + (t199 * t226 - t202 * t213 + t203 * t214 + t242 * t261 + t251 * t259 + t383) * MDP(21) + (t201 * t255 - t204 * t248 + t218 * t229 + t224 * t227 - t382) * MDP(22) + (-t197 * t255 + t211 * t230 - t212 * t229 + t327) * MDP(23) + (-t201 * t256 - t204 * t340 - t218 * t230 - t224 * t228 + t381) * MDP(24) + (t197 * t226 + t201 * t224 + t202 * t211 + t203 * t212 + t204 * t218 + t383) * MDP(25) + (-t257 * MDP(13) - t258 * MDP(14) + MDP(16) * t325 - MDP(17) * t332) * t313; -t370 * t391 + t323 * t390 + (-t276 * t310 + t375 - t338 * t313 + (t318 * t353 - t372) * qJD(3) + t349) * MDP(16) + (t278 * t310 + t366 * t313 + (qJD(3) * t353 + t284) * t320 + t334) * MDP(17) + (t248 * t260 + t331 - t376) * MDP(18) + (-t260 * t340 + t313 * t367 + t344) * MDP(19) + (-t227 * t272 - t228 * t271 - t248 * t367 + t342 - t352) * MDP(20) + (-t198 * t271 + t199 * t272 + t213 * t369 - t214 * t367 - t259 * t260) * MDP(21) + (t222 * t248 + t331 - t378) * MDP(22) + (-t227 * t265 + t228 * t266 - t248 * t368 + t343 - t352) * MDP(23) + (t222 * t340 - t313 * t368 + t312 + t345) * MDP(24) + (t197 * t265 + t198 * t266 - t211 * t369 - t212 * t368 - t218 * t222) * MDP(25) + t333 + (MDP(9) * t319 * t323 + MDP(10) * t370) * pkin(1); (-t313 * t339 + t326 + t375) * MDP(16) + (t313 * t350 + t334 - t389) * MDP(17) + (-t248 * t386 + t336 - t376) * MDP(18) + (t340 * t386 + t344 + t379) * MDP(19) + (-t380 - t216 * t248 + (-t227 * t316 - t228 * t317) * pkin(3) + t342) * MDP(20) + (t213 * t215 - t214 * t216 + (-t198 * t317 + t199 * t316 + t259 * t278) * pkin(3)) * MDP(21) + (t223 * t248 + t336 - t378) * MDP(22) + (-t227 * t305 + t228 * t306 + t248 * t360 + t343 - t380) * MDP(23) + (t223 * t340 + 0.2e1 * t312 + t345 - t379) * MDP(24) + (t197 * t305 + t198 * t306 - t211 * t215 + t212 * t360 - t218 * t223) * MDP(25) + t333; (t213 * t340 - t214 * t248 + t242) * MDP(21) + (-t211 * t340 - t212 * t248 + t201) * MDP(25) + (MDP(20) + MDP(23)) * (-t248 ^ 2 - t393) + (MDP(18) + MDP(22)) * (t313 * t340 + t227) + (MDP(19) - MDP(24)) * (t228 + t394); -t340 * t248 * MDP(22) + (t228 - t394) * MDP(23) + (-t313 ^ 2 - t393) * MDP(24) + (-t212 * t313 + t198 + t378) * MDP(25);];
tauc = t1;
