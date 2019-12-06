% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:48:17
% EndTime: 2019-12-05 18:48:22
% DurationCPUTime: 1.55s
% Computational Cost: add. (1917->219), mult. (3256->287), div. (0->0), fcn. (2038->6), ass. (0->136)
t307 = cos(qJ(3));
t300 = qJD(1) + qJD(2);
t305 = sin(qJ(2));
t379 = pkin(1) * qJD(1);
t345 = t305 * t379;
t391 = pkin(7) + pkin(8);
t335 = t391 * t300 + t345;
t247 = t335 * t307;
t303 = sin(qJ(4));
t240 = t303 * t247;
t304 = sin(qJ(3));
t246 = t335 * t304;
t243 = qJD(3) * pkin(3) - t246;
t306 = cos(qJ(4));
t331 = t306 * t243 - t240;
t272 = t303 * t307 + t304 * t306;
t255 = t272 * t300;
t374 = t255 * qJ(5);
t392 = t374 - t331;
t353 = qJD(4) * t306;
t356 = qJD(3) * t307;
t390 = -t306 * t356 - t307 * t353;
t389 = MDP(7) * t304;
t325 = qJD(3) * t335;
t308 = cos(qJ(2));
t378 = pkin(1) * qJD(2);
t344 = qJD(1) * t378;
t328 = t308 * t344;
t228 = -t304 * t325 + t307 * t328;
t388 = (qJD(4) * t243 + t228) * t306;
t387 = (t304 ^ 2 - t307 ^ 2) * MDP(8);
t342 = qJD(3) * t391;
t273 = t304 * t342;
t274 = t307 * t342;
t283 = t391 * t304;
t297 = t307 * pkin(8);
t284 = pkin(7) * t307 + t297;
t322 = t283 * t303 - t284 * t306;
t347 = t308 * t379;
t386 = qJD(4) * t322 + t272 * t347 + t303 * t273 - t306 * t274;
t369 = t306 * t307;
t371 = t303 * t304;
t271 = -t369 + t371;
t354 = qJD(4) * t303;
t385 = -t271 * t347 + t306 * t273 + t303 * t274 + t283 * t353 + t284 * t354;
t384 = t304 * MDP(12) + t307 * MDP(13);
t299 = qJD(3) + qJD(4);
t383 = t255 ^ 2;
t381 = pkin(1) * t308;
t291 = pkin(1) * t305 + pkin(7);
t380 = -pkin(8) - t291;
t343 = t300 * t371;
t253 = -t300 * t369 + t343;
t377 = qJ(5) * t253;
t376 = qJ(5) * t272;
t294 = -pkin(3) * t307 - pkin(2);
t257 = t294 * t300 - t347;
t222 = pkin(4) * t253 + qJD(5) + t257;
t375 = t222 * t255;
t373 = t257 * t255;
t372 = t300 * t304;
t309 = qJD(3) ^ 2;
t370 = t304 * t309;
t242 = t306 * t247;
t368 = t307 * t309;
t232 = t299 * t272;
t363 = -t232 * qJ(5) - t271 * qJD(5);
t367 = t363 - t385;
t326 = t299 * t371;
t231 = t326 + t390;
t321 = qJ(5) * t231 - qJD(5) * t272;
t366 = t321 + t386;
t207 = pkin(4) * t299 - t392;
t365 = t207 + t392;
t288 = t305 * t344;
t357 = qJD(3) * t304;
t341 = t300 * t357;
t259 = pkin(3) * t341 + t288;
t364 = t257 * t232 + t259 * t271;
t362 = -t257 * t231 + t259 * t272;
t361 = -t306 * t246 - t240;
t277 = -pkin(2) * t300 - t347;
t360 = t277 * t356 + t304 * t288;
t359 = t390 * t300;
t355 = qJD(3) * t308;
t351 = t307 * MDP(12);
t349 = -qJD(2) + t300;
t348 = pkin(3) * t372;
t346 = t308 * t378;
t295 = pkin(3) * t357;
t340 = t300 * t356;
t337 = -pkin(3) * t299 - t243;
t336 = pkin(4) * t232 + t295;
t334 = qJD(3) * t380;
t229 = -t304 * t328 - t307 * t325;
t333 = -t303 * t228 + t306 * t229;
t332 = t303 * t229 - t247 * t354;
t330 = t246 * t303 - t242;
t221 = t232 * t300;
t198 = -qJ(5) * t221 - qJD(5) * t253 + t332 + t388;
t220 = t300 * t326 + t359;
t324 = -t303 * t243 - t242;
t312 = qJD(4) * t324 + t333;
t199 = qJ(5) * t220 - qJD(5) * t255 + t312;
t209 = -t324 - t377;
t327 = -t198 * t271 - t199 * t272 + t207 * t231 - t209 * t232;
t215 = pkin(4) * t221 + t259;
t267 = t380 * t304;
t268 = t291 * t307 + t297;
t323 = -t267 * t303 - t268 * t306;
t320 = pkin(4) * t271 + t294;
t319 = t257 * t253 - t332;
t252 = t253 ^ 2;
t318 = t255 * t253 * MDP(14) + (-t359 + (t253 - t343) * t299) * MDP(16) + (-t252 + t383) * MDP(15);
t244 = t304 * t334 + t307 * t346;
t245 = -t304 * t346 + t307 * t334;
t317 = t306 * t244 + t303 * t245 + t267 * t353 - t268 * t354;
t315 = -t345 + t295;
t313 = -MDP(10) * t370 + (t220 * t271 - t221 * t272 + t231 * t253 - t232 * t255) * MDP(15) + (-t220 * t272 - t231 * t255) * MDP(14) - 0.2e1 * t300 * qJD(3) * t387 + 0.2e1 * t340 * t389 + MDP(9) * t368 + (-t231 * MDP(16) - t232 * MDP(17)) * t299;
t311 = qJD(4) * t323 - t244 * t303 + t306 * t245;
t296 = t305 * t378;
t293 = -pkin(2) - t381;
t292 = pkin(3) * t306 + pkin(4);
t280 = t294 - t381;
t275 = t296 + t295;
t266 = t271 * qJ(5);
t260 = t277 * t357;
t224 = -t266 - t322;
t223 = -t283 * t306 - t284 * t303 - t376;
t217 = -t266 - t323;
t216 = t267 * t306 - t268 * t303 - t376;
t211 = -t374 + t361;
t210 = t330 + t377;
t201 = t311 + t321;
t200 = t317 + t363;
t1 = [(((-qJD(1) - t300) * MDP(6) - t384 * qJD(3)) * t308 + (-qJD(1) * t351 + (t304 * MDP(13) - MDP(5) - t351) * t300) * t305) * t378 + (t198 * t217 + t209 * t200 + t199 * t216 + t207 * t201 + t215 * (t320 - t381) + t222 * (t296 + t336)) * MDP(22) + t313 + (t280 * t221 + t275 * t253 + t299 * t311 + t364) * MDP(19) + (-t291 * t368 + t293 * t341 + t260) * MDP(12) + (t291 * t370 + t293 * t340 + t360) * MDP(13) + (-t280 * t220 + t275 * t255 - t299 * t317 + t362) * MDP(20) + (-t200 * t253 - t201 * t255 + t216 * t220 - t217 * t221 + t327) * MDP(21) - t288 * MDP(5); t349 * MDP(6) * t347 + (t300 * t345 - t288) * MDP(5) + t313 + (t294 * t221 + t315 * t253 + t386 * t299 + t364) * MDP(19) + (t198 * t224 + t199 * t223 + t215 * t320 + (t336 - t345) * t222 + t367 * t209 + t366 * t207) * MDP(22) + (-pkin(2) * t341 - pkin(7) * t368 + t260 + (t305 * t307 * t349 + t304 * t355) * t379) * MDP(12) + (-pkin(2) * t340 + pkin(7) * t370 + (-t305 * t372 + t307 * t355) * t379 + t360) * MDP(13) + (-t294 * t220 + t315 * t255 + t385 * t299 + t362) * MDP(20) + (t220 * t223 - t221 * t224 - t253 * t367 - t255 * t366 + t327) * MDP(21); (-t253 * t348 - t373 - t330 * t299 + (t303 * t337 - t242) * qJD(4) + t333) * MDP(19) + (-t255 * t348 + t361 * t299 + (qJD(4) * t337 - t228) * t306 + t319) * MDP(20) + (t220 * t292 + (t209 + t210) * t255 + (-t207 + t211) * t253 + (-t221 * t303 + (-t253 * t306 + t255 * t303) * qJD(4)) * pkin(3)) * MDP(21) + (-pkin(4) * t375 + t199 * t292 - t207 * t210 - t209 * t211 + (-t222 * t372 + t198 * t303 + (-t207 * t303 + t209 * t306) * qJD(4)) * pkin(3)) * MDP(22) + t318 + t384 * (-t277 * t300 - t328) + (-t307 * t389 + t387) * t300 ^ 2; (-t299 * t324 + t312 - t373) * MDP(19) + (t299 * t331 + t319 - t388) * MDP(20) + (pkin(4) * t220 - t253 * t365) * MDP(21) + (t365 * t209 + (t199 - t375) * pkin(4)) * MDP(22) + t318; (-t252 - t383) * MDP(21) + (t207 * t255 + t209 * t253 + t215) * MDP(22);];
tauc = t1;
