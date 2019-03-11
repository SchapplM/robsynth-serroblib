% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:07:51
% EndTime: 2019-03-08 20:07:58
% DurationCPUTime: 3.16s
% Computational Cost: add. (2921->322), mult. (7511->446), div. (0->0), fcn. (5886->10), ass. (0->142)
t317 = sin(pkin(11));
t393 = pkin(8) + qJ(3);
t302 = t393 * t317;
t319 = cos(pkin(11));
t303 = t393 * t319;
t322 = sin(qJ(4));
t325 = cos(qJ(4));
t268 = t302 * t325 + t322 * t303;
t298 = t317 * t322 - t325 * t319;
t242 = -t298 * qJD(3) - t268 * qJD(4);
t318 = sin(pkin(6));
t326 = cos(qJ(2));
t379 = t318 * t326;
t329 = t298 * t379;
t271 = qJD(1) * t329;
t407 = t242 + t271;
t293 = t298 * qJD(4);
t299 = t317 * t325 + t319 * t322;
t294 = t299 * qJD(4);
t323 = sin(qJ(2));
t363 = qJD(1) * t318;
t351 = t323 * t363;
t406 = pkin(4) * t294 + pkin(9) * t293 - t351;
t301 = qJD(2) * qJ(3) + t351;
t320 = cos(pkin(6));
t362 = qJD(1) * t320;
t309 = t319 * t362;
t391 = pkin(8) * qJD(2);
t265 = t309 + (-t301 - t391) * t317;
t278 = t319 * t301 + t317 * t362;
t266 = t319 * t391 + t278;
t229 = t322 * t265 + t325 * t266;
t405 = qJD(4) * t229;
t292 = qJD(2) * t299;
t324 = cos(qJ(5));
t359 = qJD(2) * t325;
t310 = t319 * t359;
t307 = qJD(4) * t310;
t361 = qJD(2) * t322;
t349 = t317 * t361;
t337 = qJD(4) * t349 - t307;
t356 = t324 * qJD(4);
t321 = sin(qJ(5));
t358 = qJD(5) * t321;
t239 = -qJD(5) * t356 + t292 * t358 + t324 * t337;
t275 = qJD(4) * t321 + t292 * t324;
t283 = qJD(2) * t294;
t225 = qJD(4) * pkin(9) + t229;
t312 = -pkin(3) * t319 - pkin(2);
t350 = t326 * t363;
t342 = qJD(3) - t350;
t284 = t312 * qJD(2) + t342;
t290 = -t310 + t349;
t232 = pkin(4) * t290 - pkin(9) * t292 + t284;
t215 = t225 * t324 + t232 * t321;
t396 = t265 * t325 - t322 * t266;
t297 = (qJD(3) + t350) * qJD(2);
t399 = t298 * t297;
t218 = qJD(4) * t396 - t399;
t360 = qJD(2) * t323;
t348 = t318 * t360;
t306 = qJD(1) * t348;
t238 = t283 * pkin(4) + t337 * pkin(9) + t306;
t234 = t324 * t238;
t328 = -t215 * qJD(5) - t218 * t321 + t234;
t203 = pkin(5) * t283 + qJ(6) * t239 - qJD(6) * t275 + t328;
t273 = t292 * t321 - t356;
t210 = -qJ(6) * t273 + t215;
t285 = qJD(5) + t290;
t404 = t285 * t210 + t203;
t240 = t275 * qJD(5) - t321 * t337;
t357 = qJD(5) * t324;
t332 = t324 * t218 - t225 * t358 + t232 * t357 + t321 * t238;
t204 = -qJ(6) * t240 - qJD(6) * t273 + t332;
t214 = -t225 * t321 + t324 * t232;
t209 = -qJ(6) * t275 + t214;
t207 = pkin(5) * t285 + t209;
t403 = -t285 * t207 + t204;
t392 = -qJ(6) - pkin(9);
t402 = -qJ(6) * t290 + qJD(5) * t392;
t364 = t317 ^ 2 + t319 ^ 2;
t398 = t364 * MDP(7);
t397 = -t271 * t321 + t324 * t406;
t258 = pkin(4) * t298 - pkin(9) * t299 + t312;
t395 = t258 * t357 + t406 * t321 + t324 * t407;
t394 = t275 ^ 2;
t390 = qJD(2) * pkin(2);
t388 = t239 * t321;
t386 = t273 * t290;
t385 = t275 * t285;
t384 = t299 * t321;
t383 = t299 * t324;
t380 = t318 * t323;
t376 = t321 * t283;
t375 = t321 * t285;
t269 = -t302 * t322 + t303 * t325;
t264 = t324 * t269;
t280 = t324 * t283;
t338 = qJ(6) * t293 - qJD(6) * t299;
t374 = pkin(5) * t294 - t242 * t321 + t338 * t324 + (-t264 + (qJ(6) * t299 - t258) * t321) * qJD(5) + t397;
t347 = t299 * t357;
t373 = -qJ(6) * t347 + (-qJD(5) * t269 + t338) * t321 + t395;
t372 = t207 - t209;
t256 = pkin(4) * t292 + pkin(9) * t290;
t371 = t321 * t256 + t324 * t396;
t370 = -t321 * t240 - t273 * t357;
t330 = t299 * t379;
t369 = -qJD(1) * t330 + t299 * qJD(3) + t269 * qJD(4);
t368 = -t290 * t375 + t280;
t367 = t321 * t258 + t264;
t366 = qJD(6) * t324 + t402 * t321 - t371;
t251 = t324 * t256;
t365 = -pkin(5) * t292 - t251 + t402 * t324 + (-qJD(6) + t396) * t321;
t345 = t364 * t297;
t344 = t285 * t324;
t219 = t299 * t297 + t405;
t340 = (-t301 * t317 + t309) * t317 - t278 * t319;
t288 = -t317 * t380 + t319 * t320;
t289 = t317 * t320 + t319 * t380;
t339 = t288 * t325 - t289 * t322;
t247 = t288 * t322 + t289 * t325;
t224 = -qJD(4) * pkin(4) - t396;
t235 = -t247 * t321 - t324 * t379;
t336 = -t247 * t324 + t321 * t379;
t335 = -t293 * t321 + t347;
t334 = -t293 * t324 - t299 * t358;
t208 = pkin(5) * t240 + t219;
t333 = -pkin(9) * t283 + t285 * t224;
t327 = qJD(2) ^ 2;
t305 = t392 * t324;
t304 = t392 * t321;
t300 = t342 - t390;
t272 = t273 ^ 2;
t254 = t324 * t258;
t227 = qJD(2) * t330 + t247 * qJD(4);
t226 = -qJD(2) * t329 + t339 * qJD(4);
t222 = -qJ(6) * t384 + t367;
t221 = pkin(5) * t273 + qJD(6) + t224;
t220 = pkin(5) * t298 - qJ(6) * t383 - t269 * t321 + t254;
t213 = t336 * qJD(5) - t226 * t321 + t324 * t348;
t212 = t235 * qJD(5) + t226 * t324 + t321 * t348;
t1 = [((-t288 * t317 + t289 * t319) * t297 + (t300 * t323 + (-t340 - t351) * t326) * t318 * qJD(2)) * MDP(8) + (-qJD(4) * t227 + (-t283 * t326 + t290 * t360) * t318) * MDP(14) + (-t226 * qJD(4) + (t292 * t360 + t326 * t337) * t318) * MDP(15) + (t213 * t285 + t227 * t273 + t235 * t283 - t240 * t339) * MDP(21) + (-t212 * t285 + t227 * t275 + t239 * t339 + t283 * t336) * MDP(22) + (-t212 * t273 - t213 * t275 + t235 * t239 + t240 * t336) * MDP(23) + (t203 * t235 - t204 * t336 + t207 * t213 - t208 * t339 + t210 * t212 + t221 * t227) * MDP(24) + ((-MDP(4) + t398) * t326 + (-t319 * MDP(5) + t317 * MDP(6) - MDP(3)) * t323) * t318 * t327; (t342 * qJD(2) * t364 + t345) * MDP(7) + (-t340 * qJD(3) + qJ(3) * t345 + (t340 * t326 + (-t300 - t390) * t323) * t363) * MDP(8) + (-t292 * t293 - t337 * t299) * MDP(9) + (-t299 * t283 + t293 * t290 - t292 * t294 + t337 * t298) * MDP(10) + (t283 * t312 + t284 * t294 + (qJD(2) * t298 - t290) * t351) * MDP(14) + (-t284 * t293 + t312 * t307) * MDP(15) + (-t239 * t383 + t334 * t275) * MDP(16) + (-(-t273 * t324 - t275 * t321) * t293 + (t388 - t240 * t324 + (t273 * t321 - t275 * t324) * qJD(5)) * t299) * MDP(17) + (-t239 * t298 + t275 * t294 + t299 * t280 + t334 * t285) * MDP(18) + (-t240 * t298 - t273 * t294 - t335 * t285 - t299 * t376) * MDP(19) + (t283 * t298 + t285 * t294) * MDP(20) + (t254 * t283 + (-t225 * t357 + t234) * t298 + t214 * t294 + t268 * t240 + t224 * t347 + (-t269 * t357 + t397) * t285 + t369 * t273 + ((-qJD(5) * t258 - t242) * t285 - t269 * t283 + (-qJD(5) * t232 - t218) * t298 + t219 * t299 - t224 * t293) * t321) * MDP(21) + (-t367 * t283 - t332 * t298 - t215 * t294 - t268 * t239 + t219 * t383 + (t269 * t358 - t395) * t285 + t369 * t275 + t334 * t224) * MDP(22) + (t220 * t239 - t222 * t240 - (-t207 * t324 - t210 * t321) * t293 - t374 * t275 - t373 * t273 + (-t203 * t324 - t204 * t321 + (t207 * t321 - t210 * t324) * qJD(5)) * t299) * MDP(23) + (t204 * t222 + t203 * t220 + t208 * (pkin(5) * t384 + t268) + (t335 * pkin(5) + t369) * t221 + t373 * t210 + t374 * t207) * MDP(24) + (-t293 * MDP(11) - t294 * MDP(12) - t369 * MDP(14) + (-t312 * t349 - t407) * MDP(15)) * qJD(4); (t340 * qJD(2) + t306) * MDP(8) + t307 * MDP(15) + t368 * MDP(21) + t370 * MDP(23) - t327 * t398 + (-t273 * MDP(21) - t275 * MDP(22) - t221 * MDP(24)) * t292 + (-qJD(5) * t285 * MDP(21) - t283 * MDP(22) + MDP(23) * t385 + t403 * MDP(24)) * t321 + ((t239 - t386) * MDP(23) + t404 * MDP(24) - t285 ^ 2 * MDP(22)) * t324 + ((t317 * t359 + t319 * t361 + t292) * MDP(14) + (-t290 - t349) * MDP(15)) * qJD(4); -t290 ^ 2 * MDP(10) + (t307 + (t290 - t349) * qJD(4)) * MDP(11) + (-t219 + t405) * MDP(14) + (t284 * t290 + t399) * MDP(15) + (t275 * t344 - t388) * MDP(16) + ((-t239 - t386) * t324 - t275 * t375 + t370) * MDP(17) + (t285 * t344 + t376) * MDP(18) + (-t285 * t358 + t368) * MDP(19) + (-pkin(4) * t240 - t219 * t324 - t229 * t273 + (-pkin(9) * t357 - t251) * t285 + (t285 * t396 + t333) * t321) * MDP(21) + (pkin(4) * t239 + t219 * t321 - t229 * t275 + (pkin(9) * t358 + t371) * t285 + t333 * t324) * MDP(22) + (t239 * t304 + t240 * t305 - t366 * t273 - t365 * t275 - t404 * t321 + t403 * t324) * MDP(23) + (-t204 * t305 + t203 * t304 + t208 * (-pkin(5) * t324 - pkin(4)) + (pkin(5) * t375 - t229) * t221 + t366 * t210 + t365 * t207) * MDP(24) + (MDP(10) * t292 - t284 * MDP(14) - t275 * MDP(18) + t273 * MDP(19) - t285 * MDP(20) - t214 * MDP(21) + t215 * MDP(22) + t290 * MDP(9)) * t292; t275 * t273 * MDP(16) + (-t272 + t394) * MDP(17) + (t273 * t285 - t239) * MDP(18) + (-t240 + t385) * MDP(19) + t283 * MDP(20) + (t215 * t285 - t224 * t275 + t328) * MDP(21) + (t214 * t285 + t224 * t273 - t332) * MDP(22) + (pkin(5) * t239 - t372 * t273) * MDP(23) + (t372 * t210 + (-t221 * t275 + t203) * pkin(5)) * MDP(24); (-t272 - t394) * MDP(23) + (t207 * t275 + t210 * t273 + t208) * MDP(24);];
tauc  = t1;
