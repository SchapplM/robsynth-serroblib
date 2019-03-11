% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:19:51
% EndTime: 2019-03-09 03:19:57
% DurationCPUTime: 3.48s
% Computational Cost: add. (3105->348), mult. (7796->444), div. (0->0), fcn. (5616->6), ass. (0->147)
t333 = cos(pkin(9));
t404 = cos(qJ(3));
t372 = t404 * t333;
t357 = qJD(1) * t372;
t332 = sin(pkin(9));
t335 = sin(qJ(3));
t395 = t332 * t335;
t371 = qJD(1) * t395;
t299 = -t357 + t371;
t334 = sin(qJ(5));
t336 = cos(qJ(5));
t279 = qJD(3) * t336 + t299 * t334;
t310 = t332 * t404 + t335 * t333;
t413 = t310 * qJD(1);
t418 = qJD(5) + t413;
t414 = t336 * t418;
t423 = t279 * t414;
t327 = -pkin(2) * t333 - pkin(1);
t313 = qJD(1) * t327 + qJD(2);
t341 = -qJ(4) * t413 + t313;
t405 = pkin(3) + pkin(8);
t237 = t299 * t405 + t341;
t403 = pkin(7) + qJ(2);
t316 = t403 * t332;
t311 = qJD(1) * t316;
t317 = t403 * t333;
t312 = qJD(1) * t317;
t271 = t404 * t311 + t335 * t312;
t412 = qJD(4) + t271;
t351 = pkin(4) * t413 + t412;
t243 = -qJD(3) * t405 + t351;
t223 = t237 * t336 + t243 * t334;
t380 = qJD(5) * t334;
t304 = t310 * qJD(3);
t339 = qJD(1) * t304;
t379 = qJD(5) * t336;
t388 = t299 * t379 + t334 * t339;
t257 = qJD(3) * t380 - t388;
t322 = qJD(3) * t357;
t382 = qJD(3) * t335;
t370 = t332 * t382;
t286 = qJD(1) * t370 - t322;
t340 = t413 * qJD(3);
t242 = pkin(3) * t340 + t286 * qJ(4) - qJD(4) * t413;
t228 = pkin(8) * t339 + t242;
t369 = qJD(2) * t404;
t356 = qJD(1) * t369;
t374 = qJD(1) * qJD(2);
t367 = t335 * t374;
t368 = qJD(3) * t404;
t247 = -t311 * t382 + t312 * t368 + t332 * t356 + t333 * t367;
t232 = -pkin(4) * t286 + t247;
t364 = -t228 * t334 + t336 * t232;
t212 = -pkin(5) * t286 + qJ(6) * t257 - qJD(5) * t223 - qJD(6) * t279 + t364;
t277 = qJD(3) * t334 - t336 * t299;
t218 = -qJ(6) * t277 + t223;
t422 = t218 * t418 + t212;
t283 = t336 * t339;
t381 = qJD(5) * t279;
t258 = -t283 + t381;
t373 = -t336 * t228 - t334 * t232 - t243 * t379;
t347 = -t237 * t380 - t373;
t213 = -qJ(6) * t258 - qJD(6) * t277 + t347;
t222 = -t237 * t334 + t336 * t243;
t217 = -qJ(6) * t279 + t222;
t216 = pkin(5) * t418 + t217;
t421 = -t216 * t418 + t213;
t359 = t418 * t277;
t420 = t257 - t359;
t282 = t336 * t286;
t362 = t334 * t418;
t419 = -t362 * t418 - t282;
t309 = -t372 + t395;
t350 = -qJ(4) * t310 + t327;
t254 = t309 * t405 + t350;
t274 = t316 * t404 + t335 * t317;
t265 = t310 * pkin(4) + t274;
t389 = t336 * t254 + t334 * t265;
t275 = t335 * t316 - t404 * t317;
t411 = (t332 ^ 2 + t333 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t272 = -t335 * t311 + t404 * t312;
t256 = -pkin(4) * t299 + t272;
t331 = qJD(3) * qJ(4);
t248 = t256 + t331;
t409 = t248 * t418 + t286 * t405;
t261 = (qJD(2) * t332 + qJD(3) * t317) * t335 + t316 * t368 - t333 * t369;
t408 = t279 ^ 2;
t407 = t299 ^ 2;
t406 = t413 ^ 2;
t401 = qJ(4) * t299;
t400 = qJ(6) * t336;
t399 = t257 * t336;
t398 = t277 * t299;
t397 = t279 * t299;
t396 = t309 * t334;
t394 = t334 * t286;
t392 = qJ(6) + t405;
t391 = t216 - t217;
t251 = t405 * t413 + t401;
t390 = t336 * t251 + t334 * t256;
t250 = t336 * t256;
t377 = qJD(6) * t336;
t387 = t392 * t380 - t377 + pkin(5) * t299 - t250 - (-qJ(6) * t413 - t251) * t334;
t315 = t392 * t336;
t386 = -qJD(5) * t315 - qJD(6) * t334 - t400 * t413 - t390;
t383 = qJD(3) * t261;
t378 = qJD(5) * t405;
t262 = t310 * qJD(2) - t275 * qJD(3);
t376 = t262 * qJD(3);
t365 = -qJ(6) * t309 - t254;
t363 = t418 ^ 2;
t358 = -t311 * t368 - t312 * t382 - t332 * t367 + t333 * t356;
t330 = qJD(3) * qJD(4);
t244 = -t330 - t358;
t303 = -t333 * t368 + t370;
t354 = qJ(4) * t303 - qJD(4) * t310;
t349 = t334 * t304 + t309 * t379;
t348 = -t304 * t336 + t309 * t380;
t234 = t304 * t405 + t354;
t239 = -t303 * pkin(4) + t262;
t346 = t336 * t234 + t334 * t239 - t254 * t380 + t265 * t379;
t345 = -qJD(3) * t271 - t358;
t344 = qJD(3) * t272 - t247;
t343 = -t414 * t418 + t394;
t238 = -pkin(4) * t304 - t261;
t229 = -pkin(4) * t339 - t244;
t224 = t258 * pkin(5) + t229;
t314 = t392 * t334;
t289 = qJD(3) * t299;
t276 = t277 ^ 2;
t273 = t286 * t310;
t270 = pkin(3) * t309 + t350;
t269 = pkin(3) * t413 + t401;
t268 = -t331 - t272;
t267 = -qJD(3) * pkin(3) + t412;
t266 = -pkin(4) * t309 - t275;
t264 = t336 * t265;
t260 = pkin(3) * t299 + t341;
t253 = pkin(3) * t304 + t354;
t252 = t336 * t258;
t236 = t336 * t239;
t227 = pkin(5) * t277 + qJD(6) + t248;
t225 = t309 * t400 + t389;
t220 = pkin(5) * t310 + t334 * t365 + t264;
t215 = -qJ(6) * t348 + t309 * t377 + t346;
t214 = -pkin(5) * t303 + t236 + t365 * t379 + (-qJ(6) * t304 - qJD(5) * t265 - qJD(6) * t309 - t234) * t334;
t1 = [(-t303 * t413 - t273) * MDP(8) + (t286 * t309 + t303 * t299 - t304 * t413 - t310 * t340) * MDP(9) + (t313 * t304 + t327 * t339 - t376) * MDP(13) + (-t286 * t327 - t303 * t313 + t383) * MDP(14) + (t244 * t309 + t247 * t310 + t261 * t299 + t262 * t413 - t267 * t303 + t268 * t304 - t274 * t286 + t275 * t339) * MDP(15) + (-t242 * t309 - t253 * t299 - t260 * t304 - t270 * t340 + t376) * MDP(16) + (-t242 * t310 - t253 * t413 + t260 * t303 + t270 * t286 - t383) * MDP(17) + (t242 * t270 + t244 * t275 + t247 * t274 + t253 * t260 + t261 * t268 + t262 * t267) * MDP(18) + (-t257 * t396 + t279 * t349) * MDP(19) + ((-t277 * t334 + t279 * t336) * t304 + (-t399 - t258 * t334 + (-t277 * t336 - t279 * t334) * qJD(5)) * t309) * MDP(20) + (-t257 * t310 - t279 * t303 - t309 * t394 + t349 * t418) * MDP(21) + (-t258 * t310 + t277 * t303 - t282 * t309 - t348 * t418) * MDP(22) + (-t303 * t418 - t273) * MDP(23) + ((-t234 * t334 + t236) * t418 - (-t254 * t334 + t264) * t286 + t364 * t310 - t222 * t303 + t238 * t277 + t266 * t258 + (-t229 * t309 - t248 * t304) * t336 + (-t223 * t310 + t248 * t396 - t389 * t418) * qJD(5)) * MDP(24) + (t223 * t303 + t229 * t396 + t238 * t279 + t248 * t349 - t266 * t257 + t286 * t389 - t310 * t347 - t346 * t418) * MDP(25) + (-t214 * t279 - t215 * t277 + t220 * t257 - t225 * t258 + (-t216 * t334 + t218 * t336) * t304 + (-t212 * t334 + t213 * t336 + (-t216 * t336 - t218 * t334) * qJD(5)) * t309) * MDP(26) + (t213 * t225 + t218 * t215 + t212 * t220 + t216 * t214 + t224 * ((-pkin(5) * t336 - pkin(4)) * t309 - t275) + t227 * (pkin(5) * t348 + t238)) * MDP(27) + 0.2e1 * t374 * t411 + (-MDP(10) * t303 - MDP(11) * t304) * qJD(3); t322 * MDP(14) + (-t406 - t407) * MDP(15) + (t286 + t289) * MDP(17) + (-t267 * t413 - t268 * t299 + t242) * MDP(18) + (t343 + t398) * MDP(24) + (t397 - t419) * MDP(25) + (-t334 * t420 - t252 + t423) * MDP(26) + (t227 * t299 - t334 * t422 + t336 * t421) * MDP(27) - qJD(1) ^ 2 * t411 + ((-t299 - t371) * MDP(14) + (0.2e1 * MDP(13) - 0.2e1 * MDP(16)) * t413) * qJD(3); (t406 - t407) * MDP(9) + (-qJD(3) * t371 + t322) * MDP(10) + (-t313 * t413 + t344) * MDP(13) + t345 * MDP(14) + (pkin(3) * t286 - qJ(4) * t339 + (-t268 - t272) * t413) * MDP(15) + (t260 * t413 - t344) * MDP(16) + (t269 * t413 + 0.2e1 * t330 - t345) * MDP(17) + (-pkin(3) * t247 - qJ(4) * t244 - t260 * t269 - t267 * t272 - t412 * t268) * MDP(18) + (-t279 * t362 - t399) * MDP(19) + (-t252 - t423 + (t257 + t359) * t334) * MDP(20) + (t397 + t419) * MDP(21) + (t343 - t398) * MDP(22) + (qJ(4) * t258 + t229 * t334 + (-t250 + (t251 + t378) * t334) * t418 + t351 * t277 + t409 * t336) * MDP(24) + (-qJ(4) * t257 + t229 * t336 + (t336 * t378 + t390) * t418 + t351 * t279 - t409 * t334) * MDP(25) + (-t257 * t315 + t258 * t314 - t386 * t277 - t387 * t279 - t334 * t421 - t336 * t422) * MDP(26) + (-t213 * t314 - t212 * t315 + t224 * (pkin(5) * t334 + qJ(4)) + (pkin(5) * t414 + t351) * t227 + t386 * t218 + t387 * t216) * MDP(27) + (t413 * MDP(8) + qJD(3) * MDP(10) + t313 * MDP(14) + (t267 - t412) * MDP(15) + t269 * MDP(16) - t260 * MDP(17) + t418 * MDP(23) + t222 * MDP(24) - t223 * MDP(25)) * t299; (t289 + t322) * MDP(15) + t247 * MDP(18) - t282 * MDP(24) + (-t299 * MDP(16) - MDP(17) * t413 + t260 * MDP(18)) * t413 + (-MDP(15) * t371 - qJD(3) * MDP(17) + t268 * MDP(18) - t277 * MDP(24) - t279 * MDP(25) - MDP(27) * t227) * qJD(3) + (-MDP(25) * t363 + MDP(26) * t420 + MDP(27) * t422) * t336 + (t286 * MDP(25) + (t279 * t413 - t258 + t381) * MDP(26) + t421 * MDP(27) - MDP(24) * t363) * t334; (-t276 + t408) * MDP(20) + t388 * MDP(21) + (t279 * t418 + t283) * MDP(22) - t286 * MDP(23) + (t223 * t418 - t248 * t279 + t364) * MDP(24) + (t222 * t418 + t373) * MDP(25) + t391 * MDP(27) * t218 + (t257 * MDP(26) + (-t227 * t279 + t212) * MDP(27)) * pkin(5) + (t279 * MDP(19) + MDP(21) * t418 + t248 * MDP(25) - MDP(26) * t391) * t277 + ((-MDP(22) * qJD(3) - MDP(24) * t237) * t336 + (-MDP(21) * qJD(3) - t299 * MDP(22) - MDP(24) * t243 + MDP(25) * t237) * t334) * qJD(5); (-t276 - t408) * MDP(26) + (t216 * t279 + t218 * t277 + t224) * MDP(27);];
tauc  = t1;
