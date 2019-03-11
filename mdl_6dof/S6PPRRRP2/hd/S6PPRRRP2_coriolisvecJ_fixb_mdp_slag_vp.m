% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPRRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:58:16
% EndTime: 2019-03-08 18:58:23
% DurationCPUTime: 3.49s
% Computational Cost: add. (3246->363), mult. (8516->516), div. (0->0), fcn. (7128->12), ass. (0->158)
t331 = cos(pkin(6));
t318 = qJD(1) * t331 + qJD(2);
t326 = sin(pkin(12));
t328 = sin(pkin(6));
t334 = sin(qJ(3));
t337 = cos(qJ(3));
t329 = cos(pkin(12));
t330 = cos(pkin(7));
t413 = t329 * t330;
t344 = (t326 * t337 + t334 * t413) * t328;
t327 = sin(pkin(7));
t415 = t327 * t334;
t279 = qJD(1) * t344 + t318 * t415;
t333 = sin(qJ(4));
t336 = cos(qJ(4));
t364 = pkin(4) * t333 - pkin(10) * t336;
t313 = t364 * qJD(4);
t442 = t279 - t313;
t414 = t327 * t337;
t437 = (-t326 * t334 + t337 * t413) * t328;
t441 = t331 * t414 + t437;
t403 = MDP(12) * t336;
t398 = qJD(3) * t336;
t319 = -qJD(5) + t398;
t440 = qJD(1) * t437 + t318 * t414;
t439 = MDP(6) * t333;
t324 = t333 ^ 2;
t438 = MDP(7) * (-t336 ^ 2 + t324);
t272 = t440 * qJD(3);
t276 = qJD(3) * pkin(9) + t279;
t402 = qJD(1) * t328;
t376 = t329 * t402;
t295 = t318 * t330 - t327 * t376;
t433 = -t333 * t276 + t295 * t336;
t234 = qJD(4) * t433 + t272 * t336;
t255 = t336 * t276 + t333 * t295;
t253 = qJD(4) * pkin(10) + t255;
t365 = t330 * t376;
t399 = qJD(3) * t334;
t375 = t327 * t399;
t377 = t326 * t402;
t397 = qJD(3) * t337;
t273 = t318 * t375 + t365 * t399 + t377 * t397;
t261 = qJD(3) * t313 + t273;
t278 = -t334 * t377 + t337 * (t318 * t327 + t365);
t316 = -pkin(4) * t336 - pkin(10) * t333 - pkin(3);
t264 = qJD(3) * t316 - t278;
t332 = sin(qJ(5));
t335 = cos(qJ(5));
t392 = qJD(5) * t335;
t393 = qJD(5) * t332;
t366 = -t332 * t234 - t253 * t392 + t335 * t261 - t264 * t393;
t387 = qJD(3) * qJD(4);
t369 = t333 * t387;
t224 = -pkin(5) * t369 - t366;
t231 = t253 * t335 + t264 * t332;
t229 = -qJ(6) * t319 + t231;
t436 = t229 * t319 + t224;
t412 = t332 * t336;
t435 = -t278 * t412 + t316 * t393 + t442 * t335;
t411 = t335 * t336;
t434 = -t278 * t411 + t316 * t392 - t442 * t332;
t432 = qJD(3) * t279 - t273;
t431 = MDP(18) + MDP(20);
t430 = MDP(19) - MDP(22);
t429 = qJD(3) * pkin(3);
t394 = qJD(4) * t336;
t395 = qJD(4) * t333;
t235 = t333 * t272 + t276 * t394 + t295 * t395;
t368 = t336 * t387;
t371 = t333 * t393;
t386 = qJD(4) * qJD(5);
t293 = qJD(3) * t371 + (-t368 - t386) * t335;
t367 = t332 * t386;
t370 = t333 * t392;
t294 = t367 + (t332 * t394 + t370) * qJD(3);
t396 = qJD(4) * t332;
t400 = qJD(3) * t333;
t309 = t335 * t400 + t396;
t227 = pkin(5) * t294 + qJ(6) * t293 - qJD(6) * t309 + t235;
t428 = t227 * t332;
t427 = t227 * t335;
t425 = t235 * t332;
t424 = t235 * t335;
t389 = t335 * qJD(4);
t307 = t332 * t400 - t389;
t423 = t278 * t307;
t422 = t278 * t309;
t421 = t293 * t332;
t419 = t307 * t319;
t418 = t316 * t335;
t417 = t319 * t332;
t416 = t319 * t335;
t360 = pkin(5) * t332 - qJ(6) * t335;
t410 = qJD(6) * t332 + t319 * t360 + t255;
t391 = qJD(5) * t336;
t409 = -pkin(5) * t395 + (-t332 * t395 + t335 * t391) * pkin(9) + t435;
t408 = -qJ(6) * t395 + qJD(6) * t336 - (-t332 * t391 - t333 * t389) * pkin(9) - t434;
t312 = t364 * qJD(3);
t407 = t332 * t312 + t335 * t433;
t405 = pkin(9) * t411 + t332 * t316;
t390 = qJD(6) * t319;
t230 = -t253 * t332 + t264 * t335;
t388 = qJD(6) - t230;
t383 = pkin(10) * t417;
t382 = pkin(10) * t416;
t381 = pkin(10) * t395;
t380 = pkin(10) * t389;
t374 = t327 * t397;
t373 = t319 * t392;
t372 = t319 * t393;
t361 = pkin(5) * t335 + qJ(6) * t332;
t228 = pkin(5) * t319 + t388;
t359 = t228 * t335 - t229 * t332;
t358 = t312 * t335 - t332 * t433;
t283 = t331 * t415 + t344;
t298 = -t327 * t328 * t329 + t330 * t331;
t263 = t283 * t336 + t298 * t333;
t243 = t263 * t335 - t332 * t441;
t242 = t263 * t332 + t335 * t441;
t262 = t283 * t333 - t298 * t336;
t356 = qJD(3) * t324 - t319 * t336;
t355 = pkin(9) + t360;
t252 = -qJD(4) * pkin(4) - t433;
t338 = qJD(4) ^ 2;
t354 = pkin(9) * t338 - t432;
t275 = -t278 - t429;
t353 = qJD(4) * (t275 + t278 - t429);
t300 = t330 * t333 + t336 * t415;
t287 = t300 * t332 + t335 * t414;
t288 = t300 * t335 - t332 * t414;
t299 = -t330 * t336 + t333 * t415;
t350 = -MDP(11) * t336 + MDP(12) * t333 - MDP(4);
t348 = -t335 * t234 + t253 * t393 - t332 * t261 - t264 * t392;
t341 = -t230 * t319 + t348;
t339 = qJD(3) ^ 2;
t315 = -pkin(4) - t361;
t296 = t355 * t333;
t292 = -t418 + (pkin(9) * t332 + pkin(5)) * t336;
t291 = -qJ(6) * t336 + t405;
t286 = qJD(4) * t300 + t333 * t374;
t285 = -qJD(4) * t299 + t336 * t374;
t284 = pkin(5) * t309 + qJ(6) * t307;
t281 = t283 * qJD(3);
t280 = t441 * qJD(3);
t271 = -t293 - t419;
t266 = (qJD(5) * t361 - qJD(6) * t335) * t333 + t355 * t394;
t251 = qJD(5) * t288 + t285 * t332 - t335 * t375;
t250 = -qJD(5) * t287 + t285 * t335 + t332 * t375;
t241 = -pkin(5) * t400 - t358;
t240 = qJ(6) * t400 + t407;
t239 = -qJD(4) * t262 + t280 * t336;
t238 = qJD(4) * t263 + t280 * t333;
t237 = pkin(5) * t307 - qJ(6) * t309 + t252;
t226 = -qJD(5) * t242 + t239 * t335 + t281 * t332;
t225 = qJD(5) * t243 + t239 * t332 - t281 * t335;
t223 = qJ(6) * t369 - t348 - t390;
t1 = [(t225 * t309 - t226 * t307 - t242 * t293 - t243 * t294) * MDP(21) + (t223 * t243 + t224 * t242 + t225 * t228 + t226 * t229 + t227 * t262 + t237 * t238) * MDP(23) + (-MDP(11) * t238 - MDP(12) * t239) * qJD(4) + (-t280 * MDP(5) + t350 * t281 + (-t441 * t403 + (-MDP(11) * t441 - t242 * t431 - t243 * t430) * t333) * qJD(4)) * qJD(3) + t431 * (t225 * t319 + t238 * t307 + t262 * t294) + t430 * (t226 * t319 + t238 * t309 - t262 * t293); (-t250 * t307 + t251 * t309 - t287 * t293 - t288 * t294) * MDP(21) + (t223 * t288 + t224 * t287 + t227 * t299 + t228 * t251 + t229 * t250 + t237 * t286) * MDP(23) + (-t286 * MDP(11) - t285 * MDP(12) + (-t287 * t431 - t288 * t430) * t400) * qJD(4) + ((-MDP(11) * t333 - t403) * t337 * t387 + (-t337 * MDP(5) + t334 * t350) * t339) * t327 + t431 * (t251 * t319 + t286 * t307 + t299 * t294) + t430 * (t250 * t319 + t286 * t309 - t293 * t299); t432 * MDP(4) + (t278 - t440) * qJD(3) * MDP(5) + 0.2e1 * t368 * t439 - 0.2e1 * t387 * t438 + (t333 * t353 - t336 * t354) * MDP(11) + (t333 * t354 + t336 * t353) * MDP(12) + (-t293 * t333 * t335 + (t336 * t389 - t371) * t309) * MDP(13) + ((-t307 * t335 - t309 * t332) * t394 + (t421 - t294 * t335 + (t307 * t332 - t309 * t335) * qJD(5)) * t333) * MDP(14) + (t319 * t371 + t293 * t336 + (t309 * t333 + t335 * t356) * qJD(4)) * MDP(15) + (t319 * t370 + t294 * t336 + (-t307 * t333 - t332 * t356) * qJD(4)) * MDP(16) + (-t319 - t398) * MDP(17) * t395 + (t435 * t319 + (t252 * t396 + (qJD(4) * t307 + t373) * pkin(9) - t366) * t336 + (t252 * t392 + pkin(9) * t294 + t425 - t423 + (-pkin(9) * t417 + (-pkin(9) * t412 + t418) * qJD(3) + t230) * qJD(4)) * t333) * MDP(18) + (t434 * t319 + (t252 * t389 + (qJD(4) * t309 - t372) * pkin(9) - t348) * t336 + (-t252 * t393 - pkin(9) * t293 + t424 - t422 + (-pkin(9) * t416 - qJD(3) * t405 - t231) * qJD(4)) * t333) * MDP(19) + (t266 * t307 + t294 * t296 + (t237 * t396 + t224) * t336 + t409 * t319 + (t237 * t392 + t428 - t423 + (-qJD(3) * t292 - t228) * qJD(4)) * t333) * MDP(20) + (-t291 * t294 - t292 * t293 + t409 * t309 + t408 * t307 + t359 * t394 + (-t223 * t332 + t224 * t335 + (-t228 * t332 - t229 * t335) * qJD(5)) * t333) * MDP(21) + (-t266 * t309 + t293 * t296 + (-t237 * t389 - t223) * t336 + t408 * t319 + (t237 * t393 - t427 + t422 + (qJD(3) * t291 + t229) * qJD(4)) * t333) * MDP(22) + (t223 * t291 + t224 * t292 + t227 * t296 + (-t278 * t333 + t266) * t237 - t408 * t229 + t409 * t228) * MDP(23) + (MDP(8) * t336 - MDP(9) * t333) * t338; (qJD(4) * t255 - t275 * t400 - t235) * MDP(11) + (-qJD(3) * t275 - t272) * t403 + (-t309 * t416 - t421) * MDP(13) + ((-t293 + t419) * t335 + (t309 * t319 - t294) * t332) * MDP(14) + (-t373 + (t319 * t411 + (-t309 + t396) * t333) * qJD(3)) * MDP(15) + (t372 + (-t319 * t412 + (t307 + t389) * t333) * qJD(3)) * MDP(16) + t319 * MDP(17) * t400 + (-pkin(4) * t294 - t424 + t358 * t319 - t255 * t307 + (t252 * t332 + t382) * qJD(5) + (-t230 * t333 + (-t252 * t336 - t381) * t332) * qJD(3)) * MDP(18) + (pkin(4) * t293 + t425 - t407 * t319 - t255 * t309 + (t252 * t335 - t383) * qJD(5) + (-t252 * t411 + (t231 - t380) * t333) * qJD(3)) * MDP(19) + (-t427 - t241 * t319 + t294 * t315 - t410 * t307 + (t237 * t332 + t382) * qJD(5) + (t228 * t333 + (-t237 * t336 - t381) * t332) * qJD(3)) * MDP(20) + (t240 * t307 - t241 * t309 + (t223 - t319 * t228 + (qJD(5) * t309 - t294) * pkin(10)) * t335 + ((qJD(5) * t307 - t293) * pkin(10) + t436) * t332) * MDP(21) + (-t428 + t240 * t319 + t293 * t315 + t410 * t309 + (-t237 * t335 + t383) * qJD(5) + (t237 * t411 + (-t229 + t380) * t333) * qJD(3)) * MDP(22) + (t227 * t315 - t228 * t241 - t229 * t240 - t410 * t237 + (qJD(5) * t359 + t223 * t335 + t224 * t332) * pkin(10)) * MDP(23) + (-t336 * t439 + t438) * t339; t271 * MDP(15) - MDP(16) * t367 + t341 * MDP(19) + (pkin(5) * t293 - qJ(6) * t294) * MDP(21) + (-t341 - 0.2e1 * t390) * MDP(22) + (-pkin(5) * t224 + qJ(6) * t223 - t228 * t231 + t229 * t388 - t237 * t284) * MDP(23) + (-t319 * MDP(16) - t252 * MDP(18) - t237 * MDP(20) + (t229 - t231) * MDP(21) + t284 * MDP(22) + MDP(14) * t309) * t309 + (-MDP(16) * t370 + (-MDP(16) * t412 + (0.2e1 * pkin(5) * MDP(20) + 0.2e1 * qJ(6) * MDP(22) + MDP(17)) * t333) * qJD(4)) * qJD(3) + (t309 * MDP(13) + t252 * MDP(19) - t284 * MDP(20) + (t228 - t388) * MDP(21) - t237 * MDP(22) - MDP(14) * t307) * t307 + t431 * (-t231 * t319 + t366); (t307 * t309 - t369) * MDP(20) + t271 * MDP(21) + (-t309 ^ 2 - t319 ^ 2) * MDP(22) + (t237 * t309 + t436) * MDP(23);];
tauc  = t1;
