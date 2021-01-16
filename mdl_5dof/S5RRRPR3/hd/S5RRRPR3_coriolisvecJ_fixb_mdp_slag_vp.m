% Calculate Coriolis joint torque vector for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:59:43
% EndTime: 2021-01-15 22:59:51
% DurationCPUTime: 2.28s
% Computational Cost: add. (2268->268), mult. (3972->359), div. (0->0), fcn. (2720->8), ass. (0->149)
t372 = qJD(1) + qJD(2);
t376 = cos(pkin(9));
t381 = cos(qJ(3));
t429 = t376 * t381;
t403 = t372 * t429;
t375 = sin(pkin(9));
t378 = sin(qJ(3));
t433 = t375 * t378;
t324 = t372 * t433 - t403;
t380 = cos(qJ(5));
t310 = t380 * t324;
t344 = t375 * t381 + t376 * t378;
t336 = t344 * qJD(3);
t316 = t372 * t336;
t416 = qJD(3) * t378;
t402 = t372 * t416;
t352 = t375 * t402;
t415 = qJD(3) * t381;
t401 = t372 * t415;
t317 = t376 * t401 - t352;
t326 = t344 * t372;
t377 = sin(qJ(5));
t413 = qJD(5) * t377;
t235 = -qJD(5) * t310 - t377 * t316 + t380 * t317 - t326 * t413;
t267 = -t326 * t377 - t310;
t389 = t324 * t377 - t380 * t326;
t384 = qJD(5) * t389 - t380 * t316 - t317 * t377;
t371 = qJD(3) + qJD(5);
t435 = t267 * t371;
t436 = t389 * t371;
t455 = t267 * MDP(18) * t389 + (-t267 ^ 2 + t389 ^ 2) * MDP(19) + (t235 - t435) * MDP(20) + (t384 - t436) * MDP(21);
t454 = qJ(4) + pkin(7);
t379 = sin(qJ(2));
t438 = pkin(1) * qJD(1);
t407 = t379 * t438;
t399 = t454 * t372 + t407;
t312 = t399 * t378;
t303 = qJD(3) * pkin(3) - t312;
t313 = t399 * t381;
t431 = t376 * t313;
t259 = t375 * t303 + t431;
t443 = pkin(8) * t324;
t245 = t259 - t443;
t365 = -pkin(3) * t381 - pkin(2);
t382 = cos(qJ(2));
t406 = t382 * t438;
t323 = t365 * t372 + qJD(4) - t406;
t274 = pkin(4) * t324 + t323;
t453 = t245 * t413 - t274 * t267;
t437 = pkin(1) * qJD(2);
t404 = qJD(1) * t437;
t395 = t382 * t404;
t386 = qJD(4) * t372 + t395;
t388 = qJD(3) * t399;
t272 = -t378 * t388 + t381 * t386;
t273 = -t378 * t386 - t381 * t388;
t237 = -t272 * t375 + t376 * t273;
t229 = -pkin(8) * t317 + t237;
t238 = t376 * t272 + t375 * t273;
t230 = -pkin(8) * t316 + t238;
t452 = t380 * t229 - t377 * t230 + t274 * t389;
t450 = MDP(7) * t378;
t449 = MDP(8) * (t378 ^ 2 - t381 ^ 2);
t368 = t381 * qJD(4);
t400 = qJD(3) * t454;
t333 = -t378 * t400 + t368;
t334 = -qJD(4) * t378 - t381 * t400;
t423 = -t333 * t375 + t376 * t334 + t344 * t406;
t343 = -t429 + t433;
t422 = t376 * t333 + t375 * t334 + t343 * t406;
t448 = t378 * MDP(12) + t381 * MDP(13);
t447 = qJD(5) - t371;
t446 = pkin(1) * t382;
t445 = pkin(3) * t375;
t444 = pkin(3) * t378;
t442 = pkin(8) * t326;
t337 = t343 * qJD(3);
t441 = pkin(8) * t337;
t440 = pkin(8) * t344;
t434 = t372 * t378;
t297 = t375 * t313;
t383 = qJD(3) ^ 2;
t428 = t378 * t383;
t427 = t381 * t383;
t363 = pkin(1) * t379 + pkin(7);
t426 = -qJ(4) - t363;
t284 = t343 * t380 + t344 * t377;
t250 = -qJD(5) * t284 - t336 * t377 - t337 * t380;
t359 = t379 * t404;
t335 = pkin(3) * t402 + t359;
t279 = pkin(4) * t316 + t335;
t285 = -t343 * t377 + t344 * t380;
t425 = t274 * t250 + t279 * t285;
t251 = qJD(5) * t285 + t380 * t336 - t337 * t377;
t424 = t274 * t251 + t279 * t284;
t421 = t323 * t336 + t335 * t343;
t420 = -t323 * t337 + t335 * t344;
t396 = qJD(3) * t426;
t405 = t382 * t437;
t292 = t378 * t396 + t381 * t405 + t368;
t293 = (-qJD(4) - t405) * t378 + t381 * t396;
t255 = t376 * t292 + t375 * t293;
t261 = -t376 * t312 - t297;
t341 = t426 * t378;
t369 = t381 * qJ(4);
t342 = t363 * t381 + t369;
t283 = t375 * t341 + t376 * t342;
t351 = -pkin(2) * t372 - t406;
t419 = t351 * t415 + t378 * t359;
t356 = t454 * t378;
t357 = pkin(7) * t381 + t369;
t302 = -t375 * t356 + t376 * t357;
t414 = qJD(3) * t382;
t411 = t381 * MDP(12);
t409 = -qJD(2) + t372;
t408 = pkin(3) * t434;
t366 = pkin(3) * t416;
t309 = pkin(4) * t336 + t366;
t254 = -t292 * t375 + t376 * t293;
t258 = t376 * t303 - t297;
t260 = t312 * t375 - t431;
t282 = t376 * t341 - t342 * t375;
t301 = -t376 * t356 - t357 * t375;
t394 = -t237 * t344 - t238 * t343 + t258 * t337 - t259 * t336;
t393 = t309 - t407;
t340 = t343 * pkin(8);
t392 = qJD(5) * (-t340 + t302) - t441 - t423;
t332 = t336 * pkin(8);
t391 = -qJD(5) * (t301 - t440) + t332 - t422;
t242 = qJD(3) * pkin(4) + t258 - t442;
t390 = -t377 * t242 - t380 * t245;
t319 = pkin(4) * t343 + t365;
t385 = -MDP(10) * t428 + (-t235 * t284 + t250 * t267 + t251 * t389 + t285 * t384) * MDP(19) + (t235 * t285 - t250 * t389) * MDP(18) - 0.2e1 * t372 * qJD(3) * t449 + 0.2e1 * t401 * t450 + MDP(9) * t427 + (t250 * MDP(20) - t251 * MDP(21)) * t371;
t367 = t379 * t437;
t364 = -pkin(2) - t446;
t361 = pkin(3) * t376 + pkin(4);
t353 = t365 - t446;
t349 = t367 + t366;
t338 = t351 * t416;
t308 = t319 - t446;
t296 = t309 + t367;
t291 = pkin(4) * t326 + t408;
t263 = -t340 + t283;
t262 = t282 - t440;
t247 = t261 - t442;
t246 = t260 + t443;
t240 = -t332 + t255;
t239 = t254 + t441;
t1 = [t385 + (-t296 * t267 - t308 * t384 + (t239 * t380 - t240 * t377 + (-t262 * t377 - t263 * t380) * qJD(5)) * t371 + t424) * MDP(23) + (-t296 * t389 + t308 * t235 - (t239 * t377 + t240 * t380 + (t262 * t380 - t263 * t377) * qJD(5)) * t371 + t425) * MDP(24) + (qJD(3) * t254 + t316 * t353 + t324 * t349 + t421) * MDP(14) + (t237 * t282 + t238 * t283 + t254 * t258 + t255 * t259 + t323 * t349 + t335 * t353) * MDP(17) + (((-qJD(1) - t372) * MDP(6) - t448 * qJD(3)) * t382 + (-qJD(1) * t411 + (t378 * MDP(13) - MDP(5) - t411) * t372) * t379) * t437 + (-t363 * t427 + t364 * t402 + t338) * MDP(12) + (t363 * t428 + t364 * t401 + t419) * MDP(13) + (-qJD(3) * t255 + t317 * t353 + t326 * t349 + t420) * MDP(15) + (-t254 * t326 - t255 * t324 - t282 * t317 - t283 * t316 + t394) * MDP(16) - t359 * MDP(5); t409 * MDP(6) * t406 + (t319 * t235 + (t377 * t392 + t380 * t391) * t371 - t393 * t389 + t425) * MDP(24) + (-t319 * t384 + (t377 * t391 - t380 * t392) * t371 - t393 * t267 + t424) * MDP(23) + (-pkin(2) * t401 + pkin(7) * t428 + (-t379 * t434 + t381 * t414) * t438 + t419) * MDP(13) + t385 + (-pkin(2) * t402 - pkin(7) * t427 + t338 + (t379 * t381 * t409 + t378 * t414) * t438) * MDP(12) + (t372 * t407 - t359) * MDP(5) + (t237 * t301 + t238 * t302 + t335 * t365 + (t366 - t407) * t323 + t422 * t259 + t423 * t258) * MDP(17) + (-t326 * t407 + t317 * t365 + (t326 * t444 - t422) * qJD(3) + t420) * MDP(15) + (-t301 * t317 - t302 * t316 - t324 * t422 - t326 * t423 + t394) * MDP(16) + (-t324 * t407 + t316 * t365 + (t324 * t444 + t423) * qJD(3) + t421) * MDP(14); (-qJD(3) * t260 - t323 * t326 - t324 * t408 + t237) * MDP(14) + (qJD(3) * t261 + t323 * t324 - t326 * t408 - t238) * MDP(15) + ((t259 + t260) * t326 + (-t258 + t261) * t324 + (-t316 * t375 - t317 * t376) * pkin(3)) * MDP(16) + (-t258 * t260 - t259 * t261 + (t237 * t376 + t238 * t375 - t323 * t434) * pkin(3)) * MDP(17) + (t291 * t267 - (t246 * t380 - t247 * t377) * t371 + ((-t361 * t377 - t380 * t445) * t371 + t390) * qJD(5) + t452) * MDP(23) + (-t380 * t230 - t377 * t229 + t291 * t389 + (t246 * t377 + t247 * t380) * t371 + (-(t361 * t380 - t377 * t445) * t371 - t380 * t242) * qJD(5) + t453) * MDP(24) + t448 * (-t351 * t372 - t395) + (-t381 * t450 + t449) * t372 ^ 2 + t455; -t352 * MDP(15) + (-t324 ^ 2 - t326 ^ 2) * MDP(16) + (t258 * t326 + t259 * t324 + t335) * MDP(17) + (-t384 - t436) * MDP(23) + (t235 + t435) * MDP(24) + (0.2e1 * t326 * MDP(14) + (-t324 + t403) * MDP(15)) * qJD(3); (t447 * t390 + t452) * MDP(23) + ((-t245 * t371 - t229) * t377 + (-t447 * t242 - t230) * t380 + t453) * MDP(24) + t455;];
tauc = t1;
