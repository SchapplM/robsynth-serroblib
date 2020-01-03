% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:20
% EndTime: 2019-12-31 20:18:27
% DurationCPUTime: 3.60s
% Computational Cost: add. (3193->295), mult. (8368->408), div. (0->0), fcn. (6437->8), ass. (0->152)
t370 = cos(qJ(5));
t409 = qJD(5) * t370;
t365 = sin(pkin(9));
t366 = cos(pkin(9));
t369 = sin(qJ(2));
t372 = cos(qJ(2));
t345 = -t365 * t369 + t366 * t372;
t335 = t345 * qJD(1);
t371 = cos(qJ(4));
t323 = t371 * t335;
t346 = t365 * t372 + t366 * t369;
t337 = t346 * qJD(1);
t368 = sin(qJ(4));
t294 = -t337 * t368 + t323;
t454 = t294 * t370;
t458 = t409 - t454;
t384 = t335 * t368 + t337 * t371;
t367 = sin(qJ(5));
t410 = qJD(5) * t367;
t336 = t346 * qJD(2);
t326 = qJD(1) * t336;
t405 = qJD(1) * qJD(2);
t400 = t372 * t405;
t401 = t369 * t405;
t327 = -t365 * t401 + t366 * t400;
t411 = qJD(4) * t368;
t262 = qJD(4) * t323 - t326 * t368 + t327 * t371 - t337 * t411;
t362 = qJD(2) + qJD(4);
t415 = t262 * t370 + t362 * t409;
t244 = -t384 * t410 + t415;
t243 = t244 * t370;
t282 = t362 * t367 + t370 * t384;
t433 = t262 * t367;
t245 = qJD(5) * t282 + t433;
t425 = t384 * t367;
t279 = -t362 * t370 + t425;
t457 = -t367 * t245 - t279 * t458 + t243;
t242 = t244 * t367;
t263 = qJD(4) * t384 + t326 * t371 + t327 * t368;
t406 = -qJD(5) + t294;
t259 = t367 * t263;
t416 = -t406 * t409 + t259;
t422 = t370 * t406;
t427 = t294 * t362;
t429 = t384 * t362;
t431 = t282 * t384;
t456 = (-t263 + t429) * MDP(16) - t294 ^ 2 * MDP(14) + (t294 * t422 + t416 - t431) * MDP(22) + (-MDP(13) * t294 + MDP(14) * t384 + MDP(24) * t406) * t384 + (t262 - t427) * MDP(15) + (t282 * t458 + t242) * MDP(20);
t455 = t294 * t367;
t265 = pkin(4) * t384 - pkin(8) * t294;
t437 = -qJ(3) - pkin(6);
t399 = qJD(2) * t437;
t332 = qJD(3) * t372 + t369 * t399;
t314 = t332 * qJD(1);
t333 = -qJD(3) * t369 + t372 * t399;
t315 = t333 * qJD(1);
t278 = -t314 * t365 + t315 * t366;
t269 = -pkin(7) * t327 + t278;
t281 = t314 * t366 + t315 * t365;
t270 = -pkin(7) * t326 + t281;
t354 = t437 * t372;
t351 = qJD(1) * t354;
t340 = t365 * t351;
t353 = t437 * t369;
t350 = qJD(1) * t353;
t436 = qJD(2) * pkin(2);
t344 = t350 + t436;
t296 = t344 * t366 + t340;
t438 = pkin(7) * t337;
t274 = qJD(2) * pkin(3) + t296 - t438;
t424 = t366 * t351;
t297 = t344 * t365 - t424;
t439 = pkin(7) * t335;
t277 = t297 + t439;
t232 = (qJD(4) * t274 + t270) * t371 + t269 * t368 - t277 * t411;
t403 = -pkin(2) * t372 - pkin(1);
t388 = t403 * qJD(1);
t352 = qJD(3) + t388;
t302 = -pkin(3) * t335 + t352;
t453 = -t294 * t302 - t232;
t450 = -0.2e1 * t405;
t448 = MDP(4) * t369;
t447 = MDP(5) * (t369 ^ 2 - t372 ^ 2);
t432 = t279 * t384;
t261 = t370 * t263;
t445 = -t406 * t410 - t261;
t251 = t274 * t368 + t277 * t371;
t233 = qJD(4) * t251 - t269 * t371 + t270 * t368;
t250 = t274 * t371 - t277 * t368;
t248 = -pkin(4) * t362 - t250;
t249 = pkin(8) * t362 + t251;
t252 = -pkin(4) * t294 - pkin(8) * t384 + t302;
t386 = t249 * t367 - t252 * t370;
t444 = -t233 * t370 + t248 * t410 + t384 * t386;
t235 = t249 * t370 + t252 * t367;
t443 = t233 * t367 + t235 * t384 + t248 * t409;
t442 = -t302 * t384 - t233;
t286 = -t332 * t365 + t333 * t366;
t339 = t345 * qJD(2);
t271 = -pkin(7) * t339 + t286;
t287 = t332 * t366 + t333 * t365;
t272 = -pkin(7) * t336 + t287;
t304 = t353 * t366 + t354 * t365;
t288 = -pkin(7) * t346 + t304;
t305 = t353 * t365 - t354 * t366;
t289 = pkin(7) * t345 + t305;
t385 = t288 * t371 - t289 * t368;
t236 = qJD(4) * t385 + t271 * t368 + t272 * t371;
t299 = t345 * t368 + t346 * t371;
t317 = -pkin(3) * t345 + t403;
t383 = t345 * t371 - t346 * t368;
t256 = -pkin(4) * t383 - pkin(8) * t299 + t317;
t258 = t288 * t368 + t289 * t371;
t266 = qJD(4) * t383 - t336 * t368 + t339 * t371;
t441 = t233 * t299 + t248 * t266 - t258 * t263 + (qJD(5) * t256 + t236) * t406 + (qJD(5) * t252 + t232) * t383;
t440 = pkin(2) * t365;
t435 = t248 * t299;
t434 = t256 * t263;
t430 = t282 * t367;
t373 = qJD(2) ^ 2;
t423 = t369 * t373;
t421 = t372 * t373;
t374 = qJD(1) ^ 2;
t420 = t372 * t374;
t300 = -t350 * t365 + t424;
t283 = t300 - t439;
t301 = t350 * t366 + t340;
t284 = t301 - t438;
t359 = pkin(2) * t366 + pkin(3);
t380 = t359 * t371 - t368 * t440;
t417 = -t380 * qJD(4) + t283 * t368 + t284 * t371;
t381 = t359 * t368 + t371 * t440;
t414 = qJD(4) * t381 + t283 * t371 - t284 * t368;
t407 = t369 * qJD(1);
t361 = t369 * t436;
t357 = pkin(2) * t401;
t303 = pkin(3) * t326 + t357;
t310 = pkin(3) * t336 + t361;
t309 = pkin(2) * t407 + pkin(3) * t337;
t398 = pkin(1) * t450;
t393 = t406 * t367;
t331 = pkin(8) + t381;
t389 = qJD(5) * t331 + t265 + t309;
t387 = -t248 * t294 - t263 * t331;
t382 = -t406 * t455 - t445;
t379 = t266 * t370 - t299 * t410;
t330 = -pkin(4) - t380;
t267 = qJD(4) * t299 + t336 * t371 + t339 * t368;
t240 = pkin(4) * t267 - pkin(8) * t266 + t310;
t239 = pkin(4) * t263 - pkin(8) * t262 + t303;
t238 = t370 * t239;
t237 = qJD(4) * t258 - t271 * t371 + t272 * t368;
t1 = [0.2e1 * t400 * t448 + t447 * t450 + MDP(6) * t421 - MDP(7) * t423 + (-pkin(6) * t421 + t369 * t398) * MDP(9) + (pkin(6) * t423 + t372 * t398) * MDP(10) + (-t278 * t346 + t281 * t345 - t286 * t337 + t287 * t335 - t296 * t339 - t297 * t336 - t304 * t327 - t305 * t326) * MDP(11) + (t278 * t304 + t281 * t305 + t296 * t286 + t297 * t287 + (t352 + t388) * t361) * MDP(12) + (t262 * t299 + t266 * t384) * MDP(13) + (t262 * t383 - t263 * t299 + t266 * t294 - t267 * t384) * MDP(14) + (t263 * t317 + t267 * t302 - t294 * t310 - t303 * t383) * MDP(18) + (t262 * t317 + t266 * t302 + t299 * t303 + t310 * t384) * MDP(19) + (t243 * t299 + t282 * t379) * MDP(20) + ((-t279 * t370 - t430) * t266 + (-t242 - t245 * t370 + (t279 * t367 - t282 * t370) * qJD(5)) * t299) * MDP(21) + (-t244 * t383 + t261 * t299 + t267 * t282 - t379 * t406) * MDP(22) + (-t299 * t259 + t245 * t383 - t267 * t279 - (-t266 * t367 - t299 * t409) * t406) * MDP(23) + (-t263 * t383 - t267 * t406) * MDP(24) + (-t386 * t267 + t237 * t279 - t238 * t383 - t385 * t245 + (-t240 * t406 + t434 + (t249 * t383 + t258 * t406 + t435) * qJD(5)) * t370 + t441 * t367) * MDP(25) + (-t235 * t267 + t237 * t282 - t385 * t244 + ((-qJD(5) * t258 + t240) * t406 - t434 + (-qJD(5) * t249 + t239) * t383 - qJD(5) * t435) * t367 + t441 * t370) * MDP(26) + (MDP(15) * t266 - MDP(16) * t267 - MDP(18) * t237 - MDP(19) * t236) * t362; -t420 * t448 + t374 * t447 + ((t297 + t300) * t337 + (t296 - t301) * t335 + (-t326 * t365 - t327 * t366) * pkin(2)) * MDP(11) + (-t296 * t300 - t297 * t301 + (t278 * t366 + t281 * t365 - t352 * t407) * pkin(2)) * MDP(12) + (t294 * t309 - t362 * t414 + t442) * MDP(18) + (-t309 * t384 + t362 * t417 + t453) * MDP(19) + (t406 * t430 + t457) * MDP(21) + (t382 + t432) * MDP(23) + (t330 * t245 + t387 * t367 + t414 * t279 - (t367 * t417 - t370 * t389) * t406 + t444) * MDP(25) + (t330 * t244 + t387 * t370 + t414 * t282 - (t367 * t389 + t370 * t417) * t406 + t443) * MDP(26) + (MDP(9) * t369 * t374 + MDP(10) * t420) * pkin(1) + t456; (-t335 ^ 2 - t337 ^ 2) * MDP(11) + (t296 * t337 - t297 * t335 + t357) * MDP(12) + (t263 + t429) * MDP(18) + (t262 + t427) * MDP(19) + (t382 - t432) * MDP(25) + (-t406 * t422 - t259 - t431) * MDP(26); (t251 * t362 + t442) * MDP(18) + (t250 * t362 + t453) * MDP(19) + (t282 * t393 + t457) * MDP(21) + (-t393 * t406 + t261 + t432) * MDP(23) + (-pkin(4) * t245 + (-t250 * t367 + t265 * t370) * t406 - t251 * t279 - t248 * t455 - t416 * pkin(8) + t444) * MDP(25) + (-pkin(4) * t244 - (t250 * t370 + t265 * t367) * t406 - t251 * t282 - t248 * t454 + t445 * pkin(8) + t443) * MDP(26) + t456; t282 * t279 * MDP(20) + (-t279 ^ 2 + t282 ^ 2) * MDP(21) + (-t279 * t406 + t415) * MDP(22) + (-t282 * t406 - t433) * MDP(23) + t263 * MDP(24) + (-t232 * t367 - t235 * t406 - t248 * t282 + t238) * MDP(25) + (-t232 * t370 - t239 * t367 + t248 * t279 + t386 * t406) * MDP(26) + (-MDP(22) * t425 - MDP(23) * t282 - MDP(25) * t235 + MDP(26) * t386) * qJD(5);];
tauc = t1;
