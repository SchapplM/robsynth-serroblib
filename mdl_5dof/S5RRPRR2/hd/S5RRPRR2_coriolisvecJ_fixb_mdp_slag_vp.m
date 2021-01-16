% Calculate Coriolis joint torque vector for
% S5RRPRR2
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
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:23:44
% EndTime: 2021-01-15 21:23:57
% DurationCPUTime: 4.30s
% Computational Cost: add. (2983->279), mult. (8028->375), div. (0->0), fcn. (6191->8), ass. (0->141)
t375 = sin(pkin(9));
t376 = cos(pkin(9));
t379 = sin(qJ(2));
t382 = cos(qJ(2));
t354 = t375 * t382 + t376 * t379;
t344 = t354 * qJD(2);
t333 = qJD(1) * t344;
t416 = qJD(2) * t379;
t407 = qJD(1) * t416;
t363 = t375 * t407;
t417 = qJD(1) * t382;
t408 = t376 * t417;
t334 = qJD(2) * t408 - t363;
t418 = qJD(1) * t379;
t342 = t375 * t418 - t408;
t378 = sin(qJ(4));
t381 = cos(qJ(4));
t414 = qJD(4) * t381;
t415 = qJD(4) * t378;
t420 = qJD(1) * t354;
t253 = -t378 * t333 + t381 * t334 - t342 * t414 - t415 * t420;
t301 = -t381 * t342 - t378 * t420;
t393 = t342 * t378 - t381 * t420;
t385 = t393 * qJD(4) - t381 * t333 - t334 * t378;
t372 = qJD(2) + qJD(4);
t429 = t301 * t372;
t430 = t393 * t372;
t380 = cos(qJ(5));
t292 = t301 * t380;
t377 = sin(qJ(5));
t413 = qJD(5) * t377;
t236 = qJD(5) * t292 + t380 * t253 + t377 * t385 + t393 * t413;
t256 = t301 * t377 - t380 * t393;
t260 = t377 * t393 + t292;
t386 = -qJD(5) * t256 - t253 * t377 + t380 * t385;
t371 = qJD(5) + t372;
t431 = t260 * t371;
t432 = t256 * t371;
t444 = (t236 - t431) * MDP(24) + (t386 + t432) * MDP(25) + (t256 ^ 2 - t260 ^ 2) * MDP(23) - t260 * t256 * MDP(22);
t457 = t444 + t301 * t393 * MDP(15) + (-t301 ^ 2 + t393 ^ 2) * MDP(16) + (t253 - t429) * MDP(17) + (t385 - t430) * MDP(18);
t433 = -qJ(3) - pkin(6);
t362 = t433 * t382;
t359 = qJD(1) * t362;
t348 = t375 * t359;
t361 = t433 * t379;
t358 = qJD(1) * t361;
t352 = qJD(2) * pkin(2) + t358;
t303 = t376 * t352 + t348;
t434 = pkin(7) * t420;
t279 = qJD(2) * pkin(3) + t303 - t434;
t428 = t376 * t359;
t304 = t375 * t352 - t428;
t435 = pkin(7) * t342;
t284 = t304 - t435;
t396 = -t279 * t378 - t284 * t381;
t452 = pkin(8) * t301;
t243 = -t396 + t452;
t368 = -pkin(2) * t382 - pkin(1);
t419 = qJD(1) * t368;
t360 = qJD(3) + t419;
t309 = pkin(3) * t342 + t360;
t270 = -pkin(4) * t301 + t309;
t456 = t243 * t413 - t270 * t260;
t405 = qJD(2) * t433;
t339 = qJD(3) * t382 + t379 * t405;
t321 = t339 * qJD(1);
t340 = -qJD(3) * t379 + t382 * t405;
t322 = t340 * qJD(1);
t285 = -t321 * t375 + t376 * t322;
t268 = -pkin(7) * t334 + t285;
t286 = t376 * t321 + t375 * t322;
t269 = -pkin(7) * t333 + t286;
t389 = -(qJD(4) * t279 + t269) * t381 - t378 * t268 + t284 * t415;
t232 = pkin(8) * t385 - t389;
t388 = t396 * qJD(4) + t381 * t268 - t378 * t269;
t233 = -pkin(8) * t253 + t388;
t448 = -t377 * t232 + t380 * t233 - t270 * t256;
t436 = pkin(4) * t393;
t451 = pkin(8) * t393;
t449 = (-t243 * t371 - t233) * t377 + t456;
t447 = -t309 * t301 + t389;
t446 = t309 * t393 + t388;
t307 = -t358 * t375 + t428;
t287 = t307 + t435;
t308 = t376 * t358 + t348;
t288 = t308 - t434;
t367 = pkin(2) * t376 + pkin(3);
t437 = pkin(2) * t375;
t392 = t367 * t381 - t378 * t437;
t443 = -t392 * qJD(4) + t378 * t287 + t381 * t288;
t338 = t367 * t378 + t381 * t437;
t442 = -t338 * qJD(4) - t381 * t287 + t288 * t378;
t441 = t382 * MDP(4) - pkin(1) * MDP(9);
t440 = qJD(5) - t371;
t439 = pkin(1) * t382 * MDP(10) + (t379 ^ 2 - t382 ^ 2) * MDP(5);
t426 = t380 * t243;
t424 = t442 + t452;
t423 = t443 + t451;
t295 = t376 * t339 + t375 * t340;
t312 = t375 * t361 - t376 * t362;
t410 = 0.2e1 * qJD(1);
t369 = pkin(2) * t418;
t365 = pkin(2) * t407;
t310 = pkin(3) * t333 + t365;
t317 = pkin(2) * t416 + pkin(3) * t344;
t316 = pkin(3) * t420 + t369;
t401 = t381 * t279 - t284 * t378;
t242 = t401 + t451;
t240 = pkin(4) * t372 + t242;
t406 = -pkin(4) * t371 - t240;
t294 = -t339 * t375 + t376 * t340;
t311 = t376 * t361 + t362 * t375;
t397 = -t377 * t240 - t426;
t296 = -pkin(7) * t354 + t311;
t353 = t375 * t379 - t376 * t382;
t297 = -pkin(7) * t353 + t312;
t395 = -t296 * t378 - t297 * t381;
t305 = t381 * t353 + t354 * t378;
t306 = -t353 * t378 + t354 * t381;
t265 = t305 * t380 + t306 * t377;
t266 = -t305 * t377 + t306 * t380;
t324 = pkin(3) * t353 + t368;
t347 = t353 * qJD(2);
t275 = pkin(7) * t347 + t294;
t276 = -pkin(7) * t344 + t295;
t390 = -t378 * t275 - t381 * t276 - t296 * t414 + t297 * t415;
t387 = t395 * qJD(4) + t381 * t275 - t276 * t378;
t337 = pkin(4) + t392;
t280 = pkin(4) * t305 + t324;
t271 = t316 - t436;
t264 = t306 * qJD(4) + t381 * t344 - t347 * t378;
t263 = -t305 * qJD(4) - t344 * t378 - t347 * t381;
t249 = pkin(4) * t264 + t317;
t248 = -pkin(4) * t385 + t310;
t247 = -pkin(8) * t305 - t395;
t246 = -pkin(8) * t306 + t296 * t381 - t297 * t378;
t239 = t266 * qJD(5) + t263 * t377 + t380 * t264;
t238 = -t265 * qJD(5) + t263 * t380 - t264 * t377;
t235 = -pkin(8) * t263 + t387;
t234 = -pkin(8) * t264 - t390;
t1 = [(t333 * t368 + t344 * t360) * MDP(11) + (t334 * t368 - t347 * t360) * MDP(12) + (-t285 * t354 - t286 * t353 - t294 * t420 - t295 * t342 + t303 * t347 - t304 * t344 - t311 * t334 - t312 * t333) * MDP(13) + (t285 * t311 + t286 * t312 + t294 * t303 + t295 * t304) * MDP(14) + (t253 * t306 - t263 * t393) * MDP(15) + (-t253 * t305 + t263 * t301 + t264 * t393 + t306 * t385) * MDP(16) + (t309 * t264 - t301 * t317 + t310 * t305 - t324 * t385) * MDP(20) + (t324 * t253 + t309 * t263 + t310 * t306 - t317 * t393) * MDP(21) + (t236 * t266 + t238 * t256) * MDP(22) + (-t236 * t265 + t238 * t260 - t239 * t256 + t266 * t386) * MDP(23) + (t270 * t239 + t248 * t265 - t249 * t260 - t280 * t386) * MDP(27) + (t280 * t236 + t270 * t238 + t248 * t266 + t249 * t256) * MDP(28) + (t263 * MDP(17) - t264 * MDP(18) + t387 * MDP(20) + t390 * MDP(21)) * t372 + (t238 * MDP(24) - t239 * MDP(25) + (-t234 * t377 + t235 * t380) * MDP(27) + (-t234 * t380 - t235 * t377) * MDP(28) + ((-t246 * t377 - t247 * t380) * MDP(27) + (-t246 * t380 + t247 * t377) * MDP(28)) * qJD(5)) * t371 + (t294 * MDP(11) - t295 * MDP(12) - t439 * t410 + (t441 * t410 + ((qJD(1) * t353 + t342) * MDP(11) + 0.2e1 * t420 * MDP(12) + (t360 + t419) * MDP(14)) * pkin(2)) * t379 + (MDP(6) * t382 - MDP(7) * t379 + (MDP(10) * t379 - MDP(9) * t382) * pkin(6)) * qJD(2)) * qJD(2); (-qJD(2) * t307 - t342 * t369 - t360 * t420 + t285) * MDP(11) + (qJD(2) * t308 + t342 * t360 - t369 * t420 - t286) * MDP(12) + ((t304 + t307) * t420 + (-t303 + t308) * t342 + (-t333 * t375 - t334 * t376) * pkin(2)) * MDP(13) + (-t303 * t307 - t304 * t308 + (t285 * t376 + t286 * t375 - t360 * t418) * pkin(2)) * MDP(14) + (t301 * t316 + t372 * t442 + t446) * MDP(20) + (t316 * t393 + t372 * t443 + t447) * MDP(21) + (t271 * t260 + (t423 * t377 + t424 * t380) * t371 + ((-t337 * t377 - t338 * t380) * t371 + t397) * qJD(5) + t448) * MDP(27) + (-t271 * t256 + (-t233 + (qJD(5) * t338 - t424) * t371) * t377 + (-qJD(5) * t240 - t232 + (-qJD(5) * t337 + t423) * t371) * t380 + t456) * MDP(28) + (-t379 * t441 + t439) * qJD(1) ^ 2 + t457; -t363 * MDP(12) + (-t342 ^ 2 - t420 ^ 2) * MDP(13) + (t303 * t420 + t304 * t342 + t365) * MDP(14) + (-t385 - t430) * MDP(20) + (t253 + t429) * MDP(21) + (-t386 + t432) * MDP(27) + (t236 + t431) * MDP(28) + ((t375 * t417 + t376 * t418 + t420) * MDP(11) + (-t342 + t408) * MDP(12)) * qJD(2); (-t396 * t372 + t446) * MDP(20) + (t401 * t372 + t447) * MDP(21) + (-t260 * t436 - (-t242 * t377 - t426) * t371 + (t406 * t377 - t426) * qJD(5) + t448) * MDP(27) + (t256 * t436 + (t406 * qJD(5) + t242 * t371 - t232) * t380 + t449) * MDP(28) + t457; (t397 * t440 + t448) * MDP(27) + ((-t240 * t440 - t232) * t380 + t449) * MDP(28) + t444;];
tauc = t1;
