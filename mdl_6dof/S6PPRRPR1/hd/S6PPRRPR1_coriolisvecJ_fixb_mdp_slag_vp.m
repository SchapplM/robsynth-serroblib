% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPRRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:47:35
% EndTime: 2019-03-08 18:47:42
% DurationCPUTime: 3.95s
% Computational Cost: add. (2807->365), mult. (7685->550), div. (0->0), fcn. (6662->14), ass. (0->171)
t367 = cos(pkin(6));
t352 = qJD(1) * t367 + qJD(2);
t361 = sin(pkin(12));
t363 = sin(pkin(6));
t370 = sin(qJ(3));
t373 = cos(qJ(3));
t365 = cos(pkin(12));
t366 = cos(pkin(7));
t445 = t365 * t366;
t380 = (t361 * t373 + t370 * t445) * t363;
t362 = sin(pkin(7));
t449 = t362 * t370;
t294 = qJD(1) * t380 + t352 * t449;
t369 = sin(qJ(4));
t372 = cos(qJ(4));
t404 = pkin(4) * t369 - qJ(5) * t372;
t322 = t404 * qJD(4) - qJD(5) * t369;
t464 = t294 - t322;
t448 = t362 * t373;
t461 = (-t361 * t370 + t373 * t445) * t363;
t463 = t367 * t448 + t461;
t462 = qJD(1) * t461 + t352 * t448;
t460 = (t369 ^ 2 - t372 ^ 2) * MDP(7);
t436 = qJD(1) * t363;
t416 = t365 * t436;
t407 = t366 * t416;
t417 = t361 * t436;
t292 = -t370 * t417 + (t352 * t362 + t407) * t373;
t360 = sin(pkin(13));
t429 = qJD(4) * t369;
t421 = pkin(9) * t429;
t351 = t360 * t421;
t364 = cos(pkin(13));
t450 = t360 * t372;
t443 = -t292 * t450 + t364 * t464 - t351;
t291 = qJD(3) * pkin(9) + t294;
t312 = t352 * t366 - t362 * t416;
t459 = -t369 * t291 + t312 * t372;
t433 = qJD(3) * t370;
t414 = t362 * t433;
t431 = qJD(3) * t373;
t287 = t352 * t414 + t407 * t433 + t417 * t431;
t458 = qJD(3) * t294 - t287;
t446 = t364 * t372;
t457 = -t292 * t446 - t360 * t464;
t456 = pkin(10) + qJ(5);
t455 = qJD(3) * pkin(3);
t286 = t462 * qJD(3);
t428 = qJD(4) * t372;
t250 = t369 * t286 + t291 * t428 + t312 * t429;
t454 = t250 * t360;
t453 = t250 * t364;
t368 = sin(qJ(6));
t451 = t360 * t368;
t447 = t364 * t369;
t249 = t286 * t372 + (qJD(5) + t459) * qJD(4);
t271 = t322 * qJD(3) + t287;
t240 = t364 * t249 + t360 * t271;
t394 = pkin(5) * t369 - pkin(10) * t446;
t385 = t394 * qJD(4);
t444 = -t385 + t443;
t442 = (-pkin(9) * t447 - pkin(10) * t450) * qJD(4) + t457;
t408 = t364 * t421;
t441 = t408 - t457;
t268 = t372 * t291 + t369 * t312;
t266 = qJD(4) * qJ(5) + t268;
t347 = -pkin(4) * t372 - qJ(5) * t369 - pkin(3);
t278 = t347 * qJD(3) - t292;
t246 = t364 * t266 + t360 * t278;
t342 = t404 * qJD(3);
t258 = t360 * t342 + t364 * t459;
t371 = cos(qJ(6));
t340 = -t371 * t364 + t451;
t386 = t372 * t340;
t440 = qJD(3) * t386 - t340 * qJD(6);
t341 = t360 * t371 + t364 * t368;
t387 = t372 * t341;
t439 = -qJD(3) * t387 + t341 * qJD(6);
t424 = t364 * qJD(4);
t434 = qJD(3) * t369;
t336 = t360 * t434 - t424;
t422 = qJD(3) * qJD(4);
t411 = t372 * t422;
t406 = t364 * t411;
t425 = qJD(6) * t371;
t438 = -t336 * t425 + t371 * t406;
t314 = pkin(9) * t446 + t360 * t347;
t432 = qJD(3) * t372;
t430 = qJD(4) * t360;
t242 = -pkin(10) * t336 + t246;
t427 = qJD(6) * t242;
t426 = qJD(6) * t369;
t423 = MDP(12) * qJD(4);
t418 = pkin(5) * t360 + pkin(9);
t415 = t360 * t432;
t413 = t362 * t431;
t412 = MDP(21) * t429;
t410 = -qJD(4) * pkin(4) + qJD(5);
t239 = -t249 * t360 + t364 * t271;
t245 = -t266 * t360 + t364 * t278;
t257 = t364 * t342 - t360 * t459;
t405 = t360 * t411;
t238 = -pkin(10) * t405 + t240;
t338 = t364 * t434 + t430;
t241 = -pkin(5) * t432 - pkin(10) * t338 + t245;
t409 = -qJD(6) * t241 - t238;
t234 = t241 * t371 - t242 * t368;
t235 = t241 * t368 + t242 * t371;
t302 = t367 * t449 + t380;
t323 = -t362 * t363 * t365 + t366 * t367;
t277 = t302 * t372 + t323 * t369;
t255 = -t277 * t360 - t364 * t463;
t256 = t277 * t364 - t360 * t463;
t403 = t255 * t371 - t256 * t368;
t402 = t255 * t368 + t256 * t371;
t335 = t364 * t347;
t300 = -pkin(10) * t447 + t335 + (-pkin(9) * t360 - pkin(5)) * t372;
t309 = -pkin(10) * t360 * t369 + t314;
t400 = t300 * t371 - t309 * t368;
t399 = t300 * t368 + t309 * t371;
t276 = t302 * t369 - t323 * t372;
t328 = t366 * t369 + t372 * t449;
t305 = -t328 * t360 - t364 * t448;
t306 = t328 * t364 - t360 * t448;
t398 = t305 * t371 - t306 * t368;
t397 = t305 * t368 + t306 * t371;
t396 = t336 * t368 - t338 * t371;
t395 = t360 * MDP(13) + t364 * MDP(14);
t374 = qJD(4) ^ 2;
t393 = pkin(9) * t374 - t458;
t290 = -t292 - t455;
t392 = qJD(4) * (t290 + t292 - t455);
t327 = -t366 * t372 + t369 * t449;
t390 = -t372 * MDP(11) + t369 * MDP(12) - MDP(4);
t350 = t456 * t364;
t389 = t394 * qJD(3) + qJD(5) * t360 + qJD(6) * t350 + t257;
t349 = t456 * t360;
t388 = pkin(10) * t415 + qJD(5) * t364 - qJD(6) * t349 - t258;
t264 = -t459 + t410;
t383 = -qJD(6) * t338 - t405;
t379 = qJD(4) * t387;
t377 = -qJ(5) * t429 + (-t264 + t410) * t372;
t324 = t371 * t336;
t297 = t338 * t368 + t324;
t376 = -MDP(11) * qJD(4) + t336 * MDP(13) + t338 * MDP(14) + t264 * MDP(16) + t297 * MDP(22) - MDP(23) * t396;
t273 = qJD(3) * t379 - t396 * qJD(6);
t375 = qJD(3) ^ 2;
t356 = -pkin(5) * t364 - pkin(4);
t354 = -qJD(6) + t432;
t343 = t418 * t369;
t332 = t418 * t428;
t320 = t340 * t369;
t319 = t341 * t369;
t313 = -pkin(9) * t450 + t335;
t307 = -t327 * qJD(4) + t372 * t413;
t296 = t302 * qJD(3);
t295 = t463 * qJD(3);
t289 = t425 * t447 - t426 * t451 + t379;
t288 = -qJD(4) * t386 - t341 * t426;
t285 = t307 * t364 + t360 * t414;
t284 = -t307 * t360 + t364 * t414;
t272 = t383 * t368 + t438;
t262 = pkin(5) * t415 + t268;
t259 = pkin(5) * t336 + t264;
t254 = -t276 * qJD(4) + t295 * t372;
t247 = pkin(5) * t405 + t250;
t244 = t254 * t364 + t296 * t360;
t243 = -t254 * t360 + t296 * t364;
t237 = qJD(3) * t385 + t239;
t236 = t371 * t237;
t1 = [-t254 * t423 + (-t243 * t338 - t244 * t336) * MDP(15) + (t239 * t255 + t240 * t256 + t243 * t245 + t244 * t246 + t250 * t276) * MDP(16) + (-(-t402 * qJD(6) + t243 * t371 - t244 * t368) * t354 + t276 * t273) * MDP(22) + ((t403 * qJD(6) + t243 * t368 + t244 * t371) * t354 + t276 * t272) * MDP(23) + t376 * (t277 * qJD(4) + t295 * t369) + (-t295 * MDP(5) + (-MDP(13) * t243 + MDP(14) * t244) * t372 + t390 * t296 + ((-t463 * MDP(12) + (-t255 * t364 - t256 * t360) * MDP(15) + t395 * t276) * t372 + (-MDP(11) * t463 + t255 * MDP(13) - t256 * MDP(14) + t403 * MDP(22) - t402 * MDP(23)) * t369) * qJD(4)) * qJD(3); -t307 * t423 + (-t284 * t338 - t285 * t336) * MDP(15) + (t239 * t305 + t240 * t306 + t245 * t284 + t246 * t285 + t250 * t327) * MDP(16) + (-(-t397 * qJD(6) + t284 * t371 - t285 * t368) * t354 + t327 * t273) * MDP(22) + ((t398 * qJD(6) + t284 * t368 + t285 * t371) * t354 + t327 * t272) * MDP(23) + (-MDP(5) * t373 + t390 * t370) * t375 * t362 + t376 * (t328 * qJD(4) + t369 * t413) + ((-MDP(13) * t284 + MDP(14) * t285) * t372 + ((-MDP(12) * t448 + (-t305 * t364 - t306 * t360) * MDP(15) + t395 * t327) * t372 + (-MDP(11) * t448 + t305 * MDP(13) - t306 * MDP(14) + t398 * MDP(22) - t397 * MDP(23)) * t369) * qJD(4)) * qJD(3); t458 * MDP(4) + (t292 - t462) * qJD(3) * MDP(5) - 0.2e1 * t422 * t460 + (t443 * t338 + t441 * t336 + (-t245 * t364 - t246 * t360 + (-t313 * t364 - t314 * t360) * qJD(3)) * t428) * MDP(15) + (t264 * t428 * pkin(9) + t239 * t313 + t240 * t314 - t443 * t245 - t441 * t246) * MDP(16) + (-t272 * t320 - t288 * t396) * MDP(17) + (-t272 * t319 + t273 * t320 - t288 * t297 + t289 * t396) * MDP(18) - t432 * t412 + (t247 * t319 + t259 * t289 + t343 * t273 + t332 * t297) * MDP(22) + (-t247 * t320 + t259 * t288 + t343 * t272 - t332 * t396) * MDP(23) + ((-qJD(3) * t320 - t396) * MDP(19) + (-qJD(3) * t319 - t297) * MDP(20)) * t429 + (-t288 * MDP(19) + t289 * MDP(20) - t412 + (t399 * qJD(6) + t442 * t368 + t444 * t371) * MDP(22) + (t400 * qJD(6) - t444 * t368 + t442 * t371) * MDP(23)) * t354 + (t374 * MDP(8) - t393 * MDP(11) + t392 * MDP(12) + (-t239 + (pkin(9) * t336 + t264 * t360) * qJD(4) + (t351 + t443) * qJD(3)) * MDP(13) + (t240 + (pkin(9) * t338 + t264 * t364) * qJD(4) + (t408 - t441) * qJD(3)) * MDP(14) - t272 * MDP(19) + t273 * MDP(20) + (t235 * qJD(6) + t238 * t368 - t236) * MDP(22) + (t234 * qJD(6) + t237 * t368 + t238 * t371) * MDP(23)) * t372 + (0.2e1 * MDP(6) * t411 - t374 * MDP(9) + t392 * MDP(11) + t393 * MDP(12) + (t454 - t292 * t336 + (qJD(3) * t313 + t245) * qJD(4)) * MDP(13) + (t453 - t292 * t338 + (-qJD(3) * t314 - t246) * qJD(4)) * MDP(14) + (-t239 * t364 - t240 * t360) * MDP(15) + (t250 * pkin(9) - t264 * t292) * MDP(16) + (-t292 * t297 + (t400 * qJD(3) + t234) * qJD(4)) * MDP(22) + (t292 * t396 + (-t399 * qJD(3) - t235) * qJD(4)) * MDP(23)) * t369; (qJD(4) * t268 - t290 * t434 - t250) * MDP(11) + (-qJD(3) * t290 - t286) * t372 * MDP(12) + (-t453 - t268 * t336 + (-t245 * t369 + t257 * t372 + t377 * t360) * qJD(3)) * MDP(13) + (t454 - t268 * t338 + (t246 * t369 - t258 * t372 + t377 * t364) * qJD(3)) * MDP(14) + (t257 * t338 + t258 * t336 + (-qJD(5) * t336 + t245 * t432 + t240) * t364 + (qJD(5) * t338 + t246 * t432 - t239) * t360) * MDP(15) + (-pkin(4) * t250 - t245 * t257 - t246 * t258 - t264 * t268 + (-t245 * t360 + t246 * t364) * qJD(5) + (-t239 * t360 + t240 * t364) * qJ(5)) * MDP(16) + (t272 * t341 - t396 * t440) * MDP(17) + (-t272 * t340 - t273 * t341 - t440 * t297 + t396 * t439) * MDP(18) + (-t440 * t354 + (qJD(4) * t341 + t396) * t434) * MDP(19) + (t439 * t354 + (-qJD(4) * t340 + t297) * t434) * MDP(20) + t354 * MDP(21) * t434 + (t247 * t340 - t262 * t297 + t356 * t273 + (t388 * t368 + t389 * t371) * t354 + t439 * t259 + ((-t349 * t371 - t350 * t368) * qJD(4) - t234) * t434) * MDP(22) + (t247 * t341 + t262 * t396 + t356 * t272 + (-t389 * t368 + t388 * t371) * t354 + t440 * t259 + (-(-t349 * t368 + t350 * t371) * qJD(4) + t235) * t434) * MDP(23) + (-t372 * t369 * MDP(6) + t460) * t375; (-t336 ^ 2 - t338 ^ 2) * MDP(15) + (t245 * t338 + t246 * t336 + t250) * MDP(16) + (t396 * t354 + t273) * MDP(22) + (t324 * t354 + (-t405 + (-qJD(6) + t354) * t338) * t368 + t438) * MDP(23) + ((-t338 + t430) * MDP(13) + (t336 + t424) * MDP(14)) * t432; -t297 ^ 2 * MDP(18) + (-t297 * t354 + t438) * MDP(19) + qJD(3) * t412 + (-t235 * t354 + t236) * MDP(22) + (-t234 * t354 + t259 * t297) * MDP(23) - (MDP(17) * t297 - MDP(18) * t396 - t354 * MDP(20) - t259 * MDP(22)) * t396 + (t383 * MDP(20) - MDP(22) * t427 + t409 * MDP(23)) * t371 + (t383 * MDP(19) + (qJD(6) * t336 - t406) * MDP(20) + t409 * MDP(22) + (-t237 + t427) * MDP(23)) * t368;];
tauc  = t1;
