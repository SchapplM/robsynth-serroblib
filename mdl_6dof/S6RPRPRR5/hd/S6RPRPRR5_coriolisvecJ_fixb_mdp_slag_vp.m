% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRPRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:49:42
% EndTime: 2019-03-09 03:49:51
% DurationCPUTime: 4.92s
% Computational Cost: add. (3829->362), mult. (9831->465), div. (0->0), fcn. (7662->8), ass. (0->155)
t415 = sin(qJ(6));
t464 = qJD(6) * t415;
t409 = qJD(3) - qJD(5);
t418 = cos(qJ(6));
t413 = cos(pkin(10));
t498 = cos(qJ(3));
t458 = t498 * t413;
t439 = qJD(1) * t458;
t397 = qJD(3) * t439;
t412 = sin(pkin(10));
t417 = sin(qJ(3));
t467 = qJD(3) * t417;
t456 = t412 * t467;
t352 = qJD(1) * t456 - t397;
t386 = t412 * t498 + t417 * t413;
t377 = t386 * qJD(3);
t353 = qJD(1) * t377;
t480 = t412 * t417;
t457 = qJD(1) * t480;
t372 = -t439 + t457;
t374 = t386 * qJD(1);
t416 = sin(qJ(5));
t419 = cos(qJ(5));
t465 = qJD(5) * t419;
t466 = qJD(5) * t416;
t440 = t419 * t352 - t416 * t353 - t372 * t465 + t374 * t466;
t463 = qJD(6) * t418;
t474 = -t409 * t463 - t418 * t440;
t505 = t372 * t416 + t419 * t374;
t273 = -t464 * t505 + t474;
t317 = -t409 * t415 + t418 * t505;
t492 = t440 * t415;
t274 = t317 * qJD(6) - t492;
t287 = qJD(5) * t505 - t352 * t416 - t419 * t353;
t477 = t418 * t287;
t506 = -t419 * t372 + t374 * t416;
t511 = qJD(6) + t506;
t450 = t464 * t511 - t477;
t436 = t415 * t506 * t511 + t450;
t485 = t506 * t409;
t486 = t505 * t409;
t488 = t317 * t505;
t489 = t317 * t511;
t484 = t505 * t415;
t315 = t418 * t409 + t484;
t490 = t315 * t505;
t491 = t315 * t511;
t495 = t273 * t415;
t479 = t415 * t287;
t521 = t511 * t418;
t517 = t511 * t521 + t479;
t529 = -((t274 + t489) * t415 - (t273 - t491) * t418) * MDP(27) + (t317 * t521 + t495) * MDP(26) + (-t488 + t517) * MDP(28) - (t440 + t485) * MDP(21) - (t287 + t486) * MDP(22) - (t436 - t490) * MDP(29) + (t505 ^ 2 - t506 ^ 2) * MDP(20) + (MDP(19) * t506 - t511 * MDP(30)) * t505;
t497 = pkin(7) + qJ(2);
t393 = t497 * t413;
t388 = qJD(1) * t393;
t369 = t417 * t388;
t392 = t497 * t412;
t387 = qJD(1) * t392;
t338 = -t498 * t387 - t369;
t461 = qJD(4) - t338;
t515 = pkin(5) * t505;
t404 = -t413 * pkin(2) - pkin(1);
t391 = t404 * qJD(1) + qJD(2);
t318 = t372 * pkin(3) - t374 * qJ(4) + t391;
t299 = -pkin(4) * t372 - t318;
t272 = pkin(5) * t506 - pkin(9) * t505 + t299;
t420 = -pkin(3) - pkin(4);
t512 = -pkin(8) * t374 + t461;
t306 = qJD(3) * t420 + t512;
t339 = -t417 * t387 + t498 * t388;
t314 = pkin(8) * t372 + t339;
t411 = qJD(3) * qJ(4);
t309 = t314 + t411;
t280 = t306 * t416 + t309 * t419;
t277 = -pkin(9) * t409 + t280;
t265 = t272 * t418 - t277 * t415;
t514 = t265 * t505;
t266 = t272 * t415 + t277 * t418;
t513 = t266 * t505;
t410 = qJD(3) * qJD(4);
t459 = qJD(1) * qJD(2);
t455 = qJD(2) * t498;
t437 = qJD(1) * t455;
t454 = qJD(3) * t498;
t472 = -t387 * t454 + t413 * t437;
t307 = t410 + (-qJD(3) * t388 - t412 * t459) * t417 + t472;
t293 = pkin(8) * t353 + t307;
t453 = t417 * t459;
t308 = -t387 * t467 + t388 * t454 + t412 * t437 + t413 * t453;
t297 = pkin(8) * t352 + t308;
t263 = t419 * t293 + t416 * t297 + t306 * t465 - t309 * t466;
t509 = t299 * t506 - t263;
t264 = t416 * t293 - t297 * t419 + t306 * t466 + t309 * t465;
t508 = -t299 * t505 - t264;
t341 = -t417 * t392 + t498 * t393;
t504 = (t412 ^ 2 + t413 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t334 = t374 * pkin(3) + t372 * qJ(4);
t310 = -pkin(4) * t374 - t334;
t435 = qJ(4) * t419 + t416 * t420;
t390 = -pkin(9) + t435;
t503 = t511 * (-pkin(9) * t506 + qJD(6) * t390 + t310 - t515) - t264;
t502 = t511 * (t511 * pkin(9) + t515) + t264;
t501 = qJD(3) * (t338 + t369) + t412 * t453 - t472;
t319 = -t392 * t454 + t413 * t455 + (-qJD(2) * t412 - qJD(3) * t393) * t417;
t300 = pkin(8) * t377 + t319;
t320 = t386 * qJD(2) + qJD(3) * t341;
t376 = -t413 * t454 + t456;
t301 = t376 * pkin(8) + t320;
t340 = t392 * t498 + t417 * t393;
t322 = -t386 * pkin(8) + t340;
t385 = -t458 + t480;
t323 = pkin(8) * t385 + t341;
t433 = t322 * t419 - t323 * t416;
t269 = qJD(5) * t433 + t300 * t419 + t301 * t416;
t279 = t306 * t419 - t309 * t416;
t276 = pkin(5) * t409 - t279;
t335 = t385 * pkin(3) - t386 * qJ(4) + t404;
t312 = -pkin(4) * t385 - t335;
t337 = t385 * t416 + t386 * t419;
t431 = t419 * t385 - t386 * t416;
t278 = -pkin(5) * t431 - pkin(9) * t337 + t312;
t284 = t322 * t416 + t323 * t419;
t295 = qJD(5) * t431 - t376 * t419 + t377 * t416;
t500 = t264 * t337 + t276 * t295 - t284 * t287 - (qJD(6) * t278 + t269) * t511 + (qJD(6) * t272 + t263) * t431;
t499 = t374 ^ 2;
t494 = t276 * t337;
t493 = t278 * t287;
t487 = t318 * t374;
t483 = t372 * t374;
t434 = -qJ(4) * t416 + t419 * t420;
t475 = -qJD(5) * t434 + t314 * t416 - t512 * t419;
t473 = qJD(5) * t435 + t314 * t419 + t512 * t416;
t469 = qJD(3) * t319;
t468 = qJD(3) * t320;
t304 = t353 * pkin(3) + t352 * qJ(4) - t374 * qJD(4);
t311 = t377 * pkin(3) + t376 * qJ(4) - t386 * qJD(4);
t442 = t409 ^ 2;
t441 = t419 * t409;
t292 = -pkin(4) * t353 - t304;
t298 = -pkin(4) * t377 - t311;
t429 = t295 * t418 - t337 * t464;
t428 = qJD(3) * t339 - t308;
t425 = -pkin(9) * t287 + (t276 + t279) * t511;
t423 = -t390 * t287 + (-t276 + t475) * t511;
t389 = pkin(5) - t434;
t366 = t372 ^ 2;
t333 = t411 + t339;
t332 = t397 + (t372 - t457) * qJD(3);
t330 = -qJD(3) * pkin(3) + t461;
t296 = qJD(5) * t337 - t376 * t416 - t419 * t377;
t271 = pkin(5) * t296 - pkin(9) * t295 + t298;
t270 = qJD(5) * t284 + t300 * t416 - t301 * t419;
t268 = pkin(5) * t287 + pkin(9) * t440 + t292;
t267 = t418 * t268;
t1 = [(-t352 * t386 - t374 * t376) * MDP(8) + (-t307 * t385 + t308 * t386 - t319 * t372 + t320 * t374 - t330 * t376 - t333 * t377 - t340 * t352 - t341 * t353) * MDP(16) + (t352 * t385 - t353 * t386 + t372 * t376 - t374 * t377) * MDP(9) + (t273 * t337 * t418 + t317 * t429) * MDP(26) + (-t287 * t431 + t296 * t511) * MDP(30) + (-t287 * t337 - t295 * t506 - t296 * t505 - t431 * t440) * MDP(20) + (t295 * t505 - t337 * t440) * MDP(19) + (t304 * t335 + t307 * t341 + t308 * t340 + t311 * t318 + t319 * t333 + t320 * t330) * MDP(18) + (t287 * t312 - t292 * t431 + t296 * t299 + t298 * t506) * MDP(24) + (t292 * t337 + t295 * t299 + t298 * t505 - t312 * t440) * MDP(25) + (t353 * t404 + t377 * t391 - t468) * MDP(13) + (t304 * t385 + t311 * t372 + t318 * t377 + t335 * t353 - t468) * MDP(15) + (-t352 * t404 - t376 * t391 - t469) * MDP(14) + (-t304 * t386 - t311 * t374 + t318 * t376 + t335 * t352 + t469) * MDP(17) + (-t273 * t431 + t296 * t317 + t337 * t477 + t429 * t511) * MDP(28) + (-t337 * t479 + t274 * t431 - t296 * t315 + (-t295 * t415 - t337 * t463) * t511) * MDP(29) + (-t266 * t296 + t270 * t317 - t433 * t273 + (-(-qJD(6) * t284 + t271) * t511 - t493 + (-qJD(6) * t277 + t268) * t431 - qJD(6) * t494) * t415 + t500 * t418) * MDP(32) + (t265 * t296 - t267 * t431 + t270 * t315 - t433 * t274 + (t271 * t511 + t493 + (t277 * t431 - t284 * t511 + t494) * qJD(6)) * t418 + t500 * t415) * MDP(31) + ((-t315 * t418 - t317 * t415) * t295 + (-t495 - t274 * t418 + (t315 * t415 - t317 * t418) * qJD(6)) * t337) * MDP(27) + 0.2e1 * t459 * t504 + (-MDP(21) * t295 + MDP(22) * t296 + MDP(24) * t270 + MDP(25) * t269) * t409 + (-MDP(10) * t376 - MDP(11) * t377) * qJD(3); (-t366 - t499) * MDP(16) + (-t330 * t374 + t333 * t372 + t304) * MDP(18) + (-t287 + t486) * MDP(24) + (t440 - t485) * MDP(25) + (t436 + t490) * MDP(31) + (t488 + t517) * MDP(32) - qJD(1) ^ 2 * t504 + 0.2e1 * (MDP(13) + MDP(15)) * qJD(3) * t374 + (-MDP(14) + MDP(17)) * (-t397 + (t372 + t457) * qJD(3)); (-t310 * t506 + t409 * t473 - t508) * MDP(24) + (-t310 * t505 - t409 * t475 - t509) * MDP(25) + (-t318 * t372 + t334 * t374 + 0.2e1 * t410 - t501) * MDP(17) + (t372 * t391 + t501) * MDP(14) + (-t334 * t372 + t428 - t487) * MDP(15) + (-t374 * t391 + t428) * MDP(13) + (pkin(3) * t352 - qJ(4) * t353 + (t333 - t339) * t374 + (t330 - t461) * t372) * MDP(16) + t332 * MDP(10) + (t389 * t273 + t473 * t317 + t415 * t503 + t423 * t418 - t513) * MDP(32) + (t389 * t274 + t473 * t315 + t423 * t415 - t418 * t503 + t514) * MDP(31) + (-pkin(3) * t308 + qJ(4) * t307 - t318 * t334 - t330 * t339 + t333 * t461) * MDP(18) + MDP(8) * t483 + (-t366 + t499) * MDP(9) - t529; MDP(15) * t483 + t332 * MDP(16) + (-qJD(3) ^ 2 - t499) * MDP(17) + (-qJD(3) * t333 + t308 + t487) * MDP(18) + (-t374 * t506 - t416 * t442) * MDP(24) + (-t374 * t505 - t419 * t442) * MDP(25) + (-t419 * t274 + (-t418 * t374 + t415 * t441) * t511 + (-t315 * t409 - t463 * t511 - t479) * t416) * MDP(31) + (-t419 * t273 + (t415 * t374 + t418 * t441) * t511 + (-t317 * t409 + t450) * t416) * MDP(32); (-t280 * t409 + t508) * MDP(24) + (-t279 * t409 + t509) * MDP(25) + (-pkin(5) * t274 - t280 * t315 + t425 * t415 - t418 * t502 - t514) * MDP(31) + (-pkin(5) * t273 - t280 * t317 + t415 * t502 + t425 * t418 + t513) * MDP(32) + t529; t317 * t315 * MDP(26) + (-t315 ^ 2 + t317 ^ 2) * MDP(27) + (t474 + t491) * MDP(28) + (t489 + t492) * MDP(29) + t287 * MDP(30) + (-t263 * t415 + t266 * t511 - t276 * t317 + t267) * MDP(31) + (-t263 * t418 + t265 * t511 - t268 * t415 + t276 * t315) * MDP(32) + (-MDP(28) * t484 - MDP(29) * t317 - MDP(31) * t266 - MDP(32) * t265) * qJD(6);];
tauc  = t1;
