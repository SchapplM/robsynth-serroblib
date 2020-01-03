% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRP9
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
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:06:14
% EndTime: 2019-12-31 22:06:24
% DurationCPUTime: 4.88s
% Computational Cost: add. (4025->429), mult. (9764->563), div. (0->0), fcn. (6413->6), ass. (0->181)
t410 = sin(qJ(3));
t411 = sin(qJ(2));
t477 = qJD(1) * t411;
t452 = t410 * t477;
t412 = cos(qJ(3));
t471 = qJD(2) * t412;
t358 = -t452 + t471;
t473 = qJD(2) * t410;
t359 = t412 * t477 + t473;
t409 = sin(qJ(4));
t512 = cos(qJ(4));
t317 = -t512 * t358 + t359 * t409;
t426 = t409 * t358 + t359 * t512;
t536 = t317 * t426;
t413 = cos(qJ(2));
t495 = t412 * t413;
t432 = pkin(3) * t411 - pkin(8) * t495;
t513 = -pkin(8) - pkin(7);
t458 = qJD(3) * t513;
t434 = pkin(2) * t411 - pkin(7) * t413;
t362 = t434 * qJD(1);
t480 = pkin(6) * t452 + t412 * t362;
t535 = qJD(1) * t432 - t412 * t458 + t480;
t348 = t410 * t362;
t498 = t411 * t412;
t499 = t410 * t413;
t534 = -t348 - (-pkin(6) * t498 - pkin(8) * t499) * qJD(1) + t410 * t458;
t468 = qJD(3) * t411;
t447 = qJD(1) * t468;
t462 = qJD(1) * qJD(2);
t448 = t413 * t462;
t530 = qJD(2) * qJD(3) + t448;
t481 = t530 * t412;
t420 = t410 * t447 - t481;
t450 = t512 * qJD(4);
t459 = t530 * t410 + t412 * t447;
t466 = qJD(4) * t409;
t282 = -t358 * t450 + t359 * t466 + t409 * t459 + t512 * t420;
t476 = qJD(1) * t413;
t394 = -qJD(3) + t476;
t380 = -qJD(4) + t394;
t269 = -t317 * t380 - t282;
t283 = qJD(4) * t426 - t409 * t420 + t512 * t459;
t472 = qJD(2) * t411;
t446 = MDP(22) * t472;
t514 = t426 ^ 2;
t533 = qJD(1) * t446 + (-t380 * t426 - t283) * MDP(21) + MDP(18) * t536 + (-t317 ^ 2 + t514) * MDP(19) + t269 * MDP(20);
t373 = -qJD(2) * pkin(2) + pkin(6) * t477;
t335 = -pkin(3) * t358 + t373;
t279 = pkin(4) * t317 - qJ(5) * t426 + t335;
t532 = t279 * t317;
t531 = t317 * t335;
t502 = t409 * t410;
t424 = t412 * t512 - t502;
t517 = qJD(3) + qJD(4);
t518 = t512 * qJD(3) + t450;
t486 = -t518 * t412 + t424 * t476 + t517 * t502;
t361 = t409 * t412 + t410 * t512;
t328 = t517 * t361;
t485 = -t361 * t476 + t328;
t467 = qJD(3) * t412;
t453 = t411 * t467;
t470 = qJD(2) * t413;
t456 = t410 * t470;
t529 = t453 + t456;
t292 = pkin(4) * t426 + qJ(5) * t317;
t526 = -0.2e1 * t462;
t524 = MDP(4) * t411;
t407 = t411 ^ 2;
t523 = MDP(5) * (-t413 ^ 2 + t407);
t508 = t279 * t426;
t522 = t335 * t426;
t375 = t513 * t410;
t376 = t513 * t412;
t425 = t375 * t512 + t409 * t376;
t521 = qJD(4) * t425 - t535 * t409 + t534 * t512;
t334 = t409 * t375 - t376 * t512;
t520 = qJD(4) * t334 + t534 * t409 + t535 * t512;
t368 = -pkin(2) * t413 - pkin(7) * t411 - pkin(1);
t357 = t412 * t368;
t510 = pkin(6) * t410;
t323 = -pkin(8) * t498 + t357 + (-pkin(3) - t510) * t413;
t396 = pkin(6) * t495;
t479 = t410 * t368 + t396;
t500 = t410 * t411;
t329 = -pkin(8) * t500 + t479;
t519 = t409 * t323 + t512 * t329;
t404 = pkin(6) * t476;
t469 = qJD(3) * t410;
t436 = -t404 + (-t410 * t476 + t469) * pkin(3);
t365 = t434 * qJD(2);
t482 = t412 * t365 + t472 * t510;
t293 = t432 * qJD(2) + (-t396 + (pkin(8) * t411 - t368) * t410) * qJD(3) + t482;
t454 = t413 * t469;
t483 = t410 * t365 + t368 * t467;
t295 = -t529 * pkin(8) + (-t411 * t471 - t454) * pkin(6) + t483;
t516 = -qJD(4) * t519 + t293 * t512 - t409 * t295;
t511 = pkin(6) * t394;
t509 = qJD(2) * pkin(7);
t507 = t358 * t394;
t506 = t373 * t410;
t505 = t373 * t412;
t504 = t394 * t412;
t374 = t404 + t509;
t496 = t412 * t374;
t352 = t368 * qJD(1);
t501 = t410 * t352;
t325 = t496 + t501;
t308 = pkin(8) * t358 + t325;
t503 = t409 * t308;
t414 = qJD(2) ^ 2;
t497 = t411 * t414;
t494 = t413 * t414;
t415 = qJD(1) ^ 2;
t493 = t413 * t415;
t324 = t412 * t352 - t374 * t410;
t307 = -pkin(8) * t359 + t324;
t275 = t307 * t512 - t503;
t492 = -pkin(3) * t450 - qJD(5) + t275;
t491 = t485 * pkin(4) + t486 * qJ(5) - qJD(5) * t361 + t436;
t490 = qJ(5) * t477 - t521;
t489 = pkin(4) * t477 + t520;
t353 = qJD(1) * t365;
t449 = t411 * t462;
t440 = pkin(6) * t449;
t484 = t412 * t353 + t410 * t440;
t366 = pkin(3) * t500 + t411 * pkin(6);
t475 = qJD(2) * t425;
t474 = qJD(2) * t334;
t464 = t411 * MDP(15);
t303 = -pkin(3) * t394 + t307;
t272 = t303 * t512 - t503;
t463 = qJD(5) - t272;
t336 = t529 * pkin(3) + pkin(6) * t470;
t402 = -pkin(3) * t412 - pkin(2);
t457 = t512 * t308;
t455 = t410 * t468;
t445 = qJD(2) * t464;
t444 = pkin(1) * t526;
t443 = -t358 + t471;
t442 = -t359 + t473;
t441 = pkin(4) * t449;
t278 = -t481 * pkin(8) + pkin(3) * t449 + (-t496 + (pkin(8) * t477 - t352) * t410) * qJD(3) + t484;
t431 = t352 * t467 + t410 * t353 - t374 * t469;
t417 = -t412 * t440 + t431;
t284 = -pkin(8) * t459 + t417;
t439 = t409 * t278 + t512 * t284 + t303 * t450 - t308 * t466;
t438 = -t512 * t278 + t409 * t284 + t303 * t466 + t308 * t450;
t437 = t512 * t470;
t274 = t409 * t307 + t457;
t435 = pkin(3) * t466 - t274;
t433 = qJD(1) * t407 - t394 * t413;
t315 = pkin(3) * t459 + pkin(6) * t448;
t367 = t380 * qJD(5);
t389 = qJ(5) * t449;
t263 = t389 - t367 + t439;
t429 = t323 * t512 - t409 * t329;
t273 = t409 * t303 + t457;
t423 = -t272 * t380 - t439;
t422 = -t273 * t380 - t438;
t421 = t409 * t293 + t512 * t295 + t323 * t450 - t329 * t466;
t419 = t420 * t411;
t264 = t438 - t441;
t401 = -pkin(3) * t512 - pkin(4);
t397 = pkin(3) * t409 + qJ(5);
t344 = t424 * t411;
t343 = t361 * t411;
t314 = -pkin(4) * t424 - qJ(5) * t361 + t402;
t304 = pkin(4) * t343 - qJ(5) * t344 + t366;
t297 = t410 * t437 - t409 * t455 - t466 * t500 + (t409 * t470 + t518 * t411) * t412;
t296 = t328 * t411 + t409 * t456 - t412 * t437;
t291 = t413 * pkin(4) - t429;
t290 = -qJ(5) * t413 + t519;
t286 = pkin(3) * t359 + t292;
t271 = -t380 * qJ(5) + t273;
t270 = t380 * pkin(4) + t463;
t268 = pkin(4) * t297 + qJ(5) * t296 - qJD(5) * t344 + t336;
t267 = -pkin(4) * t472 - t516;
t266 = qJ(5) * t472 - qJD(5) * t413 + t421;
t265 = t283 * pkin(4) + t282 * qJ(5) - qJD(5) * t426 + t315;
t1 = [t523 * t526 + (-pkin(6) * t494 + t411 * t444) * MDP(9) + (pkin(6) * t497 + t413 * t444) * MDP(10) + (-t412 * t419 + (t412 * t470 - t455) * t359) * MDP(11) + ((t358 * t412 - t359 * t410) * t470 + (-t412 * t459 - t481 * t410 + (-t359 * t412 + (-t358 + t452) * t410) * qJD(3)) * t411) * MDP(12) + (t394 * t455 + t420 * t413 + (t359 * t411 + t412 * t433) * qJD(2)) * MDP(13) + (t394 * t453 + t459 * t413 + (t358 * t411 - t410 * t433) * qJD(2)) * MDP(14) + (-t394 - t476) * t445 + (-(-t368 * t469 + t482) * t394 + (pkin(6) * t459 + t373 * t467 + (t357 * qJD(1) + t324) * qJD(2)) * t411 + ((-pkin(6) * t358 + t506) * qJD(2) + (t501 + (t374 + t511) * t412) * qJD(3) - t484) * t413) * MDP(16) + (-pkin(6) * t419 - t373 * t455 + (-pkin(6) * t454 + t483) * t394 + t431 * t413 + ((pkin(6) * t359 + t505) * t413 + (-pkin(6) * t504 - qJD(1) * t479 - t325) * t411) * qJD(2)) * MDP(17) + (-t282 * t344 - t296 * t426) * MDP(18) + (t282 * t343 - t283 * t344 + t296 * t317 - t297 * t426) * MDP(19) + (t282 * t413 + t296 * t380 + (qJD(1) * t344 + t426) * t472) * MDP(20) + (t283 * t413 + t297 * t380 + (-qJD(1) * t343 - t317) * t472) * MDP(21) + (-t380 - t476) * t446 + (t272 * t472 + t366 * t283 + t335 * t297 + t315 * t343 + t336 * t317 - t516 * t380 + t438 * t413 + t429 * t449) * MDP(23) + (t421 * t380 + t439 * t413 + t336 * t426 - t366 * t282 + t315 * t344 - t335 * t296 + (-qJD(1) * t519 - t273) * t472) * MDP(24) + (t264 * t413 + t265 * t343 + t267 * t380 + t268 * t317 + t279 * t297 + t283 * t304 + (-qJD(1) * t291 - t270) * t472) * MDP(25) + (-t263 * t343 + t264 * t344 - t266 * t317 + t267 * t426 - t270 * t296 - t271 * t297 - t282 * t291 - t283 * t290) * MDP(26) + (-t263 * t413 - t265 * t344 - t266 * t380 - t268 * t426 + t279 * t296 + t282 * t304 + (qJD(1) * t290 + t271) * t472) * MDP(27) + (t263 * t290 + t264 * t291 + t265 * t304 + t266 * t271 + t267 * t270 + t268 * t279) * MDP(28) + 0.2e1 * t448 * t524 - MDP(7) * t497 + MDP(6) * t494; -t493 * t524 + t415 * t523 + (-t359 * t504 - t410 * t420) * MDP(11) + ((t481 - t507) * t412 + (-t359 * qJD(3) + (t413 * t359 - t453) * qJD(1) - t459) * t410) * MDP(12) + (-t394 * t467 + (t394 * t495 + t411 * t442) * qJD(1)) * MDP(13) + (t394 * t469 + (-t394 * t499 + t411 * t443) * qJD(1)) * MDP(14) + t394 * qJD(1) * t464 + (-pkin(2) * t459 + t480 * t394 + (pkin(7) * t504 + t506) * qJD(3) + ((-pkin(7) * t473 - t324) * t411 + (-pkin(6) * t443 - t506) * t413) * qJD(1)) * MDP(16) + (-pkin(2) * t481 - t348 * t394 + (-pkin(7) * t394 * t410 + t505) * qJD(3) + ((pkin(6) * t442 - t505) * t413 + (pkin(2) * t469 + t325 + (-t509 + t511) * t412) * t411) * qJD(1)) * MDP(17) + (-t282 * t361 - t426 * t486) * MDP(18) + (-t282 * t424 - t283 * t361 + t317 * t486 - t426 * t485) * MDP(19) + (t486 * t380 + (qJD(2) * t361 - t426) * t477) * MDP(20) + (t485 * t380 + (qJD(2) * t424 + t317) * t477) * MDP(21) + t380 * MDP(22) * t477 + (t402 * t283 - t315 * t424 + t520 * t380 + t485 * t335 + t436 * t317 + (-t272 + t475) * t477) * MDP(23) + (-t402 * t282 + t315 * t361 + t521 * t380 - t486 * t335 + t436 * t426 + (t273 - t474) * t477) * MDP(24) + (-t265 * t424 + t283 * t314 + t489 * t380 + t491 * t317 + t485 * t279 + (t270 + t475) * t477) * MDP(25) + (t263 * t424 + t264 * t361 - t270 * t486 - t271 * t485 + t282 * t425 - t283 * t334 + t317 * t490 + t426 * t489) * MDP(26) + (-t265 * t361 + t282 * t314 + t490 * t380 - t491 * t426 + t486 * t279 + (-t271 + t474) * t477) * MDP(27) + (t263 * t334 - t264 * t425 + t265 * t314 + t270 * t489 - t271 * t490 + t279 * t491) * MDP(28) + (MDP(9) * t411 * t415 + MDP(10) * t493) * pkin(1); -t359 * t358 * MDP(11) + (-t358 ^ 2 + t359 ^ 2) * MDP(12) + (-t420 + t507) * MDP(13) + (-t359 * t394 - t459) * MDP(14) + qJD(1) * t445 + (-t359 * t373 + t484 + (-qJD(3) - t394) * t325) * MDP(16) + (-t324 * t394 - t358 * t373 - t417) * MDP(17) + (-t274 * t380 - t522 + (-t317 * t359 + t380 * t466 + t449 * t512) * pkin(3) - t438) * MDP(23) + (-t275 * t380 + t531 + (-t359 * t426 + t380 * t450 - t409 * t449) * pkin(3) - t439) * MDP(24) + (-t508 - t286 * t317 + t435 * t380 + (pkin(4) - t401) * t449 - t438) * MDP(25) + (-t282 * t401 - t283 * t397 + (t271 + t435) * t426 + (t270 + t492) * t317) * MDP(26) + (t286 * t426 + t380 * t492 + t397 * t449 + t263 - t532) * MDP(27) + (t263 * t397 + t264 * t401 + t270 * t435 - t271 * t492 - t279 * t286) * MDP(28) + t533; (t422 - t522) * MDP(23) + (t423 + t531) * MDP(24) + (-t292 * t317 + t422 + 0.2e1 * t441 - t508) * MDP(25) + (pkin(4) * t282 - qJ(5) * t283 + (t271 - t273) * t426 + (t270 - t463) * t317) * MDP(26) + (t292 * t426 - 0.2e1 * t367 + 0.2e1 * t389 - t423 - t532) * MDP(27) + (-pkin(4) * t264 + qJ(5) * t263 - t270 * t273 + t271 * t463 - t279 * t292) * MDP(28) + t533; (-t449 + t536) * MDP(25) + t269 * MDP(26) + (-t380 ^ 2 - t514) * MDP(27) + (t271 * t380 + t264 + t508) * MDP(28);];
tauc = t1;
