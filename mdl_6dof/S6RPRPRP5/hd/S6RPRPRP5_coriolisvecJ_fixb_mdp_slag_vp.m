% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:16:30
% EndTime: 2019-03-09 03:16:39
% DurationCPUTime: 5.92s
% Computational Cost: add. (6572->424), mult. (17069->541), div. (0->0), fcn. (13185->8), ass. (0->162)
t440 = sin(pkin(10));
t442 = cos(pkin(10));
t441 = sin(pkin(9));
t443 = cos(pkin(9));
t445 = sin(qJ(3));
t521 = cos(qJ(3));
t419 = t521 * t441 + t445 * t443;
t526 = t419 * qJD(1);
t387 = qJD(3) * t440 + t442 * t526;
t444 = sin(qJ(5));
t446 = cos(qJ(5));
t395 = t440 * t526;
t529 = qJD(3) * t442 - t395;
t461 = t446 * t529;
t341 = -t387 * t444 + t461;
t546 = t341 ^ 2;
t478 = t521 * t443;
t432 = qJD(1) * t478;
t500 = t445 * t441;
t477 = qJD(1) * t500;
t406 = -t432 + t477;
t401 = qJD(5) + t406;
t545 = t341 * t401;
t544 = MDP(24) + MDP(26);
t543 = -MDP(25) + MDP(28);
t534 = t446 * t387 + t444 * t529;
t522 = t534 ^ 2;
t412 = t419 * qJD(3);
t397 = qJD(1) * t412;
t498 = t446 * t442;
t416 = t440 * t444 - t498;
t418 = t440 * t446 + t442 * t444;
t410 = t418 * qJD(5);
t490 = t418 * t406 + t410;
t469 = -t416 * t397 - t401 * t490;
t485 = qJD(5) * t446;
t486 = qJD(5) * t444;
t528 = -t440 * t486 + t442 * t485;
t491 = -t416 * t406 + t528;
t468 = t418 * t397 + t401 * t491;
t536 = t419 * qJD(2);
t535 = qJD(3) * t526;
t429 = qJD(3) * t432;
t396 = -qJD(3) * t477 + t429;
t533 = t396 * t418;
t437 = -pkin(2) * t443 - pkin(1);
t457 = t478 - t500;
t372 = -pkin(3) * t457 - qJ(4) * t419 + t437;
t518 = pkin(7) + qJ(2);
t424 = t518 * t441;
t426 = t518 * t443;
t385 = -t445 * t424 + t426 * t521;
t328 = t442 * t372 - t385 * t440;
t519 = pkin(8) * t442;
t315 = -pkin(4) * t457 - t419 * t519 + t328;
t329 = t440 * t372 + t442 * t385;
t503 = t419 * t440;
t321 = -pkin(8) * t503 + t329;
t532 = t444 * t315 + t446 * t321;
t370 = pkin(3) * t526 + qJ(4) * t406;
t420 = qJD(1) * t424;
t421 = qJD(1) * t426;
t458 = t521 * t420 + t445 * t421;
t325 = t442 * t370 + t440 * t458;
t308 = pkin(4) * t526 + t406 * t519 + t325;
t326 = t440 * t370 - t442 * t458;
t505 = t406 * t440;
t317 = pkin(8) * t505 + t326;
t517 = pkin(8) + qJ(4);
t423 = t517 * t440;
t425 = t517 * t442;
t462 = -t423 * t446 - t425 * t444;
t531 = qJD(4) * t416 - qJD(5) * t462 + t444 * t308 + t446 * t317;
t384 = -t423 * t444 + t425 * t446;
t530 = -qJD(4) * t418 - qJD(5) * t384 - t308 * t446 + t317 * t444;
t527 = -t521 * t424 - t445 * t426;
t525 = (t441 ^ 2 + t443 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t506 = t396 * t440;
t309 = t444 * (qJD(5) * t387 + t506) - qJD(5) * t461 - t396 * t498;
t524 = t309 * t416 - t490 * t534;
t411 = t457 * qJD(3);
t348 = pkin(3) * t412 - qJ(4) * t411 - qJD(4) * t419;
t451 = t457 * qJD(2);
t355 = qJD(3) * t527 + t451;
t313 = t442 * t348 - t355 * t440;
t294 = pkin(4) * t412 - t411 * t519 + t313;
t314 = t440 * t348 + t442 * t355;
t504 = t411 * t440;
t300 = -pkin(8) * t504 + t314;
t523 = -qJD(5) * t532 + t294 * t446 - t300 * t444;
t402 = t406 ^ 2;
t520 = pkin(5) * t397;
t515 = qJ(6) * t397;
t422 = qJD(1) * t437 + qJD(2);
t352 = pkin(3) * t406 - qJ(4) * t526 + t422;
t376 = -t445 * t420 + t521 * t421;
t369 = qJD(3) * qJ(4) + t376;
t319 = t442 * t352 - t369 * t440;
t297 = pkin(4) * t406 - pkin(8) * t387 + t319;
t320 = t440 * t352 + t442 * t369;
t307 = pkin(8) * t529 + t320;
t279 = t297 * t444 + t307 * t446;
t514 = t279 * t401;
t512 = t341 * t526;
t511 = t534 * t341;
t510 = t534 * t526;
t508 = t462 * t397;
t507 = t384 * t397;
t501 = t442 * t396;
t497 = qJ(6) * t526 + t531;
t496 = -pkin(5) * t526 + t530;
t346 = -pkin(4) * t505 + t376;
t495 = -t490 * pkin(5) + t491 * qJ(6) + qJD(6) * t418 + t346;
t335 = pkin(3) * t397 - qJ(4) * t396 - qJD(4) * t526;
t342 = qJD(1) * t451 + (qJD(4) - t458) * qJD(3);
t302 = t440 * t335 + t442 * t342;
t487 = qJD(3) * t445;
t278 = t297 * t446 - t307 * t444;
t484 = qJD(6) - t278;
t482 = qJD(1) * qJD(2);
t436 = -pkin(4) * t442 - pkin(3);
t474 = qJD(3) * t521;
t301 = t442 * t335 - t342 * t440;
t288 = pkin(4) * t397 - pkin(8) * t501 + t301;
t291 = -pkin(8) * t506 + t302;
t471 = -t446 * t288 + t444 * t291 + t297 * t486 + t307 * t485;
t345 = t536 * qJD(1) - t420 * t487 + t421 * t474;
t356 = -t424 * t487 + t426 * t474 + t536;
t310 = qJD(5) * t534 + t533;
t470 = -t418 * t310 + t341 * t491;
t324 = pkin(4) * t506 + t345;
t330 = pkin(4) * t504 + t356;
t465 = t315 * t446 - t321 * t444;
t463 = -t319 * t440 + t320 * t442;
t357 = pkin(4) * t503 - t527;
t367 = -qJD(3) * pkin(3) + qJD(4) + t458;
t331 = -pkin(4) * t529 + t367;
t285 = -pkin(5) * t341 - qJ(6) * t534 + t331;
t456 = t285 * t534 + t471;
t455 = t444 * t288 + t446 * t291 + t297 * t485 - t307 * t486;
t454 = t444 * t294 + t446 * t300 + t315 * t485 - t321 * t486;
t453 = t345 * t419 + t367 * t411 - t396 * t527;
t450 = -pkin(3) * t396 - qJ(4) * t397 + (-qJD(4) + t367) * t406;
t274 = pkin(5) * t310 + qJ(6) * t309 - qJD(6) * t534 + t324;
t371 = pkin(5) * t416 - qJ(6) * t418 + t436;
t364 = t416 * t419;
t363 = t418 * t419;
t323 = t411 * t418 + t419 * t528;
t322 = t410 * t419 + t411 * t416;
t305 = pkin(5) * t534 - qJ(6) * t341;
t299 = t363 * pkin(5) + t364 * qJ(6) + t357;
t284 = -t309 - t545;
t283 = pkin(5) * t457 - t465;
t282 = -qJ(6) * t457 + t532;
t277 = qJ(6) * t401 + t279;
t276 = -pkin(5) * t401 + t484;
t275 = pkin(5) * t323 + qJ(6) * t322 + qJD(6) * t364 + t330;
t273 = -pkin(5) * t412 - t523;
t272 = qJ(6) * t412 - qJD(6) * t457 + t454;
t271 = t471 - t520;
t270 = qJD(6) * t401 + t455 + t515;
t1 = [(t396 * t419 + t411 * t526) * MDP(8) + (t396 * t457 - t397 * t419 - t406 * t411 - t412 * t526) * MDP(9) + t411 * qJD(3) * MDP(10) - t412 * qJD(3) * MDP(11) + (-qJD(3) * t356 + t397 * t437 + t412 * t422) * MDP(13) + (-qJD(3) * t355 + t396 * t437 + t411 * t422) * MDP(14) + (-t301 * t457 + t313 * t406 + t319 * t412 + t328 * t397 - t356 * t529 + t440 * t453) * MDP(15) + (t302 * t457 - t314 * t406 - t320 * t412 - t329 * t397 + t356 * t387 + t442 * t453) * MDP(16) + (-t313 * t387 - t314 * t395 + (qJD(3) * t314 - t301 * t419 - t319 * t411 - t328 * t396) * t442 + (-t302 * t419 - t320 * t411 - t329 * t396) * t440) * MDP(17) + (t301 * t328 + t302 * t329 + t313 * t319 + t314 * t320 - t345 * t527 + t356 * t367) * MDP(18) + (t309 * t364 - t322 * t534) * MDP(19) + (t309 * t363 + t310 * t364 - t322 * t341 - t323 * t534) * MDP(20) + (t309 * t457 - t322 * t401 - t364 * t397 + t412 * t534) * MDP(21) + (t310 * t457 - t323 * t401 + t341 * t412 - t363 * t397) * MDP(22) + (-t397 * t457 + t401 * t412) * MDP(23) + (t278 * t412 + t357 * t310 + t331 * t323 + t324 * t363 - t330 * t341 + t465 * t397 + t401 * t523 + t457 * t471) * MDP(24) + (-t279 * t412 - t357 * t309 - t331 * t322 - t324 * t364 + t330 * t534 - t397 * t532 - t401 * t454 + t455 * t457) * MDP(25) + (t271 * t457 - t273 * t401 + t274 * t363 - t275 * t341 - t276 * t412 - t283 * t397 + t285 * t323 + t299 * t310) * MDP(26) + (-t270 * t363 - t271 * t364 + t272 * t341 + t273 * t534 - t276 * t322 - t277 * t323 - t282 * t310 - t283 * t309) * MDP(27) + (-t270 * t457 + t272 * t401 + t274 * t364 - t275 * t534 + t277 * t412 + t282 * t397 + t285 * t322 + t299 * t309) * MDP(28) + (t270 * t282 + t271 * t283 + t272 * t277 + t273 * t276 + t274 * t299 + t275 * t285) * MDP(29) + 0.2e1 * t482 * t525; 0.2e1 * MDP(13) * t535 + (t429 + (-t406 - t477) * qJD(3)) * MDP(14) + (t442 * t397 - t402 * t440 + t526 * t529) * MDP(15) + (-t387 * t526 - t397 * t440 - t402 * t442) * MDP(16) + ((t440 * t387 + t442 * t529) * t406 + (-t440 ^ 2 - t442 ^ 2) * t396) * MDP(17) + (t301 * t442 + t302 * t440 - t367 * t526 + t406 * t463) * MDP(18) + (t470 - t524) * MDP(27) + (t270 * t418 + t271 * t416 + t276 * t490 + t277 * t491 - t285 * t526) * MDP(29) - qJD(1) ^ 2 * t525 + t543 * (t468 + t510) + t544 * (t469 + t512); -t402 * MDP(9) + (t429 + (t406 - t477) * qJD(3)) * MDP(10) + (qJD(3) * t376 - t345) * MDP(13) + (t422 * t406 - t457 * t482) * MDP(14) + (-t325 * t406 - t345 * t442 + t376 * t529 + t440 * t450) * MDP(15) + (t326 * t406 + t345 * t440 - t376 * t387 + t442 * t450) * MDP(16) + (t325 * t387 + t326 * t395 + (-qJD(4) * t395 - t319 * t406 + t302 + (qJD(4) * t442 - t326) * qJD(3)) * t442 + (qJD(4) * t387 - t320 * t406 - t301) * t440) * MDP(17) + (-pkin(3) * t345 - t319 * t325 - t320 * t326 - t367 * t376 + t463 * qJD(4) + (-t301 * t440 + t302 * t442) * qJ(4)) * MDP(18) + (-t309 * t418 + t491 * t534) * MDP(19) + (t470 + t524) * MDP(20) + (t468 - t510) * MDP(21) + (t469 - t512) * MDP(22) + (t436 * t310 + t324 * t416 + t490 * t331 + t341 * t346 + t401 * t530 + t508) * MDP(24) + (-t436 * t309 + t324 * t418 + t491 * t331 - t346 * t534 + t401 * t531 - t507) * MDP(25) + (t274 * t416 + t285 * t490 + t310 * t371 + t341 * t495 + t401 * t496 + t508) * MDP(26) + (-t270 * t416 + t271 * t418 + t276 * t491 - t277 * t490 + t309 * t462 - t310 * t384 - t341 * t497 - t496 * t534) * MDP(27) + (-t274 * t418 - t285 * t491 + t309 * t371 - t401 * t497 + t495 * t534 + t507) * MDP(28) + (t270 * t384 - t271 * t462 + t274 * t371 - t276 * t496 - t277 * t497 - t285 * t495) * MDP(29) + (-t422 * MDP(13) - t319 * MDP(15) + t320 * MDP(16) - t401 * MDP(23) - t278 * MDP(24) + t279 * MDP(25) + t276 * MDP(26) - t277 * MDP(28) + MDP(8) * t406 + MDP(9) * t526) * t526; (t387 * t406 + t506) * MDP(15) + (t406 * t529 + t501) * MDP(16) + (-t387 ^ 2 - t529 ^ 2) * MDP(17) + (t319 * t387 - t320 * t529 + t345) * MDP(18) + (-t522 - t546) * MDP(27) + (-t276 * t534 - t277 * t341 + t274) * MDP(29) + t544 * (t401 * t534 + t310) + t543 * (t309 - t545); -MDP(19) * t511 + (t522 - t546) * MDP(20) + t284 * MDP(21) + (-t533 + (-qJD(5) + t401) * t534) * MDP(22) + t397 * MDP(23) + (-t331 * t534 - t471 + t514) * MDP(24) + (t278 * t401 - t331 * t341 - t455) * MDP(25) + (t305 * t341 - t456 + t514 + 0.2e1 * t520) * MDP(26) + (pkin(5) * t309 - qJ(6) * t310 + (t277 - t279) * t534 - (t276 - t484) * t341) * MDP(27) + (0.2e1 * t515 + t285 * t341 + t305 * t534 + (0.2e1 * qJD(6) - t278) * t401 + t455) * MDP(28) + (-pkin(5) * t271 + qJ(6) * t270 - t276 * t279 + t277 * t484 - t285 * t305) * MDP(29); (-t511 - t535) * MDP(26) + t284 * MDP(27) + (-t401 ^ 2 - t522) * MDP(28) + (-t277 * t401 + t456 - t520) * MDP(29);];
tauc  = t1;
