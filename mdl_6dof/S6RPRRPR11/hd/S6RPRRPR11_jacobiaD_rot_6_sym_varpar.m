% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR11_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:45
% EndTime: 2019-02-26 21:06:48
% DurationCPUTime: 3.55s
% Computational Cost: add. (14207->166), mult. (41237->320), div. (726->12), fcn. (53008->17), ass. (0->146)
t495 = sin(pkin(12));
t500 = cos(pkin(6));
t455 = t500 * t495;
t498 = cos(pkin(12));
t501 = sin(qJ(1));
t503 = cos(qJ(1));
t395 = -t501 * t455 + t503 * t498;
t392 = t395 * qJD(1);
t406 = sin(qJ(3));
t407 = cos(qJ(3));
t394 = t503 * t455 + t501 * t498;
t496 = sin(pkin(7));
t497 = sin(pkin(6));
t453 = t497 * t496;
t440 = t503 * t453;
t499 = cos(pkin(7));
t456 = t500 * t498;
t510 = -t503 * t456 + t501 * t495;
t511 = t510 * t499;
t417 = t511 + t440;
t506 = t394 * t406 + t417 * t407;
t426 = t501 * t456 + t503 * t495;
t391 = t426 * qJD(1);
t437 = t501 * t453;
t516 = qJD(1) * t437 - t391 * t499;
t354 = qJD(3) * t506 - t392 * t407 - t516 * t406;
t405 = sin(qJ(4));
t454 = t499 * t497;
t438 = t501 * t454;
t427 = qJD(1) * t438 + t391 * t496;
t502 = cos(qJ(4));
t476 = t394 * t407;
t375 = t417 * t406 - t476;
t418 = -t503 * t454 + t496 * t510;
t525 = -t375 * t405 - t418 * t502;
t332 = qJD(4) * t525 + t354 * t502 - t427 * t405;
t363 = t375 * t502 - t418 * t405;
t527 = t363 * qJD(4) + t354 * t405 + t427 * t502;
t357 = t525 ^ 2;
t423 = t498 * t454 + t500 * t496;
t452 = t497 * t495;
t386 = t423 * t406 + t407 * t452;
t393 = -t498 * t453 + t500 * t499;
t443 = -t386 * t405 + t393 * t502;
t368 = 0.1e1 / t443 ^ 2;
t345 = t357 * t368 + 0.1e1;
t343 = 0.1e1 / t345;
t371 = t386 * t502 + t393 * t405;
t385 = -t406 * t452 + t423 * t407;
t380 = t385 * qJD(3);
t355 = t371 * qJD(4) + t380 * t405;
t367 = 0.1e1 / t443;
t481 = t525 * t368;
t448 = t355 * t481 - t367 * t527;
t312 = t448 * t343;
t346 = atan2(-t525, -t443);
t341 = sin(t346);
t342 = cos(t346);
t451 = t341 * t443 - t342 * t525;
t307 = t451 * t312 + t341 * t527 + t342 * t355;
t324 = -t341 * t525 - t342 * t443;
t322 = 0.1e1 / t324 ^ 2;
t524 = t307 * t322;
t521 = -t426 * t499 + t437;
t377 = t395 * t407 + t521 * t406;
t416 = t426 * t496 + t438;
t365 = t377 * t502 + t416 * t405;
t508 = t521 * t407;
t376 = t395 * t406 - t508;
t404 = pkin(13) + qJ(6);
t402 = sin(t404);
t403 = cos(t404);
t339 = t365 * t402 - t376 * t403;
t520 = 0.2e1 * t339;
t321 = 0.1e1 / t324;
t519 = t321 * t524;
t390 = t394 * qJD(1);
t475 = qJD(3) * t406;
t509 = t417 * qJD(1);
t350 = qJD(3) * t508 - t390 * t407 - t395 * t475 + t406 * t509;
t415 = qJD(1) * t418;
t328 = t365 * qJD(4) + t350 * t405 + t502 * t415;
t414 = t416 * t502;
t364 = t377 * t405 - t414;
t504 = 0.2e1 * t364;
t460 = t504 * t519;
t515 = -t322 * t328 + t460;
t514 = (t406 * t511 - t476) * qJD(3) - t392 * t406 + t516 * t407 + t440 * t475;
t358 = t364 ^ 2;
t320 = t322 * t358 + 0.1e1;
t488 = t322 * t364;
t472 = 0.2e1 * (t328 * t488 - t358 * t519) / t320 ^ 2;
t513 = t355 * t368;
t445 = -t367 * t506 + t385 * t481;
t512 = t405 * t445;
t462 = qJD(4) * t502;
t340 = t365 * t403 + t376 * t402;
t334 = 0.1e1 / t340;
t335 = 0.1e1 / t340 ^ 2;
t505 = -0.2e1 * t525;
t474 = qJD(4) * t405;
t329 = qJD(4) * t414 + t350 * t502 - t377 * t474 - t405 * t415;
t349 = t377 * qJD(3) - t390 * t406 - t407 * t509;
t316 = t340 * qJD(6) + t329 * t402 - t349 * t403;
t333 = t339 ^ 2;
t327 = t333 * t335 + 0.1e1;
t486 = t335 * t339;
t473 = qJD(6) * t339;
t317 = t329 * t403 + t349 * t402 - t473;
t490 = t317 * t334 * t335;
t493 = (t316 * t486 - t333 * t490) / t327 ^ 2;
t483 = t367 * t513;
t491 = (t357 * t483 - t481 * t527) / t345 ^ 2;
t325 = 0.1e1 / t327;
t487 = t325 * t335;
t485 = t341 * t364;
t484 = t342 * t364;
t482 = t525 * t367;
t479 = t376 * t405;
t471 = -0.2e1 * t493;
t470 = -0.2e1 * t491;
t468 = t335 * t493;
t467 = t367 * t491;
t466 = t316 * t487;
t465 = t339 * t490;
t463 = t376 * t502;
t459 = 0.2e1 * t465;
t458 = t483 * t505;
t338 = t363 * t403 - t402 * t506;
t337 = t363 * t402 + t403 * t506;
t449 = qJD(6) * t463 + t350;
t447 = -t334 * t402 + t403 * t486;
t446 = -t363 * t367 + t371 * t481;
t436 = -t341 + (-t342 * t482 + t341) * t343;
t430 = qJD(6) * t377 - t502 * t349 + t376 * t474;
t381 = t386 * qJD(3);
t356 = t443 * qJD(4) + t380 * t502;
t348 = t377 * t402 - t403 * t463;
t318 = 0.1e1 / t320;
t315 = t343 * t512;
t313 = t446 * t343;
t309 = (t341 * t506 + t342 * t385) * t405 + t451 * t315;
t308 = t451 * t313 + t341 * t363 + t342 * t371;
t305 = t446 * t470 + (-t371 * t458 - t332 * t367 + (-t355 * t363 + t356 * t525 - t371 * t527) * t368) * t343;
t304 = t470 * t512 + ((-t385 * t458 + t514 * t367 + (-t355 * t506 - t381 * t525 - t385 * t527) * t368) * t405 + t445 * t462) * t343;
t1 = [-t467 * t504 + (t328 * t367 + t364 * t513) * t343, 0, t304, t305, 0, 0; t525 * t321 * t472 + (t527 * t321 + t525 * t524 - (t436 * t328 + ((t312 * t343 * t482 + t470) * t341 + (-t467 * t505 - t312 + (t312 - t448) * t343) * t342) * t364) * t488) * t318 + (t515 * t318 + t488 * t472) * t436 * t364, 0 (t309 * t488 + t321 * t479) * t472 + ((-t349 * t405 - t376 * t462) * t321 + t515 * t309 + (t479 * t307 - (t385 * t462 - t304 * t525 + t315 * t527 - t381 * t405 + (t315 * t443 + t405 * t506) * t312) * t484 - (t506 * t462 + t304 * t443 - t315 * t355 - t514 * t405 + (t315 * t525 - t385 * t405) * t312) * t485) * t322) * t318 (t308 * t488 - t321 * t365) * t472 + (t308 * t460 + t329 * t321 + (-t365 * t307 - t308 * t328 - (-t305 * t525 + t313 * t527 + t356 + (t313 * t443 + t363) * t312) * t484 - (t305 * t443 - t313 * t355 + t332 + (t313 * t525 - t371) * t312) * t485) * t322) * t318, 0, 0; 0.2e1 * (-t334 * t337 + t338 * t486) * t493 + ((t338 * qJD(6) + t332 * t402 - t403 * t514) * t334 + t338 * t459 + (-t337 * t317 - (-t337 * qJD(6) + t332 * t403 + t402 * t514) * t339 - t338 * t316) * t335) * t325, 0 (t468 * t520 - t466) * t348 + (-t317 * t487 + t334 * t471) * (-t377 * t403 - t402 * t463) + ((t430 * t402 - t449 * t403) * t334 - (t449 * t402 + t430 * t403) * t486 + t348 * t459) * t325, t447 * t364 * t471 + (t447 * t328 + ((-qJD(6) * t334 - 0.2e1 * t465) * t403 + (t316 * t403 + (t317 - t473) * t402) * t335) * t364) * t325, 0, t471 + (t466 + (-t325 * t490 - t468) * t339) * t520;];
JaD_rot  = t1;
