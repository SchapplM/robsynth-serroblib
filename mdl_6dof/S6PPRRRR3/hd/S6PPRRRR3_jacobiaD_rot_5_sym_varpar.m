% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRRR3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobiaD_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:56
% EndTime: 2019-02-26 19:43:59
% DurationCPUTime: 2.50s
% Computational Cost: add. (15480->143), mult. (47073->274), div. (538->12), fcn. (61235->19), ass. (0->130)
t446 = sin(pkin(14));
t447 = sin(pkin(13));
t411 = t447 * t446;
t450 = cos(pkin(14));
t451 = cos(pkin(13));
t417 = t451 * t450;
t453 = cos(pkin(6));
t391 = -t453 * t417 + t411;
t448 = sin(pkin(7));
t449 = sin(pkin(6));
t415 = t449 * t448;
t452 = cos(pkin(7));
t460 = -t391 * t452 - t451 * t415;
t413 = t447 * t450;
t416 = t451 * t446;
t360 = -t453 * t413 - t416;
t412 = t447 * t449;
t459 = t452 * t360 + t448 * t412;
t418 = t452 * t449;
t458 = t450 * t418 + t453 * t448;
t368 = sin(qJ(4));
t369 = sin(qJ(3));
t392 = t453 * t416 + t413;
t455 = cos(qJ(3));
t381 = -t392 * t369 + t455 * t460;
t454 = cos(qJ(4));
t378 = t381 * t454;
t348 = t369 * t460 + t392 * t455;
t380 = qJD(3) * t348;
t457 = qJD(4) * t378 - t368 * t380;
t366 = cos(pkin(8));
t365 = sin(pkin(8));
t383 = (t391 * t448 - t451 * t418) * t365;
t325 = t348 * t454 + (t381 * t366 + t383) * t368;
t345 = t381 * qJD(3);
t425 = t366 * t454;
t299 = t325 * qJD(4) + t345 * t368 + t380 * t425;
t382 = t454 * t383;
t436 = t348 * t368;
t323 = -t366 * t378 - t382 + t436;
t321 = t323 ^ 2;
t414 = t449 * t446;
t356 = t369 * t458 + t455 * t414;
t355 = -t369 * t414 + t455 * t458;
t424 = t454 * t355;
t435 = (-t450 * t415 + t453 * t452) * t365;
t393 = -t356 * t368 + t366 * t424 + t454 * t435;
t334 = 0.1e1 / t393 ^ 2;
t313 = t321 * t334 + 0.1e1;
t437 = t323 * t334;
t337 = t356 * t454 + (t355 * t366 + t435) * t368;
t353 = t355 * qJD(3);
t354 = t356 * qJD(3);
t319 = t337 * qJD(4) + t353 * t368 + t354 * t425;
t333 = 0.1e1 / t393;
t438 = t319 * t333 * t334;
t456 = -0.2e1 * (t299 * t437 + t321 * t438) / t313 ^ 2;
t361 = -t453 * t411 + t417;
t349 = -t361 * t369 + t455 * t459;
t350 = t361 * t455 + t369 * t459;
t314 = atan2(-t323, -t393);
t309 = sin(t314);
t310 = cos(t314);
t293 = -t309 * t323 - t310 * t393;
t290 = 0.1e1 / t293;
t395 = -t360 * t448 + t452 * t412;
t389 = t395 * t365;
t327 = t350 * t454 + (t349 * t366 + t389) * t368;
t338 = -t349 * t365 + t395 * t366;
t367 = sin(qJ(5));
t370 = cos(qJ(5));
t308 = t327 * t370 + t338 * t367;
t304 = 0.1e1 / t308;
t291 = 0.1e1 / t293 ^ 2;
t305 = 0.1e1 / t308 ^ 2;
t326 = -t349 * t425 + t350 * t368 - t454 * t389;
t322 = t326 ^ 2;
t289 = t291 * t322 + 0.1e1;
t346 = t349 * qJD(3);
t347 = t350 * qJD(3);
t301 = t327 * qJD(4) + t346 * t368 + t347 * t425;
t442 = t291 * t326;
t311 = 0.1e1 / t313;
t283 = (t299 * t333 + t319 * t437) * t311;
t410 = t309 * t393 - t310 * t323;
t279 = t410 * t283 - t309 * t299 + t310 * t319;
t444 = t279 * t290 * t291;
t445 = (t301 * t442 - t322 * t444) / t289 ^ 2;
t287 = 0.1e1 / t289;
t443 = t287 * t291;
t431 = t366 * t368;
t302 = -t326 * qJD(4) + t346 * t454 - t347 * t431;
t307 = t327 * t367 - t338 * t370;
t430 = qJD(5) * t307;
t433 = t365 * t367;
t295 = t302 * t370 + t347 * t433 - t430;
t441 = t295 * t304 * t305;
t432 = t365 * t370;
t294 = t308 * qJD(5) + t302 * t367 - t347 * t432;
t303 = t307 ^ 2;
t298 = t303 * t305 + 0.1e1;
t439 = t305 * t307;
t440 = 0.1e1 / t298 ^ 2 * (t294 * t439 - t303 * t441);
t429 = 0.2e1 * t445;
t428 = -0.2e1 * t440;
t427 = t307 * t441;
t423 = qJD(4) * t436;
t421 = 0.2e1 * t326 * t444;
t420 = 0.2e1 * t323 * t438;
t407 = -t304 * t367 + t370 * t439;
t406 = t325 * t333 + t337 * t437;
t329 = t348 * t425 + t381 * t368;
t339 = t355 * t368 + t356 * t425;
t405 = t329 * t333 + t339 * t437;
t331 = t349 * t454 - t350 * t431;
t404 = -t331 * t367 + t350 * t432;
t318 = t331 * t370 + t350 * t433;
t398 = -t349 * t368 - t350 * t425;
t328 = t353 * t425 - t354 * t368 + (-t356 * t431 + t424) * qJD(4);
t320 = t393 * qJD(4) + t353 * t454 - t354 * t431;
t316 = t398 * qJD(4) - t346 * t431 - t347 * t454;
t315 = t345 * t425 - t366 * t423 + t457;
t300 = qJD(4) * t382 + t345 * t454 + t457 * t366 - t423;
t296 = 0.1e1 / t298;
t285 = t405 * t311;
t284 = t406 * t311;
t280 = t410 * t284 - t309 * t325 + t310 * t337;
t278 = t405 * t456 + (t339 * t420 + t315 * t333 + (t299 * t339 + t319 * t329 + t323 * t328) * t334) * t311;
t276 = t406 * t456 + (t337 * t420 + t300 * t333 + (t299 * t337 + t319 * t325 + t320 * t323) * t334) * t311;
t1 = [0, 0, t278, t276, 0, 0; 0, 0 -(-t279 * t443 - 0.2e1 * t290 * t445) * t398 + ((t331 * qJD(4) + t346 * t425 - t347 * t368) * t290 - ((-t278 * t323 - t285 * t299 + t328 + (t285 * t393 - t329) * t283) * t310 + (t278 * t393 - t285 * t319 - t315 + (t285 * t323 - t339) * t283) * t309) * t442) * t287 + (t287 * t421 - t301 * t443 + t442 * t429) * (t410 * t285 - t309 * t329 + t310 * t339) (t280 * t442 - t290 * t327) * t429 + (t280 * t421 + t302 * t290 + (-t327 * t279 - t280 * t301 + (-(-t276 * t323 - t284 * t299 + t320 + (t284 * t393 - t325) * t283) * t310 - (t276 * t393 - t284 * t319 - t300 + (t284 * t323 - t337) * t283) * t309) * t326) * t291) * t287, 0, 0; 0, 0, 0.2e1 * (t304 * t404 + t318 * t439) * t440 + ((t318 * qJD(5) + t316 * t367 - t346 * t432) * t304 + 0.2e1 * t318 * t427 + (t404 * t295 - (t404 * qJD(5) + t316 * t370 + t346 * t433) * t307 - t318 * t294) * t305) * t296, t407 * t326 * t428 + (t407 * t301 + ((-qJD(5) * t304 - 0.2e1 * t427) * t370 + (t294 * t370 + (t295 - t430) * t367) * t305) * t326) * t296, t428 + 0.2e1 * (t294 * t305 * t296 + (-t296 * t441 - t305 * t440) * t307) * t307, 0;];
JaD_rot  = t1;
