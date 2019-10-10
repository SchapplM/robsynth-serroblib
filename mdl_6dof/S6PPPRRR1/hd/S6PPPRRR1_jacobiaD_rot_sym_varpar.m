% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPPRRR1
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
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6PPPRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPPRRR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPPRRR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobiaD_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:22
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (431->14), mult. (1286->36), div. (18->4), fcn. (1680->14), ass. (0->28)
	t139 = sin(pkin(13));
	t145 = cos(pkin(13));
	t146 = cos(pkin(12));
	t140 = sin(pkin(12));
	t156 = t140 * cos(pkin(6));
	t137 = -t139 * t156 + t146 * t145;
	t138 = sin(pkin(14));
	t144 = cos(pkin(14));
	t136 = -t146 * t139 - t145 * t156;
	t142 = sin(pkin(7));
	t148 = cos(pkin(7));
	t157 = t140 * sin(pkin(6));
	t154 = t136 * t148 + t142 * t157;
	t133 = t137 * t144 + t154 * t138;
	t150 = sin(qJ(4));
	t151 = cos(qJ(4));
	t155 = (-t136 * t142 + t148 * t157) * sin(pkin(8)) + (-t137 * t138 + t154 * t144) * cos(pkin(8));
	t129 = t133 * t150 - t155 * t151;
	t130 = t133 * t151 + t155 * t150;
	t127 = 0.1e1 / t130 ^ 2;
	t164 = qJD(4) * t127;
	t126 = t129 ^ 2;
	t123 = t126 * t127 + 0.1e1;
	t161 = t130 * t164;
	t162 = t129 / t130 * t164;
	t163 = (t126 * t162 + t129 * t161) / t123 ^ 2;
	t121 = 0.1e1 / t123;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, -0.2e1 * t163 + 0.2e1 * (t121 * t161 + (t121 * t162 - t127 * t163) * t129) * t129, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:22
	% EndTime: 2019-10-10 08:49:24
	% DurationCPUTime: 1.69s
	% Computational Cost: add. (7555->76), mult. (22578->165), div. (281->12), fcn. (29955->19), ass. (0->93)
	t289 = sin(pkin(14));
	t339 = sin(pkin(13));
	t340 = sin(pkin(12));
	t318 = t340 * t339;
	t345 = cos(pkin(13));
	t346 = cos(pkin(12));
	t324 = t346 * t345;
	t349 = cos(pkin(6));
	t308 = t324 * t349 - t318;
	t320 = t340 * t345;
	t323 = t346 * t339;
	t309 = t323 * t349 + t320;
	t343 = sin(pkin(6));
	t348 = cos(pkin(7));
	t325 = t348 * t343;
	t341 = sin(pkin(8));
	t342 = sin(pkin(7));
	t344 = cos(pkin(14));
	t347 = cos(pkin(8));
	t322 = t343 * t342;
	t354 = t308 * t348 - t346 * t322;
	t356 = (-t309 * t289 + t354 * t344) * t347 + (-t308 * t342 - t325 * t346) * t341;
	t311 = -t318 * t349 + t324;
	t310 = -t320 * t349 - t323;
	t319 = t340 * t343;
	t353 = t310 * t348 + t342 * t319;
	t301 = t311 * t289 - t353 * t344;
	t304 = -t310 * t342 + t319 * t348;
	t355 = t301 * t347 - t304 * t341;
	t352 = t345 * t325 + t349 * t342;
	t321 = t343 * t339;
	t351 = t347 * (-t289 * t321 + t352 * t344) + t341 * (-t345 * t322 + t349 * t348);
	t282 = t354 * t289 + t309 * t344;
	t291 = sin(qJ(4));
	t350 = cos(qJ(4));
	t266 = t282 * t291 - t356 * t350;
	t287 = t352 * t289 + t344 * t321;
	t307 = -t287 * t291 + t351 * t350;
	t259 = atan2(-t266, -t307);
	t254 = sin(t259);
	t255 = cos(t259);
	t242 = -t254 * t266 - t255 * t307;
	t239 = 0.1e1 / t242;
	t283 = t353 * t289 + t311 * t344;
	t270 = t283 * t350 - t355 * t291;
	t279 = t301 * t341 + t304 * t347;
	t290 = sin(qJ(5));
	t292 = cos(qJ(5));
	t253 = t270 * t292 + t279 * t290;
	t249 = 0.1e1 / t253;
	t274 = 0.1e1 / t307;
	t240 = 0.1e1 / t242 ^ 2;
	t250 = 0.1e1 / t253 ^ 2;
	t275 = 0.1e1 / t307 ^ 2;
	t264 = t266 ^ 2;
	t258 = t264 * t275 + 0.1e1;
	t256 = 0.1e1 / t258;
	t268 = t282 * t350 + t356 * t291;
	t261 = t268 * qJD(4);
	t278 = t287 * t350 + t351 * t291;
	t272 = t278 * qJD(4);
	t333 = t266 * t275;
	t233 = (t261 * t274 + t272 * t333) * t256;
	t317 = t254 * t307 - t255 * t266;
	t230 = t233 * t317 - t254 * t261 + t255 * t272;
	t338 = t230 * t239 * t240;
	t252 = t270 * t290 - t279 * t292;
	t248 = t252 ^ 2;
	t245 = t248 * t250 + 0.1e1;
	t269 = t283 * t291 + t355 * t350;
	t262 = t269 * qJD(4);
	t246 = qJD(5) * t253 - t262 * t290;
	t334 = t250 * t252;
	t330 = qJD(5) * t252;
	t247 = -t262 * t292 - t330;
	t335 = t247 * t249 * t250;
	t337 = (t246 * t334 - t248 * t335) / t245 ^ 2;
	t336 = t240 * t269;
	t332 = t266 * t278;
	t331 = t272 * t274 * t275;
	t329 = -0.2e1 * t337;
	t316 = -t249 * t290 + t292 * t334;
	t315 = t268 * t274 + t275 * t332;
	t271 = t307 * qJD(4);
	t265 = t269 ^ 2;
	t263 = t270 * qJD(4);
	t260 = t266 * qJD(4);
	t243 = 0.1e1 / t245;
	t237 = t265 * t240 + 0.1e1;
	t234 = t315 * t256;
	t231 = t234 * t317 - t254 * t268 + t255 * t278;
	t229 = -0.2e1 * t315 / t258 ^ 2 * (t261 * t333 + t264 * t331) + (0.2e1 * t331 * t332 - t260 * t274 + (t261 * t278 + t266 * t271 + t268 * t272) * t275) * t256;
	t1 = [0, 0, 0, t229, 0, 0; 0, 0, 0, 0.2e1 * (t231 * t336 - t239 * t270) / t237 ^ 2 * (t263 * t336 - t265 * t338) + (-t262 * t239 + (-t270 * t230 - t231 * t263) * t240 + (0.2e1 * t231 * t338 + (-(-t229 * t266 - t234 * t261 + t271 + (t234 * t307 - t268) * t233) * t255 - (t229 * t307 - t234 * t272 + t260 + (t234 * t266 - t278) * t233) * t254) * t240) * t269) / t237, 0, 0; 0, 0, 0, t316 * t269 * t329 + (t316 * t263 + ((-qJD(5) * t249 - 0.2e1 * t252 * t335) * t292 + (t246 * t292 + (t247 - t330) * t290) * t250) * t269) * t243, t329 + 0.2e1 * (t243 * t246 * t250 + (-t243 * t335 - t250 * t337) * t252) * t252, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:23
	% EndTime: 2019-10-10 08:49:28
	% DurationCPUTime: 4.98s
	% Computational Cost: add. (25805->134), mult. (75366->274), div. (559->12), fcn. (100504->21), ass. (0->138)
	t476 = sin(pkin(12));
	t481 = cos(pkin(13));
	t442 = t476 * t481;
	t475 = sin(pkin(13));
	t482 = cos(pkin(12));
	t447 = t482 * t475;
	t485 = cos(pkin(6));
	t430 = -t485 * t442 - t447;
	t484 = cos(pkin(7));
	t427 = t430 * t484;
	t441 = t476 * t475;
	t448 = t482 * t481;
	t431 = -t485 * t441 + t448;
	t474 = sin(pkin(14));
	t478 = sin(pkin(7));
	t443 = t478 * t474;
	t479 = sin(pkin(6));
	t434 = t479 * t443;
	t480 = cos(pkin(14));
	t398 = t474 * t427 + t431 * t480 + t476 * t434;
	t406 = sin(qJ(4));
	t409 = cos(qJ(4));
	t444 = t478 * t480;
	t435 = t479 * t444;
	t422 = -t480 * t427 + t431 * t474 - t476 * t435;
	t449 = t484 * t479;
	t425 = -t430 * t478 + t476 * t449;
	t477 = sin(pkin(8));
	t483 = cos(pkin(8));
	t487 = t422 * t483 - t425 * t477;
	t388 = t398 * t406 + t487 * t409;
	t428 = t485 * t448 - t441;
	t426 = t428 * t484;
	t429 = t485 * t447 + t442;
	t397 = t474 * t426 + t429 * t480 - t482 * t434;
	t421 = -t480 * t426 + t429 * t474 + t482 * t435;
	t424 = -t428 * t478 - t482 * t449;
	t417 = -t421 * t483 + t424 * t477;
	t387 = t397 * t409 + t417 * t406;
	t405 = sin(qJ(5));
	t408 = cos(qJ(5));
	t416 = t421 * t477 + t424 * t483;
	t373 = t387 * t408 + t416 * t405;
	t386 = -t397 * t406 + t417 * t409;
	t376 = t386 * qJD(4);
	t349 = t373 * qJD(5) + t376 * t405;
	t371 = t387 * t405 - t416 * t408;
	t369 = t371 ^ 2;
	t446 = t479 * t481;
	t436 = t484 * t446;
	t445 = t479 * t475;
	t402 = t474 * t436 + t485 * t443 + t480 * t445;
	t401 = t480 * t436 + t485 * t444 - t474 * t445;
	t403 = -t478 * t446 + t485 * t484;
	t433 = t483 * t401 + t477 * t403;
	t395 = t402 * t409 + t433 * t406;
	t399 = -t401 * t477 + t403 * t483;
	t384 = t395 * t405 - t399 * t408;
	t382 = 0.1e1 / t384 ^ 2;
	t363 = t369 * t382 + 0.1e1;
	t361 = 0.1e1 / t363;
	t385 = t395 * t408 + t399 * t405;
	t394 = -t402 * t406 + t433 * t409;
	t390 = t394 * qJD(4);
	t367 = t385 * qJD(5) + t390 * t405;
	t381 = 0.1e1 / t384;
	t463 = t371 * t382;
	t333 = (-t349 * t381 + t367 * t463) * t361;
	t364 = atan2(-t371, t384);
	t359 = sin(t364);
	t360 = cos(t364);
	t440 = -t359 * t384 - t360 * t371;
	t329 = t440 * t333 - t359 * t349 + t360 * t367;
	t343 = -t359 * t371 + t360 * t384;
	t340 = 0.1e1 / t343;
	t341 = 0.1e1 / t343 ^ 2;
	t490 = t329 * t340 * t341;
	t389 = t398 * t409 - t406 * t487;
	t418 = t422 * t477 + t425 * t483;
	t374 = t389 * t405 - t418 * t408;
	t489 = 0.2e1 * t374 * t490;
	t437 = -t381 * t386 + t394 * t463;
	t488 = t405 * t437;
	t464 = t367 * t381 * t382;
	t486 = -0.2e1 * (t349 * t463 - t369 * t464) / t363 ^ 2;
	t375 = t389 * t408 + t418 * t405;
	t404 = sin(qJ(6));
	t407 = cos(qJ(6));
	t358 = t375 * t407 + t388 * t404;
	t354 = 0.1e1 / t358;
	t355 = 0.1e1 / t358 ^ 2;
	t378 = t388 * qJD(4);
	t352 = -t374 * qJD(5) - t378 * t408;
	t379 = t389 * qJD(4);
	t344 = t358 * qJD(6) + t352 * t404 - t379 * t407;
	t357 = t375 * t404 - t388 * t407;
	t353 = t357 ^ 2;
	t348 = t353 * t355 + 0.1e1;
	t468 = t355 * t357;
	t456 = qJD(6) * t357;
	t345 = t352 * t407 + t379 * t404 - t456;
	t471 = t345 * t354 * t355;
	t473 = (t344 * t468 - t353 * t471) / t348 ^ 2;
	t472 = t341 * t374;
	t351 = t375 * qJD(5) - t378 * t405;
	t470 = t351 * t341;
	t469 = t354 * t404;
	t467 = t357 * t407;
	t466 = t359 * t374;
	t465 = t360 * t374;
	t462 = t388 * t405;
	t461 = t388 * t408;
	t457 = qJD(5) * t408;
	t370 = t374 ^ 2;
	t339 = t370 * t341 + 0.1e1;
	t455 = 0.2e1 * (-t370 * t490 + t374 * t470) / t339 ^ 2;
	t454 = -0.2e1 * t473;
	t452 = t357 * t471;
	t451 = -0.2e1 * t371 * t464;
	t450 = qJD(6) * t461 - t378;
	t439 = t355 * t467 - t469;
	t438 = -t373 * t381 + t385 * t463;
	t432 = qJD(5) * t462 + qJD(6) * t389 - t379 * t408;
	t391 = t395 * qJD(4);
	t377 = t387 * qJD(4);
	t368 = -t384 * qJD(5) + t390 * t408;
	t366 = t389 * t404 - t407 * t461;
	t365 = -t389 * t407 - t404 * t461;
	t350 = -t371 * qJD(5) + t376 * t408;
	t346 = 0.1e1 / t348;
	t337 = 0.1e1 / t339;
	t335 = t361 * t488;
	t334 = t438 * t361;
	t331 = (-t359 * t386 + t360 * t394) * t405 + t440 * t335;
	t330 = t440 * t334 - t359 * t373 + t360 * t385;
	t327 = t438 * t486 + (t385 * t451 - t350 * t381 + (t349 * t385 + t367 * t373 + t368 * t371) * t382) * t361;
	t326 = t486 * t488 + (t437 * t457 + (t394 * t451 + t377 * t381 + (t349 * t394 + t367 * t386 - t371 * t391) * t382) * t405) * t361;
	t1 = [0, 0, 0, t326, t327, 0; 0, 0, 0, (t331 * t472 + t340 * t462) * t455 + ((-t379 * t405 - t388 * t457) * t340 + (-t470 + t489) * t331 + (t462 * t329 - (t394 * t457 - t326 * t371 - t335 * t349 - t391 * t405 + (-t335 * t384 - t386 * t405) * t333) * t465 - (-t386 * t457 - t326 * t384 - t335 * t367 + t377 * t405 + (t335 * t371 - t394 * t405) * t333) * t466) * t341) * t337, (t330 * t472 - t340 * t375) * t455 + (t330 * t489 + t352 * t340 + (-t375 * t329 - t330 * t351 - (-t327 * t371 - t334 * t349 + t368 + (-t334 * t384 - t373) * t333) * t465 - (-t327 * t384 - t334 * t367 - t350 + (t334 * t371 - t385) * t333) * t466) * t341) * t337, 0; 0, 0, 0, 0.2e1 * (-t354 * t365 + t366 * t468) * t473 + (0.2e1 * t366 * t452 - t450 * t354 * t407 + t432 * t469 + (-t450 * t357 * t404 - t366 * t344 - t365 * t345 - t432 * t467) * t355) * t346, t439 * t374 * t454 + (t439 * t351 + ((-qJD(6) * t354 - 0.2e1 * t452) * t407 + (t344 * t407 + (t345 - t456) * t404) * t355) * t374) * t346, t454 + 0.2e1 * (t344 * t355 * t346 + (-t346 * t471 - t355 * t473) * t357) * t357;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end