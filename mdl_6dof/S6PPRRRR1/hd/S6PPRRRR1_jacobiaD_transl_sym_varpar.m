% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PPRRRR1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:18
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:21
	% EndTime: 2019-10-09 21:18:21
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (23->19), mult. (86->45), div. (0->0), fcn. (88->10), ass. (0->22)
	t111 = sin(pkin(12));
	t117 = cos(pkin(6));
	t126 = t111 * t117;
	t112 = sin(pkin(7));
	t113 = sin(pkin(6));
	t125 = t112 * t113;
	t115 = cos(pkin(12));
	t124 = t115 * t117;
	t116 = cos(pkin(7));
	t118 = sin(qJ(3));
	t123 = t116 * t118;
	t119 = cos(qJ(3));
	t122 = t116 * t119;
	t121 = t118 * t125;
	t120 = t119 * t125;
	t114 = cos(pkin(13));
	t110 = sin(pkin(13));
	t109 = -t110 * t126 + t115 * t114;
	t108 = -t115 * t110 - t114 * t126;
	t107 = t110 * t124 + t111 * t114;
	t106 = -t111 * t110 + t114 * t124;
	t1 = [0, 0, ((-t108 * t123 - t109 * t119 - t111 * t121) * r_i_i_C(1) + (-t108 * t122 + t109 * t118 - t111 * t120) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, 0, ((-t106 * t123 - t107 * t119 + t115 * t121) * r_i_i_C(1) + (-t106 * t122 + t107 * t118 + t115 * t120) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, 0, ((-r_i_i_C(1) * t118 - r_i_i_C(2) * t119) * t117 * t112 + ((-t110 * t119 - t114 * t123) * r_i_i_C(1) + (t110 * t118 - t114 * t122) * r_i_i_C(2)) * t113) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:22
	% EndTime: 2019-10-09 21:18:22
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (153->45), mult. (509->95), div. (0->0), fcn. (560->12), ass. (0->45)
	t281 = sin(pkin(13));
	t284 = sin(pkin(6));
	t290 = sin(qJ(3));
	t292 = cos(qJ(3));
	t285 = cos(pkin(13));
	t287 = cos(pkin(7));
	t306 = t285 * t287;
	t283 = sin(pkin(7));
	t288 = cos(pkin(6));
	t309 = t283 * t288;
	t268 = (t281 * t292 + t290 * t306) * t284 + t290 * t309;
	t286 = cos(pkin(12));
	t282 = sin(pkin(12));
	t311 = t282 * t288;
	t277 = -t281 * t311 + t286 * t285;
	t276 = -t286 * t281 - t285 * t311;
	t310 = t283 * t284;
	t296 = t276 * t287 + t282 * t310;
	t264 = t277 * t292 + t296 * t290;
	t305 = t286 * t288;
	t275 = t281 * t305 + t282 * t285;
	t274 = -t282 * t281 + t285 * t305;
	t302 = t286 * t310;
	t297 = -t274 * t287 + t302;
	t316 = -t275 * t292 + t297 * t290;
	t315 = -pkin(9) - r_i_i_C(3);
	t314 = t275 * t290;
	t308 = t284 * t285;
	t307 = t284 * t287;
	t304 = qJD(3) * t290;
	t303 = qJD(3) * t292;
	t300 = t283 * t303;
	t299 = t287 * t303;
	t289 = sin(qJ(4));
	t291 = cos(qJ(4));
	t298 = t289 * r_i_i_C(1) + t291 * r_i_i_C(2);
	t294 = qJD(4) * t298;
	t293 = qJD(3) * (t291 * r_i_i_C(1) - t289 * r_i_i_C(2) + pkin(3));
	t273 = -t283 * t308 + t288 * t287;
	t270 = -t276 * t283 + t282 * t307;
	t269 = -t274 * t283 - t286 * t307;
	t265 = t284 * t281 * t304 - t288 * t300 - t299 * t308;
	t259 = -t282 * t284 * t300 - t276 * t299 + t277 * t304;
	t257 = -t274 * t299 + (t292 * t302 + t314) * qJD(3);
	t1 = [0, 0, t315 * t259 - (-t277 * t290 + t296 * t292) * t294 - t264 * t293, t298 * t259 + ((-t264 * t291 - t270 * t289) * r_i_i_C(1) + (t264 * t289 - t270 * t291) * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, t315 * t257 - (-t297 * t292 - t314) * t294 + t316 * t293, t298 * t257 + ((-t269 * t289 + t291 * t316) * r_i_i_C(1) + (-t269 * t291 - t289 * t316) * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, t315 * t265 - (t292 * t309 + (-t281 * t290 + t292 * t306) * t284) * t294 - t268 * t293, t298 * t265 + ((-t268 * t291 - t273 * t289) * r_i_i_C(1) + (t268 * t289 - t273 * t291) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:22
	% EndTime: 2019-10-09 21:18:22
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (349->58), mult. (870->112), div. (0->0), fcn. (970->14), ass. (0->56)
	t325 = sin(pkin(13));
	t328 = sin(pkin(6));
	t334 = sin(qJ(3));
	t336 = cos(qJ(3));
	t329 = cos(pkin(13));
	t331 = cos(pkin(7));
	t356 = t329 * t331;
	t327 = sin(pkin(7));
	t332 = cos(pkin(6));
	t359 = t327 * t332;
	t307 = (t325 * t336 + t334 * t356) * t328 + t334 * t359;
	t330 = cos(pkin(12));
	t326 = sin(pkin(12));
	t361 = t326 * t332;
	t316 = -t325 * t361 + t330 * t329;
	t315 = -t330 * t325 - t329 * t361;
	t360 = t327 * t328;
	t341 = t315 * t331 + t326 * t360;
	t303 = t316 * t336 + t341 * t334;
	t355 = t330 * t332;
	t314 = t325 * t355 + t326 * t329;
	t313 = -t326 * t325 + t329 * t355;
	t349 = t330 * t360;
	t342 = -t313 * t331 + t349;
	t368 = -t314 * t336 + t342 * t334;
	t367 = -r_i_i_C(3) - pkin(10) - pkin(9);
	t366 = t314 * t334;
	t324 = qJ(4) + qJ(5);
	t321 = sin(t324);
	t323 = qJD(4) + qJD(5);
	t363 = t321 * t323;
	t322 = cos(t324);
	t362 = t322 * t323;
	t358 = t328 * t329;
	t357 = t328 * t331;
	t350 = qJD(3) * t336;
	t346 = t331 * t350;
	t296 = -t313 * t346 + (t336 * t349 + t366) * qJD(3);
	t308 = -t313 * t327 - t330 * t357;
	t345 = -t308 * t323 + t296;
	t354 = (t345 * t321 + t362 * t368) * r_i_i_C(1) + (t345 * t322 - t363 * t368) * r_i_i_C(2);
	t347 = t327 * t350;
	t351 = qJD(3) * t334;
	t298 = -t326 * t328 * t347 - t315 * t346 + t316 * t351;
	t309 = -t315 * t327 + t326 * t357;
	t344 = -t309 * t323 + t298;
	t353 = (-t303 * t362 + t344 * t321) * r_i_i_C(1) + (t303 * t363 + t344 * t322) * r_i_i_C(2);
	t304 = t328 * t325 * t351 - t332 * t347 - t346 * t358;
	t312 = -t327 * t358 + t332 * t331;
	t343 = -t312 * t323 + t304;
	t352 = (-t307 * t362 + t343 * t321) * r_i_i_C(1) + (t307 * t363 + t343 * t322) * r_i_i_C(2);
	t335 = cos(qJ(4));
	t339 = qJD(3) * (t335 * pkin(4) + r_i_i_C(1) * t322 - r_i_i_C(2) * t321 + pkin(3));
	t333 = sin(qJ(4));
	t338 = -pkin(4) * qJD(4) * t333 + (-r_i_i_C(1) * t321 - r_i_i_C(2) * t322) * t323;
	t1 = [0, 0, t367 * t298 + t338 * (-t316 * t334 + t341 * t336) - t303 * t339, (t298 * t333 + (-t303 * t335 - t309 * t333) * qJD(4)) * pkin(4) + t353, t353, 0; 0, 0, t367 * t296 + t338 * (-t342 * t336 - t366) + t368 * t339, (t296 * t333 + (-t308 * t333 + t335 * t368) * qJD(4)) * pkin(4) + t354, t354, 0; 0, 0, t367 * t304 + t338 * (t336 * t359 + (-t325 * t334 + t336 * t356) * t328) - t307 * t339, (t304 * t333 + (-t307 * t335 - t312 * t333) * qJD(4)) * pkin(4) + t352, t352, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:25
	% EndTime: 2019-10-09 21:18:26
	% DurationCPUTime: 0.89s
	% Computational Cost: add. (1128->102), mult. (2701->183), div. (0->0), fcn. (3130->16), ass. (0->76)
	t529 = sin(pkin(13));
	t530 = sin(pkin(12));
	t510 = t530 * t529;
	t533 = cos(pkin(13));
	t534 = cos(pkin(12));
	t516 = t534 * t533;
	t536 = cos(pkin(6));
	t498 = -t536 * t516 + t510;
	t535 = cos(pkin(7));
	t496 = t498 * t535;
	t531 = sin(pkin(7));
	t532 = sin(pkin(6));
	t514 = t532 * t531;
	t506 = t534 * t514;
	t545 = t496 + t506;
	t512 = t530 * t533;
	t515 = t534 * t529;
	t499 = t536 * t512 + t515;
	t511 = t530 * t532;
	t544 = t499 * t535 - t531 * t511;
	t517 = t535 * t532;
	t543 = t533 * t517 + t536 * t531;
	t542 = pkin(11) + r_i_i_C(3);
	t487 = cos(qJ(6));
	t524 = qJD(6) * t487;
	t484 = sin(qJ(6));
	t525 = qJD(6) * t484;
	t541 = -r_i_i_C(1) * t525 - t524 * r_i_i_C(2);
	t509 = t487 * r_i_i_C(1) - t484 * r_i_i_C(2) + pkin(5);
	t470 = t536 * t515 + t512;
	t486 = sin(qJ(3));
	t538 = cos(qJ(3));
	t454 = t470 * t486 + t545 * t538;
	t540 = t544 * t538;
	t513 = t532 * t529;
	t461 = t486 * t513 - t543 * t538;
	t483 = qJ(4) + qJ(5);
	t480 = sin(t483);
	t482 = qJD(4) + qJD(5);
	t528 = t480 * t482;
	t481 = cos(t483);
	t527 = t481 * t482;
	t526 = qJD(3) * t486;
	t522 = t470 * t538;
	t448 = t454 * qJD(3);
	t463 = t498 * t531 - t534 * t517;
	t521 = t463 * t482 - t448;
	t471 = -t536 * t510 + t516;
	t450 = t540 * qJD(3) + t471 * t526;
	t464 = t499 * t531 + t535 * t511;
	t520 = t464 * t482 - t450;
	t459 = t461 * qJD(3);
	t469 = -t533 * t514 + t536 * t535;
	t519 = t469 * t482 - t459;
	t488 = cos(qJ(4));
	t501 = -t488 * pkin(4) - t480 * t542 - t509 * t481 - pkin(3);
	t455 = -t545 * t486 + t522;
	t434 = -t455 * t528 + t521 * t481;
	t495 = t541 * (-t455 * t480 + t463 * t481) + t542 * t434 + t509 * (-t455 * t527 - t521 * t480);
	t457 = t471 * t538 - t544 * t486;
	t436 = -t457 * t528 + t520 * t481;
	t494 = t541 * (-t457 * t480 + t464 * t481) + t542 * t436 + t509 * (-t457 * t527 - t520 * t480);
	t462 = t543 * t486 + t538 * t513;
	t441 = -t462 * t528 + t519 * t481;
	t493 = t541 * (-t462 * t480 + t469 * t481) + t542 * t441 + t509 * (-t462 * t527 - t519 * t480);
	t485 = sin(qJ(4));
	t490 = qJD(4) * t485 * pkin(4) + (t484 * r_i_i_C(1) + t487 * r_i_i_C(2)) * t481 * qJD(6) + (t509 * t480 - t481 * t542) * t482;
	t489 = -pkin(10) - pkin(9);
	t460 = t462 * qJD(3);
	t456 = t471 * t486 + t540;
	t453 = t462 * t481 + t469 * t480;
	t451 = t457 * qJD(3);
	t449 = -t506 * t526 + (-t486 * t496 + t522) * qJD(3);
	t445 = t457 * t481 + t464 * t480;
	t443 = t455 * t481 + t463 * t480;
	t1 = [0, 0, (-t450 * t484 + t457 * t524) * r_i_i_C(1) + (-t450 * t487 - t457 * t525) * r_i_i_C(2) + t450 * t489 + t501 * t451 + t490 * t456, (t450 * t485 + (-t457 * t488 - t464 * t485) * qJD(4)) * pkin(4) + t494, t494, (-t436 * t484 + t451 * t487) * r_i_i_C(1) + (-t436 * t487 - t451 * t484) * r_i_i_C(2) + ((-t445 * t487 - t456 * t484) * r_i_i_C(1) + (t445 * t484 - t456 * t487) * r_i_i_C(2)) * qJD(6); 0, 0, (-t448 * t484 + t455 * t524) * r_i_i_C(1) + (-t448 * t487 - t455 * t525) * r_i_i_C(2) + t448 * t489 + t501 * t449 + t490 * t454, (t448 * t485 + (-t455 * t488 - t463 * t485) * qJD(4)) * pkin(4) + t495, t495, (-t434 * t484 + t449 * t487) * r_i_i_C(1) + (-t434 * t487 - t449 * t484) * r_i_i_C(2) + ((-t443 * t487 - t454 * t484) * r_i_i_C(1) + (t443 * t484 - t454 * t487) * r_i_i_C(2)) * qJD(6); 0, 0, (-t459 * t484 + t462 * t524) * r_i_i_C(1) + (-t459 * t487 - t462 * t525) * r_i_i_C(2) + t459 * t489 + t501 * t460 + t490 * t461, (t459 * t485 + (-t462 * t488 - t469 * t485) * qJD(4)) * pkin(4) + t493, t493, (-t441 * t484 + t460 * t487) * r_i_i_C(1) + (-t441 * t487 - t460 * t484) * r_i_i_C(2) + ((-t453 * t487 - t461 * t484) * r_i_i_C(1) + (t453 * t484 - t461 * t487) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end