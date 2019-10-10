% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR9
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:31
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR9_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR9_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR9_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:25
	% EndTime: 2019-10-10 11:31:25
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (35->18), mult. (110->35), div. (0->0), fcn. (94->6), ass. (0->20)
	t136 = sin(pkin(6));
	t151 = t136 * (pkin(8) + r_i_i_C(3));
	t138 = sin(qJ(2));
	t139 = sin(qJ(1));
	t149 = t138 * t139;
	t141 = cos(qJ(1));
	t148 = t138 * t141;
	t140 = cos(qJ(2));
	t147 = t139 * t140;
	t146 = t140 * t141;
	t137 = cos(pkin(6));
	t145 = -t137 * t146 + t149;
	t144 = t137 * t147 + t148;
	t143 = t137 * t148 + t147;
	t142 = t137 * t149 - t146;
	t135 = t142 * qJD(1) + t145 * qJD(2);
	t134 = t144 * qJD(1) + t143 * qJD(2);
	t133 = t143 * qJD(1) + t144 * qJD(2);
	t132 = t145 * qJD(1) + t142 * qJD(2);
	t1 = [t135 * r_i_i_C(1) + t134 * r_i_i_C(2) + (-pkin(1) * t141 - t139 * t151) * qJD(1), t132 * r_i_i_C(1) + t133 * r_i_i_C(2), 0, 0, 0, 0; -t133 * r_i_i_C(1) + t132 * r_i_i_C(2) + (-pkin(1) * t139 + t141 * t151) * qJD(1), -t134 * r_i_i_C(1) + t135 * r_i_i_C(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t138 - r_i_i_C(2) * t140) * t136 * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:26
	% EndTime: 2019-10-10 11:31:26
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (145->51), mult. (441->91), div. (0->0), fcn. (412->8), ass. (0->38)
	t242 = cos(pkin(6));
	t244 = sin(qJ(2));
	t245 = sin(qJ(1));
	t264 = t245 * t244;
	t258 = t242 * t264;
	t247 = cos(qJ(2));
	t248 = cos(qJ(1));
	t261 = t248 * t247;
	t233 = -qJD(1) * t258 - qJD(2) * t264 + (qJD(2) * t242 + qJD(1)) * t261;
	t241 = sin(pkin(6));
	t265 = t241 * t248;
	t256 = qJD(3) * t265;
	t269 = t233 - t256;
	t268 = -r_i_i_C(3) - pkin(9);
	t267 = t241 * t245;
	t246 = cos(qJ(3));
	t266 = t241 * t246;
	t263 = t245 * t247;
	t262 = t248 * t244;
	t260 = qJD(1) * t241;
	t235 = t242 * t262 + t263;
	t259 = qJD(3) * t235;
	t257 = t248 * t260;
	t243 = sin(qJ(3));
	t255 = -r_i_i_C(1) * t243 - r_i_i_C(2) * t246;
	t254 = t246 * r_i_i_C(1) - t243 * r_i_i_C(2) + pkin(2);
	t253 = t242 * t261 - t264;
	t252 = t242 * t263 + t262;
	t251 = t258 - t261;
	t250 = qJD(3) * t255;
	t249 = t245 * t260 - t259;
	t238 = t246 * t256;
	t232 = t252 * qJD(1) + t235 * qJD(2);
	t231 = t235 * qJD(1) + t252 * qJD(2);
	t230 = -t253 * qJD(1) + t251 * qJD(2);
	t229 = t243 * t257 - t231 * t246 + (t243 * t251 + t245 * t266) * qJD(3);
	t228 = t246 * t257 + t231 * t243 + (-t243 * t267 + t246 * t251) * qJD(3);
	t1 = [(-t233 * t246 + t243 * t259 + t238) * r_i_i_C(1) + (t269 * t243 + t246 * t259) * r_i_i_C(2) - t233 * pkin(2) + t268 * t232 + (-t248 * pkin(1) + (-pkin(8) + t255) * t267) * qJD(1), t254 * t230 + t268 * t231 - t252 * t250, t228 * r_i_i_C(1) - t229 * r_i_i_C(2), 0, 0, 0; -t231 * pkin(2) + t229 * r_i_i_C(1) + t228 * r_i_i_C(2) + t268 * t230 + (-pkin(1) * t245 + pkin(8) * t265) * qJD(1), -t254 * t232 - t268 * t233 + t253 * t250, t238 * r_i_i_C(2) + (t249 * r_i_i_C(1) - t233 * r_i_i_C(2)) * t246 + (-t269 * r_i_i_C(1) - t249 * r_i_i_C(2)) * t243, 0, 0, 0; 0, (t247 * t250 + (-t254 * t244 - t268 * t247) * qJD(2)) * t241, t255 * t247 * t241 * qJD(2) + ((-t242 * t243 - t244 * t266) * r_i_i_C(1) + (t241 * t243 * t244 - t242 * t246) * r_i_i_C(2)) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:27
	% EndTime: 2019-10-10 11:31:27
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (344->68), mult. (1036->114), div. (0->0), fcn. (1016->10), ass. (0->47)
	t325 = cos(pkin(6));
	t341 = qJD(2) * t325 + qJD(1);
	t327 = sin(qJ(2));
	t328 = sin(qJ(1));
	t355 = t328 * t327;
	t348 = t325 * t355;
	t330 = cos(qJ(2));
	t331 = cos(qJ(1));
	t352 = t331 * t330;
	t307 = -qJD(1) * t348 - qJD(2) * t355 + t341 * t352;
	t323 = sin(pkin(6));
	t356 = t323 * t331;
	t363 = -qJD(3) * t356 + t307;
	t353 = t331 * t327;
	t354 = t328 * t330;
	t311 = t325 * t353 + t354;
	t351 = qJD(1) * t323;
	t362 = -qJD(3) * t311 + t328 * t351;
	t326 = sin(qJ(3));
	t329 = cos(qJ(3));
	t361 = t362 * t326 + t363 * t329;
	t322 = sin(pkin(11));
	t324 = cos(pkin(11));
	t340 = t324 * r_i_i_C(1) - t322 * r_i_i_C(2) + pkin(3);
	t359 = r_i_i_C(3) + qJ(4);
	t360 = t359 * t326 + t340 * t329 + pkin(2);
	t358 = t323 * t328;
	t357 = t323 * t329;
	t350 = qJD(2) * t330;
	t347 = t325 * t352;
	t345 = t331 * t351;
	t344 = t323 * t350;
	t339 = t322 * r_i_i_C(1) + t324 * r_i_i_C(2) + pkin(9);
	t313 = -t348 + t352;
	t338 = -t313 * t326 + t328 * t357;
	t337 = t313 * t329 + t326 * t358;
	t336 = t325 * t326 + t327 * t357;
	t335 = t325 * t354 + t353;
	t300 = t363 * t326 - t362 * t329;
	t332 = t326 * qJD(4) + (-t340 * t326 + t359 * t329) * qJD(3);
	t308 = t336 * qJD(3) + t326 * t344;
	t306 = t335 * qJD(1) + t311 * qJD(2);
	t305 = t311 * qJD(1) + t335 * qJD(2);
	t304 = -qJD(1) * t347 - t331 * t350 + t341 * t355;
	t299 = t338 * qJD(3) - t305 * t329 + t326 * t345;
	t298 = t337 * qJD(3) - t305 * t326 - t329 * t345;
	t1 = [(-t306 * t322 - t324 * t361) * r_i_i_C(1) + (-t306 * t324 + t322 * t361) * r_i_i_C(2) - t361 * pkin(3) - (t311 * t326 + t329 * t356) * qJD(4) - t307 * pkin(2) - t306 * pkin(9) - t359 * t300 + (-t331 * pkin(1) - pkin(8) * t358) * qJD(1), t360 * t304 - t339 * t305 - t332 * t335, t337 * qJD(4) - t340 * t298 + t359 * t299, t298, 0, 0; (t299 * t324 - t304 * t322) * r_i_i_C(1) + (-t299 * t322 - t304 * t324) * r_i_i_C(2) + t299 * pkin(3) - t338 * qJD(4) - t305 * pkin(2) - t304 * pkin(9) + t359 * t298 + (-t328 * pkin(1) + pkin(8) * t356) * qJD(1), t339 * t307 + t332 * (t347 - t355) - t360 * t306, -(-t311 * t329 + t326 * t356) * qJD(4) + t359 * t361 - t340 * t300, t300, 0, 0; 0, (t332 * t330 + (-t327 * t360 + t339 * t330) * qJD(2)) * t323, t336 * qJD(4) + t359 * (t329 * t344 + (-t323 * t326 * t327 + t325 * t329) * qJD(3)) - t340 * t308, t308, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:28
	% EndTime: 2019-10-10 11:31:28
	% DurationCPUTime: 0.68s
	% Computational Cost: add. (529->96), mult. (1587->159), div. (0->0), fcn. (1585->10), ass. (0->62)
	t370 = sin(qJ(3));
	t373 = cos(qJ(3));
	t369 = cos(pkin(6));
	t371 = sin(qJ(2));
	t375 = cos(qJ(1));
	t405 = t375 * t371;
	t372 = sin(qJ(1));
	t374 = cos(qJ(2));
	t406 = t372 * t374;
	t354 = t369 * t405 + t406;
	t367 = sin(pkin(6));
	t403 = qJD(1) * t367;
	t418 = -qJD(3) * t354 + t372 * t403;
	t386 = qJD(2) * t369 + qJD(1);
	t407 = t372 * t371;
	t396 = t369 * t407;
	t402 = qJD(2) * t371;
	t404 = t375 * t374;
	t348 = -qJD(1) * t396 - t372 * t402 + t386 * t404;
	t409 = t367 * t375;
	t419 = -qJD(3) * t409 + t348;
	t340 = t418 * t370 + t419 * t373;
	t355 = t369 * t406 + t405;
	t347 = t355 * qJD(1) + t354 * qJD(2);
	t366 = sin(pkin(11));
	t368 = cos(pkin(11));
	t420 = t340 * t366 - t347 * t368;
	t398 = qJD(5) * t366;
	t413 = r_i_i_C(2) + qJ(4);
	t376 = (-pkin(3) * t370 + t413 * t373) * qJD(3) + t370 * qJD(4) + t373 * t398;
	t417 = t373 * pkin(3) + t413 * t370 + pkin(2);
	t415 = pkin(4) + r_i_i_C(1);
	t412 = r_i_i_C(3) + qJ(5);
	t410 = t367 * t372;
	t408 = t371 * t373;
	t401 = qJD(2) * t374;
	t399 = qJD(3) * t370;
	t397 = t368 * qJD(5);
	t395 = t369 * t404;
	t392 = t375 * t403;
	t391 = t367 * t401;
	t389 = t374 * t399;
	t387 = t354 * t373 - t370 * t409;
	t385 = t354 * t370 + t373 * t409;
	t356 = -t396 + t404;
	t351 = -t356 * t370 + t373 * t410;
	t352 = t356 * t373 + t370 * t410;
	t384 = t367 * t408 + t369 * t370;
	t383 = -t367 * t371 * t370 + t369 * t373;
	t345 = -qJD(1) * t395 - t375 * t401 + t386 * t407;
	t382 = t345 * t373 + t355 * t399;
	t353 = t395 - t407;
	t381 = -t347 * t373 - t353 * t399;
	t377 = -t412 * t366 - t415 * t368 - pkin(3);
	t339 = t419 * t370 - t418 * t373;
	t350 = t383 * qJD(3) + t373 * t391;
	t349 = t384 * qJD(3) + t370 * t391;
	t346 = t354 * qJD(1) + t355 * qJD(2);
	t338 = t351 * qJD(3) - t346 * t373 + t370 * t392;
	t337 = t352 * qJD(3) - t346 * t370 - t373 * t392;
	t329 = t338 * t366 + t345 * t368;
	t1 = [-(t353 * t368 + t387 * t366) * qJD(5) - t340 * pkin(3) - t385 * qJD(4) - t348 * pkin(2) - t347 * pkin(9) - t413 * t339 + t415 * (-t340 * t368 - t347 * t366) - t412 * t420 + (-t375 * pkin(1) - pkin(8) * t410) * qJD(1), -t356 * t397 - t346 * pkin(9) + t415 * (-t346 * t366 + t382 * t368) + t412 * (t346 * t368 + t382 * t366) - t376 * t355 + t417 * t345, t352 * qJD(4) + t377 * t337 + t413 * t338 + t351 * t398, t337, t329, 0; -(-t352 * t366 + t355 * t368) * qJD(5) + t338 * pkin(3) - t351 * qJD(4) - t346 * pkin(2) - t345 * pkin(9) + t413 * t337 + t415 * (t338 * t368 - t345 * t366) + t412 * t329 + (-t372 * pkin(1) + pkin(8) * t409) * qJD(1), -t354 * t397 + t348 * pkin(9) + t415 * (t348 * t366 + t381 * t368) + t412 * (-t348 * t368 + t381 * t366) + t376 * t353 - t417 * t347, t387 * qJD(4) + t377 * t339 + t413 * t340 - t385 * t398, t339, t420, 0; 0, (t415 * (-t368 * t389 + (t366 * t374 - t368 * t408) * qJD(2)) - t412 * (t366 * t389 + (t366 * t408 + t368 * t374) * qJD(2)) - t371 * t397 + t376 * t374 + (t374 * pkin(9) - t371 * t417) * qJD(2)) * t367, t384 * qJD(4) + t377 * t349 + t413 * t350 + t383 * t398, t349, -t367 * t368 * t402 + t350 * t366, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:28
	% EndTime: 2019-10-10 11:31:29
	% DurationCPUTime: 1.50s
	% Computational Cost: add. (1049->158), mult. (3138->268), div. (0->0), fcn. (3263->12), ass. (0->94)
	t471 = cos(pkin(6));
	t474 = sin(qJ(2));
	t523 = cos(qJ(1));
	t499 = t523 * t474;
	t475 = sin(qJ(1));
	t478 = cos(qJ(2));
	t510 = t475 * t478;
	t454 = t471 * t499 + t510;
	t473 = sin(qJ(3));
	t477 = cos(qJ(3));
	t469 = sin(pkin(6));
	t500 = t469 * t523;
	t443 = t454 * t477 - t473 * t500;
	t498 = t523 * t478;
	t489 = t471 * t498;
	t511 = t475 * t474;
	t453 = -t489 + t511;
	t468 = sin(pkin(11));
	t470 = cos(pkin(11));
	t420 = t443 * t468 - t453 * t470;
	t421 = t443 * t470 + t453 * t468;
	t472 = sin(qJ(6));
	t476 = cos(qJ(6));
	t530 = -t420 * t476 + t421 * t472;
	t529 = t420 * t472 + t421 * t476;
	t492 = t523 * qJD(2);
	t493 = t523 * qJD(1);
	t502 = t471 * t511;
	t508 = qJD(2) * t474;
	t437 = -qJD(1) * t502 - t475 * t508 + (t471 * t492 + t493) * t478;
	t528 = -qJD(3) * t500 + t437;
	t504 = -r_i_i_C(3) - pkin(10) + qJ(4);
	t527 = (-pkin(3) * t473 + t504 * t477) * qJD(3) + t473 * qJD(4);
	t526 = t477 * pkin(3) + t504 * t473 + pkin(2);
	t524 = pkin(4) + pkin(5);
	t455 = t471 * t510 + t499;
	t436 = t455 * qJD(1) + t454 * qJD(2);
	t521 = t436 * t468;
	t520 = t436 * t470;
	t517 = t468 * t477;
	t516 = t469 * t475;
	t515 = t469 * t478;
	t514 = t470 * t477;
	t513 = t470 * t478;
	t512 = t474 * t477;
	t509 = t477 * t478;
	t507 = qJD(3) * t473;
	t506 = qJD(3) * t477;
	t503 = t469 * t474 * t473;
	t497 = qJD(1) * t516;
	t496 = t469 * t508;
	t495 = qJD(2) * t515;
	t494 = t478 * t507;
	t491 = t528 * t477;
	t490 = t477 * t500;
	t487 = t472 * r_i_i_C(1) + t476 * r_i_i_C(2) + qJ(5);
	t486 = t476 * r_i_i_C(1) - t472 * r_i_i_C(2) + t524;
	t456 = t498 - t502;
	t446 = t456 * t477 + t473 * t516;
	t452 = t469 * t512 + t471 * t473;
	t434 = -qJD(1) * t489 - t478 * t492 + (qJD(2) * t471 + qJD(1)) * t511;
	t485 = t434 * t477 + t455 * t507;
	t484 = -t436 * t477 + t453 * t507;
	t483 = t454 * t473 + t490;
	t416 = t454 * t506 + t528 * t473 - t477 * t497;
	t480 = -t487 * t468 - t486 * t470 - pkin(3);
	t479 = t468 * qJD(5) + ((t468 * t476 - t470 * t472) * r_i_i_C(1) + (-t468 * t472 - t470 * t476) * r_i_i_C(2)) * qJD(6);
	t448 = (t468 * t474 + t470 * t509) * t469;
	t447 = (t468 * t509 - t470 * t474) * t469;
	t445 = -t456 * t473 + t477 * t516;
	t441 = -qJD(3) * t503 + (qJD(3) * t471 + t495) * t477;
	t440 = t452 * qJD(3) + t473 * t495;
	t439 = t452 * t470 - t468 * t515;
	t438 = t452 * t468 + t469 * t513;
	t435 = t454 * qJD(1) + t455 * qJD(2);
	t431 = -t455 * t514 + t456 * t468;
	t430 = -t455 * t517 - t456 * t470;
	t429 = -t453 * t514 + t454 * t468;
	t428 = -t453 * t517 - t454 * t470;
	t427 = t441 * t470 + t468 * t496;
	t426 = t441 * t468 - t470 * t496;
	t425 = t446 * t470 + t455 * t468;
	t424 = t446 * t468 - t455 * t470;
	t419 = (qJD(3) * t454 - t497) * t473 - t491;
	t417 = -t454 * t507 + t473 * t497 + t491;
	t415 = -t435 * t477 - t456 * t507 + (t473 * t493 + t475 * t506) * t469;
	t414 = -qJD(1) * t490 + t446 * qJD(3) - t435 * t473;
	t407 = t417 * t470 + t521;
	t406 = t417 * t468 - t520;
	t405 = t415 * t470 - t434 * t468;
	t404 = t415 * t468 + t434 * t470;
	t403 = t404 * t472 + t405 * t476 + (t424 * t476 - t425 * t472) * qJD(6);
	t402 = t404 * t476 - t405 * t472 + (-t424 * t472 - t425 * t476) * qJD(6);
	t1 = [-t420 * qJD(5) + t419 * pkin(3) - t483 * qJD(4) - t437 * pkin(2) - t436 * pkin(9) + t487 * (t419 * t468 + t520) + t486 * (t419 * t470 - t521) + (t530 * r_i_i_C(1) + t529 * r_i_i_C(2)) * qJD(6) + (-t523 * pkin(1) - pkin(8) * t516) * qJD(1) - t504 * t416, -t435 * pkin(9) + t430 * qJD(5) + t487 * (t435 * t470 + t485 * t468) + t486 * (-t435 * t468 + t485 * t470) + ((t430 * t476 - t431 * t472) * r_i_i_C(1) + (-t430 * t472 - t431 * t476) * r_i_i_C(2)) * qJD(6) - t527 * t455 + t526 * t434, t446 * qJD(4) + t480 * t414 + t504 * t415 + t479 * t445, t414, t404, t402 * r_i_i_C(1) - t403 * r_i_i_C(2); -t435 * pkin(2) + t415 * pkin(3) - t434 * pkin(9) + t403 * r_i_i_C(1) + t402 * r_i_i_C(2) + t404 * qJ(5) - t445 * qJD(4) + t424 * qJD(5) + t524 * t405 + (-pkin(1) * t475 + pkin(8) * t500) * qJD(1) + t504 * t414, t437 * pkin(9) + t428 * qJD(5) + t487 * (-t437 * t470 + t484 * t468) + t486 * (t437 * t468 + t484 * t470) + ((t428 * t476 - t429 * t472) * r_i_i_C(1) + (-t428 * t472 - t429 * t476) * r_i_i_C(2)) * qJD(6) - t527 * t453 - t526 * t436, t443 * qJD(4) + t480 * t416 + t504 * t417 - t479 * t483, t416, t406, (t406 * t476 - t407 * t472) * r_i_i_C(1) + (-t406 * t472 - t407 * t476) * r_i_i_C(2) + (-t529 * r_i_i_C(1) + t530 * r_i_i_C(2)) * qJD(6); 0, t447 * qJD(5) + ((t447 * t476 - t448 * t472) * r_i_i_C(1) + (-t447 * t472 - t448 * t476) * r_i_i_C(2)) * qJD(6) + (-t487 * (t468 * t494 + (t468 * t512 + t513) * qJD(2)) + t486 * (-t470 * t494 + (t468 * t478 - t470 * t512) * qJD(2)) + t527 * t478 + (pkin(9) * t478 - t474 * t526) * qJD(2)) * t469, t452 * qJD(4) + t504 * t441 + t479 * (t471 * t477 - t503) + t480 * t440, t440, t426, (t426 * t476 - t427 * t472) * r_i_i_C(1) + (-t426 * t472 - t427 * t476) * r_i_i_C(2) + ((-t438 * t472 - t439 * t476) * r_i_i_C(1) + (-t438 * t476 + t439 * t472) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end