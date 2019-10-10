% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:49
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:49:37
	% EndTime: 2019-10-10 11:49:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:49:37
	% EndTime: 2019-10-10 11:49:37
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
	% StartTime: 2019-10-10 11:49:38
	% EndTime: 2019-10-10 11:49:38
	% DurationCPUTime: 0.07s
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
	% StartTime: 2019-10-10 11:49:39
	% EndTime: 2019-10-10 11:49:39
	% DurationCPUTime: 0.23s
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
	% StartTime: 2019-10-10 11:49:40
	% EndTime: 2019-10-10 11:49:40
	% DurationCPUTime: 0.42s
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
	% StartTime: 2019-10-10 11:49:40
	% EndTime: 2019-10-10 11:49:41
	% DurationCPUTime: 0.81s
	% Computational Cost: add. (600->95), mult. (1491->158), div. (0->0), fcn. (1502->12), ass. (0->66)
	t386 = sin(qJ(1));
	t382 = cos(pkin(6));
	t401 = qJD(2) * t382 + qJD(1);
	t385 = sin(qJ(2));
	t419 = t386 * t385;
	t410 = t382 * t419;
	t414 = qJD(2) * t385;
	t388 = cos(qJ(2));
	t389 = cos(qJ(1));
	t416 = t389 * t388;
	t352 = -qJD(1) * t410 - t386 * t414 + t401 * t416;
	t381 = sin(pkin(6));
	t420 = t381 * t389;
	t433 = -qJD(3) * t420 + t352;
	t417 = t389 * t385;
	t418 = t386 * t388;
	t363 = t382 * t417 + t418;
	t415 = qJD(1) * t381;
	t432 = -qJD(3) * t363 + t386 * t415;
	t384 = sin(qJ(3));
	t387 = cos(qJ(3));
	t356 = t363 * t387 - t384 * t420;
	t409 = t382 * t416;
	t362 = -t409 + t419;
	t379 = pkin(11) + qJ(5);
	t377 = sin(t379);
	t378 = cos(t379);
	t431 = t356 * t377 - t362 * t378;
	t430 = t356 * t378 + t362 * t377;
	t399 = t377 * r_i_i_C(1) + t378 * r_i_i_C(2);
	t395 = qJD(5) * t399;
	t376 = cos(pkin(11)) * pkin(4) + pkin(3);
	t400 = t378 * r_i_i_C(1) - t377 * r_i_i_C(2);
	t398 = t376 + t400;
	t426 = r_i_i_C(3) + pkin(10) + qJ(4);
	t390 = (t398 * t384 - t426 * t387) * qJD(3) - t384 * qJD(4) + t387 * t395;
	t346 = t432 * t384 + t433 * t387;
	t427 = t426 * t384 + t398 * t387 + pkin(2);
	t423 = t381 * t386;
	t422 = t381 * t387;
	t421 = t381 * t388;
	t413 = qJD(2) * t388;
	t408 = -sin(pkin(11)) * pkin(4) - pkin(9);
	t406 = t389 * t415;
	t405 = t381 * t414;
	t404 = t381 * t413;
	t397 = t363 * t384 + t387 * t420;
	t365 = -t410 + t416;
	t358 = -t365 * t384 + t386 * t422;
	t359 = t365 * t387 + t384 * t423;
	t361 = t382 * t384 + t385 * t422;
	t396 = -t381 * t385 * t384 + t382 * t387;
	t364 = t382 * t418 + t417;
	t394 = t400 * qJD(5);
	t392 = t399 - t408;
	t345 = t433 * t384 - t432 * t387;
	t354 = t396 * qJD(3) + t387 * t404;
	t353 = t361 * qJD(3) + t384 * t404;
	t351 = t364 * qJD(1) + t363 * qJD(2);
	t350 = t363 * qJD(1) + t364 * qJD(2);
	t349 = -qJD(1) * t409 - t389 * t413 + t401 * t419;
	t344 = t358 * qJD(3) - t350 * t387 + t384 * t406;
	t343 = t359 * qJD(3) - t350 * t384 - t387 * t406;
	t342 = t344 * t378 - t349 * t377 + (-t359 * t377 + t364 * t378) * qJD(5);
	t341 = -t344 * t377 - t349 * t378 + (-t359 * t378 - t364 * t377) * qJD(5);
	t1 = [-t397 * qJD(4) - t352 * pkin(2) - t398 * t346 - t392 * t351 - t426 * t345 + (t431 * r_i_i_C(1) + t430 * r_i_i_C(2)) * qJD(5) + (-t389 * pkin(1) - pkin(8) * t423) * qJD(1), t427 * t349 - t392 * t350 + t390 * t364 + t365 * t394, t359 * qJD(4) - t398 * t343 + t426 * t344 - t358 * t395, t343, t341 * r_i_i_C(1) - t342 * r_i_i_C(2), 0; -t350 * pkin(2) + t342 * r_i_i_C(1) + t341 * r_i_i_C(2) - t358 * qJD(4) + t344 * t376 + t408 * t349 + t426 * t343 + (-pkin(1) * t386 + pkin(8) * t420) * qJD(1), -t351 * t427 + t392 * t352 + t390 * t362 + t363 * t394, t356 * qJD(4) - t398 * t345 + t426 * t346 + t397 * t395, t345, (-t346 * t377 + t351 * t378) * r_i_i_C(1) + (-t346 * t378 - t351 * t377) * r_i_i_C(2) + (-t430 * r_i_i_C(1) + t431 * r_i_i_C(2)) * qJD(5), 0; 0, ((-qJD(2) * t427 + t394) * t385 + (t392 * qJD(2) - t390) * t388) * t381, t361 * qJD(4) - t398 * t353 + t426 * t354 - t396 * t395, t353, (-t354 * t377 + t378 * t405) * r_i_i_C(1) + (-t354 * t378 - t377 * t405) * r_i_i_C(2) + ((-t361 * t378 + t377 * t421) * r_i_i_C(1) + (t361 * t377 + t378 * t421) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:49:42
	% EndTime: 2019-10-10 11:49:43
	% DurationCPUTime: 0.87s
	% Computational Cost: add. (1062->122), mult. (2499->192), div. (0->0), fcn. (2581->12), ass. (0->77)
	t472 = sin(qJ(3));
	t475 = cos(qJ(3));
	t470 = cos(pkin(6));
	t473 = sin(qJ(2));
	t477 = cos(qJ(1));
	t516 = t477 * t473;
	t474 = sin(qJ(1));
	t476 = cos(qJ(2));
	t517 = t474 * t476;
	t451 = t470 * t516 + t517;
	t469 = sin(pkin(6));
	t514 = qJD(1) * t469;
	t539 = -qJD(3) * t451 + t474 * t514;
	t518 = t474 * t473;
	t506 = t470 * t518;
	t513 = qJD(2) * t473;
	t515 = t477 * t476;
	t440 = -qJD(1) * t506 - t474 * t513 + (qJD(2) * t470 + qJD(1)) * t515;
	t519 = t469 * t477;
	t540 = -qJD(3) * t519 + t440;
	t430 = t539 * t472 + t540 * t475;
	t452 = t470 * t517 + t516;
	t439 = t452 * qJD(1) + t451 * qJD(2);
	t467 = pkin(11) + qJ(5);
	t465 = sin(t467);
	t466 = cos(t467);
	t444 = t451 * t475 - t472 * t519;
	t450 = -t470 * t515 + t518;
	t493 = t444 * t466 + t450 * t465;
	t417 = t493 * qJD(5) + t430 * t465 - t439 * t466;
	t492 = t444 * t465 - t450 * t466;
	t541 = -t492 * qJD(5) + t430 * t466 + t439 * t465;
	t464 = cos(pkin(11)) * pkin(4) + pkin(3);
	t509 = qJD(6) * t465;
	t530 = r_i_i_C(2) + pkin(10) + qJ(4);
	t536 = (-t464 * t472 + t530 * t475) * qJD(3) + t472 * qJD(4) + t475 * t509;
	t521 = t469 * t475;
	t449 = t470 * t472 + t473 * t521;
	t520 = t469 * t476;
	t535 = -t449 * t466 + t465 * t520;
	t534 = t475 * t464 + t530 * t472 + pkin(2);
	t511 = qJD(3) * t472;
	t533 = (qJD(2) * t475 - qJD(5)) * t473 + t476 * t511;
	t531 = r_i_i_C(1) + pkin(5);
	t529 = r_i_i_C(3) + qJ(6);
	t522 = t469 * t474;
	t510 = qJD(5) * t475;
	t508 = t466 * qJD(6);
	t505 = sin(pkin(11)) * pkin(4) + pkin(9);
	t503 = t477 * t514;
	t502 = t469 * t513;
	t501 = qJD(2) * t520;
	t438 = t451 * qJD(1) + t452 * qJD(2);
	t495 = t452 * t510 - t438;
	t494 = t450 * t510 + t440;
	t486 = t506 - t515;
	t447 = t472 * t522 - t475 * t486;
	t491 = t447 * t466 + t452 * t465;
	t490 = -t447 * t465 + t452 * t466;
	t489 = (qJD(2) - t510) * t476;
	t488 = t451 * t472 + t475 * t519;
	t446 = t472 * t486 + t474 * t521;
	t487 = -t469 * t473 * t472 + t470 * t475;
	t429 = t540 * t472 - t539 * t475;
	t482 = -t529 * t465 - t531 * t466 - t464;
	t437 = t450 * qJD(1) + t486 * qJD(2);
	t481 = -qJD(5) * t486 + t437 * t475 + t452 * t511;
	t480 = qJD(5) * t451 - t439 * t475 + t450 * t511;
	t479 = t509 + (-t531 * t465 + t529 * t466) * qJD(5);
	t442 = t487 * qJD(3) + t475 * t501;
	t441 = t449 * qJD(3) + t472 * t501;
	t428 = t446 * qJD(3) - t438 * t475 + t472 * t503;
	t427 = t447 * qJD(3) - t438 * t472 - t475 * t503;
	t425 = -t535 * qJD(5) + t442 * t465 - t466 * t502;
	t416 = t490 * qJD(5) + t428 * t466 - t437 * t465;
	t415 = t491 * qJD(5) + t428 * t465 + t437 * t466;
	t1 = [-t492 * qJD(6) - t430 * t464 - t488 * qJD(4) - t440 * pkin(2) - t505 * t439 - t530 * t429 - t531 * t541 - t529 * t417 + (-t477 * pkin(1) - pkin(8) * t522) * qJD(1), t486 * t508 - t505 * t438 + t531 * (t495 * t465 + t481 * t466) + t529 * (t481 * t465 - t495 * t466) - t536 * t452 + t534 * t437, t447 * qJD(4) + t482 * t427 + t530 * t428 + t479 * t446, t427, t491 * qJD(6) - t531 * t415 + t529 * t416, t415; -t490 * qJD(6) + t428 * t464 - t446 * qJD(4) - t438 * pkin(2) - t505 * t437 + t530 * t427 + t531 * t416 + t529 * t415 + (-t474 * pkin(1) + pkin(8) * t519) * qJD(1), -t451 * t508 + t505 * t440 + t531 * (t494 * t465 + t480 * t466) + t529 * (t480 * t465 - t494 * t466) - t536 * t450 - t534 * t439, t444 * qJD(4) + t482 * t429 + t530 * t430 - t479 * t488, t429, t493 * qJD(6) - t531 * t417 + t529 * t541, t417; 0, (t531 * (t465 * t489 - t533 * t466) - t529 * (t533 * t465 + t466 * t489) - t473 * t508 + t536 * t476 + (-t473 * t534 + t505 * t476) * qJD(2)) * t469, t449 * qJD(4) + t482 * t441 + t530 * t442 + t479 * t487, t441, -t535 * qJD(6) + t529 * (t465 * t502 + t442 * t466 + (-t449 * t465 - t466 * t520) * qJD(5)) - t531 * t425, t425;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end