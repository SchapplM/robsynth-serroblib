% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S7RRRRRRR1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% JaD_transl [3x7]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 17:10
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S7RRRRRRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),uint8(0),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_sym_varpar: qJD has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S7RRRRRRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_sym_varpar: pkin has to be [4x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (15->13), mult. (56->29), div. (0->0), fcn. (36->4), ass. (0->12)
	t16 = sin(qJ(1));
	t25 = qJD(1) * t16;
	t18 = cos(qJ(1));
	t24 = qJD(1) * t18;
	t23 = qJD(2) * t16;
	t22 = qJD(2) * t18;
	t15 = sin(qJ(2));
	t17 = cos(qJ(2));
	t21 = -r_i_i_C(1) * t17 + r_i_i_C(2) * t15;
	t20 = r_i_i_C(1) * t15 + r_i_i_C(2) * t17;
	t19 = t20 * qJD(2);
	t1 = [t20 * t23 + (-r_i_i_C(3) * t16 + t21 * t18) * qJD(1), (t15 * t22 + t17 * t25) * r_i_i_C(2) + (t15 * t25 - t17 * t22) * r_i_i_C(1), 0, 0, 0, 0, 0; -t18 * t19 + (r_i_i_C(3) * t18 + t21 * t16) * qJD(1), (t15 * t23 - t17 * t24) * r_i_i_C(2) + (-t15 * t24 - t17 * t23) * r_i_i_C(1), 0, 0, 0, 0, 0; 0, -t19, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:05
	% EndTime: 2019-10-10 17:10:05
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (70->32), mult. (236->67), div. (0->0), fcn. (190->6), ass. (0->28)
	t217 = pkin(2) + r_i_i_C(3);
	t196 = cos(qJ(3));
	t198 = cos(qJ(1));
	t216 = t196 * t198;
	t195 = sin(qJ(1));
	t215 = qJD(1) * t195;
	t214 = qJD(1) * t198;
	t194 = sin(qJ(2));
	t213 = qJD(2) * t194;
	t197 = cos(qJ(2));
	t212 = qJD(2) * t197;
	t211 = qJD(2) * t198;
	t210 = qJD(3) * t194;
	t209 = qJD(3) * t197;
	t208 = qJD(1) + t209;
	t207 = qJD(1) * t197 + qJD(3);
	t193 = sin(qJ(3));
	t206 = r_i_i_C(1) * t196 - r_i_i_C(2) * t193;
	t205 = r_i_i_C(1) * t193 + r_i_i_C(2) * t196;
	t204 = t208 * t193;
	t203 = qJD(2) * t206;
	t200 = t194 * t211 + t207 * t195;
	t199 = t217 * qJD(2) + t205 * qJD(3);
	t192 = -t207 * t216 + (t196 * t213 + t204) * t195;
	t191 = t208 * t196 * t195 + (-t195 * t213 + t207 * t198) * t193;
	t190 = t200 * t196 + t198 * t204;
	t189 = t200 * t193 - t208 * t216;
	t1 = [t192 * r_i_i_C(1) + t191 * r_i_i_C(2) + t217 * (t194 * t214 + t195 * t212), (-t198 * t203 + t217 * t215) * t197 + (t199 * t198 + t206 * t215) * t194, r_i_i_C(1) * t189 + r_i_i_C(2) * t190, 0, 0, 0, 0; -t190 * r_i_i_C(1) + t189 * r_i_i_C(2) + t217 * (t194 * t215 - t197 * t211), (-t195 * t203 - t217 * t214) * t197 + (t199 * t195 - t206 * t214) * t194, -r_i_i_C(1) * t191 + r_i_i_C(2) * t192, 0, 0, 0, 0; 0, -t205 * t209 + (-t206 * t194 - t217 * t197) * qJD(2), (t193 * t210 - t196 * t212) * r_i_i_C(2) + (-t193 * t212 - t196 * t210) * r_i_i_C(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:06
	% EndTime: 2019-10-10 17:10:07
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (190->67), mult. (618->126), div. (0->0), fcn. (577->8), ass. (0->57)
	t301 = cos(qJ(3));
	t348 = -qJD(4) * t301 + qJD(2);
	t297 = sin(qJ(3));
	t303 = cos(qJ(1));
	t338 = t303 * t297;
	t299 = sin(qJ(1));
	t302 = cos(qJ(2));
	t339 = t299 * t302;
	t289 = t301 * t339 + t338;
	t298 = sin(qJ(2));
	t334 = qJD(2) * t302;
	t335 = qJD(1) * t303;
	t313 = t298 * t335 + t299 * t334;
	t308 = -qJD(4) * t289 + t313;
	t324 = qJD(1) * t302 + qJD(3);
	t326 = t299 * qJD(2) * t298;
	t327 = t297 * t339;
	t336 = qJD(1) * t299;
	t337 = t303 * t301;
	t287 = -qJD(3) * t327 - t297 * t336 - t301 * t326 + t324 * t337;
	t330 = qJD(4) * t298;
	t320 = t299 * t330 + t287;
	t347 = t308 * r_i_i_C(1) - t320 * r_i_i_C(2);
	t296 = sin(qJ(4));
	t300 = cos(qJ(4));
	t332 = qJD(3) * t297;
	t309 = t296 * t332 + t348 * t300;
	t310 = -t348 * t296 + t300 * t332;
	t341 = qJD(3) * r_i_i_C(3);
	t346 = qJD(2) * pkin(2) + t310 * r_i_i_C(1) - t309 * r_i_i_C(2) + t301 * t341;
	t340 = t298 * t301;
	t342 = t297 * r_i_i_C(3);
	t345 = (-t296 * t302 + t300 * t340) * r_i_i_C(1) - (t296 * t340 + t300 * t302) * r_i_i_C(2) + t302 * pkin(2) - t298 * t342;
	t344 = t296 * r_i_i_C(1);
	t343 = t296 * r_i_i_C(2);
	t333 = qJD(2) * t303;
	t331 = qJD(4) * t296;
	t329 = qJD(4) * t300;
	t325 = qJD(3) * t302 + qJD(1);
	t323 = qJD(2) * t301 - qJD(4);
	t322 = -t300 * r_i_i_C(1) + t343;
	t311 = t298 * t333 + t324 * t299;
	t285 = t311 * t301 + t325 * t338;
	t321 = t303 * t330 - t285;
	t319 = t323 * t300;
	t314 = t299 * t301 + t302 * t338;
	t312 = t298 * t336 - t302 * t333;
	t307 = -qJD(4) * (-t299 * t297 + t302 * t337) - t312;
	t306 = -r_i_i_C(1) * t319 + qJD(2) * t342 + t323 * t343;
	t305 = -t320 * r_i_i_C(1) - t308 * r_i_i_C(2);
	t304 = t346 * t298 + t306 * t302;
	t288 = -t327 + t337;
	t286 = t314 * qJD(1) + t289 * qJD(3) - t297 * t326;
	t284 = t311 * t297 - t325 * t337;
	t283 = t307 * t296 + t321 * t300;
	t282 = -t321 * t296 + t307 * t300;
	t1 = [t313 * pkin(2) + t286 * r_i_i_C(3) - t347 * t296 + t305 * t300, t304 * t303 + t345 * t336, t285 * r_i_i_C(3) + (-t284 * t296 + t314 * t329) * r_i_i_C(2) + (t284 * t300 + t314 * t331) * r_i_i_C(1), t282 * r_i_i_C(1) - t283 * r_i_i_C(2), 0, 0, 0; t312 * pkin(2) + t283 * r_i_i_C(1) + t282 * r_i_i_C(2) + t284 * r_i_i_C(3), t304 * t299 - t345 * t335, -t287 * r_i_i_C(3) + (t286 * t296 - t288 * t329) * r_i_i_C(2) + (-t286 * t300 - t288 * t331) * r_i_i_C(1), t305 * t296 + t347 * t300, 0, 0, 0; 0, t306 * t298 - t346 * t302, (t322 * t298 * qJD(3) - r_i_i_C(3) * t334) * t301 + (t322 * t334 + (t341 + (t300 * r_i_i_C(2) + t344) * qJD(4)) * t298) * t297, (-r_i_i_C(2) * t319 - t323 * t344) * t302 + (t309 * r_i_i_C(1) + t310 * r_i_i_C(2)) * t298, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:09
	% EndTime: 2019-10-10 17:10:10
	% DurationCPUTime: 1.56s
	% Computational Cost: add. (553->140), mult. (1736->256), div. (0->0), fcn. (1741->10), ass. (0->97)
	t460 = sin(qJ(4));
	t465 = cos(qJ(4));
	t463 = sin(qJ(1));
	t466 = cos(qJ(3));
	t467 = cos(qJ(2));
	t522 = t466 * t467;
	t461 = sin(qJ(3));
	t468 = cos(qJ(1));
	t527 = t461 * t468;
	t448 = t463 * t522 + t527;
	t462 = sin(qJ(2));
	t516 = qJD(2) * t467;
	t501 = t463 * t516;
	t518 = qJD(1) * t468;
	t479 = t462 * t518 + t501;
	t475 = -qJD(4) * t448 + t479;
	t519 = qJD(1) * t467;
	t495 = qJD(3) + t519;
	t512 = qJD(3) * t467;
	t497 = t461 * t512;
	t517 = qJD(2) * t462;
	t502 = t463 * t517;
	t520 = qJD(1) * t463;
	t503 = t461 * t520;
	t521 = t468 * t466;
	t437 = -t463 * t497 - t466 * t502 + t495 * t521 - t503;
	t511 = qJD(4) * t462;
	t487 = t463 * t511 + t437;
	t423 = t475 * t460 + t487 * t465;
	t480 = t463 * t466 + t467 * t527;
	t436 = t480 * qJD(1) + t448 * qJD(3) - t461 * t502;
	t459 = sin(qJ(5));
	t464 = cos(qJ(5));
	t551 = t423 * t459 + t436 * t464;
	t550 = -t423 * t464 + t436 * t459;
	t494 = qJD(4) * t466 - qJD(2);
	t514 = qJD(3) * t461;
	t526 = t462 * t465;
	t493 = qJD(2) * t466 - qJD(4);
	t541 = t493 * t467;
	t549 = (-t462 * t514 + t541) * t460 + t494 * t526;
	t489 = t464 * r_i_i_C(1) - t459 * r_i_i_C(2);
	t536 = r_i_i_C(3) + pkin(3);
	t471 = (t489 * t460 - t536 * t465) * qJD(4);
	t488 = t459 * r_i_i_C(1) + t464 * r_i_i_C(2);
	t507 = qJD(5) * t465;
	t548 = t488 * t507 + t471;
	t505 = t536 * t460;
	t547 = t489 * t465 + t505;
	t529 = t460 * t462;
	t439 = t448 * t465 + t463 * t529;
	t523 = t463 * t461;
	t447 = t467 * t523 - t521;
	t545 = t439 * t459 + t447 * t464;
	t544 = t439 * t464 - t447 * t459;
	t470 = (t494 * t460 + t465 * t514) * t462 - t465 * t541;
	t515 = qJD(2) * t468;
	t500 = t467 * t515;
	t478 = t462 * t520 - t500;
	t538 = t493 * t462 + t497;
	t528 = t461 * t465;
	t525 = t462 * t466;
	t524 = t462 * t468;
	t513 = qJD(3) * t466;
	t510 = qJD(5) * t459;
	t509 = qJD(5) * t461;
	t508 = qJD(5) * t464;
	t506 = pkin(2) * t519;
	t496 = qJD(1) + t512;
	t491 = t461 * t501;
	t490 = t461 * t500;
	t484 = t494 * t467;
	t486 = t460 * t484 + t538 * t465 + t467 * t509;
	t485 = t462 * t509 + t470;
	t451 = t467 * t521 - t523;
	t442 = t451 * t465 + t460 * t524;
	t446 = -t460 * t467 + t465 * t525;
	t481 = t460 * t525 + t465 * t467;
	t444 = t446 * t468;
	t477 = qJD(1) * t481;
	t476 = t462 * t515 + t495 * t463;
	t474 = qJD(5) * t446 + t461 * t516 + t462 * t513;
	t473 = -qJD(5) * (t465 * t522 + t529) + t461 * t517 - t466 * t512;
	t422 = -t487 * t460 + t475 * t465;
	t469 = (t459 * t513 + t461 * t508) * r_i_i_C(1) + (-t459 * t509 + t464 * t513) * r_i_i_C(2) + qJD(2) * pkin(2);
	t443 = t446 * t463;
	t441 = -t451 * t460 + t465 * t524;
	t438 = -t448 * t460 + t463 * t526;
	t435 = t476 * t466 + t496 * t527;
	t434 = t476 * t461 - t496 * t521;
	t429 = -qJD(1) * t444 + t470 * t463;
	t427 = t446 * t520 + t470 * t468;
	t421 = (t468 * t511 - t435) * t465 + (-qJD(4) * t451 - t478) * t460;
	t420 = t442 * qJD(4) - t435 * t460 + t478 * t465;
	t419 = t421 * t464 + t434 * t459 + (-t442 * t459 - t464 * t480) * qJD(5);
	t418 = -t421 * t459 + t434 * t464 + (-t442 * t464 + t459 * t480) * qJD(5);
	t1 = [t550 * r_i_i_C(1) + t551 * r_i_i_C(2) + t536 * t422 + (t545 * r_i_i_C(1) + t544 * r_i_i_C(2)) * qJD(5) + t479 * pkin(2), (t427 * t464 + t444 * t510 + t459 * t490) * r_i_i_C(1) + (-t427 * t459 + t444 * t508 + t464 * t490) * r_i_i_C(2) + t463 * t506 + t536 * (t463 * t477 - t468 * t549) + (t469 * t468 - t488 * t503) * t462, (t435 * t459 - t451 * t508) * r_i_i_C(1) + (t435 * t464 + t451 * t510) * r_i_i_C(2) + t547 * t434 + t548 * t480, t536 * t421 + (t420 * t459 - t441 * t508) * r_i_i_C(2) + (-t420 * t464 - t441 * t510) * r_i_i_C(1), r_i_i_C(1) * t418 - r_i_i_C(2) * t419, 0, 0; t478 * pkin(2) + t419 * r_i_i_C(1) + t418 * r_i_i_C(2) + t536 * t420, (t429 * t464 + t443 * t510 + t459 * t491) * r_i_i_C(1) + (-t429 * t459 + t443 * t508 + t464 * t491) * r_i_i_C(2) - t468 * t506 - t536 * (t549 * t463 + t468 * t477) + (t488 * t461 * t518 + t469 * t463) * t462, (-t437 * t459 - t448 * t508) * r_i_i_C(1) + (-t437 * t464 + t448 * t510) * r_i_i_C(2) - t547 * t436 + t548 * t447, t536 * t423 + (-t422 * t459 - t438 * t508) * r_i_i_C(2) + (t422 * t464 - t438 * t510) * r_i_i_C(1), -t551 * r_i_i_C(1) + t550 * r_i_i_C(2) + (-t544 * r_i_i_C(1) + t545 * r_i_i_C(2)) * qJD(5), 0, 0; 0, -pkin(2) * t516 - t536 * (t538 * t460 - t465 * t484) + (-t486 * r_i_i_C(1) + t473 * r_i_i_C(2)) * t464 + (t473 * r_i_i_C(1) + t486 * r_i_i_C(2)) * t459, ((-t459 * t466 - t464 * t528) * r_i_i_C(1) + (t459 * t528 - t464 * t466) * r_i_i_C(2) - t461 * t505) * t516 + ((-qJD(3) * t547 - t489 * qJD(5)) * t466 + (t471 + t488 * (qJD(3) + t507)) * t461) * t462, -t536 * t470 + (t459 * t549 + t481 * t508) * r_i_i_C(2) + (-t464 * t549 + t481 * t510) * r_i_i_C(1), (-t474 * r_i_i_C(1) + t485 * r_i_i_C(2)) * t464 + (t485 * r_i_i_C(1) + t474 * r_i_i_C(2)) * t459, 0, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:13
	% EndTime: 2019-10-10 17:10:18
	% DurationCPUTime: 4.60s
	% Computational Cost: add. (1171->258), mult. (3601->450), div. (0->0), fcn. (3789->12), ass. (0->148)
	t692 = sin(qJ(3));
	t694 = sin(qJ(1));
	t698 = cos(qJ(3));
	t699 = cos(qJ(2));
	t765 = qJD(1) * t699;
	t727 = qJD(3) + t765;
	t758 = qJD(3) * t699;
	t734 = t692 * t758;
	t693 = sin(qJ(2));
	t763 = qJD(2) * t693;
	t740 = t694 * t763;
	t766 = qJD(1) * t694;
	t700 = cos(qJ(1));
	t767 = t700 * t698;
	t643 = -t692 * t766 - t694 * t734 - t698 * t740 + t727 * t767;
	t768 = t698 * t699;
	t777 = t692 * t700;
	t665 = t694 * t768 + t777;
	t691 = sin(qJ(4));
	t697 = cos(qJ(4));
	t762 = qJD(2) * t699;
	t764 = qJD(1) * t700;
	t713 = t693 * t764 + t694 * t762;
	t756 = qJD(4) * t693;
	t622 = (-qJD(4) * t665 + t713) * t691 + (t694 * t756 + t643) * t697;
	t715 = t694 * t698 + t699 * t777;
	t642 = t715 * qJD(1) + t665 * qJD(3) - t692 * t740;
	t690 = sin(qJ(5));
	t696 = cos(qJ(5));
	t781 = t691 * t693;
	t651 = t665 * t697 + t694 * t781;
	t773 = t694 * t692;
	t664 = t699 * t773 - t767;
	t793 = t651 * t690 + t664 * t696;
	t611 = t793 * qJD(5) - t622 * t696 + t642 * t690;
	t621 = t651 * qJD(4) + t643 * t691 - t713 * t697;
	t689 = sin(qJ(6));
	t695 = cos(qJ(6));
	t805 = t611 * t689 + t621 * t695;
	t804 = t611 * t695 - t621 * t689;
	t631 = t651 * t696 - t664 * t690;
	t776 = t693 * t697;
	t650 = t665 * t691 - t694 * t776;
	t803 = t631 * t689 - t650 * t695;
	t802 = t631 * t695 + t650 * t689;
	t609 = -t631 * qJD(5) - t622 * t690 - t642 * t696;
	t726 = qJD(4) * t698 - qJD(2);
	t718 = t726 * t691;
	t760 = qJD(3) * t697;
	t799 = (t692 * t760 + t718) * t693;
	t761 = qJD(2) * t700;
	t712 = t693 * t766 - t699 * t761;
	t779 = t692 * t693;
	t736 = qJD(3) * t779;
	t757 = qJD(4) * t691;
	t792 = -t691 * t736 - t697 * t763 - t699 * t757;
	t749 = qJD(6) * t695;
	t750 = qJD(6) * t689;
	t755 = qJD(4) * t697;
	t702 = -(t689 * t755 + t691 * t749) * r_i_i_C(1) - (-t691 * t750 + t695 * t755) * r_i_i_C(2) - pkin(3) * t755;
	t751 = qJD(5) * t697;
	t711 = t690 * t757 - t696 * t751;
	t791 = t711 * r_i_i_C(3) + t702;
	t790 = r_i_i_C(3) * t690;
	t789 = t691 * pkin(3);
	t784 = t689 * t691;
	t783 = t689 * t696;
	t782 = t690 * t697;
	t780 = t691 * t695;
	t778 = t692 * t696;
	t775 = t693 * t698;
	t774 = t693 * t700;
	t772 = t695 * t696;
	t771 = t696 * t697;
	t770 = t696 * t698;
	t769 = t697 * t699;
	t759 = qJD(3) * t698;
	t754 = qJD(5) * t690;
	t753 = qJD(5) * t692;
	t752 = qJD(5) * t696;
	t748 = qJD(6) * t696;
	t746 = t690 * t779;
	t745 = t697 * t768;
	t742 = t699 * t764;
	t741 = t692 * t762;
	t733 = t693 * t759;
	t732 = t693 * t755;
	t731 = t693 * t753;
	t729 = r_i_i_C(3) * qJD(1) * t778;
	t728 = qJD(1) + t758;
	t725 = qJD(2) * t698 - qJD(4);
	t724 = t696 * t741;
	t720 = (-t725 * t693 - t734) * t697 + (-t753 - t718) * t699;
	t639 = t725 * t769 - t799;
	t719 = -t639 + t731;
	t669 = t699 * t767 - t773;
	t656 = t669 * t697 + t691 * t774;
	t634 = t656 * t696 - t690 * t715;
	t633 = -t656 * t690 - t696 * t715;
	t717 = t725 * t699;
	t716 = -t690 * t698 - t692 * t771;
	t663 = -t691 * t699 + t697 * t775;
	t662 = t691 * t775 + t769;
	t714 = (-t689 * r_i_i_C(1) - t695 * r_i_i_C(2) - pkin(3)) * t691;
	t661 = t663 * t700;
	t708 = t693 * t761 + t727 * t694;
	t706 = -qJD(5) * t663 - t733 - t741;
	t667 = t745 + t781;
	t705 = -qJD(5) * t667 + t692 * t763 - t698 * t758;
	t704 = (t690 * t753 - t696 * t759) * r_i_i_C(3) + qJD(2) * pkin(2);
	t703 = -t697 * t717 + t799;
	t701 = (t689 * t748 + t695 * t754) * r_i_i_C(1) + (-t689 * t754 + t695 * t748) * r_i_i_C(2) - r_i_i_C(3) * t752;
	t666 = t691 * t768 - t776;
	t660 = t662 * t700;
	t659 = t663 * t694;
	t658 = t662 * t694;
	t657 = t716 * t693;
	t655 = t669 * t691 - t697 * t774;
	t654 = -t690 * t692 * t699 + t667 * t696;
	t649 = t663 * t696 - t746;
	t648 = -t663 * t690 - t693 * t778;
	t647 = -t661 * t696 + t700 * t746;
	t646 = -t659 * t696 + t694 * t746;
	t645 = -t669 * t690 - t715 * t771;
	t644 = -t664 * t771 - t665 * t690;
	t641 = t708 * t698 + t728 * t777;
	t640 = t708 * t692 - t728 * t767;
	t638 = (t691 * t762 + t732) * t698 + t792;
	t636 = t662 * qJD(2) - qJD(4) * t745 + (t734 - t756) * t691;
	t629 = -qJD(1) * t661 + t703 * t694;
	t628 = t713 * t691 * t698 + t697 * t742 + (t732 * t698 + t792) * t694;
	t627 = t663 * t766 + t703 * t700;
	t626 = t662 * t766 + (-t726 * t776 + (-t717 + t736) * t691) * t700;
	t625 = t716 * t762 + ((-qJD(5) - t760) * t770 + (t696 * t757 + (qJD(3) + t751) * t690) * t692) * t693;
	t620 = (t700 * t756 - t641) * t697 + (-qJD(4) * t669 - t712) * t691;
	t619 = t656 * qJD(4) - t641 * t691 + t697 * t712;
	t618 = t706 * t690 - t719 * t696;
	t617 = t719 * t690 + t706 * t696;
	t616 = t705 * t690 + t720 * t696;
	t615 = (t694 * t731 + t629) * t696 + (qJD(5) * t659 + t713 * t692 + t694 * t733) * t690;
	t614 = (t700 * t731 + t627) * t696 + (qJD(5) * t661 - t712 * t692 + t700 * t733) * t690;
	t613 = (t664 * t751 - t643) * t690 + (-qJD(5) * t665 - t642 * t697 + t664 * t757) * t696;
	t612 = (t715 * t751 + t641) * t690 + (-qJD(5) * t669 + t640 * t697 + t715 * t757) * t696;
	t608 = t633 * qJD(5) + t620 * t696 + t640 * t690;
	t607 = t634 * qJD(5) + t620 * t690 - t640 * t696;
	t606 = t608 * t695 + t619 * t689 + (-t634 * t689 + t655 * t695) * qJD(6);
	t605 = -t608 * t689 + t619 * t695 + (-t634 * t695 - t655 * t689) * qJD(6);
	t1 = [t804 * r_i_i_C(1) - t805 * r_i_i_C(2) + t609 * r_i_i_C(3) - t621 * pkin(3) + (t803 * r_i_i_C(1) + t802 * r_i_i_C(2)) * qJD(6) + t713 * pkin(2), (t614 * t695 + t626 * t689) * r_i_i_C(1) + (-t614 * t689 + t626 * t695) * r_i_i_C(2) + (t627 * t690 - t661 * t752 - t700 * t724) * r_i_i_C(3) + t626 * pkin(3) + t694 * pkin(2) * t765 + ((-t647 * t689 - t660 * t695) * r_i_i_C(1) + (-t647 * t695 + t660 * t689) * r_i_i_C(2)) * qJD(6) + (t694 * t729 + t704 * t700) * t693, (t612 * t695 + t640 * t784 - t645 * t750) * r_i_i_C(1) + (-t612 * t689 + t640 * t780 - t645 * t749) * r_i_i_C(2) + (t640 * t782 - t641 * t696 - t669 * t754) * r_i_i_C(3) + t640 * t789 + t791 * t715, (-t619 * t772 + t620 * t689 + t656 * t749) * r_i_i_C(1) + (t619 * t783 + t620 * t695 - t656 * t750) * r_i_i_C(2) - t619 * t790 + t620 * pkin(3) + t701 * t655, r_i_i_C(3) * t608 + (t607 * t689 - t633 * t749) * r_i_i_C(2) + (-t607 * t695 - t633 * t750) * r_i_i_C(1), r_i_i_C(1) * t605 - r_i_i_C(2) * t606, 0; t712 * pkin(2) + t619 * pkin(3) + t606 * r_i_i_C(1) + t605 * r_i_i_C(2) + t607 * r_i_i_C(3), (t615 * t695 - t628 * t689) * r_i_i_C(1) + (-t615 * t689 - t628 * t695) * r_i_i_C(2) + (t629 * t690 - t659 * t752 - t694 * t724) * r_i_i_C(3) - t628 * pkin(3) - pkin(2) * t742 + ((-t646 * t689 - t658 * t695) * r_i_i_C(1) + (-t646 * t695 + t658 * t689) * r_i_i_C(2)) * qJD(6) + (t704 * t694 - t700 * t729) * t693, (t613 * t695 - t642 * t784 - t644 * t750) * r_i_i_C(1) + (-t613 * t689 - t642 * t780 - t644 * t749) * r_i_i_C(2) + (-t642 * t782 + t643 * t696 - t665 * t754) * r_i_i_C(3) - t642 * t789 + t791 * t664, (-t621 * t772 + t622 * t689 + t651 * t749) * r_i_i_C(1) + (t621 * t783 + t622 * t695 - t651 * t750) * r_i_i_C(2) - t621 * t790 + t622 * pkin(3) + t701 * t650, -r_i_i_C(3) * t611 + (-t609 * t689 + t749 * t793) * r_i_i_C(2) + (t609 * t695 + t750 * t793) * r_i_i_C(1), t805 * r_i_i_C(1) + t804 * r_i_i_C(2) + (-t802 * r_i_i_C(1) + t803 * r_i_i_C(2)) * qJD(6), 0; 0, (t616 * t695 - t636 * t689) * r_i_i_C(1) + (-t616 * t689 - t636 * t695) * r_i_i_C(2) - t636 * pkin(3) - pkin(2) * t762 + t720 * t790 - t705 * r_i_i_C(3) * t696 + ((-t654 * t689 + t666 * t695) * r_i_i_C(1) + (-t654 * t695 - t666 * t689) * r_i_i_C(2)) * qJD(6), (t625 * t695 - t657 * t750) * r_i_i_C(1) + (-t625 * t689 - t657 * t749) * r_i_i_C(2) + ((-t692 * t782 + t770) * r_i_i_C(3) + t692 * t714) * t762 + ((-r_i_i_C(3) * t754 + (-r_i_i_C(3) * t782 + t714) * qJD(3)) * t698 + ((-qJD(3) * t696 + t711) * r_i_i_C(3) + t702) * t692) * t693, (-t638 * t772 + t639 * t689 + t663 * t749) * r_i_i_C(1) + (t638 * t783 + t639 * t695 - t663 * t750) * r_i_i_C(2) - t638 * t790 + t639 * pkin(3) + t701 * t662, r_i_i_C(3) * t618 + (-t617 * t689 - t648 * t749) * r_i_i_C(2) + (t617 * t695 - t648 * t750) * r_i_i_C(1), (-t618 * t689 + t638 * t695) * r_i_i_C(1) + (-t618 * t695 - t638 * t689) * r_i_i_C(2) + ((-t649 * t695 - t662 * t689) * r_i_i_C(1) + (t649 * t689 - t662 * t695) * r_i_i_C(2)) * qJD(6), 0;];
	JaD_transl = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:26
	% EndTime: 2019-10-10 17:10:35
	% DurationCPUTime: 8.95s
	% Computational Cost: add. (2918->345), mult. (8814->601), div. (0->0), fcn. (9592->14), ass. (0->214)
	t978 = cos(qJ(2));
	t1072 = qJD(2) * t978;
	t979 = cos(qJ(1));
	t1074 = qJD(1) * t979;
	t971 = sin(qJ(2));
	t972 = sin(qJ(1));
	t1003 = t972 * t1072 + t971 * t1074;
	t1066 = qJD(4) * t971;
	t1075 = qJD(1) * t978;
	t1033 = qJD(3) + t1075;
	t1068 = qJD(3) * t978;
	t970 = sin(qJ(3));
	t1039 = t970 * t1068;
	t1073 = qJD(2) * t971;
	t1046 = t972 * t1073;
	t1076 = qJD(1) * t972;
	t977 = cos(qJ(3));
	t1077 = t979 * t977;
	t908 = t1033 * t1077 - t972 * t1039 - t977 * t1046 - t970 * t1076;
	t1078 = t977 * t978;
	t1086 = t970 * t979;
	t941 = t972 * t1078 + t1086;
	t969 = sin(qJ(4));
	t976 = cos(qJ(4));
	t868 = (t972 * t1066 + t908) * t976 + (-qJD(4) * t941 + t1003) * t969;
	t1090 = t969 * t971;
	t923 = t972 * t1090 + t941 * t976;
	t1082 = t972 * t970;
	t940 = t978 * t1082 - t1077;
	t968 = sin(qJ(5));
	t975 = cos(qJ(5));
	t893 = t923 * t975 - t940 * t968;
	t1006 = t978 * t1086 + t972 * t977;
	t907 = t1006 * qJD(1) + t941 * qJD(3) - t970 * t1046;
	t845 = t893 * qJD(5) + t868 * t968 + t907 * t975;
	t966 = sin(qJ(7));
	t1118 = t845 * t966;
	t973 = cos(qJ(7));
	t1117 = t845 * t973;
	t1085 = t971 * t976;
	t922 = -t972 * t1085 + t941 * t969;
	t967 = sin(qJ(6));
	t974 = cos(qJ(6));
	t872 = t893 * t974 + t922 * t967;
	t892 = t923 * t968 + t940 * t975;
	t1116 = t872 * t966 + t892 * t973;
	t1115 = t872 * t973 - t892 * t966;
	t848 = t892 * qJD(5) - t868 * t975 + t907 * t968;
	t867 = t923 * qJD(4) - t1003 * t976 + t908 * t969;
	t1114 = t872 * qJD(6) - t848 * t967 - t867 * t974;
	t1111 = t893 * t967 - t922 * t974;
	t1110 = t867 * t967;
	t1032 = qJD(4) * t977 - qJD(2);
	t1014 = t1032 * t969;
	t1069 = qJD(3) * t976;
	t1108 = (t970 * t1069 + t1014) * t971;
	t1031 = qJD(2) * t977 - qJD(4);
	t1013 = t1031 * t978;
	t1063 = qJD(5) * t970;
	t1036 = t971 * t1063;
	t1103 = -t976 * t1013 + t1036 + t1108;
	t1079 = t976 * t978;
	t904 = t1031 * t1079 - t1108;
	t1101 = t904 - t1036;
	t1071 = qJD(2) * t979;
	t1002 = -t978 * t1071 + t971 * t1076;
	t1070 = qJD(3) * t971;
	t1041 = t970 * t1070;
	t1067 = qJD(4) * t969;
	t1100 = -t969 * t1041 - t978 * t1067 - t976 * t1073;
	t1038 = t977 * t1070;
	t1084 = t971 * t977;
	t939 = t976 * t1084 - t969 * t978;
	t935 = t939 * t979;
	t1099 = -qJD(5) * t935 + t1002 * t970 - t979 * t1038;
	t933 = t939 * t972;
	t1098 = qJD(5) * t933 + t1003 * t970 + t972 * t1038;
	t1097 = r_i_i_C(3) + pkin(4);
	t1092 = t967 * t969;
	t1091 = t968 * t976;
	t1089 = t970 * t971;
	t1088 = t970 * t975;
	t1087 = t970 * t978;
	t1083 = t971 * t979;
	t1081 = t974 * t975;
	t1080 = t975 * t976;
	t1065 = qJD(4) * t976;
	t1064 = qJD(5) * t968;
	t1062 = qJD(5) * t975;
	t1061 = qJD(5) * t976;
	t1060 = qJD(6) * t967;
	t1059 = qJD(6) * t969;
	t1058 = qJD(6) * t975;
	t1057 = qJD(7) * t966;
	t1056 = qJD(7) * t968;
	t1055 = qJD(7) * t973;
	t1054 = t968 * t1089;
	t1053 = t969 * t1089;
	t1052 = t971 * t1088;
	t1051 = t970 * t1083;
	t1050 = t976 * t1078;
	t1047 = t978 * t1074;
	t1044 = t971 * t1071;
	t1037 = t971 * t1065;
	t1034 = qJD(1) + t1068;
	t1030 = qJD(3) + t1061;
	t1028 = r_i_i_C(1) * t973 - r_i_i_C(2) * t966;
	t1025 = (-t1031 * t971 - t1039) * t976 + (-t1063 - t1014) * t978;
	t988 = t1033 * t972 + t1044;
	t906 = t1034 * t1086 + t988 * t977;
	t1024 = -t1006 * t1061 - t906;
	t1023 = t940 * t1061 - t908;
	t905 = -t1034 * t1077 + t988 * t970;
	t945 = t978 * t1077 - t1082;
	t998 = -qJD(5) * t945 + t1006 * t1067 + t905 * t976;
	t1022 = -t1006 * t1059 - t1024 * t968 + t998 * t975;
	t997 = -qJD(5) * t941 + t940 * t1067 - t907 * t976;
	t1021 = -t1023 * t968 + t940 * t1059 - t997 * t975;
	t866 = (t979 * t1066 - t906) * t976 + (-qJD(4) * t945 - t1002) * t969;
	t928 = -t976 * t1083 + t945 * t969;
	t1020 = t928 * t1058 + t866;
	t1019 = t922 * t1058 + t868;
	t938 = t969 * t1084 + t1079;
	t1018 = t938 * t1058 + t904;
	t929 = t969 * t1083 + t945 * t976;
	t865 = t929 * qJD(4) + t1002 * t976 - t906 * t969;
	t994 = qJD(6) * t929 + t928 * t1064 - t865 * t975;
	t1017 = t1020 * t967 + t928 * t1056 + t994 * t974;
	t993 = qJD(6) * t923 + t922 * t1064 - t867 * t975;
	t1016 = t1019 * t967 + t922 * t1056 + t993 * t974;
	t903 = (t969 * t1072 + t1037) * t977 + t1100;
	t992 = qJD(6) * t939 + t938 * t1064 - t903 * t975;
	t1015 = t1018 * t967 + t938 * t1056 + t992 * t974;
	t898 = -t1006 * t968 + t929 * t975;
	t875 = t898 * t974 + t928 * t967;
	t916 = t972 * t1054 - t933 * t975;
	t932 = t938 * t972;
	t886 = t916 * t974 - t932 * t967;
	t918 = t968 * t1051 - t935 * t975;
	t934 = t938 * t979;
	t887 = t918 * t974 - t934 * t967;
	t921 = t939 * t975 - t1054;
	t891 = t921 * t974 + t938 * t967;
	t943 = t1050 + t1090;
	t927 = -t968 * t1087 + t943 * t975;
	t942 = t969 * t1078 - t1085;
	t896 = t927 * t974 + t942 * t967;
	t897 = t1006 * t975 + t929 * t968;
	t1012 = (-qJD(5) - t1069) * t977;
	t1011 = t939 * t1076 + t1103 * t979;
	t1010 = -qJD(1) * t935 + t1103 * t972;
	t1008 = -t970 * t1080 - t968 * t977;
	t1009 = qJD(6) * t1053 - t1008 * t1072 - (t975 * t1012 + (t1030 * t968 + t975 * t1067) * t970) * t971;
	t1007 = -t970 * t1091 + t975 * t977;
	t1005 = -t1006 * t1065 + t905 * t969;
	t1004 = -t940 * t1065 - t907 * t969;
	t1001 = -t970 * t1072 - t1038;
	t913 = -t1006 * t1080 - t945 * t968;
	t996 = -qJD(6) * t913 + t1005;
	t911 = -t940 * t1080 - t941 * t968;
	t995 = -qJD(6) * t911 + t1004;
	t991 = -qJD(7) * (-t928 * t1081 + t929 * t967) + t928 * t1062 + t865 * t968;
	t990 = -qJD(7) * (-t922 * t1081 + t923 * t967) + t922 * t1062 + t867 * t968;
	t989 = -qJD(7) * (-t938 * t1081 + t939 * t967) + t938 * t1062 + t903 * t968;
	t986 = qJD(5) * t939 - t1001;
	t985 = -qJD(5) * t943 - t977 * t1068 + t970 * t1073;
	t984 = -t1028 * t974 + t1097 * t967;
	t983 = t1001 * t969 - t970 * t1037;
	t931 = t1008 * t971;
	t982 = -qJD(6) * t931 + t983;
	t980 = (r_i_i_C(1) * t966 + r_i_i_C(2) * t973) * t974 * qJD(7) + (t1028 * t967 + t1097 * t974) * qJD(6);
	t930 = t1007 * t971;
	t926 = t975 * t1087 + t943 * t968;
	t920 = t939 * t968 + t1052;
	t917 = -t975 * t1051 - t935 * t968;
	t915 = -t972 * t1052 - t933 * t968;
	t914 = -t967 * t1053 + t931 * t974;
	t912 = -t1006 * t1091 + t945 * t975;
	t910 = -t940 * t1091 + t941 * t975;
	t901 = t938 * qJD(2) - qJD(4) * t1050 + (t1039 - t1066) * t969;
	t890 = -t921 * t967 + t938 * t974;
	t885 = -t1006 * t1092 + t913 * t974;
	t884 = -t940 * t1092 + t911 * t974;
	t880 = t1003 * t969 * t977 + t976 * t1047 + (t1037 * t977 + t1100) * t972;
	t878 = t938 * t1076 + (-t1032 * t1085 + (-t1013 + t1041) * t969) * t979;
	t876 = t1007 * t1072 + (-t1030 * t1088 + (t970 * t1067 + t1012) * t968) * t971;
	t874 = -t898 * t967 + t928 * t974;
	t864 = t1101 * t975 - t986 * t968;
	t863 = t1101 * t968 + t975 * t986;
	t862 = t1025 * t975 + t985 * t968;
	t861 = t1025 * t968 - t985 * t975;
	t860 = t1010 * t975 + t1098 * t968;
	t859 = t1010 * t968 - t1098 * t975;
	t858 = t1011 * t975 - t1099 * t968;
	t857 = t1011 * t968 + t1099 * t975;
	t856 = -t1009 * t974 + t982 * t967;
	t853 = -t1023 * t975 + t997 * t968;
	t851 = t1024 * t975 + t998 * t968;
	t844 = -t897 * qJD(5) + t866 * t975 + t905 * t968;
	t843 = t898 * qJD(5) + t866 * t968 - t905 * t975;
	t842 = t903 * t967 - t921 * t1060 + (qJD(6) * t938 + t864) * t974;
	t841 = -t891 * qJD(6) - t864 * t967 + t903 * t974;
	t840 = t862 * t974 - t901 * t967 + (-t927 * t967 + t942 * t974) * qJD(6);
	t838 = -t1021 * t974 + t995 * t967;
	t836 = t1022 * t974 + t996 * t967;
	t830 = t860 * t974 - t880 * t967 + (-t916 * t967 - t932 * t974) * qJD(6);
	t828 = t858 * t974 + t878 * t967 + (-t918 * t967 - t934 * t974) * qJD(6);
	t826 = t1111 * qJD(6) + t848 * t974 - t1110;
	t824 = t1110 - t893 * t1060 + (qJD(6) * t922 - t848) * t974;
	t822 = t865 * t967 - t898 * t1060 + (qJD(6) * t928 + t844) * t974;
	t821 = -t875 * qJD(6) - t844 * t967 + t865 * t974;
	t820 = t822 * t973 - t843 * t966 + (-t875 * t966 - t897 * t973) * qJD(7);
	t819 = -t822 * t966 - t843 * t973 + (-t875 * t973 + t897 * t966) * qJD(7);
	t1 = [(t826 * t973 + t1118) * r_i_i_C(1) + (-t826 * t966 + t1117) * r_i_i_C(2) - t867 * pkin(3) + t1097 * t1114 + (t1116 * r_i_i_C(1) + t1115 * r_i_i_C(2)) * qJD(7) + t1003 * pkin(2), (t828 * t973 - t857 * t966) * r_i_i_C(1) + (-t828 * t966 - t857 * t973) * r_i_i_C(2) + t878 * pkin(3) + t1097 * (-t887 * qJD(6) - t858 * t967 + t878 * t974) + ((-t887 * t966 - t917 * t973) * r_i_i_C(1) + (-t887 * t973 + t917 * t966) * r_i_i_C(2)) * qJD(7) + (t972 * t1075 + t1044) * pkin(2), (t836 * t973 - t851 * t966) * r_i_i_C(1) + (-t836 * t966 - t851 * t973) * r_i_i_C(2) + t1097 * (-t1022 * t967 + t996 * t974) + ((-t885 * t966 - t912 * t973) * r_i_i_C(1) + (-t885 * t973 + t912 * t966) * r_i_i_C(2)) * qJD(7) + t1005 * pkin(3), t866 * pkin(3) + t1097 * (t1020 * t974 - t994 * t967) + (t1017 * r_i_i_C(1) + t991 * r_i_i_C(2)) * t973 + (t991 * r_i_i_C(1) - t1017 * r_i_i_C(2)) * t966, (-t1055 * t898 - t844 * t966) * r_i_i_C(1) + (t1057 * t898 - t844 * t973) * r_i_i_C(2) + t984 * t843 + t980 * t897, -t1097 * t822 + (-t1055 * t874 - t821 * t966) * r_i_i_C(2) + (-t1057 * t874 + t821 * t973) * r_i_i_C(1), r_i_i_C(1) * t819 - r_i_i_C(2) * t820; t1002 * pkin(2) + t865 * pkin(3) + t820 * r_i_i_C(1) + t819 * r_i_i_C(2) + t1097 * t821, (t830 * t973 - t859 * t966) * r_i_i_C(1) + (-t830 * t966 - t859 * t973) * r_i_i_C(2) - t880 * pkin(3) + t1097 * (-t886 * qJD(6) - t860 * t967 - t880 * t974) + ((-t886 * t966 - t915 * t973) * r_i_i_C(1) + (-t886 * t973 + t915 * t966) * r_i_i_C(2)) * qJD(7) + (t1046 - t1047) * pkin(2), (t838 * t973 - t853 * t966) * r_i_i_C(1) + (-t838 * t966 - t853 * t973) * r_i_i_C(2) + t1097 * (t1021 * t967 + t995 * t974) + ((-t884 * t966 - t910 * t973) * r_i_i_C(1) + (-t884 * t973 + t910 * t966) * r_i_i_C(2)) * qJD(7) + t1004 * pkin(3), t868 * pkin(3) + t1097 * (t1019 * t974 - t993 * t967) + (t1016 * r_i_i_C(1) + t990 * r_i_i_C(2)) * t973 + (t990 * r_i_i_C(1) - t1016 * r_i_i_C(2)) * t966, (-t1055 * t893 + t848 * t966) * r_i_i_C(1) + (t1057 * t893 + t848 * t973) * r_i_i_C(2) + t984 * t845 + t980 * t892, -t1097 * t824 + (t1055 * t1111 + t1114 * t966) * r_i_i_C(2) + (t1057 * t1111 - t1114 * t973) * r_i_i_C(1), (-t824 * t966 - t1117) * r_i_i_C(1) + (-t824 * t973 + t1118) * r_i_i_C(2) + (-t1115 * r_i_i_C(1) + t1116 * r_i_i_C(2)) * qJD(7); 0, (t840 * t973 - t861 * t966) * r_i_i_C(1) + (-t840 * t966 - t861 * t973) * r_i_i_C(2) - t901 * pkin(3) - pkin(2) * t1072 + t1097 * (-t896 * qJD(6) - t862 * t967 - t901 * t974) + ((-t896 * t966 - t926 * t973) * r_i_i_C(1) + (-t896 * t973 + t926 * t966) * r_i_i_C(2)) * qJD(7), (t856 * t973 - t876 * t966) * r_i_i_C(1) + (-t856 * t966 - t876 * t973) * r_i_i_C(2) + t1097 * (t1009 * t967 + t982 * t974) + ((-t914 * t966 - t930 * t973) * r_i_i_C(1) + (-t914 * t973 + t930 * t966) * r_i_i_C(2)) * qJD(7) + t983 * pkin(3), t904 * pkin(3) + t1097 * (t1018 * t974 - t992 * t967) + (t1015 * r_i_i_C(1) + t989 * r_i_i_C(2)) * t973 + (t989 * r_i_i_C(1) - t1015 * r_i_i_C(2)) * t966, (-t1055 * t921 - t864 * t966) * r_i_i_C(1) + (t1057 * t921 - t864 * t973) * r_i_i_C(2) + t984 * t863 + t980 * t920, -t1097 * t842 + (-t1055 * t890 - t841 * t966) * r_i_i_C(2) + (-t1057 * t890 + t841 * t973) * r_i_i_C(1), (-t842 * t966 - t863 * t973) * r_i_i_C(1) + (-t842 * t973 + t863 * t966) * r_i_i_C(2) + ((-t891 * t973 + t920 * t966) * r_i_i_C(1) + (t891 * t966 + t920 * t973) * r_i_i_C(2)) * qJD(7);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,7);
end