% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14V3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14V3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14V3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobiaD_transl_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
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
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:23
	% DurationCPUTime: 0.10s
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
	t1 = [t20 * t23 + (-r_i_i_C(3) * t16 + t21 * t18) * qJD(1), (t15 * t22 + t17 * t25) * r_i_i_C(2) + (t15 * t25 - t17 * t22) * r_i_i_C(1), 0, 0, 0, 0; -t18 * t19 + (r_i_i_C(3) * t18 + t21 * t16) * qJD(1), (t15 * t23 - t17 * t24) * r_i_i_C(2) + (-t15 * t24 - t17 * t23) * r_i_i_C(1), 0, 0, 0, 0; 0, -t19, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:23
	% EndTime: 2019-10-10 11:13:24
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (31->15), mult. (100->32), div. (0->0), fcn. (71->4), ass. (0->14)
	t132 = sin(qJ(2));
	t134 = cos(qJ(2));
	t145 = r_i_i_C(3) + qJ(3);
	t139 = -r_i_i_C(1) * t132 + t145 * t134;
	t136 = t139 * qJD(2) + t132 * qJD(3);
	t148 = r_i_i_C(2) * qJD(1) + t136;
	t133 = sin(qJ(1));
	t143 = qJD(1) * t133;
	t135 = cos(qJ(1));
	t142 = qJD(1) * t135;
	t141 = qJD(2) * t135;
	t138 = -r_i_i_C(1) * t134 - t145 * t132;
	t137 = qJD(1) * t138;
	t1 = [-t148 * t133 + t135 * t137, (r_i_i_C(1) * t143 - t145 * t141) * t132 + ((-r_i_i_C(1) * qJD(2) + qJD(3)) * t135 - t145 * t143) * t134, -t132 * t143 + t134 * t141, 0, 0, 0; t133 * t137 + t148 * t135, t139 * t142 + (t138 * qJD(2) + qJD(3) * t134) * t133, t133 * qJD(2) * t134 + t132 * t142, 0, 0, 0; 0, t136, qJD(2) * t132, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:24
	% EndTime: 2019-10-10 11:13:24
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (77->35), mult. (254->70), div. (0->0), fcn. (208->6), ass. (0->31)
	t197 = sin(qJ(2));
	t196 = sin(qJ(4));
	t199 = cos(qJ(4));
	t207 = r_i_i_C(1) * t199 - r_i_i_C(2) * t196;
	t200 = cos(qJ(2));
	t222 = r_i_i_C(3) + qJ(3);
	t223 = t222 * t200;
	t203 = -t207 * t197 + t223;
	t201 = cos(qJ(1));
	t221 = t199 * t201;
	t198 = sin(qJ(1));
	t220 = qJD(1) * t198;
	t219 = qJD(1) * t201;
	t218 = qJD(2) * t197;
	t217 = qJD(2) * t198;
	t216 = qJD(2) * t200;
	t215 = qJD(2) * t201;
	t214 = qJD(4) * t197;
	t213 = qJD(4) * t200;
	t210 = qJD(1) * t222;
	t209 = -qJD(1) + t213;
	t208 = qJD(1) * t200 - qJD(4);
	t206 = r_i_i_C(1) * t196 + r_i_i_C(2) * t199;
	t205 = t209 * t196;
	t204 = t197 * t215 + t208 * t198;
	t202 = qJD(3) * t200 + t206 * t214 + (-t222 * t197 - t207 * t200) * qJD(2);
	t195 = -t208 * t221 + (t199 * t218 + t205) * t198;
	t194 = t209 * t199 * t198 + (-t197 * t217 + t208 * t201) * t196;
	t193 = t204 * t199 + t201 * t205;
	t192 = t204 * t196 - t209 * t221;
	t1 = [t195 * r_i_i_C(1) + t194 * r_i_i_C(2) - t223 * t217 + (-qJD(3) * t198 - t201 * t210) * t197, t202 * t201 - t203 * t220, -t197 * t220 + t200 * t215, t192 * r_i_i_C(1) + t193 * r_i_i_C(2), 0, 0; -t193 * r_i_i_C(1) + t192 * r_i_i_C(2) + t223 * t215 + (qJD(3) * t201 - t198 * t210) * t197, t202 * t198 + t203 * t219, t197 * t219 + t198 * t216, -t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0, 0; 0, t203 * qJD(2) + t197 * qJD(3) - t206 * t213, t218, (t196 * t214 - t199 * t216) * r_i_i_C(2) + (-t196 * t216 - t199 * t214) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:25
	% EndTime: 2019-10-10 11:13:26
	% DurationCPUTime: 0.48s
	% Computational Cost: add. (197->79), mult. (636->146), div. (0->0), fcn. (595->8), ass. (0->61)
	t311 = cos(qJ(4));
	t361 = -qJD(5) * t311 + qJD(2);
	t308 = sin(qJ(2));
	t309 = sin(qJ(1));
	t312 = cos(qJ(2));
	t332 = qJD(1) * t312 - qJD(4);
	t313 = cos(qJ(1));
	t344 = qJD(2) * t313;
	t360 = t308 * t344 + t332 * t309;
	t306 = sin(qJ(5));
	t310 = cos(qJ(5));
	t307 = sin(qJ(4));
	t343 = qJD(4) * t307;
	t318 = t306 * t343 + t361 * t310;
	t319 = -t361 * t306 + t310 * t343;
	t342 = qJD(4) * t311;
	t359 = t319 * r_i_i_C(1) - t318 * r_i_i_C(2) - r_i_i_C(3) * t342 - qJD(2) * qJ(3);
	t354 = t308 * t311;
	t355 = t307 * r_i_i_C(3);
	t358 = (-t306 * t312 + t310 * t354) * r_i_i_C(1) - (t306 * t354 + t310 * t312) * r_i_i_C(2) - t312 * qJ(3) + t308 * t355;
	t357 = t306 * r_i_i_C(1);
	t356 = t306 * r_i_i_C(2);
	t353 = t308 * t313;
	t352 = t309 * t311;
	t351 = t309 * t312;
	t350 = t313 * t307;
	t349 = t313 * t311;
	t348 = qJD(1) * t309;
	t347 = qJD(1) * t313;
	t346 = qJD(2) * t308;
	t345 = qJD(2) * t312;
	t341 = qJD(5) * t306;
	t340 = qJD(5) * t309;
	t339 = qJD(5) * t310;
	t337 = t307 * t351;
	t336 = t309 * t346;
	t335 = t309 * t345;
	t333 = -qJD(4) * t312 + qJD(1);
	t331 = qJD(2) * t311 - qJD(5);
	t330 = -t310 * r_i_i_C(1) + t356;
	t327 = t333 * t313;
	t295 = t307 * t327 - t360 * t311;
	t329 = qJD(5) * t353 + t295;
	t301 = t309 * t307 + t312 * t349;
	t297 = t301 * qJD(1) - qJD(4) * t337 - t311 * t336 - t313 * t342;
	t328 = -t308 * t340 - t297;
	t326 = t331 * t310;
	t321 = t308 * t347 + t335;
	t320 = -t308 * t348 + t312 * t344;
	t299 = t311 * t351 - t350;
	t317 = -qJD(5) * t299 + t321;
	t316 = -qJD(5) * t301 + t320;
	t315 = -r_i_i_C(1) * t326 - qJD(2) * t355 + t331 * t356 + qJD(3);
	t314 = t359 * t308 + t315 * t312;
	t300 = -t312 * t350 + t352;
	t298 = -t337 - t349;
	t296 = t333 * t352 + (-t332 * t313 + t336) * t307;
	t294 = t360 * t307 + t311 * t327;
	t293 = t316 * t306 + t329 * t310;
	t292 = -t329 * t306 + t316 * t310;
	t1 = [(-t297 * t310 + t299 * t341 - t306 * t335) * r_i_i_C(1) + (t297 * t306 + t299 * t339 - t310 * t335) * r_i_i_C(2) + t296 * r_i_i_C(3) - qJ(3) * t335 + ((-t306 * t347 - t309 * t339) * r_i_i_C(1) + (t306 * t340 - t310 * t347) * r_i_i_C(2) - qJ(3) * t347 - t309 * qJD(3)) * t308, t314 * t313 + t358 * t348, t320, t295 * r_i_i_C(3) + (-t294 * t306 - t300 * t339) * r_i_i_C(2) + (t294 * t310 - t300 * t341) * r_i_i_C(1), t292 * r_i_i_C(1) - t293 * r_i_i_C(2), 0; t293 * r_i_i_C(1) + t292 * r_i_i_C(2) - t294 * r_i_i_C(3) + t320 * qJ(3) + qJD(3) * t353, t314 * t309 - t358 * t347, t321, t297 * r_i_i_C(3) + (-t296 * t306 - t298 * t339) * r_i_i_C(2) + (t296 * t310 - t298 * t341) * r_i_i_C(1), (t317 * r_i_i_C(1) + t328 * r_i_i_C(2)) * t310 + (t328 * r_i_i_C(1) - t317 * r_i_i_C(2)) * t306, 0; 0, t315 * t308 - t359 * t312, t346, (t330 * t308 * qJD(4) + r_i_i_C(3) * t345) * t311 + (t330 * t345 + (-qJD(4) * r_i_i_C(3) + (t310 * r_i_i_C(2) + t357) * qJD(5)) * t308) * t307, (-r_i_i_C(2) * t326 - t331 * t357) * t312 + (t318 * r_i_i_C(1) + t319 * r_i_i_C(2)) * t308, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:13:27
	% EndTime: 2019-10-10 11:13:28
	% DurationCPUTime: 1.46s
	% Computational Cost: add. (479->149), mult. (1504->264), div. (0->0), fcn. (1518->10), ass. (0->95)
	t459 = sin(qJ(4));
	t460 = sin(qJ(2));
	t458 = sin(qJ(5));
	t464 = cos(qJ(4));
	t538 = -qJD(5) * t464 + qJD(2);
	t483 = t538 * t458;
	t463 = cos(qJ(5));
	t511 = qJD(4) * t463;
	t546 = (t459 * t511 - t483) * t460;
	t465 = cos(qJ(2));
	t466 = cos(qJ(1));
	t521 = t464 * t466;
	t461 = sin(qJ(1));
	t525 = t461 * t459;
	t443 = t465 * t521 + t525;
	t508 = qJD(4) * t466;
	t490 = t464 * t508;
	t512 = qJD(4) * t459;
	t492 = t461 * t512;
	t515 = qJD(2) * t460;
	t496 = t461 * t515;
	t429 = t443 * qJD(1) - t464 * t496 - t465 * t492 - t490;
	t519 = t466 * t459;
	t522 = t464 * t465;
	t440 = t461 * t522 - t519;
	t514 = qJD(2) * t465;
	t517 = qJD(1) * t466;
	t479 = t460 * t517 + t461 * t514;
	t506 = qJD(5) * t460;
	t419 = (-qJD(5) * t440 + t479) * t458 + (t461 * t506 + t429) * t463;
	t510 = qJD(4) * t464;
	t518 = qJD(1) * t461;
	t428 = t461 * t510 * t465 - t464 * t518 + (t517 * t465 - t496 - t508) * t459;
	t457 = sin(qJ(6));
	t462 = cos(qJ(6));
	t545 = t419 * t457 - t428 * t462;
	t544 = -t419 * t462 - t428 * t457;
	t523 = t463 * t465;
	t527 = t460 * t464;
	t480 = t458 * t527 + t523;
	t481 = (r_i_i_C(1) * t457 + r_i_i_C(2) * t462) * t459;
	t542 = -t480 * r_i_i_C(3) + t465 * qJ(3) - t460 * t481;
	t529 = t458 * t460;
	t431 = t440 * t463 + t461 * t529;
	t439 = t465 * t525 + t521;
	t541 = t431 * t457 - t439 * t462;
	t540 = t431 * t462 + t439 * t457;
	t513 = qJD(2) * t466;
	t478 = -t460 * t518 + t465 * t513;
	t502 = qJD(6) * t462;
	t503 = qJD(6) * t459;
	t536 = (t457 * t510 + t459 * t502) * r_i_i_C(1) + (-t457 * t503 + t462 * t510) * r_i_i_C(2) - (t458 * t512 + t538 * t463) * r_i_i_C(3) + qJD(2) * qJ(3);
	t535 = r_i_i_C(3) * t458;
	t530 = t457 * t463;
	t528 = t460 * t463;
	t526 = t460 * t466;
	t524 = t462 * t463;
	t509 = qJD(4) * t465;
	t507 = qJD(5) * t458;
	t504 = qJD(6) * t457;
	t501 = qJD(6) * t463;
	t500 = t460 * qJD(3);
	t499 = qJD(5) * t463 * r_i_i_C(3);
	t491 = t460 * t513;
	t488 = qJD(2) * t464 - qJD(5);
	t425 = t488 * t523 - t546;
	t484 = -t460 * t503 - t425;
	t482 = t488 * t465;
	t434 = t443 * t463 + t458 * t526;
	t438 = -t458 * t465 + t463 * t527;
	t475 = t457 * t501 + t462 * t507;
	t474 = -t457 * t507 + t462 * t501;
	t436 = t438 * t466;
	t471 = -qJD(6) * t438 + t459 * t514 + t460 * t510;
	t418 = -t431 * qJD(5) - t429 * t458 + t479 * t463;
	t470 = r_i_i_C(3) * t507 + qJD(3) + (-t464 * t535 - t481) * qJD(2);
	t469 = -t463 * t482 + t546;
	t468 = t475 * r_i_i_C(1) + t474 * r_i_i_C(2) - t499;
	t467 = -t536 * t460 + t470 * t465;
	t442 = -t461 * t464 + t465 * t519;
	t441 = t463 * t522 + t529;
	t435 = t438 * t461;
	t433 = -t443 * t458 + t463 * t526;
	t430 = -t440 * t458 + t461 * t528;
	t427 = (qJD(1) - t509) * t519 + (-t491 + (-qJD(1) * t465 + qJD(4)) * t461) * t464;
	t426 = t439 * qJD(1) + t459 * t491 - t465 * t490 - t492;
	t424 = t538 * t528 + (t460 * t512 - t482) * t458;
	t423 = t465 * t483 + (-t459 * t509 - t488 * t460) * t463;
	t422 = -qJD(1) * t436 + t469 * t461;
	t421 = t438 * t518 + t469 * t466;
	t417 = (t466 * t506 + t427) * t463 + (-qJD(5) * t443 + t478) * t458;
	t416 = t434 * qJD(5) + t427 * t458 - t478 * t463;
	t415 = t417 * t462 - t426 * t457 + (-t434 * t457 + t442 * t462) * qJD(6);
	t414 = -t417 * t457 - t426 * t462 + (-t434 * t462 - t442 * t457) * qJD(6);
	t1 = [t544 * r_i_i_C(1) + t545 * r_i_i_C(2) + t418 * r_i_i_C(3) - t461 * t500 + (t541 * r_i_i_C(1) + t540 * r_i_i_C(2)) * qJD(6) - t479 * qJ(3), (t421 * t462 + t436 * t504) * r_i_i_C(1) + (-t421 * t457 + t436 * t502) * r_i_i_C(2) - t542 * t518 + t467 * t466, t478, (t426 * t524 + t427 * t457 + t443 * t502) * r_i_i_C(1) + (-t426 * t530 + t427 * t462 - t443 * t504) * r_i_i_C(2) + t426 * t535 + t468 * t442, t417 * r_i_i_C(3) + (t416 * t457 - t433 * t502) * r_i_i_C(2) + (-t416 * t462 - t433 * t504) * r_i_i_C(1), r_i_i_C(1) * t414 - r_i_i_C(2) * t415; t415 * r_i_i_C(1) + t414 * r_i_i_C(2) + t416 * r_i_i_C(3) + t478 * qJ(3) + t466 * t500, (t422 * t462 + t435 * t504) * r_i_i_C(1) + (-t422 * t457 + t435 * t502) * r_i_i_C(2) + t542 * t517 + t467 * t461, t479, (-t428 * t524 + t429 * t457 + t440 * t502) * r_i_i_C(1) + (t428 * t530 + t429 * t462 - t440 * t504) * r_i_i_C(2) - t428 * t535 + t468 * t439, t419 * r_i_i_C(3) + (-t418 * t457 - t430 * t502) * r_i_i_C(2) + (t418 * t462 - t430 * t504) * r_i_i_C(1), -t545 * r_i_i_C(1) + t544 * r_i_i_C(2) + (-t540 * r_i_i_C(1) + t541 * r_i_i_C(2)) * qJD(6); 0, (t423 * t462 - t441 * t504) * r_i_i_C(1) + (-t423 * t457 - t441 * t502) * r_i_i_C(2) + t470 * t460 + t536 * t465, t515, ((t457 * t464 - t459 * t524) * r_i_i_C(1) + (t459 * t530 + t462 * t464) * r_i_i_C(2) - t459 * t535) * t514 + ((-qJD(4) * t535 + (-r_i_i_C(1) * t462 + r_i_i_C(2) * t457) * (-qJD(6) + t511)) * t464 + ((-qJD(4) * t457 + t475) * r_i_i_C(1) + (-qJD(4) * t462 + t474) * r_i_i_C(2) - t499) * t459) * t460, t425 * r_i_i_C(3) + (-t424 * t457 + t480 * t502) * r_i_i_C(2) + (t424 * t462 + t480 * t504) * r_i_i_C(1), (t471 * r_i_i_C(1) + t484 * r_i_i_C(2)) * t462 + (t484 * r_i_i_C(1) - t471 * r_i_i_C(2)) * t457;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end