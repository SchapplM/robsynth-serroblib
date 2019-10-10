% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:57
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
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
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:46
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (90->45), mult. (273->84), div. (0->0), fcn. (260->8), ass. (0->34)
	t190 = cos(pkin(6));
	t187 = sin(pkin(12));
	t189 = cos(pkin(12));
	t191 = sin(qJ(2));
	t193 = cos(qJ(2));
	t196 = t193 * t187 + t191 * t189;
	t175 = t196 * t190;
	t209 = t191 * pkin(2);
	t192 = sin(qJ(1));
	t208 = t191 * t192;
	t194 = cos(qJ(1));
	t207 = t191 * t194;
	t206 = t192 * t193;
	t205 = t193 * t194;
	t188 = sin(pkin(6));
	t204 = qJD(1) * t188;
	t203 = qJD(2) * t191;
	t202 = qJD(2) * t193;
	t201 = pkin(2) * t203;
	t200 = t190 * t202;
	t180 = t191 * t187 - t193 * t189;
	t199 = -t193 * pkin(2) + t180 * r_i_i_C(1) - pkin(1);
	t172 = t190 * t187 * t203 - t189 * t200;
	t178 = -t187 * t202 - t189 * t203;
	t198 = t194 * t172 - t192 * t178;
	t197 = t192 * t172 + t194 * t178;
	t195 = t175 * r_i_i_C(1) + t190 * t209 + (-r_i_i_C(3) - pkin(8) - qJ(3)) * t188;
	t179 = pkin(2) * t200 - t188 * qJD(3);
	t177 = t180 * qJD(2);
	t174 = t180 * t190;
	t173 = qJD(2) * t175;
	t171 = -t194 * t173 + t192 * t177 + (t174 * t192 - t194 * t196) * qJD(1);
	t170 = t192 * t173 + t194 * t177 + (t174 * t194 + t192 * t196) * qJD(1);
	t1 = [t198 * r_i_i_C(1) - t171 * r_i_i_C(2) + t192 * t201 - t194 * t179 + (t195 * t192 + t199 * t194) * qJD(1), t170 * r_i_i_C(1) + ((t175 * t194 - t180 * t192) * qJD(1) - t197) * r_i_i_C(2) + ((t190 * t208 - t205) * qJD(2) + (-t190 * t205 + t208) * qJD(1)) * pkin(2), t194 * t204, 0, 0, 0; t197 * r_i_i_C(1) + t170 * r_i_i_C(2) - t194 * t201 - t192 * t179 + (t199 * t192 - t195 * t194) * qJD(1), t171 * r_i_i_C(1) + ((t175 * t192 + t180 * t194) * qJD(1) + t198) * r_i_i_C(2) + ((-t190 * t207 - t206) * qJD(2) + (-t190 * t206 - t207) * qJD(1)) * pkin(2), t192 * t204, 0, 0, 0; 0, (-t196 * r_i_i_C(1) + t180 * r_i_i_C(2) - t209) * t188 * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:46
	% EndTime: 2019-10-10 10:57:47
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (297->77), mult. (888->135), div. (0->0), fcn. (922->10), ass. (0->59)
	t316 = cos(pkin(6));
	t313 = sin(pkin(12));
	t315 = cos(pkin(12));
	t321 = cos(qJ(2));
	t340 = qJD(2) * t321;
	t318 = sin(qJ(2));
	t341 = qJD(2) * t318;
	t352 = t313 * t341 - t315 * t340;
	t294 = t352 * t316;
	t328 = t321 * t313 + t318 * t315;
	t300 = t328 * t316;
	t304 = t318 * t313 - t321 * t315;
	t322 = cos(qJ(1));
	t319 = sin(qJ(1));
	t342 = qJD(1) * t319;
	t302 = -t313 * t340 - t315 * t341;
	t345 = t319 * t302;
	t286 = t345 - t300 * t342 + (-qJD(1) * t304 - t294) * t322;
	t314 = sin(pkin(6));
	t348 = t314 * t322;
	t333 = qJD(4) * t348;
	t353 = t286 - t333;
	t351 = r_i_i_C(3) + pkin(9);
	t350 = pkin(2) * t316;
	t349 = t314 * t319;
	t347 = t318 * t319;
	t346 = t318 * t322;
	t344 = t319 * t321;
	t343 = t321 * t322;
	t288 = t322 * t300 - t319 * t304;
	t339 = qJD(4) * t288;
	t338 = pkin(2) * t341;
	t337 = t314 * t342;
	t336 = qJD(1) * t348;
	t317 = sin(qJ(4));
	t320 = cos(qJ(4));
	t332 = t317 * r_i_i_C(1) + t320 * r_i_i_C(2);
	t299 = t304 * t316;
	t331 = -t322 * t299 - t319 * t328;
	t330 = t319 * t299 - t322 * t328;
	t329 = t319 * t300 + t322 * t304;
	t327 = t320 * r_i_i_C(1) - t317 * r_i_i_C(2) + pkin(3);
	t326 = qJD(4) * t332;
	t325 = t337 - t339;
	t324 = qJD(2) * t328;
	t323 = t304 * qJD(2);
	t283 = -t288 * qJD(1) + t319 * t294 + t322 * t302;
	t312 = t321 * pkin(2) + pkin(1);
	t308 = t320 * t333;
	t303 = -t314 * qJD(3) + t340 * t350;
	t301 = t318 * t350 + (-pkin(8) - qJ(3)) * t314;
	t298 = t328 * t314;
	t295 = t316 * t324;
	t292 = t352 * t314;
	t285 = t330 * qJD(1) - t322 * t295 + t319 * t323;
	t282 = t331 * qJD(1) - t319 * t295 - t322 * t323;
	t280 = t317 * t336 + t283 * t320 + (t317 * t329 + t320 * t349) * qJD(4);
	t279 = t320 * t336 - t283 * t317 + (-t317 * t349 + t320 * t329) * qJD(4);
	t1 = [(-t286 * t320 + t317 * t339 + t308) * r_i_i_C(1) + (t353 * t317 + t320 * t339) * r_i_i_C(2) - t286 * pkin(3) + t319 * t338 - t322 * t303 + t351 * t285 + (-t322 * t312 + (-t332 * t314 + t301) * t319) * qJD(1), t351 * t283 - t330 * t326 - t327 * t282 + ((t316 * t347 - t343) * qJD(2) + (-t316 * t343 + t347) * qJD(1)) * pkin(2), t336, t279 * r_i_i_C(1) - t280 * r_i_i_C(2), 0, 0; -t322 * t338 + t283 * pkin(3) + t280 * r_i_i_C(1) + t279 * r_i_i_C(2) - t319 * t303 + t351 * t282 + (-t301 * t322 - t312 * t319) * qJD(1), -t351 * (qJD(1) * t329 + t322 * t294 - t345) - t331 * t326 + t327 * t285 + ((-t316 * t346 - t344) * qJD(2) + (-t316 * t344 - t346) * qJD(1)) * pkin(2), t337, t308 * r_i_i_C(2) + (r_i_i_C(1) * t325 - t286 * r_i_i_C(2)) * t320 + (-t353 * r_i_i_C(1) - t325 * r_i_i_C(2)) * t317, 0, 0; 0, -t351 * t292 + (t304 * t326 - t324 * t327 - t338) * t314, 0, t332 * t292 + ((-t298 * t320 - t316 * t317) * r_i_i_C(1) + (t298 * t317 - t316 * t320) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:48
	% EndTime: 2019-10-10 10:57:49
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (845->121), mult. (2486->210), div. (0->0), fcn. (2730->12), ass. (0->77)
	t453 = cos(pkin(6));
	t451 = sin(pkin(12));
	t506 = cos(pkin(12));
	t508 = cos(qJ(2));
	t480 = t508 * t506;
	t456 = sin(qJ(2));
	t494 = qJD(2) * t456;
	t512 = -qJD(2) * t480 + t451 * t494;
	t431 = t512 * t453;
	t483 = t456 * t506;
	t469 = t508 * t451 + t483;
	t436 = t469 * t453;
	t484 = t508 * qJD(2);
	t438 = -qJD(2) * t483 - t451 * t484;
	t457 = sin(qJ(1));
	t460 = cos(qJ(1));
	t468 = -t456 * t451 + t480;
	t496 = qJD(1) * t457;
	t407 = t436 * t496 - t457 * t438 + (-qJD(1) * t468 + t431) * t460;
	t455 = sin(qJ(4));
	t459 = cos(qJ(4));
	t477 = t460 * t436 + t457 * t468;
	t452 = sin(pkin(6));
	t487 = t452 * t496;
	t500 = t452 * t460;
	t490 = t459 * t500;
	t401 = (-qJD(4) * t477 + t487) * t455 - qJD(4) * t490 - t407 * t459;
	t432 = qJD(2) * t436;
	t495 = qJD(1) * t460;
	t467 = t468 * t453;
	t513 = qJD(1) * t467 + t468 * qJD(2);
	t408 = -t460 * t432 - t513 * t457 - t469 * t495;
	t454 = sin(qJ(5));
	t458 = cos(qJ(5));
	t519 = t401 * t454 + t408 * t458;
	t518 = -t401 * t458 + t408 * t454;
	t414 = t455 * t500 - t459 * t477;
	t418 = -t457 * t469 + t460 * t467;
	t517 = -t414 * t454 + t418 * t458;
	t516 = t414 * t458 + t418 * t454;
	t472 = qJD(5) * (t454 * r_i_i_C(1) + t458 * r_i_i_C(2));
	t475 = t458 * r_i_i_C(1) - t454 * r_i_i_C(2) + pkin(4);
	t509 = r_i_i_C(3) + pkin(10);
	t515 = (t475 * t455 - t509 * t459) * qJD(4) + t459 * t472;
	t510 = t509 * t455 + t475 * t459 + pkin(3);
	t507 = pkin(2) * t453;
	t501 = t452 * t457;
	t498 = t456 * t457;
	t497 = t456 * t460;
	t493 = qJD(5) * t454;
	t492 = qJD(5) * t458;
	t491 = pkin(2) * t494;
	t489 = t508 * t457;
	t488 = t508 * t460;
	t486 = t452 * t495;
	t435 = t469 * t452;
	t424 = t435 * t459 + t453 * t455;
	t478 = -t435 * t455 + t453 * t459;
	t476 = -t457 * t436 + t460 * t468;
	t473 = -t455 * t476 + t459 * t501;
	t416 = t455 * t501 + t459 * t476;
	t463 = -t477 * qJD(1) + t457 * t431 + t460 * t438;
	t462 = t414 * qJD(4) + t407 * t455 + t459 * t487;
	t450 = t508 * pkin(2) + pkin(1);
	t439 = -t452 * qJD(3) + t484 * t507;
	t437 = t456 * t507 + (-pkin(8) - qJ(3)) * t452;
	t434 = t468 * t452;
	t430 = qJD(2) * t435;
	t429 = t512 * t452;
	t421 = -t457 * t467 - t460 * t469;
	t411 = t478 * qJD(4) - t429 * t459;
	t405 = -t457 * t432 + t513 * t460 - t469 * t496;
	t399 = t473 * qJD(4) + t455 * t486 + t459 * t463;
	t398 = t416 * qJD(4) + t455 * t463 - t459 * t486;
	t397 = t399 * t458 + t405 * t454 + (-t416 * t454 - t421 * t458) * qJD(5);
	t396 = -t399 * t454 + t405 * t458 + (-t416 * t458 + t421 * t454) * qJD(5);
	t1 = [t518 * r_i_i_C(1) + t519 * r_i_i_C(2) - t401 * pkin(4) + t407 * pkin(3) + t408 * pkin(9) + t457 * t491 - t460 * t439 + t509 * t462 + (t517 * r_i_i_C(1) - t516 * r_i_i_C(2)) * qJD(5) + (t457 * t437 - t460 * t450) * qJD(1), (t454 * t463 + t476 * t492) * r_i_i_C(1) + (t458 * t463 - t476 * t493) * r_i_i_C(2) + t463 * pkin(9) - t510 * t405 + ((t453 * t498 - t488) * qJD(2) + (-t453 * t488 + t498) * qJD(1)) * pkin(2) - t515 * t421, t486, -t475 * t398 + t509 * t399 - t473 * t472, t396 * r_i_i_C(1) - t397 * r_i_i_C(2), 0; -t460 * t491 + t463 * pkin(3) + t399 * pkin(4) + t405 * pkin(9) + t397 * r_i_i_C(1) + t396 * r_i_i_C(2) - t457 * t439 + t509 * t398 + (-t437 * t460 - t450 * t457) * qJD(1), (-t407 * t454 + t477 * t492) * r_i_i_C(1) + (-t407 * t458 - t477 * t493) * r_i_i_C(2) - t407 * pkin(9) + t510 * t408 + ((-t453 * t497 - t489) * qJD(2) + (-t453 * t489 - t497) * qJD(1)) * pkin(2) - t515 * t418, t487, t509 * t401 - (-t455 * t477 - t490) * t472 + t475 * t462, -t519 * r_i_i_C(1) + t518 * r_i_i_C(2) + (t516 * r_i_i_C(1) + t517 * r_i_i_C(2)) * qJD(5), 0; 0, (-t429 * t454 + t435 * t492) * r_i_i_C(1) + (-t429 * t458 - t435 * t493) * r_i_i_C(2) - t429 * pkin(9) - t452 * t491 - t510 * t430 - t515 * t434, 0, t509 * t411 - t478 * t472 + t475 * (-t424 * qJD(4) + t429 * t455), (-t411 * t454 + t430 * t458) * r_i_i_C(1) + (-t411 * t458 - t430 * t454) * r_i_i_C(2) + ((-t424 * t458 + t434 * t454) * r_i_i_C(1) + (t424 * t454 + t434 * t458) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:48
	% EndTime: 2019-10-10 10:57:50
	% DurationCPUTime: 1.43s
	% Computational Cost: add. (1327->143), mult. (3394->238), div. (0->0), fcn. (3744->14), ass. (0->88)
	t491 = cos(pkin(6));
	t489 = sin(pkin(12));
	t555 = cos(pkin(12));
	t558 = cos(qJ(2));
	t519 = t558 * t555;
	t494 = sin(qJ(2));
	t541 = qJD(2) * t494;
	t562 = -qJD(2) * t519 + t489 * t541;
	t464 = t562 * t491;
	t530 = t494 * t555;
	t509 = t558 * t489 + t530;
	t469 = t509 * t491;
	t531 = qJD(2) * t558;
	t471 = -qJD(2) * t530 - t489 * t531;
	t495 = sin(qJ(1));
	t498 = cos(qJ(1));
	t508 = -t494 * t489 + t519;
	t543 = qJD(1) * t495;
	t440 = t469 * t543 - t495 * t471 + (-qJD(1) * t508 + t464) * t498;
	t493 = sin(qJ(4));
	t497 = cos(qJ(4));
	t516 = t469 * t498 + t495 * t508;
	t490 = sin(pkin(6));
	t534 = t490 * t543;
	t550 = t490 * t498;
	t537 = t497 * t550;
	t434 = (-qJD(4) * t516 + t534) * t493 - qJD(4) * t537 - t440 * t497;
	t507 = t491 * t508;
	t451 = -t495 * t509 + t498 * t507;
	t487 = qJD(5) + qJD(6);
	t567 = t451 * t487 - t434;
	t496 = cos(qJ(5));
	t483 = pkin(5) * t496 + pkin(4);
	t488 = qJ(5) + qJ(6);
	t485 = sin(t488);
	t486 = cos(t488);
	t514 = r_i_i_C(1) * t486 - r_i_i_C(2) * t485 + t483;
	t556 = r_i_i_C(3) + pkin(11) + pkin(10);
	t492 = sin(qJ(5));
	t561 = qJD(5) * t492 * pkin(5) + (r_i_i_C(1) * t485 + r_i_i_C(2) * t486) * t487;
	t566 = t561 * t497 + (t514 * t493 - t556 * t497) * qJD(4);
	t465 = qJD(2) * t469;
	t542 = qJD(1) * t498;
	t563 = qJD(1) * t507 + qJD(2) * t508;
	t441 = -t498 * t465 - t563 * t495 - t509 * t542;
	t447 = t493 * t550 - t497 * t516;
	t564 = -t447 * t487 + t441;
	t559 = t556 * t493 + t514 * t497 + pkin(3);
	t557 = pkin(2) * t491;
	t553 = t485 * t487;
	t552 = t486 * t487;
	t551 = t490 * t495;
	t548 = t494 * t495;
	t547 = t494 * t498;
	t438 = -t495 * t465 + t563 * t498 - t509 * t543;
	t515 = -t469 * t495 + t498 * t508;
	t449 = t493 * t551 + t497 * t515;
	t526 = -t449 * t487 + t438;
	t502 = -t516 * qJD(1) + t495 * t464 + t498 * t471;
	t512 = -t493 * t515 + t497 * t551;
	t533 = t490 * t542;
	t432 = t512 * qJD(4) + t493 * t533 + t497 * t502;
	t454 = -t495 * t507 - t498 * t509;
	t529 = t454 * t487 - t432;
	t427 = t529 * t485 + t526 * t486;
	t428 = t526 * t485 - t529 * t486;
	t546 = t427 * r_i_i_C(1) - t428 * r_i_i_C(2);
	t545 = (t485 * t567 - t486 * t564) * r_i_i_C(1) + (t485 * t564 + t486 * t567) * r_i_i_C(2);
	t468 = t509 * t490;
	t457 = t468 * t497 + t491 * t493;
	t463 = qJD(2) * t468;
	t521 = t457 * t487 - t463;
	t462 = t562 * t490;
	t517 = -t468 * t493 + t491 * t497;
	t444 = t517 * qJD(4) - t462 * t497;
	t467 = t508 * t490;
	t522 = t467 * t487 - t444;
	t544 = (t522 * t485 - t521 * t486) * r_i_i_C(1) + (t521 * t485 + t522 * t486) * r_i_i_C(2);
	t540 = qJD(5) * t496;
	t539 = pkin(2) * t541;
	t536 = t558 * t495;
	t535 = t558 * t498;
	t501 = t447 * qJD(4) + t440 * t493 + t497 * t534;
	t484 = t558 * pkin(2) + pkin(1);
	t472 = -t490 * qJD(3) + t531 * t557;
	t470 = t494 * t557 + (-pkin(8) - qJ(3)) * t490;
	t431 = t449 * qJD(4) + t493 * t502 - t497 * t533;
	t1 = [t495 * t539 + t440 * pkin(3) + t441 * pkin(9) - t434 * t483 - t498 * t472 + (r_i_i_C(1) * t567 + r_i_i_C(2) * t564) * t486 + (r_i_i_C(1) * t564 - r_i_i_C(2) * t567) * t485 + t556 * t501 + (t470 * t495 - t484 * t498) * qJD(1) + (t441 * t492 + (-t447 * t492 + t451 * t496) * qJD(5)) * pkin(5), (t485 * t502 + t515 * t552) * r_i_i_C(1) + (t486 * t502 - t515 * t553) * r_i_i_C(2) + t502 * pkin(9) + (t492 * t502 + t515 * t540) * pkin(5) - t559 * t438 + ((t491 * t548 - t535) * qJD(2) + (-t491 * t535 + t548) * qJD(1)) * pkin(2) - t566 * t454, t533, -t514 * t431 + t556 * t432 - t512 * t561, (-t432 * t492 + t438 * t496 + (-t449 * t496 + t454 * t492) * qJD(5)) * pkin(5) + t546, t546; -t498 * t539 + t502 * pkin(3) + t438 * pkin(9) + t428 * r_i_i_C(1) + t427 * r_i_i_C(2) + t432 * t483 - t495 * t472 + t556 * t431 + (-t470 * t498 - t484 * t495) * qJD(1) + (t438 * t492 + (-t449 * t492 - t454 * t496) * qJD(5)) * pkin(5), (-t440 * t485 + t516 * t552) * r_i_i_C(1) + (-t440 * t486 - t516 * t553) * r_i_i_C(2) - t440 * pkin(9) + (-t440 * t492 + t516 * t540) * pkin(5) + t559 * t441 + ((-t491 * t547 - t536) * qJD(2) + (-t491 * t536 - t547) * qJD(1)) * pkin(2) - t566 * t451, t534, t556 * t434 - t561 * (-t493 * t516 - t537) + t514 * t501, (-t434 * t492 - t441 * t496 + (t447 * t496 + t451 * t492) * qJD(5)) * pkin(5) + t545, t545; 0, (-t462 * t485 + t468 * t552) * r_i_i_C(1) + (-t462 * t486 - t468 * t553) * r_i_i_C(2) - t462 * pkin(9) - t490 * t539 + (-t462 * t492 + t468 * t540) * pkin(5) - t559 * t463 - t566 * t467, 0, t556 * t444 - t561 * t517 + t514 * (-t457 * qJD(4) + t462 * t493), (-t444 * t492 + t463 * t496 + (-t457 * t496 + t467 * t492) * qJD(5)) * pkin(5) + t544, t544;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end