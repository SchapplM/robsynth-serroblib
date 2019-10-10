% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:37
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:18
	% EndTime: 2019-10-10 10:37:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:19
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
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:19
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
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:19
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (90->45), mult. (273->84), div. (0->0), fcn. (260->8), ass. (0->34)
	t190 = cos(pkin(6));
	t187 = sin(pkin(11));
	t189 = cos(pkin(11));
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
	% StartTime: 2019-10-10 10:37:20
	% EndTime: 2019-10-10 10:37:21
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (297->77), mult. (888->135), div. (0->0), fcn. (922->10), ass. (0->59)
	t316 = cos(pkin(6));
	t313 = sin(pkin(11));
	t315 = cos(pkin(11));
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
	t1 = [(-t286 * t320 + t317 * t339 + t308) * r_i_i_C(1) + (t353 * t317 + t320 * t339) * r_i_i_C(2) - t286 * pkin(3) + t319 * t338 - t322 * t303 + t351 * t285 + (-t322 * t312 + (-t332 * t314 + t301) * t319) * qJD(1), t351 * t283 - t330 * t326 - t327 * t282 + ((t316 * t347 - t343) * qJD(2) + (-t316 * t343 + t347) * qJD(1)) * pkin(2), t336, t279 * r_i_i_C(1) - t280 * r_i_i_C(2), 0, 0; -t322 * t338 + t283 * pkin(3) + t280 * r_i_i_C(1) + t279 * r_i_i_C(2) - t319 * t303 + t351 * t282 + (-t301 * t322 - t312 * t319) * qJD(1), -t351 * (t329 * qJD(1) + t322 * t294 - t345) - t331 * t326 + t327 * t285 + ((-t316 * t346 - t344) * qJD(2) + (-t316 * t344 - t346) * qJD(1)) * pkin(2), t337, t308 * r_i_i_C(2) + (t325 * r_i_i_C(1) - t286 * r_i_i_C(2)) * t320 + (-t353 * r_i_i_C(1) - t325 * r_i_i_C(2)) * t317, 0, 0; 0, -t351 * t292 + (t304 * t326 - t324 * t327 - t338) * t314, 0, t332 * t292 + ((-t298 * t320 - t316 * t317) * r_i_i_C(1) + (t298 * t317 - t316 * t320) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:22
	% EndTime: 2019-10-10 10:37:23
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (845->121), mult. (2486->210), div. (0->0), fcn. (2730->12), ass. (0->77)
	t453 = cos(pkin(6));
	t451 = sin(pkin(11));
	t506 = cos(pkin(11));
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
	% StartTime: 2019-10-10 10:37:25
	% EndTime: 2019-10-10 10:37:27
	% DurationCPUTime: 1.53s
	% Computational Cost: add. (1486->142), mult. (4318->232), div. (0->0), fcn. (4849->12), ass. (0->90)
	t554 = cos(pkin(6));
	t552 = sin(pkin(11));
	t619 = cos(pkin(11));
	t623 = cos(qJ(2));
	t591 = t623 * t619;
	t557 = sin(qJ(2));
	t606 = qJD(2) * t557;
	t628 = -qJD(2) * t591 + t552 * t606;
	t532 = t628 * t554;
	t594 = t557 * t619;
	t573 = t623 * t552 + t594;
	t537 = t573 * t554;
	t595 = t623 * qJD(2);
	t539 = -qJD(2) * t594 - t552 * t595;
	t558 = sin(qJ(1));
	t561 = cos(qJ(1));
	t572 = -t557 * t552 + t591;
	t608 = qJD(1) * t558;
	t507 = t537 * t608 - t558 * t539 + (-qJD(1) * t572 + t532) * t561;
	t556 = sin(qJ(4));
	t560 = cos(qJ(4));
	t581 = t561 * t537 + t558 * t572;
	t553 = sin(pkin(6));
	t598 = t553 * t608;
	t613 = t553 * t561;
	t602 = t560 * t613;
	t495 = (-qJD(4) * t581 + t598) * t556 - qJD(4) * t602 - t507 * t560;
	t533 = qJD(2) * t537;
	t607 = qJD(1) * t561;
	t567 = t572 * t554;
	t629 = qJD(1) * t567 + t572 * qJD(2);
	t508 = -t561 * t533 - t629 * t558 - t573 * t607;
	t555 = sin(qJ(5));
	t559 = cos(qJ(5));
	t514 = t556 * t613 - t560 * t581;
	t518 = -t558 * t573 + t561 * t567;
	t630 = t514 * t559 + t518 * t555;
	t635 = t630 * qJD(5) - t495 * t555 - t508 * t559;
	t631 = t514 * t555 - t518 * t559;
	t634 = t631 * qJD(5) + t495 * t559 - t508 * t555;
	t624 = r_i_i_C(2) + pkin(10);
	t574 = qJD(4) * (-pkin(4) * t556 + t624 * t560);
	t627 = t560 * pkin(4) + t624 * t556 + pkin(3);
	t620 = r_i_i_C(3) + qJ(6);
	t625 = r_i_i_C(1) + pkin(5);
	t571 = t620 * t555 + t625 * t559 + pkin(4);
	t622 = pkin(2) * t554;
	t614 = t553 * t558;
	t612 = t555 * t560;
	t610 = t557 * t558;
	t609 = t557 * t561;
	t605 = qJD(4) * t556;
	t604 = qJD(5) * t560;
	t603 = pkin(2) * t606;
	t600 = t623 * t558;
	t599 = t623 * t561;
	t597 = t553 * t607;
	t521 = -t558 * t567 - t561 * t573;
	t563 = -t581 * qJD(1) + t558 * t532 + t561 * t539;
	t590 = t521 * t604 - t563;
	t589 = t518 * t604 + t507;
	t530 = t628 * t553;
	t535 = t572 * t553;
	t588 = -t535 * t604 - t530;
	t580 = -t558 * t537 + t561 * t572;
	t516 = t556 * t614 + t560 * t580;
	t585 = t516 * t559 - t521 * t555;
	t584 = -t516 * t555 - t521 * t559;
	t536 = t573 * t553;
	t524 = t536 * t560 + t554 * t556;
	t583 = t524 * t559 - t535 * t555;
	t582 = -t536 * t556 + t554 * t560;
	t578 = -t556 * t580 + t560 * t614;
	t505 = -t558 * t533 + t629 * t561 - t573 * t608;
	t570 = qJD(5) * t580 - t505 * t560 - t521 * t605;
	t569 = qJD(5) * t581 + t508 * t560 - t518 * t605;
	t531 = qJD(2) * t536;
	t568 = qJD(5) * t536 - t531 * t560 - t535 * t605;
	t564 = qJD(6) * t555 + (-t625 * t555 + t620 * t559) * qJD(5);
	t562 = t514 * qJD(4) + t507 * t556 + t560 * t598;
	t551 = t623 * pkin(2) + pkin(1);
	t540 = -t553 * qJD(3) + t595 * t622;
	t538 = t557 * t622 + (-pkin(8) - qJ(3)) * t553;
	t511 = t582 * qJD(4) - t530 * t560;
	t498 = t583 * qJD(5) + t511 * t555 - t531 * t559;
	t493 = t578 * qJD(4) + t556 * t597 + t560 * t563;
	t492 = t516 * qJD(4) + t556 * t563 - t560 * t597;
	t483 = t584 * qJD(5) + t493 * t559 + t505 * t555;
	t482 = t585 * qJD(5) + t493 * t555 - t505 * t559;
	t1 = [t631 * qJD(6) - t495 * pkin(4) + t507 * pkin(3) + t508 * pkin(9) + t558 * t603 - t561 * t540 + t624 * t562 - t625 * t634 + t620 * t635 + (t558 * t538 - t561 * t551) * qJD(1), -(-t521 * t612 + t559 * t580) * qJD(6) + t563 * pkin(9) + t625 * (-t590 * t555 + t570 * t559) + t620 * (t570 * t555 + t590 * t559) + t521 * t574 - t627 * t505 + ((t554 * t610 - t599) * qJD(2) + (-t554 * t599 + t610) * qJD(1)) * pkin(2), t597, -t571 * t492 + t624 * t493 + t564 * t578, t585 * qJD(6) - t625 * t482 + t620 * t483, t482; -t584 * qJD(6) + t493 * pkin(4) + t563 * pkin(3) + t505 * pkin(9) - t561 * t603 - t558 * t540 + t624 * t492 + t625 * t483 + t620 * t482 + (-t561 * t538 - t558 * t551) * qJD(1), -(-t518 * t612 + t559 * t581) * qJD(6) - t507 * pkin(9) + t625 * (-t589 * t555 + t569 * t559) + t620 * (t569 * t555 + t589 * t559) + t518 * t574 + t627 * t508 + ((-t554 * t609 - t600) * qJD(2) + (-t554 * t600 - t609) * qJD(1)) * pkin(2), t598, t624 * t495 + t564 * (-t556 * t581 - t602) + t571 * t562, -t630 * qJD(6) + t620 * t634 + t625 * t635, -t635; 0, -(-t535 * t612 + t536 * t559) * qJD(6) - t530 * pkin(9) - t553 * t603 + t625 * (t588 * t555 + t568 * t559) + t620 * (t568 * t555 - t588 * t559) + t535 * t574 - t627 * t531, 0, t624 * t511 + t564 * t582 + t571 * (-t524 * qJD(4) + t530 * t556), t583 * qJD(6) + t620 * (t511 * t559 + t531 * t555 + (-t524 * t555 - t535 * t559) * qJD(5)) - t625 * t498, t498;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end