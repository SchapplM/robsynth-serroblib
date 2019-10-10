% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
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
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
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
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:39
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
	% StartTime: 2019-10-10 10:11:39
	% EndTime: 2019-10-10 10:11:40
	% DurationCPUTime: 0.43s
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
	% StartTime: 2019-10-10 10:11:41
	% EndTime: 2019-10-10 10:11:41
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (661->92), mult. (1960->157), div. (0->0), fcn. (2112->12), ass. (0->65)
	t397 = sin(pkin(6));
	t405 = cos(qJ(1));
	t437 = t397 * t405;
	t399 = cos(pkin(6));
	t396 = sin(pkin(11));
	t404 = cos(qJ(2));
	t439 = cos(pkin(11));
	t422 = qJD(2) * t439;
	t401 = sin(qJ(2));
	t430 = qJD(2) * t401;
	t444 = t396 * t430 - t404 * t422;
	t373 = t444 * t399;
	t384 = -t404 * t396 - t401 * t439;
	t379 = t384 * t399;
	t429 = qJD(2) * t404;
	t381 = -t396 * t429 - t401 * t422;
	t402 = sin(qJ(1));
	t410 = -t401 * t396 + t404 * t439;
	t432 = qJD(1) * t402;
	t442 = -t379 * t432 - t402 * t381 + (-qJD(1) * t410 + t373) * t405;
	t446 = qJD(4) * t437 + t442;
	t366 = -t405 * t379 + t402 * t410;
	t426 = t397 * t432;
	t445 = -qJD(4) * t366 + t426;
	t400 = sin(qJ(4));
	t403 = cos(qJ(4));
	t443 = t445 * t400 - t403 * t446;
	t395 = sin(pkin(12));
	t398 = cos(pkin(12));
	t416 = t398 * r_i_i_C(1) - t395 * r_i_i_C(2) + pkin(4);
	t440 = r_i_i_C(3) + qJ(5);
	t407 = t440 * t400 + t416 * t403 + pkin(3);
	t441 = pkin(2) * t399;
	t438 = t397 * t402;
	t436 = t401 * t402;
	t435 = t401 * t405;
	t434 = t402 * t404;
	t433 = t404 * t405;
	t431 = qJD(1) * t405;
	t427 = pkin(2) * t430;
	t425 = t397 * t431;
	t377 = t384 * t397;
	t418 = -t377 * t403 + t399 * t400;
	t378 = t410 * t399;
	t417 = t405 * t378 + t402 * t384;
	t415 = t395 * r_i_i_C(1) + t398 * r_i_i_C(2) + pkin(9);
	t368 = t402 * t379 + t405 * t410;
	t414 = -t368 * t400 + t403 * t438;
	t413 = t368 * t403 + t400 * t438;
	t409 = qJD(2) * t384;
	t408 = t410 * qJD(2);
	t353 = -t400 * t446 - t445 * t403;
	t359 = -t366 * qJD(1) + t402 * t373 + t405 * t381;
	t406 = t400 * qJD(5) + (-t416 * t400 + t440 * t403) * qJD(4);
	t394 = t404 * pkin(2) + pkin(1);
	t382 = -t397 * qJD(3) + t429 * t441;
	t380 = t401 * t441 + (-pkin(8) - qJ(3)) * t397;
	t374 = t399 * t409;
	t371 = t444 * t397;
	t363 = t418 * qJD(4) - t371 * t400;
	t361 = t405 * t374 - t378 * t432 + t384 * t431 - t402 * t408;
	t358 = t417 * qJD(1) + t402 * t374 + t405 * t408;
	t352 = t414 * qJD(4) + t359 * t403 + t400 * t425;
	t351 = t413 * qJD(4) + t359 * t400 - t403 * t425;
	t1 = [(t361 * t395 - t398 * t443) * r_i_i_C(1) + (t361 * t398 + t395 * t443) * r_i_i_C(2) - t443 * pkin(4) - (t366 * t400 + t403 * t437) * qJD(5) + t442 * pkin(3) + t361 * pkin(9) + t402 * t427 - t405 * t382 - t440 * t353 + (t402 * t380 - t405 * t394) * qJD(1), t415 * t359 + t406 * (-t402 * t378 + t405 * t384) - t407 * t358 + ((t399 * t436 - t433) * qJD(2) + (-t399 * t433 + t436) * qJD(1)) * pkin(2), t425, t413 * qJD(5) - t416 * t351 + t440 * t352, t351, 0; (t352 * t398 + t358 * t395) * r_i_i_C(1) + (-t352 * t395 + t358 * t398) * r_i_i_C(2) + t352 * pkin(4) - t414 * qJD(5) + t359 * pkin(3) + t358 * pkin(9) - t405 * t427 - t402 * t382 + t440 * t351 + (-t405 * t380 - t402 * t394) * qJD(1), -t415 * t442 + t406 * t417 + t407 * t361 + ((-t399 * t435 - t434) * qJD(2) + (-t399 * t434 - t435) * qJD(1)) * pkin(2), t426, -(-t366 * t403 + t400 * t437) * qJD(5) + t440 * t443 - t416 * t353, t353, 0; 0, -t415 * t371 + (t406 * t410 + t407 * t409 - t427) * t397, 0, t418 * qJD(5) + t440 * (-t371 * t403 + (t377 * t400 + t399 * t403) * qJD(4)) - t416 * t363, t363, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:41
	% EndTime: 2019-10-10 10:11:42
	% DurationCPUTime: 1.10s
	% Computational Cost: add. (1046->119), mult. (2775->199), div. (0->0), fcn. (3058->14), ass. (0->83)
	t463 = sin(pkin(6));
	t470 = cos(qJ(1));
	t510 = t463 * t470;
	t464 = cos(pkin(6));
	t462 = sin(pkin(11));
	t514 = cos(pkin(11));
	t517 = cos(qJ(2));
	t489 = t517 * t514;
	t467 = sin(qJ(2));
	t504 = qJD(2) * t467;
	t521 = -qJD(2) * t489 + t462 * t504;
	t434 = t521 * t464;
	t492 = t467 * t514;
	t477 = t517 * t462 + t492;
	t439 = t477 * t464;
	t493 = t517 * qJD(2);
	t441 = -qJD(2) * t492 - t462 * t493;
	t468 = sin(qJ(1));
	t476 = -t467 * t462 + t489;
	t506 = qJD(1) * t468;
	t520 = t439 * t506 - t468 * t441 + (-qJD(1) * t476 + t434) * t470;
	t528 = qJD(4) * t510 + t520;
	t422 = t470 * t439 + t468 * t476;
	t497 = t463 * t506;
	t527 = -qJD(4) * t422 + t497;
	t466 = sin(qJ(4));
	t469 = cos(qJ(4));
	t416 = t422 * t469 - t466 * t510;
	t475 = t476 * t464;
	t421 = -t468 * t477 + t470 * t475;
	t460 = pkin(12) + qJ(6);
	t458 = sin(t460);
	t459 = cos(t460);
	t526 = t416 * t458 + t421 * t459;
	t525 = t416 * t459 - t421 * t458;
	t487 = t458 * r_i_i_C(1) + t459 * r_i_i_C(2);
	t481 = qJD(6) * t487;
	t456 = cos(pkin(12)) * pkin(5) + pkin(4);
	t488 = t459 * r_i_i_C(1) - t458 * r_i_i_C(2);
	t485 = t456 + t488;
	t515 = r_i_i_C(3) + pkin(10) + qJ(5);
	t524 = (t485 * t466 - t515 * t469) * qJD(4) - t466 * qJD(5) + t469 * t481;
	t522 = qJD(1) * t475 + t476 * qJD(2);
	t404 = t527 * t466 - t528 * t469;
	t518 = t515 * t466 + t485 * t469 + pkin(3);
	t516 = pkin(2) * t464;
	t511 = t463 * t468;
	t508 = t467 * t468;
	t507 = t467 * t470;
	t505 = qJD(1) * t470;
	t501 = pkin(2) * t504;
	t500 = sin(pkin(12)) * pkin(5) + pkin(9);
	t499 = t517 * t468;
	t498 = t517 * t470;
	t496 = t463 * t505;
	t438 = t477 * t463;
	t427 = t438 * t469 + t464 * t466;
	t486 = -t438 * t466 + t464 * t469;
	t425 = -t468 * t439 + t470 * t476;
	t483 = t422 * t466 + t469 * t510;
	t418 = -t425 * t466 + t469 * t511;
	t419 = t425 * t469 + t466 * t511;
	t482 = qJD(6) * t488;
	t478 = t487 + t500;
	t403 = -t528 * t466 - t527 * t469;
	t409 = -t422 * qJD(1) + t468 * t434 + t470 * t441;
	t457 = t517 * pkin(2) + pkin(1);
	t442 = -t463 * qJD(3) + t493 * t516;
	t440 = t467 * t516 + (-pkin(8) - qJ(3)) * t463;
	t437 = t476 * t463;
	t435 = qJD(2) * t439;
	t433 = qJD(2) * t438;
	t432 = t521 * t463;
	t424 = -t468 * t475 - t470 * t477;
	t414 = t486 * qJD(4) - t432 * t469;
	t413 = t427 * qJD(4) - t432 * t466;
	t411 = -t470 * t435 - t522 * t468 - t477 * t505;
	t408 = -t468 * t435 + t522 * t470 - t477 * t506;
	t402 = t418 * qJD(4) + t409 * t469 + t466 * t496;
	t401 = t419 * qJD(4) + t409 * t466 - t469 * t496;
	t400 = t402 * t459 + t408 * t458 + (-t419 * t458 - t424 * t459) * qJD(6);
	t399 = -t402 * t458 + t408 * t459 + (-t419 * t459 + t424 * t458) * qJD(6);
	t1 = [-t483 * qJD(5) + t520 * pkin(3) + t468 * t501 - t470 * t442 - t485 * t404 + t478 * t411 - t515 * t403 + (t526 * r_i_i_C(1) + t525 * r_i_i_C(2)) * qJD(6) + (t468 * t440 - t470 * t457) * qJD(1), t425 * t482 + t478 * t409 - t518 * t408 + ((t464 * t508 - t498) * qJD(2) + (-t464 * t498 + t508) * qJD(1)) * pkin(2) - t524 * t424, t496, t419 * qJD(5) - t485 * t401 + t515 * t402 - t418 * t481, t401, t399 * r_i_i_C(1) - t400 * r_i_i_C(2); -t470 * t501 + t409 * pkin(3) + t400 * r_i_i_C(1) + t399 * r_i_i_C(2) - t418 * qJD(5) + t402 * t456 - t468 * t442 + t500 * t408 + t515 * t401 + (-t440 * t470 - t457 * t468) * qJD(1), t422 * t482 - t478 * t520 + t518 * t411 + ((-t464 * t507 - t499) * qJD(2) + (-t464 * t499 - t507) * qJD(1)) * pkin(2) - t524 * t421, t497, t416 * qJD(5) - t485 * t403 + t515 * t404 + t483 * t481, t403, (-t404 * t458 - t411 * t459) * r_i_i_C(1) + (-t404 * t459 + t411 * t458) * r_i_i_C(2) + (-t525 * r_i_i_C(1) + t526 * r_i_i_C(2)) * qJD(6); 0, -t478 * t432 - t433 * t518 - t524 * t437 + t438 * t482 - t463 * t501, 0, t427 * qJD(5) - t485 * t413 + t515 * t414 - t486 * t481, t413, (-t414 * t458 + t433 * t459) * r_i_i_C(1) + (-t414 * t459 - t433 * t458) * r_i_i_C(2) + ((-t427 * t459 + t437 * t458) * r_i_i_C(1) + (t427 * t458 + t437 * t459) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end