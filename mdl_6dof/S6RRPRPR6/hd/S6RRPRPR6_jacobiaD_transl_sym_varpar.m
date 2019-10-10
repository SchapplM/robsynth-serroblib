% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
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
	% StartTime: 2019-10-10 10:13:33
	% EndTime: 2019-10-10 10:13:33
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
	% StartTime: 2019-10-10 10:13:33
	% EndTime: 2019-10-10 10:13:33
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
	% StartTime: 2019-10-10 10:13:34
	% EndTime: 2019-10-10 10:13:34
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
	t1 = [(-t286 * t320 + t317 * t339 + t308) * r_i_i_C(1) + (t353 * t317 + t320 * t339) * r_i_i_C(2) - t286 * pkin(3) + t319 * t338 - t322 * t303 + t351 * t285 + (-t322 * t312 + (-t332 * t314 + t301) * t319) * qJD(1), t351 * t283 - t330 * t326 - t327 * t282 + ((t316 * t347 - t343) * qJD(2) + (-t316 * t343 + t347) * qJD(1)) * pkin(2), t336, t279 * r_i_i_C(1) - t280 * r_i_i_C(2), 0, 0; -t322 * t338 + t283 * pkin(3) + t280 * r_i_i_C(1) + t279 * r_i_i_C(2) - t319 * t303 + t351 * t282 + (-t301 * t322 - t312 * t319) * qJD(1), -t351 * (t329 * qJD(1) + t322 * t294 - t345) - t331 * t326 + t327 * t285 + ((-t316 * t346 - t344) * qJD(2) + (-t316 * t344 - t346) * qJD(1)) * pkin(2), t337, t308 * r_i_i_C(2) + (r_i_i_C(1) * t325 - t286 * r_i_i_C(2)) * t320 + (-t353 * r_i_i_C(1) - t325 * r_i_i_C(2)) * t317, 0, 0; 0, -t351 * t292 + (t304 * t326 - t324 * t327 - t338) * t314, 0, t332 * t292 + ((-t298 * t320 - t316 * t317) * r_i_i_C(1) + (t298 * t317 - t316 * t320) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:35
	% EndTime: 2019-10-10 10:13:35
	% DurationCPUTime: 0.52s
	% Computational Cost: add. (538->86), mult. (1584->146), div. (0->0), fcn. (1696->10), ass. (0->64)
	t353 = sin(pkin(11));
	t360 = cos(qJ(2));
	t393 = cos(pkin(11));
	t376 = qJD(2) * t393;
	t357 = sin(qJ(2));
	t382 = qJD(2) * t357;
	t400 = t353 * t382 - t360 * t376;
	t355 = cos(pkin(6));
	t334 = t400 * t355;
	t345 = -t360 * t353 - t357 * t393;
	t340 = t345 * t355;
	t361 = cos(qJ(1));
	t368 = -t357 * t353 + t360 * t393;
	t358 = sin(qJ(1));
	t384 = qJD(1) * t358;
	t381 = qJD(2) * t360;
	t342 = -t353 * t381 - t357 * t376;
	t388 = t358 * t342;
	t323 = t388 + t340 * t384 + (qJD(1) * t368 - t334) * t361;
	t356 = sin(qJ(4));
	t359 = cos(qJ(4));
	t328 = -t361 * t340 + t358 * t368;
	t354 = sin(pkin(6));
	t391 = t354 * t361;
	t371 = t328 * t356 + t359 * t391;
	t379 = t354 * t384;
	t399 = t371 * qJD(4) - t323 * t359 - t356 * t379;
	t370 = -t328 * t359 + t356 * t391;
	t398 = t370 * qJD(4) - t323 * t356 + t359 * t379;
	t394 = r_i_i_C(3) + qJ(5);
	t397 = pkin(4) - r_i_i_C(2);
	t365 = t394 * t356 + t397 * t359 + pkin(3);
	t396 = r_i_i_C(1) + pkin(9);
	t395 = pkin(2) * t355;
	t372 = -t358 * t340 - t361 * t368;
	t392 = t372 * t356;
	t390 = t357 * t358;
	t389 = t357 * t361;
	t387 = t358 * t359;
	t386 = t358 * t360;
	t385 = t360 * t361;
	t383 = qJD(1) * t361;
	t380 = pkin(2) * t382;
	t378 = t354 * t383;
	t338 = t345 * t354;
	t374 = -t338 * t359 + t355 * t356;
	t339 = t368 * t355;
	t373 = t361 * t339 + t358 * t345;
	t369 = t358 * t354 * t356 - t359 * t372;
	t364 = qJD(2) * t345;
	t363 = t368 * qJD(2);
	t362 = qJD(5) * t356 + (-t397 * t356 + t394 * t359) * qJD(4);
	t320 = -t328 * qJD(1) + t358 * t334 + t361 * t342;
	t352 = t360 * pkin(2) + pkin(1);
	t343 = -t354 * qJD(3) + t381 * t395;
	t341 = t357 * t395 + (-pkin(8) - qJ(3)) * t354;
	t335 = t355 * t364;
	t332 = t400 * t354;
	t324 = t374 * qJD(4) - t332 * t356;
	t322 = t361 * t335 - t339 * t384 + t345 * t383 - t358 * t363;
	t319 = t373 * qJD(1) + t358 * t335 + t361 * t363;
	t313 = t320 * t359 + qJD(4) * t392 + (qJD(4) * t387 + t356 * t383) * t354;
	t312 = t369 * qJD(4) + t320 * t356 - t359 * t378;
	t1 = [-t371 * qJD(5) - t323 * pkin(3) + t358 * t380 - t361 * t343 + t396 * t322 + t397 * t399 + t394 * t398 + (t358 * t341 - t361 * t352) * qJD(1), t396 * t320 + ((t355 * t390 - t385) * qJD(2) + (-t355 * t385 + t390) * qJD(1)) * pkin(2) + t362 * (-t358 * t339 + t361 * t345) - t365 * t319, t378, t369 * qJD(5) - t397 * t312 + t394 * t313, t312, 0; -(t354 * t387 + t392) * qJD(5) + t320 * pkin(3) - t361 * t380 - t358 * t343 + t396 * t319 + t397 * t313 + t394 * t312 + (-t361 * t341 - t358 * t352) * qJD(1), -t396 * (t372 * qJD(1) + t361 * t334 - t388) + ((-t355 * t389 - t386) * qJD(2) + (-t355 * t386 - t389) * qJD(1)) * pkin(2) + t362 * t373 + t365 * t322, t379, -t370 * qJD(5) - t394 * t399 + t397 * t398, -t398, 0; 0, -t396 * t332 + (t362 * t368 + t364 * t365 - t380) * t354, 0, t374 * qJD(5) + t394 * (-t332 * t359 + (t338 * t356 + t355 * t359) * qJD(4)) - t397 * t324, t324, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:36
	% EndTime: 2019-10-10 10:13:37
	% DurationCPUTime: 1.04s
	% Computational Cost: add. (1031->116), mult. (3020->197), div. (0->0), fcn. (3324->12), ass. (0->80)
	t453 = sin(qJ(1));
	t456 = cos(qJ(1));
	t449 = cos(pkin(6));
	t447 = sin(pkin(11));
	t452 = sin(qJ(2));
	t501 = cos(pkin(11));
	t503 = cos(qJ(2));
	t476 = t503 * t501;
	t463 = -t452 * t447 + t476;
	t461 = t449 * t463;
	t478 = t452 * t501;
	t464 = t503 * t447 + t478;
	t412 = -t453 * t464 + t456 * t461;
	t450 = sin(qJ(6));
	t454 = cos(qJ(6));
	t430 = t464 * t449;
	t413 = t430 * t456 + t453 * t463;
	t451 = sin(qJ(4));
	t455 = cos(qJ(4));
	t448 = sin(pkin(6));
	t495 = t448 * t456;
	t511 = t413 * t451 + t455 * t495;
	t515 = -t412 * t454 + t450 * t511;
	t514 = -t412 * t450 - t454 * t511;
	t475 = t454 * r_i_i_C(1) - t450 * r_i_i_C(2);
	t462 = t475 * qJD(6) + qJD(5);
	t474 = -r_i_i_C(1) * t450 - r_i_i_C(2) * t454;
	t472 = qJ(5) - t474;
	t487 = r_i_i_C(3) + pkin(10) + pkin(4);
	t513 = (t487 * t451 - t472 * t455) * qJD(4) - t462 * t451;
	t510 = qJD(1) * t461 + qJD(2) * t463;
	t489 = qJD(2) * t452;
	t509 = -qJD(2) * t476 + t447 * t489;
	t425 = t509 * t449;
	t479 = t503 * qJD(2);
	t432 = -qJD(2) * t478 - t447 * t479;
	t491 = qJD(1) * t453;
	t508 = t430 * t491 - t453 * t432 + (-qJD(1) * t463 + t425) * t456;
	t482 = t448 * t491;
	t507 = qJD(4) * t511 - t451 * t482 + t455 * t508;
	t505 = t472 * t451 + t487 * t455 + pkin(3);
	t504 = pkin(5) + pkin(9);
	t502 = pkin(2) * t449;
	t416 = -t453 * t430 + t456 * t463;
	t497 = t416 * t451;
	t496 = t448 * t453;
	t493 = t452 * t453;
	t492 = t452 * t456;
	t490 = qJD(1) * t456;
	t488 = qJD(4) * t455;
	t486 = pkin(2) * t489;
	t485 = t451 * t495;
	t484 = t503 * t453;
	t483 = t503 * t456;
	t481 = t448 * t490;
	t429 = t464 * t448;
	t473 = t429 * t455 + t449 * t451;
	t417 = t429 * t451 - t449 * t455;
	t470 = t475 + t504;
	t468 = t416 * t455 + t451 * t496;
	t467 = qJD(6) * t474;
	t393 = -qJD(4) * t485 + t413 * t488 - t451 * t508 - t455 * t482;
	t399 = -t413 * qJD(1) + t453 * t425 + t432 * t456;
	t446 = t503 * pkin(2) + pkin(1);
	t433 = -t448 * qJD(3) + t479 * t502;
	t431 = t452 * t502 + (-pkin(8) - qJ(3)) * t448;
	t428 = t463 * t448;
	t426 = qJD(2) * t430;
	t424 = qJD(2) * t429;
	t423 = t509 * t448;
	t415 = -t453 * t461 - t456 * t464;
	t408 = -t455 * t496 + t497;
	t403 = t473 * qJD(4) - t423 * t451;
	t401 = -t456 * t426 - t510 * t453 - t464 * t490;
	t398 = -t453 * t426 + t510 * t456 - t464 * t491;
	t392 = t399 * t455 - qJD(4) * t497 + (t451 * t490 + t453 * t488) * t448;
	t391 = t468 * qJD(4) + t399 * t451 - t455 * t481;
	t390 = t391 * t450 + t398 * t454 + (t408 * t454 + t415 * t450) * qJD(6);
	t389 = t391 * t454 - t398 * t450 + (-t408 * t450 + t415 * t454) * qJD(6);
	t1 = [t453 * t486 + t508 * pkin(3) - t511 * qJD(5) - t456 * t433 - t472 * t393 + t470 * t401 + (t514 * r_i_i_C(1) + t515 * r_i_i_C(2)) * qJD(6) + (t431 * t453 - t446 * t456) * qJD(1) + t487 * t507, t416 * t467 + t470 * t399 + ((t449 * t493 - t483) * qJD(2) + (-t449 * t483 + t493) * qJD(1)) * pkin(2) - t505 * t398 - t513 * t415, t481, -t487 * t391 + t472 * t392 + t462 * t468, t391, r_i_i_C(1) * t389 - r_i_i_C(2) * t390; -t456 * t486 + t399 * pkin(3) + t390 * r_i_i_C(1) + t389 * r_i_i_C(2) + t391 * qJ(5) + t408 * qJD(5) - t453 * t433 + t504 * t398 + (-t431 * t456 - t446 * t453) * qJD(1) + t487 * t392, t413 * t467 - t470 * t508 + ((-t449 * t492 - t484) * qJD(2) + (-t449 * t484 - t492) * qJD(1)) * pkin(2) + t505 * t401 - t513 * t412, t482, t462 * (t413 * t455 - t485) - t472 * t507 - t487 * t393, t393, (t393 * t454 + t401 * t450) * r_i_i_C(1) + (-t393 * t450 + t401 * t454) * r_i_i_C(2) + (-t515 * r_i_i_C(1) + t514 * r_i_i_C(2)) * qJD(6); 0, -t470 * t423 - t424 * t505 - t513 * t428 + t429 * t467 - t448 * t486, 0, t462 * t473 + t472 * (-t417 * qJD(4) - t423 * t455) - t487 * t403, t403, (t403 * t454 - t424 * t450) * r_i_i_C(1) + (-t403 * t450 - t424 * t454) * r_i_i_C(2) + ((-t417 * t450 + t428 * t454) * r_i_i_C(1) + (-t417 * t454 - t428 * t450) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end