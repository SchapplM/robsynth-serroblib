% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR4
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
% Datum: 2019-10-10 10:55
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
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
	% StartTime: 2019-10-10 10:55:49
	% EndTime: 2019-10-10 10:55:49
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
	% StartTime: 2019-10-10 10:55:49
	% EndTime: 2019-10-10 10:55:49
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (90->45), mult. (273->84), div. (0->0), fcn. (260->8), ass. (0->34)
	t190 = cos(pkin(6));
	t187 = sin(pkin(12));
	t189 = cos(pkin(12));
	t191 = sin(qJ(2));
	t193 = cos(qJ(2));
	t196 = t187 * t193 + t189 * t191;
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
	t180 = t187 * t191 - t193 * t189;
	t199 = -pkin(2) * t193 + r_i_i_C(1) * t180 - pkin(1);
	t172 = t187 * t190 * t203 - t189 * t200;
	t178 = -t187 * t202 - t189 * t203;
	t198 = t194 * t172 - t192 * t178;
	t197 = t192 * t172 + t194 * t178;
	t195 = t175 * r_i_i_C(1) + t190 * t209 + (-r_i_i_C(3) - pkin(8) - qJ(3)) * t188;
	t179 = pkin(2) * t200 - qJD(3) * t188;
	t177 = t180 * qJD(2);
	t174 = t180 * t190;
	t173 = qJD(2) * t175;
	t171 = -t194 * t173 + t192 * t177 + (t174 * t192 - t194 * t196) * qJD(1);
	t170 = t192 * t173 + t194 * t177 + (t174 * t194 + t192 * t196) * qJD(1);
	t1 = [t198 * r_i_i_C(1) - t171 * r_i_i_C(2) + t192 * t201 - t194 * t179 + (t192 * t195 + t194 * t199) * qJD(1), t170 * r_i_i_C(1) + ((t175 * t194 - t180 * t192) * qJD(1) - t197) * r_i_i_C(2) + ((t190 * t208 - t205) * qJD(2) + (-t190 * t205 + t208) * qJD(1)) * pkin(2), t194 * t204, 0, 0, 0; t197 * r_i_i_C(1) + t170 * r_i_i_C(2) - t194 * t201 - t192 * t179 + (t192 * t199 - t194 * t195) * qJD(1), t171 * r_i_i_C(1) + ((t175 * t192 + t180 * t194) * qJD(1) + t198) * r_i_i_C(2) + ((-t190 * t207 - t206) * qJD(2) + (-t190 * t206 - t207) * qJD(1)) * pkin(2), t192 * t204, 0, 0, 0; 0, (-r_i_i_C(1) * t196 + r_i_i_C(2) * t180 - t209) * t188 * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:50
	% EndTime: 2019-10-10 10:55:50
	% DurationCPUTime: 0.44s
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
	% StartTime: 2019-10-10 10:55:50
	% EndTime: 2019-10-10 10:55:50
	% DurationCPUTime: 0.60s
	% Computational Cost: add. (535->94), mult. (1256->158), div. (0->0), fcn. (1310->12), ass. (0->71)
	t350 = qJ(4) + qJ(5);
	t347 = sin(t350);
	t354 = cos(pkin(6));
	t351 = sin(pkin(12));
	t353 = cos(pkin(12));
	t359 = cos(qJ(2));
	t382 = qJD(2) * t359;
	t356 = sin(qJ(2));
	t383 = qJD(2) * t356;
	t399 = t351 * t383 - t353 * t382;
	t327 = t399 * t354;
	t368 = t359 * t351 + t356 * t353;
	t333 = t368 * t354;
	t338 = t356 * t351 - t359 * t353;
	t360 = cos(qJ(1));
	t357 = sin(qJ(1));
	t384 = qJD(1) * t357;
	t335 = -t351 * t382 - t353 * t383;
	t390 = t357 * t335;
	t317 = t390 - t333 * t384 + (-qJD(1) * t338 - t327) * t360;
	t349 = qJD(4) + qJD(5);
	t352 = sin(pkin(6));
	t393 = t352 * t360;
	t400 = -t349 * t393 + t317;
	t401 = t400 * t347;
	t398 = pkin(2) * t354;
	t397 = r_i_i_C(3) + pkin(10) + pkin(9);
	t396 = t347 * t349;
	t348 = cos(t350);
	t395 = t348 * t349;
	t394 = t352 * t357;
	t392 = t356 * t357;
	t391 = t356 * t360;
	t389 = t357 * t359;
	t388 = t359 * t360;
	t369 = t357 * t333 + t360 * t338;
	t378 = qJD(1) * t393;
	t365 = t349 * t369 + t378;
	t321 = t360 * t333 - t357 * t338;
	t314 = -t321 * qJD(1) + t357 * t327 + t360 * t335;
	t372 = t349 * t394 + t314;
	t310 = -t372 * t347 + t365 * t348;
	t311 = t365 * t347 + t372 * t348;
	t387 = t310 * r_i_i_C(1) - t311 * r_i_i_C(2);
	t379 = t352 * t384;
	t366 = -t321 * t349 + t379;
	t375 = t400 * t348;
	t386 = (t366 * t348 - t401) * r_i_i_C(1) + (-t366 * t347 - t375) * r_i_i_C(2);
	t331 = t368 * t352;
	t325 = t399 * t352;
	t374 = -t349 * t354 + t325;
	t385 = (-t331 * t395 + t374 * t347) * r_i_i_C(1) + (t331 * t396 + t374 * t348) * r_i_i_C(2);
	t381 = pkin(2) * t383;
	t373 = -r_i_i_C(1) * t347 - r_i_i_C(2) * t348;
	t332 = t338 * t354;
	t371 = -t360 * t332 - t357 * t368;
	t370 = t357 * t332 - t360 * t368;
	t358 = cos(qJ(4));
	t345 = t358 * pkin(4) + pkin(3);
	t367 = t348 * r_i_i_C(1) - t347 * r_i_i_C(2) + t345;
	t364 = qJD(2) * t368;
	t363 = t338 * qJD(2);
	t355 = sin(qJ(4));
	t362 = -qJD(4) * t355 * pkin(4) + t373 * t349;
	t346 = t359 * pkin(2) + pkin(1);
	t336 = -t352 * qJD(3) + t382 * t398;
	t334 = t356 * t398 + (-pkin(8) - qJ(3)) * t352;
	t328 = t354 * t364;
	t316 = t370 * qJD(1) - t360 * t328 + t357 * t363;
	t313 = t371 * qJD(1) - t357 * t328 - t360 * t363;
	t1 = [(t321 * t396 - t375) * r_i_i_C(1) + (t321 * t395 + t401) * r_i_i_C(2) - t317 * t345 + t357 * t381 - t360 * t336 + t397 * t316 + (-t360 * t346 + (t373 * t352 + t334) * t357) * qJD(1) + (-t355 * t379 + (t321 * t355 + t358 * t393) * qJD(4)) * pkin(4), t397 * t314 + t362 * t370 - t367 * t313 + ((t354 * t392 - t388) * qJD(2) + (-t354 * t388 + t392) * qJD(1)) * pkin(2), t378, (t358 * t378 - t314 * t355 + (-t355 * t394 + t358 * t369) * qJD(4)) * pkin(4) + t387, t387, 0; -t360 * t381 + t311 * r_i_i_C(1) + t310 * r_i_i_C(2) + t314 * t345 - t357 * t336 + t397 * t313 + (-t334 * t360 - t346 * t357) * qJD(1) + (t355 * t378 + (t355 * t369 + t358 * t394) * qJD(4)) * pkin(4), -t397 * (t369 * qJD(1) + t360 * t327 - t390) + t362 * t371 + t367 * t316 + ((-t354 * t391 - t389) * qJD(2) + (-t354 * t389 - t391) * qJD(1)) * pkin(2), t379, (t358 * t379 - t317 * t355 + (-t321 * t358 + t355 * t393) * qJD(4)) * pkin(4) + t386, t386, 0; 0, -t397 * t325 + (-t338 * t362 - t364 * t367 - t381) * t352, 0, (t325 * t355 + (-t331 * t358 - t354 * t355) * qJD(4)) * pkin(4) + t385, t385, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:52
	% EndTime: 2019-10-10 10:55:54
	% DurationCPUTime: 1.45s
	% Computational Cost: add. (1457->150), mult. (3222->250), div. (0->0), fcn. (3524->14), ass. (0->91)
	t579 = pkin(11) + r_i_i_C(3);
	t520 = cos(qJ(6));
	t558 = qJD(6) * t520;
	t516 = sin(qJ(6));
	t559 = qJD(6) * t516;
	t589 = -r_i_i_C(1) * t559 - t558 * r_i_i_C(2);
	t512 = qJ(4) + qJ(5);
	t509 = sin(t512);
	t510 = cos(t512);
	t515 = cos(pkin(6));
	t513 = sin(pkin(12));
	t518 = sin(qJ(2));
	t575 = cos(pkin(12));
	t547 = t518 * t575;
	t578 = cos(qJ(2));
	t534 = t578 * t513 + t547;
	t493 = t534 * t515;
	t519 = sin(qJ(1));
	t522 = cos(qJ(1));
	t543 = t578 * t575;
	t533 = -t518 * t513 + t543;
	t540 = t522 * t493 + t519 * t533;
	t514 = sin(pkin(6));
	t567 = t514 * t522;
	t468 = t509 * t567 - t510 * t540;
	t532 = t533 * t515;
	t474 = -t519 * t534 + t522 * t532;
	t588 = -t468 * t516 + t474 * t520;
	t587 = t468 * t520 + t474 * t516;
	t511 = qJD(4) + qJD(5);
	t517 = sin(qJ(4));
	t585 = -t520 * r_i_i_C(1) + t516 * r_i_i_C(2) - pkin(5);
	t586 = (-t509 * t585 - t579 * t510) * t511 + (t516 * r_i_i_C(1) + t520 * r_i_i_C(2)) * t510 * qJD(6) + qJD(4) * t517 * pkin(4);
	t584 = qJD(1) * t532 + t533 * qJD(2);
	t561 = qJD(2) * t518;
	t583 = -qJD(2) * t543 + t513 * t561;
	t488 = t583 * t515;
	t548 = t578 * qJD(2);
	t495 = -qJD(2) * t547 - t513 * t548;
	t563 = qJD(1) * t519;
	t458 = t493 * t563 - t519 * t495 + (-qJD(1) * t533 + t488) * t522;
	t521 = cos(qJ(4));
	t507 = t521 * pkin(4) + pkin(3);
	t580 = t579 * t509 - t510 * t585 + t507;
	t577 = pkin(2) * t515;
	t489 = qJD(2) * t493;
	t562 = qJD(1) * t522;
	t459 = -t522 * t489 - t584 * t519 - t534 * t562;
	t574 = t459 * t516;
	t573 = t459 * t520;
	t570 = t509 * t511;
	t569 = t510 * t511;
	t568 = t514 * t519;
	t565 = t518 * t519;
	t564 = t518 * t522;
	t557 = pkin(2) * t561;
	t554 = t510 * t567;
	t553 = t578 * t519;
	t552 = t578 * t522;
	t551 = t514 * t563;
	t550 = t514 * t562;
	t546 = t458 * t510 + t511 * t554;
	t486 = t583 * t514;
	t545 = t511 * t515 - t486;
	t528 = -t540 * qJD(1) + t519 * t488 + t522 * t495;
	t541 = t511 * t568 + t528;
	t539 = -t519 * t493 + t522 * t533;
	t535 = -t511 * t540 + t551;
	t492 = t534 * t514;
	t451 = t535 * t510 + (t511 * t567 + t458) * t509;
	t452 = t509 * t551 - t540 * t570 - t546;
	t527 = t589 * (-t509 * t540 - t554) + t579 * t452 - t585 * t451;
	t449 = t541 * t509 - t510 * t550 + t539 * t569;
	t450 = t509 * t550 + t541 * t510 - t539 * t570;
	t526 = t589 * (-t509 * t539 + t510 * t568) + t579 * t450 + t585 * t449;
	t465 = -t492 * t570 + t545 * t510;
	t525 = t589 * (-t492 * t509 + t515 * t510) + t579 * t465 - t585 * (-t492 * t569 - t545 * t509);
	t523 = -pkin(10) - pkin(9);
	t508 = t578 * pkin(2) + pkin(1);
	t496 = -t514 * qJD(3) + t548 * t577;
	t494 = t518 * t577 + (-pkin(8) - qJ(3)) * t514;
	t491 = t533 * t514;
	t487 = qJD(2) * t492;
	t480 = t492 * t510 + t515 * t509;
	t477 = -t519 * t532 - t522 * t534;
	t470 = t509 * t568 + t510 * t539;
	t456 = -t519 * t489 + t584 * t522 - t534 * t563;
	t454 = -t535 * t509 + t546;
	t442 = t450 * t520 + t456 * t516 + (-t470 * t516 - t477 * t520) * qJD(6);
	t441 = -t450 * t516 + t456 * t520 + (-t470 * t520 + t477 * t516) * qJD(6);
	t1 = [(t454 * t520 + t574) * r_i_i_C(1) + (-t454 * t516 + t573) * r_i_i_C(2) + t454 * pkin(5) + t458 * t507 - t459 * t523 + t519 * t557 - t522 * t496 + t579 * t451 + (t588 * r_i_i_C(1) - t587 * r_i_i_C(2)) * qJD(6) + (t519 * t494 - t522 * t508) * qJD(1) + (-t517 * t551 + (t517 * t540 + t521 * t567) * qJD(4)) * pkin(4), (t516 * t528 + t539 * t558) * r_i_i_C(1) + (t520 * t528 - t539 * t559) * r_i_i_C(2) - t528 * t523 - t580 * t456 + ((t515 * t565 - t552) * qJD(2) + (-t515 * t552 + t565) * qJD(1)) * pkin(2) - t586 * t477, t550, (t521 * t550 - t528 * t517 + (-t517 * t568 - t521 * t539) * qJD(4)) * pkin(4) + t526, t526, t441 * r_i_i_C(1) - t442 * r_i_i_C(2); -t522 * t557 + t450 * pkin(5) + t442 * r_i_i_C(1) + t441 * r_i_i_C(2) - t456 * t523 + t528 * t507 - t519 * t496 + t579 * t449 + (-t494 * t522 - t508 * t519) * qJD(1) + (t517 * t550 + (-t517 * t539 + t521 * t568) * qJD(4)) * pkin(4), (-t458 * t516 + t540 * t558) * r_i_i_C(1) + (-t458 * t520 - t540 * t559) * r_i_i_C(2) + t458 * t523 + t580 * t459 + ((-t515 * t564 - t553) * qJD(2) + (-t515 * t553 - t564) * qJD(1)) * pkin(2) - t586 * t474, t551, (t521 * t551 + t458 * t517 + (t517 * t567 - t521 * t540) * qJD(4)) * pkin(4) + t527, t527, (-t452 * t516 - t573) * r_i_i_C(1) + (-t452 * t520 + t574) * r_i_i_C(2) + (t587 * r_i_i_C(1) + t588 * r_i_i_C(2)) * qJD(6); 0, (-t486 * t516 + t492 * t558) * r_i_i_C(1) + (-t486 * t520 - t492 * t559) * r_i_i_C(2) + t486 * t523 - t514 * t557 - t580 * t487 - t586 * t491, 0, (t486 * t517 + (-t492 * t521 - t515 * t517) * qJD(4)) * pkin(4) + t525, t525, (-t465 * t516 + t487 * t520) * r_i_i_C(1) + (-t465 * t520 - t487 * t516) * r_i_i_C(2) + ((-t480 * t520 + t491 * t516) * r_i_i_C(1) + (t480 * t516 + t491 * t520) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end