% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR10
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
% Datum: 2019-10-10 11:05
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:05:08
	% EndTime: 2019-10-10 11:05:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:05:08
	% EndTime: 2019-10-10 11:05:08
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
	% StartTime: 2019-10-10 11:05:08
	% EndTime: 2019-10-10 11:05:08
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
	% StartTime: 2019-10-10 11:05:09
	% EndTime: 2019-10-10 11:05:09
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (103->31), mult. (307->49), div. (0->0), fcn. (282->8), ass. (0->26)
	t199 = sin(pkin(12));
	t200 = sin(pkin(6));
	t201 = cos(pkin(12));
	t220 = t200 * (r_i_i_C(1) * t199 + r_i_i_C(2) * t201 + pkin(8));
	t219 = r_i_i_C(3) + qJ(3);
	t203 = sin(qJ(2));
	t204 = sin(qJ(1));
	t218 = t204 * t203;
	t205 = cos(qJ(2));
	t217 = t204 * t205;
	t206 = cos(qJ(1));
	t216 = t206 * t203;
	t215 = t206 * t205;
	t214 = qJD(2) * t203;
	t202 = cos(pkin(6));
	t213 = t202 * t218;
	t212 = t202 * t215;
	t211 = qJD(2) * t202 + qJD(1);
	t210 = t201 * r_i_i_C(1) - t199 * r_i_i_C(2) + pkin(2);
	t208 = t202 * t217 + t216;
	t207 = t202 * t216 + t217;
	t194 = -qJD(1) * t213 - t204 * t214 + t211 * t215;
	t193 = t208 * qJD(1) + t207 * qJD(2);
	t192 = t207 * qJD(1) + t208 * qJD(2);
	t191 = -qJD(1) * t212 - qJD(2) * t215 + t211 * t218;
	t1 = [-(-t212 + t218) * qJD(3) - t219 * t193 - t210 * t194 + (-t206 * pkin(1) - t204 * t220) * qJD(1), -(t213 - t215) * qJD(3) - t219 * t192 + t210 * t191, -t191, 0, 0, 0; t208 * qJD(3) - t219 * t191 - t210 * t192 + (-t204 * pkin(1) + t206 * t220) * qJD(1), t207 * qJD(3) - t210 * t193 + t219 * t194, t193, 0, 0, 0; 0, (t203 * qJD(3) + (-t210 * t203 + t219 * t205) * qJD(2)) * t200, t200 * t214, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:05:09
	% EndTime: 2019-10-10 11:05:09
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (222->61), mult. (494->100), div. (0->0), fcn. (467->10), ass. (0->45)
	t259 = sin(qJ(1));
	t256 = cos(pkin(6));
	t267 = qJD(2) * t256 + qJD(1);
	t258 = sin(qJ(2));
	t280 = t259 * t258;
	t272 = t256 * t280;
	t275 = qJD(2) * t258;
	t260 = cos(qJ(2));
	t261 = cos(qJ(1));
	t277 = t261 * t260;
	t239 = -qJD(1) * t272 - t259 * t275 + t267 * t277;
	t255 = sin(pkin(6));
	t281 = t255 * t261;
	t268 = qJD(4) * t281;
	t285 = t239 - t268;
	t253 = pkin(12) + qJ(4);
	t251 = sin(t253);
	t252 = cos(t253);
	t266 = r_i_i_C(1) * t251 + r_i_i_C(2) * t252;
	t263 = qJD(4) * t266;
	t284 = r_i_i_C(3) + pkin(9) + qJ(3);
	t283 = t255 * t258;
	t282 = t255 * t259;
	t279 = t259 * t260;
	t278 = t261 * t258;
	t276 = qJD(1) * t255;
	t274 = qJD(2) * t260;
	t241 = t256 * t278 + t279;
	t273 = qJD(4) * t241;
	t271 = t256 * t277;
	t270 = pkin(3) * sin(pkin(12)) + pkin(8);
	t269 = t261 * t276;
	t250 = cos(pkin(12)) * pkin(3) + pkin(2);
	t265 = t252 * r_i_i_C(1) - t251 * r_i_i_C(2) + t250;
	t264 = t256 * t279 + t278;
	t262 = t259 * t276 - t273;
	t244 = t252 * t268;
	t243 = -t272 + t277;
	t240 = -t271 + t280;
	t238 = t264 * qJD(1) + t241 * qJD(2);
	t237 = t241 * qJD(1) + t264 * qJD(2);
	t236 = -qJD(1) * t271 - t261 * t274 + t267 * t280;
	t235 = t251 * t269 - t237 * t252 + (-t243 * t251 + t252 * t282) * qJD(4);
	t234 = t252 * t269 + t237 * t251 + (-t243 * t252 - t251 * t282) * qJD(4);
	t1 = [(-t239 * t252 + t251 * t273 + t244) * r_i_i_C(1) + (t285 * t251 + t252 * t273) * r_i_i_C(2) - t239 * t250 - t240 * qJD(3) - t284 * t238 + (-t261 * pkin(1) + (-t266 - t270) * t282) * qJD(1), t243 * qJD(3) + t265 * t236 - t284 * t237 + t263 * t264, -t236, t234 * r_i_i_C(1) - t235 * r_i_i_C(2), 0, 0; t235 * r_i_i_C(1) + t234 * r_i_i_C(2) + t264 * qJD(3) - t237 * t250 - t284 * t236 + (-pkin(1) * t259 + t270 * t281) * qJD(1), t241 * qJD(3) - t265 * t238 + t284 * t239 + t240 * t263, t238, t244 * r_i_i_C(2) + (t262 * r_i_i_C(1) - t239 * r_i_i_C(2)) * t252 + (-t285 * r_i_i_C(1) - t262 * r_i_i_C(2)) * t251, 0, 0; 0, (qJD(3) * t258 - t260 * t263 + (-t265 * t258 + t284 * t260) * qJD(2)) * t255, t255 * t275, -t266 * t255 * t274 + ((-t251 * t256 - t252 * t283) * r_i_i_C(1) + (t251 * t283 - t252 * t256) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:05:10
	% EndTime: 2019-10-10 11:05:11
	% DurationCPUTime: 0.80s
	% Computational Cost: add. (656->106), mult. (1384->172), div. (0->0), fcn. (1389->12), ass. (0->69)
	t392 = sin(qJ(1));
	t388 = cos(pkin(6));
	t407 = qJD(2) * t388 + qJD(1);
	t391 = sin(qJ(2));
	t426 = t392 * t391;
	t416 = t388 * t426;
	t421 = qJD(2) * t391;
	t394 = cos(qJ(2));
	t395 = cos(qJ(1));
	t423 = t395 * t394;
	t361 = -qJD(1) * t416 - t392 * t421 + t407 * t423;
	t424 = t395 * t391;
	t425 = t392 * t394;
	t372 = t388 * t424 + t425;
	t385 = pkin(12) + qJ(4);
	t383 = sin(t385);
	t384 = cos(t385);
	t387 = sin(pkin(6));
	t422 = qJD(1) * t387;
	t412 = t392 * t422;
	t427 = t387 * t395;
	t415 = t384 * t427;
	t355 = (-qJD(4) * t372 + t412) * t383 - qJD(4) * t415 + t361 * t384;
	t373 = t388 * t425 + t424;
	t360 = t373 * qJD(1) + t372 * qJD(2);
	t390 = sin(qJ(5));
	t393 = cos(qJ(5));
	t444 = t355 * t390 - t360 * t393;
	t443 = -t355 * t393 - t360 * t390;
	t405 = t393 * r_i_i_C(1) - t390 * r_i_i_C(2);
	t403 = pkin(4) + t405;
	t437 = r_i_i_C(3) + pkin(10);
	t442 = (t403 * t383 - t437 * t384) * qJD(4);
	t366 = -t372 * t384 + t383 * t427;
	t414 = t388 * t423;
	t371 = -t414 + t426;
	t441 = -t366 * t390 - t371 * t393;
	t440 = t366 * t393 - t371 * t390;
	t404 = t390 * r_i_i_C(1) + t393 * r_i_i_C(2);
	t382 = cos(pkin(12)) * pkin(3) + pkin(2);
	t438 = t437 * t383 + t403 * t384 + t382;
	t430 = t387 * t391;
	t429 = t387 * t392;
	t428 = t387 * t394;
	t420 = qJD(2) * t394;
	t419 = qJD(5) * t384;
	t418 = qJD(5) * t390;
	t417 = qJD(5) * t393;
	t413 = pkin(3) * sin(pkin(12)) + pkin(8);
	t411 = t395 * t422;
	t410 = t387 * t420;
	t409 = t387 * t421;
	t374 = -t416 + t423;
	t402 = -t374 * t383 + t384 * t429;
	t368 = t374 * t384 + t383 * t429;
	t370 = t388 * t383 + t384 * t430;
	t401 = -t383 * t430 + t388 * t384;
	t400 = qJD(5) * t404;
	t397 = t366 * qJD(4) - t361 * t383 + t384 * t412;
	t396 = t404 * t419 + t442;
	t389 = -pkin(9) - qJ(3);
	t363 = t401 * qJD(4) + t384 * t410;
	t359 = t372 * qJD(1) + t373 * qJD(2);
	t358 = -qJD(1) * t414 - t395 * t420 + t407 * t426;
	t353 = t402 * qJD(4) - t359 * t384 + t383 * t411;
	t352 = t368 * qJD(4) - t359 * t383 - t384 * t411;
	t351 = t353 * t393 - t358 * t390 + (-t368 * t390 + t373 * t393) * qJD(5);
	t350 = -t353 * t390 - t358 * t393 + (-t368 * t393 - t373 * t390) * qJD(5);
	t1 = [t443 * r_i_i_C(1) + t444 * r_i_i_C(2) - t355 * pkin(4) - t361 * t382 + t360 * t389 - t371 * qJD(3) + t437 * t397 + (t441 * r_i_i_C(1) - t440 * r_i_i_C(2)) * qJD(5) + (-t395 * pkin(1) - t413 * t429) * qJD(1), (-t359 * t390 + t374 * t417) * r_i_i_C(1) + (-t359 * t393 - t374 * t418) * r_i_i_C(2) + t359 * t389 + t374 * qJD(3) + t438 * t358 + t396 * t373, -t358, -t403 * t352 + t437 * t353 - t402 * t400, t350 * r_i_i_C(1) - t351 * r_i_i_C(2), 0; t353 * pkin(4) + t351 * r_i_i_C(1) + t350 * r_i_i_C(2) + t373 * qJD(3) + t358 * t389 - t359 * t382 + t437 * t352 + (-pkin(1) * t392 + t413 * t427) * qJD(1), (t361 * t390 + t372 * t417) * r_i_i_C(1) + (t361 * t393 - t372 * t418) * r_i_i_C(2) - t361 * t389 + t372 * qJD(3) - t438 * t360 + t396 * t371, t360, t437 * t355 - (-t372 * t383 - t415) * t400 + t403 * t397, -t444 * r_i_i_C(1) + t443 * r_i_i_C(2) + (t440 * r_i_i_C(1) + t441 * r_i_i_C(2)) * qJD(5), 0; 0, ((-qJD(2) * t438 + t405 * qJD(5) + qJD(3)) * t391 + (-qJD(2) * t389 - t442 + t404 * (qJD(2) - t419)) * t394) * t387, t409, t437 * t363 - t401 * t400 + t403 * (-t370 * qJD(4) - t383 * t410), (-t363 * t390 + t393 * t409) * r_i_i_C(1) + (-t363 * t393 - t390 * t409) * r_i_i_C(2) + ((-t370 * t393 + t390 * t428) * r_i_i_C(1) + (t370 * t390 + t393 * t428) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:05:11
	% EndTime: 2019-10-10 11:05:12
	% DurationCPUTime: 1.10s
	% Computational Cost: add. (1063->129), mult. (1896->202), div. (0->0), fcn. (1907->14), ass. (0->80)
	t431 = sin(qJ(1));
	t484 = cos(pkin(6));
	t445 = qJD(2) * t484 + qJD(1);
	t430 = sin(qJ(2));
	t460 = t431 * t484;
	t450 = t430 * t460;
	t470 = qJD(2) * t430;
	t433 = cos(qJ(2));
	t434 = cos(qJ(1));
	t475 = t434 * t433;
	t396 = -qJD(1) * t450 - t431 * t470 + t445 * t475;
	t459 = t434 * t484;
	t407 = t430 * t459 + t431 * t433;
	t423 = pkin(12) + qJ(4);
	t419 = sin(t423);
	t420 = cos(t423);
	t427 = sin(pkin(6));
	t471 = qJD(1) * t427;
	t464 = t431 * t471;
	t477 = t427 * t434;
	t466 = t420 * t477;
	t390 = (-qJD(4) * t407 + t464) * t419 - qJD(4) * t466 + t396 * t420;
	t449 = t433 * t459;
	t476 = t431 * t430;
	t406 = -t449 + t476;
	t424 = qJD(5) + qJD(6);
	t456 = t406 * t424 + t390;
	t432 = cos(qJ(5));
	t418 = t432 * pkin(5) + pkin(4);
	t425 = qJ(5) + qJ(6);
	t421 = sin(t425);
	t422 = cos(t425);
	t448 = t422 * r_i_i_C(1) - t421 * r_i_i_C(2);
	t444 = t418 + t448;
	t485 = r_i_i_C(3) + pkin(11) + pkin(10);
	t493 = (t444 * t419 - t485 * t420) * qJD(4);
	t408 = t434 * t430 + t433 * t460;
	t395 = t408 * qJD(1) + t407 * qJD(2);
	t401 = -t407 * t420 + t419 * t477;
	t492 = -t401 * t424 - t395;
	t447 = t421 * r_i_i_C(1) + t422 * r_i_i_C(2);
	t429 = sin(qJ(5));
	t486 = t429 * pkin(5);
	t491 = qJD(5) * t486 + t447 * t424;
	t417 = cos(pkin(12)) * pkin(3) + pkin(2);
	t489 = t485 * t419 + t444 * t420 + t417;
	t482 = t421 * t424;
	t481 = t422 * t424;
	t480 = t427 * t430;
	t479 = t427 * t431;
	t478 = t427 * t433;
	t469 = qJD(2) * t433;
	t393 = -qJD(1) * t449 - t434 * t469 + t445 * t476;
	t409 = -t450 + t475;
	t403 = t409 * t420 + t419 * t479;
	t455 = -t403 * t424 - t393;
	t394 = t407 * qJD(1) + t408 * qJD(2);
	t443 = -t409 * t419 + t420 * t479;
	t463 = t434 * t471;
	t388 = t443 * qJD(4) - t394 * t420 + t419 * t463;
	t458 = t408 * t424 + t388;
	t383 = -t458 * t421 + t455 * t422;
	t384 = t455 * t421 + t458 * t422;
	t474 = t383 * r_i_i_C(1) - t384 * r_i_i_C(2);
	t473 = (-t456 * t421 - t422 * t492) * r_i_i_C(1) + (t421 * t492 - t456 * t422) * r_i_i_C(2);
	t405 = t484 * t419 + t420 * t480;
	t461 = t427 * t470;
	t442 = -t405 * t424 + t461;
	t440 = -t419 * t480 + t484 * t420;
	t462 = t427 * t469;
	t398 = t440 * qJD(4) + t420 * t462;
	t446 = t424 * t478 - t398;
	t472 = (t446 * t421 + t442 * t422) * r_i_i_C(1) + (-t442 * t421 + t446 * t422) * r_i_i_C(2);
	t468 = qJD(5) * t432;
	t465 = pkin(3) * sin(pkin(12)) + pkin(8);
	t437 = t401 * qJD(4) - t396 * t419 + t420 * t464;
	t436 = t491 * t420 + t493;
	t428 = -pkin(9) - qJ(3);
	t387 = t403 * qJD(4) - t394 * t419 - t420 * t463;
	t1 = [-t406 * qJD(3) - t390 * t418 + t395 * t428 - t396 * t417 + (-t456 * r_i_i_C(1) + r_i_i_C(2) * t492) * t422 + (r_i_i_C(1) * t492 + t456 * r_i_i_C(2)) * t421 + t485 * t437 + (-t434 * pkin(1) - t465 * t479) * qJD(1) + (-t395 * t429 + (-t401 * t429 - t406 * t432) * qJD(5)) * pkin(5), (-t394 * t421 + t409 * t481) * r_i_i_C(1) + (-t394 * t422 - t409 * t482) * r_i_i_C(2) + t394 * t428 + t409 * qJD(3) + (-t394 * t429 + t409 * t468) * pkin(5) + t489 * t393 + t436 * t408, -t393, -t444 * t387 + t485 * t388 - t443 * t491, (-t388 * t429 - t393 * t432 + (-t403 * t432 - t408 * t429) * qJD(5)) * pkin(5) + t474, t474; t384 * r_i_i_C(1) + t383 * r_i_i_C(2) + t408 * qJD(3) + t388 * t418 + t393 * t428 - t394 * t417 + t485 * t387 + (-pkin(1) * t431 + t465 * t477) * qJD(1) + (-t393 * t429 + (-t403 * t429 + t408 * t432) * qJD(5)) * pkin(5), (t396 * t421 + t407 * t481) * r_i_i_C(1) + (t396 * t422 - t407 * t482) * r_i_i_C(2) - t396 * t428 + t407 * qJD(3) + (t396 * t429 + t407 * t468) * pkin(5) - t489 * t395 + t436 * t406, t395, t485 * t390 - t491 * (-t407 * t419 - t466) + t444 * t437, (-t390 * t429 + t395 * t432 + (t401 * t432 - t406 * t429) * qJD(5)) * pkin(5) + t473, t473; 0, ((pkin(5) * t468 - qJD(2) * t489 + t448 * t424 + qJD(3)) * t430 + (-qJD(2) * t428 + (-qJD(5) * t420 + qJD(2)) * t486 - t493 + t447 * (-t420 * t424 + qJD(2))) * t433) * t427, t461, t485 * t398 - t491 * t440 + t444 * (-t405 * qJD(4) - t419 * t462), (t432 * t461 - t398 * t429 + (-t405 * t432 + t429 * t478) * qJD(5)) * pkin(5) + t472, t472;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end