% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRP9
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
% Datum: 2019-10-10 10:42
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP9_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP9_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP9_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:42:53
	% EndTime: 2019-10-10 10:42:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:42:53
	% EndTime: 2019-10-10 10:42:53
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
	% StartTime: 2019-10-10 10:42:54
	% EndTime: 2019-10-10 10:42:54
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
	% StartTime: 2019-10-10 10:42:54
	% EndTime: 2019-10-10 10:42:54
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (103->31), mult. (307->49), div. (0->0), fcn. (282->8), ass. (0->26)
	t199 = sin(pkin(11));
	t200 = sin(pkin(6));
	t201 = cos(pkin(11));
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
	% StartTime: 2019-10-10 10:42:55
	% EndTime: 2019-10-10 10:42:55
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
	t253 = pkin(11) + qJ(4);
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
	t270 = pkin(3) * sin(pkin(11)) + pkin(8);
	t269 = t261 * t276;
	t250 = cos(pkin(11)) * pkin(3) + pkin(2);
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
	% StartTime: 2019-10-10 10:42:56
	% EndTime: 2019-10-10 10:42:57
	% DurationCPUTime: 0.78s
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
	t385 = pkin(11) + qJ(4);
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
	t437 = pkin(10) + r_i_i_C(3);
	t442 = (t403 * t383 - t437 * t384) * qJD(4);
	t366 = -t372 * t384 + t383 * t427;
	t414 = t388 * t423;
	t371 = -t414 + t426;
	t441 = -t366 * t390 - t371 * t393;
	t440 = t366 * t393 - t371 * t390;
	t404 = t390 * r_i_i_C(1) + t393 * r_i_i_C(2);
	t382 = cos(pkin(11)) * pkin(3) + pkin(2);
	t438 = t437 * t383 + t403 * t384 + t382;
	t430 = t387 * t391;
	t429 = t387 * t392;
	t428 = t387 * t394;
	t420 = qJD(2) * t394;
	t419 = qJD(5) * t384;
	t418 = qJD(5) * t390;
	t417 = qJD(5) * t393;
	t413 = pkin(3) * sin(pkin(11)) + pkin(8);
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
	% StartTime: 2019-10-10 10:42:56
	% EndTime: 2019-10-10 10:42:57
	% DurationCPUTime: 0.90s
	% Computational Cost: add. (853->119), mult. (1759->184), div. (0->0), fcn. (1771->12), ass. (0->73)
	t398 = sin(qJ(1));
	t400 = cos(qJ(2));
	t443 = cos(pkin(6));
	t448 = cos(qJ(1));
	t418 = t443 * t448;
	t397 = sin(qJ(2));
	t424 = t398 * t443;
	t419 = t397 * t424;
	t425 = t448 * qJD(1);
	t436 = qJD(2) * t397;
	t363 = -qJD(1) * t419 - t398 * t436 + (qJD(2) * t418 + t425) * t400;
	t393 = sin(pkin(6));
	t430 = t393 * t448;
	t458 = -qJD(4) * t430 + t363;
	t374 = t397 * t418 + t398 * t400;
	t439 = t393 * t398;
	t457 = qJD(1) * t439 - qJD(4) * t374;
	t391 = pkin(11) + qJ(4);
	t389 = sin(t391);
	t390 = cos(t391);
	t367 = t374 * t390 - t389 * t430;
	t414 = t400 * t418;
	t437 = t398 * t397;
	t373 = -t414 + t437;
	t396 = sin(qJ(5));
	t399 = cos(qJ(5));
	t456 = t367 * t396 - t373 * t399;
	t455 = t367 * t399 + t373 * t396;
	t449 = pkin(5) + r_i_i_C(1);
	t452 = t399 * r_i_i_C(2) + t449 * t396;
	t388 = t399 * pkin(5) + pkin(4);
	t446 = t396 * r_i_i_C(2);
	t413 = t399 * r_i_i_C(1) + t388 - t446;
	t444 = r_i_i_C(3) + qJ(6) + pkin(10);
	t454 = -(t413 * t389 - t444 * t390) * qJD(4) + t389 * qJD(6);
	t453 = pkin(8) + pkin(3) * sin(pkin(11));
	t357 = t457 * t389 + t458 * t390;
	t387 = cos(pkin(11)) * pkin(3) + pkin(2);
	t450 = t444 * t389 + t413 * t390 + t387;
	t440 = t393 * t397;
	t438 = t393 * t400;
	t434 = qJD(5) * t390;
	t433 = qJD(5) * t396;
	t432 = qJD(5) * t399;
	t429 = t448 * t400;
	t427 = qJD(2) * t438;
	t426 = t393 * t436;
	t420 = t393 * t425;
	t375 = t448 * t397 + t400 * t424;
	t362 = t375 * qJD(1) + t374 * qJD(2);
	t417 = -t357 * t396 + t362 * t399;
	t372 = t443 * t389 + t390 * t440;
	t411 = -t372 * t399 + t396 * t438;
	t376 = t429 - t419;
	t369 = -t376 * t389 + t390 * t439;
	t370 = t376 * t390 + t389 * t439;
	t407 = t374 * t389 + t390 * t430;
	t404 = -t389 * t440 + t443 * t390;
	t365 = t404 * qJD(4) + t390 * t427;
	t406 = -t365 * t396 + t399 * t426;
	t405 = qJD(5) * t452;
	t356 = t458 * t389 - t457 * t390;
	t360 = -qJD(1) * t414 - qJD(2) * t429 + (qJD(2) * t443 + qJD(1)) * t437;
	t403 = -t360 * t396 + (-t370 * t396 + t375 * t399) * qJD(5);
	t361 = t374 * qJD(1) + t375 * qJD(2);
	t355 = t369 * qJD(4) - t361 * t390 + t389 * t420;
	t352 = -t355 * t396 - t360 * t399 + (-t370 * t399 - t375 * t396) * qJD(5);
	t401 = t452 * t434 - t454;
	t395 = -pkin(9) - qJ(3);
	t364 = t372 * qJD(4) + t389 * t427;
	t354 = t370 * qJD(4) - t361 * t389 - t390 * t420;
	t353 = t355 * t399 + t403;
	t1 = [-t407 * qJD(6) - t363 * t387 - t373 * qJD(3) - t413 * t357 + (t395 - t452) * t362 - t444 * t356 + (-t448 * pkin(1) - t453 * t439) * qJD(1) + (t455 * r_i_i_C(2) + t449 * t456) * qJD(5), (-t361 * t399 - t376 * t433) * r_i_i_C(2) + t361 * t395 + t376 * qJD(3) + t450 * t360 + t401 * t375 + t449 * (-t361 * t396 + t376 * t432), -t360, t370 * qJD(6) - t413 * t354 + t444 * t355 - t369 * t405, -t353 * r_i_i_C(2) + t449 * t352, t354; t353 * r_i_i_C(1) + t352 * r_i_i_C(2) + t375 * qJD(3) - t369 * qJD(6) + t355 * t388 + t360 * t395 - t361 * t387 + t444 * t354 + (-pkin(1) * t398 + t453 * t430) * qJD(1) + t403 * pkin(5), (t363 * t399 - t374 * t433) * r_i_i_C(2) - t363 * t395 + t374 * qJD(3) - t450 * t362 + t401 * t373 + t449 * (t363 * t396 + t374 * t432), t362, t367 * qJD(6) - t413 * t356 + t444 * t357 + t407 * t405, t417 * r_i_i_C(1) + (-t357 * t399 - t362 * t396) * r_i_i_C(2) + (-r_i_i_C(1) * t455 + t456 * r_i_i_C(2)) * qJD(5) + (-qJD(5) * t455 + t417) * pkin(5), t356; 0, ((qJD(3) + (t449 * t399 - t446) * qJD(5) - t450 * qJD(2)) * t397 + (-qJD(2) * t395 + t452 * (qJD(2) - t434) + t454) * t400) * t393, t426, t372 * qJD(6) - t413 * t364 + t444 * t365 - t404 * t405, t406 * r_i_i_C(1) + (-t365 * t399 - t396 * t426) * r_i_i_C(2) + (t411 * r_i_i_C(1) + (t372 * t396 + t399 * t438) * r_i_i_C(2)) * qJD(5) + (t411 * qJD(5) + t406) * pkin(5), t364;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end