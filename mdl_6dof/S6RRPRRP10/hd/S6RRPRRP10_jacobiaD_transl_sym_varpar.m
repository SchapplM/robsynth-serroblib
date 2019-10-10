% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRP10
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
% Datum: 2019-10-10 10:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
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
	% StartTime: 2019-10-10 10:44:48
	% EndTime: 2019-10-10 10:44:48
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
	% StartTime: 2019-10-10 10:44:49
	% EndTime: 2019-10-10 10:44:49
	% DurationCPUTime: 0.18s
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
	% StartTime: 2019-10-10 10:44:49
	% EndTime: 2019-10-10 10:44:49
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
	% StartTime: 2019-10-10 10:44:50
	% EndTime: 2019-10-10 10:44:51
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
	% StartTime: 2019-10-10 10:44:53
	% EndTime: 2019-10-10 10:44:54
	% DurationCPUTime: 1.16s
	% Computational Cost: add. (1151->127), mult. (2392->198), div. (0->0), fcn. (2468->12), ass. (0->78)
	t471 = cos(pkin(6));
	t475 = sin(qJ(1));
	t474 = sin(qJ(2));
	t517 = t475 * t474;
	t507 = t471 * t517;
	t512 = qJD(2) * t474;
	t477 = cos(qJ(2));
	t478 = cos(qJ(1));
	t514 = t478 * t477;
	t444 = -qJD(1) * t507 - t475 * t512 + (qJD(2) * t471 + qJD(1)) * t514;
	t515 = t478 * t474;
	t516 = t475 * t477;
	t455 = t471 * t515 + t516;
	t468 = pkin(11) + qJ(4);
	t466 = sin(t468);
	t467 = cos(t468);
	t470 = sin(pkin(6));
	t513 = qJD(1) * t470;
	t504 = t475 * t513;
	t519 = t470 * t478;
	t506 = t467 * t519;
	t432 = t466 * (-qJD(4) * t455 + t504) - qJD(4) * t506 + t444 * t467;
	t456 = t471 * t516 + t515;
	t443 = qJD(1) * t456 + qJD(2) * t455;
	t473 = sin(qJ(5));
	t476 = cos(qJ(5));
	t449 = -t455 * t467 + t466 * t519;
	t454 = -t471 * t514 + t517;
	t537 = t449 * t476 - t454 * t473;
	t542 = t537 * qJD(5) - t432 * t473 + t443 * t476;
	t538 = t449 * t473 + t454 * t476;
	t541 = t538 * qJD(5) + t432 * t476 + t443 * t473;
	t522 = t470 * t474;
	t453 = t471 * t466 + t467 * t522;
	t518 = t473 * t477;
	t536 = -t453 * t476 + t470 * t518;
	t465 = cos(pkin(11)) * pkin(3) + pkin(2);
	t531 = r_i_i_C(2) + pkin(10);
	t535 = t467 * pkin(4) + t531 * t466 + t465;
	t534 = -pkin(4) * t466 + t531 * t467;
	t510 = qJD(4) * t477;
	t533 = (qJD(2) * t467 - qJD(5)) * t474 + t466 * t510;
	t529 = r_i_i_C(3) + qJ(6);
	t532 = r_i_i_C(1) + pkin(5);
	t483 = t529 * t473 + t532 * t476 + pkin(4);
	t523 = t467 * t473;
	t521 = t470 * t475;
	t520 = t470 * t477;
	t511 = qJD(4) * t466;
	t509 = qJD(5) * t467;
	t505 = pkin(3) * sin(pkin(11)) + pkin(8);
	t503 = t478 * t513;
	t502 = qJD(2) * t520;
	t501 = t470 * t512;
	t442 = qJD(1) * t455 + qJD(2) * t456;
	t496 = t456 * t509 - t442;
	t495 = t454 * t509 + t444;
	t487 = t507 - t514;
	t451 = t466 * t521 - t467 * t487;
	t492 = t451 * t476 + t456 * t473;
	t491 = -t451 * t473 + t456 * t476;
	t490 = (qJD(2) - t509) * t477;
	t489 = t466 * t487 + t467 * t521;
	t488 = -t466 * t522 + t471 * t467;
	t484 = qJD(4) * t534;
	t441 = qJD(1) * t454 + qJD(2) * t487;
	t482 = -qJD(5) * t487 + t441 * t467 + t456 * t511;
	t481 = qJD(5) * t455 - t443 * t467 + t454 * t511;
	t480 = qJD(6) * t473 + (-t532 * t473 + t529 * t476) * qJD(5);
	t479 = qJD(4) * t449 - t444 * t466 + t467 * t504;
	t472 = -pkin(9) - qJ(3);
	t446 = qJD(4) * t488 + t467 * t502;
	t435 = -t536 * qJD(5) + t446 * t473 - t476 * t501;
	t430 = t489 * qJD(4) - t442 * t467 + t466 * t503;
	t429 = qJD(4) * t451 - t442 * t466 - t467 * t503;
	t420 = qJD(5) * t491 + t430 * t476 - t441 * t473;
	t419 = qJD(5) * t492 + t430 * t473 + t441 * t476;
	t1 = [t538 * qJD(6) - t432 * pkin(4) - t444 * t465 + t443 * t472 - t454 * qJD(3) + t531 * t479 - t532 * t541 + t529 * t542 + (-t478 * pkin(1) - t505 * t521) * qJD(1), -(t456 * t523 - t476 * t487) * qJD(6) + t442 * t472 - t487 * qJD(3) + t532 * (t473 * t496 + t482 * t476) + t529 * (t482 * t473 - t476 * t496) - t456 * t484 + t535 * t441, -t441, -t483 * t429 + t531 * t430 + t480 * t489, t492 * qJD(6) - t532 * t419 + t529 * t420, t419; -t491 * qJD(6) + t430 * pkin(4) - t442 * t465 + t441 * t472 + t456 * qJD(3) + t531 * t429 + t532 * t420 + t529 * t419 + (-t475 * pkin(1) + t505 * t519) * qJD(1), -(t454 * t523 + t455 * t476) * qJD(6) - t444 * t472 + t455 * qJD(3) + t532 * (t473 * t495 + t476 * t481) + t529 * (t473 * t481 - t476 * t495) - t454 * t484 - t535 * t443, t443, t531 * t432 + t480 * (-t455 * t466 - t506) + t483 * t479, -t537 * qJD(6) + t529 * t541 + t532 * t542, -t542; 0, (t532 * (t473 * t490 - t533 * t476) - t529 * (t533 * t473 + t476 * t490) - (-t467 * t518 + t474 * t476) * qJD(6) + t474 * qJD(3) + t534 * t510 + (-t477 * t472 - t474 * t535) * qJD(2)) * t470, t501, t531 * t446 + t480 * t488 + t483 * (-qJD(4) * t453 - t466 * t502), -t536 * qJD(6) + t529 * (t473 * t501 + t446 * t476 + (-t453 * t473 - t476 * t520) * qJD(5)) - t532 * t435, t435;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end