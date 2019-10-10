% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:29:20
	% EndTime: 2019-10-10 12:29:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:29:20
	% EndTime: 2019-10-10 12:29:20
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
	% StartTime: 2019-10-10 12:29:20
	% EndTime: 2019-10-10 12:29:20
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
	% StartTime: 2019-10-10 12:29:21
	% EndTime: 2019-10-10 12:29:21
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (145->51), mult. (441->91), div. (0->0), fcn. (412->8), ass. (0->38)
	t242 = cos(pkin(6));
	t244 = sin(qJ(2));
	t245 = sin(qJ(1));
	t264 = t245 * t244;
	t258 = t242 * t264;
	t247 = cos(qJ(2));
	t248 = cos(qJ(1));
	t261 = t248 * t247;
	t233 = -qJD(1) * t258 - qJD(2) * t264 + (qJD(2) * t242 + qJD(1)) * t261;
	t241 = sin(pkin(6));
	t265 = t241 * t248;
	t256 = qJD(3) * t265;
	t269 = t233 - t256;
	t268 = -r_i_i_C(3) - pkin(9);
	t267 = t241 * t245;
	t246 = cos(qJ(3));
	t266 = t241 * t246;
	t263 = t245 * t247;
	t262 = t248 * t244;
	t260 = qJD(1) * t241;
	t235 = t242 * t262 + t263;
	t259 = qJD(3) * t235;
	t257 = t248 * t260;
	t243 = sin(qJ(3));
	t255 = -r_i_i_C(1) * t243 - r_i_i_C(2) * t246;
	t254 = t246 * r_i_i_C(1) - t243 * r_i_i_C(2) + pkin(2);
	t253 = t242 * t261 - t264;
	t252 = t242 * t263 + t262;
	t251 = t258 - t261;
	t250 = qJD(3) * t255;
	t249 = t245 * t260 - t259;
	t238 = t246 * t256;
	t232 = t252 * qJD(1) + t235 * qJD(2);
	t231 = t235 * qJD(1) + t252 * qJD(2);
	t230 = -t253 * qJD(1) + t251 * qJD(2);
	t229 = t243 * t257 - t231 * t246 + (t243 * t251 + t245 * t266) * qJD(3);
	t228 = t246 * t257 + t231 * t243 + (-t243 * t267 + t246 * t251) * qJD(3);
	t1 = [(-t233 * t246 + t243 * t259 + t238) * r_i_i_C(1) + (t269 * t243 + t246 * t259) * r_i_i_C(2) - t233 * pkin(2) + t268 * t232 + (-t248 * pkin(1) + (-pkin(8) + t255) * t267) * qJD(1), t254 * t230 + t268 * t231 - t252 * t250, t228 * r_i_i_C(1) - t229 * r_i_i_C(2), 0, 0, 0; -t231 * pkin(2) + t229 * r_i_i_C(1) + t228 * r_i_i_C(2) + t268 * t230 + (-pkin(1) * t245 + pkin(8) * t265) * qJD(1), -t254 * t232 - t268 * t233 + t253 * t250, t238 * r_i_i_C(2) + (t249 * r_i_i_C(1) - t233 * r_i_i_C(2)) * t246 + (-t269 * r_i_i_C(1) - t249 * r_i_i_C(2)) * t243, 0, 0, 0; 0, (t247 * t250 + (-t254 * t244 - t268 * t247) * qJD(2)) * t241, t255 * t247 * t241 * qJD(2) + ((-t242 * t243 - t244 * t266) * r_i_i_C(1) + (t241 * t243 * t244 - t242 * t246) * r_i_i_C(2)) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:29:23
	% EndTime: 2019-10-10 12:29:23
	% DurationCPUTime: 0.78s
	% Computational Cost: add. (444->96), mult. (1331->167), div. (0->0), fcn. (1334->10), ass. (0->65)
	t374 = sin(qJ(1));
	t370 = cos(pkin(6));
	t390 = qJD(2) * t370 + qJD(1);
	t373 = sin(qJ(2));
	t408 = t374 * t373;
	t398 = t370 * t408;
	t403 = qJD(2) * t373;
	t377 = cos(qJ(2));
	t378 = cos(qJ(1));
	t405 = t378 * t377;
	t348 = -qJD(1) * t398 - t374 * t403 + t390 * t405;
	t406 = t378 * t373;
	t407 = t374 * t377;
	t359 = t370 * t406 + t407;
	t372 = sin(qJ(3));
	t376 = cos(qJ(3));
	t369 = sin(pkin(6));
	t404 = qJD(1) * t369;
	t395 = t374 * t404;
	t409 = t369 * t378;
	t397 = t376 * t409;
	t342 = t372 * (-qJD(3) * t359 + t395) - qJD(3) * t397 + t348 * t376;
	t360 = t370 * t407 + t406;
	t347 = qJD(1) * t360 + t359 * qJD(2);
	t371 = sin(qJ(4));
	t375 = cos(qJ(4));
	t426 = t342 * t371 - t347 * t375;
	t425 = -t342 * t375 - t347 * t371;
	t388 = t375 * r_i_i_C(1) - t371 * r_i_i_C(2);
	t386 = pkin(3) + t388;
	t419 = pkin(10) + r_i_i_C(3);
	t424 = (t386 * t372 - t419 * t376) * qJD(3);
	t353 = -t359 * t376 + t372 * t409;
	t396 = t370 * t405;
	t358 = -t396 + t408;
	t423 = -t353 * t371 - t358 * t375;
	t422 = t353 * t375 - t358 * t371;
	t387 = t371 * r_i_i_C(1) + t375 * r_i_i_C(2);
	t420 = t419 * t372 + t386 * t376 + pkin(2);
	t412 = t369 * t374;
	t411 = t369 * t376;
	t410 = t369 * t377;
	t402 = qJD(2) * t377;
	t401 = qJD(4) * t371;
	t400 = qJD(4) * t375;
	t399 = qJD(4) * t376;
	t394 = t378 * t404;
	t393 = t369 * t403;
	t392 = t369 * t402;
	t361 = -t398 + t405;
	t385 = -t361 * t372 + t374 * t411;
	t355 = t361 * t376 + t372 * t412;
	t357 = t370 * t372 + t373 * t411;
	t384 = -t369 * t373 * t372 + t370 * t376;
	t383 = qJD(4) * t387;
	t380 = t353 * qJD(3) - t348 * t372 + t376 * t395;
	t379 = t387 * t399 + t424;
	t350 = t384 * qJD(3) + t376 * t392;
	t346 = t359 * qJD(1) + qJD(2) * t360;
	t345 = -qJD(1) * t396 - t378 * t402 + t390 * t408;
	t340 = t385 * qJD(3) - t346 * t376 + t372 * t394;
	t339 = t355 * qJD(3) - t346 * t372 - t376 * t394;
	t338 = t340 * t375 - t345 * t371 + (-t355 * t371 + t360 * t375) * qJD(4);
	t337 = -t340 * t371 - t345 * t375 + (-t355 * t375 - t360 * t371) * qJD(4);
	t1 = [t425 * r_i_i_C(1) + t426 * r_i_i_C(2) - t342 * pkin(3) - t348 * pkin(2) - t347 * pkin(9) + t419 * t380 + (t423 * r_i_i_C(1) - t422 * r_i_i_C(2)) * qJD(4) + (-t378 * pkin(1) - pkin(8) * t412) * qJD(1), (-t346 * t371 + t361 * t400) * r_i_i_C(1) + (-t346 * t375 - t361 * t401) * r_i_i_C(2) - t346 * pkin(9) + t420 * t345 + t379 * t360, -t386 * t339 + t419 * t340 - t385 * t383, t337 * r_i_i_C(1) - t338 * r_i_i_C(2), 0, 0; -t346 * pkin(2) + t340 * pkin(3) - t345 * pkin(9) + t338 * r_i_i_C(1) + t337 * r_i_i_C(2) + t419 * t339 + (-pkin(1) * t374 + pkin(8) * t409) * qJD(1), (t348 * t371 + t359 * t400) * r_i_i_C(1) + (t348 * t375 - t359 * t401) * r_i_i_C(2) + t348 * pkin(9) - t420 * t347 + t379 * t358, t419 * t342 - (-t359 * t372 - t397) * t383 + t386 * t380, -t426 * r_i_i_C(1) + t425 * r_i_i_C(2) + (t422 * r_i_i_C(1) + t423 * r_i_i_C(2)) * qJD(4), 0, 0; 0, ((-qJD(2) * t420 + t388 * qJD(4)) * t373 + (qJD(2) * pkin(9) - t424 + t387 * (qJD(2) - t399)) * t377) * t369, t419 * t350 - t384 * t383 + t386 * (-t357 * qJD(3) - t372 * t392), (-t350 * t371 + t375 * t393) * r_i_i_C(1) + (-t350 * t375 - t371 * t393) * r_i_i_C(2) + ((-t357 * t375 + t371 * t410) * r_i_i_C(1) + (t357 * t371 + t375 * t410) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:29:23
	% EndTime: 2019-10-10 12:29:23
	% DurationCPUTime: 0.89s
	% Computational Cost: add. (669->126), mult. (1706->207), div. (0->0), fcn. (1716->12), ass. (0->69)
	t386 = sin(qJ(1));
	t389 = cos(qJ(2));
	t426 = cos(pkin(6));
	t432 = cos(qJ(1));
	t400 = t426 * t432;
	t385 = sin(qJ(2));
	t406 = t386 * t426;
	t401 = t385 * t406;
	t407 = t432 * qJD(1);
	t419 = qJD(2) * t385;
	t353 = -qJD(1) * t401 - t386 * t419 + (qJD(2) * t400 + t407) * t389;
	t381 = sin(pkin(6));
	t412 = t381 * t432;
	t439 = -qJD(3) * t412 + t353;
	t364 = t385 * t400 + t386 * t389;
	t423 = t381 * t386;
	t438 = qJD(1) * t423 - qJD(3) * t364;
	t384 = sin(qJ(3));
	t388 = cos(qJ(3));
	t357 = t364 * t388 - t384 * t412;
	t398 = t389 * t400;
	t420 = t385 * t386;
	t363 = -t398 + t420;
	t380 = qJ(4) + pkin(11);
	t378 = sin(t380);
	t379 = cos(t380);
	t437 = t357 * t378 - t363 * t379;
	t436 = t357 * t379 + t363 * t378;
	t387 = cos(qJ(4));
	t431 = pkin(4) * t387;
	t377 = pkin(3) + t431;
	t399 = r_i_i_C(1) * t379 - r_i_i_C(2) * t378;
	t397 = t377 + t399;
	t427 = r_i_i_C(3) + qJ(5) + pkin(10);
	t435 = -(t397 * t384 - t427 * t388) * qJD(3) + t384 * qJD(5);
	t383 = sin(qJ(4));
	t395 = t383 * pkin(4) + t378 * r_i_i_C(1) + t379 * r_i_i_C(2);
	t347 = t438 * t384 + t439 * t388;
	t433 = t427 * t384 + t397 * t388 + pkin(2);
	t422 = t381 * t388;
	t421 = t381 * t389;
	t417 = qJD(4) * t378;
	t416 = qJD(4) * t379;
	t415 = qJD(4) * t387;
	t414 = qJD(4) * t388;
	t411 = t432 * t389;
	t409 = t381 * t419;
	t408 = qJD(2) * t421;
	t403 = t381 * t407;
	t366 = t411 - t401;
	t359 = -t366 * t384 + t386 * t422;
	t360 = t366 * t388 + t384 * t423;
	t394 = t364 * t384 + t388 * t412;
	t393 = -t381 * t385 * t384 + t426 * t388;
	t362 = t426 * t384 + t385 * t422;
	t346 = t439 * t384 - t438 * t388;
	t392 = qJD(4) * t395;
	t365 = t432 * t385 + t389 * t406;
	t390 = t395 * t414 - t435;
	t355 = t393 * qJD(3) + t388 * t408;
	t354 = t362 * qJD(3) + t384 * t408;
	t352 = t365 * qJD(1) + t364 * qJD(2);
	t351 = t364 * qJD(1) + t365 * qJD(2);
	t350 = -qJD(1) * t398 - qJD(2) * t411 + (qJD(2) * t426 + qJD(1)) * t420;
	t345 = t359 * qJD(3) - t351 * t388 + t384 * t403;
	t344 = t360 * qJD(3) - t351 * t384 - t388 * t403;
	t343 = t345 * t379 - t350 * t378 + (-t360 * t378 + t365 * t379) * qJD(4);
	t342 = -t345 * t378 - t350 * t379 + (-t360 * t379 - t365 * t378) * qJD(4);
	t1 = [-t394 * qJD(5) - t353 * pkin(2) - t397 * t347 + (-pkin(9) - t395) * t352 - t427 * t346 + (-t432 * pkin(1) - pkin(8) * t423) * qJD(1) + (t437 * r_i_i_C(1) + t436 * r_i_i_C(2) + (t357 * t383 - t363 * t387) * pkin(4)) * qJD(4), (-t351 * t378 + t366 * t416) * r_i_i_C(1) + (-t351 * t379 - t366 * t417) * r_i_i_C(2) - t351 * pkin(9) + (-t351 * t383 + t366 * t415) * pkin(4) + t433 * t350 + t390 * t365, qJD(5) * t360 - t397 * t344 + t427 * t345 - t359 * t392, r_i_i_C(1) * t342 - r_i_i_C(2) * t343 + (-t345 * t383 - t350 * t387 + (-t360 * t387 - t365 * t383) * qJD(4)) * pkin(4), t344, 0; -t351 * pkin(2) - t350 * pkin(9) + t343 * r_i_i_C(1) + t342 * r_i_i_C(2) - t359 * qJD(5) + t345 * t377 + t427 * t344 + (-pkin(1) * t386 + pkin(8) * t412) * qJD(1) + (-t350 * t383 + (-t360 * t383 + t365 * t387) * qJD(4)) * pkin(4), (t353 * t378 + t364 * t416) * r_i_i_C(1) + (t353 * t379 - t364 * t417) * r_i_i_C(2) + t353 * pkin(9) + (t353 * t383 + t364 * t415) * pkin(4) - t433 * t352 + t390 * t363, qJD(5) * t357 - t397 * t346 + t427 * t347 + t394 * t392, (-t347 * t378 + t352 * t379) * r_i_i_C(1) + (-t347 * t379 - t352 * t378) * r_i_i_C(2) + (-t436 * r_i_i_C(1) + t437 * r_i_i_C(2)) * qJD(4) + (-t347 * t383 + t352 * t387 + (-t357 * t387 - t363 * t383) * qJD(4)) * pkin(4), t346, 0; 0, (((t399 + t431) * qJD(4) - t433 * qJD(2)) * t385 + (qJD(2) * pkin(9) + t395 * (qJD(2) - t414) + t435) * t389) * t381, qJD(5) * t362 - t397 * t354 + t427 * t355 - t393 * t392, (-t355 * t378 + t379 * t409) * r_i_i_C(1) + (-t355 * t379 - t378 * t409) * r_i_i_C(2) + ((-t362 * t379 + t378 * t421) * r_i_i_C(1) + (t362 * t378 + t379 * t421) * r_i_i_C(2)) * qJD(4) + (t387 * t409 - t355 * t383 + (-t362 * t387 + t383 * t421) * qJD(4)) * pkin(4), t354, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:29:25
	% EndTime: 2019-10-10 12:29:26
	% DurationCPUTime: 1.20s
	% Computational Cost: add. (1131->149), mult. (2714->231), div. (0->0), fcn. (2795->12), ass. (0->81)
	t474 = sin(qJ(3));
	t478 = cos(qJ(3));
	t475 = sin(qJ(2));
	t476 = sin(qJ(1));
	t479 = cos(qJ(2));
	t531 = cos(pkin(6));
	t535 = cos(qJ(1));
	t498 = t531 * t535;
	t452 = t475 * t498 + t476 * t479;
	t471 = sin(pkin(6));
	t525 = t471 * t476;
	t544 = qJD(1) * t525 - qJD(3) * t452;
	t505 = t476 * t531;
	t499 = t475 * t505;
	t507 = t535 * qJD(1);
	t521 = qJD(2) * t475;
	t441 = -qJD(1) * t499 - t476 * t521 + (qJD(2) * t498 + t507) * t479;
	t513 = t471 * t535;
	t545 = -qJD(3) * t513 + t441;
	t431 = t544 * t474 + t545 * t478;
	t453 = t535 * t475 + t479 * t505;
	t440 = t453 * qJD(1) + t452 * qJD(2);
	t470 = qJ(4) + pkin(11);
	t468 = sin(t470);
	t469 = cos(t470);
	t445 = t452 * t478 - t474 * t513;
	t490 = t479 * t498;
	t522 = t476 * t475;
	t451 = -t490 + t522;
	t495 = t445 * t469 + t451 * t468;
	t418 = t495 * qJD(4) + t431 * t468 - t440 * t469;
	t494 = t445 * t468 - t451 * t469;
	t546 = -t494 * qJD(4) + t431 * t469 + t440 * t468;
	t477 = cos(qJ(4));
	t467 = t477 * pkin(4) + pkin(3);
	t516 = t468 * qJD(6);
	t533 = r_i_i_C(2) + qJ(5) + pkin(10);
	t473 = sin(qJ(4));
	t534 = t473 * pkin(4);
	t480 = -(t533 * qJD(3) - qJD(4) * t534 + t516) * t478 + (qJD(3) * t467 - qJD(5)) * t474;
	t524 = t471 * t478;
	t450 = t531 * t474 + t475 * t524;
	t523 = t471 * t479;
	t540 = -t450 * t469 + t468 * t523;
	t539 = t478 * t467 + t533 * t474 + pkin(2);
	t519 = qJD(3) * t474;
	t538 = (qJD(2) * t478 - qJD(4)) * t475 + t479 * t519;
	t536 = r_i_i_C(1) + pkin(5);
	t532 = r_i_i_C(3) + qJ(6);
	t518 = qJD(4) * t477;
	t517 = qJD(4) * t478;
	t515 = t469 * qJD(6);
	t512 = t535 * t479;
	t510 = t471 * t521;
	t509 = qJD(2) * t523;
	t500 = t471 * t507;
	t439 = t452 * qJD(1) + t453 * qJD(2);
	t497 = t453 * t517 - t439;
	t496 = t451 * t517 + t441;
	t454 = t512 - t499;
	t448 = t454 * t478 + t474 * t525;
	t493 = t448 * t469 + t453 * t468;
	t492 = -t448 * t468 + t453 * t469;
	t491 = (qJD(2) - t517) * t479;
	t447 = -t454 * t474 + t476 * t524;
	t486 = t452 * t474 + t478 * t513;
	t485 = -t471 * t475 * t474 + t531 * t478;
	t430 = t545 * t474 - t544 * t478;
	t484 = -t532 * t468 - t536 * t469 - t467;
	t438 = -qJD(1) * t490 - qJD(2) * t512 + (qJD(2) * t531 + qJD(1)) * t522;
	t483 = qJD(4) * t454 + t438 * t478 + t453 * t519;
	t482 = qJD(4) * t452 - t440 * t478 + t451 * t519;
	t481 = t516 + (-t536 * t468 + t532 * t469 - t534) * qJD(4);
	t443 = t485 * qJD(3) + t478 * t509;
	t442 = t450 * qJD(3) + t474 * t509;
	t429 = t447 * qJD(3) - t439 * t478 + t474 * t500;
	t428 = t448 * qJD(3) - t439 * t474 - t478 * t500;
	t426 = -t540 * qJD(4) + t443 * t468 - t469 * t510;
	t417 = t492 * qJD(4) + t429 * t469 - t438 * t468;
	t416 = t493 * qJD(4) + t429 * t468 + t438 * t469;
	t1 = [-t494 * qJD(6) - t431 * t467 - t486 * qJD(5) - t441 * pkin(2) - t440 * pkin(9) - t533 * t430 - t536 * t546 - t532 * t418 + (-t535 * pkin(1) - pkin(8) * t525) * qJD(1) + (-t440 * t473 + (t445 * t473 - t451 * t477) * qJD(4)) * pkin(4), -t454 * t515 - t439 * pkin(9) + t536 * (t497 * t468 + t483 * t469) + t532 * (t483 * t468 - t497 * t469) + (-t439 * t473 + t454 * t518) * pkin(4) + t539 * t438 + t480 * t453, t448 * qJD(5) + t484 * t428 + t533 * t429 + t481 * t447, t493 * qJD(6) + t532 * t417 - t536 * t416 + (-t429 * t473 - t438 * t477 + (-t448 * t477 - t453 * t473) * qJD(4)) * pkin(4), t428, t416; -t492 * qJD(6) + t429 * t467 - t447 * qJD(5) - t439 * pkin(2) - t438 * pkin(9) + t533 * t428 + t536 * t417 + t532 * t416 + (-t476 * pkin(1) + pkin(8) * t513) * qJD(1) + (-t438 * t473 + (-t448 * t473 + t453 * t477) * qJD(4)) * pkin(4), -t452 * t515 + t441 * pkin(9) + t536 * (t496 * t468 + t482 * t469) + t532 * (t482 * t468 - t496 * t469) + (t441 * t473 + t452 * t518) * pkin(4) - t539 * t440 + t480 * t451, t445 * qJD(5) + t484 * t430 + t533 * t431 - t481 * t486, t495 * qJD(6) + t532 * t546 - t536 * t418 + (-t431 * t473 + t440 * t477 + (-t445 * t477 - t451 * t473) * qJD(4)) * pkin(4), t430, t418; 0, (t536 * (t468 * t491 - t538 * t469) - t532 * (t538 * t468 + t469 * t491) + (pkin(4) * t518 - qJD(2) * t539 - t515) * t475 + ((pkin(9) + t534) * qJD(2) - t480) * t479) * t471, t450 * qJD(5) + t484 * t442 + t533 * t443 + t481 * t485, -t540 * qJD(6) + t532 * (t468 * t510 + t443 * t469 + (-t450 * t468 - t469 * t523) * qJD(4)) - t536 * t426 + (t477 * t510 - t443 * t473 + (-t450 * t477 + t473 * t523) * qJD(4)) * pkin(4), t442, t426;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end