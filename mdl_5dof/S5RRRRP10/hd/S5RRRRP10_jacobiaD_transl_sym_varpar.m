% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRP10
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRP10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRP10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:44:20
	% EndTime: 2019-12-29 20:44:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:44:20
	% EndTime: 2019-12-29 20:44:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:44:21
	% EndTime: 2019-12-29 20:44:21
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (35->18), mult. (110->35), div. (0->0), fcn. (94->6), ass. (0->20)
	t136 = sin(pkin(5));
	t151 = t136 * (pkin(7) + r_i_i_C(3));
	t138 = sin(qJ(2));
	t139 = sin(qJ(1));
	t149 = t138 * t139;
	t141 = cos(qJ(1));
	t148 = t138 * t141;
	t140 = cos(qJ(2));
	t147 = t139 * t140;
	t146 = t140 * t141;
	t137 = cos(pkin(5));
	t145 = -t137 * t146 + t149;
	t144 = t137 * t147 + t148;
	t143 = t137 * t148 + t147;
	t142 = t137 * t149 - t146;
	t135 = t142 * qJD(1) + t145 * qJD(2);
	t134 = t144 * qJD(1) + t143 * qJD(2);
	t133 = t143 * qJD(1) + t144 * qJD(2);
	t132 = t145 * qJD(1) + t142 * qJD(2);
	t1 = [t135 * r_i_i_C(1) + t134 * r_i_i_C(2) + (-pkin(1) * t141 - t139 * t151) * qJD(1), t132 * r_i_i_C(1) + t133 * r_i_i_C(2), 0, 0, 0; -t133 * r_i_i_C(1) + t132 * r_i_i_C(2) + (-pkin(1) * t139 + t141 * t151) * qJD(1), -t134 * r_i_i_C(1) + t135 * r_i_i_C(2), 0, 0, 0; 0, (-r_i_i_C(1) * t138 - r_i_i_C(2) * t140) * t136 * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:44:22
	% EndTime: 2019-12-29 20:44:22
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (145->51), mult. (441->91), div. (0->0), fcn. (412->8), ass. (0->38)
	t242 = cos(pkin(5));
	t244 = sin(qJ(2));
	t245 = sin(qJ(1));
	t264 = t245 * t244;
	t258 = t242 * t264;
	t247 = cos(qJ(2));
	t248 = cos(qJ(1));
	t261 = t248 * t247;
	t233 = -qJD(1) * t258 - qJD(2) * t264 + (qJD(2) * t242 + qJD(1)) * t261;
	t241 = sin(pkin(5));
	t265 = t241 * t248;
	t256 = qJD(3) * t265;
	t269 = t233 - t256;
	t268 = -r_i_i_C(3) - pkin(8);
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
	t1 = [(-t233 * t246 + t243 * t259 + t238) * r_i_i_C(1) + (t269 * t243 + t246 * t259) * r_i_i_C(2) - t233 * pkin(2) + t268 * t232 + (-t248 * pkin(1) + (-pkin(7) + t255) * t267) * qJD(1), t254 * t230 + t268 * t231 - t252 * t250, t228 * r_i_i_C(1) - t229 * r_i_i_C(2), 0, 0; -t231 * pkin(2) + t229 * r_i_i_C(1) + t228 * r_i_i_C(2) + t268 * t230 + (-pkin(1) * t245 + pkin(7) * t265) * qJD(1), -t254 * t232 - t268 * t233 + t253 * t250, t238 * r_i_i_C(2) + (t249 * r_i_i_C(1) - t233 * r_i_i_C(2)) * t246 + (-t269 * r_i_i_C(1) - t249 * r_i_i_C(2)) * t243, 0, 0; 0, (t247 * t250 + (-t254 * t244 - t268 * t247) * qJD(2)) * t241, t255 * t247 * t241 * qJD(2) + ((-t242 * t243 - t244 * t266) * r_i_i_C(1) + (t241 * t243 * t244 - t242 * t246) * r_i_i_C(2)) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:44:24
	% EndTime: 2019-12-29 20:44:25
	% DurationCPUTime: 1.16s
	% Computational Cost: add. (444->96), mult. (1331->167), div. (0->0), fcn. (1334->10), ass. (0->65)
	t374 = sin(qJ(1));
	t370 = cos(pkin(5));
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
	t369 = sin(pkin(5));
	t404 = qJD(1) * t369;
	t395 = t374 * t404;
	t409 = t369 * t378;
	t397 = t376 * t409;
	t342 = (-qJD(3) * t359 + t395) * t372 - qJD(3) * t397 + t348 * t376;
	t360 = t370 * t407 + t406;
	t347 = t360 * qJD(1) + t359 * qJD(2);
	t371 = sin(qJ(4));
	t375 = cos(qJ(4));
	t426 = t342 * t371 - t347 * t375;
	t425 = -t342 * t375 - t347 * t371;
	t388 = t375 * r_i_i_C(1) - t371 * r_i_i_C(2);
	t386 = pkin(3) + t388;
	t419 = pkin(9) + r_i_i_C(3);
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
	t346 = t359 * qJD(1) + t360 * qJD(2);
	t345 = -qJD(1) * t396 - t378 * t402 + t390 * t408;
	t340 = t385 * qJD(3) - t346 * t376 + t372 * t394;
	t339 = t355 * qJD(3) - t346 * t372 - t376 * t394;
	t338 = t340 * t375 - t345 * t371 + (-t355 * t371 + t360 * t375) * qJD(4);
	t337 = -t340 * t371 - t345 * t375 + (-t355 * t375 - t360 * t371) * qJD(4);
	t1 = [t425 * r_i_i_C(1) + t426 * r_i_i_C(2) - t342 * pkin(3) - t348 * pkin(2) - t347 * pkin(8) + t419 * t380 + (t423 * r_i_i_C(1) - t422 * r_i_i_C(2)) * qJD(4) + (-t378 * pkin(1) - pkin(7) * t412) * qJD(1), (-t346 * t371 + t361 * t400) * r_i_i_C(1) + (-t346 * t375 - t361 * t401) * r_i_i_C(2) - t346 * pkin(8) + t420 * t345 + t379 * t360, -t386 * t339 + t419 * t340 - t385 * t383, r_i_i_C(1) * t337 - r_i_i_C(2) * t338, 0; -t346 * pkin(2) + t340 * pkin(3) - t345 * pkin(8) + t338 * r_i_i_C(1) + t337 * r_i_i_C(2) + t419 * t339 + (-pkin(1) * t374 + pkin(7) * t409) * qJD(1), (t348 * t371 + t359 * t400) * r_i_i_C(1) + (t348 * t375 - t359 * t401) * r_i_i_C(2) + t348 * pkin(8) - t420 * t347 + t379 * t358, t419 * t342 - (-t359 * t372 - t397) * t383 + t386 * t380, -t426 * r_i_i_C(1) + t425 * r_i_i_C(2) + (t422 * r_i_i_C(1) + t423 * r_i_i_C(2)) * qJD(4), 0; 0, ((-qJD(2) * t420 + t388 * qJD(4)) * t373 + (qJD(2) * pkin(8) - t424 + t387 * (qJD(2) - t399)) * t377) * t369, t419 * t350 - t384 * t383 + t386 * (-t357 * qJD(3) - t372 * t392), (-t350 * t371 + t375 * t393) * r_i_i_C(1) + (-t350 * t375 - t371 * t393) * r_i_i_C(2) + ((-t357 * t375 + t371 * t410) * r_i_i_C(1) + (t357 * t371 + t375 * t410) * r_i_i_C(2)) * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:44:24
	% EndTime: 2019-12-29 20:44:25
	% DurationCPUTime: 1.26s
	% Computational Cost: add. (587->109), mult. (1706->179), div. (0->0), fcn. (1716->10), ass. (0->69)
	t380 = sin(qJ(1));
	t383 = cos(qJ(2));
	t426 = cos(pkin(5));
	t430 = cos(qJ(1));
	t401 = t426 * t430;
	t379 = sin(qJ(2));
	t407 = t380 * t426;
	t402 = t379 * t407;
	t408 = t430 * qJD(1);
	t419 = qJD(2) * t379;
	t350 = -qJD(1) * t402 - t380 * t419 + (qJD(2) * t401 + t408) * t383;
	t375 = sin(pkin(5));
	t413 = t375 * t430;
	t439 = -qJD(3) * t413 + t350;
	t361 = t379 * t401 + t380 * t383;
	t423 = t375 * t380;
	t438 = qJD(1) * t423 - qJD(3) * t361;
	t378 = sin(qJ(3));
	t382 = cos(qJ(3));
	t354 = t361 * t382 - t378 * t413;
	t397 = t383 * t401;
	t420 = t379 * t380;
	t360 = -t397 + t420;
	t377 = sin(qJ(4));
	t381 = cos(qJ(4));
	t437 = t354 * t377 - t360 * t381;
	t436 = t354 * t381 + t360 * t377;
	t431 = pkin(4) + r_i_i_C(1);
	t434 = t381 * r_i_i_C(2) + t431 * t377;
	t374 = pkin(4) * t381 + pkin(3);
	t429 = t377 * r_i_i_C(2);
	t396 = r_i_i_C(1) * t381 + t374 - t429;
	t427 = r_i_i_C(3) + qJ(5) + pkin(9);
	t435 = -(t396 * t378 - t427 * t382) * qJD(3) + t378 * qJD(5);
	t344 = t438 * t378 + t439 * t382;
	t432 = t427 * t378 + t396 * t382 + pkin(2);
	t422 = t375 * t382;
	t421 = t375 * t383;
	t417 = qJD(4) * t377;
	t416 = qJD(4) * t381;
	t415 = qJD(4) * t382;
	t412 = t430 * t383;
	t410 = t375 * t419;
	t409 = qJD(2) * t421;
	t404 = t375 * t408;
	t362 = t430 * t379 + t383 * t407;
	t349 = t362 * qJD(1) + t361 * qJD(2);
	t400 = -t344 * t377 + t349 * t381;
	t359 = t426 * t378 + t379 * t422;
	t394 = -t359 * t381 + t377 * t421;
	t363 = t412 - t402;
	t356 = -t363 * t378 + t380 * t422;
	t357 = t363 * t382 + t378 * t423;
	t390 = t361 * t378 + t382 * t413;
	t387 = -t375 * t379 * t378 + t426 * t382;
	t352 = t387 * qJD(3) + t382 * t409;
	t389 = -t352 * t377 + t381 * t410;
	t388 = qJD(4) * t434;
	t343 = t439 * t378 - t438 * t382;
	t347 = -qJD(1) * t397 - qJD(2) * t412 + (qJD(2) * t426 + qJD(1)) * t420;
	t386 = -t347 * t377 + (-t357 * t377 + t362 * t381) * qJD(4);
	t348 = t361 * qJD(1) + t362 * qJD(2);
	t342 = t356 * qJD(3) - t348 * t382 + t378 * t404;
	t339 = -t342 * t377 - t347 * t381 + (-t357 * t381 - t362 * t377) * qJD(4);
	t384 = t434 * t415 - t435;
	t351 = t359 * qJD(3) + t378 * t409;
	t341 = t357 * qJD(3) - t348 * t378 - t382 * t404;
	t340 = t342 * t381 + t386;
	t1 = [-t390 * qJD(5) - t350 * pkin(2) - t396 * t344 + (-pkin(8) - t434) * t349 - t427 * t343 + (-t430 * pkin(1) - pkin(7) * t423) * qJD(1) + (t436 * r_i_i_C(2) + t431 * t437) * qJD(4), (-t348 * t381 - t363 * t417) * r_i_i_C(2) - t348 * pkin(8) + t432 * t347 + t384 * t362 + t431 * (-t348 * t377 + t363 * t416), qJD(5) * t357 - t396 * t341 + t427 * t342 - t356 * t388, -r_i_i_C(2) * t340 + t431 * t339, t341; -t348 * pkin(2) - t347 * pkin(8) + t340 * r_i_i_C(1) + t339 * r_i_i_C(2) - t356 * qJD(5) + t342 * t374 + t427 * t341 + (-pkin(1) * t380 + pkin(7) * t413) * qJD(1) + t386 * pkin(4), (t350 * t381 - t361 * t417) * r_i_i_C(2) + t350 * pkin(8) - t432 * t349 + t384 * t360 + t431 * (t350 * t377 + t361 * t416), qJD(5) * t354 - t396 * t343 + t427 * t344 + t390 * t388, t400 * r_i_i_C(1) + (-t344 * t381 - t349 * t377) * r_i_i_C(2) + (-r_i_i_C(1) * t436 + t437 * r_i_i_C(2)) * qJD(4) + (-qJD(4) * t436 + t400) * pkin(4), t343; 0, (((t431 * t381 - t429) * qJD(4) - t432 * qJD(2)) * t379 + (qJD(2) * pkin(8) + t434 * (qJD(2) - t415) + t435) * t383) * t375, qJD(5) * t359 - t396 * t351 + t427 * t352 - t387 * t388, t389 * r_i_i_C(1) + (-t352 * t381 - t377 * t410) * r_i_i_C(2) + (t394 * r_i_i_C(1) + (t359 * t377 + t381 * t421) * r_i_i_C(2)) * qJD(4) + (t394 * qJD(4) + t389) * pkin(4), t351;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end