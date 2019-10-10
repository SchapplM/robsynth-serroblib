% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:11:34
	% EndTime: 2019-10-10 13:11:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:11:34
	% EndTime: 2019-10-10 13:11:34
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
	% StartTime: 2019-10-10 13:11:35
	% EndTime: 2019-10-10 13:11:35
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
	% StartTime: 2019-10-10 13:11:35
	% EndTime: 2019-10-10 13:11:35
	% DurationCPUTime: 0.23s
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
	% StartTime: 2019-10-10 13:11:37
	% EndTime: 2019-10-10 13:11:37
	% DurationCPUTime: 0.79s
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
	t347 = qJD(1) * t360 + qJD(2) * t359;
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
	t350 = qJD(3) * t384 + t376 * t392;
	t346 = qJD(1) * t359 + qJD(2) * t360;
	t345 = -qJD(1) * t396 - t378 * t402 + t390 * t408;
	t340 = t385 * qJD(3) - t346 * t376 + t372 * t394;
	t339 = t355 * qJD(3) - t346 * t372 - t376 * t394;
	t338 = t340 * t375 - t345 * t371 + (-t355 * t371 + t360 * t375) * qJD(4);
	t337 = -t340 * t371 - t345 * t375 + (-t355 * t375 - t360 * t371) * qJD(4);
	t1 = [t425 * r_i_i_C(1) + t426 * r_i_i_C(2) - t342 * pkin(3) - t348 * pkin(2) - t347 * pkin(9) + t419 * t380 + (t423 * r_i_i_C(1) - t422 * r_i_i_C(2)) * qJD(4) + (-t378 * pkin(1) - pkin(8) * t412) * qJD(1), (-t346 * t371 + t361 * t400) * r_i_i_C(1) + (-t346 * t375 - t361 * t401) * r_i_i_C(2) - t346 * pkin(9) + t420 * t345 + t379 * t360, -t386 * t339 + t419 * t340 - t385 * t383, t337 * r_i_i_C(1) - t338 * r_i_i_C(2), 0, 0; -t346 * pkin(2) + t340 * pkin(3) - t345 * pkin(9) + t338 * r_i_i_C(1) + t337 * r_i_i_C(2) + t419 * t339 + (-pkin(1) * t374 + pkin(8) * t409) * qJD(1), (t348 * t371 + t359 * t400) * r_i_i_C(1) + (t348 * t375 - t359 * t401) * r_i_i_C(2) + t348 * pkin(9) - t420 * t347 + t379 * t358, t419 * t342 - (-t359 * t372 - t397) * t383 + t386 * t380, -t426 * r_i_i_C(1) + t425 * r_i_i_C(2) + (t422 * r_i_i_C(1) + t423 * r_i_i_C(2)) * qJD(4), 0, 0; 0, ((-qJD(2) * t420 + t388 * qJD(4)) * t373 + (qJD(2) * pkin(9) - t424 + t387 * (qJD(2) - t399)) * t377) * t369, t419 * t350 - t384 * t383 + t386 * (-qJD(3) * t357 - t372 * t392), (-t350 * t371 + t375 * t393) * r_i_i_C(1) + (-t350 * t375 - t371 * t393) * r_i_i_C(2) + ((-t357 * t375 + t371 * t410) * r_i_i_C(1) + (t357 * t371 + t375 * t410) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:11:37
	% EndTime: 2019-10-10 13:11:38
	% DurationCPUTime: 1.04s
	% Computational Cost: add. (787->119), mult. (1843->197), div. (0->0), fcn. (1852->12), ass. (0->76)
	t412 = sin(qJ(1));
	t465 = cos(pkin(6));
	t427 = qJD(2) * t465 + qJD(1);
	t411 = sin(qJ(2));
	t442 = t412 * t465;
	t432 = t411 * t442;
	t451 = qJD(2) * t411;
	t415 = cos(qJ(2));
	t416 = cos(qJ(1));
	t456 = t416 * t415;
	t382 = -qJD(1) * t432 - t412 * t451 + t427 * t456;
	t441 = t416 * t465;
	t393 = t411 * t441 + t412 * t415;
	t410 = sin(qJ(3));
	t414 = cos(qJ(3));
	t408 = sin(pkin(6));
	t452 = qJD(1) * t408;
	t446 = t412 * t452;
	t458 = t408 * t416;
	t447 = t414 * t458;
	t376 = (-qJD(3) * t393 + t446) * t410 - qJD(3) * t447 + t382 * t414;
	t431 = t415 * t441;
	t457 = t411 * t412;
	t392 = -t431 + t457;
	t406 = qJD(4) + qJD(5);
	t438 = t392 * t406 + t376;
	t413 = cos(qJ(4));
	t403 = pkin(4) * t413 + pkin(3);
	t407 = qJ(4) + qJ(5);
	t404 = sin(t407);
	t405 = cos(t407);
	t430 = r_i_i_C(1) * t405 - r_i_i_C(2) * t404;
	t426 = t403 + t430;
	t466 = r_i_i_C(3) + pkin(11) + pkin(10);
	t474 = (t426 * t410 - t466 * t414) * qJD(3);
	t394 = t416 * t411 + t415 * t442;
	t381 = t394 * qJD(1) + t393 * qJD(2);
	t387 = -t393 * t414 + t410 * t458;
	t473 = -t387 * t406 - t381;
	t429 = t404 * r_i_i_C(1) + t405 * r_i_i_C(2);
	t409 = sin(qJ(4));
	t467 = t409 * pkin(4);
	t472 = qJD(4) * t467 + t429 * t406;
	t470 = t466 * t410 + t426 * t414 + pkin(2);
	t463 = t404 * t406;
	t462 = t405 * t406;
	t461 = t408 * t412;
	t460 = t408 * t414;
	t459 = t408 * t415;
	t450 = qJD(2) * t415;
	t379 = -qJD(1) * t431 - t416 * t450 + t427 * t457;
	t395 = -t432 + t456;
	t389 = t395 * t414 + t410 * t461;
	t437 = -t389 * t406 - t379;
	t380 = t393 * qJD(1) + t394 * qJD(2);
	t425 = -t395 * t410 + t412 * t460;
	t445 = t416 * t452;
	t374 = t425 * qJD(3) - t380 * t414 + t410 * t445;
	t440 = t394 * t406 + t374;
	t369 = -t440 * t404 + t437 * t405;
	t370 = t437 * t404 + t440 * t405;
	t455 = t369 * r_i_i_C(1) - t370 * r_i_i_C(2);
	t454 = (-t438 * t404 - t405 * t473) * r_i_i_C(1) + (t404 * t473 - t438 * t405) * r_i_i_C(2);
	t391 = t465 * t410 + t411 * t460;
	t444 = t408 * t451;
	t424 = -t391 * t406 + t444;
	t422 = -t408 * t411 * t410 + t465 * t414;
	t443 = t408 * t450;
	t384 = t422 * qJD(3) + t414 * t443;
	t428 = t406 * t459 - t384;
	t453 = (t428 * t404 + t424 * t405) * r_i_i_C(1) + (-t424 * t404 + t428 * t405) * r_i_i_C(2);
	t449 = qJD(4) * t413;
	t419 = t387 * qJD(3) - t382 * t410 + t414 * t446;
	t418 = t472 * t414 + t474;
	t373 = t389 * qJD(3) - t380 * t410 - t414 * t445;
	t1 = [-t382 * pkin(2) - t381 * pkin(9) - t376 * t403 + (-t438 * r_i_i_C(1) + r_i_i_C(2) * t473) * t405 + (r_i_i_C(1) * t473 + t438 * r_i_i_C(2)) * t404 + t466 * t419 + (-pkin(1) * t416 - pkin(8) * t461) * qJD(1) + (-t381 * t409 + (-t387 * t409 - t392 * t413) * qJD(4)) * pkin(4), (-t380 * t404 + t395 * t462) * r_i_i_C(1) + (-t380 * t405 - t395 * t463) * r_i_i_C(2) - t380 * pkin(9) + (-t380 * t409 + t395 * t449) * pkin(4) + t470 * t379 + t418 * t394, -t426 * t373 + t466 * t374 - t425 * t472, (-t374 * t409 - t379 * t413 + (-t389 * t413 - t394 * t409) * qJD(4)) * pkin(4) + t455, t455, 0; -t380 * pkin(2) - t379 * pkin(9) + t370 * r_i_i_C(1) + t369 * r_i_i_C(2) + t374 * t403 + t466 * t373 + (-pkin(1) * t412 + pkin(8) * t458) * qJD(1) + (-t379 * t409 + (-t389 * t409 + t394 * t413) * qJD(4)) * pkin(4), (t382 * t404 + t393 * t462) * r_i_i_C(1) + (t382 * t405 - t393 * t463) * r_i_i_C(2) + t382 * pkin(9) + (t382 * t409 + t393 * t449) * pkin(4) - t470 * t381 + t418 * t392, t466 * t376 - t472 * (-t393 * t410 - t447) + t426 * t419, (-t376 * t409 + t381 * t413 + (t387 * t413 - t392 * t409) * qJD(4)) * pkin(4) + t454, t454, 0; 0, ((pkin(4) * t449 - qJD(2) * t470 + t430 * t406) * t411 + (qJD(2) * pkin(9) + (-qJD(4) * t414 + qJD(2)) * t467 - t474 + t429 * (-t406 * t414 + qJD(2))) * t415) * t408, t466 * t384 - t472 * t422 + t426 * (-t391 * qJD(3) - t410 * t443), (t413 * t444 - t384 * t409 + (-t391 * t413 + t409 * t459) * qJD(4)) * pkin(4) + t453, t453, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:11:39
	% EndTime: 2019-10-10 13:11:41
	% DurationCPUTime: 1.47s
	% Computational Cost: add. (1430->149), mult. (3116->234), div. (0->0), fcn. (3218->12), ass. (0->87)
	t518 = qJ(4) + qJ(5);
	t515 = sin(t518);
	t523 = sin(qJ(1));
	t582 = cos(pkin(6));
	t543 = qJD(2) * t582 + qJD(1);
	t522 = sin(qJ(2));
	t556 = t523 * t582;
	t548 = t522 * t556;
	t570 = qJD(2) * t522;
	t526 = cos(qJ(2));
	t527 = cos(qJ(1));
	t572 = t527 * t526;
	t488 = -qJD(1) * t548 - t523 * t570 + t543 * t572;
	t555 = t527 * t582;
	t501 = t522 * t555 + t523 * t526;
	t521 = sin(qJ(3));
	t525 = cos(qJ(3));
	t519 = sin(pkin(6));
	t571 = qJD(1) * t519;
	t563 = t523 * t571;
	t574 = t519 * t527;
	t564 = t525 * t574;
	t474 = (-qJD(3) * t501 + t563) * t521 - qJD(3) * t564 + t488 * t525;
	t547 = t526 * t555;
	t573 = t523 * t522;
	t500 = -t547 + t573;
	t517 = qJD(4) + qJD(5);
	t553 = t500 * t517 + t474;
	t592 = t553 * t515;
	t516 = cos(t518);
	t591 = t553 * t516;
	t586 = r_i_i_C(1) + pkin(5);
	t583 = r_i_i_C(3) + qJ(6);
	t524 = cos(qJ(4));
	t514 = t524 * pkin(4) + pkin(3);
	t520 = sin(qJ(4));
	t585 = t520 * pkin(4);
	t542 = -qJD(4) * t585 + t515 * qJD(6);
	t568 = qJD(3) * t521;
	t584 = r_i_i_C(2) + pkin(11) + pkin(10);
	t590 = -t514 * t568 + (t584 * qJD(3) + t542) * t525;
	t589 = t525 * t514 + t584 * t521 + pkin(2);
	t588 = (qJD(2) * t525 - t517) * t522 + t526 * t568;
	t534 = t583 * t515 + t586 * t516 + t514;
	t580 = t515 * t517;
	t579 = t516 * t517;
	t578 = t517 * t525;
	t577 = t519 * t523;
	t576 = t519 * t525;
	t575 = t519 * t526;
	t569 = qJD(2) * t526;
	t567 = qJD(4) * t524;
	t566 = t516 * qJD(6);
	t565 = t515 * t575;
	t562 = t527 * t571;
	t561 = t519 * t570;
	t560 = t519 * t569;
	t502 = t527 * t522 + t526 * t556;
	t486 = t501 * qJD(1) + t502 * qJD(2);
	t503 = -t548 + t572;
	t541 = -t503 * t521 + t523 * t576;
	t472 = t541 * qJD(3) - t486 * t525 + t521 * t562;
	t554 = t502 * t517 + t472;
	t487 = t502 * qJD(1) + t501 * qJD(2);
	t494 = -t501 * t525 + t521 * t574;
	t551 = t494 * t517 + t487;
	t546 = t502 * t578 - t486;
	t545 = t500 * t578 + t488;
	t544 = (qJD(2) - t578) * t526;
	t496 = t503 * t525 + t521 * t577;
	t538 = -t519 * t522 * t521 + t582 * t525;
	t499 = t582 * t521 + t522 * t576;
	t485 = -qJD(1) * t547 - t527 * t569 + t543 * t573;
	t456 = t485 * t516 + t496 * t579 + t554 * t515;
	t457 = -t485 * t515 - t496 * t580 + t554 * t516;
	t537 = -(-t496 * t516 - t502 * t515) * qJD(6) + t583 * t457 - t586 * t456;
	t458 = -t487 * t516 - t494 * t579 + t592;
	t536 = -(t494 * t516 - t500 * t515) * qJD(6) + t583 * (t487 * t515 + t494 * t580 + t591) - t586 * t458;
	t491 = t538 * qJD(3) + t525 * t560;
	t469 = t491 * t515 + t499 * t579 - t516 * t561 - t517 * t565;
	t535 = -(-t499 * t516 + t565) * qJD(6) + t583 * (t491 * t516 - t499 * t580 + t515 * t561 - t575 * t579) - t586 * t469;
	t533 = t485 * t525 + t502 * t568 + t503 * t517;
	t532 = -t487 * t525 + t500 * t568 + t501 * t517;
	t531 = (-t586 * t515 + t583 * t516) * t517 + t542;
	t530 = t494 * qJD(3) - t488 * t521 + t525 * t563;
	t471 = t496 * qJD(3) - t486 * t521 - t525 * t562;
	t1 = [-(-t494 * t515 - t500 * t516) * qJD(6) - t474 * t514 - t488 * pkin(2) - t487 * pkin(9) + t584 * t530 + t586 * (-t551 * t515 - t591) + t583 * (t551 * t516 - t592) + (-t527 * pkin(1) - pkin(8) * t577) * qJD(1) + (-t487 * t520 + (-t494 * t520 - t500 * t524) * qJD(4)) * pkin(4), -t503 * t566 - t486 * pkin(9) + t586 * (t546 * t515 + t533 * t516) + t583 * (t533 * t515 - t546 * t516) + (-t486 * t520 + t503 * t567) * pkin(4) - t590 * t502 + t589 * t485, -t534 * t471 + t584 * t472 + t531 * t541, (-t472 * t520 - t485 * t524 + (-t496 * t524 - t502 * t520) * qJD(4)) * pkin(4) + t537, t537, t456; -(-t496 * t515 + t502 * t516) * qJD(6) + t472 * t514 - t486 * pkin(2) - t485 * pkin(9) + t584 * t471 + t586 * t457 + t583 * t456 + (-t523 * pkin(1) + pkin(8) * t574) * qJD(1) + (-t485 * t520 + (-t496 * t520 + t502 * t524) * qJD(4)) * pkin(4), -t501 * t566 + t488 * pkin(9) + t586 * (t545 * t515 + t532 * t516) + t583 * (t532 * t515 - t545 * t516) + (t488 * t520 + t501 * t567) * pkin(4) - t590 * t500 - t589 * t487, t584 * t474 + t534 * t530 + t531 * (-t501 * t521 - t564), (-t474 * t520 + t487 * t524 + (t494 * t524 - t500 * t520) * qJD(4)) * pkin(4) + t536, t536, t458; 0, (t586 * (t515 * t544 - t588 * t516) - t583 * (t588 * t515 + t516 * t544) + (pkin(4) * t567 - t566) * t522 + t590 * t526 + ((pkin(9) + t585) * t526 - t589 * t522) * qJD(2)) * t519, t584 * t491 + t534 * (-t499 * qJD(3) - t521 * t560) + t531 * t538, (t524 * t561 - t491 * t520 + (-t499 * t524 + t520 * t575) * qJD(4)) * pkin(4) + t535, t535, t469;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end