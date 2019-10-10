% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPP9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:33
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP9_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP9_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP9_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:10
	% EndTime: 2019-10-10 12:33:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:10
	% EndTime: 2019-10-10 12:33:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:11
	% EndTime: 2019-10-10 12:33:11
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
	% StartTime: 2019-10-10 12:33:12
	% EndTime: 2019-10-10 12:33:12
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
	% StartTime: 2019-10-10 12:33:13
	% EndTime: 2019-10-10 12:33:14
	% DurationCPUTime: 0.77s
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
	t342 = (-qJD(3) * t359 + t395) * t372 - qJD(3) * t397 + t348 * t376;
	t360 = t370 * t407 + t406;
	t347 = t360 * qJD(1) + t359 * qJD(2);
	t371 = sin(qJ(4));
	t375 = cos(qJ(4));
	t426 = t342 * t371 - t347 * t375;
	t425 = -t342 * t375 - t347 * t371;
	t388 = t375 * r_i_i_C(1) - t371 * r_i_i_C(2);
	t386 = pkin(3) + t388;
	t419 = r_i_i_C(3) + pkin(10);
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
	t1 = [t425 * r_i_i_C(1) + t426 * r_i_i_C(2) - t342 * pkin(3) - t348 * pkin(2) - t347 * pkin(9) + t419 * t380 + (t423 * r_i_i_C(1) - t422 * r_i_i_C(2)) * qJD(4) + (-t378 * pkin(1) - pkin(8) * t412) * qJD(1), (-t346 * t371 + t361 * t400) * r_i_i_C(1) + (-t346 * t375 - t361 * t401) * r_i_i_C(2) - t346 * pkin(9) + t420 * t345 + t379 * t360, -t386 * t339 + t419 * t340 - t385 * t383, t337 * r_i_i_C(1) - t338 * r_i_i_C(2), 0, 0; -t346 * pkin(2) + t340 * pkin(3) - t345 * pkin(9) + t338 * r_i_i_C(1) + t337 * r_i_i_C(2) + t419 * t339 + (-pkin(1) * t374 + pkin(8) * t409) * qJD(1), (t348 * t371 + t359 * t400) * r_i_i_C(1) + (t348 * t375 - t359 * t401) * r_i_i_C(2) + t348 * pkin(9) - t420 * t347 + t379 * t358, t419 * t342 - (-t359 * t372 - t397) * t383 + t386 * t380, -t426 * r_i_i_C(1) + t425 * r_i_i_C(2) + (t422 * r_i_i_C(1) + t423 * r_i_i_C(2)) * qJD(4), 0, 0; 0, ((-qJD(2) * t420 + t388 * qJD(4)) * t373 + (qJD(2) * pkin(9) - t424 + t387 * (qJD(2) - t399)) * t377) * t369, t419 * t350 - t384 * t383 + t386 * (-t357 * qJD(3) - t372 * t392), (-t350 * t371 + t375 * t393) * r_i_i_C(1) + (-t350 * t375 - t371 * t393) * r_i_i_C(2) + ((-t357 * t375 + t371 * t410) * r_i_i_C(1) + (t357 * t371 + t375 * t410) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:15
	% EndTime: 2019-10-10 12:33:16
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (793->118), mult. (2339->192), div. (0->0), fcn. (2413->10), ass. (0->74)
	t438 = cos(pkin(6));
	t442 = sin(qJ(1));
	t441 = sin(qJ(2));
	t483 = t442 * t441;
	t473 = t438 * t483;
	t477 = qJD(2) * t441;
	t445 = cos(qJ(2));
	t446 = cos(qJ(1));
	t479 = t446 * t445;
	t417 = -qJD(1) * t473 - t442 * t477 + (qJD(2) * t438 + qJD(1)) * t479;
	t480 = t446 * t441;
	t482 = t442 * t445;
	t429 = t438 * t480 + t482;
	t440 = sin(qJ(3));
	t444 = cos(qJ(3));
	t437 = sin(pkin(6));
	t478 = qJD(1) * t437;
	t470 = t442 * t478;
	t486 = t437 * t446;
	t472 = t444 * t486;
	t406 = (-qJD(3) * t429 + t470) * t440 - qJD(3) * t472 + t417 * t444;
	t430 = t438 * t482 + t480;
	t416 = t430 * qJD(1) + t429 * qJD(2);
	t439 = sin(qJ(4));
	t443 = cos(qJ(4));
	t423 = -t429 * t444 + t440 * t486;
	t428 = -t438 * t479 + t483;
	t500 = t423 * t443 - t428 * t439;
	t505 = t500 * qJD(4) - t406 * t439 + t416 * t443;
	t501 = t423 * t439 + t428 * t443;
	t504 = t501 * qJD(4) + t406 * t443 + t416 * t439;
	t497 = r_i_i_C(1) + pkin(10);
	t499 = t444 * pkin(3) + t497 * t440 + pkin(2);
	t498 = -pkin(3) * t440 + t497 * t444;
	t494 = r_i_i_C(3) + qJ(5);
	t496 = r_i_i_C(2) - pkin(4);
	t452 = t494 * t439 - t496 * t443 + pkin(3);
	t456 = t473 - t479;
	t488 = t437 * t442;
	t425 = t440 * t488 - t444 * t456;
	t491 = t425 * t439;
	t487 = t437 * t444;
	t485 = t439 * t444;
	t484 = t439 * t445;
	t481 = t443 * t445;
	t476 = qJD(3) * t440;
	t475 = qJD(3) * t445;
	t474 = qJD(4) * t444;
	t469 = t446 * t478;
	t468 = t437 * t477;
	t467 = qJD(2) * t437 * t445;
	t465 = -qJD(2) + t474;
	t415 = t429 * qJD(1) + t430 * qJD(2);
	t464 = t430 * t474 - t415;
	t463 = t428 * t474 + t417;
	t460 = t425 * t443 + t430 * t439;
	t427 = t438 * t440 + t441 * t487;
	t459 = -t427 * t443 + t437 * t484;
	t458 = t440 * t456 + t442 * t487;
	t457 = -t437 * t441 * t440 + t438 * t444;
	t453 = qJD(3) * t498;
	t414 = t428 * qJD(1) + t456 * qJD(2);
	t451 = -qJD(4) * t456 + t414 * t444 + t430 * t476;
	t450 = qJD(4) * t429 - t416 * t444 + t428 * t476;
	t449 = t440 * t475 + (qJD(2) * t444 - qJD(4)) * t441;
	t448 = qJD(5) * t439 + (t496 * t439 + t494 * t443) * qJD(4);
	t447 = t423 * qJD(3) - t417 * t440 + t444 * t470;
	t420 = t457 * qJD(3) + t444 * t467;
	t409 = -t459 * qJD(4) + t420 * t439 - t443 * t468;
	t404 = t458 * qJD(3) - t415 * t444 + t440 * t469;
	t403 = t425 * qJD(3) - t415 * t440 - t444 * t469;
	t394 = -t414 * t439 - qJD(4) * t491 + (qJD(4) * t430 + t404) * t443;
	t393 = t460 * qJD(4) + t404 * t439 + t414 * t443;
	t1 = [t501 * qJD(5) - t406 * pkin(3) - t417 * pkin(2) - t416 * pkin(9) + t497 * t447 + t496 * t504 + t494 * t505 + (-t446 * pkin(1) - pkin(8) * t488) * qJD(1), -(t430 * t485 - t443 * t456) * qJD(5) - t415 * pkin(9) - t496 * (t464 * t439 + t451 * t443) + t494 * (t451 * t439 - t464 * t443) - t430 * t453 + t499 * t414, -t452 * t403 + t497 * t404 + t448 * t458, t460 * qJD(5) + t496 * t393 + t494 * t394, t393, 0; -(t430 * t443 - t491) * qJD(5) + t404 * pkin(3) - t415 * pkin(2) - t414 * pkin(9) + t497 * t403 - t496 * t394 + t494 * t393 + (-t442 * pkin(1) + pkin(8) * t486) * qJD(1), -(t428 * t485 + t429 * t443) * qJD(5) + t417 * pkin(9) - t496 * (t463 * t439 + t450 * t443) + t494 * (t450 * t439 - t463 * t443) - t428 * t453 - t499 * t416, t497 * t406 + t448 * (-t429 * t440 - t472) + t452 * t447, -t500 * qJD(5) + t494 * t504 - t496 * t505, -t505, 0; 0, (t496 * (t449 * t443 + t465 * t484) - t494 * (t449 * t439 - t465 * t481) - (t441 * t443 - t444 * t484) * qJD(5) + t498 * t475 + (t445 * pkin(9) - t441 * t499) * qJD(2)) * t437, t497 * t420 + t448 * t457 + t452 * (-t427 * qJD(3) - t440 * t467), -t459 * qJD(5) + t494 * (t439 * t468 + t420 * t443 + (-t427 * t439 - t437 * t481) * qJD(4)) + t496 * t409, t409, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:15
	% EndTime: 2019-10-10 12:33:16
	% DurationCPUTime: 1.42s
	% Computational Cost: add. (1056->135), mult. (3096->219), div. (0->0), fcn. (3220->10), ass. (0->80)
	t443 = sin(qJ(1));
	t442 = sin(qJ(2));
	t503 = cos(pkin(6));
	t470 = t443 * t503;
	t467 = t442 * t470;
	t487 = qJD(2) * t442;
	t446 = cos(qJ(2));
	t447 = cos(qJ(1));
	t489 = t447 * t446;
	t415 = -qJD(1) * t467 - t443 * t487 + (qJD(2) * t503 + qJD(1)) * t489;
	t469 = t447 * t503;
	t427 = t442 * t469 + t443 * t446;
	t441 = sin(qJ(3));
	t445 = cos(qJ(3));
	t439 = sin(pkin(6));
	t488 = qJD(1) * t439;
	t476 = t443 * t488;
	t494 = t439 * t447;
	t478 = t445 * t494;
	t401 = (-qJD(3) * t427 + t476) * t441 - qJD(3) * t478 + t415 * t445;
	t428 = t447 * t442 + t446 * t470;
	t414 = t428 * qJD(1) + t427 * qJD(2);
	t440 = sin(qJ(4));
	t444 = cos(qJ(4));
	t421 = -t427 * t445 + t441 * t494;
	t426 = t443 * t442 - t446 * t469;
	t462 = t421 * t444 - t426 * t440;
	t511 = t462 * qJD(4) - t401 * t440 + t414 * t444;
	t508 = -t421 * t440 - t426 * t444;
	t509 = qJD(4) * t508 - t401 * t444 - t414 * t440;
	t480 = pkin(5) + pkin(10) + r_i_i_C(1);
	t507 = t445 * pkin(3) + t480 * t441 + pkin(2);
	t506 = -pkin(3) * t441 + t480 * t445;
	t479 = -r_i_i_C(3) - qJ(6) - pkin(4);
	t504 = r_i_i_C(2) + qJ(5);
	t450 = t504 * t440 - t479 * t444 + pkin(3);
	t498 = t439 * t442;
	t497 = t439 * t443;
	t496 = t439 * t445;
	t495 = t439 * t446;
	t493 = t440 * t445;
	t492 = t444 * t445;
	t491 = t444 * t446;
	t490 = t445 * t446;
	t486 = qJD(2) * t445;
	t485 = qJD(3) * t441;
	t484 = qJD(3) * t446;
	t483 = qJD(4) * t440;
	t482 = qJD(4) * t444;
	t481 = qJD(4) * t445;
	t475 = t447 * t488;
	t474 = t439 * t487;
	t473 = qJD(2) * t495;
	t472 = t446 * t483;
	t471 = t441 * t484;
	t413 = t427 * qJD(1) + t428 * qJD(2);
	t466 = t428 * t481 - t413;
	t465 = t426 * t481 + t415;
	t454 = t467 - t489;
	t423 = t441 * t497 - t445 * t454;
	t460 = t423 * t444 + t428 * t440;
	t425 = t503 * t441 + t442 * t496;
	t459 = t425 * t440 + t439 * t491;
	t458 = t441 * t454 + t443 * t496;
	t455 = -t441 * t498 + t503 * t445;
	t453 = qJD(3) * t506;
	t412 = t426 * qJD(1) + t454 * qJD(2);
	t452 = -qJD(4) * t454 + t412 * t445 + t428 * t485;
	t451 = qJD(4) * t427 - t414 * t445 + t426 * t485;
	t449 = t421 * qJD(3) - t415 * t441 + t445 * t476;
	t448 = t440 * qJD(5) + t444 * qJD(6) + (t479 * t440 + t504 * t444) * qJD(4);
	t418 = t455 * qJD(3) + t445 * t473;
	t408 = t423 * t440 - t428 * t444;
	t405 = -t459 * qJD(4) + t418 * t444 + t440 * t474;
	t404 = t418 * t440 + t425 * t482 - t439 * t472 - t444 * t474;
	t399 = t458 * qJD(3) - t413 * t445 + t441 * t475;
	t398 = t423 * qJD(3) - t413 * t441 - t445 * t475;
	t389 = -t412 * t440 - t423 * t483 + (qJD(4) * t428 + t399) * t444;
	t388 = t460 * qJD(4) + t399 * t440 + t412 * t444;
	t1 = [t462 * qJD(6) - t508 * qJD(5) - t401 * pkin(3) - t415 * pkin(2) - t414 * pkin(9) + t504 * t511 + (-pkin(1) * t447 - pkin(8) * t497) * qJD(1) + t480 * t449 - t479 * t509, -(t428 * t492 + t440 * t454) * qJD(6) - (t428 * t493 - t444 * t454) * qJD(5) - t413 * pkin(9) + t504 * (t452 * t440 - t466 * t444) - t479 * (t466 * t440 + t452 * t444) - t428 * t453 + t507 * t412, -t450 * t398 + t480 * t399 + t448 * t458, qJD(5) * t460 - t408 * qJD(6) + t479 * t388 + t504 * t389, t388, t389; -t413 * pkin(2) + t399 * pkin(3) - t412 * pkin(9) + t408 * qJD(5) + t460 * qJD(6) + t504 * t388 + (-pkin(1) * t443 + pkin(8) * t494) * qJD(1) + t480 * t398 - t479 * t389, -(t426 * t492 - t427 * t440) * qJD(6) - (t426 * t493 + t427 * t444) * qJD(5) + t415 * pkin(9) + t504 * (t451 * t440 - t465 * t444) - t479 * (t465 * t440 + t451 * t444) - t426 * t453 - t507 * t414, t480 * t401 + t450 * t449 + t448 * (-t427 * t441 - t478), -qJD(5) * t462 - qJD(6) * t508 - t479 * t511 - t504 * t509, -t511, -t509; 0, t479 * (-t440 * t473 - t482 * t498) + (-t504 * ((qJD(2) - t481) * t491 + (t471 + (-qJD(4) + t486) * t442) * t440) + t479 * (t445 * t472 + (t442 * t486 + t471) * t444) - (-t440 * t442 - t444 * t490) * qJD(6) - (-t440 * t490 + t442 * t444) * qJD(5) + t506 * t484 + (t446 * pkin(9) - t442 * t507) * qJD(2)) * t439, t480 * t418 + t450 * (-t425 * qJD(3) - t441 * t473) + t448 * t455, -t459 * qJD(6) - (-t425 * t444 + t440 * t495) * qJD(5) + t504 * t405 + t479 * t404, t404, t405;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end