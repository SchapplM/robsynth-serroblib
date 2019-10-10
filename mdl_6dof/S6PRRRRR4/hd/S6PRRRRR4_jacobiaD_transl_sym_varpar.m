% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRRR4
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:19
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(13));
	t50 = sin(pkin(13));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:24
	% EndTime: 2019-10-09 23:19:24
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (95->48), mult. (333->96), div. (0->0), fcn. (328->10), ass. (0->36)
	t218 = sin(pkin(7));
	t248 = (pkin(9) + r_i_i_C(3)) * t218;
	t217 = sin(pkin(13));
	t220 = cos(pkin(13));
	t224 = sin(qJ(2));
	t222 = cos(pkin(6));
	t226 = cos(qJ(2));
	t241 = t222 * t226;
	t211 = -t217 * t224 + t220 * t241;
	t219 = sin(pkin(6));
	t245 = t218 * t219;
	t221 = cos(pkin(7));
	t223 = sin(qJ(3));
	t244 = t221 * t223;
	t225 = cos(qJ(3));
	t243 = t221 * t225;
	t242 = t222 * t224;
	t240 = t223 * t224;
	t239 = t223 * t226;
	t238 = t224 * t225;
	t237 = t225 * t226;
	t236 = t223 * t245;
	t235 = t225 * t245;
	t233 = r_i_i_C(1) * t223 + r_i_i_C(2) * t225;
	t232 = r_i_i_C(1) * t225 - r_i_i_C(2) * t223 + pkin(2);
	t212 = t217 * t226 + t220 * t242;
	t231 = t217 * t241 + t220 * t224;
	t230 = t217 * t242 - t220 * t226;
	t229 = t233 * t221 - t248;
	t228 = (-t221 * t238 - t239) * r_i_i_C(1) + (t221 * t240 - t237) * r_i_i_C(2);
	t227 = (-t221 * t239 - t238) * r_i_i_C(1) + (-t221 * t237 + t240) * r_i_i_C(2);
	t210 = t230 * qJD(2);
	t209 = t231 * qJD(2);
	t208 = t212 * qJD(2);
	t207 = t211 * qJD(2);
	t1 = [0, t232 * t210 + t229 * t209 + ((t223 * t231 + t230 * t243) * r_i_i_C(1) + (t225 * t231 - t230 * t244) * r_i_i_C(2)) * qJD(3), (t209 * t223 + t210 * t243) * r_i_i_C(1) + (t209 * t225 - t210 * t244) * r_i_i_C(2) + ((-t217 * t236 + t225 * t230 + t231 * t244) * r_i_i_C(1) + (-t217 * t235 - t223 * t230 + t231 * t243) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, -t232 * t208 - t229 * t207 + ((-t211 * t223 - t212 * t243) * r_i_i_C(1) + (-t211 * t225 + t212 * t244) * r_i_i_C(2)) * qJD(3), (-t207 * t223 - t208 * t243) * r_i_i_C(1) + (-t207 * t225 + t208 * t244) * r_i_i_C(2) + ((-t211 * t244 - t212 * t225 + t220 * t236) * r_i_i_C(1) + (-t211 * t243 + t212 * t223 + t220 * t235) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, (t228 * qJD(3) + (-t224 * pkin(2) + t226 * t248 + t227) * qJD(2)) * t219, -t233 * t222 * t218 * qJD(3) + (t228 * qJD(2) + t227 * qJD(3)) * t219, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:25
	% EndTime: 2019-10-09 23:19:26
	% DurationCPUTime: 0.84s
	% Computational Cost: add. (380->110), mult. (1266->209), div. (0->0), fcn. (1334->12), ass. (0->67)
	t400 = sin(qJ(3));
	t403 = cos(qJ(3));
	t393 = sin(pkin(13));
	t396 = cos(pkin(13));
	t404 = cos(qJ(2));
	t398 = cos(pkin(6));
	t401 = sin(qJ(2));
	t429 = t398 * t401;
	t412 = t393 * t429 - t396 * t404;
	t397 = cos(pkin(7));
	t428 = t398 * t404;
	t413 = t393 * t428 + t396 * t401;
	t394 = sin(pkin(7));
	t395 = sin(pkin(6));
	t436 = t394 * t395;
	t414 = t393 * t436 - t397 * t413;
	t372 = t414 * t400 - t403 * t412;
	t388 = t393 * t404 + t396 * t429;
	t419 = t396 * t428;
	t387 = -t393 * t401 + t419;
	t415 = -t387 * t397 + t396 * t436;
	t440 = -t388 * t403 + t415 * t400;
	t439 = r_i_i_C(3) + pkin(10);
	t435 = t394 * t398;
	t399 = sin(qJ(4));
	t434 = t394 * t399;
	t402 = cos(qJ(4));
	t433 = t394 * t402;
	t432 = t395 * t397;
	t431 = t397 * t400;
	t430 = t397 * t403;
	t427 = t400 * t401;
	t426 = t400 * t404;
	t425 = t401 * t403;
	t424 = t403 * t404;
	t423 = qJD(2) * t401;
	t422 = qJD(2) * t404;
	t421 = qJD(4) * t399;
	t420 = qJD(4) * t402;
	t418 = qJD(3) * t435;
	t417 = t423 * t436;
	t416 = t402 * r_i_i_C(1) - t399 * r_i_i_C(2) + pkin(3);
	t375 = t387 * t403 - t388 * t431;
	t376 = -t403 * t413 + t412 * t431;
	t411 = t397 * t424 - t427;
	t410 = -t397 * t425 - t426;
	t409 = t397 * t426 + t425;
	t408 = t397 * t427 - t424;
	t407 = qJD(4) * (-t399 * r_i_i_C(1) - t402 * r_i_i_C(2));
	t406 = -t388 * t400 - t415 * t403;
	t405 = t400 * t412 + t414 * t403;
	t386 = t398 * t397 - t404 * t436;
	t385 = t412 * qJD(2);
	t384 = t413 * qJD(2);
	t383 = t388 * qJD(2);
	t382 = -qJD(2) * t419 + t393 * t423;
	t381 = t408 * t395;
	t380 = t393 * t432 + t394 * t413;
	t379 = -t387 * t394 - t396 * t432;
	t378 = t409 * t395 + t400 * t435;
	t374 = (-t409 * qJD(2) + t410 * qJD(3)) * t395;
	t368 = t403 * t418 + (-t408 * qJD(2) + t411 * qJD(3)) * t395;
	t366 = t384 * t431 + t385 * t403 + (t400 * t413 + t412 * t430) * qJD(3);
	t364 = t382 * t431 - t383 * t403 + (-t387 * t400 - t388 * t430) * qJD(3);
	t362 = t405 * qJD(3) - t384 * t403 + t385 * t431;
	t360 = t406 * qJD(3) - t382 * t403 - t383 * t431;
	t1 = [0, (t366 * t402 - t376 * t421) * r_i_i_C(1) + (-t366 * t399 - t376 * t420) * r_i_i_C(2) + t366 * pkin(3) + t385 * pkin(2) + t439 * (t376 * qJD(3) - t384 * t430 + t385 * t400) + ((-t384 * t399 - t412 * t420) * r_i_i_C(1) + (-t384 * t402 + t412 * t421) * r_i_i_C(2) - t384 * pkin(9)) * t394, t439 * t362 + t405 * t407 + t416 * (-t372 * qJD(3) + t384 * t400 + t385 * t430), (-t362 * t399 - t385 * t433) * r_i_i_C(1) + (-t362 * t402 + t385 * t434) * r_i_i_C(2) + ((-t372 * t402 - t380 * t399) * r_i_i_C(1) + (t372 * t399 - t380 * t402) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (t364 * t402 - t375 * t421) * r_i_i_C(1) + (-t364 * t399 - t375 * t420) * r_i_i_C(2) + t364 * pkin(3) - t383 * pkin(2) + t439 * (t375 * qJD(3) - t382 * t430 - t383 * t400) + ((-t382 * t399 + t388 * t420) * r_i_i_C(1) + (-t382 * t402 - t388 * t421) * r_i_i_C(2) - t382 * pkin(9)) * t394, t439 * t360 + t406 * t407 + t416 * (t440 * qJD(3) + t382 * t400 - t383 * t430), (-t360 * t399 + t383 * t433) * r_i_i_C(1) + (-t360 * t402 - t383 * t434) * r_i_i_C(2) + ((-t379 * t399 + t402 * t440) * r_i_i_C(1) + (-t379 * t402 - t399 * t440) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (t374 * t402 + t381 * t421) * r_i_i_C(1) + (-t374 * t399 + t381 * t420) * r_i_i_C(2) + t374 * pkin(3) + (-t439 * (-t411 * qJD(2) + t408 * qJD(3)) - pkin(2) * t423 + ((t399 * t422 + t401 * t420) * r_i_i_C(1) + (-t401 * t421 + t402 * t422) * r_i_i_C(2) + pkin(9) * t422) * t394) * t395, t439 * t368 + (t411 * t395 + t403 * t435) * t407 + t416 * (-t400 * t418 + (t410 * qJD(2) - t409 * qJD(3)) * t395), (-t368 * t399 + t402 * t417) * r_i_i_C(1) + (-t368 * t402 - t399 * t417) * r_i_i_C(2) + ((-t378 * t402 - t386 * t399) * r_i_i_C(1) + (t378 * t399 - t386 * t402) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:26
	% EndTime: 2019-10-09 23:19:27
	% DurationCPUTime: 1.04s
	% Computational Cost: add. (695->129), mult. (1858->229), div. (0->0), fcn. (1972->14), ass. (0->83)
	t442 = sin(qJ(3));
	t445 = cos(qJ(3));
	t435 = sin(pkin(13));
	t438 = cos(pkin(13));
	t446 = cos(qJ(2));
	t440 = cos(pkin(6));
	t443 = sin(qJ(2));
	t481 = t440 * t443;
	t456 = t435 * t481 - t438 * t446;
	t439 = cos(pkin(7));
	t480 = t440 * t446;
	t457 = t435 * t480 + t438 * t443;
	t436 = sin(pkin(7));
	t437 = sin(pkin(6));
	t487 = t436 * t437;
	t458 = t435 * t487 - t439 * t457;
	t409 = t442 * t458 - t445 * t456;
	t425 = t435 * t446 + t438 * t481;
	t468 = t438 * t480;
	t424 = -t435 * t443 + t468;
	t459 = -t424 * t439 + t438 * t487;
	t494 = -t425 * t445 + t442 * t459;
	t493 = r_i_i_C(3) + pkin(11) + pkin(10);
	t434 = qJ(4) + qJ(5);
	t431 = sin(t434);
	t433 = qJD(4) + qJD(5);
	t490 = t431 * t433;
	t432 = cos(t434);
	t489 = t432 * t433;
	t488 = t433 * t443;
	t486 = t436 * t440;
	t444 = cos(qJ(4));
	t485 = t436 * t444;
	t484 = t437 * t439;
	t483 = t439 * t442;
	t482 = t439 * t445;
	t479 = t442 * t443;
	t478 = t442 * t446;
	t477 = t443 * t445;
	t476 = t445 * t446;
	t420 = t425 * qJD(2);
	t462 = -t420 * t436 - t433 * t494;
	t472 = qJD(2) * t443;
	t419 = -qJD(2) * t468 + t435 * t472;
	t449 = -t425 * t442 - t445 * t459;
	t397 = qJD(3) * t449 - t419 * t445 - t420 * t483;
	t416 = -t424 * t436 - t438 * t484;
	t466 = -t416 * t433 - t397;
	t475 = (t431 * t466 - t462 * t432) * r_i_i_C(1) + (t431 * t462 + t432 * t466) * r_i_i_C(2);
	t422 = t456 * qJD(2);
	t461 = t409 * t433 + t422 * t436;
	t421 = t457 * qJD(2);
	t448 = t442 * t456 + t445 * t458;
	t399 = qJD(3) * t448 - t421 * t445 + t422 * t483;
	t417 = t435 * t484 + t436 * t457;
	t465 = -t417 * t433 - t399;
	t474 = (t431 * t465 - t432 * t461) * r_i_i_C(1) + (t431 * t461 + t432 * t465) * r_i_i_C(2);
	t453 = t439 * t478 + t477;
	t415 = t437 * t453 + t442 * t486;
	t463 = t472 * t487;
	t451 = -t415 * t433 + t463;
	t452 = t439 * t479 - t476;
	t455 = t439 * t476 - t479;
	t467 = qJD(3) * t486;
	t405 = t445 * t467 + (-qJD(2) * t452 + qJD(3) * t455) * t437;
	t423 = t440 * t439 - t446 * t487;
	t464 = -t423 * t433 - t405;
	t473 = (t431 * t464 + t432 * t451) * r_i_i_C(1) + (-t431 * t451 + t432 * t464) * r_i_i_C(2);
	t471 = qJD(2) * t446;
	t470 = qJD(4) * t444;
	t441 = sin(qJ(4));
	t469 = qJD(4) * t441 * pkin(4);
	t430 = t444 * pkin(4) + pkin(3);
	t460 = t432 * r_i_i_C(1) - t431 * r_i_i_C(2) + t430;
	t412 = t424 * t445 - t425 * t483;
	t413 = -t445 * t457 + t456 * t483;
	t454 = -t439 * t477 - t478;
	t450 = -t469 + (-t431 * r_i_i_C(1) - t432 * r_i_i_C(2)) * t433;
	t418 = t452 * t437;
	t411 = (-qJD(2) * t453 + qJD(3) * t454) * t437;
	t403 = t421 * t483 + t422 * t445 + (t442 * t457 + t456 * t482) * qJD(3);
	t401 = t419 * t483 - t420 * t445 + (-t424 * t442 - t425 * t482) * qJD(3);
	t1 = [0, (t403 * t432 - t413 * t490) * r_i_i_C(1) + (-t403 * t431 - t413 * t489) * r_i_i_C(2) + t403 * t430 - t413 * t469 + t422 * pkin(2) + t493 * (qJD(3) * t413 - t421 * t482 + t422 * t442) + ((-t421 * t431 - t456 * t489) * r_i_i_C(1) + (-t421 * t432 + t456 * t490) * r_i_i_C(2) - t421 * pkin(9) + (-t421 * t441 - t456 * t470) * pkin(4)) * t436, t493 * t399 + t450 * t448 + t460 * (-t409 * qJD(3) + t421 * t442 + t422 * t482), (-t422 * t485 - t399 * t441 + (-t409 * t444 - t417 * t441) * qJD(4)) * pkin(4) + t474, t474, 0; 0, (t401 * t432 - t412 * t490) * r_i_i_C(1) + (-t401 * t431 - t412 * t489) * r_i_i_C(2) + t401 * t430 - t412 * t469 - t420 * pkin(2) + t493 * (qJD(3) * t412 - t419 * t482 - t420 * t442) + ((-t419 * t431 + t425 * t489) * r_i_i_C(1) + (-t419 * t432 - t425 * t490) * r_i_i_C(2) - t419 * pkin(9) + (-t419 * t441 + t425 * t470) * pkin(4)) * t436, t493 * t397 + t450 * t449 + t460 * (t494 * qJD(3) + t419 * t442 - t420 * t482), (t420 * t485 - t397 * t441 + (-t416 * t441 + t444 * t494) * qJD(4)) * pkin(4) + t475, t475, 0; 0, (t411 * t432 + t418 * t490) * r_i_i_C(1) + (-t411 * t431 + t418 * t489) * r_i_i_C(2) + t411 * t430 + t418 * t469 + (-t493 * (-qJD(2) * t455 + qJD(3) * t452) - pkin(2) * t472 + ((t431 * t471 + t432 * t488) * r_i_i_C(1) + (-t431 * t488 + t432 * t471) * r_i_i_C(2) + pkin(9) * t471 + (t441 * t471 + t443 * t470) * pkin(4)) * t436) * t437, t493 * t405 + t450 * (t437 * t455 + t445 * t486) + t460 * (-t442 * t467 + (qJD(2) * t454 - qJD(3) * t453) * t437), (t444 * t463 - t405 * t441 + (-t415 * t444 - t423 * t441) * qJD(4)) * pkin(4) + t473, t473, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:30
	% EndTime: 2019-10-09 23:19:32
	% DurationCPUTime: 2.00s
	% Computational Cost: add. (1934->206), mult. (4907->360), div. (0->0), fcn. (5425->16), ass. (0->123)
	t694 = r_i_i_C(3) + pkin(12);
	t623 = cos(qJ(6));
	t674 = qJD(6) * t623;
	t619 = sin(qJ(6));
	t675 = qJD(6) * t619;
	t695 = -r_i_i_C(1) * t675 - t674 * r_i_i_C(2);
	t645 = r_i_i_C(1) * t623 - r_i_i_C(2) * t619 + pkin(5);
	t622 = sin(qJ(2));
	t625 = cos(qJ(2));
	t617 = cos(pkin(13));
	t690 = cos(pkin(6));
	t661 = t617 * t690;
	t689 = sin(pkin(13));
	t598 = -t689 * t622 + t625 * t661;
	t693 = cos(qJ(3));
	t691 = pkin(4) * qJD(4);
	t593 = t598 * qJD(2);
	t615 = sin(pkin(7));
	t688 = t593 * t615;
	t649 = t690 * t689;
	t634 = t617 * t622 + t625 * t649;
	t595 = t634 * qJD(2);
	t687 = t595 * t615;
	t614 = qJ(4) + qJ(5);
	t611 = sin(t614);
	t613 = qJD(4) + qJD(5);
	t686 = t611 * t613;
	t685 = t611 * t615;
	t684 = t613 * t615;
	t616 = sin(pkin(6));
	t683 = t615 * t616;
	t620 = sin(qJ(4));
	t682 = t615 * t620;
	t624 = cos(qJ(4));
	t681 = t615 * t624;
	t680 = t615 * t625;
	t679 = t616 * t617;
	t618 = cos(pkin(7));
	t621 = sin(qJ(3));
	t678 = t618 * t621;
	t677 = t621 * t622;
	t676 = t621 * t625;
	t673 = t620 * t691;
	t671 = t615 * t679;
	t670 = t622 * t683;
	t669 = t616 * t677;
	t636 = -t622 * t661 - t689 * t625;
	t668 = t636 * t693;
	t667 = t618 * t693;
	t666 = t693 * t622;
	t665 = t693 * t625;
	t664 = qJD(2) * t683;
	t663 = t615 * t690;
	t662 = t616 * t689;
	t594 = t636 * qJD(2);
	t632 = t598 * t667 + t621 * t636 - t693 * t671;
	t550 = t632 * qJD(3) + t593 * t693 + t594 * t678;
	t583 = -t598 * t615 - t618 * t679;
	t659 = t583 * t613 + t550;
	t635 = -t617 * t625 + t622 * t649;
	t596 = t635 * qJD(2);
	t652 = t615 * t662;
	t631 = t621 * t635 - t634 * t667 + t693 * t652;
	t552 = t631 * qJD(3) - t595 * t693 + t596 * t678;
	t584 = t615 * t634 + t618 * t662;
	t658 = t584 * t613 + t552;
	t640 = -t618 * t677 + t665;
	t644 = t693 * t663;
	t656 = t618 * t665;
	t566 = qJD(3) * t644 + ((t656 - t677) * qJD(3) + t640 * qJD(2)) * t616;
	t597 = -t616 * t680 + t690 * t618;
	t657 = t597 * t613 + t566;
	t655 = t622 * t664;
	t653 = t621 * t663;
	t642 = -t598 * t621 + t636 * t667;
	t558 = t642 * qJD(3) - t593 * t678 + t594 * t693;
	t651 = -t636 * t684 + t558;
	t641 = t621 * t634 + t635 * t667;
	t560 = t641 * qJD(3) + t595 * t678 + t596 * t693;
	t650 = -t635 * t684 + t560;
	t648 = t616 * t656;
	t576 = t598 * t693 + t636 * t678;
	t647 = t576 * t613 - t688;
	t578 = -t634 * t693 + t635 * t678;
	t646 = t578 * t613 + t687;
	t638 = t618 * t666 + t676;
	t639 = t618 * t676 + t666;
	t574 = (-t639 * qJD(2) - t638 * qJD(3)) * t616;
	t643 = t613 * t670 + t574;
	t592 = t640 * t616;
	t637 = -t592 * t613 + t625 * t664;
	t610 = pkin(4) * t624 + pkin(3);
	t612 = cos(t614);
	t633 = -t694 * t611 - t645 * t612 - t610;
	t572 = -t635 * t693 + (-t618 * t634 + t652) * t621;
	t570 = -t668 + (t598 * t618 - t671) * t621;
	t535 = -t570 * t686 - t594 * t685 + t659 * t612;
	t630 = t695 * (-t570 * t611 + t583 * t612) + t694 * t535 + t645 * ((-t570 * t613 - t594 * t615) * t612 - t659 * t611);
	t537 = -t572 * t686 - t596 * t685 + t658 * t612;
	t629 = t695 * (-t572 * t611 + t584 * t612) + t694 * t537 + t645 * ((-t572 * t613 - t596 * t615) * t612 - t658 * t611);
	t582 = t639 * t616 + t653;
	t546 = -t582 * t686 + t611 * t655 + t657 * t612;
	t628 = t695 * (-t582 * t611 + t597 * t612) + t694 * t546 + t645 * ((-t582 * t613 + t655) * t612 - t657 * t611);
	t627 = t673 + (t619 * r_i_i_C(1) + t623 * r_i_i_C(2)) * t612 * qJD(6) + (t645 * t611 - t694 * t612) * t613;
	t626 = -pkin(11) - pkin(10);
	t591 = t638 * t616;
	t581 = -t644 - t648 + t669;
	t579 = t592 * t612 + t611 * t670;
	t573 = -qJD(2) * t648 - t616 * qJD(3) * t665 + (qJD(3) * t618 + qJD(2)) * t669;
	t568 = t582 * t612 + t597 * t611;
	t565 = qJD(3) * t653 + (t638 * qJD(2) + t639 * qJD(3)) * t616;
	t562 = t578 * t612 - t635 * t685;
	t561 = t576 * t612 - t636 * t685;
	t559 = t578 * qJD(3) - t595 * t667 + t596 * t621;
	t557 = t576 * qJD(3) + t593 * t667 + t594 * t621;
	t556 = t572 * t612 + t584 * t611;
	t554 = t570 * t612 + t583 * t611;
	t551 = t572 * qJD(3) - t595 * t621 - t596 * t667;
	t549 = t593 * t621 - t594 * t667 + (t598 * t678 - t621 * t671 - t668) * qJD(3);
	t548 = t637 * t611 + t643 * t612;
	t541 = -t646 * t611 + t650 * t612;
	t539 = -t647 * t611 + t651 * t612;
	t1 = [0, (t541 * t623 + t559 * t619) * r_i_i_C(1) + (-t541 * t619 + t559 * t623) * r_i_i_C(2) + t541 * pkin(5) + t560 * t610 - t559 * t626 + t596 * pkin(2) - pkin(9) * t687 + t694 * (t650 * t611 + t646 * t612) + ((-t562 * t619 - t623 * t641) * r_i_i_C(1) + (-t562 * t623 + t619 * t641) * r_i_i_C(2)) * qJD(6) + (-t595 * t682 + (-t578 * t620 - t635 * t681) * qJD(4)) * pkin(4), (t552 * t619 + t572 * t674) * r_i_i_C(1) + (t552 * t623 - t572 * t675) * r_i_i_C(2) - t552 * t626 + t633 * t551 - t627 * t631, (-t596 * t681 - t552 * t620 + (-t572 * t624 - t584 * t620) * qJD(4)) * pkin(4) + t629, t629, (-t537 * t619 + t551 * t623) * r_i_i_C(1) + (-t537 * t623 - t551 * t619) * r_i_i_C(2) + ((-t556 * t623 + t619 * t631) * r_i_i_C(1) + (t556 * t619 + t623 * t631) * r_i_i_C(2)) * qJD(6); 0, (t539 * t623 + t557 * t619) * r_i_i_C(1) + (-t539 * t619 + t557 * t623) * r_i_i_C(2) + t539 * pkin(5) + t558 * t610 - t557 * t626 + t594 * pkin(2) + pkin(9) * t688 + t694 * (t651 * t611 + t647 * t612) + ((-t561 * t619 - t623 * t642) * r_i_i_C(1) + (-t561 * t623 + t619 * t642) * r_i_i_C(2)) * qJD(6) + (t593 * t682 + (-t576 * t620 - t636 * t681) * qJD(4)) * pkin(4), (t550 * t619 + t570 * t674) * r_i_i_C(1) + (t550 * t623 - t570 * t675) * r_i_i_C(2) - t550 * t626 + t633 * t549 - t627 * t632, (-t594 * t681 - t550 * t620 + (-t570 * t624 - t583 * t620) * qJD(4)) * pkin(4) + t630, t630, (-t535 * t619 + t549 * t623) * r_i_i_C(1) + (-t535 * t623 - t549 * t619) * r_i_i_C(2) + ((-t554 * t623 + t619 * t632) * r_i_i_C(1) + (t554 * t619 + t623 * t632) * r_i_i_C(2)) * qJD(6); 0, (t548 * t623 - t573 * t619) * r_i_i_C(1) + (-t548 * t619 - t573 * t623) * r_i_i_C(2) + t548 * pkin(5) + t574 * t610 - t592 * t673 + t573 * t626 + t694 * (t643 * t611 - t637 * t612) + ((-t579 * t619 + t591 * t623) * r_i_i_C(1) + (-t579 * t623 - t591 * t619) * r_i_i_C(2)) * qJD(6) + (t622 * t681 * t691 + (-pkin(2) * t622 + (pkin(4) * t620 + pkin(9)) * t680) * qJD(2)) * t616, (t566 * t619 + t582 * t674) * r_i_i_C(1) + (t566 * t623 - t582 * t675) * r_i_i_C(2) - t566 * t626 + t633 * t565 + t627 * t581, (t624 * t655 - t566 * t620 + (-t582 * t624 - t597 * t620) * qJD(4)) * pkin(4) + t628, t628, (-t546 * t619 + t565 * t623) * r_i_i_C(1) + (-t546 * t623 - t565 * t619) * r_i_i_C(2) + ((-t568 * t623 - t581 * t619) * r_i_i_C(1) + (t568 * t619 - t581 * t623) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end