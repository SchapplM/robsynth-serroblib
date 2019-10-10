% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(12));
	t50 = sin(pkin(12));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:06
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (95->48), mult. (333->96), div. (0->0), fcn. (328->10), ass. (0->36)
	t218 = sin(pkin(7));
	t248 = (pkin(9) + r_i_i_C(3)) * t218;
	t217 = sin(pkin(12));
	t220 = cos(pkin(12));
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
	% StartTime: 2019-10-09 23:00:07
	% EndTime: 2019-10-09 23:00:08
	% DurationCPUTime: 0.84s
	% Computational Cost: add. (380->110), mult. (1266->209), div. (0->0), fcn. (1334->12), ass. (0->67)
	t400 = sin(qJ(3));
	t403 = cos(qJ(3));
	t393 = sin(pkin(12));
	t396 = cos(pkin(12));
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
	% StartTime: 2019-10-09 23:00:09
	% EndTime: 2019-10-09 23:00:09
	% DurationCPUTime: 0.82s
	% Computational Cost: add. (725->118), mult. (2362->212), div. (0->0), fcn. (2557->12), ass. (0->79)
	t437 = sin(pkin(12));
	t440 = cos(pkin(12));
	t445 = sin(qJ(2));
	t442 = cos(pkin(6));
	t448 = cos(qJ(2));
	t478 = t442 * t448;
	t431 = -t437 * t445 + t440 * t478;
	t444 = sin(qJ(3));
	t447 = cos(qJ(3));
	t479 = t442 * t445;
	t458 = t437 * t479 - t440 * t448;
	t441 = cos(pkin(7));
	t459 = t437 * t478 + t440 * t445;
	t438 = sin(pkin(7));
	t439 = sin(pkin(6));
	t487 = t438 * t439;
	t460 = t437 * t487 - t441 * t459;
	t416 = t460 * t444 - t447 * t458;
	t432 = t437 * t448 + t440 * t479;
	t461 = -t431 * t441 + t440 * t487;
	t495 = -t432 * t447 + t461 * t444;
	t494 = r_i_i_C(1) + pkin(10);
	t493 = r_i_i_C(2) - pkin(4);
	t492 = pkin(9) * t438;
	t491 = r_i_i_C(3) + qJ(5);
	t486 = t438 * t442;
	t443 = sin(qJ(4));
	t485 = t438 * t443;
	t446 = cos(qJ(4));
	t484 = t438 * t446;
	t483 = t438 * t448;
	t482 = t439 * t441;
	t481 = t441 * t444;
	t480 = t441 * t447;
	t477 = t444 * t445;
	t476 = t444 * t448;
	t475 = t445 * t447;
	t474 = t447 * t448;
	t473 = qJD(2) * t439;
	t472 = t445 * t487;
	t470 = t438 * t473;
	t469 = qJD(3) * t486;
	t468 = t445 * t470;
	t467 = t448 * t470;
	t423 = -t431 * t438 - t440 * t482;
	t466 = t423 * t443 - t446 * t495;
	t424 = t437 * t482 + t438 * t459;
	t465 = t416 * t446 + t424 * t443;
	t455 = t441 * t476 + t475;
	t422 = t455 * t439 + t444 * t486;
	t430 = -t439 * t483 + t442 * t441;
	t464 = t422 * t446 + t430 * t443;
	t419 = t431 * t447 - t432 * t481;
	t463 = -t419 * t443 + t432 * t484;
	t420 = -t447 * t459 + t458 * t481;
	t462 = -t420 * t443 - t458 * t484;
	t457 = t441 * t474 - t477;
	t456 = -t441 * t475 - t476;
	t454 = t441 * t477 - t474;
	t425 = t454 * t439;
	t453 = t425 * t443 + t446 * t472;
	t452 = t491 * t443 - t493 * t446 + pkin(3);
	t451 = -t432 * t444 - t461 * t447;
	t450 = t444 * t458 + t460 * t447;
	t449 = qJD(5) * t443 + (t493 * t443 + t491 * t446) * qJD(4);
	t429 = t458 * qJD(2);
	t428 = t459 * qJD(2);
	t427 = t432 * qJD(2);
	t426 = t431 * qJD(2);
	t418 = (-t455 * qJD(2) + t456 * qJD(3)) * t439;
	t412 = t447 * t469 + (-t454 * qJD(2) + t457 * qJD(3)) * t439;
	t410 = t428 * t481 + t429 * t447 + (t444 * t459 + t458 * t480) * qJD(3);
	t408 = -t426 * t481 - t427 * t447 + (-t431 * t444 - t432 * t480) * qJD(3);
	t404 = t450 * qJD(3) - t428 * t447 + t429 * t481;
	t402 = t451 * qJD(3) + t426 * t447 - t427 * t481;
	t399 = t464 * qJD(4) + t412 * t443 - t446 * t468;
	t393 = t465 * qJD(4) + t404 * t443 + t429 * t484;
	t391 = t466 * qJD(4) + t402 * t443 - t427 * t484;
	t1 = [0, -t462 * qJD(5) + t410 * pkin(3) + t429 * pkin(2) - t428 * t492 + t494 * (t420 * qJD(3) - t428 * t480 + t429 * t444) - t493 * (t462 * qJD(4) + t410 * t446 - t428 * t485) + t491 * (t428 * t484 + t410 * t443 + (t420 * t446 - t458 * t485) * qJD(4)), t494 * t404 + t449 * t450 + t452 * (-t416 * qJD(3) + t428 * t444 + t429 * t480), t465 * qJD(5) + t491 * (-t429 * t485 + t404 * t446 + (-t416 * t443 + t424 * t446) * qJD(4)) + t493 * t393, t393, 0; 0, -t463 * qJD(5) + t408 * pkin(3) - t427 * pkin(2) + t426 * t492 + t494 * (t419 * qJD(3) + t426 * t480 - t427 * t444) - t493 * (t463 * qJD(4) + t408 * t446 + t426 * t485) + t491 * (-t426 * t484 + t408 * t443 + (t419 * t446 + t432 * t485) * qJD(4)), t494 * t402 + t449 * t451 + t452 * (t495 * qJD(3) - t426 * t444 - t427 * t480), t466 * qJD(5) + t491 * (t427 * t485 + t402 * t446 + (t423 * t446 + t443 * t495) * qJD(4)) + t493 * t391, t391, 0; 0, -t453 * qJD(5) + t418 * pkin(3) - t494 * (-t457 * qJD(2) + t454 * qJD(3)) * t439 - t493 * (t453 * qJD(4) + t418 * t446 + t443 * t467) + t491 * (-t446 * t467 + t418 * t443 + (-t425 * t446 + t443 * t472) * qJD(4)) + (-pkin(2) * t445 + pkin(9) * t483) * t473, t494 * t412 + t449 * (t457 * t439 + t447 * t486) + t452 * (-t444 * t469 + (t456 * qJD(2) - t455 * qJD(3)) * t439), t464 * qJD(5) + t491 * (t443 * t468 + t412 * t446 + (-t422 * t443 + t430 * t446) * qJD(4)) + t493 * t399, t399, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:10
	% EndTime: 2019-10-09 23:00:12
	% DurationCPUTime: 1.89s
	% Computational Cost: add. (1414->173), mult. (4554->309), div. (0->0), fcn. (5043->14), ass. (0->115)
	t568 = sin(qJ(2));
	t571 = cos(qJ(2));
	t564 = cos(pkin(12));
	t645 = cos(pkin(6));
	t622 = t564 * t645;
	t643 = sin(pkin(12));
	t591 = t643 * t568 - t571 * t622;
	t644 = cos(pkin(7));
	t578 = t591 * t644;
	t562 = sin(pkin(7));
	t563 = sin(pkin(6));
	t639 = t563 * t564;
	t628 = t562 * t639;
	t655 = -t578 - t628;
	t609 = t645 * t643;
	t588 = t564 * t568 + t571 * t609;
	t623 = t563 * t643;
	t654 = t562 * t623 - t588 * t644;
	t647 = cos(qJ(3));
	t653 = t655 * t647;
	t652 = t654 * t647;
	t567 = sin(qJ(3));
	t589 = -t564 * t571 + t568 * t609;
	t528 = t654 * t567 - t589 * t647;
	t584 = qJD(2) * t589;
	t651 = t528 * qJD(4) + t562 * t584;
	t590 = -t568 * t622 - t643 * t571;
	t627 = t590 * t647;
	t526 = t655 * t567 - t627;
	t585 = qJD(2) * t590;
	t650 = t526 * qJD(4) + t562 * t585;
	t631 = r_i_i_C(3) + pkin(11) + pkin(4);
	t565 = sin(qJ(6));
	t569 = cos(qJ(6));
	t611 = t569 * r_i_i_C(1) - t565 * r_i_i_C(2);
	t649 = pkin(5) + pkin(10) + t611;
	t595 = t611 * qJD(6) + qJD(5);
	t646 = pkin(9) * t562;
	t566 = sin(qJ(4));
	t642 = t562 * t566;
	t570 = cos(qJ(4));
	t641 = t562 * t570;
	t640 = t562 * t571;
	t638 = t567 * t568;
	t637 = qJD(2) * t563;
	t636 = qJD(3) * t567;
	t635 = qJD(4) * t566;
	t634 = qJD(4) * t570;
	t630 = t562 * t563 * t568;
	t629 = t563 * t638;
	t626 = t647 * t571;
	t625 = t562 * t637;
	t624 = t562 * t645;
	t621 = t567 * t644;
	t619 = t644 * qJD(3);
	t618 = t568 * t625;
	t617 = t571 * t625;
	t615 = t567 * t624;
	t613 = t644 * t647;
	t612 = t567 * t619;
	t610 = -t565 * r_i_i_C(1) - t569 * r_i_i_C(2);
	t593 = t647 * t568 + t571 * t621;
	t539 = t593 * t563 + t615;
	t550 = -t563 * t640 + t645 * t644;
	t607 = t539 * t570 + t550 * t566;
	t529 = t539 * t566 - t550 * t570;
	t606 = t647 * t624;
	t605 = t571 * t613;
	t604 = qJ(5) - t610;
	t603 = qJD(3) * t613;
	t582 = t591 * t647;
	t534 = t590 * t621 - t582;
	t601 = -t534 * t566 - t590 * t641;
	t583 = t588 * t647;
	t536 = t589 * t621 - t583;
	t600 = -t536 * t566 - t589 * t641;
	t599 = qJD(6) * t610;
	t598 = t563 * t605;
	t594 = -t568 * t621 + t626;
	t547 = t594 * t563;
	t596 = -t547 * t566 + t570 * t630;
	t592 = t567 * t571 + t568 * t613;
	t587 = t588 * t567;
	t586 = t591 * t567;
	t577 = -t604 * t566 - t631 * t570 - pkin(3);
	t576 = t644 * t585;
	t575 = t644 * t584;
	t572 = -t595 * t566 + (t631 * t566 - t604 * t570) * qJD(4);
	t549 = t588 * qJD(2);
	t548 = t591 * qJD(2);
	t546 = t592 * t563;
	t541 = t588 * t562 + t644 * t623;
	t540 = t591 * t562 - t644 * t639;
	t538 = -t598 - t606 + t629;
	t535 = -t589 * t613 - t587;
	t533 = -t590 * t613 - t586;
	t532 = (-t593 * qJD(2) - t592 * qJD(3)) * t563;
	t527 = -t567 * t589 - t652;
	t525 = -t567 * t590 - t653;
	t524 = qJD(3) * t606 + ((t605 - t638) * qJD(3) + t594 * qJD(2)) * t563;
	t523 = qJD(3) * t615 + (t592 * qJD(2) + t593 * qJD(3)) * t563;
	t519 = t528 * t566 - t541 * t570;
	t517 = t526 * t566 - t540 * t570;
	t516 = qJD(3) * t587 + t549 * t621 + t647 * t584 + t589 * t603;
	t514 = qJD(3) * t586 + t548 * t621 + t647 * t585 + t590 * t603;
	t510 = t652 * qJD(3) - t549 * t647 + t567 * t575 + t589 * t636;
	t509 = t528 * qJD(3) - t549 * t567 - t647 * t575;
	t508 = t653 * qJD(3) - t548 * t647 + t567 * t576 + t590 * t636;
	t507 = -t548 * t567 - t647 * t576 - t628 * t636 + (-t567 * t578 - t627) * qJD(3);
	t505 = t607 * qJD(4) + t524 * t566 - t570 * t618;
	t503 = t549 * t641 + t516 * t566 + (t536 * t570 - t589 * t642) * qJD(4);
	t501 = t548 * t641 + t514 * t566 + (t534 * t570 - t590 * t642) * qJD(4);
	t499 = t510 * t566 + t541 * t635 + t651 * t570;
	t497 = t508 * t566 + t540 * t635 + t650 * t570;
	t1 = [0, (t503 * t565 + (-t535 * t565 - t569 * t600) * qJD(6)) * r_i_i_C(1) + (t503 * t569 + (-t535 * t569 + t565 * t600) * qJD(6)) * r_i_i_C(2) + t503 * qJ(5) - t600 * qJD(5) + t516 * pkin(3) + pkin(2) * t584 - t549 * t646 + t649 * (-qJD(3) * t583 - t549 * t613 + t567 * t584 + t589 * t612) + t631 * (t600 * qJD(4) + t516 * t570 - t549 * t642), t577 * t509 + t510 * t649 + t572 * t527 + t528 * t599, t595 * (t528 * t570 + t541 * t566) + t604 * (t510 * t570 + t541 * t634 - t651 * t566) - t631 * t499, t499, (t499 * t569 - t509 * t565) * r_i_i_C(1) + (-t499 * t565 - t509 * t569) * r_i_i_C(2) + ((-t519 * t565 - t527 * t569) * r_i_i_C(1) + (-t519 * t569 + t527 * t565) * r_i_i_C(2)) * qJD(6); 0, (t501 * t565 + (-t533 * t565 - t569 * t601) * qJD(6)) * r_i_i_C(1) + (t501 * t569 + (-t533 * t569 + t565 * t601) * qJD(6)) * r_i_i_C(2) + t501 * qJ(5) - t601 * qJD(5) + t514 * pkin(3) + pkin(2) * t585 - t548 * t646 + t649 * (-qJD(3) * t582 - t548 * t613 + t567 * t585 + t590 * t612) + t631 * (t601 * qJD(4) + t514 * t570 - t548 * t642), t577 * t507 + t508 * t649 + t572 * t525 + t526 * t599, t595 * (t526 * t570 + t540 * t566) + t604 * (t508 * t570 + t540 * t634 - t650 * t566) - t631 * t497, t497, (t497 * t569 - t507 * t565) * r_i_i_C(1) + (-t497 * t565 - t507 * t569) * r_i_i_C(2) + ((-t517 * t565 - t525 * t569) * r_i_i_C(1) + (-t517 * t569 + t525 * t565) * r_i_i_C(2)) * qJD(6); 0, t532 * pkin(3) - t596 * qJD(5) + t604 * (-t570 * t617 + t532 * t566 + (t547 * t570 + t566 * t630) * qJD(4)) - t649 * (-qJD(2) * t598 - t563 * qJD(3) * t626 + (t619 + qJD(2)) * t629) + ((-t546 * t565 - t569 * t596) * r_i_i_C(1) + (-t546 * t569 + t565 * t596) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t568 + pkin(9) * t640) * t637 + t631 * (t596 * qJD(4) + t532 * t570 + t566 * t617), t577 * t523 + t524 * t649 + t572 * t538 + t539 * t599, t595 * t607 + t604 * (-t529 * qJD(4) + t524 * t570 + t566 * t618) - t631 * t505, t505, (t505 * t569 - t523 * t565) * r_i_i_C(1) + (-t505 * t565 - t523 * t569) * r_i_i_C(2) + ((-t529 * t565 - t538 * t569) * r_i_i_C(1) + (-t529 * t569 + t538 * t565) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end