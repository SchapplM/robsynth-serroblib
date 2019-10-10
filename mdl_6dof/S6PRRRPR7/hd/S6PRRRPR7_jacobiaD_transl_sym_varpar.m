% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
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
	% StartTime: 2019-10-09 22:58:05
	% EndTime: 2019-10-09 22:58:05
	% DurationCPUTime: 0.35s
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
	% StartTime: 2019-10-09 22:58:06
	% EndTime: 2019-10-09 22:58:07
	% DurationCPUTime: 0.83s
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
	% StartTime: 2019-10-09 22:58:09
	% EndTime: 2019-10-09 22:58:10
	% DurationCPUTime: 0.93s
	% Computational Cost: add. (893->134), mult. (2926->236), div. (0->0), fcn. (3175->14), ass. (0->87)
	t485 = sin(pkin(12));
	t489 = cos(pkin(12));
	t494 = sin(qJ(2));
	t491 = cos(pkin(6));
	t497 = cos(qJ(2));
	t529 = t491 * t497;
	t475 = -t485 * t494 + t489 * t529;
	t493 = sin(qJ(3));
	t496 = cos(qJ(3));
	t530 = t491 * t494;
	t507 = t485 * t530 - t489 * t497;
	t490 = cos(pkin(7));
	t508 = t485 * t529 + t489 * t494;
	t486 = sin(pkin(7));
	t487 = sin(pkin(6));
	t538 = t486 * t487;
	t509 = t485 * t538 - t490 * t508;
	t460 = t509 * t493 - t496 * t507;
	t476 = t485 * t497 + t489 * t530;
	t510 = -t475 * t490 + t489 * t538;
	t544 = -t476 * t496 + t510 * t493;
	t543 = pkin(9) * t486;
	t542 = r_i_i_C(3) + qJ(5);
	t537 = t486 * t491;
	t492 = sin(qJ(4));
	t536 = t486 * t492;
	t495 = cos(qJ(4));
	t535 = t486 * t495;
	t534 = t486 * t497;
	t533 = t487 * t490;
	t532 = t490 * t493;
	t531 = t490 * t496;
	t528 = t493 * t494;
	t527 = t493 * t497;
	t526 = t494 * t496;
	t525 = t496 * t497;
	t524 = qJD(2) * t487;
	t523 = t494 * t538;
	t521 = t497 * t524;
	t520 = qJD(3) * t537;
	t519 = qJD(2) * t523;
	t518 = t486 * t521;
	t467 = -t475 * t486 - t489 * t533;
	t517 = t467 * t492 - t495 * t544;
	t468 = t485 * t533 + t486 * t508;
	t516 = t460 * t495 + t468 * t492;
	t504 = t490 * t527 + t526;
	t466 = t504 * t487 + t493 * t537;
	t474 = -t487 * t534 + t491 * t490;
	t515 = t466 * t495 + t474 * t492;
	t484 = sin(pkin(13));
	t488 = cos(pkin(13));
	t514 = t488 * r_i_i_C(1) - t484 * r_i_i_C(2) + pkin(4);
	t513 = t484 * r_i_i_C(1) + t488 * r_i_i_C(2) + pkin(10);
	t463 = t475 * t496 - t476 * t532;
	t512 = -t463 * t492 + t476 * t535;
	t464 = -t496 * t508 + t507 * t532;
	t511 = -t464 * t492 - t507 * t535;
	t506 = t490 * t525 - t528;
	t505 = -t490 * t526 - t527;
	t503 = -t490 * t528 + t525;
	t469 = t503 * t487;
	t502 = -t469 * t492 + t495 * t523;
	t501 = -t476 * t493 - t510 * t496;
	t500 = t493 * t507 + t509 * t496;
	t499 = t542 * t492 + t514 * t495 + pkin(3);
	t498 = t492 * qJD(5) + (-t514 * t492 + t542 * t495) * qJD(4);
	t473 = t507 * qJD(2);
	t472 = t508 * qJD(2);
	t471 = t476 * qJD(2);
	t470 = t475 * qJD(2);
	t462 = (-t504 * qJD(2) + t505 * qJD(3)) * t487;
	t461 = -t521 * t531 + (-qJD(3) * t525 + (qJD(3) * t490 + qJD(2)) * t528) * t487;
	t456 = t496 * t520 + (t503 * qJD(2) + t506 * qJD(3)) * t487;
	t454 = t472 * t532 + t473 * t496 + (t493 * t508 + t507 * t531) * qJD(3);
	t453 = t464 * qJD(3) - t472 * t531 + t473 * t493;
	t452 = -t470 * t532 - t471 * t496 + (-t475 * t493 - t476 * t531) * qJD(3);
	t451 = t463 * qJD(3) + t470 * t531 - t471 * t493;
	t450 = t502 * qJD(4) + t462 * t495 + t492 * t518;
	t448 = t500 * qJD(3) - t472 * t496 + t473 * t532;
	t446 = t501 * qJD(3) + t470 * t496 - t471 * t532;
	t443 = t515 * qJD(4) + t456 * t492 - t495 * t519;
	t442 = t511 * qJD(4) + t454 * t495 - t472 * t536;
	t440 = t512 * qJD(4) + t452 * t495 + t470 * t536;
	t437 = t516 * qJD(4) + t448 * t492 + t473 * t535;
	t435 = t517 * qJD(4) + t446 * t492 - t471 * t535;
	t1 = [0, (t442 * t488 + t453 * t484) * r_i_i_C(1) + (-t442 * t484 + t453 * t488) * r_i_i_C(2) + t442 * pkin(4) - t511 * qJD(5) + t454 * pkin(3) + t453 * pkin(10) + t473 * pkin(2) - t472 * t543 + t542 * (t472 * t535 + t454 * t492 + (t464 * t495 - t507 * t536) * qJD(4)), t513 * t448 + t498 * t500 + t499 * (-t460 * qJD(3) + t472 * t493 + t473 * t531), t516 * qJD(5) + t542 * (-t473 * t536 + t448 * t495 + (-t460 * t492 + t468 * t495) * qJD(4)) - t514 * t437, t437, 0; 0, (t440 * t488 + t451 * t484) * r_i_i_C(1) + (-t440 * t484 + t451 * t488) * r_i_i_C(2) + t440 * pkin(4) - t512 * qJD(5) + t452 * pkin(3) + t451 * pkin(10) - t471 * pkin(2) + t470 * t543 + t542 * (-t470 * t535 + t452 * t492 + (t463 * t495 + t476 * t536) * qJD(4)), t513 * t446 + t498 * t501 + t499 * (t544 * qJD(3) - t470 * t493 - t471 * t531), t517 * qJD(5) + t542 * (t471 * t536 + t446 * t495 + (t467 * t495 + t492 * t544) * qJD(4)) - t514 * t435, t435, 0; 0, (t450 * t488 - t461 * t484) * r_i_i_C(1) + (-t450 * t484 - t461 * t488) * r_i_i_C(2) + t450 * pkin(4) - t502 * qJD(5) + t462 * pkin(3) - t461 * pkin(10) + t542 * (-t495 * t518 + t462 * t492 + (t469 * t495 + t492 * t523) * qJD(4)) + (-pkin(2) * t494 + pkin(9) * t534) * t524, t513 * t456 + t498 * (t506 * t487 + t496 * t537) + t499 * (-t493 * t520 + (t505 * qJD(2) - t504 * qJD(3)) * t487), t515 * qJD(5) + t542 * (t492 * t519 + t456 * t495 + (-t466 * t492 + t474 * t495) * qJD(4)) - t514 * t443, t443, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:09
	% EndTime: 2019-10-09 22:58:10
	% DurationCPUTime: 1.37s
	% Computational Cost: add. (1403->169), mult. (4167->294), div. (0->0), fcn. (4622->16), ass. (0->106)
	t567 = sin(qJ(2));
	t569 = cos(qJ(2));
	t562 = cos(pkin(12));
	t624 = cos(pkin(6));
	t604 = t562 * t624;
	t623 = sin(pkin(12));
	t542 = -t623 * t567 + t569 * t604;
	t558 = pkin(13) + qJ(6);
	t556 = sin(t558);
	t557 = cos(t558);
	t595 = t556 * r_i_i_C(1) + t557 * r_i_i_C(2);
	t584 = qJD(6) * t595;
	t627 = cos(qJ(3));
	t560 = sin(pkin(7));
	t626 = pkin(9) * t560;
	t625 = r_i_i_C(3) + pkin(11) + qJ(5);
	t565 = sin(qJ(4));
	t622 = t560 * t565;
	t568 = cos(qJ(4));
	t621 = t560 * t568;
	t620 = t560 * t569;
	t561 = sin(pkin(6));
	t619 = t561 * t562;
	t563 = cos(pkin(7));
	t566 = sin(qJ(3));
	t618 = t563 * t566;
	t617 = t566 * t567;
	t616 = t566 * t569;
	t615 = qJD(2) * t561;
	t614 = t560 * t619;
	t613 = t560 * t561 * t567;
	t612 = t561 * t617;
	t576 = -t567 * t604 - t623 * t569;
	t611 = t576 * t627;
	t610 = t563 * t627;
	t609 = t627 * t567;
	t608 = t627 * t569;
	t607 = t560 * t615;
	t606 = t560 * t624;
	t605 = t561 * t623;
	t602 = t563 * t608;
	t601 = t567 * t607;
	t600 = t569 * t607;
	t598 = t566 * t606;
	t597 = t560 * t605;
	t596 = r_i_i_C(1) * t557 - r_i_i_C(2) * t556;
	t594 = t624 * t623;
	t593 = t561 * t602;
	t515 = -t611 + (t542 * t563 - t614) * t566;
	t529 = -t542 * t560 - t563 * t619;
	t507 = t515 * t568 + t529 * t565;
	t592 = -t515 * t565 + t529 * t568;
	t574 = t562 * t567 + t569 * t594;
	t575 = -t562 * t569 + t567 * t594;
	t517 = -t575 * t627 + (-t563 * t574 + t597) * t566;
	t530 = t560 * t574 + t563 * t605;
	t509 = t517 * t568 + t530 * t565;
	t591 = -t517 * t565 + t530 * t568;
	t579 = t563 * t616 + t609;
	t528 = t579 * t561 + t598;
	t541 = -t561 * t620 + t624 * t563;
	t519 = t528 * t568 + t541 * t565;
	t590 = -t528 * t565 + t541 * t568;
	t589 = t627 * t606;
	t588 = cos(pkin(13)) * pkin(5) + pkin(4) + t596;
	t523 = t542 * t627 + t576 * t618;
	t587 = -t523 * t565 - t576 * t621;
	t510 = t523 * t568 - t576 * t622;
	t525 = -t574 * t627 + t575 * t618;
	t586 = -t525 * t565 - t575 * t621;
	t511 = t525 * t568 - t575 * t622;
	t585 = qJD(6) * t596;
	t580 = -t563 * t617 + t608;
	t536 = t580 * t561;
	t583 = -t536 * t565 + t568 * t613;
	t526 = t536 * t568 + t565 * t613;
	t582 = -t542 * t566 + t576 * t610;
	t581 = t566 * t574 + t575 * t610;
	t578 = t563 * t609 + t616;
	t577 = sin(pkin(13)) * pkin(5) + pkin(10) + t595;
	t573 = -t625 * t565 - t588 * t568 - pkin(3);
	t572 = t542 * t610 + t566 * t576 - t627 * t614;
	t571 = t566 * t575 - t574 * t610 + t627 * t597;
	t570 = -t565 * qJD(5) + t568 * t584 + (t588 * t565 - t625 * t568) * qJD(4);
	t540 = t575 * qJD(2);
	t539 = t574 * qJD(2);
	t538 = t576 * qJD(2);
	t537 = t542 * qJD(2);
	t535 = t578 * t561;
	t527 = -t589 - t593 + t612;
	t521 = (-t579 * qJD(2) - t578 * qJD(3)) * t561;
	t513 = qJD(3) * t589 + ((t602 - t617) * qJD(3) + t580 * qJD(2)) * t561;
	t512 = qJD(3) * t598 + (t578 * qJD(2) + t579 * qJD(3)) * t561;
	t505 = t581 * qJD(3) + t539 * t618 + t540 * t627;
	t503 = t582 * qJD(3) - t537 * t618 + t538 * t627;
	t499 = t571 * qJD(3) - t539 * t627 + t540 * t618;
	t498 = t517 * qJD(3) - t539 * t566 - t540 * t610;
	t497 = t572 * qJD(3) + t537 * t627 + t538 * t618;
	t496 = t537 * t566 - t538 * t610 + (t542 * t618 - t566 * t614 - t611) * qJD(3);
	t495 = t590 * qJD(4) + t513 * t568 + t565 * t601;
	t494 = t519 * qJD(4) + t513 * t565 - t568 * t601;
	t489 = t591 * qJD(4) + t499 * t568 - t540 * t622;
	t488 = t509 * qJD(4) + t499 * t565 + t540 * t621;
	t487 = t592 * qJD(4) + t497 * t568 - t538 * t622;
	t486 = t507 * qJD(4) + t497 * t565 + t538 * t621;
	t1 = [0, -t586 * qJD(5) + t505 * pkin(3) + t540 * pkin(2) - t539 * t626 + t588 * (t586 * qJD(4) + t505 * t568 - t539 * t622) + t577 * (t525 * qJD(3) - t539 * t610 + t540 * t566) + t625 * (t511 * qJD(4) + t505 * t565 + t539 * t621) + ((-t511 * t556 - t557 * t581) * r_i_i_C(1) + (-t511 * t557 + t556 * t581) * r_i_i_C(2)) * qJD(6), t573 * t498 + t577 * t499 + t517 * t585 - t570 * t571, qJD(5) * t509 - t588 * t488 + t625 * t489 - t591 * t584, t488, (-t489 * t556 + t498 * t557) * r_i_i_C(1) + (-t489 * t557 - t498 * t556) * r_i_i_C(2) + ((-t509 * t557 + t556 * t571) * r_i_i_C(1) + (t509 * t556 + t557 * t571) * r_i_i_C(2)) * qJD(6); 0, -t587 * qJD(5) + t503 * pkin(3) + t538 * pkin(2) + t537 * t626 + t588 * (t587 * qJD(4) + t503 * t568 + t537 * t622) + t577 * (t523 * qJD(3) + t537 * t610 + t538 * t566) + t625 * (t510 * qJD(4) + t503 * t565 - t537 * t621) + ((-t510 * t556 - t557 * t582) * r_i_i_C(1) + (-t510 * t557 + t556 * t582) * r_i_i_C(2)) * qJD(6), t573 * t496 + t577 * t497 + t515 * t585 - t570 * t572, qJD(5) * t507 - t588 * t486 + t625 * t487 - t592 * t584, t486, (-t487 * t556 + t496 * t557) * r_i_i_C(1) + (-t487 * t557 - t496 * t556) * r_i_i_C(2) + ((-t507 * t557 + t556 * t572) * r_i_i_C(1) + (t507 * t556 + t557 * t572) * r_i_i_C(2)) * qJD(6); 0, -t583 * qJD(5) + t521 * pkin(3) + t588 * (t583 * qJD(4) + t521 * t568 + t565 * t600) - t577 * (-qJD(2) * t593 - t561 * qJD(3) * t608 + (qJD(3) * t563 + qJD(2)) * t612) + t625 * (t526 * qJD(4) + t521 * t565 - t568 * t600) + ((-t526 * t556 + t535 * t557) * r_i_i_C(1) + (-t526 * t557 - t535 * t556) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t567 + pkin(9) * t620) * t615, t573 * t512 + t577 * t513 + t570 * t527 + t528 * t585, qJD(5) * t519 - t588 * t494 + t625 * t495 - t590 * t584, t494, (-t495 * t556 + t512 * t557) * r_i_i_C(1) + (-t495 * t557 - t512 * t556) * r_i_i_C(2) + ((-t519 * t557 - t527 * t556) * r_i_i_C(1) + (t519 * t556 - t527 * t557) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end