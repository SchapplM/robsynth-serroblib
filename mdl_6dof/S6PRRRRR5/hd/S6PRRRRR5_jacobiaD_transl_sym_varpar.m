% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRRR5
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
% Datum: 2019-10-09 23:21
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:21:28
	% EndTime: 2019-10-09 23:21:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:21:28
	% EndTime: 2019-10-09 23:21:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:21:28
	% EndTime: 2019-10-09 23:21:28
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
	% StartTime: 2019-10-09 23:21:29
	% EndTime: 2019-10-09 23:21:30
	% DurationCPUTime: 0.35s
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
	% StartTime: 2019-10-09 23:21:31
	% EndTime: 2019-10-09 23:21:32
	% DurationCPUTime: 0.83s
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
	% StartTime: 2019-10-09 23:21:34
	% EndTime: 2019-10-09 23:21:35
	% DurationCPUTime: 1.33s
	% Computational Cost: add. (1151->179), mult. (3722->319), div. (0->0), fcn. (4116->14), ass. (0->103)
	t552 = sin(qJ(2));
	t555 = cos(qJ(2));
	t547 = cos(pkin(13));
	t606 = cos(pkin(6));
	t585 = t547 * t606;
	t605 = sin(pkin(13));
	t533 = -t605 * t552 + t555 * t585;
	t549 = sin(qJ(5));
	t553 = cos(qJ(5));
	t568 = (t549 * r_i_i_C(1) + t553 * r_i_i_C(2)) * qJD(5);
	t609 = r_i_i_C(3) + pkin(11);
	t608 = cos(qJ(3));
	t545 = sin(pkin(7));
	t607 = pkin(9) * t545;
	t550 = sin(qJ(4));
	t604 = t545 * t550;
	t554 = cos(qJ(4));
	t603 = t545 * t554;
	t602 = t545 * t555;
	t546 = sin(pkin(6));
	t601 = t546 * t547;
	t548 = cos(pkin(7));
	t551 = sin(qJ(3));
	t600 = t548 * t551;
	t599 = t551 * t552;
	t598 = t551 * t555;
	t597 = qJD(2) * t546;
	t596 = qJD(5) * t549;
	t595 = qJD(5) * t553;
	t594 = t545 * t601;
	t593 = t545 * t546 * t552;
	t592 = t546 * t599;
	t562 = -t552 * t585 - t605 * t555;
	t591 = t562 * t608;
	t590 = t548 * t608;
	t589 = t608 * t552;
	t588 = t608 * t555;
	t587 = t545 * t597;
	t586 = t546 * t605;
	t584 = t606 * t545;
	t582 = t548 * t588;
	t581 = t552 * t587;
	t580 = t555 * t587;
	t578 = t551 * t584;
	t577 = t545 * t586;
	t575 = t606 * t605;
	t574 = t546 * t582;
	t506 = -t591 + (t533 * t548 - t594) * t551;
	t520 = -t533 * t545 - t548 * t601;
	t498 = t506 * t554 + t520 * t550;
	t573 = -t506 * t550 + t520 * t554;
	t560 = t547 * t552 + t555 * t575;
	t561 = -t547 * t555 + t552 * t575;
	t508 = -t561 * t608 + (-t548 * t560 + t577) * t551;
	t521 = t545 * t560 + t548 * t586;
	t500 = t508 * t554 + t521 * t550;
	t572 = -t508 * t550 + t521 * t554;
	t564 = t548 * t598 + t589;
	t519 = t564 * t546 + t578;
	t532 = -t546 * t602 + t606 * t548;
	t510 = t519 * t554 + t532 * t550;
	t571 = -t519 * t550 + t532 * t554;
	t570 = t553 * r_i_i_C(1) - t549 * r_i_i_C(2) + pkin(4);
	t569 = t608 * t584;
	t514 = t533 * t608 + t562 * t600;
	t501 = t514 * t554 - t562 * t604;
	t516 = -t560 * t608 + t561 * t600;
	t502 = t516 * t554 - t561 * t604;
	t565 = -t548 * t599 + t588;
	t527 = t565 * t546;
	t517 = t527 * t554 + t550 * t593;
	t567 = -t533 * t551 + t562 * t590;
	t566 = t551 * t560 + t561 * t590;
	t563 = t548 * t589 + t598;
	t559 = -t609 * t550 - t570 * t554 - pkin(3);
	t558 = t533 * t590 + t551 * t562 - t608 * t594;
	t557 = t551 * t561 - t560 * t590 + t608 * t577;
	t556 = t554 * t568 + (t570 * t550 - t609 * t554) * qJD(4);
	t531 = t561 * qJD(2);
	t530 = t560 * qJD(2);
	t529 = t562 * qJD(2);
	t528 = t533 * qJD(2);
	t526 = t563 * t546;
	t518 = -t569 - t574 + t592;
	t512 = (-qJD(2) * t564 - qJD(3) * t563) * t546;
	t511 = -qJD(2) * t574 - t546 * qJD(3) * t588 + (qJD(3) * t548 + qJD(2)) * t592;
	t504 = qJD(3) * t569 + ((t582 - t599) * qJD(3) + t565 * qJD(2)) * t546;
	t503 = qJD(3) * t578 + (qJD(2) * t563 + qJD(3) * t564) * t546;
	t496 = t566 * qJD(3) + t530 * t600 + t531 * t608;
	t495 = t516 * qJD(3) - t530 * t590 + t531 * t551;
	t494 = t567 * qJD(3) - t528 * t600 + t529 * t608;
	t493 = t514 * qJD(3) + t528 * t590 + t529 * t551;
	t492 = t550 * t580 + t512 * t554 + (-t527 * t550 + t554 * t593) * qJD(4);
	t490 = t557 * qJD(3) - t530 * t608 + t531 * t600;
	t489 = t508 * qJD(3) - t530 * t551 - t531 * t590;
	t488 = t558 * qJD(3) + t528 * t608 + t529 * t600;
	t487 = t528 * t551 - t529 * t590 + (t533 * t600 - t551 * t594 - t591) * qJD(3);
	t486 = t571 * qJD(4) + t504 * t554 + t550 * t581;
	t484 = -t530 * t604 + t496 * t554 + (-t516 * t550 - t561 * t603) * qJD(4);
	t482 = t528 * t604 + t494 * t554 + (-t514 * t550 - t562 * t603) * qJD(4);
	t480 = t572 * qJD(4) + t490 * t554 - t531 * t604;
	t478 = t573 * qJD(4) + t488 * t554 - t529 * t604;
	t1 = [0, (t484 * t553 + t495 * t549) * r_i_i_C(1) + (-t484 * t549 + t495 * t553) * r_i_i_C(2) + t484 * pkin(4) + t496 * pkin(3) + t495 * pkin(10) + t531 * pkin(2) - t530 * t607 + t609 * (qJD(4) * t502 + t496 * t550 + t530 * t603) + ((-t502 * t549 - t553 * t566) * r_i_i_C(1) + (-t502 * t553 + t549 * t566) * r_i_i_C(2)) * qJD(5), (t490 * t549 + t508 * t595) * r_i_i_C(1) + (t490 * t553 - t508 * t596) * r_i_i_C(2) + t490 * pkin(10) + t559 * t489 - t556 * t557, t609 * t480 - t572 * t568 + t570 * (-t500 * qJD(4) - t490 * t550 - t531 * t603), (-t480 * t549 + t489 * t553) * r_i_i_C(1) + (-t480 * t553 - t489 * t549) * r_i_i_C(2) + ((-t500 * t553 + t549 * t557) * r_i_i_C(1) + (t500 * t549 + t553 * t557) * r_i_i_C(2)) * qJD(5), 0; 0, (t482 * t553 + t493 * t549) * r_i_i_C(1) + (-t482 * t549 + t493 * t553) * r_i_i_C(2) + t482 * pkin(4) + t494 * pkin(3) + t493 * pkin(10) + t529 * pkin(2) + t528 * t607 + t609 * (t501 * qJD(4) + t494 * t550 - t528 * t603) + ((-t501 * t549 - t553 * t567) * r_i_i_C(1) + (-t501 * t553 + t549 * t567) * r_i_i_C(2)) * qJD(5), (t488 * t549 + t506 * t595) * r_i_i_C(1) + (t488 * t553 - t506 * t596) * r_i_i_C(2) + t488 * pkin(10) + t559 * t487 - t556 * t558, t609 * t478 - t573 * t568 + t570 * (-t498 * qJD(4) - t488 * t550 - t529 * t603), (-t478 * t549 + t487 * t553) * r_i_i_C(1) + (-t478 * t553 - t487 * t549) * r_i_i_C(2) + ((-t498 * t553 + t549 * t558) * r_i_i_C(1) + (t498 * t549 + t553 * t558) * r_i_i_C(2)) * qJD(5), 0; 0, (t492 * t553 - t511 * t549) * r_i_i_C(1) + (-t492 * t549 - t511 * t553) * r_i_i_C(2) + t492 * pkin(4) + t512 * pkin(3) - t511 * pkin(10) + t609 * (qJD(4) * t517 + t512 * t550 - t554 * t580) + ((-t517 * t549 + t526 * t553) * r_i_i_C(1) + (-t517 * t553 - t526 * t549) * r_i_i_C(2)) * qJD(5) + (-pkin(2) * t552 + pkin(9) * t602) * t597, (t504 * t549 + t519 * t595) * r_i_i_C(1) + (t504 * t553 - t519 * t596) * r_i_i_C(2) + t504 * pkin(10) + t559 * t503 + t556 * t518, t609 * t486 - t571 * t568 + t570 * (-t510 * qJD(4) - t504 * t550 + t554 * t581), (-t486 * t549 + t503 * t553) * r_i_i_C(1) + (-t486 * t553 - t503 * t549) * r_i_i_C(2) + ((-t510 * t553 - t518 * t549) * r_i_i_C(1) + (t510 * t549 - t518 * t553) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:21:34
	% EndTime: 2019-10-09 23:21:36
	% DurationCPUTime: 1.84s
	% Computational Cost: add. (1773->198), mult. (5101->335), div. (0->0), fcn. (5674->16), ass. (0->124)
	t597 = sin(qJ(2));
	t600 = cos(qJ(2));
	t592 = cos(pkin(13));
	t669 = cos(pkin(6));
	t642 = t592 * t669;
	t668 = sin(pkin(13));
	t573 = -t668 * t597 + t600 * t642;
	t589 = qJ(5) + qJ(6);
	t586 = sin(t589);
	t587 = cos(t589);
	t588 = qJD(5) + qJD(6);
	t594 = sin(qJ(5));
	t673 = qJD(5) * t594 * pkin(5) + (t586 * r_i_i_C(1) + t587 * r_i_i_C(2)) * t588;
	t672 = cos(qJ(3));
	t590 = sin(pkin(7));
	t671 = pkin(9) * t590;
	t670 = r_i_i_C(3) + pkin(12) + pkin(11);
	t667 = t586 * t588;
	t666 = t587 * t588;
	t595 = sin(qJ(4));
	t665 = t590 * t595;
	t599 = cos(qJ(4));
	t664 = t590 * t599;
	t663 = t590 * t600;
	t591 = sin(pkin(6));
	t662 = t591 * t592;
	t593 = cos(pkin(7));
	t596 = sin(qJ(3));
	t661 = t593 * t596;
	t660 = t596 * t597;
	t659 = t596 * t600;
	t568 = t573 * qJD(2);
	t609 = -t597 * t642 - t668 * t600;
	t569 = t609 * qJD(2);
	t648 = t593 * t672;
	t649 = t609 * t672;
	t652 = t590 * t662;
	t527 = t568 * t596 - t569 * t648 + (t573 * t661 - t596 * t652 - t649) * qJD(3);
	t546 = -t649 + (t573 * t593 - t652) * t596;
	t560 = -t573 * t590 - t593 * t662;
	t538 = t546 * t599 + t560 * t595;
	t635 = t538 * t588 - t527;
	t604 = t573 * t648 + t596 * t609 - t672 * t652;
	t528 = t604 * qJD(3) + t568 * t672 + t569 * t661;
	t619 = -t546 * t595 + t560 * t599;
	t518 = qJD(4) * t619 + t528 * t599 - t569 * t665;
	t640 = t588 * t604 - t518;
	t658 = (t586 * t640 - t587 * t635) * r_i_i_C(1) + (t586 * t635 + t587 * t640) * r_i_i_C(2);
	t621 = t669 * t668;
	t607 = t592 * t597 + t600 * t621;
	t608 = -t592 * t600 + t597 * t621;
	t643 = t591 * t668;
	t623 = t590 * t643;
	t548 = -t608 * t672 + (-t593 * t607 + t623) * t596;
	t570 = t607 * qJD(2);
	t571 = t608 * qJD(2);
	t529 = t548 * qJD(3) - t570 * t596 - t571 * t648;
	t561 = t590 * t607 + t593 * t643;
	t540 = t548 * t599 + t561 * t595;
	t634 = t540 * t588 - t529;
	t603 = t596 * t608 - t607 * t648 + t672 * t623;
	t530 = t603 * qJD(3) - t570 * t672 + t571 * t661;
	t618 = -t548 * t595 + t561 * t599;
	t520 = qJD(4) * t618 + t530 * t599 - t571 * t665;
	t639 = t588 * t603 - t520;
	t657 = (t586 * t639 - t587 * t634) * r_i_i_C(1) + (t586 * t634 + t587 * t639) * r_i_i_C(2);
	t647 = t672 * t597;
	t610 = t593 * t647 + t659;
	t611 = t593 * t659 + t647;
	t644 = t590 * t669;
	t624 = t596 * t644;
	t543 = qJD(3) * t624 + (qJD(2) * t610 + qJD(3) * t611) * t591;
	t559 = t591 * t611 + t624;
	t572 = -t591 * t663 + t669 * t593;
	t550 = t559 * t599 + t572 * t595;
	t630 = t550 * t588 - t543;
	t646 = t672 * t600;
	t612 = -t593 * t660 + t646;
	t616 = t672 * t644;
	t628 = t593 * t646;
	t544 = qJD(3) * t616 + ((t628 - t660) * qJD(3) + t612 * qJD(2)) * t591;
	t617 = -t559 * t595 + t572 * t599;
	t655 = qJD(2) * t591;
	t645 = t590 * t655;
	t627 = t597 * t645;
	t526 = qJD(4) * t617 + t544 * t599 + t595 * t627;
	t620 = t591 * t628;
	t650 = t591 * t660;
	t558 = -t616 - t620 + t650;
	t636 = -t558 * t588 - t526;
	t656 = (t586 * t636 - t587 * t630) * r_i_i_C(1) + (t586 * t630 + t587 * t636) * r_i_i_C(2);
	t598 = cos(qJ(5));
	t654 = qJD(5) * t598;
	t651 = t590 * t591 * t597;
	t614 = -t573 * t596 + t609 * t648;
	t534 = t614 * qJD(3) - t568 * t661 + t569 * t672;
	t554 = t573 * t672 + t609 * t661;
	t522 = t568 * t665 + t534 * t599 + (-t554 * t595 - t609 * t664) * qJD(4);
	t638 = -t588 * t614 + t522;
	t613 = t596 * t607 + t608 * t648;
	t536 = t613 * qJD(3) + t570 * t661 + t571 * t672;
	t556 = -t607 * t672 + t608 * t661;
	t524 = -t570 * t665 + t536 * t599 + (-t556 * t595 - t608 * t664) * qJD(4);
	t637 = -t588 * t613 + t524;
	t552 = (-qJD(2) * t611 - qJD(3) * t610) * t591;
	t567 = t612 * t591;
	t626 = t600 * t645;
	t532 = t595 * t626 + t552 * t599 + (-t567 * t595 + t599 * t651) * qJD(4);
	t566 = t610 * t591;
	t633 = t566 * t588 + t532;
	t533 = t554 * qJD(3) + t568 * t648 + t569 * t596;
	t541 = t554 * t599 - t609 * t665;
	t632 = -t541 * t588 + t533;
	t535 = t556 * qJD(3) - t570 * t648 + t571 * t596;
	t542 = t556 * t599 - t608 * t665;
	t631 = -t542 * t588 + t535;
	t551 = -qJD(2) * t620 - t591 * qJD(3) * t646 + (qJD(3) * t593 + qJD(2)) * t650;
	t557 = t567 * t599 + t595 * t651;
	t629 = -t557 * t588 - t551;
	t585 = pkin(5) * t598 + pkin(4);
	t615 = r_i_i_C(1) * t587 - r_i_i_C(2) * t586 + t585;
	t605 = -t670 * t595 - t615 * t599 - pkin(3);
	t602 = t673 * t599 + (t615 * t595 - t670 * t599) * qJD(4);
	t1 = [0, -t570 * t671 + t571 * pkin(2) + t536 * pkin(3) + t535 * pkin(10) + t524 * t585 + (r_i_i_C(1) * t637 + r_i_i_C(2) * t631) * t587 + (r_i_i_C(1) * t631 - r_i_i_C(2) * t637) * t586 + t670 * (qJD(4) * t542 + t536 * t595 + t570 * t664) + (t535 * t594 + (-t542 * t594 - t598 * t613) * qJD(5)) * pkin(5), (t530 * t586 + t548 * t666) * r_i_i_C(1) + (t530 * t587 - t548 * t667) * r_i_i_C(2) + t530 * pkin(10) + (t530 * t594 + t548 * t654) * pkin(5) + t605 * t529 - t602 * t603, t670 * t520 - t673 * t618 + t615 * (-qJD(4) * t540 - t530 * t595 - t571 * t664), (-t520 * t594 + t529 * t598 + (-t540 * t598 + t594 * t603) * qJD(5)) * pkin(5) + t657, t657; 0, t568 * t671 + t569 * pkin(2) + t534 * pkin(3) + t533 * pkin(10) + t522 * t585 + (r_i_i_C(1) * t638 + r_i_i_C(2) * t632) * t587 + (r_i_i_C(1) * t632 - r_i_i_C(2) * t638) * t586 + t670 * (qJD(4) * t541 + t534 * t595 - t568 * t664) + (t533 * t594 + (-t541 * t594 - t598 * t614) * qJD(5)) * pkin(5), (t528 * t586 + t546 * t666) * r_i_i_C(1) + (t528 * t587 - t546 * t667) * r_i_i_C(2) + t528 * pkin(10) + (t528 * t594 + t546 * t654) * pkin(5) + t605 * t527 - t602 * t604, t670 * t518 - t673 * t619 + t615 * (-qJD(4) * t538 - t528 * t595 - t569 * t664), (-t518 * t594 + t527 * t598 + (-t538 * t598 + t594 * t604) * qJD(5)) * pkin(5) + t658, t658; 0, t552 * pkin(3) - t551 * pkin(10) + t532 * t585 + (r_i_i_C(1) * t633 + r_i_i_C(2) * t629) * t587 + (r_i_i_C(1) * t629 - r_i_i_C(2) * t633) * t586 + t670 * (qJD(4) * t557 + t552 * t595 - t599 * t626) + (-pkin(2) * t597 + pkin(9) * t663) * t655 + (-t551 * t594 + (-t557 * t594 + t566 * t598) * qJD(5)) * pkin(5), (t544 * t586 + t559 * t666) * r_i_i_C(1) + (t544 * t587 - t559 * t667) * r_i_i_C(2) + t544 * pkin(10) + (t544 * t594 + t559 * t654) * pkin(5) + t605 * t543 + t602 * t558, t670 * t526 - t673 * t617 + t615 * (-qJD(4) * t550 - t544 * t595 + t599 * t627), (-t526 * t594 + t543 * t598 + (-t550 * t598 - t558 * t594) * qJD(5)) * pkin(5) + t656, t656;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end