% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:11
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRP6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:11:38
	% EndTime: 2019-10-09 23:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:11:38
	% EndTime: 2019-10-09 23:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:11:38
	% EndTime: 2019-10-09 23:11:38
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
	% StartTime: 2019-10-09 23:11:39
	% EndTime: 2019-10-09 23:11:39
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
	% StartTime: 2019-10-09 23:11:40
	% EndTime: 2019-10-09 23:11:41
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
	% StartTime: 2019-10-09 23:11:43
	% EndTime: 2019-10-09 23:11:44
	% DurationCPUTime: 1.32s
	% Computational Cost: add. (1151->179), mult. (3722->319), div. (0->0), fcn. (4116->14), ass. (0->103)
	t552 = sin(qJ(2));
	t555 = cos(qJ(2));
	t547 = cos(pkin(12));
	t606 = cos(pkin(6));
	t585 = t547 * t606;
	t605 = sin(pkin(12));
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
	t486 = qJD(4) * t571 + t504 * t554 + t550 * t581;
	t484 = -t530 * t604 + t496 * t554 + (-t516 * t550 - t561 * t603) * qJD(4);
	t482 = t528 * t604 + t494 * t554 + (-t514 * t550 - t562 * t603) * qJD(4);
	t480 = qJD(4) * t572 + t490 * t554 - t531 * t604;
	t478 = qJD(4) * t573 + t488 * t554 - t529 * t604;
	t1 = [0, (t484 * t553 + t495 * t549) * r_i_i_C(1) + (-t484 * t549 + t495 * t553) * r_i_i_C(2) + t484 * pkin(4) + t496 * pkin(3) + t495 * pkin(10) + t531 * pkin(2) - t530 * t607 + t609 * (qJD(4) * t502 + t496 * t550 + t530 * t603) + ((-t502 * t549 - t553 * t566) * r_i_i_C(1) + (-t502 * t553 + t549 * t566) * r_i_i_C(2)) * qJD(5), (t490 * t549 + t508 * t595) * r_i_i_C(1) + (t490 * t553 - t508 * t596) * r_i_i_C(2) + t490 * pkin(10) + t559 * t489 - t556 * t557, t609 * t480 - t572 * t568 + t570 * (-qJD(4) * t500 - t490 * t550 - t531 * t603), (-t480 * t549 + t489 * t553) * r_i_i_C(1) + (-t480 * t553 - t489 * t549) * r_i_i_C(2) + ((-t500 * t553 + t549 * t557) * r_i_i_C(1) + (t500 * t549 + t553 * t557) * r_i_i_C(2)) * qJD(5), 0; 0, (t482 * t553 + t493 * t549) * r_i_i_C(1) + (-t482 * t549 + t493 * t553) * r_i_i_C(2) + t482 * pkin(4) + t494 * pkin(3) + t493 * pkin(10) + t529 * pkin(2) + t528 * t607 + t609 * (qJD(4) * t501 + t494 * t550 - t528 * t603) + ((-t501 * t549 - t553 * t567) * r_i_i_C(1) + (-t501 * t553 + t549 * t567) * r_i_i_C(2)) * qJD(5), (t488 * t549 + t506 * t595) * r_i_i_C(1) + (t488 * t553 - t506 * t596) * r_i_i_C(2) + t488 * pkin(10) + t559 * t487 - t556 * t558, t609 * t478 - t573 * t568 + t570 * (-qJD(4) * t498 - t488 * t550 - t529 * t603), (-t478 * t549 + t487 * t553) * r_i_i_C(1) + (-t478 * t553 - t487 * t549) * r_i_i_C(2) + ((-t498 * t553 + t549 * t558) * r_i_i_C(1) + (t498 * t549 + t553 * t558) * r_i_i_C(2)) * qJD(5), 0; 0, (t492 * t553 - t511 * t549) * r_i_i_C(1) + (-t492 * t549 - t511 * t553) * r_i_i_C(2) + t492 * pkin(4) + t512 * pkin(3) - t511 * pkin(10) + t609 * (qJD(4) * t517 + t512 * t550 - t554 * t580) + ((-t517 * t549 + t526 * t553) * r_i_i_C(1) + (-t517 * t553 - t526 * t549) * r_i_i_C(2)) * qJD(5) + (-pkin(2) * t552 + pkin(9) * t602) * t597, (t504 * t549 + t519 * t595) * r_i_i_C(1) + (t504 * t553 - t519 * t596) * r_i_i_C(2) + t504 * pkin(10) + t559 * t503 + t556 * t518, t609 * t486 - t571 * t568 + t570 * (-qJD(4) * t510 - t504 * t550 + t554 * t581), (-t486 * t549 + t503 * t553) * r_i_i_C(1) + (-t486 * t553 - t503 * t549) * r_i_i_C(2) + ((-t510 * t553 - t518 * t549) * r_i_i_C(1) + (t510 * t549 - t518 * t553) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:11:47
	% EndTime: 2019-10-09 23:11:48
	% DurationCPUTime: 1.70s
	% Computational Cost: add. (2052->202), mult. (6535->341), div. (0->0), fcn. (7368->14), ass. (0->121)
	t640 = sin(qJ(2));
	t643 = cos(qJ(2));
	t635 = cos(pkin(12));
	t706 = cos(pkin(6));
	t683 = t635 * t706;
	t705 = sin(pkin(12));
	t621 = -t705 * t640 + t643 * t683;
	t711 = r_i_i_C(1) + pkin(5);
	t710 = r_i_i_C(2) + pkin(11);
	t709 = cos(qJ(3));
	t633 = sin(pkin(7));
	t708 = pkin(9) * t633;
	t707 = r_i_i_C(3) + qJ(6);
	t638 = sin(qJ(4));
	t704 = t633 * t638;
	t642 = cos(qJ(4));
	t703 = t633 * t642;
	t702 = t633 * t643;
	t634 = sin(pkin(6));
	t701 = t634 * t635;
	t636 = cos(pkin(7));
	t639 = sin(qJ(3));
	t700 = t636 * t639;
	t637 = sin(qJ(5));
	t699 = t637 * t642;
	t698 = t639 * t640;
	t697 = t639 * t643;
	t696 = qJD(2) * t634;
	t695 = qJD(4) * t638;
	t694 = qJD(5) * t642;
	t693 = t633 * t701;
	t692 = t633 * t634 * t640;
	t691 = t634 * t698;
	t652 = -t640 * t683 - t705 * t643;
	t690 = t652 * t709;
	t689 = t636 * t709;
	t688 = t709 * t640;
	t687 = t709 * t643;
	t686 = t633 * t696;
	t685 = t633 * t706;
	t684 = t634 * t705;
	t681 = t636 * t687;
	t680 = t640 * t686;
	t679 = t643 * t686;
	t677 = t639 * t685;
	t676 = t633 * t684;
	t616 = t621 * qJD(2);
	t617 = t652 * qJD(2);
	t645 = t621 * t689 + t639 * t652 - t709 * t693;
	t575 = t645 * qJD(3) + t616 * t709 + t617 * t700;
	t675 = -t645 * t694 + t575;
	t672 = t706 * t705;
	t650 = t635 * t640 + t643 * t672;
	t618 = t650 * qJD(2);
	t651 = -t635 * t643 + t640 * t672;
	t619 = t651 * qJD(2);
	t644 = t639 * t651 - t650 * t689 + t709 * t676;
	t577 = t644 * qJD(3) - t618 * t709 + t619 * t700;
	t674 = -t644 * t694 + t577;
	t657 = -t636 * t698 + t687;
	t661 = t709 * t685;
	t592 = qJD(3) * t661 + ((t681 - t698) * qJD(3) + t657 * qJD(2)) * t634;
	t671 = t634 * t681;
	t606 = -t661 - t671 + t691;
	t673 = t606 * t694 + t592;
	t594 = -t690 + (t621 * t636 - t693) * t639;
	t608 = -t621 * t633 - t636 * t701;
	t585 = t594 * t642 + t608 * t638;
	t641 = cos(qJ(5));
	t670 = t585 * t641 - t637 * t645;
	t596 = -t651 * t709 + (-t636 * t650 + t676) * t639;
	t609 = t633 * t650 + t636 * t684;
	t587 = t596 * t642 + t609 * t638;
	t669 = t587 * t641 - t637 * t644;
	t602 = t621 * t709 + t652 * t700;
	t588 = t602 * t642 - t652 * t704;
	t659 = -t621 * t639 + t652 * t689;
	t668 = -t588 * t637 - t641 * t659;
	t604 = -t650 * t709 + t651 * t700;
	t589 = t604 * t642 - t651 * t704;
	t658 = t639 * t650 + t651 * t689;
	t667 = -t589 * t637 - t641 * t658;
	t666 = -t594 * t638 + t608 * t642;
	t665 = -t596 * t638 + t609 * t642;
	t656 = t636 * t697 + t688;
	t607 = t656 * t634 + t677;
	t620 = -t634 * t702 + t706 * t636;
	t598 = t607 * t642 + t620 * t638;
	t664 = t598 * t641 + t606 * t637;
	t615 = t657 * t634;
	t605 = t615 * t642 + t638 * t692;
	t655 = t636 * t688 + t697;
	t614 = t655 * t634;
	t663 = -t605 * t637 + t614 * t641;
	t662 = -t607 * t638 + t620 * t642;
	t660 = -t642 * pkin(4) - t710 * t638 - pkin(3);
	t654 = qJD(4) * (pkin(4) * t638 - t710 * t642);
	t653 = t707 * t637 + t711 * t641 + pkin(4);
	t574 = t616 * t639 - t617 * t689 + (t621 * t700 - t639 * t693 - t690) * qJD(3);
	t649 = qJD(5) * t594 - t574 * t642 - t645 * t695;
	t576 = t596 * qJD(3) - t618 * t639 - t619 * t689;
	t648 = qJD(5) * t596 - t576 * t642 - t644 * t695;
	t591 = qJD(3) * t677 + (t655 * qJD(2) + t656 * qJD(3)) * t634;
	t647 = qJD(5) * t607 - t591 * t642 + t606 * t695;
	t646 = qJD(6) * t637 + (-t711 * t637 + t707 * t641) * qJD(5);
	t600 = (-t656 * qJD(2) - t655 * qJD(3)) * t634;
	t599 = -qJD(2) * t671 - t634 * qJD(3) * t687 + (qJD(3) * t636 + qJD(2)) * t691;
	t583 = t658 * qJD(3) + t618 * t700 + t619 * t709;
	t582 = t604 * qJD(3) - t618 * t689 + t619 * t639;
	t581 = t659 * qJD(3) - t616 * t700 + t617 * t709;
	t580 = t602 * qJD(3) + t616 * t689 + t617 * t639;
	t579 = t638 * t679 + t600 * t642 + (-t615 * t638 + t642 * t692) * qJD(4);
	t571 = t662 * qJD(4) + t592 * t642 + t638 * t680;
	t569 = -t618 * t704 + t583 * t642 + (-t604 * t638 - t651 * t703) * qJD(4);
	t567 = t616 * t704 + t581 * t642 + (-t602 * t638 - t652 * t703) * qJD(4);
	t563 = t665 * qJD(4) + t577 * t642 - t619 * t704;
	t561 = t666 * qJD(4) + t575 * t642 - t617 * t704;
	t556 = t664 * qJD(5) + t571 * t637 - t591 * t641;
	t546 = t669 * qJD(5) + t563 * t637 - t576 * t641;
	t544 = t670 * qJD(5) + t561 * t637 - t574 * t641;
	t1 = [0, -t667 * qJD(6) + t569 * pkin(4) + t583 * pkin(3) + t582 * pkin(10) + t619 * pkin(2) - t618 * t708 + t710 * (t589 * qJD(4) + t583 * t638 + t618 * t703) + t711 * (t667 * qJD(5) + t569 * t641 + t582 * t637) + t707 * (t569 * t637 - t582 * t641 + (t589 * t641 - t637 * t658) * qJD(5)), -(t596 * t641 - t644 * t699) * qJD(6) + t577 * pkin(10) + t711 * (t674 * t637 + t648 * t641) + t707 * (t648 * t637 - t674 * t641) - t644 * t654 + t660 * t576, t710 * t563 + t646 * t665 + t653 * (-t587 * qJD(4) - t577 * t638 - t619 * t703), t669 * qJD(6) + t707 * (t563 * t641 + t576 * t637 + (-t587 * t637 - t641 * t644) * qJD(5)) - t711 * t546, t546; 0, -t668 * qJD(6) + t567 * pkin(4) + t581 * pkin(3) + t580 * pkin(10) + t617 * pkin(2) + t616 * t708 + t710 * (t588 * qJD(4) + t581 * t638 - t616 * t703) + t711 * (t668 * qJD(5) + t567 * t641 + t580 * t637) + t707 * (t567 * t637 - t580 * t641 + (t588 * t641 - t637 * t659) * qJD(5)), -(t594 * t641 - t645 * t699) * qJD(6) + t575 * pkin(10) + t711 * (t675 * t637 + t649 * t641) + t707 * (t649 * t637 - t675 * t641) - t645 * t654 + t660 * t574, t710 * t561 + t646 * t666 + t653 * (-t585 * qJD(4) - t575 * t638 - t617 * t703), t670 * qJD(6) + t707 * (t561 * t641 + t574 * t637 + (-t585 * t637 - t641 * t645) * qJD(5)) - t711 * t544, t544; 0, -t663 * qJD(6) + t579 * pkin(4) + t600 * pkin(3) - t599 * pkin(10) + t710 * (t605 * qJD(4) + t600 * t638 - t642 * t679) + t711 * (t663 * qJD(5) + t579 * t641 - t599 * t637) + t707 * (t579 * t637 + t599 * t641 + (t605 * t641 + t614 * t637) * qJD(5)) + (-pkin(2) * t640 + pkin(9) * t702) * t696, -(t606 * t699 + t607 * t641) * qJD(6) + t592 * pkin(10) + t711 * (t673 * t637 + t647 * t641) + t707 * (t647 * t637 - t673 * t641) + t606 * t654 + t660 * t591, t710 * t571 + t646 * t662 + t653 * (-t598 * qJD(4) - t592 * t638 + t642 * t680), t664 * qJD(6) + t707 * (t571 * t641 + t591 * t637 + (-t598 * t637 + t606 * t641) * qJD(5)) - t711 * t556, t556;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end