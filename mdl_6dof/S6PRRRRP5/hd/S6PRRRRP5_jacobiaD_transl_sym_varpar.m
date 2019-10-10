% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRRP5
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
% Datum: 2019-10-09 23:09
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRP5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:09:35
	% EndTime: 2019-10-09 23:09:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:09:35
	% EndTime: 2019-10-09 23:09:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:09:35
	% EndTime: 2019-10-09 23:09:35
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
	% StartTime: 2019-10-09 23:09:36
	% EndTime: 2019-10-09 23:09:36
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
	% StartTime: 2019-10-09 23:09:37
	% EndTime: 2019-10-09 23:09:38
	% DurationCPUTime: 0.82s
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
	% StartTime: 2019-10-09 23:09:40
	% EndTime: 2019-10-09 23:09:42
	% DurationCPUTime: 1.33s
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
	t512 = (-t564 * qJD(2) - t563 * qJD(3)) * t546;
	t511 = -qJD(2) * t574 - t546 * qJD(3) * t588 + (qJD(3) * t548 + qJD(2)) * t592;
	t504 = qJD(3) * t569 + ((t582 - t599) * qJD(3) + t565 * qJD(2)) * t546;
	t503 = qJD(3) * t578 + (t563 * qJD(2) + t564 * qJD(3)) * t546;
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
	t1 = [0, (t484 * t553 + t495 * t549) * r_i_i_C(1) + (-t484 * t549 + t495 * t553) * r_i_i_C(2) + t484 * pkin(4) + t496 * pkin(3) + t495 * pkin(10) + t531 * pkin(2) - t530 * t607 + t609 * (t502 * qJD(4) + t496 * t550 + t530 * t603) + ((-t502 * t549 - t553 * t566) * r_i_i_C(1) + (-t502 * t553 + t549 * t566) * r_i_i_C(2)) * qJD(5), (t490 * t549 + t508 * t595) * r_i_i_C(1) + (t490 * t553 - t508 * t596) * r_i_i_C(2) + t490 * pkin(10) + t559 * t489 - t556 * t557, t609 * t480 - t572 * t568 + t570 * (-t500 * qJD(4) - t490 * t550 - t531 * t603), (-t480 * t549 + t489 * t553) * r_i_i_C(1) + (-t480 * t553 - t489 * t549) * r_i_i_C(2) + ((-t500 * t553 + t549 * t557) * r_i_i_C(1) + (t500 * t549 + t553 * t557) * r_i_i_C(2)) * qJD(5), 0; 0, (t482 * t553 + t493 * t549) * r_i_i_C(1) + (-t482 * t549 + t493 * t553) * r_i_i_C(2) + t482 * pkin(4) + t494 * pkin(3) + t493 * pkin(10) + t529 * pkin(2) + t528 * t607 + t609 * (t501 * qJD(4) + t494 * t550 - t528 * t603) + ((-t501 * t549 - t553 * t567) * r_i_i_C(1) + (-t501 * t553 + t549 * t567) * r_i_i_C(2)) * qJD(5), (t488 * t549 + t506 * t595) * r_i_i_C(1) + (t488 * t553 - t506 * t596) * r_i_i_C(2) + t488 * pkin(10) + t559 * t487 - t556 * t558, t609 * t478 - t573 * t568 + t570 * (-t498 * qJD(4) - t488 * t550 - t529 * t603), (-t478 * t549 + t487 * t553) * r_i_i_C(1) + (-t478 * t553 - t487 * t549) * r_i_i_C(2) + ((-t498 * t553 + t549 * t558) * r_i_i_C(1) + (t498 * t549 + t553 * t558) * r_i_i_C(2)) * qJD(5), 0; 0, (t492 * t553 - t511 * t549) * r_i_i_C(1) + (-t492 * t549 - t511 * t553) * r_i_i_C(2) + t492 * pkin(4) + t512 * pkin(3) - t511 * pkin(10) + t609 * (t517 * qJD(4) + t512 * t550 - t554 * t580) + ((-t517 * t549 + t526 * t553) * r_i_i_C(1) + (-t517 * t553 - t526 * t549) * r_i_i_C(2)) * qJD(5) + (-pkin(2) * t552 + pkin(9) * t602) * t597, (t504 * t549 + t519 * t595) * r_i_i_C(1) + (t504 * t553 - t519 * t596) * r_i_i_C(2) + t504 * pkin(10) + t559 * t503 + t556 * t518, t609 * t486 - t571 * t568 + t570 * (-t510 * qJD(4) - t504 * t550 + t554 * t581), (-t486 * t549 + t503 * t553) * r_i_i_C(1) + (-t486 * t553 - t503 * t549) * r_i_i_C(2) + ((-t510 * t553 - t518 * t549) * r_i_i_C(1) + (t510 * t549 - t518 * t553) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:09:40
	% EndTime: 2019-10-09 23:09:42
	% DurationCPUTime: 1.85s
	% Computational Cost: add. (1491->182), mult. (4729->314), div. (0->0), fcn. (5265->14), ass. (0->115)
	t559 = sin(qJ(5));
	t563 = cos(qJ(5));
	t638 = pkin(5) + r_i_i_C(1);
	t586 = t563 * r_i_i_C(2) + t638 * t559;
	t580 = pkin(10) + t586;
	t575 = qJD(5) * t586;
	t562 = sin(qJ(2));
	t565 = cos(qJ(2));
	t556 = cos(pkin(12));
	t633 = cos(pkin(6));
	t611 = t556 * t633;
	t632 = sin(pkin(12));
	t574 = t632 * t562 - t565 * t611;
	t637 = cos(qJ(3));
	t554 = sin(pkin(7));
	t636 = pkin(9) * t554;
	t634 = r_i_i_C(3) + qJ(6) + pkin(11);
	t560 = sin(qJ(4));
	t631 = t554 * t560;
	t564 = cos(qJ(4));
	t630 = t554 * t564;
	t629 = t554 * t565;
	t555 = sin(pkin(6));
	t628 = t555 * t556;
	t557 = cos(pkin(7));
	t561 = sin(qJ(3));
	t627 = t557 * t561;
	t626 = t561 * t562;
	t625 = t561 * t565;
	t624 = qJD(2) * t555;
	t623 = qJD(5) * t559;
	t622 = qJD(5) * t563;
	t621 = t554 * t555 * t562;
	t620 = t555 * t626;
	t619 = t554 * t628;
	t573 = -t562 * t611 - t632 * t565;
	t618 = t573 * t637;
	t617 = t557 * t637;
	t616 = t637 * t562;
	t615 = t637 * t565;
	t614 = t554 * t624;
	t613 = t554 * t633;
	t612 = t555 * t632;
	t609 = t557 * t615;
	t608 = t562 * t614;
	t607 = t565 * t614;
	t605 = t561 * t613;
	t604 = t554 * t612;
	t603 = t633 * t632;
	t602 = t555 * t609;
	t568 = t574 * t637;
	t514 = t557 * t568 - t561 * t573 + t637 * t619;
	t537 = t574 * qJD(2);
	t538 = t573 * qJD(2);
	t497 = -t514 * qJD(3) - t537 * t637 + t538 * t627;
	t515 = -t618 + (-t574 * t557 - t619) * t561;
	t529 = t574 * t554 - t557 * t628;
	t594 = -t515 * t560 + t529 * t564;
	t487 = t594 * qJD(4) + t497 * t564 - t538 * t631;
	t570 = t574 * t561;
	t496 = -t537 * t561 - t538 * t617 + (-t557 * t570 - t561 * t619 - t618) * qJD(3);
	t601 = -t487 * t559 + t496 * t563;
	t571 = t556 * t562 + t565 * t603;
	t569 = t571 * t637;
	t572 = -t556 * t565 + t562 * t603;
	t516 = t557 * t569 - t561 * t572 - t637 * t604;
	t539 = t571 * qJD(2);
	t540 = t572 * qJD(2);
	t499 = -t516 * qJD(3) - t539 * t637 + t540 * t627;
	t517 = -t572 * t637 + (-t571 * t557 + t604) * t561;
	t530 = t571 * t554 + t557 * t612;
	t593 = -t517 * t560 + t530 * t564;
	t489 = t593 * qJD(4) + t499 * t564 - t540 * t631;
	t498 = t517 * qJD(3) - t539 * t561 - t540 * t617;
	t600 = -t489 * t559 + t498 * t563;
	t578 = -t557 * t626 + t615;
	t589 = t637 * t613;
	t513 = qJD(3) * t589 + ((t609 - t626) * qJD(3) + t578 * qJD(2)) * t555;
	t577 = t557 * t625 + t616;
	t528 = t577 * t555 + t605;
	t541 = -t555 * t629 + t633 * t557;
	t590 = -t528 * t560 + t541 * t564;
	t495 = t590 * qJD(4) + t513 * t564 + t560 * t608;
	t576 = t557 * t616 + t625;
	t512 = qJD(3) * t605 + (t576 * qJD(2) + t577 * qJD(3)) * t555;
	t599 = -t495 * t559 + t512 * t563;
	t507 = t515 * t564 + t529 * t560;
	t598 = -t507 * t563 - t514 * t559;
	t509 = t517 * t564 + t530 * t560;
	t597 = -t509 * t563 - t516 * t559;
	t519 = t528 * t564 + t541 * t560;
	t527 = -t589 - t602 + t620;
	t592 = -t519 * t563 - t527 * t559;
	t588 = -r_i_i_C(2) * t559 + t638 * t563 + pkin(4);
	t523 = t573 * t627 - t568;
	t585 = -t523 * t560 - t573 * t630;
	t510 = t523 * t564 - t573 * t631;
	t525 = t572 * t627 - t569;
	t584 = -t525 * t560 - t572 * t630;
	t511 = t525 * t564 - t572 * t631;
	t536 = t578 * t555;
	t579 = -t536 * t560 + t564 * t621;
	t526 = t536 * t564 + t560 * t621;
	t567 = -t634 * t560 - t588 * t564 - pkin(3);
	t524 = -t571 * t561 - t572 * t617;
	t522 = -t573 * t617 - t570;
	t566 = -t560 * qJD(6) + t564 * t575 + (t588 * t560 - t634 * t564) * qJD(4);
	t535 = t576 * t555;
	t521 = (-t577 * qJD(2) - t576 * qJD(3)) * t555;
	t505 = -t524 * qJD(3) + t539 * t627 + t540 * t637;
	t503 = -t522 * qJD(3) + t537 * t627 + t538 * t637;
	t494 = t519 * qJD(4) + t513 * t560 - t564 * t608;
	t488 = t509 * qJD(4) + t499 * t560 + t540 * t630;
	t486 = t507 * qJD(4) + t497 * t560 + t538 * t630;
	t1 = [0, -t584 * qJD(6) + t505 * pkin(3) + t540 * pkin(2) - t539 * t636 + t588 * (t584 * qJD(4) + t505 * t564 - t539 * t631) + t580 * (t525 * qJD(3) - t539 * t617 + t540 * t561) + t634 * (t511 * qJD(4) + t505 * t560 + t539 * t630) + ((-t511 * t563 - t524 * t559) * r_i_i_C(2) + t638 * (-t511 * t559 + t524 * t563)) * qJD(5), (t499 * t563 - t517 * t623) * r_i_i_C(2) + t499 * pkin(10) + t567 * t498 + t566 * t516 + t638 * (t499 * t559 + t517 * t622), qJD(6) * t509 - t588 * t488 + t634 * t489 - t593 * t575, t600 * r_i_i_C(1) + (-t489 * t563 - t498 * t559) * r_i_i_C(2) + (t597 * r_i_i_C(1) + (t509 * t559 - t516 * t563) * r_i_i_C(2)) * qJD(5) + (t597 * qJD(5) + t600) * pkin(5), t488; 0, -t585 * qJD(6) + t503 * pkin(3) + t538 * pkin(2) - t537 * t636 + t588 * (t585 * qJD(4) + t503 * t564 - t537 * t631) + t580 * (t523 * qJD(3) - t537 * t617 + t538 * t561) + t634 * (t510 * qJD(4) + t503 * t560 + t537 * t630) + ((-t510 * t563 - t522 * t559) * r_i_i_C(2) + t638 * (-t510 * t559 + t522 * t563)) * qJD(5), (t497 * t563 - t515 * t623) * r_i_i_C(2) + t497 * pkin(10) + t567 * t496 + t566 * t514 + t638 * (t497 * t559 + t515 * t622), qJD(6) * t507 - t588 * t486 + t634 * t487 - t594 * t575, t601 * r_i_i_C(1) + (-t487 * t563 - t496 * t559) * r_i_i_C(2) + (t598 * r_i_i_C(1) + (t507 * t559 - t514 * t563) * r_i_i_C(2)) * qJD(5) + (t598 * qJD(5) + t601) * pkin(5), t486; 0, -t579 * qJD(6) + t521 * pkin(3) + t588 * (t579 * qJD(4) + t521 * t564 + t560 * t607) - t580 * (-qJD(2) * t602 - t555 * qJD(3) * t615 + (qJD(3) * t557 + qJD(2)) * t620) + t634 * (t526 * qJD(4) + t521 * t560 - t564 * t607) + (-pkin(2) * t562 + pkin(9) * t629) * t624 + ((-t526 * t563 - t535 * t559) * r_i_i_C(2) + t638 * (-t526 * t559 + t535 * t563)) * qJD(5), (t513 * t563 - t528 * t623) * r_i_i_C(2) + t513 * pkin(10) + t567 * t512 + t566 * t527 + t638 * (t513 * t559 + t528 * t622), qJD(6) * t519 - t588 * t494 + t634 * t495 - t590 * t575, t599 * r_i_i_C(1) + (-t495 * t563 - t512 * t559) * r_i_i_C(2) + (t592 * r_i_i_C(1) + (t519 * t559 - t527 * t563) * r_i_i_C(2)) * qJD(5) + (t592 * qJD(5) + t599) * pkin(5), t494;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end