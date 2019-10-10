% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRRR7
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:05
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobiaD_transl_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
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
	% StartTime: 2019-10-09 22:05:03
	% EndTime: 2019-10-09 22:05:03
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (41->25), mult. (151->54), div. (0->0), fcn. (138->10), ass. (0->24)
	t193 = r_i_i_C(3) + qJ(3);
	t175 = sin(pkin(14));
	t181 = cos(pkin(7));
	t192 = t175 * t181;
	t177 = sin(pkin(7));
	t183 = sin(qJ(2));
	t191 = t177 * t183;
	t179 = cos(pkin(14));
	t190 = t179 * t181;
	t184 = cos(qJ(2));
	t189 = t181 * t184;
	t182 = cos(pkin(6));
	t188 = t182 * t183;
	t187 = t182 * t184;
	t176 = sin(pkin(13));
	t180 = cos(pkin(13));
	t186 = -t176 * t184 - t180 * t188;
	t185 = t176 * t188 - t180 * t184;
	t178 = sin(pkin(6));
	t174 = t185 * qJD(2);
	t173 = (t176 * t187 + t180 * t183) * qJD(2);
	t172 = t186 * qJD(2);
	t171 = (t176 * t183 - t180 * t187) * qJD(2);
	t1 = [0, (t173 * t192 + t174 * t179) * r_i_i_C(1) + (t173 * t190 - t174 * t175) * r_i_i_C(2) + t174 * pkin(2) + (-t185 * qJD(3) - t193 * t173) * t177, -t174 * t177, 0, 0, 0; 0, (t171 * t192 + t172 * t179) * r_i_i_C(1) + (t171 * t190 - t172 * t175) * r_i_i_C(2) + t172 * pkin(2) + (-t186 * qJD(3) - t193 * t171) * t177, -t172 * t177, 0, 0, 0; 0, (qJD(3) * t191 + ((-t175 * t189 - t179 * t183) * r_i_i_C(1) + (t175 * t183 - t179 * t189) * r_i_i_C(2) - t183 * pkin(2) + t193 * t184 * t177) * qJD(2)) * t178, t178 * qJD(2) * t191, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:04
	% EndTime: 2019-10-09 22:05:05
	% DurationCPUTime: 0.90s
	% Computational Cost: add. (271->109), mult. (972->215), div. (0->0), fcn. (1026->14), ass. (0->74)
	t373 = cos(pkin(8));
	t412 = pkin(10) + r_i_i_C(3);
	t413 = t412 * t373 + qJ(3);
	t366 = sin(pkin(14));
	t374 = cos(pkin(7));
	t411 = t366 * t374;
	t368 = sin(pkin(8));
	t376 = sin(qJ(4));
	t410 = t368 * t376;
	t378 = cos(qJ(4));
	t409 = t368 * t378;
	t369 = sin(pkin(7));
	t370 = sin(pkin(6));
	t408 = t369 * t370;
	t375 = cos(pkin(6));
	t407 = t369 * t375;
	t406 = t370 * t374;
	t371 = cos(pkin(14));
	t405 = t371 * t374;
	t377 = sin(qJ(2));
	t404 = t371 * t377;
	t403 = t373 * t376;
	t402 = t373 * t378;
	t401 = t374 * t377;
	t379 = cos(qJ(2));
	t400 = t374 * t379;
	t399 = t375 * t377;
	t398 = t375 * t379;
	t396 = qJD(2) * t377;
	t395 = qJD(2) * t379;
	t394 = qJD(4) * t376;
	t393 = qJD(4) * t378;
	t392 = t369 * t410;
	t391 = t369 * t409;
	t372 = cos(pkin(13));
	t390 = t372 * t398;
	t389 = t396 * t408;
	t387 = t368 * t389;
	t386 = t378 * r_i_i_C(1) - t376 * r_i_i_C(2) + pkin(3);
	t367 = sin(pkin(13));
	t360 = -t367 * t377 + t390;
	t385 = t360 * t374 - t372 * t408;
	t382 = t367 * t398 + t372 * t377;
	t384 = t367 * t408 - t374 * t382;
	t383 = -t366 * t377 + t371 * t400;
	t361 = t367 * t379 + t372 * t399;
	t381 = t367 * t399 - t372 * t379;
	t353 = (-t366 * t379 - t371 * t401) * t370;
	t354 = (-t366 * t401 + t371 * t379) * t370;
	t380 = (t376 * r_i_i_C(1) + t378 * r_i_i_C(2)) * t373 - t412 * t368;
	t359 = t375 * t374 - t379 * t408;
	t358 = t381 * qJD(2);
	t357 = t382 * qJD(2);
	t356 = t361 * qJD(2);
	t355 = -qJD(2) * t390 + t367 * t396;
	t352 = qJD(2) * t354;
	t351 = qJD(2) * t353;
	t348 = t367 * t406 + t369 * t382;
	t347 = -t360 * t369 - t372 * t406;
	t346 = t370 * t404 + (t370 * t400 + t407) * t366;
	t345 = t383 * t370 + t371 * t407;
	t344 = -t371 * t382 + t381 * t411;
	t343 = t366 * t382 + t381 * t405;
	t342 = t360 * t371 - t361 * t411;
	t341 = -t360 * t366 - t361 * t405;
	t338 = -t357 * t371 + t358 * t411;
	t337 = t357 * t366 + t358 * t405;
	t334 = -t355 * t371 - t356 * t411;
	t333 = t355 * t366 - t356 * t405;
	t332 = t366 * t384 - t371 * t381;
	t331 = t366 * t381 + t371 * t384;
	t330 = t361 * t371 + t366 * t385;
	t329 = -t361 * t366 + t371 * t385;
	t1 = [0, t358 * pkin(2) + t386 * (t357 * t411 + t358 * t371) + t380 * (t357 * t405 - t358 * t366) + ((t343 * t402 - t344 * t376) * r_i_i_C(1) + (-t343 * t403 - t344 * t378) * r_i_i_C(2)) * qJD(4) + (-t381 * qJD(3) - t413 * t357 + ((-t357 * t376 - t381 * t393) * r_i_i_C(1) + (-t357 * t378 + t381 * t394) * r_i_i_C(2)) * t368) * t369, -t358 * t369, (t337 * t402 - t338 * t376 - t358 * t391) * r_i_i_C(1) + (-t337 * t403 - t338 * t378 + t358 * t392) * r_i_i_C(2) + ((-t331 * t403 - t332 * t378 - t348 * t410) * r_i_i_C(1) + (-t331 * t402 + t332 * t376 - t348 * t409) * r_i_i_C(2)) * qJD(4), 0, 0; 0, -t356 * pkin(2) + t386 * (t355 * t411 - t356 * t371) + t380 * (t355 * t405 + t356 * t366) + ((t341 * t402 - t342 * t376) * r_i_i_C(1) + (-t341 * t403 - t342 * t378) * r_i_i_C(2)) * qJD(4) + (t361 * qJD(3) - t413 * t355 + ((-t355 * t376 + t361 * t393) * r_i_i_C(1) + (-t355 * t378 - t361 * t394) * r_i_i_C(2)) * t368) * t369, t356 * t369, (t333 * t402 - t334 * t376 + t356 * t391) * r_i_i_C(1) + (-t333 * t403 - t334 * t378 - t356 * t392) * r_i_i_C(2) + ((-t329 * t403 - t330 * t378 - t347 * t410) * r_i_i_C(1) + (-t329 * t402 + t330 * t376 - t347 * t409) * r_i_i_C(2)) * qJD(4), 0, 0; 0, ((t353 * t402 - t354 * t376) * r_i_i_C(1) + (-t353 * t403 - t354 * t378) * r_i_i_C(2)) * qJD(4) + (-pkin(2) * t396 + (t377 * qJD(3) + t413 * t395 + ((t376 * t395 + t377 * t393) * r_i_i_C(1) + (-t377 * t394 + t378 * t395) * r_i_i_C(2)) * t368) * t369 + (t386 * (-t366 * t400 - t404) - t380 * t383) * qJD(2)) * t370, t389, (t351 * t402 - t352 * t376 + t378 * t387) * r_i_i_C(1) + (-t351 * t403 - t352 * t378 - t376 * t387) * r_i_i_C(2) + ((-t345 * t403 - t346 * t378 - t359 * t410) * r_i_i_C(1) + (-t345 * t402 + t346 * t376 - t359 * t409) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:08
	% EndTime: 2019-10-09 22:05:09
	% DurationCPUTime: 1.30s
	% Computational Cost: add. (1016->167), mult. (3456->319), div. (0->0), fcn. (3884->16), ass. (0->120)
	t583 = sin(pkin(13));
	t588 = cos(pkin(13));
	t594 = sin(qJ(2));
	t591 = cos(pkin(6));
	t597 = cos(qJ(2));
	t627 = t591 * t597;
	t576 = -t583 * t594 + t588 * t627;
	t582 = sin(pkin(14));
	t586 = sin(pkin(6));
	t590 = cos(pkin(7));
	t629 = t590 * t597;
	t587 = cos(pkin(14));
	t631 = t587 * t594;
	t585 = sin(pkin(7));
	t636 = t585 * t591;
	t562 = t586 * t631 + (t586 * t629 + t636) * t582;
	t593 = sin(qJ(4));
	t596 = cos(qJ(4));
	t607 = -t582 * t594 + t587 * t629;
	t561 = t586 * t607 + t587 * t636;
	t634 = t585 * t597;
	t575 = -t586 * t634 + t590 * t591;
	t584 = sin(pkin(8));
	t589 = cos(pkin(8));
	t617 = t561 * t589 + t575 * t584;
	t531 = t562 * t596 + t593 * t617;
	t628 = t591 * t594;
	t605 = t583 * t628 - t588 * t597;
	t606 = t583 * t627 + t588 * t594;
	t638 = t585 * t586;
	t608 = t583 * t638 - t590 * t606;
	t545 = t582 * t608 - t587 * t605;
	t544 = t582 * t605 + t587 * t608;
	t633 = t586 * t590;
	t564 = t583 * t633 + t585 * t606;
	t618 = t544 * t589 + t564 * t584;
	t525 = t545 * t596 + t593 * t618;
	t577 = t583 * t597 + t588 * t628;
	t609 = t576 * t590 - t588 * t638;
	t543 = t577 * t587 + t582 * t609;
	t542 = -t577 * t582 + t587 * t609;
	t563 = -t576 * t585 - t588 * t633;
	t619 = t542 * t589 + t563 * t584;
	t523 = t543 * t596 + t593 * t619;
	t648 = r_i_i_C(3) + pkin(11);
	t571 = t576 * qJD(2);
	t572 = t577 * qJD(2);
	t632 = t587 * t590;
	t548 = -t571 * t632 + t572 * t582;
	t645 = t548 * t584;
	t573 = t606 * qJD(2);
	t574 = t605 * qJD(2);
	t552 = t573 * t632 - t574 * t582;
	t644 = t552 * t584;
	t626 = qJD(2) * t586;
	t565 = t607 * t626;
	t642 = t565 * t584;
	t641 = t582 * t590;
	t639 = t584 * t585;
	t637 = t585 * t589;
	t635 = t585 * t594;
	t630 = t590 * t594;
	t625 = t586 * t635;
	t623 = t585 * t626;
	t622 = pkin(10) * t589 + qJ(3);
	t621 = t597 * t623;
	t620 = t594 * t623;
	t592 = sin(qJ(5));
	t595 = cos(qJ(5));
	t616 = r_i_i_C(1) * t595 - r_i_i_C(2) * t592 + pkin(4);
	t546 = -t571 * t582 - t572 * t632;
	t615 = t546 * t589 + t572 * t639;
	t614 = -t548 * t589 - t571 * t639;
	t550 = t573 * t582 + t574 * t632;
	t613 = t550 * t589 - t574 * t639;
	t612 = -t552 * t589 + t573 * t639;
	t554 = -t576 * t582 - t577 * t632;
	t611 = t554 * t589 + t577 * t639;
	t556 = t582 * t606 + t605 * t632;
	t610 = t556 * t589 - t605 * t639;
	t604 = qJD(5) * (-r_i_i_C(1) * t592 - r_i_i_C(2) * t595);
	t569 = (-t582 * t597 - t587 * t630) * t586;
	t603 = t569 * t589 + t584 * t625;
	t570 = (-t582 * t630 + t587 * t597) * t586;
	t602 = -t565 * t589 + t584 * t621;
	t567 = qJD(2) * t569;
	t601 = t567 * t589 + t584 * t620;
	t600 = -t543 * t593 + t596 * t619;
	t599 = -t545 * t593 + t596 * t618;
	t598 = -t562 * t593 + t596 * t617;
	t555 = t576 * t587 - t577 * t641;
	t528 = t555 * t596 + t593 * t611;
	t557 = -t587 * t606 + t605 * t641;
	t529 = t557 * t596 + t593 * t610;
	t540 = t570 * t596 + t593 * t603;
	t568 = qJD(2) * t570;
	t566 = (-t582 * t629 - t631) * t626;
	t560 = -t569 * t584 + t589 * t625;
	t559 = -t567 * t584 + t589 * t620;
	t558 = t589 * t621 + t642;
	t553 = t573 * t641 + t574 * t587;
	t551 = -t573 * t587 + t574 * t641;
	t549 = -t571 * t641 - t572 * t587;
	t547 = t571 * t587 - t572 * t641;
	t541 = -t561 * t584 + t575 * t589;
	t539 = -t556 * t584 - t605 * t637;
	t538 = -t554 * t584 + t577 * t637;
	t537 = -t573 * t637 - t644;
	t536 = -t550 * t584 - t574 * t637;
	t535 = t571 * t637 - t645;
	t534 = -t546 * t584 + t572 * t637;
	t533 = -t544 * t584 + t564 * t589;
	t532 = -t542 * t584 + t563 * t589;
	t527 = t566 * t596 + t602 * t593 + (-t570 * t593 + t596 * t603) * qJD(4);
	t521 = qJD(4) * t598 + t568 * t596 + t593 * t601;
	t519 = t553 * t596 - t612 * t593 + (-t557 * t593 + t596 * t610) * qJD(4);
	t517 = t549 * t596 - t614 * t593 + (-t555 * t593 + t596 * t611) * qJD(4);
	t515 = qJD(4) * t599 + t551 * t596 + t593 * t613;
	t513 = qJD(4) * t600 + t547 * t596 + t593 * t615;
	t1 = [0, (t519 * t595 + t537 * t592) * r_i_i_C(1) + (-t519 * t592 + t537 * t595) * r_i_i_C(2) + t519 * pkin(4) + t553 * pkin(3) - pkin(10) * t644 + t574 * pkin(2) + t648 * (qJD(4) * t529 + t553 * t593 + t596 * t612) + ((-t529 * t592 + t539 * t595) * r_i_i_C(1) + (-t529 * t595 - t539 * t592) * r_i_i_C(2)) * qJD(5) + (-qJD(3) * t605 - t573 * t622) * t585, -t574 * t585, t648 * t515 + t599 * t604 + t616 * (-qJD(4) * t525 - t551 * t593 + t613 * t596), (-t515 * t592 + t536 * t595) * r_i_i_C(1) + (-t515 * t595 - t536 * t592) * r_i_i_C(2) + ((-t525 * t595 - t533 * t592) * r_i_i_C(1) + (t525 * t592 - t533 * t595) * r_i_i_C(2)) * qJD(5), 0; 0, (t517 * t595 + t535 * t592) * r_i_i_C(1) + (-t517 * t592 + t535 * t595) * r_i_i_C(2) + t517 * pkin(4) + t549 * pkin(3) - pkin(10) * t645 - t572 * pkin(2) + t648 * (qJD(4) * t528 + t549 * t593 + t596 * t614) + ((-t528 * t592 + t538 * t595) * r_i_i_C(1) + (-t528 * t595 - t538 * t592) * r_i_i_C(2)) * qJD(5) + (t577 * qJD(3) + t571 * t622) * t585, t572 * t585, t648 * t513 + t600 * t604 + t616 * (-qJD(4) * t523 - t547 * t593 + t615 * t596), (-t513 * t592 + t534 * t595) * r_i_i_C(1) + (-t513 * t595 - t534 * t592) * r_i_i_C(2) + ((-t523 * t595 - t532 * t592) * r_i_i_C(1) + (t523 * t592 - t532 * t595) * r_i_i_C(2)) * qJD(5), 0; 0, (t527 * t595 + t558 * t592) * r_i_i_C(1) + (-t527 * t592 + t558 * t595) * r_i_i_C(2) + t527 * pkin(4) + t566 * pkin(3) + pkin(10) * t642 + t648 * (qJD(4) * t540 + t566 * t593 - t596 * t602) + ((-t540 * t592 + t560 * t595) * r_i_i_C(1) + (-t540 * t595 - t560 * t592) * r_i_i_C(2)) * qJD(5) + (qJD(3) * t635 + (-pkin(2) * t594 + t622 * t634) * qJD(2)) * t586, t620, t648 * t521 + t598 * t604 + t616 * (-qJD(4) * t531 - t568 * t593 + t601 * t596), (-t521 * t592 + t559 * t595) * r_i_i_C(1) + (-t521 * t595 - t559 * t592) * r_i_i_C(2) + ((-t531 * t595 - t541 * t592) * r_i_i_C(1) + (t531 * t592 - t541 * t595) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:14
	% EndTime: 2019-10-09 22:05:17
	% DurationCPUTime: 3.10s
	% Computational Cost: add. (3033->247), mult. (9956->442), div. (0->0), fcn. (11639->18), ass. (0->155)
	t769 = sin(pkin(8));
	t774 = cos(pkin(8));
	t768 = sin(pkin(14));
	t772 = cos(pkin(14));
	t770 = sin(pkin(7));
	t775 = cos(pkin(7));
	t773 = cos(pkin(13));
	t779 = sin(qJ(2));
	t782 = cos(qJ(2));
	t862 = sin(pkin(13));
	t863 = cos(pkin(6));
	t827 = t863 * t862;
	t811 = t773 * t779 + t782 * t827;
	t771 = sin(pkin(6));
	t836 = t771 * t862;
	t795 = t770 * t836 - t811 * t775;
	t808 = -t773 * t782 + t779 * t827;
	t788 = -t768 * t808 - t795 * t772;
	t796 = t811 * t770 + t775 * t836;
	t869 = t796 * t769 - t788 * t774;
	t835 = t773 * t863;
	t810 = t862 * t779 - t782 * t835;
	t853 = t771 * t773;
	t798 = -t770 * t853 - t810 * t775;
	t809 = -t779 * t835 - t862 * t782;
	t790 = -t768 * t809 - t798 * t772;
	t799 = t810 * t770 - t775 * t853;
	t868 = t799 * t769 - t790 * t774;
	t848 = t775 * t782;
	t818 = -t768 * t779 + t772 * t848;
	t834 = t863 * t770;
	t801 = t818 * t771 + t772 * t834;
	t854 = t770 * t782;
	t814 = -t771 * t854 + t863 * t775;
	t866 = t814 * t769 + t801 * t774;
	t864 = cos(qJ(4));
	t871 = t869 * t864;
	t870 = t868 * t864;
	t867 = t866 * t864;
	t776 = sin(qJ(6));
	t780 = cos(qJ(6));
	t816 = qJD(6) * (t776 * r_i_i_C(1) + t780 * r_i_i_C(2));
	t865 = r_i_i_C(3) + pkin(12);
	t757 = t810 * qJD(2);
	t758 = t809 * qJD(2);
	t852 = t772 * t775;
	t734 = t757 * t852 - t758 * t768;
	t861 = t734 * t769;
	t759 = t811 * qJD(2);
	t760 = t808 * qJD(2);
	t737 = t759 * t852 - t760 * t768;
	t860 = t737 * t769;
	t846 = qJD(2) * t771;
	t751 = t818 * t846;
	t859 = t751 * t769;
	t849 = t775 * t779;
	t817 = t768 * t782 + t772 * t849;
	t755 = t817 * t771;
	t858 = t755 * t774;
	t857 = t768 * t775;
	t856 = t770 * t769;
	t855 = t770 * t774;
	t851 = t772 * t779;
	t847 = t779 * t770;
	t778 = sin(qJ(4));
	t845 = qJD(4) * t778;
	t844 = qJD(6) * t776;
	t843 = qJD(6) * t780;
	t842 = t778 * t856;
	t841 = t774 * t847;
	t840 = t774 * t864;
	t839 = t770 * t846;
	t838 = qJD(4) * t864;
	t837 = pkin(10) * t774 + qJ(3);
	t832 = t864 * t856;
	t831 = t782 * t839;
	t830 = t779 * t839;
	t826 = t771 * t832;
	t729 = t798 * t768 - t772 * t809;
	t704 = t729 * t864 + t868 * t778;
	t716 = t790 * t769 + t799 * t774;
	t777 = sin(qJ(5));
	t781 = cos(qJ(5));
	t694 = t704 * t781 + t716 * t777;
	t825 = -t704 * t777 + t716 * t781;
	t730 = t795 * t768 - t772 * t808;
	t706 = t730 * t864 + t869 * t778;
	t717 = t788 * t769 + t796 * t774;
	t696 = t706 * t781 + t717 * t777;
	t824 = -t706 * t777 + t717 * t781;
	t739 = t810 * t768 + t809 * t852;
	t740 = -t810 * t772 + t809 * t857;
	t710 = t740 * t864 + (t739 * t774 - t809 * t856) * t778;
	t722 = -t739 * t769 - t809 * t855;
	t697 = t710 * t781 + t722 * t777;
	t741 = t811 * t768 + t808 * t852;
	t742 = -t811 * t772 + t808 * t857;
	t712 = t742 * t864 + (t741 * t774 - t808 * t856) * t778;
	t723 = -t741 * t769 - t808 * t855;
	t698 = t712 * t781 + t723 * t777;
	t749 = t771 * t851 + (t771 * t848 + t834) * t768;
	t715 = t749 * t864 + t866 * t778;
	t728 = -t801 * t769 + t814 * t774;
	t702 = t715 * t781 + t728 * t777;
	t823 = -t715 * t777 + t728 * t781;
	t756 = (-t768 * t849 + t772 * t782) * t771;
	t725 = t756 * t864 + (t769 * t771 * t847 - t858) * t778;
	t745 = t755 * t769 + t771 * t841;
	t713 = t725 * t781 + t745 * t777;
	t822 = t780 * r_i_i_C(1) - t776 * r_i_i_C(2) + pkin(5);
	t820 = t757 * t768 + t758 * t852;
	t819 = t759 * t768 + t760 * t852;
	t815 = qJD(2) * t826;
	t813 = t820 * t774;
	t812 = t819 * t774;
	t806 = qJD(2) * t858;
	t804 = -t865 * t777 - t822 * t781 - pkin(4);
	t803 = t739 * t840 - t740 * t778 - t809 * t832;
	t802 = t741 * t840 - t742 * t778 - t808 * t832;
	t800 = -t755 * t840 - t756 * t778 + t779 * t826;
	t786 = t781 * t816 + (t822 * t777 - t865 * t781) * qJD(5);
	t753 = qJD(2) * t756;
	t752 = (-t768 * t848 - t851) * t846;
	t744 = (t817 * t769 + t841) * t846;
	t743 = t774 * t831 + t859;
	t738 = t759 * t857 + t760 * t772;
	t736 = -t759 * t772 + t760 * t857;
	t735 = t757 * t857 + t758 * t772;
	t733 = -t757 * t772 + t758 * t857;
	t721 = -t759 * t855 - t860;
	t720 = -t760 * t855 - t819 * t769;
	t719 = -t757 * t855 - t861;
	t718 = -t758 * t855 - t820 * t769;
	t714 = t749 * t778 - t867;
	t708 = t752 * t864 + (-t751 * t774 + t769 * t831) * t778 + t800 * qJD(4);
	t707 = t725 * qJD(4) + t751 * t840 + t752 * t778 - t782 * t815;
	t705 = t730 * t778 - t871;
	t703 = t729 * t778 - t870;
	t700 = -t749 * t845 + t753 * t864 + (t769 * t830 - t806) * t778 + t867 * qJD(4);
	t699 = t749 * t838 + t753 * t778 - t779 * t815 + t864 * t806 + t866 * t845;
	t692 = t738 * t864 + (t737 * t774 - t759 * t856) * t778 + t802 * qJD(4);
	t691 = t712 * qJD(4) - t737 * t840 + t738 * t778 + t759 * t832;
	t690 = t735 * t864 + (t734 * t774 - t757 * t856) * t778 + t803 * qJD(4);
	t689 = t710 * qJD(4) - t734 * t840 + t735 * t778 + t757 * t832;
	t688 = t871 * qJD(4) - t730 * t845 + t736 * t864 - t760 * t842 + t778 * t812;
	t687 = t730 * t838 + t736 * t778 + t760 * t832 - t864 * t812 + t869 * t845;
	t686 = t870 * qJD(4) - t729 * t845 + t733 * t864 - t758 * t842 + t778 * t813;
	t685 = t729 * t838 + t733 * t778 + t758 * t832 - t864 * t813 + t868 * t845;
	t684 = t708 * t781 + t743 * t777 + (-t725 * t777 + t745 * t781) * qJD(5);
	t682 = t823 * qJD(5) + t700 * t781 + t744 * t777;
	t680 = t692 * t781 + t721 * t777 + (-t712 * t777 + t723 * t781) * qJD(5);
	t678 = t690 * t781 + t719 * t777 + (-t710 * t777 + t722 * t781) * qJD(5);
	t676 = t824 * qJD(5) + t688 * t781 + t720 * t777;
	t674 = t825 * qJD(5) + t686 * t781 + t718 * t777;
	t1 = [0, (t680 * t780 + t691 * t776) * r_i_i_C(1) + (-t680 * t776 + t691 * t780) * r_i_i_C(2) + t680 * pkin(5) + t692 * pkin(4) + t691 * pkin(11) + t738 * pkin(3) - pkin(10) * t860 + t760 * pkin(2) + t865 * (t698 * qJD(5) + t692 * t777 - t721 * t781) + ((-t698 * t776 - t780 * t802) * r_i_i_C(1) + (-t698 * t780 + t776 * t802) * r_i_i_C(2)) * qJD(6) + (-qJD(3) * t808 - t837 * t759) * t770, -t760 * t770, (t688 * t776 + t706 * t843) * r_i_i_C(1) + (t688 * t780 - t706 * t844) * r_i_i_C(2) + t688 * pkin(11) + t804 * t687 + t786 * t705, t865 * t676 - t824 * t816 + t822 * (-t696 * qJD(5) - t688 * t777 + t720 * t781), (-t676 * t776 + t687 * t780) * r_i_i_C(1) + (-t676 * t780 - t687 * t776) * r_i_i_C(2) + ((-t696 * t780 - t705 * t776) * r_i_i_C(1) + (t696 * t776 - t705 * t780) * r_i_i_C(2)) * qJD(6); 0, (t678 * t780 + t689 * t776) * r_i_i_C(1) + (-t678 * t776 + t689 * t780) * r_i_i_C(2) + t678 * pkin(5) + t690 * pkin(4) + t689 * pkin(11) + t735 * pkin(3) - pkin(10) * t861 + t758 * pkin(2) + t865 * (t697 * qJD(5) + t690 * t777 - t719 * t781) + ((-t697 * t776 - t780 * t803) * r_i_i_C(1) + (-t697 * t780 + t776 * t803) * r_i_i_C(2)) * qJD(6) + (-qJD(3) * t809 - t837 * t757) * t770, -t758 * t770, (t686 * t776 + t704 * t843) * r_i_i_C(1) + (t686 * t780 - t704 * t844) * r_i_i_C(2) + t686 * pkin(11) + t804 * t685 + t786 * t703, t865 * t674 - t825 * t816 + t822 * (-t694 * qJD(5) - t686 * t777 + t718 * t781), (-t674 * t776 + t685 * t780) * r_i_i_C(1) + (-t674 * t780 - t685 * t776) * r_i_i_C(2) + ((-t694 * t780 - t703 * t776) * r_i_i_C(1) + (t694 * t776 - t703 * t780) * r_i_i_C(2)) * qJD(6); 0, (t684 * t780 + t707 * t776) * r_i_i_C(1) + (-t684 * t776 + t707 * t780) * r_i_i_C(2) + t684 * pkin(5) + t708 * pkin(4) + t707 * pkin(11) + t752 * pkin(3) + pkin(10) * t859 + t865 * (t713 * qJD(5) + t708 * t777 - t743 * t781) + ((-t713 * t776 - t780 * t800) * r_i_i_C(1) + (-t713 * t780 + t776 * t800) * r_i_i_C(2)) * qJD(6) + (qJD(3) * t847 + (-pkin(2) * t779 + t837 * t854) * qJD(2)) * t771, t830, (t700 * t776 + t715 * t843) * r_i_i_C(1) + (t700 * t780 - t715 * t844) * r_i_i_C(2) + t700 * pkin(11) + t804 * t699 + t786 * t714, t865 * t682 - t823 * t816 + t822 * (-t702 * qJD(5) - t700 * t777 + t744 * t781), (-t682 * t776 + t699 * t780) * r_i_i_C(1) + (-t682 * t780 - t699 * t776) * r_i_i_C(2) + ((-t702 * t780 - t714 * t776) * r_i_i_C(1) + (t702 * t776 - t714 * t780) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end