% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPPRRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPPRRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPPRRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobiaD_transl_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:22
	% EndTime: 2019-10-10 08:49:22
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (65->32), mult. (208->71), div. (0->0), fcn. (252->14), ass. (0->40)
	t181 = sin(pkin(12));
	t190 = cos(pkin(6));
	t204 = t181 * t190;
	t182 = sin(pkin(8));
	t191 = sin(qJ(4));
	t203 = t182 * t191;
	t192 = cos(qJ(4));
	t202 = t182 * t192;
	t183 = sin(pkin(7));
	t184 = sin(pkin(6));
	t201 = t183 * t184;
	t200 = t183 * t190;
	t189 = cos(pkin(7));
	t199 = t184 * t189;
	t186 = cos(pkin(13));
	t198 = t186 * t189;
	t187 = cos(pkin(12));
	t197 = t187 * t190;
	t188 = cos(pkin(8));
	t196 = t188 * t191;
	t195 = t188 * t192;
	t180 = sin(pkin(13));
	t175 = -t181 * t180 + t186 * t197;
	t194 = t175 * t189 - t187 * t201;
	t177 = -t187 * t180 - t186 * t204;
	t193 = t177 * t189 + t181 * t201;
	t185 = cos(pkin(14));
	t179 = sin(pkin(14));
	t178 = -t180 * t204 + t187 * t186;
	t176 = t180 * t197 + t181 * t186;
	t174 = -t186 * t201 + t190 * t189;
	t173 = -t177 * t183 + t181 * t199;
	t172 = -t175 * t183 - t187 * t199;
	t171 = t184 * t180 * t185 + (t184 * t198 + t200) * t179;
	t170 = t185 * t200 + (-t179 * t180 + t185 * t198) * t184;
	t169 = t178 * t185 + t193 * t179;
	t168 = -t178 * t179 + t193 * t185;
	t167 = t176 * t185 + t194 * t179;
	t166 = -t176 * t179 + t194 * t185;
	t1 = [0, 0, 0, ((-t168 * t196 - t169 * t192 - t173 * t203) * r_i_i_C(1) + (-t168 * t195 + t169 * t191 - t173 * t202) * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, 0, ((-t166 * t196 - t167 * t192 - t172 * t203) * r_i_i_C(1) + (-t166 * t195 + t167 * t191 - t172 * t202) * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, 0, ((-t170 * t196 - t171 * t192 - t174 * t203) * r_i_i_C(1) + (-t170 * t195 + t171 * t191 - t174 * t202) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:24
	% EndTime: 2019-10-10 08:49:24
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (410->57), mult. (1254->116), div. (0->0), fcn. (1562->16), ass. (0->58)
	t391 = sin(pkin(14));
	t392 = sin(pkin(13));
	t396 = sin(pkin(6));
	t397 = cos(pkin(14));
	t398 = cos(pkin(13));
	t401 = cos(pkin(7));
	t420 = t398 * t401;
	t395 = sin(pkin(7));
	t402 = cos(pkin(6));
	t422 = t395 * t402;
	t382 = t396 * t392 * t397 + (t396 * t420 + t422) * t391;
	t404 = sin(qJ(4));
	t406 = cos(qJ(4));
	t381 = t397 * t422 + (-t391 * t392 + t397 * t420) * t396;
	t423 = t395 * t396;
	t386 = -t398 * t423 + t402 * t401;
	t394 = sin(pkin(8));
	t400 = cos(pkin(8));
	t411 = t381 * t400 + t386 * t394;
	t368 = t382 * t406 + t411 * t404;
	t399 = cos(pkin(12));
	t393 = sin(pkin(12));
	t424 = t393 * t402;
	t390 = -t392 * t424 + t399 * t398;
	t389 = -t399 * t392 - t398 * t424;
	t409 = t389 * t401 + t393 * t423;
	t376 = t390 * t397 + t409 * t391;
	t375 = -t390 * t391 + t409 * t397;
	t421 = t396 * t401;
	t384 = -t389 * t395 + t393 * t421;
	t412 = t375 * t400 + t384 * t394;
	t364 = t376 * t406 + t412 * t404;
	t419 = t399 * t402;
	t388 = t392 * t419 + t393 * t398;
	t387 = -t393 * t392 + t398 * t419;
	t410 = t387 * t401 - t399 * t423;
	t374 = t388 * t397 + t410 * t391;
	t373 = -t388 * t391 + t410 * t397;
	t383 = -t387 * t395 - t399 * t421;
	t413 = t373 * t400 + t383 * t394;
	t362 = t374 * t406 + t413 * t404;
	t428 = -pkin(10) - r_i_i_C(3);
	t418 = qJD(4) * t404;
	t417 = qJD(4) * t406;
	t416 = t394 * t417;
	t415 = t400 * t417;
	t403 = sin(qJ(5));
	t405 = cos(qJ(5));
	t414 = t403 * r_i_i_C(1) + t405 * r_i_i_C(2);
	t408 = qJD(5) * t414;
	t407 = qJD(4) * (t405 * r_i_i_C(1) - t403 * r_i_i_C(2) + pkin(4));
	t377 = -t381 * t394 + t386 * t400;
	t370 = -t375 * t394 + t384 * t400;
	t369 = -t373 * t394 + t383 * t400;
	t365 = -t381 * t415 + t382 * t418 - t386 * t416;
	t359 = -t375 * t415 + t376 * t418 - t384 * t416;
	t357 = -t373 * t415 + t374 * t418 - t383 * t416;
	t1 = [0, 0, 0, t428 * t359 - (-t376 * t404 + t412 * t406) * t408 - t364 * t407, t414 * t359 + ((-t364 * t405 - t370 * t403) * r_i_i_C(1) + (t364 * t403 - t370 * t405) * r_i_i_C(2)) * qJD(5), 0; 0, 0, 0, t428 * t357 - (-t374 * t404 + t413 * t406) * t408 - t362 * t407, t414 * t357 + ((-t362 * t405 - t369 * t403) * r_i_i_C(1) + (t362 * t403 - t369 * t405) * r_i_i_C(2)) * qJD(5), 0; 0, 0, 0, t428 * t365 - (-t382 * t404 + t411 * t406) * t408 - t368 * t407, t414 * t365 + ((-t368 * t405 - t377 * t403) * r_i_i_C(1) + (t368 * t403 - t377 * t405) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:28
	% EndTime: 2019-10-10 08:49:29
	% DurationCPUTime: 1.03s
	% Computational Cost: add. (1582->99), mult. (4747->188), div. (0->0), fcn. (6032->18), ass. (0->84)
	t586 = sin(pkin(13));
	t587 = sin(pkin(12));
	t572 = t587 * t586;
	t592 = cos(pkin(13));
	t593 = cos(pkin(12));
	t579 = t593 * t592;
	t596 = cos(pkin(6));
	t559 = t596 * t579 - t572;
	t595 = cos(pkin(7));
	t556 = t559 * t595;
	t573 = t587 * t592;
	t578 = t593 * t586;
	t560 = t596 * t578 + t573;
	t589 = sin(pkin(7));
	t591 = cos(pkin(14));
	t575 = t589 * t591;
	t590 = sin(pkin(6));
	t565 = t590 * t575;
	t585 = sin(pkin(14));
	t542 = -t591 * t556 + t560 * t585 + t593 * t565;
	t580 = t595 * t590;
	t551 = -t559 * t589 - t593 * t580;
	t588 = sin(pkin(8));
	t594 = cos(pkin(8));
	t604 = t542 * t594 - t551 * t588;
	t561 = -t596 * t573 - t578;
	t557 = t561 * t595;
	t562 = -t596 * t572 + t579;
	t543 = -t591 * t557 + t562 * t585 - t587 * t565;
	t552 = -t561 * t589 + t587 * t580;
	t603 = t543 * t594 - t552 * t588;
	t577 = t590 * t592;
	t566 = t595 * t577;
	t576 = t590 * t586;
	t550 = -t591 * t566 - t596 * t575 + t585 * t576;
	t558 = -t589 * t577 + t596 * t595;
	t602 = t550 * t594 - t558 * t588;
	t597 = cos(qJ(4));
	t601 = t604 * t597;
	t600 = t603 * t597;
	t599 = t602 * t597;
	t533 = sin(qJ(6));
	t536 = cos(qJ(6));
	t567 = qJD(6) * (t533 * r_i_i_C(1) + t536 * r_i_i_C(2));
	t598 = pkin(11) + r_i_i_C(3);
	t535 = sin(qJ(4));
	t584 = qJD(4) * t535;
	t583 = qJD(6) * t533;
	t582 = qJD(6) * t536;
	t574 = t589 * t585;
	t564 = t590 * t574;
	t521 = t585 * t556 + t560 * t591 - t593 * t564;
	t508 = t521 * t597 - t604 * t535;
	t515 = t542 * t588 + t551 * t594;
	t534 = sin(qJ(5));
	t537 = cos(qJ(5));
	t498 = t508 * t537 + t515 * t534;
	t571 = -t508 * t534 + t515 * t537;
	t522 = t585 * t557 + t562 * t591 + t587 * t564;
	t510 = t522 * t597 - t603 * t535;
	t516 = t543 * t588 + t552 * t594;
	t500 = t510 * t537 + t516 * t534;
	t570 = -t510 * t534 + t516 * t537;
	t530 = t585 * t566 + t596 * t574 + t591 * t576;
	t514 = t530 * t597 - t602 * t535;
	t523 = t550 * t588 + t558 * t594;
	t506 = t514 * t537 + t523 * t534;
	t569 = -t514 * t534 + t523 * t537;
	t568 = t536 * r_i_i_C(1) - t533 * r_i_i_C(2) + pkin(5);
	t563 = -t598 * t534 - t568 * t537 - pkin(4);
	t553 = t537 * t567 + (t568 * t534 - t598 * t537) * qJD(5);
	t513 = t530 * t535 + t599;
	t512 = t514 * qJD(4);
	t511 = t599 * qJD(4) + t530 * t584;
	t509 = t522 * t535 + t600;
	t507 = t521 * t535 + t601;
	t504 = t510 * qJD(4);
	t503 = t600 * qJD(4) + t522 * t584;
	t502 = t508 * qJD(4);
	t501 = t601 * qJD(4) + t521 * t584;
	t496 = t569 * qJD(5) - t511 * t537;
	t494 = t570 * qJD(5) - t503 * t537;
	t492 = t571 * qJD(5) - t501 * t537;
	t1 = [0, 0, 0, (-t503 * t533 + t510 * t582) * r_i_i_C(1) + (-t503 * t536 - t510 * t583) * r_i_i_C(2) - t503 * pkin(10) + t563 * t504 + t553 * t509, t598 * t494 - t570 * t567 + t568 * (-t500 * qJD(5) + t503 * t534), (-t494 * t533 + t504 * t536) * r_i_i_C(1) + (-t494 * t536 - t504 * t533) * r_i_i_C(2) + ((-t500 * t536 - t509 * t533) * r_i_i_C(1) + (t500 * t533 - t509 * t536) * r_i_i_C(2)) * qJD(6); 0, 0, 0, (-t501 * t533 + t508 * t582) * r_i_i_C(1) + (-t501 * t536 - t508 * t583) * r_i_i_C(2) - t501 * pkin(10) + t563 * t502 + t553 * t507, t598 * t492 - t571 * t567 + t568 * (-t498 * qJD(5) + t501 * t534), (-t492 * t533 + t502 * t536) * r_i_i_C(1) + (-t492 * t536 - t502 * t533) * r_i_i_C(2) + ((-t498 * t536 - t507 * t533) * r_i_i_C(1) + (t498 * t533 - t507 * t536) * r_i_i_C(2)) * qJD(6); 0, 0, 0, (-t511 * t533 + t514 * t582) * r_i_i_C(1) + (-t511 * t536 - t514 * t583) * r_i_i_C(2) - t511 * pkin(10) + t563 * t512 + t553 * t513, t598 * t496 - t569 * t567 + t568 * (-t506 * qJD(5) + t511 * t534), (-t496 * t533 + t512 * t536) * r_i_i_C(1) + (-t496 * t536 - t512 * t533) * r_i_i_C(2) + ((-t506 * t536 - t513 * t533) * r_i_i_C(1) + (t506 * t533 - t513 * t536) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end