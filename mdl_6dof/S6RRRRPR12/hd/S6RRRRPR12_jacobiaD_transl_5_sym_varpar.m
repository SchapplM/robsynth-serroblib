% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR12_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR12_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:37:02
% EndTime: 2019-02-26 22:37:04
% DurationCPUTime: 1.25s
% Computational Cost: add. (1013->180), mult. (2885->302), div. (0->0), fcn. (3026->14), ass. (0->99)
t521 = sin(qJ(1));
t520 = sin(qJ(2));
t577 = cos(pkin(6));
t550 = t521 * t577;
t543 = t520 * t550;
t561 = qJD(2) * t520;
t524 = cos(qJ(2));
t525 = cos(qJ(1));
t563 = t525 * t524;
t489 = -qJD(1) * t543 - t521 * t561 + (qJD(2) * t577 + qJD(1)) * t563;
t514 = sin(pkin(7));
t515 = sin(pkin(6));
t570 = t515 * t525;
t554 = t514 * t570;
t588 = -qJD(3) * t554 + t489;
t549 = t525 * t577;
t498 = t520 * t549 + t521 * t524;
t529 = t525 * t520 + t524 * t550;
t488 = t529 * qJD(1) + t498 * qJD(2);
t516 = cos(pkin(7));
t519 = sin(qJ(3));
t523 = cos(qJ(3));
t562 = qJD(1) * t515;
t552 = t521 * t562;
t547 = t514 * t552;
t497 = t521 * t520 - t524 * t549;
t574 = t497 * t516;
t463 = (-qJD(3) * t498 - t488 * t516 + t547) * t519 + (-qJD(3) * t574 + t588) * t523;
t575 = t488 * t514;
t479 = t516 * t552 + t575;
t513 = qJ(4) + pkin(13);
t511 = sin(t513);
t512 = cos(t513);
t587 = t463 * t511 - t479 * t512;
t586 = -t463 * t512 - t479 * t511;
t569 = t516 * t519;
t538 = t497 * t569 - t498 * t523;
t474 = t519 * t554 + t538;
t492 = -t497 * t514 + t516 * t570;
t583 = t474 * t512 + t492 * t511;
t582 = -t474 * t511 + t492 * t512;
t581 = (t554 + t574) * t523 + t498 * t519;
t568 = t516 * t523;
t580 = t538 * qJD(3) - t488 * t568 - t588 * t519 + t523 * t547;
t518 = sin(qJ(4));
t579 = t518 * pkin(4);
t578 = r_i_i_C(3) + qJ(5) + pkin(11);
t530 = t543 - t563;
t486 = t497 * qJD(1) + t530 * qJD(2);
t576 = t486 * t514;
t572 = t514 * t515;
t571 = t515 * t521;
t567 = t519 * t520;
t566 = t519 * t524;
t565 = t520 * t523;
t564 = t523 * t524;
t560 = qJD(2) * t524;
t559 = qJD(4) * t511;
t558 = qJD(4) * t512;
t557 = qJD(4) * t520;
t522 = cos(qJ(4));
t556 = qJD(4) * t522;
t555 = qJD(4) * t579;
t553 = pkin(10) * t516 + pkin(9);
t551 = t525 * t562;
t548 = t577 * t514;
t546 = t514 * t551;
t545 = t561 * t572;
t542 = qJD(3) * t548;
t510 = t522 * pkin(4) + pkin(3);
t540 = -t512 * r_i_i_C(1) + t511 * r_i_i_C(2) - t510;
t539 = t497 * t519 - t498 * t568;
t484 = -t497 * t523 - t498 * t569;
t536 = t519 * t529 + t530 * t568;
t485 = -t523 * t529 + t530 * t569;
t535 = t514 * t571 - t516 * t529;
t534 = t516 * t564 - t567;
t533 = t516 * t565 + t566;
t532 = t516 * t566 + t565;
t531 = t516 * t567 - t564;
t527 = qJD(4) * (-t511 * r_i_i_C(1) - t512 * r_i_i_C(2) - t579);
t475 = t519 * t530 + t535 * t523;
t476 = t535 * t519 - t523 * t530;
t496 = t577 * t516 - t524 * t572;
t495 = t531 * t515;
t494 = t514 * t529 + t516 * t571;
t491 = t532 * t515 + t519 * t548;
t487 = t498 * qJD(1) + t529 * qJD(2);
t481 = (-t532 * qJD(2) - t533 * qJD(3)) * t515;
t477 = t516 * t551 - t576;
t471 = t523 * t542 + (-t531 * qJD(2) + t534 * qJD(3)) * t515;
t470 = t519 * t542 + (t533 * qJD(2) + t532 * qJD(3)) * t515;
t469 = t539 * qJD(3) - t488 * t523 - t489 * t569;
t467 = t536 * qJD(3) + t486 * t523 + t487 * t569;
t461 = -t487 * t523 + (t486 * t516 + t546) * t519 + t475 * qJD(3);
t460 = t476 * qJD(3) - t486 * t568 - t487 * t519 - t523 * t546;
t459 = t461 * t512 + t477 * t511 + (-t476 * t511 + t494 * t512) * qJD(4);
t458 = -t461 * t511 + t477 * t512 + (-t476 * t512 - t494 * t511) * qJD(4);
t1 = [t586 * r_i_i_C(1) + t587 * r_i_i_C(2) - t463 * t510 - t581 * qJD(5) - t479 * t579 - t489 * pkin(2) - pkin(10) * t575 + t578 * t580 + (-t525 * pkin(1) - t553 * t571) * qJD(1) + (t582 * r_i_i_C(1) - t583 * r_i_i_C(2) + (-t474 * t518 + t492 * t522) * pkin(4)) * qJD(4) (t467 * t512 - t485 * t559) * r_i_i_C(1) + (-t467 * t511 - t485 * t558) * r_i_i_C(2) + t467 * t510 - t485 * t555 - t536 * qJD(5) + t486 * pkin(2) + t578 * (t485 * qJD(3) + t486 * t519 - t487 * t568) + ((-t487 * t511 - t530 * t558) * r_i_i_C(1) + (-t487 * t512 + t530 * t559) * r_i_i_C(2) - t487 * pkin(10) + (-t487 * t518 - t530 * t556) * pkin(4)) * t514, t476 * qJD(5) + t540 * t460 + t578 * t461 + t475 * t527, t458 * r_i_i_C(1) - t459 * r_i_i_C(2) + (-t461 * t518 + t477 * t522 + (-t476 * t522 - t494 * t518) * qJD(4)) * pkin(4), t460, 0; -pkin(10) * t576 - t487 * pkin(2) + t459 * r_i_i_C(1) + t458 * r_i_i_C(2) - t475 * qJD(5) + t461 * t510 + t578 * t460 + (-pkin(1) * t521 + t553 * t570) * qJD(1) + (t477 * t518 + (-t476 * t518 + t494 * t522) * qJD(4)) * pkin(4) (t469 * t512 - t484 * t559) * r_i_i_C(1) + (-t469 * t511 - t484 * t558) * r_i_i_C(2) + t469 * t510 - t484 * t555 - t539 * qJD(5) - t488 * pkin(2) + t578 * (t484 * qJD(3) - t488 * t519 + t489 * t568) + ((t489 * t511 + t498 * t558) * r_i_i_C(1) + (t489 * t512 - t498 * t559) * r_i_i_C(2) + t489 * pkin(10) + (t489 * t518 + t498 * t556) * pkin(4)) * t514, -qJD(5) * t474 + t578 * t463 - t581 * t527 - t540 * t580, -t587 * r_i_i_C(1) + t586 * r_i_i_C(2) + (t583 * r_i_i_C(1) + t582 * r_i_i_C(2)) * qJD(4) + (-t463 * t518 + t479 * t522 + (t474 * t522 + t492 * t518) * qJD(4)) * pkin(4), -t580, 0; 0 (t481 * t512 + t495 * t559) * r_i_i_C(1) + (-t481 * t511 + t495 * t558) * r_i_i_C(2) + t481 * t510 + t495 * t555 + (-t578 * (-t534 * qJD(2) + t531 * qJD(3)) + t533 * qJD(5) - pkin(2) * t561 + ((t511 * t560 + t512 * t557) * r_i_i_C(1) + (-t511 * t557 + t512 * t560) * r_i_i_C(2) + pkin(10) * t560 + (t518 * t560 + t520 * t556) * pkin(4)) * t514) * t515, t491 * qJD(5) + t578 * t471 + t540 * t470 + (t534 * t515 + t523 * t548) * t527 (-t471 * t511 + t512 * t545) * r_i_i_C(1) + (-t471 * t512 - t511 * t545) * r_i_i_C(2) + ((-t491 * t512 - t496 * t511) * r_i_i_C(1) + (t491 * t511 - t496 * t512) * r_i_i_C(2)) * qJD(4) + (t522 * t545 - t471 * t518 + (-t491 * t522 - t496 * t518) * qJD(4)) * pkin(4), t470, 0;];
JaD_transl  = t1;
