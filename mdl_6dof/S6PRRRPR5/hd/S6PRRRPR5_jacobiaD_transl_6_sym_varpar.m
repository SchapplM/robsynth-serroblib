% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:58
% EndTime: 2019-02-26 20:12:59
% DurationCPUTime: 1.15s
% Computational Cost: add. (1504->210), mult. (4137->363), div. (0->0), fcn. (4566->16), ass. (0->112)
t569 = sin(qJ(2));
t572 = cos(qJ(2));
t563 = cos(pkin(12));
t625 = cos(pkin(6));
t602 = t563 * t625;
t624 = sin(pkin(12));
t545 = -t624 * t569 + t572 * t602;
t566 = sin(qJ(6));
t570 = cos(qJ(6));
t585 = (t566 * r_i_i_C(1) + t570 * r_i_i_C(2)) * qJD(6);
t630 = r_i_i_C(3) + pkin(11);
t629 = cos(qJ(3));
t561 = sin(pkin(7));
t628 = t561 * pkin(9);
t567 = sin(qJ(4));
t627 = t567 * pkin(4);
t626 = pkin(4) * qJD(4);
t560 = qJ(4) + pkin(13);
t558 = sin(t560);
t623 = t558 * t561;
t559 = cos(t560);
t622 = t559 * t561;
t562 = sin(pkin(6));
t621 = t561 * t562;
t620 = t561 * t567;
t571 = cos(qJ(4));
t619 = t561 * t571;
t618 = t561 * t572;
t617 = t562 * t563;
t564 = cos(pkin(7));
t568 = sin(qJ(3));
t616 = t564 * t568;
t615 = t568 * t569;
t614 = t568 * t572;
t613 = qJD(6) * t566;
t612 = qJD(6) * t570;
t611 = t561 * t617;
t610 = t569 * t621;
t609 = t562 * t615;
t579 = -t569 * t602 - t572 * t624;
t608 = t579 * t629;
t607 = t564 * t629;
t606 = t629 * t569;
t605 = t629 * t572;
t604 = qJD(2) * t621;
t603 = t562 * t624;
t601 = t625 * t561;
t599 = t564 * t605;
t598 = t569 * t604;
t597 = t572 * t604;
t595 = t568 * t601;
t594 = t561 * t603;
t592 = t625 * t624;
t591 = t562 * t599;
t520 = -t608 + (t545 * t564 - t611) * t568;
t532 = -t545 * t561 - t564 * t617;
t506 = t520 * t559 + t532 * t558;
t590 = -t520 * t558 + t532 * t559;
t577 = t563 * t569 + t572 * t592;
t578 = -t563 * t572 + t569 * t592;
t522 = -t578 * t629 + (-t564 * t577 + t594) * t568;
t533 = t561 * t577 + t564 * t603;
t508 = t522 * t559 + t533 * t558;
t589 = -t522 * t558 + t533 * t559;
t581 = t564 * t614 + t606;
t531 = t562 * t581 + t595;
t544 = -t562 * t618 + t564 * t625;
t516 = t531 * t559 + t544 * t558;
t588 = -t531 * t558 + t544 * t559;
t587 = r_i_i_C(1) * t570 - r_i_i_C(2) * t566 + pkin(5);
t586 = t629 * t601;
t526 = t545 * t629 + t579 * t616;
t513 = t526 * t559 - t579 * t623;
t528 = -t577 * t629 + t578 * t616;
t514 = t528 * t559 - t578 * t623;
t582 = -t564 * t615 + t605;
t539 = t582 * t562;
t529 = t539 * t559 + t558 * t610;
t584 = -t545 * t568 + t579 * t607;
t583 = t568 * t577 + t578 * t607;
t580 = t564 * t606 + t614;
t557 = pkin(4) * t571 + pkin(3);
t576 = -t558 * t630 - t587 * t559 - t557;
t575 = t545 * t607 + t568 * t579 - t611 * t629;
t574 = t568 * t578 - t577 * t607 + t594 * t629;
t573 = t559 * t585 + (t587 * t558 - t559 * t630 + t627) * qJD(4);
t565 = -qJ(5) - pkin(10);
t543 = t578 * qJD(2);
t542 = t577 * qJD(2);
t541 = t579 * qJD(2);
t540 = t545 * qJD(2);
t538 = t580 * t562;
t530 = -t586 - t591 + t609;
t524 = (-qJD(2) * t581 - qJD(3) * t580) * t562;
t523 = -qJD(2) * t591 - t562 * qJD(3) * t605 + (qJD(3) * t564 + qJD(2)) * t609;
t518 = qJD(3) * t586 + ((t599 - t615) * qJD(3) + t582 * qJD(2)) * t562;
t517 = qJD(3) * t595 + (qJD(2) * t580 + qJD(3) * t581) * t562;
t512 = qJD(3) * t583 + t542 * t616 + t543 * t629;
t511 = qJD(3) * t528 - t542 * t607 + t543 * t568;
t510 = qJD(3) * t584 - t540 * t616 + t541 * t629;
t509 = qJD(3) * t526 + t540 * t607 + t541 * t568;
t504 = qJD(3) * t574 - t542 * t629 + t543 * t616;
t503 = qJD(3) * t522 - t542 * t568 - t543 * t607;
t502 = qJD(3) * t575 + t540 * t629 + t541 * t616;
t501 = t540 * t568 - t541 * t607 + (t545 * t616 - t568 * t611 - t608) * qJD(3);
t500 = t558 * t597 + t524 * t559 + (-t539 * t558 + t559 * t610) * qJD(4);
t498 = qJD(4) * t588 + t518 * t559 + t558 * t598;
t496 = -t542 * t623 + t512 * t559 + (-t528 * t558 - t578 * t622) * qJD(4);
t494 = t540 * t623 + t510 * t559 + (-t526 * t558 - t579 * t622) * qJD(4);
t492 = qJD(4) * t589 + t504 * t559 - t543 * t623;
t490 = qJD(4) * t590 + t502 * t559 - t541 * t623;
t1 = [0 (t496 * t570 + t511 * t566) * r_i_i_C(1) + (-t496 * t566 + t511 * t570) * r_i_i_C(2) + t496 * pkin(5) + t512 * t557 - t511 * t565 - t583 * qJD(5) + t543 * pkin(2) - t542 * t628 + t630 * (qJD(4) * t514 + t512 * t558 + t542 * t622) + ((-t514 * t566 - t570 * t583) * r_i_i_C(1) + (-t514 * t570 + t566 * t583) * r_i_i_C(2)) * qJD(6) + (-t542 * t620 + (-t528 * t567 - t578 * t619) * qJD(4)) * pkin(4) (t504 * t566 + t522 * t612) * r_i_i_C(1) + (t504 * t570 - t522 * t613) * r_i_i_C(2) - t504 * t565 + t522 * qJD(5) + t576 * t503 - t573 * t574, t630 * t492 - t589 * t585 + t587 * (-qJD(4) * t508 - t504 * t558 - t543 * t622) + (-t543 * t619 - t504 * t567 + (-t522 * t571 - t533 * t567) * qJD(4)) * pkin(4), t503 (-t492 * t566 + t503 * t570) * r_i_i_C(1) + (-t492 * t570 - t503 * t566) * r_i_i_C(2) + ((-t508 * t570 + t566 * t574) * r_i_i_C(1) + (t508 * t566 + t570 * t574) * r_i_i_C(2)) * qJD(6); 0 (t494 * t570 + t509 * t566) * r_i_i_C(1) + (-t494 * t566 + t509 * t570) * r_i_i_C(2) + t494 * pkin(5) + t510 * t557 - t509 * t565 - t584 * qJD(5) + t541 * pkin(2) + t540 * t628 + t630 * (qJD(4) * t513 + t510 * t558 - t540 * t622) + ((-t513 * t566 - t570 * t584) * r_i_i_C(1) + (-t513 * t570 + t566 * t584) * r_i_i_C(2)) * qJD(6) + (t540 * t620 + (-t526 * t567 - t579 * t619) * qJD(4)) * pkin(4) (t502 * t566 + t520 * t612) * r_i_i_C(1) + (t502 * t570 - t520 * t613) * r_i_i_C(2) - t502 * t565 + t520 * qJD(5) + t576 * t501 - t573 * t575, t630 * t490 - t590 * t585 + t587 * (-qJD(4) * t506 - t502 * t558 - t541 * t622) + (-t541 * t619 - t502 * t567 + (-t520 * t571 - t532 * t567) * qJD(4)) * pkin(4), t501 (-t490 * t566 + t501 * t570) * r_i_i_C(1) + (-t490 * t570 - t501 * t566) * r_i_i_C(2) + ((-t506 * t570 + t566 * t575) * r_i_i_C(1) + (t506 * t566 + t570 * t575) * r_i_i_C(2)) * qJD(6); 0 (t500 * t570 - t523 * t566) * r_i_i_C(1) + (-t500 * t566 - t523 * t570) * r_i_i_C(2) + t500 * pkin(5) + t524 * t557 - t539 * t567 * t626 + t523 * t565 + t538 * qJD(5) + t630 * (qJD(4) * t529 + t524 * t558 - t559 * t597) + ((-t529 * t566 + t538 * t570) * r_i_i_C(1) + (-t529 * t570 - t538 * t566) * r_i_i_C(2)) * qJD(6) + (t569 * t619 * t626 + (-pkin(2) * t569 + (pkin(9) + t627) * t618) * qJD(2)) * t562 (t518 * t566 + t531 * t612) * r_i_i_C(1) + (t518 * t570 - t531 * t613) * r_i_i_C(2) - t518 * t565 + t531 * qJD(5) + t576 * t517 + t573 * t530, t630 * t498 - t588 * t585 + t587 * (-qJD(4) * t516 - t518 * t558 + t559 * t598) + (t571 * t598 - t518 * t567 + (-t531 * t571 - t544 * t567) * qJD(4)) * pkin(4), t517 (-t498 * t566 + t517 * t570) * r_i_i_C(1) + (-t498 * t570 - t517 * t566) * r_i_i_C(2) + ((-t516 * t570 - t530 * t566) * r_i_i_C(1) + (t516 * t566 - t530 * t570) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
