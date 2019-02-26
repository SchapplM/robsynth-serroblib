% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:10
% EndTime: 2019-02-26 20:07:11
% DurationCPUTime: 0.99s
% Computational Cost: add. (1427->189), mult. (3882->327), div. (0->0), fcn. (4294->16), ass. (0->107)
t568 = sin(qJ(2));
t570 = cos(qJ(2));
t563 = cos(pkin(12));
t623 = cos(pkin(6));
t601 = t563 * t623;
t622 = sin(pkin(12));
t544 = -t622 * t568 + t570 * t601;
t566 = sin(qJ(6));
t569 = cos(qJ(6));
t583 = (t566 * r_i_i_C(1) + t569 * r_i_i_C(2)) * qJD(6);
t625 = r_i_i_C(3) + pkin(11);
t624 = cos(qJ(3));
t559 = pkin(13) + qJ(5);
t557 = sin(t559);
t561 = sin(pkin(7));
t621 = t557 * t561;
t558 = cos(t559);
t620 = t558 * t561;
t619 = t561 * t570;
t562 = sin(pkin(6));
t618 = t562 * t563;
t564 = cos(pkin(7));
t567 = sin(qJ(3));
t617 = t564 * t567;
t616 = t567 * t568;
t615 = t567 * t570;
t614 = qJD(2) * t562;
t613 = qJD(6) * t566;
t612 = qJD(6) * t569;
t611 = t561 * t618;
t610 = t561 * t562 * t568;
t609 = t562 * t616;
t608 = pkin(4) * sin(pkin(13)) + pkin(9);
t577 = -t568 * t601 - t622 * t570;
t607 = t577 * t624;
t606 = t564 * t624;
t605 = t624 * t568;
t604 = t624 * t570;
t603 = t561 * t614;
t602 = t562 * t622;
t600 = t623 * t561;
t598 = t564 * t604;
t597 = t568 * t603;
t596 = t570 * t603;
t594 = t567 * t600;
t593 = t561 * t602;
t592 = t608 * t561;
t590 = t623 * t622;
t589 = t562 * t598;
t519 = -t607 + (t544 * t564 - t611) * t567;
t531 = -t544 * t561 - t564 * t618;
t505 = t519 * t558 + t531 * t557;
t588 = -t519 * t557 + t531 * t558;
t575 = t563 * t568 + t570 * t590;
t576 = -t563 * t570 + t568 * t590;
t521 = -t576 * t624 + (-t564 * t575 + t593) * t567;
t532 = t561 * t575 + t564 * t602;
t507 = t521 * t558 + t532 * t557;
t587 = -t521 * t557 + t532 * t558;
t579 = t564 * t615 + t605;
t530 = t579 * t562 + t594;
t543 = -t562 * t619 + t623 * t564;
t515 = t530 * t558 + t543 * t557;
t586 = -t530 * t557 + t543 * t558;
t585 = t569 * r_i_i_C(1) - t566 * r_i_i_C(2) + pkin(5);
t584 = t624 * t600;
t525 = t544 * t624 + t577 * t617;
t512 = t525 * t558 - t577 * t621;
t527 = -t575 * t624 + t576 * t617;
t513 = t527 * t558 - t576 * t621;
t580 = -t564 * t616 + t604;
t538 = t580 * t562;
t528 = t538 * t558 + t557 * t610;
t582 = -t544 * t567 + t577 * t606;
t581 = t567 * t575 + t576 * t606;
t578 = t564 * t605 + t615;
t556 = cos(pkin(13)) * pkin(4) + pkin(3);
t574 = -t625 * t557 - t585 * t558 - t556;
t573 = t544 * t606 + t567 * t577 - t624 * t611;
t572 = t567 * t576 - t575 * t606 + t624 * t593;
t571 = t558 * t583 + (t585 * t557 - t625 * t558) * qJD(5);
t565 = -pkin(10) - qJ(4);
t542 = t576 * qJD(2);
t541 = t575 * qJD(2);
t540 = t577 * qJD(2);
t539 = t544 * qJD(2);
t537 = t578 * t562;
t529 = -t584 - t589 + t609;
t523 = (-t579 * qJD(2) - t578 * qJD(3)) * t562;
t522 = -qJD(2) * t589 - t562 * qJD(3) * t604 + (qJD(3) * t564 + qJD(2)) * t609;
t517 = qJD(3) * t584 + ((t598 - t616) * qJD(3) + t580 * qJD(2)) * t562;
t516 = qJD(3) * t594 + (t578 * qJD(2) + t579 * qJD(3)) * t562;
t511 = t581 * qJD(3) + t541 * t617 + t542 * t624;
t510 = t527 * qJD(3) - t541 * t606 + t542 * t567;
t509 = t582 * qJD(3) - t539 * t617 + t540 * t624;
t508 = t525 * qJD(3) + t539 * t606 + t540 * t567;
t503 = t572 * qJD(3) - t541 * t624 + t542 * t617;
t502 = t521 * qJD(3) - t541 * t567 - t542 * t606;
t501 = t573 * qJD(3) + t539 * t624 + t540 * t617;
t500 = t539 * t567 - t540 * t606 + (t544 * t617 - t567 * t611 - t607) * qJD(3);
t499 = t557 * t596 + t523 * t558 + (-t538 * t557 + t558 * t610) * qJD(5);
t497 = t586 * qJD(5) + t517 * t558 + t557 * t597;
t495 = -t541 * t621 + t511 * t558 + (-t527 * t557 - t576 * t620) * qJD(5);
t493 = t539 * t621 + t509 * t558 + (-t525 * t557 - t577 * t620) * qJD(5);
t491 = t587 * qJD(5) + t503 * t558 - t542 * t621;
t489 = t588 * qJD(5) + t501 * t558 - t540 * t621;
t1 = [0 (t495 * t569 + t510 * t566) * r_i_i_C(1) + (-t495 * t566 + t510 * t569) * r_i_i_C(2) + t495 * pkin(5) + t511 * t556 - t510 * t565 - t581 * qJD(4) + t542 * pkin(2) - t541 * t592 + t625 * (t513 * qJD(5) + t511 * t557 + t541 * t620) + ((-t513 * t566 - t569 * t581) * r_i_i_C(1) + (-t513 * t569 + t566 * t581) * r_i_i_C(2)) * qJD(6) (t503 * t566 + t521 * t612) * r_i_i_C(1) + (t503 * t569 - t521 * t613) * r_i_i_C(2) - t503 * t565 + t521 * qJD(4) + t574 * t502 - t571 * t572, t502, t625 * t491 - t587 * t583 + t585 * (-t507 * qJD(5) - t503 * t557 - t542 * t620) (-t491 * t566 + t502 * t569) * r_i_i_C(1) + (-t491 * t569 - t502 * t566) * r_i_i_C(2) + ((-t507 * t569 + t566 * t572) * r_i_i_C(1) + (t507 * t566 + t569 * t572) * r_i_i_C(2)) * qJD(6); 0 (t493 * t569 + t508 * t566) * r_i_i_C(1) + (-t493 * t566 + t508 * t569) * r_i_i_C(2) + t493 * pkin(5) + t509 * t556 - t508 * t565 - t582 * qJD(4) + t540 * pkin(2) + t539 * t592 + t625 * (t512 * qJD(5) + t509 * t557 - t539 * t620) + ((-t512 * t566 - t569 * t582) * r_i_i_C(1) + (-t512 * t569 + t566 * t582) * r_i_i_C(2)) * qJD(6) (t501 * t566 + t519 * t612) * r_i_i_C(1) + (t501 * t569 - t519 * t613) * r_i_i_C(2) - t501 * t565 + t519 * qJD(4) + t574 * t500 - t571 * t573, t500, t625 * t489 - t588 * t583 + t585 * (-t505 * qJD(5) - t501 * t557 - t540 * t620) (-t489 * t566 + t500 * t569) * r_i_i_C(1) + (-t489 * t569 - t500 * t566) * r_i_i_C(2) + ((-t505 * t569 + t566 * t573) * r_i_i_C(1) + (t505 * t566 + t569 * t573) * r_i_i_C(2)) * qJD(6); 0 (t499 * t569 - t522 * t566) * r_i_i_C(1) + (-t499 * t566 - t522 * t569) * r_i_i_C(2) + t499 * pkin(5) + t523 * t556 + t522 * t565 + t537 * qJD(4) + t625 * (t528 * qJD(5) + t523 * t557 - t558 * t596) + ((-t528 * t566 + t537 * t569) * r_i_i_C(1) + (-t528 * t569 - t537 * t566) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t568 + t608 * t619) * t614 (t517 * t566 + t530 * t612) * r_i_i_C(1) + (t517 * t569 - t530 * t613) * r_i_i_C(2) - t517 * t565 + t530 * qJD(4) + t574 * t516 + t571 * t529, t516, t625 * t497 - t586 * t583 + t585 * (-t515 * qJD(5) - t517 * t557 + t558 * t597) (-t497 * t566 + t516 * t569) * r_i_i_C(1) + (-t497 * t569 - t516 * t566) * r_i_i_C(2) + ((-t515 * t569 - t529 * t566) * r_i_i_C(1) + (t515 * t566 - t529 * t569) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
