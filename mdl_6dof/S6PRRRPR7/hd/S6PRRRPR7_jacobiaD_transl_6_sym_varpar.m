% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR7
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
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:07
% EndTime: 2019-02-26 20:14:08
% DurationCPUTime: 0.92s
% Computational Cost: add. (1403->169), mult. (4167->294), div. (0->0), fcn. (4622->16), ass. (0->106)
t568 = sin(qJ(2));
t570 = cos(qJ(2));
t563 = cos(pkin(12));
t625 = cos(pkin(6));
t605 = t563 * t625;
t624 = sin(pkin(12));
t543 = -t624 * t568 + t570 * t605;
t559 = pkin(13) + qJ(6);
t557 = sin(t559);
t558 = cos(t559);
t596 = t557 * r_i_i_C(1) + t558 * r_i_i_C(2);
t585 = qJD(6) * t596;
t628 = cos(qJ(3));
t561 = sin(pkin(7));
t627 = pkin(9) * t561;
t626 = r_i_i_C(3) + pkin(11) + qJ(5);
t566 = sin(qJ(4));
t623 = t561 * t566;
t569 = cos(qJ(4));
t622 = t561 * t569;
t621 = t561 * t570;
t562 = sin(pkin(6));
t620 = t562 * t563;
t564 = cos(pkin(7));
t567 = sin(qJ(3));
t619 = t564 * t567;
t618 = t567 * t568;
t617 = t567 * t570;
t616 = qJD(2) * t562;
t615 = t561 * t620;
t614 = t561 * t562 * t568;
t613 = t562 * t618;
t577 = -t568 * t605 - t624 * t570;
t612 = t577 * t628;
t611 = t564 * t628;
t610 = t628 * t568;
t609 = t628 * t570;
t608 = t561 * t616;
t607 = t561 * t625;
t606 = t562 * t624;
t603 = t564 * t609;
t602 = t568 * t608;
t601 = t570 * t608;
t599 = t567 * t607;
t598 = t561 * t606;
t597 = t558 * r_i_i_C(1) - t557 * r_i_i_C(2);
t595 = t625 * t624;
t594 = t562 * t603;
t516 = -t612 + (t543 * t564 - t615) * t567;
t530 = -t543 * t561 - t564 * t620;
t508 = t516 * t569 + t530 * t566;
t593 = -t516 * t566 + t530 * t569;
t575 = t563 * t568 + t570 * t595;
t576 = -t563 * t570 + t568 * t595;
t518 = -t576 * t628 + (-t564 * t575 + t598) * t567;
t531 = t561 * t575 + t564 * t606;
t510 = t518 * t569 + t531 * t566;
t592 = -t518 * t566 + t531 * t569;
t580 = t564 * t617 + t610;
t529 = t580 * t562 + t599;
t542 = -t562 * t621 + t625 * t564;
t520 = t529 * t569 + t542 * t566;
t591 = -t529 * t566 + t542 * t569;
t590 = t628 * t607;
t589 = cos(pkin(13)) * pkin(5) + pkin(4) + t597;
t524 = t543 * t628 + t577 * t619;
t588 = -t524 * t566 - t577 * t622;
t511 = t524 * t569 - t577 * t623;
t526 = -t575 * t628 + t576 * t619;
t587 = -t526 * t566 - t576 * t622;
t512 = t526 * t569 - t576 * t623;
t586 = qJD(6) * t597;
t581 = -t564 * t618 + t609;
t537 = t581 * t562;
t584 = -t537 * t566 + t569 * t614;
t527 = t537 * t569 + t566 * t614;
t583 = -t543 * t567 + t577 * t611;
t582 = t567 * t575 + t576 * t611;
t579 = t564 * t610 + t617;
t578 = sin(pkin(13)) * pkin(5) + pkin(10) + t596;
t574 = -t626 * t566 - t589 * t569 - pkin(3);
t573 = t543 * t611 + t567 * t577 - t628 * t615;
t572 = t567 * t576 - t575 * t611 + t628 * t598;
t571 = -t566 * qJD(5) + t569 * t585 + (t589 * t566 - t626 * t569) * qJD(4);
t541 = t576 * qJD(2);
t540 = t575 * qJD(2);
t539 = t577 * qJD(2);
t538 = t543 * qJD(2);
t536 = t579 * t562;
t528 = -t590 - t594 + t613;
t522 = (-t580 * qJD(2) - t579 * qJD(3)) * t562;
t514 = qJD(3) * t590 + ((t603 - t618) * qJD(3) + t581 * qJD(2)) * t562;
t513 = qJD(3) * t599 + (t579 * qJD(2) + t580 * qJD(3)) * t562;
t506 = t582 * qJD(3) + t540 * t619 + t541 * t628;
t504 = t583 * qJD(3) - t538 * t619 + t539 * t628;
t500 = t572 * qJD(3) - t540 * t628 + t541 * t619;
t499 = t518 * qJD(3) - t540 * t567 - t541 * t611;
t498 = t573 * qJD(3) + t538 * t628 + t539 * t619;
t497 = t538 * t567 - t539 * t611 + (t543 * t619 - t567 * t615 - t612) * qJD(3);
t496 = t591 * qJD(4) + t514 * t569 + t566 * t602;
t495 = t520 * qJD(4) + t514 * t566 - t569 * t602;
t490 = t592 * qJD(4) + t500 * t569 - t541 * t623;
t489 = t510 * qJD(4) + t500 * t566 + t541 * t622;
t488 = t593 * qJD(4) + t498 * t569 - t539 * t623;
t487 = t508 * qJD(4) + t498 * t566 + t539 * t622;
t1 = [0, -t587 * qJD(5) + t506 * pkin(3) + t541 * pkin(2) - t540 * t627 + t589 * (t587 * qJD(4) + t506 * t569 - t540 * t623) + t578 * (t526 * qJD(3) - t540 * t611 + t541 * t567) + t626 * (t512 * qJD(4) + t506 * t566 + t540 * t622) + ((-t512 * t557 - t558 * t582) * r_i_i_C(1) + (-t512 * t558 + t557 * t582) * r_i_i_C(2)) * qJD(6), t574 * t499 + t578 * t500 + t518 * t586 - t571 * t572, t510 * qJD(5) - t589 * t489 + t626 * t490 - t592 * t585, t489 (-t490 * t557 + t499 * t558) * r_i_i_C(1) + (-t490 * t558 - t499 * t557) * r_i_i_C(2) + ((-t510 * t558 + t557 * t572) * r_i_i_C(1) + (t510 * t557 + t558 * t572) * r_i_i_C(2)) * qJD(6); 0, -t588 * qJD(5) + t504 * pkin(3) + t539 * pkin(2) + t538 * t627 + t589 * (t588 * qJD(4) + t504 * t569 + t538 * t623) + t578 * (t524 * qJD(3) + t538 * t611 + t539 * t567) + t626 * (t511 * qJD(4) + t504 * t566 - t538 * t622) + ((-t511 * t557 - t558 * t583) * r_i_i_C(1) + (-t511 * t558 + t557 * t583) * r_i_i_C(2)) * qJD(6), t574 * t497 + t578 * t498 + t516 * t586 - t571 * t573, t508 * qJD(5) - t589 * t487 + t626 * t488 - t593 * t585, t487 (-t488 * t557 + t497 * t558) * r_i_i_C(1) + (-t488 * t558 - t497 * t557) * r_i_i_C(2) + ((-t508 * t558 + t557 * t573) * r_i_i_C(1) + (t508 * t557 + t558 * t573) * r_i_i_C(2)) * qJD(6); 0, -t584 * qJD(5) + t522 * pkin(3) + t589 * (t584 * qJD(4) + t522 * t569 + t566 * t601) - t578 * (-qJD(2) * t594 - t562 * qJD(3) * t609 + (qJD(3) * t564 + qJD(2)) * t613) + t626 * (t527 * qJD(4) + t522 * t566 - t569 * t601) + ((-t527 * t557 + t536 * t558) * r_i_i_C(1) + (-t527 * t558 - t536 * t557) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t568 + pkin(9) * t621) * t616, t574 * t513 + t578 * t514 + t571 * t528 + t529 * t586, t520 * qJD(5) - t589 * t495 + t626 * t496 - t591 * t585, t495 (-t496 * t557 + t513 * t558) * r_i_i_C(1) + (-t496 * t558 - t513 * t557) * r_i_i_C(2) + ((-t520 * t558 - t528 * t557) * r_i_i_C(1) + (t520 * t557 - t528 * t558) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
