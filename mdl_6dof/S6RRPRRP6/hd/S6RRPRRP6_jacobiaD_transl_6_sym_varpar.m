% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:54
% EndTime: 2019-02-26 21:48:55
% DurationCPUTime: 1.06s
% Computational Cost: add. (1486->142), mult. (4318->232), div. (0->0), fcn. (4849->12), ass. (0->90)
t554 = cos(pkin(6));
t552 = sin(pkin(11));
t618 = cos(pkin(11));
t622 = cos(qJ(2));
t590 = t622 * t618;
t557 = sin(qJ(2));
t605 = qJD(2) * t557;
t627 = -qJD(2) * t590 + t552 * t605;
t532 = t627 * t554;
t593 = t557 * t618;
t542 = -t622 * t552 - t593;
t537 = t542 * t554;
t594 = qJD(2) * t622;
t539 = -qJD(2) * t593 - t552 * t594;
t558 = sin(qJ(1));
t561 = cos(qJ(1));
t572 = -t557 * t552 + t590;
t607 = qJD(1) * t558;
t507 = -t537 * t607 - t558 * t539 + (-qJD(1) * t572 + t532) * t561;
t556 = sin(qJ(4));
t560 = cos(qJ(4));
t580 = -t537 * t561 + t558 * t572;
t553 = sin(pkin(6));
t597 = t553 * t607;
t612 = t553 * t561;
t601 = t560 * t612;
t495 = (-qJD(4) * t580 + t597) * t556 - qJD(4) * t601 - t507 * t560;
t533 = qJD(2) * t537;
t606 = qJD(1) * t561;
t567 = t554 * t572;
t628 = qJD(1) * t567 + qJD(2) * t572;
t508 = t561 * t533 + t542 * t606 - t628 * t558;
t555 = sin(qJ(5));
t559 = cos(qJ(5));
t514 = t556 * t612 - t560 * t580;
t518 = t558 * t542 + t561 * t567;
t629 = t514 * t559 + t518 * t555;
t634 = t629 * qJD(5) - t495 * t555 - t508 * t559;
t630 = t514 * t555 - t518 * t559;
t633 = t630 * qJD(5) + t495 * t559 - t508 * t555;
t623 = r_i_i_C(2) + pkin(10);
t573 = qJD(4) * (-pkin(4) * t556 + t623 * t560);
t626 = t560 * pkin(4) + t623 * t556 + pkin(3);
t619 = r_i_i_C(3) + qJ(6);
t624 = r_i_i_C(1) + pkin(5);
t571 = t619 * t555 + t624 * t559 + pkin(4);
t621 = pkin(2) * t554;
t613 = t553 * t558;
t611 = t555 * t560;
t609 = t557 * t558;
t608 = t557 * t561;
t604 = qJD(4) * t556;
t603 = qJD(5) * t560;
t602 = pkin(2) * t605;
t599 = t622 * t558;
t598 = t622 * t561;
t596 = t553 * t606;
t521 = t542 * t561 - t558 * t567;
t563 = -t580 * qJD(1) + t558 * t532 + t561 * t539;
t589 = t521 * t603 - t563;
t588 = t518 * t603 + t507;
t530 = t627 * t553;
t535 = t572 * t553;
t587 = -t535 * t603 - t530;
t579 = t558 * t537 + t561 * t572;
t516 = t556 * t613 + t560 * t579;
t584 = t516 * t559 - t521 * t555;
t583 = -t516 * t555 - t521 * t559;
t536 = t542 * t553;
t524 = -t536 * t560 + t554 * t556;
t582 = t524 * t559 - t535 * t555;
t581 = t536 * t556 + t554 * t560;
t577 = -t556 * t579 + t560 * t613;
t505 = t558 * t533 + t542 * t607 + t628 * t561;
t570 = qJD(5) * t579 - t505 * t560 - t521 * t604;
t569 = qJD(5) * t580 + t508 * t560 - t518 * t604;
t531 = qJD(2) * t536;
t568 = -qJD(5) * t536 + t531 * t560 - t535 * t604;
t564 = qJD(6) * t555 + (-t624 * t555 + t619 * t559) * qJD(5);
t562 = t514 * qJD(4) + t507 * t556 + t560 * t597;
t551 = t622 * pkin(2) + pkin(1);
t540 = -t553 * qJD(3) + t594 * t621;
t538 = t557 * t621 + (-pkin(8) - qJ(3)) * t553;
t511 = t581 * qJD(4) - t530 * t560;
t498 = t582 * qJD(5) + t511 * t555 + t531 * t559;
t493 = t577 * qJD(4) + t556 * t596 + t560 * t563;
t492 = t516 * qJD(4) + t556 * t563 - t560 * t596;
t483 = t583 * qJD(5) + t493 * t559 + t505 * t555;
t482 = t584 * qJD(5) + t493 * t555 - t505 * t559;
t1 = [t630 * qJD(6) - t495 * pkin(4) + t507 * pkin(3) + t508 * pkin(9) + t558 * t602 - t561 * t540 + t623 * t562 - t624 * t633 + t619 * t634 + (t538 * t558 - t551 * t561) * qJD(1) -(-t521 * t611 + t559 * t579) * qJD(6) + t563 * pkin(9) + t624 * (-t589 * t555 + t570 * t559) + t619 * (t570 * t555 + t589 * t559) + t521 * t573 - t626 * t505 + ((t554 * t609 - t598) * qJD(2) + (-t554 * t598 + t609) * qJD(1)) * pkin(2), t596, -t571 * t492 + t623 * t493 + t564 * t577, t584 * qJD(6) - t624 * t482 + t619 * t483, t482; -t583 * qJD(6) + t493 * pkin(4) + t563 * pkin(3) + t505 * pkin(9) - t561 * t602 - t558 * t540 + t623 * t492 + t624 * t483 + t619 * t482 + (-t538 * t561 - t551 * t558) * qJD(1) -(-t518 * t611 + t559 * t580) * qJD(6) - t507 * pkin(9) + t624 * (-t588 * t555 + t569 * t559) + t619 * (t569 * t555 + t588 * t559) + t518 * t573 + t626 * t508 + ((-t554 * t608 - t599) * qJD(2) + (-t554 * t599 - t608) * qJD(1)) * pkin(2), t597, t623 * t495 + t564 * (-t556 * t580 - t601) + t571 * t562, -t629 * qJD(6) + t619 * t633 + t624 * t634, -t634; 0 -(-t535 * t611 - t536 * t559) * qJD(6) - t530 * pkin(9) - t553 * t602 + t624 * (t587 * t555 + t568 * t559) + t619 * (t568 * t555 - t587 * t559) + t535 * t573 + t626 * t531, 0, t623 * t511 + t564 * t581 + t571 * (-t524 * qJD(4) + t530 * t556) t582 * qJD(6) + t619 * (t511 * t559 - t531 * t555 + (-t524 * t555 - t535 * t559) * qJD(5)) - t624 * t498, t498;];
JaD_transl  = t1;
