% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRP5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:17:32
% EndTime: 2019-02-26 20:17:33
% DurationCPUTime: 1.28s
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
t1 = [0, -t584 * qJD(6) + t505 * pkin(3) + t540 * pkin(2) - t539 * t636 + t588 * (t584 * qJD(4) + t505 * t564 - t539 * t631) + t580 * (t525 * qJD(3) - t539 * t617 + t540 * t561) + t634 * (t511 * qJD(4) + t505 * t560 + t539 * t630) + ((-t511 * t563 - t524 * t559) * r_i_i_C(2) + t638 * (-t511 * t559 + t524 * t563)) * qJD(5) (t499 * t563 - t517 * t623) * r_i_i_C(2) + t499 * pkin(10) + t567 * t498 + t566 * t516 + t638 * (t499 * t559 + t517 * t622) qJD(6) * t509 - t588 * t488 + t634 * t489 - t593 * t575, t600 * r_i_i_C(1) + (-t489 * t563 - t498 * t559) * r_i_i_C(2) + (t597 * r_i_i_C(1) + (t509 * t559 - t516 * t563) * r_i_i_C(2)) * qJD(5) + (t597 * qJD(5) + t600) * pkin(5), t488; 0, -t585 * qJD(6) + t503 * pkin(3) + t538 * pkin(2) - t537 * t636 + t588 * (t585 * qJD(4) + t503 * t564 - t537 * t631) + t580 * (t523 * qJD(3) - t537 * t617 + t538 * t561) + t634 * (t510 * qJD(4) + t503 * t560 + t537 * t630) + ((-t510 * t563 - t522 * t559) * r_i_i_C(2) + t638 * (-t510 * t559 + t522 * t563)) * qJD(5) (t497 * t563 - t515 * t623) * r_i_i_C(2) + t497 * pkin(10) + t567 * t496 + t566 * t514 + t638 * (t497 * t559 + t515 * t622) qJD(6) * t507 - t588 * t486 + t634 * t487 - t594 * t575, t601 * r_i_i_C(1) + (-t487 * t563 - t496 * t559) * r_i_i_C(2) + (t598 * r_i_i_C(1) + (t507 * t559 - t514 * t563) * r_i_i_C(2)) * qJD(5) + (t598 * qJD(5) + t601) * pkin(5), t486; 0, -t579 * qJD(6) + t521 * pkin(3) + t588 * (t579 * qJD(4) + t521 * t564 + t560 * t607) - t580 * (-qJD(2) * t602 - t555 * qJD(3) * t615 + (qJD(3) * t557 + qJD(2)) * t620) + t634 * (t526 * qJD(4) + t521 * t560 - t564 * t607) + (-pkin(2) * t562 + pkin(9) * t629) * t624 + ((-t526 * t563 - t535 * t559) * r_i_i_C(2) + t638 * (-t526 * t559 + t535 * t563)) * qJD(5) (t513 * t563 - t528 * t623) * r_i_i_C(2) + t513 * pkin(10) + t567 * t512 + t566 * t527 + t638 * (t513 * t559 + t528 * t622) qJD(6) * t519 - t588 * t494 + t634 * t495 - t590 * t575, t599 * r_i_i_C(1) + (-t495 * t563 - t512 * t559) * r_i_i_C(2) + (t592 * r_i_i_C(1) + (t519 * t559 - t527 * t563) * r_i_i_C(2)) * qJD(5) + (t592 * qJD(5) + t599) * pkin(5), t494;];
JaD_transl  = t1;
