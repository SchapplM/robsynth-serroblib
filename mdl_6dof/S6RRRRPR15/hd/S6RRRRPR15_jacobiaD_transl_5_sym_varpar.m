% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR15_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR15_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:39:02
% EndTime: 2019-02-26 22:39:03
% DurationCPUTime: 1.18s
% Computational Cost: add. (1291->149), mult. (3990->253), div. (0->0), fcn. (4251->12), ass. (0->97)
t563 = sin(qJ(2));
t564 = sin(qJ(1));
t626 = cos(pkin(6));
t599 = t564 * t626;
t593 = t563 * t599;
t567 = cos(qJ(2));
t568 = cos(qJ(1));
t608 = t568 * t567;
t610 = t564 * t563;
t540 = -qJD(1) * t593 - qJD(2) * t610 + (qJD(2) * t626 + qJD(1)) * t608;
t598 = t568 * t626;
t549 = t563 * t598 + t564 * t567;
t562 = sin(qJ(3));
t566 = cos(qJ(3));
t574 = t568 * t563 + t567 * t599;
t539 = t574 * qJD(1) + t549 * qJD(2);
t558 = sin(pkin(7));
t560 = cos(pkin(7));
t559 = sin(pkin(6));
t607 = qJD(1) * t559;
t602 = t564 * t607;
t576 = -t539 * t560 + t558 * t602;
t548 = -t567 * t598 + t610;
t616 = t559 * t568;
t584 = t548 * t560 + t558 * t616;
t510 = (-qJD(3) * t549 + t576) * t562 + (-t584 * qJD(3) + t540) * t566;
t623 = t539 * t558;
t531 = t560 * t602 + t623;
t561 = sin(qJ(4));
t565 = cos(qJ(4));
t527 = -t549 * t566 + t584 * t562;
t543 = -t548 * t558 + t560 * t616;
t635 = t527 * t565 + t543 * t561;
t640 = t635 * qJD(4) - t510 * t561 + t531 * t565;
t636 = t527 * t561 - t543 * t565;
t639 = t636 * qJD(4) + t510 * t565 + t531 * t561;
t627 = r_i_i_C(3) + qJ(5);
t629 = r_i_i_C(2) - pkin(4);
t573 = t627 * t561 - t629 * t565 + pkin(3);
t630 = r_i_i_C(1) + pkin(11);
t628 = pkin(10) * t558;
t575 = t593 - t608;
t617 = t559 * t564;
t583 = t558 * t617 - t560 * t574;
t529 = t583 * t562 - t566 * t575;
t625 = t529 * t561;
t537 = t548 * qJD(1) + t575 * qJD(2);
t624 = t537 * t558;
t620 = t558 * t561;
t619 = t558 * t565;
t618 = t558 * t567;
t615 = t560 * t562;
t614 = t560 * t566;
t613 = t562 * t563;
t612 = t562 * t567;
t611 = t563 * t566;
t609 = t566 * t567;
t606 = qJD(2) * t559;
t605 = t558 * t559 * t563;
t603 = pkin(10) * t560 + pkin(9);
t601 = t568 * t607;
t600 = t558 * t606;
t597 = t626 * t558;
t596 = t558 * t601;
t595 = t563 * t600;
t594 = t567 * t600;
t592 = qJD(3) * t597;
t545 = t558 * t574 + t560 * t617;
t588 = t529 * t565 + t545 * t561;
t580 = t560 * t612 + t611;
t542 = t580 * t559 + t562 * t597;
t547 = -t559 * t618 + t626 * t560;
t587 = t542 * t565 + t547 * t561;
t535 = -t548 * t566 - t549 * t615;
t586 = -t535 * t561 + t549 * t619;
t536 = -t566 * t574 + t575 * t615;
t585 = -t536 * t561 - t575 * t619;
t582 = t560 * t609 - t613;
t581 = -t560 * t611 - t612;
t579 = t560 * t613 - t609;
t546 = t579 * t559;
t578 = t546 * t561 + t565 * t605;
t577 = t560 * t601 - t624;
t571 = t562 * t575 + t583 * t566;
t570 = qJD(5) * t561 + (t629 * t561 + t627 * t565) * qJD(4);
t569 = t527 * qJD(3) - t540 * t562 + t576 * t566;
t538 = t549 * qJD(1) + t574 * qJD(2);
t533 = (-t580 * qJD(2) + t581 * qJD(3)) * t559;
t522 = t566 * t592 + (-t579 * qJD(2) + t582 * qJD(3)) * t559;
t518 = -t540 * t615 - t539 * t566 + (t548 * t562 - t549 * t614) * qJD(3);
t516 = t538 * t615 + t537 * t566 + (t562 * t574 + t575 * t614) * qJD(3);
t513 = t587 * qJD(4) + t522 * t561 - t565 * t595;
t508 = -t538 * t566 + (t537 * t560 + t596) * t562 + t571 * qJD(3);
t507 = t529 * qJD(3) - t537 * t614 - t538 * t562 - t566 * t596;
t498 = -qJD(4) * t625 + t577 * t561 + (t545 * qJD(4) + t508) * t565;
t497 = t588 * qJD(4) + t508 * t561 - t577 * t565;
t1 = [t636 * qJD(5) - t510 * pkin(3) - t540 * pkin(2) - pkin(10) * t623 + t630 * t569 + t629 * t639 + t627 * t640 + (-t568 * pkin(1) - t603 * t617) * qJD(1), -t585 * qJD(5) + t516 * pkin(3) + t537 * pkin(2) - t538 * t628 + t630 * (t536 * qJD(3) + t537 * t562 - t538 * t614) - t629 * (t585 * qJD(4) + t516 * t565 - t538 * t620) + t627 * (t538 * t619 + t516 * t561 + (t536 * t565 - t575 * t620) * qJD(4)) -t507 * t573 + t630 * t508 + t570 * t571, t588 * qJD(5) + t629 * t497 + t627 * t498, t497, 0; -(t545 * t565 - t625) * qJD(5) + t508 * pkin(3) - t538 * pkin(2) - pkin(10) * t624 + t630 * t507 - t629 * t498 + t627 * t497 + (-t564 * pkin(1) + t603 * t616) * qJD(1), -t586 * qJD(5) + t518 * pkin(3) - t539 * pkin(2) + t540 * t628 + t630 * (t535 * qJD(3) - t539 * t562 + t540 * t614) - t629 * (t586 * qJD(4) + t518 * t565 + t540 * t620) + t627 * (-t540 * t619 + t518 * t561 + (t535 * t565 + t549 * t620) * qJD(4)) t630 * t510 + t570 * (-t549 * t562 - t584 * t566) + t573 * t569, -t635 * qJD(5) + t627 * t639 - t629 * t640, -t640, 0; 0, -t578 * qJD(5) + t533 * pkin(3) - t630 * (-t582 * qJD(2) + t579 * qJD(3)) * t559 - t629 * (t578 * qJD(4) + t533 * t565 + t561 * t594) + t627 * (-t565 * t594 + t533 * t561 + (-t546 * t565 + t561 * t605) * qJD(4)) + (-pkin(2) * t563 + pkin(10) * t618) * t606, t630 * t522 + t570 * (t582 * t559 + t566 * t597) + t573 * (-t562 * t592 + (t581 * qJD(2) - t580 * qJD(3)) * t559) t587 * qJD(5) + t627 * (t561 * t595 + t522 * t565 + (-t542 * t561 + t547 * t565) * qJD(4)) + t629 * t513, t513, 0;];
JaD_transl  = t1;
