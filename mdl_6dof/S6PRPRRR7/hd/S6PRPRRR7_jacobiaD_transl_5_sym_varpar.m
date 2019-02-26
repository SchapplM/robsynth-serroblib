% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR7
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobiaD_transl_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:35
% EndTime: 2019-02-26 19:57:36
% DurationCPUTime: 0.94s
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
t561 = t607 * t586 + t587 * t636;
t634 = t585 * t597;
t575 = -t586 * t634 + t591 * t590;
t584 = sin(pkin(8));
t589 = cos(pkin(8));
t617 = t561 * t589 + t575 * t584;
t531 = t562 * t596 + t617 * t593;
t628 = t591 * t594;
t605 = t583 * t628 - t588 * t597;
t606 = t583 * t627 + t588 * t594;
t638 = t585 * t586;
t608 = t583 * t638 - t590 * t606;
t545 = t608 * t582 - t587 * t605;
t544 = t582 * t605 + t608 * t587;
t633 = t586 * t590;
t564 = t583 * t633 + t585 * t606;
t618 = t544 * t589 + t564 * t584;
t525 = t545 * t596 + t618 * t593;
t577 = t583 * t597 + t588 * t628;
t609 = t576 * t590 - t588 * t638;
t543 = t577 * t587 + t609 * t582;
t542 = -t577 * t582 + t609 * t587;
t563 = -t576 * t585 - t588 * t633;
t619 = t542 * t589 + t563 * t584;
t523 = t543 * t596 + t619 * t593;
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
t616 = t595 * r_i_i_C(1) - t592 * r_i_i_C(2) + pkin(4);
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
t604 = qJD(5) * (-t592 * r_i_i_C(1) - t595 * r_i_i_C(2));
t569 = (-t582 * t597 - t587 * t630) * t586;
t603 = t569 * t589 + t584 * t625;
t570 = (-t582 * t630 + t587 * t597) * t586;
t602 = -t565 * t589 + t584 * t621;
t567 = qJD(2) * t569;
t601 = t567 * t589 + t584 * t620;
t600 = -t543 * t593 + t619 * t596;
t599 = -t545 * t593 + t618 * t596;
t598 = -t562 * t593 + t617 * t596;
t555 = t576 * t587 - t577 * t641;
t528 = t555 * t596 + t611 * t593;
t557 = -t587 * t606 + t605 * t641;
t529 = t557 * t596 + t610 * t593;
t540 = t570 * t596 + t603 * t593;
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
t527 = t566 * t596 + t602 * t593 + (-t570 * t593 + t603 * t596) * qJD(4);
t521 = t598 * qJD(4) + t568 * t596 + t601 * t593;
t519 = t553 * t596 - t612 * t593 + (-t557 * t593 + t610 * t596) * qJD(4);
t517 = t549 * t596 - t614 * t593 + (-t555 * t593 + t611 * t596) * qJD(4);
t515 = t599 * qJD(4) + t551 * t596 + t613 * t593;
t513 = t600 * qJD(4) + t547 * t596 + t615 * t593;
t1 = [0 (t519 * t595 + t537 * t592) * r_i_i_C(1) + (-t519 * t592 + t537 * t595) * r_i_i_C(2) + t519 * pkin(4) + t553 * pkin(3) - pkin(10) * t644 + t574 * pkin(2) + t648 * (t529 * qJD(4) + t553 * t593 + t612 * t596) + ((-t529 * t592 + t539 * t595) * r_i_i_C(1) + (-t529 * t595 - t539 * t592) * r_i_i_C(2)) * qJD(5) + (-qJD(3) * t605 - t622 * t573) * t585, -t574 * t585, t648 * t515 + t599 * t604 + t616 * (-t525 * qJD(4) - t551 * t593 + t613 * t596) (-t515 * t592 + t536 * t595) * r_i_i_C(1) + (-t515 * t595 - t536 * t592) * r_i_i_C(2) + ((-t525 * t595 - t533 * t592) * r_i_i_C(1) + (t525 * t592 - t533 * t595) * r_i_i_C(2)) * qJD(5), 0; 0 (t517 * t595 + t535 * t592) * r_i_i_C(1) + (-t517 * t592 + t535 * t595) * r_i_i_C(2) + t517 * pkin(4) + t549 * pkin(3) - pkin(10) * t645 - t572 * pkin(2) + t648 * (t528 * qJD(4) + t549 * t593 + t614 * t596) + ((-t528 * t592 + t538 * t595) * r_i_i_C(1) + (-t528 * t595 - t538 * t592) * r_i_i_C(2)) * qJD(5) + (t577 * qJD(3) + t622 * t571) * t585, t572 * t585, t648 * t513 + t600 * t604 + t616 * (-t523 * qJD(4) - t547 * t593 + t615 * t596) (-t513 * t592 + t534 * t595) * r_i_i_C(1) + (-t513 * t595 - t534 * t592) * r_i_i_C(2) + ((-t523 * t595 - t532 * t592) * r_i_i_C(1) + (t523 * t592 - t532 * t595) * r_i_i_C(2)) * qJD(5), 0; 0 (t527 * t595 + t558 * t592) * r_i_i_C(1) + (-t527 * t592 + t558 * t595) * r_i_i_C(2) + t527 * pkin(4) + t566 * pkin(3) + pkin(10) * t642 + t648 * (t540 * qJD(4) + t566 * t593 - t602 * t596) + ((-t540 * t592 + t560 * t595) * r_i_i_C(1) + (-t540 * t595 - t560 * t592) * r_i_i_C(2)) * qJD(5) + (qJD(3) * t635 + (-pkin(2) * t594 + t622 * t634) * qJD(2)) * t586, t620, t648 * t521 + t598 * t604 + t616 * (-t531 * qJD(4) - t568 * t593 + t601 * t596) (-t521 * t592 + t559 * t595) * r_i_i_C(1) + (-t521 * t595 - t559 * t592) * r_i_i_C(2) + ((-t531 * t595 - t541 * t592) * r_i_i_C(1) + (t531 * t592 - t541 * t595) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
