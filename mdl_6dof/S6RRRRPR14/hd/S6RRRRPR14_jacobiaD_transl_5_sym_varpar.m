% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR14
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
% Datum: 2019-02-26 22:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR14_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR14_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:18
% EndTime: 2019-02-26 22:38:19
% DurationCPUTime: 1.52s
% Computational Cost: add. (1597->175), mult. (4956->293), div. (0->0), fcn. (5291->14), ass. (0->112)
t609 = sin(qJ(2));
t610 = sin(qJ(1));
t680 = cos(pkin(6));
t650 = t610 * t680;
t644 = t609 * t650;
t613 = cos(qJ(2));
t614 = cos(qJ(1));
t661 = t614 * t613;
t663 = t610 * t609;
t581 = -qJD(1) * t644 - qJD(2) * t663 + (qJD(2) * t680 + qJD(1)) * t661;
t649 = t614 * t680;
t589 = -t613 * t649 + t663;
t590 = t609 * t649 + t610 * t613;
t608 = sin(qJ(3));
t612 = cos(qJ(3));
t623 = t614 * t609 + t613 * t650;
t580 = t623 * qJD(1) + t590 * qJD(2);
t603 = sin(pkin(7));
t606 = cos(pkin(7));
t604 = sin(pkin(6));
t660 = qJD(1) * t610;
t655 = t604 * t660;
t625 = -t580 * t606 + t603 * t655;
t657 = qJD(3) * t612;
t652 = t604 * t657;
t645 = t603 * t652;
t658 = qJD(3) * t606;
t553 = (-qJD(3) * t590 + t625) * t608 - (t589 * t658 - t581) * t612 - t614 * t645;
t607 = sin(qJ(4));
t611 = cos(qJ(4));
t678 = t580 * t603;
t626 = t606 * t655 + t678;
t669 = t604 * t614;
t633 = t589 * t606 + t603 * t669;
t570 = -t590 * t612 + t633 * t608;
t584 = -t589 * t603 + t606 * t669;
t690 = t570 * t611 + t584 * t607;
t695 = t690 * qJD(4) - t553 * t607 + t626 * t611;
t691 = t570 * t607 - t584 * t611;
t694 = t691 * qJD(4) + t553 * t611 + t626 * t607;
t681 = r_i_i_C(3) + qJ(5);
t685 = pkin(10) * t606 + pkin(9);
t624 = t644 - t661;
t670 = t604 * t610;
t675 = t623 * t606;
t632 = t603 * t670 - t675;
t572 = t632 * t608 - t612 * t624;
t586 = t603 * t623 + t606 * t670;
t684 = -t572 * t607 + t586 * t611;
t602 = sin(pkin(13));
t605 = cos(pkin(13));
t637 = t605 * r_i_i_C(1) - t602 * r_i_i_C(2) + pkin(4);
t620 = t681 * t607 + t637 * t611 + pkin(3);
t683 = pkin(10) * t603;
t674 = t624 * t608;
t673 = t603 * t607;
t672 = t603 * t611;
t671 = t603 * t613;
t668 = t606 * t608;
t667 = t606 * t612;
t666 = t608 * t609;
t665 = t608 * t613;
t664 = t609 * t612;
t662 = t612 * t613;
t659 = qJD(2) * t604;
t656 = t603 * t604 * t609;
t654 = qJD(1) * t669;
t653 = t613 * t659;
t651 = t603 * t680;
t648 = t603 * t654;
t647 = qJD(2) * t656;
t646 = t603 * t653;
t643 = qJD(3) * t651;
t639 = t572 * t611 + t586 * t607;
t629 = t606 * t665 + t664;
t583 = t629 * t604 + t608 * t651;
t588 = -t604 * t671 + t680 * t606;
t638 = t583 * t611 + t588 * t607;
t636 = t602 * r_i_i_C(1) + t605 * r_i_i_C(2) + pkin(11);
t577 = -t589 * t612 - t590 * t668;
t635 = -t577 * t607 + t590 * t672;
t578 = -t612 * t623 + t624 * t668;
t634 = -t578 * t607 - t624 * t672;
t631 = t606 * t662 - t666;
t630 = -t606 * t664 - t665;
t628 = -t606 * t666 + t662;
t587 = t628 * t604;
t627 = -t587 * t607 + t611 * t656;
t622 = t624 * qJD(2);
t619 = t607 * qJD(5) + (-t637 * t607 + t681 * t611) * qJD(4);
t618 = t589 * qJD(1) + t622;
t617 = t618 * t608;
t616 = t618 * t612;
t554 = t570 * qJD(3) - t581 * t608 + t625 * t612;
t615 = t584 * qJD(1) - t603 * t622;
t579 = t590 * qJD(1) + t623 * qJD(2);
t575 = (-t629 * qJD(2) + t630 * qJD(3)) * t604;
t574 = -t653 * t667 - t613 * t652 + (qJD(2) + t658) * t604 * t666;
t565 = t612 * t643 + (t628 * qJD(2) + t631 * qJD(3)) * t604;
t563 = t627 * qJD(4) + t575 * t611 + t607 * t646;
t561 = -t581 * t668 - t580 * t612 + (t589 * t608 - t590 * t667) * qJD(3);
t560 = t577 * qJD(3) - t580 * t608 + t581 * t667;
t559 = t616 + t579 * t668 + (t608 * t623 + t624 * t667) * qJD(3);
t558 = t578 * qJD(3) - t579 * t667 + t617;
t556 = t638 * qJD(4) + t565 * t607 - t611 * t647;
t551 = qJD(3) * t674 - t579 * t612 + t606 * t617 + t608 * t648 + t610 * t645 - t657 * t675;
t550 = t572 * qJD(3) - t579 * t608 - t606 * t616 - t612 * t648;
t549 = t635 * qJD(4) + t561 * t611 + t581 * t673;
t547 = t634 * qJD(4) + t559 * t611 - t579 * t673;
t541 = t684 * qJD(4) + t551 * t611 + t615 * t607;
t540 = t639 * qJD(4) + t551 * t607 - t615 * t611;
t1 = [(t554 * t602 - t605 * t694) * r_i_i_C(1) + (t554 * t605 + t602 * t694) * r_i_i_C(2) - t694 * pkin(4) + t691 * qJD(5) - t553 * pkin(3) + t554 * pkin(11) - t581 * pkin(2) - pkin(10) * t678 + t681 * t695 + (-t614 * pkin(1) - t685 * t670) * qJD(1) (t547 * t605 + t558 * t602) * r_i_i_C(1) + (-t547 * t602 + t558 * t605) * r_i_i_C(2) + t547 * pkin(4) - t634 * qJD(5) + t559 * pkin(3) + t558 * pkin(11) + t618 * pkin(2) - t579 * t683 + t681 * (t579 * t672 + t559 * t607 + (t578 * t611 - t624 * t673) * qJD(4)) t636 * t551 + t619 * (t632 * t612 + t674) - t620 * t550, t639 * qJD(5) - t637 * t540 + t681 * t541, t540, 0; (t541 * t605 + t550 * t602) * r_i_i_C(1) + (-t541 * t602 + t550 * t605) * r_i_i_C(2) + t541 * pkin(4) - t684 * qJD(5) + t551 * pkin(3) + t550 * pkin(11) - t579 * pkin(2) - t618 * t683 - pkin(1) * t660 + t685 * t654 + t681 * t540 (t549 * t605 + t560 * t602) * r_i_i_C(1) + (-t549 * t602 + t560 * t605) * r_i_i_C(2) + t549 * pkin(4) - t635 * qJD(5) + t561 * pkin(3) + t560 * pkin(11) - t580 * pkin(2) + t581 * t683 + t681 * (-t581 * t672 + t561 * t607 + (t577 * t611 + t590 * t673) * qJD(4)) t636 * t553 + t619 * (-t590 * t608 - t633 * t612) + t620 * t554, -t690 * qJD(5) + t637 * t695 + t681 * t694, -t695, 0; 0 (t563 * t605 - t574 * t602) * r_i_i_C(1) + (-t563 * t602 - t574 * t605) * r_i_i_C(2) + t563 * pkin(4) - t627 * qJD(5) + t575 * pkin(3) - t574 * pkin(11) + t681 * (-t611 * t646 + t575 * t607 + (t587 * t611 + t607 * t656) * qJD(4)) + (-pkin(2) * t609 + pkin(10) * t671) * t659, t636 * t565 + t619 * (t631 * t604 + t612 * t651) + t620 * (-t608 * t643 + (t630 * qJD(2) - t629 * qJD(3)) * t604) t638 * qJD(5) + t681 * (t607 * t647 + t565 * t611 + (-t583 * t607 + t588 * t611) * qJD(4)) - t637 * t556, t556, 0;];
JaD_transl  = t1;
