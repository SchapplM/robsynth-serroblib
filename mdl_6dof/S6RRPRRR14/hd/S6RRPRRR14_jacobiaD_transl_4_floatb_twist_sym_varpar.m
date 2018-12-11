% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14_jacobiaD_transl_4_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_4_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_4_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobiaD_transl_4_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_transl_4_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:25
% EndTime: 2018-12-10 18:38:26
% DurationCPUTime: 0.88s
% Computational Cost: add. (2167->165), mult. (2764->247), div. (0->0), fcn. (2168->26), ass. (0->112)
t663 = cos(qJ(2));
t664 = cos(qJ(1));
t693 = pkin(6) + qJ(2);
t678 = sin(t693) / 0.2e1;
t694 = pkin(6) - qJ(2);
t685 = sin(t694);
t625 = t678 - t685 / 0.2e1;
t666 = t625 * qJD(2);
t660 = sin(qJ(2));
t661 = sin(qJ(1));
t681 = cos(t693) / 0.2e1;
t687 = cos(t694);
t673 = t687 / 0.2e1 + t681;
t667 = t664 * t660 + t661 * t673;
t698 = qJD(2) * t661;
t589 = qJD(1) * t667 + t663 * t698 + t664 * t666;
t615 = t673 * qJD(2);
t636 = t660 * t698;
t699 = qJD(1) * t661;
t590 = -t625 * t699 - t636 + (qJD(1) * t663 + t615) * t664;
t649 = pkin(7) + pkin(14);
t637 = sin(t649) / 0.2e1;
t650 = pkin(7) - pkin(14);
t647 = sin(t650);
t618 = t637 + t647 / 0.2e1;
t638 = cos(t650) / 0.2e1;
t648 = cos(t649);
t620 = t638 + t648 / 0.2e1;
t651 = sin(pkin(14));
t654 = sin(pkin(6));
t690 = t654 * t699;
t563 = -t589 * t620 - t590 * t651 + t618 * t690;
t619 = t637 - t647 / 0.2e1;
t621 = t638 - t648 / 0.2e1;
t655 = cos(pkin(14));
t564 = -t589 * t619 + t590 * t655 + t621 * t690;
t670 = t664 * t673;
t602 = t660 * t661 - t670;
t603 = t664 * t625 + t661 * t663;
t700 = t654 * t664;
t573 = t602 * t620 + t603 * t651 + t618 * t700;
t574 = t602 * t619 - t603 * t655 + t621 * t700;
t653 = sin(pkin(7));
t657 = cos(pkin(7));
t579 = t589 * t653 + t657 * t690;
t598 = -t602 * t653 + t657 * t700;
t691 = pkin(8) + qJ(4);
t683 = sin(t691);
t676 = t683 / 0.2e1;
t692 = pkin(8) - qJ(4);
t684 = sin(t692);
t677 = t684 / 0.2e1;
t622 = t676 + t677;
t610 = t622 * qJD(4);
t679 = cos(t691) / 0.2e1;
t686 = cos(t692);
t627 = t686 / 0.2e1 + t679;
t612 = t627 * qJD(4);
t623 = t676 - t684 / 0.2e1;
t626 = t679 - t686 / 0.2e1;
t662 = cos(qJ(4));
t659 = sin(qJ(4));
t696 = qJD(4) * t659;
t707 = -t563 * t623 - t564 * t662 + t573 * t612 - t574 * t696 + t579 * t626 + t598 * t610;
t611 = (t677 - t683 / 0.2e1) * qJD(4);
t613 = t626 * qJD(4);
t695 = qJD(4) * t662;
t706 = t563 * t627 - t564 * t659 - t573 * t611 + t574 * t695 + t579 * t622 - t598 * t613;
t705 = r_i_i_C(3) + pkin(11);
t628 = t681 - t687 / 0.2e1;
t697 = qJD(2) * t664;
t587 = qJD(1) * t603 + t661 * t615 + t660 * t697;
t616 = t628 * qJD(2);
t703 = t616 * t653;
t702 = t653 * qJ(3);
t701 = t654 * t661;
t689 = qJD(1) * t700;
t688 = qJ(3) * t657 + pkin(10);
t607 = t661 * t625 - t664 * t663;
t672 = -t610 * r_i_i_C(1) - t613 * r_i_i_C(2) - qJD(3);
t600 = t653 * t667 + t657 * t701;
t652 = sin(pkin(8));
t669 = t623 * r_i_i_C(1) + t627 * r_i_i_C(2) - t652 * t705;
t624 = t678 + t685 / 0.2e1;
t656 = cos(pkin(8));
t668 = t626 * r_i_i_C(1) - t622 * r_i_i_C(2) - t656 * t705 - qJ(3);
t586 = -qJD(1) * t670 + t660 * t699 + t661 * t666 - t663 * t697;
t561 = t586 * t620 + t587 * t651 + t618 * t689;
t562 = t586 * t619 - t587 * t655 + t621 * t689;
t575 = t607 * t651 + t618 * t701 - t620 * t667;
t576 = -t607 * t655 - t619 * t667 + t621 * t701;
t577 = -t586 * t653 + t657 * t689;
t665 = -t561 * t623 - t562 * t662 - t575 * t612 + t576 * t696 + t577 * t626 - t600 * t610;
t658 = cos(pkin(6));
t614 = t624 * qJD(2);
t601 = -t624 * t653 + t657 * t658;
t597 = t619 * t628 + t624 * t655;
t596 = t620 * t628 - t624 * t651;
t595 = -t614 * t619 + t616 * t655;
t593 = t614 * t655 + t616 * t619;
t592 = -t614 * t651 + t616 * t620;
t591 = qJD(1) * t607 - t664 * t615 + t636;
t585 = t619 * t624 + t621 * t658 - t628 * t655;
t584 = t618 * t658 + t620 * t624 + t628 * t651;
t583 = t607 * t619 - t655 * t667;
t582 = t607 * t620 + t651 * t667;
t581 = -t602 * t655 - t603 * t619;
t580 = t602 * t651 - t603 * t620;
t570 = -t589 * t655 + t591 * t619;
t568 = t586 * t655 + t587 * t619;
t560 = t561 * t627 - t562 * t659 + t575 * t611 - t576 * t695 + t577 * t622 + t600 * t613;
t1 = [t707 * r_i_i_C(1) - t706 * r_i_i_C(2) - t564 * pkin(3) - t590 * pkin(2) - t589 * t702 + t598 * qJD(3) + (-t664 * pkin(1) - t688 * t701) * qJD(1) + t705 * (t563 * t652 - t579 * t656) (t568 * t662 + t582 * t612 - t583 * t696) * r_i_i_C(1) + (-t568 * t659 + t582 * t611 - t583 * t695) * r_i_i_C(2) + t568 * pkin(3) + t586 * pkin(2) + t669 * (-t586 * t651 + t587 * t620) + (t587 * t668 + t607 * t672) * t653, t577, t560 * r_i_i_C(1) + r_i_i_C(2) * t665, 0, 0; -t665 * r_i_i_C(1) + t560 * r_i_i_C(2) + t562 * pkin(3) - t587 * pkin(2) - t586 * t702 + t600 * qJD(3) + (-t661 * pkin(1) + t688 * t700) * qJD(1) + t705 * (-t561 * t652 + t577 * t656) (t570 * t662 + t580 * t612 - t581 * t696) * r_i_i_C(1) + (-t570 * t659 + t580 * t611 - t581 * t695) * r_i_i_C(2) + t570 * pkin(3) - t589 * pkin(2) + t669 * (t589 * t651 + t591 * t620) + (t591 * t668 - t603 * t672) * t653, t579, t706 * r_i_i_C(1) + r_i_i_C(2) * t707, 0, 0; 0 (t595 * t662 + t596 * t612 - t597 * t696) * r_i_i_C(1) + (-t595 * t659 + t596 * t611 - t597 * t695) * r_i_i_C(2) + t595 * pkin(3) + t616 * pkin(2) + t669 * (-t614 * t620 - t616 * t651) + (-t614 * t668 + t628 * t672) * t653, -t703 (t584 * t611 - t585 * t695 + t592 * t627 - t593 * t659 + t601 * t613 - t622 * t703) * r_i_i_C(1) + (-t584 * t612 + t585 * t696 - t592 * t623 - t593 * t662 - t601 * t610 - t626 * t703) * r_i_i_C(2), 0, 0;];
JaD_transl  = t1;
