% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP12_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP12_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:14:25
% EndTime: 2019-02-26 21:14:26
% DurationCPUTime: 1.55s
% Computational Cost: add. (2090->149), mult. (6608->243), div. (0->0), fcn. (7481->14), ass. (0->98)
t719 = cos(pkin(12));
t721 = cos(pkin(6));
t696 = t721 * t719;
t717 = sin(pkin(12));
t724 = sin(qJ(1));
t726 = cos(qJ(1));
t662 = t724 * t696 + t726 * t717;
t631 = t662 * qJD(1);
t694 = t721 * t717;
t636 = -t724 * t694 + t726 * t719;
t632 = t636 * qJD(1);
t635 = t694 * t726 + t719 * t724;
t653 = sin(qJ(3));
t650 = sin(pkin(6));
t718 = sin(pkin(7));
t707 = t650 * t718;
t681 = t724 * t707;
t720 = cos(pkin(7));
t725 = cos(qJ(3));
t634 = -t696 * t726 + t717 * t724;
t683 = t726 * t707;
t674 = t725 * t683;
t703 = t720 * t725;
t730 = -t634 * t703 - t674;
t597 = t653 * (qJD(1) * t681 - qJD(3) * t635 - t631 * t720) + t632 * t725 + t730 * qJD(3);
t702 = t720 * t724;
t682 = t650 * t702;
t708 = t631 * t718;
t620 = qJD(1) * t682 + t708;
t652 = sin(qJ(4));
t655 = cos(qJ(4));
t706 = t720 * t634;
t711 = t635 * t725;
t612 = t653 * (t683 + t706) - t711;
t710 = t650 * t726;
t624 = t634 * t718 - t720 * t710;
t687 = t612 * t652 + t624 * t655;
t587 = -qJD(4) * t687 - t597 * t655 - t620 * t652;
t673 = t725 * t681;
t680 = t653 * t683;
t598 = qJD(1) * t673 - t631 * t703 - t632 * t653 + (t653 * t706 + t680 - t711) * qJD(3);
t651 = sin(qJ(5));
t654 = cos(qJ(5));
t604 = t612 * t655 - t624 * t652;
t609 = t635 * t653 - t730;
t743 = t604 * t654 - t609 * t651;
t748 = t743 * qJD(5) + t587 * t651 - t598 * t654;
t744 = t604 * t651 + t609 * t654;
t747 = t744 * qJD(5) - t587 * t654 - t598 * t651;
t742 = qJD(4) * t604 - t597 * t652 + t620 * t655;
t728 = r_i_i_C(1) + pkin(5);
t727 = r_i_i_C(2) + pkin(11);
t723 = r_i_i_C(3) + qJ(6);
t660 = t662 * t720;
t729 = -t636 * t653 - t725 * t660 + t673;
t693 = t720 * t719;
t695 = t721 * t718;
t621 = -t725 * t695 + (t653 * t717 - t693 * t725) * t650;
t656 = t624 * qJD(1);
t671 = t723 * t651 + t654 * t728 + pkin(4);
t625 = t662 * t718 + t682;
t716 = t625 * t655;
t714 = t651 * t655;
t713 = qJD(4) * t652;
t712 = qJD(5) * t655;
t709 = qJD(1) * t724;
t705 = qJD(1) * t710;
t630 = t635 * qJD(1);
t658 = qJD(1) * t706;
t595 = qJD(1) * t680 + t729 * qJD(3) - t630 * t725 + t653 * t658;
t699 = -t712 * t729 + t595;
t698 = t609 * t712 + t597;
t617 = t621 * qJD(3);
t697 = t621 * t712 - t617;
t614 = t636 * t725 + (-t660 + t681) * t653;
t606 = t614 * t655 + t625 * t652;
t690 = t606 * t654 - t651 * t729;
t689 = -t606 * t651 - t654 * t729;
t622 = t653 * t695 + (t653 * t693 + t717 * t725) * t650;
t633 = -t707 * t719 + t720 * t721;
t608 = t622 * t655 + t633 * t652;
t688 = t608 * t654 + t621 * t651;
t686 = -t622 * t652 + t633 * t655;
t679 = -t655 * pkin(4) - t652 * t727 - pkin(3);
t675 = qJD(4) * (pkin(4) * t652 - t655 * t727);
t594 = -qJD(1) * t674 + qJD(3) * t614 - t630 * t653 - t658 * t725;
t669 = qJD(5) * t614 - t594 * t655 - t713 * t729;
t668 = -qJD(5) * t612 + t598 * t655 + t609 * t713;
t618 = t622 * qJD(3);
t667 = qJD(5) * t622 - t618 * t655 + t621 * t713;
t663 = qJD(6) * t651 + (-t651 * t728 + t723 * t654) * qJD(5);
t601 = qJD(4) * t686 - t617 * t655;
t588 = qJD(5) * t688 + t601 * t651 - t618 * t654;
t583 = qJD(4) * t716 + t595 * t655 - t614 * t713 - t652 * t656;
t582 = qJD(4) * t606 + t595 * t652 + t655 * t656;
t573 = qJD(5) * t689 + t583 * t654 + t594 * t651;
t572 = qJD(5) * t690 + t583 * t651 - t594 * t654;
t1 = [t744 * qJD(6) + t587 * pkin(4) - t597 * pkin(3) + t598 * pkin(10) - t632 * pkin(2) - pkin(9) * t708 + qJD(2) * t710 + t727 * t742 - t728 * t747 + t723 * t748 + (-t726 * pkin(1) + (-pkin(9) * t702 - qJ(2) * t724) * t650) * qJD(1), t705 -(t614 * t654 - t714 * t729) * qJD(6) + t595 * pkin(10) + t728 * (t651 * t699 + t654 * t669) + t723 * (t651 * t669 - t654 * t699) - t729 * t675 + t679 * t594, t727 * t583 + t663 * (-t614 * t652 + t716) - t671 * t582, t690 * qJD(6) - t572 * t728 + t723 * t573, t572; t724 * t650 * qJD(2) - pkin(1) * t709 - t630 * pkin(2) + t595 * pkin(3) + t583 * pkin(4) - pkin(9) * t656 + t594 * pkin(10) + qJ(2) * t705 - t689 * qJD(6) + t723 * t572 + t728 * t573 + t727 * t582, t650 * t709 -(t609 * t714 - t612 * t654) * qJD(6) + t597 * pkin(10) + t728 * (t651 * t698 + t654 * t668) + t723 * (t651 * t668 - t654 * t698) + t609 * t675 - t679 * t598, -t587 * t727 + t663 * t687 + t671 * t742, -t743 * qJD(6) + t723 * t747 + t728 * t748, -t748; 0, 0 -(t621 * t714 + t622 * t654) * qJD(6) - t617 * pkin(10) + t728 * (t651 * t697 + t654 * t667) + t723 * (t651 * t667 - t654 * t697) + t621 * t675 + t679 * t618, t727 * t601 + t663 * t686 + t671 * (-qJD(4) * t608 + t617 * t652) t688 * qJD(6) + t723 * (t601 * t654 + t618 * t651 + (-t608 * t651 + t621 * t654) * qJD(5)) - t728 * t588, t588;];
JaD_transl  = t1;
