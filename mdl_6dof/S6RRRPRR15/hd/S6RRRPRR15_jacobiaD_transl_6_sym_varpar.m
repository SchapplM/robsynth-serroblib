% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR15_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR15_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:21
% EndTime: 2019-02-26 22:24:23
% DurationCPUTime: 1.70s
% Computational Cost: add. (2098->210), mult. (6474->344), div. (0->0), fcn. (7012->14), ass. (0->124)
t690 = cos(qJ(2));
t691 = cos(qJ(1));
t752 = cos(pkin(6));
t723 = t691 * t752;
t686 = sin(qJ(2));
t753 = sin(qJ(1));
t727 = t753 * t686;
t663 = -t690 * t723 + t727;
t664 = t686 * t723 + t690 * t753;
t685 = sin(qJ(3));
t689 = cos(qJ(3));
t680 = sin(pkin(7));
t681 = sin(pkin(6));
t743 = t681 * t691;
t732 = t680 * t743;
t682 = cos(pkin(7));
t742 = t682 * t685;
t633 = t663 * t742 - t664 * t689 + t685 * t732;
t716 = t752 * t753;
t694 = t691 * t686 + t690 * t716;
t650 = qJD(1) * t694 + qJD(2) * t664;
t707 = t686 * t716;
t736 = t691 * t690;
t651 = -qJD(1) * t707 - qJD(2) * t727 + (qJD(2) * t752 + qJD(1)) * t736;
t729 = t681 * t753;
t717 = qJD(1) * t729;
t709 = t680 * t717;
t741 = t682 * t689;
t609 = t633 * qJD(3) - t650 * t741 - t651 * t685 + t689 * t709;
t750 = t650 * t680;
t638 = t682 * t717 + t750;
t684 = sin(qJ(5));
t688 = cos(qJ(5));
t767 = t609 * t684 - t638 * t688;
t766 = -t609 * t688 - t638 * t684;
t748 = t664 * t685;
t757 = t689 * (t663 * t682 + t732) + t748;
t765 = t757 * qJD(3) - (-t650 * t682 + t709) * t685 - t651 * t689;
t683 = sin(qJ(6));
t764 = t633 * t683;
t687 = cos(qJ(6));
t763 = t633 * t687;
t654 = -t663 * t680 + t682 * t743;
t760 = t654 * t684;
t759 = t654 * t688;
t734 = (pkin(4) + pkin(10)) * t680;
t756 = pkin(3) + pkin(11);
t754 = r_i_i_C(3) + pkin(12);
t695 = t707 - t736;
t648 = qJD(1) * t663 + qJD(2) * t695;
t751 = t648 * t680;
t747 = t695 * t685;
t746 = t680 * t681;
t745 = t680 * t684;
t744 = t680 * t688;
t740 = t685 * t686;
t739 = t685 * t690;
t738 = t686 * t689;
t737 = t689 * t690;
t735 = qJD(2) * t681;
t733 = t686 * t746;
t731 = t681 * t740;
t730 = t682 * t737;
t728 = t682 * t753;
t726 = qJD(1) * t743;
t725 = t680 * t735;
t724 = t680 * t752;
t722 = t681 * t730;
t721 = t680 * t729;
t720 = t680 * t726;
t719 = t686 * t725;
t718 = t690 * t725;
t715 = qJD(3) * t724;
t714 = -t687 * r_i_i_C(1) + t683 * r_i_i_C(2);
t713 = -t683 * r_i_i_C(1) - t687 * r_i_i_C(2);
t630 = t663 * t741 + t689 * t732 + t748;
t712 = t630 * t688 + t760;
t621 = t630 * t684 - t759;
t619 = -t684 * t757 + t759;
t634 = -t689 * t721 + t694 * t741 - t747;
t656 = t680 * t694 + t681 * t728;
t711 = t634 * t688 - t656 * t684;
t623 = t634 * t684 + t656 * t688;
t652 = -t689 * t724 - t722 + t731;
t662 = t682 * t752 - t690 * t746;
t710 = t652 * t688 - t662 * t684;
t629 = t652 * t684 + t662 * t688;
t708 = pkin(5) - t714;
t706 = -t713 + t756;
t643 = -t663 * t685 + t664 * t741;
t624 = t643 * t684 + t664 * t744;
t645 = -t685 * t694 - t695 * t741;
t625 = t645 * t684 - t695 * t744;
t644 = -t663 * t689 - t664 * t742;
t646 = -t689 * t694 + t695 * t742;
t703 = t682 * t738 + t739;
t702 = t682 * t739 + t738;
t701 = -t682 * t740 + t737;
t700 = qJD(6) * t714;
t699 = qJD(6) * t713;
t659 = t703 * t681;
t647 = t659 * t684 + t688 * t733;
t698 = -t682 * t694 + t721;
t635 = t685 * t698 - t689 * t695;
t693 = t684 * t708 - t688 * t754 + qJ(4);
t692 = qJD(4) + t684 * t699 + (t684 * t754 + t688 * t708) * qJD(5);
t660 = t701 * t681;
t653 = t681 * t702 + t685 * t724;
t649 = qJD(1) * t664 + qJD(2) * t694;
t639 = -qJD(2) * t722 - t681 * qJD(3) * t737 + (qJD(3) * t682 + qJD(2)) * t731;
t636 = t682 * t726 - t751;
t627 = t689 * t715 + ((t730 - t740) * qJD(3) + t701 * qJD(2)) * t681;
t626 = t685 * t715 + (qJD(2) * t703 + qJD(3) * t702) * t681;
t615 = qJD(3) * t644 - t650 * t685 + t651 * t741;
t613 = qJD(3) * t646 + t648 * t685 - t649 * t741;
t612 = qJD(5) * t710 + t626 * t684 + t688 * t719;
t606 = qJD(3) * t747 + (t648 * t682 + t720) * t685 + (qJD(3) * t698 - t649) * t689;
t605 = qJD(3) * t635 - t648 * t741 - t649 * t685 - t689 * t720;
t599 = qJD(5) * t712 - t767;
t595 = qJD(5) * t711 + t605 * t684 + t636 * t688;
t594 = qJD(5) * t623 - t605 * t688 + t636 * t684;
t593 = t595 * t687 + t606 * t683 + (-t623 * t683 + t635 * t687) * qJD(6);
t592 = -t595 * t683 + t606 * t687 + (-t623 * t687 - t635 * t683) * qJD(6);
t1 = [-pkin(10) * t750 - t651 * pkin(2) - t638 * pkin(4) + t609 * qJ(4) - t757 * qJD(4) + t708 * ((-t688 * t757 - t760) * qJD(5) + t767) + t706 * t765 + t754 * (qJD(5) * t619 + t766) + ((-t619 * t683 + t763) * r_i_i_C(1) + (-t619 * t687 - t764) * r_i_i_C(2)) * qJD(6) + (-t691 * pkin(1) + (-pkin(9) * t753 - pkin(10) * t728) * t681) * qJD(1), t648 * pkin(2) + t613 * qJ(4) + t645 * qJD(4) - t649 * t734 + t708 * (-t649 * t744 + t613 * t684 + (t645 * t688 + t695 * t745) * qJD(5)) + t706 * (-qJD(3) * t645 + t648 * t689 + t649 * t742) - t754 * (-qJD(5) * t625 + t613 * t688 + t649 * t745) + ((-t625 * t683 + t646 * t687) * r_i_i_C(1) + (-t625 * t687 - t646 * t683) * r_i_i_C(2)) * qJD(6), -t605 * t706 + t606 * t693 + t634 * t700 + t635 * t692, t605, -t594 * t708 + t595 * t754 + t699 * t711, r_i_i_C(1) * t592 - t593 * r_i_i_C(2); -pkin(10) * t751 - t649 * pkin(2) + t636 * pkin(4) + t595 * pkin(5) + t593 * r_i_i_C(1) + t592 * r_i_i_C(2) + t605 * qJ(4) + t634 * qJD(4) + t756 * t606 + t754 * t594 + (-t753 * pkin(1) + (pkin(10) * t682 + pkin(9)) * t743) * qJD(1), -t650 * pkin(2) + t615 * qJ(4) + t643 * qJD(4) + t651 * t734 + t708 * (t651 * t744 + t615 * t684 + (t643 * t688 - t664 * t745) * qJD(5)) + t706 * (-qJD(3) * t643 - t650 * t689 - t651 * t742) - t754 * (-qJD(5) * t624 + t615 * t688 - t651 * t745) + ((-t624 * t683 + t644 * t687) * r_i_i_C(1) + (-t624 * t687 - t644 * t683) * r_i_i_C(2)) * qJD(6), t609 * t706 + t630 * t700 - t633 * t692 - t693 * t765, -t609, t754 * t599 + t712 * t699 + t708 * (-qJD(5) * t621 + t766) (-t599 * t683 - t687 * t765) * r_i_i_C(1) + (-t599 * t687 + t683 * t765) * r_i_i_C(2) + ((-t621 * t687 + t764) * r_i_i_C(1) + (t621 * t683 + t763) * r_i_i_C(2)) * qJD(6); 0, -t639 * qJ(4) + t659 * qJD(4) + t708 * (t688 * t718 - t639 * t684 + (t659 * t688 - t684 * t733) * qJD(5)) - t706 * (qJD(2) * t702 + qJD(3) * t703) * t681 + t754 * (qJD(5) * t647 + t639 * t688 + t684 * t718) + ((-t647 * t683 + t660 * t687) * r_i_i_C(1) + (-t647 * t687 - t660 * t683) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t686 + t690 * t734) * t735, -t626 * t706 + t627 * t693 + t652 * t700 + t653 * t692, t626, t754 * t612 + t710 * t699 + t708 * (-qJD(5) * t629 + t626 * t688 - t684 * t719) (-t612 * t683 + t627 * t687) * r_i_i_C(1) + (-t612 * t687 - t627 * t683) * r_i_i_C(2) + ((-t629 * t687 - t653 * t683) * r_i_i_C(1) + (t629 * t683 - t653 * t687) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
