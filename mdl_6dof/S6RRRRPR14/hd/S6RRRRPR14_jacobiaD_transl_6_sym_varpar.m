% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JaD_transl = S6RRRRPR14_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR14_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:19
% EndTime: 2019-02-26 22:38:21
% DurationCPUTime: 2.14s
% Computational Cost: add. (2277->209), mult. (6623->346), div. (0->0), fcn. (7210->16), ass. (0->131)
t688 = cos(qJ(2));
t689 = cos(qJ(1));
t765 = cos(pkin(6));
t736 = t689 * t765;
t686 = sin(qJ(2));
t769 = sin(qJ(1));
t742 = t769 * t686;
t657 = -t688 * t736 + t742;
t682 = cos(pkin(7));
t685 = sin(qJ(3));
t658 = t686 * t736 + t688 * t769;
t770 = cos(qJ(3));
t747 = t658 * t770;
t680 = sin(pkin(7));
t681 = sin(pkin(6));
t756 = t681 * t689;
t749 = t680 * t756;
t631 = t685 * (t657 * t682 + t749) - t747;
t649 = -t657 * t680 + t682 * t756;
t684 = sin(qJ(4));
t687 = cos(qJ(4));
t617 = t631 * t687 + t649 * t684;
t723 = t770 * t749;
t746 = t682 * t770;
t628 = t657 * t746 + t658 * t685 + t723;
t678 = pkin(13) + qJ(6);
t676 = sin(t678);
t677 = cos(t678);
t784 = -t617 * t676 - t628 * t677;
t783 = t617 * t677 - t628 * t676;
t727 = t765 * t769;
t700 = t689 * t686 + t688 * t727;
t645 = qJD(1) * t700 + qJD(2) * t658;
t718 = t686 * t727;
t752 = t689 * t688;
t646 = -qJD(1) * t718 - qJD(2) * t742 + (qJD(2) * t765 + qJD(1)) * t752;
t738 = qJD(3) * t770;
t729 = t682 * t738;
t739 = qJD(1) * t769;
t730 = t681 * t739;
t604 = t685 * (-qJD(3) * t658 - t645 * t682 + t680 * t730) - qJD(3) * t723 + t646 * t770 - t657 * t729;
t763 = t645 * t680;
t702 = t682 * t730 + t763;
t782 = qJD(4) * t617 - t604 * t684 + t687 * t702;
t722 = t631 * t684 - t649 * t687;
t594 = qJD(4) * t722 + t604 * t687 + t684 * t702;
t773 = pkin(11) + pkin(5) * sin(pkin(13));
t766 = r_i_i_C(3) + pkin(12) + qJ(5);
t701 = t718 - t752;
t760 = t680 * t681;
t733 = t769 * t760;
t711 = t770 * t733;
t772 = t685 * t701 + t711;
t633 = -t701 * t770 + (-t682 * t700 + t733) * t685;
t745 = t682 * t769;
t651 = t680 * t700 + t681 * t745;
t618 = -t633 * t684 + t651 * t687;
t725 = r_i_i_C(1) * t676 + r_i_i_C(2) * t677;
t712 = qJD(6) * t725;
t703 = t725 + t773;
t735 = t685 * t749;
t755 = t682 * t685;
t771 = qJD(1) * t711 - t645 * t746 - t646 * t685 + (t657 * t755 + t735 - t747) * qJD(3);
t767 = pkin(10) * t680;
t759 = t680 * t684;
t758 = t680 * t687;
t757 = t681 * t688;
t754 = t685 * t686;
t753 = t685 * t688;
t751 = qJD(2) * t681;
t750 = t686 * t760;
t748 = t681 * t754;
t744 = t770 * t686;
t743 = t770 * t688;
t740 = t680 * t751;
t737 = t680 * t765;
t734 = t682 * t743;
t732 = t686 * t740;
t731 = t688 * t740;
t728 = t685 * t737;
t726 = r_i_i_C(1) * t677 - r_i_i_C(2) * t676;
t724 = t681 * t734;
t619 = t633 * t687 + t651 * t684;
t705 = t682 * t753 + t744;
t648 = t681 * t705 + t728;
t656 = -t680 * t757 + t682 * t765;
t625 = t648 * t687 + t656 * t684;
t720 = -t648 * t684 + t656 * t687;
t719 = t770 * t737;
t675 = cos(pkin(13)) * pkin(5) + pkin(4);
t717 = t675 + t726;
t640 = -t657 * t770 - t658 * t755;
t716 = -t640 * t684 + t658 * t758;
t620 = t640 * t687 + t658 * t759;
t642 = -t700 * t770 + t701 * t755;
t715 = -t642 * t684 - t701 * t758;
t621 = t642 * t687 - t701 * t759;
t713 = qJD(6) * t726;
t706 = -t682 * t754 + t743;
t655 = t706 * t681;
t710 = -t655 * t684 + t687 * t750;
t643 = t655 * t687 + t684 * t750;
t708 = t657 * t685 - t658 * t746;
t707 = t685 * t700 + t701 * t746;
t704 = t682 * t744 + t753;
t698 = t701 * qJD(2);
t695 = -t684 * t766 - t687 * t717 - pkin(3);
t694 = qJD(1) * t657 + t698;
t693 = t694 * t685;
t692 = -t684 * qJD(5) + t687 * t712 + (t684 * t717 - t687 * t766) * qJD(4);
t691 = t694 * t770;
t690 = qJD(1) * t649 - t680 * t698;
t654 = t704 * t681;
t647 = -t719 - t724 + t748;
t644 = qJD(1) * t658 + qJD(2) * t700;
t636 = (-qJD(2) * t705 - qJD(3) * t704) * t681;
t632 = t700 * t746 - t772;
t623 = qJD(3) * t719 + ((t734 - t754) * qJD(3) + t706 * qJD(2)) * t681;
t622 = qJD(3) * t728 + (qJD(2) * t704 + qJD(3) * t705) * t681;
t612 = t708 * qJD(3) - t645 * t770 - t646 * t755;
t610 = qJD(3) * t707 + t644 * t755 + t691;
t608 = qJD(4) * t720 + t623 * t687 + t684 * t732;
t607 = qJD(4) * t625 + t623 * t684 - t687 * t732;
t602 = qJD(1) * t735 + t772 * qJD(3) - t644 * t770 + t682 * t693 - t700 * t729;
t601 = -qJD(1) * t723 + qJD(3) * t633 - t644 * t685 - t682 * t691;
t598 = qJD(4) * t715 + t610 * t687 - t644 * t759;
t592 = t618 * qJD(4) + t602 * t687 + t690 * t684;
t591 = qJD(4) * t619 + t602 * t684 - t687 * t690;
t590 = t592 * t677 + t601 * t676 + (-t619 * t676 + t632 * t677) * qJD(6);
t589 = -t592 * t676 + t601 * t677 + (-t619 * t677 - t632 * t676) * qJD(6);
t1 = [t722 * qJD(5) - t604 * pkin(3) - t646 * pkin(2) - pkin(10) * t763 - t717 * t594 + t703 * t771 + t766 * t782 + (t784 * r_i_i_C(1) - t783 * r_i_i_C(2)) * qJD(6) + (-t689 * pkin(1) + (-pkin(9) * t769 - pkin(10) * t745) * t681) * qJD(1) (t598 * t677 + (-t621 * t676 - t677 * t707) * qJD(6)) * r_i_i_C(1) + (-t598 * t676 + (-t621 * t677 + t676 * t707) * qJD(6)) * r_i_i_C(2) + t598 * t675 - t715 * qJD(5) + t610 * pkin(3) + t694 * pkin(2) - t644 * t767 + t703 * (qJD(3) * t642 - t644 * t746 + t693) + t766 * (qJD(4) * t621 + t610 * t684 + t644 * t758) t601 * t695 + t602 * t703 + t632 * t692 + t633 * t713, t619 * qJD(5) - t591 * t717 + t592 * t766 - t618 * t712, t591, r_i_i_C(1) * t589 - t590 * r_i_i_C(2); -pkin(1) * t739 - t644 * pkin(2) + t602 * pkin(3) + t590 * r_i_i_C(1) + t589 * r_i_i_C(2) - t618 * qJD(5) + t592 * t675 - t694 * t767 + (pkin(10) * t682 + pkin(9)) * qJD(1) * t756 + t773 * t601 + t766 * t591, -t716 * qJD(5) + t612 * pkin(3) - t645 * pkin(2) + t646 * t767 + t717 * (qJD(4) * t716 + t612 * t687 + t646 * t759) + t703 * (qJD(3) * t640 - t645 * t685 + t646 * t746) + t766 * (qJD(4) * t620 + t612 * t684 - t646 * t758) + ((-t620 * t676 - t677 * t708) * r_i_i_C(1) + (-t620 * t677 + t676 * t708) * r_i_i_C(2)) * qJD(6), t604 * t703 + t628 * t692 - t631 * t713 - t695 * t771, -qJD(5) * t617 + t594 * t766 - t712 * t722 + t717 * t782, -t782 (-t594 * t676 - t677 * t771) * r_i_i_C(1) + (-t594 * t677 + t676 * t771) * r_i_i_C(2) + (t783 * r_i_i_C(1) + t784 * r_i_i_C(2)) * qJD(6); 0, -t710 * qJD(5) + t636 * pkin(3) + t717 * (qJD(4) * t710 + t636 * t687 + t684 * t731) - t703 * (-qJD(2) * t724 - t738 * t757 + (qJD(3) * t682 + qJD(2)) * t748) + t766 * (qJD(4) * t643 + t636 * t684 - t687 * t731) + ((-t643 * t676 + t654 * t677) * r_i_i_C(1) + (-t643 * t677 - t654 * t676) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t686 + t688 * t767) * t751, t622 * t695 + t623 * t703 + t647 * t692 + t648 * t713, t625 * qJD(5) - t607 * t717 + t608 * t766 - t712 * t720, t607 (-t608 * t676 + t622 * t677) * r_i_i_C(1) + (-t608 * t677 - t622 * t676) * r_i_i_C(2) + ((-t625 * t677 - t647 * t676) * r_i_i_C(1) + (t625 * t676 - t647 * t677) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
