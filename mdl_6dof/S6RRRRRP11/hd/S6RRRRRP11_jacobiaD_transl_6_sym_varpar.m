% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP11_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP11_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:45
% EndTime: 2019-02-26 22:45:48
% DurationCPUTime: 2.49s
% Computational Cost: add. (2401->227), mult. (7333->372), div. (0->0), fcn. (8019->14), ass. (0->137)
t682 = cos(qJ(2));
t683 = cos(qJ(1));
t774 = cos(pkin(6));
t741 = t683 * t774;
t679 = sin(qJ(2));
t778 = sin(qJ(1));
t745 = t778 * t679;
t655 = -t682 * t741 + t745;
t674 = cos(pkin(7));
t678 = sin(qJ(3));
t656 = t679 * t741 + t778 * t682;
t779 = cos(qJ(3));
t750 = t656 * t779;
t672 = sin(pkin(7));
t673 = sin(pkin(6));
t763 = t673 * t683;
t752 = t672 * t763;
t629 = (t655 * t674 + t752) * t678 - t750;
t647 = -t655 * t672 + t674 * t763;
t677 = sin(qJ(4));
t681 = cos(qJ(4));
t615 = t629 * t681 + t647 * t677;
t729 = t779 * t752;
t749 = t674 * t779;
t783 = -t655 * t749 - t729;
t626 = t656 * t678 - t783;
t676 = sin(qJ(5));
t680 = cos(qJ(5));
t793 = -t615 * t676 - t626 * t680;
t726 = t615 * t680 - t626 * t676;
t731 = t774 * t778;
t696 = t683 * t679 + t682 * t731;
t643 = t696 * qJD(1) + t656 * qJD(2);
t717 = t679 * t731;
t759 = t683 * t682;
t644 = -qJD(1) * t717 - qJD(2) * t745 + (qJD(2) * t774 + qJD(1)) * t759;
t742 = t778 * qJD(1);
t733 = t673 * t742;
t602 = (-qJD(3) * t656 - t643 * t674 + t672 * t733) * t678 + t644 * t779 + t783 * qJD(3);
t770 = t643 * t672;
t698 = t674 * t733 + t770;
t792 = t615 * qJD(4) - t602 * t677 + t698 * t681;
t722 = t629 * t677 - t647 * t681;
t592 = t722 * qJD(4) + t602 * t681 + t698 * t677;
t780 = pkin(5) + r_i_i_C(1);
t714 = t680 * r_i_i_C(2) + t780 * t676;
t706 = pkin(11) + t714;
t699 = qJD(5) * t714;
t775 = r_i_i_C(3) + qJ(6) + pkin(12);
t697 = t717 - t759;
t767 = t672 * t673;
t736 = t778 * t767;
t631 = -t697 * t779 + (-t696 * t674 + t736) * t678;
t748 = t674 * t778;
t649 = t696 * t672 + t673 * t748;
t616 = -t631 * t677 + t649 * t681;
t691 = t696 * t779;
t710 = t779 * t736;
t782 = -t674 * t691 + t678 * t697 + t710;
t739 = t678 * t752;
t762 = t674 * t678;
t781 = qJD(1) * t710 - t643 * t749 - t644 * t678 + (t655 * t762 + t739 - t750) * qJD(3);
t777 = pkin(10) * t672;
t642 = t656 * qJD(1) + t696 * qJD(2);
t694 = t697 * qJD(2);
t688 = t655 * qJD(1) + t694;
t686 = t688 * t779;
t599 = -qJD(1) * t729 + t631 * qJD(3) - t642 * t678 - t674 * t686;
t773 = t599 * t676;
t640 = t697 * t762 - t691;
t687 = t688 * t678;
t607 = t640 * qJD(3) - t642 * t749 + t687;
t772 = t607 * t676;
t766 = t672 * t677;
t765 = t672 * t681;
t764 = t672 * t682;
t761 = t678 * t679;
t760 = t678 * t682;
t758 = qJD(2) * t673;
t757 = qJD(5) * t676;
t756 = qJD(5) * t680;
t755 = pkin(5) * t757;
t754 = pkin(5) * t756;
t753 = t679 * t767;
t751 = t673 * t761;
t747 = t779 * t679;
t746 = t779 * t682;
t743 = t672 * t758;
t740 = t774 * t672;
t737 = t674 * t746;
t735 = t679 * t743;
t734 = t682 * t743;
t732 = t678 * t740;
t730 = t673 * t737;
t728 = -t592 * t676 - t680 * t781;
t702 = -t674 * t761 + t746;
t718 = t779 * t740;
t621 = qJD(3) * t718 + ((t737 - t761) * qJD(3) + t702 * qJD(2)) * t673;
t701 = t674 * t760 + t747;
t646 = t701 * t673 + t732;
t654 = -t673 * t764 + t774 * t674;
t719 = -t646 * t677 + t654 * t681;
t606 = t719 * qJD(4) + t621 * t681 + t677 * t735;
t700 = t674 * t747 + t760;
t620 = qJD(3) * t732 + (t700 * qJD(2) + t701 * qJD(3)) * t673;
t727 = -t606 * t676 + t620 * t680;
t623 = t646 * t681 + t654 * t677;
t645 = -t718 - t730 + t751;
t723 = -t623 * t680 - t645 * t676;
t617 = t631 * t681 + t649 * t677;
t671 = t680 * pkin(5) + pkin(4);
t716 = t680 * r_i_i_C(1) - t676 * r_i_i_C(2) + t671;
t638 = -t655 * t779 - t656 * t762;
t713 = -t638 * t677 + t656 * t765;
t618 = t638 * t681 + t656 * t766;
t712 = -t640 * t677 - t697 * t765;
t619 = t640 * t681 - t697 * t766;
t653 = t702 * t673;
t705 = -t653 * t677 + t681 * t753;
t641 = t653 * t681 + t677 * t753;
t703 = t655 * t678 - t656 * t749;
t690 = -t775 * t677 - t716 * t681 - pkin(3);
t600 = qJD(1) * t739 + qJD(3) * t782 - t642 * t779 + t674 * t687;
t684 = t647 * qJD(1) - t672 * t694;
t590 = qJD(4) * t616 + t600 * t681 + t684 * t677;
t587 = -t590 * t676 + t599 * t680 + (-t617 * t680 + t676 * t782) * qJD(5);
t639 = -t696 * t678 - t697 * t749;
t685 = -t677 * qJD(6) + t681 * t699 + (t716 * t677 - t775 * t681) * qJD(4);
t652 = t700 * t673;
t634 = (-t701 * qJD(2) - t700 * qJD(3)) * t673;
t610 = t703 * qJD(3) - t643 * t779 - t644 * t762;
t608 = -t639 * qJD(3) + t642 * t762 + t686;
t605 = t623 * qJD(4) + t621 * t677 - t681 * t735;
t596 = t712 * qJD(4) + t608 * t681 - t642 * t766;
t589 = t617 * qJD(4) + t600 * t677 - t684 * t681;
t588 = t590 * t680 + t773 + (-t617 * t676 - t680 * t782) * qJD(5);
t1 = [t722 * qJD(6) - t602 * pkin(3) - t644 * pkin(2) - pkin(10) * t770 - t716 * t592 + t706 * t781 + t775 * t792 + (-t683 * pkin(1) + (-t778 * pkin(9) - pkin(10) * t748) * t673) * qJD(1) + (-t726 * r_i_i_C(2) + t780 * t793) * qJD(5) (t596 * t680 + t772 + (-t619 * t676 + t639 * t680) * qJD(5)) * r_i_i_C(1) + (-t596 * t676 + t607 * t680 + (-t619 * t680 - t639 * t676) * qJD(5)) * r_i_i_C(2) + t596 * t671 - t619 * t755 - t712 * qJD(6) + pkin(5) * t772 + t639 * t754 + t608 * pkin(3) + t607 * pkin(11) + t688 * pkin(2) - t642 * t777 + t775 * (t619 * qJD(4) + t608 * t677 + t642 * t765) (t600 * t680 - t631 * t757) * r_i_i_C(2) + t600 * pkin(11) + t690 * t599 - t685 * t782 + t780 * (t600 * t676 + t631 * t756) t617 * qJD(6) - t716 * t589 + t775 * t590 - t616 * t699, -t588 * r_i_i_C(2) + t587 * t780, t589; -pkin(1) * t742 - t642 * pkin(2) + t600 * pkin(3) + pkin(5) * t773 + t599 * pkin(11) + t588 * r_i_i_C(1) + t587 * r_i_i_C(2) - t616 * qJD(6) + t590 * t671 - t617 * t755 - t782 * t754 - t688 * t777 + (pkin(10) * t674 + pkin(9)) * qJD(1) * t763 + t775 * t589, -t713 * qJD(6) + t610 * pkin(3) - t643 * pkin(2) + t644 * t777 + t716 * (t713 * qJD(4) + t610 * t681 + t644 * t766) + t706 * (t638 * qJD(3) - t643 * t678 + t644 * t749) + t775 * (t618 * qJD(4) + t610 * t677 - t644 * t765) + ((-t618 * t680 + t676 * t703) * r_i_i_C(2) + t780 * (-t618 * t676 - t680 * t703)) * qJD(5) (t602 * t680 + t629 * t757) * r_i_i_C(2) + t602 * pkin(11) - t690 * t781 + t685 * t626 + t780 * (t602 * t676 - t629 * t756) -qJD(6) * t615 + t775 * t592 - t722 * t699 + t716 * t792, t728 * r_i_i_C(1) + (-t592 * t680 + t676 * t781) * r_i_i_C(2) + (t726 * r_i_i_C(1) + t793 * r_i_i_C(2)) * qJD(5) + (t726 * qJD(5) + t728) * pkin(5), -t792; 0, -t705 * qJD(6) + t634 * pkin(3) + t716 * (t705 * qJD(4) + t634 * t681 + t677 * t734) - t706 * (-qJD(2) * t730 - t673 * qJD(3) * t746 + (qJD(3) * t674 + qJD(2)) * t751) + t775 * (t641 * qJD(4) + t634 * t677 - t681 * t734) + (-pkin(2) * t679 + pkin(10) * t764) * t758 + ((-t641 * t680 - t652 * t676) * r_i_i_C(2) + t780 * (-t641 * t676 + t652 * t680)) * qJD(5) (t621 * t680 - t646 * t757) * r_i_i_C(2) + t621 * pkin(11) + t690 * t620 + t685 * t645 + t780 * (t621 * t676 + t646 * t756) t623 * qJD(6) - t716 * t605 + t775 * t606 - t719 * t699, t727 * r_i_i_C(1) + (-t606 * t680 - t620 * t676) * r_i_i_C(2) + (t723 * r_i_i_C(1) + (t623 * t676 - t645 * t680) * r_i_i_C(2)) * qJD(5) + (t723 * qJD(5) + t727) * pkin(5), t605;];
JaD_transl  = t1;
