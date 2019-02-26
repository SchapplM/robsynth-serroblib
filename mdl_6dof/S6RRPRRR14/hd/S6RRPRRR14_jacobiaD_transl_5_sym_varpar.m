% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_transl_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:55:52
% EndTime: 2019-02-26 22:55:56
% DurationCPUTime: 2.84s
% Computational Cost: add. (1958->200), mult. (6206->363), div. (0->0), fcn. (6786->16), ass. (0->140)
t739 = sin(qJ(2));
t809 = cos(pkin(6));
t810 = sin(qJ(1));
t772 = t809 * t810;
t764 = t739 * t772;
t781 = t810 * t739;
t742 = cos(qJ(2));
t743 = cos(qJ(1));
t786 = t743 * t742;
t709 = -qJD(1) * t764 - qJD(2) * t781 + (qJD(2) * t809 + qJD(1)) * t786;
t730 = sin(pkin(14));
t734 = cos(pkin(14));
t777 = t743 * t809;
t723 = t739 * t777 + t810 * t742;
t749 = t743 * t739 + t742 * t772;
t708 = t749 * qJD(1) + t723 * qJD(2);
t732 = sin(pkin(7));
t736 = cos(pkin(7));
t733 = sin(pkin(6));
t783 = t733 * t810;
t773 = qJD(1) * t783;
t751 = -t708 * t736 + t732 * t773;
t675 = t709 * t734 + t751 * t730;
t738 = sin(qJ(4));
t741 = cos(qJ(4));
t674 = -t709 * t730 + t751 * t734;
t803 = t708 * t732;
t698 = t736 * t773 + t803;
t731 = sin(pkin(8));
t735 = cos(pkin(8));
t771 = t674 * t735 + t698 * t731;
t722 = -t742 * t777 + t781;
t791 = t733 * t743;
t759 = t722 * t736 + t732 * t791;
t693 = -t723 * t734 + t759 * t730;
t692 = t723 * t730 + t759 * t734;
t712 = -t722 * t732 + t736 * t791;
t768 = t692 * t735 + t712 * t731;
t822 = -t693 * t738 + t768 * t741;
t649 = t822 * qJD(4) - t675 * t741 - t771 * t738;
t658 = t674 * t731 - t698 * t735;
t737 = sin(qJ(5));
t740 = cos(qJ(5));
t832 = t649 * t737 - t658 * t740;
t831 = t649 * t740 + t658 * t737;
t661 = t693 * t741 + t768 * t738;
t678 = t692 * t731 - t712 * t735;
t830 = -t661 * t737 - t678 * t740;
t829 = t661 * t740 - t678 * t737;
t828 = t661 * qJD(4) - t675 * t738 + t771 * t741;
t780 = qJD(1) * t791;
t750 = t764 - t786;
t706 = t722 * qJD(1) + t750 * qJD(2);
t804 = t706 * t732;
t696 = t736 * t780 - t804;
t707 = t723 * qJD(1) + t749 * qJD(2);
t754 = t706 * t736 + t732 * t780;
t745 = t707 * t730 + t754 * t734;
t812 = t696 * t731 + t745 * t735;
t776 = t809 * t732;
t787 = t736 * t742;
t789 = t734 * t739;
t711 = t733 * t789 + (t733 * t787 + t776) * t730;
t758 = -t730 * t739 + t734 * t787;
t710 = t758 * t733 + t734 * t776;
t792 = t732 * t742;
t721 = -t733 * t792 + t809 * t736;
t766 = t710 * t735 + t721 * t731;
t672 = t711 * t741 + t766 * t738;
t811 = r_i_i_C(3) + pkin(12);
t790 = t734 * t736;
t681 = -t706 * t730 + t707 * t790;
t808 = t681 * t731;
t683 = t708 * t730 - t709 * t790;
t807 = t683 * t731;
t785 = qJD(2) * t733;
t715 = t758 * t785;
t799 = t715 * t731;
t796 = t730 * t736;
t795 = t731 * t732;
t794 = t732 * t735;
t793 = t732 * t739;
t788 = t736 * t739;
t784 = t733 * t793;
t782 = t736 * t810;
t779 = t732 * t785;
t778 = pkin(11) * t735 + qJ(3);
t775 = t742 * t779;
t774 = t739 * t779;
t755 = t732 * t783 - t736 * t749;
t694 = t730 * t750 + t755 * t734;
t714 = t732 * t749 + t733 * t782;
t767 = t694 * t735 + t714 * t731;
t765 = r_i_i_C(1) * t740 - r_i_i_C(2) * t737 + pkin(4);
t763 = -t681 * t735 + t707 * t795;
t762 = t683 * t735 + t709 * t795;
t701 = t722 * t730 - t723 * t790;
t761 = t701 * t735 + t723 * t795;
t703 = t730 * t749 + t750 * t790;
t760 = t703 * t735 - t750 * t795;
t757 = qJD(5) * (-r_i_i_C(1) * t737 - r_i_i_C(2) * t740);
t719 = (-t730 * t742 - t734 * t788) * t733;
t756 = t719 * t735 + t731 * t784;
t720 = (-t730 * t788 + t734 * t742) * t733;
t753 = -t715 * t735 + t731 * t775;
t717 = qJD(2) * t719;
t752 = t717 * t735 + t731 * t774;
t695 = t755 * t730 - t734 * t750;
t747 = -t695 * t738 + t767 * t741;
t663 = t695 * t741 + t767 * t738;
t746 = -t711 * t738 + t766 * t741;
t702 = -t722 * t734 - t723 * t796;
t668 = t702 * t741 + t761 * t738;
t704 = -t734 * t749 + t750 * t796;
t669 = t704 * t741 + t760 * t738;
t685 = t720 * t741 + t756 * t738;
t656 = t696 * t735 - t745 * t731;
t718 = qJD(2) * t720;
t716 = (-t730 * t787 - t789) * t785;
t705 = -t719 * t731 + t735 * t784;
t700 = -t717 * t731 + t735 * t774;
t699 = t735 * t775 + t799;
t689 = -t710 * t731 + t721 * t735;
t687 = -t703 * t731 - t750 * t794;
t686 = -t701 * t731 + t723 * t794;
t684 = -t708 * t734 - t709 * t796;
t682 = t706 * t734 + t707 * t796;
t680 = -t694 * t731 + t714 * t735;
t673 = -t707 * t734 + t754 * t730;
t667 = t709 * t794 - t807;
t666 = -t707 * t794 - t808;
t665 = t716 * t741 + t753 * t738 + (-t720 * t738 + t756 * t741) * qJD(4);
t655 = t746 * qJD(4) + t718 * t741 + t752 * t738;
t653 = t684 * t741 + t762 * t738 + (-t702 * t738 + t761 * t741) * qJD(4);
t651 = t682 * t741 - t763 * t738 + (-t704 * t738 + t760 * t741) * qJD(4);
t645 = t747 * qJD(4) + t673 * t741 + t812 * t738;
t644 = t663 * qJD(4) + t673 * t738 - t812 * t741;
t643 = t645 * t740 + t656 * t737 + (-t663 * t737 + t680 * t740) * qJD(5);
t642 = -t645 * t737 + t656 * t740 + (-t663 * t740 - t680 * t737) * qJD(5);
t1 = [t831 * r_i_i_C(1) - t832 * r_i_i_C(2) + t649 * pkin(4) - t675 * pkin(3) - t709 * pkin(2) - qJ(3) * t803 + t811 * t828 + (r_i_i_C(1) * t830 - t829 * r_i_i_C(2)) * qJD(5) + t712 * qJD(3) + t658 * pkin(11) + (-t743 * pkin(1) + (-t810 * pkin(10) - qJ(3) * t782) * t733) * qJD(1) (t651 * t740 + t666 * t737) * r_i_i_C(1) + (-t651 * t737 + t666 * t740) * r_i_i_C(2) + t651 * pkin(4) + t682 * pkin(3) - pkin(11) * t808 + t706 * pkin(2) + t811 * (t669 * qJD(4) + t682 * t738 + t763 * t741) + ((-t669 * t737 + t687 * t740) * r_i_i_C(1) + (-t669 * t740 - t687 * t737) * r_i_i_C(2)) * qJD(5) + (-qJD(3) * t750 - t778 * t707) * t732, t696, -t765 * t644 + t811 * t645 + t747 * t757, r_i_i_C(1) * t642 - r_i_i_C(2) * t643, 0; -qJ(3) * t804 - t707 * pkin(2) + t673 * pkin(3) + t645 * pkin(4) + t643 * r_i_i_C(1) + t642 * r_i_i_C(2) + t811 * t644 + t714 * qJD(3) + (-t810 * pkin(1) + (qJ(3) * t736 + pkin(10)) * t791) * qJD(1) + t656 * pkin(11) (t653 * t740 + t667 * t737) * r_i_i_C(1) + (-t653 * t737 + t667 * t740) * r_i_i_C(2) + t653 * pkin(4) + t684 * pkin(3) - pkin(11) * t807 - t708 * pkin(2) + t811 * (t668 * qJD(4) + t684 * t738 - t762 * t741) + ((-t668 * t737 + t686 * t740) * r_i_i_C(1) + (-t668 * t740 - t686 * t737) * r_i_i_C(2)) * qJD(5) + (t723 * qJD(3) + t778 * t709) * t732, t698, -t649 * t811 - t822 * t757 + t765 * t828, t832 * r_i_i_C(1) + t831 * r_i_i_C(2) + (t829 * r_i_i_C(1) + r_i_i_C(2) * t830) * qJD(5), 0; 0 (t665 * t740 + t699 * t737) * r_i_i_C(1) + (-t665 * t737 + t699 * t740) * r_i_i_C(2) + t665 * pkin(4) + t716 * pkin(3) + pkin(11) * t799 + t811 * (t685 * qJD(4) + t716 * t738 - t753 * t741) + ((-t685 * t737 + t705 * t740) * r_i_i_C(1) + (-t685 * t740 - t705 * t737) * r_i_i_C(2)) * qJD(5) + (qJD(3) * t793 + (-pkin(2) * t739 + t778 * t792) * qJD(2)) * t733, t774, t811 * t655 + t746 * t757 + t765 * (-t672 * qJD(4) - t718 * t738 + t752 * t741) (-t655 * t737 + t700 * t740) * r_i_i_C(1) + (-t655 * t740 - t700 * t737) * r_i_i_C(2) + ((-t672 * t740 - t689 * t737) * r_i_i_C(1) + (t672 * t737 - t689 * t740) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
