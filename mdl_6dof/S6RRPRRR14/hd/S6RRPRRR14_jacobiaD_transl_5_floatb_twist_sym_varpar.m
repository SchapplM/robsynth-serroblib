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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:25
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14_jacobiaD_transl_5_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_5_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_5_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobiaD_transl_5_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_transl_5_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:25:40
% EndTime: 2019-01-03 10:25:42
% DurationCPUTime: 2.70s
% Computational Cost: add. (1958->200), mult. (6206->363), div. (0->0), fcn. (6786->16), ass. (0->140)
t740 = sin(qJ(2));
t810 = cos(pkin(6));
t811 = sin(qJ(1));
t773 = t810 * t811;
t765 = t740 * t773;
t782 = t811 * t740;
t743 = cos(qJ(2));
t744 = cos(qJ(1));
t787 = t744 * t743;
t710 = -qJD(1) * t765 - qJD(2) * t782 + (qJD(2) * t810 + qJD(1)) * t787;
t731 = sin(pkin(14));
t735 = cos(pkin(14));
t778 = t744 * t810;
t724 = t740 * t778 + t743 * t811;
t750 = t744 * t740 + t743 * t773;
t709 = qJD(1) * t750 + qJD(2) * t724;
t733 = sin(pkin(7));
t737 = cos(pkin(7));
t734 = sin(pkin(6));
t784 = t734 * t811;
t774 = qJD(1) * t784;
t752 = -t709 * t737 + t733 * t774;
t676 = t710 * t735 + t731 * t752;
t739 = sin(qJ(4));
t742 = cos(qJ(4));
t675 = -t710 * t731 + t735 * t752;
t804 = t709 * t733;
t699 = t737 * t774 + t804;
t732 = sin(pkin(8));
t736 = cos(pkin(8));
t772 = t675 * t736 + t699 * t732;
t723 = -t743 * t778 + t782;
t792 = t734 * t744;
t760 = t723 * t737 + t733 * t792;
t694 = -t724 * t735 + t731 * t760;
t693 = t724 * t731 + t735 * t760;
t713 = -t723 * t733 + t737 * t792;
t769 = t693 * t736 + t713 * t732;
t823 = -t694 * t739 + t742 * t769;
t650 = t823 * qJD(4) - t676 * t742 - t739 * t772;
t659 = t675 * t732 - t699 * t736;
t738 = sin(qJ(5));
t741 = cos(qJ(5));
t833 = t650 * t738 - t659 * t741;
t832 = t650 * t741 + t659 * t738;
t662 = t694 * t742 + t739 * t769;
t679 = t693 * t732 - t713 * t736;
t831 = -t662 * t738 - t679 * t741;
t830 = t662 * t741 - t679 * t738;
t829 = qJD(4) * t662 - t676 * t739 + t772 * t742;
t781 = qJD(1) * t792;
t751 = t765 - t787;
t707 = qJD(1) * t723 + qJD(2) * t751;
t805 = t707 * t733;
t697 = t737 * t781 - t805;
t708 = qJD(1) * t724 + qJD(2) * t750;
t755 = t707 * t737 + t733 * t781;
t746 = t708 * t731 + t735 * t755;
t813 = t697 * t732 + t746 * t736;
t777 = t810 * t733;
t788 = t737 * t743;
t790 = t735 * t740;
t712 = t734 * t790 + (t734 * t788 + t777) * t731;
t759 = -t731 * t740 + t735 * t788;
t711 = t734 * t759 + t735 * t777;
t793 = t733 * t743;
t722 = -t734 * t793 + t737 * t810;
t767 = t711 * t736 + t722 * t732;
t673 = t712 * t742 + t739 * t767;
t812 = r_i_i_C(3) + pkin(12);
t791 = t735 * t737;
t682 = -t707 * t731 + t708 * t791;
t809 = t682 * t732;
t684 = t709 * t731 - t710 * t791;
t808 = t684 * t732;
t786 = qJD(2) * t734;
t716 = t759 * t786;
t800 = t716 * t732;
t797 = t731 * t737;
t796 = t732 * t733;
t795 = t733 * t736;
t794 = t733 * t740;
t789 = t737 * t740;
t785 = t734 * t794;
t783 = t737 * t811;
t780 = t733 * t786;
t779 = pkin(11) * t736 + qJ(3);
t776 = t743 * t780;
t775 = t740 * t780;
t756 = t733 * t784 - t737 * t750;
t695 = t731 * t751 + t735 * t756;
t715 = t733 * t750 + t734 * t783;
t768 = t695 * t736 + t715 * t732;
t766 = r_i_i_C(1) * t741 - r_i_i_C(2) * t738 + pkin(4);
t764 = -t682 * t736 + t708 * t796;
t763 = t684 * t736 + t710 * t796;
t702 = t723 * t731 - t724 * t791;
t762 = t702 * t736 + t724 * t796;
t704 = t731 * t750 + t751 * t791;
t761 = t704 * t736 - t751 * t796;
t758 = qJD(5) * (-r_i_i_C(1) * t738 - r_i_i_C(2) * t741);
t720 = (-t731 * t743 - t735 * t789) * t734;
t757 = t720 * t736 + t732 * t785;
t721 = (-t731 * t789 + t735 * t743) * t734;
t754 = -t716 * t736 + t732 * t776;
t718 = qJD(2) * t720;
t753 = t718 * t736 + t732 * t775;
t696 = t731 * t756 - t735 * t751;
t748 = -t696 * t739 + t742 * t768;
t664 = t696 * t742 + t739 * t768;
t747 = -t712 * t739 + t742 * t767;
t703 = -t723 * t735 - t724 * t797;
t669 = t703 * t742 + t739 * t762;
t705 = -t735 * t750 + t751 * t797;
t670 = t705 * t742 + t739 * t761;
t686 = t721 * t742 + t739 * t757;
t657 = t697 * t736 - t732 * t746;
t719 = qJD(2) * t721;
t717 = (-t731 * t788 - t790) * t786;
t706 = -t720 * t732 + t736 * t785;
t701 = -t718 * t732 + t736 * t775;
t700 = t736 * t776 + t800;
t690 = -t711 * t732 + t722 * t736;
t688 = -t704 * t732 - t751 * t795;
t687 = -t702 * t732 + t724 * t795;
t685 = -t709 * t735 - t710 * t797;
t683 = t707 * t735 + t708 * t797;
t681 = -t695 * t732 + t715 * t736;
t674 = -t708 * t735 + t731 * t755;
t668 = t710 * t795 - t808;
t667 = -t708 * t795 - t809;
t666 = t717 * t742 + t754 * t739 + (-t721 * t739 + t742 * t757) * qJD(4);
t656 = qJD(4) * t747 + t719 * t742 + t739 * t753;
t654 = t685 * t742 + t763 * t739 + (-t703 * t739 + t742 * t762) * qJD(4);
t652 = t683 * t742 - t764 * t739 + (-t705 * t739 + t742 * t761) * qJD(4);
t646 = t748 * qJD(4) + t674 * t742 + t813 * t739;
t645 = t664 * qJD(4) + t674 * t739 - t813 * t742;
t644 = t646 * t741 + t657 * t738 + (-t664 * t738 + t681 * t741) * qJD(5);
t643 = -t646 * t738 + t657 * t741 + (-t664 * t741 - t681 * t738) * qJD(5);
t1 = [t832 * r_i_i_C(1) - t833 * r_i_i_C(2) + t650 * pkin(4) - t676 * pkin(3) - t710 * pkin(2) - qJ(3) * t804 + t812 * t829 + (t831 * r_i_i_C(1) - t830 * r_i_i_C(2)) * qJD(5) + t713 * qJD(3) + t659 * pkin(11) + (-t744 * pkin(1) + (-pkin(10) * t811 - qJ(3) * t783) * t734) * qJD(1) (t652 * t741 + t667 * t738) * r_i_i_C(1) + (-t652 * t738 + t667 * t741) * r_i_i_C(2) + t652 * pkin(4) + t683 * pkin(3) - pkin(11) * t809 + t707 * pkin(2) + t812 * (qJD(4) * t670 + t683 * t739 + t742 * t764) + ((-t670 * t738 + t688 * t741) * r_i_i_C(1) + (-t670 * t741 - t688 * t738) * r_i_i_C(2)) * qJD(5) + (-qJD(3) * t751 - t708 * t779) * t733, t697, -t766 * t645 + t646 * t812 + t748 * t758, r_i_i_C(1) * t643 - r_i_i_C(2) * t644, 0; -qJ(3) * t805 - t708 * pkin(2) + t674 * pkin(3) + t646 * pkin(4) + t644 * r_i_i_C(1) + t643 * r_i_i_C(2) + t812 * t645 + t715 * qJD(3) + (-t811 * pkin(1) + (qJ(3) * t737 + pkin(10)) * t792) * qJD(1) + t657 * pkin(11) (t654 * t741 + t668 * t738) * r_i_i_C(1) + (-t654 * t738 + t668 * t741) * r_i_i_C(2) + t654 * pkin(4) + t685 * pkin(3) - pkin(11) * t808 - t709 * pkin(2) + t812 * (qJD(4) * t669 + t685 * t739 - t742 * t763) + ((-t669 * t738 + t687 * t741) * r_i_i_C(1) + (-t669 * t741 - t687 * t738) * r_i_i_C(2)) * qJD(5) + (t724 * qJD(3) + t710 * t779) * t733, t699, -t650 * t812 - t823 * t758 + t766 * t829, t833 * r_i_i_C(1) + t832 * r_i_i_C(2) + (t830 * r_i_i_C(1) + t831 * r_i_i_C(2)) * qJD(5), 0; 0 (t666 * t741 + t700 * t738) * r_i_i_C(1) + (-t666 * t738 + t700 * t741) * r_i_i_C(2) + t666 * pkin(4) + t717 * pkin(3) + pkin(11) * t800 + t812 * (qJD(4) * t686 + t717 * t739 - t742 * t754) + ((-t686 * t738 + t706 * t741) * r_i_i_C(1) + (-t686 * t741 - t706 * t738) * r_i_i_C(2)) * qJD(5) + (qJD(3) * t794 + (-pkin(2) * t740 + t779 * t793) * qJD(2)) * t734, t775, t812 * t656 + t747 * t758 + t766 * (-t673 * qJD(4) - t719 * t739 + t753 * t742) (-t656 * t738 + t701 * t741) * r_i_i_C(1) + (-t656 * t741 - t701 * t738) * r_i_i_C(2) + ((-t673 * t741 - t690 * t738) * r_i_i_C(1) + (t673 * t738 - t690 * t741) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
