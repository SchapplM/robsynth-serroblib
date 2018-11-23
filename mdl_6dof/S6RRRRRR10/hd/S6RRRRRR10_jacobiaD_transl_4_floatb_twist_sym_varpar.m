% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR10_jacobiaD_transl_4_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_transl_4_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_transl_4_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10_jacobiaD_transl_4_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiaD_transl_4_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:23
% EndTime: 2018-11-23 11:27:24
% DurationCPUTime: 1.38s
% Computational Cost: add. (3837->213), mult. (4743->315), div. (0->0), fcn. (3702->26), ass. (0->121)
t781 = cos(qJ(2));
t782 = cos(qJ(1));
t777 = sin(qJ(2));
t778 = sin(qJ(1));
t825 = pkin(6) + qJ(2);
t809 = cos(t825) / 0.2e1;
t826 = pkin(6) - qJ(2);
t818 = cos(t826);
t790 = t818 / 0.2e1 + t809;
t784 = t782 * t777 + t778 * t790;
t804 = sin(t825) / 0.2e1;
t815 = sin(t826);
t740 = t804 - t815 / 0.2e1;
t788 = qJD(2) * t740;
t832 = qJD(2) * t778;
t706 = t784 * qJD(1) + t781 * t832 + t782 * t788;
t732 = t790 * qJD(2);
t757 = t777 * t832;
t834 = qJD(1) * t778;
t707 = -t740 * t834 - t757 + (qJD(1) * t781 + t732) * t782;
t789 = t782 * t790;
t715 = t778 * t777 - t789;
t716 = t740 * t782 + t778 * t781;
t824 = pkin(7) - qJ(3);
t814 = sin(t824);
t803 = t814 / 0.2e1;
t823 = pkin(7) + qJ(3);
t813 = sin(t823);
t728 = (t803 - t813 / 0.2e1) * qJD(3);
t807 = cos(t823) / 0.2e1;
t817 = cos(t824);
t743 = t807 - t817 / 0.2e1;
t730 = t743 * qJD(3);
t802 = t813 / 0.2e1;
t737 = t802 + t803;
t744 = t817 / 0.2e1 + t807;
t771 = sin(pkin(6));
t776 = sin(qJ(3));
t780 = cos(qJ(3));
t829 = qJD(3) * t780;
t671 = -t706 * t744 - t707 * t776 - t715 * t728 - t716 * t829 + (-t730 * t782 + t737 * t834) * t771;
t727 = t737 * qJD(3);
t729 = t744 * qJD(3);
t738 = t802 - t814 / 0.2e1;
t830 = qJD(3) * t776;
t675 = t706 * t738 - t707 * t780 + t715 * t729 + t716 * t830 + (t727 * t782 + t743 * t834) * t771;
t835 = t771 * t782;
t686 = -t715 * t738 + t716 * t780 + t743 * t835;
t688 = t715 * t744 + t716 * t776 + t737 * t835;
t773 = cos(pkin(7));
t819 = qJD(1) * t771 * t773;
t770 = sin(pkin(7));
t838 = t706 * t770;
t695 = t778 * t819 + t838;
t711 = -t715 * t770 + t773 * t835;
t821 = pkin(8) + qJ(4);
t811 = sin(t821);
t800 = t811 / 0.2e1;
t822 = pkin(8) - qJ(4);
t812 = sin(t822);
t801 = t812 / 0.2e1;
t735 = t800 + t801;
t723 = t735 * qJD(4);
t805 = cos(t821) / 0.2e1;
t816 = cos(t822);
t742 = t816 / 0.2e1 + t805;
t725 = t742 * qJD(4);
t736 = t800 - t812 / 0.2e1;
t741 = t805 - t816 / 0.2e1;
t779 = cos(qJ(4));
t775 = sin(qJ(4));
t828 = qJD(4) * t775;
t846 = -t671 * t736 + t675 * t779 + t686 * t828 + t688 * t725 + t695 * t741 + t711 * t723;
t724 = (t801 - t811 / 0.2e1) * qJD(4);
t726 = t741 * qJD(4);
t827 = qJD(4) * t779;
t845 = t671 * t742 + t675 * t775 - t686 * t827 - t688 * t724 + t695 * t735 - t711 * t726;
t842 = pkin(12) + r_i_i_C(3);
t745 = t809 - t818 / 0.2e1;
t831 = qJD(2) * t782;
t703 = -qJD(1) * t789 + t777 * t834 + t778 * t788 - t781 * t831;
t704 = t716 * qJD(1) + t778 * t732 + t777 * t831;
t720 = t778 * t740 - t781 * t782;
t833 = qJD(1) * t782;
t841 = -t703 * t738 + t704 * t780 + t784 * t729 - t720 * t830 + (-t727 * t778 + t743 * t833) * t771;
t839 = t703 * t770;
t733 = t745 * qJD(2);
t837 = t733 * t770;
t836 = t771 * t778;
t820 = pkin(11) * t773 + pkin(10);
t799 = -r_i_i_C(1) * t723 - r_i_i_C(2) * t726;
t739 = t804 + t815 / 0.2e1;
t774 = cos(pkin(6));
t701 = t738 * t739 - t743 * t774 - t745 * t780;
t769 = sin(pkin(8));
t792 = t736 * r_i_i_C(1) + t742 * r_i_i_C(2) - t769 * t842;
t692 = t720 * t780 + t738 * t784 + t743 * t836;
t772 = cos(pkin(8));
t791 = t741 * r_i_i_C(1) - t735 * r_i_i_C(2) - t772 * t842 - pkin(11);
t731 = t739 * qJD(2);
t681 = t727 * t774 + t729 * t739 + t731 * t780 + t733 * t738 + t745 * t830;
t668 = t720 * t829 + t703 * t744 + t704 * t776 - t784 * t728 + (t730 * t778 + t737 * t833) * t771;
t690 = t720 * t776 + t737 * t836 - t744 * t784;
t693 = t782 * t819 - t839;
t713 = t770 * t784 + t773 * t836;
t783 = -t668 * t736 - t690 * t725 - t692 * t828 + t693 * t741 - t713 * t723 + t779 * t841;
t714 = -t739 * t770 + t773 * t774;
t710 = t738 * t745 + t739 * t780;
t709 = -t739 * t776 + t744 * t745;
t708 = t720 * qJD(1) - t782 * t732 + t757;
t700 = t737 * t774 + t739 * t744 + t745 * t776;
t699 = t720 * t738 - t780 * t784;
t698 = t720 * t744 + t776 * t784;
t697 = -t715 * t780 - t716 * t738;
t696 = t715 * t776 - t716 * t744;
t684 = t729 * t745 - t731 * t738 + t733 * t780 - t739 * t830;
t680 = t728 * t739 + t730 * t774 - t731 * t776 + t733 * t744 + t745 * t829;
t679 = -t706 * t780 + t708 * t738 + t715 * t830 - t716 * t729;
t677 = t703 * t780 + t704 * t738 + t720 * t729 + t784 * t830;
t667 = t668 * t742 + t690 * t724 + t692 * t827 + t693 * t735 + t713 * t726 + t775 * t841;
t1 = [t846 * r_i_i_C(1) - t845 * r_i_i_C(2) + t675 * pkin(3) - t707 * pkin(2) - pkin(11) * t838 + (-t782 * pkin(1) - t820 * t836) * qJD(1) + t842 * (t671 * t769 - t695 * t772) (t677 * t779 + t698 * t725 - t699 * t828) * r_i_i_C(1) + (-t677 * t775 + t698 * t724 - t699 * t827) * r_i_i_C(2) + t677 * pkin(3) + t703 * pkin(2) + t792 * (-t703 * t776 + t704 * t744 + t720 * t728 + t784 * t829) + (t791 * t704 + t799 * t720) * t770 (t668 * t779 - t690 * t828 + t692 * t725) * r_i_i_C(1) + (-t668 * t775 - t690 * t827 + t692 * t724) * r_i_i_C(2) + t668 * pkin(3) + t792 * t841, t667 * r_i_i_C(1) + t783 * r_i_i_C(2), 0, 0; -t783 * r_i_i_C(1) + t667 * r_i_i_C(2) - t841 * pkin(3) - t704 * pkin(2) - pkin(11) * t839 + (-t778 * pkin(1) + t820 * t835) * qJD(1) + t842 * (-t668 * t769 + t693 * t772) (t679 * t779 + t696 * t725 - t697 * t828) * r_i_i_C(1) + (-t679 * t775 + t696 * t724 - t697 * t827) * r_i_i_C(2) + t679 * pkin(3) - t706 * pkin(2) + t792 * (t706 * t776 + t708 * t744 + t715 * t829 - t716 * t728) + (t791 * t708 - t716 * t799) * t770 (t671 * t779 - t686 * t725 + t688 * t828) * r_i_i_C(1) + (-t671 * t775 - t686 * t724 + t688 * t827) * r_i_i_C(2) + t671 * pkin(3) + t792 * t675, t845 * r_i_i_C(1) + t846 * r_i_i_C(2), 0, 0; 0 (t684 * t779 + t709 * t725 - t710 * t828) * r_i_i_C(1) + (-t684 * t775 + t709 * t724 - t710 * t827) * r_i_i_C(2) + t684 * pkin(3) + t733 * pkin(2) + t792 * (t728 * t745 - t731 * t744 - t733 * t776 - t739 * t829) + (-t791 * t731 + t799 * t745) * t770 (t680 * t779 - t700 * t828 - t701 * t725) * r_i_i_C(1) + (-t680 * t775 - t700 * t827 - t701 * t724) * r_i_i_C(2) + t680 * pkin(3) - t792 * t681 (t680 * t742 - t681 * t775 + t700 * t724 - t701 * t827 + t714 * t726 - t735 * t837) * r_i_i_C(1) + (-t680 * t736 - t681 * t779 - t700 * t725 + t701 * t828 - t714 * t723 - t741 * t837) * r_i_i_C(2), 0, 0;];
JaD_transl  = t1;
