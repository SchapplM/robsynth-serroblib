% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:14
% EndTime: 2019-02-26 22:52:17
% DurationCPUTime: 2.59s
% Computational Cost: add. (2727->240), mult. (7763->391), div. (0->0), fcn. (8472->16), ass. (0->150)
t721 = sin(qJ(2));
t724 = cos(qJ(2));
t725 = cos(qJ(1));
t816 = cos(pkin(6));
t779 = t725 * t816;
t820 = sin(qJ(1));
t694 = t721 * t779 + t820 * t724;
t756 = t816 * t820;
t738 = t725 * t721 + t724 * t756;
t681 = t738 * qJD(1) + t694 * qJD(2);
t749 = t721 * t756;
t785 = t820 * t721;
t801 = t725 * t724;
t682 = -qJD(1) * t749 - qJD(2) * t785 + (qJD(2) * t816 + qJD(1)) * t801;
t693 = -t724 * t779 + t785;
t715 = sin(pkin(7));
t717 = cos(pkin(7));
t720 = sin(qJ(3));
t716 = sin(pkin(6));
t805 = t716 * t725;
t792 = t715 * t805;
t821 = cos(qJ(3));
t753 = t821 * t792;
t781 = qJD(3) * t821;
t758 = t717 * t781;
t782 = qJD(1) * t820;
t759 = t716 * t782;
t640 = (-qJD(3) * t694 - t681 * t717 + t715 * t759) * t720 - qJD(3) * t753 + t682 * t821 - t693 * t758;
t814 = t681 * t715;
t670 = t717 * t759 + t814;
t719 = sin(qJ(4));
t723 = cos(qJ(4));
t790 = t694 * t821;
t666 = (t693 * t717 + t792) * t720 - t790;
t685 = -t693 * t715 + t717 * t805;
t752 = t666 * t719 - t685 * t723;
t632 = -t752 * qJD(4) - t640 * t723 - t670 * t719;
t789 = t717 * t821;
t663 = t693 * t789 + t694 * t720 + t753;
t713 = qJD(5) + qJD(6);
t836 = -t663 * t713 + t632;
t809 = t715 * t716;
t762 = t820 * t809;
t746 = t821 * t762;
t764 = t720 * t792;
t804 = t717 * t720;
t641 = qJD(1) * t746 - t681 * t789 - t682 * t720 + (t693 * t804 + t764 - t790) * qJD(3);
t653 = t666 * t723 + t685 * t719;
t835 = -t653 * t713 + t641;
t834 = t653 * qJD(4) - t640 * t719 + t670 * t723;
t817 = r_i_i_C(3) + pkin(13) + pkin(12);
t718 = sin(qJ(5));
t819 = pkin(5) * t718;
t825 = t819 + pkin(11);
t739 = t749 - t801;
t824 = t720 * t739 + t746;
t668 = -t739 * t821 + (-t717 * t738 + t762) * t720;
t788 = t717 * t820;
t687 = t715 * t738 + t716 * t788;
t823 = -t668 * t719 + t687 * t723;
t714 = qJ(5) + qJ(6);
t711 = sin(t714);
t712 = cos(t714);
t794 = qJD(5) * t819;
t822 = (t711 * r_i_i_C(1) + t712 * r_i_i_C(2)) * t713 + t794;
t818 = pkin(10) * t715;
t811 = t711 * t713;
t810 = t712 * t713;
t808 = t715 * t719;
t807 = t715 * t723;
t806 = t716 * t724;
t803 = t720 * t721;
t802 = t720 * t724;
t680 = t694 * qJD(1) + t738 * qJD(2);
t735 = t739 * qJD(2);
t731 = t693 * qJD(1) + t735;
t729 = t731 * t821;
t637 = -qJD(1) * t753 + t668 * qJD(3) - t680 * t720 - t717 * t729;
t655 = t668 * t723 + t687 * t719;
t773 = -t655 * t713 + t637;
t730 = t731 * t720;
t638 = qJD(1) * t764 + qJD(3) * t824 - t680 * t821 + t717 * t730 - t738 * t758;
t727 = t685 * qJD(1) - t715 * t735;
t628 = qJD(4) * t823 + t638 * t723 + t727 * t719;
t667 = t738 * t789 - t824;
t778 = t667 * t713 + t628;
t623 = -t778 * t711 + t773 * t712;
t624 = t773 * t711 + t778 * t712;
t800 = t623 * r_i_i_C(1) - t624 * r_i_i_C(2);
t799 = (t711 * t836 - t712 * t835) * r_i_i_C(1) + (t711 * t835 + t712 * t836) * r_i_i_C(2);
t787 = t821 * t721;
t740 = t717 * t787 + t802;
t741 = t717 * t802 + t787;
t780 = t715 * t816;
t757 = t720 * t780;
t658 = qJD(3) * t757 + (t740 * qJD(2) + t741 * qJD(3)) * t716;
t684 = t741 * t716 + t757;
t692 = -t715 * t806 + t816 * t717;
t661 = t684 * t723 + t692 * t719;
t766 = t661 * t713 - t658;
t786 = t821 * t724;
t742 = -t717 * t803 + t786;
t750 = t821 * t780;
t763 = t717 * t786;
t659 = qJD(3) * t750 + ((t763 - t803) * qJD(3) + t742 * qJD(2)) * t716;
t751 = -t684 * t719 + t692 * t723;
t797 = qJD(2) * t716;
t783 = t715 * t797;
t761 = t721 * t783;
t644 = t751 * qJD(4) + t659 * t723 + t719 * t761;
t754 = t716 * t763;
t791 = t716 * t803;
t683 = -t750 - t754 + t791;
t770 = -t683 * t713 - t644;
t798 = (t770 * t711 - t766 * t712) * r_i_i_C(1) + (t766 * t711 + t770 * t712) * r_i_i_C(2);
t722 = cos(qJ(5));
t796 = qJD(5) * t722;
t795 = pkin(5) * t796;
t793 = t721 * t809;
t743 = t720 * t738 + t739 * t789;
t646 = t743 * qJD(3) + t680 * t804 + t729;
t678 = -t738 * t821 + t739 * t804;
t634 = -t680 * t808 + t646 * t723 + (-t678 * t719 - t739 * t807) * qJD(4);
t775 = -t713 * t743 + t634;
t744 = t693 * t720 - t694 * t789;
t648 = t744 * qJD(3) - t681 * t821 - t682 * t804;
t676 = -t693 * t821 - t694 * t804;
t636 = t682 * t808 + t648 * t723 + (-t676 * t719 + t694 * t807) * qJD(4);
t774 = -t713 * t744 + t636;
t645 = t678 * qJD(3) - t680 * t789 + t730;
t657 = t678 * t723 - t739 * t808;
t769 = -t657 * t713 + t645;
t647 = t676 * qJD(3) - t681 * t720 + t682 * t789;
t656 = t676 * t723 + t694 * t808;
t768 = -t656 * t713 + t647;
t672 = (-t741 * qJD(2) - t740 * qJD(3)) * t716;
t691 = t742 * t716;
t760 = t724 * t783;
t650 = t719 * t760 + t672 * t723 + (-t691 * t719 + t723 * t793) * qJD(4);
t690 = t740 * t716;
t767 = t690 * t713 + t650;
t671 = -qJD(2) * t754 - t781 * t806 + (qJD(3) * t717 + qJD(2)) * t791;
t679 = t691 * t723 + t719 * t793;
t765 = -t679 * t713 - t671;
t710 = t722 * pkin(5) + pkin(4);
t748 = t712 * r_i_i_C(1) - t711 * r_i_i_C(2) + t710;
t732 = -t817 * t719 - t748 * t723 - pkin(3);
t728 = t822 * t723 + (t748 * t719 - t817 * t723) * qJD(4);
t627 = t655 * qJD(4) + t638 * t719 - t727 * t723;
t1 = [-pkin(10) * t814 - t682 * pkin(2) - t640 * pkin(3) + t641 * pkin(11) + t632 * t710 + (r_i_i_C(1) * t836 + r_i_i_C(2) * t835) * t712 + (r_i_i_C(1) * t835 - r_i_i_C(2) * t836) * t711 + t817 * t834 + (-t725 * pkin(1) + (-t820 * pkin(9) - pkin(10) * t788) * t716) * qJD(1) + (t641 * t718 + (-t653 * t718 - t663 * t722) * qJD(5)) * pkin(5) (t769 * t711 + t775 * t712) * r_i_i_C(1) + (-t775 * t711 + t769 * t712) * r_i_i_C(2) + t634 * t710 - t657 * t794 - t743 * t795 + t646 * pkin(3) + t731 * pkin(2) - t680 * t818 + t825 * t645 + t817 * (t657 * qJD(4) + t646 * t719 + t680 * t807) (t638 * t711 + t668 * t810) * r_i_i_C(1) + (t638 * t712 - t668 * t811) * r_i_i_C(2) + t638 * pkin(11) + (t638 * t718 + t668 * t796) * pkin(5) + t732 * t637 + t728 * t667, -t748 * t627 + t817 * t628 - t822 * t823 (-t628 * t718 + t637 * t722 + (-t655 * t722 - t667 * t718) * qJD(5)) * pkin(5) + t800, t800; -pkin(1) * t782 - t680 * pkin(2) + t638 * pkin(3) + t624 * r_i_i_C(1) + t623 * r_i_i_C(2) + t628 * t710 - t655 * t794 + t667 * t795 - t731 * t818 + (pkin(10) * t717 + pkin(9)) * qJD(1) * t805 + t825 * t637 + t817 * t627, t682 * t818 - t681 * pkin(2) + t648 * pkin(3) + t647 * pkin(11) + t636 * t710 + (t774 * r_i_i_C(1) + t768 * r_i_i_C(2)) * t712 + (t768 * r_i_i_C(1) - t774 * r_i_i_C(2)) * t711 + t817 * (t656 * qJD(4) + t648 * t719 - t682 * t807) + (t647 * t718 + (-t656 * t718 - t722 * t744) * qJD(5)) * pkin(5) (t640 * t711 - t666 * t810) * r_i_i_C(1) + (t640 * t712 + t666 * t811) * r_i_i_C(2) + t640 * pkin(11) + (t640 * t718 - t666 * t796) * pkin(5) - t732 * t641 + t728 * t663, -t632 * t817 + t748 * t834 - t822 * t752 (t632 * t718 - t641 * t722 + (t653 * t722 - t663 * t718) * qJD(5)) * pkin(5) + t799, t799; 0, t672 * pkin(3) - t671 * pkin(11) + t650 * t710 + (t767 * r_i_i_C(1) + t765 * r_i_i_C(2)) * t712 + (t765 * r_i_i_C(1) - t767 * r_i_i_C(2)) * t711 + t817 * (t679 * qJD(4) + t672 * t719 - t723 * t760) + (-pkin(2) * t721 + t724 * t818) * t797 + (-t671 * t718 + (-t679 * t718 + t690 * t722) * qJD(5)) * pkin(5) (t659 * t711 + t684 * t810) * r_i_i_C(1) + (t659 * t712 - t684 * t811) * r_i_i_C(2) + t659 * pkin(11) + (t659 * t718 + t684 * t796) * pkin(5) + t732 * t658 + t728 * t683, t817 * t644 - t822 * t751 + t748 * (-qJD(4) * t661 - t659 * t719 + t723 * t761) (-t644 * t718 + t658 * t722 + (-t661 * t722 - t683 * t718) * qJD(5)) * pkin(5) + t798, t798;];
JaD_transl  = t1;
