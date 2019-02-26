% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR12
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
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR12_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR12_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:37:06
% EndTime: 2019-02-26 22:37:09
% DurationCPUTime: 2.74s
% Computational Cost: add. (2398->252), mult. (6559->414), div. (0->0), fcn. (7118->16), ass. (0->139)
t697 = sin(qJ(2));
t700 = cos(qJ(2));
t701 = cos(qJ(1));
t773 = cos(pkin(6));
t741 = t701 * t773;
t779 = sin(qJ(1));
t670 = t697 * t741 + t779 * t700;
t732 = t773 * t779;
t715 = t701 * t697 + t700 * t732;
t657 = t715 * qJD(1) + t670 * qJD(2);
t724 = t697 * t732;
t746 = t779 * t697;
t761 = t701 * t700;
t658 = -qJD(1) * t724 - qJD(2) * t746 + (qJD(2) * t773 + qJD(1)) * t761;
t690 = sin(pkin(7));
t692 = cos(pkin(7));
t696 = sin(qJ(3));
t691 = sin(pkin(6));
t743 = t779 * qJD(1);
t734 = t691 * t743;
t780 = cos(qJ(3));
t669 = -t700 * t741 + t746;
t765 = t691 * t701;
t751 = t690 * t765;
t729 = t780 * t751;
t750 = t692 * t780;
t785 = -t669 * t750 - t729;
t616 = (-qJD(3) * t670 - t657 * t692 + t690 * t734) * t696 + t658 * t780 + t785 * qJD(3);
t772 = t657 * t690;
t646 = t692 * t734 + t772;
t689 = qJ(4) + pkin(13);
t687 = sin(t689);
t688 = cos(t689);
t661 = -t669 * t690 + t692 * t765;
t679 = t696 * t751;
t764 = t692 * t696;
t793 = t669 * t764 - t670 * t780 + t679;
t784 = t661 * t688 - t687 * t793;
t608 = t784 * qJD(4) - t616 * t688 - t646 * t687;
t766 = t691 * t690;
t737 = t779 * t766;
t722 = t780 * t737;
t617 = qJD(1) * t722 + qJD(3) * t793 - t657 * t750 - t658 * t696;
t694 = sin(qJ(6));
t698 = cos(qJ(6));
t800 = t608 * t694 - t617 * t698;
t799 = t608 * t698 + t617 * t694;
t628 = -t661 * t687 - t688 * t793;
t639 = t670 * t696 - t785;
t798 = t628 * t694 - t639 * t698;
t797 = t628 * t698 + t639 * t694;
t796 = -t628 * qJD(4) - t616 * t687 + t646 * t688;
t716 = t724 - t761;
t644 = -t716 * t780 + (-t715 * t692 + t737) * t696;
t713 = t716 * qJD(2);
t703 = t661 * qJD(1) - t690 * t713;
t792 = -t644 * qJD(4) + t703;
t781 = r_i_i_C(3) + pkin(12);
t783 = t792 * pkin(4);
t723 = qJD(6) * (r_i_i_C(1) * t694 + r_i_i_C(2) * t698);
t726 = t698 * r_i_i_C(1) - t694 * r_i_i_C(2) + pkin(5);
t710 = t715 * t780;
t782 = -t692 * t710 + t696 * t716 + t722;
t778 = t690 * pkin(10);
t695 = sin(qJ(4));
t776 = t695 * pkin(4);
t774 = pkin(4) * qJD(4);
t770 = t687 * t690;
t769 = t688 * t690;
t768 = t690 * t695;
t699 = cos(qJ(4));
t767 = t690 * t699;
t763 = t696 * t697;
t762 = t696 * t700;
t760 = qJD(6) * t694;
t759 = qJD(6) * t698;
t749 = t692 * t779;
t663 = t715 * t690 + t691 * t749;
t757 = t663 * qJD(4);
t756 = t695 * t774;
t754 = pkin(4) * t757;
t753 = t697 * t766;
t752 = t691 * t763;
t748 = t780 * t697;
t747 = t780 * t700;
t744 = qJD(2) * t766;
t742 = t690 * t773;
t740 = t767 * t774;
t738 = t692 * t747;
t736 = t697 * t744;
t735 = t700 * t744;
t733 = t696 * t742;
t730 = t691 * t738;
t631 = t644 * t688 + t663 * t687;
t718 = t692 * t762 + t748;
t660 = t718 * t691 + t733;
t668 = t773 * t692 - t700 * t766;
t635 = t660 * t688 + t668 * t687;
t727 = -t660 * t687 + t668 * t688;
t725 = t780 * t742;
t653 = -t669 * t780 - t670 * t764;
t632 = t653 * t688 + t670 * t770;
t655 = t716 * t764 - t710;
t633 = t655 * t688 - t716 * t770;
t719 = -t692 * t763 + t747;
t667 = t719 * t691;
t651 = t667 * t688 + t687 * t753;
t720 = t669 * t696 - t670 * t750;
t717 = t692 * t748 + t762;
t686 = pkin(4) * t699 + pkin(3);
t709 = -t781 * t687 - t726 * t688 - t686;
t654 = -t715 * t696 - t716 * t750;
t707 = t669 * qJD(1) + t713;
t706 = t688 * t723 + (t726 * t687 - t781 * t688 + t776) * qJD(4);
t705 = t707 * t696;
t704 = t707 * t780;
t693 = -qJ(5) - pkin(11);
t666 = t717 * t691;
t659 = -t725 - t730 + t752;
t656 = t670 * qJD(1) + t715 * qJD(2);
t648 = (-t718 * qJD(2) - t717 * qJD(3)) * t691;
t647 = -qJD(2) * t730 - t691 * qJD(3) * t747 + (qJD(3) * t692 + qJD(2)) * t752;
t637 = qJD(3) * t725 + ((t738 - t763) * qJD(3) + t719 * qJD(2)) * t691;
t636 = qJD(3) * t733 + (t717 * qJD(2) + t718 * qJD(3)) * t691;
t626 = t687 * t735 + t648 * t688 + (-t667 * t687 + t688 * t753) * qJD(4);
t624 = t720 * qJD(3) - t657 * t780 - t658 * t764;
t623 = t653 * qJD(3) - t657 * t696 + t658 * t750;
t622 = -t654 * qJD(3) + t656 * t764 + t704;
t621 = t655 * qJD(3) - t656 * t750 + t705;
t620 = t727 * qJD(4) + t637 * t688 + t687 * t736;
t614 = qJD(1) * t679 + t782 * qJD(3) - t656 * t780 + t692 * t705;
t613 = -qJD(1) * t729 + t644 * qJD(3) - t656 * t696 - t692 * t704;
t612 = t658 * t770 + t624 * t688 + (-t653 * t687 + t670 * t769) * qJD(4);
t610 = -t656 * t770 + t622 * t688 + (-t655 * t687 - t716 * t769) * qJD(4);
t604 = (t614 + t757) * t688 + t792 * t687;
t603 = t631 * qJD(4) + t614 * t687 - t703 * t688;
t602 = t604 * t698 + t613 * t694 + (-t631 * t694 - t698 * t782) * qJD(6);
t601 = -t604 * t694 + t613 * t698 + (-t631 * t698 + t694 * t782) * qJD(6);
t1 = [t799 * r_i_i_C(1) - t800 * r_i_i_C(2) + t608 * pkin(5) - t616 * t686 - t617 * t693 - t639 * qJD(5) - t658 * pkin(2) - pkin(10) * t772 + t781 * t796 + (t798 * r_i_i_C(1) + t797 * r_i_i_C(2)) * qJD(6) + (-t701 * pkin(1) + (-t779 * pkin(9) - pkin(10) * t749) * t691) * qJD(1) + (-t646 * t695 + (t661 * t699 - t695 * t793) * qJD(4)) * pkin(4) (t610 * t698 + t621 * t694 + (-t633 * t694 + t654 * t698) * qJD(6)) * r_i_i_C(1) + (-t610 * t694 + t621 * t698 + (-t633 * t698 - t654 * t694) * qJD(6)) * r_i_i_C(2) + t610 * pkin(5) + t622 * t686 - t655 * t756 - t621 * t693 + t654 * qJD(5) - t716 * t740 + t707 * pkin(2) + (-t768 * pkin(4) - t778) * t656 + t781 * (t633 * qJD(4) + t622 * t687 + t656 * t769) (t614 * t694 + t644 * t759) * r_i_i_C(1) + (t614 * t698 - t644 * t760) * r_i_i_C(2) - t614 * t693 + t644 * qJD(5) + t709 * t613 - t706 * t782, -t614 * t776 - t695 * t754 + t783 * t699 + (-r_i_i_C(1) * t760 - r_i_i_C(2) * t759) * (-t644 * t687 + t663 * t688) + t781 * t604 - t726 * t603, t613, r_i_i_C(1) * t601 - r_i_i_C(2) * t602; -pkin(1) * t743 - t656 * pkin(2) + t604 * pkin(5) + t602 * r_i_i_C(1) + t601 * r_i_i_C(2) - t782 * qJD(5) - t613 * t693 + t614 * t686 + t699 * t754 - t707 * t778 + (pkin(10) * t692 + pkin(9)) * qJD(1) * t765 + t783 * t695 + t781 * t603 (t612 * t698 + t623 * t694) * r_i_i_C(1) + (-t612 * t694 + t623 * t698) * r_i_i_C(2) + t612 * pkin(5) + t624 * t686 - t623 * t693 - t720 * qJD(5) - t657 * pkin(2) + t658 * t778 + t781 * (t632 * qJD(4) + t624 * t687 - t658 * t769) + ((-t632 * t694 - t698 * t720) * r_i_i_C(1) + (-t632 * t698 + t694 * t720) * r_i_i_C(2)) * qJD(6) + (t658 * t768 + (-t653 * t695 + t670 * t767) * qJD(4)) * pkin(4) (t616 * t694 - t759 * t793) * r_i_i_C(1) + (t616 * t698 + t760 * t793) * r_i_i_C(2) - t616 * t693 - t793 * qJD(5) - t709 * t617 + t706 * t639, -t781 * t608 + t784 * t723 + t726 * t796 + (-t616 * t695 + t646 * t699 + (t661 * t695 + t699 * t793) * qJD(4)) * pkin(4), -t617, t800 * r_i_i_C(1) + t799 * r_i_i_C(2) + (-t797 * r_i_i_C(1) + t798 * r_i_i_C(2)) * qJD(6); 0 (t626 * t698 - t647 * t694) * r_i_i_C(1) + (-t626 * t694 - t647 * t698) * r_i_i_C(2) + t626 * pkin(5) + t648 * t686 - t667 * t756 + t647 * t693 + t666 * qJD(5) + t781 * (t651 * qJD(4) + t648 * t687 - t688 * t735) + ((-t651 * t694 + t666 * t698) * r_i_i_C(1) + (-t651 * t698 - t666 * t694) * r_i_i_C(2)) * qJD(6) + (t697 * t740 + (-pkin(2) * t697 + (pkin(10) + t776) * t700 * t690) * qJD(2)) * t691 (t637 * t694 + t660 * t759) * r_i_i_C(1) + (t637 * t698 - t660 * t760) * r_i_i_C(2) - t637 * t693 + t660 * qJD(5) + t709 * t636 + t706 * t659, t781 * t620 - t727 * t723 + t726 * (-t635 * qJD(4) - t637 * t687 + t688 * t736) + (t699 * t736 - t637 * t695 + (-t660 * t699 - t668 * t695) * qJD(4)) * pkin(4), t636 (-t620 * t694 + t636 * t698) * r_i_i_C(1) + (-t620 * t698 - t636 * t694) * r_i_i_C(2) + ((-t635 * t698 - t659 * t694) * r_i_i_C(1) + (t635 * t694 - t659 * t698) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
