% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:44
% EndTime: 2019-02-26 22:20:47
% DurationCPUTime: 2.69s
% Computational Cost: add. (2842->270), mult. (8652->451), div. (0->0), fcn. (9703->16), ass. (0->142)
t712 = sin(pkin(7));
t711 = sin(pkin(13));
t783 = cos(pkin(13));
t788 = cos(qJ(3));
t750 = t788 * t783;
t717 = sin(qJ(3));
t771 = qJD(3) * t717;
t790 = -qJD(3) * t750 + t711 * t771;
t680 = t790 * t712;
t758 = t717 * t783;
t700 = -t711 * t788 - t758;
t685 = t700 * t712;
t713 = sin(pkin(6));
t723 = cos(qJ(1));
t718 = sin(qJ(2));
t719 = sin(qJ(1));
t722 = cos(qJ(2));
t784 = cos(pkin(6));
t756 = t723 * t784;
t694 = t718 * t756 + t719 * t722;
t757 = t719 * t784;
t733 = t718 * t723 + t722 * t757;
t671 = qJD(1) * t733 + qJD(2) * t694;
t751 = t718 * t757;
t772 = qJD(2) * t718;
t775 = t723 * t722;
t672 = -qJD(1) * t751 - t719 * t772 + (qJD(2) * t784 + qJD(1)) * t775;
t714 = cos(pkin(7));
t682 = t790 * t714;
t687 = t700 * t714;
t759 = t788 * qJD(3);
t692 = -qJD(3) * t758 - t711 * t759;
t693 = t718 * t719 - t722 * t756;
t732 = -t711 * t717 + t750;
t728 = -t671 * t687 - t672 * t732 - t693 * t682 - t694 * t692;
t774 = qJD(1) * t719;
t625 = (t680 * t723 - t685 * t774) * t713 - t728;
t762 = t713 * t774;
t666 = t671 * t712 + t714 * t762;
t716 = sin(qJ(5));
t721 = cos(qJ(5));
t777 = t713 * t723;
t652 = -t685 * t777 - t687 * t693 - t694 * t732;
t766 = t714 * t777;
t673 = -t693 * t712 + t766;
t749 = t652 * t716 - t673 * t721;
t611 = qJD(5) * t749 + t625 * t721 + t666 * t716;
t730 = qJD(3) * t700;
t681 = t712 * t730;
t683 = t714 * t730;
t684 = t732 * t712;
t686 = t732 * t714;
t691 = t732 * qJD(3);
t624 = -t671 * t686 + t672 * t700 - t693 * t683 - t694 * t691 + (-t681 * t723 + t684 * t774) * t713;
t715 = sin(qJ(6));
t720 = cos(qJ(6));
t803 = t611 * t715 + t624 * t720;
t802 = -t611 * t720 + t624 * t715;
t636 = t652 * t721 + t673 * t716;
t740 = -t684 * t777 - t686 * t693 + t694 * t700;
t801 = -t636 * t715 + t720 * t740;
t800 = t636 * t720 + t715 * t740;
t799 = qJD(5) * t636 - t625 * t716 + t666 * t721;
t734 = t751 - t775;
t669 = qJD(1) * t693 + qJD(2) * t734;
t670 = qJD(1) * t694 + qJD(2) * t733;
t773 = qJD(1) * t723;
t622 = -t669 * t687 - t670 * t732 + t733 * t682 - t734 * t692 + (-t680 * t719 - t685 * t773) * t713;
t789 = r_i_i_C(3) + pkin(12);
t787 = pkin(3) * t717;
t785 = pkin(10) + qJ(4);
t786 = t712 * t787 + t714 * t785 + pkin(9);
t782 = t692 * t718;
t781 = t712 * t713;
t780 = t712 * t716;
t779 = t712 * t721;
t778 = t713 * t719;
t776 = t717 * t722;
t770 = qJD(6) * t715;
t769 = qJD(6) * t720;
t768 = pkin(3) * t771;
t767 = t718 * t781;
t765 = t722 * t781;
t764 = t714 * t788;
t763 = t788 * t718;
t761 = t713 * t772;
t755 = t784 * t680;
t754 = pkin(3) * t759;
t753 = qJD(2) * t765;
t752 = t712 * t761;
t675 = t712 * t733 + t714 * t778;
t731 = -t685 * t778 + t687 * t733 - t732 * t734;
t638 = t675 * t716 + t721 * t731;
t748 = t675 * t721 - t716 * t731;
t690 = t714 * t784 - t765;
t745 = -t687 * t722 + t718 * t732;
t726 = -t685 * t784 + t713 * t745;
t645 = t690 * t716 + t721 * t726;
t747 = t690 * t721 - t716 * t726;
t746 = t686 * t722 + t700 * t718;
t744 = t687 * t718 + t722 * t732;
t743 = qJD(1) * t788 * t781;
t742 = r_i_i_C(1) * t720 - r_i_i_C(2) * t715 + pkin(5);
t660 = t687 * t694 - t693 * t732;
t646 = t660 * t721 + t694 * t780;
t662 = -t687 * t734 - t732 * t733;
t647 = t662 * t721 - t734 * t780;
t739 = qJD(6) * (-r_i_i_C(1) * t715 - r_i_i_C(2) * t720);
t668 = t744 * t713;
t663 = t668 * t721 + t716 * t767;
t735 = qJD(1) * t766 - t669 * t712;
t725 = t716 * t789 + t721 * t742 + pkin(4);
t724 = t721 * t739 + (-t716 * t742 + t721 * t789) * qJD(5);
t710 = pkin(3) * t788 + pkin(2);
t698 = -qJD(4) * t712 + t714 * t754;
t697 = qJD(4) * t714 + t712 * t754;
t689 = -t712 * t785 + t714 * t787;
t667 = (t686 * t718 - t700 * t722) * t713;
t661 = -t686 * t734 + t700 * t733;
t659 = t686 * t694 + t693 * t700;
t657 = t684 * t784 + t713 * t746;
t654 = t684 * t778 - t686 * t733 - t700 * t734;
t643 = (-qJD(2) * t745 + t682 * t718 + t692 * t722) * t713;
t642 = (qJD(2) * t746 + t683 * t718 + t691 * t722) * t713;
t641 = -t755 + (qJD(2) * t744 - t682 * t722 + t782) * t713;
t640 = t784 * t681 - t686 * t761 + (-t691 * t718 + (qJD(2) * t700 + t683) * t722) * t713;
t639 = t755 - t687 * t761 + (-t782 + (-qJD(2) * t732 + t682) * t722) * t713;
t633 = -t671 * t732 + t672 * t687 + t682 * t694 - t692 * t693;
t632 = -t671 * t700 - t672 * t686 - t683 * t694 + t691 * t693;
t631 = t669 * t732 - t670 * t687 - t682 * t734 - t692 * t733;
t630 = t669 * t700 + t670 * t686 + t683 * t734 + t691 * t733;
t629 = t716 * t753 + t643 * t721 + (-t668 * t716 + t721 * t767) * qJD(5);
t623 = -t680 * t777 + t685 * t762 + t728;
t621 = t669 * t686 - t670 * t700 - t733 * t683 + t734 * t691 + (t681 * t719 + t684 * t773) * t713;
t619 = qJD(5) * t747 + t641 * t721 + t716 * t752;
t617 = t672 * t780 + t633 * t721 + (-t660 * t716 + t694 * t779) * qJD(5);
t615 = -t670 * t780 + t631 * t721 + (-t662 * t716 - t734 * t779) * qJD(5);
t609 = qJD(5) * t748 + t622 * t721 + t716 * t735;
t608 = qJD(5) * t638 + t622 * t716 - t721 * t735;
t607 = t609 * t720 - t621 * t715 + (-t638 * t715 - t654 * t720) * qJD(6);
t606 = -t609 * t715 - t621 * t720 + (-t638 * t720 + t654 * t715) * qJD(6);
t1 = [t802 * r_i_i_C(1) + t803 * r_i_i_C(2) - t611 * pkin(5) - t625 * pkin(4) + t624 * pkin(11) - t672 * t710 + t694 * t768 + t671 * t689 + t693 * t698 + t697 * t777 + t789 * t799 + (r_i_i_C(1) * t801 - r_i_i_C(2) * t800) * qJD(6) + (-pkin(1) * t723 - t778 * t786) * qJD(1) (t615 * t720 - t630 * t715) * r_i_i_C(1) + (-t615 * t715 - t630 * t720) * r_i_i_C(2) + t615 * pkin(5) + t631 * pkin(4) - t630 * pkin(11) + t669 * t710 + t733 * t768 + t670 * t689 + t734 * t698 + t789 * (qJD(5) * t647 + t631 * t716 + t670 * t779) + ((-t647 * t715 + t661 * t720) * r_i_i_C(1) + (-t647 * t720 - t661 * t715) * r_i_i_C(2)) * qJD(6) (t622 * t715 + t731 * t769) * r_i_i_C(1) + (t622 * t720 - t731 * t770) * r_i_i_C(2) + t622 * pkin(11) + t725 * t621 + t724 * t654 + (t723 * t743 + t669 * t764 + t670 * t717 + (t788 * t734 + (-t712 * t778 + t714 * t733) * t717) * qJD(3)) * pkin(3), t735, -t608 * t742 + t609 * t789 + t739 * t748, r_i_i_C(1) * t606 - r_i_i_C(2) * t607; t734 * t768 + t697 * t778 + t622 * pkin(4) + t609 * pkin(5) - t621 * pkin(11) + t607 * r_i_i_C(1) + t606 * r_i_i_C(2) + t669 * t689 - t670 * t710 - t733 * t698 + t789 * t608 + (-pkin(1) * t719 + t777 * t786) * qJD(1) (t617 * t720 - t632 * t715) * r_i_i_C(1) + (-t617 * t715 - t632 * t720) * r_i_i_C(2) + t617 * pkin(5) + t633 * pkin(4) - t632 * pkin(11) - t671 * t710 + t693 * t768 - t672 * t689 - t694 * t698 + t789 * (qJD(5) * t646 + t633 * t716 - t672 * t779) + ((-t646 * t715 + t659 * t720) * r_i_i_C(1) + (-t646 * t720 - t659 * t715) * r_i_i_C(2)) * qJD(6) (-t623 * t715 - t652 * t769) * r_i_i_C(1) + (-t623 * t720 + t652 * t770) * r_i_i_C(2) - t623 * pkin(11) + t725 * t624 + t724 * t740 + (t719 * t743 - t671 * t764 - t672 * t717 + (-t788 * t694 + (t693 * t714 + t712 * t777) * t717) * qJD(3)) * pkin(3), t666, t611 * t789 + t739 * t749 + t742 * t799, -t803 * r_i_i_C(1) + t802 * r_i_i_C(2) + (r_i_i_C(1) * t800 + r_i_i_C(2) * t801) * qJD(6); 0 (t629 * t720 + t642 * t715) * r_i_i_C(1) + (-t629 * t715 + t642 * t720) * r_i_i_C(2) + t629 * pkin(5) + t643 * pkin(4) + t642 * pkin(11) + t789 * (qJD(5) * t663 + t643 * t716 - t721 * t753) + ((-t663 * t715 + t667 * t720) * r_i_i_C(1) + (-t663 * t720 - t667 * t715) * r_i_i_C(2)) * qJD(6) + (-t722 * t768 - t718 * t698 + (-t689 * t722 - t710 * t718) * qJD(2)) * t713 (-t639 * t715 + t726 * t769) * r_i_i_C(1) + (-t639 * t720 - t726 * t770) * r_i_i_C(2) - t639 * pkin(11) + t725 * t640 + t724 * t657 + (-t784 * t712 * t771 + ((-t714 * t776 - t763) * qJD(3) + (-t714 * t763 - t776) * qJD(2)) * t713) * pkin(3), t752, t789 * t619 + t747 * t739 + t742 * (-qJD(5) * t645 - t641 * t716 + t721 * t752) (-t619 * t715 - t640 * t720) * r_i_i_C(1) + (-t619 * t720 + t640 * t715) * r_i_i_C(2) + ((-t645 * t720 + t657 * t715) * r_i_i_C(1) + (t645 * t715 + t657 * t720) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
