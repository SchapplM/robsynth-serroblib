% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR13
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
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR13_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR13_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:12
% EndTime: 2019-02-26 22:23:14
% DurationCPUTime: 2.35s
% Computational Cost: add. (2289->229), mult. (6212->378), div. (0->0), fcn. (6746->16), ass. (0->133)
t697 = sin(qJ(2));
t699 = cos(qJ(2));
t700 = cos(qJ(1));
t769 = cos(pkin(6));
t738 = t700 * t769;
t771 = sin(qJ(1));
t669 = t697 * t738 + t771 * t699;
t729 = t769 * t771;
t711 = t700 * t697 + t699 * t729;
t656 = qJD(1) * t711 + qJD(2) * t669;
t721 = t697 * t729;
t744 = t771 * t697;
t755 = t700 * t699;
t657 = -qJD(1) * t721 - qJD(2) * t744 + (qJD(2) * t769 + qJD(1)) * t755;
t668 = -t699 * t738 + t744;
t691 = sin(pkin(7));
t693 = cos(pkin(7));
t696 = sin(qJ(3));
t692 = sin(pkin(6));
t756 = t700 * t692;
t749 = t691 * t756;
t772 = cos(qJ(3));
t726 = t772 * t749;
t741 = qJD(3) * t772;
t732 = t693 * t741;
t740 = t771 * qJD(1);
t733 = t692 * t740;
t615 = t696 * (-qJD(3) * t669 - t656 * t693 + t691 * t733) - qJD(3) * t726 + t657 * t772 - t668 * t732;
t766 = t656 * t691;
t645 = t693 * t733 + t766;
t689 = pkin(13) + qJ(5);
t687 = sin(t689);
t688 = cos(t689);
t660 = -t668 * t691 + t693 * t756;
t679 = t696 * t749;
t759 = t693 * t696;
t784 = t668 * t759 - t669 * t772 + t679;
t777 = t660 * t688 - t687 * t784;
t607 = qJD(5) * t777 - t615 * t688 - t645 * t687;
t761 = t692 * t691;
t736 = t771 * t761;
t719 = t772 * t736;
t748 = t693 * t772;
t616 = qJD(1) * t719 + t784 * qJD(3) - t656 * t748 - t657 * t696;
t695 = sin(qJ(6));
t698 = cos(qJ(6));
t791 = t607 * t695 - t616 * t698;
t790 = t607 * t698 + t616 * t695;
t627 = -t660 * t687 - t688 * t784;
t638 = t668 * t748 + t669 * t696 + t726;
t789 = t627 * t695 - t638 * t698;
t788 = t627 * t698 + t638 * t695;
t787 = -qJD(5) * t627 - t615 * t687 + t645 * t688;
t770 = pkin(4) * sin(pkin(13));
t730 = (-t770 - pkin(10)) * t691;
t773 = r_i_i_C(3) + pkin(12);
t712 = t721 - t755;
t775 = t696 * t712 + t719;
t643 = -t712 * t772 + (-t693 * t711 + t736) * t696;
t747 = t693 * t771;
t662 = t691 * t711 + t692 * t747;
t774 = -t643 * t687 + t662 * t688;
t720 = qJD(6) * (t695 * r_i_i_C(1) + t698 * r_i_i_C(2));
t763 = t687 * t691;
t762 = t688 * t691;
t760 = t692 * t699;
t758 = t696 * t697;
t757 = t696 * t699;
t754 = qJD(2) * t692;
t753 = qJD(6) * t695;
t752 = qJD(6) * t698;
t751 = t697 * t761;
t750 = t692 * t758;
t746 = t772 * t697;
t745 = t772 * t699;
t743 = qJD(1) * t756;
t742 = t691 * t754;
t739 = t691 * t769;
t737 = t693 * t745;
t735 = t697 * t742;
t734 = t699 * t742;
t731 = t696 * t739;
t727 = t692 * t737;
t630 = t643 * t688 + t662 * t687;
t714 = t693 * t757 + t746;
t659 = t714 * t692 + t731;
t667 = -t691 * t760 + t769 * t693;
t634 = t659 * t688 + t667 * t687;
t724 = -t659 * t687 + t667 * t688;
t723 = t698 * r_i_i_C(1) - t695 * r_i_i_C(2) + pkin(5);
t722 = t772 * t739;
t652 = -t668 * t772 - t669 * t759;
t631 = t652 * t688 + t669 * t763;
t654 = -t711 * t772 + t712 * t759;
t632 = t654 * t688 - t712 * t763;
t715 = -t693 * t758 + t745;
t666 = t715 * t692;
t650 = t666 * t688 + t687 * t751;
t717 = t668 * t696 - t669 * t748;
t716 = t696 * t711 + t712 * t748;
t713 = t693 * t746 + t757;
t709 = t712 * qJD(2);
t686 = cos(pkin(13)) * pkin(4) + pkin(3);
t706 = -t773 * t687 - t723 * t688 - t686;
t705 = t688 * t720 + (t723 * t687 - t773 * t688) * qJD(5);
t704 = qJD(1) * t668 + t709;
t703 = t704 * t696;
t702 = t704 * t772;
t701 = qJD(1) * t660 - t691 * t709;
t694 = -pkin(11) - qJ(4);
t665 = t713 * t692;
t658 = -t722 - t727 + t750;
t655 = qJD(1) * t669 + qJD(2) * t711;
t647 = (-qJD(2) * t714 - qJD(3) * t713) * t692;
t646 = -qJD(2) * t727 - t741 * t760 + (qJD(3) * t693 + qJD(2)) * t750;
t642 = t711 * t748 - t775;
t636 = qJD(3) * t722 + ((t737 - t758) * qJD(3) + t715 * qJD(2)) * t692;
t635 = qJD(3) * t731 + (qJD(2) * t713 + qJD(3) * t714) * t692;
t625 = t687 * t734 + t647 * t688 + (-t666 * t687 + t688 * t751) * qJD(5);
t623 = t717 * qJD(3) - t656 * t772 - t657 * t759;
t622 = t652 * qJD(3) - t656 * t696 + t657 * t748;
t621 = qJD(3) * t716 + t655 * t759 + t702;
t620 = t654 * qJD(3) - t655 * t748 + t703;
t619 = qJD(5) * t724 + t636 * t688 + t687 * t735;
t613 = qJD(1) * t679 + qJD(3) * t775 - t655 * t772 + t693 * t703 - t711 * t732;
t612 = -qJD(1) * t726 + t643 * qJD(3) - t655 * t696 - t693 * t702;
t611 = t657 * t763 + t623 * t688 + (-t652 * t687 + t669 * t762) * qJD(5);
t609 = -t655 * t763 + t621 * t688 + (-t654 * t687 - t712 * t762) * qJD(5);
t603 = qJD(5) * t774 + t613 * t688 + t687 * t701;
t602 = qJD(5) * t630 + t613 * t687 - t688 * t701;
t601 = t603 * t698 + t612 * t695 + (-t630 * t695 + t642 * t698) * qJD(6);
t600 = -t603 * t695 + t612 * t698 + (-t630 * t698 - t642 * t695) * qJD(6);
t1 = [t790 * r_i_i_C(1) - t791 * r_i_i_C(2) + t607 * pkin(5) - t615 * t686 - t616 * t694 - t638 * qJD(4) - t645 * t770 - t657 * pkin(2) - pkin(10) * t766 + t773 * t787 + (t789 * r_i_i_C(1) + t788 * r_i_i_C(2)) * qJD(6) + (-t700 * pkin(1) + (-t771 * pkin(9) - pkin(10) * t747) * t692) * qJD(1) (t609 * t698 + t620 * t695 + (-t632 * t695 - t698 * t716) * qJD(6)) * r_i_i_C(1) + (-t609 * t695 + t620 * t698 + (-t632 * t698 + t695 * t716) * qJD(6)) * r_i_i_C(2) + t609 * pkin(5) + t621 * t686 - t620 * t694 - t716 * qJD(4) + t704 * pkin(2) + t655 * t730 + t773 * (qJD(5) * t632 + t621 * t687 + t655 * t762) (t613 * t695 + t643 * t752) * r_i_i_C(1) + (t613 * t698 - t643 * t753) * r_i_i_C(2) - t613 * t694 + t643 * qJD(4) + t706 * t612 + t705 * t642, t612, -t723 * t602 + t773 * t603 - t720 * t774, r_i_i_C(1) * t600 - t601 * r_i_i_C(2); -pkin(1) * t740 - t655 * pkin(2) + t603 * pkin(5) + pkin(9) * t743 + t601 * r_i_i_C(1) + t600 * r_i_i_C(2) + t642 * qJD(4) - t612 * t694 + t613 * t686 + t701 * t770 + t773 * t602 + (-t691 * t704 + t693 * t743) * pkin(10) (t611 * t698 + t622 * t695) * r_i_i_C(1) + (-t611 * t695 + t622 * t698) * r_i_i_C(2) + t611 * pkin(5) + t623 * t686 - t622 * t694 - t717 * qJD(4) - t656 * pkin(2) - t657 * t730 + t773 * (qJD(5) * t631 + t623 * t687 - t657 * t762) + ((-t631 * t695 - t698 * t717) * r_i_i_C(1) + (-t631 * t698 + t695 * t717) * r_i_i_C(2)) * qJD(6) (t615 * t695 - t752 * t784) * r_i_i_C(1) + (t615 * t698 + t753 * t784) * r_i_i_C(2) - t615 * t694 - t784 * qJD(4) - t706 * t616 + t705 * t638, -t616, -t607 * t773 + t777 * t720 + t723 * t787, t791 * r_i_i_C(1) + t790 * r_i_i_C(2) + (-t788 * r_i_i_C(1) + t789 * r_i_i_C(2)) * qJD(6); 0 (t625 * t698 - t646 * t695) * r_i_i_C(1) + (-t625 * t695 - t646 * t698) * r_i_i_C(2) + t625 * pkin(5) + t647 * t686 + t646 * t694 + t665 * qJD(4) + t773 * (qJD(5) * t650 + t647 * t687 - t688 * t734) + ((-t650 * t695 + t665 * t698) * r_i_i_C(1) + (-t650 * t698 - t665 * t695) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t697 - t699 * t730) * t754 (t636 * t695 + t659 * t752) * r_i_i_C(1) + (t636 * t698 - t659 * t753) * r_i_i_C(2) - t636 * t694 + t659 * qJD(4) + t706 * t635 + t705 * t658, t635, t773 * t619 - t724 * t720 + t723 * (-qJD(5) * t634 - t636 * t687 + t688 * t735) (-t619 * t695 + t635 * t698) * r_i_i_C(1) + (-t619 * t698 - t635 * t695) * r_i_i_C(2) + ((-t634 * t698 - t658 * t695) * r_i_i_C(1) + (t634 * t695 - t658 * t698) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
