% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JaD_transl [3x7]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S7RRRRRRR1_jacobiaD_transl_6_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_6_floatb_twist_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_6_floatb_twist_sym_varpar: qJD has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobiaD_transl_6_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_6_floatb_twist_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 21:20:59
% EndTime: 2018-11-26 21:21:02
% DurationCPUTime: 3.14s
% Computational Cost: add. (1171->258), mult. (3601->450), div. (0->0), fcn. (3789->12), ass. (0->148)
t692 = sin(qJ(3));
t694 = sin(qJ(1));
t698 = cos(qJ(3));
t699 = cos(qJ(2));
t765 = qJD(1) * t699;
t727 = qJD(3) + t765;
t758 = qJD(3) * t699;
t735 = t692 * t758;
t693 = sin(qJ(2));
t763 = qJD(2) * t693;
t740 = t694 * t763;
t766 = qJD(1) * t694;
t700 = cos(qJ(1));
t767 = t700 * t698;
t643 = -t692 * t766 - t694 * t735 - t698 * t740 + t727 * t767;
t768 = t700 * t692;
t774 = t694 * t699;
t665 = t698 * t774 + t768;
t691 = sin(qJ(4));
t697 = cos(qJ(4));
t762 = qJD(2) * t699;
t764 = qJD(1) * t700;
t713 = t693 * t764 + t694 * t762;
t756 = qJD(4) * t693;
t622 = (-qJD(4) * t665 + t713) * t691 + (t694 * t756 + t643) * t697;
t715 = t694 * t698 + t699 * t768;
t642 = t715 * qJD(1) + t665 * qJD(3) - t692 * t740;
t690 = sin(qJ(5));
t696 = cos(qJ(5));
t778 = t693 * t691;
t651 = t665 * t697 + t694 * t778;
t664 = t692 * t774 - t767;
t793 = t651 * t690 + t664 * t696;
t611 = t793 * qJD(5) - t622 * t696 + t642 * t690;
t621 = t651 * qJD(4) + t643 * t691 - t713 * t697;
t689 = sin(qJ(6));
t695 = cos(qJ(6));
t805 = t611 * t689 + t621 * t695;
t804 = t611 * t695 - t621 * t689;
t631 = t651 * t696 - t664 * t690;
t777 = t693 * t697;
t650 = t665 * t691 - t694 * t777;
t803 = t631 * t689 - t650 * t695;
t802 = t631 * t695 + t650 * t689;
t609 = -t631 * qJD(5) - t622 * t690 - t642 * t696;
t726 = qJD(4) * t698 - qJD(2);
t718 = t726 * t691;
t760 = qJD(3) * t697;
t799 = (t692 * t760 + t718) * t693;
t761 = qJD(2) * t700;
t712 = t693 * t766 - t699 * t761;
t780 = t692 * t693;
t737 = qJD(3) * t780;
t757 = qJD(4) * t691;
t792 = -t691 * t737 - t697 * t763 - t699 * t757;
t749 = qJD(6) * t695;
t750 = qJD(6) * t689;
t755 = qJD(4) * t697;
t702 = -(t689 * t755 + t691 * t749) * r_i_i_C(1) - (-t691 * t750 + t695 * t755) * r_i_i_C(2) - pkin(3) * t755;
t751 = qJD(5) * t697;
t711 = t690 * t757 - t696 * t751;
t791 = t711 * r_i_i_C(3) + t702;
t790 = r_i_i_C(3) * t690;
t789 = t691 * pkin(3);
t784 = t689 * t691;
t783 = t689 * t696;
t782 = t690 * t697;
t781 = t691 * t695;
t779 = t692 * t696;
t776 = t693 * t698;
t775 = t693 * t700;
t773 = t695 * t696;
t772 = t696 * t697;
t771 = t696 * t698;
t770 = t699 * t691;
t769 = t699 * t697;
t759 = qJD(3) * t698;
t754 = qJD(5) * t690;
t753 = qJD(5) * t692;
t752 = qJD(5) * t696;
t748 = qJD(6) * t696;
t746 = t690 * t780;
t745 = t698 * t769;
t743 = t699 * t764;
t741 = t692 * t762;
t734 = t693 * t759;
t733 = t693 * t753;
t732 = t693 * t755;
t729 = r_i_i_C(3) * qJD(1) * t779;
t728 = qJD(1) + t758;
t725 = qJD(2) * t698 - qJD(4);
t724 = t696 * t741;
t720 = (-t725 * t693 - t735) * t697 + (-t753 - t718) * t699;
t639 = t725 * t769 - t799;
t719 = -t639 + t733;
t669 = -t694 * t692 + t699 * t767;
t656 = t669 * t697 + t691 * t775;
t634 = t656 * t696 - t690 * t715;
t633 = -t656 * t690 - t696 * t715;
t717 = t725 * t699;
t716 = -t690 * t698 - t692 * t772;
t663 = t697 * t776 - t770;
t662 = t691 * t776 + t769;
t714 = (-t689 * r_i_i_C(1) - t695 * r_i_i_C(2) - pkin(3)) * t691;
t661 = t663 * t700;
t708 = t693 * t761 + t727 * t694;
t706 = -qJD(5) * t663 - t734 - t741;
t667 = t745 + t778;
t705 = -qJD(5) * t667 + t692 * t763 - t698 * t758;
t704 = (t690 * t753 - t696 * t759) * r_i_i_C(3) + qJD(2) * pkin(2);
t703 = -t697 * t717 + t799;
t701 = (t689 * t748 + t695 * t754) * r_i_i_C(1) + (-t689 * t754 + t695 * t748) * r_i_i_C(2) - r_i_i_C(3) * t752;
t666 = t698 * t770 - t777;
t660 = t662 * t700;
t659 = t663 * t694;
t658 = t662 * t694;
t657 = t716 * t693;
t655 = t669 * t691 - t697 * t775;
t654 = -t699 * t692 * t690 + t667 * t696;
t649 = t663 * t696 - t746;
t648 = -t663 * t690 - t693 * t779;
t647 = -t661 * t696 + t700 * t746;
t646 = -t659 * t696 + t694 * t746;
t645 = -t669 * t690 - t715 * t772;
t644 = -t664 * t772 - t665 * t690;
t641 = t708 * t698 + t728 * t768;
t640 = t708 * t692 - t728 * t767;
t638 = (t691 * t762 + t732) * t698 + t792;
t636 = t662 * qJD(2) - qJD(4) * t745 + (t735 - t756) * t691;
t629 = -qJD(1) * t661 + t703 * t694;
t628 = t713 * t691 * t698 + t697 * t743 + (t732 * t698 + t792) * t694;
t627 = t663 * t766 + t703 * t700;
t626 = t662 * t766 + (-t726 * t777 + (-t717 + t737) * t691) * t700;
t625 = t716 * t762 + ((-qJD(5) - t760) * t771 + (t696 * t757 + (qJD(3) + t751) * t690) * t692) * t693;
t620 = (t700 * t756 - t641) * t697 + (-qJD(4) * t669 - t712) * t691;
t619 = t656 * qJD(4) - t641 * t691 + t697 * t712;
t618 = t706 * t690 - t719 * t696;
t617 = t719 * t690 + t706 * t696;
t616 = t705 * t690 + t720 * t696;
t615 = (t694 * t733 + t629) * t696 + (qJD(5) * t659 + t713 * t692 + t694 * t734) * t690;
t614 = (t700 * t733 + t627) * t696 + (qJD(5) * t661 - t712 * t692 + t700 * t734) * t690;
t613 = (t664 * t751 - t643) * t690 + (-qJD(5) * t665 - t642 * t697 + t664 * t757) * t696;
t612 = (t715 * t751 + t641) * t690 + (-qJD(5) * t669 + t640 * t697 + t715 * t757) * t696;
t608 = t633 * qJD(5) + t620 * t696 + t640 * t690;
t607 = t634 * qJD(5) + t620 * t690 - t640 * t696;
t606 = t608 * t695 + t619 * t689 + (-t634 * t689 + t655 * t695) * qJD(6);
t605 = -t608 * t689 + t619 * t695 + (-t634 * t695 - t655 * t689) * qJD(6);
t1 = [t804 * r_i_i_C(1) - t805 * r_i_i_C(2) + t609 * r_i_i_C(3) - t621 * pkin(3) + (t803 * r_i_i_C(1) + t802 * r_i_i_C(2)) * qJD(6) + t713 * pkin(2) (t614 * t695 + t626 * t689) * r_i_i_C(1) + (-t614 * t689 + t626 * t695) * r_i_i_C(2) + (t627 * t690 - t661 * t752 - t700 * t724) * r_i_i_C(3) + t626 * pkin(3) + t694 * pkin(2) * t765 + ((-t647 * t689 - t660 * t695) * r_i_i_C(1) + (-t647 * t695 + t660 * t689) * r_i_i_C(2)) * qJD(6) + (t694 * t729 + t704 * t700) * t693 (t612 * t695 + t640 * t784 - t645 * t750) * r_i_i_C(1) + (-t612 * t689 + t640 * t781 - t645 * t749) * r_i_i_C(2) + (t640 * t782 - t641 * t696 - t669 * t754) * r_i_i_C(3) + t640 * t789 + t791 * t715 (-t619 * t773 + t620 * t689 + t656 * t749) * r_i_i_C(1) + (t619 * t783 + t620 * t695 - t656 * t750) * r_i_i_C(2) - t619 * t790 + t620 * pkin(3) + t701 * t655, t608 * r_i_i_C(3) + (t607 * t689 - t633 * t749) * r_i_i_C(2) + (-t607 * t695 - t633 * t750) * r_i_i_C(1), t605 * r_i_i_C(1) - t606 * r_i_i_C(2), 0; t712 * pkin(2) + t619 * pkin(3) + t606 * r_i_i_C(1) + t605 * r_i_i_C(2) + t607 * r_i_i_C(3) (t615 * t695 - t628 * t689) * r_i_i_C(1) + (-t615 * t689 - t628 * t695) * r_i_i_C(2) + (t629 * t690 - t659 * t752 - t694 * t724) * r_i_i_C(3) - t628 * pkin(3) - pkin(2) * t743 + ((-t646 * t689 - t658 * t695) * r_i_i_C(1) + (-t646 * t695 + t658 * t689) * r_i_i_C(2)) * qJD(6) + (t704 * t694 - t700 * t729) * t693 (t613 * t695 - t642 * t784 - t644 * t750) * r_i_i_C(1) + (-t613 * t689 - t642 * t781 - t644 * t749) * r_i_i_C(2) + (-t642 * t782 + t643 * t696 - t665 * t754) * r_i_i_C(3) - t642 * t789 + t791 * t664 (-t621 * t773 + t622 * t689 + t651 * t749) * r_i_i_C(1) + (t621 * t783 + t622 * t695 - t651 * t750) * r_i_i_C(2) - t621 * t790 + t622 * pkin(3) + t701 * t650, -t611 * r_i_i_C(3) + (-t609 * t689 + t749 * t793) * r_i_i_C(2) + (t609 * t695 + t750 * t793) * r_i_i_C(1), t805 * r_i_i_C(1) + t804 * r_i_i_C(2) + (-t802 * r_i_i_C(1) + t803 * r_i_i_C(2)) * qJD(6), 0; 0 (t616 * t695 - t636 * t689) * r_i_i_C(1) + (-t616 * t689 - t636 * t695) * r_i_i_C(2) - t636 * pkin(3) - pkin(2) * t762 + t720 * t790 - t705 * r_i_i_C(3) * t696 + ((-t654 * t689 + t666 * t695) * r_i_i_C(1) + (-t654 * t695 - t666 * t689) * r_i_i_C(2)) * qJD(6) (t625 * t695 - t657 * t750) * r_i_i_C(1) + (-t625 * t689 - t657 * t749) * r_i_i_C(2) + ((-t692 * t782 + t771) * r_i_i_C(3) + t692 * t714) * t762 + ((-r_i_i_C(3) * t754 + (-r_i_i_C(3) * t782 + t714) * qJD(3)) * t698 + ((-qJD(3) * t696 + t711) * r_i_i_C(3) + t702) * t692) * t693 (-t638 * t773 + t639 * t689 + t663 * t749) * r_i_i_C(1) + (t638 * t783 + t639 * t695 - t663 * t750) * r_i_i_C(2) - t638 * t790 + t639 * pkin(3) + t701 * t662, t618 * r_i_i_C(3) + (-t617 * t689 - t648 * t749) * r_i_i_C(2) + (t617 * t695 - t648 * t750) * r_i_i_C(1) (-t618 * t689 + t638 * t695) * r_i_i_C(1) + (-t618 * t695 - t638 * t689) * r_i_i_C(2) + ((-t649 * t695 - t662 * t689) * r_i_i_C(1) + (t649 * t689 - t662 * t695) * r_i_i_C(2)) * qJD(6), 0;];
JaD_transl  = t1;
