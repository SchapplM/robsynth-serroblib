% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR15_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR15_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:53
% EndTime: 2019-02-26 22:38:55
% DurationCPUTime: 2.13s
% Computational Cost: add. (2360->202), mult. (7270->337), div. (0->0), fcn. (7905->14), ass. (0->124)
t670 = cos(qJ(2));
t671 = cos(qJ(1));
t748 = cos(pkin(6));
t719 = t671 * t748;
t667 = sin(qJ(2));
t750 = sin(qJ(1));
t724 = t750 * t667;
t645 = -t670 * t719 + t724;
t663 = cos(pkin(7));
t666 = sin(qJ(3));
t646 = t667 * t719 + t750 * t670;
t751 = cos(qJ(3));
t729 = t646 * t751;
t661 = sin(pkin(7));
t662 = sin(pkin(6));
t739 = t662 * t671;
t731 = t661 * t739;
t619 = (t645 * t663 + t731) * t666 - t729;
t637 = -t645 * t661 + t663 * t739;
t665 = sin(qJ(4));
t669 = cos(qJ(4));
t604 = t619 * t665 - t637 * t669;
t706 = t751 * t731;
t728 = t663 * t751;
t756 = -t645 * t728 - t706;
t616 = t646 * t666 - t756;
t664 = sin(qJ(6));
t668 = cos(qJ(6));
t769 = t604 * t664 - t616 * t668;
t768 = t604 * t668 + t616 * t664;
t710 = t748 * t750;
t684 = t671 * t667 + t670 * t710;
t633 = t684 * qJD(1) + t646 * qJD(2);
t701 = t667 * t710;
t735 = t671 * t670;
t634 = -qJD(1) * t701 - qJD(2) * t724 + (qJD(2) * t748 + qJD(1)) * t735;
t721 = qJD(1) * t750;
t712 = t662 * t721;
t591 = (-qJD(3) * t646 - t633 * t663 + t661 * t712) * t666 + t634 * t751 + t756 * qJD(3);
t746 = t633 * t661;
t686 = t663 * t712 + t746;
t767 = t604 * qJD(4) + t591 * t669 + t686 * t665;
t763 = t619 * t669 + t637 * t665;
t766 = t763 * qJD(4) - t591 * t665 + t686 * t669;
t755 = pkin(5) + pkin(11);
t685 = t701 - t735;
t743 = t661 * t662;
t715 = t750 * t743;
t621 = -t685 * t751 + (-t684 * t663 + t715) * t666;
t727 = t663 * t750;
t639 = t684 * t661 + t662 * t727;
t754 = -t621 * t665 + t639 * t669;
t733 = r_i_i_C(3) + pkin(12) + pkin(4);
t709 = r_i_i_C(1) * t668 - r_i_i_C(2) * t664;
t699 = t709 + t755;
t679 = t684 * t751;
t694 = t751 * t715;
t753 = -t663 * t679 + t666 * t685 + t694;
t718 = t666 * t731;
t738 = t663 * t666;
t752 = qJD(1) * t694 - t633 * t728 - t634 * t666 + (t645 * t738 + t718 - t729) * qJD(3);
t687 = t709 * qJD(6) + qJD(5);
t749 = pkin(10) * t661;
t742 = t661 * t665;
t741 = t661 * t669;
t740 = t661 * t670;
t737 = t666 * t667;
t736 = t666 * t670;
t734 = qJD(2) * t662;
t732 = t667 * t743;
t730 = t662 * t737;
t726 = t751 * t667;
t725 = t751 * t670;
t722 = t661 * t734;
t720 = t661 * t748;
t716 = t663 * t725;
t714 = t667 * t722;
t713 = t670 * t722;
t711 = t666 * t720;
t708 = -r_i_i_C(1) * t664 - r_i_i_C(2) * t668;
t707 = t662 * t716;
t704 = t621 * t669 + t639 * t665;
t689 = t663 * t736 + t726;
t636 = t689 * t662 + t711;
t644 = -t662 * t740 + t748 * t663;
t703 = t636 * t669 + t644 * t665;
t612 = t636 * t665 - t644 * t669;
t702 = t751 * t720;
t700 = qJ(5) - t708;
t628 = -t645 * t751 - t646 * t738;
t698 = -t628 * t665 + t646 * t741;
t630 = t685 * t738 - t679;
t697 = -t630 * t665 - t685 * t741;
t695 = qJD(6) * t708;
t690 = -t663 * t737 + t725;
t643 = t690 * t662;
t693 = -t643 * t665 + t669 * t732;
t691 = t645 * t666 - t646 * t728;
t688 = t663 * t726 + t736;
t682 = t685 * qJD(2);
t677 = -t700 * t665 - t733 * t669 - pkin(3);
t629 = -t684 * t666 - t685 * t728;
t676 = t645 * qJD(1) + t682;
t675 = t676 * t666;
t674 = -t687 * t665 + (t733 * t665 - t700 * t669) * qJD(4);
t673 = t676 * t751;
t672 = t637 * qJD(1) - t661 * t682;
t642 = t688 * t662;
t635 = -t702 - t707 + t730;
t632 = t646 * qJD(1) + t684 * qJD(2);
t624 = (-t689 * qJD(2) - t688 * qJD(3)) * t662;
t610 = qJD(3) * t702 + ((t716 - t737) * qJD(3) + t690 * qJD(2)) * t662;
t609 = qJD(3) * t711 + (t688 * qJD(2) + t689 * qJD(3)) * t662;
t599 = t691 * qJD(3) - t633 * t751 - t634 * t738;
t597 = -t629 * qJD(3) + t632 * t738 + t673;
t594 = t703 * qJD(4) + t610 * t665 - t669 * t714;
t589 = qJD(1) * t718 + t753 * qJD(3) - t632 * t751 + t663 * t675;
t588 = -qJD(1) * t706 + t621 * qJD(3) - t632 * t666 - t663 * t673;
t584 = t632 * t741 + t597 * t665 + (t630 * t669 - t685 * t742) * qJD(4);
t579 = t754 * qJD(4) + t589 * t669 + t672 * t665;
t578 = t704 * qJD(4) + t589 * t665 - t672 * t669;
t577 = t578 * t664 + t588 * t668 + (t664 * t753 - t668 * t754) * qJD(6);
t576 = t578 * t668 - t588 * t664 + (t664 * t754 + t668 * t753) * qJD(6);
t1 = [-pkin(10) * t746 - t634 * pkin(2) - t591 * pkin(3) + t604 * qJD(5) + t700 * t766 + t699 * t752 + (t768 * r_i_i_C(1) - r_i_i_C(2) * t769) * qJD(6) - t733 * t767 + (-t671 * pkin(1) + (-t750 * pkin(9) - pkin(10) * t727) * t662) * qJD(1) (t584 * t664 + (-t629 * t664 - t668 * t697) * qJD(6)) * r_i_i_C(1) + (t584 * t668 + (-t629 * t668 + t664 * t697) * qJD(6)) * r_i_i_C(2) + t584 * qJ(5) - t697 * qJD(5) + t597 * pkin(3) + t676 * pkin(2) - t632 * t749 + t699 * (t630 * qJD(3) - t632 * t728 + t675) + t733 * (t697 * qJD(4) + t597 * t669 - t632 * t742) t677 * t588 + t699 * t589 + t621 * t695 - t674 * t753, -t733 * t578 + t700 * t579 + t687 * t704, t578, r_i_i_C(1) * t576 - t577 * r_i_i_C(2); -pkin(1) * t721 - t632 * pkin(2) + t589 * pkin(3) + t577 * r_i_i_C(1) + t576 * r_i_i_C(2) + t578 * qJ(5) - t754 * qJD(5) - t676 * t749 + (pkin(10) * t663 + pkin(9)) * qJD(1) * t739 + t755 * t588 + t733 * t579, t634 * t749 - t633 * pkin(2) + t599 * pkin(3) - t698 * qJD(5) + t700 * (-t634 * t741 + t599 * t665 + (t628 * t669 + t646 * t742) * qJD(4)) + t699 * (t628 * qJD(3) - t633 * t666 + t634 * t728) + ((t664 * t691 - t668 * t698) * r_i_i_C(1) + (t664 * t698 + t668 * t691) * r_i_i_C(2)) * qJD(6) + t733 * (t698 * qJD(4) + t599 * t669 + t634 * t742) t699 * t591 + t674 * t616 - t619 * t695 - t677 * t752, -t687 * t763 + t700 * t767 + t733 * t766, -t766 (t664 * t752 - t668 * t766) * r_i_i_C(1) + (t664 * t766 + t668 * t752) * r_i_i_C(2) + (r_i_i_C(1) * t769 + t768 * r_i_i_C(2)) * qJD(6); 0, t624 * pkin(3) - t693 * qJD(5) + t700 * (-t669 * t713 + t624 * t665 + (t643 * t669 + t665 * t732) * qJD(4)) - t699 * (-qJD(2) * t707 - t662 * qJD(3) * t725 + (qJD(3) * t663 + qJD(2)) * t730) + ((-t642 * t664 - t668 * t693) * r_i_i_C(1) + (-t642 * t668 + t664 * t693) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t667 + pkin(10) * t740) * t734 + t733 * (t693 * qJD(4) + t624 * t669 + t665 * t713) t677 * t609 + t699 * t610 + t674 * t635 + t636 * t695, t687 * t703 + t700 * (-t612 * qJD(4) + t610 * t669 + t665 * t714) - t733 * t594, t594 (t594 * t668 - t609 * t664) * r_i_i_C(1) + (-t594 * t664 - t609 * t668) * r_i_i_C(2) + ((-t612 * t664 - t635 * t668) * r_i_i_C(1) + (-t612 * t668 + t635 * t664) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
