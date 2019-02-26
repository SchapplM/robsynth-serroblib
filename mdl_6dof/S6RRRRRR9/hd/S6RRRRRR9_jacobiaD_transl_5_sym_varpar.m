% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function JaD_transl = S6RRRRRR9_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR9_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:14
% EndTime: 2019-02-26 22:52:16
% DurationCPUTime: 2.25s
% Computational Cost: add. (1927->216), mult. (5954->367), div. (0->0), fcn. (6464->14), ass. (0->129)
t677 = sin(qJ(2));
t680 = cos(qJ(2));
t681 = cos(qJ(1));
t751 = cos(pkin(6));
t720 = t681 * t751;
t753 = sin(qJ(1));
t655 = t677 * t720 + t753 * t680;
t711 = t751 * t753;
t692 = t681 * t677 + t680 * t711;
t642 = t692 * qJD(1) + t655 * qJD(2);
t703 = t677 * t711;
t726 = t753 * t677;
t738 = t681 * t680;
t643 = -qJD(1) * t703 - qJD(2) * t726 + (qJD(2) * t751 + qJD(1)) * t738;
t654 = -t680 * t720 + t726;
t671 = sin(pkin(7));
t673 = cos(pkin(7));
t676 = sin(qJ(3));
t672 = sin(pkin(6));
t742 = t672 * t681;
t732 = t671 * t742;
t754 = cos(qJ(3));
t708 = t754 * t732;
t723 = qJD(3) * t754;
t713 = t673 * t723;
t722 = t753 * qJD(1);
t714 = t672 * t722;
t601 = (-qJD(3) * t655 - t642 * t673 + t671 * t714) * t676 - qJD(3) * t708 + t643 * t754 - t654 * t713;
t749 = t642 * t671;
t631 = t673 * t714 + t749;
t675 = sin(qJ(4));
t679 = cos(qJ(4));
t731 = t655 * t754;
t627 = (t654 * t673 + t732) * t676 - t731;
t646 = -t654 * t671 + t673 * t742;
t707 = t627 * t675 - t646 * t679;
t593 = -t707 * qJD(4) - t601 * t679 - t631 * t675;
t744 = t672 * t671;
t717 = t753 * t744;
t700 = t754 * t717;
t719 = t676 * t732;
t730 = t673 * t754;
t741 = t673 * t676;
t602 = qJD(1) * t700 - t642 * t730 - t643 * t676 + (t654 * t741 + t719 - t731) * qJD(3);
t674 = sin(qJ(5));
t678 = cos(qJ(5));
t773 = t593 * t674 - t602 * t678;
t772 = t593 * t678 + t602 * t674;
t614 = t627 * t679 + t646 * t675;
t624 = t654 * t730 + t655 * t676 + t708;
t771 = -t614 * t674 - t624 * t678;
t770 = t614 * t678 - t624 * t674;
t769 = t614 * qJD(4) - t601 * t675 + t631 * t679;
t755 = r_i_i_C(3) + pkin(12);
t693 = t703 - t738;
t757 = t676 * t693 + t700;
t629 = -t693 * t754 + (-t673 * t692 + t717) * t676;
t729 = t673 * t753;
t648 = t671 * t692 + t672 * t729;
t756 = -t629 * t675 + t648 * t679;
t701 = qJD(5) * (t674 * r_i_i_C(1) + t678 * r_i_i_C(2));
t752 = pkin(10) * t671;
t746 = t671 * t675;
t745 = t671 * t679;
t743 = t672 * t680;
t740 = t676 * t677;
t739 = t676 * t680;
t737 = qJD(2) * t672;
t736 = qJD(5) * t674;
t735 = qJD(5) * t678;
t734 = t677 * t744;
t733 = t672 * t740;
t728 = t754 * t677;
t727 = t754 * t680;
t724 = t671 * t737;
t721 = t671 * t751;
t718 = t673 * t727;
t716 = t677 * t724;
t715 = t680 * t724;
t712 = t676 * t721;
t709 = t672 * t718;
t616 = t629 * t679 + t648 * t675;
t695 = t673 * t739 + t728;
t645 = t695 * t672 + t712;
t653 = -t671 * t743 + t751 * t673;
t622 = t645 * t679 + t653 * t675;
t706 = -t645 * t675 + t653 * t679;
t705 = t678 * r_i_i_C(1) - t674 * r_i_i_C(2) + pkin(4);
t704 = t754 * t721;
t637 = -t654 * t754 - t655 * t741;
t617 = t637 * t679 + t655 * t746;
t639 = -t692 * t754 + t693 * t741;
t618 = t639 * t679 - t693 * t746;
t696 = -t673 * t740 + t727;
t652 = t696 * t672;
t640 = t652 * t679 + t675 * t734;
t698 = t654 * t676 - t655 * t730;
t697 = t676 * t692 + t693 * t730;
t694 = t673 * t728 + t739;
t690 = t693 * qJD(2);
t687 = -t755 * t675 - t705 * t679 - pkin(3);
t686 = t679 * t701 + (t705 * t675 - t755 * t679) * qJD(4);
t685 = t654 * qJD(1) + t690;
t684 = t685 * t676;
t683 = t685 * t754;
t682 = t646 * qJD(1) - t671 * t690;
t651 = t694 * t672;
t644 = -t704 - t709 + t733;
t641 = t655 * qJD(1) + t692 * qJD(2);
t633 = (-t695 * qJD(2) - t694 * qJD(3)) * t672;
t632 = -qJD(2) * t709 - t723 * t743 + (qJD(3) * t673 + qJD(2)) * t733;
t628 = t692 * t730 - t757;
t620 = qJD(3) * t704 + ((t718 - t740) * qJD(3) + t696 * qJD(2)) * t672;
t619 = qJD(3) * t712 + (t694 * qJD(2) + t695 * qJD(3)) * t672;
t611 = t675 * t715 + t633 * t679 + (-t652 * t675 + t679 * t734) * qJD(4);
t609 = t698 * qJD(3) - t642 * t754 - t643 * t741;
t608 = t637 * qJD(3) - t642 * t676 + t643 * t730;
t607 = t697 * qJD(3) + t641 * t741 + t683;
t606 = t639 * qJD(3) - t641 * t730 + t684;
t605 = t706 * qJD(4) + t620 * t679 + t675 * t716;
t599 = qJD(1) * t719 + qJD(3) * t757 - t641 * t754 + t673 * t684 - t692 * t713;
t598 = -qJD(1) * t708 + t629 * qJD(3) - t641 * t676 - t673 * t683;
t597 = t643 * t746 + t609 * t679 + (-t637 * t675 + t655 * t745) * qJD(4);
t595 = -t641 * t746 + t607 * t679 + (-t639 * t675 - t693 * t745) * qJD(4);
t589 = qJD(4) * t756 + t599 * t679 + t682 * t675;
t588 = t616 * qJD(4) + t599 * t675 - t682 * t679;
t587 = t589 * t678 + t598 * t674 + (-t616 * t674 + t628 * t678) * qJD(5);
t586 = -t589 * t674 + t598 * t678 + (-t616 * t678 - t628 * t674) * qJD(5);
t1 = [t772 * r_i_i_C(1) - t773 * r_i_i_C(2) + t593 * pkin(4) - t601 * pkin(3) + t602 * pkin(11) - t643 * pkin(2) - pkin(10) * t749 + t755 * t769 + (t771 * r_i_i_C(1) - t770 * r_i_i_C(2)) * qJD(5) + (-t681 * pkin(1) + (-t753 * pkin(9) - pkin(10) * t729) * t672) * qJD(1) (t595 * t678 + t606 * t674 + (-t618 * t674 - t678 * t697) * qJD(5)) * r_i_i_C(1) + (-t595 * t674 + t606 * t678 + (-t618 * t678 + t674 * t697) * qJD(5)) * r_i_i_C(2) + t595 * pkin(4) + t607 * pkin(3) + t606 * pkin(11) + t685 * pkin(2) - t641 * t752 + t755 * (t618 * qJD(4) + t607 * t675 + t641 * t745) (t599 * t674 + t629 * t735) * r_i_i_C(1) + (t599 * t678 - t629 * t736) * r_i_i_C(2) + t599 * pkin(11) + t687 * t598 + t686 * t628, -t705 * t588 + t755 * t589 - t701 * t756, r_i_i_C(1) * t586 - r_i_i_C(2) * t587, 0; -pkin(1) * t722 - t641 * pkin(2) + t599 * pkin(3) + t589 * pkin(4) + t598 * pkin(11) + t587 * r_i_i_C(1) + t586 * r_i_i_C(2) - t685 * t752 + (pkin(10) * t673 + pkin(9)) * qJD(1) * t742 + t755 * t588 (t597 * t678 + t608 * t674) * r_i_i_C(1) + (-t597 * t674 + t608 * t678) * r_i_i_C(2) + t597 * pkin(4) + t609 * pkin(3) + t608 * pkin(11) - t642 * pkin(2) + t643 * t752 + t755 * (t617 * qJD(4) + t609 * t675 - t643 * t745) + ((-t617 * t674 - t678 * t698) * r_i_i_C(1) + (-t617 * t678 + t674 * t698) * r_i_i_C(2)) * qJD(5) (t601 * t674 - t627 * t735) * r_i_i_C(1) + (t601 * t678 + t627 * t736) * r_i_i_C(2) + t601 * pkin(11) - t687 * t602 + t686 * t624, -t593 * t755 - t707 * t701 + t705 * t769, t773 * r_i_i_C(1) + t772 * r_i_i_C(2) + (t770 * r_i_i_C(1) + t771 * r_i_i_C(2)) * qJD(5), 0; 0 (t611 * t678 - t632 * t674) * r_i_i_C(1) + (-t611 * t674 - t632 * t678) * r_i_i_C(2) + t611 * pkin(4) + t633 * pkin(3) - t632 * pkin(11) + t755 * (t640 * qJD(4) + t633 * t675 - t679 * t715) + ((-t640 * t674 + t651 * t678) * r_i_i_C(1) + (-t640 * t678 - t651 * t674) * r_i_i_C(2)) * qJD(5) + (-pkin(2) * t677 + t680 * t752) * t737 (t620 * t674 + t645 * t735) * r_i_i_C(1) + (t620 * t678 - t645 * t736) * r_i_i_C(2) + t620 * pkin(11) + t687 * t619 + t686 * t644, t755 * t605 - t706 * t701 + t705 * (-t622 * qJD(4) - t620 * t675 + t679 * t716) (-t605 * t674 + t619 * t678) * r_i_i_C(1) + (-t605 * t678 - t619 * t674) * r_i_i_C(2) + ((-t622 * t678 - t644 * t674) * r_i_i_C(1) + (t622 * t674 - t644 * t678) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
