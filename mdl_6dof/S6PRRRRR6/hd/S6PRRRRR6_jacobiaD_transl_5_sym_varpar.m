% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobiaD_transl_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:22:06
% EndTime: 2019-02-26 20:22:07
% DurationCPUTime: 1.36s
% Computational Cost: add. (1882->235), mult. (6202->424), div. (0->0), fcn. (6949->16), ass. (0->141)
t674 = sin(pkin(6));
t677 = cos(pkin(7));
t682 = sin(qJ(2));
t685 = cos(qJ(3));
t729 = t682 * t685;
t681 = sin(qJ(3));
t686 = cos(qJ(2));
t730 = t681 * t686;
t695 = t677 * t730 + t729;
t673 = sin(pkin(7));
t678 = cos(pkin(6));
t741 = t673 * t678;
t718 = t681 * t741;
t649 = t695 * t674 + t718;
t680 = sin(qJ(4));
t684 = cos(qJ(4));
t728 = t685 * t686;
t731 = t681 * t682;
t697 = t677 * t728 - t731;
t648 = t697 * t674 + t685 * t741;
t738 = t674 * t686;
t660 = -t673 * t738 + t678 * t677;
t672 = sin(pkin(8));
t676 = cos(pkin(8));
t711 = t648 * t676 + t660 * t672;
t621 = t649 * t684 + t711 * t680;
t671 = sin(pkin(14));
t675 = cos(pkin(14));
t733 = t678 * t682;
t698 = t671 * t733 - t675 * t686;
t732 = t678 * t686;
t699 = t671 * t732 + t675 * t682;
t743 = t673 * t674;
t700 = t671 * t743 - t677 * t699;
t640 = t700 * t681 - t685 * t698;
t639 = t681 * t698 + t700 * t685;
t740 = t674 * t677;
t651 = t671 * t740 + t673 * t699;
t712 = t639 * t676 + t651 * t672;
t607 = t640 * t684 + t712 * t680;
t719 = t675 * t732;
t661 = -t671 * t682 + t719;
t721 = t675 * t743;
t735 = t677 * t681;
t662 = t671 * t686 + t675 * t733;
t746 = t662 * t685;
t638 = t661 * t735 - t681 * t721 + t746;
t702 = -t661 * t677 + t721;
t637 = -t662 * t681 - t702 * t685;
t650 = -t661 * t673 - t675 * t740;
t713 = t637 * t676 + t650 * t672;
t605 = t638 * t684 + t713 * t680;
t753 = r_i_i_C(3) + pkin(12);
t752 = t673 * pkin(10);
t694 = t677 * t731 - t728;
t641 = (-t697 * qJD(2) + t694 * qJD(3)) * t674;
t749 = t641 * t672;
t657 = t662 * qJD(2);
t747 = t657 * t681;
t744 = t672 * t673;
t742 = t673 * t676;
t739 = t674 * t682;
t737 = t676 * t680;
t736 = t676 * t684;
t734 = t677 * t685;
t727 = qJD(2) * t674;
t726 = qJD(2) * t682;
t725 = qJD(3) * t681;
t724 = qJD(3) * t685;
t679 = sin(qJ(5));
t723 = qJD(5) * t679;
t683 = cos(qJ(5));
t722 = qJD(5) * t683;
t720 = t673 * t739;
t717 = t674 * t726;
t716 = t686 * t727;
t715 = t673 * t717;
t714 = t673 * t716;
t710 = t683 * r_i_i_C(1) - t679 * r_i_i_C(2) + pkin(4);
t656 = -qJD(2) * t719 + t671 * t726;
t616 = -t657 * t734 + t656 * t681 + (t702 * t681 - t746) * qJD(3);
t709 = t616 * t676 + t657 * t744;
t658 = t699 * qJD(2);
t659 = t698 * qJD(2);
t618 = -t640 * qJD(3) + t658 * t681 + t659 * t734;
t708 = t618 * t676 - t659 * t744;
t703 = -t661 * t685 + t662 * t735;
t622 = t703 * qJD(3) + t656 * t734 + t747;
t610 = -t622 * t672 - t656 * t742;
t707 = -t622 * t676 + t656 * t744;
t701 = t685 * t699 - t698 * t735;
t624 = t701 * qJD(3) + t658 * t734 - t659 * t681;
t611 = -t624 * t672 - t658 * t742;
t706 = -t624 * t676 + t658 * t744;
t614 = t637 * t684 - t638 * t737;
t615 = t639 * t684 - t640 * t737;
t643 = -t661 * t681 - t662 * t734;
t705 = t643 * t676 + t662 * t744;
t645 = t681 * t699 + t698 * t734;
t704 = t645 * t676 - t698 * t744;
t629 = t648 * t684 - t649 * t737;
t696 = -t677 * t729 - t730;
t693 = qJD(5) * (-t679 * r_i_i_C(1) - t683 * r_i_i_C(2));
t652 = t696 * t674;
t692 = t652 * t676 + t672 * t720;
t634 = -qJD(3) * t718 + (t696 * qJD(2) - t695 * qJD(3)) * t674;
t691 = t634 * t676 + t672 * t715;
t690 = t641 * t676 + t672 * t714;
t689 = -t638 * t680 + t713 * t684;
t688 = -t640 * t680 + t712 * t684;
t687 = -t649 * t680 + t711 * t684;
t612 = t705 * t680 - t684 * t703;
t613 = t704 * t680 - t684 * t701;
t653 = t694 * t674;
t632 = -t653 * t684 + t692 * t680;
t647 = -t652 * t672 + t676 * t720;
t642 = (-t695 * qJD(2) + t696 * qJD(3)) * t674;
t636 = -t648 * t672 + t660 * t676;
t635 = -t717 * t735 - t725 * t739 + (t716 + (t677 * t738 + t741) * qJD(3)) * t685;
t633 = t676 * t714 - t749;
t631 = -t645 * t672 - t698 * t742;
t630 = -t643 * t672 + t662 * t742;
t628 = -t634 * t672 + t676 * t715;
t627 = -t639 * t672 + t651 * t676;
t626 = -t637 * t672 + t650 * t676;
t625 = t645 * qJD(3) + t658 * t735 + t659 * t685;
t623 = t643 * qJD(3) + t656 * t735 - t657 * t685;
t619 = t659 * t735 + t698 * t725 + (t700 * qJD(3) - t658) * t685;
t617 = -t656 * t685 - t662 * t725 - t721 * t724 + (t661 * t724 - t747) * t677;
t609 = -t618 * t672 - t659 * t742;
t608 = -t616 * t672 + t657 * t742;
t603 = t642 * t684 + t690 * t680 + (t653 * t680 + t692 * t684) * qJD(4);
t601 = -t635 * t737 + t634 * t684 + (-t648 * t680 - t649 * t736) * qJD(4);
t599 = t687 * qJD(4) + t635 * t684 + t691 * t680;
t597 = -t619 * t737 + t618 * t684 + (-t639 * t680 - t640 * t736) * qJD(4);
t595 = -t617 * t737 + t616 * t684 + (-t637 * t680 - t638 * t736) * qJD(4);
t593 = t625 * t684 - t706 * t680 + (t680 * t701 + t704 * t684) * qJD(4);
t591 = t623 * t684 - t707 * t680 + (t680 * t703 + t705 * t684) * qJD(4);
t589 = t688 * qJD(4) + t619 * t684 + t708 * t680;
t587 = t689 * qJD(4) + t617 * t684 + t709 * t680;
t1 = [0 (t593 * t683 + t611 * t679) * r_i_i_C(1) + (-t593 * t679 + t611 * t683) * r_i_i_C(2) + t593 * pkin(4) + t625 * pkin(3) + t659 * pkin(2) - t658 * t752 + t753 * (t613 * qJD(4) + t625 * t680 + t706 * t684) + ((-t613 * t679 + t631 * t683) * r_i_i_C(1) + (-t613 * t683 - t631 * t679) * r_i_i_C(2)) * qJD(5) + t611 * pkin(11) (t597 * t683 - t615 * t723) * r_i_i_C(1) + (-t597 * t679 - t615 * t722) * r_i_i_C(2) + t597 * pkin(4) + t618 * pkin(3) + t753 * (t615 * qJD(4) + t618 * t680 + t619 * t736) + ((t619 * t679 + t640 * t722) * r_i_i_C(1) + (t619 * t683 - t640 * t723) * r_i_i_C(2) + t619 * pkin(11)) * t672, t753 * t589 + t688 * t693 + t710 * (-t607 * qJD(4) - t619 * t680 + t708 * t684) (-t589 * t679 + t609 * t683) * r_i_i_C(1) + (-t589 * t683 - t609 * t679) * r_i_i_C(2) + ((-t607 * t683 - t627 * t679) * r_i_i_C(1) + (t607 * t679 - t627 * t683) * r_i_i_C(2)) * qJD(5), 0; 0 (t591 * t683 + t610 * t679) * r_i_i_C(1) + (-t591 * t679 + t610 * t683) * r_i_i_C(2) + t591 * pkin(4) + t623 * pkin(3) - t657 * pkin(2) - t656 * t752 + t753 * (t612 * qJD(4) + t623 * t680 + t707 * t684) + ((-t612 * t679 + t630 * t683) * r_i_i_C(1) + (-t612 * t683 - t630 * t679) * r_i_i_C(2)) * qJD(5) + t610 * pkin(11) (t595 * t683 - t614 * t723) * r_i_i_C(1) + (-t595 * t679 - t614 * t722) * r_i_i_C(2) + t595 * pkin(4) + t616 * pkin(3) + t753 * (t614 * qJD(4) + t616 * t680 + t617 * t736) + ((t617 * t679 + t638 * t722) * r_i_i_C(1) + (t617 * t683 - t638 * t723) * r_i_i_C(2) + t617 * pkin(11)) * t672, t753 * t587 + t689 * t693 + t710 * (-t605 * qJD(4) - t617 * t680 + t709 * t684) (-t587 * t679 + t608 * t683) * r_i_i_C(1) + (-t587 * t683 - t608 * t679) * r_i_i_C(2) + ((-t605 * t683 - t626 * t679) * r_i_i_C(1) + (t605 * t679 - t626 * t683) * r_i_i_C(2)) * qJD(5), 0; 0 (t603 * t683 + t633 * t679) * r_i_i_C(1) + (-t603 * t679 + t633 * t683) * r_i_i_C(2) + t603 * pkin(4) + t642 * pkin(3) - pkin(11) * t749 + t753 * (t632 * qJD(4) + t642 * t680 - t690 * t684) + ((-t632 * t679 + t647 * t683) * r_i_i_C(1) + (-t632 * t683 - t647 * t679) * r_i_i_C(2)) * qJD(5) + (-pkin(2) * t682 + (pkin(11) * t676 + pkin(10)) * t686 * t673) * t727 (t601 * t683 - t629 * t723) * r_i_i_C(1) + (-t601 * t679 - t629 * t722) * r_i_i_C(2) + t601 * pkin(4) + t634 * pkin(3) + t753 * (t629 * qJD(4) + t634 * t680 + t635 * t736) + ((t635 * t679 + t649 * t722) * r_i_i_C(1) + (t635 * t683 - t649 * t723) * r_i_i_C(2) + t635 * pkin(11)) * t672, t753 * t599 + t687 * t693 + t710 * (-t621 * qJD(4) - t635 * t680 + t691 * t684) (-t599 * t679 + t628 * t683) * r_i_i_C(1) + (-t599 * t683 - t628 * t679) * r_i_i_C(2) + ((-t621 * t683 - t636 * t679) * r_i_i_C(1) + (t621 * t679 - t636 * t683) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
