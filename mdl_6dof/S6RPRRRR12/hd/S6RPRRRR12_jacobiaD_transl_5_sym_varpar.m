% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR12_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR12_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobiaD_transl_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:16
% EndTime: 2019-02-26 21:21:18
% DurationCPUTime: 1.89s
% Computational Cost: add. (1915->180), mult. (6243->309), div. (0->0), fcn. (7022->16), ass. (0->115)
t699 = cos(pkin(6));
t696 = cos(pkin(14));
t707 = cos(qJ(1));
t744 = t707 * t696;
t692 = sin(pkin(14));
t703 = sin(qJ(1));
t747 = t703 * t692;
t683 = -t699 * t744 + t747;
t702 = sin(qJ(3));
t694 = sin(pkin(7));
t695 = sin(pkin(6));
t752 = t695 * t707;
t737 = t694 * t752;
t698 = cos(pkin(7));
t748 = t698 * t702;
t745 = t707 * t692;
t746 = t703 * t696;
t684 = t699 * t745 + t746;
t706 = cos(qJ(3));
t755 = t684 * t706;
t663 = t683 * t748 + t702 * t737 - t755;
t701 = sin(qJ(4));
t705 = cos(qJ(4));
t757 = t683 * t698;
t722 = t737 + t757;
t662 = t684 * t702 + t706 * t722;
t673 = -t683 * t694 + t698 * t752;
t693 = sin(pkin(8));
t697 = cos(pkin(8));
t726 = t662 * t697 + t673 * t693;
t641 = t663 * t705 + t701 * t726;
t655 = t662 * t693 - t673 * t697;
t700 = sin(qJ(5));
t704 = cos(qJ(5));
t782 = -t641 * t700 - t655 * t704;
t781 = t641 * t704 - t655 * t700;
t734 = t699 * t747;
t743 = qJD(1) * t707;
t681 = -qJD(1) * t734 + t696 * t743;
t719 = t699 * t746 + t745;
t680 = t719 * qJD(1);
t753 = t695 * t703;
t733 = qJD(1) * t753;
t730 = t694 * t733;
t714 = -t680 * t698 + t730;
t647 = (t702 * t722 - t755) * qJD(3) - t681 * t702 + t714 * t706;
t759 = t680 * t694;
t670 = t698 * t733 + t759;
t729 = t647 * t697 + t670 * t693;
t780 = qJD(4) * t641 + t729 * t705;
t772 = -t663 * t701 + t705 * t726;
t779 = -t772 * qJD(4) + t701 * t729;
t636 = t647 * t693 - t670 * t697;
t778 = t636 * t700;
t777 = t636 * t704;
t754 = t694 * t699;
t672 = t695 * (t692 * t706 + t696 * t748) + t702 * t754;
t735 = t706 * t754;
t751 = t696 * t698;
t671 = t735 + (-t692 * t702 + t706 * t751) * t695;
t682 = -t695 * t696 * t694 + t699 * t698;
t724 = t671 * t697 + t682 * t693;
t654 = t672 * t705 + t701 * t724;
t718 = t734 - t744;
t721 = t694 * t753 - t698 * t719;
t665 = t702 * t721 - t718 * t706;
t763 = r_i_i_C(3) + pkin(12);
t668 = t672 * qJD(3);
t761 = t668 * t693;
t750 = t697 * t701;
t749 = t697 * t705;
t742 = qJD(3) * t702;
t741 = qJD(3) * t706;
t740 = qJD(5) * t700;
t739 = qJD(5) * t704;
t738 = t695 * qJD(2);
t732 = t695 * t743;
t731 = t695 * t741;
t713 = t718 * t702;
t664 = t706 * t721 + t713;
t675 = t694 * t719 + t698 * t753;
t725 = t664 * t697 + t675 * t693;
t723 = t704 * r_i_i_C(1) - t700 * r_i_i_C(2) + pkin(4);
t651 = -t662 * t705 + t663 * t750;
t652 = t664 * t705 - t665 * t750;
t658 = t671 * t705 - t672 * t750;
t717 = qJD(5) * (-t700 * r_i_i_C(1) - t704 * r_i_i_C(2));
t678 = t683 * qJD(1);
t716 = -t678 * t694 + t698 * t732;
t715 = t678 * t698 + t694 * t732;
t711 = t716 * t693;
t709 = -t665 * t701 + t705 * t725;
t643 = t665 * t705 + t701 * t725;
t708 = -t672 * t701 + t705 * t724;
t679 = t684 * qJD(1);
t645 = -t665 * qJD(3) + t679 * t702 + t715 * t706;
t634 = -t645 * t693 + t697 * t716;
t687 = t707 * t694 * t731;
t667 = t692 * t695 * t742 - qJD(3) * t735 - t731 * t751;
t659 = -t671 * t693 + t682 * t697;
t657 = -t664 * t693 + t675 * t697;
t650 = t687 + (qJD(3) * t757 - t681) * t706 + (qJD(3) * t684 - t714) * t702;
t648 = t702 * t730 + t681 * t706 - t684 * t742 - t687 + (-t680 * t702 - t683 * t741) * t698;
t646 = qJD(3) * t713 + t715 * t702 + (qJD(3) * t721 - t679) * t706;
t638 = t667 * t750 - t668 * t705 + (-t671 * t701 - t672 * t749) * qJD(4);
t633 = qJD(4) * t708 - t667 * t705 - t668 * t750;
t631 = -t648 * t750 + t647 * t705 + (t662 * t701 + t663 * t749) * qJD(4);
t629 = -t646 * t750 + t645 * t705 + (-t664 * t701 - t665 * t749) * qJD(4);
t627 = t650 * t705 - t779;
t625 = t648 * t705 + t779;
t623 = t646 * t705 + (t645 * t697 + t711) * t701 + t709 * qJD(4);
t622 = qJD(4) * t643 - t645 * t749 + t646 * t701 - t705 * t711;
t621 = t623 * t704 + t634 * t700 + (-t643 * t700 + t657 * t704) * qJD(5);
t620 = -t623 * t700 + t634 * t704 + (-t643 * t704 - t657 * t700) * qJD(5);
t1 = [(t627 * t704 + t778) * r_i_i_C(1) + (-t627 * t700 + t777) * r_i_i_C(2) + t627 * pkin(4) + t650 * pkin(3) - t681 * pkin(2) - pkin(10) * t759 + t707 * t738 + t763 * (t650 * t701 + t780) + (t782 * r_i_i_C(1) - t781 * r_i_i_C(2)) * qJD(5) + t636 * pkin(11) + (-t707 * pkin(1) + (-pkin(10) * t698 - qJ(2)) * t753) * qJD(1), t732 (t629 * t704 - t652 * t740) * r_i_i_C(1) + (-t629 * t700 - t652 * t739) * r_i_i_C(2) + t629 * pkin(4) + t645 * pkin(3) + t763 * (qJD(4) * t652 + t645 * t701 + t646 * t749) + ((t646 * t700 + t665 * t739) * r_i_i_C(1) + (t646 * t704 - t665 * t740) * r_i_i_C(2) + t646 * pkin(11)) * t693, -t723 * t622 + t623 * t763 + t709 * t717, t620 * r_i_i_C(1) - t621 * r_i_i_C(2), 0; t703 * t738 - t679 * pkin(2) + t646 * pkin(3) + t623 * pkin(4) + t621 * r_i_i_C(1) + t620 * r_i_i_C(2) + t763 * t622 + (-t703 * pkin(1) + qJ(2) * t752) * qJD(1) + t634 * pkin(11) + t716 * pkin(10), t733 (t631 * t704 - t651 * t740) * r_i_i_C(1) + (-t631 * t700 - t651 * t739) * r_i_i_C(2) + t631 * pkin(4) + t647 * pkin(3) + t763 * (qJD(4) * t651 + t647 * t701 + t648 * t749) + ((t648 * t700 - t663 * t739) * r_i_i_C(1) + (t648 * t704 + t663 * t740) * r_i_i_C(2) + t648 * pkin(11)) * t693, t763 * t625 - t772 * t717 + t723 * (-t648 * t701 + t780) (-t625 * t700 - t777) * r_i_i_C(1) + (-t625 * t704 + t778) * r_i_i_C(2) + (t781 * r_i_i_C(1) + t782 * r_i_i_C(2)) * qJD(5), 0; 0, 0 (t638 * t704 - t658 * t740) * r_i_i_C(1) + (-t638 * t700 - t658 * t739) * r_i_i_C(2) + t638 * pkin(4) - t668 * pkin(3) + t763 * (qJD(4) * t658 - t667 * t749 - t668 * t701) + ((-t667 * t700 + t672 * t739) * r_i_i_C(1) + (-t667 * t704 - t672 * t740) * r_i_i_C(2) - t667 * pkin(11)) * t693, t763 * t633 + t708 * t717 + t723 * (-t654 * qJD(4) + t667 * t701 - t668 * t749) (-t633 * t700 + t704 * t761) * r_i_i_C(1) + (-t633 * t704 - t700 * t761) * r_i_i_C(2) + ((-t654 * t704 - t659 * t700) * r_i_i_C(1) + (t654 * t700 - t659 * t704) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
