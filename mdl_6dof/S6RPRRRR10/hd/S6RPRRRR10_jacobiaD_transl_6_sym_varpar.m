% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR10_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR10_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:55
% EndTime: 2019-02-26 21:19:57
% DurationCPUTime: 1.40s
% Computational Cost: add. (1920->158), mult. (4975->261), div. (0->0), fcn. (5524->16), ass. (0->105)
t639 = qJ(4) + qJ(5);
t637 = cos(t639);
t646 = cos(qJ(1));
t713 = cos(pkin(13));
t715 = cos(pkin(6));
t681 = t715 * t713;
t711 = sin(pkin(13));
t719 = sin(qJ(1));
t662 = t646 * t711 + t719 * t681;
t616 = t662 * qJD(1);
t679 = t715 * t711;
t621 = t646 * t713 - t719 * t679;
t617 = t621 * qJD(1);
t643 = sin(qJ(3));
t660 = -t646 * t679 - t719 * t713;
t640 = sin(pkin(6));
t712 = sin(pkin(7));
t691 = t712 * t640;
t673 = t719 * t691;
t714 = cos(pkin(7));
t720 = cos(qJ(3));
t661 = t646 * t681 - t719 * t711;
t684 = t646 * t691;
t672 = t720 * t684;
t683 = t714 * t720;
t727 = t661 * t683 - t672;
t577 = (qJD(1) * t673 + qJD(3) * t660 - t714 * t616) * t643 + t617 * t720 + t727 * qJD(3);
t638 = qJD(4) + qJD(5);
t706 = t640 * t646;
t725 = t661 * t712 + t714 * t706;
t690 = -t638 * t725 + t577;
t739 = t690 * t637;
t692 = t714 * t661;
t699 = t660 * t720;
t597 = (t684 - t692) * t643 + t699;
t636 = sin(t639);
t587 = t597 * t637 + t636 * t725;
t594 = -t643 * t660 - t727;
t641 = sin(qJ(6));
t644 = cos(qJ(6));
t738 = -t587 * t641 - t594 * t644;
t737 = t587 * t644 - t594 * t641;
t682 = t714 * t719;
t674 = t640 * t682;
t695 = t616 * t712;
t605 = qJD(1) * t674 + t695;
t688 = t597 * t638 + t605;
t736 = -t690 * t636 + t688 * t637;
t668 = t720 * t673;
t676 = t643 * t684;
t578 = qJD(1) * t668 - t616 * t683 - t617 * t643 + (-t643 * t692 + t676 + t699) * qJD(3);
t735 = t578 * t641;
t734 = t578 * t644;
t721 = r_i_i_C(3) + pkin(12);
t732 = -t644 * r_i_i_C(1) - pkin(5);
t705 = qJD(6) * t641;
t702 = r_i_i_C(1) * t705;
t704 = qJD(6) * t644;
t726 = -r_i_i_C(2) * t704 - t702;
t677 = -t641 * r_i_i_C(2) - t732;
t656 = t662 * t714;
t599 = t621 * t720 + (-t656 + t673) * t643;
t650 = t725 * qJD(1);
t717 = pkin(4) * qJD(4);
t724 = pkin(4) * t650 - t599 * t717;
t723 = -t621 * t643 - t720 * t656 + t668;
t678 = t714 * t713;
t680 = t715 * t712;
t606 = -t720 * t680 + (t643 * t711 - t678 * t720) * t640;
t615 = t660 * qJD(1);
t654 = qJD(1) * t692;
t575 = qJD(1) * t676 + qJD(3) * t723 + t615 * t720 - t643 * t654;
t610 = t662 * t712 + t674;
t707 = t637 * t638;
t568 = t599 * t707 - t637 * t650 + (t610 * t638 + t575) * t636;
t708 = t636 * t638;
t569 = t575 * t637 - t599 * t708 + t610 * t707 + t636 * t650;
t588 = -t599 * t636 + t610 * t637;
t722 = (t568 * t641 - t588 * t704) * r_i_i_C(2) - t588 * t702 + t721 * t569 + t732 * t568;
t703 = t640 * qJD(2);
t700 = t610 * t717;
t698 = qJD(1) * t706;
t696 = qJD(1) * t719;
t602 = t606 * qJD(3);
t618 = -t713 * t691 + t715 * t714;
t687 = t618 * t638 - t602;
t645 = cos(qJ(4));
t635 = t645 * pkin(4) + pkin(3);
t659 = -t721 * t636 - t677 * t637 - t635;
t571 = t597 * t708 + t605 * t636 + t739;
t652 = t726 * (t597 * t636 - t637 * t725) + t721 * t571 + t677 * t736;
t607 = t643 * t680 + (t643 * t678 + t720 * t711) * t640;
t584 = -t607 * t708 + t687 * t637;
t651 = t726 * (-t607 * t636 + t618 * t637) + t721 * t584 + t677 * (-t607 * t707 - t687 * t636);
t642 = sin(qJ(4));
t648 = t642 * t717 + (t641 * r_i_i_C(1) + t644 * r_i_i_C(2)) * t637 * qJD(6) + (t677 * t636 - t721 * t637) * t638;
t647 = -pkin(11) - pkin(10);
t603 = t607 * qJD(3);
t593 = t607 * t637 + t618 * t636;
t589 = t599 * t637 + t610 * t636;
t574 = -qJD(1) * t672 + t599 * qJD(3) + t615 * t643 + t720 * t654;
t573 = -t688 * t636 - t739;
t561 = t569 * t644 + t574 * t641 + (-t589 * t641 - t644 * t723) * qJD(6);
t560 = -t569 * t641 + t574 * t644 + (-t589 * t644 + t641 * t723) * qJD(6);
t1 = [(t573 * t644 + t735) * r_i_i_C(1) + (-t573 * t641 + t734) * r_i_i_C(2) + t573 * pkin(5) - t577 * t635 - t578 * t647 - t617 * pkin(2) - pkin(9) * t695 + t646 * t703 + t721 * t736 + (t738 * r_i_i_C(1) - t737 * r_i_i_C(2)) * qJD(6) + (-t646 * pkin(1) + (-pkin(9) * t682 - t719 * qJ(2)) * t640) * qJD(1) + (-t605 * t642 + (-t597 * t642 + t645 * t725) * qJD(4)) * pkin(4), t698 (t575 * t641 + t599 * t704) * r_i_i_C(1) + (t575 * t644 - t599 * t705) * r_i_i_C(2) - t575 * t647 + t659 * t574 - t648 * t723, t724 * t645 + (-pkin(4) * t575 - t700) * t642 + t722, t722, t560 * r_i_i_C(1) - t561 * r_i_i_C(2); -pkin(1) * t696 + t615 * pkin(2) + t569 * pkin(5) + pkin(9) * t650 + t561 * r_i_i_C(1) + t560 * r_i_i_C(2) + qJ(2) * t698 + t568 * t721 - t574 * t647 + t575 * t635 + t724 * t642 + t645 * t700 + t719 * t703, t640 * t696 (t577 * t641 - t597 * t704) * r_i_i_C(1) + (t577 * t644 + t597 * t705) * r_i_i_C(2) - t577 * t647 - t659 * t578 + t648 * t594 (-t577 * t642 + t605 * t645 + (t597 * t645 + t642 * t725) * qJD(4)) * pkin(4) + t652, t652 (-t571 * t641 - t734) * r_i_i_C(1) + (-t571 * t644 + t735) * r_i_i_C(2) + (t737 * r_i_i_C(1) + t738 * r_i_i_C(2)) * qJD(6); 0, 0 (-t602 * t641 + t607 * t704) * r_i_i_C(1) + (-t602 * t644 - t607 * t705) * r_i_i_C(2) + t602 * t647 + t659 * t603 + t648 * t606 (t602 * t642 + (-t607 * t645 - t618 * t642) * qJD(4)) * pkin(4) + t651, t651 (-t584 * t641 + t603 * t644) * r_i_i_C(1) + (-t584 * t644 - t603 * t641) * r_i_i_C(2) + ((-t593 * t644 - t606 * t641) * r_i_i_C(1) + (t593 * t641 - t606 * t644) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
