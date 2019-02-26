% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRRR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobiaD_transl_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:02
% EndTime: 2019-02-26 19:44:03
% DurationCPUTime: 1.18s
% Computational Cost: add. (3022->198), mult. (9505->345), div. (0->0), fcn. (11498->18), ass. (0->111)
t677 = sin(qJ(6));
t681 = cos(qJ(6));
t697 = qJD(6) * (r_i_i_C(1) * t677 + r_i_i_C(2) * t681);
t669 = sin(pkin(14));
t672 = sin(pkin(6));
t673 = cos(pkin(14));
t676 = cos(pkin(6));
t680 = sin(qJ(3));
t683 = cos(qJ(3));
t675 = cos(pkin(7));
t717 = t675 * t680;
t726 = sin(pkin(7));
t657 = (t669 * t683 + t673 * t717) * t672 + t676 * t726 * t680;
t727 = cos(pkin(13));
t707 = t727 * t673;
t670 = sin(pkin(13));
t722 = t670 * t676;
t664 = -t669 * t722 + t707;
t708 = t727 * t669;
t663 = -t673 * t722 - t708;
t710 = t672 * t726;
t691 = t663 * t675 + t670 * t710;
t648 = t664 * t683 + t680 * t691;
t730 = r_i_i_C(3) + pkin(12);
t729 = cos(qJ(4));
t671 = sin(pkin(8));
t728 = pkin(10) * t671;
t662 = t670 * t673 + t676 * t708;
t725 = t662 * t680;
t724 = t662 * t683;
t678 = sin(qJ(5));
t721 = t671 * t678;
t682 = cos(qJ(5));
t720 = t671 * t682;
t719 = t672 * t675;
t674 = cos(pkin(8));
t679 = sin(qJ(4));
t718 = t674 * t679;
t716 = t675 * t683;
t715 = qJD(3) * t680;
t714 = qJD(6) * t677;
t713 = qJD(6) * t681;
t712 = t674 * t729;
t711 = qJD(3) * t716;
t709 = t683 * t726;
t705 = qJD(3) * t709;
t661 = -t669 * t670 + t676 * t707;
t699 = t727 * t710;
t689 = -t661 * t675 + t699;
t645 = -t683 * t689 - t725;
t646 = t661 * t717 - t680 * t699 + t724;
t690 = -t661 * t726 - t719 * t727;
t686 = t690 * t671;
t619 = t646 * t729 + (t645 * t674 + t686) * t679;
t633 = -t645 * t671 + t674 * t690;
t607 = t619 * t682 + t633 * t678;
t703 = -t619 * t678 + t633 * t682;
t647 = -t664 * t680 + t683 * t691;
t693 = -t663 * t726 + t670 * t719;
t688 = t693 * t671;
t621 = t648 * t729 + (t647 * t674 + t688) * t679;
t634 = -t647 * t671 + t674 * t693;
t609 = t621 * t682 + t634 * t678;
t702 = -t621 * t678 + t634 * t682;
t656 = t676 * t709 + (-t669 * t680 + t673 * t716) * t672;
t692 = -t673 * t710 + t675 * t676;
t687 = t692 * t671;
t632 = t657 * t729 + (t656 * t674 + t687) * t679;
t649 = -t656 * t671 + t674 * t692;
t623 = t632 * t682 + t649 * t678;
t701 = -t632 * t678 + t649 * t682;
t700 = r_i_i_C(1) * t681 - r_i_i_C(2) * t677 + pkin(5);
t627 = t645 * t729 - t646 * t718;
t614 = t627 * t682 + t646 * t721;
t629 = t647 * t729 - t648 * t718;
t615 = t629 * t682 + t648 * t721;
t636 = t656 * t729 - t657 * t718;
t630 = t636 * t682 + t657 * t721;
t696 = -t645 * t679 - t646 * t712;
t695 = -t647 * t679 - t648 * t712;
t694 = -t656 * t679 - t657 * t712;
t685 = -t678 * t730 - t682 * t700 - pkin(4);
t620 = -t647 * t712 + t648 * t679 - t688 * t729;
t631 = -t656 * t712 + t657 * t679 - t687 * t729;
t618 = -t645 * t712 + t646 * t679 - t686 * t729;
t684 = t682 * t697 + (t678 * t700 - t682 * t730) * qJD(5);
t655 = t657 * qJD(3);
t654 = -t676 * t705 + (t669 * t715 - t673 * t711) * t672;
t644 = t648 * qJD(3);
t643 = -t670 * t672 * t705 - t663 * t711 + t664 * t715;
t642 = (t680 * t689 - t724) * qJD(3);
t641 = -t661 * t711 + (t683 * t699 + t725) * qJD(3);
t625 = qJD(4) * t694 + t654 * t718 - t655 * t729;
t624 = qJD(4) * t636 - t654 * t712 - t655 * t679;
t617 = -qJD(4) * t631 - t654 * t729 - t655 * t718;
t616 = qJD(4) * t632 - t654 * t679 + t655 * t712;
t613 = qJD(4) * t695 + t643 * t718 - t644 * t729;
t612 = qJD(4) * t629 - t643 * t712 - t644 * t679;
t611 = qJD(4) * t696 + t641 * t718 + t642 * t729;
t610 = qJD(4) * t627 - t641 * t712 + t642 * t679;
t605 = -qJD(4) * t620 - t643 * t729 - t644 * t718;
t604 = qJD(4) * t621 - t643 * t679 + t644 * t712;
t603 = -qJD(4) * t618 - t641 * t729 + t642 * t718;
t602 = qJD(4) * t619 - t641 * t679 - t642 * t712;
t601 = -t654 * t721 + t625 * t682 + (-t636 * t678 + t657 * t720) * qJD(5);
t599 = qJD(5) * t701 + t617 * t682 + t655 * t721;
t597 = -t643 * t721 + t613 * t682 + (-t629 * t678 + t648 * t720) * qJD(5);
t595 = -t641 * t721 + t611 * t682 + (-t627 * t678 + t646 * t720) * qJD(5);
t593 = qJD(5) * t702 + t605 * t682 + t644 * t721;
t591 = qJD(5) * t703 + t603 * t682 - t642 * t721;
t1 = [0, 0 (t597 * t681 + t612 * t677) * r_i_i_C(1) + (-t597 * t677 + t612 * t681) * r_i_i_C(2) + t597 * pkin(5) + t613 * pkin(4) + t612 * pkin(11) - t644 * pkin(3) - t643 * t728 + t730 * (qJD(5) * t615 + t613 * t678 + t643 * t720) + ((-t615 * t677 - t681 * t695) * r_i_i_C(1) + (-t615 * t681 + t677 * t695) * r_i_i_C(2)) * qJD(6) (t605 * t677 + t621 * t713) * r_i_i_C(1) + (t605 * t681 - t621 * t714) * r_i_i_C(2) + t605 * pkin(11) + t685 * t604 + t684 * t620, t730 * t593 - t702 * t697 + t700 * (-qJD(5) * t609 - t605 * t678 + t644 * t720) (-t593 * t677 + t604 * t681) * r_i_i_C(1) + (-t593 * t681 - t604 * t677) * r_i_i_C(2) + ((-t609 * t681 - t620 * t677) * r_i_i_C(1) + (t609 * t677 - t620 * t681) * r_i_i_C(2)) * qJD(6); 0, 0 (t595 * t681 + t610 * t677) * r_i_i_C(1) + (-t595 * t677 + t610 * t681) * r_i_i_C(2) + t595 * pkin(5) + t611 * pkin(4) + t610 * pkin(11) + t642 * pkin(3) - t641 * t728 + t730 * (qJD(5) * t614 + t611 * t678 + t641 * t720) + ((-t614 * t677 - t681 * t696) * r_i_i_C(1) + (-t614 * t681 + t677 * t696) * r_i_i_C(2)) * qJD(6) (t603 * t677 + t619 * t713) * r_i_i_C(1) + (t603 * t681 - t619 * t714) * r_i_i_C(2) + t603 * pkin(11) + t685 * t602 + t684 * t618, t730 * t591 - t703 * t697 + t700 * (-qJD(5) * t607 - t603 * t678 - t642 * t720) (-t591 * t677 + t602 * t681) * r_i_i_C(1) + (-t591 * t681 - t602 * t677) * r_i_i_C(2) + ((-t607 * t681 - t618 * t677) * r_i_i_C(1) + (t607 * t677 - t618 * t681) * r_i_i_C(2)) * qJD(6); 0, 0 (t601 * t681 + t624 * t677) * r_i_i_C(1) + (-t601 * t677 + t624 * t681) * r_i_i_C(2) + t601 * pkin(5) + t625 * pkin(4) + t624 * pkin(11) - t655 * pkin(3) - t654 * t728 + t730 * (qJD(5) * t630 + t625 * t678 + t654 * t720) + ((-t630 * t677 - t681 * t694) * r_i_i_C(1) + (-t630 * t681 + t677 * t694) * r_i_i_C(2)) * qJD(6) (t617 * t677 + t632 * t713) * r_i_i_C(1) + (t617 * t681 - t632 * t714) * r_i_i_C(2) + t617 * pkin(11) + t685 * t616 + t684 * t631, t730 * t599 - t701 * t697 + t700 * (-qJD(5) * t623 - t617 * t678 + t655 * t720) (-t599 * t677 + t616 * t681) * r_i_i_C(1) + (-t599 * t681 - t616 * t677) * r_i_i_C(2) + ((-t623 * t681 - t631 * t677) * r_i_i_C(1) + (t623 * t677 - t631 * t681) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
