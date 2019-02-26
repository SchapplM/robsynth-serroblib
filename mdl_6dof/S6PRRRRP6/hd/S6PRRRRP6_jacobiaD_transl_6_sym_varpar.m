% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRP6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:08
% EndTime: 2019-02-26 20:18:09
% DurationCPUTime: 1.06s
% Computational Cost: add. (2052->202), mult. (6535->341), div. (0->0), fcn. (7368->14), ass. (0->121)
t639 = sin(qJ(2));
t642 = cos(qJ(2));
t634 = cos(pkin(12));
t705 = cos(pkin(6));
t682 = t634 * t705;
t704 = sin(pkin(12));
t620 = -t704 * t639 + t642 * t682;
t710 = r_i_i_C(1) + pkin(5);
t709 = r_i_i_C(2) + pkin(11);
t708 = cos(qJ(3));
t632 = sin(pkin(7));
t707 = pkin(9) * t632;
t706 = r_i_i_C(3) + qJ(6);
t637 = sin(qJ(4));
t703 = t632 * t637;
t641 = cos(qJ(4));
t702 = t632 * t641;
t701 = t632 * t642;
t633 = sin(pkin(6));
t700 = t633 * t634;
t635 = cos(pkin(7));
t638 = sin(qJ(3));
t699 = t635 * t638;
t636 = sin(qJ(5));
t698 = t636 * t641;
t697 = t638 * t639;
t696 = t638 * t642;
t695 = qJD(2) * t633;
t694 = qJD(4) * t637;
t693 = qJD(5) * t641;
t692 = t632 * t700;
t691 = t632 * t633 * t639;
t690 = t633 * t697;
t651 = -t639 * t682 - t704 * t642;
t689 = t651 * t708;
t688 = t635 * t708;
t687 = t708 * t639;
t686 = t708 * t642;
t685 = t632 * t695;
t684 = t632 * t705;
t683 = t633 * t704;
t680 = t635 * t686;
t679 = t639 * t685;
t678 = t642 * t685;
t676 = t638 * t684;
t675 = t632 * t683;
t615 = t620 * qJD(2);
t616 = t651 * qJD(2);
t644 = t620 * t688 + t638 * t651 - t708 * t692;
t574 = t644 * qJD(3) + t615 * t708 + t616 * t699;
t674 = -t644 * t693 + t574;
t671 = t705 * t704;
t649 = t634 * t639 + t642 * t671;
t617 = t649 * qJD(2);
t650 = -t634 * t642 + t639 * t671;
t618 = t650 * qJD(2);
t643 = t638 * t650 - t649 * t688 + t708 * t675;
t576 = t643 * qJD(3) - t617 * t708 + t618 * t699;
t673 = -t643 * t693 + t576;
t656 = -t635 * t697 + t686;
t660 = t708 * t684;
t591 = qJD(3) * t660 + ((t680 - t697) * qJD(3) + t656 * qJD(2)) * t633;
t670 = t633 * t680;
t605 = -t660 - t670 + t690;
t672 = t605 * t693 + t591;
t593 = -t689 + (t620 * t635 - t692) * t638;
t607 = -t620 * t632 - t635 * t700;
t584 = t593 * t641 + t607 * t637;
t640 = cos(qJ(5));
t669 = t584 * t640 - t636 * t644;
t595 = -t650 * t708 + (-t635 * t649 + t675) * t638;
t608 = t632 * t649 + t635 * t683;
t586 = t595 * t641 + t608 * t637;
t668 = t586 * t640 - t636 * t643;
t601 = t620 * t708 + t651 * t699;
t587 = t601 * t641 - t651 * t703;
t658 = -t620 * t638 + t651 * t688;
t667 = -t587 * t636 - t640 * t658;
t603 = -t649 * t708 + t650 * t699;
t588 = t603 * t641 - t650 * t703;
t657 = t638 * t649 + t650 * t688;
t666 = -t588 * t636 - t640 * t657;
t665 = -t593 * t637 + t607 * t641;
t664 = -t595 * t637 + t608 * t641;
t655 = t635 * t696 + t687;
t606 = t655 * t633 + t676;
t619 = -t633 * t701 + t705 * t635;
t597 = t606 * t641 + t619 * t637;
t663 = t597 * t640 + t605 * t636;
t614 = t656 * t633;
t604 = t614 * t641 + t637 * t691;
t654 = t635 * t687 + t696;
t613 = t654 * t633;
t662 = -t604 * t636 + t613 * t640;
t661 = -t606 * t637 + t619 * t641;
t659 = -t641 * pkin(4) - t709 * t637 - pkin(3);
t653 = qJD(4) * (pkin(4) * t637 - t709 * t641);
t652 = t706 * t636 + t710 * t640 + pkin(4);
t573 = t615 * t638 - t616 * t688 + (t620 * t699 - t638 * t692 - t689) * qJD(3);
t648 = qJD(5) * t593 - t573 * t641 - t644 * t694;
t575 = t595 * qJD(3) - t617 * t638 - t618 * t688;
t647 = qJD(5) * t595 - t575 * t641 - t643 * t694;
t590 = qJD(3) * t676 + (t654 * qJD(2) + t655 * qJD(3)) * t633;
t646 = qJD(5) * t606 - t590 * t641 + t605 * t694;
t645 = qJD(6) * t636 + (-t710 * t636 + t706 * t640) * qJD(5);
t599 = (-t655 * qJD(2) - t654 * qJD(3)) * t633;
t598 = -qJD(2) * t670 - t633 * qJD(3) * t686 + (qJD(3) * t635 + qJD(2)) * t690;
t582 = t657 * qJD(3) + t617 * t699 + t618 * t708;
t581 = t603 * qJD(3) - t617 * t688 + t618 * t638;
t580 = t658 * qJD(3) - t615 * t699 + t616 * t708;
t579 = t601 * qJD(3) + t615 * t688 + t616 * t638;
t578 = t637 * t678 + t599 * t641 + (-t614 * t637 + t641 * t691) * qJD(4);
t570 = t661 * qJD(4) + t591 * t641 + t637 * t679;
t568 = -t617 * t703 + t582 * t641 + (-t603 * t637 - t650 * t702) * qJD(4);
t566 = t615 * t703 + t580 * t641 + (-t601 * t637 - t651 * t702) * qJD(4);
t562 = t664 * qJD(4) + t576 * t641 - t618 * t703;
t560 = t665 * qJD(4) + t574 * t641 - t616 * t703;
t555 = t663 * qJD(5) + t570 * t636 - t590 * t640;
t545 = t668 * qJD(5) + t562 * t636 - t575 * t640;
t543 = t669 * qJD(5) + t560 * t636 - t573 * t640;
t1 = [0, -t666 * qJD(6) + t568 * pkin(4) + t582 * pkin(3) + t581 * pkin(10) + t618 * pkin(2) - t617 * t707 + t709 * (t588 * qJD(4) + t582 * t637 + t617 * t702) + t710 * (t666 * qJD(5) + t568 * t640 + t581 * t636) + t706 * (t568 * t636 - t581 * t640 + (t588 * t640 - t636 * t657) * qJD(5)) -(t595 * t640 - t643 * t698) * qJD(6) + t576 * pkin(10) + t710 * (t673 * t636 + t647 * t640) + t706 * (t647 * t636 - t673 * t640) - t643 * t653 + t659 * t575, t709 * t562 + t645 * t664 + t652 * (-t586 * qJD(4) - t576 * t637 - t618 * t702) t668 * qJD(6) + t706 * (t562 * t640 + t575 * t636 + (-t586 * t636 - t640 * t643) * qJD(5)) - t710 * t545, t545; 0, -t667 * qJD(6) + t566 * pkin(4) + t580 * pkin(3) + t579 * pkin(10) + t616 * pkin(2) + t615 * t707 + t709 * (t587 * qJD(4) + t580 * t637 - t615 * t702) + t710 * (t667 * qJD(5) + t566 * t640 + t579 * t636) + t706 * (t566 * t636 - t579 * t640 + (t587 * t640 - t636 * t658) * qJD(5)) -(t593 * t640 - t644 * t698) * qJD(6) + t574 * pkin(10) + t710 * (t674 * t636 + t648 * t640) + t706 * (t648 * t636 - t674 * t640) - t644 * t653 + t659 * t573, t709 * t560 + t645 * t665 + t652 * (-t584 * qJD(4) - t574 * t637 - t616 * t702) t669 * qJD(6) + t706 * (t560 * t640 + t573 * t636 + (-t584 * t636 - t640 * t644) * qJD(5)) - t710 * t543, t543; 0, -t662 * qJD(6) + t578 * pkin(4) + t599 * pkin(3) - t598 * pkin(10) + t709 * (t604 * qJD(4) + t599 * t637 - t641 * t678) + t710 * (t662 * qJD(5) + t578 * t640 - t598 * t636) + t706 * (t578 * t636 + t598 * t640 + (t604 * t640 + t613 * t636) * qJD(5)) + (-pkin(2) * t639 + pkin(9) * t701) * t695 -(t605 * t698 + t606 * t640) * qJD(6) + t591 * pkin(10) + t710 * (t672 * t636 + t646 * t640) + t706 * (t646 * t636 - t672 * t640) + t605 * t653 + t659 * t590, t709 * t570 + t645 * t661 + t652 * (-t597 * qJD(4) - t591 * t637 + t641 * t679) t663 * qJD(6) + t706 * (t570 * t640 + t590 * t636 + (-t597 * t636 + t605 * t640) * qJD(5)) - t710 * t555, t555;];
JaD_transl  = t1;
