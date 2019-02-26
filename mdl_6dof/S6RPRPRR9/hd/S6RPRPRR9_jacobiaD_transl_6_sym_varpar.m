% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:29
% EndTime: 2019-02-26 20:53:30
% DurationCPUTime: 1.48s
% Computational Cost: add. (1858->170), mult. (5833->294), div. (0->0), fcn. (6652->16), ass. (0->106)
t631 = sin(pkin(7));
t629 = sin(pkin(13));
t642 = cos(qJ(3));
t690 = cos(pkin(13));
t666 = qJD(3) * t690;
t638 = sin(qJ(3));
t677 = qJD(3) * t638;
t693 = t629 * t677 - t642 * t666;
t596 = t693 * t631;
t620 = -t642 * t629 - t638 * t690;
t601 = t620 * t631;
t632 = sin(pkin(6));
t643 = cos(qJ(1));
t634 = cos(pkin(7));
t598 = t693 * t634;
t603 = t620 * t634;
t635 = cos(pkin(6));
t630 = sin(pkin(12));
t682 = t643 * t630;
t633 = cos(pkin(12));
t639 = sin(qJ(1));
t683 = t639 * t633;
t659 = t635 * t683 + t682;
t606 = t659 * qJD(1);
t684 = t639 * t630;
t670 = t635 * t684;
t678 = qJD(1) * t643;
t607 = -qJD(1) * t670 + t633 * t678;
t681 = t643 * t633;
t611 = -t635 * t681 + t684;
t612 = t635 * t682 + t683;
t676 = qJD(3) * t642;
t616 = -t629 * t676 - t638 * t666;
t653 = -t638 * t629 + t642 * t690;
t647 = t611 * t598 + t606 * t603 + t607 * t653 + t612 * t616;
t679 = qJD(1) * t639;
t562 = t632 * (t596 * t643 - t601 * t679) + t647;
t669 = t632 * t679;
t590 = t606 * t631 + t634 * t669;
t637 = sin(qJ(5));
t641 = cos(qJ(5));
t687 = t632 * t643;
t578 = -t601 * t687 - t611 * t603 - t612 * t653;
t591 = -t611 * t631 + t634 * t687;
t664 = t578 * t637 - t591 * t641;
t554 = qJD(5) * t664 + t562 * t641 + t590 * t637;
t651 = qJD(3) * t620;
t597 = t631 * t651;
t599 = t634 * t651;
t600 = t653 * t631;
t602 = t653 * t634;
t615 = t653 * qJD(3);
t561 = -t611 * t599 - t606 * t602 + t607 * t620 - t612 * t615 + t632 * (-t597 * t643 + t600 * t679);
t636 = sin(qJ(6));
t640 = cos(qJ(6));
t706 = t554 * t636 + t561 * t640;
t705 = -t554 * t640 + t561 * t636;
t569 = t578 * t641 + t591 * t637;
t658 = -t600 * t687 - t611 * t602 + t612 * t620;
t704 = -t569 * t636 + t640 * t658;
t703 = t569 * t640 + t636 * t658;
t702 = qJD(5) * t569 - t562 * t637 + t590 * t641;
t604 = t611 * qJD(1);
t605 = t612 * qJD(1);
t614 = -t670 + t681;
t559 = t659 * t598 - t604 * t603 - t605 * t653 + t614 * t616 + t632 * (-t596 * t639 - t601 * t678);
t584 = t635 * t596 + (t598 * t633 - t616 * t630) * t632;
t692 = pkin(11) + r_i_i_C(3);
t691 = pkin(9) + qJ(4);
t689 = t631 * t638;
t688 = t632 * t639;
t686 = t634 * t638;
t680 = pkin(3) * t689 + t634 * t691 + qJ(2);
t675 = qJD(6) * t636;
t674 = qJD(6) * t640;
t671 = pkin(3) * t676;
t673 = t634 * qJD(4) + t631 * t671 + qJD(2);
t672 = pkin(3) * t677;
t668 = t632 * t678;
t593 = t631 * t659 + t634 * t688;
t650 = -t601 * t688 + t603 * t659 + t614 * t653;
t571 = t593 * t637 + t641 * t650;
t663 = t593 * t641 - t637 * t650;
t610 = -t632 * t633 * t631 + t635 * t634;
t649 = -t635 * t601 + (-t603 * t633 + t630 * t653) * t632;
t573 = t610 * t637 + t641 * t649;
t662 = t610 * t641 - t637 * t649;
t660 = t640 * r_i_i_C(1) - t636 * r_i_i_C(2) + pkin(5);
t657 = qJD(6) * (-t636 * r_i_i_C(1) - t640 * r_i_i_C(2));
t652 = -t604 * t631 + t634 * t668;
t645 = t637 * t692 + t660 * t641 + pkin(4);
t644 = t641 * t657 + (-t660 * t637 + t641 * t692) * qJD(5);
t628 = t642 * pkin(3) + pkin(2);
t618 = -t631 * qJD(4) + t634 * t671;
t609 = pkin(3) * t686 - t631 * t691;
t586 = t635 * t600 + (t602 * t633 + t620 * t630) * t632;
t582 = t635 * t597 + (t599 * t633 - t615 * t630) * t632;
t580 = t600 * t688 - t602 * t659 + t614 * t620;
t566 = qJD(5) * t662 - t584 * t641;
t560 = -t596 * t687 + t601 * t669 - t647;
t558 = -t659 * t599 + t604 * t602 - t605 * t620 - t614 * t615 + (t597 * t639 + t600 * t678) * t632;
t552 = qJD(5) * t663 + t559 * t641 + t637 * t652;
t551 = qJD(5) * t571 + t559 * t637 - t641 * t652;
t550 = t552 * t640 - t558 * t636 + (-t571 * t636 - t580 * t640) * qJD(6);
t549 = -t552 * t636 - t558 * t640 + (-t571 * t640 + t580 * t636) * qJD(6);
t1 = [t705 * r_i_i_C(1) + t706 * r_i_i_C(2) - t554 * pkin(5) - t562 * pkin(4) + t561 * pkin(10) - t607 * t628 + t612 * t672 + t606 * t609 + t611 * t618 - pkin(1) * t678 + t692 * t702 + (t704 * r_i_i_C(1) - t703 * r_i_i_C(2)) * qJD(6) + (t643 * t673 - t679 * t680) * t632, t668 (t559 * t636 + t650 * t674) * r_i_i_C(1) + (t559 * t640 - t650 * t675) * r_i_i_C(2) + t559 * pkin(10) + t645 * t558 + t644 * t580 + (t605 * t638 + (t604 * t634 + t631 * t668) * t642 + (-t614 * t642 + (-t631 * t688 + t634 * t659) * t638) * qJD(3)) * pkin(3), t652, -t660 * t551 + t552 * t692 + t663 * t657, t549 * r_i_i_C(1) - t550 * r_i_i_C(2); -t614 * t672 - pkin(1) * t679 + t559 * pkin(4) + t552 * pkin(5) - t558 * pkin(10) + t550 * r_i_i_C(1) + t549 * r_i_i_C(2) + t604 * t609 - t605 * t628 - t659 * t618 + t692 * t551 + (t639 * t673 + t678 * t680) * t632, t669 (-t560 * t636 - t578 * t674) * r_i_i_C(1) + (-t560 * t640 + t578 * t675) * r_i_i_C(2) - t560 * pkin(10) + t645 * t561 + t644 * t658 + (-t607 * t638 + (-t606 * t634 + t631 * t669) * t642 + (-t612 * t642 + (t611 * t634 + t631 * t687) * t638) * qJD(3)) * pkin(3), t590, t692 * t554 + t664 * t657 + t660 * t702, -t706 * r_i_i_C(1) + t705 * r_i_i_C(2) + (t703 * r_i_i_C(1) + t704 * r_i_i_C(2)) * qJD(6); 0, 0 (-t584 * t636 + t649 * t674) * r_i_i_C(1) + (-t584 * t640 - t649 * t675) * r_i_i_C(2) - t584 * pkin(10) + (-t635 * t689 + (-t630 * t642 - t633 * t686) * t632) * qJD(3) * pkin(3) + t645 * t582 + t644 * t586, 0, t692 * t566 + t662 * t657 + t660 * (-qJD(5) * t573 + t584 * t637) (-t566 * t636 - t582 * t640) * r_i_i_C(1) + (-t566 * t640 + t582 * t636) * r_i_i_C(2) + ((-t573 * t640 + t586 * t636) * r_i_i_C(1) + (t573 * t636 + t586 * t640) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
