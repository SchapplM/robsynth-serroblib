% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:20:40
% EndTime: 2019-02-26 20:20:42
% DurationCPUTime: 1.34s
% Computational Cost: add. (1934->206), mult. (4907->360), div. (0->0), fcn. (5425->16), ass. (0->123)
t693 = r_i_i_C(3) + pkin(12);
t622 = cos(qJ(6));
t673 = qJD(6) * t622;
t618 = sin(qJ(6));
t674 = qJD(6) * t618;
t694 = -r_i_i_C(1) * t674 - t673 * r_i_i_C(2);
t644 = r_i_i_C(1) * t622 - r_i_i_C(2) * t618 + pkin(5);
t621 = sin(qJ(2));
t624 = cos(qJ(2));
t616 = cos(pkin(13));
t689 = cos(pkin(6));
t660 = t616 * t689;
t688 = sin(pkin(13));
t597 = -t688 * t621 + t624 * t660;
t692 = cos(qJ(3));
t690 = pkin(4) * qJD(4);
t592 = t597 * qJD(2);
t614 = sin(pkin(7));
t687 = t592 * t614;
t648 = t689 * t688;
t633 = t616 * t621 + t624 * t648;
t594 = t633 * qJD(2);
t686 = t594 * t614;
t613 = qJ(4) + qJ(5);
t610 = sin(t613);
t612 = qJD(4) + qJD(5);
t685 = t610 * t612;
t684 = t610 * t614;
t683 = t612 * t614;
t615 = sin(pkin(6));
t682 = t614 * t615;
t619 = sin(qJ(4));
t681 = t614 * t619;
t623 = cos(qJ(4));
t680 = t614 * t623;
t679 = t614 * t624;
t678 = t615 * t616;
t617 = cos(pkin(7));
t620 = sin(qJ(3));
t677 = t617 * t620;
t676 = t620 * t621;
t675 = t620 * t624;
t671 = t619 * t690;
t670 = t614 * t678;
t669 = t621 * t682;
t668 = t615 * t676;
t635 = -t621 * t660 - t688 * t624;
t667 = t635 * t692;
t666 = t617 * t692;
t665 = t692 * t621;
t664 = t692 * t624;
t663 = qJD(2) * t682;
t662 = t614 * t689;
t661 = t615 * t688;
t593 = t635 * qJD(2);
t631 = t597 * t666 + t620 * t635 - t692 * t670;
t549 = t631 * qJD(3) + t592 * t692 + t593 * t677;
t582 = -t597 * t614 - t617 * t678;
t658 = t582 * t612 + t549;
t634 = -t616 * t624 + t621 * t648;
t595 = t634 * qJD(2);
t651 = t614 * t661;
t630 = t620 * t634 - t633 * t666 + t692 * t651;
t551 = t630 * qJD(3) - t594 * t692 + t595 * t677;
t583 = t614 * t633 + t617 * t661;
t657 = t583 * t612 + t551;
t639 = -t617 * t676 + t664;
t643 = t692 * t662;
t655 = t617 * t664;
t565 = qJD(3) * t643 + ((t655 - t676) * qJD(3) + t639 * qJD(2)) * t615;
t596 = -t615 * t679 + t689 * t617;
t656 = t596 * t612 + t565;
t654 = t621 * t663;
t652 = t620 * t662;
t641 = -t597 * t620 + t635 * t666;
t557 = t641 * qJD(3) - t592 * t677 + t593 * t692;
t650 = -t635 * t683 + t557;
t640 = t620 * t633 + t634 * t666;
t559 = t640 * qJD(3) + t594 * t677 + t595 * t692;
t649 = -t634 * t683 + t559;
t647 = t615 * t655;
t575 = t597 * t692 + t635 * t677;
t646 = t575 * t612 - t687;
t577 = -t633 * t692 + t634 * t677;
t645 = t577 * t612 + t686;
t637 = t617 * t665 + t675;
t638 = t617 * t675 + t665;
t573 = (-t638 * qJD(2) - t637 * qJD(3)) * t615;
t642 = t612 * t669 + t573;
t591 = t639 * t615;
t636 = -t591 * t612 + t624 * t663;
t609 = pkin(4) * t623 + pkin(3);
t611 = cos(t613);
t632 = -t693 * t610 - t644 * t611 - t609;
t571 = -t634 * t692 + (-t617 * t633 + t651) * t620;
t569 = -t667 + (t597 * t617 - t670) * t620;
t534 = -t569 * t685 - t593 * t684 + t658 * t611;
t629 = t694 * (-t569 * t610 + t582 * t611) + t693 * t534 + t644 * ((-t569 * t612 - t593 * t614) * t611 - t658 * t610);
t536 = -t571 * t685 - t595 * t684 + t657 * t611;
t628 = t694 * (-t571 * t610 + t583 * t611) + t693 * t536 + t644 * ((-t571 * t612 - t595 * t614) * t611 - t657 * t610);
t581 = t638 * t615 + t652;
t545 = -t581 * t685 + t610 * t654 + t656 * t611;
t627 = t694 * (-t581 * t610 + t596 * t611) + t693 * t545 + t644 * ((-t581 * t612 + t654) * t611 - t656 * t610);
t626 = t671 + (r_i_i_C(1) * t618 + r_i_i_C(2) * t622) * t611 * qJD(6) + (t644 * t610 - t693 * t611) * t612;
t625 = -pkin(11) - pkin(10);
t590 = t637 * t615;
t580 = -t643 - t647 + t668;
t578 = t591 * t611 + t610 * t669;
t572 = -qJD(2) * t647 - t615 * qJD(3) * t664 + (qJD(3) * t617 + qJD(2)) * t668;
t567 = t581 * t611 + t596 * t610;
t564 = qJD(3) * t652 + (t637 * qJD(2) + t638 * qJD(3)) * t615;
t561 = t577 * t611 - t634 * t684;
t560 = t575 * t611 - t635 * t684;
t558 = t577 * qJD(3) - t594 * t666 + t595 * t620;
t556 = t575 * qJD(3) + t592 * t666 + t593 * t620;
t555 = t571 * t611 + t583 * t610;
t553 = t569 * t611 + t582 * t610;
t550 = t571 * qJD(3) - t594 * t620 - t595 * t666;
t548 = t592 * t620 - t593 * t666 + (t597 * t677 - t620 * t670 - t667) * qJD(3);
t547 = t636 * t610 + t642 * t611;
t540 = -t645 * t610 + t649 * t611;
t538 = -t646 * t610 + t650 * t611;
t1 = [0 (t540 * t622 + t558 * t618) * r_i_i_C(1) + (-t540 * t618 + t558 * t622) * r_i_i_C(2) + t540 * pkin(5) + t559 * t609 - t558 * t625 + t595 * pkin(2) - pkin(9) * t686 + t693 * (t649 * t610 + t645 * t611) + ((-t561 * t618 - t622 * t640) * r_i_i_C(1) + (-t561 * t622 + t618 * t640) * r_i_i_C(2)) * qJD(6) + (-t594 * t681 + (-t577 * t619 - t634 * t680) * qJD(4)) * pkin(4) (t551 * t618 + t571 * t673) * r_i_i_C(1) + (t551 * t622 - t571 * t674) * r_i_i_C(2) - t551 * t625 + t632 * t550 - t626 * t630 (-t595 * t680 - t551 * t619 + (-t571 * t623 - t583 * t619) * qJD(4)) * pkin(4) + t628, t628 (-t536 * t618 + t550 * t622) * r_i_i_C(1) + (-t536 * t622 - t550 * t618) * r_i_i_C(2) + ((-t555 * t622 + t618 * t630) * r_i_i_C(1) + (t555 * t618 + t622 * t630) * r_i_i_C(2)) * qJD(6); 0 (t538 * t622 + t556 * t618) * r_i_i_C(1) + (-t538 * t618 + t556 * t622) * r_i_i_C(2) + t538 * pkin(5) + t557 * t609 - t556 * t625 + t593 * pkin(2) + pkin(9) * t687 + t693 * (t650 * t610 + t646 * t611) + ((-t560 * t618 - t622 * t641) * r_i_i_C(1) + (-t560 * t622 + t618 * t641) * r_i_i_C(2)) * qJD(6) + (t592 * t681 + (-t575 * t619 - t635 * t680) * qJD(4)) * pkin(4) (t549 * t618 + t569 * t673) * r_i_i_C(1) + (t549 * t622 - t569 * t674) * r_i_i_C(2) - t549 * t625 + t632 * t548 - t626 * t631 (-t593 * t680 - t549 * t619 + (-t569 * t623 - t582 * t619) * qJD(4)) * pkin(4) + t629, t629 (-t534 * t618 + t548 * t622) * r_i_i_C(1) + (-t534 * t622 - t548 * t618) * r_i_i_C(2) + ((-t553 * t622 + t618 * t631) * r_i_i_C(1) + (t553 * t618 + t622 * t631) * r_i_i_C(2)) * qJD(6); 0 (t547 * t622 - t572 * t618) * r_i_i_C(1) + (-t547 * t618 - t572 * t622) * r_i_i_C(2) + t547 * pkin(5) + t573 * t609 - t591 * t671 + t572 * t625 + t693 * (t642 * t610 - t636 * t611) + ((-t578 * t618 + t590 * t622) * r_i_i_C(1) + (-t578 * t622 - t590 * t618) * r_i_i_C(2)) * qJD(6) + (t621 * t680 * t690 + (-pkin(2) * t621 + (pkin(4) * t619 + pkin(9)) * t679) * qJD(2)) * t615 (t565 * t618 + t581 * t673) * r_i_i_C(1) + (t565 * t622 - t581 * t674) * r_i_i_C(2) - t565 * t625 + t632 * t564 + t626 * t580 (t623 * t654 - t565 * t619 + (-t581 * t623 - t596 * t619) * qJD(4)) * pkin(4) + t627, t627 (-t545 * t618 + t564 * t622) * r_i_i_C(1) + (-t545 * t622 - t564 * t618) * r_i_i_C(2) + ((-t567 * t622 - t580 * t618) * r_i_i_C(1) + (t567 * t618 - t580 * t622) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
