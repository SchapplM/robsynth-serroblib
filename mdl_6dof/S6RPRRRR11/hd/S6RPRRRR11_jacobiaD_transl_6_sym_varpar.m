% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR11
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

function JaD_transl = S6RPRRRR11_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR11_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:20:30
% EndTime: 2019-02-26 21:20:32
% DurationCPUTime: 1.37s
% Computational Cost: add. (1789->149), mult. (5163->246), div. (0->0), fcn. (5766->16), ass. (0->99)
t623 = cos(qJ(1));
t694 = cos(pkin(13));
t696 = cos(pkin(6));
t657 = t696 * t694;
t692 = sin(pkin(13));
t700 = sin(qJ(1));
t636 = t623 * t692 + t700 * t657;
t593 = t636 * qJD(1);
t655 = t696 * t692;
t598 = t623 * t694 - t700 * t655;
t594 = t598 * qJD(1);
t620 = sin(qJ(3));
t634 = -t623 * t655 - t700 * t694;
t617 = sin(pkin(6));
t693 = sin(pkin(7));
t672 = t693 * t617;
t648 = t700 * t672;
t695 = cos(pkin(7));
t701 = cos(qJ(3));
t635 = t623 * t657 - t700 * t692;
t661 = t623 * t672;
t646 = t701 * t661;
t660 = t695 * t701;
t706 = t635 * t660 - t646;
t560 = (qJD(1) * t648 + qJD(3) * t634 - t695 * t593) * t620 + t594 * t701 + t706 * qJD(3);
t659 = t695 * t700;
t649 = t617 * t659;
t676 = t593 * t693;
t582 = qJD(1) * t649 + t676;
t619 = sin(qJ(4));
t622 = cos(qJ(4));
t673 = t695 * t635;
t679 = t634 * t701;
t575 = (t661 - t673) * t620 + t679;
t686 = t617 * t623;
t705 = t635 * t693 + t695 * t686;
t653 = t575 * t619 - t622 * t705;
t554 = -t653 * qJD(4) - t560 * t622 - t582 * t619;
t572 = -t620 * t634 - t706;
t615 = qJD(5) + qJD(6);
t717 = -t572 * t615 + t554;
t642 = t701 * t648;
t651 = t620 * t661;
t561 = qJD(1) * t642 - t593 * t660 - t594 * t620 + (-t620 * t673 + t651 + t679) * qJD(3);
t567 = t575 * t622 + t619 * t705;
t716 = -t567 * t615 + t561;
t715 = t567 * qJD(4) - t560 * t619 + t582 * t622;
t698 = r_i_i_C(3) + pkin(12) + pkin(11);
t630 = t636 * t695;
t577 = t598 * t701 + (-t630 + t648) * t620;
t587 = t636 * t693 + t649;
t704 = -t577 * t619 + t587 * t622;
t703 = -t598 * t620 - t701 * t630 + t642;
t654 = t695 * t694;
t656 = t696 * t693;
t583 = -t701 * t656 + (t620 * t692 - t654 * t701) * t617;
t616 = qJ(5) + qJ(6);
t613 = sin(t616);
t614 = cos(t616);
t618 = sin(qJ(5));
t699 = t618 * pkin(5);
t680 = qJD(5) * t699;
t702 = (t613 * r_i_i_C(1) + t614 * r_i_i_C(2)) * t615 + t680;
t626 = t705 * qJD(1);
t688 = t613 * t615;
t687 = t614 * t615;
t592 = t634 * qJD(1);
t628 = qJD(1) * t673;
t557 = -qJD(1) * t646 + t577 * qJD(3) + t592 * t620 + t701 * t628;
t569 = t577 * t622 + t587 * t619;
t668 = -t569 * t615 + t557;
t558 = qJD(1) * t651 + qJD(3) * t703 + t592 * t701 - t620 * t628;
t550 = qJD(4) * t704 + t558 * t622 + t619 * t626;
t671 = -t615 * t703 + t550;
t547 = -t671 * t613 + t668 * t614;
t548 = t668 * t613 + t671 * t614;
t685 = t547 * r_i_i_C(1) - t548 * r_i_i_C(2);
t684 = (t613 * t717 - t614 * t716) * r_i_i_C(1) + (t613 * t716 + t614 * t717) * r_i_i_C(2);
t584 = t620 * t656 + (t620 * t654 + t701 * t692) * t617;
t595 = -t694 * t672 + t696 * t695;
t571 = t584 * t622 + t595 * t619;
t580 = t584 * qJD(3);
t664 = t571 * t615 - t580;
t579 = t583 * qJD(3);
t652 = -t584 * t619 + t595 * t622;
t564 = t652 * qJD(4) - t579 * t622;
t665 = -t583 * t615 - t564;
t683 = (t665 * t613 - t664 * t614) * r_i_i_C(1) + (t664 * t613 + t665 * t614) * r_i_i_C(2);
t621 = cos(qJ(5));
t682 = qJD(5) * t621;
t681 = t617 * qJD(2);
t678 = qJD(1) * t686;
t677 = qJD(1) * t700;
t612 = t621 * pkin(5) + pkin(4);
t647 = r_i_i_C(1) * t614 - r_i_i_C(2) * t613 + t612;
t632 = -t698 * t619 - t647 * t622 - pkin(3);
t625 = t702 * t622 + (t647 * t619 - t698 * t622) * qJD(4);
t549 = t569 * qJD(4) + t558 * t619 - t622 * t626;
t1 = [-pkin(9) * t676 + t623 * t681 - t594 * pkin(2) - t560 * pkin(3) + t561 * pkin(10) + t554 * t612 + (r_i_i_C(1) * t717 + r_i_i_C(2) * t716) * t614 + (r_i_i_C(1) * t716 - r_i_i_C(2) * t717) * t613 + t698 * t715 + (-t623 * pkin(1) + (-pkin(9) * t659 - t700 * qJ(2)) * t617) * qJD(1) + (t561 * t618 + (-t567 * t618 - t572 * t621) * qJD(5)) * pkin(5), t678 (t558 * t613 + t577 * t687) * r_i_i_C(1) + (t558 * t614 - t577 * t688) * r_i_i_C(2) + t558 * pkin(10) + (t558 * t618 + t577 * t682) * pkin(5) + t632 * t557 - t625 * t703, -t647 * t549 + t698 * t550 - t702 * t704 (-t550 * t618 + t557 * t621 + (-t569 * t621 + t618 * t703) * qJD(5)) * pkin(5) + t685, t685; -t703 * pkin(5) * t682 - pkin(1) * t677 + t592 * pkin(2) + t558 * pkin(3) + t548 * r_i_i_C(1) + t547 * r_i_i_C(2) + qJ(2) * t678 + t550 * t612 - t569 * t680 + t700 * t681 + pkin(9) * t626 + (pkin(10) + t699) * t557 + t698 * t549, t617 * t677 (t560 * t613 - t575 * t687) * r_i_i_C(1) + (t560 * t614 + t575 * t688) * r_i_i_C(2) + t560 * pkin(10) + (t560 * t618 - t575 * t682) * pkin(5) - t632 * t561 + t625 * t572, -t554 * t698 + t647 * t715 - t702 * t653 (t554 * t618 - t561 * t621 + (t567 * t621 - t572 * t618) * qJD(5)) * pkin(5) + t684, t684; 0, 0 (-t579 * t613 + t584 * t687) * r_i_i_C(1) + (-t579 * t614 - t584 * t688) * r_i_i_C(2) - t579 * pkin(10) + (-t579 * t618 + t584 * t682) * pkin(5) + t632 * t580 + t625 * t583, t698 * t564 - t702 * t652 + t647 * (-t571 * qJD(4) + t579 * t619) (-t564 * t618 + t580 * t621 + (-t571 * t621 - t583 * t618) * qJD(5)) * pkin(5) + t683, t683;];
JaD_transl  = t1;
