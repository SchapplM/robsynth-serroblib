% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:18
% EndTime: 2019-02-26 20:05:19
% DurationCPUTime: 1.22s
% Computational Cost: add. (1830->227), mult. (5726->394), div. (0->0), fcn. (6521->16), ass. (0->121)
t606 = sin(pkin(13));
t667 = cos(pkin(13));
t669 = cos(qJ(3));
t639 = t669 * t667;
t615 = sin(qJ(3));
t654 = qJD(3) * t615;
t671 = -qJD(3) * t639 + t606 * t654;
t670 = r_i_i_C(3) + pkin(11);
t611 = cos(pkin(7));
t668 = pkin(3) * t611;
t642 = t615 * t667;
t643 = t669 * qJD(3);
t594 = -qJD(3) * t642 - t606 * t643;
t616 = sin(qJ(2));
t666 = t594 * t616;
t607 = sin(pkin(12));
t609 = sin(pkin(6));
t665 = t607 * t609;
t608 = sin(pkin(7));
t664 = t608 * t609;
t614 = sin(qJ(5));
t663 = t608 * t614;
t618 = cos(qJ(5));
t662 = t608 * t618;
t610 = cos(pkin(12));
t661 = t609 * t610;
t660 = t609 * t611;
t575 = t671 * t608;
t612 = cos(pkin(6));
t659 = t612 * t575;
t658 = t612 * t616;
t619 = cos(qJ(2));
t657 = t612 * t619;
t656 = t615 * t619;
t655 = qJD(2) * t616;
t613 = sin(qJ(6));
t653 = qJD(6) * t613;
t617 = cos(qJ(6));
t652 = qJD(6) * t617;
t651 = pkin(3) * t654;
t650 = t616 * t664;
t649 = t619 * t664;
t648 = t610 * t657;
t647 = t611 * t669;
t646 = t669 * t616;
t645 = t609 * t655;
t641 = qJD(2) * t649;
t640 = t608 * t645;
t589 = -t607 * t616 + t648;
t568 = -t589 * t608 - t610 * t660;
t627 = t606 * t669 + t642;
t580 = t627 * t608;
t582 = t627 * t611;
t590 = t607 * t619 + t610 * t658;
t626 = -t615 * t606 + t639;
t625 = -t580 * t661 + t582 * t589 + t590 * t626;
t540 = t568 * t614 + t618 * t625;
t638 = t568 * t618 - t614 * t625;
t630 = t607 * t657 + t610 * t616;
t569 = t607 * t660 + t608 * t630;
t629 = t607 * t658 - t610 * t619;
t624 = t580 * t665 - t582 * t630 - t626 * t629;
t542 = t569 * t614 + t618 * t624;
t637 = t569 * t618 - t614 * t624;
t588 = t612 * t611 - t649;
t634 = t582 * t619 + t616 * t626;
t622 = t612 * t580 + t609 * t634;
t557 = t588 * t614 + t618 * t622;
t636 = t588 * t618 - t614 * t622;
t581 = t626 * t611;
t635 = t581 * t619 - t616 * t627;
t633 = -t582 * t616 + t619 * t626;
t632 = t617 * r_i_i_C(1) - t613 * r_i_i_C(2) + pkin(5);
t562 = -t582 * t590 + t589 * t626;
t548 = t562 * t618 + t590 * t663;
t564 = t582 * t629 - t626 * t630;
t549 = t564 * t618 - t629 * t663;
t628 = qJD(6) * (-t613 * r_i_i_C(1) - t617 * r_i_i_C(2));
t567 = t633 * t609;
t565 = t567 * t618 + t614 * t650;
t623 = qJD(3) * t627;
t621 = t614 * t670 + t632 * t618 + pkin(4);
t577 = t671 * t611;
t583 = -qJD(2) * t648 + t607 * t655;
t584 = t590 * qJD(2);
t531 = t575 * t661 - t577 * t589 - t582 * t584 - t583 * t626 + t590 * t594;
t585 = t630 * qJD(2);
t586 = t629 * qJD(2);
t532 = t575 * t665 - t577 * t630 - t582 * t586 + t585 * t626 + t594 * t629;
t620 = t618 * t628 + (-t632 * t614 + t618 * t670) * qJD(5);
t605 = pkin(3) * t669 + pkin(2);
t595 = -t608 * qJD(4) + t643 * t668;
t593 = t626 * qJD(3);
t587 = t615 * t668 + (-pkin(9) - qJ(4)) * t608;
t579 = t626 * t608;
t578 = t611 * t623;
t576 = t608 * t623;
t566 = (t581 * t616 + t619 * t627) * t609;
t563 = -t581 * t629 - t627 * t630;
t561 = t581 * t590 + t589 * t627;
t559 = t612 * t579 + t609 * t635;
t554 = t579 * t665 - t581 * t630 + t627 * t629;
t551 = -t579 * t661 + t581 * t589 - t590 * t627;
t547 = (-qJD(2) * t634 + t577 * t616 + t594 * t619) * t609;
t546 = (qJD(2) * t635 - t578 * t616 + t593 * t619) * t609;
t545 = -t659 + (qJD(2) * t633 - t577 * t619 + t666) * t609;
t544 = -t612 * t576 - t581 * t645 + (-t593 * t616 + (-qJD(2) * t627 - t578) * t619) * t609;
t543 = t659 + t582 * t645 + (-t666 + (-qJD(2) * t626 + t577) * t619) * t609;
t538 = -t577 * t629 + t582 * t585 + t586 * t626 - t594 * t630;
t537 = -t578 * t629 + t581 * t585 - t586 * t627 + t593 * t630;
t536 = t577 * t590 + t582 * t583 - t584 * t626 + t589 * t594;
t535 = t578 * t590 + t581 * t583 + t584 * t627 - t589 * t593;
t533 = -t576 * t665 + t578 * t630 + t581 * t586 + t585 * t627 + t593 * t629;
t530 = t576 * t661 - t578 * t589 - t581 * t584 + t583 * t627 - t590 * t593;
t528 = t614 * t641 + t547 * t618 + (-t567 * t614 + t618 * t650) * qJD(5);
t526 = qJD(5) * t636 + t545 * t618 + t614 * t640;
t524 = -t585 * t663 + t538 * t618 + (-t564 * t614 - t629 * t662) * qJD(5);
t522 = -t583 * t663 + t536 * t618 + (-t562 * t614 + t590 * t662) * qJD(5);
t520 = qJD(5) * t637 - t532 * t618 - t586 * t663;
t518 = qJD(5) * t638 + t531 * t618 + t584 * t663;
t1 = [0 (t524 * t617 - t537 * t613) * r_i_i_C(1) + (-t524 * t613 - t537 * t617) * r_i_i_C(2) + t524 * pkin(5) + t538 * pkin(4) - t537 * pkin(10) + t586 * t605 + t630 * t651 + t585 * t587 + t629 * t595 + t670 * (qJD(5) * t549 + t538 * t614 + t585 * t662) + ((-t549 * t613 + t563 * t617) * r_i_i_C(1) + (-t549 * t617 - t563 * t613) * r_i_i_C(2)) * qJD(6) (-t532 * t613 + t624 * t652) * r_i_i_C(1) + (-t532 * t617 - t624 * t653) * r_i_i_C(2) - t532 * pkin(10) + t621 * t533 + t620 * t554 + (t586 * t647 + t585 * t615 + (t669 * t629 + (-t607 * t664 + t611 * t630) * t615) * qJD(3)) * pkin(3), -t586 * t608, t670 * t520 + t637 * t628 + t632 * (-qJD(5) * t542 + t532 * t614 - t586 * t662) (-t520 * t613 - t533 * t617) * r_i_i_C(1) + (-t520 * t617 + t533 * t613) * r_i_i_C(2) + ((-t542 * t617 + t554 * t613) * r_i_i_C(1) + (t542 * t613 + t554 * t617) * r_i_i_C(2)) * qJD(6); 0 (t522 * t617 - t535 * t613) * r_i_i_C(1) + (-t522 * t613 - t535 * t617) * r_i_i_C(2) + t522 * pkin(5) + t536 * pkin(4) - t535 * pkin(10) - t584 * t605 - t589 * t651 + t583 * t587 - t590 * t595 + t670 * (qJD(5) * t548 + t536 * t614 + t583 * t662) + ((-t548 * t613 + t561 * t617) * r_i_i_C(1) + (-t548 * t617 - t561 * t613) * r_i_i_C(2)) * qJD(6) (t531 * t613 + t625 * t652) * r_i_i_C(1) + (t531 * t617 - t625 * t653) * r_i_i_C(2) + t531 * pkin(10) + t621 * t530 + t620 * t551 + (-t584 * t647 + t583 * t615 + (-t669 * t590 + (-t589 * t611 + t608 * t661) * t615) * qJD(3)) * pkin(3), t584 * t608, t670 * t518 + t638 * t628 + t632 * (-qJD(5) * t540 - t531 * t614 + t584 * t662) (-t518 * t613 - t530 * t617) * r_i_i_C(1) + (-t518 * t617 + t530 * t613) * r_i_i_C(2) + ((-t540 * t617 + t551 * t613) * r_i_i_C(1) + (t540 * t613 + t551 * t617) * r_i_i_C(2)) * qJD(6); 0 (t528 * t617 + t546 * t613) * r_i_i_C(1) + (-t528 * t613 + t546 * t617) * r_i_i_C(2) + t528 * pkin(5) + t547 * pkin(4) + t546 * pkin(10) + t670 * (qJD(5) * t565 + t547 * t614 - t618 * t641) + ((-t565 * t613 + t566 * t617) * r_i_i_C(1) + (-t565 * t617 - t566 * t613) * r_i_i_C(2)) * qJD(6) + (-t619 * t651 - t595 * t616 + (-t587 * t619 - t605 * t616) * qJD(2)) * t609 (-t543 * t613 + t622 * t652) * r_i_i_C(1) + (-t543 * t617 - t622 * t653) * r_i_i_C(2) - t543 * pkin(10) + t621 * t544 + t620 * t559 + (-t612 * t608 * t654 + ((-t611 * t656 - t646) * qJD(3) + (-t611 * t646 - t656) * qJD(2)) * t609) * pkin(3), t640, t670 * t526 + t636 * t628 + t632 * (-qJD(5) * t557 - t545 * t614 + t618 * t640) (-t526 * t613 - t544 * t617) * r_i_i_C(1) + (-t526 * t617 + t544 * t613) * r_i_i_C(2) + ((-t557 * t617 + t559 * t613) * r_i_i_C(1) + (t557 * t613 + t559 * t617) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
