% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR5
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
% Datum: 2019-02-26 20:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:14
% EndTime: 2019-02-26 20:21:15
% DurationCPUTime: 1.23s
% Computational Cost: add. (1773->198), mult. (5101->335), div. (0->0), fcn. (5674->16), ass. (0->124)
t596 = sin(qJ(2));
t599 = cos(qJ(2));
t591 = cos(pkin(13));
t668 = cos(pkin(6));
t642 = t591 * t668;
t667 = sin(pkin(13));
t572 = -t667 * t596 + t599 * t642;
t588 = qJ(5) + qJ(6);
t585 = sin(t588);
t586 = cos(t588);
t587 = qJD(5) + qJD(6);
t593 = sin(qJ(5));
t672 = qJD(5) * t593 * pkin(5) + (t585 * r_i_i_C(1) + t586 * r_i_i_C(2)) * t587;
t671 = cos(qJ(3));
t589 = sin(pkin(7));
t670 = pkin(9) * t589;
t669 = r_i_i_C(3) + pkin(12) + pkin(11);
t666 = t585 * t587;
t665 = t586 * t587;
t594 = sin(qJ(4));
t664 = t589 * t594;
t598 = cos(qJ(4));
t663 = t589 * t598;
t662 = t589 * t599;
t590 = sin(pkin(6));
t661 = t590 * t591;
t592 = cos(pkin(7));
t595 = sin(qJ(3));
t660 = t592 * t595;
t659 = t595 * t596;
t658 = t595 * t599;
t567 = t572 * qJD(2);
t608 = -t596 * t642 - t667 * t599;
t568 = t608 * qJD(2);
t647 = t592 * t671;
t648 = t608 * t671;
t651 = t589 * t661;
t526 = t567 * t595 - t568 * t647 + (t572 * t660 - t595 * t651 - t648) * qJD(3);
t545 = -t648 + (t572 * t592 - t651) * t595;
t559 = -t572 * t589 - t592 * t661;
t537 = t545 * t598 + t559 * t594;
t634 = t537 * t587 - t526;
t603 = t572 * t647 + t595 * t608 - t671 * t651;
t527 = t603 * qJD(3) + t567 * t671 + t568 * t660;
t618 = -t545 * t594 + t559 * t598;
t517 = t618 * qJD(4) + t527 * t598 - t568 * t664;
t639 = t587 * t603 - t517;
t657 = (t639 * t585 - t634 * t586) * r_i_i_C(1) + (t634 * t585 + t639 * t586) * r_i_i_C(2);
t620 = t668 * t667;
t606 = t591 * t596 + t599 * t620;
t607 = -t591 * t599 + t596 * t620;
t643 = t590 * t667;
t622 = t589 * t643;
t547 = -t607 * t671 + (-t592 * t606 + t622) * t595;
t569 = t606 * qJD(2);
t570 = t607 * qJD(2);
t528 = t547 * qJD(3) - t569 * t595 - t570 * t647;
t560 = t589 * t606 + t592 * t643;
t539 = t547 * t598 + t560 * t594;
t633 = t539 * t587 - t528;
t602 = t595 * t607 - t606 * t647 + t671 * t622;
t529 = t602 * qJD(3) - t569 * t671 + t570 * t660;
t617 = -t547 * t594 + t560 * t598;
t519 = t617 * qJD(4) + t529 * t598 - t570 * t664;
t638 = t587 * t602 - t519;
t656 = (t638 * t585 - t633 * t586) * r_i_i_C(1) + (t633 * t585 + t638 * t586) * r_i_i_C(2);
t646 = t671 * t596;
t609 = t592 * t646 + t658;
t610 = t592 * t658 + t646;
t641 = t668 * t589;
t623 = t595 * t641;
t542 = qJD(3) * t623 + (t609 * qJD(2) + t610 * qJD(3)) * t590;
t558 = t610 * t590 + t623;
t571 = -t590 * t662 + t668 * t592;
t549 = t558 * t598 + t571 * t594;
t629 = t549 * t587 - t542;
t645 = t671 * t599;
t611 = -t592 * t659 + t645;
t615 = t671 * t641;
t627 = t592 * t645;
t543 = qJD(3) * t615 + ((t627 - t659) * qJD(3) + t611 * qJD(2)) * t590;
t616 = -t558 * t594 + t571 * t598;
t654 = qJD(2) * t590;
t644 = t589 * t654;
t626 = t596 * t644;
t525 = t616 * qJD(4) + t543 * t598 + t594 * t626;
t619 = t590 * t627;
t649 = t590 * t659;
t557 = -t615 - t619 + t649;
t635 = -t557 * t587 - t525;
t655 = (t635 * t585 - t629 * t586) * r_i_i_C(1) + (t629 * t585 + t635 * t586) * r_i_i_C(2);
t597 = cos(qJ(5));
t653 = qJD(5) * t597;
t650 = t589 * t590 * t596;
t613 = -t572 * t595 + t608 * t647;
t533 = t613 * qJD(3) - t567 * t660 + t568 * t671;
t553 = t572 * t671 + t608 * t660;
t521 = t567 * t664 + t533 * t598 + (-t553 * t594 - t608 * t663) * qJD(4);
t637 = -t587 * t613 + t521;
t612 = t595 * t606 + t607 * t647;
t535 = t612 * qJD(3) + t569 * t660 + t570 * t671;
t555 = -t606 * t671 + t607 * t660;
t523 = -t569 * t664 + t535 * t598 + (-t555 * t594 - t607 * t663) * qJD(4);
t636 = -t587 * t612 + t523;
t551 = (-t610 * qJD(2) - t609 * qJD(3)) * t590;
t566 = t611 * t590;
t625 = t599 * t644;
t531 = t594 * t625 + t551 * t598 + (-t566 * t594 + t598 * t650) * qJD(4);
t565 = t609 * t590;
t632 = t565 * t587 + t531;
t532 = t553 * qJD(3) + t567 * t647 + t568 * t595;
t540 = t553 * t598 - t608 * t664;
t631 = -t540 * t587 + t532;
t534 = t555 * qJD(3) - t569 * t647 + t570 * t595;
t541 = t555 * t598 - t607 * t664;
t630 = -t541 * t587 + t534;
t550 = -qJD(2) * t619 - t590 * qJD(3) * t645 + (qJD(3) * t592 + qJD(2)) * t649;
t556 = t566 * t598 + t594 * t650;
t628 = -t556 * t587 - t550;
t584 = t597 * pkin(5) + pkin(4);
t614 = t586 * r_i_i_C(1) - t585 * r_i_i_C(2) + t584;
t604 = -t669 * t594 - t614 * t598 - pkin(3);
t601 = t672 * t598 + (t614 * t594 - t669 * t598) * qJD(4);
t1 = [0, -t569 * t670 + t570 * pkin(2) + t535 * pkin(3) + t534 * pkin(10) + t523 * t584 + (t636 * r_i_i_C(1) + t630 * r_i_i_C(2)) * t586 + (t630 * r_i_i_C(1) - t636 * r_i_i_C(2)) * t585 + t669 * (t541 * qJD(4) + t535 * t594 + t569 * t663) + (t534 * t593 + (-t541 * t593 - t597 * t612) * qJD(5)) * pkin(5) (t529 * t585 + t547 * t665) * r_i_i_C(1) + (t529 * t586 - t547 * t666) * r_i_i_C(2) + t529 * pkin(10) + (t529 * t593 + t547 * t653) * pkin(5) + t604 * t528 - t601 * t602, t669 * t519 - t672 * t617 + t614 * (-t539 * qJD(4) - t529 * t594 - t570 * t663) (-t519 * t593 + t528 * t597 + (-t539 * t597 + t593 * t602) * qJD(5)) * pkin(5) + t656, t656; 0, t567 * t670 + t568 * pkin(2) + t533 * pkin(3) + t532 * pkin(10) + t521 * t584 + (t637 * r_i_i_C(1) + t631 * r_i_i_C(2)) * t586 + (t631 * r_i_i_C(1) - t637 * r_i_i_C(2)) * t585 + t669 * (t540 * qJD(4) + t533 * t594 - t567 * t663) + (t532 * t593 + (-t540 * t593 - t597 * t613) * qJD(5)) * pkin(5) (t527 * t585 + t545 * t665) * r_i_i_C(1) + (t527 * t586 - t545 * t666) * r_i_i_C(2) + t527 * pkin(10) + (t527 * t593 + t545 * t653) * pkin(5) + t604 * t526 - t601 * t603, t669 * t517 - t672 * t618 + t614 * (-t537 * qJD(4) - t527 * t594 - t568 * t663) (-t517 * t593 + t526 * t597 + (-t537 * t597 + t593 * t603) * qJD(5)) * pkin(5) + t657, t657; 0, t551 * pkin(3) - t550 * pkin(10) + t531 * t584 + (t632 * r_i_i_C(1) + t628 * r_i_i_C(2)) * t586 + (t628 * r_i_i_C(1) - t632 * r_i_i_C(2)) * t585 + t669 * (t556 * qJD(4) + t551 * t594 - t598 * t625) + (-pkin(2) * t596 + pkin(9) * t662) * t654 + (-t550 * t593 + (-t556 * t593 + t565 * t597) * qJD(5)) * pkin(5) (t543 * t585 + t558 * t665) * r_i_i_C(1) + (t543 * t586 - t558 * t666) * r_i_i_C(2) + t543 * pkin(10) + (t543 * t593 + t558 * t653) * pkin(5) + t604 * t542 + t601 * t557, t669 * t525 - t672 * t616 + t614 * (-t549 * qJD(4) - t543 * t594 + t598 * t626) (-t525 * t593 + t542 * t597 + (-t549 * t597 - t557 * t593) * qJD(5)) * pkin(5) + t655, t655;];
JaD_transl  = t1;
