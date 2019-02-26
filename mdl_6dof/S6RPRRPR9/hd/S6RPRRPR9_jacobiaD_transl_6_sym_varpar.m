% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:33
% EndTime: 2019-02-26 21:05:35
% DurationCPUTime: 1.34s
% Computational Cost: add. (1500->152), mult. (4191->252), div. (0->0), fcn. (4652->16), ass. (0->94)
t600 = cos(qJ(1));
t656 = cos(pkin(12));
t658 = cos(pkin(6));
t632 = t658 * t656;
t654 = sin(pkin(12));
t663 = sin(qJ(1));
t613 = t600 * t654 + t663 * t632;
t569 = t613 * qJD(1);
t631 = t658 * t654;
t574 = t600 * t656 - t663 * t631;
t570 = t574 * qJD(1);
t597 = sin(qJ(3));
t611 = -t600 * t631 - t663 * t656;
t593 = sin(pkin(6));
t655 = sin(pkin(7));
t641 = t593 * t655;
t623 = t663 * t641;
t657 = cos(pkin(7));
t664 = cos(qJ(3));
t612 = t600 * t632 - t663 * t654;
t636 = t600 * t641;
t621 = t664 * t636;
t635 = t657 * t664;
t670 = t612 * t635 - t621;
t536 = (qJD(1) * t623 + qJD(3) * t611 - t657 * t569) * t597 + t570 * t664 + t670 * qJD(3);
t634 = t657 * t663;
t624 = t593 * t634;
t642 = t569 * t655;
t558 = qJD(1) * t624 + t642;
t592 = qJ(4) + pkin(13);
t590 = sin(t592);
t591 = cos(t592);
t652 = t593 * t600;
t669 = t612 * t655 + t657 * t652;
t582 = t597 * t636;
t671 = t612 * t657;
t679 = -t597 * t671 + t611 * t664 + t582;
t668 = -t590 * t679 + t591 * t669;
t532 = qJD(4) * t668 - t536 * t591 - t558 * t590;
t617 = t664 * t623;
t537 = qJD(1) * t617 + t679 * qJD(3) - t569 * t635 - t570 * t597;
t595 = sin(qJ(6));
t598 = cos(qJ(6));
t686 = t532 * t595 - t537 * t598;
t685 = t532 * t598 + t537 * t595;
t542 = -t590 * t669 - t591 * t679;
t548 = -t597 * t611 - t670;
t684 = t542 * t595 - t548 * t598;
t683 = t542 * t598 + t548 * t595;
t682 = -t542 * qJD(4) - t536 * t590 + t558 * t591;
t607 = t613 * t657;
t553 = t574 * t664 + (-t607 + t623) * t597;
t603 = t669 * qJD(1);
t678 = -t553 * qJD(4) + t603;
t665 = r_i_i_C(3) + pkin(11);
t667 = t678 * pkin(4);
t622 = qJD(6) * (t595 * r_i_i_C(1) + t598 * r_i_i_C(2));
t626 = t598 * r_i_i_C(1) - t595 * r_i_i_C(2) + pkin(5);
t666 = -t574 * t597 - t664 * t607 + t617;
t629 = t655 * t658;
t630 = t657 * t656;
t559 = -t664 * t629 + (t597 * t654 - t630 * t664) * t593;
t596 = sin(qJ(4));
t661 = t596 * pkin(4);
t651 = qJD(6) * t595;
t650 = qJD(6) * t598;
t563 = t613 * t655 + t624;
t648 = t563 * qJD(4);
t647 = t593 * qJD(2);
t645 = pkin(4) * t648;
t644 = qJD(1) * t652;
t643 = qJD(1) * t663;
t545 = t553 * t591 + t563 * t590;
t560 = t597 * t629 + (t597 * t630 + t664 * t654) * t593;
t571 = -t656 * t641 + t658 * t657;
t547 = t560 * t591 + t571 * t590;
t627 = -t560 * t590 + t571 * t591;
t599 = cos(qJ(4));
t589 = t599 * pkin(4) + pkin(3);
t610 = -t665 * t590 - t626 * t591 - t589;
t605 = qJD(1) * t671;
t601 = t591 * t622 + (t626 * t590 - t665 * t591 + t661) * qJD(4);
t594 = -qJ(5) - pkin(10);
t568 = t611 * qJD(1);
t556 = t560 * qJD(3);
t555 = t559 * qJD(3);
t540 = t627 * qJD(4) - t555 * t591;
t534 = qJD(1) * t582 + t666 * qJD(3) + t568 * t664 - t597 * t605;
t533 = -qJD(1) * t621 + t553 * qJD(3) + t568 * t597 + t664 * t605;
t528 = (t534 + t648) * t591 + t678 * t590;
t527 = t545 * qJD(4) + t534 * t590 - t591 * t603;
t526 = t528 * t598 + t533 * t595 + (-t545 * t595 - t598 * t666) * qJD(6);
t525 = -t528 * t595 + t533 * t598 + (-t545 * t598 + t595 * t666) * qJD(6);
t1 = [t685 * r_i_i_C(1) - t686 * r_i_i_C(2) + t532 * pkin(5) - t536 * t589 - t537 * t594 - t548 * qJD(5) - t570 * pkin(2) - pkin(9) * t642 + t600 * t647 + t665 * t682 + (t684 * r_i_i_C(1) + t683 * r_i_i_C(2)) * qJD(6) + (-t600 * pkin(1) + (-pkin(9) * t634 - t663 * qJ(2)) * t593) * qJD(1) + (-t558 * t596 + (-t596 * t679 + t599 * t669) * qJD(4)) * pkin(4), t644 (t534 * t595 + t553 * t650) * r_i_i_C(1) + (t534 * t598 - t553 * t651) * r_i_i_C(2) - t534 * t594 + t553 * qJD(5) + t610 * t533 - t601 * t666, -t534 * t661 - t596 * t645 + t667 * t599 + (-r_i_i_C(1) * t651 - r_i_i_C(2) * t650) * (-t553 * t590 + t563 * t591) + t665 * t528 - t626 * t527, t533, t525 * r_i_i_C(1) - t526 * r_i_i_C(2); -pkin(1) * t643 + t568 * pkin(2) + t528 * pkin(5) + pkin(9) * t603 + t526 * r_i_i_C(1) + t525 * r_i_i_C(2) + qJ(2) * t644 - qJD(5) * t666 + t665 * t527 - t533 * t594 + t534 * t589 + t667 * t596 + t599 * t645 + t663 * t647, t593 * t643 (t536 * t595 - t650 * t679) * r_i_i_C(1) + (t536 * t598 + t651 * t679) * r_i_i_C(2) - t536 * t594 - t679 * qJD(5) - t610 * t537 + t601 * t548, -t665 * t532 + t668 * t622 + t626 * t682 + (-t536 * t596 + t558 * t599 + (t596 * t669 + t599 * t679) * qJD(4)) * pkin(4), -t537, t686 * r_i_i_C(1) + t685 * r_i_i_C(2) + (-t683 * r_i_i_C(1) + t684 * r_i_i_C(2)) * qJD(6); 0, 0 (-t555 * t595 + t560 * t650) * r_i_i_C(1) + (-t555 * t598 - t560 * t651) * r_i_i_C(2) + t555 * t594 + t560 * qJD(5) + t610 * t556 + t601 * t559, t665 * t540 - t627 * t622 + t626 * (-t547 * qJD(4) + t555 * t590) + (t555 * t596 + (-t560 * t599 - t571 * t596) * qJD(4)) * pkin(4), t556 (-t540 * t595 + t556 * t598) * r_i_i_C(1) + (-t540 * t598 - t556 * t595) * r_i_i_C(2) + ((-t547 * t598 - t559 * t595) * r_i_i_C(1) + (t547 * t595 - t559 * t598) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
