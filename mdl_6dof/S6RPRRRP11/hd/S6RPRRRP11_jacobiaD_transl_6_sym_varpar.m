% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP11_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP11_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:36
% EndTime: 2019-02-26 21:13:37
% DurationCPUTime: 1.28s
% Computational Cost: add. (1519->141), mult. (4786->227), div. (0->0), fcn. (5354->14), ass. (0->91)
t592 = sin(qJ(3));
t662 = sin(pkin(7));
t663 = sin(pkin(6));
t642 = t663 * t662;
t672 = cos(qJ(1));
t619 = t672 * t642;
t664 = cos(pkin(12));
t666 = cos(pkin(6));
t646 = t666 * t664;
t661 = sin(pkin(12));
t670 = sin(qJ(1));
t572 = -t646 * t672 + t661 * t670;
t665 = cos(pkin(7));
t652 = t665 * t572;
t644 = t666 * t661;
t573 = t644 * t672 + t664 * t670;
t671 = cos(qJ(3));
t653 = t573 * t671;
t552 = t592 * (t619 + t652) - t653;
t643 = t665 * t663;
t563 = t572 * t662 - t672 * t643;
t591 = sin(qJ(4));
t594 = cos(qJ(4));
t544 = t552 * t594 - t563 * t591;
t609 = t671 * t619;
t678 = -t671 * t652 - t609;
t549 = t573 * t592 - t678;
t590 = sin(qJ(5));
t593 = cos(qJ(5));
t687 = -t544 * t590 - t549 * t593;
t635 = t544 * t593 - t549 * t590;
t574 = -t670 * t644 + t672 * t664;
t570 = t574 * qJD(1);
t615 = t670 * t642;
t602 = t670 * t646 + t672 * t661;
t656 = t602 * qJD(1);
t676 = t656 * t665;
t539 = -t678 * qJD(3) - t570 * t671 + (-qJD(1) * t615 + t573 * qJD(3) + t676) * t592;
t616 = t670 * t643;
t638 = t656 * t662;
t601 = qJD(1) * t616 + t638;
t686 = qJD(4) * t544 + t539 * t591 + t594 * t601;
t632 = t552 * t591 + t563 * t594;
t531 = qJD(4) * t632 - t539 * t594 + t591 * t601;
t683 = t664 * t643 + t666 * t662;
t673 = pkin(5) + r_i_i_C(1);
t625 = t593 * r_i_i_C(2) + t590 * t673;
t613 = qJD(5) * t625;
t608 = t671 * t615;
t612 = t592 * t619;
t536 = t671 * t676 - qJD(1) * t608 + t570 * t592 - (t592 * t652 + t612 - t653) * qJD(3);
t668 = r_i_i_C(3) + qJ(6) + pkin(11);
t600 = t602 * t665;
t554 = t574 * t671 + (-t600 + t615) * t592;
t564 = t602 * t662 + t616;
t545 = -t554 * t591 + t564 * t594;
t675 = -t574 * t592 - t671 * t600 + t608;
t641 = t663 * t661;
t560 = t592 * t641 - t671 * t683;
t596 = t563 * qJD(1);
t569 = t573 * qJD(1);
t598 = qJD(1) * t652;
t534 = -qJD(1) * t609 + qJD(3) * t554 - t569 * t592 - t598 * t671;
t660 = t534 * t590;
t655 = qJD(5) * t590;
t654 = qJD(5) * t593;
t650 = t672 * t663;
t648 = t670 * t663;
t637 = -t531 * t590 + t536 * t593;
t557 = t560 * qJD(3);
t561 = t592 * t683 + t671 * t641;
t571 = -t642 * t664 + t665 * t666;
t630 = -t561 * t591 + t571 * t594;
t541 = qJD(4) * t630 - t557 * t594;
t558 = t561 * qJD(3);
t636 = -t541 * t590 + t558 * t593;
t548 = t561 * t594 + t571 * t591;
t633 = -t548 * t593 - t560 * t590;
t546 = t554 * t594 + t564 * t591;
t588 = t593 * pkin(5) + pkin(4);
t628 = t593 * r_i_i_C(1) - t590 * r_i_i_C(2) + t588;
t627 = qJD(1) * t650;
t603 = -t591 * t668 - t594 * t628 - pkin(3);
t535 = qJD(1) * t612 + t675 * qJD(3) - t569 * t671 + t592 * t598;
t529 = t545 * qJD(4) + t535 * t594 - t591 * t596;
t526 = -t529 * t590 + t534 * t593 + (-t546 * t593 + t590 * t675) * qJD(5);
t595 = -t591 * qJD(6) + t594 * t613 + (t591 * t628 - t594 * t668) * qJD(4);
t540 = qJD(4) * t548 - t557 * t591;
t528 = qJD(4) * t546 + t535 * t591 + t594 * t596;
t527 = t529 * t593 + t660 + (-t546 * t590 - t593 * t675) * qJD(5);
t1 = [t632 * qJD(6) + t539 * pkin(3) - t570 * pkin(2) - pkin(9) * t638 + qJD(2) * t650 - t628 * t531 - (pkin(10) + t625) * t536 + t668 * t686 + (-pkin(1) * t672 - pkin(9) * t616 - qJ(2) * t648) * qJD(1) + (-t635 * r_i_i_C(2) + t673 * t687) * qJD(5), t627 (t535 * t593 - t554 * t655) * r_i_i_C(2) + t535 * pkin(10) + t603 * t534 - t595 * t675 + t673 * (t535 * t590 + t554 * t654) t546 * qJD(6) - t528 * t628 + t529 * t668 - t545 * t613, -t527 * r_i_i_C(2) + t673 * t526, t528; -qJD(1) * t670 * pkin(1) - t569 * pkin(2) + t535 * pkin(3) + t534 * pkin(10) + t527 * r_i_i_C(1) + t526 * r_i_i_C(2) + qJ(2) * t627 + qJD(2) * t648 - t545 * qJD(6) + t529 * t588 - pkin(9) * t596 + t668 * t528 + (-t546 * t655 - t654 * t675 + t660) * pkin(5), qJD(1) * t648 (-t539 * t593 + t552 * t655) * r_i_i_C(2) - t539 * pkin(10) + t603 * t536 + t595 * t549 + t673 * (-t539 * t590 - t552 * t654) -qJD(6) * t544 + t531 * t668 - t613 * t632 + t628 * t686, t637 * r_i_i_C(1) + (-t531 * t593 - t536 * t590) * r_i_i_C(2) + (t635 * r_i_i_C(1) + t687 * r_i_i_C(2)) * qJD(5) + (qJD(5) * t635 + t637) * pkin(5), -t686; 0, 0 (-t557 * t593 - t561 * t655) * r_i_i_C(2) - t557 * pkin(10) + t603 * t558 + t595 * t560 + t673 * (-t557 * t590 + t561 * t654) t548 * qJD(6) - t540 * t628 + t541 * t668 - t613 * t630, t636 * r_i_i_C(1) + (-t541 * t593 - t558 * t590) * r_i_i_C(2) + (t633 * r_i_i_C(1) + (t548 * t590 - t560 * t593) * r_i_i_C(2)) * qJD(5) + (qJD(5) * t633 + t636) * pkin(5), t540;];
JaD_transl  = t1;
