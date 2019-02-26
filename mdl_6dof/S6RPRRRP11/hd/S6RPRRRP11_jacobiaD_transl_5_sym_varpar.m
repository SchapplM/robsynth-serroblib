% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function JaD_transl = S6RPRRRP11_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP11_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:46
% EndTime: 2019-02-26 21:13:47
% DurationCPUTime: 1.12s
% Computational Cost: add. (1177->127), mult. (3776->220), div. (0->0), fcn. (4194->14), ass. (0->86)
t584 = cos(qJ(1));
t641 = cos(pkin(12));
t643 = cos(pkin(6));
t617 = t643 * t641;
t639 = sin(pkin(12));
t645 = sin(qJ(1));
t596 = t584 * t639 + t645 * t617;
t559 = t596 * qJD(1);
t615 = t643 * t639;
t564 = t584 * t641 - t645 * t615;
t560 = t564 * qJD(1);
t581 = sin(qJ(3));
t594 = -t584 * t615 - t641 * t645;
t578 = sin(pkin(6));
t640 = sin(pkin(7));
t627 = t578 * t640;
t607 = t645 * t627;
t642 = cos(pkin(7));
t646 = cos(qJ(3));
t595 = t584 * t617 - t639 * t645;
t621 = t584 * t627;
t605 = t646 * t621;
t620 = t642 * t646;
t651 = t595 * t620 - t605;
t526 = t581 * (qJD(1) * t607 + qJD(3) * t594 - t559 * t642) + t560 * t646 + t651 * qJD(3);
t619 = t642 * t645;
t608 = t578 * t619;
t628 = t559 * t640;
t548 = qJD(1) * t608 + t628;
t580 = sin(qJ(4));
t583 = cos(qJ(4));
t624 = t642 * t595;
t631 = t594 * t646;
t541 = t581 * (t621 - t624) + t631;
t635 = t578 * t584;
t650 = t595 * t640 + t642 * t635;
t613 = t541 * t580 - t583 * t650;
t522 = -qJD(4) * t613 - t526 * t583 - t548 * t580;
t601 = t646 * t607;
t610 = t581 * t621;
t527 = qJD(1) * t601 - t559 * t620 - t560 * t581 + (-t581 * t624 + t610 + t631) * qJD(3);
t579 = sin(qJ(5));
t582 = cos(qJ(5));
t667 = t522 * t579 - t527 * t582;
t666 = t522 * t582 + t527 * t579;
t533 = t541 * t583 + t580 * t650;
t538 = -t581 * t594 - t651;
t665 = -t533 * t579 - t538 * t582;
t664 = t533 * t582 - t538 * t579;
t663 = qJD(4) * t533 - t526 * t580 + t548 * t583;
t647 = r_i_i_C(3) + pkin(11);
t590 = t596 * t642;
t543 = t564 * t646 + (-t590 + t607) * t581;
t553 = t596 * t640 + t608;
t649 = -t543 * t580 + t553 * t583;
t606 = qJD(5) * (t579 * r_i_i_C(1) + t582 * r_i_i_C(2));
t648 = -t564 * t581 - t646 * t590 + t601;
t614 = t642 * t641;
t616 = t643 * t640;
t549 = -t646 * t616 + (t581 * t639 - t614 * t646) * t578;
t585 = qJD(1) * t650;
t634 = qJD(5) * t579;
t633 = qJD(5) * t582;
t632 = t578 * qJD(2);
t630 = qJD(1) * t635;
t629 = qJD(1) * t645;
t535 = t543 * t583 + t553 * t580;
t550 = t581 * t616 + (t581 * t614 + t639 * t646) * t578;
t561 = -t627 * t641 + t642 * t643;
t537 = t550 * t583 + t561 * t580;
t612 = -t550 * t580 + t561 * t583;
t611 = r_i_i_C(1) * t582 - r_i_i_C(2) * t579 + pkin(4);
t593 = -t580 * t647 - t611 * t583 - pkin(3);
t588 = qJD(1) * t624;
t586 = t583 * t606 + (t611 * t580 - t583 * t647) * qJD(4);
t558 = t594 * qJD(1);
t546 = t550 * qJD(3);
t545 = t549 * qJD(3);
t530 = qJD(4) * t612 - t545 * t583;
t524 = qJD(1) * t610 + t648 * qJD(3) + t558 * t646 - t581 * t588;
t523 = -qJD(1) * t605 + qJD(3) * t543 + t558 * t581 + t588 * t646;
t518 = t649 * qJD(4) + t524 * t583 + t580 * t585;
t517 = qJD(4) * t535 + t524 * t580 - t583 * t585;
t516 = t518 * t582 + t523 * t579 + (-t535 * t579 - t582 * t648) * qJD(5);
t515 = -t518 * t579 + t523 * t582 + (-t535 * t582 + t579 * t648) * qJD(5);
t1 = [t666 * r_i_i_C(1) - t667 * r_i_i_C(2) + t522 * pkin(4) - t526 * pkin(3) + t527 * pkin(10) - t560 * pkin(2) - pkin(9) * t628 + t584 * t632 + t647 * t663 + (t665 * r_i_i_C(1) - t664 * r_i_i_C(2)) * qJD(5) + (-t584 * pkin(1) + (-pkin(9) * t619 - qJ(2) * t645) * t578) * qJD(1), t630 (t524 * t579 + t543 * t633) * r_i_i_C(1) + (t524 * t582 - t543 * t634) * r_i_i_C(2) + t524 * pkin(10) + t593 * t523 - t586 * t648, -t611 * t517 + t647 * t518 - t649 * t606, r_i_i_C(1) * t515 - r_i_i_C(2) * t516, 0; -pkin(1) * t629 + t558 * pkin(2) + t524 * pkin(3) + t518 * pkin(4) + pkin(9) * t585 + t523 * pkin(10) + t516 * r_i_i_C(1) + t515 * r_i_i_C(2) + qJ(2) * t630 + t647 * t517 + t645 * t632, t578 * t629 (t526 * t579 - t541 * t633) * r_i_i_C(1) + (t526 * t582 + t541 * t634) * r_i_i_C(2) + t526 * pkin(10) - t593 * t527 + t586 * t538, -t522 * t647 - t613 * t606 + t611 * t663, t667 * r_i_i_C(1) + t666 * r_i_i_C(2) + (t664 * r_i_i_C(1) + t665 * r_i_i_C(2)) * qJD(5), 0; 0, 0 (-t545 * t579 + t550 * t633) * r_i_i_C(1) + (-t545 * t582 - t550 * t634) * r_i_i_C(2) - t545 * pkin(10) + t593 * t546 + t586 * t549, t647 * t530 - t612 * t606 + t611 * (-qJD(4) * t537 + t545 * t580) (-t530 * t579 + t546 * t582) * r_i_i_C(1) + (-t530 * t582 - t546 * t579) * r_i_i_C(2) + ((-t537 * t582 - t549 * t579) * r_i_i_C(1) + (t537 * t579 - t549 * t582) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
