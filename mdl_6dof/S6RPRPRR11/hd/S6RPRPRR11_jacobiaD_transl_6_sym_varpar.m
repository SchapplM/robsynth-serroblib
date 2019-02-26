% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR11
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
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR11_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR11_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:36
% EndTime: 2019-02-26 20:54:37
% DurationCPUTime: 1.17s
% Computational Cost: add. (1421->136), mult. (3937->227), div. (0->0), fcn. (4375->16), ass. (0->89)
t601 = cos(qJ(1));
t654 = cos(pkin(12));
t656 = cos(pkin(6));
t632 = t656 * t654;
t652 = sin(pkin(12));
t659 = sin(qJ(1));
t613 = t601 * t652 + t659 * t632;
t571 = t613 * qJD(1);
t630 = t656 * t652;
t576 = t601 * t654 - t659 * t630;
t572 = t576 * qJD(1);
t599 = sin(qJ(3));
t611 = -t601 * t630 - t659 * t654;
t596 = sin(pkin(6));
t653 = sin(pkin(7));
t641 = t596 * t653;
t623 = t659 * t641;
t655 = cos(pkin(7));
t660 = cos(qJ(3));
t612 = t601 * t632 - t659 * t652;
t638 = t601 * t641;
t621 = t660 * t638;
t635 = t655 * t660;
t666 = t612 * t635 - t621;
t538 = (qJD(1) * t623 + qJD(3) * t611 - t655 * t571) * t599 + t572 * t660 + t666 * qJD(3);
t634 = t655 * t659;
t624 = t596 * t634;
t642 = t571 * t653;
t560 = qJD(1) * t624 + t642;
t594 = pkin(13) + qJ(5);
t592 = sin(t594);
t593 = cos(t594);
t648 = t596 * t601;
t665 = t612 * t653 + t655 * t648;
t584 = t599 * t638;
t667 = t612 * t655;
t674 = -t599 * t667 + t611 * t660 + t584;
t664 = -t592 * t674 + t593 * t665;
t534 = qJD(5) * t664 - t538 * t593 - t560 * t592;
t617 = t660 * t623;
t539 = qJD(1) * t617 + t674 * qJD(3) - t571 * t635 - t572 * t599;
t598 = sin(qJ(6));
t600 = cos(qJ(6));
t681 = t534 * t598 - t539 * t600;
t680 = t534 * t600 + t539 * t598;
t544 = -t592 * t665 - t593 * t674;
t550 = -t599 * t611 - t666;
t679 = t544 * t598 - t550 * t600;
t678 = t544 * t600 + t550 * t598;
t677 = -t544 * qJD(5) - t538 * t592 + t560 * t593;
t661 = r_i_i_C(3) + pkin(11);
t607 = t613 * t655;
t555 = t576 * t660 + (-t607 + t623) * t599;
t565 = t613 * t653 + t624;
t663 = -t555 * t592 + t565 * t593;
t622 = qJD(6) * (t598 * r_i_i_C(1) + t600 * r_i_i_C(2));
t662 = -t576 * t599 - t660 * t607 + t617;
t629 = t655 * t654;
t631 = t656 * t653;
t561 = -t660 * t631 + (t599 * t652 - t629 * t660) * t596;
t602 = qJD(1) * t665;
t658 = pkin(4) * sin(pkin(13));
t647 = qJD(2) * t596;
t646 = qJD(6) * t598;
t645 = qJD(6) * t600;
t644 = qJD(1) * t648;
t643 = qJD(1) * t659;
t547 = t555 * t593 + t565 * t592;
t562 = t599 * t631 + (t599 * t629 + t660 * t652) * t596;
t573 = -t654 * t641 + t656 * t655;
t549 = t562 * t593 + t573 * t592;
t627 = -t562 * t592 + t573 * t593;
t626 = t600 * r_i_i_C(1) - t598 * r_i_i_C(2) + pkin(5);
t591 = cos(pkin(13)) * pkin(4) + pkin(3);
t610 = -t661 * t592 - t626 * t593 - t591;
t605 = qJD(1) * t667;
t603 = t593 * t622 + (t626 * t592 - t661 * t593) * qJD(5);
t597 = -pkin(10) - qJ(4);
t570 = t611 * qJD(1);
t558 = t562 * qJD(3);
t557 = t561 * qJD(3);
t542 = t627 * qJD(5) - t557 * t593;
t536 = qJD(1) * t584 + t662 * qJD(3) + t570 * t660 - t599 * t605;
t535 = -qJD(1) * t621 + t555 * qJD(3) + t570 * t599 + t660 * t605;
t530 = t663 * qJD(5) + t536 * t593 + t592 * t602;
t529 = t547 * qJD(5) + t536 * t592 - t593 * t602;
t528 = t530 * t600 + t535 * t598 + (-t547 * t598 - t600 * t662) * qJD(6);
t527 = -t530 * t598 + t535 * t600 + (-t547 * t600 + t598 * t662) * qJD(6);
t1 = [t680 * r_i_i_C(1) - t681 * r_i_i_C(2) + t534 * pkin(5) - t538 * t591 - t539 * t597 - t550 * qJD(4) - t560 * t658 - t572 * pkin(2) - pkin(9) * t642 + t601 * t647 + t661 * t677 + (t679 * r_i_i_C(1) + t678 * r_i_i_C(2)) * qJD(6) + (-t601 * pkin(1) + (-pkin(9) * t634 - t659 * qJ(2)) * t596) * qJD(1), t644 (t536 * t598 + t555 * t645) * r_i_i_C(1) + (t536 * t600 - t555 * t646) * r_i_i_C(2) - t536 * t597 + t555 * qJD(4) + t610 * t535 - t603 * t662, t535, -t626 * t529 + t661 * t530 - t663 * t622, t527 * r_i_i_C(1) - t528 * r_i_i_C(2); -pkin(1) * t643 + t570 * pkin(2) + t530 * pkin(5) + t528 * r_i_i_C(1) + t527 * r_i_i_C(2) + qJ(2) * t644 - t662 * qJD(4) + t661 * t529 - t535 * t597 + t536 * t591 + t659 * t647 + (pkin(9) + t658) * t602, t596 * t643 (t538 * t598 - t645 * t674) * r_i_i_C(1) + (t538 * t600 + t646 * t674) * r_i_i_C(2) - t538 * t597 - t674 * qJD(4) - t610 * t539 + t603 * t550, -t539, -t534 * t661 + t664 * t622 + t626 * t677, t681 * r_i_i_C(1) + t680 * r_i_i_C(2) + (-t678 * r_i_i_C(1) + t679 * r_i_i_C(2)) * qJD(6); 0, 0 (-t557 * t598 + t562 * t645) * r_i_i_C(1) + (-t557 * t600 - t562 * t646) * r_i_i_C(2) + t557 * t597 + t562 * qJD(4) + t610 * t558 + t603 * t561, t558, t661 * t542 - t627 * t622 + t626 * (-t549 * qJD(5) + t557 * t592) (-t542 * t598 + t558 * t600) * r_i_i_C(1) + (-t542 * t600 - t558 * t598) * r_i_i_C(2) + ((-t549 * t600 - t561 * t598) * r_i_i_C(1) + (t549 * t598 - t561 * t600) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
