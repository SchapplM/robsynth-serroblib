% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR11
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR11_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR11_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:55
% EndTime: 2019-02-26 22:21:56
% DurationCPUTime: 1.18s
% Computational Cost: add. (1501->145), mult. (4424->242), div. (0->0), fcn. (4722->12), ass. (0->84)
t625 = sin(qJ(2));
t626 = sin(qJ(1));
t630 = cos(qJ(2));
t676 = cos(pkin(6));
t677 = cos(qJ(1));
t654 = t676 * t677;
t607 = t625 * t654 + t626 * t630;
t624 = sin(qJ(3));
t629 = cos(qJ(3));
t621 = sin(pkin(6));
t665 = t621 * t677;
t597 = t607 * t624 + t629 * t665;
t598 = t607 * t629 - t624 * t665;
t623 = sin(qJ(5));
t628 = cos(qJ(5));
t581 = t597 * t623 + t598 * t628;
t606 = t626 * t625 - t630 * t654;
t622 = sin(qJ(6));
t627 = cos(qJ(6));
t705 = t581 * t622 + t606 * t627;
t704 = t581 * t627 - t606 * t622;
t672 = t621 * t626;
t698 = qJD(1) * t672 - qJD(3) * t607;
t659 = t626 * t676;
t655 = t625 * t659;
t661 = t677 * qJD(1);
t669 = qJD(2) * t625;
t594 = -qJD(1) * t655 - t626 * t669 + (qJD(2) * t654 + t661) * t630;
t699 = -qJD(3) * t665 + t594;
t573 = t699 * t624 - t698 * t629;
t574 = t698 * t624 + t699 * t629;
t703 = t581 * qJD(5) - t573 * t628 + t574 * t623;
t646 = t597 * t628 - t598 * t623;
t559 = t646 * qJD(5) + t573 * t623 + t574 * t628;
t643 = t623 * t629 - t624 * t628;
t686 = qJD(5) - qJD(3);
t632 = t686 * t643;
t634 = t677 * t625 + t630 * t659;
t592 = t607 * qJD(1) + t634 * qJD(2);
t635 = -t677 * t630 + t655;
t602 = t624 * t672 - t629 * t635;
t656 = t621 * t661;
t571 = t602 * qJD(3) - t592 * t624 - t629 * t656;
t671 = t621 * t629;
t639 = t624 * t635 + t626 * t671;
t572 = t639 * qJD(3) - t592 * t629 + t624 * t656;
t585 = t602 * t628 - t623 * t639;
t554 = t585 * qJD(5) - t571 * t628 + t572 * t623;
t645 = -t602 * t623 - t628 * t639;
t555 = t645 * qJD(5) + t571 * t623 + t572 * t628;
t652 = -t622 * r_i_i_C(1) - t627 * r_i_i_C(2);
t638 = qJD(6) * t652;
t653 = -t627 * r_i_i_C(1) + t622 * r_i_i_C(2);
t641 = pkin(5) - t653;
t678 = -r_i_i_C(3) - pkin(11);
t697 = -t641 * t554 - t678 * t555 + t645 * t638;
t680 = pkin(3) + pkin(4);
t695 = (-qJ(4) * t629 + t680 * t624) * qJD(3) - t624 * qJD(4);
t605 = t676 * t624 + t625 * t671;
t670 = t621 * t630;
t662 = qJD(2) * t670;
t595 = t605 * qJD(3) + t624 * t662;
t604 = t621 * t625 * t624 - t676 * t629;
t596 = -t604 * qJD(3) + t629 * t662;
t644 = t604 * t628 - t605 * t623;
t569 = t644 * qJD(5) + t595 * t623 + t596 * t628;
t590 = t604 * t623 + t605 * t628;
t694 = -t641 * (qJD(5) * t590 - t595 * t628 + t596 * t623) - t678 * t569 + t644 * t638;
t693 = -t678 * t559 + t646 * t638 - t641 * t703;
t642 = t623 * t624 + t628 * t629;
t692 = t642 * t634;
t684 = t624 * qJ(4) + t680 * t629 + pkin(2);
t679 = -pkin(9) + pkin(10);
t667 = qJD(6) * t642 * t670;
t663 = t621 * t669;
t640 = -t652 + t679;
t631 = t686 * t642;
t593 = t634 * qJD(1) + t607 * qJD(2);
t591 = t606 * qJD(1) + t635 * qJD(2);
t586 = t642 * t606;
t578 = (-t632 * t630 - t642 * t669) * t621;
t551 = t555 * t627 + t591 * t622 + (-t585 * t622 - t627 * t634) * qJD(6);
t550 = -t555 * t622 + t591 * t627 + (-t585 * t627 + t622 * t634) * qJD(6);
t1 = [-t594 * pkin(2) - t573 * qJ(4) - t597 * qJD(4) - t641 * t559 + t640 * t593 - t680 * t574 + t678 * t703 + (t705 * r_i_i_C(1) + t704 * r_i_i_C(2)) * qJD(6) + (-t677 * pkin(1) - pkin(8) * t672) * qJD(1), t641 * (t642 * t591 + t632 * t634) + t640 * t592 - t678 * (t643 * t591 - t686 * t692) + ((t622 * t692 + t627 * t635) * r_i_i_C(1) + (-t622 * t635 + t627 * t692) * r_i_i_C(2)) * qJD(6) + t695 * t634 + t684 * t591, t572 * qJ(4) + t602 * qJD(4) - t680 * t571 - t697, t571, t697, t550 * r_i_i_C(1) - t551 * r_i_i_C(2); -t592 * pkin(2) + t555 * pkin(5) + t551 * r_i_i_C(1) + t550 * r_i_i_C(2) + t571 * qJ(4) - t639 * qJD(4) + t679 * t591 + t680 * t572 - t678 * t554 + (-pkin(1) * t626 + pkin(8) * t665) * qJD(1), t641 * (-t642 * t593 + t606 * t632) - t640 * t594 - t678 * (-t643 * t593 - t631 * t606) + ((t586 * t622 - t607 * t627) * r_i_i_C(1) + (t586 * t627 + t607 * t622) * r_i_i_C(2)) * qJD(6) + t695 * t606 - t684 * t593, t574 * qJ(4) + t598 * qJD(4) - t680 * t573 - t693, t573, t693 (-t559 * t622 - t593 * t627) * r_i_i_C(1) + (-t559 * t627 + t593 * t622) * r_i_i_C(2) + (-t704 * r_i_i_C(1) + t705 * r_i_i_C(2)) * qJD(6); 0 (t578 * t627 - t622 * t667) * r_i_i_C(1) + (-t578 * t622 - t627 * t667) * r_i_i_C(2) + t578 * pkin(5) + (t678 * (-t630 * t631 + t643 * t669) + t653 * t625 * qJD(6) - t695 * t630 + (-t625 * t684 - t630 * t640) * qJD(2)) * t621, t596 * qJ(4) + t605 * qJD(4) - t680 * t595 - t694, t595, t694 (-t569 * t622 - t627 * t663) * r_i_i_C(1) + (-t569 * t627 + t622 * t663) * r_i_i_C(2) + ((-t590 * t627 - t622 * t670) * r_i_i_C(1) + (t590 * t622 - t627 * t670) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
