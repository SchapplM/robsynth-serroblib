% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR14_jacobiaD_rot_5_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_5_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_5_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_rot_5_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:23
% EndTime: 2018-12-10 18:38:30
% DurationCPUTime: 6.86s
% Computational Cost: add. (93051->265), mult. (94962->446), div. (705->12), fcn. (88776->31), ass. (0->214)
t719 = pkin(7) + pkin(14);
t686 = sin(t719) / 0.2e1;
t720 = pkin(7) - pkin(14);
t701 = sin(t720);
t609 = t686 - t701 / 0.2e1;
t631 = cos(pkin(14));
t728 = pkin(6) + qJ(2);
t710 = sin(t728);
t692 = t710 / 0.2e1;
t729 = pkin(6) - qJ(2);
t711 = sin(t729);
t693 = -t711 / 0.2e1;
t671 = t692 + t693;
t754 = sin(qJ(1));
t755 = cos(qJ(2));
t698 = t754 * t755;
t756 = cos(qJ(1));
t652 = t671 * t756 + t698;
t687 = cos(t720) / 0.2e1;
t702 = cos(t719);
t610 = t687 - t702 / 0.2e1;
t751 = sin(pkin(6));
t716 = t610 * t751;
t650 = -t631 * t652 + t756 * t716;
t753 = sin(qJ(2));
t629 = t754 * t753;
t696 = cos(t728) / 0.2e1;
t713 = cos(t729);
t673 = t713 / 0.2e1 + t696;
t657 = -t756 * t673 + t629;
t573 = t609 * t657 + t650;
t635 = sin(qJ(4));
t726 = pkin(8) + qJ(4);
t708 = sin(t726);
t690 = t708 / 0.2e1;
t727 = pkin(8) - qJ(4);
t709 = sin(t727);
t691 = t709 / 0.2e1;
t670 = t690 + t691;
t694 = cos(t726) / 0.2e1;
t712 = cos(t727);
t672 = t712 / 0.2e1 + t694;
t667 = t686 + t701 / 0.2e1;
t660 = t667 * t751;
t656 = t756 * t660;
t668 = t687 + t702 / 0.2e1;
t749 = sin(pkin(14));
t758 = t652 * t749 + t657 * t668 + t656;
t752 = cos(pkin(7));
t688 = t752 * t751;
t678 = t756 * t688;
t750 = sin(pkin(7));
t759 = -t657 * t750 + t678;
t535 = t573 * t635 - t670 * t759 - t672 * t758;
t611 = t690 - t709 / 0.2e1;
t613 = t694 - t712 / 0.2e1;
t637 = cos(qJ(4));
t536 = t573 * t637 + t611 * t758 - t613 * t759;
t625 = -t710 / 0.2e1;
t626 = t711 / 0.2e1;
t733 = t626 + t625;
t606 = t733 * qJD(2);
t699 = t756 * t753;
t653 = t673 * t754 + t699;
t588 = qJD(1) * t653 + qJD(2) * t698 - t606 * t756;
t700 = t756 * t755;
t654 = -t754 * t671 + t700;
t607 = t673 * qJD(2);
t689 = qJD(2) * t629 - t607 * t756;
t589 = qJD(1) * t654 - t689;
t683 = t754 * t716;
t555 = qJD(1) * t683 - t588 * t609 + t589 * t631;
t655 = t754 * t660;
t556 = -qJD(1) * t655 + t588 * t668 + t589 * t749;
t677 = t754 * t688;
t658 = qJD(1) * t677 + t588 * t750;
t675 = t613 * qJD(4);
t676 = (t691 - t708 / 0.2e1) * qJD(4);
t731 = qJD(4) * t637;
t766 = -t555 * t635 - t556 * t672 + t573 * t731 + t658 * t670 - t675 * t759 - t676 * t758;
t603 = t670 * qJD(4);
t604 = t672 * qJD(4);
t732 = qJD(4) * t635;
t507 = t555 * t637 - t556 * t611 + t573 * t732 - t603 * t759 - t604 * t758 - t613 * t658;
t530 = t535 ^ 2;
t612 = t692 + t626;
t614 = t696 - t713 / 0.2e1;
t633 = cos(pkin(6));
t582 = t612 * t668 + t614 * t749 + t633 * t667;
t583 = t609 * t612 + t610 * t633 - t614 * t631;
t599 = -t612 * t750 + t633 * t752;
t549 = -t582 * t672 + t583 * t635 - t599 * t670;
t547 = 0.1e1 / t549 ^ 2;
t515 = t530 * t547 + 0.1e1;
t513 = 0.1e1 / t515;
t605 = (t625 + t693) * qJD(2);
t608 = t614 * qJD(2);
t591 = t605 * t749 + t608 * t668;
t592 = -t605 * t631 + t608 * t609;
t661 = t670 * t750;
t525 = -t582 * t676 + t583 * t731 - t591 * t672 + t592 * t635 - t599 * t675 + t608 * t661;
t546 = 0.1e1 / t549;
t735 = t535 * t547;
t682 = -t525 * t735 + t546 * t766;
t488 = t682 * t513;
t516 = atan2(t535, t549);
t511 = sin(t516);
t512 = cos(t516);
t685 = -t511 * t549 + t512 * t535;
t483 = t488 * t685 + t511 * t766 + t512 * t525;
t500 = t511 * t535 + t512 * t549;
t498 = 0.1e1 / t500 ^ 2;
t761 = t483 * t498;
t760 = t525 * t547;
t497 = 0.1e1 / t500;
t630 = sin(pkin(8));
t632 = cos(pkin(8));
t645 = t653 * t668 + t654 * t749 - t655;
t649 = -t653 * t750 - t677;
t561 = t630 * t645 - t632 * t649;
t634 = sin(qJ(5));
t636 = cos(qJ(5));
t574 = -t609 * t653 + t631 * t654 + t683;
t644 = t574 * t637 - t611 * t645 + t613 * t649;
t524 = t561 * t634 + t636 * t644;
t518 = 0.1e1 / t524;
t519 = 0.1e1 / t524 ^ 2;
t757 = 0.2e1 * t535;
t537 = t574 * t635 + t645 * t672 + t649 * t670;
t531 = t537 ^ 2;
t494 = t498 * t531 + 0.1e1;
t586 = qJD(1) * t657 - qJD(2) * t700 - t754 * t606;
t669 = -qJD(2) * t699 - t607 * t754;
t554 = qJD(1) * t650 + t586 * t609 + t631 * t669;
t646 = (-qJD(1) * t652 + t669) * t749 - t586 * t668 - qJD(1) * t656;
t659 = -qJD(1) * t678 + t586 * t750;
t504 = t554 * t635 + t574 * t731 + t645 * t676 + t646 * t672 + t649 * t675 + t659 * t670;
t741 = t498 * t537;
t747 = t497 * t761;
t748 = (t504 * t741 - t531 * t747) / t494 ^ 2;
t505 = t554 * t637 - t574 * t732 - t603 * t649 - t604 * t645 - t611 * t646 + t613 * t659;
t540 = t630 * t646 - t632 * t659;
t495 = qJD(5) * t524 + t505 * t634 - t540 * t636;
t523 = -t561 * t636 + t634 * t644;
t517 = t523 ^ 2;
t503 = t517 * t519 + 0.1e1;
t738 = t519 * t523;
t730 = qJD(5) * t523;
t496 = t505 * t636 + t540 * t634 - t730;
t743 = t496 * t518 * t519;
t745 = (t495 * t738 - t517 * t743) / t503 ^ 2;
t737 = t546 * t760;
t744 = (-t530 * t737 + t735 * t766) / t515 ^ 2;
t742 = t498 * t504;
t740 = t511 * t537;
t739 = t512 * t537;
t736 = t535 * t546;
t725 = 0.2e1 * t748;
t724 = 0.2e1 * t747;
t723 = -0.2e1 * t745;
t722 = 0.2e1 * t745;
t721 = -0.2e1 * t744;
t718 = t546 * t744;
t717 = t523 * t743;
t715 = t613 * t750;
t714 = t632 * t750;
t707 = -0.2e1 * t497 * t748;
t706 = t498 * t725;
t705 = t537 * t724;
t704 = 0.2e1 * t717;
t703 = t737 * t757;
t560 = -t630 * t758 + t632 * t759;
t522 = t536 * t636 + t560 * t634;
t521 = t536 * t634 - t560 * t636;
t601 = -t733 * t754 - t700;
t580 = t601 * t668 + t653 * t749;
t581 = t601 * t609 - t631 * t653;
t545 = t580 * t611 + t581 * t637 + t601 * t715;
t566 = -t580 * t630 - t601 * t714;
t529 = t545 * t636 + t566 * t634;
t528 = t545 * t634 - t566 * t636;
t681 = -t518 * t634 + t636 * t738;
t550 = t582 * t611 + t583 * t637 - t599 * t613;
t680 = t536 * t546 - t550 * t735;
t666 = -t733 * t756 + t698;
t579 = -t609 * t666 - t631 * t657;
t648 = -t657 * t749 + t666 * t668;
t543 = t579 * t635 + t648 * t672 - t661 * t666;
t596 = -t612 * t749 + t614 * t668;
t597 = t609 * t614 + t612 * t631;
t559 = -t596 * t672 + t597 * t635 + t614 * t661;
t679 = -t543 * t546 - t559 * t735;
t674 = -t511 + (-t512 * t736 + t511) * t513;
t665 = t675 * t750;
t590 = qJD(1) * t601 + t689;
t587 = qJD(1) * t666 - t669;
t563 = t586 * t631 + t587 * t609;
t562 = -t586 * t749 + t587 * t668;
t544 = -t580 * t672 + t581 * t635 + t601 * t661;
t542 = -t562 * t630 - t587 * t714;
t541 = -t556 * t630 - t632 * t658;
t527 = (t605 * t609 + t608 * t631) * t635 + t597 * t731 - (t605 * t668 - t608 * t749) * t672 - t596 * t676 + t605 * t661 + t614 * t665;
t526 = t582 * t604 - t583 * t732 + t591 * t611 + t592 * t637 + t599 * t603 + t608 * t715;
t510 = (-t588 * t631 + t590 * t609) * t635 + t579 * t731 - (t588 * t749 + t590 * t668) * t672 + t648 * t676 + t590 * t661 - t666 * t665;
t509 = -t601 * t603 * t750 + t562 * t611 + t563 * t637 + t580 * t604 - t581 * t732 + t587 * t715;
t501 = 0.1e1 / t503;
t492 = 0.1e1 / t494;
t491 = t679 * t513;
t489 = t680 * t513;
t487 = t674 * t537;
t485 = t491 * t685 - t511 * t543 + t512 * t559;
t482 = t679 * t721 + (t559 * t703 - t510 * t546 + (t525 * t543 - t527 * t535 - t559 * t766) * t547) * t513;
t480 = t680 * t721 + (t550 * t703 - t507 * t546 + (-t525 * t536 - t526 * t535 - t550 * t766) * t547) * t513;
t1 = [0.2e1 * t537 * t718 + (-t504 * t546 + t537 * t760) * t513, t482, 0, t480, 0, 0; t535 * t707 + (t766 * t497 + (-t483 * t535 - t487 * t504) * t498) * t492 + (t487 * t706 + (t487 * t724 - t674 * t742 + (-(t488 * t513 * t736 + t721) * t740 - (t718 * t757 - t488 + (t488 - t682) * t513) * t739) * t498) * t492) * t537 (t485 * t741 - t497 * t544) * t725 + ((-t562 * t672 + t563 * t635 - t580 * t676 + t581 * t731 + t587 * t661 + t601 * t665) * t497 + t485 * t705 + (-t544 * t483 - t485 * t504 - (t482 * t535 + t491 * t766 + t527 + (-t491 * t549 - t543) * t488) * t739 - (-t482 * t549 - t491 * t525 - t510 + (-t491 * t535 - t559) * t488) * t740) * t498) * t492, 0, t644 * t707 + (t505 * t497 - t644 * t761 - ((t480 * t535 + t489 * t766 + t526 + (-t489 * t549 + t536) * t488) * t512 + (-t480 * t549 - t489 * t525 - t507 + (-t489 * t535 - t550) * t488) * t511) * t741) * t492 + (t537 * t706 + (-t742 + t705) * t492) * (t489 * t685 + t511 * t536 + t512 * t550) 0, 0; (-t518 * t521 + t522 * t738) * t722 + ((qJD(5) * t522 - t507 * t634 - t541 * t636) * t518 + t522 * t704 + (-t521 * t496 - (-qJD(5) * t521 - t507 * t636 + t541 * t634) * t523 - t522 * t495) * t519) * t501 (-t518 * t528 + t529 * t738) * t722 + ((qJD(5) * t529 + t509 * t634 - t542 * t636) * t518 + t529 * t704 + (-t528 * t496 - (-qJD(5) * t528 + t509 * t636 + t542 * t634) * t523 - t529 * t495) * t519) * t501, 0, t681 * t537 * t723 + (t681 * t504 + ((-qJD(5) * t518 - 0.2e1 * t717) * t636 + (t495 * t636 + (t496 - t730) * t634) * t519) * t537) * t501, t723 + 0.2e1 * (t495 * t501 * t519 + (-t501 * t743 - t519 * t745) * t523) * t523, 0;];
JaD_rot  = t1;
