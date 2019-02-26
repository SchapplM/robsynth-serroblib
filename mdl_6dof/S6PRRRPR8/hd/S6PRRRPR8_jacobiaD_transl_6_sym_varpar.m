% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:40
% EndTime: 2019-02-26 20:14:42
% DurationCPUTime: 1.30s
% Computational Cost: add. (1414->171), mult. (4554->305), div. (0->0), fcn. (5043->14), ass. (0->115)
t570 = sin(qJ(2));
t573 = cos(qJ(2));
t566 = cos(pkin(12));
t647 = cos(pkin(6));
t625 = t566 * t647;
t645 = sin(pkin(12));
t593 = t645 * t570 - t573 * t625;
t646 = cos(pkin(7));
t580 = t593 * t646;
t564 = sin(pkin(7));
t565 = sin(pkin(6));
t641 = t565 * t566;
t630 = t564 * t641;
t659 = -t580 - t630;
t611 = t647 * t645;
t590 = t566 * t570 + t573 * t611;
t626 = t565 * t645;
t658 = t564 * t626 - t590 * t646;
t550 = t593 * qJD(2);
t569 = sin(qJ(3));
t592 = -t570 * t625 - t645 * t573;
t587 = t592 * qJD(2);
t578 = t646 * t587;
t638 = qJD(3) * t569;
t649 = cos(qJ(3));
t655 = t659 * t649;
t510 = t655 * qJD(3) - t550 * t649 + t569 * t578 + t592 * t638;
t542 = t593 * t564 - t646 * t641;
t657 = t542 * qJD(4) + t510;
t551 = t590 * qJD(2);
t591 = -t566 * t573 + t570 * t611;
t586 = t591 * qJD(2);
t577 = t646 * t586;
t654 = t658 * t649;
t512 = t654 * qJD(3) - t551 * t649 + t569 * t577 + t591 * t638;
t543 = t590 * t564 + t646 * t626;
t656 = t543 * qJD(4) + t512;
t530 = t658 * t569 - t591 * t649;
t653 = t530 * qJD(4) + t564 * t586;
t629 = t592 * t649;
t528 = t659 * t569 - t629;
t652 = t528 * qJD(4) + t564 * t587;
t633 = r_i_i_C(3) + pkin(11) + pkin(4);
t567 = sin(qJ(6));
t571 = cos(qJ(6));
t613 = t571 * r_i_i_C(1) - t567 * r_i_i_C(2);
t651 = pkin(5) + pkin(10) + t613;
t597 = t613 * qJD(6) + qJD(5);
t648 = pkin(9) * t564;
t568 = sin(qJ(4));
t644 = t564 * t568;
t572 = cos(qJ(4));
t643 = t564 * t572;
t642 = t565 * t564;
t640 = t569 * t570;
t639 = qJD(2) * t565;
t632 = t570 * t642;
t631 = t565 * t640;
t628 = t649 * t573;
t627 = t564 * t639;
t624 = t569 * t646;
t623 = t647 * t564;
t621 = t646 * qJD(3);
t620 = t570 * t627;
t619 = t573 * t627;
t617 = t569 * t623;
t615 = t646 * t649;
t614 = t569 * t621;
t612 = -t567 * r_i_i_C(1) - t571 * r_i_i_C(2);
t595 = t649 * t570 + t573 * t624;
t541 = t595 * t565 + t617;
t552 = -t573 * t642 + t647 * t646;
t609 = t541 * t572 + t552 * t568;
t531 = t541 * t568 - t552 * t572;
t608 = t649 * t623;
t607 = t573 * t615;
t606 = qJ(5) - t612;
t605 = qJD(3) * t615;
t584 = t593 * t649;
t536 = t592 * t624 - t584;
t603 = -t536 * t568 - t592 * t643;
t585 = t590 * t649;
t538 = t591 * t624 - t585;
t602 = -t538 * t568 - t591 * t643;
t601 = qJD(6) * t612;
t600 = t565 * t607;
t596 = -t570 * t624 + t628;
t549 = t596 * t565;
t598 = -t549 * t568 + t572 * t632;
t594 = t569 * t573 + t570 * t615;
t589 = t590 * t569;
t588 = t593 * t569;
t579 = -t606 * t568 - t633 * t572 - pkin(3);
t574 = -t597 * t568 + (t633 * t568 - t606 * t572) * qJD(4);
t548 = t594 * t565;
t540 = -t600 - t608 + t631;
t537 = -t591 * t615 - t589;
t535 = -t592 * t615 - t588;
t534 = (-t595 * qJD(2) - t594 * qJD(3)) * t565;
t529 = -t569 * t591 - t654;
t527 = -t569 * t592 - t655;
t526 = qJD(3) * t608 + ((t607 - t640) * qJD(3) + t596 * qJD(2)) * t565;
t525 = qJD(3) * t617 + (t594 * qJD(2) + t595 * qJD(3)) * t565;
t521 = t530 * t568 - t543 * t572;
t519 = t528 * t568 - t542 * t572;
t518 = qJD(3) * t589 + t551 * t624 + t649 * t586 + t591 * t605;
t516 = qJD(3) * t588 + t550 * t624 + t649 * t587 + t592 * t605;
t511 = t530 * qJD(3) - t551 * t569 - t649 * t577;
t509 = -t550 * t569 - t649 * t578 - t630 * t638 + (-t569 * t580 - t629) * qJD(3);
t507 = t609 * qJD(4) + t526 * t568 - t572 * t620;
t505 = t551 * t643 + t518 * t568 + (t538 * t572 - t591 * t644) * qJD(4);
t503 = t550 * t643 + t516 * t568 + (t536 * t572 - t592 * t644) * qJD(4);
t501 = t656 * t568 + t653 * t572;
t499 = t657 * t568 + t652 * t572;
t1 = [0 (t505 * t567 + (-t537 * t567 - t571 * t602) * qJD(6)) * r_i_i_C(1) + (t505 * t571 + (-t537 * t571 + t567 * t602) * qJD(6)) * r_i_i_C(2) + t505 * qJ(5) - t602 * qJD(5) + t518 * pkin(3) + pkin(2) * t586 - t551 * t648 + t651 * (-qJD(3) * t585 - t551 * t615 + t569 * t586 + t591 * t614) + t633 * (t602 * qJD(4) + t518 * t572 - t551 * t644) t579 * t511 + t512 * t651 + t574 * t529 + t530 * t601, t597 * (t530 * t572 + t543 * t568) + t606 * (-t653 * t568 + t656 * t572) - t633 * t501, t501 (t501 * t571 - t511 * t567) * r_i_i_C(1) + (-t501 * t567 - t511 * t571) * r_i_i_C(2) + ((-t521 * t567 - t529 * t571) * r_i_i_C(1) + (-t521 * t571 + t529 * t567) * r_i_i_C(2)) * qJD(6); 0 (t503 * t567 + (-t535 * t567 - t571 * t603) * qJD(6)) * r_i_i_C(1) + (t503 * t571 + (-t535 * t571 + t567 * t603) * qJD(6)) * r_i_i_C(2) + t503 * qJ(5) - t603 * qJD(5) + t516 * pkin(3) + pkin(2) * t587 - t550 * t648 + t651 * (-qJD(3) * t584 - t550 * t615 + t569 * t587 + t592 * t614) + t633 * (t603 * qJD(4) + t516 * t572 - t550 * t644) t579 * t509 + t510 * t651 + t574 * t527 + t528 * t601, t597 * (t528 * t572 + t542 * t568) + t606 * (-t652 * t568 + t657 * t572) - t633 * t499, t499 (t499 * t571 - t509 * t567) * r_i_i_C(1) + (-t499 * t567 - t509 * t571) * r_i_i_C(2) + ((-t519 * t567 - t527 * t571) * r_i_i_C(1) + (-t519 * t571 + t527 * t567) * r_i_i_C(2)) * qJD(6); 0, t534 * pkin(3) - t598 * qJD(5) + t606 * (-t572 * t619 + t534 * t568 + (t549 * t572 + t568 * t632) * qJD(4)) - t651 * (-qJD(2) * t600 - t565 * qJD(3) * t628 + (t621 + qJD(2)) * t631) + ((-t548 * t567 - t571 * t598) * r_i_i_C(1) + (-t548 * t571 + t567 * t598) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t570 + t573 * t648) * t639 + t633 * (t598 * qJD(4) + t534 * t572 + t568 * t619) t579 * t525 + t526 * t651 + t574 * t540 + t541 * t601, t597 * t609 + t606 * (-t531 * qJD(4) + t526 * t572 + t568 * t620) - t633 * t507, t507 (t507 * t571 - t525 * t567) * r_i_i_C(1) + (-t507 * t567 - t525 * t571) * r_i_i_C(2) + ((-t531 * t567 - t540 * t571) * r_i_i_C(1) + (-t531 * t571 + t540 * t567) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
