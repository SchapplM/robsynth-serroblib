% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRR10
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR10_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiaD_transl_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:53:01
% EndTime: 2019-02-26 22:53:02
% DurationCPUTime: 1.56s
% Computational Cost: add. (1104->185), mult. (3507->321), div. (0->0), fcn. (3702->14), ass. (0->115)
t575 = sin(qJ(1));
t648 = cos(pkin(6));
t650 = sin(qJ(2));
t606 = t648 * t650;
t598 = t575 * t606;
t613 = t650 * qJD(2);
t578 = cos(qJ(2));
t579 = cos(qJ(1));
t626 = t579 * t578;
t545 = -qJD(1) * t598 - t575 * t613 + (qJD(2) * t648 + qJD(1)) * t626;
t574 = sin(qJ(3));
t577 = cos(qJ(3));
t557 = t575 * t578 + t579 * t606;
t612 = t578 * t648;
t582 = t575 * t612 + t579 * t650;
t544 = qJD(1) * t582 + qJD(2) * t557;
t572 = cos(pkin(7));
t569 = sin(pkin(7));
t570 = sin(pkin(6));
t625 = qJD(1) * t570;
t616 = t575 * t625;
t610 = t569 * t616;
t585 = -t544 * t572 + t610;
t633 = t570 * t579;
t620 = t569 * t633;
t556 = t575 * t650 - t579 * t612;
t642 = t556 * t572;
t594 = t620 + t642;
t640 = t557 * t577;
t517 = (t574 * t594 - t640) * qJD(3) - t545 * t574 + t585 * t577;
t645 = t544 * t569;
t535 = t572 * t616 + t645;
t571 = cos(pkin(8));
t576 = cos(qJ(4));
t631 = t571 * t576;
t568 = sin(pkin(8));
t637 = t568 * t576;
t662 = t517 * t631 + t535 * t637;
t573 = sin(qJ(4));
t632 = t571 * t573;
t638 = t568 * t573;
t661 = -t517 * t632 - t535 * t638;
t529 = t557 * t574 + t577 * t594;
t630 = t572 * t574;
t530 = t556 * t630 + t574 * t620 - t640;
t548 = -t556 * t569 + t572 * t633;
t660 = t529 * t631 - t530 * t573 + t548 * t637;
t659 = t529 * t632 + t530 * t576 + t548 * t638;
t654 = pkin(12) + r_i_i_C(3);
t581 = t598 - t626;
t635 = t570 * t575;
t592 = t569 * t635 - t572 * t582;
t532 = t574 * t592 - t577 * t581;
t649 = t569 * pkin(11);
t542 = qJD(1) * t556 + qJD(2) * t581;
t646 = t542 * t569;
t644 = t544 * t574;
t636 = t569 * t571;
t634 = t570 * t578;
t629 = t572 * t577;
t628 = t574 * t578;
t627 = t577 * t578;
t624 = qJD(3) * t574;
t623 = qJD(3) * t577;
t622 = t569 * t638;
t621 = t569 * t637;
t619 = pkin(11) * t572 + pkin(10);
t618 = t650 * t574;
t617 = t650 * t577;
t615 = t579 * t625;
t614 = qJD(2) * t634;
t611 = t648 * t569;
t609 = t569 * t614;
t608 = t570 * t613;
t607 = t574 * t611;
t605 = t650 * t568 * t569 * t570;
t543 = qJD(1) * t557 + qJD(2) * t582;
t586 = t542 * t572 + t569 * t615;
t515 = -t532 * qJD(3) + t543 * t574 + t586 * t577;
t533 = t572 * t615 - t646;
t603 = t515 * t571 + t533 * t568;
t531 = t574 * t581 + t577 * t592;
t600 = t531 * t571 + (t569 * t582 + t572 * t635) * t568;
t599 = t576 * r_i_i_C(1) - t573 * r_i_i_C(2) + pkin(3);
t538 = t556 * t574 - t557 * t629;
t595 = t556 * t577 + t557 * t630;
t540 = t574 * t582 + t581 * t629;
t593 = t577 * t582 - t581 * t630;
t591 = qJD(2) * t605;
t590 = t572 * t628 + t617;
t589 = t572 * t627 - t618;
t588 = t572 * t618 - t627;
t587 = -t572 * t617 - t628;
t584 = t587 * t570 * t571 + t605;
t536 = (-qJD(2) * t589 + qJD(3) * t588) * t570;
t583 = t536 * t571 + t568 * t609;
t580 = (-t573 * r_i_i_C(1) - t576 * r_i_i_C(2)) * t571 + t654 * t568;
t561 = t620 * t623;
t555 = -t569 * t634 + t572 * t648;
t552 = t588 * t570;
t547 = t570 * t590 + t607;
t546 = t570 * t589 + t577 * t611;
t537 = (-qJD(2) * t590 + qJD(3) * t587) * t570;
t526 = -t608 * t630 - t570 * qJD(3) * t618 + (t614 + (t572 * t634 + t611) * qJD(3)) * t577;
t525 = -qJD(3) * t607 + (qJD(2) * t587 - qJD(3) * t590) * t570;
t524 = qJD(3) * t538 - t544 * t577 - t545 * t630;
t523 = qJD(3) * t595 - t545 * t629 + t644;
t522 = qJD(3) * t540 + t542 * t577 + t543 * t630;
t521 = qJD(3) * t593 - t542 * t574 + t543 * t629;
t520 = t561 + (qJD(3) * t642 - t545) * t577 + (qJD(3) * t557 - t585) * t574;
t518 = t574 * t610 + t545 * t577 - t557 * t624 - t561 + (-t556 * t623 - t644) * t572;
t516 = t581 * t624 + t586 * t574 + (qJD(3) * t592 - t543) * t577;
t514 = t516 * t576 + t603 * t573 + (-t532 * t573 + t576 * t600) * qJD(4);
t513 = -t516 * t573 + t603 * t576 + (-t532 * t576 - t573 * t600) * qJD(4);
t1 = [(t520 * t576 + t661) * r_i_i_C(1) + (-t520 * t573 - t662) * r_i_i_C(2) + t520 * pkin(3) - t545 * pkin(2) - pkin(11) * t645 + (-t579 * pkin(1) - t619 * t635) * qJD(1) + (t660 * r_i_i_C(1) - t659 * r_i_i_C(2)) * qJD(4) + t654 * (t517 * t568 - t535 * t571) (t521 * t632 + t522 * t576 - t543 * t622) * r_i_i_C(1) + (t521 * t631 - t522 * t573 - t543 * t621) * r_i_i_C(2) + t522 * pkin(3) + t542 * pkin(2) - t543 * t649 + ((t540 * t631 + t573 * t593 - t581 * t621) * r_i_i_C(1) + (-t540 * t632 + t576 * t593 + t581 * t622) * r_i_i_C(2)) * qJD(4) + t654 * (-t521 * t568 - t543 * t636) t599 * t515 + t580 * t516 + ((-t531 * t573 - t532 * t631) * r_i_i_C(1) + (-t531 * t576 + t532 * t632) * r_i_i_C(2)) * qJD(4), t513 * r_i_i_C(1) - t514 * r_i_i_C(2), 0, 0; t514 * r_i_i_C(1) + t513 * r_i_i_C(2) + t516 * pkin(3) - t543 * pkin(2) - pkin(11) * t646 + (-t575 * pkin(1) + t619 * t633) * qJD(1) + t654 * (-t515 * t568 + t533 * t571) (t523 * t632 + t524 * t576 + t545 * t622) * r_i_i_C(1) + (t523 * t631 - t524 * t573 + t545 * t621) * r_i_i_C(2) + t524 * pkin(3) - t544 * pkin(2) + t545 * t649 + ((t538 * t631 + t557 * t621 + t573 * t595) * r_i_i_C(1) + (-t538 * t632 - t557 * t622 + t576 * t595) * r_i_i_C(2)) * qJD(4) + t654 * (-t523 * t568 + t545 * t636) t599 * t517 + t580 * t518 + ((t529 * t573 + t530 * t631) * r_i_i_C(1) + (t529 * t576 - t530 * t632) * r_i_i_C(2)) * qJD(4) (-t518 * t573 + t662) * r_i_i_C(1) + (-t518 * t576 + t661) * r_i_i_C(2) + (t659 * r_i_i_C(1) + t660 * r_i_i_C(2)) * qJD(4), 0, 0; 0 (t537 * t576 + t583 * t573 + (t552 * t573 + t576 * t584) * qJD(4)) * r_i_i_C(1) + (-t537 * t573 + t583 * t576 + (t552 * t576 - t573 * t584) * qJD(4)) * r_i_i_C(2) + t537 * pkin(3) - pkin(2) * t608 + pkin(11) * t609 + t654 * (-t536 * t568 + t571 * t609) t599 * t525 + t580 * t526 + ((-t546 * t573 - t547 * t631) * r_i_i_C(1) + (-t546 * t576 + t547 * t632) * r_i_i_C(2)) * qJD(4) (t525 * t631 - t526 * t573 + t576 * t591) * r_i_i_C(1) + (-t525 * t632 - t526 * t576 - t573 * t591) * r_i_i_C(2) + ((-t546 * t632 - t547 * t576 - t555 * t638) * r_i_i_C(1) + (-t546 * t631 + t547 * t573 - t555 * t637) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
