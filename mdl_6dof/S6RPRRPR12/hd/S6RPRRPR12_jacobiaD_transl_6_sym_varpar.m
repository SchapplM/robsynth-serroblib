% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR12_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR12_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:20
% EndTime: 2019-02-26 21:07:21
% DurationCPUTime: 1.10s
% Computational Cost: add. (1445->121), mult. (4618->199), div. (0->0), fcn. (5138->14), ass. (0->86)
t572 = sin(qJ(3));
t635 = sin(pkin(7));
t636 = sin(pkin(6));
t615 = t636 * t635;
t643 = cos(qJ(1));
t599 = t643 * t615;
t637 = cos(pkin(12));
t639 = cos(pkin(6));
t619 = t639 * t637;
t634 = sin(pkin(12));
t641 = sin(qJ(1));
t554 = -t643 * t619 + t641 * t634;
t638 = cos(pkin(7));
t627 = t638 * t554;
t617 = t639 * t634;
t555 = t643 * t617 + t641 * t637;
t642 = cos(qJ(3));
t628 = t555 * t642;
t534 = (t599 + t627) * t572 - t628;
t616 = t638 * t636;
t545 = t554 * t635 - t643 * t616;
t571 = sin(qJ(4));
t574 = cos(qJ(4));
t525 = t534 * t571 + t545 * t574;
t589 = t642 * t599;
t650 = -t642 * t627 - t589;
t531 = t555 * t572 - t650;
t570 = sin(qJ(6));
t573 = cos(qJ(6));
t662 = t525 * t570 - t531 * t573;
t661 = t525 * t573 + t531 * t570;
t556 = -t641 * t617 + t643 * t637;
t552 = t556 * qJD(1);
t595 = t641 * t615;
t583 = t641 * t619 + t643 * t634;
t630 = t583 * qJD(1);
t647 = t630 * t638;
t520 = -t650 * qJD(3) - t552 * t642 + (-qJD(1) * t595 + t555 * qJD(3) + t647) * t572;
t596 = t641 * t616;
t611 = t630 * t635;
t581 = qJD(1) * t596 + t611;
t660 = t525 * qJD(4) - t520 * t574 + t581 * t571;
t655 = t534 * t574 - t545 * t571;
t659 = qJD(4) * t655 + t520 * t571 + t581 * t574;
t656 = t637 * t616 + t639 * t635;
t588 = t642 * t595;
t593 = t572 * t599;
t517 = t642 * t647 - qJD(1) * t588 + t552 * t572 - (t572 * t627 + t593 - t628) * qJD(3);
t648 = pkin(5) + pkin(10);
t580 = t583 * t638;
t536 = t556 * t642 + (-t580 + t595) * t572;
t546 = t583 * t635 + t596;
t646 = -t536 * t571 + t546 * t574;
t629 = pkin(4) + pkin(11) + r_i_i_C(3);
t645 = -t556 * t572 - t642 * t580 + t588;
t614 = t636 * t634;
t542 = t572 * t614 - t642 * t656;
t575 = t545 * qJD(1);
t622 = t573 * r_i_i_C(1) - t570 * r_i_i_C(2);
t592 = t622 * qJD(6) + qJD(5);
t625 = t643 * t636;
t623 = t641 * t636;
t621 = -t570 * r_i_i_C(1) - t573 * r_i_i_C(2);
t609 = t536 * t574 + t546 * t571;
t543 = t572 * t656 + t642 * t614;
t553 = -t637 * t615 + t639 * t638;
t608 = t543 * t574 + t553 * t571;
t529 = t543 * t571 - t553 * t574;
t606 = qJ(5) - t621;
t605 = qJD(1) * t625;
t603 = t622 + t648;
t602 = qJD(6) * t621;
t582 = -t606 * t571 - t629 * t574 - pkin(3);
t578 = qJD(1) * t627;
t576 = -t592 * t571 + (t629 * t571 - t606 * t574) * qJD(4);
t551 = t555 * qJD(1);
t540 = t543 * qJD(3);
t539 = t542 * qJD(3);
t521 = t608 * qJD(4) - t539 * t571;
t516 = qJD(1) * t593 + t645 * qJD(3) - t551 * t642 + t572 * t578;
t515 = -qJD(1) * t589 + t536 * qJD(3) - t551 * t572 - t642 * t578;
t510 = t646 * qJD(4) + t516 * t574 - t571 * t575;
t509 = t609 * qJD(4) + t516 * t571 + t574 * t575;
t508 = t509 * t570 + t515 * t573 + (t570 * t645 - t573 * t646) * qJD(6);
t507 = t509 * t573 - t515 * t570 + (t570 * t646 + t573 * t645) * qJD(6);
t1 = [t525 * qJD(5) + t520 * pkin(3) - t552 * pkin(2) - pkin(9) * t611 + qJD(2) * t625 + t606 * t659 - t603 * t517 + (t661 * r_i_i_C(1) - t662 * r_i_i_C(2)) * qJD(6) - t629 * t660 + (-t643 * pkin(1) - pkin(9) * t596 - qJ(2) * t623) * qJD(1), t605, t582 * t515 + t603 * t516 + t536 * t602 - t576 * t645, -t629 * t509 + t606 * t510 + t592 * t609, t509, t507 * r_i_i_C(1) - t508 * r_i_i_C(2); -qJD(1) * t641 * pkin(1) - t551 * pkin(2) + t516 * pkin(3) - pkin(9) * t575 + t508 * r_i_i_C(1) + t507 * r_i_i_C(2) + qJ(2) * t605 + t509 * qJ(5) + qJD(2) * t623 - qJD(5) * t646 + t629 * t510 + t648 * t515, qJD(1) * t623, t582 * t517 - t520 * t603 + t576 * t531 - t534 * t602, -t592 * t655 + t606 * t660 + t629 * t659, -t659 (-t517 * t570 - t573 * t659) * r_i_i_C(1) + (-t517 * t573 + t570 * t659) * r_i_i_C(2) + (t662 * r_i_i_C(1) + t661 * r_i_i_C(2)) * qJD(6); 0, 0, -t603 * t539 + t582 * t540 + t576 * t542 + t543 * t602, t592 * t608 + t606 * (-t529 * qJD(4) - t539 * t574) - t629 * t521, t521 (t521 * t573 - t540 * t570) * r_i_i_C(1) + (-t521 * t570 - t540 * t573) * r_i_i_C(2) + ((-t529 * t570 - t542 * t573) * r_i_i_C(1) + (-t529 * t573 + t542 * t570) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
