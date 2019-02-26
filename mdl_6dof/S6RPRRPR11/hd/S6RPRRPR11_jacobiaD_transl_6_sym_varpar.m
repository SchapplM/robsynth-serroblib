% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR11
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
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR11_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR11_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:48
% EndTime: 2019-02-26 21:06:49
% DurationCPUTime: 1.08s
% Computational Cost: add. (1421->127), mult. (4224->209), div. (0->0), fcn. (4706->16), ass. (0->89)
t598 = sin(qJ(3));
t595 = sin(pkin(6));
t658 = sin(pkin(7));
t649 = t595 * t658;
t667 = cos(qJ(1));
t627 = t667 * t649;
t659 = cos(pkin(12));
t661 = cos(pkin(6));
t639 = t661 * t659;
t657 = sin(pkin(12));
t665 = sin(qJ(1));
t574 = -t667 * t639 + t665 * t657;
t660 = cos(pkin(7));
t648 = t660 * t574;
t637 = t661 * t657;
t575 = t667 * t637 + t665 * t659;
t666 = cos(qJ(3));
t652 = t575 * t666;
t554 = (t627 + t648) * t598 - t652;
t651 = t595 * t667;
t565 = t574 * t658 - t660 * t651;
t597 = sin(qJ(4));
t599 = cos(qJ(4));
t546 = t554 * t599 - t565 * t597;
t615 = t666 * t627;
t673 = -t666 * t648 - t615;
t551 = t575 * t598 - t673;
t593 = pkin(13) + qJ(6);
t591 = sin(t593);
t592 = cos(t593);
t682 = -t546 * t591 - t551 * t592;
t681 = t546 * t592 - t551 * t591;
t576 = -t665 * t637 + t667 * t659;
t572 = t576 * qJD(1);
t625 = t665 * t649;
t607 = t665 * t639 + t667 * t657;
t653 = t607 * qJD(1);
t670 = t653 * t660;
t541 = -qJD(3) * t673 - t572 * t666 + t598 * (-qJD(1) * t625 + t575 * qJD(3) + t670);
t645 = t660 * t665;
t626 = t595 * t645;
t633 = t653 * t658;
t606 = qJD(1) * t626 + t633;
t680 = t546 * qJD(4) + t541 * t597 + t606 * t599;
t632 = t554 * t597 + t565 * t599;
t533 = t632 * qJD(4) - t541 * t599 + t606 * t597;
t614 = t666 * t625;
t621 = t598 * t627;
t538 = t666 * t670 - qJD(1) * t614 + t572 * t598 - (t598 * t648 + t621 - t652) * qJD(3);
t671 = pkin(10) + sin(pkin(13)) * pkin(5);
t663 = r_i_i_C(3) + pkin(11) + qJ(5);
t605 = t607 * t660;
t556 = t576 * t666 + (-t605 + t625) * t598;
t566 = t607 * t658 + t626;
t547 = -t556 * t597 + t566 * t599;
t641 = t591 * r_i_i_C(1) + t592 * r_i_i_C(2);
t622 = qJD(6) * t641;
t669 = -t576 * t598 - t666 * t605 + t614;
t636 = t660 * t659;
t638 = t661 * t658;
t562 = -t666 * t638 + (t598 * t657 - t636 * t666) * t595;
t601 = t565 * qJD(1);
t650 = qJD(1) * t665;
t647 = qJD(1) * t651;
t642 = t592 * r_i_i_C(1) - t591 * r_i_i_C(2);
t548 = t556 * t599 + t566 * t597;
t563 = t598 * t638 + (t598 * t636 + t657 * t666) * t595;
t573 = -t659 * t649 + t661 * t660;
t550 = t563 * t599 + t573 * t597;
t630 = -t563 * t597 + t573 * t599;
t590 = cos(pkin(13)) * pkin(5) + pkin(4);
t624 = t590 + t642;
t623 = qJD(6) * t642;
t616 = t641 + t671;
t608 = -t663 * t597 - t624 * t599 - pkin(3);
t603 = qJD(1) * t648;
t600 = -t597 * qJD(5) + t599 * t622 + (t624 * t597 - t663 * t599) * qJD(4);
t571 = t575 * qJD(1);
t560 = t563 * qJD(3);
t559 = t562 * qJD(3);
t543 = t630 * qJD(4) - t559 * t599;
t542 = t550 * qJD(4) - t559 * t597;
t537 = qJD(1) * t621 + t669 * qJD(3) - t571 * t666 + t598 * t603;
t536 = -qJD(1) * t615 + t556 * qJD(3) - t571 * t598 - t666 * t603;
t531 = t547 * qJD(4) + t537 * t599 - t597 * t601;
t530 = t548 * qJD(4) + t537 * t597 + t599 * t601;
t529 = t531 * t592 + t536 * t591 + (-t548 * t591 - t592 * t669) * qJD(6);
t528 = -t531 * t591 + t536 * t592 + (-t548 * t592 + t591 * t669) * qJD(6);
t1 = [t632 * qJD(5) + t541 * pkin(3) - t572 * pkin(2) - pkin(9) * t633 + qJD(2) * t651 - t624 * t533 - t616 * t538 + t663 * t680 + (r_i_i_C(1) * t682 - t681 * r_i_i_C(2)) * qJD(6) + (-t667 * pkin(1) + (-pkin(9) * t645 - t665 * qJ(2)) * t595) * qJD(1), t647, t608 * t536 + t616 * t537 + t556 * t623 - t600 * t669, t548 * qJD(5) - t624 * t530 + t663 * t531 - t547 * t622, t530, t528 * r_i_i_C(1) - t529 * r_i_i_C(2); t665 * t595 * qJD(2) - pkin(1) * t650 - t571 * pkin(2) + t537 * pkin(3) - pkin(9) * t601 + t529 * r_i_i_C(1) + t528 * r_i_i_C(2) + qJ(2) * t647 - t547 * qJD(5) + t663 * t530 + t531 * t590 + t671 * t536, t595 * t650, t608 * t538 - t541 * t616 + t600 * t551 - t554 * t623, -qJD(5) * t546 + t663 * t533 - t632 * t622 + t624 * t680, -t680 (-t533 * t591 + t538 * t592) * r_i_i_C(1) + (-t533 * t592 - t538 * t591) * r_i_i_C(2) + (t681 * r_i_i_C(1) + r_i_i_C(2) * t682) * qJD(6); 0, 0, -t616 * t559 + t608 * t560 + t600 * t562 + t563 * t623, t550 * qJD(5) - t624 * t542 + t663 * t543 - t630 * t622, t542 (-t543 * t591 + t560 * t592) * r_i_i_C(1) + (-t543 * t592 - t560 * t591) * r_i_i_C(2) + ((-t550 * t592 - t562 * t591) * r_i_i_C(1) + (t550 * t591 - t562 * t592) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
