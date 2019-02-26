% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR10_jacobigD_rot_6_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_6_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_6_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobigD_rot_6_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:22
% EndTime: 2018-11-23 11:27:22
% DurationCPUTime: 0.41s
% Computational Cost: add. (1304->139), mult. (1583->216), div. (0->0), fcn. (1294->28), ass. (0->121)
t666 = qJD(2) / 0.2e1;
t665 = qJD(3) / 0.2e1;
t664 = qJD(4) / 0.2e1;
t627 = pkin(6) + qJ(2);
t621 = cos(t627);
t604 = t621 * t666;
t628 = pkin(6) - qJ(2);
t622 = cos(t628);
t657 = qJD(2) * t622;
t586 = t604 - t657 / 0.2e1;
t630 = sin(pkin(7));
t663 = t586 * t630;
t631 = sin(pkin(6));
t639 = sin(qJ(1));
t662 = t631 * t639;
t644 = cos(qJ(1));
t661 = t631 * t644;
t660 = qJD(1) * t639;
t659 = qJD(1) * t644;
t615 = sin(t627);
t658 = qJD(2) * t615;
t656 = qJD(2) * t639;
t655 = qJD(2) * t644;
t625 = pkin(7) + qJ(3);
t613 = sin(t625);
t654 = qJD(3) * t613;
t626 = pkin(7) - qJ(3);
t620 = cos(t626);
t653 = qJD(3) * t620;
t637 = sin(qJ(3));
t652 = qJD(3) * t637;
t642 = cos(qJ(3));
t651 = qJD(3) * t642;
t623 = pkin(8) + qJ(4);
t611 = sin(t623);
t650 = qJD(4) * t611;
t624 = pkin(8) - qJ(4);
t618 = cos(t624);
t649 = qJD(4) * t618;
t636 = sin(qJ(4));
t648 = qJD(4) * t636;
t641 = cos(qJ(4));
t647 = qJD(4) * t641;
t646 = t631 * t660;
t645 = t631 * t659;
t607 = t615 / 0.2e1;
t616 = sin(t628);
t592 = t607 - t616 / 0.2e1;
t643 = cos(qJ(2));
t572 = t644 * t592 + t639 * t643;
t574 = -t639 * t592 + t644 * t643;
t610 = t622 / 0.2e1;
t597 = t610 + t621 / 0.2e1;
t638 = sin(qJ(2));
t571 = t644 * t597 - t639 * t638;
t573 = -t639 * t597 - t644 * t638;
t640 = cos(qJ(5));
t635 = sin(qJ(5));
t634 = cos(pkin(6));
t633 = cos(pkin(7));
t632 = cos(pkin(8));
t629 = sin(pkin(8));
t619 = cos(t625);
t617 = cos(t623);
t614 = sin(t626);
t612 = sin(t624);
t609 = t620 / 0.2e1;
t608 = t618 / 0.2e1;
t606 = t613 / 0.2e1;
t605 = t611 / 0.2e1;
t603 = t616 * t666;
t602 = t619 * t665;
t601 = t614 * t665;
t600 = t617 * t664;
t599 = t612 * t664;
t598 = t610 - t621 / 0.2e1;
t596 = t609 - t619 / 0.2e1;
t595 = t609 + t619 / 0.2e1;
t594 = t608 - t617 / 0.2e1;
t593 = t608 + t617 / 0.2e1;
t591 = t607 + t616 / 0.2e1;
t590 = t606 - t614 / 0.2e1;
t589 = t606 + t614 / 0.2e1;
t588 = t605 - t612 / 0.2e1;
t587 = t605 + t612 / 0.2e1;
t585 = t604 + t657 / 0.2e1;
t584 = t603 - t658 / 0.2e1;
t583 = t603 + t658 / 0.2e1;
t582 = t602 - t653 / 0.2e1;
t581 = t602 + t653 / 0.2e1;
t580 = t601 - t654 / 0.2e1;
t579 = t601 + t654 / 0.2e1;
t578 = t600 - t649 / 0.2e1;
t577 = t600 + t649 / 0.2e1;
t576 = t599 - t650 / 0.2e1;
t575 = t599 + t650 / 0.2e1;
t570 = -t591 * t630 + t634 * t633;
t569 = -t573 * t630 + t633 * t662;
t568 = -t571 * t630 - t633 * t661;
t567 = qJD(1) * t574 + t644 * t585 - t638 * t656;
t566 = qJD(1) * t573 + t644 * t584 - t643 * t656;
t565 = -qJD(1) * t572 - t639 * t585 - t638 * t655;
t564 = -qJD(1) * t571 - t639 * t584 - t643 * t655;
t563 = t591 * t590 + t634 * t596 + t598 * t642;
t562 = t634 * t589 + t591 * t595 - t598 * t637;
t561 = -t566 * t630 + t633 * t646;
t560 = -t564 * t630 + t633 * t645;
t559 = t573 * t590 + t574 * t642 + t596 * t662;
t558 = t573 * t595 - t574 * t637 + t589 * t662;
t557 = t571 * t590 + t572 * t642 - t596 * t661;
t556 = t571 * t595 - t572 * t637 - t589 * t661;
t555 = t634 * t579 + t591 * t581 + t583 * t642 + t586 * t590 - t598 * t652;
t554 = t591 * t580 + t634 * t582 - t583 * t637 + t586 * t595 - t598 * t651;
t553 = -t554 * t629 - t632 * t663;
t552 = -t572 * t652 + t566 * t590 + t567 * t642 + t571 * t581 + (-t579 * t644 + t596 * t660) * t631;
t551 = -t572 * t651 + t566 * t595 - t567 * t637 + t571 * t580 + (-t582 * t644 + t589 * t660) * t631;
t550 = -t574 * t652 + t564 * t590 + t565 * t642 + t573 * t581 + (t579 * t639 + t596 * t659) * t631;
t549 = -t574 * t651 + t564 * t595 - t565 * t637 + t573 * t580 + (t582 * t639 + t589 * t659) * t631;
t548 = -t551 * t629 + t561 * t632;
t547 = -t549 * t629 + t560 * t632;
t1 = [0, t645, t560, t547, -t549 * t593 + t550 * t636 - t558 * t576 + t559 * t647 - t560 * t587 - t569 * t578 (t549 * t588 + t550 * t641 + t558 * t577 - t559 * t648 + t560 * t594 + t569 * t575) * t635 - t547 * t640 + ((t558 * t588 + t559 * t641 + t569 * t594) * t640 + (-t558 * t629 + t569 * t632) * t635) * qJD(5); 0, t646, t561, t548, -t551 * t593 + t552 * t636 - t556 * t576 + t557 * t647 - t561 * t587 - t568 * t578 (t551 * t588 + t552 * t641 + t556 * t577 - t557 * t648 + t561 * t594 + t568 * t575) * t635 - t548 * t640 + ((t556 * t588 + t557 * t641 + t568 * t594) * t640 + (-t556 * t629 + t568 * t632) * t635) * qJD(5); 0, 0, -t663, t553, -t554 * t593 + t555 * t636 - t562 * t576 + t563 * t647 - t570 * t578 + t587 * t663 (t554 * t588 + t555 * t641 + t562 * t577 - t563 * t648 + t570 * t575 - t594 * t663) * t635 - t553 * t640 + ((t562 * t588 + t563 * t641 + t570 * t594) * t640 + (-t562 * t629 + t570 * t632) * t635) * qJD(5);];
JgD_rot  = t1;
