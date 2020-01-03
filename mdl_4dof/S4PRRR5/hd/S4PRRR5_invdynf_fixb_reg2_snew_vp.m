% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRRR5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR5_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:54
% EndTime: 2019-12-31 16:33:55
% DurationCPUTime: 0.80s
% Computational Cost: add. (1619->96), mult. (2364->137), div. (0->0), fcn. (1646->8), ass. (0->78)
t633 = qJD(2) + qJD(3);
t631 = t633 ^ 2;
t632 = qJDD(2) + qJDD(3);
t640 = sin(qJ(3));
t643 = cos(qJ(3));
t611 = t640 * t631 - t643 * t632;
t641 = sin(qJ(2));
t644 = cos(qJ(2));
t648 = -t643 * t631 - t640 * t632;
t656 = t641 * t611 + t644 * t648;
t594 = t644 * t611 - t641 * t648;
t637 = sin(pkin(7));
t638 = cos(pkin(7));
t620 = t637 * g(1) - t638 * g(2);
t653 = t637 * t620;
t639 = sin(qJ(4));
t652 = t639 * t632;
t642 = cos(qJ(4));
t651 = t642 * t632;
t621 = -t638 * g(1) - t637 * g(2);
t636 = -g(3) + qJDD(1);
t605 = t644 * t621 + t641 * t636;
t646 = qJD(2) ^ 2;
t603 = -t646 * pkin(2) + t605;
t604 = -t641 * t621 + t644 * t636;
t647 = qJDD(2) * pkin(2) + t604;
t583 = t643 * t603 + t640 * t647;
t634 = t639 ^ 2;
t635 = t642 ^ 2;
t650 = t634 + t635;
t649 = qJD(4) * t633;
t582 = -t640 * t603 + t643 * t647;
t645 = qJD(4) ^ 2;
t624 = t639 * t631 * t642;
t623 = -t635 * t631 - t645;
t622 = -t634 * t631 - t645;
t619 = t644 * qJDD(2) - t641 * t646;
t618 = -t641 * qJDD(2) - t644 * t646;
t617 = -qJDD(4) + t624;
t616 = qJDD(4) + t624;
t614 = t638 * t620;
t613 = t650 * t631;
t608 = t650 * t632;
t607 = -0.2e1 * t639 * t649 + t651;
t606 = 0.2e1 * t642 * t649 + t652;
t599 = t642 * t617 - t639 * t622;
t598 = -t639 * t616 + t642 * t623;
t597 = t639 * t617 + t642 * t622;
t596 = t642 * t616 + t639 * t623;
t593 = t643 * t608 - t640 * t613;
t590 = t640 * t608 + t643 * t613;
t589 = -t641 * t604 + t644 * t605;
t588 = t644 * t604 + t641 * t605;
t587 = t643 * t599 + t640 * t606;
t586 = t643 * t598 - t640 * t607;
t585 = t640 * t599 - t643 * t606;
t584 = t640 * t598 + t643 * t607;
t581 = -t631 * pkin(3) + t632 * pkin(6) + t583;
t580 = -t632 * pkin(3) - t631 * pkin(6) - t582;
t579 = -t641 * t590 + t644 * t593;
t578 = t644 * t590 + t641 * t593;
t577 = t642 * t581 - t639 * t620;
t576 = -t639 * t581 - t642 * t620;
t575 = -t641 * t585 + t644 * t587;
t574 = -t641 * t584 + t644 * t586;
t573 = t644 * t585 + t641 * t587;
t572 = t644 * t584 + t641 * t586;
t571 = -t640 * t582 + t643 * t583;
t570 = t643 * t582 + t640 * t583;
t569 = -t639 * t576 + t642 * t577;
t568 = t642 * t576 + t639 * t577;
t567 = -t641 * t570 + t644 * t571;
t566 = t644 * t570 + t641 * t571;
t565 = t643 * t569 + t640 * t580;
t564 = t640 * t569 - t643 * t580;
t563 = -t641 * t564 + t644 * t565;
t562 = t644 * t564 + t641 * t565;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t638 * t621 - t653, 0, 0, 0, 0, 0, 0, t638 * t618, -t638 * t619, 0, t638 * t589 - t653, 0, 0, 0, 0, 0, 0, t638 * t656, t638 * t594, 0, t638 * t567 - t653, 0, 0, 0, 0, 0, 0, t638 * t574 + t637 * t596, t638 * t575 + t637 * t597, t638 * t579, t638 * t563 + t637 * t568; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t637 * t621 + t614, 0, 0, 0, 0, 0, 0, t637 * t618, -t637 * t619, 0, t637 * t589 + t614, 0, 0, 0, 0, 0, 0, t637 * t656, t637 * t594, 0, t637 * t567 + t614, 0, 0, 0, 0, 0, 0, t637 * t574 - t638 * t596, t637 * t575 - t638 * t597, t637 * t579, t637 * t563 - t638 * t568; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t636, 0, 0, 0, 0, 0, 0, t619, t618, 0, t588, 0, 0, 0, 0, 0, 0, -t594, t656, 0, t566, 0, 0, 0, 0, 0, 0, t572, t573, t578, t562; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t621, 0, 0, 0, 0, 0, 0, t618, -t619, 0, t589, 0, 0, 0, 0, 0, 0, t656, t594, 0, t567, 0, 0, 0, 0, 0, 0, t574, t575, t579, t563; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t620, 0, 0, 0, 0, 0, 0, 0, 0, 0, t620, 0, 0, 0, 0, 0, 0, 0, 0, 0, t620, 0, 0, 0, 0, 0, 0, -t596, -t597, 0, -t568; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t636, 0, 0, 0, 0, 0, 0, t619, t618, 0, t588, 0, 0, 0, 0, 0, 0, -t594, t656, 0, t566, 0, 0, 0, 0, 0, 0, t572, t573, t578, t562; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t646, -qJDD(2), 0, t605, 0, 0, 0, 0, 0, 0, t648, t611, 0, t571, 0, 0, 0, 0, 0, 0, t586, t587, t593, t565; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t646, 0, t604, 0, 0, 0, 0, 0, 0, -t611, t648, 0, t570, 0, 0, 0, 0, 0, 0, t584, t585, t590, t564; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t620, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t620, 0, 0, 0, 0, 0, 0, t596, t597, 0, t568; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t631, -t632, 0, t583, 0, 0, 0, 0, 0, 0, t598, t599, t608, t569; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t632, -t631, 0, t582, 0, 0, 0, 0, 0, 0, t607, -t606, t613, -t580; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t620, 0, 0, 0, 0, 0, 0, t596, t597, 0, t568; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t623, t617, t651, t577; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t616, t622, -t652, t576; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t607, t606, -t613, t580;];
f_new_reg = t1;
