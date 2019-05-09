% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRPR1
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRPR1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR1_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:21:48
% EndTime: 2019-05-04 19:21:49
% DurationCPUTime: 0.94s
% Computational Cost: add. (3310->89), mult. (4606->105), div. (0->0), fcn. (2880->8), ass. (0->63)
t604 = qJD(1) + qJD(2);
t600 = qJD(4) + t604;
t598 = t600 ^ 2;
t603 = qJDD(1) + qJDD(2);
t599 = qJDD(4) + t603;
t608 = sin(qJ(4));
t611 = cos(qJ(4));
t579 = t608 * t598 - t611 * t599;
t606 = sin(pkin(7));
t607 = cos(pkin(7));
t617 = -t611 * t598 - t608 * t599;
t564 = t607 * t579 - t606 * t617;
t609 = sin(qJ(2));
t612 = cos(qJ(2));
t627 = t606 * t579 + t607 * t617;
t554 = t612 * t564 - t609 * t627;
t610 = sin(qJ(1));
t613 = cos(qJ(1));
t634 = t609 * t564 + t612 * t627;
t638 = t610 * t554 + t613 * t634;
t637 = t613 * t554 - t610 * t634;
t602 = t604 ^ 2;
t585 = t606 * t602 - t607 * t603;
t616 = -t607 * t602 - t606 * t603;
t571 = t612 * t585 - t609 * t616;
t626 = t609 * t585 + t612 * t616;
t633 = t610 * t571 + t613 * t626;
t632 = t613 * t571 - t610 * t626;
t590 = t609 * t602 - t612 * t603;
t615 = -t612 * t602 - t609 * t603;
t625 = t610 * t590 + t613 * t615;
t624 = t613 * t590 - t610 * t615;
t596 = t610 * g(1) - t613 * g(2);
t592 = qJDD(1) * pkin(1) + t596;
t597 = -t613 * g(1) - t610 * g(2);
t614 = qJD(1) ^ 2;
t593 = -t614 * pkin(1) + t597;
t575 = t612 * t592 - t609 * t593;
t573 = t603 * pkin(2) + t575;
t576 = t609 * t592 + t612 * t593;
t574 = -t602 * pkin(2) + t576;
t559 = t606 * t573 + t607 * t574;
t558 = t607 * t573 - t606 * t574;
t605 = -g(3) + qJDD(3);
t595 = -t610 * qJDD(1) - t613 * t614;
t594 = t613 * qJDD(1) - t610 * t614;
t561 = -t609 * t575 + t612 * t576;
t560 = t612 * t575 + t609 * t576;
t557 = -t602 * pkin(3) + t559;
t556 = t603 * pkin(3) + t558;
t551 = -t606 * t558 + t607 * t559;
t550 = t607 * t558 + t606 * t559;
t549 = t608 * t556 + t611 * t557;
t548 = t611 * t556 - t608 * t557;
t547 = -t609 * t550 + t612 * t551;
t546 = t612 * t550 + t609 * t551;
t545 = -t608 * t548 + t611 * t549;
t544 = t611 * t548 + t608 * t549;
t543 = -t606 * t544 + t607 * t545;
t542 = t607 * t544 + t606 * t545;
t541 = -t609 * t542 + t612 * t543;
t540 = t612 * t542 + t609 * t543;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t595, -t594, 0, -t610 * t596 + t613 * t597, 0, 0, 0, 0, 0, 0, t625, t624, 0, -t610 * t560 + t613 * t561, 0, 0, 0, 0, 0, 0, t633, t632, 0, -t610 * t546 + t613 * t547, 0, 0, 0, 0, 0, 0, t638, t637, 0, -t610 * t540 + t613 * t541; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t594, t595, 0, t613 * t596 + t610 * t597, 0, 0, 0, 0, 0, 0, -t624, t625, 0, t613 * t560 + t610 * t561, 0, 0, 0, 0, 0, 0, -t632, t633, 0, t613 * t546 + t610 * t547, 0, 0, 0, 0, 0, 0, -t637, t638, 0, t613 * t540 + t610 * t541; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t605, 0, 0, 0, 0, 0, 0, 0, 0, 0, t605; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t614, -qJDD(1), 0, t597, 0, 0, 0, 0, 0, 0, t615, t590, 0, t561, 0, 0, 0, 0, 0, 0, t626, t571, 0, t547, 0, 0, 0, 0, 0, 0, t634, t554, 0, t541; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t614, 0, t596, 0, 0, 0, 0, 0, 0, -t590, t615, 0, t560, 0, 0, 0, 0, 0, 0, -t571, t626, 0, t546, 0, 0, 0, 0, 0, 0, -t554, t634, 0, t540; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t605, 0, 0, 0, 0, 0, 0, 0, 0, 0, t605; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t602, -t603, 0, t576, 0, 0, 0, 0, 0, 0, t616, t585, 0, t551, 0, 0, 0, 0, 0, 0, t627, t564, 0, t543; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t603, -t602, 0, t575, 0, 0, 0, 0, 0, 0, -t585, t616, 0, t550, 0, 0, 0, 0, 0, 0, -t564, t627, 0, t542; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t605, 0, 0, 0, 0, 0, 0, 0, 0, 0, t605; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t602, -t603, 0, t559, 0, 0, 0, 0, 0, 0, t617, t579, 0, t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t603, -t602, 0, t558, 0, 0, 0, 0, 0, 0, -t579, t617, 0, t544; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t605, 0, 0, 0, 0, 0, 0, 0, 0, 0, t605; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t598, -t599, 0, t549; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t599, -t598, 0, t548; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t605;];
f_new_reg  = t1;
