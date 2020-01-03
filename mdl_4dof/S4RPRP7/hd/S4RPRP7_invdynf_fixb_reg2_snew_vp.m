% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPRP7
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPRP7_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP7_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:23
% EndTime: 2019-12-31 16:47:24
% DurationCPUTime: 0.72s
% Computational Cost: add. (628->122), mult. (1297->99), div. (0->0), fcn. (612->4), ass. (0->61)
t591 = sin(qJ(3));
t593 = cos(qJ(3));
t596 = qJD(1) ^ 2;
t603 = t591 * t596 * t593;
t574 = qJDD(3) + t603;
t590 = t593 ^ 2;
t595 = qJD(3) ^ 2;
t579 = t590 * t596 + t595;
t557 = t591 * t574 + t593 * t579;
t607 = qJD(1) * qJD(3);
t602 = t591 * t607;
t605 = t593 * qJDD(1);
t569 = -0.2e1 * t602 + t605;
t592 = sin(qJ(1));
t594 = cos(qJ(1));
t617 = t592 * t557 - t594 * t569;
t616 = t594 * t557 + t592 * t569;
t615 = -2 * qJD(2);
t614 = 2 * qJD(4);
t613 = -pkin(5) - pkin(1);
t612 = t591 * g(3);
t589 = t591 ^ 2;
t609 = t589 + t590;
t600 = t591 * pkin(3) - t593 * qJ(4);
t608 = t596 * t600;
t606 = t591 * qJDD(1);
t604 = qJD(1) * t615;
t601 = t593 * t607;
t576 = t592 * g(1) - t594 * g(2);
t598 = -t596 * qJ(2) + qJDD(2) - t576;
t564 = t613 * qJDD(1) + t598;
t553 = -t593 * g(3) + t591 * t564;
t577 = -t594 * g(1) - t592 * g(2);
t560 = t593 * t574 - t591 * t579;
t599 = -qJDD(1) * qJ(2) - t577;
t597 = -t613 * t596 + t599;
t578 = -t589 * t596 - t595;
t575 = qJDD(3) - t603;
t573 = t609 * t596;
t572 = t592 * qJDD(1) + t594 * t596;
t571 = t594 * qJDD(1) - t592 * t596;
t570 = t609 * qJDD(1);
t568 = 0.2e1 * t601 + t606;
t566 = qJDD(1) * pkin(1) - t598;
t565 = t596 * pkin(1) + t599 + t604;
t563 = t597 + t604;
t559 = -t591 * t575 + t593 * t578;
t556 = t593 * t575 + t591 * t578;
t555 = -t592 * t570 - t594 * t573;
t554 = t594 * t570 - t592 * t573;
t552 = t593 * t564 + t612;
t551 = t592 * t556 + t594 * t568;
t550 = -t594 * t556 + t592 * t568;
t549 = qJDD(3) * pkin(3) + t612 + t595 * qJ(4) - qJDD(4) + (t564 - t608) * t593;
t548 = -t595 * pkin(3) + qJDD(3) * qJ(4) + (qJD(3) * t614) - t591 * t608 + t553;
t547 = -t591 * t552 + t593 * t553;
t546 = t593 * t552 + t591 * t553;
t545 = -qJ(4) * t602 - pkin(3) * t601 - t600 * qJDD(1) + (t593 * t614 + t615 + (-pkin(3) * t593 - qJ(4) * t591) * qJD(3)) * qJD(1) + t597;
t544 = t593 * t548 - t591 * t549;
t543 = t591 * t548 + t593 * t549;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t572, -t571, 0, -t592 * t576 + t594 * t577, 0, 0, 0, 0, 0, 0, 0, t572, t571, -t594 * t565 - t592 * t566, 0, 0, 0, 0, 0, 0, t551, -t617, t555, t592 * t546 - t594 * t563, 0, 0, 0, 0, 0, 0, t551, t555, t617, t592 * t543 - t594 * t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t571, -t572, 0, t594 * t576 + t592 * t577, 0, 0, 0, 0, 0, 0, 0, -t571, t572, -t592 * t565 + t594 * t566, 0, 0, 0, 0, 0, 0, t550, t616, t554, -t594 * t546 - t592 * t563, 0, 0, 0, 0, 0, 0, t550, t554, -t616, -t594 * t543 - t592 * t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t559, -t560, 0, t547, 0, 0, 0, 0, 0, 0, t559, 0, t560, t544; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t596, -qJDD(1), 0, t577, 0, 0, 0, 0, 0, 0, 0, t596, qJDD(1), -t565, 0, 0, 0, 0, 0, 0, t568, t569, -t573, -t563, 0, 0, 0, 0, 0, 0, t568, -t573, -t569, -t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t596, 0, t576, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), t596, t566, 0, 0, 0, 0, 0, 0, -t556, t557, t570, -t546, 0, 0, 0, 0, 0, 0, -t556, t570, -t557, -t543; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t559, -t560, 0, t547, 0, 0, 0, 0, 0, 0, t559, 0, t560, t544; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t559, -t560, 0, t547, 0, 0, 0, 0, 0, 0, t559, 0, t560, t544; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t596, -qJDD(1), t565, 0, 0, 0, 0, 0, 0, -t568, -t569, t573, t563, 0, 0, 0, 0, 0, 0, -t568, t573, t569, t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t596, -t566, 0, 0, 0, 0, 0, 0, t556, -t557, -t570, t546, 0, 0, 0, 0, 0, 0, t556, -t570, t557, t543; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t578, -t574, -t606, t553, 0, 0, 0, 0, 0, 0, t578, -t606, t574, t548; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t575, -t579, -t605, t552, 0, 0, 0, 0, 0, 0, t575, -t605, t579, t549; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568, t569, -t573, -t563, 0, 0, 0, 0, 0, 0, t568, -t573, -t569, -t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t578, -t606, t574, t548; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t568, -t573, -t569, -t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t575, t605, -t579, -t549;];
f_new_reg = t1;
