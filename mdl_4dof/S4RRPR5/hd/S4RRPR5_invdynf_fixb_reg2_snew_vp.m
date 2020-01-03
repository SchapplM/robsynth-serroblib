% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRPR5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRPR5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR5_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:40
% EndTime: 2019-12-31 17:03:40
% DurationCPUTime: 0.73s
% Computational Cost: add. (1113->107), mult. (1502->115), div. (0->0), fcn. (848->6), ass. (0->63)
t590 = (qJD(1) + qJD(2));
t588 = t590 ^ 2;
t589 = qJDD(1) + qJDD(2);
t594 = sin(qJ(2));
t597 = cos(qJ(2));
t568 = t597 * t588 + t594 * t589;
t570 = t594 * t588 - t597 * t589;
t595 = sin(qJ(1));
t598 = cos(qJ(1));
t602 = t598 * t568 - t595 * t570;
t604 = t595 * t568 + t598 * t570;
t593 = sin(qJ(4));
t609 = t593 * t589;
t596 = cos(qJ(4));
t608 = t596 * t589;
t591 = t593 ^ 2;
t592 = t596 ^ 2;
t607 = t591 + t592;
t606 = qJD(4) * t590;
t605 = t593 * t588 * t596;
t581 = t595 * g(1) - t598 * g(2);
t573 = qJDD(1) * pkin(1) + t581;
t582 = -t598 * g(1) - t595 * g(2);
t600 = qJD(1) ^ 2;
t574 = -t600 * pkin(1) + t582;
t556 = t597 * t573 - t594 * t574;
t557 = t594 * t573 + t597 * t574;
t553 = -t589 * pkin(2) - t588 * qJ(3) + qJDD(3) - t556;
t552 = -t588 * pkin(2) + t589 * qJ(3) + (2 * qJD(3) * t590) + t557;
t599 = qJD(4) ^ 2;
t580 = -t592 * t588 - t599;
t579 = -t591 * t588 - t599;
t578 = -t595 * qJDD(1) - t598 * t600;
t577 = t598 * qJDD(1) - t595 * t600;
t576 = -qJDD(4) - t605;
t575 = qJDD(4) - t605;
t572 = t607 * t588;
t567 = t607 * t589;
t565 = -0.2e1 * t593 * t606 + t608;
t564 = 0.2e1 * t596 * t606 + t609;
t561 = t596 * t576 - t593 * t580;
t560 = -t593 * t575 + t596 * t579;
t559 = t593 * t576 + t596 * t580;
t558 = t596 * t575 + t593 * t579;
t555 = -t594 * t567 - t597 * t572;
t554 = t597 * t567 - t594 * t572;
t551 = -t589 * pkin(6) + t553;
t550 = t594 * t559 + t597 * t565;
t549 = t594 * t558 + t597 * t564;
t548 = -t597 * t559 + t594 * t565;
t547 = -t597 * t558 + t594 * t564;
t546 = -t588 * pkin(6) + t552;
t545 = -t596 * g(3) + t593 * t551;
t544 = t593 * g(3) + t596 * t551;
t543 = -t594 * t556 + t597 * t557;
t542 = t597 * t556 + t594 * t557;
t541 = t597 * t552 + t594 * t553;
t540 = t594 * t552 - t597 * t553;
t539 = -t593 * t544 + t596 * t545;
t538 = t596 * t544 + t593 * t545;
t537 = t594 * t538 + t597 * t546;
t536 = -t597 * t538 + t594 * t546;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t578, -t577, 0, -t595 * t581 + t598 * t582, 0, 0, 0, 0, 0, 0, -t602, t604, 0, -t595 * t542 + t598 * t543, 0, 0, 0, 0, 0, 0, 0, t602, -t604, -t595 * t540 + t598 * t541, 0, 0, 0, 0, 0, 0, -t595 * t547 + t598 * t549, -t595 * t548 + t598 * t550, -t595 * t554 + t598 * t555, -t595 * t536 + t598 * t537; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t577, t578, 0, t598 * t581 + t595 * t582, 0, 0, 0, 0, 0, 0, -t604, -t602, 0, t598 * t542 + t595 * t543, 0, 0, 0, 0, 0, 0, 0, t604, t602, t598 * t540 + t595 * t541, 0, 0, 0, 0, 0, 0, t598 * t547 + t595 * t549, t598 * t548 + t595 * t550, t598 * t554 + t595 * t555, t598 * t536 + t595 * t537; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t560, t561, 0, t539; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t600, -qJDD(1), 0, t582, 0, 0, 0, 0, 0, 0, -t568, t570, 0, t543, 0, 0, 0, 0, 0, 0, 0, t568, -t570, t541, 0, 0, 0, 0, 0, 0, t549, t550, t555, t537; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t600, 0, t581, 0, 0, 0, 0, 0, 0, -t570, -t568, 0, t542, 0, 0, 0, 0, 0, 0, 0, t570, t568, t540, 0, 0, 0, 0, 0, 0, t547, t548, t554, t536; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t560, t561, 0, t539; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t588, -t589, 0, t557, 0, 0, 0, 0, 0, 0, 0, t588, t589, t552, 0, 0, 0, 0, 0, 0, t564, t565, -t572, t546; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t589, -t588, 0, t556, 0, 0, 0, 0, 0, 0, 0, -t589, t588, -t553, 0, 0, 0, 0, 0, 0, -t558, -t559, t567, -t538; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t560, t561, 0, t539; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t560, t561, 0, t539; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t588, -t589, -t552, 0, 0, 0, 0, 0, 0, -t564, -t565, t572, -t546; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t589, -t588, t553, 0, 0, 0, 0, 0, 0, t558, t559, -t567, t538; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t579, t576, -t609, t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t575, t580, -t608, t544; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t564, t565, -t572, t546;];
f_new_reg = t1;
