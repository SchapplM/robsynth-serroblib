% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPPR4
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPPR4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR4_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:00
% EndTime: 2019-12-31 16:39:01
% DurationCPUTime: 0.67s
% Computational Cost: add. (804->99), mult. (1502->114), div. (0->0), fcn. (848->6), ass. (0->61)
t571 = sin(pkin(6));
t572 = cos(pkin(6));
t579 = qJD(1) ^ 2;
t551 = t571 * qJDD(1) + t572 * t579;
t552 = -t572 * qJDD(1) + t571 * t579;
t575 = sin(qJ(1));
t577 = cos(qJ(1));
t581 = t575 * t551 + t577 * t552;
t583 = t577 * t551 - t575 * t552;
t574 = sin(qJ(4));
t567 = t574 ^ 2;
t576 = cos(qJ(4));
t568 = t576 ^ 2;
t588 = t567 + t568;
t587 = qJD(1) * qJD(4);
t586 = t574 * qJDD(1);
t585 = t576 * qJDD(1);
t584 = t574 * t579 * t576;
t559 = t575 * g(1) - t577 * g(2);
t547 = qJDD(1) * pkin(1) + t559;
t560 = -t577 * g(1) - t575 * g(2);
t548 = -t579 * pkin(1) + t560;
t534 = t572 * t547 - t571 * t548;
t535 = t571 * t547 + t572 * t548;
t533 = -qJDD(1) * pkin(2) - t579 * qJ(3) + qJDD(3) - t534;
t532 = -t579 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t535;
t578 = qJD(4) ^ 2;
t569 = -g(3) + qJDD(2);
t562 = -t568 * t579 - t578;
t561 = -t567 * t579 - t578;
t558 = -qJDD(4) - t584;
t557 = qJDD(4) - t584;
t556 = t588 * t579;
t555 = -t575 * qJDD(1) - t577 * t579;
t554 = t577 * qJDD(1) - t575 * t579;
t553 = t588 * qJDD(1);
t550 = -0.2e1 * t574 * t587 + t585;
t549 = 0.2e1 * t576 * t587 + t586;
t541 = t576 * t558 - t574 * t562;
t540 = -t574 * t557 + t576 * t561;
t539 = t574 * t558 + t576 * t562;
t538 = t576 * t557 + t574 * t561;
t537 = -t571 * t553 - t572 * t556;
t536 = t572 * t553 - t571 * t556;
t531 = -qJDD(1) * pkin(5) + t533;
t530 = t571 * t539 + t572 * t550;
t529 = t571 * t538 + t572 * t549;
t528 = -t572 * t539 + t571 * t550;
t527 = -t572 * t538 + t571 * t549;
t526 = -t579 * pkin(5) + t532;
t525 = t574 * t531 + t576 * t569;
t524 = t576 * t531 - t574 * t569;
t523 = -t571 * t534 + t572 * t535;
t522 = t572 * t534 + t571 * t535;
t521 = t572 * t532 + t571 * t533;
t520 = t571 * t532 - t572 * t533;
t519 = -t574 * t524 + t576 * t525;
t518 = t576 * t524 + t574 * t525;
t517 = t571 * t518 + t572 * t526;
t516 = -t572 * t518 + t571 * t526;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t555, -t554, 0, -t575 * t559 + t577 * t560, 0, 0, 0, 0, 0, 0, -t583, t581, 0, -t575 * t522 + t577 * t523, 0, 0, 0, 0, 0, 0, 0, t583, -t581, -t575 * t520 + t577 * t521, 0, 0, 0, 0, 0, 0, -t575 * t527 + t577 * t529, -t575 * t528 + t577 * t530, -t575 * t536 + t577 * t537, -t575 * t516 + t577 * t517; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t554, t555, 0, t577 * t559 + t575 * t560, 0, 0, 0, 0, 0, 0, -t581, -t583, 0, t577 * t522 + t575 * t523, 0, 0, 0, 0, 0, 0, 0, t581, t583, t577 * t520 + t575 * t521, 0, 0, 0, 0, 0, 0, t577 * t527 + t575 * t529, t577 * t528 + t575 * t530, t577 * t536 + t575 * t537, t577 * t516 + t575 * t517; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t569, 0, 0, 0, 0, 0, 0, 0, 0, 0, t569, 0, 0, 0, 0, 0, 0, t540, t541, 0, t519; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t579, -qJDD(1), 0, t560, 0, 0, 0, 0, 0, 0, -t551, t552, 0, t523, 0, 0, 0, 0, 0, 0, 0, t551, -t552, t521, 0, 0, 0, 0, 0, 0, t529, t530, t537, t517; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t579, 0, t559, 0, 0, 0, 0, 0, 0, -t552, -t551, 0, t522, 0, 0, 0, 0, 0, 0, 0, t552, t551, t520, 0, 0, 0, 0, 0, 0, t527, t528, t536, t516; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t569, 0, 0, 0, 0, 0, 0, 0, 0, 0, t569, 0, 0, 0, 0, 0, 0, t540, t541, 0, t519; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t579, -qJDD(1), 0, t535, 0, 0, 0, 0, 0, 0, 0, t579, qJDD(1), t532, 0, 0, 0, 0, 0, 0, t549, t550, -t556, t526; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t579, 0, t534, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), t579, -t533, 0, 0, 0, 0, 0, 0, -t538, -t539, t553, -t518; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t569, 0, 0, 0, 0, 0, 0, 0, 0, 0, t569, 0, 0, 0, 0, 0, 0, t540, t541, 0, t519; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t569, 0, 0, 0, 0, 0, 0, t540, t541, 0, t519; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t579, -qJDD(1), -t532, 0, 0, 0, 0, 0, 0, -t549, -t550, t556, -t526; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t579, t533, 0, 0, 0, 0, 0, 0, t538, t539, -t553, t518; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t561, t558, -t586, t525; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t557, t562, -t585, t524; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t549, t550, -t556, t526;];
f_new_reg = t1;
