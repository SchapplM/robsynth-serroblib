% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPRP6
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPRP6_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP6_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:19
% EndTime: 2019-12-31 16:46:20
% DurationCPUTime: 0.59s
% Computational Cost: add. (641->108), mult. (1363->97), div. (0->0), fcn. (645->4), ass. (0->60)
t589 = -2 * qJD(1);
t588 = -pkin(5) - pkin(1);
t569 = sin(qJ(3));
t567 = t569 ^ 2;
t574 = qJD(1) ^ 2;
t587 = t567 * t574;
t586 = t569 * t574;
t570 = sin(qJ(1));
t572 = cos(qJ(1));
t556 = t570 * g(1) - t572 * g(2);
t576 = -(t574 * qJ(2)) + qJDD(2) - t556;
t543 = t588 * qJDD(1) + t576;
t571 = cos(qJ(3));
t532 = t569 * g(3) + t571 * t543;
t568 = t571 ^ 2;
t585 = t567 + t568;
t584 = qJD(1) * t571;
t583 = qJD(1) * qJD(3);
t582 = t569 * qJDD(1);
t581 = t571 * qJDD(1);
t580 = qJD(4) * t589;
t579 = t571 * t586;
t578 = t571 * t583;
t533 = -t571 * g(3) + t569 * t543;
t557 = -t572 * g(1) - t570 * g(2);
t577 = -(qJD(3) * pkin(3)) + qJ(4) * t584;
t575 = -qJDD(1) * qJ(2) + (qJD(2) * t589) - t557;
t573 = qJD(3) ^ 2;
t559 = -t568 * t574 - t573;
t558 = -t573 - t587;
t555 = -qJDD(3) - t579;
t554 = qJDD(3) - t579;
t553 = t585 * t574;
t552 = t570 * qJDD(1) + t572 * t574;
t551 = t572 * qJDD(1) - t570 * t574;
t550 = t585 * qJDD(1);
t549 = -0.2e1 * t569 * t583 + t581;
t548 = -t578 - t582;
t547 = 0.2e1 * t578 + t582;
t545 = qJDD(1) * pkin(1) - t576;
t544 = t574 * pkin(1) + t575;
t542 = -t588 * t574 + t575;
t539 = t571 * t555 - t569 * t559;
t538 = -t569 * t554 + t571 * t558;
t537 = t569 * t555 + t571 * t559;
t536 = t571 * t554 + t569 * t558;
t535 = -t570 * t550 - t572 * t553;
t534 = t572 * t550 - t570 * t553;
t531 = t570 * t537 + t572 * t549;
t530 = t570 * t536 + t572 * t547;
t529 = -t572 * t537 + t570 * t549;
t528 = -t572 * t536 + t570 * t547;
t527 = t548 * pkin(3) - qJDD(4) + t577 * t584 + (qJ(4) * t567 - t588) * t574 + t575;
t526 = -t569 * t532 + t571 * t533;
t525 = t571 * t532 + t569 * t533;
t524 = -pkin(3) * t587 + t548 * qJ(4) + qJD(3) * t577 + t569 * t580 + t533;
t523 = qJDD(3) * pkin(3) + (-pkin(3) * t586 - qJDD(1) * qJ(4) + t580) * t571 + t532;
t522 = -t569 * t523 + t571 * t524;
t521 = t571 * t523 + t569 * t524;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t552, -t551, 0, -t570 * t556 + t572 * t557, 0, 0, 0, 0, 0, 0, 0, t552, t551, -t572 * t544 - t570 * t545, 0, 0, 0, 0, 0, 0, t530, t531, t535, t570 * t525 - t572 * t542, 0, 0, 0, 0, 0, 0, t530, t531, t535, t570 * t521 - t572 * t527; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t551, -t552, 0, t572 * t556 + t570 * t557, 0, 0, 0, 0, 0, 0, 0, -t551, t552, -t570 * t544 + t572 * t545, 0, 0, 0, 0, 0, 0, t528, t529, t534, -t572 * t525 - t570 * t542, 0, 0, 0, 0, 0, 0, t528, t529, t534, -t572 * t521 - t570 * t527; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t538, t539, 0, t526, 0, 0, 0, 0, 0, 0, t538, t539, 0, t522; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t574, -qJDD(1), 0, t557, 0, 0, 0, 0, 0, 0, 0, t574, qJDD(1), -t544, 0, 0, 0, 0, 0, 0, t547, t549, -t553, -t542, 0, 0, 0, 0, 0, 0, t547, t549, -t553, -t527; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t574, 0, t556, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), t574, t545, 0, 0, 0, 0, 0, 0, -t536, -t537, t550, -t525, 0, 0, 0, 0, 0, 0, -t536, -t537, t550, -t521; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t538, t539, 0, t526, 0, 0, 0, 0, 0, 0, t538, t539, 0, t522; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t538, t539, 0, t526, 0, 0, 0, 0, 0, 0, t538, t539, 0, t522; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t574, -qJDD(1), t544, 0, 0, 0, 0, 0, 0, -t547, -t549, t553, t542, 0, 0, 0, 0, 0, 0, -t547, -t549, t553, t527; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t574, -t545, 0, 0, 0, 0, 0, 0, t536, t537, -t550, t525, 0, 0, 0, 0, 0, 0, t536, t537, -t550, t521; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t558, t555, -t582, t533, 0, 0, 0, 0, 0, 0, t558, t555, -t582, t524; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t554, t559, -t581, t532, 0, 0, 0, 0, 0, 0, t554, t559, -t581, t523; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t547, t549, -t553, -t542, 0, 0, 0, 0, 0, 0, t547, t549, -t553, -t527; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t558, t555, -t582, t524; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t554, t559, -t581, t523; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t547, t549, -t553, -t527;];
f_new_reg = t1;
