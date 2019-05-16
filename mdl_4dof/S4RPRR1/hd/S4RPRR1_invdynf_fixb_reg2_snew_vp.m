% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPRR1
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPRR1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR1_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_invdynf_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:16:52
% EndTime: 2019-05-04 19:16:53
% DurationCPUTime: 0.90s
% Computational Cost: add. (2975->86), mult. (4606->105), div. (0->0), fcn. (2880->8), ass. (0->63)
t597 = qJD(1) + qJD(3);
t593 = qJD(4) + t597;
t590 = t593 ^ 2;
t596 = qJDD(1) + qJDD(3);
t591 = qJDD(4) + t596;
t601 = sin(qJ(4));
t604 = cos(qJ(4));
t572 = t601 * t590 - t604 * t591;
t602 = sin(qJ(3));
t605 = cos(qJ(3));
t609 = -t604 * t590 - t601 * t591;
t557 = t605 * t572 - t602 * t609;
t599 = sin(pkin(7));
t600 = cos(pkin(7));
t617 = t602 * t572 + t605 * t609;
t547 = t600 * t557 - t599 * t617;
t603 = sin(qJ(1));
t606 = cos(qJ(1));
t624 = t599 * t557 + t600 * t617;
t628 = t603 * t547 + t606 * t624;
t627 = t606 * t547 - t603 * t624;
t595 = t597 ^ 2;
t579 = t602 * t595 - t605 * t596;
t608 = -t605 * t595 - t602 * t596;
t564 = t600 * t579 - t599 * t608;
t616 = t599 * t579 + t600 * t608;
t623 = t603 * t564 + t606 * t616;
t622 = t606 * t564 - t603 * t616;
t588 = t603 * g(1) - t606 * g(2);
t582 = qJDD(1) * pkin(1) + t588;
t589 = -t606 * g(1) - t603 * g(2);
t607 = qJD(1) ^ 2;
t583 = -t607 * pkin(1) + t589;
t568 = t600 * t582 - t599 * t583;
t566 = qJDD(1) * pkin(2) + t568;
t569 = t599 * t582 + t600 * t583;
t567 = -t607 * pkin(2) + t569;
t552 = t602 * t566 + t605 * t567;
t551 = t605 * t566 - t602 * t567;
t584 = -t599 * qJDD(1) - t600 * t607;
t585 = t600 * qJDD(1) - t599 * t607;
t611 = t606 * t584 - t603 * t585;
t610 = t603 * t584 + t606 * t585;
t598 = -g(3) + qJDD(2);
t587 = -t603 * qJDD(1) - t606 * t607;
t586 = t606 * qJDD(1) - t603 * t607;
t554 = -t599 * t568 + t600 * t569;
t553 = t600 * t568 + t599 * t569;
t550 = -t595 * pkin(3) + t552;
t549 = t596 * pkin(3) + t551;
t544 = -t602 * t551 + t605 * t552;
t543 = t605 * t551 + t602 * t552;
t542 = t601 * t549 + t604 * t550;
t541 = t604 * t549 - t601 * t550;
t540 = -t599 * t543 + t600 * t544;
t539 = t600 * t543 + t599 * t544;
t538 = -t601 * t541 + t604 * t542;
t537 = t604 * t541 + t601 * t542;
t536 = -t602 * t537 + t605 * t538;
t535 = t605 * t537 + t602 * t538;
t534 = -t599 * t535 + t600 * t536;
t533 = t600 * t535 + t599 * t536;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t587, -t586, 0, -t603 * t588 + t606 * t589, 0, 0, 0, 0, 0, 0, t611, -t610, 0, -t603 * t553 + t606 * t554, 0, 0, 0, 0, 0, 0, t623, t622, 0, -t603 * t539 + t606 * t540, 0, 0, 0, 0, 0, 0, t628, t627, 0, -t603 * t533 + t606 * t534; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t586, t587, 0, t606 * t588 + t603 * t589, 0, 0, 0, 0, 0, 0, t610, t611, 0, t606 * t553 + t603 * t554, 0, 0, 0, 0, 0, 0, -t622, t623, 0, t606 * t539 + t603 * t540, 0, 0, 0, 0, 0, 0, -t627, t628, 0, t606 * t533 + t603 * t534; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t598, 0, 0, 0, 0, 0, 0, 0, 0, 0, t598, 0, 0, 0, 0, 0, 0, 0, 0, 0, t598; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t607, -qJDD(1), 0, t589, 0, 0, 0, 0, 0, 0, t584, -t585, 0, t554, 0, 0, 0, 0, 0, 0, t616, t564, 0, t540, 0, 0, 0, 0, 0, 0, t624, t547, 0, t534; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t607, 0, t588, 0, 0, 0, 0, 0, 0, t585, t584, 0, t553, 0, 0, 0, 0, 0, 0, -t564, t616, 0, t539, 0, 0, 0, 0, 0, 0, -t547, t624, 0, t533; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t598, 0, 0, 0, 0, 0, 0, 0, 0, 0, t598, 0, 0, 0, 0, 0, 0, 0, 0, 0, t598; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t607, -qJDD(1), 0, t569, 0, 0, 0, 0, 0, 0, t608, t579, 0, t544, 0, 0, 0, 0, 0, 0, t617, t557, 0, t536; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t607, 0, t568, 0, 0, 0, 0, 0, 0, -t579, t608, 0, t543, 0, 0, 0, 0, 0, 0, -t557, t617, 0, t535; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t598, 0, 0, 0, 0, 0, 0, 0, 0, 0, t598, 0, 0, 0, 0, 0, 0, 0, 0, 0, t598; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t595, -t596, 0, t552, 0, 0, 0, 0, 0, 0, t609, t572, 0, t538; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t596, -t595, 0, t551, 0, 0, 0, 0, 0, 0, -t572, t609, 0, t537; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t598, 0, 0, 0, 0, 0, 0, 0, 0, 0, t598; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t590, -t591, 0, t542; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t591, -t590, 0, t541; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t598;];
f_new_reg  = t1;
