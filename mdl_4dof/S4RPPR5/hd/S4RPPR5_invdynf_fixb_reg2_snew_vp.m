% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPPR5
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPPR5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR5_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:54
% EndTime: 2019-12-31 16:39:55
% DurationCPUTime: 0.60s
% Computational Cost: add. (1025->114), mult. (1830->112), div. (0->0), fcn. (826->6), ass. (0->62)
t490 = -pkin(1) - pkin(2);
t480 = qJD(1) ^ 2;
t476 = sin(qJ(1));
t478 = cos(qJ(1));
t461 = -g(1) * t478 - t476 * g(2);
t483 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t461;
t446 = t490 * t480 + t483;
t472 = sin(pkin(6));
t473 = cos(pkin(6));
t460 = t476 * g(1) - g(2) * t478;
t482 = -t480 * qJ(2) + qJDD(2) - t460;
t481 = t490 * qJDD(1) + t482;
t432 = t473 * t446 + t472 * t481;
t475 = sin(qJ(4));
t469 = t475 ^ 2;
t477 = cos(qJ(4));
t470 = t477 ^ 2;
t489 = t469 + t470;
t488 = qJD(1) * qJD(4);
t487 = t475 * qJDD(1);
t486 = t477 * qJDD(1);
t431 = -t446 * t472 + t473 * t481;
t452 = -t472 * qJDD(1) + t480 * t473;
t453 = qJDD(1) * t473 + t480 * t472;
t485 = -t452 * t476 + t478 * t453;
t484 = t452 * t478 + t476 * t453;
t479 = qJD(4) ^ 2;
t471 = g(3) + qJDD(3);
t464 = t477 * t480 * t475;
t463 = -t480 * t470 - t479;
t462 = -t480 * t469 - t479;
t459 = -qJDD(4) + t464;
t458 = qJDD(4) + t464;
t457 = t489 * t480;
t456 = qJDD(1) * t478 - t476 * t480;
t455 = t476 * qJDD(1) + t480 * t478;
t454 = t489 * qJDD(1);
t451 = -0.2e1 * t475 * t488 + t486;
t450 = 0.2e1 * t477 * t488 + t487;
t448 = qJDD(1) * pkin(1) - t482;
t447 = -t480 * pkin(1) + t483;
t443 = t459 * t477 - t462 * t475;
t442 = -t458 * t475 + t463 * t477;
t441 = t459 * t475 + t462 * t477;
t440 = t458 * t477 + t463 * t475;
t438 = -t454 * t473 - t457 * t472;
t437 = -t454 * t472 + t457 * t473;
t436 = t443 * t473 - t450 * t472;
t435 = t442 * t473 + t451 * t472;
t434 = t443 * t472 + t450 * t473;
t433 = t442 * t472 - t451 * t473;
t430 = -t480 * pkin(3) - qJDD(1) * pkin(5) + t432;
t429 = qJDD(1) * pkin(3) - t480 * pkin(5) - t431;
t428 = t430 * t477 + t471 * t475;
t427 = -t430 * t475 + t471 * t477;
t426 = -t431 * t472 + t432 * t473;
t425 = t431 * t473 + t432 * t472;
t424 = -t427 * t475 + t428 * t477;
t423 = t427 * t477 + t428 * t475;
t422 = t424 * t473 + t429 * t472;
t421 = t424 * t472 - t429 * t473;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t455, -t456, 0, -t476 * t460 + t461 * t478, 0, 0, 0, 0, 0, 0, -t455, 0, t456, t447 * t478 - t476 * t448, 0, 0, 0, 0, 0, 0, -t484, t485, 0, t476 * t425 + t426 * t478, 0, 0, 0, 0, 0, 0, t476 * t433 + t435 * t478, t476 * t434 + t436 * t478, t476 * t437 + t438 * t478, t476 * t421 + t422 * t478; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t456, -t455, 0, t460 * t478 + t476 * t461, 0, 0, 0, 0, 0, 0, t456, 0, t455, t476 * t447 + t448 * t478, 0, 0, 0, 0, 0, 0, t485, t484, 0, -t425 * t478 + t476 * t426, 0, 0, 0, 0, 0, 0, -t433 * t478 + t476 * t435, -t434 * t478 + t476 * t436, -t437 * t478 + t476 * t438, -t421 * t478 + t476 * t422; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t471, 0, 0, 0, 0, 0, 0, -t440, -t441, 0, -t423; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t480, -qJDD(1), 0, t461, 0, 0, 0, 0, 0, 0, -t480, 0, qJDD(1), t447, 0, 0, 0, 0, 0, 0, -t452, t453, 0, t426, 0, 0, 0, 0, 0, 0, t435, t436, t438, t422; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t480, 0, t460, 0, 0, 0, 0, 0, 0, qJDD(1), 0, t480, t448, 0, 0, 0, 0, 0, 0, t453, t452, 0, -t425, 0, 0, 0, 0, 0, 0, -t433, -t434, -t437, -t421; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t471, 0, 0, 0, 0, 0, 0, -t440, -t441, 0, -t423; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t480, 0, qJDD(1), t447, 0, 0, 0, 0, 0, 0, -t452, t453, 0, t426, 0, 0, 0, 0, 0, 0, t435, t436, t438, t422; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t471, 0, 0, 0, 0, 0, 0, -t440, -t441, 0, -t423; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t480, -t448, 0, 0, 0, 0, 0, 0, -t453, -t452, 0, t425, 0, 0, 0, 0, 0, 0, t433, t434, t437, t421; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t480, qJDD(1), 0, t432, 0, 0, 0, 0, 0, 0, t442, t443, -t454, t424; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), -t480, 0, t431, 0, 0, 0, 0, 0, 0, -t451, t450, t457, -t429; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t471, 0, 0, 0, 0, 0, 0, t440, t441, 0, t423; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t463, t459, -t486, t428; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t458, t462, t487, t427; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t451, -t450, -t457, t429;];
f_new_reg = t1;
