% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRPR7
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRPR7_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR7_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:00
% EndTime: 2019-12-31 16:26:01
% DurationCPUTime: 0.52s
% Computational Cost: add. (620->93), mult. (1146->102), div. (0->0), fcn. (722->6), ass. (0->61)
t507 = sin(qJ(2));
t509 = cos(qJ(2));
t511 = qJD(2) ^ 2;
t487 = t507 * qJDD(2) + t509 * t511;
t503 = sin(pkin(6));
t524 = t503 * t487;
t488 = -t509 * qJDD(2) + t507 * t511;
t523 = t503 * t488;
t504 = cos(pkin(6));
t490 = t503 * g(1) - t504 * g(2);
t522 = t503 * t490;
t521 = t504 * t487;
t520 = t504 * t488;
t506 = sin(qJ(4));
t499 = t506 ^ 2;
t508 = cos(qJ(4));
t500 = t508 ^ 2;
t519 = t499 + t500;
t518 = qJD(2) * qJD(4);
t517 = t506 * qJDD(2);
t516 = t508 * qJDD(2);
t515 = t508 * t511 * t506;
t491 = -t504 * g(1) - t503 * g(2);
t501 = -g(3) + qJDD(1);
t479 = -t507 * t491 + t509 * t501;
t510 = qJD(4) ^ 2;
t514 = -t500 * t511 - t510;
t513 = qJDD(4) - t515;
t472 = -qJDD(2) * pkin(2) - t511 * qJ(3) + qJDD(3) - t479;
t470 = -qJDD(2) * pkin(5) + t472;
t512 = t508 * t470 + t506 * t490;
t480 = t509 * t491 + t507 * t501;
t471 = -t511 * pkin(2) + qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t480;
t493 = -t499 * t511 - t510;
t492 = -qJDD(4) - t515;
t489 = t519 * t511;
t486 = t519 * qJDD(2);
t484 = -0.2e1 * t506 * t518 + t516;
t483 = 0.2e1 * t508 * t518 + t517;
t481 = t504 * t490;
t478 = t508 * t492 - t506 * t514;
t477 = t508 * t493 - t506 * t513;
t476 = t506 * t492 + t508 * t514;
t475 = t506 * t493 + t508 * t513;
t474 = -t507 * t486 - t509 * t489;
t473 = t509 * t486 - t507 * t489;
t469 = -t511 * pkin(5) + t471;
t468 = t507 * t476 + t509 * t484;
t467 = t507 * t475 + t509 * t483;
t466 = -t509 * t476 + t507 * t484;
t465 = -t509 * t475 + t507 * t483;
t464 = -t507 * t479 + t509 * t480;
t463 = t509 * t479 + t507 * t480;
t462 = t506 * t470 - t508 * t490;
t460 = t509 * t471 + t507 * t472;
t459 = t507 * t471 - t509 * t472;
t458 = t508 * t462 - t506 * t512;
t457 = t506 * t462 + t508 * t512;
t456 = t507 * t457 + t509 * t469;
t455 = -t509 * t457 + t507 * t469;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t504 * t491 - t522, 0, 0, 0, 0, 0, 0, -t521, t520, 0, t504 * t464 - t522, 0, 0, 0, 0, 0, 0, 0, t521, -t520, t504 * t460 - t522, 0, 0, 0, 0, 0, 0, t504 * t467 + t503 * t477, t504 * t468 + t503 * t478, t504 * t474, t504 * t456 + t503 * t458; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t503 * t491 + t481, 0, 0, 0, 0, 0, 0, -t524, t523, 0, t503 * t464 + t481, 0, 0, 0, 0, 0, 0, 0, t524, -t523, t503 * t460 + t481, 0, 0, 0, 0, 0, 0, t503 * t467 - t504 * t477, t503 * t468 - t504 * t478, t503 * t474, t503 * t456 - t504 * t458; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t501, 0, 0, 0, 0, 0, 0, -t488, -t487, 0, t463, 0, 0, 0, 0, 0, 0, 0, t488, t487, t459, 0, 0, 0, 0, 0, 0, t465, t466, t473, t455; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t491, 0, 0, 0, 0, 0, 0, -t487, t488, 0, t464, 0, 0, 0, 0, 0, 0, 0, t487, -t488, t460, 0, 0, 0, 0, 0, 0, t467, t468, t474, t456; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t490, 0, 0, 0, 0, 0, 0, 0, 0, 0, t490, 0, 0, 0, 0, 0, 0, 0, 0, 0, t490, 0, 0, 0, 0, 0, 0, -t477, -t478, 0, -t458; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t501, 0, 0, 0, 0, 0, 0, -t488, -t487, 0, t463, 0, 0, 0, 0, 0, 0, 0, t488, t487, t459, 0, 0, 0, 0, 0, 0, t465, t466, t473, t455; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t511, -qJDD(2), 0, t480, 0, 0, 0, 0, 0, 0, 0, t511, qJDD(2), t471, 0, 0, 0, 0, 0, 0, t483, t484, -t489, t469; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t511, 0, t479, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), t511, -t472, 0, 0, 0, 0, 0, 0, -t475, -t476, t486, -t457; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t490, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t490, 0, 0, 0, 0, 0, 0, t477, t478, 0, t458; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t490, 0, 0, 0, 0, 0, 0, t477, t478, 0, t458; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t511, -qJDD(2), -t471, 0, 0, 0, 0, 0, 0, -t483, -t484, t489, -t469; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t511, t472, 0, 0, 0, 0, 0, 0, t475, t476, -t486, t457; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t493, t492, -t517, t462; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t513, t514, -t516, t512; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t483, t484, -t489, t469;];
f_new_reg = t1;
