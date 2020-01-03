% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPRR5
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPRR5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR5_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:46
% EndTime: 2019-12-31 16:51:46
% DurationCPUTime: 0.64s
% Computational Cost: add. (1243->113), mult. (1830->113), div. (0->0), fcn. (826->6), ass. (0->64)
t514 = sin(qJ(1));
t517 = cos(qJ(1));
t508 = -qJD(1) + qJD(3);
t506 = t508 ^ 2;
t507 = qJDD(1) - qJDD(3);
t513 = sin(qJ(3));
t516 = cos(qJ(3));
t523 = t513 * t506 + t516 * t507;
t524 = -t516 * t506 + t513 * t507;
t531 = t514 * t524 + t517 * t523;
t530 = t514 * t523 - t517 * t524;
t529 = -pkin(1) - pkin(2);
t512 = sin(qJ(4));
t528 = t512 * t507;
t515 = cos(qJ(4));
t527 = t515 * t507;
t519 = qJD(1) ^ 2;
t500 = -t517 * g(1) - t514 * g(2);
t522 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t500;
t481 = t519 * t529 + t522;
t499 = t514 * g(1) - t517 * g(2);
t521 = -t519 * qJ(2) + qJDD(2) - t499;
t520 = qJDD(1) * t529 + t521;
t471 = t516 * t481 + t513 * t520;
t510 = t512 ^ 2;
t511 = t515 ^ 2;
t526 = t510 + t511;
t525 = qJD(4) * t508;
t470 = -t513 * t481 + t516 * t520;
t518 = qJD(4) ^ 2;
t498 = t515 * t506 * t512;
t497 = -t511 * t506 - t518;
t496 = -t510 * t506 - t518;
t495 = t517 * qJDD(1) - t514 * t519;
t494 = t514 * qJDD(1) + t517 * t519;
t493 = -qJDD(4) + t498;
t492 = qJDD(4) + t498;
t491 = t526 * t506;
t486 = t526 * t507;
t485 = -0.2e1 * t512 * t525 - t527;
t484 = 0.2e1 * t515 * t525 - t528;
t483 = qJDD(1) * pkin(1) - t521;
t482 = -t519 * pkin(1) + t522;
t477 = t515 * t493 - t512 * t496;
t476 = -t512 * t492 + t515 * t497;
t475 = t512 * t493 + t515 * t496;
t474 = t515 * t492 + t512 * t497;
t473 = -t516 * t486 - t513 * t491;
t472 = -t513 * t486 + t516 * t491;
t469 = t516 * t477 + t513 * t484;
t468 = t516 * t476 - t513 * t485;
t467 = t513 * t477 - t516 * t484;
t466 = t513 * t476 + t516 * t485;
t465 = -t506 * pkin(3) - t507 * pkin(6) + t471;
t464 = t507 * pkin(3) - t506 * pkin(6) - t470;
t463 = t512 * g(3) + t515 * t465;
t462 = t515 * g(3) - t512 * t465;
t461 = -t513 * t470 + t516 * t471;
t460 = t516 * t470 + t513 * t471;
t459 = -t512 * t462 + t515 * t463;
t458 = t515 * t462 + t512 * t463;
t457 = t516 * t459 + t513 * t464;
t456 = t513 * t459 - t516 * t464;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t494, -t495, 0, -t514 * t499 + t517 * t500, 0, 0, 0, 0, 0, 0, -t494, 0, t495, t517 * t482 - t514 * t483, 0, 0, 0, 0, 0, 0, -t530, t531, 0, t514 * t460 + t517 * t461, 0, 0, 0, 0, 0, 0, t514 * t466 + t517 * t468, t514 * t467 + t517 * t469, t514 * t472 + t517 * t473, t514 * t456 + t517 * t457; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t495, -t494, 0, t517 * t499 + t514 * t500, 0, 0, 0, 0, 0, 0, t495, 0, t494, t514 * t482 + t517 * t483, 0, 0, 0, 0, 0, 0, t531, t530, 0, -t517 * t460 + t514 * t461, 0, 0, 0, 0, 0, 0, -t517 * t466 + t514 * t468, -t517 * t467 + t514 * t469, -t517 * t472 + t514 * t473, -t517 * t456 + t514 * t457; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t474, -t475, 0, -t458; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t519, -qJDD(1), 0, t500, 0, 0, 0, 0, 0, 0, -t519, 0, qJDD(1), t482, 0, 0, 0, 0, 0, 0, t524, t523, 0, t461, 0, 0, 0, 0, 0, 0, t468, t469, t473, t457; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t519, 0, t499, 0, 0, 0, 0, 0, 0, qJDD(1), 0, t519, t483, 0, 0, 0, 0, 0, 0, t523, -t524, 0, -t460, 0, 0, 0, 0, 0, 0, -t466, -t467, -t472, -t456; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t474, -t475, 0, -t458; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t519, 0, qJDD(1), t482, 0, 0, 0, 0, 0, 0, t524, t523, 0, t461, 0, 0, 0, 0, 0, 0, t468, t469, t473, t457; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t474, -t475, 0, -t458; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t519, -t483, 0, 0, 0, 0, 0, 0, -t523, t524, 0, t460, 0, 0, 0, 0, 0, 0, t466, t467, t472, t456; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t506, t507, 0, t471, 0, 0, 0, 0, 0, 0, t476, t477, -t486, t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t507, -t506, 0, t470, 0, 0, 0, 0, 0, 0, t485, -t484, t491, -t464; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, t474, t475, 0, t458; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t497, t493, -t527, t463; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t492, t496, t528, t462; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t485, t484, -t491, t464;];
f_new_reg = t1;
