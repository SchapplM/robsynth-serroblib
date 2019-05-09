% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRRP1
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRRP1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP1_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:24:29
% EndTime: 2019-05-04 19:24:29
% DurationCPUTime: 0.69s
% Computational Cost: add. (1568->83), mult. (2018->81), div. (0->0), fcn. (1288->6), ass. (0->49)
t524 = qJD(1) + qJD(2);
t520 = qJD(3) + t524;
t518 = t520 ^ 2;
t523 = qJDD(1) + qJDD(2);
t519 = qJDD(3) + t523;
t526 = sin(qJ(3));
t529 = cos(qJ(3));
t503 = t526 * t518 - t529 * t519;
t527 = sin(qJ(2));
t530 = cos(qJ(2));
t534 = -t529 * t518 - t526 * t519;
t492 = t530 * t503 - t527 * t534;
t528 = sin(qJ(1));
t531 = cos(qJ(1));
t541 = t527 * t503 + t530 * t534;
t544 = t528 * t492 + t531 * t541;
t482 = t531 * t492 - t528 * t541;
t522 = t524 ^ 2;
t510 = t527 * t522 - t530 * t523;
t533 = -t530 * t522 - t527 * t523;
t540 = t528 * t510 + t531 * t533;
t539 = t531 * t510 - t528 * t533;
t516 = t528 * g(1) - t531 * g(2);
t512 = qJDD(1) * pkin(1) + t516;
t517 = -t531 * g(1) - t528 * g(2);
t532 = qJD(1) ^ 2;
t513 = -t532 * pkin(1) + t517;
t499 = t530 * t512 - t527 * t513;
t497 = t523 * pkin(2) + t499;
t500 = t527 * t512 + t530 * t513;
t498 = -t522 * pkin(2) + t500;
t487 = t526 * t497 + t529 * t498;
t486 = t529 * t497 - t526 * t498;
t525 = -g(3) + qJDD(4);
t515 = -t528 * qJDD(1) - t531 * t532;
t514 = t531 * qJDD(1) - t528 * t532;
t489 = -t527 * t499 + t530 * t500;
t488 = t530 * t499 + t527 * t500;
t485 = -t518 * pkin(3) + t487;
t484 = t519 * pkin(3) + t486;
t479 = -t526 * t486 + t529 * t487;
t478 = t529 * t486 + t526 * t487;
t477 = -t526 * t484 + t529 * t485;
t476 = t529 * t484 + t526 * t485;
t475 = -t527 * t478 + t530 * t479;
t474 = t530 * t478 + t527 * t479;
t473 = -t527 * t476 + t530 * t477;
t472 = t530 * t476 + t527 * t477;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t515, -t514, 0, -t528 * t516 + t531 * t517, 0, 0, 0, 0, 0, 0, t540, t539, 0, -t528 * t488 + t531 * t489, 0, 0, 0, 0, 0, 0, t544, t482, 0, -t528 * t474 + t531 * t475, 0, 0, 0, 0, 0, 0, t544, t482, 0, -t528 * t472 + t531 * t473; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t514, t515, 0, t531 * t516 + t528 * t517, 0, 0, 0, 0, 0, 0, -t539, t540, 0, t531 * t488 + t528 * t489, 0, 0, 0, 0, 0, 0, -t482, t544, 0, t531 * t474 + t528 * t475, 0, 0, 0, 0, 0, 0, -t482, t544, 0, t531 * t472 + t528 * t473; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t525; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t532, -qJDD(1), 0, t517, 0, 0, 0, 0, 0, 0, t533, t510, 0, t489, 0, 0, 0, 0, 0, 0, t541, t492, 0, t475, 0, 0, 0, 0, 0, 0, t541, t492, 0, t473; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t532, 0, t516, 0, 0, 0, 0, 0, 0, -t510, t533, 0, t488, 0, 0, 0, 0, 0, 0, -t492, t541, 0, t474, 0, 0, 0, 0, 0, 0, -t492, t541, 0, t472; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t525; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t522, -t523, 0, t500, 0, 0, 0, 0, 0, 0, t534, t503, 0, t479, 0, 0, 0, 0, 0, 0, t534, t503, 0, t477; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t523, -t522, 0, t499, 0, 0, 0, 0, 0, 0, -t503, t534, 0, t478, 0, 0, 0, 0, 0, 0, -t503, t534, 0, t476; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t525; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t518, -t519, 0, t487, 0, 0, 0, 0, 0, 0, -t518, -t519, 0, t485; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t519, -t518, 0, t486, 0, 0, 0, 0, 0, 0, t519, -t518, 0, t484; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t525; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t518, -t519, 0, t485; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t519, -t518, 0, t484; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t525;];
f_new_reg  = t1;
