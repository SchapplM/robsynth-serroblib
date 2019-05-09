% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRPR1
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
% Datum: 2019-05-04 18:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRPR1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR1_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:57:17
% EndTime: 2019-05-04 18:57:17
% DurationCPUTime: 0.53s
% Computational Cost: add. (870->84), mult. (1414->74), div. (0->0), fcn. (912->6), ass. (0->41)
t447 = sin(qJ(2));
t449 = cos(qJ(2));
t439 = -qJD(2) + qJD(4);
t437 = t439 ^ 2;
t438 = qJDD(2) - qJDD(4);
t446 = sin(qJ(4));
t448 = cos(qJ(4));
t451 = t446 * t437 + t448 * t438;
t453 = -t448 * t437 + t446 * t438;
t412 = t447 * t451 - t449 * t453;
t444 = sin(pkin(6));
t445 = cos(pkin(6));
t456 = t447 * t453 + t449 * t451;
t461 = -t444 * t412 + t445 * t456;
t460 = t445 * t412 + t444 * t456;
t450 = qJD(2) ^ 2;
t430 = t447 * qJDD(2) + t449 * t450;
t431 = t449 * qJDD(2) - t447 * t450;
t417 = -t444 * t430 + t445 * t431;
t457 = t445 * t430 + t444 * t431;
t432 = t444 * g(1) - t445 * g(2);
t433 = -t445 * g(1) - t444 * g(2);
t420 = t447 * t432 + t449 * t433;
t419 = t449 * t432 - t447 * t433;
t452 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t420;
t416 = -qJDD(2) * pkin(2) - t450 * qJ(3) + qJDD(3) - t419;
t442 = g(3) - qJDD(1);
t411 = -t450 * pkin(2) + t452;
t410 = -qJDD(2) * pkin(3) + t416;
t409 = (-pkin(3) - pkin(2)) * t450 + t452;
t408 = -t447 * t419 + t449 * t420;
t407 = t449 * t419 + t447 * t420;
t406 = t449 * t411 + t447 * t416;
t405 = t447 * t411 - t449 * t416;
t404 = t448 * t409 + t446 * t410;
t403 = -t446 * t409 + t448 * t410;
t402 = -t446 * t403 + t448 * t404;
t401 = t448 * t403 + t446 * t404;
t400 = t447 * t401 + t449 * t402;
t399 = -t449 * t401 + t447 * t402;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t444 * t432 + t445 * t433, 0, 0, 0, 0, 0, 0, -t457, -t417, 0, -t444 * t407 + t445 * t408, 0, 0, 0, 0, 0, 0, -t457, 0, t417, -t444 * t405 + t445 * t406, 0, 0, 0, 0, 0, 0, -t460, t461, 0, -t444 * t399 + t445 * t400; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t445 * t432 + t444 * t433, 0, 0, 0, 0, 0, 0, t417, -t457, 0, t445 * t407 + t444 * t408, 0, 0, 0, 0, 0, 0, t417, 0, t457, t445 * t405 + t444 * t406, 0, 0, 0, 0, 0, 0, t461, t460, 0, t445 * t399 + t444 * t400; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t442, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t442, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t442, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t442; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t433, 0, 0, 0, 0, 0, 0, -t430, -t431, 0, t408, 0, 0, 0, 0, 0, 0, -t430, 0, t431, t406, 0, 0, 0, 0, 0, 0, -t412, t456, 0, t400; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t432, 0, 0, 0, 0, 0, 0, t431, -t430, 0, t407, 0, 0, 0, 0, 0, 0, t431, 0, t430, t405, 0, 0, 0, 0, 0, 0, t456, t412, 0, t399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t442, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t442, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t442, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t442; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t450, -qJDD(2), 0, t420, 0, 0, 0, 0, 0, 0, -t450, 0, qJDD(2), t411, 0, 0, 0, 0, 0, 0, t453, t451, 0, t402; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t450, 0, t419, 0, 0, 0, 0, 0, 0, qJDD(2), 0, t450, -t416, 0, 0, 0, 0, 0, 0, t451, -t453, 0, -t401; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t442, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t442, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t442; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t450, 0, qJDD(2), t411, 0, 0, 0, 0, 0, 0, t453, t451, 0, t402; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t442, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t442; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), 0, -t450, t416, 0, 0, 0, 0, 0, 0, -t451, t453, 0, t401; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t437, t438, 0, t404; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t438, -t437, 0, t403; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t442;];
f_new_reg  = t1;
