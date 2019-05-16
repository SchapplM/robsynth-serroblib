% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRRP1
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRRP1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP1_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:02:32
% EndTime: 2019-05-04 19:02:33
% DurationCPUTime: 0.63s
% Computational Cost: add. (1089->79), mult. (1670->78), div. (0->0), fcn. (1280->6), ass. (0->42)
t458 = (qJD(2) + qJD(3));
t456 = t458 ^ 2;
t457 = qJDD(2) + qJDD(3);
t462 = sin(qJ(3));
t464 = cos(qJ(3));
t442 = t464 * t456 + t462 * t457;
t445 = t462 * t456 - t464 * t457;
t463 = sin(qJ(2));
t465 = cos(qJ(2));
t429 = t465 * t442 - t463 * t445;
t433 = t463 * t442 + t465 * t445;
t460 = sin(pkin(6));
t461 = cos(pkin(6));
t473 = t461 * t429 - t460 * t433;
t474 = t460 * t429 + t461 * t433;
t451 = t460 * g(1) - t461 * g(2);
t452 = -t461 * g(1) - t460 * g(2);
t440 = t463 * t451 + t465 * t452;
t466 = qJD(2) ^ 2;
t438 = -t466 * pkin(2) + t440;
t439 = t465 * t451 - t463 * t452;
t467 = qJDD(2) * pkin(2) + t439;
t426 = t464 * t438 + t462 * t467;
t425 = -t462 * t438 + t464 * t467;
t449 = t465 * qJDD(2) - t463 * t466;
t450 = -t463 * qJDD(2) - t465 * t466;
t469 = -t460 * t449 + t461 * t450;
t468 = t461 * t449 + t460 * t450;
t459 = -g(3) + qJDD(1);
t428 = -t463 * t439 + t465 * t440;
t427 = t465 * t439 + t463 * t440;
t424 = -t457 * pkin(3) - t456 * qJ(4) + qJDD(4) - t425;
t421 = -t456 * pkin(3) + t457 * qJ(4) + (2 * qJD(4) * t458) + t426;
t420 = -t462 * t425 + t464 * t426;
t419 = t464 * t425 + t462 * t426;
t418 = t464 * t421 + t462 * t424;
t417 = t462 * t421 - t464 * t424;
t416 = -t463 * t419 + t465 * t420;
t415 = t465 * t419 + t463 * t420;
t414 = -t463 * t417 + t465 * t418;
t413 = t465 * t417 + t463 * t418;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t460 * t451 + t461 * t452, 0, 0, 0, 0, 0, 0, t469, -t468, 0, -t460 * t427 + t461 * t428, 0, 0, 0, 0, 0, 0, -t473, t474, 0, -t460 * t415 + t461 * t416, 0, 0, 0, 0, 0, 0, -t473, 0, -t474, -t460 * t413 + t461 * t414; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t461 * t451 + t460 * t452, 0, 0, 0, 0, 0, 0, t468, t469, 0, t461 * t427 + t460 * t428, 0, 0, 0, 0, 0, 0, -t474, -t473, 0, t461 * t415 + t460 * t416, 0, 0, 0, 0, 0, 0, -t474, 0, t473, t461 * t413 + t460 * t414; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t459, 0, 0, 0, 0, 0, 0, 0, 0, 0, t459, 0, 0, 0, 0, 0, 0, 0, 0, 0, t459, 0, 0, 0, 0, 0, 0, 0, 0, 0, t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t452, 0, 0, 0, 0, 0, 0, t450, -t449, 0, t428, 0, 0, 0, 0, 0, 0, -t429, t433, 0, t416, 0, 0, 0, 0, 0, 0, -t429, 0, -t433, t414; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t451, 0, 0, 0, 0, 0, 0, t449, t450, 0, t427, 0, 0, 0, 0, 0, 0, -t433, -t429, 0, t415, 0, 0, 0, 0, 0, 0, -t433, 0, t429, t413; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t459, 0, 0, 0, 0, 0, 0, 0, 0, 0, t459, 0, 0, 0, 0, 0, 0, 0, 0, 0, t459, 0, 0, 0, 0, 0, 0, 0, 0, 0, t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t466, -qJDD(2), 0, t440, 0, 0, 0, 0, 0, 0, -t442, t445, 0, t420, 0, 0, 0, 0, 0, 0, -t442, 0, -t445, t418; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t466, 0, t439, 0, 0, 0, 0, 0, 0, -t445, -t442, 0, t419, 0, 0, 0, 0, 0, 0, -t445, 0, t442, t417; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t459, 0, 0, 0, 0, 0, 0, 0, 0, 0, t459, 0, 0, 0, 0, 0, 0, 0, 0, 0, t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t456, -t457, 0, t426, 0, 0, 0, 0, 0, 0, -t456, 0, t457, t421; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t457, -t456, 0, t425, 0, 0, 0, 0, 0, 0, t457, 0, t456, -t424; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t459, 0, 0, 0, 0, 0, 0, 0, 0, 0, t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t456, 0, t457, t421; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t457, 0, -t456, t424;];
f_new_reg  = t1;
