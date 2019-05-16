% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPPR1
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
% Datum: 2019-05-04 19:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPPR1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR1_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:08:36
% EndTime: 2019-05-04 19:08:37
% DurationCPUTime: 0.56s
% Computational Cost: add. (1073->92), mult. (1722->80), div. (0->0), fcn. (920->6), ass. (0->45)
t465 = sin(pkin(6));
t466 = cos(pkin(6));
t460 = -qJD(1) + qJD(4);
t458 = t460 ^ 2;
t459 = qJDD(1) - qJDD(4);
t467 = sin(qJ(4));
t469 = cos(qJ(4));
t473 = t467 * t458 + t469 * t459;
t475 = -t469 * t458 + t467 * t459;
t430 = t465 * t473 - t466 * t475;
t468 = sin(qJ(1));
t470 = cos(qJ(1));
t478 = t465 * t475 + t466 * t473;
t483 = -t468 * t430 + t470 * t478;
t482 = t470 * t430 + t468 * t478;
t471 = qJD(1) ^ 2;
t448 = t465 * qJDD(1) + t466 * t471;
t449 = t466 * qJDD(1) - t465 * t471;
t479 = t470 * t448 + t468 * t449;
t436 = -t468 * t448 + t470 * t449;
t453 = -t470 * g(1) - t468 * g(2);
t447 = -t471 * pkin(1) + t453;
t452 = t468 * g(1) - t470 * g(2);
t472 = qJDD(1) * pkin(1) + t452;
t435 = t466 * t447 + t465 * t472;
t434 = -t465 * t447 + t466 * t472;
t474 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t435;
t429 = -qJDD(1) * pkin(2) - t471 * qJ(3) + qJDD(3) - t434;
t463 = g(3) - qJDD(2);
t451 = -t468 * qJDD(1) - t470 * t471;
t450 = t470 * qJDD(1) - t468 * t471;
t428 = -t471 * pkin(2) + t474;
t427 = -qJDD(1) * pkin(3) + t429;
t426 = (-pkin(2) - pkin(3)) * t471 + t474;
t425 = -t465 * t434 + t466 * t435;
t424 = t466 * t434 + t465 * t435;
t423 = t466 * t428 + t465 * t429;
t422 = t465 * t428 - t466 * t429;
t421 = t469 * t426 + t467 * t427;
t420 = -t467 * t426 + t469 * t427;
t419 = -t467 * t420 + t469 * t421;
t418 = t469 * t420 + t467 * t421;
t417 = t465 * t418 + t466 * t419;
t416 = -t466 * t418 + t465 * t419;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t451, -t450, 0, -t468 * t452 + t470 * t453, 0, 0, 0, 0, 0, 0, -t479, -t436, 0, -t468 * t424 + t470 * t425, 0, 0, 0, 0, 0, 0, -t479, 0, t436, -t468 * t422 + t470 * t423, 0, 0, 0, 0, 0, 0, -t482, t483, 0, -t468 * t416 + t470 * t417; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t450, t451, 0, t470 * t452 + t468 * t453, 0, 0, 0, 0, 0, 0, t436, -t479, 0, t470 * t424 + t468 * t425, 0, 0, 0, 0, 0, 0, t436, 0, t479, t470 * t422 + t468 * t423, 0, 0, 0, 0, 0, 0, t483, t482, 0, t470 * t416 + t468 * t417; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t463, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t463, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t463; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t471, -qJDD(1), 0, t453, 0, 0, 0, 0, 0, 0, -t448, -t449, 0, t425, 0, 0, 0, 0, 0, 0, -t448, 0, t449, t423, 0, 0, 0, 0, 0, 0, -t430, t478, 0, t417; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t471, 0, t452, 0, 0, 0, 0, 0, 0, t449, -t448, 0, t424, 0, 0, 0, 0, 0, 0, t449, 0, t448, t422, 0, 0, 0, 0, 0, 0, t478, t430, 0, t416; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t463, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t463, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t463; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t471, -qJDD(1), 0, t435, 0, 0, 0, 0, 0, 0, -t471, 0, qJDD(1), t428, 0, 0, 0, 0, 0, 0, t475, t473, 0, t419; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t471, 0, t434, 0, 0, 0, 0, 0, 0, qJDD(1), 0, t471, -t429, 0, 0, 0, 0, 0, 0, t473, -t475, 0, -t418; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t463, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t463, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t463; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t471, 0, qJDD(1), t428, 0, 0, 0, 0, 0, 0, t475, t473, 0, t419; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t463, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t463; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t471, t429, 0, 0, 0, 0, 0, 0, -t473, t475, 0, t418; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t458, t459, 0, t421; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t459, -t458, 0, t420; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t463;];
f_new_reg  = t1;
