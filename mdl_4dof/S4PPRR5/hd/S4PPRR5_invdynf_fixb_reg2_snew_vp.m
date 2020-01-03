% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PPRR5
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PPRR5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR5_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:57
% EndTime: 2019-12-31 16:19:57
% DurationCPUTime: 0.45s
% Computational Cost: add. (571->84), mult. (1051->94), div. (0->0), fcn. (724->6), ass. (0->54)
t438 = sin(pkin(6));
t439 = cos(pkin(6));
t424 = t438 * g(1) - t439 * g(2);
t419 = -qJDD(2) + t424;
t436 = -g(3) + qJDD(1);
t442 = sin(qJ(3));
t444 = cos(qJ(3));
t412 = -t442 * t419 + t444 * t436;
t441 = sin(qJ(4));
t434 = t441 ^ 2;
t443 = cos(qJ(4));
t435 = t443 ^ 2;
t450 = t434 + t435;
t449 = qJD(3) * qJD(4);
t448 = t441 * qJDD(3);
t447 = t443 * qJDD(3);
t411 = -t444 * t419 - t442 * t436;
t446 = qJD(3) ^ 2;
t445 = qJD(4) ^ 2;
t430 = t441 * t446 * t443;
t429 = -t435 * t446 - t445;
t428 = -t434 * t446 - t445;
t427 = -qJDD(4) + t430;
t426 = qJDD(4) + t430;
t425 = t439 * g(1) + t438 * g(2);
t423 = t450 * t446;
t422 = t444 * qJDD(3) - t442 * t446;
t421 = t442 * qJDD(3) + t444 * t446;
t420 = t450 * qJDD(3);
t418 = -0.2e1 * t441 * t449 + t447;
t417 = 0.2e1 * t443 * t449 + t448;
t416 = t439 * t425;
t415 = t438 * t425;
t410 = t443 * t427 - t441 * t428;
t409 = -t441 * t426 + t443 * t429;
t408 = t441 * t427 + t443 * t428;
t407 = t443 * t426 + t441 * t429;
t406 = t444 * t420 - t442 * t423;
t405 = t442 * t420 + t444 * t423;
t404 = -t446 * pkin(3) + qJDD(3) * pkin(5) + t412;
t403 = -qJDD(3) * pkin(3) - t446 * pkin(5) - t411;
t402 = t444 * t410 + t442 * t417;
t401 = t444 * t409 - t442 * t418;
t400 = t442 * t410 - t444 * t417;
t399 = t442 * t409 + t444 * t418;
t398 = t443 * t404 - t441 * t425;
t397 = -t441 * t404 - t443 * t425;
t396 = -t442 * t411 + t444 * t412;
t395 = t444 * t411 + t442 * t412;
t394 = -t441 * t397 + t443 * t398;
t393 = t443 * t397 + t441 * t398;
t392 = t444 * t394 + t442 * t403;
t391 = t442 * t394 - t444 * t403;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t438 * t424 - t416, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t438 * t419 - t416, 0, 0, 0, 0, 0, 0, t438 * t422, -t438 * t421, 0, t438 * t395 - t416, 0, 0, 0, 0, 0, 0, t438 * t399 + t439 * t407, t438 * t400 + t439 * t408, t438 * t405, t438 * t391 + t439 * t393; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t439 * t424 - t415, 0, 0, 0, 0, 0, 0, 0, 0, 0, t439 * t419 - t415, 0, 0, 0, 0, 0, 0, -t439 * t422, t439 * t421, 0, -t439 * t395 - t415, 0, 0, 0, 0, 0, 0, -t439 * t399 + t438 * t407, -t439 * t400 + t438 * t408, -t439 * t405, -t439 * t391 + t438 * t393; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t436, 0, 0, 0, 0, 0, 0, 0, 0, 0, t436, 0, 0, 0, 0, 0, 0, -t421, -t422, 0, t396, 0, 0, 0, 0, 0, 0, t401, t402, t406, t392; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t425, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t425, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t425, 0, 0, 0, 0, 0, 0, t407, t408, 0, t393; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t424, 0, 0, 0, 0, 0, 0, 0, 0, 0, t419, 0, 0, 0, 0, 0, 0, -t422, t421, 0, -t395, 0, 0, 0, 0, 0, 0, -t399, -t400, -t405, -t391; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t436, 0, 0, 0, 0, 0, 0, 0, 0, 0, t436, 0, 0, 0, 0, 0, 0, -t421, -t422, 0, t396, 0, 0, 0, 0, 0, 0, t401, t402, t406, t392; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t436, 0, 0, 0, 0, 0, 0, -t421, -t422, 0, t396, 0, 0, 0, 0, 0, 0, t401, t402, t406, t392; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t425, 0, 0, 0, 0, 0, 0, 0, 0, 0, t425, 0, 0, 0, 0, 0, 0, -t407, -t408, 0, -t393; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t419, 0, 0, 0, 0, 0, 0, t422, -t421, 0, t395, 0, 0, 0, 0, 0, 0, t399, t400, t405, t391; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t446, -qJDD(3), 0, t412, 0, 0, 0, 0, 0, 0, t409, t410, t420, t394; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t446, 0, t411, 0, 0, 0, 0, 0, 0, t418, -t417, t423, -t403; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t425, 0, 0, 0, 0, 0, 0, t407, t408, 0, t393; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t429, t427, t447, t398; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t426, t428, -t448, t397; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t418, t417, -t423, t403;];
f_new_reg = t1;
