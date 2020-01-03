% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PPRR3
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PPRR3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR3_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:32
% EndTime: 2019-12-31 16:17:32
% DurationCPUTime: 0.50s
% Computational Cost: add. (598->92), mult. (1152->98), div. (0->0), fcn. (810->6), ass. (0->56)
t428 = sin(pkin(6));
t429 = cos(pkin(6));
t416 = t428 * g(1) - t429 * g(2);
t411 = -qJDD(2) + t416;
t417 = -t429 * g(1) - t428 * g(2);
t432 = sin(qJ(3));
t434 = cos(qJ(3));
t396 = -t432 * t411 + t434 * t417;
t431 = sin(qJ(4));
t424 = t431 ^ 2;
t433 = cos(qJ(4));
t425 = t433 ^ 2;
t442 = t424 + t425;
t441 = qJD(3) * qJD(4);
t440 = t431 * qJDD(3);
t439 = t433 * qJDD(3);
t395 = -t434 * t411 - t432 * t417;
t436 = qJD(3) ^ 2;
t413 = t432 * qJDD(3) + t434 * t436;
t414 = -t434 * qJDD(3) + t432 * t436;
t438 = -t428 * t413 + t429 * t414;
t437 = t429 * t413 + t428 * t414;
t435 = qJD(4) ^ 2;
t426 = g(3) - qJDD(1);
t422 = t431 * t436 * t433;
t421 = -t425 * t436 - t435;
t420 = -t424 * t436 - t435;
t419 = -qJDD(4) + t422;
t418 = qJDD(4) + t422;
t415 = t442 * t436;
t412 = t442 * qJDD(3);
t410 = -0.2e1 * t431 * t441 + t439;
t409 = 0.2e1 * t433 * t441 + t440;
t407 = t429 * t417;
t406 = t428 * t417;
t402 = t433 * t419 - t431 * t420;
t401 = -t431 * t418 + t433 * t421;
t400 = t431 * t419 + t433 * t420;
t399 = t433 * t418 + t431 * t421;
t398 = t434 * t412 - t432 * t415;
t397 = t432 * t412 + t434 * t415;
t394 = -t436 * pkin(3) + qJDD(3) * pkin(5) + t396;
t393 = -qJDD(3) * pkin(3) - t436 * pkin(5) - t395;
t392 = t434 * t402 + t432 * t409;
t391 = t434 * t401 - t432 * t410;
t390 = t432 * t402 - t434 * t409;
t389 = t432 * t401 + t434 * t410;
t388 = t433 * t394 + t431 * t426;
t387 = -t431 * t394 + t433 * t426;
t386 = -t432 * t395 + t434 * t396;
t385 = t434 * t395 + t432 * t396;
t384 = -t431 * t387 + t433 * t388;
t383 = t433 * t387 + t431 * t388;
t382 = t434 * t384 + t432 * t393;
t381 = t432 * t384 - t434 * t393;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t428 * t416 + t407, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t428 * t411 + t407, 0, 0, 0, 0, 0, 0, -t437, t438, 0, t428 * t385 + t429 * t386, 0, 0, 0, 0, 0, 0, t428 * t389 + t429 * t391, t428 * t390 + t429 * t392, t428 * t397 + t429 * t398, t428 * t381 + t429 * t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t429 * t416 + t406, 0, 0, 0, 0, 0, 0, 0, 0, 0, t429 * t411 + t406, 0, 0, 0, 0, 0, 0, t438, t437, 0, -t429 * t385 + t428 * t386, 0, 0, 0, 0, 0, 0, -t429 * t389 + t428 * t391, -t429 * t390 + t428 * t392, -t429 * t397 + t428 * t398, -t429 * t381 + t428 * t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426, 0, 0, 0, 0, 0, 0, -t399, -t400, 0, -t383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t417, 0, 0, 0, 0, 0, 0, 0, 0, 0, t417, 0, 0, 0, 0, 0, 0, -t413, t414, 0, t386, 0, 0, 0, 0, 0, 0, t391, t392, t398, t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t416, 0, 0, 0, 0, 0, 0, 0, 0, 0, t411, 0, 0, 0, 0, 0, 0, t414, t413, 0, -t385, 0, 0, 0, 0, 0, 0, -t389, -t390, -t397, -t381; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426, 0, 0, 0, 0, 0, 0, -t399, -t400, 0, -t383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t417, 0, 0, 0, 0, 0, 0, -t413, t414, 0, t386, 0, 0, 0, 0, 0, 0, t391, t392, t398, t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426, 0, 0, 0, 0, 0, 0, -t399, -t400, 0, -t383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t411, 0, 0, 0, 0, 0, 0, -t414, -t413, 0, t385, 0, 0, 0, 0, 0, 0, t389, t390, t397, t381; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t436, -qJDD(3), 0, t396, 0, 0, 0, 0, 0, 0, t401, t402, t412, t384; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t436, 0, t395, 0, 0, 0, 0, 0, 0, t410, -t409, t415, -t393; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t426, 0, 0, 0, 0, 0, 0, t399, t400, 0, t383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t421, t419, t439, t388; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t418, t420, -t440, t387; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t410, t409, -t415, t393;];
f_new_reg = t1;
