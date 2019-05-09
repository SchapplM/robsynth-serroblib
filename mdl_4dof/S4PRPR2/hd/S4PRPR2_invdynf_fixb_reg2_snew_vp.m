% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRPR2
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
%   pkin=[a2,a3,a4,d2,d4,theta3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:00
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRPR2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR2_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:59:58
% EndTime: 2019-05-04 18:59:59
% DurationCPUTime: 0.58s
% Computational Cost: add. (1237->71), mult. (1678->66), div. (0->0), fcn. (1120->6), ass. (0->47)
t420 = qJD(2) + qJD(4);
t418 = t420 ^ 2;
t419 = qJDD(2) + qJDD(4);
t425 = sin(qJ(4));
t427 = cos(qJ(4));
t402 = t418 * t425 - t419 * t427;
t423 = sin(pkin(6));
t424 = cos(pkin(6));
t430 = -t418 * t427 - t419 * t425;
t387 = t402 * t424 - t423 * t430;
t426 = sin(qJ(2));
t428 = cos(qJ(2));
t434 = t402 * t423 + t424 * t430;
t437 = t387 * t426 + t428 * t434;
t379 = t387 * t428 - t426 * t434;
t422 = -g(2) + qJDD(1);
t407 = g(1) * t426 + t422 * t428;
t404 = qJDD(2) * pkin(2) + t407;
t408 = -g(1) * t428 + t422 * t426;
t429 = qJD(2) ^ 2;
t405 = -pkin(2) * t429 + t408;
t390 = t404 * t423 + t405 * t424;
t389 = t404 * t424 - t405 * t423;
t409 = qJDD(2) * t424 - t423 * t429;
t410 = -qJDD(2) * t423 - t424 * t429;
t431 = -t409 * t426 + t410 * t428;
t394 = t409 * t428 + t410 * t426;
t421 = g(3) - qJDD(3);
t412 = qJDD(2) * t428 - t426 * t429;
t411 = -qJDD(2) * t426 - t428 * t429;
t392 = -t407 * t426 + t408 * t428;
t391 = t407 * t428 + t408 * t426;
t384 = -pkin(3) * t429 + t390;
t383 = qJDD(2) * pkin(3) + t389;
t382 = -t389 * t423 + t390 * t424;
t381 = t389 * t424 + t390 * t423;
t376 = t383 * t425 + t384 * t427;
t375 = t383 * t427 - t384 * t425;
t374 = -t381 * t426 + t382 * t428;
t373 = t381 * t428 + t382 * t426;
t372 = -t375 * t425 + t376 * t427;
t371 = t375 * t427 + t376 * t425;
t370 = -t371 * t423 + t372 * t424;
t369 = t371 * t424 + t372 * t423;
t368 = -t369 * t426 + t428 * t370;
t367 = t369 * t428 + t370 * t426;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t411, -t412, 0, t392, 0, 0, 0, 0, 0, 0, t431, -t394, 0, t374, 0, 0, 0, 0, 0, 0, t437, t379, 0, t368; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t422, 0, 0, 0, 0, 0, 0, t412, t411, 0, t391, 0, 0, 0, 0, 0, 0, t394, t431, 0, t373, 0, 0, 0, 0, 0, 0, -t379, t437, 0, t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t421, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t421; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t411, -t412, 0, t392, 0, 0, 0, 0, 0, 0, t431, -t394, 0, t374, 0, 0, 0, 0, 0, 0, t437, t379, 0, t368; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t421, 0, 0, 0, 0, 0, 0, 0, 0, 0, t421; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t422, 0, 0, 0, 0, 0, 0, t412, t411, 0, t391, 0, 0, 0, 0, 0, 0, t394, t431, 0, t373, 0, 0, 0, 0, 0, 0, -t379, t437, 0, t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t429, -qJDD(2), 0, t408, 0, 0, 0, 0, 0, 0, t410, -t409, 0, t382, 0, 0, 0, 0, 0, 0, t434, t387, 0, t370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t429, 0, t407, 0, 0, 0, 0, 0, 0, t409, t410, 0, t381, 0, 0, 0, 0, 0, 0, -t387, t434, 0, t369; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t421, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t421; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t429, -qJDD(2), 0, t390, 0, 0, 0, 0, 0, 0, t430, t402, 0, t372; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t429, 0, t389, 0, 0, 0, 0, 0, 0, -t402, t430, 0, t371; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t421, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t421; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t418, -t419, 0, t376; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t419, -t418, 0, t375; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t421;];
f_new_reg  = t1;
