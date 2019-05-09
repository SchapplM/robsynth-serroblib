% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRPP1
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRPP1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP1_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:53:31
% EndTime: 2019-05-04 18:53:31
% DurationCPUTime: 0.43s
% Computational Cost: add. (393->64), mult. (677->55), div. (0->0), fcn. (448->4), ass. (0->25)
t401 = sin(qJ(2));
t402 = cos(qJ(2));
t403 = qJD(2) ^ 2;
t387 = qJDD(2) * t402 - t401 * t403;
t388 = qJDD(2) * t401 + t402 * t403;
t398 = sin(pkin(5));
t399 = cos(pkin(5));
t376 = t387 * t398 + t388 * t399;
t405 = t387 * t399 - t388 * t398;
t389 = g(1) * t398 - g(2) * t399;
t390 = -g(1) * t399 - g(2) * t398;
t380 = t389 * t402 - t390 * t401;
t381 = t389 * t401 + t390 * t402;
t375 = -qJDD(2) * pkin(2) - qJ(3) * t403 + qJDD(3) - t380;
t374 = -pkin(2) * t403 + qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t381;
t396 = -g(3) + qJDD(1);
t373 = -qJ(4) * t403 + qJDD(4) + t374;
t372 = -qJDD(2) * qJ(4) - (2 * qJD(4) * qJD(2)) + t375;
t371 = -t380 * t401 + t381 * t402;
t370 = t380 * t402 + t381 * t401;
t369 = t374 * t402 + t375 * t401;
t368 = t374 * t401 - t375 * t402;
t367 = t372 * t401 + t373 * t402;
t366 = -t372 * t402 + t373 * t401;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t389 * t398 + t390 * t399, 0, 0, 0, 0, 0, 0, -t376, -t405, 0, -t370 * t398 + t371 * t399, 0, 0, 0, 0, 0, 0, 0, t376, t405, -t368 * t398 + t369 * t399, 0, 0, 0, 0, 0, 0, 0, t405, -t376, -t366 * t398 + t367 * t399; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t389 * t399 + t390 * t398, 0, 0, 0, 0, 0, 0, t405, -t376, 0, t370 * t399 + t371 * t398, 0, 0, 0, 0, 0, 0, 0, -t405, t376, t368 * t399 + t369 * t398, 0, 0, 0, 0, 0, 0, 0, t376, t405, t366 * t399 + t367 * t398; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t396, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t390, 0, 0, 0, 0, 0, 0, -t388, -t387, 0, t371, 0, 0, 0, 0, 0, 0, 0, t388, t387, t369, 0, 0, 0, 0, 0, 0, 0, t387, -t388, t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t389, 0, 0, 0, 0, 0, 0, t387, -t388, 0, t370, 0, 0, 0, 0, 0, 0, 0, -t387, t388, t368, 0, 0, 0, 0, 0, 0, 0, t388, t387, t366; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t403, -qJDD(2), 0, t381, 0, 0, 0, 0, 0, 0, 0, t403, qJDD(2), t374, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t403, t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t403, 0, t380, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), t403, -t375, 0, 0, 0, 0, 0, 0, 0, t403, qJDD(2), -t372; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t403, -qJDD(2), -t374, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), t403, -t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t403, t375, 0, 0, 0, 0, 0, 0, 0, -t403, -qJDD(2), t372; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t403, -qJDD(2), t372; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t403, t373;];
f_new_reg  = t1;
