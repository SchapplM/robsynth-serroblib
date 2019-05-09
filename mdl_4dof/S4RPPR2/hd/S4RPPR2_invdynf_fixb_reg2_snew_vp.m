% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPPR2
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPPR2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR2_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:11:22
% EndTime: 2019-05-04 19:11:23
% DurationCPUTime: 0.59s
% Computational Cost: add. (1417->89), mult. (2170->76), div. (0->0), fcn. (920->6), ass. (0->45)
t389 = -qJD(1) + qJD(4);
t387 = t389 ^ 2;
t388 = qJDD(1) - qJDD(4);
t394 = sin(qJ(4));
t396 = cos(qJ(4));
t375 = t394 * t387 + t396 * t388;
t392 = sin(pkin(6));
t393 = cos(pkin(6));
t401 = -t396 * t387 + t394 * t388;
t364 = t393 * t375 - t392 * t401;
t395 = sin(qJ(1));
t397 = cos(qJ(1));
t407 = t392 * t375 + t393 * t401;
t411 = t395 * t364 - t397 * t407;
t410 = t397 * t364 + t395 * t407;
t404 = -pkin(1) - pkin(2);
t398 = qJD(1) ^ 2;
t383 = -t397 * g(1) - t395 * g(2);
t400 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t383;
t369 = t404 * t398 + t400;
t382 = t395 * g(1) - t397 * g(2);
t399 = -t398 * qJ(2) + qJDD(2) - t382;
t370 = t404 * qJDD(1) + t399;
t361 = t393 * t369 + t392 * t370;
t360 = -t392 * t369 + t393 * t370;
t378 = -t392 * qJDD(1) + t393 * t398;
t379 = t393 * qJDD(1) + t392 * t398;
t403 = -t395 * t378 + t397 * t379;
t402 = t397 * t378 + t395 * t379;
t391 = g(3) + qJDD(3);
t381 = t397 * qJDD(1) - t395 * t398;
t380 = t395 * qJDD(1) + t397 * t398;
t372 = qJDD(1) * pkin(1) - t399;
t371 = -t398 * pkin(1) + t400;
t359 = -t398 * pkin(3) + t361;
t358 = -qJDD(1) * pkin(3) + t360;
t357 = -t392 * t360 + t393 * t361;
t356 = t393 * t360 + t392 * t361;
t355 = t394 * t358 + t396 * t359;
t354 = t396 * t358 - t394 * t359;
t353 = -t394 * t354 + t396 * t355;
t352 = t396 * t354 + t394 * t355;
t351 = -t392 * t352 + t393 * t353;
t350 = t393 * t352 + t392 * t353;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t380, -t381, 0, -t395 * t382 + t397 * t383, 0, 0, 0, 0, 0, 0, -t380, 0, t381, t397 * t371 - t395 * t372, 0, 0, 0, 0, 0, 0, -t402, t403, 0, t395 * t356 + t397 * t357, 0, 0, 0, 0, 0, 0, -t411, t410, 0, t395 * t350 + t397 * t351; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t381, -t380, 0, t397 * t382 + t395 * t383, 0, 0, 0, 0, 0, 0, t381, 0, t380, t395 * t371 + t397 * t372, 0, 0, 0, 0, 0, 0, t403, t402, 0, -t397 * t356 + t395 * t357, 0, 0, 0, 0, 0, 0, t410, t411, 0, -t397 * t350 + t395 * t351; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t391, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t391; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t398, -qJDD(1), 0, t383, 0, 0, 0, 0, 0, 0, -t398, 0, qJDD(1), t371, 0, 0, 0, 0, 0, 0, -t378, t379, 0, t357, 0, 0, 0, 0, 0, 0, t407, t364, 0, t351; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t398, 0, t382, 0, 0, 0, 0, 0, 0, qJDD(1), 0, t398, t372, 0, 0, 0, 0, 0, 0, t379, t378, 0, -t356, 0, 0, 0, 0, 0, 0, t364, -t407, 0, -t350; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t391, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t391; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t398, 0, qJDD(1), t371, 0, 0, 0, 0, 0, 0, -t378, t379, 0, t357, 0, 0, 0, 0, 0, 0, t407, t364, 0, t351; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t391, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t391; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t398, -t372, 0, 0, 0, 0, 0, 0, -t379, -t378, 0, t356, 0, 0, 0, 0, 0, 0, -t364, t407, 0, t350; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t398, qJDD(1), 0, t361, 0, 0, 0, 0, 0, 0, t401, t375, 0, t353; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), -t398, 0, t360, 0, 0, 0, 0, 0, 0, -t375, t401, 0, t352; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t391, 0, 0, 0, 0, 0, 0, 0, 0, 0, t391; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t387, t388, 0, t355; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t388, -t387, 0, t354; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t391;];
f_new_reg  = t1;
