% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PPRR2
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
%   pkin=[a2,a3,a4,d3,d4,theta2]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PPRR2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR2_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:50:59
% EndTime: 2019-05-04 18:51:00
% DurationCPUTime: 0.54s
% Computational Cost: add. (999->62), mult. (1310->60), div. (0->0), fcn. (1104->6), ass. (0->43)
t391 = qJD(3) + qJD(4);
t389 = t391 ^ 2;
t390 = qJDD(3) + qJDD(4);
t396 = sin(qJ(4));
t398 = cos(qJ(4));
t381 = t396 * t389 - t398 * t390;
t397 = sin(qJ(3));
t399 = cos(qJ(3));
t401 = -t398 * t389 - t396 * t390;
t364 = t399 * t381 - t397 * t401;
t394 = sin(pkin(6));
t395 = cos(pkin(6));
t405 = t397 * t381 + t399 * t401;
t408 = t394 * t364 + t395 * t405;
t356 = t395 * t364 - t394 * t405;
t393 = -g(2) + qJDD(1);
t384 = t394 * g(1) + t395 * t393;
t385 = -t395 * g(1) + t394 * t393;
t371 = t397 * t384 + t399 * t385;
t370 = t399 * t384 - t397 * t385;
t400 = qJD(3) ^ 2;
t386 = t399 * qJDD(3) - t397 * t400;
t387 = -t397 * qJDD(3) - t399 * t400;
t402 = -t394 * t386 + t395 * t387;
t373 = t395 * t386 + t394 * t387;
t392 = g(3) - qJDD(2);
t369 = -t394 * t384 + t395 * t385;
t368 = t395 * t384 + t394 * t385;
t367 = -t400 * pkin(3) + t371;
t366 = qJDD(3) * pkin(3) + t370;
t361 = -t397 * t370 + t399 * t371;
t360 = t399 * t370 + t397 * t371;
t359 = t396 * t366 + t398 * t367;
t358 = t398 * t366 - t396 * t367;
t353 = -t394 * t360 + t395 * t361;
t352 = t395 * t360 + t394 * t361;
t351 = -t396 * t358 + t398 * t359;
t350 = t398 * t358 + t396 * t359;
t349 = -t397 * t350 + t399 * t351;
t348 = t399 * t350 + t397 * t351;
t347 = -t394 * t348 + t395 * t349;
t346 = t395 * t348 + t394 * t349;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t369, 0, 0, 0, 0, 0, 0, t402, -t373, 0, t353, 0, 0, 0, 0, 0, 0, t408, t356, 0, t347; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t393, 0, 0, 0, 0, 0, 0, 0, 0, 0, t368, 0, 0, 0, 0, 0, 0, t373, t402, 0, t352, 0, 0, 0, 0, 0, 0, -t356, t408, 0, t346; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t392, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t392, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t392; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t369, 0, 0, 0, 0, 0, 0, t402, -t373, 0, t353, 0, 0, 0, 0, 0, 0, t408, t356, 0, t347; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t392, 0, 0, 0, 0, 0, 0, 0, 0, 0, t392, 0, 0, 0, 0, 0, 0, 0, 0, 0, t392; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t393, 0, 0, 0, 0, 0, 0, 0, 0, 0, t368, 0, 0, 0, 0, 0, 0, t373, t402, 0, t352, 0, 0, 0, 0, 0, 0, -t356, t408, 0, t346; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t385, 0, 0, 0, 0, 0, 0, t387, -t386, 0, t361, 0, 0, 0, 0, 0, 0, t405, t364, 0, t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t384, 0, 0, 0, 0, 0, 0, t386, t387, 0, t360, 0, 0, 0, 0, 0, 0, -t364, t405, 0, t348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t392, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t392, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t392; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t400, -qJDD(3), 0, t371, 0, 0, 0, 0, 0, 0, t401, t381, 0, t351; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t400, 0, t370, 0, 0, 0, 0, 0, 0, -t381, t401, 0, t350; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t392, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t392; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t389, -t390, 0, t359; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t390, -t389, 0, t358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t392;];
f_new_reg  = t1;
