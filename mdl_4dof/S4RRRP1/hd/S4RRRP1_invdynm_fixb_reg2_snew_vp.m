% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
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
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RRRP1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP1_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_invdynm_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:24:27
% EndTime: 2019-05-04 19:24:29
% DurationCPUTime: 2.02s
% Computational Cost: add. (6935->184), mult. (8999->187), div. (0->0), fcn. (5500->6), ass. (0->119)
t377 = qJD(1) + qJD(2);
t373 = qJD(3) + t377;
t371 = t373 ^ 2;
t378 = g(3) - qJDD(4);
t365 = -qJ(4) * t371 + t378;
t382 = cos(qJ(3));
t376 = qJDD(1) + qJDD(2);
t372 = qJDD(3) + t376;
t379 = sin(qJ(3));
t409 = t379 * t372;
t351 = t382 * t371 + t409;
t435 = pkin(6) * t351;
t320 = -qJ(4) * t409 + t382 * t365 - t435;
t380 = sin(qJ(2));
t383 = cos(qJ(2));
t405 = t382 * t372;
t354 = t379 * t371 - t405;
t434 = pkin(6) * t354;
t439 = -qJ(4) * t405 - t379 * t365 + t434;
t324 = t383 * t351 - t380 * t354;
t458 = pkin(5) * t324;
t285 = t383 * t320 + t380 * t439 - t458;
t328 = t380 * t351 + t383 * t354;
t381 = sin(qJ(1));
t384 = cos(qJ(1));
t297 = t384 * t324 - t381 * t328;
t296 = pkin(4) * t297;
t457 = pkin(5) * t328;
t460 = -t380 * t320 + t383 * t439 + t457;
t471 = t384 * t285 + t381 * t460 - t296;
t339 = t382 * g(3) - t435;
t438 = -t379 * g(3) + t434;
t290 = t383 * t339 + t380 * t438 - t458;
t459 = -t380 * t339 + t383 * t438 + t457;
t470 = t384 * t290 + t381 * t459 - t296;
t447 = t381 * t324 + t384 * t328;
t467 = pkin(4) * t447;
t469 = -t381 * t285 + t384 * t460 + t467;
t468 = -t381 * t290 + t384 * t459 + t467;
t366 = t381 * g(1) - t384 * g(2);
t361 = qJDD(1) * pkin(1) + t366;
t367 = t384 * g(1) + t381 * g(2);
t386 = qJD(1) ^ 2;
t362 = -t386 * pkin(1) - t367;
t336 = t380 * t361 + t383 * t362;
t375 = t377 ^ 2;
t331 = -t375 * pkin(2) + t336;
t392 = t383 * t361 - t380 * t362;
t387 = t376 * pkin(2) + t392;
t307 = t379 * t331 - t382 * t387;
t308 = t382 * t331 + t379 * t387;
t400 = t379 * t307 + t382 * t308;
t280 = t382 * t307 - t379 * t308;
t404 = t383 * t280;
t269 = -t380 * t400 + t404;
t408 = t380 * t280;
t440 = t383 * t400 + t408;
t462 = t381 * t269 + t384 * t440;
t461 = t384 * t269 - t381 * t440;
t305 = -t371 * pkin(3) + t308;
t369 = t372 * pkin(3);
t304 = t307 - t369;
t406 = t382 * t304;
t274 = -t379 * t305 + t406;
t410 = t379 * t304;
t401 = t382 * t305 + t410;
t265 = t383 * t274 - t380 * t401;
t432 = t380 * t274 + t383 * t401;
t450 = t381 * t265 + t384 * t432;
t449 = t384 * t265 - t381 * t432;
t356 = t383 * t375 + t380 * t376;
t344 = pkin(5) * t356 - t383 * g(3);
t359 = t380 * t375 - t383 * t376;
t393 = t384 * t356 - t381 * t359;
t433 = pkin(5) * t359 - t380 * g(3);
t448 = pkin(4) * t393 + t384 * t344 - t381 * t433;
t430 = t381 * t356 + t384 * t359;
t446 = pkin(4) * t430 + t381 * t344 + t384 * t433;
t399 = t383 * t336 - t380 * t392;
t313 = -t380 * t336 - t383 * t392;
t403 = t384 * t313;
t442 = -t381 * t399 + t403;
t407 = t381 * t313;
t441 = t384 * t399 + t407;
t411 = qJ(4) * t372;
t303 = pkin(3) * t304;
t402 = -pkin(2) * t274 - t303;
t397 = -t381 * t366 - t384 * t367;
t364 = t384 * qJDD(1) - t381 * t386;
t396 = -pkin(4) * t364 - t381 * g(3);
t349 = pkin(2) * t354;
t395 = -t307 - t349;
t394 = -t307 + 0.2e1 * t369;
t391 = t384 * t366 - t381 * t367;
t390 = -t349 + t394;
t289 = -pkin(2) * t351 - t308;
t385 = pkin(1) * g(3);
t363 = t381 * qJDD(1) + t384 * t386;
t355 = -pkin(4) * t363 + t384 * g(3);
t319 = pkin(1) * t328;
t316 = -pkin(1) * t359 + t392;
t315 = -pkin(1) * t356 - t336;
t310 = pkin(1) * t313;
t309 = pkin(5) * t399 + t385;
t294 = pkin(3) * t378 + qJ(4) * t305;
t284 = -t319 + t395;
t283 = -pkin(1) * t324 + t289;
t282 = -t319 + t390;
t277 = pkin(2) * t280;
t276 = pkin(2) * g(3) + pkin(6) * t400;
t262 = pkin(6) * t274 + qJ(4) * t406 - t379 * t294;
t261 = pkin(2) * t378 + pkin(6) * t401 + qJ(4) * t410 + t382 * t294;
t260 = -pkin(1) * t269 - t277;
t259 = -pkin(1) * t265 + t402;
t258 = pkin(5) * t269 + pkin(6) * t404 - t380 * t276;
t257 = pkin(5) * t440 + pkin(6) * t408 + t383 * t276 + t385;
t256 = pkin(5) * t265 - t380 * t261 + t383 * t262;
t255 = pkin(1) * t378 + pkin(5) * t432 + t383 * t261 + t380 * t262;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t364, 0, -t363, 0, t396, -t355, -t391, -pkin(4) * t391, 0, 0, -t430, 0, -t393, 0, t446, t448, t442, pkin(4) * t442 + pkin(5) * t403 - t381 * t309, 0, 0, -t447, 0, -t297, 0, t468, -t470, t461, pkin(4) * t461 - t381 * t257 + t384 * t258, 0, 0, -t447, 0, -t297, 0, t469, -t471, t449, pkin(4) * t449 - t381 * t255 + t384 * t256; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t363, 0, t364, 0, t355, t396, t397, pkin(4) * t397, 0, 0, t393, 0, -t430, 0, -t448, t446, t441, pkin(4) * t441 + pkin(5) * t407 + t384 * t309, 0, 0, t297, 0, -t447, 0, t470, t468, t462, pkin(4) * t462 + t384 * t257 + t381 * t258, 0, 0, t297, 0, -t447, 0, t471, t469, t450, pkin(4) * t450 + t384 * t255 + t381 * t256; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t366, t367, 0, 0, 0, 0, 0, 0, 0, t376, t316, t315, 0, -t310, 0, 0, 0, 0, 0, t372, t284, t283, 0, t260, 0, 0, 0, 0, 0, t372, t282, t283, 0, t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t386, 0, 0, -g(3), -t366, 0, 0, 0, -t359, 0, -t356, 0, t433, t344, t313, pkin(5) * t313, 0, 0, -t328, 0, -t324, 0, t459, -t290, t269, t258, 0, 0, -t328, 0, -t324, 0, t460, -t285, t265, t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t386, 0, qJDD(1), 0, g(3), 0, -t367, 0, 0, 0, t356, 0, -t359, 0, -t344, t433, t399, t309, 0, 0, t324, 0, -t328, 0, t290, t459, t440, t257, 0, 0, t324, 0, -t328, 0, t285, t460, t432, t255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t366, t367, 0, 0, 0, 0, 0, 0, 0, t376, t316, t315, 0, -t310, 0, 0, 0, 0, 0, t372, t284, t283, 0, t260, 0, 0, 0, 0, 0, t372, t282, t283, 0, t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t376, 0, -t375, 0, 0, -g(3), -t392, 0, 0, 0, -t354, 0, -t351, 0, t438, -t339, t280, pkin(6) * t280, 0, 0, -t354, 0, -t351, 0, t439, -t320, t274, t262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t375, 0, t376, 0, g(3), 0, t336, 0, 0, 0, t351, 0, -t354, 0, t339, t438, t400, t276, 0, 0, t351, 0, -t354, 0, t320, t439, t401, t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t376, t392, -t336, 0, 0, 0, 0, 0, 0, 0, t372, t395, t289, 0, -t277, 0, 0, 0, 0, 0, t372, t390, t289, 0, t402; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t372, 0, -t371, 0, 0, -g(3), t307, 0, 0, 0, t372, 0, -t371, 0, -t411, -t365, t304, qJ(4) * t304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t371, 0, t372, 0, g(3), 0, t308, 0, 0, 0, t371, 0, t372, 0, t365, -t411, t305, t294; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t372, -t307, -t308, 0, 0, 0, 0, 0, 0, 0, t372, t394, -t308, 0, -t303; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t372, 0, -t371, 0, 0, -t378, t304, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t371, 0, t372, 0, t378, 0, t305, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t372, -t304, -t305, 0, 0;];
m_new_reg  = t1;
