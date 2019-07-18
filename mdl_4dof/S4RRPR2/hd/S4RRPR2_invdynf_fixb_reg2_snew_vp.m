% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRPR2
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
%   pkin=[a2,a3,a4,d1,d2]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRPR2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR2_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:40
% EndTime: 2019-07-18 18:16:40
% DurationCPUTime: 0.60s
% Computational Cost: add. (1476->97), mult. (1722->81), div. (0->0), fcn. (920->6), ass. (0->47)
t480 = sin(qJ(2));
t483 = cos(qJ(2));
t478 = (qJD(1) + qJD(2));
t474 = qJD(4) - t478;
t470 = t474 ^ 2;
t477 = qJDD(1) + qJDD(2);
t471 = -qJDD(4) + t477;
t479 = sin(qJ(4));
t482 = cos(qJ(4));
t487 = t479 * t470 + t482 * t471;
t490 = -t482 * t470 + t479 * t471;
t442 = t480 * t487 - t483 * t490;
t481 = sin(qJ(1));
t484 = cos(qJ(1));
t491 = t480 * t490 + t483 * t487;
t496 = -t481 * t442 + t484 * t491;
t495 = t484 * t442 + t481 * t491;
t476 = t478 ^ 2;
t458 = t483 * t476 + t480 * t477;
t461 = t480 * t476 - t483 * t477;
t488 = t484 * t458 - t481 * t461;
t492 = t481 * t458 + t484 * t461;
t466 = -t484 * g(1) - t481 * g(2);
t485 = qJD(1) ^ 2;
t462 = -t485 * pkin(1) + t466;
t465 = t481 * g(1) - t484 * g(2);
t486 = qJDD(1) * pkin(1) + t465;
t449 = t483 * t462 + t480 * t486;
t448 = -t480 * t462 + t483 * t486;
t489 = t477 * qJ(3) + (2 * qJD(3) * t478) + t449;
t441 = -t477 * pkin(2) - t476 * qJ(3) + qJDD(3) - t448;
t464 = -t481 * qJDD(1) - t484 * t485;
t463 = t484 * qJDD(1) - t481 * t485;
t440 = -t476 * pkin(2) + t489;
t439 = -t477 * pkin(3) + t441;
t438 = (-pkin(2) - pkin(3)) * t476 + t489;
t437 = -t480 * t448 + t483 * t449;
t436 = t483 * t448 + t480 * t449;
t435 = t483 * t440 + t480 * t441;
t434 = t480 * t440 - t483 * t441;
t433 = t482 * t438 + t479 * t439;
t432 = -t479 * t438 + t482 * t439;
t431 = -t479 * t432 + t482 * t433;
t430 = t482 * t432 + t479 * t433;
t429 = t480 * t430 + t483 * t431;
t428 = -t483 * t430 + t480 * t431;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t464, -t463, 0, -t481 * t465 + t484 * t466, 0, 0, 0, 0, 0, 0, -t488, t492, 0, -t481 * t436 + t484 * t437, 0, 0, 0, 0, 0, 0, -t488, 0, -t492, -t481 * t434 + t484 * t435, 0, 0, 0, 0, 0, 0, -t495, t496, 0, -t481 * t428 + t484 * t429; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t463, t464, 0, t484 * t465 + t481 * t466, 0, 0, 0, 0, 0, 0, -t492, -t488, 0, t484 * t436 + t481 * t437, 0, 0, 0, 0, 0, 0, -t492, 0, t488, t484 * t434 + t481 * t435, 0, 0, 0, 0, 0, 0, t496, t495, 0, t484 * t428 + t481 * t429; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t485, -qJDD(1), 0, t466, 0, 0, 0, 0, 0, 0, -t458, t461, 0, t437, 0, 0, 0, 0, 0, 0, -t458, 0, -t461, t435, 0, 0, 0, 0, 0, 0, -t442, t491, 0, t429; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t485, 0, t465, 0, 0, 0, 0, 0, 0, -t461, -t458, 0, t436, 0, 0, 0, 0, 0, 0, -t461, 0, t458, t434, 0, 0, 0, 0, 0, 0, t491, t442, 0, t428; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t476, -t477, 0, t449, 0, 0, 0, 0, 0, 0, -t476, 0, t477, t440, 0, 0, 0, 0, 0, 0, t490, t487, 0, t431; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t477, -t476, 0, t448, 0, 0, 0, 0, 0, 0, t477, 0, t476, -t441, 0, 0, 0, 0, 0, 0, t487, -t490, 0, -t430; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t476, 0, t477, t440, 0, 0, 0, 0, 0, 0, t490, t487, 0, t431; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t477, 0, -t476, t441, 0, 0, 0, 0, 0, 0, -t487, t490, 0, t430; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t470, t471, 0, t433; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t471, -t470, 0, t432; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3);];
f_new_reg  = t1;
