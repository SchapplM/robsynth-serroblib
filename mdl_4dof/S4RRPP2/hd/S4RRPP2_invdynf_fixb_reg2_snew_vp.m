% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRPP2
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
% Datum: 2019-05-04 19:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRPP2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP2_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:20:33
% EndTime: 2019-05-04 19:20:34
% DurationCPUTime: 0.51s
% Computational Cost: add. (673->84), mult. (817->60), div. (0->0), fcn. (456->4), ass. (0->33)
t453 = sin(qJ(1));
t455 = cos(qJ(1));
t448 = qJDD(1) + qJDD(2);
t452 = sin(qJ(2));
t454 = cos(qJ(2));
t449 = (qJD(1) + qJD(2));
t464 = t449 ^ 2;
t465 = t452 * t448 + t454 * t464;
t466 = t454 * t448 - t452 * t464;
t470 = t453 * t466 + t455 * t465;
t469 = -t453 * t465 + t455 * t466;
t439 = -t455 * g(1) - t453 * g(2);
t456 = qJD(1) ^ 2;
t435 = -t456 * pkin(1) + t439;
t438 = t453 * g(1) - t455 * g(2);
t457 = qJDD(1) * pkin(1) + t438;
t423 = t454 * t435 + t452 * t457;
t422 = -t452 * t435 + t454 * t457;
t460 = t448 * qJ(3) + (2 * qJD(3) * t449) + t423;
t419 = -t448 * pkin(2) - qJ(3) * t464 + qJDD(3) - t422;
t451 = g(3) + qJDD(4);
t437 = -t453 * qJDD(1) - t455 * t456;
t436 = t455 * qJDD(1) - t453 * t456;
t418 = -pkin(2) * t464 + t460;
t417 = -t448 * pkin(3) + t419;
t416 = (-pkin(2) - pkin(3)) * t464 + t460;
t415 = -t452 * t422 + t454 * t423;
t414 = t454 * t422 + t452 * t423;
t413 = t454 * t418 + t452 * t419;
t412 = t452 * t418 - t454 * t419;
t411 = t454 * t416 + t452 * t417;
t410 = t452 * t416 - t454 * t417;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t437, -t436, 0, -t453 * t438 + t455 * t439, 0, 0, 0, 0, 0, 0, -t470, -t469, 0, -t453 * t414 + t455 * t415, 0, 0, 0, 0, 0, 0, -t470, 0, t469, -t453 * t412 + t455 * t413, 0, 0, 0, 0, 0, 0, -t470, t469, 0, -t453 * t410 + t455 * t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t436, t437, 0, t455 * t438 + t453 * t439, 0, 0, 0, 0, 0, 0, t469, -t470, 0, t455 * t414 + t453 * t415, 0, 0, 0, 0, 0, 0, t469, 0, t470, t455 * t412 + t453 * t413, 0, 0, 0, 0, 0, 0, t469, t470, 0, t455 * t410 + t453 * t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t451; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t456, -qJDD(1), 0, t439, 0, 0, 0, 0, 0, 0, -t465, -t466, 0, t415, 0, 0, 0, 0, 0, 0, -t465, 0, t466, t413, 0, 0, 0, 0, 0, 0, -t465, t466, 0, t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t456, 0, t438, 0, 0, 0, 0, 0, 0, t466, -t465, 0, t414, 0, 0, 0, 0, 0, 0, t466, 0, t465, t412, 0, 0, 0, 0, 0, 0, t466, t465, 0, t410; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t451; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t464, -t448, 0, t423, 0, 0, 0, 0, 0, 0, -t464, 0, t448, t418, 0, 0, 0, 0, 0, 0, -t464, t448, 0, t416; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t448, -t464, 0, t422, 0, 0, 0, 0, 0, 0, t448, 0, t464, -t419, 0, 0, 0, 0, 0, 0, t448, t464, 0, -t417; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t451; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t464, 0, t448, t418, 0, 0, 0, 0, 0, 0, -t464, t448, 0, t416; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t451; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t448, 0, -t464, t419, 0, 0, 0, 0, 0, 0, -t448, -t464, 0, t417; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t464, t448, 0, t416; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t448, -t464, 0, t417; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t451;];
f_new_reg  = t1;
