% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RRPP1
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RRPP1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP1_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:19:18
% EndTime: 2019-05-04 19:19:18
% DurationCPUTime: 0.69s
% Computational Cost: add. (1523->92), mult. (2098->84), div. (0->0), fcn. (1288->6), ass. (0->46)
t486 = (qJD(1) + qJD(2));
t484 = t486 ^ 2;
t485 = qJDD(1) + qJDD(2);
t488 = sin(pkin(6));
t489 = cos(pkin(6));
t466 = t484 * t489 + t485 * t488;
t469 = t484 * t488 - t485 * t489;
t490 = sin(qJ(2));
t492 = cos(qJ(2));
t455 = t466 * t492 - t469 * t490;
t459 = t466 * t490 + t469 * t492;
t491 = sin(qJ(1));
t493 = cos(qJ(1));
t504 = t455 * t493 - t459 * t491;
t505 = t455 * t491 + t459 * t493;
t473 = t484 * t490 - t485 * t492;
t496 = -t484 * t492 - t485 * t490;
t503 = t473 * t491 + t493 * t496;
t502 = t473 * t493 - t491 * t496;
t479 = t491 * g(1) - g(2) * t493;
t475 = qJDD(1) * pkin(1) + t479;
t480 = -g(1) * t493 - g(2) * t491;
t494 = qJD(1) ^ 2;
t476 = -pkin(1) * t494 + t480;
t463 = t490 * t475 + t492 * t476;
t461 = -pkin(2) * t484 + t463;
t462 = t475 * t492 - t476 * t490;
t495 = pkin(2) * t485 + t462;
t449 = t489 * t461 + t488 * t495;
t448 = -t461 * t488 + t489 * t495;
t487 = -g(3) + qJDD(3);
t478 = -qJDD(1) * t491 - t493 * t494;
t477 = qJDD(1) * t493 - t491 * t494;
t451 = -t462 * t490 + t463 * t492;
t450 = t462 * t492 + t463 * t490;
t445 = -t485 * pkin(3) - t484 * qJ(4) + qJDD(4) - t448;
t444 = -pkin(3) * t484 + qJ(4) * t485 + (2 * qJD(4) * t486) + t449;
t443 = -t448 * t488 + t449 * t489;
t442 = t448 * t489 + t449 * t488;
t441 = t444 * t489 + t445 * t488;
t440 = t444 * t488 - t445 * t489;
t439 = -t442 * t490 + t443 * t492;
t438 = t442 * t492 + t443 * t490;
t437 = -t440 * t490 + t441 * t492;
t436 = t440 * t492 + t441 * t490;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t478, -t477, 0, -t479 * t491 + t480 * t493, 0, 0, 0, 0, 0, 0, t503, t502, 0, -t450 * t491 + t451 * t493, 0, 0, 0, 0, 0, 0, -t504, t505, 0, -t438 * t491 + t439 * t493, 0, 0, 0, 0, 0, 0, -t504, 0, -t505, -t436 * t491 + t437 * t493; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t477, t478, 0, t479 * t493 + t480 * t491, 0, 0, 0, 0, 0, 0, -t502, t503, 0, t450 * t493 + t451 * t491, 0, 0, 0, 0, 0, 0, -t505, -t504, 0, t438 * t493 + t439 * t491, 0, 0, 0, 0, 0, 0, -t505, 0, t504, t436 * t493 + t437 * t491; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t494, -qJDD(1), 0, t480, 0, 0, 0, 0, 0, 0, t496, t473, 0, t451, 0, 0, 0, 0, 0, 0, -t455, t459, 0, t439, 0, 0, 0, 0, 0, 0, -t455, 0, -t459, t437; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t494, 0, t479, 0, 0, 0, 0, 0, 0, -t473, t496, 0, t450, 0, 0, 0, 0, 0, 0, -t459, -t455, 0, t438, 0, 0, 0, 0, 0, 0, -t459, 0, t455, t436; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t484, -t485, 0, t463, 0, 0, 0, 0, 0, 0, -t466, t469, 0, t443, 0, 0, 0, 0, 0, 0, -t466, 0, -t469, t441; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t485, -t484, 0, t462, 0, 0, 0, 0, 0, 0, -t469, -t466, 0, t442, 0, 0, 0, 0, 0, 0, -t469, 0, t466, t440; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t484, -t485, 0, t449, 0, 0, 0, 0, 0, 0, -t484, 0, t485, t444; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t485, -t484, 0, t448, 0, 0, 0, 0, 0, 0, t485, 0, t484, -t445; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t484, 0, t485, t444; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t485, 0, -t484, t445;];
f_new_reg  = t1;
