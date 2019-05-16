% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPRP1
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:14
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPRP1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP1_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:14:10
% EndTime: 2019-05-04 19:14:11
% DurationCPUTime: 0.65s
% Computational Cost: add. (1372->89), mult. (2098->84), div. (0->0), fcn. (1288->6), ass. (0->46)
t486 = (qJD(1) + qJD(3));
t484 = t486 ^ 2;
t485 = qJDD(1) + qJDD(3);
t490 = sin(qJ(3));
t492 = cos(qJ(3));
t467 = t492 * t484 + t490 * t485;
t470 = t490 * t484 - t492 * t485;
t488 = sin(pkin(6));
t489 = cos(pkin(6));
t455 = t489 * t467 - t488 * t470;
t459 = t488 * t467 + t489 * t470;
t491 = sin(qJ(1));
t493 = cos(qJ(1));
t501 = t493 * t455 - t491 * t459;
t502 = t491 * t455 + t493 * t459;
t478 = t491 * g(1) - t493 * g(2);
t472 = qJDD(1) * pkin(1) + t478;
t479 = -t493 * g(1) - t491 * g(2);
t494 = qJD(1) ^ 2;
t473 = -t494 * pkin(1) + t479;
t463 = t488 * t472 + t489 * t473;
t461 = -t494 * pkin(2) + t463;
t462 = t489 * t472 - t488 * t473;
t495 = qJDD(1) * pkin(2) + t462;
t449 = t492 * t461 + t490 * t495;
t448 = -t490 * t461 + t492 * t495;
t474 = -t488 * qJDD(1) - t489 * t494;
t475 = t489 * qJDD(1) - t488 * t494;
t497 = t493 * t474 - t491 * t475;
t496 = t491 * t474 + t493 * t475;
t487 = -g(3) + qJDD(2);
t477 = -t491 * qJDD(1) - t493 * t494;
t476 = t493 * qJDD(1) - t491 * t494;
t451 = -t488 * t462 + t489 * t463;
t450 = t489 * t462 + t488 * t463;
t445 = -t485 * pkin(3) - t484 * qJ(4) + qJDD(4) - t448;
t444 = -t484 * pkin(3) + t485 * qJ(4) + (2 * qJD(4) * t486) + t449;
t443 = -t490 * t448 + t492 * t449;
t442 = t492 * t448 + t490 * t449;
t441 = t492 * t444 + t490 * t445;
t440 = t490 * t444 - t492 * t445;
t439 = -t488 * t442 + t489 * t443;
t438 = t489 * t442 + t488 * t443;
t437 = -t488 * t440 + t489 * t441;
t436 = t489 * t440 + t488 * t441;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t477, -t476, 0, -t491 * t478 + t493 * t479, 0, 0, 0, 0, 0, 0, t497, -t496, 0, -t491 * t450 + t493 * t451, 0, 0, 0, 0, 0, 0, -t501, t502, 0, -t491 * t438 + t493 * t439, 0, 0, 0, 0, 0, 0, -t501, 0, -t502, -t491 * t436 + t493 * t437; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t476, t477, 0, t493 * t478 + t491 * t479, 0, 0, 0, 0, 0, 0, t496, t497, 0, t493 * t450 + t491 * t451, 0, 0, 0, 0, 0, 0, -t502, -t501, 0, t493 * t438 + t491 * t439, 0, 0, 0, 0, 0, 0, -t502, 0, t501, t493 * t436 + t491 * t437; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t494, -qJDD(1), 0, t479, 0, 0, 0, 0, 0, 0, t474, -t475, 0, t451, 0, 0, 0, 0, 0, 0, -t455, t459, 0, t439, 0, 0, 0, 0, 0, 0, -t455, 0, -t459, t437; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t494, 0, t478, 0, 0, 0, 0, 0, 0, t475, t474, 0, t450, 0, 0, 0, 0, 0, 0, -t459, -t455, 0, t438, 0, 0, 0, 0, 0, 0, -t459, 0, t455, t436; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t494, -qJDD(1), 0, t463, 0, 0, 0, 0, 0, 0, -t467, t470, 0, t443, 0, 0, 0, 0, 0, 0, -t467, 0, -t470, t441; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t494, 0, t462, 0, 0, 0, 0, 0, 0, -t470, -t467, 0, t442, 0, 0, 0, 0, 0, 0, -t470, 0, t467, t440; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t484, -t485, 0, t449, 0, 0, 0, 0, 0, 0, -t484, 0, t485, t444; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t485, -t484, 0, t448, 0, 0, 0, 0, 0, 0, t485, 0, t484, -t445; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t484, 0, t485, t444; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t485, 0, -t484, t445;];
f_new_reg  = t1;
