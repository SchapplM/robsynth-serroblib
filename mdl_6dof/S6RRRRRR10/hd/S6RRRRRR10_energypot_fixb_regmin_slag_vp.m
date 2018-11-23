% Calculate minimal parameter regressor of potential energy for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x38]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_energypot_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:16:23
% EndTime: 2018-11-23 11:16:23
% DurationCPUTime: 0.25s
% Computational Cost: add. (1446->86), mult. (1444->134), div. (0->0), fcn. (1452->30), ass. (0->76)
t446 = sin(pkin(6));
t455 = sin(qJ(1));
t463 = t446 * t455;
t461 = cos(qJ(1));
t462 = t446 * t461;
t460 = cos(qJ(2));
t459 = cos(qJ(3));
t458 = cos(qJ(4));
t457 = cos(qJ(5));
t456 = cos(qJ(6));
t454 = sin(qJ(2));
t453 = sin(qJ(3));
t452 = sin(qJ(4));
t451 = sin(qJ(5));
t450 = sin(qJ(6));
t449 = cos(pkin(6));
t448 = cos(pkin(7));
t447 = cos(pkin(8));
t445 = sin(pkin(7));
t444 = sin(pkin(8));
t443 = pkin(6) - qJ(2);
t442 = pkin(6) + qJ(2);
t441 = pkin(7) - qJ(3);
t440 = pkin(7) + qJ(3);
t439 = pkin(8) - qJ(4);
t438 = pkin(8) + qJ(4);
t437 = cos(t442);
t436 = cos(t440);
t435 = cos(t438);
t434 = sin(t443);
t433 = sin(t441);
t432 = sin(t439);
t431 = cos(t443) / 0.2e1;
t430 = cos(t441) / 0.2e1;
t429 = cos(t439) / 0.2e1;
t428 = sin(t442) / 0.2e1;
t427 = sin(t440) / 0.2e1;
t426 = sin(t438) / 0.2e1;
t425 = t431 - t437 / 0.2e1;
t424 = t431 + t437 / 0.2e1;
t423 = t430 - t436 / 0.2e1;
t422 = t430 + t436 / 0.2e1;
t421 = t429 - t435 / 0.2e1;
t420 = t429 + t435 / 0.2e1;
t419 = t428 - t434 / 0.2e1;
t418 = t428 + t434 / 0.2e1;
t417 = t427 - t433 / 0.2e1;
t416 = t427 + t433 / 0.2e1;
t415 = t426 - t432 / 0.2e1;
t414 = t426 + t432 / 0.2e1;
t413 = -t455 * t419 + t461 * t460;
t412 = -t455 * t424 - t461 * t454;
t411 = t461 * t419 + t455 * t460;
t410 = t461 * t424 - t455 * t454;
t409 = -t418 * t445 + t449 * t448;
t408 = -t412 * t445 + t448 * t463;
t407 = -t410 * t445 - t448 * t462;
t406 = t418 * t417 + t449 * t423 + t425 * t459;
t405 = t449 * t416 + t418 * t422 - t425 * t453;
t404 = t412 * t417 + t413 * t459 + t423 * t463;
t403 = t412 * t422 - t413 * t453 + t416 * t463;
t402 = t410 * t417 + t411 * t459 - t423 * t462;
t401 = t410 * t422 - t411 * t453 - t416 * t462;
t400 = -t405 * t444 + t409 * t447;
t399 = -t403 * t444 + t408 * t447;
t398 = -t401 * t444 + t407 * t447;
t397 = t405 * t415 + t406 * t458 + t409 * t421;
t396 = -t405 * t420 + t406 * t452 - t409 * t414;
t395 = t403 * t415 + t404 * t458 + t408 * t421;
t394 = -t403 * t420 + t404 * t452 - t408 * t414;
t393 = t401 * t415 + t402 * t458 + t407 * t421;
t392 = -t401 * t420 + t402 * t452 - t407 * t414;
t391 = t397 * t457 + t400 * t451;
t390 = t395 * t457 + t399 * t451;
t389 = t393 * t457 + t398 * t451;
t1 = [0, -g(1) * t461 - g(2) * t455, g(1) * t455 - g(2) * t461, 0, 0, 0, 0, 0, -g(1) * t413 - g(2) * t411 - g(3) * t425, -g(1) * t412 - g(2) * t410 - g(3) * t418, 0, 0, 0, 0, 0, -g(1) * t404 - g(2) * t402 - g(3) * t406, -g(1) * t403 - g(2) * t401 - g(3) * t405, 0, 0, 0, 0, 0, -g(1) * t395 - g(2) * t393 - g(3) * t397, g(1) * t394 + g(2) * t392 + g(3) * t396, 0, 0, 0, 0, 0, -g(1) * t390 - g(2) * t389 - g(3) * t391, -g(1) * (-t395 * t451 + t399 * t457) - g(2) * (-t393 * t451 + t398 * t457) - g(3) * (-t397 * t451 + t400 * t457) 0, 0, 0, 0, 0, -g(1) * (t390 * t456 + t394 * t450) - g(2) * (t389 * t456 + t392 * t450) - g(3) * (t391 * t456 + t396 * t450) -g(1) * (-t390 * t450 + t394 * t456) - g(2) * (-t389 * t450 + t392 * t456) - g(3) * (-t391 * t450 + t396 * t456);];
U_reg  = t1;
