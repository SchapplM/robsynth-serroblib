% Calculate minimal parameter regressor of potential energy for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR14_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energypot_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:32:32
% EndTime: 2018-12-10 18:32:32
% DurationCPUTime: 0.35s
% Computational Cost: add. (1498->97), mult. (1512->149), div. (0->0), fcn. (1517->30), ass. (0->77)
t450 = pkin(6) - qJ(2);
t438 = cos(t450) / 0.2e1;
t449 = pkin(6) + qJ(2);
t444 = cos(t449);
t431 = t438 + t444 / 0.2e1;
t462 = sin(qJ(2));
t463 = sin(qJ(1));
t468 = cos(qJ(1));
t417 = t468 * t431 - t463 * t462;
t453 = sin(pkin(7));
t457 = cos(pkin(7));
t454 = sin(pkin(6));
t469 = t454 * t468;
t414 = -t417 * t453 - t457 * t469;
t419 = -t463 * t431 - t468 * t462;
t470 = t454 * t463;
t415 = -t419 * t453 + t457 * t470;
t436 = sin(t449) / 0.2e1;
t442 = sin(t450);
t427 = t436 + t442 / 0.2e1;
t458 = cos(pkin(6));
t416 = -t427 * t453 + t458 * t457;
t474 = -g(1) * t415 - g(2) * t414 - g(3) * t416;
t467 = cos(qJ(2));
t466 = cos(qJ(4));
t465 = cos(qJ(5));
t464 = cos(qJ(6));
t461 = sin(qJ(4));
t460 = sin(qJ(5));
t459 = sin(qJ(6));
t456 = cos(pkin(8));
t455 = cos(pkin(14));
t452 = sin(pkin(8));
t451 = sin(pkin(14));
t448 = pkin(8) - qJ(4);
t447 = pkin(8) + qJ(4);
t446 = pkin(7) - pkin(14);
t445 = pkin(7) + pkin(14);
t443 = cos(t447);
t441 = sin(t448);
t440 = cos(t445);
t439 = sin(t446);
t437 = cos(t448) / 0.2e1;
t435 = sin(t447) / 0.2e1;
t434 = cos(t446) / 0.2e1;
t433 = sin(t445) / 0.2e1;
t432 = t438 - t444 / 0.2e1;
t430 = t437 - t443 / 0.2e1;
t429 = t437 + t443 / 0.2e1;
t428 = t436 - t442 / 0.2e1;
t426 = t435 - t441 / 0.2e1;
t425 = t435 + t441 / 0.2e1;
t424 = t434 - t440 / 0.2e1;
t423 = t434 + t440 / 0.2e1;
t422 = t433 - t439 / 0.2e1;
t421 = t433 + t439 / 0.2e1;
t420 = -t463 * t428 + t468 * t467;
t418 = t468 * t428 + t463 * t467;
t413 = t427 * t422 + t458 * t424 + t432 * t455;
t412 = t458 * t421 + t427 * t423 - t432 * t451;
t411 = t419 * t422 + t420 * t455 + t424 * t470;
t410 = t419 * t423 - t420 * t451 + t421 * t470;
t409 = t417 * t422 + t418 * t455 - t424 * t469;
t408 = t417 * t423 - t418 * t451 - t421 * t469;
t407 = -t412 * t452 + t416 * t456;
t406 = -t410 * t452 + t415 * t456;
t405 = -t408 * t452 + t414 * t456;
t404 = t412 * t426 + t413 * t466 + t416 * t430;
t403 = -t412 * t429 + t413 * t461 - t416 * t425;
t402 = t410 * t426 + t411 * t466 + t415 * t430;
t401 = -t410 * t429 + t411 * t461 - t415 * t425;
t400 = t408 * t426 + t409 * t466 + t414 * t430;
t399 = -t408 * t429 + t409 * t461 - t414 * t425;
t398 = t404 * t465 + t407 * t460;
t397 = t402 * t465 + t406 * t460;
t396 = t400 * t465 + t405 * t460;
t1 = [0, -g(1) * t468 - g(2) * t463, g(1) * t463 - g(2) * t468, 0, 0, 0, 0, 0, -g(1) * t420 - g(2) * t418 - g(3) * t432, -g(1) * t419 - g(2) * t417 - g(3) * t427, -g(1) * t411 - g(2) * t409 - g(3) * t413, -g(1) * t410 - g(2) * t408 - g(3) * t412, t474, -g(1) * (t468 * pkin(1) + t420 * pkin(2) + pkin(10) * t470) - g(2) * (t463 * pkin(1) + t418 * pkin(2) - pkin(10) * t469) - g(3) * (t432 * pkin(2) + t458 * pkin(10) + pkin(9)) + t474 * qJ(3), 0, 0, 0, 0, 0, -g(1) * t402 - g(2) * t400 - g(3) * t404, g(1) * t401 + g(2) * t399 + g(3) * t403, 0, 0, 0, 0, 0, -g(1) * t397 - g(2) * t396 - g(3) * t398, -g(1) * (-t402 * t460 + t406 * t465) - g(2) * (-t400 * t460 + t405 * t465) - g(3) * (-t404 * t460 + t407 * t465) 0, 0, 0, 0, 0, -g(1) * (t397 * t464 + t401 * t459) - g(2) * (t396 * t464 + t399 * t459) - g(3) * (t398 * t464 + t403 * t459) -g(1) * (-t397 * t459 + t401 * t464) - g(2) * (-t396 * t459 + t399 * t464) - g(3) * (-t398 * t459 + t403 * t464);];
U_reg  = t1;
