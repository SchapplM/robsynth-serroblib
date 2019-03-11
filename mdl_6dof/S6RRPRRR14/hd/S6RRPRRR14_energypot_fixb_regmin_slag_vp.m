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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
% StartTime: 2019-03-09 15:08:38
% EndTime: 2019-03-09 15:08:38
% DurationCPUTime: 0.26s
% Computational Cost: add. (409->74), mult. (1149->136), div. (0->0), fcn. (1517->18), ass. (0->60)
t410 = cos(pkin(6));
t419 = cos(qJ(2));
t420 = cos(qJ(1));
t426 = t420 * t419;
t414 = sin(qJ(2));
t415 = sin(qJ(1));
t429 = t415 * t414;
t399 = t410 * t426 - t429;
t405 = sin(pkin(7));
t409 = cos(pkin(7));
t406 = sin(pkin(6));
t431 = t406 * t420;
t396 = -t399 * t405 - t409 * t431;
t427 = t420 * t414;
t428 = t415 * t419;
t401 = -t410 * t428 - t427;
t433 = t406 * t415;
t397 = -t401 * t405 + t409 * t433;
t432 = t406 * t419;
t398 = -t405 * t432 + t409 * t410;
t439 = -g(1) * t397 - g(2) * t396 - g(3) * t398;
t435 = t405 * t410;
t434 = t406 * t414;
t430 = t409 * t419;
t400 = t410 * t427 + t428;
t403 = sin(pkin(14));
t407 = cos(pkin(14));
t422 = t399 * t409 - t405 * t431;
t390 = -t400 * t403 + t407 * t422;
t404 = sin(pkin(8));
t408 = cos(pkin(8));
t425 = t390 * t408 + t396 * t404;
t402 = -t410 * t429 + t426;
t421 = t401 * t409 + t405 * t433;
t392 = -t402 * t403 + t407 * t421;
t424 = t392 * t408 + t397 * t404;
t394 = t407 * t435 + (-t403 * t414 + t407 * t430) * t406;
t423 = t394 * t408 + t398 * t404;
t418 = cos(qJ(4));
t417 = cos(qJ(5));
t416 = cos(qJ(6));
t413 = sin(qJ(4));
t412 = sin(qJ(5));
t411 = sin(qJ(6));
t395 = t407 * t434 + (t406 * t430 + t435) * t403;
t393 = t402 * t407 + t403 * t421;
t391 = t400 * t407 + t403 * t422;
t389 = -t394 * t404 + t398 * t408;
t388 = -t392 * t404 + t397 * t408;
t387 = -t390 * t404 + t396 * t408;
t386 = t395 * t418 + t413 * t423;
t385 = t395 * t413 - t418 * t423;
t384 = t393 * t418 + t413 * t424;
t383 = t393 * t413 - t418 * t424;
t382 = t391 * t418 + t413 * t425;
t381 = t391 * t413 - t418 * t425;
t380 = t386 * t417 + t389 * t412;
t379 = t384 * t417 + t388 * t412;
t378 = t382 * t417 + t387 * t412;
t1 = [0, -g(1) * t420 - g(2) * t415, g(1) * t415 - g(2) * t420, 0, 0, 0, 0, 0, -g(1) * t402 - g(2) * t400 - g(3) * t434, -g(1) * t401 - g(2) * t399 - g(3) * t432, -g(1) * t393 - g(2) * t391 - g(3) * t395, -g(1) * t392 - g(2) * t390 - g(3) * t394, t439, -g(1) * (t420 * pkin(1) + t402 * pkin(2) + pkin(10) * t433) - g(2) * (t415 * pkin(1) + t400 * pkin(2) - pkin(10) * t431) - g(3) * (pkin(2) * t434 + t410 * pkin(10) + pkin(9)) + t439 * qJ(3), 0, 0, 0, 0, 0, -g(1) * t384 - g(2) * t382 - g(3) * t386, g(1) * t383 + g(2) * t381 + g(3) * t385, 0, 0, 0, 0, 0, -g(1) * t379 - g(2) * t378 - g(3) * t380, -g(1) * (-t384 * t412 + t388 * t417) - g(2) * (-t382 * t412 + t387 * t417) - g(3) * (-t386 * t412 + t389 * t417) 0, 0, 0, 0, 0, -g(1) * (t379 * t416 + t383 * t411) - g(2) * (t378 * t416 + t381 * t411) - g(3) * (t380 * t416 + t385 * t411) -g(1) * (-t379 * t411 + t383 * t416) - g(2) * (-t378 * t411 + t381 * t416) - g(3) * (-t380 * t411 + t385 * t416);];
U_reg  = t1;
