% Calculate minimal parameter regressor of potential energy for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_energypot_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:00:41
% EndTime: 2019-03-09 08:00:41
% DurationCPUTime: 0.20s
% Computational Cost: add. (391->69), mult. (1105->132), div. (0->0), fcn. (1464->18), ass. (0->59)
t403 = sin(pkin(7));
t408 = cos(pkin(6));
t433 = t403 * t408;
t404 = sin(pkin(6));
t405 = cos(pkin(14));
t432 = t404 * t405;
t413 = sin(qJ(1));
t431 = t404 * t413;
t418 = cos(qJ(1));
t430 = t404 * t418;
t407 = cos(pkin(7));
t429 = t405 * t407;
t401 = sin(pkin(14));
t428 = t413 * t401;
t427 = t413 * t405;
t426 = t418 * t401;
t425 = t418 * t405;
t424 = g(1) * t413 - g(2) * t418;
t398 = t408 * t426 + t427;
t412 = sin(qJ(3));
t417 = cos(qJ(3));
t397 = t408 * t425 - t428;
t420 = t397 * t407 - t403 * t430;
t388 = -t398 * t412 + t420 * t417;
t394 = -t397 * t403 - t407 * t430;
t402 = sin(pkin(8));
t406 = cos(pkin(8));
t423 = t388 * t406 + t394 * t402;
t400 = -t408 * t428 + t425;
t399 = -t408 * t427 - t426;
t419 = t399 * t407 + t403 * t431;
t390 = -t400 * t412 + t419 * t417;
t395 = -t399 * t403 + t407 * t431;
t422 = t390 * t406 + t395 * t402;
t392 = t417 * t433 + (-t401 * t412 + t417 * t429) * t404;
t396 = -t403 * t432 + t408 * t407;
t421 = t392 * t406 + t396 * t402;
t416 = cos(qJ(4));
t415 = cos(qJ(5));
t414 = cos(qJ(6));
t411 = sin(qJ(4));
t410 = sin(qJ(5));
t409 = sin(qJ(6));
t393 = t412 * t433 + (t401 * t417 + t412 * t429) * t404;
t391 = t400 * t417 + t419 * t412;
t389 = t398 * t417 + t420 * t412;
t387 = -t392 * t402 + t396 * t406;
t386 = -t390 * t402 + t395 * t406;
t385 = -t388 * t402 + t394 * t406;
t384 = t393 * t416 + t421 * t411;
t383 = t393 * t411 - t421 * t416;
t382 = t391 * t416 + t422 * t411;
t381 = t391 * t411 - t422 * t416;
t380 = t389 * t416 + t423 * t411;
t379 = t389 * t411 - t423 * t416;
t378 = t384 * t415 + t387 * t410;
t377 = t382 * t415 + t386 * t410;
t376 = t380 * t415 + t385 * t410;
t1 = [0, -g(1) * t418 - g(2) * t413, t424, -g(3) * t404 * t401 - g(1) * t400 - g(2) * t398, -g(1) * t399 - g(2) * t397 - g(3) * t432, -g(3) * t408 - t424 * t404, -g(1) * (t418 * pkin(1) + qJ(2) * t431) - g(2) * (t413 * pkin(1) - qJ(2) * t430) - g(3) * (t408 * qJ(2) + pkin(9)) 0, 0, 0, 0, 0, -g(1) * t391 - g(2) * t389 - g(3) * t393, -g(1) * t390 - g(2) * t388 - g(3) * t392, 0, 0, 0, 0, 0, -g(1) * t382 - g(2) * t380 - g(3) * t384, g(1) * t381 + g(2) * t379 + g(3) * t383, 0, 0, 0, 0, 0, -g(1) * t377 - g(2) * t376 - g(3) * t378, -g(1) * (-t382 * t410 + t386 * t415) - g(2) * (-t380 * t410 + t385 * t415) - g(3) * (-t384 * t410 + t387 * t415) 0, 0, 0, 0, 0, -g(1) * (t377 * t414 + t381 * t409) - g(2) * (t376 * t414 + t379 * t409) - g(3) * (t378 * t414 + t383 * t409) -g(1) * (-t377 * t409 + t381 * t414) - g(2) * (-t376 * t409 + t379 * t414) - g(3) * (-t378 * t409 + t383 * t414);];
U_reg  = t1;
