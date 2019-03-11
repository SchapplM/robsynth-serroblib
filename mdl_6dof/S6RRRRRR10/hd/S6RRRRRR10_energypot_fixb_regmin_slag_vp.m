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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
% StartTime: 2019-03-10 06:16:48
% EndTime: 2019-03-10 06:16:48
% DurationCPUTime: 0.18s
% Computational Cost: add. (384->63), mult. (1090->122), div. (0->0), fcn. (1452->18), ass. (0->58)
t397 = sin(pkin(7));
t401 = cos(pkin(6));
t427 = t397 * t401;
t398 = sin(pkin(6));
t407 = sin(qJ(1));
t426 = t398 * t407;
t412 = cos(qJ(2));
t425 = t398 * t412;
t413 = cos(qJ(1));
t424 = t398 * t413;
t400 = cos(pkin(7));
t423 = t400 * t412;
t406 = sin(qJ(2));
t422 = t407 * t406;
t421 = t407 * t412;
t420 = t413 * t406;
t419 = t413 * t412;
t393 = t401 * t420 + t421;
t405 = sin(qJ(3));
t411 = cos(qJ(3));
t392 = t401 * t419 - t422;
t415 = t392 * t400 - t397 * t424;
t383 = -t393 * t405 + t415 * t411;
t389 = -t392 * t397 - t400 * t424;
t396 = sin(pkin(8));
t399 = cos(pkin(8));
t418 = t383 * t399 + t389 * t396;
t395 = -t401 * t422 + t419;
t394 = -t401 * t421 - t420;
t414 = t394 * t400 + t397 * t426;
t385 = -t395 * t405 + t414 * t411;
t390 = -t394 * t397 + t400 * t426;
t417 = t385 * t399 + t390 * t396;
t387 = t411 * t427 + (-t405 * t406 + t411 * t423) * t398;
t391 = -t397 * t425 + t401 * t400;
t416 = t387 * t399 + t391 * t396;
t410 = cos(qJ(4));
t409 = cos(qJ(5));
t408 = cos(qJ(6));
t404 = sin(qJ(4));
t403 = sin(qJ(5));
t402 = sin(qJ(6));
t388 = t405 * t427 + (t405 * t423 + t406 * t411) * t398;
t386 = t395 * t411 + t414 * t405;
t384 = t393 * t411 + t415 * t405;
t382 = -t387 * t396 + t391 * t399;
t381 = -t385 * t396 + t390 * t399;
t380 = -t383 * t396 + t389 * t399;
t379 = t388 * t410 + t416 * t404;
t378 = t388 * t404 - t416 * t410;
t377 = t386 * t410 + t417 * t404;
t376 = t386 * t404 - t417 * t410;
t375 = t384 * t410 + t418 * t404;
t374 = t384 * t404 - t418 * t410;
t373 = t379 * t409 + t382 * t403;
t372 = t377 * t409 + t381 * t403;
t371 = t375 * t409 + t380 * t403;
t1 = [0, -g(1) * t413 - g(2) * t407, g(1) * t407 - g(2) * t413, 0, 0, 0, 0, 0, -g(3) * t398 * t406 - g(1) * t395 - g(2) * t393, -g(1) * t394 - g(2) * t392 - g(3) * t425, 0, 0, 0, 0, 0, -g(1) * t386 - g(2) * t384 - g(3) * t388, -g(1) * t385 - g(2) * t383 - g(3) * t387, 0, 0, 0, 0, 0, -g(1) * t377 - g(2) * t375 - g(3) * t379, g(1) * t376 + g(2) * t374 + g(3) * t378, 0, 0, 0, 0, 0, -g(1) * t372 - g(2) * t371 - g(3) * t373, -g(1) * (-t377 * t403 + t381 * t409) - g(2) * (-t375 * t403 + t380 * t409) - g(3) * (-t379 * t403 + t382 * t409) 0, 0, 0, 0, 0, -g(1) * (t372 * t408 + t376 * t402) - g(2) * (t371 * t408 + t374 * t402) - g(3) * (t373 * t408 + t378 * t402) -g(1) * (-t372 * t402 + t376 * t408) - g(2) * (-t371 * t402 + t374 * t408) - g(3) * (-t373 * t402 + t378 * t408);];
U_reg  = t1;
