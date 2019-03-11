% Calculate minimal parameter regressor of potential energy for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_energypot_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:55:53
% EndTime: 2019-03-08 20:55:54
% DurationCPUTime: 0.25s
% Computational Cost: add. (408->73), mult. (1146->136), div. (0->0), fcn. (1513->18), ass. (0->59)
t344 = sin(pkin(13));
t349 = cos(pkin(13));
t356 = sin(qJ(2));
t352 = cos(pkin(6));
t360 = cos(qJ(2));
t366 = t352 * t360;
t339 = -t344 * t356 + t349 * t366;
t346 = sin(pkin(7));
t347 = sin(pkin(6));
t351 = cos(pkin(7));
t371 = t347 * t351;
t336 = -t339 * t346 - t349 * t371;
t341 = -t344 * t366 - t349 * t356;
t337 = -t341 * t346 + t344 * t371;
t369 = t347 * t360;
t338 = -t346 * t369 + t352 * t351;
t378 = -g(1) * t337 - g(2) * t336 - g(3) * t338;
t374 = t344 * t347;
t373 = t346 * t352;
t372 = t347 * t349;
t370 = t347 * t356;
t368 = t351 * t360;
t367 = t352 * t356;
t340 = t344 * t360 + t349 * t367;
t343 = sin(pkin(14));
t348 = cos(pkin(14));
t362 = t339 * t351 - t346 * t372;
t330 = -t340 * t343 + t362 * t348;
t345 = sin(pkin(8));
t350 = cos(pkin(8));
t365 = t330 * t350 + t336 * t345;
t342 = -t344 * t367 + t349 * t360;
t361 = t341 * t351 + t346 * t374;
t332 = -t342 * t343 + t361 * t348;
t364 = t332 * t350 + t337 * t345;
t334 = t348 * t373 + (-t343 * t356 + t348 * t368) * t347;
t363 = t334 * t350 + t338 * t345;
t359 = cos(qJ(4));
t358 = cos(qJ(5));
t357 = cos(qJ(6));
t355 = sin(qJ(4));
t354 = sin(qJ(5));
t353 = sin(qJ(6));
t335 = t348 * t370 + (t347 * t368 + t373) * t343;
t333 = t342 * t348 + t361 * t343;
t331 = t340 * t348 + t362 * t343;
t329 = -t334 * t345 + t338 * t350;
t328 = -t332 * t345 + t337 * t350;
t327 = -t330 * t345 + t336 * t350;
t326 = t335 * t359 + t363 * t355;
t325 = t335 * t355 - t363 * t359;
t324 = t333 * t359 + t364 * t355;
t323 = t333 * t355 - t364 * t359;
t322 = t331 * t359 + t365 * t355;
t321 = t331 * t355 - t365 * t359;
t320 = t326 * t358 + t329 * t354;
t319 = t324 * t358 + t328 * t354;
t318 = t322 * t358 + t327 * t354;
t1 = [-g(3) * qJ(1), 0, -g(1) * t342 - g(2) * t340 - g(3) * t370, -g(1) * t341 - g(2) * t339 - g(3) * t369, -g(1) * t333 - g(2) * t331 - g(3) * t335, -g(1) * t332 - g(2) * t330 - g(3) * t334, t378, -g(1) * (t349 * pkin(1) + t342 * pkin(2) + pkin(9) * t374) - g(2) * (t344 * pkin(1) + t340 * pkin(2) - pkin(9) * t372) - g(3) * (pkin(2) * t370 + t352 * pkin(9) + qJ(1)) + t378 * qJ(3), 0, 0, 0, 0, 0, -g(1) * t324 - g(2) * t322 - g(3) * t326, g(1) * t323 + g(2) * t321 + g(3) * t325, 0, 0, 0, 0, 0, -g(1) * t319 - g(2) * t318 - g(3) * t320, -g(1) * (-t324 * t354 + t328 * t358) - g(2) * (-t322 * t354 + t327 * t358) - g(3) * (-t326 * t354 + t329 * t358) 0, 0, 0, 0, 0, -g(1) * (t319 * t357 + t323 * t353) - g(2) * (t318 * t357 + t321 * t353) - g(3) * (t320 * t357 + t325 * t353) -g(1) * (-t319 * t353 + t323 * t357) - g(2) * (-t318 * t353 + t321 * t357) - g(3) * (-t320 * t353 + t325 * t357);];
U_reg  = t1;
