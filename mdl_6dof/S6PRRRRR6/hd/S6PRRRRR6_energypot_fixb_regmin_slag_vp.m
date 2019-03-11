% Calculate minimal parameter regressor of potential energy for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_energypot_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:25:03
% EndTime: 2019-03-09 01:25:03
% DurationCPUTime: 0.17s
% Computational Cost: add. (383->62), mult. (1087->121), div. (0->0), fcn. (1448->18), ass. (0->56)
t346 = sin(pkin(7));
t347 = sin(pkin(6));
t373 = t346 * t347;
t351 = cos(pkin(6));
t372 = t346 * t351;
t350 = cos(pkin(7));
t371 = t347 * t350;
t361 = cos(qJ(2));
t370 = t347 * t361;
t369 = t350 * t361;
t356 = sin(qJ(2));
t368 = t351 * t356;
t367 = t351 * t361;
t344 = sin(pkin(14));
t348 = cos(pkin(14));
t341 = t344 * t361 + t348 * t368;
t355 = sin(qJ(3));
t360 = cos(qJ(3));
t340 = -t344 * t356 + t348 * t367;
t363 = t340 * t350 - t348 * t373;
t331 = -t341 * t355 + t363 * t360;
t337 = -t340 * t346 - t348 * t371;
t345 = sin(pkin(8));
t349 = cos(pkin(8));
t366 = t331 * t349 + t337 * t345;
t343 = -t344 * t368 + t348 * t361;
t342 = -t344 * t367 - t348 * t356;
t362 = t342 * t350 + t344 * t373;
t333 = -t343 * t355 + t362 * t360;
t338 = -t342 * t346 + t344 * t371;
t365 = t333 * t349 + t338 * t345;
t335 = t360 * t372 + (-t355 * t356 + t360 * t369) * t347;
t339 = -t346 * t370 + t351 * t350;
t364 = t335 * t349 + t339 * t345;
t359 = cos(qJ(4));
t358 = cos(qJ(5));
t357 = cos(qJ(6));
t354 = sin(qJ(4));
t353 = sin(qJ(5));
t352 = sin(qJ(6));
t336 = t355 * t372 + (t355 * t369 + t356 * t360) * t347;
t334 = t343 * t360 + t362 * t355;
t332 = t341 * t360 + t363 * t355;
t330 = -t335 * t345 + t339 * t349;
t329 = -t333 * t345 + t338 * t349;
t328 = -t331 * t345 + t337 * t349;
t327 = t336 * t359 + t364 * t354;
t326 = t336 * t354 - t364 * t359;
t325 = t334 * t359 + t365 * t354;
t324 = t334 * t354 - t365 * t359;
t323 = t332 * t359 + t366 * t354;
t322 = t332 * t354 - t366 * t359;
t321 = t327 * t358 + t330 * t353;
t320 = t325 * t358 + t329 * t353;
t319 = t323 * t358 + t328 * t353;
t1 = [-g(3) * qJ(1), 0, -g(3) * t347 * t356 - g(1) * t343 - g(2) * t341, -g(1) * t342 - g(2) * t340 - g(3) * t370, 0, 0, 0, 0, 0, -g(1) * t334 - g(2) * t332 - g(3) * t336, -g(1) * t333 - g(2) * t331 - g(3) * t335, 0, 0, 0, 0, 0, -g(1) * t325 - g(2) * t323 - g(3) * t327, g(1) * t324 + g(2) * t322 + g(3) * t326, 0, 0, 0, 0, 0, -g(1) * t320 - g(2) * t319 - g(3) * t321, -g(1) * (-t325 * t353 + t329 * t358) - g(2) * (-t323 * t353 + t328 * t358) - g(3) * (-t327 * t353 + t330 * t358) 0, 0, 0, 0, 0, -g(1) * (t320 * t357 + t324 * t352) - g(2) * (t319 * t357 + t322 * t352) - g(3) * (t321 * t357 + t326 * t352) -g(1) * (-t320 * t352 + t324 * t357) - g(2) * (-t319 * t352 + t322 * t357) - g(3) * (-t321 * t352 + t326 * t357);];
U_reg  = t1;
