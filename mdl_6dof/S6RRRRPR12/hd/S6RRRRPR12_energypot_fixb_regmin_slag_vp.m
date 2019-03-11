% Calculate minimal parameter regressor of potential energy for
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:44:39
% EndTime: 2019-03-09 23:44:39
% DurationCPUTime: 0.20s
% Computational Cost: add. (222->76), mult. (553->125), div. (0->0), fcn. (711->16), ass. (0->53)
t343 = cos(pkin(6));
t353 = cos(qJ(2));
t354 = cos(qJ(1));
t358 = t354 * t353;
t348 = sin(qJ(2));
t349 = sin(qJ(1));
t361 = t349 * t348;
t329 = t343 * t358 - t361;
t340 = sin(pkin(7));
t342 = cos(pkin(7));
t341 = sin(pkin(6));
t364 = t341 * t354;
t324 = -t329 * t340 - t342 * t364;
t346 = sin(qJ(4));
t371 = t324 * t346;
t359 = t354 * t348;
t360 = t349 * t353;
t331 = -t343 * t360 - t359;
t366 = t341 * t349;
t325 = -t331 * t340 + t342 * t366;
t370 = t325 * t346;
t365 = t341 * t353;
t328 = -t340 * t365 + t343 * t342;
t369 = t328 * t346;
t368 = t340 * t343;
t367 = t341 * t348;
t352 = cos(qJ(3));
t363 = t342 * t352;
t362 = t342 * t353;
t357 = t340 * t366;
t356 = t340 * t364;
t330 = t343 * t359 + t360;
t347 = sin(qJ(3));
t318 = -t329 * t363 + t330 * t347 + t352 * t356;
t332 = -t343 * t361 + t358;
t320 = -t331 * t363 + t332 * t347 - t352 * t357;
t322 = t347 * t367 + (-t362 * t341 - t368) * t352;
t355 = g(1) * t320 + g(2) * t318 + g(3) * t322;
t351 = cos(qJ(4));
t350 = cos(qJ(6));
t345 = sin(qJ(6));
t344 = -qJ(5) - pkin(11);
t339 = qJ(4) + pkin(13);
t338 = cos(t339);
t337 = sin(t339);
t336 = t351 * pkin(4) + pkin(3);
t323 = t347 * t368 + (t347 * t362 + t348 * t352) * t341;
t321 = t332 * t352 + (t331 * t342 + t357) * t347;
t319 = t330 * t352 + (t329 * t342 - t356) * t347;
t317 = t323 * t338 + t328 * t337;
t316 = t321 * t338 + t325 * t337;
t315 = t319 * t338 + t324 * t337;
t1 = [0, -g(1) * t354 - g(2) * t349, g(1) * t349 - g(2) * t354, 0, 0, 0, 0, 0, -g(1) * t332 - g(2) * t330 - g(3) * t367, -g(1) * t331 - g(2) * t329 - g(3) * t365, 0, 0, 0, 0, 0, -g(1) * t321 - g(2) * t319 - g(3) * t323, t355, 0, 0, 0, 0, 0, -g(1) * (t321 * t351 + t370) - g(2) * (t319 * t351 + t371) - g(3) * (t323 * t351 + t369) -g(1) * (-t321 * t346 + t325 * t351) - g(2) * (-t319 * t346 + t324 * t351) - g(3) * (-t323 * t346 + t328 * t351) -t355, -g(1) * (t354 * pkin(1) + t332 * pkin(2) + pkin(4) * t370 + pkin(9) * t366 - t320 * t344 + t321 * t336) - g(2) * (t349 * pkin(1) + t330 * pkin(2) + pkin(4) * t371 - pkin(9) * t364 - t318 * t344 + t319 * t336) - g(3) * (pkin(2) * t367 + pkin(4) * t369 + t343 * pkin(9) - t322 * t344 + t323 * t336 + pkin(8)) + (-g(1) * t325 - g(2) * t324 - g(3) * t328) * pkin(10), 0, 0, 0, 0, 0, -g(1) * (t316 * t350 + t320 * t345) - g(2) * (t315 * t350 + t318 * t345) - g(3) * (t317 * t350 + t322 * t345) -g(1) * (-t316 * t345 + t320 * t350) - g(2) * (-t315 * t345 + t318 * t350) - g(3) * (-t317 * t345 + t322 * t350);];
U_reg  = t1;
