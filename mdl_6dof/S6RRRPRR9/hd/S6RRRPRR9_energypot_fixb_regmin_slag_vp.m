% Calculate minimal parameter regressor of potential energy for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:10:16
% EndTime: 2019-03-09 19:10:16
% DurationCPUTime: 0.24s
% Computational Cost: add. (234->73), mult. (610->129), div. (0->0), fcn. (792->16), ass. (0->54)
t339 = sin(pkin(7));
t342 = cos(pkin(7));
t364 = pkin(10) + qJ(4);
t346 = sin(qJ(3));
t368 = pkin(3) * t346;
t370 = t339 * t368 + t364 * t342 + pkin(9);
t348 = sin(qJ(1));
t353 = cos(qJ(1));
t369 = -g(1) * t348 + g(2) * t353;
t340 = sin(pkin(6));
t363 = t340 * t348;
t352 = cos(qJ(2));
t362 = t340 * t352;
t361 = t340 * t353;
t347 = sin(qJ(2));
t360 = t347 * t353;
t359 = t348 * t347;
t358 = t348 * t352;
t357 = t352 * t353;
t338 = sin(pkin(13));
t341 = cos(pkin(13));
t351 = cos(qJ(3));
t356 = t338 * t351 + t341 * t346;
t336 = -t338 * t346 + t341 * t351;
t343 = cos(pkin(6));
t332 = t343 * t360 + t358;
t334 = -t343 * t359 + t357;
t355 = g(3) * t340 * t347 + g(1) * t334 + g(2) * t332;
t331 = t343 * t357 - t359;
t333 = -t343 * t358 - t360;
t354 = -g(1) * (t333 * t342 + t339 * t363) - g(2) * (t331 * t342 - t339 * t361) - g(3) * (t339 * t343 + t342 * t362);
t350 = cos(qJ(5));
t349 = cos(qJ(6));
t345 = sin(qJ(5));
t344 = sin(qJ(6));
t337 = pkin(3) * t351 + pkin(2);
t330 = -t339 * t362 + t342 * t343;
t329 = -t364 * t339 + t342 * t368;
t327 = t356 * t342;
t326 = t336 * t342;
t325 = t356 * t339;
t324 = t336 * t339;
t323 = -t333 * t339 + t342 * t363;
t322 = -t331 * t339 - t342 * t361;
t321 = t325 * t343 + (t327 * t352 + t336 * t347) * t340;
t320 = -t324 * t343 + (-t326 * t352 + t347 * t356) * t340;
t319 = t325 * t363 + t327 * t333 + t334 * t336;
t318 = -t324 * t363 - t326 * t333 + t334 * t356;
t317 = -t325 * t361 + t331 * t327 + t332 * t336;
t316 = t324 * t361 - t331 * t326 + t332 * t356;
t315 = t321 * t350 + t330 * t345;
t314 = t319 * t350 + t323 * t345;
t313 = t317 * t350 + t322 * t345;
t1 = [0, -g(1) * t353 - g(2) * t348, -t369, 0, 0, 0, 0, 0, -t355, -g(1) * t333 - g(2) * t331 - g(3) * t362, 0, 0, 0, 0, 0, t354 * t346 - t355 * t351, t355 * t346 + t354 * t351, -g(1) * t323 - g(2) * t322 - g(3) * t330, -g(1) * (pkin(1) * t353 + t333 * t329 + t334 * t337) - g(2) * (t348 * pkin(1) + t331 * t329 + t332 * t337) - g(3) * (t370 * t343 + pkin(8)) + (-g(3) * (t329 * t352 + t337 * t347) + t369 * t370) * t340, 0, 0, 0, 0, 0, -g(1) * t314 - g(2) * t313 - g(3) * t315, -g(1) * (-t319 * t345 + t323 * t350) - g(2) * (-t317 * t345 + t322 * t350) - g(3) * (-t321 * t345 + t330 * t350) 0, 0, 0, 0, 0, -g(1) * (t314 * t349 + t318 * t344) - g(2) * (t313 * t349 + t316 * t344) - g(3) * (t315 * t349 + t320 * t344) -g(1) * (-t314 * t344 + t318 * t349) - g(2) * (-t313 * t344 + t316 * t349) - g(3) * (-t315 * t344 + t320 * t349);];
U_reg  = t1;
