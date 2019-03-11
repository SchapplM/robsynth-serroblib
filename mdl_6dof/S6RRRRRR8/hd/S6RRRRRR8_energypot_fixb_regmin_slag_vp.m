% Calculate minimal parameter regressor of potential energy for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x38]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:06:43
% EndTime: 2019-03-10 05:06:44
% DurationCPUTime: 0.15s
% Computational Cost: add. (204->58), mult. (496->110), div. (0->0), fcn. (654->16), ass. (0->45)
t330 = sin(pkin(7));
t333 = cos(pkin(6));
t354 = t330 * t333;
t331 = sin(pkin(6));
t338 = sin(qJ(1));
t353 = t331 * t338;
t342 = cos(qJ(2));
t352 = t331 * t342;
t343 = cos(qJ(1));
t351 = t331 * t343;
t332 = cos(pkin(7));
t350 = t332 * t342;
t337 = sin(qJ(2));
t349 = t338 * t337;
t348 = t338 * t342;
t347 = t343 * t337;
t346 = t343 * t342;
t323 = t333 * t346 - t349;
t345 = -t323 * t332 + t330 * t351;
t325 = -t333 * t348 - t347;
t344 = t325 * t332 + t330 * t353;
t341 = cos(qJ(3));
t340 = cos(qJ(4));
t339 = cos(qJ(6));
t336 = sin(qJ(3));
t335 = sin(qJ(4));
t334 = sin(qJ(6));
t329 = qJ(4) + qJ(5);
t328 = cos(t329);
t327 = sin(t329);
t326 = -t333 * t349 + t346;
t324 = t333 * t347 + t348;
t322 = -t330 * t352 + t333 * t332;
t321 = -t325 * t330 + t332 * t353;
t320 = -t323 * t330 - t332 * t351;
t319 = t336 * t354 + (t336 * t350 + t337 * t341) * t331;
t318 = -t341 * t354 + (t336 * t337 - t341 * t350) * t331;
t317 = t326 * t341 + t344 * t336;
t316 = t326 * t336 - t344 * t341;
t315 = t324 * t341 - t345 * t336;
t314 = t324 * t336 + t345 * t341;
t313 = t319 * t328 + t322 * t327;
t312 = t317 * t328 + t321 * t327;
t311 = t315 * t328 + t320 * t327;
t1 = [0, -g(1) * t343 - g(2) * t338, g(1) * t338 - g(2) * t343, 0, 0, 0, 0, 0, -g(3) * t331 * t337 - g(1) * t326 - g(2) * t324, -g(1) * t325 - g(2) * t323 - g(3) * t352, 0, 0, 0, 0, 0, -g(1) * t317 - g(2) * t315 - g(3) * t319, g(1) * t316 + g(2) * t314 + g(3) * t318, 0, 0, 0, 0, 0, -g(1) * (t317 * t340 + t321 * t335) - g(2) * (t315 * t340 + t320 * t335) - g(3) * (t319 * t340 + t322 * t335) -g(1) * (-t317 * t335 + t321 * t340) - g(2) * (-t315 * t335 + t320 * t340) - g(3) * (-t319 * t335 + t322 * t340) 0, 0, 0, 0, 0, -g(1) * t312 - g(2) * t311 - g(3) * t313, -g(1) * (-t317 * t327 + t321 * t328) - g(2) * (-t315 * t327 + t320 * t328) - g(3) * (-t319 * t327 + t322 * t328) 0, 0, 0, 0, 0, -g(1) * (t312 * t339 + t316 * t334) - g(2) * (t311 * t339 + t314 * t334) - g(3) * (t313 * t339 + t318 * t334) -g(1) * (-t312 * t334 + t316 * t339) - g(2) * (-t311 * t334 + t314 * t339) - g(3) * (-t313 * t334 + t318 * t339);];
U_reg  = t1;
