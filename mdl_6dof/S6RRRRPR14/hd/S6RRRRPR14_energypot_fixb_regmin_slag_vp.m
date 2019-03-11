% Calculate minimal parameter regressor of potential energy for
% S6RRRRPR14
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
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR14_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR14_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:25:29
% EndTime: 2019-03-10 00:25:29
% DurationCPUTime: 0.20s
% Computational Cost: add. (324->84), mult. (844->139), div. (0->0), fcn. (1103->16), ass. (0->51)
t334 = sin(pkin(7));
t338 = cos(pkin(6));
t363 = t334 * t338;
t335 = sin(pkin(6));
t341 = sin(qJ(2));
t362 = t335 * t341;
t342 = sin(qJ(1));
t361 = t335 * t342;
t345 = cos(qJ(2));
t360 = t335 * t345;
t346 = cos(qJ(1));
t359 = t335 * t346;
t337 = cos(pkin(7));
t344 = cos(qJ(3));
t358 = t337 * t344;
t357 = t337 * t345;
t356 = t342 * t341;
t355 = t342 * t345;
t354 = t346 * t341;
t353 = t346 * t345;
t352 = t334 * t361;
t351 = t334 * t359;
t323 = t338 * t353 - t356;
t350 = t323 * t334 + t337 * t359;
t325 = -t338 * t355 - t354;
t349 = -t325 * t334 + t337 * t361;
t348 = t334 * t360 - t338 * t337;
t324 = t338 * t354 + t355;
t340 = sin(qJ(3));
t313 = t324 * t344 + (t323 * t337 - t351) * t340;
t339 = sin(qJ(4));
t343 = cos(qJ(4));
t306 = t313 * t339 + t350 * t343;
t326 = -t338 * t356 + t353;
t315 = t326 * t344 + (t325 * t337 + t352) * t340;
t308 = t315 * t339 - t349 * t343;
t319 = t340 * t363 + (t340 * t357 + t341 * t344) * t335;
t310 = t319 * t339 + t348 * t343;
t347 = g(1) * t308 + g(2) * t306 + g(3) * t310;
t336 = cos(pkin(13));
t333 = sin(pkin(13));
t332 = pkin(13) + qJ(6);
t331 = cos(t332);
t330 = sin(t332);
t318 = t340 * t362 + (-t357 * t335 - t363) * t344;
t314 = -t325 * t358 + t326 * t340 - t344 * t352;
t312 = -t323 * t358 + t324 * t340 + t344 * t351;
t311 = t319 * t343 - t348 * t339;
t309 = t315 * t343 + t349 * t339;
t307 = t313 * t343 - t350 * t339;
t1 = [0, -g(1) * t346 - g(2) * t342, g(1) * t342 - g(2) * t346, 0, 0, 0, 0, 0, -g(1) * t326 - g(2) * t324 - g(3) * t362, -g(1) * t325 - g(2) * t323 - g(3) * t360, 0, 0, 0, 0, 0, -g(1) * t315 - g(2) * t313 - g(3) * t319, g(1) * t314 + g(2) * t312 + g(3) * t318, 0, 0, 0, 0, 0, -g(1) * t309 - g(2) * t307 - g(3) * t311, t347, -g(1) * (t309 * t336 + t314 * t333) - g(2) * (t307 * t336 + t312 * t333) - g(3) * (t311 * t336 + t318 * t333) -g(1) * (-t309 * t333 + t314 * t336) - g(2) * (-t307 * t333 + t312 * t336) - g(3) * (-t311 * t333 + t318 * t336) -t347, -g(1) * (t346 * pkin(1) + t326 * pkin(2) + t315 * pkin(3) + t309 * pkin(4) + pkin(9) * t361 + t314 * pkin(11) + t308 * qJ(5)) - g(2) * (t342 * pkin(1) + t324 * pkin(2) + t313 * pkin(3) + t307 * pkin(4) - pkin(9) * t359 + t312 * pkin(11) + t306 * qJ(5)) - g(3) * (pkin(2) * t362 + t319 * pkin(3) + t311 * pkin(4) + t338 * pkin(9) + t318 * pkin(11) + t310 * qJ(5) + pkin(8)) + (-g(1) * t349 + g(2) * t350 + g(3) * t348) * pkin(10), 0, 0, 0, 0, 0, -g(1) * (t309 * t331 + t314 * t330) - g(2) * (t307 * t331 + t312 * t330) - g(3) * (t311 * t331 + t318 * t330) -g(1) * (-t309 * t330 + t314 * t331) - g(2) * (-t307 * t330 + t312 * t331) - g(3) * (-t311 * t330 + t318 * t331);];
U_reg  = t1;
