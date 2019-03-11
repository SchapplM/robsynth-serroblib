% Calculate minimal parameter regressor of potential energy for
% S6RRRPRR13
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
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR13_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR13_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:05:08
% EndTime: 2019-03-09 20:05:08
% DurationCPUTime: 0.20s
% Computational Cost: add. (260->78), mult. (639->133), div. (0->0), fcn. (828->16), ass. (0->48)
t335 = sin(pkin(7));
t339 = cos(pkin(6));
t361 = t335 * t339;
t336 = sin(pkin(6));
t342 = sin(qJ(2));
t360 = t336 * t342;
t343 = sin(qJ(1));
t359 = t336 * t343;
t346 = cos(qJ(2));
t358 = t336 * t346;
t347 = cos(qJ(1));
t357 = t336 * t347;
t338 = cos(pkin(7));
t345 = cos(qJ(3));
t356 = t338 * t345;
t355 = t338 * t346;
t354 = t343 * t342;
t353 = t343 * t346;
t352 = t347 * t342;
t351 = t347 * t346;
t350 = t335 * t359;
t349 = t335 * t357;
t324 = t339 * t351 - t354;
t319 = -t324 * t335 - t338 * t357;
t326 = -t339 * t353 - t352;
t320 = -t326 * t335 + t338 * t359;
t323 = -t335 * t358 + t338 * t339;
t325 = t339 * t352 + t353;
t341 = sin(qJ(3));
t313 = -t324 * t356 + t325 * t341 + t345 * t349;
t327 = -t339 * t354 + t351;
t315 = -t326 * t356 + t327 * t341 - t345 * t350;
t317 = t341 * t360 + (-t336 * t355 - t361) * t345;
t348 = g(1) * t315 + g(2) * t313 + g(3) * t317;
t344 = cos(qJ(6));
t340 = sin(qJ(6));
t337 = cos(pkin(13));
t334 = sin(pkin(13));
t333 = pkin(13) + qJ(5);
t332 = cos(t333);
t331 = sin(t333);
t318 = t341 * t361 + (t341 * t355 + t342 * t345) * t336;
t316 = t327 * t345 + (t326 * t338 + t350) * t341;
t314 = t325 * t345 + (t324 * t338 - t349) * t341;
t312 = t318 * t332 + t323 * t331;
t311 = t316 * t332 + t320 * t331;
t310 = t314 * t332 + t319 * t331;
t1 = [0, -g(1) * t347 - g(2) * t343, g(1) * t343 - g(2) * t347, 0, 0, 0, 0, 0, -g(1) * t327 - g(2) * t325 - g(3) * t360, -g(1) * t326 - g(2) * t324 - g(3) * t358, 0, 0, 0, 0, 0, -g(1) * t316 - g(2) * t314 - g(3) * t318, t348, -g(1) * (t316 * t337 + t320 * t334) - g(2) * (t314 * t337 + t319 * t334) - g(3) * (t318 * t337 + t323 * t334) -g(1) * (-t316 * t334 + t320 * t337) - g(2) * (-t314 * t334 + t319 * t337) - g(3) * (-t318 * t334 + t323 * t337) -t348, -g(1) * (t347 * pkin(1) + t327 * pkin(2) + t316 * pkin(3) + pkin(9) * t359 + t315 * qJ(4)) - g(2) * (t343 * pkin(1) + t325 * pkin(2) + t314 * pkin(3) - pkin(9) * t357 + t313 * qJ(4)) - g(3) * (pkin(2) * t360 + t318 * pkin(3) + t339 * pkin(9) + t317 * qJ(4) + pkin(8)) + (-g(1) * t320 - g(2) * t319 - g(3) * t323) * pkin(10), 0, 0, 0, 0, 0, -g(1) * t311 - g(2) * t310 - g(3) * t312, -g(1) * (-t316 * t331 + t320 * t332) - g(2) * (-t314 * t331 + t319 * t332) - g(3) * (-t318 * t331 + t323 * t332) 0, 0, 0, 0, 0, -g(1) * (t311 * t344 + t315 * t340) - g(2) * (t310 * t344 + t313 * t340) - g(3) * (t312 * t344 + t317 * t340) -g(1) * (-t311 * t340 + t315 * t344) - g(2) * (-t310 * t340 + t313 * t344) - g(3) * (-t312 * t340 + t317 * t344);];
U_reg  = t1;
