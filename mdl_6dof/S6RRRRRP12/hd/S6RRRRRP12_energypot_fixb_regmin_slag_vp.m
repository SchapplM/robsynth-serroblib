% Calculate minimal parameter regressor of potential energy for
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:20:37
% EndTime: 2019-03-10 03:20:37
% DurationCPUTime: 0.19s
% Computational Cost: add. (380->80), mult. (1028->127), div. (0->0), fcn. (1349->14), ass. (0->56)
t331 = sin(pkin(7));
t334 = cos(pkin(6));
t362 = t331 * t334;
t332 = sin(pkin(6));
t338 = sin(qJ(2));
t361 = t332 * t338;
t339 = sin(qJ(1));
t360 = t332 * t339;
t343 = cos(qJ(2));
t359 = t332 * t343;
t344 = cos(qJ(1));
t358 = t332 * t344;
t333 = cos(pkin(7));
t342 = cos(qJ(3));
t357 = t333 * t342;
t356 = t333 * t343;
t355 = t339 * t338;
t354 = t339 * t343;
t353 = t344 * t338;
t352 = t344 * t343;
t351 = t331 * t360;
t350 = t331 * t358;
t324 = t334 * t352 - t355;
t349 = t324 * t331 + t333 * t358;
t326 = -t334 * t354 - t353;
t348 = -t326 * t331 + t333 * t360;
t347 = t331 * t359 - t333 * t334;
t325 = t334 * t353 + t354;
t337 = sin(qJ(3));
t313 = t325 * t342 + (t324 * t333 - t350) * t337;
t336 = sin(qJ(4));
t341 = cos(qJ(4));
t305 = t313 * t341 - t336 * t349;
t312 = -t324 * t357 + t325 * t337 + t342 * t350;
t335 = sin(qJ(5));
t340 = cos(qJ(5));
t298 = t305 * t335 - t312 * t340;
t327 = -t334 * t355 + t352;
t315 = t327 * t342 + (t326 * t333 + t351) * t337;
t307 = t315 * t341 + t336 * t348;
t314 = -t326 * t357 + t327 * t337 - t342 * t351;
t300 = t307 * t335 - t314 * t340;
t320 = t337 * t362 + (t337 * t356 + t338 * t342) * t332;
t311 = t320 * t341 - t336 * t347;
t319 = t337 * t361 + (-t332 * t356 - t362) * t342;
t302 = t311 * t335 - t319 * t340;
t346 = g(1) * t300 + g(2) * t298 + g(3) * t302;
t304 = t313 * t336 + t341 * t349;
t306 = t315 * t336 - t341 * t348;
t310 = t320 * t336 + t341 * t347;
t345 = g(1) * t306 + g(2) * t304 + g(3) * t310;
t303 = t311 * t340 + t319 * t335;
t301 = t307 * t340 + t314 * t335;
t299 = t305 * t340 + t312 * t335;
t297 = -g(1) * t301 - g(2) * t299 - g(3) * t303;
t1 = [0, -g(1) * t344 - g(2) * t339, g(1) * t339 - g(2) * t344, 0, 0, 0, 0, 0, -g(1) * t327 - g(2) * t325 - g(3) * t361, -g(1) * t326 - g(2) * t324 - g(3) * t359, 0, 0, 0, 0, 0, -g(1) * t315 - g(2) * t313 - g(3) * t320, g(1) * t314 + g(2) * t312 + g(3) * t319, 0, 0, 0, 0, 0, -g(1) * t307 - g(2) * t305 - g(3) * t311, t345, 0, 0, 0, 0, 0, t297, t346, t297, -t345, -t346, -g(1) * (t344 * pkin(1) + t327 * pkin(2) + t315 * pkin(3) + t307 * pkin(4) + t301 * pkin(5) + pkin(9) * t360 + t314 * pkin(11) + t306 * pkin(12) + t300 * qJ(6)) - g(2) * (t339 * pkin(1) + t325 * pkin(2) + t313 * pkin(3) + t305 * pkin(4) + t299 * pkin(5) - pkin(9) * t358 + t312 * pkin(11) + t304 * pkin(12) + t298 * qJ(6)) - g(3) * (pkin(2) * t361 + t320 * pkin(3) + t311 * pkin(4) + t303 * pkin(5) + t334 * pkin(9) + t319 * pkin(11) + t310 * pkin(12) + t302 * qJ(6) + pkin(8)) + (-g(1) * t348 + g(2) * t349 + g(3) * t347) * pkin(10);];
U_reg  = t1;
