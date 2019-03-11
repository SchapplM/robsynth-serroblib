% Calculate minimal parameter regressor of potential energy for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:41:12
% EndTime: 2019-03-09 06:41:12
% DurationCPUTime: 0.22s
% Computational Cost: add. (272->81), mult. (715->131), div. (0->0), fcn. (919->14), ass. (0->55)
t322 = cos(pkin(6));
t353 = t322 * qJ(2) + pkin(8);
t320 = cos(pkin(12));
t317 = sin(pkin(12));
t327 = sin(qJ(1));
t342 = t327 * t317;
t331 = cos(qJ(1));
t343 = t322 * t331;
t305 = t320 * t343 - t342;
t341 = t327 * t320;
t306 = t317 * t343 + t341;
t326 = sin(qJ(3));
t330 = cos(qJ(3));
t318 = sin(pkin(7));
t319 = sin(pkin(6));
t345 = t319 * t331;
t338 = t318 * t345;
t321 = cos(pkin(7));
t344 = t321 * t330;
t294 = -t305 * t344 + t306 * t326 + t330 * t338;
t324 = sin(qJ(5));
t352 = t294 * t324;
t307 = -t317 * t331 - t322 * t341;
t308 = t320 * t331 - t322 * t342;
t346 = t319 * t327;
t339 = t318 * t346;
t296 = -t307 * t344 + t308 * t326 - t330 * t339;
t351 = t296 * t324;
t347 = t319 * t320;
t348 = t318 * t322;
t349 = t317 * t319;
t300 = t326 * t349 - t330 * t348 - t344 * t347;
t350 = t300 * t324;
t340 = t331 * pkin(1) + qJ(2) * t346;
t337 = g(1) * t327 - g(2) * t331;
t336 = t327 * pkin(1) - qJ(2) * t345;
t335 = t305 * t318 + t321 * t345;
t334 = -t307 * t318 + t321 * t346;
t333 = t318 * t347 - t321 * t322;
t295 = t306 * t330 + (t305 * t321 - t338) * t326;
t325 = sin(qJ(4));
t329 = cos(qJ(4));
t288 = t295 * t325 + t335 * t329;
t297 = t308 * t330 + (t307 * t321 + t339) * t326;
t290 = t297 * t325 - t334 * t329;
t301 = t326 * t348 + (t320 * t321 * t326 + t317 * t330) * t319;
t292 = t301 * t325 + t333 * t329;
t332 = g(1) * t290 + g(2) * t288 + g(3) * t292;
t328 = cos(qJ(5));
t323 = -qJ(6) - pkin(11);
t313 = pkin(5) * t328 + pkin(4);
t293 = t301 * t329 - t333 * t325;
t291 = t297 * t329 + t334 * t325;
t289 = t295 * t329 - t335 * t325;
t1 = [0, -g(1) * t331 - g(2) * t327, t337, -g(1) * t308 - g(2) * t306 - g(3) * t349, -g(1) * t307 - g(2) * t305 - g(3) * t347, -g(3) * t322 - t337 * t319, -g(1) * t340 - g(2) * t336 - g(3) * t353, 0, 0, 0, 0, 0, -g(1) * t297 - g(2) * t295 - g(3) * t301, g(1) * t296 + g(2) * t294 + g(3) * t300, 0, 0, 0, 0, 0, -g(1) * t291 - g(2) * t289 - g(3) * t293, t332, 0, 0, 0, 0, 0, -g(1) * (t291 * t328 + t351) - g(2) * (t289 * t328 + t352) - g(3) * (t293 * t328 + t350) -g(1) * (-t291 * t324 + t296 * t328) - g(2) * (-t289 * t324 + t294 * t328) - g(3) * (-t293 * t324 + t300 * t328) -t332, -g(1) * (t308 * pkin(2) + t297 * pkin(3) + pkin(5) * t351 + t296 * pkin(10) - t290 * t323 + t291 * t313 + t340) - g(2) * (t306 * pkin(2) + t295 * pkin(3) + pkin(5) * t352 + t294 * pkin(10) - t288 * t323 + t289 * t313 + t336) - g(3) * (pkin(2) * t349 + t301 * pkin(3) + pkin(5) * t350 + t300 * pkin(10) - t292 * t323 + t293 * t313 + t353) + (-g(1) * t334 + g(2) * t335 + g(3) * t333) * pkin(9);];
U_reg  = t1;
