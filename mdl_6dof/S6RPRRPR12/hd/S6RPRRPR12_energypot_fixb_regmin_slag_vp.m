% Calculate minimal parameter regressor of potential energy for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:52:56
% EndTime: 2019-03-09 05:52:57
% DurationCPUTime: 0.18s
% Computational Cost: add. (285->78), mult. (764->127), div. (0->0), fcn. (986->14), ass. (0->52)
t323 = cos(pkin(6));
t352 = t323 * qJ(2) + pkin(8);
t318 = sin(pkin(12));
t320 = sin(pkin(6));
t351 = t318 * t320;
t319 = sin(pkin(7));
t350 = t319 * t323;
t321 = cos(pkin(12));
t349 = t320 * t321;
t327 = sin(qJ(1));
t348 = t320 * t327;
t331 = cos(qJ(1));
t347 = t320 * t331;
t322 = cos(pkin(7));
t330 = cos(qJ(3));
t346 = t322 * t330;
t345 = t323 * t331;
t344 = t327 * t318;
t343 = t327 * t321;
t342 = t331 * pkin(1) + qJ(2) * t348;
t341 = t319 * t348;
t340 = t319 * t347;
t339 = g(1) * t327 - g(2) * t331;
t338 = t327 * pkin(1) - qJ(2) * t347;
t307 = t321 * t345 - t344;
t337 = t307 * t319 + t322 * t347;
t309 = -t318 * t331 - t323 * t343;
t336 = -t309 * t319 + t322 * t348;
t335 = t319 * t349 - t322 * t323;
t308 = t318 * t345 + t343;
t326 = sin(qJ(3));
t297 = t308 * t330 + (t307 * t322 - t340) * t326;
t325 = sin(qJ(4));
t329 = cos(qJ(4));
t290 = t297 * t325 + t337 * t329;
t310 = t321 * t331 - t323 * t344;
t299 = t310 * t330 + (t309 * t322 + t341) * t326;
t292 = t299 * t325 - t336 * t329;
t303 = t326 * t350 + (t321 * t322 * t326 + t318 * t330) * t320;
t294 = t303 * t325 + t335 * t329;
t334 = g(1) * t292 + g(2) * t290 + g(3) * t294;
t291 = t297 * t329 - t337 * t325;
t293 = t299 * t329 + t336 * t325;
t295 = t303 * t329 - t335 * t325;
t333 = g(1) * t293 + g(2) * t291 + g(3) * t295;
t296 = -t307 * t346 + t308 * t326 + t330 * t340;
t298 = -t309 * t346 + t310 * t326 - t330 * t341;
t302 = t326 * t351 - t330 * t350 - t346 * t349;
t332 = g(1) * t298 + g(2) * t296 + g(3) * t302;
t328 = cos(qJ(6));
t324 = sin(qJ(6));
t1 = [0, -g(1) * t331 - g(2) * t327, t339, -g(1) * t310 - g(2) * t308 - g(3) * t351, -g(1) * t309 - g(2) * t307 - g(3) * t349, -g(3) * t323 - t339 * t320, -g(1) * t342 - g(2) * t338 - g(3) * t352, 0, 0, 0, 0, 0, -g(1) * t299 - g(2) * t297 - g(3) * t303, t332, 0, 0, 0, 0, 0, -t333, t334, -t332, t333, -t334, -g(1) * (pkin(2) * t310 + pkin(3) * t299 + pkin(4) * t293 + pkin(10) * t298 + qJ(5) * t292 + t342) - g(2) * (t308 * pkin(2) + t297 * pkin(3) + t291 * pkin(4) + t296 * pkin(10) + t290 * qJ(5) + t338) - g(3) * (pkin(2) * t351 + pkin(3) * t303 + pkin(4) * t295 + pkin(10) * t302 + qJ(5) * t294 + t352) + (-g(1) * t336 + g(2) * t337 + g(3) * t335) * pkin(9), 0, 0, 0, 0, 0, -g(1) * (t292 * t324 + t298 * t328) - g(2) * (t290 * t324 + t296 * t328) - g(3) * (t294 * t324 + t302 * t328) -g(1) * (t292 * t328 - t298 * t324) - g(2) * (t290 * t328 - t296 * t324) - g(3) * (t294 * t328 - t302 * t324);];
U_reg  = t1;
