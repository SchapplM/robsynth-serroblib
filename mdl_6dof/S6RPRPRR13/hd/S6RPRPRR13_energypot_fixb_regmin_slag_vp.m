% Calculate minimal parameter regressor of potential energy for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR13_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:25:20
% EndTime: 2019-03-09 04:25:20
% DurationCPUTime: 0.21s
% Computational Cost: add. (222->71), mult. (596->120), div. (0->0), fcn. (761->14), ass. (0->50)
t325 = cos(pkin(6));
t323 = cos(pkin(12));
t333 = cos(qJ(1));
t341 = t333 * t323;
t320 = sin(pkin(12));
t329 = sin(qJ(1));
t344 = t329 * t320;
t309 = t325 * t341 - t344;
t321 = sin(pkin(7));
t324 = cos(pkin(7));
t322 = sin(pkin(6));
t346 = t322 * t333;
t304 = -t309 * t321 - t324 * t346;
t342 = t333 * t320;
t343 = t329 * t323;
t311 = -t325 * t343 - t342;
t347 = t322 * t329;
t305 = -t311 * t321 + t324 * t347;
t348 = t322 * t323;
t308 = -t321 * t348 + t325 * t324;
t355 = -g(1) * t305 - g(2) * t304 - g(3) * t308;
t351 = t325 * qJ(2) + pkin(8);
t350 = t320 * t322;
t349 = t321 * t325;
t332 = cos(qJ(3));
t345 = t324 * t332;
t340 = t333 * pkin(1) + qJ(2) * t347;
t339 = t321 * t347;
t338 = t321 * t346;
t337 = g(1) * t329 - g(2) * t333;
t336 = t329 * pkin(1) - qJ(2) * t346;
t310 = t325 * t342 + t343;
t328 = sin(qJ(3));
t298 = -t309 * t345 + t310 * t328 + t332 * t338;
t312 = -t325 * t344 + t341;
t300 = -t311 * t345 + t312 * t328 - t332 * t339;
t302 = t328 * t350 - t332 * t349 - t345 * t348;
t335 = g(1) * t300 + g(2) * t298 + g(3) * t302;
t299 = t310 * t332 + (t309 * t324 - t338) * t328;
t301 = t312 * t332 + (t311 * t324 + t339) * t328;
t303 = t328 * t349 + (t323 * t324 * t328 + t320 * t332) * t322;
t334 = g(1) * t301 + g(2) * t299 + g(3) * t303;
t331 = cos(qJ(5));
t330 = cos(qJ(6));
t327 = sin(qJ(5));
t326 = sin(qJ(6));
t297 = t302 * t327 + t308 * t331;
t296 = t300 * t327 + t305 * t331;
t295 = t298 * t327 + t304 * t331;
t1 = [0, -g(1) * t333 - g(2) * t329, t337, -g(1) * t312 - g(2) * t310 - g(3) * t350, -g(1) * t311 - g(2) * t309 - g(3) * t348, -g(3) * t325 - t337 * t322, -g(1) * t340 - g(2) * t336 - g(3) * t351, 0, 0, 0, 0, 0, -t334, t335, t355, t334, -t335, -g(1) * (t312 * pkin(2) + t301 * pkin(3) + t300 * qJ(4) + t340) - g(2) * (t310 * pkin(2) + t299 * pkin(3) + t298 * qJ(4) + t336) - g(3) * (pkin(2) * t350 + t303 * pkin(3) + t302 * qJ(4) + t351) + t355 * pkin(9), 0, 0, 0, 0, 0, -g(1) * t296 - g(2) * t295 - g(3) * t297, -g(1) * (t300 * t331 - t305 * t327) - g(2) * (t298 * t331 - t304 * t327) - g(3) * (t302 * t331 - t308 * t327) 0, 0, 0, 0, 0, -g(1) * (t296 * t330 + t301 * t326) - g(2) * (t295 * t330 + t299 * t326) - g(3) * (t297 * t330 + t303 * t326) -g(1) * (-t296 * t326 + t301 * t330) - g(2) * (-t295 * t326 + t299 * t330) - g(3) * (-t297 * t326 + t303 * t330);];
U_reg  = t1;
