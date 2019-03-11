% Calculate minimal parameter regressor of potential energy for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:46:04
% EndTime: 2019-03-09 07:46:04
% DurationCPUTime: 0.14s
% Computational Cost: add. (225->64), mult. (585->120), div. (0->0), fcn. (766->16), ass. (0->46)
t326 = sin(pkin(7));
t330 = cos(pkin(6));
t350 = t326 * t330;
t327 = sin(pkin(6));
t328 = cos(pkin(13));
t349 = t327 * t328;
t334 = sin(qJ(1));
t348 = t327 * t334;
t338 = cos(qJ(1));
t347 = t327 * t338;
t329 = cos(pkin(7));
t346 = t328 * t329;
t325 = sin(pkin(13));
t345 = t334 * t325;
t344 = t334 * t328;
t343 = t338 * t325;
t342 = t338 * t328;
t341 = g(1) * t334 - g(2) * t338;
t318 = t330 * t342 - t345;
t340 = -t318 * t329 + t326 * t347;
t320 = -t330 * t344 - t343;
t339 = t320 * t329 + t326 * t348;
t337 = cos(qJ(3));
t336 = cos(qJ(4));
t335 = cos(qJ(5));
t333 = sin(qJ(3));
t332 = sin(qJ(4));
t331 = sin(qJ(5));
t324 = qJ(5) + qJ(6);
t323 = cos(t324);
t322 = sin(t324);
t321 = -t330 * t345 + t342;
t319 = t330 * t343 + t344;
t317 = -t326 * t349 + t330 * t329;
t316 = -t320 * t326 + t329 * t348;
t315 = -t318 * t326 - t329 * t347;
t314 = t333 * t350 + (t325 * t337 + t333 * t346) * t327;
t313 = -t337 * t350 + (t325 * t333 - t337 * t346) * t327;
t312 = t321 * t337 + t339 * t333;
t311 = t321 * t333 - t339 * t337;
t310 = t319 * t337 - t340 * t333;
t309 = t319 * t333 + t340 * t337;
t308 = t314 * t336 + t317 * t332;
t307 = t312 * t336 + t316 * t332;
t306 = t310 * t336 + t315 * t332;
t1 = [0, -g(1) * t338 - g(2) * t334, t341, -g(3) * t327 * t325 - g(1) * t321 - g(2) * t319, -g(1) * t320 - g(2) * t318 - g(3) * t349, -g(3) * t330 - t341 * t327, -g(1) * (t338 * pkin(1) + qJ(2) * t348) - g(2) * (t334 * pkin(1) - qJ(2) * t347) - g(3) * (t330 * qJ(2) + pkin(8)) 0, 0, 0, 0, 0, -g(1) * t312 - g(2) * t310 - g(3) * t314, g(1) * t311 + g(2) * t309 + g(3) * t313, 0, 0, 0, 0, 0, -g(1) * t307 - g(2) * t306 - g(3) * t308, -g(1) * (-t312 * t332 + t316 * t336) - g(2) * (-t310 * t332 + t315 * t336) - g(3) * (-t314 * t332 + t317 * t336) 0, 0, 0, 0, 0, -g(1) * (t307 * t335 + t311 * t331) - g(2) * (t306 * t335 + t309 * t331) - g(3) * (t308 * t335 + t313 * t331) -g(1) * (-t307 * t331 + t311 * t335) - g(2) * (-t306 * t331 + t309 * t335) - g(3) * (-t308 * t331 + t313 * t335) 0, 0, 0, 0, 0, -g(1) * (t307 * t323 + t311 * t322) - g(2) * (t306 * t323 + t309 * t322) - g(3) * (t308 * t323 + t313 * t322) -g(1) * (-t307 * t322 + t311 * t323) - g(2) * (-t306 * t322 + t309 * t323) - g(3) * (-t308 * t322 + t313 * t323);];
U_reg  = t1;
