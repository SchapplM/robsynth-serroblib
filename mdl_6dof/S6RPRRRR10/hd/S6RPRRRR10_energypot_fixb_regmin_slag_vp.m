% Calculate minimal parameter regressor of potential energy for
% S6RPRRRR10
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
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:34:08
% EndTime: 2019-03-09 07:34:09
% DurationCPUTime: 0.18s
% Computational Cost: add. (211->64), mult. (511->121), div. (0->0), fcn. (666->16), ass. (0->45)
t328 = sin(pkin(7));
t332 = cos(pkin(6));
t351 = t328 * t332;
t329 = sin(pkin(6));
t330 = cos(pkin(13));
t350 = t329 * t330;
t336 = sin(qJ(1));
t349 = t329 * t336;
t340 = cos(qJ(1));
t348 = t329 * t340;
t331 = cos(pkin(7));
t347 = t330 * t331;
t346 = t332 * t340;
t327 = sin(pkin(13));
t345 = t336 * t327;
t344 = t336 * t330;
t343 = g(1) * t336 - g(2) * t340;
t320 = t330 * t346 - t345;
t342 = -t320 * t331 + t328 * t348;
t322 = -t327 * t340 - t332 * t344;
t341 = t322 * t331 + t328 * t349;
t339 = cos(qJ(3));
t338 = cos(qJ(4));
t337 = cos(qJ(6));
t335 = sin(qJ(3));
t334 = sin(qJ(4));
t333 = sin(qJ(6));
t326 = qJ(4) + qJ(5);
t325 = cos(t326);
t324 = sin(t326);
t323 = t330 * t340 - t332 * t345;
t321 = t327 * t346 + t344;
t319 = -t328 * t350 + t331 * t332;
t318 = -t322 * t328 + t331 * t349;
t317 = -t320 * t328 - t331 * t348;
t316 = t335 * t351 + (t327 * t339 + t335 * t347) * t329;
t315 = -t339 * t351 + (t327 * t335 - t339 * t347) * t329;
t314 = t323 * t339 + t341 * t335;
t313 = t323 * t335 - t341 * t339;
t312 = t321 * t339 - t342 * t335;
t311 = t321 * t335 + t342 * t339;
t310 = t316 * t325 + t319 * t324;
t309 = t314 * t325 + t318 * t324;
t308 = t312 * t325 + t317 * t324;
t1 = [0, -g(1) * t340 - g(2) * t336, t343, -g(3) * t327 * t329 - g(1) * t323 - g(2) * t321, -g(1) * t322 - g(2) * t320 - g(3) * t350, -g(3) * t332 - t343 * t329, -g(1) * (pkin(1) * t340 + qJ(2) * t349) - g(2) * (t336 * pkin(1) - qJ(2) * t348) - g(3) * (qJ(2) * t332 + pkin(8)) 0, 0, 0, 0, 0, -g(1) * t314 - g(2) * t312 - g(3) * t316, g(1) * t313 + g(2) * t311 + g(3) * t315, 0, 0, 0, 0, 0, -g(1) * (t314 * t338 + t318 * t334) - g(2) * (t312 * t338 + t317 * t334) - g(3) * (t316 * t338 + t319 * t334) -g(1) * (-t314 * t334 + t318 * t338) - g(2) * (-t312 * t334 + t317 * t338) - g(3) * (-t316 * t334 + t319 * t338) 0, 0, 0, 0, 0, -g(1) * t309 - g(2) * t308 - g(3) * t310, -g(1) * (-t314 * t324 + t318 * t325) - g(2) * (-t312 * t324 + t317 * t325) - g(3) * (-t316 * t324 + t319 * t325) 0, 0, 0, 0, 0, -g(1) * (t309 * t337 + t313 * t333) - g(2) * (t308 * t337 + t311 * t333) - g(3) * (t310 * t337 + t315 * t333) -g(1) * (-t309 * t333 + t313 * t337) - g(2) * (-t308 * t333 + t311 * t337) - g(3) * (-t310 * t333 + t315 * t337);];
U_reg  = t1;
