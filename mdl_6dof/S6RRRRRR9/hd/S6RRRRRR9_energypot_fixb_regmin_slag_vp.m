% Calculate minimal parameter regressor of potential energy for
% S6RRRRRR9
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
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:39:12
% EndTime: 2019-03-10 05:39:13
% DurationCPUTime: 0.12s
% Computational Cost: add. (218->58), mult. (570->110), div. (0->0), fcn. (754->16), ass. (0->45)
t327 = sin(pkin(7));
t330 = cos(pkin(6));
t351 = t327 * t330;
t328 = sin(pkin(6));
t335 = sin(qJ(1));
t350 = t328 * t335;
t339 = cos(qJ(2));
t349 = t328 * t339;
t340 = cos(qJ(1));
t348 = t328 * t340;
t329 = cos(pkin(7));
t347 = t329 * t339;
t334 = sin(qJ(2));
t346 = t335 * t334;
t345 = t335 * t339;
t344 = t340 * t334;
t343 = t340 * t339;
t320 = t330 * t343 - t346;
t342 = -t320 * t329 + t327 * t348;
t322 = -t330 * t345 - t344;
t341 = t322 * t329 + t327 * t350;
t338 = cos(qJ(3));
t337 = cos(qJ(4));
t336 = cos(qJ(5));
t333 = sin(qJ(3));
t332 = sin(qJ(4));
t331 = sin(qJ(5));
t326 = qJ(5) + qJ(6);
t325 = cos(t326);
t324 = sin(t326);
t323 = -t330 * t346 + t343;
t321 = t330 * t344 + t345;
t319 = -t327 * t349 + t330 * t329;
t318 = -t322 * t327 + t329 * t350;
t317 = -t320 * t327 - t329 * t348;
t316 = t333 * t351 + (t333 * t347 + t334 * t338) * t328;
t315 = -t338 * t351 + (t333 * t334 - t338 * t347) * t328;
t314 = t323 * t338 + t341 * t333;
t313 = t323 * t333 - t341 * t338;
t312 = t321 * t338 - t342 * t333;
t311 = t321 * t333 + t342 * t338;
t310 = t316 * t337 + t319 * t332;
t309 = t314 * t337 + t318 * t332;
t308 = t312 * t337 + t317 * t332;
t1 = [0, -g(1) * t340 - g(2) * t335, g(1) * t335 - g(2) * t340, 0, 0, 0, 0, 0, -g(3) * t328 * t334 - g(1) * t323 - g(2) * t321, -g(1) * t322 - g(2) * t320 - g(3) * t349, 0, 0, 0, 0, 0, -g(1) * t314 - g(2) * t312 - g(3) * t316, g(1) * t313 + g(2) * t311 + g(3) * t315, 0, 0, 0, 0, 0, -g(1) * t309 - g(2) * t308 - g(3) * t310, -g(1) * (-t314 * t332 + t318 * t337) - g(2) * (-t312 * t332 + t317 * t337) - g(3) * (-t316 * t332 + t319 * t337) 0, 0, 0, 0, 0, -g(1) * (t309 * t336 + t313 * t331) - g(2) * (t308 * t336 + t311 * t331) - g(3) * (t310 * t336 + t315 * t331) -g(1) * (-t309 * t331 + t313 * t336) - g(2) * (-t308 * t331 + t311 * t336) - g(3) * (-t310 * t331 + t315 * t336) 0, 0, 0, 0, 0, -g(1) * (t309 * t325 + t313 * t324) - g(2) * (t308 * t325 + t311 * t324) - g(3) * (t310 * t325 + t315 * t324) -g(1) * (-t309 * t324 + t313 * t325) - g(2) * (-t308 * t324 + t311 * t325) - g(3) * (-t310 * t324 + t315 * t325);];
U_reg  = t1;
