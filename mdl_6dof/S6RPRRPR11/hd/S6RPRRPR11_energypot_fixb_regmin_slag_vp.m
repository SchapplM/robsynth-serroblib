% Calculate minimal parameter regressor of potential energy for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:44:31
% EndTime: 2019-03-09 05:44:31
% DurationCPUTime: 0.24s
% Computational Cost: add. (331->87), mult. (859->144), div. (0->0), fcn. (1115->16), ass. (0->54)
t341 = cos(pkin(6));
t367 = t341 * qJ(2) + pkin(8);
t335 = sin(pkin(12));
t337 = sin(pkin(6));
t366 = t335 * t337;
t336 = sin(pkin(7));
t365 = t336 * t341;
t339 = cos(pkin(12));
t364 = t337 * t339;
t344 = sin(qJ(1));
t363 = t337 * t344;
t347 = cos(qJ(1));
t362 = t337 * t347;
t340 = cos(pkin(7));
t346 = cos(qJ(3));
t361 = t340 * t346;
t360 = t344 * t335;
t359 = t344 * t339;
t358 = t347 * t335;
t357 = t347 * t339;
t356 = t347 * pkin(1) + qJ(2) * t363;
t355 = t336 * t363;
t354 = t336 * t362;
t353 = g(1) * t344 - g(2) * t347;
t352 = t344 * pkin(1) - qJ(2) * t362;
t320 = t341 * t357 - t360;
t351 = t320 * t336 + t340 * t362;
t322 = -t341 * t359 - t358;
t350 = -t322 * t336 + t340 * t363;
t349 = t336 * t364 - t341 * t340;
t321 = t341 * t358 + t359;
t343 = sin(qJ(3));
t310 = t321 * t346 + (t320 * t340 - t354) * t343;
t342 = sin(qJ(4));
t345 = cos(qJ(4));
t303 = t310 * t342 + t351 * t345;
t323 = -t341 * t360 + t357;
t312 = t323 * t346 + (t322 * t340 + t355) * t343;
t305 = t312 * t342 - t350 * t345;
t316 = t343 * t365 + (t339 * t340 * t343 + t335 * t346) * t337;
t307 = t316 * t342 + t349 * t345;
t348 = g(1) * t305 + g(2) * t303 + g(3) * t307;
t338 = cos(pkin(13));
t334 = sin(pkin(13));
t333 = pkin(13) + qJ(6);
t329 = cos(t333);
t328 = sin(t333);
t315 = t343 * t366 - t346 * t365 - t361 * t364;
t311 = -t322 * t361 + t323 * t343 - t346 * t355;
t309 = -t320 * t361 + t321 * t343 + t346 * t354;
t308 = t316 * t345 - t349 * t342;
t306 = t312 * t345 + t350 * t342;
t304 = t310 * t345 - t351 * t342;
t1 = [0, -g(1) * t347 - g(2) * t344, t353, -g(1) * t323 - g(2) * t321 - g(3) * t366, -g(1) * t322 - g(2) * t320 - g(3) * t364, -g(3) * t341 - t353 * t337, -g(1) * t356 - g(2) * t352 - g(3) * t367, 0, 0, 0, 0, 0, -g(1) * t312 - g(2) * t310 - g(3) * t316, g(1) * t311 + g(2) * t309 + g(3) * t315, 0, 0, 0, 0, 0, -g(1) * t306 - g(2) * t304 - g(3) * t308, t348, -g(1) * (t306 * t338 + t311 * t334) - g(2) * (t304 * t338 + t309 * t334) - g(3) * (t308 * t338 + t315 * t334) -g(1) * (-t306 * t334 + t311 * t338) - g(2) * (-t304 * t334 + t309 * t338) - g(3) * (-t308 * t334 + t315 * t338) -t348, -g(1) * (t323 * pkin(2) + t312 * pkin(3) + t306 * pkin(4) + t311 * pkin(10) + t305 * qJ(5) + t356) - g(2) * (t321 * pkin(2) + t310 * pkin(3) + t304 * pkin(4) + t309 * pkin(10) + t303 * qJ(5) + t352) - g(3) * (pkin(2) * t366 + t316 * pkin(3) + t308 * pkin(4) + t315 * pkin(10) + t307 * qJ(5) + t367) + (-g(1) * t350 + g(2) * t351 + g(3) * t349) * pkin(9), 0, 0, 0, 0, 0, -g(1) * (t306 * t329 + t311 * t328) - g(2) * (t304 * t329 + t309 * t328) - g(3) * (t308 * t329 + t315 * t328) -g(1) * (-t306 * t328 + t311 * t329) - g(2) * (-t304 * t328 + t309 * t329) - g(3) * (-t308 * t328 + t315 * t329);];
U_reg  = t1;
