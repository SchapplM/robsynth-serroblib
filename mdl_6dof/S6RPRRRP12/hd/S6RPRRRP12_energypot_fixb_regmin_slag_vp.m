% Calculate minimal parameter regressor of potential energy for
% S6RPRRRP12
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
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:50:38
% EndTime: 2019-03-09 06:50:38
% DurationCPUTime: 0.21s
% Computational Cost: add. (387->83), mult. (1043->132), div. (0->0), fcn. (1361->14), ass. (0->59)
t336 = cos(pkin(6));
t365 = qJ(2) * t336 + pkin(8);
t331 = sin(pkin(12));
t333 = sin(pkin(6));
t364 = t331 * t333;
t332 = sin(pkin(7));
t363 = t332 * t336;
t334 = cos(pkin(12));
t362 = t333 * t334;
t340 = sin(qJ(1));
t361 = t333 * t340;
t344 = cos(qJ(1));
t360 = t333 * t344;
t335 = cos(pkin(7));
t343 = cos(qJ(3));
t359 = t335 * t343;
t358 = t340 * t331;
t357 = t340 * t334;
t356 = t344 * t331;
t355 = t344 * t334;
t354 = pkin(1) * t344 + qJ(2) * t361;
t353 = t332 * t361;
t352 = t332 * t360;
t351 = g(1) * t340 - g(2) * t344;
t350 = pkin(1) * t340 - qJ(2) * t360;
t320 = t336 * t355 - t358;
t349 = t320 * t332 + t335 * t360;
t322 = -t336 * t357 - t356;
t348 = -t322 * t332 + t335 * t361;
t347 = t332 * t362 - t335 * t336;
t321 = t336 * t356 + t357;
t339 = sin(qJ(3));
t309 = t321 * t343 + (t320 * t335 - t352) * t339;
t338 = sin(qJ(4));
t342 = cos(qJ(4));
t301 = t309 * t342 - t338 * t349;
t308 = -t320 * t359 + t321 * t339 + t343 * t352;
t337 = sin(qJ(5));
t341 = cos(qJ(5));
t294 = t301 * t337 - t308 * t341;
t323 = -t336 * t358 + t355;
t311 = t323 * t343 + (t322 * t335 + t353) * t339;
t303 = t311 * t342 + t338 * t348;
t310 = -t322 * t359 + t323 * t339 - t343 * t353;
t296 = t303 * t337 - t310 * t341;
t316 = t339 * t363 + (t334 * t335 * t339 + t331 * t343) * t333;
t307 = t316 * t342 - t338 * t347;
t315 = t339 * t364 - t343 * t363 - t359 * t362;
t298 = t307 * t337 - t315 * t341;
t346 = g(1) * t296 + g(2) * t294 + g(3) * t298;
t300 = t309 * t338 + t342 * t349;
t302 = t311 * t338 - t342 * t348;
t306 = t316 * t338 + t342 * t347;
t345 = g(1) * t302 + g(2) * t300 + g(3) * t306;
t299 = t307 * t341 + t315 * t337;
t297 = t303 * t341 + t310 * t337;
t295 = t301 * t341 + t308 * t337;
t293 = -g(1) * t297 - g(2) * t295 - g(3) * t299;
t1 = [0, -g(1) * t344 - g(2) * t340, t351, -g(1) * t323 - g(2) * t321 - g(3) * t364, -g(1) * t322 - g(2) * t320 - g(3) * t362, -g(3) * t336 - t333 * t351, -g(1) * t354 - g(2) * t350 - g(3) * t365, 0, 0, 0, 0, 0, -g(1) * t311 - g(2) * t309 - g(3) * t316, g(1) * t310 + g(2) * t308 + g(3) * t315, 0, 0, 0, 0, 0, -g(1) * t303 - g(2) * t301 - g(3) * t307, t345, 0, 0, 0, 0, 0, t293, t346, t293, -t345, -t346, -g(1) * (t323 * pkin(2) + t311 * pkin(3) + t303 * pkin(4) + t297 * pkin(5) + t310 * pkin(10) + t302 * pkin(11) + t296 * qJ(6) + t354) - g(2) * (t321 * pkin(2) + t309 * pkin(3) + t301 * pkin(4) + t295 * pkin(5) + t308 * pkin(10) + t300 * pkin(11) + t294 * qJ(6) + t350) - g(3) * (pkin(2) * t364 + t316 * pkin(3) + t307 * pkin(4) + t299 * pkin(5) + t315 * pkin(10) + t306 * pkin(11) + t298 * qJ(6) + t365) + (-g(1) * t348 + g(2) * t349 + g(3) * t347) * pkin(9);];
U_reg  = t1;
