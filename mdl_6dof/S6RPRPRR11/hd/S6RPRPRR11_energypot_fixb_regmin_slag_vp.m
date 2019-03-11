% Calculate minimal parameter regressor of potential energy for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:15:11
% EndTime: 2019-03-09 04:15:12
% DurationCPUTime: 0.20s
% Computational Cost: add. (267->81), mult. (654->138), div. (0->0), fcn. (840->16), ass. (0->51)
t346 = cos(pkin(6));
t369 = t346 * qJ(2) + pkin(8);
t340 = sin(pkin(12));
t342 = sin(pkin(6));
t368 = t340 * t342;
t341 = sin(pkin(7));
t367 = t341 * t346;
t344 = cos(pkin(12));
t366 = t342 * t344;
t349 = sin(qJ(1));
t365 = t342 * t349;
t352 = cos(qJ(1));
t364 = t342 * t352;
t345 = cos(pkin(7));
t351 = cos(qJ(3));
t363 = t345 * t351;
t362 = t349 * t340;
t361 = t349 * t344;
t360 = t352 * t340;
t359 = t352 * t344;
t358 = t352 * pkin(1) + qJ(2) * t365;
t357 = t341 * t365;
t356 = t341 * t364;
t355 = g(1) * t349 - g(2) * t352;
t354 = t349 * pkin(1) - qJ(2) * t364;
t325 = t346 * t359 - t362;
t320 = -t325 * t341 - t345 * t364;
t327 = -t346 * t361 - t360;
t321 = -t327 * t341 + t345 * t365;
t324 = -t341 * t366 + t346 * t345;
t326 = t346 * t360 + t361;
t348 = sin(qJ(3));
t314 = -t325 * t363 + t326 * t348 + t351 * t356;
t328 = -t346 * t362 + t359;
t316 = -t327 * t363 + t328 * t348 - t351 * t357;
t318 = t348 * t368 - t351 * t367 - t363 * t366;
t353 = g(1) * t316 + g(2) * t314 + g(3) * t318;
t350 = cos(qJ(6));
t347 = sin(qJ(6));
t343 = cos(pkin(13));
t339 = sin(pkin(13));
t338 = pkin(13) + qJ(5);
t334 = cos(t338);
t333 = sin(t338);
t319 = t348 * t367 + (t344 * t345 * t348 + t340 * t351) * t342;
t317 = t328 * t351 + (t327 * t345 + t357) * t348;
t315 = t326 * t351 + (t325 * t345 - t356) * t348;
t313 = t319 * t334 + t324 * t333;
t312 = t317 * t334 + t321 * t333;
t311 = t315 * t334 + t320 * t333;
t1 = [0, -g(1) * t352 - g(2) * t349, t355, -g(1) * t328 - g(2) * t326 - g(3) * t368, -g(1) * t327 - g(2) * t325 - g(3) * t366, -g(3) * t346 - t342 * t355, -g(1) * t358 - g(2) * t354 - g(3) * t369, 0, 0, 0, 0, 0, -g(1) * t317 - g(2) * t315 - g(3) * t319, t353, -g(1) * (t317 * t343 + t321 * t339) - g(2) * (t315 * t343 + t320 * t339) - g(3) * (t319 * t343 + t324 * t339) -g(1) * (-t317 * t339 + t321 * t343) - g(2) * (-t315 * t339 + t320 * t343) - g(3) * (-t319 * t339 + t324 * t343) -t353, -g(1) * (t328 * pkin(2) + t317 * pkin(3) + t316 * qJ(4) + t358) - g(2) * (t326 * pkin(2) + t315 * pkin(3) + t314 * qJ(4) + t354) - g(3) * (pkin(2) * t368 + t319 * pkin(3) + t318 * qJ(4) + t369) + (-g(1) * t321 - g(2) * t320 - g(3) * t324) * pkin(9), 0, 0, 0, 0, 0, -g(1) * t312 - g(2) * t311 - g(3) * t313, -g(1) * (-t317 * t333 + t321 * t334) - g(2) * (-t315 * t333 + t320 * t334) - g(3) * (-t319 * t333 + t324 * t334) 0, 0, 0, 0, 0, -g(1) * (t312 * t350 + t316 * t347) - g(2) * (t311 * t350 + t314 * t347) - g(3) * (t313 * t350 + t318 * t347) -g(1) * (-t312 * t347 + t316 * t350) - g(2) * (-t311 * t347 + t314 * t350) - g(3) * (-t313 * t347 + t318 * t350);];
U_reg  = t1;
