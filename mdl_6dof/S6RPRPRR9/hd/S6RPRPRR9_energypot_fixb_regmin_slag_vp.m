% Calculate minimal parameter regressor of potential energy for
% S6RPRPRR9
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
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:05:05
% EndTime: 2019-03-09 04:05:05
% DurationCPUTime: 0.21s
% Computational Cost: add. (241->79), mult. (625->138), div. (0->0), fcn. (804->16), ass. (0->59)
t352 = sin(qJ(3));
t374 = pkin(3) * t352;
t353 = sin(qJ(1));
t373 = g(1) * t353;
t357 = cos(qJ(1));
t372 = g(2) * t357;
t371 = pkin(9) + qJ(4);
t349 = cos(pkin(6));
t370 = t349 * qJ(2) + pkin(8);
t345 = sin(pkin(6));
t347 = cos(pkin(12));
t369 = t345 * t347;
t368 = t345 * t353;
t367 = t345 * t357;
t343 = sin(pkin(12));
t366 = t353 * t343;
t365 = t353 * t347;
t364 = t357 * t343;
t363 = t357 * t347;
t362 = t357 * pkin(1) + qJ(2) * t368;
t361 = -t372 + t373;
t342 = sin(pkin(13));
t346 = cos(pkin(13));
t356 = cos(qJ(3));
t360 = t356 * t342 + t352 * t346;
t336 = -t352 * t342 + t356 * t346;
t332 = t349 * t364 + t365;
t334 = -t349 * t366 + t363;
t359 = g(3) * t345 * t343 + g(1) * t334 + g(2) * t332;
t331 = t349 * t363 - t366;
t333 = -t349 * t365 - t364;
t344 = sin(pkin(7));
t348 = cos(pkin(7));
t358 = -g(1) * (t333 * t348 + t344 * t368) - g(2) * (t331 * t348 - t344 * t367) - g(3) * (t344 * t349 + t348 * t369);
t355 = cos(qJ(5));
t354 = cos(qJ(6));
t351 = sin(qJ(5));
t350 = sin(qJ(6));
t340 = t353 * pkin(1);
t338 = t356 * pkin(3) + pkin(2);
t330 = -t344 * t369 + t349 * t348;
t329 = -t371 * t344 + t348 * t374;
t328 = t344 * t374 + t371 * t348;
t327 = t360 * t348;
t326 = t336 * t348;
t325 = t360 * t344;
t324 = t336 * t344;
t323 = -t333 * t344 + t348 * t368;
t322 = -t331 * t344 - t348 * t367;
t321 = t349 * t325 + (t327 * t347 + t336 * t343) * t345;
t320 = -t349 * t324 + (-t326 * t347 + t343 * t360) * t345;
t319 = t325 * t368 + t333 * t327 + t334 * t336;
t318 = -t324 * t368 - t333 * t326 + t334 * t360;
t317 = -t325 * t367 + t331 * t327 + t332 * t336;
t316 = t324 * t367 - t331 * t326 + t332 * t360;
t315 = t321 * t355 + t330 * t351;
t314 = t319 * t355 + t323 * t351;
t313 = t317 * t355 + t322 * t351;
t1 = [0, -g(1) * t357 - g(2) * t353, t361, -t359, -g(1) * t333 - g(2) * t331 - g(3) * t369, -g(3) * t349 - t361 * t345, -g(1) * t362 - g(2) * (-qJ(2) * t367 + t340) - g(3) * t370, 0, 0, 0, 0, 0, t358 * t352 - t359 * t356, t359 * t352 + t358 * t356, -g(1) * t323 - g(2) * t322 - g(3) * t330, -g(1) * (t333 * t329 + t334 * t338 + t362) - g(2) * (t331 * t329 + t332 * t338 + t340) - g(3) * (t349 * t328 + t370) + (-t328 * t373 - g(3) * (t329 * t347 + t338 * t343) - (-qJ(2) - t328) * t372) * t345, 0, 0, 0, 0, 0, -g(1) * t314 - g(2) * t313 - g(3) * t315, -g(1) * (-t319 * t351 + t323 * t355) - g(2) * (-t317 * t351 + t322 * t355) - g(3) * (-t321 * t351 + t330 * t355) 0, 0, 0, 0, 0, -g(1) * (t314 * t354 + t318 * t350) - g(2) * (t313 * t354 + t316 * t350) - g(3) * (t315 * t354 + t320 * t350) -g(1) * (-t314 * t350 + t318 * t354) - g(2) * (-t313 * t350 + t316 * t354) - g(3) * (-t315 * t350 + t320 * t354);];
U_reg  = t1;
