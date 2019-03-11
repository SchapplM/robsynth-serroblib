% Calculate minimal parameter regressor of potential energy for
% S6RPRRPR9
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
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:31:33
% EndTime: 2019-03-09 05:31:33
% DurationCPUTime: 0.19s
% Computational Cost: add. (229->79), mult. (568->130), div. (0->0), fcn. (723->16), ass. (0->56)
t349 = cos(pkin(6));
t378 = t349 * qJ(2) + pkin(8);
t347 = cos(pkin(12));
t358 = cos(qJ(1));
t365 = t358 * t347;
t344 = sin(pkin(12));
t354 = sin(qJ(1));
t368 = t354 * t344;
t329 = t349 * t365 - t368;
t345 = sin(pkin(7));
t348 = cos(pkin(7));
t346 = sin(pkin(6));
t370 = t346 * t358;
t324 = -t329 * t345 - t348 * t370;
t352 = sin(qJ(4));
t377 = t324 * t352;
t366 = t358 * t344;
t367 = t354 * t347;
t331 = -t349 * t367 - t366;
t371 = t346 * t354;
t325 = -t331 * t345 + t348 * t371;
t376 = t325 * t352;
t372 = t346 * t347;
t328 = -t345 * t372 + t349 * t348;
t375 = t328 * t352;
t374 = t344 * t346;
t373 = t345 * t349;
t357 = cos(qJ(3));
t369 = t348 * t357;
t364 = t358 * pkin(1) + qJ(2) * t371;
t363 = t345 * t371;
t362 = t345 * t370;
t361 = g(1) * t354 - g(2) * t358;
t360 = t354 * pkin(1) - qJ(2) * t370;
t330 = t349 * t366 + t367;
t353 = sin(qJ(3));
t318 = -t329 * t369 + t330 * t353 + t357 * t362;
t332 = -t349 * t368 + t365;
t320 = -t331 * t369 + t332 * t353 - t357 * t363;
t322 = t353 * t374 - t357 * t373 - t369 * t372;
t359 = g(1) * t320 + g(2) * t318 + g(3) * t322;
t356 = cos(qJ(4));
t355 = cos(qJ(6));
t351 = sin(qJ(6));
t350 = -qJ(5) - pkin(10);
t343 = qJ(4) + pkin(13);
t339 = cos(t343);
t338 = sin(t343);
t337 = t356 * pkin(4) + pkin(3);
t323 = t353 * t373 + (t347 * t348 * t353 + t344 * t357) * t346;
t321 = t332 * t357 + (t331 * t348 + t363) * t353;
t319 = t330 * t357 + (t329 * t348 - t362) * t353;
t317 = t323 * t339 + t328 * t338;
t316 = t321 * t339 + t325 * t338;
t315 = t319 * t339 + t324 * t338;
t1 = [0, -g(1) * t358 - g(2) * t354, t361, -g(1) * t332 - g(2) * t330 - g(3) * t374, -g(1) * t331 - g(2) * t329 - g(3) * t372, -g(3) * t349 - t361 * t346, -g(1) * t364 - g(2) * t360 - g(3) * t378, 0, 0, 0, 0, 0, -g(1) * t321 - g(2) * t319 - g(3) * t323, t359, 0, 0, 0, 0, 0, -g(1) * (t321 * t356 + t376) - g(2) * (t319 * t356 + t377) - g(3) * (t323 * t356 + t375) -g(1) * (-t321 * t352 + t325 * t356) - g(2) * (-t319 * t352 + t324 * t356) - g(3) * (-t323 * t352 + t328 * t356) -t359, -g(1) * (t332 * pkin(2) + pkin(4) * t376 - t320 * t350 + t321 * t337 + t364) - g(2) * (t330 * pkin(2) + pkin(4) * t377 - t318 * t350 + t319 * t337 + t360) - g(3) * (pkin(2) * t374 + pkin(4) * t375 - t322 * t350 + t323 * t337 + t378) + (-g(1) * t325 - g(2) * t324 - g(3) * t328) * pkin(9), 0, 0, 0, 0, 0, -g(1) * (t316 * t355 + t320 * t351) - g(2) * (t315 * t355 + t318 * t351) - g(3) * (t317 * t355 + t322 * t351) -g(1) * (-t316 * t351 + t320 * t355) - g(2) * (-t315 * t351 + t318 * t355) - g(3) * (-t317 * t351 + t322 * t355);];
U_reg  = t1;
