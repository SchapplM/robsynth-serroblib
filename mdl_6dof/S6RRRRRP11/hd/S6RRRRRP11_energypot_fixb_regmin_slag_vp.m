% Calculate minimal parameter regressor of potential energy for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:53:37
% EndTime: 2019-03-10 02:53:37
% DurationCPUTime: 0.21s
% Computational Cost: add. (265->78), mult. (700->125), div. (0->0), fcn. (907->14), ass. (0->53)
t318 = cos(pkin(6));
t328 = cos(qJ(2));
t329 = cos(qJ(1));
t336 = t329 * t328;
t323 = sin(qJ(2));
t324 = sin(qJ(1));
t339 = t324 * t323;
t307 = t318 * t336 - t339;
t337 = t329 * t323;
t338 = t324 * t328;
t308 = t318 * t337 + t338;
t322 = sin(qJ(3));
t327 = cos(qJ(3));
t315 = sin(pkin(7));
t316 = sin(pkin(6));
t342 = t316 * t329;
t334 = t315 * t342;
t317 = cos(pkin(7));
t341 = t317 * t327;
t296 = -t307 * t341 + t308 * t322 + t327 * t334;
t320 = sin(qJ(5));
t349 = t296 * t320;
t309 = -t318 * t338 - t337;
t310 = -t318 * t339 + t336;
t344 = t316 * t324;
t335 = t315 * t344;
t298 = -t309 * t341 + t310 * t322 - t327 * t335;
t348 = t298 * t320;
t340 = t317 * t328;
t345 = t316 * t323;
t346 = t315 * t318;
t302 = t322 * t345 + (-t340 * t316 - t346) * t327;
t347 = t302 * t320;
t343 = t316 * t328;
t333 = t307 * t315 + t317 * t342;
t332 = -t309 * t315 + t317 * t344;
t331 = t315 * t343 - t318 * t317;
t297 = t308 * t327 + (t307 * t317 - t334) * t322;
t321 = sin(qJ(4));
t326 = cos(qJ(4));
t290 = t297 * t321 + t333 * t326;
t299 = t310 * t327 + (t309 * t317 + t335) * t322;
t292 = t299 * t321 - t332 * t326;
t303 = t322 * t346 + (t322 * t340 + t323 * t327) * t316;
t294 = t303 * t321 + t331 * t326;
t330 = g(1) * t292 + g(2) * t290 + g(3) * t294;
t325 = cos(qJ(5));
t319 = -qJ(6) - pkin(12);
t314 = t325 * pkin(5) + pkin(4);
t295 = t303 * t326 - t331 * t321;
t293 = t299 * t326 + t332 * t321;
t291 = t297 * t326 - t333 * t321;
t1 = [0, -g(1) * t329 - g(2) * t324, g(1) * t324 - g(2) * t329, 0, 0, 0, 0, 0, -g(1) * t310 - g(2) * t308 - g(3) * t345, -g(1) * t309 - g(2) * t307 - g(3) * t343, 0, 0, 0, 0, 0, -g(1) * t299 - g(2) * t297 - g(3) * t303, g(1) * t298 + g(2) * t296 + g(3) * t302, 0, 0, 0, 0, 0, -g(1) * t293 - g(2) * t291 - g(3) * t295, t330, 0, 0, 0, 0, 0, -g(1) * (t293 * t325 + t348) - g(2) * (t291 * t325 + t349) - g(3) * (t295 * t325 + t347) -g(1) * (-t293 * t320 + t298 * t325) - g(2) * (-t291 * t320 + t296 * t325) - g(3) * (-t295 * t320 + t302 * t325) -t330, -g(1) * (t329 * pkin(1) + t310 * pkin(2) + t299 * pkin(3) + pkin(5) * t348 + pkin(9) * t344 + t298 * pkin(11) - t292 * t319 + t293 * t314) - g(2) * (t324 * pkin(1) + t308 * pkin(2) + t297 * pkin(3) + pkin(5) * t349 - pkin(9) * t342 + t296 * pkin(11) - t290 * t319 + t291 * t314) - g(3) * (pkin(2) * t345 + t303 * pkin(3) + pkin(5) * t347 + t318 * pkin(9) + t302 * pkin(11) - t294 * t319 + t295 * t314 + pkin(8)) + (-g(1) * t332 + g(2) * t333 + g(3) * t331) * pkin(10);];
U_reg  = t1;
