% Calculate minimal parameter regressor of potential energy for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR15_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:38:39
% EndTime: 2019-03-09 20:38:39
% DurationCPUTime: 0.19s
% Computational Cost: add. (215->68), mult. (581->115), div. (0->0), fcn. (749->14), ass. (0->47)
t318 = cos(pkin(6));
t327 = cos(qJ(2));
t328 = cos(qJ(1));
t333 = t328 * t327;
t322 = sin(qJ(2));
t323 = sin(qJ(1));
t336 = t323 * t322;
t308 = t318 * t333 - t336;
t315 = sin(pkin(7));
t317 = cos(pkin(7));
t316 = sin(pkin(6));
t339 = t316 * t328;
t303 = -t308 * t315 - t317 * t339;
t334 = t328 * t322;
t335 = t323 * t327;
t310 = -t318 * t335 - t334;
t341 = t316 * t323;
t304 = -t310 * t315 + t317 * t341;
t340 = t316 * t327;
t307 = -t315 * t340 + t318 * t317;
t347 = -g(1) * t304 - g(2) * t303 - g(3) * t307;
t343 = t315 * t318;
t342 = t316 * t322;
t326 = cos(qJ(3));
t338 = t317 * t326;
t337 = t317 * t327;
t332 = t315 * t341;
t331 = t315 * t339;
t309 = t318 * t334 + t335;
t321 = sin(qJ(3));
t297 = -t308 * t338 + t309 * t321 + t326 * t331;
t311 = -t318 * t336 + t333;
t299 = -t310 * t338 + t311 * t321 - t326 * t332;
t301 = t321 * t342 + (-t337 * t316 - t343) * t326;
t330 = g(1) * t299 + g(2) * t297 + g(3) * t301;
t298 = t309 * t326 + (t308 * t317 - t331) * t321;
t300 = t311 * t326 + (t310 * t317 + t332) * t321;
t302 = t321 * t343 + (t321 * t337 + t322 * t326) * t316;
t329 = g(1) * t300 + g(2) * t298 + g(3) * t302;
t325 = cos(qJ(5));
t324 = cos(qJ(6));
t320 = sin(qJ(5));
t319 = sin(qJ(6));
t296 = t301 * t320 + t307 * t325;
t295 = t299 * t320 + t304 * t325;
t294 = t297 * t320 + t303 * t325;
t1 = [0, -g(1) * t328 - g(2) * t323, g(1) * t323 - g(2) * t328, 0, 0, 0, 0, 0, -g(1) * t311 - g(2) * t309 - g(3) * t342, -g(1) * t310 - g(2) * t308 - g(3) * t340, 0, 0, 0, 0, 0, -t329, t330, t347, t329, -t330, -g(1) * (t328 * pkin(1) + t311 * pkin(2) + t300 * pkin(3) + pkin(9) * t341 + t299 * qJ(4)) - g(2) * (t323 * pkin(1) + t309 * pkin(2) + t298 * pkin(3) - pkin(9) * t339 + t297 * qJ(4)) - g(3) * (pkin(2) * t342 + t302 * pkin(3) + t318 * pkin(9) + t301 * qJ(4) + pkin(8)) + t347 * pkin(10), 0, 0, 0, 0, 0, -g(1) * t295 - g(2) * t294 - g(3) * t296, -g(1) * (t299 * t325 - t304 * t320) - g(2) * (t297 * t325 - t303 * t320) - g(3) * (t301 * t325 - t307 * t320) 0, 0, 0, 0, 0, -g(1) * (t295 * t324 + t300 * t319) - g(2) * (t294 * t324 + t298 * t319) - g(3) * (t296 * t324 + t302 * t319) -g(1) * (-t295 * t319 + t300 * t324) - g(2) * (-t294 * t319 + t298 * t324) - g(3) * (-t296 * t319 + t302 * t324);];
U_reg  = t1;
