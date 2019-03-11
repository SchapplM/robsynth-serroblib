% Calculate minimal parameter regressor of potential energy for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR15_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:50:12
% EndTime: 2019-03-10 00:50:12
% DurationCPUTime: 0.20s
% Computational Cost: add. (278->75), mult. (749->121), div. (0->0), fcn. (974->14), ass. (0->50)
t316 = sin(pkin(7));
t319 = cos(pkin(6));
t348 = t316 * t319;
t317 = sin(pkin(6));
t323 = sin(qJ(2));
t347 = t317 * t323;
t324 = sin(qJ(1));
t346 = t317 * t324;
t328 = cos(qJ(2));
t345 = t317 * t328;
t329 = cos(qJ(1));
t344 = t317 * t329;
t318 = cos(pkin(7));
t327 = cos(qJ(3));
t343 = t318 * t327;
t342 = t318 * t328;
t341 = t324 * t323;
t340 = t324 * t328;
t339 = t329 * t323;
t338 = t329 * t328;
t337 = t316 * t346;
t336 = t316 * t344;
t309 = t319 * t338 - t341;
t335 = t309 * t316 + t318 * t344;
t311 = -t319 * t340 - t339;
t334 = -t311 * t316 + t318 * t346;
t333 = t316 * t345 - t319 * t318;
t310 = t319 * t339 + t340;
t322 = sin(qJ(3));
t299 = t310 * t327 + (t309 * t318 - t336) * t322;
t321 = sin(qJ(4));
t326 = cos(qJ(4));
t292 = t299 * t321 + t335 * t326;
t312 = -t319 * t341 + t338;
t301 = t312 * t327 + (t311 * t318 + t337) * t322;
t294 = t301 * t321 - t334 * t326;
t305 = t322 * t348 + (t322 * t342 + t323 * t327) * t317;
t296 = t305 * t321 + t333 * t326;
t332 = g(1) * t294 + g(2) * t292 + g(3) * t296;
t293 = t299 * t326 - t335 * t321;
t295 = t301 * t326 + t334 * t321;
t297 = t305 * t326 - t333 * t321;
t331 = g(1) * t295 + g(2) * t293 + g(3) * t297;
t298 = -t309 * t343 + t310 * t322 + t327 * t336;
t300 = -t311 * t343 + t312 * t322 - t327 * t337;
t304 = t322 * t347 + (-t342 * t317 - t348) * t327;
t330 = g(1) * t300 + g(2) * t298 + g(3) * t304;
t325 = cos(qJ(6));
t320 = sin(qJ(6));
t1 = [0, -g(1) * t329 - g(2) * t324, g(1) * t324 - g(2) * t329, 0, 0, 0, 0, 0, -g(1) * t312 - g(2) * t310 - g(3) * t347, -g(1) * t311 - g(2) * t309 - g(3) * t345, 0, 0, 0, 0, 0, -g(1) * t301 - g(2) * t299 - g(3) * t305, t330, 0, 0, 0, 0, 0, -t331, t332, -t330, t331, -t332, -g(1) * (t329 * pkin(1) + t312 * pkin(2) + t301 * pkin(3) + t295 * pkin(4) + pkin(9) * t346 + t300 * pkin(11) + t294 * qJ(5)) - g(2) * (t324 * pkin(1) + t310 * pkin(2) + t299 * pkin(3) + t293 * pkin(4) - pkin(9) * t344 + t298 * pkin(11) + t292 * qJ(5)) - g(3) * (pkin(2) * t347 + t305 * pkin(3) + t297 * pkin(4) + t319 * pkin(9) + t304 * pkin(11) + t296 * qJ(5) + pkin(8)) + (-g(1) * t334 + g(2) * t335 + g(3) * t333) * pkin(10), 0, 0, 0, 0, 0, -g(1) * (t294 * t320 + t300 * t325) - g(2) * (t292 * t320 + t298 * t325) - g(3) * (t296 * t320 + t304 * t325) -g(1) * (t294 * t325 - t300 * t320) - g(2) * (t292 * t325 - t298 * t320) - g(3) * (t296 * t325 - t304 * t320);];
U_reg  = t1;
