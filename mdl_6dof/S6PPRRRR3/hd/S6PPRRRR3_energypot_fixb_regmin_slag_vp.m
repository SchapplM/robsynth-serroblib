% Calculate minimal parameter regressor of potential energy for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRRRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_energypot_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:11:40
% EndTime: 2019-03-08 19:11:40
% DurationCPUTime: 0.18s
% Computational Cost: add. (380->63), mult. (1077->123), div. (0->0), fcn. (1431->18), ass. (0->56)
t303 = sin(pkin(13));
t311 = cos(pkin(6));
t331 = t303 * t311;
t305 = sin(pkin(7));
t306 = sin(pkin(6));
t330 = t305 * t306;
t329 = t305 * t311;
t308 = cos(pkin(13));
t328 = t306 * t308;
t310 = cos(pkin(7));
t327 = t306 * t310;
t307 = cos(pkin(14));
t326 = t307 * t310;
t325 = t308 * t311;
t302 = sin(pkin(14));
t299 = t302 * t325 + t303 * t307;
t315 = sin(qJ(3));
t319 = cos(qJ(3));
t298 = -t303 * t302 + t307 * t325;
t321 = t298 * t310 - t305 * t328;
t288 = -t299 * t315 + t321 * t319;
t295 = -t298 * t305 - t308 * t327;
t304 = sin(pkin(8));
t309 = cos(pkin(8));
t324 = t288 * t309 + t295 * t304;
t301 = -t302 * t331 + t308 * t307;
t300 = -t308 * t302 - t307 * t331;
t320 = t300 * t310 + t303 * t330;
t290 = -t301 * t315 + t320 * t319;
t296 = -t300 * t305 + t303 * t327;
t323 = t290 * t309 + t296 * t304;
t293 = t319 * t329 + (-t302 * t315 + t319 * t326) * t306;
t297 = -t307 * t330 + t311 * t310;
t322 = t293 * t309 + t297 * t304;
t318 = cos(qJ(4));
t317 = cos(qJ(5));
t316 = cos(qJ(6));
t314 = sin(qJ(4));
t313 = sin(qJ(5));
t312 = sin(qJ(6));
t294 = t315 * t329 + (t302 * t319 + t315 * t326) * t306;
t292 = -t293 * t304 + t297 * t309;
t291 = t301 * t319 + t320 * t315;
t289 = t299 * t319 + t321 * t315;
t287 = -t290 * t304 + t296 * t309;
t286 = -t288 * t304 + t295 * t309;
t285 = t294 * t318 + t322 * t314;
t284 = t294 * t314 - t322 * t318;
t283 = t285 * t317 + t292 * t313;
t282 = t291 * t318 + t323 * t314;
t281 = t291 * t314 - t323 * t318;
t280 = t289 * t318 + t324 * t314;
t279 = t289 * t314 - t324 * t318;
t278 = t282 * t317 + t287 * t313;
t277 = t280 * t317 + t286 * t313;
t1 = [-g(3) * qJ(1), -g(1) * (t303 * t306 * qJ(2) + t308 * pkin(1)) - g(2) * (t303 * pkin(1) - qJ(2) * t328) - g(3) * (t311 * qJ(2) + qJ(1)) 0, -g(1) * t291 - g(2) * t289 - g(3) * t294, -g(1) * t290 - g(2) * t288 - g(3) * t293, 0, 0, 0, 0, 0, -g(1) * t282 - g(2) * t280 - g(3) * t285, g(1) * t281 + g(2) * t279 + g(3) * t284, 0, 0, 0, 0, 0, -g(1) * t278 - g(2) * t277 - g(3) * t283, -g(1) * (-t282 * t313 + t287 * t317) - g(2) * (-t280 * t313 + t286 * t317) - g(3) * (-t285 * t313 + t292 * t317) 0, 0, 0, 0, 0, -g(1) * (t278 * t316 + t281 * t312) - g(2) * (t277 * t316 + t279 * t312) - g(3) * (t283 * t316 + t284 * t312) -g(1) * (-t278 * t312 + t281 * t316) - g(2) * (-t277 * t312 + t279 * t316) - g(3) * (-t283 * t312 + t284 * t316);];
U_reg  = t1;
