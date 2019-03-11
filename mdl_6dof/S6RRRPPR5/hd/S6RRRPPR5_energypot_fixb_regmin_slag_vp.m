% Calculate minimal parameter regressor of potential energy for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:44:12
% EndTime: 2019-03-09 15:44:13
% DurationCPUTime: 0.14s
% Computational Cost: add. (201->75), mult. (330->119), div. (0->0), fcn. (403->14), ass. (0->44)
t287 = sin(pkin(6));
t292 = sin(qJ(2));
t311 = t287 * t292;
t293 = sin(qJ(1));
t310 = t287 * t293;
t294 = cos(qJ(3));
t309 = t287 * t294;
t295 = cos(qJ(2));
t308 = t287 * t295;
t296 = cos(qJ(1));
t307 = t287 * t296;
t289 = cos(pkin(6));
t291 = sin(qJ(3));
t306 = t289 * t291;
t305 = t293 * t292;
t304 = t293 * t295;
t303 = t296 * t292;
t302 = t296 * t295;
t301 = t291 * t310;
t276 = t294 * pkin(3) + pkin(2);
t290 = -qJ(4) - pkin(9);
t300 = pkin(3) * t306 + t289 * pkin(8) + t276 * t311 + t290 * t308 + pkin(7);
t268 = t289 * t304 + t303;
t269 = -t289 * t305 + t302;
t299 = t296 * pkin(1) + pkin(3) * t301 + pkin(8) * t310 - t268 * t290 + t269 * t276;
t266 = -t289 * t302 + t305;
t298 = -g(1) * t268 - g(2) * t266 + g(3) * t308;
t267 = t289 * t303 + t304;
t297 = t267 * t276 - t266 * t290 + t293 * pkin(1) + (-pkin(3) * t291 - pkin(8)) * t307;
t288 = cos(pkin(12));
t286 = sin(pkin(12));
t285 = qJ(3) + pkin(11);
t284 = pkin(12) + qJ(6);
t280 = cos(t285);
t279 = cos(t284);
t278 = sin(t285);
t277 = sin(t284);
t263 = t289 * t278 + t280 * t311;
t262 = t278 * t311 - t289 * t280;
t259 = t269 * t280 + t278 * t310;
t258 = t269 * t278 - t280 * t310;
t257 = t267 * t280 - t278 * t307;
t256 = t267 * t278 + t280 * t307;
t1 = [0, -g(1) * t296 - g(2) * t293, g(1) * t293 - g(2) * t296, 0, 0, 0, 0, 0, -g(1) * t269 - g(2) * t267 - g(3) * t311, -t298, 0, 0, 0, 0, 0, -g(1) * (t269 * t294 + t301) - g(2) * (t267 * t294 - t291 * t307) - g(3) * (t292 * t309 + t306) -g(1) * (-t269 * t291 + t293 * t309) - g(2) * (-t267 * t291 - t294 * t307) - g(3) * (t289 * t294 - t291 * t311) t298, -g(1) * t299 - g(2) * t297 - g(3) * t300, -g(1) * (t259 * t288 + t268 * t286) - g(2) * (t257 * t288 + t266 * t286) - g(3) * (t263 * t288 - t286 * t308) -g(1) * (-t259 * t286 + t268 * t288) - g(2) * (-t257 * t286 + t266 * t288) - g(3) * (-t263 * t286 - t288 * t308) -g(1) * t258 - g(2) * t256 - g(3) * t262, -g(1) * (t259 * pkin(4) + t258 * qJ(5) + t299) - g(2) * (t257 * pkin(4) + t256 * qJ(5) + t297) - g(3) * (t263 * pkin(4) + t262 * qJ(5) + t300) 0, 0, 0, 0, 0, -g(1) * (t259 * t279 + t268 * t277) - g(2) * (t257 * t279 + t266 * t277) - g(3) * (t263 * t279 - t277 * t308) -g(1) * (-t259 * t277 + t268 * t279) - g(2) * (-t257 * t277 + t266 * t279) - g(3) * (-t263 * t277 - t279 * t308);];
U_reg  = t1;
