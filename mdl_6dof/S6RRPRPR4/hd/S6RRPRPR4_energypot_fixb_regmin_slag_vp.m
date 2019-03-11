% Calculate minimal parameter regressor of potential energy for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:26:16
% EndTime: 2019-03-09 10:26:17
% DurationCPUTime: 0.19s
% Computational Cost: add. (156->66), mult. (327->113), div. (0->0), fcn. (406->14), ass. (0->47)
t316 = pkin(8) + qJ(3);
t289 = sin(pkin(6));
t295 = sin(qJ(2));
t315 = t289 * t295;
t296 = sin(qJ(1));
t314 = t289 * t296;
t300 = cos(qJ(1));
t313 = t289 * t300;
t291 = cos(pkin(6));
t294 = sin(qJ(4));
t312 = t291 * t294;
t311 = t296 * t295;
t299 = cos(qJ(2));
t310 = t296 * t299;
t309 = t300 * t295;
t308 = t300 * t299;
t275 = t291 * t295 * pkin(2) - t316 * t289;
t282 = t299 * pkin(2) + pkin(1);
t307 = t300 * t275 + t296 * t282;
t306 = t294 * t313;
t305 = pkin(2) * t315 + t316 * t291 + pkin(7);
t304 = g(1) * t296 - g(2) * t300;
t288 = sin(pkin(11));
t290 = cos(pkin(11));
t303 = t299 * t288 + t295 * t290;
t302 = t295 * t288 - t299 * t290;
t301 = t302 * t291;
t298 = cos(qJ(4));
t297 = cos(qJ(6));
t293 = sin(qJ(6));
t292 = -qJ(5) - pkin(9);
t287 = qJ(4) + pkin(12);
t284 = cos(t287);
t283 = sin(t287);
t281 = t298 * pkin(4) + pkin(3);
t279 = t300 * t282;
t274 = t303 * t291;
t273 = t303 * t289;
t272 = t302 * t289;
t269 = t273 * t284 + t291 * t283;
t268 = -t296 * t274 - t300 * t302;
t267 = t296 * t301 - t300 * t303;
t266 = t300 * t274 - t296 * t302;
t265 = -t296 * t303 - t300 * t301;
t264 = t268 * t284 + t283 * t314;
t263 = t266 * t284 - t283 * t313;
t1 = [0, -g(1) * t300 - g(2) * t296, t304, 0, 0, 0, 0, 0, -g(1) * (-t291 * t311 + t308) - g(2) * (t291 * t309 + t310) - g(3) * t315, -g(1) * (-t291 * t310 - t309) - g(2) * (t291 * t308 - t311) - g(3) * t289 * t299, -g(3) * t291 - t304 * t289, -g(1) * (-t296 * t275 + t279) - g(2) * t307 - g(3) * t305, 0, 0, 0, 0, 0, -g(1) * (t268 * t298 + t294 * t314) - g(2) * (t266 * t298 - t306) - g(3) * (t273 * t298 + t312) -g(1) * (-t268 * t294 + t298 * t314) - g(2) * (-t266 * t294 - t298 * t313) - g(3) * (-t273 * t294 + t291 * t298) g(1) * t267 + g(2) * t265 - g(3) * t272, -g(1) * (t267 * t292 + t268 * t281 + t279 + (pkin(4) * t289 * t294 - t275) * t296) - g(2) * (-pkin(4) * t306 + t265 * t292 + t266 * t281 + t307) - g(3) * (pkin(4) * t312 - t272 * t292 + t273 * t281 + t305) 0, 0, 0, 0, 0, -g(1) * (t264 * t297 - t267 * t293) - g(2) * (t263 * t297 - t265 * t293) - g(3) * (t269 * t297 + t272 * t293) -g(1) * (-t264 * t293 - t267 * t297) - g(2) * (-t263 * t293 - t265 * t297) - g(3) * (-t269 * t293 + t272 * t297);];
U_reg  = t1;
