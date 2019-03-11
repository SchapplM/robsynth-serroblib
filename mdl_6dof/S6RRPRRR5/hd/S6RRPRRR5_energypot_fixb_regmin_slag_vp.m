% Calculate minimal parameter regressor of potential energy for
% S6RRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:46:47
% EndTime: 2019-03-09 13:46:47
% DurationCPUTime: 0.16s
% Computational Cost: add. (149->57), mult. (340->107), div. (0->0), fcn. (441->14), ass. (0->40)
t298 = pkin(8) + qJ(3);
t277 = sin(pkin(6));
t282 = sin(qJ(2));
t297 = t277 * t282;
t283 = sin(qJ(1));
t296 = t277 * t283;
t287 = cos(qJ(1));
t295 = t277 * t287;
t294 = t283 * t282;
t286 = cos(qJ(2));
t293 = t283 * t286;
t292 = t287 * t282;
t291 = t287 * t286;
t290 = g(1) * t283 - g(2) * t287;
t276 = sin(pkin(12));
t278 = cos(pkin(12));
t289 = t286 * t276 + t282 * t278;
t288 = t282 * t276 - t286 * t278;
t285 = cos(qJ(4));
t284 = cos(qJ(5));
t281 = sin(qJ(4));
t280 = sin(qJ(5));
t279 = cos(pkin(6));
t275 = qJ(5) + qJ(6);
t274 = cos(t275);
t273 = sin(t275);
t272 = t286 * pkin(2) + pkin(1);
t269 = t279 * t282 * pkin(2) - t298 * t277;
t268 = t289 * t279;
t267 = t288 * t279;
t266 = t289 * t277;
t265 = t288 * t277;
t264 = t266 * t285 + t279 * t281;
t263 = -t283 * t268 - t287 * t288;
t262 = -t283 * t267 + t287 * t289;
t261 = t287 * t268 - t283 * t288;
t260 = t287 * t267 + t283 * t289;
t259 = t263 * t285 + t281 * t296;
t258 = t261 * t285 - t281 * t295;
t1 = [0, -g(1) * t287 - g(2) * t283, t290, 0, 0, 0, 0, 0, -g(1) * (-t279 * t294 + t291) - g(2) * (t279 * t292 + t293) - g(3) * t297, -g(1) * (-t279 * t293 - t292) - g(2) * (t279 * t291 - t294) - g(3) * t277 * t286, -g(3) * t279 - t290 * t277, -g(1) * (-t283 * t269 + t287 * t272) - g(2) * (t287 * t269 + t283 * t272) - g(3) * (pkin(2) * t297 + t298 * t279 + pkin(7)) 0, 0, 0, 0, 0, -g(1) * t259 - g(2) * t258 - g(3) * t264, -g(1) * (-t263 * t281 + t285 * t296) - g(2) * (-t261 * t281 - t285 * t295) - g(3) * (-t266 * t281 + t279 * t285) 0, 0, 0, 0, 0, -g(1) * (t259 * t284 + t262 * t280) - g(2) * (t258 * t284 + t260 * t280) - g(3) * (t264 * t284 + t265 * t280) -g(1) * (-t259 * t280 + t262 * t284) - g(2) * (-t258 * t280 + t260 * t284) - g(3) * (-t264 * t280 + t265 * t284) 0, 0, 0, 0, 0, -g(1) * (t259 * t274 + t262 * t273) - g(2) * (t258 * t274 + t260 * t273) - g(3) * (t264 * t274 + t265 * t273) -g(1) * (-t259 * t273 + t262 * t274) - g(2) * (-t258 * t273 + t260 * t274) - g(3) * (-t264 * t273 + t265 * t274);];
U_reg  = t1;
