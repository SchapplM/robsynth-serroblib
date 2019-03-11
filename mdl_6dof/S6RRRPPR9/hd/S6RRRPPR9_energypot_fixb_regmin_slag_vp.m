% Calculate minimal parameter regressor of potential energy for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:17:54
% EndTime: 2019-03-09 16:17:54
% DurationCPUTime: 0.14s
% Computational Cost: add. (224->68), mult. (544->106), div. (0->0), fcn. (696->12), ass. (0->44)
t279 = sin(pkin(6));
t284 = sin(qJ(2));
t303 = t279 * t284;
t285 = sin(qJ(1));
t302 = t279 * t285;
t287 = cos(qJ(3));
t301 = t279 * t287;
t288 = cos(qJ(2));
t300 = t279 * t288;
t289 = cos(qJ(1));
t299 = t279 * t289;
t298 = t285 * t284;
t297 = t285 * t288;
t296 = t289 * t284;
t295 = t289 * t288;
t281 = cos(pkin(6));
t267 = t281 * t296 + t297;
t283 = sin(qJ(3));
t256 = t267 * t287 - t283 * t299;
t266 = -t281 * t295 + t298;
t278 = sin(pkin(11));
t280 = cos(pkin(11));
t247 = t256 * t278 - t266 * t280;
t269 = -t281 * t298 + t295;
t258 = t269 * t287 + t283 * t302;
t268 = t281 * t297 + t296;
t249 = t258 * t278 - t268 * t280;
t265 = t281 * t283 + t284 * t301;
t253 = t265 * t278 + t280 * t300;
t294 = g(1) * t249 + g(2) * t247 + g(3) * t253;
t255 = t267 * t283 + t287 * t299;
t257 = t269 * t283 - t285 * t301;
t264 = -t281 * t287 + t283 * t303;
t293 = g(1) * t257 + g(2) * t255 + g(3) * t264;
t292 = t289 * pkin(1) + t269 * pkin(2) + t258 * pkin(3) + pkin(8) * t302 + t268 * pkin(9) + t257 * qJ(4);
t291 = pkin(2) * t303 + t265 * pkin(3) + t281 * pkin(8) - pkin(9) * t300 + t264 * qJ(4) + pkin(7);
t290 = t285 * pkin(1) + t267 * pkin(2) + t256 * pkin(3) - pkin(8) * t299 + t266 * pkin(9) + t255 * qJ(4);
t286 = cos(qJ(6));
t282 = sin(qJ(6));
t254 = t265 * t280 - t278 * t300;
t250 = t258 * t280 + t268 * t278;
t248 = t256 * t280 + t266 * t278;
t245 = -g(1) * t250 - g(2) * t248 - g(3) * t254;
t1 = [0, -g(1) * t289 - g(2) * t285, g(1) * t285 - g(2) * t289, 0, 0, 0, 0, 0, -g(1) * t269 - g(2) * t267 - g(3) * t303, g(1) * t268 + g(2) * t266 - g(3) * t300, 0, 0, 0, 0, 0, -g(1) * t258 - g(2) * t256 - g(3) * t265, t293, t245, t294, -t293, -g(1) * t292 - g(2) * t290 - g(3) * t291, t245, -t293, -t294, -g(1) * (t250 * pkin(4) + t249 * qJ(5) + t292) - g(2) * (t248 * pkin(4) + t247 * qJ(5) + t290) - g(3) * (t254 * pkin(4) + t253 * qJ(5) + t291) 0, 0, 0, 0, 0, -g(1) * (t249 * t282 + t250 * t286) - g(2) * (t247 * t282 + t248 * t286) - g(3) * (t253 * t282 + t254 * t286) -g(1) * (t249 * t286 - t250 * t282) - g(2) * (t247 * t286 - t248 * t282) - g(3) * (t253 * t286 - t254 * t282);];
U_reg  = t1;
