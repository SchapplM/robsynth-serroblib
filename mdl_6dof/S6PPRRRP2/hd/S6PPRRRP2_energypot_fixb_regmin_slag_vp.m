% Calculate minimal parameter regressor of potential energy for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRRRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:58:08
% EndTime: 2019-03-08 18:58:08
% DurationCPUTime: 0.16s
% Computational Cost: add. (376->77), mult. (1015->126), div. (0->0), fcn. (1328->14), ass. (0->56)
t262 = sin(pkin(12));
t265 = sin(pkin(6));
t293 = t262 * t265;
t263 = sin(pkin(11));
t292 = t263 * t265;
t269 = cos(pkin(6));
t291 = t263 * t269;
t264 = sin(pkin(7));
t275 = cos(qJ(3));
t290 = t264 * t275;
t266 = cos(pkin(12));
t289 = t265 * t266;
t267 = cos(pkin(11));
t288 = t265 * t267;
t268 = cos(pkin(7));
t287 = t265 * t268;
t286 = t267 * t269;
t285 = t268 * t275;
t284 = t269 * qJ(2) + qJ(1);
t283 = t267 * pkin(1) + qJ(2) * t292;
t282 = t265 * t290;
t281 = t263 * pkin(1) - qJ(2) * t288;
t251 = -t263 * t262 + t266 * t286;
t280 = t251 * t264 + t267 * t287;
t253 = -t267 * t262 - t266 * t291;
t279 = -t253 * t264 + t263 * t287;
t278 = t264 * t289 - t269 * t268;
t252 = t262 * t286 + t263 * t266;
t272 = sin(qJ(3));
t238 = t252 * t275 + (t251 * t268 - t264 * t288) * t272;
t271 = sin(qJ(4));
t274 = cos(qJ(4));
t230 = t238 * t274 - t280 * t271;
t237 = -t251 * t285 + t252 * t272 + t267 * t282;
t270 = sin(qJ(5));
t273 = cos(qJ(5));
t225 = t230 * t270 - t237 * t273;
t254 = -t262 * t291 + t267 * t266;
t240 = t254 * t275 + (t253 * t268 + t264 * t292) * t272;
t232 = t240 * t274 + t279 * t271;
t239 = -t253 * t285 + t254 * t272 - t263 * t282;
t227 = t232 * t270 - t239 * t273;
t247 = t269 * t264 * t272 + (t266 * t268 * t272 + t262 * t275) * t265;
t242 = t247 * t274 - t278 * t271;
t246 = -t269 * t290 + t272 * t293 - t285 * t289;
t233 = t242 * t270 - t246 * t273;
t277 = g(1) * t227 + g(2) * t225 + g(3) * t233;
t229 = t238 * t271 + t280 * t274;
t231 = t240 * t271 - t279 * t274;
t241 = t247 * t271 + t278 * t274;
t276 = g(1) * t231 + g(2) * t229 + g(3) * t241;
t234 = t242 * t273 + t246 * t270;
t228 = t232 * t273 + t239 * t270;
t226 = t230 * t273 + t237 * t270;
t224 = -g(1) * t228 - g(2) * t226 - g(3) * t234;
t1 = [-g(3) * qJ(1), -g(1) * t283 - g(2) * t281 - g(3) * t284, 0, -g(1) * t240 - g(2) * t238 - g(3) * t247, g(1) * t239 + g(2) * t237 + g(3) * t246, 0, 0, 0, 0, 0, -g(1) * t232 - g(2) * t230 - g(3) * t242, t276, 0, 0, 0, 0, 0, t224, t277, t224, -t276, -t277, -g(1) * (t254 * pkin(2) + t240 * pkin(3) + t232 * pkin(4) + t228 * pkin(5) + t239 * pkin(9) + t231 * pkin(10) + t227 * qJ(6) + t283) - g(2) * (t252 * pkin(2) + t238 * pkin(3) + t230 * pkin(4) + t226 * pkin(5) + t237 * pkin(9) + t229 * pkin(10) + t225 * qJ(6) + t281) - g(3) * (pkin(2) * t293 + t247 * pkin(3) + t242 * pkin(4) + t234 * pkin(5) + t246 * pkin(9) + t241 * pkin(10) + t233 * qJ(6) + t284) + (-g(1) * t279 + g(2) * t280 + g(3) * t278) * pkin(8);];
U_reg  = t1;
