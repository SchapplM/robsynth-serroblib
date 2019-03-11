% Calculate minimal parameter regressor of potential energy for
% S6RRPRPR9
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
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:01:08
% EndTime: 2019-03-09 11:01:08
% DurationCPUTime: 0.16s
% Computational Cost: add. (218->85), mult. (355->131), div. (0->0), fcn. (438->14), ass. (0->46)
t275 = sin(pkin(11));
t300 = pkin(3) * t275;
t282 = sin(qJ(1));
t299 = g(1) * t282;
t284 = cos(qJ(1));
t298 = g(2) * t284;
t279 = cos(pkin(6));
t297 = t279 * pkin(8) + pkin(7);
t276 = sin(pkin(6));
t281 = sin(qJ(2));
t296 = t276 * t281;
t295 = t276 * t282;
t283 = cos(qJ(2));
t294 = t276 * t283;
t293 = t276 * t284;
t292 = t279 * t275;
t291 = t282 * t281;
t290 = t282 * t283;
t289 = t284 * t281;
t288 = t284 * t283;
t287 = t284 * pkin(1) + pkin(8) * t295;
t257 = t279 * t289 + t290;
t273 = pkin(11) + qJ(4);
t266 = sin(t273);
t268 = cos(t273);
t250 = t257 * t266 + t268 * t293;
t259 = -t279 * t291 + t288;
t252 = t259 * t266 - t268 * t295;
t254 = t266 * t296 - t279 * t268;
t286 = g(1) * t252 + g(2) * t250 + g(3) * t254;
t256 = -t279 * t288 + t291;
t258 = t279 * t290 + t289;
t285 = -g(1) * t258 - g(2) * t256 + g(3) * t294;
t280 = -pkin(9) - qJ(3);
t278 = cos(pkin(11));
t277 = cos(pkin(12));
t274 = sin(pkin(12));
t272 = pkin(12) + qJ(6);
t270 = t282 * pkin(1);
t267 = cos(t272);
t265 = sin(t272);
t264 = t278 * pkin(3) + pkin(2);
t255 = t279 * t266 + t268 * t296;
t253 = t259 * t268 + t266 * t295;
t251 = t257 * t268 - t266 * t293;
t1 = [0, -g(1) * t284 - g(2) * t282, -t298 + t299, 0, 0, 0, 0, 0, -g(1) * t259 - g(2) * t257 - g(3) * t296, -t285, -g(1) * (t259 * t278 + t275 * t295) - g(2) * (t257 * t278 - t275 * t293) - g(3) * (t278 * t296 + t292) -g(1) * (-t259 * t275 + t278 * t295) - g(2) * (-t257 * t275 - t278 * t293) - g(3) * (-t275 * t296 + t279 * t278) t285, -g(1) * (t259 * pkin(2) + t258 * qJ(3) + t287) - g(2) * (t257 * pkin(2) - pkin(8) * t293 + t256 * qJ(3) + t270) - g(3) * ((pkin(2) * t281 - qJ(3) * t283) * t276 + t297) 0, 0, 0, 0, 0, -g(1) * t253 - g(2) * t251 - g(3) * t255, t286, -g(1) * (t253 * t277 + t258 * t274) - g(2) * (t251 * t277 + t256 * t274) - g(3) * (t255 * t277 - t274 * t294) -g(1) * (-t253 * t274 + t258 * t277) - g(2) * (-t251 * t274 + t256 * t277) - g(3) * (-t255 * t274 - t277 * t294) -t286, -g(1) * (t253 * pkin(4) + t252 * qJ(5) - t258 * t280 + t259 * t264 + t287) - g(2) * (t251 * pkin(4) + t250 * qJ(5) - t256 * t280 + t257 * t264 + t270) - g(3) * (pkin(3) * t292 + t255 * pkin(4) + t254 * qJ(5) + t297) + (-t299 * t300 - g(3) * (t264 * t281 + t280 * t283) - (-pkin(8) - t300) * t298) * t276, 0, 0, 0, 0, 0, -g(1) * (t253 * t267 + t258 * t265) - g(2) * (t251 * t267 + t256 * t265) - g(3) * (t255 * t267 - t265 * t294) -g(1) * (-t253 * t265 + t258 * t267) - g(2) * (-t251 * t265 + t256 * t267) - g(3) * (-t255 * t265 - t267 * t294);];
U_reg  = t1;
