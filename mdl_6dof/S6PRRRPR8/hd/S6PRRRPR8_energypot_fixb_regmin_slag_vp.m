% Calculate minimal parameter regressor of potential energy for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:53:20
% EndTime: 2019-03-08 23:53:20
% DurationCPUTime: 0.19s
% Computational Cost: add. (277->74), mult. (746->123), div. (0->0), fcn. (970->14), ass. (0->49)
t271 = sin(pkin(12));
t273 = sin(pkin(6));
t302 = t271 * t273;
t272 = sin(pkin(7));
t276 = cos(pkin(6));
t301 = t272 * t276;
t274 = cos(pkin(12));
t300 = t273 * t274;
t275 = cos(pkin(7));
t299 = t273 * t275;
t280 = sin(qJ(2));
t298 = t273 * t280;
t283 = cos(qJ(3));
t297 = t273 * t283;
t284 = cos(qJ(2));
t296 = t273 * t284;
t295 = t275 * t283;
t294 = t275 * t284;
t293 = t276 * t280;
t292 = t276 * t284;
t291 = t272 * t297;
t264 = -t271 * t280 + t274 * t292;
t290 = t264 * t272 + t274 * t299;
t266 = -t271 * t292 - t274 * t280;
t289 = -t266 * t272 + t271 * t299;
t288 = t272 * t296 - t276 * t275;
t265 = t271 * t284 + t274 * t293;
t279 = sin(qJ(3));
t252 = t265 * t283 + (t264 * t275 - t272 * t300) * t279;
t278 = sin(qJ(4));
t282 = cos(qJ(4));
t247 = t252 * t278 + t290 * t282;
t267 = -t271 * t293 + t274 * t284;
t254 = t267 * t283 + (t266 * t275 + t272 * t302) * t279;
t249 = t254 * t278 - t289 * t282;
t260 = t279 * t301 + (t279 * t294 + t280 * t283) * t273;
t255 = t260 * t278 + t282 * t288;
t287 = g(1) * t249 + g(2) * t247 + g(3) * t255;
t248 = t252 * t282 - t290 * t278;
t250 = t254 * t282 + t289 * t278;
t256 = t260 * t282 - t278 * t288;
t286 = g(1) * t250 + g(2) * t248 + g(3) * t256;
t251 = -t264 * t295 + t265 * t279 + t274 * t291;
t253 = -t266 * t295 + t267 * t279 - t271 * t291;
t259 = t279 * t298 - t283 * t301 - t294 * t297;
t285 = g(1) * t253 + g(2) * t251 + g(3) * t259;
t281 = cos(qJ(6));
t277 = sin(qJ(6));
t1 = [-g(3) * qJ(1), 0, -g(1) * t267 - g(2) * t265 - g(3) * t298, -g(1) * t266 - g(2) * t264 - g(3) * t296, 0, 0, 0, 0, 0, -g(1) * t254 - g(2) * t252 - g(3) * t260, t285, 0, 0, 0, 0, 0, -t286, t287, -t285, t286, -t287, -g(1) * (t274 * pkin(1) + t267 * pkin(2) + t254 * pkin(3) + t250 * pkin(4) + pkin(8) * t302 + t253 * pkin(10) + t249 * qJ(5)) - g(2) * (t271 * pkin(1) + t265 * pkin(2) + t252 * pkin(3) + t248 * pkin(4) - pkin(8) * t300 + t251 * pkin(10) + t247 * qJ(5)) - g(3) * (pkin(2) * t298 + t260 * pkin(3) + t256 * pkin(4) + t276 * pkin(8) + t259 * pkin(10) + t255 * qJ(5) + qJ(1)) + (-g(1) * t289 + g(2) * t290 + g(3) * t288) * pkin(9), 0, 0, 0, 0, 0, -g(1) * (t249 * t277 + t253 * t281) - g(2) * (t247 * t277 + t251 * t281) - g(3) * (t255 * t277 + t259 * t281) -g(1) * (t249 * t281 - t253 * t277) - g(2) * (t247 * t281 - t251 * t277) - g(3) * (t255 * t281 - t259 * t277);];
U_reg  = t1;
