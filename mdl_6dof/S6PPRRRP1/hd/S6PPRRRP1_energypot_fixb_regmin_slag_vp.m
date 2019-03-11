% Calculate minimal parameter regressor of potential energy for
% S6PPRRRP1
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
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:54:24
% EndTime: 2019-03-08 18:54:24
% DurationCPUTime: 0.15s
% Computational Cost: add. (261->75), mult. (687->124), div. (0->0), fcn. (886->14), ass. (0->53)
t246 = sin(pkin(12));
t247 = sin(pkin(11));
t250 = cos(pkin(12));
t251 = cos(pkin(11));
t253 = cos(pkin(6));
t270 = t251 * t253;
t234 = -t247 * t246 + t250 * t270;
t235 = t246 * t270 + t247 * t250;
t257 = sin(qJ(3));
t249 = sin(pkin(6));
t248 = sin(pkin(7));
t260 = cos(qJ(3));
t274 = t248 * t260;
t266 = t249 * t274;
t252 = cos(pkin(7));
t269 = t252 * t260;
t221 = -t234 * t269 + t235 * t257 + t251 * t266;
t255 = sin(qJ(5));
t280 = t221 * t255;
t275 = t247 * t253;
t236 = -t251 * t246 - t250 * t275;
t237 = -t246 * t275 + t251 * t250;
t223 = -t236 * t269 + t237 * t257 - t247 * t266;
t279 = t223 * t255;
t273 = t249 * t250;
t277 = t246 * t249;
t229 = -t253 * t274 + t257 * t277 - t269 * t273;
t278 = t229 * t255;
t276 = t247 * t249;
t272 = t249 * t251;
t271 = t249 * t252;
t268 = t253 * qJ(2) + qJ(1);
t267 = t251 * pkin(1) + qJ(2) * t276;
t265 = t247 * pkin(1) - qJ(2) * t272;
t264 = t234 * t248 + t251 * t271;
t263 = -t236 * t248 + t247 * t271;
t262 = t248 * t273 - t253 * t252;
t222 = t235 * t260 + (t234 * t252 - t248 * t272) * t257;
t256 = sin(qJ(4));
t259 = cos(qJ(4));
t217 = t222 * t256 + t264 * t259;
t224 = t237 * t260 + (t236 * t252 + t248 * t276) * t257;
t219 = t224 * t256 - t263 * t259;
t230 = t253 * t248 * t257 + (t250 * t252 * t257 + t246 * t260) * t249;
t225 = t230 * t256 + t262 * t259;
t261 = g(1) * t219 + g(2) * t217 + g(3) * t225;
t258 = cos(qJ(5));
t254 = -qJ(6) - pkin(10);
t242 = t258 * pkin(5) + pkin(4);
t226 = t230 * t259 - t262 * t256;
t220 = t224 * t259 + t263 * t256;
t218 = t222 * t259 - t264 * t256;
t1 = [-g(3) * qJ(1), -g(1) * t267 - g(2) * t265 - g(3) * t268, 0, -g(1) * t224 - g(2) * t222 - g(3) * t230, g(1) * t223 + g(2) * t221 + g(3) * t229, 0, 0, 0, 0, 0, -g(1) * t220 - g(2) * t218 - g(3) * t226, t261, 0, 0, 0, 0, 0, -g(1) * (t220 * t258 + t279) - g(2) * (t218 * t258 + t280) - g(3) * (t226 * t258 + t278) -g(1) * (-t220 * t255 + t223 * t258) - g(2) * (-t218 * t255 + t221 * t258) - g(3) * (-t226 * t255 + t229 * t258) -t261, -g(1) * (t237 * pkin(2) + t224 * pkin(3) + pkin(5) * t279 + t223 * pkin(9) - t219 * t254 + t220 * t242 + t267) - g(2) * (t235 * pkin(2) + t222 * pkin(3) + pkin(5) * t280 + t221 * pkin(9) - t217 * t254 + t218 * t242 + t265) - g(3) * (pkin(2) * t277 + t230 * pkin(3) + pkin(5) * t278 + t229 * pkin(9) - t225 * t254 + t226 * t242 + t268) + (-g(1) * t263 + g(2) * t264 + g(3) * t262) * pkin(8);];
U_reg  = t1;
