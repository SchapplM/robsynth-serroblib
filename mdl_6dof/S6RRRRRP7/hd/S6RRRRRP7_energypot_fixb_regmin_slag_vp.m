% Calculate minimal parameter regressor of potential energy for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:43:09
% EndTime: 2019-03-10 01:43:09
% DurationCPUTime: 0.16s
% Computational Cost: add. (162->68), mult. (269->105), div. (0->0), fcn. (331->12), ass. (0->42)
t262 = sin(qJ(1));
t266 = cos(qJ(1));
t284 = -g(1) * t262 + g(2) * t266;
t257 = cos(pkin(6));
t265 = cos(qJ(2));
t270 = t266 * t265;
t261 = sin(qJ(2));
t273 = t262 * t261;
t244 = -t257 * t270 + t273;
t259 = sin(qJ(5));
t281 = t244 * t259;
t271 = t266 * t261;
t272 = t262 * t265;
t246 = t257 * t272 + t271;
t280 = t246 * t259;
t256 = sin(pkin(6));
t279 = t256 * t261;
t278 = t256 * t262;
t264 = cos(qJ(3));
t277 = t256 * t264;
t276 = t256 * t265;
t275 = t256 * t266;
t260 = sin(qJ(3));
t274 = t257 * t260;
t245 = t257 * t271 + t272;
t255 = qJ(3) + qJ(4);
t253 = sin(t255);
t254 = cos(t255);
t238 = t245 * t253 + t254 * t275;
t247 = -t257 * t273 + t270;
t240 = t247 * t253 - t254 * t278;
t242 = t253 * t279 - t257 * t254;
t268 = g(1) * t240 + g(2) * t238 + g(3) * t242;
t267 = -pkin(10) - pkin(9);
t263 = cos(qJ(5));
t258 = -qJ(6) - pkin(11);
t252 = t264 * pkin(3) + pkin(2);
t251 = t263 * pkin(5) + pkin(4);
t243 = t257 * t253 + t254 * t279;
t241 = t247 * t254 + t253 * t278;
t239 = t245 * t254 - t253 * t275;
t1 = [0, -g(1) * t266 - g(2) * t262, -t284, 0, 0, 0, 0, 0, -g(1) * t247 - g(2) * t245 - g(3) * t279, g(1) * t246 + g(2) * t244 - g(3) * t276, 0, 0, 0, 0, 0, -g(1) * (t247 * t264 + t260 * t278) - g(2) * (t245 * t264 - t260 * t275) - g(3) * (t261 * t277 + t274) -g(1) * (-t247 * t260 + t262 * t277) - g(2) * (-t245 * t260 - t264 * t275) - g(3) * (t257 * t264 - t260 * t279) 0, 0, 0, 0, 0, -g(1) * t241 - g(2) * t239 - g(3) * t243, t268, 0, 0, 0, 0, 0, -g(1) * (t241 * t263 + t280) - g(2) * (t239 * t263 + t281) - g(3) * (t243 * t263 - t259 * t276) -g(1) * (-t241 * t259 + t246 * t263) - g(2) * (-t239 * t259 + t244 * t263) - g(3) * (-t243 * t259 - t263 * t276) -t268, -g(1) * (t266 * pkin(1) + pkin(5) * t280 - t240 * t258 + t241 * t251 - t246 * t267 + t247 * t252) - g(2) * (t262 * pkin(1) + pkin(5) * t281 - t238 * t258 + t239 * t251 - t244 * t267 + t245 * t252) - g(3) * (pkin(3) * t274 + t257 * pkin(8) - t242 * t258 + t243 * t251 + pkin(7)) + (-g(3) * (t252 * t261 + (-pkin(5) * t259 + t267) * t265) + t284 * (pkin(3) * t260 + pkin(8))) * t256;];
U_reg  = t1;
