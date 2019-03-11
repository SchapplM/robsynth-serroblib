% Calculate minimal parameter regressor of potential energy for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:46:28
% EndTime: 2019-03-09 19:46:28
% DurationCPUTime: 0.14s
% Computational Cost: add. (170->70), mult. (330->116), div. (0->0), fcn. (420->14), ass. (0->37)
t270 = sin(pkin(6));
t274 = sin(qJ(2));
t288 = t270 * t274;
t275 = sin(qJ(1));
t287 = t270 * t275;
t276 = cos(qJ(3));
t286 = t270 * t276;
t277 = cos(qJ(2));
t285 = t270 * t277;
t278 = cos(qJ(1));
t284 = t270 * t278;
t283 = t275 * t274;
t282 = t275 * t277;
t281 = t278 * t274;
t280 = t278 * t277;
t268 = pkin(12) + qJ(5);
t272 = cos(pkin(6));
t257 = t272 * t281 + t282;
t273 = sin(qJ(3));
t250 = t257 * t273 + t276 * t284;
t259 = -t272 * t283 + t280;
t252 = t259 * t273 - t275 * t286;
t254 = -t272 * t276 + t273 * t288;
t279 = g(1) * t252 + g(2) * t250 + g(3) * t254;
t271 = cos(pkin(12));
t269 = sin(pkin(12));
t267 = qJ(6) + t268;
t266 = cos(t268);
t265 = sin(t268);
t264 = cos(t267);
t263 = sin(t267);
t258 = t272 * t282 + t281;
t256 = -t272 * t280 + t283;
t255 = t272 * t273 + t274 * t286;
t253 = t259 * t276 + t273 * t287;
t251 = t257 * t276 - t273 * t284;
t1 = [0, -g(1) * t278 - g(2) * t275, g(1) * t275 - g(2) * t278, 0, 0, 0, 0, 0, -g(1) * t259 - g(2) * t257 - g(3) * t288, g(1) * t258 + g(2) * t256 - g(3) * t285, 0, 0, 0, 0, 0, -g(1) * t253 - g(2) * t251 - g(3) * t255, t279, -g(1) * (t253 * t271 + t258 * t269) - g(2) * (t251 * t271 + t256 * t269) - g(3) * (t255 * t271 - t269 * t285) -g(1) * (-t253 * t269 + t258 * t271) - g(2) * (-t251 * t269 + t256 * t271) - g(3) * (-t255 * t269 - t271 * t285) -t279, -g(1) * (t278 * pkin(1) + t259 * pkin(2) + t253 * pkin(3) + pkin(8) * t287 + t258 * pkin(9) + t252 * qJ(4)) - g(2) * (t275 * pkin(1) + t257 * pkin(2) + t251 * pkin(3) - pkin(8) * t284 + t256 * pkin(9) + t250 * qJ(4)) - g(3) * (t255 * pkin(3) + t272 * pkin(8) + t254 * qJ(4) + pkin(7) + (pkin(2) * t274 - pkin(9) * t277) * t270) 0, 0, 0, 0, 0, -g(1) * (t253 * t266 + t258 * t265) - g(2) * (t251 * t266 + t256 * t265) - g(3) * (t255 * t266 - t265 * t285) -g(1) * (-t253 * t265 + t258 * t266) - g(2) * (-t251 * t265 + t256 * t266) - g(3) * (-t255 * t265 - t266 * t285) 0, 0, 0, 0, 0, -g(1) * (t253 * t264 + t258 * t263) - g(2) * (t251 * t264 + t256 * t263) - g(3) * (t255 * t264 - t263 * t285) -g(1) * (-t253 * t263 + t258 * t264) - g(2) * (-t251 * t263 + t256 * t264) - g(3) * (-t255 * t263 - t264 * t285);];
U_reg  = t1;
