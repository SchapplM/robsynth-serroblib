% Calculate minimal parameter regressor of potential energy for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:08:30
% EndTime: 2019-03-09 23:08:31
% DurationCPUTime: 0.15s
% Computational Cost: add. (168->65), mult. (281->101), div. (0->0), fcn. (348->12), ass. (0->40)
t259 = sin(qJ(1));
t263 = cos(qJ(1));
t281 = -g(1) * t259 + g(2) * t263;
t254 = sin(pkin(6));
t258 = sin(qJ(2));
t278 = t254 * t258;
t277 = t254 * t259;
t261 = cos(qJ(3));
t276 = t254 * t261;
t262 = cos(qJ(2));
t275 = t254 * t262;
t274 = t254 * t263;
t255 = cos(pkin(6));
t257 = sin(qJ(3));
t273 = t255 * t257;
t272 = t259 * t258;
t271 = t259 * t262;
t270 = t263 * t258;
t269 = t263 * t262;
t244 = t255 * t270 + t271;
t253 = qJ(3) + qJ(4);
t251 = sin(t253);
t252 = cos(t253);
t237 = t244 * t251 + t252 * t274;
t246 = -t255 * t272 + t269;
t239 = t246 * t251 - t252 * t277;
t241 = t251 * t278 - t255 * t252;
t267 = g(1) * t239 + g(2) * t237 + g(3) * t241;
t238 = t244 * t252 - t251 * t274;
t240 = t246 * t252 + t251 * t277;
t242 = t255 * t251 + t252 * t278;
t266 = g(1) * t240 + g(2) * t238 + g(3) * t242;
t243 = -t255 * t269 + t272;
t245 = t255 * t271 + t270;
t265 = -g(1) * t245 - g(2) * t243 + g(3) * t275;
t264 = -pkin(10) - pkin(9);
t260 = cos(qJ(6));
t256 = sin(qJ(6));
t250 = t261 * pkin(3) + pkin(2);
t1 = [0, -g(1) * t263 - g(2) * t259, -t281, 0, 0, 0, 0, 0, -g(1) * t246 - g(2) * t244 - g(3) * t278, -t265, 0, 0, 0, 0, 0, -g(1) * (t246 * t261 + t257 * t277) - g(2) * (t244 * t261 - t257 * t274) - g(3) * (t258 * t276 + t273) -g(1) * (-t246 * t257 + t259 * t276) - g(2) * (-t244 * t257 - t261 * t274) - g(3) * (t255 * t261 - t257 * t278) 0, 0, 0, 0, 0, -t266, t267, t265, t266, -t267, -g(1) * (t263 * pkin(1) + t240 * pkin(4) + t239 * qJ(5) - t245 * t264 + t246 * t250) - g(2) * (t259 * pkin(1) + t238 * pkin(4) + t237 * qJ(5) - t243 * t264 + t244 * t250) - g(3) * (pkin(3) * t273 + t242 * pkin(4) + t255 * pkin(8) + t241 * qJ(5) + pkin(7)) + (-g(3) * (t250 * t258 + t262 * t264) + t281 * (pkin(3) * t257 + pkin(8))) * t254, 0, 0, 0, 0, 0, -g(1) * (t239 * t256 + t245 * t260) - g(2) * (t237 * t256 + t243 * t260) - g(3) * (t241 * t256 - t260 * t275) -g(1) * (t239 * t260 - t245 * t256) - g(2) * (t237 * t260 - t243 * t256) - g(3) * (t241 * t260 + t256 * t275);];
U_reg  = t1;
