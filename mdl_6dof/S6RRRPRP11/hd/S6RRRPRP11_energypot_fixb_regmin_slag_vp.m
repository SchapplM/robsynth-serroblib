% Calculate minimal parameter regressor of potential energy for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:49:02
% EndTime: 2019-03-09 17:49:02
% DurationCPUTime: 0.18s
% Computational Cost: add. (162->66), mult. (370->93), div. (0->0), fcn. (455->10), ass. (0->39)
t252 = cos(qJ(5));
t272 = pkin(5) * t252 + pkin(4) + pkin(9);
t245 = sin(pkin(6));
t250 = sin(qJ(2));
t271 = t245 * t250;
t251 = sin(qJ(1));
t270 = t245 * t251;
t253 = cos(qJ(3));
t269 = t245 * t253;
t254 = cos(qJ(2));
t268 = t245 * t254;
t255 = cos(qJ(1));
t267 = t245 * t255;
t266 = t250 * t251;
t265 = t250 * t255;
t264 = t251 * t254;
t263 = t254 * t255;
t248 = sin(qJ(5));
t262 = pkin(5) * t248 + qJ(4);
t246 = cos(pkin(6));
t249 = sin(qJ(3));
t230 = t246 * t249 + t250 * t269;
t261 = pkin(2) * t271 + t230 * pkin(3) + t246 * pkin(8) + pkin(7);
t234 = -t246 * t266 + t263;
t225 = t234 * t253 + t249 * t270;
t260 = t255 * pkin(1) + t234 * pkin(2) + t225 * pkin(3) + pkin(8) * t270;
t232 = t246 * t265 + t264;
t223 = t232 * t253 - t249 * t267;
t259 = t251 * pkin(1) + t232 * pkin(2) + t223 * pkin(3) - pkin(8) * t267;
t222 = t232 * t249 + t253 * t267;
t224 = t234 * t249 - t251 * t269;
t229 = -t246 * t253 + t249 * t271;
t258 = g(1) * t224 + g(2) * t222 + g(3) * t229;
t257 = g(1) * t225 + g(2) * t223 + g(3) * t230;
t231 = -t246 * t263 + t266;
t233 = t246 * t264 + t265;
t256 = -g(1) * t233 - g(2) * t231 + g(3) * t268;
t247 = -qJ(6) - pkin(10);
t1 = [0, -g(1) * t255 - g(2) * t251, g(1) * t251 - g(2) * t255, 0, 0, 0, 0, 0, -g(1) * t234 - g(2) * t232 - g(3) * t271, -t256, 0, 0, 0, 0, 0, -t257, t258, t256, t257, -t258, -g(1) * (pkin(9) * t233 + qJ(4) * t224 + t260) - g(2) * (t231 * pkin(9) + t222 * qJ(4) + t259) - g(3) * (-pkin(9) * t268 + qJ(4) * t229 + t261) 0, 0, 0, 0, 0, -g(1) * (t224 * t248 + t233 * t252) - g(2) * (t222 * t248 + t231 * t252) - g(3) * (t229 * t248 - t252 * t268) -g(1) * (t224 * t252 - t233 * t248) - g(2) * (t222 * t252 - t231 * t248) - g(3) * (t229 * t252 + t248 * t268) -t257, -g(1) * (t262 * t224 - t225 * t247 + t272 * t233 + t260) - g(2) * (t262 * t222 - t223 * t247 + t272 * t231 + t259) - g(3) * (t262 * t229 - t230 * t247 - t272 * t268 + t261);];
U_reg  = t1;
