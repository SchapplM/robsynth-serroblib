% Calculate minimal parameter regressor of potential energy for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:32:14
% EndTime: 2019-03-08 21:32:14
% DurationCPUTime: 0.15s
% Computational Cost: add. (226->70), mult. (389->107), div. (0->0), fcn. (481->12), ass. (0->48)
t246 = sin(pkin(10));
t247 = sin(pkin(6));
t271 = t246 * t247;
t248 = cos(pkin(10));
t270 = t247 * t248;
t252 = sin(qJ(3));
t269 = t247 * t252;
t253 = sin(qJ(2));
t268 = t247 * t253;
t255 = cos(qJ(3));
t267 = t247 * t255;
t256 = cos(qJ(2));
t266 = t247 * t256;
t249 = cos(pkin(6));
t265 = t249 * t252;
t264 = t249 * t253;
t263 = t249 * t256;
t262 = t246 * t269;
t231 = t246 * t263 + t248 * t253;
t232 = -t246 * t264 + t248 * t256;
t239 = t255 * pkin(3) + pkin(2);
t250 = -qJ(4) - pkin(8);
t261 = t248 * pkin(1) + pkin(3) * t262 + pkin(7) * t271 - t231 * t250 + t232 * t239;
t260 = pkin(3) * t265 + t249 * pkin(7) + t239 * t268 + t250 * t266 + qJ(1);
t230 = t246 * t256 + t248 * t264;
t245 = qJ(3) + pkin(11);
t240 = sin(t245);
t241 = cos(t245);
t216 = t230 * t241 - t240 * t270;
t229 = t246 * t253 - t248 * t263;
t251 = sin(qJ(5));
t254 = cos(qJ(5));
t211 = t216 * t251 - t229 * t254;
t218 = t232 * t241 + t240 * t271;
t213 = t218 * t251 - t231 * t254;
t224 = t249 * t240 + t241 * t268;
t219 = t224 * t251 + t254 * t266;
t259 = g(1) * t213 + g(2) * t211 + g(3) * t219;
t258 = -g(1) * t231 - g(2) * t229 + g(3) * t266;
t257 = t230 * t239 - t229 * t250 + t246 * pkin(1) + (-pkin(3) * t252 - pkin(7)) * t270;
t223 = t240 * t268 - t249 * t241;
t220 = t224 * t254 - t251 * t266;
t217 = t232 * t240 - t241 * t271;
t215 = t230 * t240 + t241 * t270;
t214 = t218 * t254 + t231 * t251;
t212 = t216 * t254 + t229 * t251;
t210 = -g(1) * t214 - g(2) * t212 - g(3) * t220;
t1 = [-g(3) * qJ(1), 0, -g(1) * t232 - g(2) * t230 - g(3) * t268, -t258, 0, 0, 0, 0, 0, -g(1) * (t232 * t255 + t262) - g(2) * (t230 * t255 - t248 * t269) - g(3) * (t253 * t267 + t265) -g(1) * (-t232 * t252 + t246 * t267) - g(2) * (-t230 * t252 - t248 * t267) - g(3) * (t249 * t255 - t252 * t268) t258, -g(1) * t261 - g(2) * t257 - g(3) * t260, 0, 0, 0, 0, 0, t210, t259, t210, -g(1) * t217 - g(2) * t215 - g(3) * t223, -t259, -g(1) * (t218 * pkin(4) + t214 * pkin(5) + t217 * pkin(9) + t213 * qJ(6) + t261) - g(2) * (t216 * pkin(4) + t212 * pkin(5) + t215 * pkin(9) + t211 * qJ(6) + t257) - g(3) * (t224 * pkin(4) + t220 * pkin(5) + t223 * pkin(9) + t219 * qJ(6) + t260);];
U_reg  = t1;
