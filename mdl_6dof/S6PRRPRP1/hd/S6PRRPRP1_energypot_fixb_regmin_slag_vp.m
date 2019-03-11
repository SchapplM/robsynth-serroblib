% Calculate minimal parameter regressor of potential energy for
% S6PRRPRP1
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
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:26:34
% EndTime: 2019-03-08 21:26:34
% DurationCPUTime: 0.14s
% Computational Cost: add. (163->68), mult. (281->105), div. (0->0), fcn. (335->12), ass. (0->45)
t231 = sin(pkin(10));
t233 = cos(pkin(10));
t239 = sin(qJ(2));
t234 = cos(pkin(6));
t242 = cos(qJ(2));
t249 = t234 * t242;
t213 = t231 * t239 - t233 * t249;
t237 = sin(qJ(5));
t259 = t213 * t237;
t215 = t231 * t249 + t233 * t239;
t258 = t215 * t237;
t232 = sin(pkin(6));
t257 = t231 * t232;
t256 = t232 * t233;
t238 = sin(qJ(3));
t255 = t232 * t238;
t254 = t232 * t239;
t241 = cos(qJ(3));
t253 = t232 * t241;
t252 = t232 * t242;
t251 = t234 * t238;
t250 = t234 * t239;
t248 = t231 * t255;
t247 = t237 * t252;
t216 = -t231 * t250 + t233 * t242;
t224 = t241 * pkin(3) + pkin(2);
t236 = -qJ(4) - pkin(8);
t246 = t233 * pkin(1) + pkin(3) * t248 + pkin(7) * t257 - t215 * t236 + t216 * t224;
t245 = pkin(3) * t251 + t234 * pkin(7) + t224 * t254 + t236 * t252 + qJ(1);
t244 = -g(1) * t215 - g(2) * t213 + g(3) * t252;
t214 = t231 * t242 + t233 * t250;
t243 = t214 * t224 - t213 * t236 + t231 * pkin(1) + (-pkin(3) * t238 - pkin(7)) * t256;
t240 = cos(qJ(5));
t235 = -qJ(6) - pkin(9);
t230 = qJ(3) + pkin(11);
t226 = cos(t230);
t225 = sin(t230);
t223 = t240 * pkin(5) + pkin(4);
t210 = t234 * t225 + t226 * t254;
t209 = t225 * t254 - t234 * t226;
t206 = t216 * t226 + t225 * t257;
t205 = t216 * t225 - t226 * t257;
t204 = t214 * t226 - t225 * t256;
t203 = t214 * t225 + t226 * t256;
t1 = [-g(3) * qJ(1), 0, -g(1) * t216 - g(2) * t214 - g(3) * t254, -t244, 0, 0, 0, 0, 0, -g(1) * (t216 * t241 + t248) - g(2) * (t214 * t241 - t233 * t255) - g(3) * (t239 * t253 + t251) -g(1) * (-t216 * t238 + t231 * t253) - g(2) * (-t214 * t238 - t233 * t253) - g(3) * (t234 * t241 - t238 * t254) t244, -g(1) * t246 - g(2) * t243 - g(3) * t245, 0, 0, 0, 0, 0, -g(1) * (t206 * t240 + t258) - g(2) * (t204 * t240 + t259) - g(3) * (t210 * t240 - t247) -g(1) * (-t206 * t237 + t215 * t240) - g(2) * (-t204 * t237 + t213 * t240) - g(3) * (-t210 * t237 - t240 * t252) -g(1) * t205 - g(2) * t203 - g(3) * t209, -g(1) * (pkin(5) * t258 - t205 * t235 + t206 * t223 + t246) - g(2) * (pkin(5) * t259 - t203 * t235 + t204 * t223 + t243) - g(3) * (-pkin(5) * t247 - t209 * t235 + t210 * t223 + t245);];
U_reg  = t1;
