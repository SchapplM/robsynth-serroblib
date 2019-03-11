% Calculate minimal parameter regressor of potential energy for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:27:56
% EndTime: 2019-03-08 19:27:56
% DurationCPUTime: 0.21s
% Computational Cost: add. (153->64), mult. (319->109), div. (0->0), fcn. (397->14), ass. (0->45)
t262 = pkin(7) + qJ(3);
t235 = qJ(4) + pkin(12);
t231 = sin(t235);
t238 = sin(pkin(6));
t261 = t231 * t238;
t244 = sin(qJ(4));
t260 = t238 * t244;
t245 = sin(qJ(2));
t259 = t238 * t245;
t247 = cos(qJ(4));
t258 = t238 * t247;
t241 = cos(pkin(6));
t257 = t241 * t244;
t256 = t241 * t245;
t248 = cos(qJ(2));
t255 = t241 * t248;
t223 = pkin(2) * t256 - t262 * t238;
t230 = pkin(2) * t248 + pkin(1);
t237 = sin(pkin(10));
t240 = cos(pkin(10));
t254 = t240 * t223 + t237 * t230;
t253 = t240 * t260;
t252 = pkin(2) * t259 + t262 * t241 + qJ(1);
t236 = sin(pkin(11));
t239 = cos(pkin(11));
t251 = t236 * t248 + t245 * t239;
t250 = t245 * t236 - t239 * t248;
t249 = t250 * t241;
t246 = cos(qJ(6));
t243 = sin(qJ(6));
t242 = -qJ(5) - pkin(8);
t232 = cos(t235);
t229 = pkin(4) * t247 + pkin(3);
t227 = t240 * t230;
t222 = t251 * t241;
t221 = t251 * t238;
t220 = t250 * t238;
t217 = t221 * t232 + t231 * t241;
t216 = -t222 * t237 - t240 * t250;
t215 = t237 * t249 - t240 * t251;
t214 = t222 * t240 - t237 * t250;
t213 = -t237 * t251 - t240 * t249;
t212 = t216 * t232 + t237 * t261;
t211 = t214 * t232 - t240 * t261;
t1 = [-g(3) * qJ(1), 0, -g(1) * (-t237 * t256 + t240 * t248) - g(2) * (t237 * t248 + t240 * t256) - g(3) * t259, -g(1) * (-t237 * t255 - t240 * t245) - g(2) * (-t237 * t245 + t240 * t255) - g(3) * t238 * t248, -g(1) * (-t223 * t237 + t227) - g(2) * t254 - g(3) * t252, 0, 0, 0, 0, 0, -g(1) * (t216 * t247 + t237 * t260) - g(2) * (t214 * t247 - t253) - g(3) * (t221 * t247 + t257) -g(1) * (-t216 * t244 + t237 * t258) - g(2) * (-t214 * t244 - t240 * t258) - g(3) * (-t221 * t244 + t241 * t247) g(1) * t215 + g(2) * t213 - g(3) * t220, -g(1) * (t215 * t242 + t216 * t229 + t227 + (pkin(4) * t260 - t223) * t237) - g(2) * (-pkin(4) * t253 + t213 * t242 + t214 * t229 + t254) - g(3) * (pkin(4) * t257 - t220 * t242 + t221 * t229 + t252) 0, 0, 0, 0, 0, -g(1) * (t212 * t246 - t215 * t243) - g(2) * (t211 * t246 - t213 * t243) - g(3) * (t217 * t246 + t220 * t243) -g(1) * (-t212 * t243 - t215 * t246) - g(2) * (-t211 * t243 - t213 * t246) - g(3) * (-t217 * t243 + t220 * t246);];
U_reg  = t1;
