% Calculate minimal parameter regressor of potential energy for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:16:25
% EndTime: 2019-03-08 21:16:26
% DurationCPUTime: 0.16s
% Computational Cost: add. (223->67), mult. (541->105), div. (0->0), fcn. (692->12), ass. (0->42)
t242 = sin(pkin(6));
t263 = pkin(7) * t242;
t247 = sin(qJ(3));
t262 = t242 * t247;
t248 = sin(qJ(2));
t261 = t242 * t248;
t250 = cos(qJ(3));
t260 = t242 * t250;
t251 = cos(qJ(2));
t259 = t242 * t251;
t245 = cos(pkin(6));
t258 = t245 * t248;
t257 = t245 * t251;
t241 = sin(pkin(10));
t244 = cos(pkin(10));
t227 = t241 * t251 + t244 * t258;
t218 = t227 * t250 - t244 * t262;
t226 = t241 * t248 - t244 * t257;
t240 = sin(pkin(11));
t243 = cos(pkin(11));
t209 = t218 * t240 - t226 * t243;
t229 = -t241 * t258 + t244 * t251;
t220 = t229 * t250 + t241 * t262;
t228 = t241 * t257 + t244 * t248;
t211 = t220 * t240 - t228 * t243;
t231 = t245 * t247 + t248 * t260;
t215 = t231 * t240 + t243 * t259;
t256 = g(1) * t211 + g(2) * t209 + g(3) * t215;
t217 = t227 * t247 + t244 * t260;
t219 = t229 * t247 - t241 * t260;
t230 = -t245 * t250 + t247 * t261;
t255 = g(1) * t219 + g(2) * t217 + g(3) * t230;
t254 = t244 * pkin(1) + t229 * pkin(2) + t220 * pkin(3) + t228 * pkin(8) + t219 * qJ(4) + t241 * t263;
t253 = pkin(2) * t261 + t231 * pkin(3) + t245 * pkin(7) - pkin(8) * t259 + t230 * qJ(4) + qJ(1);
t252 = t241 * pkin(1) + t227 * pkin(2) + t218 * pkin(3) + t226 * pkin(8) + t217 * qJ(4) - t244 * t263;
t249 = cos(qJ(6));
t246 = sin(qJ(6));
t216 = t231 * t243 - t240 * t259;
t212 = t220 * t243 + t228 * t240;
t210 = t218 * t243 + t226 * t240;
t207 = -g(1) * t212 - g(2) * t210 - g(3) * t216;
t1 = [-g(3) * qJ(1), 0, -g(1) * t229 - g(2) * t227 - g(3) * t261, g(1) * t228 + g(2) * t226 - g(3) * t259, 0, 0, 0, 0, 0, -g(1) * t220 - g(2) * t218 - g(3) * t231, t255, t207, t256, -t255, -g(1) * t254 - g(2) * t252 - g(3) * t253, t207, -t255, -t256, -g(1) * (t212 * pkin(4) + t211 * qJ(5) + t254) - g(2) * (t210 * pkin(4) + t209 * qJ(5) + t252) - g(3) * (t216 * pkin(4) + t215 * qJ(5) + t253) 0, 0, 0, 0, 0, -g(1) * (t211 * t246 + t212 * t249) - g(2) * (t209 * t246 + t210 * t249) - g(3) * (t215 * t246 + t216 * t249) -g(1) * (t211 * t249 - t212 * t246) - g(2) * (t209 * t249 - t210 * t246) - g(3) * (t215 * t249 - t216 * t246);];
U_reg  = t1;
