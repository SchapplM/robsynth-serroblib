% Calculate minimal parameter regressor of potential energy for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRPRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:43:37
% EndTime: 2019-03-08 18:43:38
% DurationCPUTime: 0.21s
% Computational Cost: add. (223->72), mult. (579->132), div. (0->0), fcn. (748->16), ass. (0->55)
t248 = sin(qJ(3));
t264 = pkin(3) * t248;
t263 = pkin(8) + qJ(4);
t238 = sin(pkin(11));
t240 = sin(pkin(6));
t262 = t238 * t240;
t245 = cos(pkin(6));
t261 = t238 * t245;
t239 = sin(pkin(7));
t260 = t239 * t240;
t243 = cos(pkin(11));
t259 = t240 * t243;
t244 = cos(pkin(7));
t258 = t240 * t244;
t257 = t243 * t245;
t256 = t245 * qJ(2) + qJ(1);
t255 = t243 * pkin(1) + qJ(2) * t262;
t236 = sin(pkin(13));
t241 = cos(pkin(13));
t251 = cos(qJ(3));
t254 = t251 * t236 + t248 * t241;
t230 = -t248 * t236 + t251 * t241;
t237 = sin(pkin(12));
t242 = cos(pkin(12));
t226 = t237 * t257 + t238 * t242;
t228 = -t237 * t261 + t243 * t242;
t253 = g(3) * t237 * t240 + g(1) * t228 + g(2) * t226;
t225 = -t238 * t237 + t242 * t257;
t227 = -t243 * t237 - t242 * t261;
t252 = -g(1) * (t227 * t244 + t238 * t260) - g(2) * (t225 * t244 - t239 * t259) - g(3) * (t239 * t245 + t242 * t258);
t250 = cos(qJ(5));
t249 = cos(qJ(6));
t247 = sin(qJ(5));
t246 = sin(qJ(6));
t234 = t238 * pkin(1);
t232 = t251 * pkin(3) + pkin(2);
t224 = -t242 * t260 + t245 * t244;
t223 = -t263 * t239 + t244 * t264;
t222 = t239 * t264 + t263 * t244;
t221 = t254 * t244;
t220 = t230 * t244;
t219 = t254 * t239;
t218 = t230 * t239;
t217 = -t227 * t239 + t238 * t258;
t216 = -t225 * t239 - t243 * t258;
t215 = t245 * t219 + (t221 * t242 + t230 * t237) * t240;
t214 = -t245 * t218 + (-t220 * t242 + t237 * t254) * t240;
t213 = t215 * t250 + t224 * t247;
t212 = t219 * t262 + t227 * t221 + t228 * t230;
t211 = -t218 * t262 - t227 * t220 + t228 * t254;
t210 = -t219 * t259 + t225 * t221 + t226 * t230;
t209 = t218 * t259 - t225 * t220 + t226 * t254;
t208 = t212 * t250 + t217 * t247;
t207 = t210 * t250 + t216 * t247;
t1 = [-g(3) * qJ(1), -g(1) * t255 - g(2) * (-qJ(2) * t259 + t234) - g(3) * t256, 0, t252 * t248 - t253 * t251, t253 * t248 + t252 * t251, -g(1) * (t227 * t223 + t228 * t232 + t255) - g(2) * (t225 * t223 + t226 * t232 + t234) - g(3) * (t245 * t222 + t256) + (-g(1) * t222 * t238 - g(3) * (t223 * t242 + t232 * t237) - g(2) * (-qJ(2) - t222) * t243) * t240, 0, 0, 0, 0, 0, -g(1) * t208 - g(2) * t207 - g(3) * t213, -g(1) * (-t212 * t247 + t217 * t250) - g(2) * (-t210 * t247 + t216 * t250) - g(3) * (-t215 * t247 + t224 * t250) 0, 0, 0, 0, 0, -g(1) * (t208 * t249 + t211 * t246) - g(2) * (t207 * t249 + t209 * t246) - g(3) * (t213 * t249 + t214 * t246) -g(1) * (-t208 * t246 + t211 * t249) - g(2) * (-t207 * t246 + t209 * t249) - g(3) * (-t213 * t246 + t214 * t249);];
U_reg  = t1;
