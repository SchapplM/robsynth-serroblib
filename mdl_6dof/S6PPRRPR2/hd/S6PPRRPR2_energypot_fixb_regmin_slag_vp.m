% Calculate minimal parameter regressor of potential energy for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRRPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:50:54
% EndTime: 2019-03-08 18:50:54
% DurationCPUTime: 0.17s
% Computational Cost: add. (274->72), mult. (736->120), div. (0->0), fcn. (953->14), ass. (0->50)
t243 = sin(pkin(12));
t246 = sin(pkin(6));
t275 = t243 * t246;
t244 = sin(pkin(11));
t274 = t244 * t246;
t250 = cos(pkin(6));
t273 = t244 * t250;
t245 = sin(pkin(7));
t256 = cos(qJ(3));
t272 = t245 * t256;
t247 = cos(pkin(12));
t271 = t246 * t247;
t248 = cos(pkin(11));
t270 = t246 * t248;
t249 = cos(pkin(7));
t269 = t246 * t249;
t268 = t248 * t250;
t267 = t249 * t256;
t266 = qJ(2) * t250 + qJ(1);
t265 = pkin(1) * t248 + qJ(2) * t274;
t264 = t246 * t272;
t263 = pkin(1) * t244 - qJ(2) * t270;
t232 = -t243 * t244 + t247 * t268;
t262 = t232 * t245 + t248 * t269;
t234 = -t243 * t248 - t247 * t273;
t261 = -t234 * t245 + t244 * t269;
t260 = t245 * t271 - t249 * t250;
t233 = t243 * t268 + t244 * t247;
t253 = sin(qJ(3));
t220 = t233 * t256 + (t232 * t249 - t245 * t270) * t253;
t252 = sin(qJ(4));
t255 = cos(qJ(4));
t215 = t220 * t252 + t255 * t262;
t235 = -t243 * t273 + t247 * t248;
t222 = t235 * t256 + (t234 * t249 + t245 * t274) * t253;
t217 = t222 * t252 - t255 * t261;
t228 = t250 * t245 * t253 + (t247 * t249 * t253 + t243 * t256) * t246;
t223 = t228 * t252 + t255 * t260;
t259 = g(1) * t217 + g(2) * t215 + g(3) * t223;
t216 = t220 * t255 - t252 * t262;
t218 = t222 * t255 + t252 * t261;
t224 = t228 * t255 - t252 * t260;
t258 = g(1) * t218 + g(2) * t216 + g(3) * t224;
t219 = -t232 * t267 + t233 * t253 + t248 * t264;
t221 = -t234 * t267 + t235 * t253 - t244 * t264;
t227 = -t250 * t272 + t253 * t275 - t267 * t271;
t257 = g(1) * t221 + g(2) * t219 + g(3) * t227;
t254 = cos(qJ(6));
t251 = sin(qJ(6));
t1 = [-g(3) * qJ(1), -g(1) * t265 - g(2) * t263 - g(3) * t266, 0, -g(1) * t222 - g(2) * t220 - g(3) * t228, t257, 0, 0, 0, 0, 0, -t258, t259, -t257, t258, -t259, -g(1) * (t235 * pkin(2) + t222 * pkin(3) + t218 * pkin(4) + t221 * pkin(9) + t217 * qJ(5) + t265) - g(2) * (t233 * pkin(2) + t220 * pkin(3) + t216 * pkin(4) + t219 * pkin(9) + t215 * qJ(5) + t263) - g(3) * (pkin(2) * t275 + t228 * pkin(3) + t224 * pkin(4) + t227 * pkin(9) + t223 * qJ(5) + t266) + (-g(1) * t261 + g(2) * t262 + g(3) * t260) * pkin(8), 0, 0, 0, 0, 0, -g(1) * (t217 * t251 + t221 * t254) - g(2) * (t215 * t251 + t219 * t254) - g(3) * (t223 * t251 + t227 * t254) -g(1) * (t217 * t254 - t221 * t251) - g(2) * (t215 * t254 - t219 * t251) - g(3) * (t223 * t254 - t227 * t251);];
U_reg  = t1;
