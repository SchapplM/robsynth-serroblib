% Calculate minimal parameter regressor of potential energy for
% S6PRRPRP3
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
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:38:20
% EndTime: 2019-03-08 21:38:20
% DurationCPUTime: 0.18s
% Computational Cost: add. (244->80), mult. (489->116), div. (0->0), fcn. (616->12), ass. (0->46)
t244 = sin(pkin(6));
t265 = pkin(7) * t244;
t249 = sin(qJ(3));
t264 = t244 * t249;
t250 = sin(qJ(2));
t263 = t244 * t250;
t251 = cos(qJ(3));
t262 = t244 * t251;
t252 = cos(qJ(2));
t261 = t244 * t252;
t247 = cos(pkin(6));
t260 = t247 * t250;
t259 = t247 * t252;
t258 = pkin(2) * t263 + t247 * pkin(7) + qJ(1);
t243 = sin(pkin(10));
t246 = cos(pkin(10));
t227 = -t243 * t260 + t246 * t252;
t257 = t246 * pkin(1) + t227 * pkin(2) + t243 * t265;
t242 = sin(pkin(11));
t256 = pkin(4) * t242 + pkin(8);
t225 = t243 * t252 + t246 * t260;
t255 = t243 * pkin(1) + t225 * pkin(2) - t246 * t265;
t217 = t225 * t251 - t246 * t264;
t224 = t243 * t250 - t246 * t259;
t241 = pkin(11) + qJ(5);
t236 = sin(t241);
t237 = cos(t241);
t210 = t217 * t236 - t224 * t237;
t219 = t227 * t251 + t243 * t264;
t226 = t243 * t259 + t246 * t250;
t212 = t219 * t236 - t226 * t237;
t229 = t247 * t249 + t250 * t262;
t214 = t229 * t236 + t237 * t261;
t254 = g(1) * t212 + g(2) * t210 + g(3) * t214;
t216 = t225 * t249 + t246 * t262;
t218 = t227 * t249 - t243 * t262;
t228 = -t247 * t251 + t249 * t263;
t253 = g(1) * t218 + g(2) * t216 + g(3) * t228;
t248 = -pkin(9) - qJ(4);
t245 = cos(pkin(11));
t235 = t245 * pkin(4) + pkin(3);
t215 = t229 * t237 - t236 * t261;
t213 = t219 * t237 + t226 * t236;
t211 = t217 * t237 + t224 * t236;
t208 = -g(1) * t213 - g(2) * t211 - g(3) * t215;
t1 = [-g(3) * qJ(1), 0, -g(1) * t227 - g(2) * t225 - g(3) * t263, g(1) * t226 + g(2) * t224 - g(3) * t261, 0, 0, 0, 0, 0, -g(1) * t219 - g(2) * t217 - g(3) * t229, t253, -g(1) * (t219 * t245 + t226 * t242) - g(2) * (t217 * t245 + t224 * t242) - g(3) * (t229 * t245 - t242 * t261) -g(1) * (-t219 * t242 + t226 * t245) - g(2) * (-t217 * t242 + t224 * t245) - g(3) * (-t229 * t242 - t245 * t261) -t253, -g(1) * (t219 * pkin(3) + t226 * pkin(8) + t218 * qJ(4) + t257) - g(2) * (t217 * pkin(3) + t224 * pkin(8) + t216 * qJ(4) + t255) - g(3) * (t229 * pkin(3) - pkin(8) * t261 + t228 * qJ(4) + t258) 0, 0, 0, 0, 0, t208, t254, t208, -t253, -t254, -g(1) * (t213 * pkin(5) + t212 * qJ(6) - t218 * t248 + t219 * t235 + t256 * t226 + t257) - g(2) * (t211 * pkin(5) + t210 * qJ(6) - t216 * t248 + t217 * t235 + t256 * t224 + t255) - g(3) * (t215 * pkin(5) + t214 * qJ(6) - t228 * t248 + t229 * t235 - t256 * t261 + t258);];
U_reg  = t1;
