% Calculate minimal parameter regressor of potential energy for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:39:34
% EndTime: 2019-03-09 02:39:34
% DurationCPUTime: 0.13s
% Computational Cost: add. (125->44), mult. (100->66), div. (0->0), fcn. (99->12), ass. (0->33)
t208 = qJ(3) + pkin(10);
t199 = sin(t208);
t231 = g(3) * t199;
t230 = qJ(2) + pkin(6);
t209 = qJ(1) + pkin(9);
t200 = sin(t209);
t202 = cos(t208);
t229 = t200 * t202;
t210 = sin(pkin(11));
t228 = t200 * t210;
t211 = cos(pkin(11));
t227 = t200 * t211;
t207 = pkin(11) + qJ(6);
t198 = sin(t207);
t203 = cos(t209);
t226 = t203 * t198;
t201 = cos(t207);
t225 = t203 * t201;
t224 = t203 * t210;
t223 = t203 * t211;
t213 = sin(qJ(3));
t222 = t213 * pkin(3) + t230;
t215 = cos(qJ(3));
t197 = t215 * pkin(3) + pkin(2);
t212 = -qJ(4) - pkin(7);
t214 = sin(qJ(1));
t221 = t214 * pkin(1) + t200 * t197 + t203 * t212;
t220 = g(1) * t203 + g(2) * t200;
t216 = cos(qJ(1));
t219 = -g(1) * t216 - g(2) * t214;
t218 = t216 * pkin(1) + t203 * t197 - t200 * t212;
t217 = pkin(4) * t202 + qJ(5) * t199;
t1 = [0, t219, g(1) * t214 - g(2) * t216, pkin(1) * t219 - g(3) * t230, 0, 0, 0, 0, 0, -g(3) * t213 - t215 * t220, -g(3) * t215 + t213 * t220, -g(1) * t200 + g(2) * t203, -g(1) * t218 - g(2) * t221 - g(3) * t222, -g(1) * (t202 * t223 + t228) - g(2) * (t202 * t227 - t224) - t211 * t231, -g(1) * (-t202 * t224 + t227) - g(2) * (-t202 * t228 - t223) + t210 * t231, g(3) * t202 - t199 * t220, -g(1) * (t203 * t217 + t218) - g(2) * (t200 * t217 + t221) - g(3) * (t199 * pkin(4) - t202 * qJ(5) + t222) 0, 0, 0, 0, 0, -g(1) * (t200 * t198 + t202 * t225) - g(2) * (t201 * t229 - t226) - t201 * t231, -g(1) * (t200 * t201 - t202 * t226) - g(2) * (-t198 * t229 - t225) + t198 * t231;];
U_reg  = t1;
