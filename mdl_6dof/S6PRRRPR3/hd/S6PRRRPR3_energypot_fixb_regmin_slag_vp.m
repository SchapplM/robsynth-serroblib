% Calculate minimal parameter regressor of potential energy for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:14:00
% EndTime: 2019-03-08 23:14:00
% DurationCPUTime: 0.13s
% Computational Cost: add. (167->64), mult. (278->103), div. (0->0), fcn. (344->12), ass. (0->38)
t214 = sin(pkin(11));
t215 = sin(pkin(6));
t237 = t214 * t215;
t216 = cos(pkin(11));
t236 = t215 * t216;
t219 = sin(qJ(3));
t235 = t215 * t219;
t220 = sin(qJ(2));
t234 = t215 * t220;
t222 = cos(qJ(3));
t233 = t215 * t222;
t223 = cos(qJ(2));
t232 = t215 * t223;
t217 = cos(pkin(6));
t231 = t217 * t219;
t230 = t217 * t220;
t229 = t217 * t223;
t204 = t214 * t223 + t216 * t230;
t213 = qJ(3) + qJ(4);
t211 = sin(t213);
t212 = cos(t213);
t197 = t204 * t211 + t212 * t236;
t206 = -t214 * t230 + t216 * t223;
t199 = t206 * t211 - t212 * t237;
t201 = t211 * t234 - t217 * t212;
t227 = g(1) * t199 + g(2) * t197 + g(3) * t201;
t198 = t204 * t212 - t211 * t236;
t200 = t206 * t212 + t211 * t237;
t202 = t217 * t211 + t212 * t234;
t226 = g(1) * t200 + g(2) * t198 + g(3) * t202;
t203 = t214 * t220 - t216 * t229;
t205 = t214 * t229 + t216 * t220;
t225 = -g(1) * t205 - g(2) * t203 + g(3) * t232;
t224 = -pkin(9) - pkin(8);
t221 = cos(qJ(6));
t218 = sin(qJ(6));
t210 = t222 * pkin(3) + pkin(2);
t1 = [-g(3) * qJ(1), 0, -g(1) * t206 - g(2) * t204 - g(3) * t234, -t225, 0, 0, 0, 0, 0, -g(1) * (t206 * t222 + t214 * t235) - g(2) * (t204 * t222 - t216 * t235) - g(3) * (t220 * t233 + t231) -g(1) * (-t206 * t219 + t214 * t233) - g(2) * (-t204 * t219 - t216 * t233) - g(3) * (t217 * t222 - t219 * t234) 0, 0, 0, 0, 0, -t226, t227, t225, t226, -t227, -g(1) * (t216 * pkin(1) + t200 * pkin(4) + t199 * qJ(5) - t205 * t224 + t206 * t210) - g(2) * (t214 * pkin(1) + t198 * pkin(4) + t197 * qJ(5) - t203 * t224 + t204 * t210) - g(3) * (pkin(3) * t231 + t202 * pkin(4) + t217 * pkin(7) + t201 * qJ(5) + qJ(1)) + (-g(3) * (t210 * t220 + t223 * t224) + (-g(1) * t214 + g(2) * t216) * (pkin(3) * t219 + pkin(7))) * t215, 0, 0, 0, 0, 0, -g(1) * (t199 * t218 + t205 * t221) - g(2) * (t197 * t218 + t203 * t221) - g(3) * (t201 * t218 - t221 * t232) -g(1) * (t199 * t221 - t205 * t218) - g(2) * (t197 * t221 - t203 * t218) - g(3) * (t201 * t221 + t218 * t232);];
U_reg  = t1;
