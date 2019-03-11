% Calculate minimal parameter regressor of potential energy for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:49:33
% EndTime: 2019-03-08 19:49:33
% DurationCPUTime: 0.16s
% Computational Cost: add. (152->70), mult. (330->110), div. (0->0), fcn. (406->12), ass. (0->37)
t219 = sin(pkin(6));
t239 = pkin(7) * t219;
t223 = sin(qJ(4));
t238 = t219 * t223;
t224 = sin(qJ(2));
t237 = t219 * t224;
t225 = cos(qJ(4));
t236 = t219 * t225;
t226 = cos(qJ(2));
t235 = t219 * t226;
t222 = cos(pkin(6));
t234 = t222 * t224;
t233 = t222 * t226;
t232 = pkin(2) * t237 + pkin(7) * t222 + qJ(1);
t218 = sin(pkin(10));
t221 = cos(pkin(10));
t201 = t218 * t224 - t221 * t233;
t202 = t218 * t226 + t221 * t234;
t231 = pkin(1) * t218 + pkin(2) * t202 + t201 * qJ(3);
t203 = t218 * t233 + t221 * t224;
t204 = -t218 * t234 + t221 * t226;
t230 = pkin(1) * t221 + pkin(2) * t204 + qJ(3) * t203 + t218 * t239;
t194 = -t203 * t225 + t218 * t238;
t196 = t201 * t225 + t221 * t238;
t205 = t222 * t223 + t225 * t235;
t229 = g(1) * t194 - g(2) * t196 + g(3) * t205;
t228 = -g(1) * t203 - g(2) * t201 + g(3) * t235;
t227 = g(1) * t204 + g(2) * t202 + g(3) * t237;
t220 = cos(pkin(11));
t217 = sin(pkin(11));
t216 = pkin(11) + qJ(6);
t212 = cos(t216);
t211 = sin(t216);
t206 = t222 * t225 - t223 * t235;
t197 = t201 * t223 - t221 * t236;
t195 = t203 * t223 + t218 * t236;
t1 = [-g(3) * qJ(1), 0, -t227, -t228, t227, t228, -g(1) * t230 - g(2) * (-t221 * t239 + t231) - g(3) * (-qJ(3) * t235 + t232) 0, 0, 0, 0, 0, -g(1) * t195 - g(2) * t197 - g(3) * t206, t229, -g(1) * (t195 * t220 + t204 * t217) - g(2) * (t197 * t220 + t202 * t217) - g(3) * (t206 * t220 + t217 * t237) -g(1) * (-t195 * t217 + t204 * t220) - g(2) * (-t197 * t217 + t202 * t220) - g(3) * (-t206 * t217 + t220 * t237) -t229, -g(1) * (pkin(4) * t195 + pkin(8) * t204 + qJ(5) * t194 + t230) - g(2) * (t197 * pkin(4) + t202 * pkin(8) - t196 * qJ(5) + t231) - g(3) * (t222 * pkin(3) + t206 * pkin(4) + t205 * qJ(5) + t232) + (-g(1) * pkin(3) * t218 - g(3) * (pkin(8) * t224 - qJ(3) * t226) - g(2) * (-pkin(3) - pkin(7)) * t221) * t219, 0, 0, 0, 0, 0, -g(1) * (t195 * t212 + t204 * t211) - g(2) * (t197 * t212 + t202 * t211) - g(3) * (t206 * t212 + t211 * t237) -g(1) * (-t195 * t211 + t204 * t212) - g(2) * (-t197 * t211 + t202 * t212) - g(3) * (-t206 * t211 + t212 * t237);];
U_reg  = t1;
