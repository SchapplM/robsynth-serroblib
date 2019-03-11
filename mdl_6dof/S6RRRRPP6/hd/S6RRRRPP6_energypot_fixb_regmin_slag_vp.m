% Calculate minimal parameter regressor of potential energy for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:15:03
% EndTime: 2019-03-09 21:15:03
% DurationCPUTime: 0.11s
% Computational Cost: add. (161->58), mult. (207->74), div. (0->0), fcn. (227->8), ass. (0->33)
t218 = sin(qJ(1));
t221 = cos(qJ(1));
t228 = g(1) * t221 + g(2) * t218;
t217 = sin(qJ(2));
t239 = g(3) * t217;
t215 = qJ(3) + qJ(4);
t209 = sin(t215);
t237 = t209 * t217;
t210 = cos(t215);
t236 = t210 * t217;
t222 = -pkin(9) - pkin(8);
t235 = t217 * t222;
t216 = sin(qJ(3));
t234 = t218 * t216;
t220 = cos(qJ(2));
t233 = t218 * t220;
t232 = t221 * t209;
t231 = t221 * t210;
t230 = t221 * t216;
t219 = cos(qJ(3));
t229 = t221 * t219;
t207 = t219 * pkin(3) + pkin(2);
t227 = pkin(4) * t236 + qJ(5) * t237 + t217 * t207 + t220 * t222 + pkin(6);
t194 = t209 * t233 + t231;
t195 = t210 * t233 - t232;
t226 = t218 * pkin(1) + t195 * pkin(4) + t194 * qJ(5) + t207 * t233;
t196 = -t218 * t210 + t220 * t232;
t197 = t218 * t209 + t220 * t231;
t225 = pkin(3) * t234 + t197 * pkin(4) + t218 * pkin(7) + t196 * qJ(5) + (t207 * t220 + pkin(1)) * t221;
t224 = g(1) * t196 + g(2) * t194 + g(3) * t237;
t223 = g(1) * t197 + g(2) * t195 + g(3) * t236;
t198 = -g(3) * t220 + t228 * t217;
t1 = [0, -t228, g(1) * t218 - g(2) * t221, 0, 0, 0, 0, 0, -t228 * t220 - t239, t198, 0, 0, 0, 0, 0, -g(1) * (t220 * t229 + t234) - g(2) * (t219 * t233 - t230) - t219 * t239, -g(1) * (t218 * t219 - t220 * t230) - g(2) * (-t216 * t233 - t229) + t216 * t239, 0, 0, 0, 0, 0, -t223, t224, -t198, t223, -t224, -g(1) * (-t221 * t235 + t225) - g(2) * (-t218 * t235 + (-pkin(3) * t216 - pkin(7)) * t221 + t226) - g(3) * t227, -t198, -t224, -t223, -g(1) * (t197 * qJ(6) + t225) - g(2) * (-pkin(3) * t230 - t221 * pkin(7) + t195 * qJ(6) + t226) - g(3) * (-t220 * pkin(5) + t227) + (-g(3) * qJ(6) * t210 - t228 * (pkin(5) - t222)) * t217;];
U_reg  = t1;
