% Calculate minimal parameter regressor of potential energy for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:19:43
% EndTime: 2019-03-09 08:19:43
% DurationCPUTime: 0.12s
% Computational Cost: add. (99->56), mult. (207->72), div. (0->0), fcn. (224->8), ass. (0->34)
t217 = sin(qJ(2));
t238 = qJ(3) * t217 + pkin(1);
t220 = cos(qJ(2));
t237 = g(3) * t220;
t236 = t217 * pkin(2) + pkin(6);
t214 = sin(pkin(9));
t218 = sin(qJ(1));
t234 = t218 * t214;
t215 = cos(pkin(9));
t233 = t218 * t215;
t232 = t218 * t220;
t221 = cos(qJ(1));
t231 = t220 * t221;
t230 = t221 * t214;
t229 = t221 * t215;
t228 = pkin(2) * t232 + t238 * t218;
t227 = pkin(2) * t231 + t218 * pkin(7) + t238 * t221;
t226 = -t220 * qJ(3) + t236;
t225 = g(1) * t221 + g(2) * t218;
t224 = t218 * pkin(3) + qJ(4) * t231 + t227;
t223 = qJ(4) * t232 + (-pkin(3) - pkin(7)) * t221 + t228;
t196 = -t217 * t229 + t234;
t198 = t217 * t233 + t230;
t222 = g(1) * t196 - g(2) * t198 + t215 * t237;
t219 = cos(qJ(6));
t216 = sin(qJ(6));
t208 = t217 * qJ(4);
t200 = g(1) * t218 - g(2) * t221;
t199 = t217 * t234 - t229;
t197 = t217 * t230 + t233;
t195 = g(3) * t217 + t225 * t220;
t194 = t225 * t217 - t237;
t193 = -g(1) * t197 - g(2) * t199 + t214 * t237;
t1 = [0, -t225, t200, 0, 0, 0, 0, 0, -t195, t194, -t200, t195, -t194, -g(1) * t227 - g(2) * (-t221 * pkin(7) + t228) - g(3) * t226, t193, t222, -t195, -g(1) * t224 - g(2) * t223 - g(3) * (t208 + t226) t193, -t195, -t222, -g(1) * (t197 * pkin(4) + t196 * qJ(5) + t224) - g(2) * (t199 * pkin(4) - t198 * qJ(5) + t223) - g(3) * (t208 + t236) - (-pkin(4) * t214 + qJ(5) * t215 - qJ(3)) * t237, 0, 0, 0, 0, 0, -g(1) * (t196 * t216 + t197 * t219) - g(2) * (-t198 * t216 + t199 * t219) - (-t214 * t219 + t215 * t216) * t237, -g(1) * (t196 * t219 - t197 * t216) - g(2) * (-t198 * t219 - t199 * t216) - (t214 * t216 + t215 * t219) * t237;];
U_reg  = t1;
