% Calculate minimal parameter regressor of potential energy for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:57:08
% EndTime: 2019-03-09 09:57:08
% DurationCPUTime: 0.12s
% Computational Cost: add. (173->67), mult. (229->84), div. (0->0), fcn. (246->8), ass. (0->37)
t217 = sin(qJ(1));
t219 = cos(qJ(1));
t226 = g(1) * t219 + g(2) * t217;
t216 = sin(qJ(2));
t239 = g(3) * t216;
t238 = t219 * pkin(7);
t212 = pkin(9) + qJ(4);
t206 = sin(t212);
t236 = t206 * t216;
t207 = cos(t212);
t235 = t207 * t216;
t215 = -pkin(8) - qJ(3);
t234 = t215 * t216;
t213 = sin(pkin(9));
t233 = t217 * t213;
t218 = cos(qJ(2));
t232 = t217 * t218;
t231 = t219 * t206;
t230 = t219 * t207;
t229 = t219 * t213;
t214 = cos(pkin(9));
t228 = t219 * t214;
t227 = t219 * pkin(1) + t217 * pkin(7);
t204 = t214 * pkin(3) + pkin(2);
t225 = pkin(4) * t235 + qJ(5) * t236 + t216 * t204 + t218 * t215 + pkin(6);
t224 = pkin(2) * t218 + qJ(3) * t216;
t191 = t206 * t232 + t230;
t192 = t207 * t232 - t231;
t209 = t217 * pkin(1);
t223 = t192 * pkin(4) + t191 * qJ(5) + t204 * t232 + t209;
t193 = -t217 * t207 + t218 * t231;
t194 = t217 * t206 + t218 * t230;
t222 = t219 * t218 * t204 + pkin(3) * t233 + t194 * pkin(4) + t193 * qJ(5) + t227;
t221 = g(1) * t193 + g(2) * t191 + g(3) * t236;
t220 = g(1) * t194 + g(2) * t192 + g(3) * t235;
t195 = -g(3) * t218 + t226 * t216;
t1 = [0, -t226, g(1) * t217 - g(2) * t219, 0, 0, 0, 0, 0, -t226 * t218 - t239, t195, -g(1) * (t218 * t228 + t233) - g(2) * (t214 * t232 - t229) - t214 * t239, -g(1) * (t217 * t214 - t218 * t229) - g(2) * (-t213 * t232 - t228) + t213 * t239, -t195, -g(1) * (t224 * t219 + t227) - g(2) * (t224 * t217 + t209 - t238) - g(3) * (t216 * pkin(2) - t218 * qJ(3) + pkin(6)) 0, 0, 0, 0, 0, -t220, t221, -t195, t220, -t221, -g(1) * (-t219 * t234 + t222) - g(2) * (-t217 * t234 + (-pkin(3) * t213 - pkin(7)) * t219 + t223) - g(3) * t225, -t195, -t221, -t220, -g(1) * (t194 * qJ(6) + t222) - g(2) * (-pkin(3) * t229 + t192 * qJ(6) + t223 - t238) - g(3) * (-t218 * pkin(5) + t225) + (-g(3) * qJ(6) * t207 - t226 * (pkin(5) - t215)) * t216;];
U_reg  = t1;
