% Calculate minimal parameter regressor of potential energy for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:01:39
% EndTime: 2019-03-09 10:01:39
% DurationCPUTime: 0.13s
% Computational Cost: add. (117->62), mult. (181->76), div. (0->0), fcn. (184->8), ass. (0->42)
t219 = sin(qJ(2));
t246 = qJ(3) * t219 + pkin(1);
t222 = cos(qJ(2));
t245 = g(3) * t222;
t244 = t219 * pkin(2) + pkin(6);
t216 = qJ(4) + pkin(9);
t210 = sin(t216);
t220 = sin(qJ(1));
t242 = t220 * t210;
t211 = cos(t216);
t241 = t220 * t211;
t218 = sin(qJ(4));
t240 = t220 * t218;
t221 = cos(qJ(4));
t239 = t220 * t221;
t238 = t220 * t222;
t223 = cos(qJ(1));
t237 = t222 * t223;
t236 = t223 * t210;
t235 = t223 * t211;
t234 = t223 * t218;
t233 = t223 * t221;
t232 = t219 * t240;
t231 = t219 * t234;
t230 = pkin(2) * t238 + t246 * t220;
t229 = -pkin(4) * t218 - qJ(3);
t228 = pkin(2) * t237 + t220 * pkin(7) + t246 * t223;
t217 = -qJ(5) - pkin(8);
t227 = -t219 * t217 + t244;
t226 = g(1) * t223 + g(2) * t220;
t225 = -t223 * pkin(7) + t230;
t209 = t221 * pkin(4) + pkin(3);
t224 = pkin(4) * t231 + t220 * t209 + t228;
t203 = g(1) * t220 - g(2) * t223;
t201 = pkin(4) * t232;
t200 = g(3) * t219 + t226 * t222;
t199 = t226 * t219 - t245;
t198 = t219 * t242 - t235;
t197 = t219 * t241 + t236;
t196 = t219 * t236 + t241;
t195 = -t219 * t235 + t242;
t1 = [0, -t226, t203, 0, 0, 0, 0, 0, -t200, t199, -t203, t200, -t199, -g(1) * t228 - g(2) * t225 - g(3) * (-t222 * qJ(3) + t244) 0, 0, 0, 0, 0, -g(1) * (t231 + t239) - g(2) * (t232 - t233) + t218 * t245, -g(1) * (t219 * t233 - t240) - g(2) * (t219 * t239 + t234) + t221 * t245, -t200, -g(1) * (-t217 * t237 + t224) - g(2) * (-t217 * t238 + t201 + (-pkin(7) - t209) * t223 + t230) - g(3) * (t229 * t222 + t227) -g(1) * t196 - g(2) * t198 + t210 * t245, -t200, -g(1) * t195 + g(2) * t197 - t211 * t245, -g(1) * (t196 * pkin(5) + t195 * qJ(6) + t224) - g(2) * (t198 * pkin(5) - t197 * qJ(6) - t223 * t209 + t201 + t225) - g(3) * t227 + (-g(3) * (-pkin(5) * t210 + qJ(6) * t211 + t229) + t226 * t217) * t222;];
U_reg  = t1;
