% Calculate minimal parameter regressor of potential energy for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:54:06
% EndTime: 2019-03-08 21:54:06
% DurationCPUTime: 0.15s
% Computational Cost: add. (137->56), mult. (206->97), div. (0->0), fcn. (252->12), ass. (0->33)
t227 = sin(pkin(11));
t228 = sin(pkin(6));
t248 = t227 * t228;
t229 = cos(pkin(11));
t247 = t228 * t229;
t233 = sin(qJ(3));
t246 = t228 * t233;
t234 = sin(qJ(2));
t245 = t228 * t234;
t236 = cos(qJ(3));
t244 = t228 * t236;
t237 = cos(qJ(2));
t243 = t228 * t237;
t230 = cos(pkin(6));
t242 = t230 * t233;
t241 = t230 * t234;
t240 = t230 * t237;
t218 = t227 * t234 - t229 * t240;
t220 = t227 * t240 + t229 * t234;
t238 = -g(1) * t220 - g(2) * t218 + g(3) * t243;
t235 = cos(qJ(6));
t232 = sin(qJ(6));
t231 = -qJ(4) - pkin(8);
t226 = qJ(3) + pkin(12) + qJ(5);
t225 = t236 * pkin(3) + pkin(2);
t224 = cos(t226);
t223 = sin(t226);
t221 = -t227 * t241 + t229 * t237;
t219 = t227 * t237 + t229 * t241;
t217 = t230 * t223 + t224 * t245;
t216 = t221 * t224 + t223 * t248;
t215 = t219 * t224 - t223 * t247;
t1 = [-g(3) * qJ(1), 0, -g(1) * t221 - g(2) * t219 - g(3) * t245, -t238, 0, 0, 0, 0, 0, -g(1) * (t221 * t236 + t227 * t246) - g(2) * (t219 * t236 - t229 * t246) - g(3) * (t234 * t244 + t242) -g(1) * (-t221 * t233 + t227 * t244) - g(2) * (-t219 * t233 - t229 * t244) - g(3) * (t230 * t236 - t233 * t245) t238, -g(1) * (t229 * pkin(1) - t220 * t231 + t221 * t225) - g(2) * (t227 * pkin(1) - t218 * t231 + t219 * t225) - g(3) * (pkin(3) * t242 + t230 * pkin(7) + qJ(1)) + (-g(3) * (t225 * t234 + t231 * t237) + (-g(1) * t227 + g(2) * t229) * (pkin(3) * t233 + pkin(7))) * t228, 0, 0, 0, 0, 0, -g(1) * t216 - g(2) * t215 - g(3) * t217, -g(1) * (-t221 * t223 + t224 * t248) - g(2) * (-t219 * t223 - t224 * t247) - g(3) * (-t223 * t245 + t230 * t224) 0, 0, 0, 0, 0, -g(1) * (t216 * t235 + t220 * t232) - g(2) * (t215 * t235 + t218 * t232) - g(3) * (t217 * t235 - t232 * t243) -g(1) * (-t216 * t232 + t220 * t235) - g(2) * (-t215 * t232 + t218 * t235) - g(3) * (-t217 * t232 - t235 * t243);];
U_reg  = t1;
