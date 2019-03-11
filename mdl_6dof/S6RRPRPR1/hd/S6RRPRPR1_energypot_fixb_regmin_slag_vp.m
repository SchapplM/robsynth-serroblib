% Calculate minimal parameter regressor of potential energy for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:09:56
% EndTime: 2019-03-09 10:09:56
% DurationCPUTime: 0.08s
% Computational Cost: add. (125->46), mult. (104->62), div. (0->0), fcn. (106->12), ass. (0->31)
t227 = qJ(2) + pkin(10);
t222 = qJ(4) + t227;
t217 = sin(t222);
t246 = g(3) * t217;
t230 = -pkin(7) - qJ(3);
t231 = sin(qJ(2));
t245 = t231 * pkin(2) + pkin(6);
t233 = cos(qJ(2));
t219 = t233 * pkin(2) + pkin(1);
t226 = pkin(11) + qJ(6);
t220 = sin(t226);
t232 = sin(qJ(1));
t244 = t232 * t220;
t221 = cos(t226);
t243 = t232 * t221;
t228 = sin(pkin(11));
t242 = t232 * t228;
t229 = cos(pkin(11));
t241 = t232 * t229;
t234 = cos(qJ(1));
t240 = t234 * t220;
t239 = t234 * t221;
t238 = t234 * t228;
t237 = t234 * t229;
t236 = g(1) * t234 + g(2) * t232;
t218 = cos(t222);
t235 = pkin(4) * t218 + qJ(5) * t217 + pkin(3) * cos(t227) + t219;
t225 = -pkin(8) + t230;
t216 = g(1) * t232 - g(2) * t234;
t214 = -g(3) * t218 + t236 * t217;
t1 = [0, -t236, t216, 0, 0, 0, 0, 0, -g(3) * t231 - t236 * t233, -g(3) * t233 + t236 * t231, -t216, -g(1) * (t234 * t219 - t232 * t230) - g(2) * (t232 * t219 + t234 * t230) - g(3) * t245, 0, 0, 0, 0, 0, -t236 * t218 - t246, t214, -g(1) * (t218 * t237 + t242) - g(2) * (t218 * t241 - t238) - t229 * t246, -g(1) * (-t218 * t238 + t241) - g(2) * (-t218 * t242 - t237) + t228 * t246, -t214, -g(3) * (t217 * pkin(4) - t218 * qJ(5) + pkin(3) * sin(t227) + t245) + (-g(1) * t235 - g(2) * t225) * t234 + (g(1) * t225 - g(2) * t235) * t232, 0, 0, 0, 0, 0, -g(1) * (t218 * t239 + t244) - g(2) * (t218 * t243 - t240) - t221 * t246, -g(1) * (-t218 * t240 + t243) - g(2) * (-t218 * t244 - t239) + t220 * t246;];
U_reg  = t1;
