% Calculate minimal parameter regressor of potential energy for
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:05:29
% EndTime: 2019-03-09 14:05:29
% DurationCPUTime: 0.13s
% Computational Cost: add. (110->49), mult. (116->76), div. (0->0), fcn. (129->12), ass. (0->29)
t231 = sin(qJ(2));
t246 = g(3) * t231;
t232 = sin(qJ(1));
t233 = cos(qJ(2));
t245 = t232 * t233;
t228 = pkin(11) + qJ(4);
t227 = qJ(5) + t228;
t224 = qJ(6) + t227;
t220 = sin(t224);
t234 = cos(qJ(1));
t244 = t234 * t220;
t221 = cos(t224);
t243 = t234 * t221;
t222 = sin(t227);
t242 = t234 * t222;
t223 = cos(t227);
t241 = t234 * t223;
t225 = sin(t228);
t240 = t234 * t225;
t226 = cos(t228);
t239 = t234 * t226;
t229 = sin(pkin(11));
t238 = t234 * t229;
t230 = cos(pkin(11));
t237 = t234 * t230;
t236 = g(1) * t234 + g(2) * t232;
t235 = pkin(2) * t233 + qJ(3) * t231 + pkin(1);
t219 = -g(3) * t233 + t236 * t231;
t1 = [0, -t236, g(1) * t232 - g(2) * t234, 0, 0, 0, 0, 0, -t236 * t233 - t246, t219, -g(1) * (t232 * t229 + t233 * t237) - g(2) * (t230 * t245 - t238) - t230 * t246, -g(1) * (t232 * t230 - t233 * t238) - g(2) * (-t229 * t245 - t237) + t229 * t246, -t219, -g(3) * (t231 * pkin(2) - t233 * qJ(3) + pkin(6)) + (g(2) * pkin(7) - g(1) * t235) * t234 + (-g(1) * pkin(7) - g(2) * t235) * t232, 0, 0, 0, 0, 0, -g(1) * (t232 * t225 + t233 * t239) - g(2) * (t226 * t245 - t240) - t226 * t246, -g(1) * (t232 * t226 - t233 * t240) - g(2) * (-t225 * t245 - t239) + t225 * t246, 0, 0, 0, 0, 0, -g(1) * (t232 * t222 + t233 * t241) - g(2) * (t223 * t245 - t242) - t223 * t246, -g(1) * (t232 * t223 - t233 * t242) - g(2) * (-t222 * t245 - t241) + t222 * t246, 0, 0, 0, 0, 0, -g(1) * (t232 * t220 + t233 * t243) - g(2) * (t221 * t245 - t244) - t221 * t246, -g(1) * (t232 * t221 - t233 * t244) - g(2) * (-t220 * t245 - t243) + t220 * t246;];
U_reg  = t1;
