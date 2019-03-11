% Calculate minimal parameter regressor of potential energy for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:17:45
% EndTime: 2019-03-09 18:17:45
% DurationCPUTime: 0.09s
% Computational Cost: add. (114->46), mult. (109->67), div. (0->0), fcn. (118->12), ass. (0->33)
t224 = qJ(2) + qJ(3);
t221 = sin(t224);
t246 = g(3) * t221;
t223 = pkin(11) + qJ(5);
t220 = qJ(6) + t223;
t215 = sin(t220);
t228 = sin(qJ(1));
t245 = t228 * t215;
t216 = cos(t220);
t244 = t228 * t216;
t218 = sin(t223);
t243 = t228 * t218;
t219 = cos(t223);
t242 = t228 * t219;
t225 = sin(pkin(11));
t241 = t228 * t225;
t226 = cos(pkin(11));
t240 = t228 * t226;
t230 = cos(qJ(1));
t239 = t230 * t215;
t238 = t230 * t216;
t237 = t230 * t218;
t236 = t230 * t219;
t235 = t230 * t225;
t234 = t230 * t226;
t233 = g(1) * t230 + g(2) * t228;
t222 = cos(t224);
t229 = cos(qJ(2));
t232 = t229 * pkin(2) + pkin(3) * t222 + qJ(4) * t221 + pkin(1);
t231 = -pkin(8) - pkin(7);
t227 = sin(qJ(2));
t214 = -g(3) * t222 + t233 * t221;
t1 = [0, -t233, g(1) * t228 - g(2) * t230, 0, 0, 0, 0, 0, -g(3) * t227 - t233 * t229, -g(3) * t229 + t233 * t227, 0, 0, 0, 0, 0, -t233 * t222 - t246, t214, -g(1) * (t222 * t234 + t241) - g(2) * (t222 * t240 - t235) - t226 * t246, -g(1) * (-t222 * t235 + t240) - g(2) * (-t222 * t241 - t234) + t225 * t246, -t214, -g(3) * (t227 * pkin(2) + t221 * pkin(3) - t222 * qJ(4) + pkin(6)) + (-g(1) * t232 - g(2) * t231) * t230 + (g(1) * t231 - g(2) * t232) * t228, 0, 0, 0, 0, 0, -g(1) * (t222 * t236 + t243) - g(2) * (t222 * t242 - t237) - t219 * t246, -g(1) * (-t222 * t237 + t242) - g(2) * (-t222 * t243 - t236) + t218 * t246, 0, 0, 0, 0, 0, -g(1) * (t222 * t238 + t245) - g(2) * (t222 * t244 - t239) - t216 * t246, -g(1) * (-t222 * t239 + t244) - g(2) * (-t222 * t245 - t238) + t215 * t246;];
U_reg  = t1;
