% Calculate minimal parameter regressor of potential energy for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x38]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:42:37
% EndTime: 2019-03-10 03:42:37
% DurationCPUTime: 0.07s
% Computational Cost: add. (88->34), mult. (84->54), div. (0->0), fcn. (96->12), ass. (0->30)
t224 = qJ(2) + qJ(3);
t219 = sin(t224);
t244 = g(3) * t219;
t223 = qJ(4) + qJ(5);
t222 = qJ(6) + t223;
t216 = sin(t222);
t227 = sin(qJ(1));
t243 = t227 * t216;
t217 = cos(t222);
t242 = t227 * t217;
t218 = sin(t223);
t241 = t227 * t218;
t220 = cos(t223);
t240 = t227 * t220;
t225 = sin(qJ(4));
t239 = t227 * t225;
t228 = cos(qJ(4));
t238 = t227 * t228;
t230 = cos(qJ(1));
t237 = t230 * t216;
t236 = t230 * t217;
t235 = t230 * t218;
t234 = t230 * t220;
t233 = t230 * t225;
t232 = t230 * t228;
t231 = g(1) * t230 + g(2) * t227;
t229 = cos(qJ(2));
t226 = sin(qJ(2));
t221 = cos(t224);
t1 = [0, -t231, g(1) * t227 - g(2) * t230, 0, 0, 0, 0, 0, -g(3) * t226 - t231 * t229, -g(3) * t229 + t231 * t226, 0, 0, 0, 0, 0, -t231 * t221 - t244, -g(3) * t221 + t231 * t219, 0, 0, 0, 0, 0, -g(1) * (t221 * t232 + t239) - g(2) * (t221 * t238 - t233) - t228 * t244, -g(1) * (-t221 * t233 + t238) - g(2) * (-t221 * t239 - t232) + t225 * t244, 0, 0, 0, 0, 0, -g(1) * (t221 * t234 + t241) - g(2) * (t221 * t240 - t235) - t220 * t244, -g(1) * (-t221 * t235 + t240) - g(2) * (-t221 * t241 - t234) + t218 * t244, 0, 0, 0, 0, 0, -g(1) * (t221 * t236 + t243) - g(2) * (t221 * t242 - t237) - t217 * t244, -g(1) * (-t221 * t237 + t242) - g(2) * (-t221 * t243 - t236) + t216 * t244;];
U_reg  = t1;
