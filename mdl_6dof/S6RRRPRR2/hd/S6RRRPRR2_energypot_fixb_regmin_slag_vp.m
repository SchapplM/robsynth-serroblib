% Calculate minimal parameter regressor of potential energy for
% S6RRRPRR2
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
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:08:39
% EndTime: 2019-03-09 18:08:39
% DurationCPUTime: 0.06s
% Computational Cost: add. (84->36), mult. (79->52), div. (0->0), fcn. (84->12), ass. (0->28)
t228 = qJ(2) + qJ(3);
t221 = pkin(11) + t228;
t244 = g(3) * sin(t221);
t227 = qJ(5) + qJ(6);
t222 = sin(t227);
t231 = sin(qJ(1));
t243 = t231 * t222;
t224 = cos(t227);
t242 = t231 * t224;
t229 = sin(qJ(5));
t241 = t231 * t229;
t232 = cos(qJ(5));
t240 = t231 * t232;
t234 = cos(qJ(1));
t239 = t234 * t222;
t238 = t234 * t224;
t237 = t234 * t229;
t236 = t234 * t232;
t235 = g(1) * t234 + g(2) * t231;
t233 = cos(qJ(2));
t230 = sin(qJ(2));
t226 = -qJ(4) - pkin(8) - pkin(7);
t225 = cos(t228);
t223 = sin(t228);
t220 = cos(t221);
t218 = g(1) * t231 - g(2) * t234;
t217 = t233 * pkin(2) + pkin(3) * t225 + pkin(1);
t1 = [0, -t235, t218, 0, 0, 0, 0, 0, -g(3) * t230 - t235 * t233, -g(3) * t233 + t235 * t230, 0, 0, 0, 0, 0, -g(3) * t223 - t235 * t225, -g(3) * t225 + t235 * t223, -t218, -g(1) * (t234 * t217 - t231 * t226) - g(2) * (t231 * t217 + t234 * t226) - g(3) * (t230 * pkin(2) + pkin(3) * t223 + pkin(6)) 0, 0, 0, 0, 0, -g(1) * (t220 * t236 + t241) - g(2) * (t220 * t240 - t237) - t232 * t244, -g(1) * (-t220 * t237 + t240) - g(2) * (-t220 * t241 - t236) + t229 * t244, 0, 0, 0, 0, 0, -g(1) * (t220 * t238 + t243) - g(2) * (t220 * t242 - t239) - t224 * t244, -g(1) * (-t220 * t239 + t242) - g(2) * (-t220 * t243 - t238) + t222 * t244;];
U_reg  = t1;
