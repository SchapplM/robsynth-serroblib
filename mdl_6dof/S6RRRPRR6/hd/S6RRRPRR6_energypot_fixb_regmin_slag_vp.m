% Calculate minimal parameter regressor of potential energy for
% S6RRRPRR6
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
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:29:29
% EndTime: 2019-03-09 18:29:29
% DurationCPUTime: 0.10s
% Computational Cost: add. (100->43), mult. (103->65), div. (0->0), fcn. (112->10), ass. (0->22)
t230 = sin(qJ(2));
t240 = g(3) * t230;
t231 = sin(qJ(1));
t233 = cos(qJ(2));
t239 = t231 * t233;
t234 = cos(qJ(1));
t238 = t233 * t234;
t229 = sin(qJ(3));
t237 = pkin(3) * t229 + pkin(7);
t227 = qJ(3) + pkin(11) + qJ(5);
t236 = g(1) * t234 + g(2) * t231;
t232 = cos(qJ(3));
t226 = pkin(3) * t232 + pkin(2);
t228 = -qJ(4) - pkin(8);
t235 = t226 * t233 - t228 * t230 + pkin(1);
t225 = qJ(6) + t227;
t224 = cos(t227);
t223 = sin(t227);
t222 = cos(t225);
t221 = sin(t225);
t220 = -g(3) * t233 + t236 * t230;
t1 = [0, -t236, g(1) * t231 - g(2) * t234, 0, 0, 0, 0, 0, -t236 * t233 - t240, t220, 0, 0, 0, 0, 0, -g(1) * (t231 * t229 + t232 * t238) - g(2) * (-t229 * t234 + t232 * t239) - t232 * t240, -g(1) * (-t229 * t238 + t231 * t232) - g(2) * (-t229 * t239 - t232 * t234) + t229 * t240, -t220, -g(3) * (t226 * t230 + t228 * t233 + pkin(6)) + (-g(1) * t235 + g(2) * t237) * t234 + (-g(1) * t237 - g(2) * t235) * t231, 0, 0, 0, 0, 0, -g(1) * (t231 * t223 + t224 * t238) - g(2) * (-t223 * t234 + t224 * t239) - t224 * t240, -g(1) * (-t223 * t238 + t231 * t224) - g(2) * (-t223 * t239 - t224 * t234) + t223 * t240, 0, 0, 0, 0, 0, -g(1) * (t231 * t221 + t222 * t238) - g(2) * (-t221 * t234 + t222 * t239) - t222 * t240, -g(1) * (-t221 * t238 + t231 * t222) - g(2) * (-t221 * t239 - t222 * t234) + t221 * t240;];
U_reg  = t1;
