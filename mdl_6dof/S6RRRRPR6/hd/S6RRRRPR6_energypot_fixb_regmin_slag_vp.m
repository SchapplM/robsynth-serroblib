% Calculate minimal parameter regressor of potential energy for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:23:11
% EndTime: 2019-03-09 22:23:11
% DurationCPUTime: 0.10s
% Computational Cost: add. (103->45), mult. (108->67), div. (0->0), fcn. (117->10), ass. (0->22)
t230 = sin(qJ(2));
t240 = g(3) * t230;
t228 = qJ(3) + qJ(4);
t225 = sin(t228);
t229 = sin(qJ(3));
t239 = pkin(3) * t229 + pkin(4) * t225 + pkin(7);
t231 = sin(qJ(1));
t233 = cos(qJ(2));
t238 = t231 * t233;
t234 = cos(qJ(1));
t237 = t233 * t234;
t236 = g(1) * t234 + g(2) * t231;
t226 = cos(t228);
t232 = cos(qJ(3));
t220 = pkin(3) * t232 + pkin(4) * t226 + pkin(2);
t227 = -qJ(5) - pkin(9) - pkin(8);
t235 = t220 * t233 - t227 * t230 + pkin(1);
t224 = pkin(11) + qJ(6) + t228;
t223 = cos(t224);
t222 = sin(t224);
t219 = -g(3) * t233 + t236 * t230;
t1 = [0, -t236, g(1) * t231 - g(2) * t234, 0, 0, 0, 0, 0, -t236 * t233 - t240, t219, 0, 0, 0, 0, 0, -g(1) * (t231 * t229 + t232 * t237) - g(2) * (-t229 * t234 + t232 * t238) - t232 * t240, -g(1) * (-t229 * t237 + t231 * t232) - g(2) * (-t229 * t238 - t232 * t234) + t229 * t240, 0, 0, 0, 0, 0, -g(1) * (t231 * t225 + t226 * t237) - g(2) * (-t225 * t234 + t226 * t238) - t226 * t240, -g(1) * (-t225 * t237 + t231 * t226) - g(2) * (-t225 * t238 - t226 * t234) + t225 * t240, -t219, -g(3) * (t220 * t230 + t227 * t233 + pkin(6)) + (-g(1) * t235 + g(2) * t239) * t234 + (-g(1) * t239 - g(2) * t235) * t231, 0, 0, 0, 0, 0, -g(1) * (t231 * t222 + t223 * t237) - g(2) * (-t222 * t234 + t223 * t238) - t223 * t240, -g(1) * (-t222 * t237 + t231 * t223) - g(2) * (-t222 * t238 - t223 * t234) + t222 * t240;];
U_reg  = t1;
