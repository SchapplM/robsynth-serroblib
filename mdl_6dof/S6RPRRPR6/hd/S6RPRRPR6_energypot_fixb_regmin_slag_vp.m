% Calculate minimal parameter regressor of potential energy for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:17:42
% EndTime: 2019-03-09 05:17:42
% DurationCPUTime: 0.10s
% Computational Cost: add. (103->45), mult. (105->62), div. (0->0), fcn. (107->10), ass. (0->29)
t221 = pkin(10) + qJ(3);
t218 = sin(t221);
t241 = g(3) * t218;
t220 = qJ(4) + pkin(11) + qJ(6);
t214 = sin(t220);
t227 = sin(qJ(1));
t240 = t227 * t214;
t215 = cos(t220);
t239 = t227 * t215;
t226 = sin(qJ(4));
t238 = t227 * t226;
t228 = cos(qJ(4));
t237 = t227 * t228;
t229 = cos(qJ(1));
t236 = t229 * t214;
t235 = t229 * t215;
t234 = t229 * t226;
t233 = t229 * t228;
t232 = pkin(4) * t226 + pkin(7) + qJ(2);
t231 = g(1) * t229 + g(2) * t227;
t217 = t228 * pkin(4) + pkin(3);
t219 = cos(t221);
t223 = cos(pkin(10));
t224 = -qJ(5) - pkin(8);
t230 = t223 * pkin(2) + t217 * t219 - t218 * t224 + pkin(1);
t222 = sin(pkin(10));
t213 = g(1) * t227 - g(2) * t229;
t212 = -g(3) * t219 + t231 * t218;
t1 = [0, -t231, t213, -g(3) * t222 - t231 * t223, -g(3) * t223 + t231 * t222, -t213, -g(1) * (t229 * pkin(1) + t227 * qJ(2)) - g(2) * (t227 * pkin(1) - t229 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t231 * t219 - t241, t212, 0, 0, 0, 0, 0, -g(1) * (t219 * t233 + t238) - g(2) * (t219 * t237 - t234) - t228 * t241, -g(1) * (-t219 * t234 + t237) - g(2) * (-t219 * t238 - t233) + t226 * t241, -t212, -g(3) * (t222 * pkin(2) + t218 * t217 + t219 * t224 + pkin(6)) + (-g(1) * t230 + g(2) * t232) * t229 + (-g(1) * t232 - g(2) * t230) * t227, 0, 0, 0, 0, 0, -g(1) * (t219 * t235 + t240) - g(2) * (t219 * t239 - t236) - t215 * t241, -g(1) * (-t219 * t236 + t239) - g(2) * (-t219 * t240 - t235) + t214 * t241;];
U_reg  = t1;
