% Calculate minimal parameter regressor of potential energy for
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:31:26
% EndTime: 2019-03-10 01:31:26
% DurationCPUTime: 0.12s
% Computational Cost: add. (151->52), mult. (148->74), div. (0->0), fcn. (165->10), ass. (0->33)
t225 = sin(qJ(2));
t241 = g(3) * t225;
t223 = qJ(3) + qJ(4);
t219 = sin(t223);
t224 = sin(qJ(3));
t240 = t224 * pkin(3) + pkin(4) * t219 + pkin(7);
t226 = sin(qJ(1));
t228 = cos(qJ(2));
t239 = t226 * t228;
t221 = qJ(5) + t223;
t217 = sin(t221);
t229 = cos(qJ(1));
t238 = t229 * t217;
t218 = cos(t221);
t237 = t229 * t218;
t236 = t229 * t219;
t220 = cos(t223);
t235 = t229 * t220;
t234 = t229 * t224;
t227 = cos(qJ(3));
t233 = t229 * t227;
t232 = g(1) * t229 + g(2) * t226;
t214 = t227 * pkin(3) + pkin(4) * t220 + pkin(2);
t222 = -pkin(10) - pkin(9) - pkin(8);
t231 = t214 * t228 - t222 * t225 + pkin(1);
t209 = t217 * t239 + t237;
t211 = -t226 * t218 + t228 * t238;
t230 = g(1) * t211 + g(2) * t209 + t217 * t241;
t213 = -g(3) * t228 + t232 * t225;
t212 = t226 * t217 + t228 * t237;
t210 = t218 * t239 - t238;
t208 = -g(1) * t212 - g(2) * t210 - t218 * t241;
t1 = [0, -t232, g(1) * t226 - g(2) * t229, 0, 0, 0, 0, 0, -t232 * t228 - t241, t213, 0, 0, 0, 0, 0, -g(1) * (t226 * t224 + t228 * t233) - g(2) * (t227 * t239 - t234) - t227 * t241, -g(1) * (t226 * t227 - t228 * t234) - g(2) * (-t224 * t239 - t233) + t224 * t241, 0, 0, 0, 0, 0, -g(1) * (t226 * t219 + t228 * t235) - g(2) * (t220 * t239 - t236) - t220 * t241, -g(1) * (t226 * t220 - t228 * t236) - g(2) * (-t219 * t239 - t235) + t219 * t241, 0, 0, 0, 0, 0, t208, t230, t208, -t213, -t230, -g(1) * (t212 * pkin(5) + t211 * qJ(6)) - g(2) * (t210 * pkin(5) + t209 * qJ(6)) - g(3) * (t228 * t222 + pkin(6)) - (pkin(5) * t218 + qJ(6) * t217 + t214) * t241 + (-g(1) * t231 + g(2) * t240) * t229 + (-g(1) * t240 - g(2) * t231) * t226;];
U_reg  = t1;
