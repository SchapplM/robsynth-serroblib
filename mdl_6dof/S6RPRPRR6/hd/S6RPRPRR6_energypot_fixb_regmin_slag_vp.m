% Calculate minimal parameter regressor of potential energy for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:52:52
% EndTime: 2019-03-09 03:52:52
% DurationCPUTime: 0.09s
% Computational Cost: add. (119->51), mult. (118->74), div. (0->0), fcn. (124->12), ass. (0->34)
t223 = pkin(10) + qJ(3);
t218 = sin(t223);
t245 = g(3) * t218;
t222 = pkin(11) + qJ(5);
t221 = qJ(6) + t222;
t214 = sin(t221);
t229 = sin(qJ(1));
t244 = t229 * t214;
t215 = cos(t221);
t243 = t229 * t215;
t217 = sin(t222);
t242 = t229 * t217;
t219 = cos(t222);
t241 = t229 * t219;
t224 = sin(pkin(11));
t240 = t229 * t224;
t226 = cos(pkin(11));
t239 = t229 * t226;
t230 = cos(qJ(1));
t238 = t230 * t214;
t237 = t230 * t215;
t236 = t230 * t217;
t235 = t230 * t219;
t234 = t230 * t224;
t233 = t230 * t226;
t232 = g(1) * t230 + g(2) * t229;
t220 = cos(t223);
t227 = cos(pkin(10));
t231 = t227 * pkin(2) + pkin(3) * t220 + qJ(4) * t218 + pkin(1);
t228 = -pkin(7) - qJ(2);
t225 = sin(pkin(10));
t213 = g(1) * t229 - g(2) * t230;
t212 = -g(3) * t220 + t232 * t218;
t1 = [0, -t232, t213, -g(3) * t225 - t232 * t227, -g(3) * t227 + t232 * t225, -t213, -g(1) * (t230 * pkin(1) + t229 * qJ(2)) - g(2) * (t229 * pkin(1) - t230 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t232 * t220 - t245, t212, -g(1) * (t220 * t233 + t240) - g(2) * (t220 * t239 - t234) - t226 * t245, -g(1) * (-t220 * t234 + t239) - g(2) * (-t220 * t240 - t233) + t224 * t245, -t212, -g(3) * (t225 * pkin(2) + t218 * pkin(3) - t220 * qJ(4) + pkin(6)) + (-g(1) * t231 - g(2) * t228) * t230 + (g(1) * t228 - g(2) * t231) * t229, 0, 0, 0, 0, 0, -g(1) * (t220 * t235 + t242) - g(2) * (t220 * t241 - t236) - t219 * t245, -g(1) * (-t220 * t236 + t241) - g(2) * (-t220 * t242 - t235) + t217 * t245, 0, 0, 0, 0, 0, -g(1) * (t220 * t237 + t244) - g(2) * (t220 * t243 - t238) - t215 * t245, -g(1) * (-t220 * t238 + t243) - g(2) * (-t220 * t244 - t237) + t214 * t245;];
U_reg  = t1;
