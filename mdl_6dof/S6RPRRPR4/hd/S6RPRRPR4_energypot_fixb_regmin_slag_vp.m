% Calculate minimal parameter regressor of potential energy for
% S6RPRRPR4
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
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:10:15
% EndTime: 2019-03-09 05:10:16
% DurationCPUTime: 0.08s
% Computational Cost: add. (130->47), mult. (111->66), div. (0->0), fcn. (113->12), ass. (0->30)
t220 = pkin(10) + qJ(3);
t217 = qJ(4) + t220;
t211 = sin(t217);
t237 = g(3) * t211;
t219 = pkin(11) + qJ(6);
t213 = sin(t219);
t225 = sin(qJ(1));
t236 = t225 * t213;
t215 = cos(t219);
t235 = t225 * t215;
t221 = sin(pkin(11));
t234 = t225 * t221;
t223 = cos(pkin(11));
t233 = t225 * t223;
t226 = cos(qJ(1));
t232 = t226 * t213;
t231 = t226 * t215;
t230 = t226 * t221;
t229 = t226 * t223;
t228 = g(1) * t226 + g(2) * t225;
t212 = cos(t217);
t216 = cos(t220);
t224 = cos(pkin(10));
t227 = t224 * pkin(2) + pkin(3) * t216 + pkin(4) * t212 + qJ(5) * t211 + pkin(1);
t222 = sin(pkin(10));
t218 = -pkin(8) - pkin(7) - qJ(2);
t214 = sin(t220);
t210 = g(1) * t225 - g(2) * t226;
t208 = -g(3) * t212 + t228 * t211;
t1 = [0, -t228, t210, -g(3) * t222 - t228 * t224, -g(3) * t224 + t228 * t222, -t210, -g(1) * (t226 * pkin(1) + t225 * qJ(2)) - g(2) * (t225 * pkin(1) - t226 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t214 - t228 * t216, -g(3) * t216 + t228 * t214, 0, 0, 0, 0, 0, -t228 * t212 - t237, t208, -g(1) * (t212 * t229 + t234) - g(2) * (t212 * t233 - t230) - t223 * t237, -g(1) * (-t212 * t230 + t233) - g(2) * (-t212 * t234 - t229) + t221 * t237, -t208, -g(3) * (t222 * pkin(2) + pkin(3) * t214 + t211 * pkin(4) - t212 * qJ(5) + pkin(6)) + (-g(1) * t227 - g(2) * t218) * t226 + (g(1) * t218 - g(2) * t227) * t225, 0, 0, 0, 0, 0, -g(1) * (t212 * t231 + t236) - g(2) * (t212 * t235 - t232) - t215 * t237, -g(1) * (-t212 * t232 + t235) - g(2) * (-t212 * t236 - t231) + t213 * t237;];
U_reg  = t1;
