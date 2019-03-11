% Calculate minimal parameter regressor of potential energy for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:02:29
% EndTime: 2019-03-09 07:02:29
% DurationCPUTime: 0.11s
% Computational Cost: add. (91->34), mult. (79->60), div. (0->0), fcn. (88->12), ass. (0->25)
t212 = sin(qJ(3));
t225 = g(3) * t212;
t209 = qJ(1) + pkin(11);
t204 = sin(t209);
t215 = cos(qJ(3));
t224 = t204 * t215;
t205 = cos(t209);
t223 = t205 * t215;
t210 = qJ(4) + qJ(5);
t206 = sin(t210);
t222 = t206 * t215;
t207 = cos(t210);
t221 = t207 * t215;
t211 = sin(qJ(4));
t220 = t211 * t215;
t214 = cos(qJ(4));
t219 = t214 * t215;
t218 = g(1) * t205 + g(2) * t204;
t213 = sin(qJ(1));
t216 = cos(qJ(1));
t217 = -g(1) * t216 - g(2) * t213;
t208 = qJ(6) + t210;
t203 = cos(t208);
t202 = sin(t208);
t1 = [0, t217, g(1) * t213 - g(2) * t216, -g(3) * (qJ(2) + pkin(6)) + t217 * pkin(1), 0, 0, 0, 0, 0, -t218 * t215 - t225, -g(3) * t215 + t218 * t212, 0, 0, 0, 0, 0, -g(1) * (t204 * t211 + t205 * t219) - g(2) * (t204 * t219 - t205 * t211) - t214 * t225, -g(1) * (t204 * t214 - t205 * t220) - g(2) * (-t204 * t220 - t205 * t214) + t211 * t225, 0, 0, 0, 0, 0, -g(1) * (t204 * t206 + t205 * t221) - g(2) * (t204 * t221 - t205 * t206) - t207 * t225, -g(1) * (t204 * t207 - t205 * t222) - g(2) * (-t204 * t222 - t205 * t207) + t206 * t225, 0, 0, 0, 0, 0, -g(1) * (t204 * t202 + t203 * t223) - g(2) * (-t205 * t202 + t203 * t224) - t203 * t225, -g(1) * (-t202 * t223 + t204 * t203) - g(2) * (-t202 * t224 - t205 * t203) + t202 * t225;];
U_reg  = t1;
