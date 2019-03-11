% Calculate minimal parameter regressor of potential energy for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:45:55
% EndTime: 2019-03-09 16:45:55
% DurationCPUTime: 0.09s
% Computational Cost: add. (126->45), mult. (148->54), div. (0->0), fcn. (154->8), ass. (0->30)
t205 = qJ(2) + qJ(3);
t202 = sin(t205);
t203 = cos(t205);
t210 = cos(qJ(2));
t225 = t210 * pkin(2) + pkin(3) * t203 + qJ(4) * t202 + pkin(1);
t223 = g(3) * t203;
t206 = sin(qJ(5));
t208 = sin(qJ(1));
t221 = t208 * t206;
t209 = cos(qJ(5));
t220 = t208 * t209;
t211 = cos(qJ(1));
t219 = t211 * t206;
t218 = t211 * t209;
t207 = sin(qJ(2));
t217 = t207 * pkin(2) + t202 * pkin(3) + pkin(6);
t212 = -pkin(8) - pkin(7);
t216 = t225 * t208 + t211 * t212;
t215 = g(1) * t211 + g(2) * t208;
t214 = -t208 * t212 + t225 * t211;
t187 = -t202 * t218 + t221;
t189 = t202 * t220 + t219;
t213 = g(1) * t187 - g(2) * t189 + t209 * t223;
t192 = g(1) * t208 - g(2) * t211;
t190 = t202 * t221 - t218;
t188 = t202 * t219 + t220;
t186 = g(3) * t202 + t215 * t203;
t185 = t215 * t202 - t223;
t184 = -g(1) * t188 - g(2) * t190 + t206 * t223;
t1 = [0, -t215, t192, 0, 0, 0, 0, 0, -g(3) * t207 - t215 * t210, -g(3) * t210 + t215 * t207, 0, 0, 0, 0, 0, -t186, t185, -t192, t186, -t185, -g(1) * t214 - g(2) * t216 - g(3) * (-t203 * qJ(4) + t217) 0, 0, 0, 0, 0, t184, t213, t184, -t186, -t213, -g(1) * (t208 * pkin(4) + t188 * pkin(5) + t187 * qJ(6) + t214) - g(2) * (-t211 * pkin(4) + t190 * pkin(5) - t189 * qJ(6) + t216) - g(3) * (t202 * pkin(9) + t217) + (-g(3) * (-pkin(5) * t206 + qJ(6) * t209 - qJ(4)) - t215 * pkin(9)) * t203;];
U_reg  = t1;
