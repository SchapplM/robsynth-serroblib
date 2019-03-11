% Calculate minimal parameter regressor of potential energy for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:44:41
% EndTime: 2019-03-09 04:44:41
% DurationCPUTime: 0.10s
% Computational Cost: add. (161->51), mult. (198->66), div. (0->0), fcn. (211->8), ass. (0->32)
t206 = pkin(9) + qJ(3);
t204 = cos(t206);
t208 = cos(pkin(9));
t228 = t208 * pkin(2) + pkin(3) * t204 + pkin(1);
t203 = sin(t206);
t210 = sin(qJ(4));
t226 = t203 * t210;
t211 = sin(qJ(1));
t225 = t203 * t211;
t212 = cos(qJ(4));
t224 = t203 * t212;
t213 = cos(qJ(1));
t223 = t203 * t213;
t222 = t211 * t210;
t221 = t211 * t212;
t220 = t213 * t210;
t219 = t213 * t212;
t218 = g(1) * t213 + g(2) * t211;
t207 = sin(pkin(9));
t217 = t207 * pkin(2) + t203 * pkin(3) + pkin(4) * t224 + qJ(5) * t226 + pkin(6);
t186 = t204 * t222 + t219;
t187 = t204 * t221 - t220;
t209 = -pkin(7) - qJ(2);
t216 = t187 * pkin(4) + pkin(8) * t225 + t186 * qJ(5) + t213 * t209 + t228 * t211;
t188 = t204 * t220 - t221;
t215 = g(1) * t188 + g(2) * t186 + g(3) * t226;
t189 = t204 * t219 + t222;
t214 = t189 * pkin(4) + pkin(8) * t223 + t188 * qJ(5) - t211 * t209 + t228 * t213;
t193 = g(1) * t211 - g(2) * t213;
t183 = -g(3) * t204 + t218 * t203;
t182 = -g(1) * t189 - g(2) * t187 - g(3) * t224;
t1 = [0, -t218, t193, -g(3) * t207 - t218 * t208, -g(3) * t208 + t218 * t207, -t193, -g(1) * (t213 * pkin(1) + t211 * qJ(2)) - g(2) * (t211 * pkin(1) - t213 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t203 - t218 * t204, t183, 0, 0, 0, 0, 0, t182, t215, t182, -t183, -t215, -g(1) * t214 - g(2) * t216 - g(3) * (-t204 * pkin(8) + t217) t182, -t215, t183, -g(1) * (t189 * pkin(5) - qJ(6) * t223 + t214) - g(2) * (t187 * pkin(5) - qJ(6) * t225 + t216) - g(3) * (pkin(5) * t224 + (-pkin(8) + qJ(6)) * t204 + t217);];
U_reg  = t1;
