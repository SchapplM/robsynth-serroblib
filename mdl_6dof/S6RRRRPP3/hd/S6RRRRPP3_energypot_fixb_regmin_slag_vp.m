% Calculate minimal parameter regressor of potential energy for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:56:28
% EndTime: 2019-03-09 20:56:29
% DurationCPUTime: 0.10s
% Computational Cost: add. (156->49), mult. (189->59), div. (0->0), fcn. (205->8), ass. (0->31)
t207 = qJ(2) + qJ(3);
t205 = cos(t207);
t212 = cos(qJ(2));
t230 = t212 * pkin(2) + pkin(3) * t205 + pkin(1);
t204 = sin(t207);
t208 = sin(qJ(4));
t228 = t204 * t208;
t210 = sin(qJ(1));
t227 = t204 * t210;
t211 = cos(qJ(4));
t226 = t204 * t211;
t213 = cos(qJ(1));
t225 = t204 * t213;
t224 = t210 * t208;
t223 = t210 * t211;
t222 = t213 * t208;
t221 = t213 * t211;
t220 = g(1) * t213 + g(2) * t210;
t209 = sin(qJ(2));
t219 = t209 * pkin(2) + t204 * pkin(3) + pkin(4) * t226 + qJ(5) * t228 + pkin(6);
t186 = t205 * t224 + t221;
t187 = t205 * t223 - t222;
t214 = -pkin(8) - pkin(7);
t218 = t187 * pkin(4) + pkin(9) * t227 + t186 * qJ(5) + t230 * t210 + t213 * t214;
t188 = t205 * t222 - t223;
t217 = g(1) * t188 + g(2) * t186 + g(3) * t228;
t189 = t205 * t221 + t224;
t216 = g(1) * t189 + g(2) * t187 + g(3) * t226;
t215 = t189 * pkin(4) + pkin(9) * t225 + t188 * qJ(5) - t210 * t214 + t230 * t213;
t183 = -g(3) * t205 + t220 * t204;
t1 = [0, -t220, g(1) * t210 - g(2) * t213, 0, 0, 0, 0, 0, -g(3) * t209 - t220 * t212, -g(3) * t212 + t220 * t209, 0, 0, 0, 0, 0, -g(3) * t204 - t220 * t205, t183, 0, 0, 0, 0, 0, -t216, t217, -t183, t216, -t217, -g(1) * t215 - g(2) * t218 - g(3) * (-t205 * pkin(9) + t219) -t183, -t217, -t216, -g(1) * (pkin(5) * t225 + t189 * qJ(6) + t215) - g(2) * (pkin(5) * t227 + t187 * qJ(6) + t218) - g(3) * (qJ(6) * t226 + (-pkin(5) - pkin(9)) * t205 + t219);];
U_reg  = t1;
