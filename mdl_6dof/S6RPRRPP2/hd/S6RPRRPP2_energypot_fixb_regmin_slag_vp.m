% Calculate minimal parameter regressor of potential energy for
% S6RPRRPP2
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
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:33:20
% EndTime: 2019-03-09 04:33:20
% DurationCPUTime: 0.11s
% Computational Cost: add. (161->45), mult. (182->61), div. (0->0), fcn. (195->8), ass. (0->30)
t195 = cos(qJ(3));
t211 = pkin(3) * t195 + pkin(2);
t209 = qJ(2) + pkin(6);
t190 = qJ(1) + pkin(9);
t185 = sin(t190);
t192 = sin(qJ(3));
t208 = t185 * t192;
t186 = cos(t190);
t207 = t186 * t192;
t191 = sin(qJ(4));
t206 = t191 * t192;
t205 = t191 * t195;
t194 = cos(qJ(4));
t204 = t192 * t194;
t203 = t194 * t195;
t202 = g(1) * t186 + g(2) * t185;
t193 = sin(qJ(1));
t196 = cos(qJ(1));
t201 = -g(1) * t196 - g(2) * t193;
t200 = t192 * pkin(3) + pkin(4) * t204 + qJ(5) * t206 + t209;
t171 = t185 * t205 + t186 * t194;
t173 = -t185 * t194 + t186 * t205;
t199 = g(1) * t173 + g(2) * t171 + g(3) * t206;
t174 = t185 * t191 + t186 * t203;
t198 = t196 * pkin(1) + t174 * pkin(4) + t185 * pkin(7) + pkin(8) * t207 + t173 * qJ(5) + t211 * t186;
t172 = t185 * t203 - t186 * t191;
t197 = t193 * pkin(1) + t172 * pkin(4) - t186 * pkin(7) + pkin(8) * t208 + t171 * qJ(5) + t211 * t185;
t168 = -g(3) * t195 + t202 * t192;
t167 = -g(1) * t174 - g(2) * t172 - g(3) * t204;
t1 = [0, t201, g(1) * t193 - g(2) * t196, t201 * pkin(1) - g(3) * t209, 0, 0, 0, 0, 0, -g(3) * t192 - t202 * t195, t168, 0, 0, 0, 0, 0, t167, t199, t167, -t168, -t199, -g(1) * t198 - g(2) * t197 - g(3) * (-t195 * pkin(8) + t200) t167, -t199, t168, -g(1) * (t174 * pkin(5) - qJ(6) * t207 + t198) - g(2) * (t172 * pkin(5) - qJ(6) * t208 + t197) - g(3) * (pkin(5) * t204 + (-pkin(8) + qJ(6)) * t195 + t200);];
U_reg  = t1;
