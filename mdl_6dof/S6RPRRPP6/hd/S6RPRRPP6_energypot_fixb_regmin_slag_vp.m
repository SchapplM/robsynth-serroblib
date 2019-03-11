% Calculate minimal parameter regressor of potential energy for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:48:27
% EndTime: 2019-03-09 04:48:27
% DurationCPUTime: 0.13s
% Computational Cost: add. (104->52), mult. (147->67), div. (0->0), fcn. (150->8), ass. (0->36)
t194 = cos(qJ(4));
t181 = t194 * pkin(4) + pkin(3);
t190 = -qJ(5) - pkin(8);
t192 = sin(qJ(3));
t195 = cos(qJ(3));
t214 = t181 * t192 + t190 * t195;
t196 = cos(qJ(1));
t213 = g(2) * t196;
t212 = g(3) * t195;
t189 = qJ(4) + pkin(9);
t182 = sin(t189);
t193 = sin(qJ(1));
t209 = t193 * t182;
t183 = cos(t189);
t208 = t193 * t183;
t191 = sin(qJ(4));
t207 = t193 * t191;
t206 = t193 * t194;
t205 = t196 * t182;
t204 = t196 * t183;
t203 = t196 * t191;
t202 = t196 * t194;
t201 = t196 * pkin(1) + t193 * qJ(2);
t186 = t193 * pkin(1);
t200 = pkin(4) * t207 + t193 * pkin(7) + t186;
t175 = g(1) * t193 - t213;
t199 = t195 * t181 - t192 * t190 + pkin(2) + pkin(6);
t198 = pkin(4) * t203 + t196 * pkin(7) + t214 * t193 + t201;
t197 = (-qJ(2) - t214) * t213;
t176 = g(1) * t196 + g(2) * t193;
t173 = -g(3) * t192 + t175 * t195;
t172 = -t192 * t204 + t209;
t171 = t192 * t205 + t208;
t170 = t192 * t208 + t205;
t169 = t192 * t209 - t204;
t1 = [0, -t176, t175, t176, -t175, -g(1) * t201 - g(2) * (-t196 * qJ(2) + t186) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t175 * t192 - t212, -t173, 0, 0, 0, 0, 0, -g(1) * (t192 * t206 + t203) - g(2) * (-t192 * t202 + t207) - t194 * t212, -g(1) * (-t192 * t207 + t202) - g(2) * (t192 * t203 + t206) + t191 * t212, t173, -g(1) * t198 - g(2) * t200 - g(3) * t199 - t197, -g(1) * t170 - g(2) * t172 - t183 * t212, t173, -g(1) * t169 + g(2) * t171 - t182 * t212, -g(1) * (t170 * pkin(5) + t169 * qJ(6) + t198) - g(2) * (t172 * pkin(5) - t171 * qJ(6) + t200) - g(3) * ((pkin(5) * t183 + qJ(6) * t182) * t195 + t199) - t197;];
U_reg  = t1;
