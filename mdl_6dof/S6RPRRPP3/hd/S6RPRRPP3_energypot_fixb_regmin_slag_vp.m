% Calculate minimal parameter regressor of potential energy for
% S6RPRRPP3
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
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:36:58
% EndTime: 2019-03-09 04:36:58
% DurationCPUTime: 0.09s
% Computational Cost: add. (161->48), mult. (182->61), div. (0->0), fcn. (195->8), ass. (0->30)
t193 = cos(qJ(3));
t210 = pkin(3) * t193 + pkin(2);
t208 = qJ(2) + pkin(6);
t188 = qJ(1) + pkin(9);
t182 = sin(t188);
t190 = sin(qJ(3));
t207 = t182 * t190;
t183 = cos(t188);
t206 = t183 * t190;
t189 = sin(qJ(4));
t205 = t189 * t190;
t204 = t189 * t193;
t192 = cos(qJ(4));
t203 = t190 * t192;
t202 = t192 * t193;
t201 = g(1) * t183 + g(2) * t182;
t191 = sin(qJ(1));
t194 = cos(qJ(1));
t200 = -g(1) * t194 - g(2) * t191;
t199 = t190 * pkin(3) + pkin(4) * t203 + qJ(5) * t205 + t208;
t167 = t182 * t204 + t183 * t192;
t169 = -t182 * t192 + t183 * t204;
t198 = g(1) * t169 + g(2) * t167 + g(3) * t205;
t168 = t182 * t202 - t183 * t189;
t170 = t182 * t189 + t183 * t202;
t197 = g(1) * t170 + g(2) * t168 + g(3) * t203;
t196 = t194 * pkin(1) + t170 * pkin(4) + t182 * pkin(7) + pkin(8) * t206 + t169 * qJ(5) + t210 * t183;
t195 = t191 * pkin(1) + t168 * pkin(4) - t183 * pkin(7) + pkin(8) * t207 + t167 * qJ(5) + t210 * t182;
t164 = -g(3) * t193 + t201 * t190;
t1 = [0, t200, g(1) * t191 - g(2) * t194, t200 * pkin(1) - g(3) * t208, 0, 0, 0, 0, 0, -g(3) * t190 - t201 * t193, t164, 0, 0, 0, 0, 0, -t197, t198, -t164, t197, -t198, -g(1) * t196 - g(2) * t195 - g(3) * (-t193 * pkin(8) + t199) -t164, -t198, -t197, -g(1) * (pkin(5) * t206 + t170 * qJ(6) + t196) - g(2) * (pkin(5) * t207 + t168 * qJ(6) + t195) - g(3) * (qJ(6) * t203 + (-pkin(5) - pkin(8)) * t193 + t199);];
U_reg  = t1;
