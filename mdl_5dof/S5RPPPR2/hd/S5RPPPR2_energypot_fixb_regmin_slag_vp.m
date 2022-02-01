% Calculate minimal parameter regressor of potential energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:32
% EndTime: 2022-01-23 08:59:32
% DurationCPUTime: 0.13s
% Computational Cost: add. (93->64), mult. (192->104), div. (0->0), fcn. (219->10), ass. (0->41)
t190 = sin(pkin(7));
t213 = g(3) * t190;
t212 = t190 * qJ(3) + pkin(1);
t189 = sin(pkin(8));
t194 = sin(qJ(5));
t211 = t189 * t194;
t196 = cos(qJ(5));
t210 = t189 * t196;
t192 = cos(pkin(8));
t209 = t190 * t192;
t195 = sin(qJ(1));
t208 = t190 * t195;
t197 = cos(qJ(1));
t207 = t190 * t197;
t191 = cos(pkin(9));
t193 = cos(pkin(7));
t206 = t193 * t191;
t205 = t195 * t189;
t204 = t195 * t192;
t203 = t197 * qJ(2);
t202 = t197 * t189;
t201 = t197 * t192;
t200 = -t193 * qJ(3) + pkin(5);
t199 = g(1) * t197 + g(2) * t195;
t178 = t193 * t205 + t201;
t198 = g(1) * (t193 * t202 - t204) + g(2) * t178 + t189 * t213;
t188 = sin(pkin(9));
t187 = t195 * qJ(2);
t185 = g(1) * t195 - g(2) * t197;
t184 = t192 * pkin(3) + t189 * qJ(4) + pkin(2);
t183 = pkin(2) * t193 + t212;
t182 = -t189 * pkin(3) + qJ(4) * t192 - qJ(2);
t181 = t193 * t201 + t205;
t179 = t193 * t204 - t202;
t177 = t191 * t211 + t196 * t192;
t176 = t190 * t188 + t192 * t206;
t175 = -t193 * t188 + t191 * t209;
t174 = -g(3) * t193 + t199 * t190;
t173 = t184 * t193 + t212;
t172 = -t176 * t194 + t193 * t210;
t1 = [0, -t199, t185, -t199 * t193 - t213, t174, -t185, -g(1) * (t197 * pkin(1) + t187) - g(2) * (t195 * pkin(1) - t203) - g(3) * pkin(5), -g(1) * t181 - g(2) * t179 - g(3) * t209, t198, -t174, -g(1) * (t183 * t197 + t187) - g(2) * (t183 * t195 - t203) - g(3) * (t190 * pkin(2) + t200), -g(1) * (t181 * t191 + t188 * t207) - g(2) * (t179 * t191 + t188 * t208) - g(3) * t175, -g(1) * (-t181 * t188 + t191 * t207) - g(2) * (-t179 * t188 + t191 * t208) - g(3) * (-t188 * t209 - t206), -t198, -g(1) * (t173 * t197 - t182 * t195) - g(2) * (t173 * t195 + t182 * t197) - g(3) * (t184 * t190 + t200), 0, 0, 0, 0, 0, -g(1) * ((t176 * t196 + t193 * t211) * t197 + t195 * (t191 * t210 - t194 * t192)) - g(2) * ((t176 * t195 - t191 * t202) * t196 + t178 * t194) - g(3) * (t175 * t196 + t190 * t211), -g(1) * (t172 * t197 - t195 * t177) - g(2) * (t172 * t195 + t197 * t177) - g(3) * (-t175 * t194 + t190 * t210);];
U_reg = t1;
