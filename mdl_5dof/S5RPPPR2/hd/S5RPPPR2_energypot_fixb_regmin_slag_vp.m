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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:22:54
% EndTime: 2020-01-03 11:22:54
% DurationCPUTime: 0.12s
% Computational Cost: add. (93->58), mult. (210->91), div. (0->0), fcn. (241->10), ass. (0->36)
t190 = cos(pkin(7));
t209 = pkin(2) * t190;
t186 = sin(pkin(8));
t187 = sin(pkin(7));
t208 = t186 * t187;
t189 = cos(pkin(8));
t207 = t187 * t189;
t192 = sin(qJ(1));
t206 = t187 * t192;
t194 = cos(qJ(1));
t205 = t187 * t194;
t204 = t192 * qJ(2);
t203 = t192 * t186;
t202 = t192 * t189;
t201 = t194 * t186;
t200 = t194 * t189;
t184 = t192 * pkin(1);
t199 = qJ(3) * t206 + t192 * t209 + t184;
t198 = t187 * pkin(2) - t190 * qJ(3) + pkin(5);
t197 = -g(2) * t192 + g(3) * t194;
t175 = t190 * t203 + t200;
t177 = t190 * t201 - t202;
t196 = g(1) * t208 + g(2) * t175 - g(3) * t177;
t195 = (g(2) * qJ(2) - g(3) * (-qJ(3) * t187 - pkin(1) - t209)) * t194;
t193 = cos(qJ(5));
t191 = sin(qJ(5));
t188 = cos(pkin(9));
t185 = sin(pkin(9));
t179 = g(2) * t194 + g(3) * t192;
t178 = -t190 * t200 - t203;
t176 = t190 * t202 - t201;
t174 = -t190 * t185 + t188 * t207;
t173 = g(1) * t190 + t197 * t187;
t172 = t178 * t188 - t185 * t205;
t171 = t176 * t188 + t185 * t206;
t1 = [0, t197, -t179, -g(1) * t187 + t197 * t190, -t173, t179, -g(1) * pkin(5) - g(2) * (-t194 * qJ(2) + t184) - g(3) * (-t194 * pkin(1) - t204), -g(1) * t207 - g(2) * t176 - g(3) * t178, t196, t173, -g(1) * t198 - g(2) * t199 + g(3) * t204 + t195, -g(1) * t174 - g(2) * t171 - g(3) * t172, -g(1) * (-t185 * t207 - t190 * t188) - g(2) * (-t176 * t185 + t188 * t206) - g(3) * (-t178 * t185 - t188 * t205), -t196, -g(1) * ((pkin(3) * t189 + qJ(4) * t186) * t187 + t198) - g(2) * (t176 * pkin(3) + t175 * qJ(4) + t199) - g(3) * (t178 * pkin(3) - t177 * qJ(4) - t204) + t195, 0, 0, 0, 0, 0, -g(1) * (t174 * t193 + t191 * t208) - g(2) * (t171 * t193 + t175 * t191) - g(3) * (t172 * t193 - t177 * t191), -g(1) * (-t174 * t191 + t193 * t208) - g(2) * (-t171 * t191 + t175 * t193) - g(3) * (-t172 * t191 - t177 * t193);];
U_reg = t1;
