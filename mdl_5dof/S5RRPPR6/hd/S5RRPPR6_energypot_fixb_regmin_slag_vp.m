% Calculate minimal parameter regressor of potential energy for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:04
% EndTime: 2019-12-31 19:33:04
% DurationCPUTime: 0.08s
% Computational Cost: add. (80->39), mult. (91->57), div. (0->0), fcn. (93->10), ass. (0->30)
t184 = qJ(2) + pkin(8);
t179 = sin(t184);
t205 = g(3) * t179;
t188 = sin(qJ(2));
t204 = t188 * pkin(2) + pkin(5);
t183 = pkin(9) + qJ(5);
t178 = sin(t183);
t189 = sin(qJ(1));
t203 = t189 * t178;
t180 = cos(t183);
t202 = t189 * t180;
t185 = sin(pkin(9));
t201 = t189 * t185;
t186 = cos(pkin(9));
t200 = t189 * t186;
t191 = cos(qJ(1));
t199 = t191 * t178;
t198 = t191 * t180;
t197 = t191 * t185;
t196 = t191 * t186;
t190 = cos(qJ(2));
t177 = t190 * pkin(2) + pkin(1);
t187 = -pkin(6) - qJ(3);
t195 = t189 * t177 + t191 * t187;
t194 = t191 * t177 - t189 * t187;
t193 = g(1) * t191 + g(2) * t189;
t181 = cos(t184);
t192 = pkin(3) * t181 + qJ(4) * t179;
t173 = g(1) * t189 - g(2) * t191;
t1 = [0, -t193, t173, 0, 0, 0, 0, 0, -g(3) * t188 - t193 * t190, -g(3) * t190 + t193 * t188, -t173, -g(1) * t194 - g(2) * t195 - g(3) * t204, -g(1) * (t181 * t196 + t201) - g(2) * (t181 * t200 - t197) - t186 * t205, -g(1) * (-t181 * t197 + t200) - g(2) * (-t181 * t201 - t196) + t185 * t205, g(3) * t181 - t193 * t179, -g(1) * (t192 * t191 + t194) - g(2) * (t192 * t189 + t195) - g(3) * (t179 * pkin(3) - t181 * qJ(4) + t204), 0, 0, 0, 0, 0, -g(1) * (t181 * t198 + t203) - g(2) * (t181 * t202 - t199) - t180 * t205, -g(1) * (-t181 * t199 + t202) - g(2) * (-t181 * t203 - t198) + t178 * t205;];
U_reg = t1;
