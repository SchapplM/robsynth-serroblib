% Calculate minimal parameter regressor of potential energy for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:17:16
% EndTime: 2019-12-31 21:17:16
% DurationCPUTime: 0.08s
% Computational Cost: add. (80->37), mult. (89->53), div. (0->0), fcn. (94->10), ass. (0->26)
t174 = qJ(2) + qJ(3);
t171 = sin(t174);
t192 = g(3) * t171;
t173 = pkin(9) + qJ(5);
t169 = sin(t173);
t178 = sin(qJ(1));
t191 = t178 * t169;
t170 = cos(t173);
t190 = t178 * t170;
t175 = sin(pkin(9));
t189 = t178 * t175;
t176 = cos(pkin(9));
t188 = t178 * t176;
t180 = cos(qJ(1));
t187 = t180 * t169;
t186 = t180 * t170;
t185 = t180 * t175;
t184 = t180 * t176;
t183 = g(1) * t180 + g(2) * t178;
t172 = cos(t174);
t179 = cos(qJ(2));
t182 = t179 * pkin(2) + pkin(3) * t172 + qJ(4) * t171 + pkin(1);
t181 = -pkin(7) - pkin(6);
t177 = sin(qJ(2));
t167 = -g(3) * t172 + t183 * t171;
t1 = [0, -t183, g(1) * t178 - g(2) * t180, 0, 0, 0, 0, 0, -g(3) * t177 - t183 * t179, -g(3) * t179 + t183 * t177, 0, 0, 0, 0, 0, -t183 * t172 - t192, t167, -g(1) * (t172 * t184 + t189) - g(2) * (t172 * t188 - t185) - t176 * t192, -g(1) * (-t172 * t185 + t188) - g(2) * (-t172 * t189 - t184) + t175 * t192, -t167, -g(3) * (t177 * pkin(2) + t171 * pkin(3) - t172 * qJ(4) + pkin(5)) + (-g(1) * t182 - g(2) * t181) * t180 + (g(1) * t181 - g(2) * t182) * t178, 0, 0, 0, 0, 0, -g(1) * (t172 * t186 + t191) - g(2) * (t172 * t190 - t187) - t170 * t192, -g(1) * (-t172 * t187 + t190) - g(2) * (-t172 * t191 - t186) + t169 * t192;];
U_reg = t1;
