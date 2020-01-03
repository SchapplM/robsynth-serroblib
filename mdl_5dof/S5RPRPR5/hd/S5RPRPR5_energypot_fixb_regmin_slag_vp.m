% Calculate minimal parameter regressor of potential energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:43:01
% EndTime: 2020-01-03 11:43:01
% DurationCPUTime: 0.11s
% Computational Cost: add. (67->39), mult. (92->57), div. (0->0), fcn. (94->8), ass. (0->27)
t174 = sin(pkin(8));
t192 = g(1) * t174;
t172 = qJ(3) + pkin(9) + qJ(5);
t169 = sin(t172);
t178 = sin(qJ(1));
t191 = t178 * t169;
t170 = cos(t172);
t190 = t178 * t170;
t177 = sin(qJ(3));
t189 = t178 * t177;
t179 = cos(qJ(3));
t188 = t178 * t179;
t180 = cos(qJ(1));
t187 = t180 * t169;
t186 = t180 * t170;
t185 = t180 * t177;
t184 = t180 * t179;
t183 = pkin(3) * t177 + qJ(2);
t182 = -g(2) * t178 + g(3) * t180;
t171 = t179 * pkin(3) + pkin(2);
t175 = cos(pkin(8));
t176 = -qJ(4) - pkin(6);
t181 = t171 * t175 - t174 * t176;
t173 = t178 * pkin(1);
t168 = g(2) * t180 + g(3) * t178;
t167 = g(1) * t175 + t182 * t174;
t1 = [0, t182, -t168, t182 * t175 - t192, -t167, t168, -g(1) * pkin(5) - g(2) * (-t180 * qJ(2) + t173) - g(3) * (-t180 * pkin(1) - t178 * qJ(2)), 0, 0, 0, 0, 0, -t179 * t192 - g(2) * (t175 * t188 - t185) - g(3) * (-t175 * t184 - t189), t177 * t192 - g(2) * (-t175 * t189 - t184) - g(3) * (t175 * t185 - t188), t167, -g(1) * (t174 * t171 + t175 * t176 + pkin(5)) - g(2) * t173 + (-g(2) * t181 + g(3) * t183) * t178 + (g(2) * t183 - g(3) * (-pkin(1) - t181)) * t180, 0, 0, 0, 0, 0, -t170 * t192 - g(2) * (t175 * t190 - t187) - g(3) * (-t175 * t186 - t191), t169 * t192 - g(2) * (-t175 * t191 - t186) - g(3) * (t175 * t187 - t190);];
U_reg = t1;
