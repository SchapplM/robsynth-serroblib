% Calculate minimal parameter regressor of potential energy for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR13_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:33:51
% EndTime: 2019-12-31 20:33:51
% DurationCPUTime: 0.09s
% Computational Cost: add. (72->40), mult. (96->62), div. (0->0), fcn. (105->10), ass. (0->24)
t183 = sin(qJ(2));
t196 = g(3) * t183;
t184 = sin(qJ(1));
t185 = cos(qJ(2));
t195 = t184 * t185;
t180 = pkin(9) + qJ(4);
t179 = qJ(5) + t180;
t175 = sin(t179);
t186 = cos(qJ(1));
t194 = t186 * t175;
t176 = cos(t179);
t193 = t186 * t176;
t177 = sin(t180);
t192 = t186 * t177;
t178 = cos(t180);
t191 = t186 * t178;
t181 = sin(pkin(9));
t190 = t186 * t181;
t182 = cos(pkin(9));
t189 = t186 * t182;
t188 = g(1) * t186 + g(2) * t184;
t187 = pkin(2) * t185 + qJ(3) * t183 + pkin(1);
t174 = -g(3) * t185 + t188 * t183;
t1 = [0, -t188, g(1) * t184 - g(2) * t186, 0, 0, 0, 0, 0, -t188 * t185 - t196, t174, -g(1) * (t184 * t181 + t185 * t189) - g(2) * (t182 * t195 - t190) - t182 * t196, -g(1) * (t184 * t182 - t185 * t190) - g(2) * (-t181 * t195 - t189) + t181 * t196, -t174, -g(3) * (t183 * pkin(2) - t185 * qJ(3) + pkin(5)) + (g(2) * pkin(6) - g(1) * t187) * t186 + (-g(1) * pkin(6) - g(2) * t187) * t184, 0, 0, 0, 0, 0, -g(1) * (t184 * t177 + t185 * t191) - g(2) * (t178 * t195 - t192) - t178 * t196, -g(1) * (t184 * t178 - t185 * t192) - g(2) * (-t177 * t195 - t191) + t177 * t196, 0, 0, 0, 0, 0, -g(1) * (t184 * t175 + t185 * t193) - g(2) * (t176 * t195 - t194) - t176 * t196, -g(1) * (t184 * t176 - t185 * t194) - g(2) * (-t175 * t195 - t193) + t175 * t196;];
U_reg = t1;
