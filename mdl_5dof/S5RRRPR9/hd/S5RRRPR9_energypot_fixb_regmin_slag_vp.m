% Calculate minimal parameter regressor of potential energy for
% S5RRRPR9
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
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:24:37
% EndTime: 2019-12-31 21:24:38
% DurationCPUTime: 0.07s
% Computational Cost: add. (62->34), mult. (83->50), div. (0->0), fcn. (88->8), ass. (0->22)
t180 = sin(qJ(2));
t193 = g(3) * t180;
t181 = sin(qJ(1));
t183 = cos(qJ(2));
t192 = t181 * t183;
t177 = qJ(3) + pkin(9) + qJ(5);
t174 = sin(t177);
t184 = cos(qJ(1));
t191 = t184 * t174;
t175 = cos(t177);
t190 = t184 * t175;
t179 = sin(qJ(3));
t189 = t184 * t179;
t182 = cos(qJ(3));
t188 = t184 * t182;
t187 = pkin(3) * t179 + pkin(6);
t186 = g(1) * t184 + g(2) * t181;
t176 = t182 * pkin(3) + pkin(2);
t178 = -qJ(4) - pkin(7);
t185 = t176 * t183 - t178 * t180 + pkin(1);
t173 = -g(3) * t183 + t186 * t180;
t1 = [0, -t186, g(1) * t181 - g(2) * t184, 0, 0, 0, 0, 0, -t186 * t183 - t193, t173, 0, 0, 0, 0, 0, -g(1) * (t181 * t179 + t183 * t188) - g(2) * (t182 * t192 - t189) - t182 * t193, -g(1) * (t181 * t182 - t183 * t189) - g(2) * (-t179 * t192 - t188) + t179 * t193, -t173, -g(3) * (t180 * t176 + t183 * t178 + pkin(5)) + (-g(1) * t185 + g(2) * t187) * t184 + (-g(1) * t187 - g(2) * t185) * t181, 0, 0, 0, 0, 0, -g(1) * (t181 * t174 + t183 * t190) - g(2) * (t175 * t192 - t191) - t175 * t193, -g(1) * (t181 * t175 - t183 * t191) - g(2) * (-t174 * t192 - t190) + t174 * t193;];
U_reg = t1;
