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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:57:11
% EndTime: 2019-12-05 17:57:11
% DurationCPUTime: 0.09s
% Computational Cost: add. (67->39), mult. (92->57), div. (0->0), fcn. (94->8), ass. (0->27)
t186 = cos(qJ(3));
t176 = t186 * pkin(3) + pkin(2);
t181 = sin(pkin(8));
t182 = cos(pkin(8));
t183 = -qJ(4) - pkin(6);
t201 = t176 * t182 - t181 * t183;
t200 = g(1) * t181;
t177 = qJ(3) + pkin(9) + qJ(5);
t174 = sin(t177);
t185 = sin(qJ(1));
t197 = t185 * t174;
t175 = cos(t177);
t196 = t185 * t175;
t184 = sin(qJ(3));
t195 = t185 * t184;
t194 = t185 * t186;
t187 = cos(qJ(1));
t193 = t187 * t174;
t192 = t187 * t175;
t191 = t187 * t184;
t190 = t187 * t186;
t189 = t187 * pkin(1) + t185 * qJ(2);
t188 = g(2) * t185 - g(3) * t187;
t179 = t187 * qJ(2);
t173 = g(2) * t187 + g(3) * t185;
t172 = g(1) * t182 + t188 * t181;
t1 = [0, t188, t173, t188 * t182 - t200, -t172, -t173, -g(1) * pkin(5) - g(2) * (-t185 * pkin(1) + t179) - g(3) * t189, 0, 0, 0, 0, 0, -t186 * t200 - g(2) * (-t182 * t194 + t191) - g(3) * (t182 * t190 + t195), t184 * t200 - g(2) * (t182 * t195 + t190) - g(3) * (-t182 * t191 + t194), t172, -g(1) * (t181 * t176 + t182 * t183 + pkin(5)) - g(2) * (pkin(3) * t191 + t179) - g(3) * (t201 * t187 + t189) + (-g(2) * (-pkin(1) - t201) - g(3) * pkin(3) * t184) * t185, 0, 0, 0, 0, 0, -t175 * t200 - g(2) * (-t182 * t196 + t193) - g(3) * (t182 * t192 + t197), t174 * t200 - g(2) * (t182 * t197 + t192) - g(3) * (-t182 * t193 + t196);];
U_reg = t1;
