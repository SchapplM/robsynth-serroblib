% Calculate minimal parameter regressor of potential energy for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x34]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:25:18
% EndTime: 2019-03-09 07:25:18
% DurationCPUTime: 0.07s
% Computational Cost: add. (66->36), mult. (85->57), div. (0->0), fcn. (94->10), ass. (0->28)
t186 = cos(qJ(3));
t200 = g(3) * t186;
t181 = qJ(4) + qJ(5);
t180 = qJ(6) + t181;
t176 = sin(t180);
t184 = sin(qJ(1));
t199 = t184 * t176;
t177 = cos(t180);
t198 = t184 * t177;
t178 = sin(t181);
t197 = t184 * t178;
t179 = cos(t181);
t196 = t184 * t179;
t182 = sin(qJ(4));
t195 = t184 * t182;
t185 = cos(qJ(4));
t194 = t184 * t185;
t187 = cos(qJ(1));
t193 = t187 * t176;
t192 = t187 * t177;
t191 = t187 * t178;
t190 = t187 * t179;
t189 = t187 * t182;
t188 = t187 * t185;
t174 = g(1) * t184 - g(2) * t187;
t183 = sin(qJ(3));
t175 = g(1) * t187 + g(2) * t184;
t1 = [0, -t175, t174, t175, -t174, -g(1) * (t187 * pkin(1) + t184 * qJ(2)) - g(2) * (t184 * pkin(1) - t187 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t174 * t183 - t200, g(3) * t183 - t174 * t186, 0, 0, 0, 0, 0, -g(1) * (t183 * t194 + t189) - g(2) * (-t183 * t188 + t195) - t185 * t200, -g(1) * (-t183 * t195 + t188) - g(2) * (t183 * t189 + t194) + t182 * t200, 0, 0, 0, 0, 0, -g(1) * (t183 * t196 + t191) - g(2) * (-t183 * t190 + t197) - t179 * t200, -g(1) * (-t183 * t197 + t190) - g(2) * (t183 * t191 + t196) + t178 * t200, 0, 0, 0, 0, 0, -g(1) * (t183 * t198 + t193) - g(2) * (-t183 * t192 + t199) - t177 * t200, -g(1) * (-t183 * t199 + t192) - g(2) * (t183 * t193 + t198) + t176 * t200;];
U_reg  = t1;
