% Calculate minimal parameter regressor of potential energy for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:36:19
% EndTime: 2019-03-09 02:36:19
% DurationCPUTime: 0.09s
% Computational Cost: add. (68->36), mult. (86->52), div. (0->0), fcn. (88->10), ass. (0->26)
t182 = pkin(10) + qJ(4);
t176 = cos(t182);
t200 = g(3) * t176;
t183 = qJ(5) + qJ(6);
t177 = sin(t183);
t187 = sin(qJ(1));
t199 = t187 * t177;
t178 = cos(t183);
t198 = t187 * t178;
t186 = sin(qJ(5));
t197 = t187 * t186;
t188 = cos(qJ(5));
t196 = t187 * t188;
t189 = cos(qJ(1));
t195 = t189 * t177;
t194 = t189 * t178;
t193 = t189 * t186;
t192 = t189 * t188;
t191 = pkin(1) * t189 + qJ(2) * t187;
t190 = pkin(1) * t187 - qJ(2) * t189;
t173 = g(1) * t187 - g(2) * t189;
t185 = cos(pkin(10));
t184 = sin(pkin(10));
t175 = sin(t182);
t174 = g(1) * t189 + g(2) * t187;
t1 = [0, -t174, t173, t174, -t173, -g(3) * pkin(6) - g(1) * t191 - g(2) * t190, -g(3) * t185 - t173 * t184, g(3) * t184 - t173 * t185, -t174, -g(1) * (qJ(3) * t189 + t191) - g(2) * (qJ(3) * t187 + t190) - g(3) * (pkin(2) + pkin(6)) 0, 0, 0, 0, 0, -t173 * t175 - t200, g(3) * t175 - t173 * t176, 0, 0, 0, 0, 0, -g(1) * (t175 * t196 + t193) - g(2) * (-t175 * t192 + t197) - t188 * t200, -g(1) * (-t175 * t197 + t192) - g(2) * (t175 * t193 + t196) + t186 * t200, 0, 0, 0, 0, 0, -g(1) * (t175 * t198 + t195) - g(2) * (-t175 * t194 + t199) - t178 * t200, -g(1) * (-t175 * t199 + t194) - g(2) * (t175 * t195 + t198) + t177 * t200;];
U_reg  = t1;
