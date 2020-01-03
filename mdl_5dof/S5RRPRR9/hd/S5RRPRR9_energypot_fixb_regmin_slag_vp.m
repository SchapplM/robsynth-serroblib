% Calculate minimal parameter regressor of potential energy for
% S5RRPRR9
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
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:21:57
% EndTime: 2019-12-31 20:21:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (54->31), mult. (66->46), div. (0->0), fcn. (71->10), ass. (0->25)
t180 = qJ(2) + pkin(9);
t198 = g(3) * sin(t180);
t181 = qJ(4) + qJ(5);
t178 = sin(t181);
t185 = sin(qJ(1));
t197 = t185 * t178;
t179 = cos(t181);
t196 = t185 * t179;
t183 = sin(qJ(4));
t195 = t185 * t183;
t186 = cos(qJ(4));
t194 = t185 * t186;
t188 = cos(qJ(1));
t193 = t188 * t178;
t192 = t188 * t179;
t191 = t188 * t183;
t190 = t188 * t186;
t189 = g(1) * t188 + g(2) * t185;
t187 = cos(qJ(2));
t184 = sin(qJ(2));
t182 = -pkin(6) - qJ(3);
t177 = cos(t180);
t175 = t187 * pkin(2) + pkin(1);
t174 = g(1) * t185 - g(2) * t188;
t1 = [0, -t189, t174, 0, 0, 0, 0, 0, -g(3) * t184 - t189 * t187, -g(3) * t187 + t189 * t184, -t174, -g(1) * (t188 * t175 - t185 * t182) - g(2) * (t185 * t175 + t188 * t182) - g(3) * (t184 * pkin(2) + pkin(5)), 0, 0, 0, 0, 0, -g(1) * (t177 * t190 + t195) - g(2) * (t177 * t194 - t191) - t186 * t198, -g(1) * (-t177 * t191 + t194) - g(2) * (-t177 * t195 - t190) + t183 * t198, 0, 0, 0, 0, 0, -g(1) * (t177 * t192 + t197) - g(2) * (t177 * t196 - t193) - t179 * t198, -g(1) * (-t177 * t193 + t196) - g(2) * (-t177 * t197 - t192) + t178 * t198;];
U_reg = t1;
