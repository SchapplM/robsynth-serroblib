% Calculate minimal parameter regressor of potential energy for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:21:10
% EndTime: 2019-12-05 17:21:10
% DurationCPUTime: 0.08s
% Computational Cost: add. (83->40), mult. (181->75), div. (0->0), fcn. (234->12), ass. (0->27)
t183 = sin(pkin(5));
t187 = sin(qJ(3));
t197 = t183 * t187;
t188 = sin(qJ(2));
t196 = t183 * t188;
t190 = cos(qJ(3));
t195 = t183 * t190;
t191 = cos(qJ(2));
t194 = t183 * t191;
t185 = cos(pkin(5));
t193 = t185 * t188;
t192 = t185 * t191;
t189 = cos(qJ(4));
t186 = sin(qJ(4));
t184 = cos(pkin(10));
t182 = sin(pkin(10));
t181 = qJ(4) + qJ(5);
t180 = cos(t181);
t179 = sin(t181);
t178 = t185 * t187 + t188 * t195;
t177 = -t182 * t193 + t184 * t191;
t176 = t182 * t192 + t184 * t188;
t175 = t182 * t191 + t184 * t193;
t174 = t182 * t188 - t184 * t192;
t173 = t177 * t190 + t182 * t197;
t172 = t175 * t190 - t184 * t197;
t1 = [-g(3) * qJ(1), 0, -g(1) * t177 - g(2) * t175 - g(3) * t196, g(1) * t176 + g(2) * t174 - g(3) * t194, 0, 0, 0, 0, 0, -g(1) * t173 - g(2) * t172 - g(3) * t178, -g(1) * (-t177 * t187 + t182 * t195) - g(2) * (-t175 * t187 - t184 * t195) - g(3) * (t185 * t190 - t187 * t196), 0, 0, 0, 0, 0, -g(1) * (t173 * t189 + t176 * t186) - g(2) * (t172 * t189 + t174 * t186) - g(3) * (t178 * t189 - t186 * t194), -g(1) * (-t173 * t186 + t176 * t189) - g(2) * (-t172 * t186 + t174 * t189) - g(3) * (-t178 * t186 - t189 * t194), 0, 0, 0, 0, 0, -g(1) * (t173 * t180 + t176 * t179) - g(2) * (t172 * t180 + t174 * t179) - g(3) * (t178 * t180 - t179 * t194), -g(1) * (-t173 * t179 + t176 * t180) - g(2) * (-t172 * t179 + t174 * t180) - g(3) * (-t178 * t179 - t180 * t194);];
U_reg = t1;
