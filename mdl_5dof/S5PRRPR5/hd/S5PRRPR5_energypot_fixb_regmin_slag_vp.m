% Calculate minimal parameter regressor of potential energy for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:27:53
% EndTime: 2019-12-05 16:27:53
% DurationCPUTime: 0.12s
% Computational Cost: add. (87->49), mult. (170->84), div. (0->0), fcn. (206->12), ass. (0->32)
t183 = qJ(3) + pkin(10);
t181 = sin(t183);
t185 = sin(pkin(5));
t204 = t181 * t185;
t190 = sin(qJ(3));
t203 = t185 * t190;
t191 = sin(qJ(2));
t202 = t185 * t191;
t193 = cos(qJ(3));
t201 = t185 * t193;
t194 = cos(qJ(2));
t200 = t185 * t194;
t187 = cos(pkin(5));
t199 = t187 * t190;
t198 = t187 * t191;
t197 = t187 * t194;
t184 = sin(pkin(9));
t186 = cos(pkin(9));
t175 = t184 * t191 - t186 * t197;
t177 = t184 * t197 + t186 * t191;
t195 = -g(1) * t177 - g(2) * t175 + g(3) * t200;
t192 = cos(qJ(5));
t189 = sin(qJ(5));
t188 = -qJ(4) - pkin(7);
t182 = cos(t183);
t180 = t193 * pkin(3) + pkin(2);
t178 = -t184 * t198 + t186 * t194;
t176 = t184 * t194 + t186 * t198;
t174 = t187 * t181 + t182 * t202;
t173 = t178 * t182 + t184 * t204;
t172 = t176 * t182 - t186 * t204;
t1 = [-g(3) * qJ(1), 0, -g(1) * t178 - g(2) * t176 - g(3) * t202, -t195, 0, 0, 0, 0, 0, -g(1) * (t178 * t193 + t184 * t203) - g(2) * (t176 * t193 - t186 * t203) - g(3) * (t191 * t201 + t199), -g(1) * (-t178 * t190 + t184 * t201) - g(2) * (-t176 * t190 - t186 * t201) - g(3) * (t187 * t193 - t190 * t202), t195, -g(1) * (t186 * pkin(1) - t177 * t188 + t178 * t180) - g(2) * (t184 * pkin(1) - t175 * t188 + t176 * t180) - g(3) * (pkin(3) * t199 + t187 * pkin(6) + qJ(1)) + (-g(3) * (t180 * t191 + t188 * t194) + (-g(1) * t184 + g(2) * t186) * (pkin(3) * t190 + pkin(6))) * t185, 0, 0, 0, 0, 0, -g(1) * (t173 * t192 + t177 * t189) - g(2) * (t172 * t192 + t175 * t189) - g(3) * (t174 * t192 - t189 * t200), -g(1) * (-t173 * t189 + t177 * t192) - g(2) * (-t172 * t189 + t175 * t192) - g(3) * (-t174 * t189 - t192 * t200);];
U_reg = t1;
