% Calculate minimal parameter regressor of potential energy for
% S5PRRPR6
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
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:46
% EndTime: 2019-12-05 16:32:46
% DurationCPUTime: 0.12s
% Computational Cost: add. (121->58), mult. (265->97), div. (0->0), fcn. (334->12), ass. (0->32)
t187 = sin(pkin(5));
t202 = pkin(6) * t187;
t191 = sin(qJ(3));
t201 = t187 * t191;
t192 = sin(qJ(2));
t200 = t187 * t192;
t193 = cos(qJ(3));
t199 = t187 * t193;
t194 = cos(qJ(2));
t198 = t187 * t194;
t190 = cos(pkin(5));
t197 = t190 * t192;
t196 = t190 * t194;
t186 = sin(pkin(9));
t189 = cos(pkin(9));
t174 = t186 * t194 + t189 * t197;
t169 = t174 * t191 + t189 * t199;
t176 = -t186 * t197 + t189 * t194;
t171 = t176 * t191 - t186 * t199;
t177 = -t190 * t193 + t191 * t200;
t195 = g(1) * t171 + g(2) * t169 + g(3) * t177;
t188 = cos(pkin(10));
t185 = sin(pkin(10));
t184 = pkin(10) + qJ(5);
t183 = cos(t184);
t182 = sin(t184);
t178 = t190 * t191 + t192 * t199;
t175 = t186 * t196 + t189 * t192;
t173 = t186 * t192 - t189 * t196;
t172 = t176 * t193 + t186 * t201;
t170 = t174 * t193 - t189 * t201;
t1 = [-g(3) * qJ(1), 0, -g(1) * t176 - g(2) * t174 - g(3) * t200, g(1) * t175 + g(2) * t173 - g(3) * t198, 0, 0, 0, 0, 0, -g(1) * t172 - g(2) * t170 - g(3) * t178, t195, -g(1) * (t172 * t188 + t175 * t185) - g(2) * (t170 * t188 + t173 * t185) - g(3) * (t178 * t188 - t185 * t198), -g(1) * (-t172 * t185 + t175 * t188) - g(2) * (-t170 * t185 + t173 * t188) - g(3) * (-t178 * t185 - t188 * t198), -t195, -g(1) * (t189 * pkin(1) + t176 * pkin(2) + t172 * pkin(3) + t175 * pkin(7) + t171 * qJ(4) + t186 * t202) - g(2) * (t186 * pkin(1) + t174 * pkin(2) + t170 * pkin(3) + t173 * pkin(7) + t169 * qJ(4) - t189 * t202) - g(3) * (t178 * pkin(3) + t190 * pkin(6) + t177 * qJ(4) + qJ(1) + (pkin(2) * t192 - pkin(7) * t194) * t187), 0, 0, 0, 0, 0, -g(1) * (t172 * t183 + t175 * t182) - g(2) * (t170 * t183 + t173 * t182) - g(3) * (t178 * t183 - t182 * t198), -g(1) * (-t172 * t182 + t175 * t183) - g(2) * (-t170 * t182 + t173 * t183) - g(3) * (-t178 * t182 - t183 * t198);];
U_reg = t1;
