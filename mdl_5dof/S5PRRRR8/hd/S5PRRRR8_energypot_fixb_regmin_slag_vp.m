% Calculate minimal parameter regressor of potential energy for
% S5PRRRR8
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:39
% EndTime: 2019-12-05 17:16:39
% DurationCPUTime: 0.09s
% Computational Cost: add. (85->40), mult. (155->77), div. (0->0), fcn. (198->12), ass. (0->29)
t175 = sin(pkin(10));
t176 = sin(pkin(5));
t192 = t175 * t176;
t177 = cos(pkin(10));
t191 = t176 * t177;
t180 = sin(qJ(3));
t190 = t176 * t180;
t181 = sin(qJ(2));
t189 = t176 * t181;
t183 = cos(qJ(3));
t188 = t176 * t183;
t184 = cos(qJ(2));
t187 = t176 * t184;
t178 = cos(pkin(5));
t186 = t178 * t181;
t185 = t178 * t184;
t182 = cos(qJ(5));
t179 = sin(qJ(5));
t174 = qJ(3) + qJ(4);
t173 = cos(t174);
t172 = sin(t174);
t171 = -t175 * t186 + t177 * t184;
t170 = t175 * t185 + t177 * t181;
t169 = t175 * t184 + t177 * t186;
t168 = t175 * t181 - t177 * t185;
t167 = t178 * t172 + t173 * t189;
t166 = t171 * t173 + t172 * t192;
t165 = t169 * t173 - t172 * t191;
t1 = [-g(3) * qJ(1), 0, -g(1) * t171 - g(2) * t169 - g(3) * t189, g(1) * t170 + g(2) * t168 - g(3) * t187, 0, 0, 0, 0, 0, -g(1) * (t171 * t183 + t175 * t190) - g(2) * (t169 * t183 - t177 * t190) - g(3) * (t178 * t180 + t181 * t188), -g(1) * (-t171 * t180 + t175 * t188) - g(2) * (-t169 * t180 - t177 * t188) - g(3) * (t178 * t183 - t180 * t189), 0, 0, 0, 0, 0, -g(1) * t166 - g(2) * t165 - g(3) * t167, -g(1) * (-t171 * t172 + t173 * t192) - g(2) * (-t169 * t172 - t173 * t191) - g(3) * (-t172 * t189 + t178 * t173), 0, 0, 0, 0, 0, -g(1) * (t166 * t182 + t170 * t179) - g(2) * (t165 * t182 + t168 * t179) - g(3) * (t167 * t182 - t179 * t187), -g(1) * (-t166 * t179 + t170 * t182) - g(2) * (-t165 * t179 + t168 * t182) - g(3) * (-t167 * t179 - t182 * t187);];
U_reg = t1;
