% Calculate minimal parameter regressor of potential energy for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:58:43
% EndTime: 2019-03-09 01:58:43
% DurationCPUTime: 0.07s
% Computational Cost: add. (112->43), mult. (94->57), div. (0->0), fcn. (89->10), ass. (0->29)
t174 = pkin(10) + qJ(4);
t168 = sin(t174);
t194 = g(3) * t168;
t179 = qJ(2) + pkin(6);
t193 = g(3) * t179;
t175 = qJ(1) + pkin(9);
t169 = sin(t175);
t181 = sin(qJ(5));
t192 = t169 * t181;
t183 = cos(qJ(5));
t191 = t169 * t183;
t171 = cos(t175);
t190 = t171 * t181;
t189 = t171 * t183;
t188 = pkin(5) * t181 + pkin(7) + qJ(3);
t187 = g(1) * t171 + g(2) * t169;
t182 = sin(qJ(1));
t184 = cos(qJ(1));
t186 = -g(1) * t184 - g(2) * t182;
t167 = t183 * pkin(5) + pkin(4);
t170 = cos(t174);
t177 = cos(pkin(10));
t178 = -qJ(6) - pkin(8);
t185 = t177 * pkin(3) + t167 * t170 - t168 * t178 + pkin(2);
t176 = sin(pkin(10));
t173 = t184 * pkin(1);
t172 = t182 * pkin(1);
t165 = -g(3) * t170 + t187 * t168;
t1 = [0, t186, g(1) * t182 - g(2) * t184, t186 * pkin(1) - t193, -g(3) * t176 - t187 * t177, -g(3) * t177 + t187 * t176, -g(1) * t169 + g(2) * t171, -g(1) * (t171 * pkin(2) + t169 * qJ(3) + t173) - g(2) * (t169 * pkin(2) - t171 * qJ(3) + t172) - t193, 0, 0, 0, 0, 0, -t187 * t170 - t194, t165, 0, 0, 0, 0, 0, -g(1) * (t170 * t189 + t192) - g(2) * (t170 * t191 - t190) - t183 * t194, -g(1) * (-t170 * t190 + t191) - g(2) * (-t170 * t192 - t189) + t181 * t194, -t165, -g(1) * t173 - g(2) * t172 - g(3) * (t176 * pkin(3) + t168 * t167 + t170 * t178 + t179) + (-g(1) * t185 + g(2) * t188) * t171 + (-g(1) * t188 - g(2) * t185) * t169;];
U_reg  = t1;
