% Calculate minimal parameter regressor of potential energy for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:25
% EndTime: 2019-12-05 16:56:25
% DurationCPUTime: 0.10s
% Computational Cost: add. (96->52), mult. (219->84), div. (0->0), fcn. (270->10), ass. (0->32)
t166 = sin(pkin(9));
t168 = cos(pkin(9));
t173 = sin(qJ(2));
t169 = cos(pkin(5));
t176 = cos(qJ(2));
t178 = t169 * t176;
t156 = t166 * t173 - t168 * t178;
t171 = sin(qJ(4));
t185 = t156 * t171;
t158 = t166 * t178 + t168 * t173;
t184 = t158 * t171;
t167 = sin(pkin(5));
t172 = sin(qJ(3));
t183 = t167 * t172;
t182 = t167 * t173;
t175 = cos(qJ(3));
t181 = t167 * t175;
t180 = t167 * t176;
t179 = t169 * t173;
t157 = t166 * t176 + t168 * t179;
t152 = t157 * t172 + t168 * t181;
t159 = -t166 * t179 + t168 * t176;
t154 = t159 * t172 - t166 * t181;
t160 = -t169 * t175 + t172 * t182;
t177 = g(1) * t154 + g(2) * t152 + g(3) * t160;
t174 = cos(qJ(4));
t170 = -qJ(5) - pkin(8);
t165 = t174 * pkin(4) + pkin(3);
t161 = t169 * t172 + t173 * t181;
t155 = t159 * t175 + t166 * t183;
t153 = t157 * t175 - t168 * t183;
t1 = [-g(3) * qJ(1), 0, -g(1) * t159 - g(2) * t157 - g(3) * t182, g(1) * t158 + g(2) * t156 - g(3) * t180, 0, 0, 0, 0, 0, -g(1) * t155 - g(2) * t153 - g(3) * t161, t177, 0, 0, 0, 0, 0, -g(1) * (t155 * t174 + t184) - g(2) * (t153 * t174 + t185) - g(3) * (t161 * t174 - t171 * t180), -g(1) * (-t155 * t171 + t158 * t174) - g(2) * (-t153 * t171 + t156 * t174) - g(3) * (-t161 * t171 - t174 * t180), -t177, -g(1) * (t168 * pkin(1) + t159 * pkin(2) + pkin(4) * t184 + t158 * pkin(7) - t154 * t170 + t155 * t165) - g(2) * (t166 * pkin(1) + t157 * pkin(2) + pkin(4) * t185 + t156 * pkin(7) - t152 * t170 + t153 * t165) - g(3) * (t169 * pkin(6) - t160 * t170 + t161 * t165 + qJ(1)) + (-g(3) * (pkin(2) * t173 + (-pkin(4) * t171 - pkin(7)) * t176) + (-g(1) * t166 + g(2) * t168) * pkin(6)) * t167;];
U_reg = t1;
