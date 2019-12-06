% Calculate minimal parameter regressor of potential energy for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:31
% EndTime: 2019-12-05 16:13:31
% DurationCPUTime: 0.11s
% Computational Cost: add. (105->48), mult. (233->72), div. (0->0), fcn. (272->8), ass. (0->33)
t174 = cos(qJ(2));
t188 = pkin(2) * t174 + pkin(1);
t168 = sin(pkin(7));
t172 = sin(qJ(2));
t186 = t168 * t172;
t170 = cos(pkin(7));
t185 = t170 * t172;
t171 = sin(qJ(3));
t184 = t171 * t172;
t183 = t171 * t174;
t173 = cos(qJ(3));
t182 = t172 * t173;
t181 = t173 * t174;
t180 = g(1) * t170 + g(2) * t168;
t179 = t172 * pkin(2) + pkin(3) * t182 - t174 * pkin(6) + qJ(4) * t184 + qJ(1);
t149 = t168 * t181 - t170 * t171;
t167 = sin(pkin(8));
t169 = cos(pkin(8));
t142 = t149 * t167 - t169 * t186;
t153 = t168 * t171 + t170 * t181;
t144 = t153 * t167 - t169 * t185;
t150 = t167 * t182 + t174 * t169;
t178 = g(1) * t144 + g(2) * t142 + g(3) * t150;
t152 = -t168 * t173 + t170 * t183;
t177 = t153 * pkin(3) + t168 * pkin(5) + pkin(6) * t185 + t152 * qJ(4) + t188 * t170;
t148 = t168 * t183 + t170 * t173;
t176 = g(1) * t152 + g(2) * t148 + g(3) * t184;
t175 = t149 * pkin(3) - t170 * pkin(5) + pkin(6) * t186 + t148 * qJ(4) + t188 * t168;
t151 = -t174 * t167 + t169 * t182;
t145 = t153 * t169 + t167 * t185;
t143 = t149 * t169 + t167 * t186;
t140 = -g(1) * t145 - g(2) * t143 - g(3) * t151;
t1 = [-g(3) * qJ(1), 0, -g(3) * t172 - t180 * t174, -g(3) * t174 + t180 * t172, 0, 0, 0, 0, 0, -g(1) * t153 - g(2) * t149 - g(3) * t182, t176, t140, t178, -t176, -g(1) * t177 - g(2) * t175 - g(3) * t179, t140, -t176, -t178, -g(1) * (t145 * pkin(4) + t144 * qJ(5) + t177) - g(2) * (t143 * pkin(4) + t142 * qJ(5) + t175) - g(3) * (t151 * pkin(4) + t150 * qJ(5) + t179);];
U_reg = t1;
