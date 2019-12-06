% Calculate minimal parameter regressor of potential energy for
% S5PRRRP8
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
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:35
% EndTime: 2019-12-05 17:00:35
% DurationCPUTime: 0.09s
% Computational Cost: add. (135->54), mult. (327->85), div. (0->0), fcn. (416->10), ass. (0->37)
t181 = sin(pkin(5));
t198 = pkin(6) * t181;
t185 = sin(qJ(3));
t197 = t181 * t185;
t186 = sin(qJ(2));
t196 = t181 * t186;
t188 = cos(qJ(3));
t195 = t181 * t188;
t189 = cos(qJ(2));
t194 = t181 * t189;
t183 = cos(pkin(5));
t193 = t183 * t186;
t192 = t183 * t189;
t180 = sin(pkin(9));
t182 = cos(pkin(9));
t172 = t180 * t189 + t182 * t193;
t164 = t172 * t188 - t182 * t197;
t171 = t180 * t186 - t182 * t192;
t184 = sin(qJ(4));
t187 = cos(qJ(4));
t159 = t164 * t184 - t171 * t187;
t174 = -t180 * t193 + t182 * t189;
t166 = t174 * t188 + t180 * t197;
t173 = t180 * t192 + t182 * t186;
t161 = t166 * t184 - t173 * t187;
t176 = t183 * t185 + t186 * t195;
t167 = t176 * t184 + t187 * t194;
t191 = g(1) * t161 + g(2) * t159 + g(3) * t167;
t163 = t172 * t185 + t182 * t195;
t165 = t174 * t185 - t180 * t195;
t175 = -t183 * t188 + t185 * t196;
t190 = g(1) * t165 + g(2) * t163 + g(3) * t175;
t168 = t176 * t187 - t184 * t194;
t162 = t166 * t187 + t173 * t184;
t160 = t164 * t187 + t171 * t184;
t158 = -g(1) * t162 - g(2) * t160 - g(3) * t168;
t1 = [-g(3) * qJ(1), 0, -g(1) * t174 - g(2) * t172 - g(3) * t196, g(1) * t173 + g(2) * t171 - g(3) * t194, 0, 0, 0, 0, 0, -g(1) * t166 - g(2) * t164 - g(3) * t176, t190, 0, 0, 0, 0, 0, t158, t191, t158, -t190, -t191, -g(1) * (t182 * pkin(1) + t174 * pkin(2) + t166 * pkin(3) + t162 * pkin(4) + t173 * pkin(7) + t165 * pkin(8) + t161 * qJ(5) + t180 * t198) - g(2) * (t180 * pkin(1) + t172 * pkin(2) + t164 * pkin(3) + t160 * pkin(4) + t171 * pkin(7) + t163 * pkin(8) + t159 * qJ(5) - t182 * t198) - g(3) * (t176 * pkin(3) + t168 * pkin(4) + t183 * pkin(6) + t175 * pkin(8) + t167 * qJ(5) + qJ(1) + (pkin(2) * t186 - pkin(7) * t189) * t181);];
U_reg = t1;
