% Calculate minimal parameter regressor of potential energy for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:47:58
% EndTime: 2019-12-31 19:47:58
% DurationCPUTime: 0.13s
% Computational Cost: add. (64->43), mult. (112->58), div. (0->0), fcn. (114->8), ass. (0->29)
t172 = sin(qJ(2));
t192 = qJ(3) * t172 + pkin(1);
t174 = cos(qJ(2));
t191 = g(3) * t174;
t169 = pkin(8) + qJ(5);
t163 = sin(t169);
t173 = sin(qJ(1));
t189 = t173 * t163;
t164 = cos(t169);
t188 = t173 * t164;
t170 = sin(pkin(8));
t187 = t173 * t170;
t171 = cos(pkin(8));
t186 = t173 * t171;
t185 = t173 * t174;
t175 = cos(qJ(1));
t184 = t174 * t175;
t183 = t175 * t163;
t182 = t175 * t164;
t181 = t175 * t170;
t180 = t175 * t171;
t179 = pkin(2) * t185 + t173 * t192;
t178 = pkin(2) * t184 + t173 * pkin(6) + t175 * t192;
t177 = t172 * pkin(2) - qJ(3) * t174 + pkin(5);
t176 = g(1) * t175 + g(2) * t173;
t158 = g(1) * t173 - g(2) * t175;
t157 = g(3) * t172 + t174 * t176;
t156 = t172 * t176 - t191;
t1 = [0, -t176, t158, 0, 0, 0, 0, 0, -t157, t156, -t158, t157, -t156, -g(1) * t178 - g(2) * (-pkin(6) * t175 + t179) - g(3) * t177, -g(1) * (t172 * t181 + t186) - g(2) * (t172 * t187 - t180) + t170 * t191, -g(1) * (t172 * t180 - t187) - g(2) * (t172 * t186 + t181) + t171 * t191, -t157, -g(1) * (pkin(3) * t173 + qJ(4) * t184 + t178) - g(2) * (qJ(4) * t185 + (-pkin(3) - pkin(6)) * t175 + t179) - g(3) * (qJ(4) * t172 + t177), 0, 0, 0, 0, 0, -g(1) * (t172 * t183 + t188) - g(2) * (t172 * t189 - t182) + t163 * t191, -g(1) * (t172 * t182 - t189) - g(2) * (t172 * t188 + t183) + t164 * t191;];
U_reg = t1;
