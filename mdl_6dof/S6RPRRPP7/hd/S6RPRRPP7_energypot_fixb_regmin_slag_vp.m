% Calculate minimal parameter regressor of potential energy for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:52:11
% EndTime: 2019-03-09 04:52:11
% DurationCPUTime: 0.08s
% Computational Cost: add. (94->51), mult. (188->62), div. (0->0), fcn. (201->6), ass. (0->32)
t172 = sin(qJ(3));
t192 = pkin(3) * t172;
t176 = cos(qJ(1));
t191 = g(2) * t176;
t171 = sin(qJ(4));
t175 = cos(qJ(3));
t190 = t171 * t175;
t173 = sin(qJ(1));
t189 = t173 * t171;
t174 = cos(qJ(4));
t188 = t173 * t174;
t187 = t173 * t175;
t186 = t174 * t175;
t185 = t175 * t176;
t184 = t176 * t171;
t183 = t176 * t174;
t182 = t176 * pkin(1) + t173 * qJ(2);
t181 = -qJ(2) - t192;
t155 = g(1) * t173 - t191;
t180 = t175 * pkin(3) + pkin(4) * t186 + t172 * pkin(8) + qJ(5) * t190 + pkin(2) + pkin(6);
t153 = t172 * t184 + t188;
t154 = -t172 * t183 + t189;
t166 = t173 * pkin(1);
t179 = t154 * pkin(4) + t173 * pkin(7) + pkin(8) * t185 - t153 * qJ(5) + t166;
t151 = t172 * t189 - t183;
t152 = t172 * t188 + t184;
t178 = t152 * pkin(4) + t176 * pkin(7) + t151 * qJ(5) + t173 * t192 + t182;
t177 = g(1) * t151 - g(2) * t153 + g(3) * t190;
t156 = g(1) * t176 + g(2) * t173;
t148 = g(1) * t187 - g(2) * t185 - g(3) * t172;
t147 = -g(1) * t152 - g(2) * t154 - g(3) * t186;
t1 = [0, -t156, t155, t156, -t155, -g(1) * t182 - g(2) * (-t176 * qJ(2) + t166) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t175 - t155 * t172, -t148, 0, 0, 0, 0, 0, t147, t177, t147, t148, -t177, -g(1) * (-pkin(8) * t187 + t178) - g(2) * (t181 * t176 + t179) - g(3) * t180, t147, -t177, -t148, -g(1) * (t152 * pkin(5) + (-pkin(8) + qJ(6)) * t187 + t178) - g(2) * (t154 * pkin(5) + t179) - g(3) * (pkin(5) * t186 - t172 * qJ(6) + t180) - (-qJ(6) * t175 + t181) * t191;];
U_reg  = t1;
