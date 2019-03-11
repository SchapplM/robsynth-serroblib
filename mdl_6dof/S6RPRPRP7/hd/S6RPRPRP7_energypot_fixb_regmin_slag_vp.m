% Calculate minimal parameter regressor of potential energy for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:22:45
% EndTime: 2019-03-09 03:22:45
% DurationCPUTime: 0.08s
% Computational Cost: add. (76->44), mult. (97->56), div. (0->0), fcn. (92->8), ass. (0->28)
t174 = sin(qJ(3));
t190 = pkin(3) * t174;
t170 = qJ(3) + pkin(9);
t165 = cos(t170);
t189 = g(3) * t165;
t173 = sin(qJ(5));
t175 = sin(qJ(1));
t188 = t175 * t173;
t176 = cos(qJ(5));
t187 = t175 * t176;
t178 = cos(qJ(1));
t186 = t178 * t173;
t185 = t178 * t176;
t184 = t178 * pkin(1) + t175 * qJ(2);
t177 = cos(qJ(3));
t183 = t177 * pkin(3) + pkin(2) + pkin(6);
t182 = t175 * t190 + t184;
t181 = -qJ(2) - t190;
t172 = -qJ(4) - pkin(7);
t180 = pkin(5) * t173 - t172;
t160 = g(1) * t175 - g(2) * t178;
t163 = t176 * pkin(5) + pkin(4);
t164 = sin(t170);
t171 = -qJ(6) - pkin(8);
t179 = t163 * t164 + t165 * t171;
t167 = t175 * pkin(1);
t161 = g(1) * t178 + g(2) * t175;
t1 = [0, -t161, t160, t161, -t160, -g(1) * t184 - g(2) * (-t178 * qJ(2) + t167) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t177 - t160 * t174, g(3) * t174 - t160 * t177, -t161, -g(1) * (-t178 * t172 + t182) - g(2) * (-t175 * t172 + t181 * t178 + t167) - g(3) * t183, 0, 0, 0, 0, 0, -g(1) * (t164 * t187 + t186) - g(2) * (-t164 * t185 + t188) - t176 * t189, -g(1) * (-t164 * t188 + t185) - g(2) * (t164 * t186 + t187) + t173 * t189, -g(3) * t164 + t160 * t165, -g(1) * t182 - g(2) * t167 - g(3) * (t165 * t163 - t164 * t171 + t183) + (-g(1) * t179 - g(2) * t180) * t175 + (-g(1) * t180 - g(2) * (-t179 + t181)) * t178;];
U_reg  = t1;
