% Calculate minimal parameter regressor of potential energy for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:44:57
% EndTime: 2019-03-09 01:44:58
% DurationCPUTime: 0.06s
% Computational Cost: add. (78->35), mult. (70->46), div. (0->0), fcn. (65->10), ass. (0->26)
t177 = sin(qJ(4));
t192 = pkin(4) * t177;
t172 = qJ(4) + pkin(10);
t191 = g(3) * cos(t172);
t175 = qJ(2) + pkin(6);
t190 = g(3) * t175;
t173 = qJ(1) + pkin(9);
t167 = sin(t173);
t176 = sin(qJ(6));
t189 = t167 * t176;
t179 = cos(qJ(6));
t188 = t167 * t179;
t169 = cos(t173);
t187 = t169 * t176;
t186 = t169 * t179;
t178 = sin(qJ(1));
t185 = t178 * pkin(1) + t167 * pkin(2);
t181 = cos(qJ(1));
t184 = t181 * pkin(1) + t169 * pkin(2) + t167 * qJ(3);
t183 = -g(1) * t167 + g(2) * t169;
t182 = -g(1) * t181 - g(2) * t178;
t180 = cos(qJ(4));
t174 = -qJ(5) - pkin(7);
t166 = sin(t172);
t162 = g(1) * t169 + g(2) * t167;
t1 = [0, t182, g(1) * t178 - g(2) * t181, t182 * pkin(1) - t190, t162, t183, -g(1) * t184 - g(2) * (-t169 * qJ(3) + t185) - t190, 0, 0, 0, 0, 0, -g(3) * t180 + t183 * t177, g(3) * t177 + t183 * t180, -t162, -g(1) * (t167 * t192 - t169 * t174 + t184) - g(2) * (-t167 * t174 + (-qJ(3) - t192) * t169 + t185) - g(3) * (t180 * pkin(4) + pkin(3) + t175) 0, 0, 0, 0, 0, -g(1) * (t166 * t188 + t187) - g(2) * (-t166 * t186 + t189) - t179 * t191, -g(1) * (-t166 * t189 + t186) - g(2) * (t166 * t187 + t188) + t176 * t191;];
U_reg  = t1;
