% Calculate minimal parameter regressor of potential energy for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:47
% EndTime: 2019-03-09 02:33:47
% DurationCPUTime: 0.06s
% Computational Cost: add. (66->30), mult. (76->42), div. (0->0), fcn. (74->10), ass. (0->22)
t177 = pkin(10) + qJ(4);
t173 = qJ(5) + t177;
t170 = cos(t173);
t190 = g(3) * t170;
t180 = sin(qJ(6));
t181 = sin(qJ(1));
t189 = t181 * t180;
t182 = cos(qJ(6));
t188 = t181 * t182;
t183 = cos(qJ(1));
t187 = t183 * t180;
t186 = t183 * t182;
t185 = t183 * pkin(1) + t181 * qJ(2);
t184 = t181 * pkin(1) - t183 * qJ(2);
t167 = g(1) * t181 - g(2) * t183;
t179 = cos(pkin(10));
t178 = sin(pkin(10));
t172 = cos(t177);
t171 = sin(t177);
t169 = sin(t173);
t168 = g(1) * t183 + g(2) * t181;
t1 = [0, -t168, t167, t168, -t167, -g(3) * pkin(6) - g(1) * t185 - g(2) * t184, -g(3) * t179 - t167 * t178, g(3) * t178 - t167 * t179, -t168, -g(1) * (t183 * qJ(3) + t185) - g(2) * (t181 * qJ(3) + t184) - g(3) * (pkin(2) + pkin(6)) 0, 0, 0, 0, 0, -g(3) * t172 - t167 * t171, g(3) * t171 - t167 * t172, 0, 0, 0, 0, 0, -t167 * t169 - t190, g(3) * t169 - t167 * t170, 0, 0, 0, 0, 0, -g(1) * (t169 * t188 + t187) - g(2) * (-t169 * t186 + t189) - t182 * t190, -g(1) * (-t169 * t189 + t186) - g(2) * (t169 * t187 + t188) + t180 * t190;];
U_reg  = t1;
